// lqcore -- lqmap (see lqmap.hpp)
#include "lqmap.hpp"

#include "sam.hpp"
#include "scan.hpp"

#include <condition_variable>
#include <deque>
#include <mutex>
#include <sys/stat.h>
#include <thread>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace lq {

namespace {
struct CodeT {
    uint8_t t[256];
    CodeT() {
        for (int i = 0; i < 256; ++i) t[i] = 4;
        t['A'] = t['a'] = 0; t['C'] = t['c'] = 1; t['G'] = t['g'] = 2; t['T'] = t['t'] = 3;
    }
};
const CodeT CT;
const char DEC[5] = {'A', 'C', 'G', 'T', 'N'};

inline void put_uint(std::string& o, unsigned long long v) {
    char b[24];
    int n = 0;
    do { b[n++] = (char)('0' + v % 10); v /= 10; } while (v);
    while (n) o.push_back(b[--n]);
}
inline void put_int(std::string& o, long long v) {
    if (v < 0) { o.push_back('-'); v = -v; }
    put_uint(o, (unsigned long long)v);
}
}  // namespace

// ================================================================== index
void MapIndex::build(const std::string& fasta, int kk, int threads) {
    const double t0 = wall();
    k = kk;
    names.clear(); lens.clear(); tstart.assign(1, 0); seq.clear();
    {
        LineReader rd(fasta);
        std::string_view l;
        bool in = false;
        while (rd.next(l)) {
            if (!l.empty() && l[0] == '>') {
                if (in) { lens.push_back((long long)(seq.size() - tstart.back())); tstart.push_back(seq.size()); }
                auto w = split_ws(l.substr(1));
                names.emplace_back(w.empty() ? std::string() : std::string(w[0]));
                in = true;
                continue;
            }
            if (!in) continue;
            for (char c : l)
                if (!isspace((unsigned char)c)) seq.push_back(CT.t[(uint8_t)c]);
        }
        if (in) { lens.push_back((long long)(seq.size() - tstart.back())); tstart.push_back(seq.size()); }
    }
    if (names.empty()) die("no sequences in " + fasta);
    if (seq.size() >= 0xFFFFFFF0ull) die("transcriptome too large for lqmap (32-bit positions)");
    const size_t R = names.size();
    if (k == 0) {   // sequences only (sa)
        logm(sfmt("[INFO] transcripts: %s (%s bp), %.0fs", commas((long long)R).c_str(), commas((long long)seq.size()).c_str(), wall() - t0));
        return;
    }
    if (k < 12 || k > 29) die("--map-k must be 12..29");
    const uint64_t kmask = (1ULL << (2 * k)) - 1;
    b = 10;
    while (b < 28 && b < 2 * k && (1ULL << b) < seq.size() / 8) ++b;
    lowbits = 2 * k - b;
    if (lowbits > 32) { lowbits = 32; b = 2 * k - 32; }
    const size_t NB = (size_t)1 << b;
    std::vector<std::atomic<uint32_t>> cnt(NB);
    for (auto& x : cnt) x.store(0, std::memory_order_relaxed);
    auto windows = [&](size_t r, auto&& f) {
        uint64_t code = 0;
        int run = 0;
        for (uint64_t i = tstart[r]; i < tstart[r + 1]; ++i) {
            uint8_t c = seq[i];
            if (c < 4) { code = ((code << 2) | c) & kmask; ++run; }
            else run = 0;
            if (run >= k) f(code, (uint32_t)(i + 1 - k));
        }
    };
#pragma omp parallel for schedule(dynamic, 256) num_threads(threads)
    for (int64_t r = 0; r < (int64_t)R; ++r)
        windows((size_t)r, [&](uint64_t code, uint32_t) { cnt[code >> lowbits].fetch_add(1, std::memory_order_relaxed); });
    off.assign(NB + 1, 0);
    uint64_t tot = 0;
    for (size_t i = 0; i < NB; ++i) { off[i] = (uint32_t)tot; tot += cnt[i].load(std::memory_order_relaxed); }
    if (tot >= 0xFFFFFFF0ull) die("too many k-mer positions for lqmap");
    off[NB] = (uint32_t)tot;
    for (size_t i = 0; i < NB; ++i) cnt[i].store(off[i], std::memory_order_relaxed);
    low.resize(tot);
    pos.resize(tot);
    const uint64_t lmask = lowbits == 32 ? 0xFFFFFFFFull : ((1ULL << lowbits) - 1);
#pragma omp parallel for schedule(dynamic, 256) num_threads(threads)
    for (int64_t r = 0; r < (int64_t)R; ++r)
        windows((size_t)r, [&](uint64_t code, uint32_t p) {
            uint32_t s = cnt[code >> lowbits].fetch_add(1, std::memory_order_relaxed);
            low[s] = (uint32_t)(code & lmask);
            pos[s] = p;
        });
    std::vector<std::atomic<uint32_t>>().swap(cnt);
    // sort each bucket by (low, pos): deterministic whatever the fill order
#pragma omp parallel num_threads(threads)
    {
        std::vector<uint64_t> tmp;
#pragma omp for schedule(dynamic, 4096)
        for (int64_t bk = 0; bk < (int64_t)NB; ++bk) {
            const uint32_t a = off[bk], z = off[bk + 1];
            if (z - a < 2) continue;
            tmp.clear();
            for (uint32_t i = a; i < z; ++i) tmp.push_back(((uint64_t)low[i] << 32) | pos[i]);
            std::sort(tmp.begin(), tmp.end());
            for (uint32_t i = a; i < z; ++i) { low[i] = (uint32_t)(tmp[i - a] >> 32); pos[i] = (uint32_t)tmp[i - a]; }
        }
    }
    logm(sfmt("[INFO] lqmap index: %s transcripts, %s bp, %s %d-mer positions, %s buckets (%.1f GB), %.0fs",
              commas((long long)R).c_str(), commas((long long)seq.size()).c_str(), commas((long long)tot).c_str(), k,
              commas((long long)NB).c_str(), (seq.size() + 8.0 * tot + 4.0 * NB) / 1e9, wall() - t0));
}

namespace {
const char MX_MAGIC[8] = {'L', 'Q', 'M', 'A', 'P', 'I', 'X', '1'};
}

void MapIndex::save(const std::string& path, const std::string& sig) const {
    std::string tmp = path + sfmt(".tmp%d", (int)getpid());
    FILE* f = fopen(tmp.c_str(), "wb");
    if (!f) { logm("[WARN] cannot write the lqmap index cache " + path); return; }
    std::vector<char> buf(1 << 24);
    setvbuf(f, buf.data(), _IOFBF, buf.size());
    auto wr = [&](const void* p, size_t n) { return n == 0 || fwrite(p, 1, n, f) == n; };
    auto wv = [&](const auto& v) { uint64_t n = v.size(); return wr(&n, 8) && wr(v.data(), n * sizeof(v[0])); };
    uint64_t sl = sig.size(), nn = names.size();
    int32_t hd[3] = {k, lowbits, b};
    bool ok = wr(MX_MAGIC, 8) && wr(&sl, 8) && wr(sig.data(), sl) && wr(hd, sizeof hd) && wr(&nn, 8);
    for (auto& n : names) { uint64_t l = n.size(); ok = ok && wr(&l, 8) && wr(n.data(), l); }
    ok = ok && wv(lens) && wv(tstart) && wv(seq) && wv(off) && wv(low) && wv(pos);
    ok = (fclose(f) == 0) && ok;
    if (ok && rename(tmp.c_str(), path.c_str()) == 0) logm("[INFO] lqmap index cached: " + path);
    else { unlink(tmp.c_str()); logm("[WARN] cannot write the lqmap index cache " + path); }
}

bool MapIndex::load(const std::string& path, const std::string& sig) {
    FILE* f = fopen(path.c_str(), "rb");
    if (!f) return false;
    std::vector<char> buf(1 << 24);
    setvbuf(f, buf.data(), _IOFBF, buf.size());
    auto rd = [&](void* p, size_t n) { return n == 0 || fread(p, 1, n, f) == n; };
    auto rv = [&](auto& v) {
        uint64_t n = 0;
        if (!rd(&n, 8) || n > (1ULL << 36)) return false;
        v.resize(n);
        return rd(v.data(), n * sizeof(v[0]));
    };
    char m[8];
    uint64_t sl = 0, nn = 0;
    int32_t hd[3];
    bool ok = rd(m, 8) && memcmp(m, MX_MAGIC, 8) == 0 && rd(&sl, 8) && sl < (1u << 20);
    std::string s(ok ? sl : 0, '\0');
    ok = ok && rd(&s[0], sl) && s == sig && rd(hd, sizeof hd) && hd[0] == k && rd(&nn, 8) && nn < (1u << 30);
    if (ok) {
        lowbits = hd[1]; b = hd[2];
        names.resize(nn);
        for (auto& n : names) {
            uint64_t l = 0;
            ok = ok && rd(&l, 8) && l < 65536;
            if (!ok) break;
            n.resize(l);
            ok = rd(&n[0], l);
        }
    }
    ok = ok && rv(lens) && rv(tstart) && rv(seq) && rv(off) && rv(low) && rv(pos);
    fclose(f);
    if (!ok) { names.clear(); lens.clear(); tstart.clear(); seq.clear(); off.clear(); low.clear(); pos.clear(); }
    return ok;
}

// ================================================================== mapper
void Mapper::encode(std::string_view s, int m) {
    const size_t L = s.size();
    auto& f = q_[m][0];
    auto& r = q_[m][1];
    f.resize(L);
    r.resize(L);
    for (size_t i = 0; i < L; ++i) {
        uint8_t c = CT.t[(uint8_t)s[i]];
        f[i] = c;
        r[L - 1 - i] = c < 4 ? (uint8_t)(3 - c) : (uint8_t)4;
    }
}

// mismatches: mm = both ACGT and different (scored), nm = any difference or ambiguous base (the NM tag)
bool Mapper::verify(const uint8_t* q, int L, uint64_t g, int maxmm, int& mm, int& nm) const {
    const uint8_t* r = X_.seq.data() + g;
    int a = 0, n = 0, amb = 0;
    const int nceil = (int)(o_.n_ceil_frac * L);
    int i = 0;
    for (; i + 16 <= L; i += 16) {
        for (int j = 0; j < 16; ++j) {
            const uint8_t x = q[i + j], y = r[i + j];
            const int ab = (x == 4) | (y == 4);
            amb += ab;
            n += (x != y) | ab;
            a += (x != y) & !ab;
        }
        if (a > maxmm || amb > nceil) return false;
    }
    for (; i < L; ++i) {
        const uint8_t x = q[i], y = r[i];
        const int ab = (x == 4) | (y == 4);
        amb += ab;
        n += (x != y) | ab;
        a += (x != y) & !ab;
    }
    if (a > maxmm || amb > nceil) return false;
    mm = a;
    nm = n;
    return true;
}

void Mapper::seeds(int m, int L, int maxmm, bool mm1) {
    auto& H = h_[m];
    if (!mm1) H.clear();
    const size_t had = H.size();
    const int K = X_.k;
    if (L < K) return;
    for (int pass = 0; pass < 2; ++pass) {
        const uint32_t cap = pass == 0 ? (uint32_t)o_.max_occ : (uint32_t)o_.max_occ * 25u;
        diag_.clear();
        int skipped = 0;
        for (int sd = 0; sd < 2; ++sd) {
            const uint8_t* q = q_[m][sd].data();
            int last = -1;
            for (int o = 0;; o += o_.stride) {
                if (o > L - K) o = L - K;
                if (o == last) break;
                last = o;
                uint64_t code = 0;
                bool ok = true;
                for (int j = 0; j < K; ++j) {
                    if (q[o + j] > 3) { ok = false; break; }
                    code = (code << 2) | q[o + j];
                }
                if (ok) {
                    auto add = [&](uint64_t c) {
                        uint32_t lo, hi;
                        X_.find(c, lo, hi);
                        if (hi - lo > cap) { ++skipped; return; }
                        for (uint32_t i = lo; i < hi; ++i) {
                            const uint32_t g = X_.pos[i];
                            if (g >= (uint32_t)o) diag_.push_back(((uint64_t)(g - o) << 1) | (uint64_t)sd);
                        }
                    };
                    if (!mm1) add(code);
                    else
                        for (int j = 0; j < K; ++j) {          // every seed with exactly one substitution
                            const int sh = 2 * (K - 1 - j);
                            const uint64_t cur = (code >> sh) & 3;
                            for (uint64_t v = 0; v < 4; ++v)
                                if (v != cur) add((code & ~(3ULL << sh)) | (v << sh));
                        }
                }
                if (o == L - K) break;
            }
        }
        if (!diag_.empty() || skipped == 0) break;
    }
    std::sort(diag_.begin(), diag_.end());
    diag_.erase(std::unique(diag_.begin(), diag_.end()), diag_.end());
    for (uint64_t d : diag_) {
        const int sd = (int)(d & 1);
        const uint64_t g0 = d >> 1;
        const uint32_t t = X_.tid_of(g0);
        if (g0 + (uint64_t)L > X_.tstart[t + 1]) continue;       // runs off the transcript end
        int mm, nm;
        if (!verify(q_[m][sd].data(), L, g0, maxmm, mm, nm)) continue;
        H.push_back({t, (int32_t)(g0 - X_.tstart[t]), (uint8_t)sd, (int16_t)mm, (int16_t)nm});
    }
    if (mm1 && H.size() > had) {
        std::sort(H.begin(), H.end(), [](const Hit& a, const Hit& b) {
            return a.t != b.t ? a.t < b.t : a.st != b.st ? a.st < b.st : a.rev < b.rev;
        });
        H.erase(std::unique(H.begin(), H.end(), [](const Hit& a, const Hit& b) { return a.t == b.t && a.st == b.st && a.rev == b.rev; }),
                H.end());
    }
}

// scan the fragment window for the other mate (ungapped), next to anchor a
void Mapper::rescue(int have, const Hit& a, int Lm, int maxmm, std::vector<Hit>& add) {
    const int other = 1 - have;
    const int La = (int)q_[have][0].size();
    const long long tl = X_.lens[a.t];
    const int sd = a.rev ? 0 : 1;           // the other mate lies on the opposite strand
    long long s0, s1;
    if (!a.rev) { s0 = a.st; s1 = (long long)a.st + o_.max_frag - Lm; }
    else { s0 = (long long)a.st + La - o_.max_frag; s1 = a.st; }
    s0 = std::max(0LL, s0);
    s1 = std::min(s1, tl - Lm);
    if (s1 < s0) return;
    const uint8_t* q = q_[other][sd].data();
    int best = 1 << 30;
    size_t first = add.size();
    for (long long s = s0; s <= s1; ++s) {
        int mm, nm;
        if (!verify(q, Lm, X_.tstart[a.t] + (uint64_t)s, maxmm, mm, nm)) continue;
        if (nm > best + o_.max_extra) continue;
        best = std::min(best, nm);
        add.push_back({a.t, (int32_t)s, (uint8_t)sd, (int16_t)mm, (int16_t)nm});
    }
    // keep only those within best + max_extra
    size_t w = first;
    for (size_t i = first; i < add.size(); ++i) if (add[i].nm <= best + o_.max_extra) add[w++] = add[i];
    add.resize(w);
}

void Mapper::map_pair(std::string_view name, std::string_view r1, std::string_view r2, std::string& out, std::vector<uint32_t>& ls) {
    encode(r1, 0);
    encode(r2, 1);
    const int L[2] = {(int)r1.size(), (int)r2.size()};
    const int maxmm[2] = {(int)(o_.mm_frac * L[0] + 1e-9), (int)(o_.mm_frac * L[1] + 1e-9)};
    for (int m = 0; m < 2; ++m) seeds(m, L[m], maxmm[m]);
    finish(name, L, maxmm, out, ls, o_.rescue);
}

void Mapper::map_hits(std::string_view name, std::string_view r1, std::string_view r2,
                      const std::vector<std::array<int64_t, 3>>* given, std::string& out, std::vector<uint32_t>& ls) {
    encode(r1, 0);
    encode(r2, 1);
    const int L[2] = {(int)r1.size(), (int)r2.size()};
    const int maxmm[2] = {(int)(o_.mm_frac * L[0] + 1e-9), (int)(o_.mm_frac * L[1] + 1e-9)};
    for (int m = 0; m < 2; ++m) {
        auto& H = h_[m];
        H.clear();
        for (auto& x : given[m]) {          // (transcript, 0-based start, reverse)
            const uint32_t t = (uint32_t)x[0];
            const int64_t s = x[1];
            if (s < 0 || s + L[m] > X_.lens[t]) continue;
            int mm, nm;
            if (!verify(q_[m][x[2] ? 1 : 0].data(), L[m], X_.tstart[t] + (uint64_t)s, maxmm[m], mm, nm)) continue;
            H.push_back({t, (int32_t)s, (uint8_t)(x[2] ? 1 : 0), (int16_t)mm, (int16_t)nm});
        }
        std::sort(H.begin(), H.end(), [](const Hit& a, const Hit& b) {
            return a.t != b.t ? a.t < b.t : a.st != b.st ? a.st < b.st : a.rev < b.rev;
        });
        H.erase(std::unique(H.begin(), H.end(), [](const Hit& a, const Hit& b) { return a.t == b.t && a.st == b.st && a.rev == b.rev; }),
                H.end());
    }
    finish(name, L, maxmm, out, ls, false);
}

void Mapper::finish(std::string_view name, const int* L, const int* maxmm, std::string& out, std::vector<uint32_t>& ls, bool resc) {
    st_.pairs++;
    auto bytid = [](const Hit& a, const Hit& b) { return a.t != b.t ? a.t < b.t : a.st < b.st; };
    std::sort(h_[0].begin(), h_[0].end(), bytid);
    std::sort(h_[1].begin(), h_[1].end(), bytid);
    std::vector<PC>& pc = pc_;
    auto pairup = [&]() {
        pc.clear();
        const auto &A = h_[0], &B = h_[1];
        size_t j0 = 0;
        for (size_t i = 0; i < A.size();) {
            const uint32_t t = A[i].t;
            size_t i2 = i;
            while (i2 < A.size() && A[i2].t == t) ++i2;
            while (j0 < B.size() && B[j0].t < t) ++j0;
            size_t j2 = j0;
            while (j2 < B.size() && B[j2].t == t) ++j2;
            PC best{t, 0, 0, 1 << 30, 0};
            for (size_t x = i; x < i2; ++x)
                for (size_t y = j0; y < j2; ++y) {
                    const Hit &a = A[x], &bb = B[y];
                    if (a.rev == bb.rev) continue;
                    const Hit& f = a.rev ? bb : a;
                    const Hit& r = a.rev ? a : bb;
                    const int Lf = a.rev ? L[1] : L[0], Lr = a.rev ? L[0] : L[1];
                    if (f.st > r.st) continue;                                   // dovetail
                    const int fr = std::max(f.st + Lf, r.st + Lr) - f.st;
                    if (fr < o_.min_frag || fr > o_.max_frag) continue;
                    const int nm = a.nm + bb.nm;
                    if (nm < best.nm || (nm == best.nm && (fr < best.frag))) best = {t, (uint32_t)x, (uint32_t)y, nm, fr};
                }
            if (best.nm < (1 << 30)) pc.push_back(best);
            i = i2;
        }
    };
    pairup();
    int best = 1 << 30;
    for (auto& p : pc) best = std::min(best, p.nm);
    // divergent fragment (no pair, or best pair at least half the mismatch budget): seeds with one substitution
    if (resc && (pc.empty() || 2 * best >= maxmm[0] + maxmm[1])) {
        st_.reseeded++;
        for (int m = 0; m < 2; ++m) seeds(m, L[m], maxmm[m], true);
        std::sort(h_[0].begin(), h_[0].end(), bytid);
        std::sort(h_[1].begin(), h_[1].end(), bytid);
        pairup();
        best = 1 << 30;
        for (auto& p : pc) best = std::min(best, p.nm);
    }
    if (resc) {
        // mate rescue: a transcript where one mate aligns but no concordant pair was found can still hold a candidate
        // if that mate alone is within best pair + max_extra; scan its fragment window for the other mate
        const int lim = pc.empty() ? (1 << 30) : best + o_.max_extra;
        struct An { int nm, m; size_t i; };
        std::vector<An> an;
        for (int m = 0; m < 2; ++m)
            for (size_t i = 0; i < h_[m].size(); ++i) {
                const Hit& h = h_[m][i];
                if (h.nm > lim) continue;
                bool paired = false;
                for (auto& p : pc) if (p.t == h.t) { paired = true; break; }
                if (!paired) an.push_back({h.nm, m, i});
            }
        if (!an.empty()) {
            std::sort(an.begin(), an.end(), [](const An& x, const An& y) { return x.nm != y.nm ? x.nm < y.nm : x.m != y.m ? x.m < y.m : x.i < y.i; });
            if (an.size() > 64) an.resize(64);
            std::vector<Hit> add[2];
            for (auto& x : an) {
                const Hit h = h_[x.m][x.i];
                const int room = lim == (1 << 30) ? maxmm[1 - x.m] : std::min(maxmm[1 - x.m], lim);
                (void)room;
                rescue(x.m, h, L[1 - x.m], maxmm[1 - x.m], add[1 - x.m]);
            }
            for (int m = 0; m < 2; ++m) {
                h_[m].insert(h_[m].end(), add[m].begin(), add[m].end());
                std::sort(h_[m].begin(), h_[m].end(), [](const Hit& a, const Hit& b) {
                    return a.t != b.t ? a.t < b.t : a.st != b.st ? a.st < b.st : a.rev < b.rev;
                });
                h_[m].erase(std::unique(h_[m].begin(), h_[m].end(), [](const Hit& a, const Hit& b) {
                    return a.t == b.t && a.st == b.st && a.rev == b.rev; }), h_[m].end());
            }
            size_t before = pc.size();
            pairup();
            if (pc.size() > before) st_.rescued++;
            best = 1 << 30;
            for (auto& p : pc) best = std::min(best, p.nm);
        }
    }
    // keep candidates within best + max_extra, best first (then transcript order)
    size_t w = 0;
    for (auto& p : pc) if (p.nm <= best + o_.max_extra) pc[w++] = p;
    pc.resize(w);
    std::sort(pc.begin(), pc.end(), [](const PC& a, const PC& b) { return a.nm != b.nm ? a.nm < b.nm : a.t < b.t; });
    if ((int)pc.size() > o_.max_cand) pc.resize(o_.max_cand);
    auto seqstr = [&](int m, int sd) {
        const auto& q = q_[m][sd];
        for (uint8_t c : q) out.push_back(DEC[c]);
    };
    if (pc.empty()) {
        st_.conc0++;
        for (int m = 0; m < 2; ++m) {
            ls.push_back((uint32_t)out.size());
            out.append(name);
            out.append(m == 0 ? "\t77\t*\t0\t0\t*\t*\t0\t0\t" : "\t141\t*\t0\t0\t*\t*\t0\t0\t");
            seqstr(m, 0);
            out.append("\t*\tYT:Z:UP\n");
        }
        return;
    }
    (pc.size() == 1 ? st_.conc1 : st_.concm)++;
    st_.cands += (long long)pc.size();
    for (size_t j = 0; j < pc.size(); ++j) {
        const Hit* hh[2] = {&h_[0][pc[j].i0], &h_[1][pc[j].i1]};
        const int left = std::min(hh[0]->st, hh[1]->st);
        for (int m = 0; m < 2; ++m) {
            const Hit& h = *hh[m];
            const Hit& o = *hh[1 - m];
            int flag = 1 | 2 | (h.rev ? 16 : 0) | (o.rev ? 32 : 0) | (m == 0 ? 64 : 128) | (j ? 256 : 0);
            long long tlen = pc[j].frag;
            if (!(h.st == left && (o.st != left || m == 0))) tlen = -tlen;
            ls.push_back((uint32_t)out.size());
            out.append(name);
            out.push_back('\t'); put_uint(out, (unsigned)flag);
            out.push_back('\t'); out.append(X_.names[h.t]);
            out.push_back('\t'); put_uint(out, (unsigned long long)h.st + 1);
            out.append("\t255\t"); put_uint(out, (unsigned)L[m]); out.push_back('M');
            out.append("\t=\t"); put_uint(out, (unsigned long long)o.st + 1);
            out.push_back('\t'); put_int(out, tlen);
            out.push_back('\t'); seqstr(m, h.rev);
            out.append("\t*\tNM:i:"); put_uint(out, (unsigned)h.nm);
            out.push_back('\n');
        }
    }
}

// ================================================================== FASTQ pairs -> batches
namespace {
struct RecBlock {
    std::string data;                 // name\0seq\0 per record
    std::vector<uint32_t> off;        // name start, seq start per record (+ end)
    size_t n() const { return off.size() / 2; }
};

class FqStream {
public:
    FqStream(const std::string& path, size_t per) : path_(path), per_(per) {}
    void run() {
        LineReader rd(path_);
        std::string_view l;
        auto b = std::make_unique<RecBlock>();
        long long li = 0;
        while (rd.next(l)) {
            int m = (int)(li++ % 4);
            if (m == 0) {
                auto w = split_ws(l.size() ? l.substr(1) : l);
                std::string_view nm = w.empty() ? std::string_view() : w[0];
                if (nm.size() > 2 && nm[nm.size() - 2] == '/' && (nm.back() == '1' || nm.back() == '2')) nm.remove_suffix(2);
                b->off.push_back((uint32_t)b->data.size());
                b->data.append(nm);
                b->data.push_back('\0');
            } else if (m == 1) {
                b->off.push_back((uint32_t)b->data.size());
                size_t s = 0, e = l.size();
                while (s < e && isspace((unsigned char)l[s])) ++s;
                while (e > s && isspace((unsigned char)l[e - 1])) --e;
                b->data.append(l.substr(s, e - s));
                b->data.push_back('\0');
            } else if (m == 3 && b->n() >= per_) {
                if (!put(std::move(b))) break;
                b = std::make_unique<RecBlock>();
            }
        }
        if (b) {
            if (li % 4 == 1) b->off.pop_back();                                   // drop a truncated last record
            else if (li % 4 == 2 || li % 4 == 3) { b->off.pop_back(); b->off.pop_back(); }
            if (b->n()) put(std::move(b));
        }
        std::lock_guard<std::mutex> g(mu_);
        done_ = true;
        cv_.notify_all();
    }
    std::unique_ptr<RecBlock> take() {
        std::unique_lock<std::mutex> g(mu_);
        cv_.wait(g, [&] { return !q_.empty() || done_; });
        if (q_.empty()) return nullptr;
        auto b = std::move(q_.front());
        q_.pop_front();
        cv_.notify_all();
        return b;
    }
    void stop() { std::lock_guard<std::mutex> g(mu_); stop_ = true; cv_.notify_all(); }

private:
    std::string path_;
    size_t per_;
    std::deque<std::unique_ptr<RecBlock>> q_;
    std::mutex mu_;
    std::condition_variable cv_;
    bool done_ = false, stop_ = false;
    bool put(std::unique_ptr<RecBlock> b) {
        std::unique_lock<std::mutex> g(mu_);
        cv_.wait(g, [&] { return q_.size() < 6 || stop_; });
        if (stop_) return false;
        q_.push_back(std::move(b));
        cv_.notify_all();
        return true;
    }
};
}  // namespace

// reads R1/R2 (two decompression threads) and emits batches of read pairs (raw = name\0r1\0r2\0)
long long produce_pairs(const std::string& r1, const std::string& r2, int chunk_i, int chunk_n, long long block,
                        size_t per_batch, const ScanEmit& emit) {
    const size_t PER = 65536;
    FqStream A(r1, PER), B(r2, PER);
    std::thread ta([&] { A.run(); }), tb([&] { B.run(); });
    std::unique_ptr<RecBlock> a, b;
    size_t ia = 0, ib = 0;
    long long pidx = 0, kept = 0;
    bool mismatch = false;
    auto batch = std::make_unique<ScanBatch>();
    for (;;) {
        if (!a || ia >= a->n()) { a = A.take(); ia = 0; }
        if (!b || ib >= b->n()) { b = B.take(); ib = 0; }
        if (!a || !b) { mismatch = (bool)a != (bool)b; break; }
        const size_t n = std::min(a->n() - ia, b->n() - ib);
        for (size_t i = 0; i < n; ++i, ++ia, ++ib, ++pidx) {
            if (chunk_n > 1 && (pidx / block) % chunk_n != chunk_i - 1) continue;
            const char* na = a->data.data() + a->off[2 * ia];
            const char* sa = a->data.data() + a->off[2 * ia + 1];
            const char* sb = b->data.data() + b->off[2 * ib + 1];
            batch->roff.push_back((uint32_t)batch->raw.size());
            batch->raw.append(na); batch->raw.push_back('\0');
            batch->raw.append(sa); batch->raw.push_back('\0');
            batch->raw.append(sb); batch->raw.push_back('\0');
            ++kept;
            if (batch->roff.size() >= per_batch) {
                batch->roff.push_back((uint32_t)batch->raw.size());
                emit(std::move(batch));
                batch = std::make_unique<ScanBatch>();
            }
        }
    }
    if (!batch->roff.empty()) { batch->roff.push_back((uint32_t)batch->raw.size()); emit(std::move(batch)); }
    A.stop(); B.stop();
    while (A.take()) {}
    while (B.take()) {}
    ta.join(); tb.join();
    if (mismatch) die("R1 and R2 contain different numbers of reads");
    return kept;
}

void map_batch(Mapper& mp, ScanBatch& b) {
    b.data.clear(); b.ls.clear(); b.fs.clear();
    const size_t n = b.roff.empty() ? 0 : b.roff.size() - 1;
    b.data.reserve(n * 1500);
    for (size_t i = 0; i < n; ++i) {
        const char* p = b.raw.data() + b.roff[i];
        std::string_view nm(p);
        std::string_view s1(p + nm.size() + 1);
        std::string_view s2(p + nm.size() + 1 + s1.size() + 1);
        b.fs.push_back((uint32_t)b.ls.size());
        mp.map_pair(nm, s1, s2, b.data, b.ls);
    }
    b.ls.push_back((uint32_t)b.data.size());
    b.fs.push_back((uint32_t)(b.ls.size() - 1));
    std::string().swap(b.raw);
}

void log_map_stats(const MapStats& s, const char* who) {
    const long long n = s.pairs.load(), c0 = s.conc0.load(), c1 = s.conc1.load(), cm = s.concm.load();
    auto pc = [&](long long x) { return n ? 100.0 * (double)x / (double)n : 0.0; };
    logm(sfmt("[INFO] %s: %s pairs; concordant 0 times %s (%.2f%%), exactly 1 time %s (%.2f%%), >1 times %s (%.2f%%); "
              "overall %.2f%%; %.2f candidates per aligned pair; re-seeded %s, rescued %s", who, commas(n).c_str(), commas(c0).c_str(), pc(c0),
              commas(c1).c_str(), pc(c1), commas(cm).c_str(), pc(cm), pc(c1 + cm),
              (c1 + cm) ? (double)s.cands.load() / (double)(c1 + cm) : 0.0, commas(s.reseeded.load()).c_str(), commas(s.rescued.load()).c_str()));
}

MapOpts map_opts_from(const Opts& a) {
    MapOpts m;
    m.stride = std::max(1, a.map_stride);
    m.max_extra = a.map_max_extra;
    m.max_occ = a.map_max_occ;
    m.max_cand = a.map_max_cand;
    m.max_frag = a.map_max_frag;
    m.mm_frac = a.map_mm_frac;
    m.rescue = !a.map_no_rescue;
    return m;
}

void load_map_index(MapIndex& X, const Opts& a, int k, int threads) {
    if (a.fasta.empty()) die("--aligner " + a.aligner + " needs --fasta (the transcriptome)");
    struct stat st;
    if (stat(a.fasta.c_str(), &st) != 0) die("cannot stat " + a.fasta);
    const std::string sig = "k=" + std::to_string(k) + "|" + a.fasta + "|" + std::to_string((long long)st.st_size) + "|" +
                            std::to_string((long long)st.st_mtime);
    X.k = k;
    const double t0 = wall();
    if (k > 0 && !a.map_cache.empty() && X.load(a.map_cache, sig)) {
        logm(sfmt("[INFO] lqmap index loaded from %s: %s transcripts, %s positions (%.0fs)", a.map_cache.c_str(),
                  commas((long long)X.names.size()).c_str(), commas((long long)X.pos.size()).c_str(), wall() - t0));
        return;
    }
    X.build(a.fasta, k, threads);
    if (k > 0 && !a.map_cache.empty()) X.save(a.map_cache, sig);
}

// ================================================================== sa: Salmon selective-alignment records -> LQ candidates
void sa_batch(Mapper& mp, const std::unordered_map<std::string, uint32_t>& name2t, ScanBatch& b, SaStats& st) {
    std::string out;
    std::vector<uint32_t> ls, fs;
    out.reserve(b.data.size());
    std::vector<std::array<int64_t, 3>> given[2];
    std::string seq[2], tmp;
    auto rcs = [](std::string_view s, std::string& o) {
        o.resize(s.size());
        for (size_t i = 0; i < s.size(); ++i) {
            char c = s[s.size() - 1 - i];
            switch (c) {
                case 'A': case 'a': c = 'T'; break;
                case 'C': case 'c': c = 'G'; break;
                case 'G': case 'g': c = 'C'; break;
                case 'T': case 't': c = 'A'; break;
                default: c = 'N';
            }
            o[i] = c;
        }
    };
    const size_t nf = b.fs.empty() ? 0 : b.fs.size() - 1;
    for (size_t f = 0; f < nf; ++f) {
        given[0].clear(); given[1].clear();
        seq[0].clear(); seq[1].clear();
        std::string_view qn;
        for (uint32_t li = b.fs[f]; li < b.fs[f + 1]; ++li) {
            std::string_view l(b.data.data() + b.ls[li], b.ls[li + 1] - b.ls[li] - 1);
            Rec r;
            if (!parse_sam(l, r)) die("malformed SAM line from salmon: " + std::string(l.substr(0, 200)));
            st.recs++;
            qn = r.qname;
            const int m = (r.flag & 0x80) ? 1 : 0;
            const bool rev = (r.flag & 16) != 0;
            if (r.seq != "*" && seq[m].empty()) {
                if (rev) { rcs(r.seq, tmp); seq[m] = tmp; }
                else seq[m].assign(r.seq);
            }
            if (r.flag & 4) continue;
            auto it = name2t.find(std::string(r.rname));
            if (it == name2t.end()) { st.unknown_ref++; continue; }
            given[m].push_back({(int64_t)it->second, (int64_t)r.pos - 1, rev ? 1 : 0});
        }
        if (qn.size() > 2 && qn[qn.size() - 2] == '/' && (qn.back() == '1' || qn.back() == '2')) qn.remove_suffix(2);
        st.frags++;
        if (seq[0].empty() || seq[1].empty()) { st.no_seq++; continue; }
        fs.push_back((uint32_t)ls.size());
        mp.map_hits(qn, seq[0], seq[1], given, out, ls);
    }
    ls.push_back((uint32_t)out.size());
    fs.push_back((uint32_t)(ls.size() - 1));
    b.data.swap(out);
    b.ls.swap(ls);
    b.fs.swap(fs);
}

// ================================================================== lqcore map: SAM to stdout (inspection, tests)
int map_main(int argc, char** argv) {
    Opts a = parse_opts(argc, argv, [] {
        std::fprintf(stderr, "usage: lqcore map --fasta TX.fa --r1 R1 --r2 R2 [--threads N] [--map-k 21] [--map-stride 8]\n"
                             "        [--map-max-extra 3] [--map-max-occ 4000] [--map-max-frag 1000] [--map-mm-frac 0.1]\n"
                             "        [--map-no-rescue] [--map-cache FILE] > out.sam\n");
    });
    a.aligner = "lqmap";
    if (a.r1.empty() || a.r2.empty()) die("--r1 and --r2 are required");
    // SAM goes to the original stdout, log lines (logm prints to stdout) to stderr
    int samfd = dup(1);
    if (samfd < 0 || dup2(2, 1) < 0) die("cannot set up output");
    FILE* so = fdopen(samfd, "w");
    if (!so) die("cannot open SAM output");
    const int T = std::max(1, a.threads);
#ifdef _OPENMP
    omp_set_num_threads(T);
#endif
    MapIndex X;
    load_map_index(X, a, a.map_k, T);
    MapOpts mo = map_opts_from(a);
    MapStats st;
    std::vector<std::unique_ptr<Mapper>> mps;
    for (int i = 0; i < T; ++i) mps.emplace_back(new Mapper(X, mo, st));
    std::vector<char> obuf(1 << 24);
    setvbuf(so, obuf.data(), _IOFBF, obuf.size());
    std::fprintf(so, "@HD\tVN:1.0\tSO:unsorted\n");
    for (size_t i = 0; i < X.names.size(); ++i) std::fprintf(so, "@SQ\tSN:%s\tLN:%lld\n", X.names[i].c_str(), X.lens[i]);
    std::fprintf(so, "@PG\tID:lqmap\tPN:lqcore map\n");
    std::vector<std::unique_ptr<ScanBatch>> pend;
    auto flush = [&]() {
        const int n = (int)pend.size();
#pragma omp parallel for schedule(dynamic, 1) num_threads(T)
        for (int i = 0; i < n; ++i) {
#ifdef _OPENMP
            map_batch(*mps[omp_get_thread_num()], *pend[i]);
#else
            map_batch(*mps[0], *pend[i]);
#endif
        }
        for (auto& b : pend)
            if (fwrite(b->data.data(), 1, b->data.size(), so) != b->data.size()) die("error writing SAM");
        pend.clear();
    };
    produce_pairs(a.r1, a.r2, a.chunk_i, a.chunk_n, a.block_pairs, 4096, [&](std::unique_ptr<ScanBatch> b) {
        pend.push_back(std::move(b));
        if ((int)pend.size() >= 2 * T) flush();
    });
    flush();
    if (fclose(so) != 0) die("error writing SAM");
    log_map_stats(st, "lqmap");
    return 0;
}

}  // namespace lq
