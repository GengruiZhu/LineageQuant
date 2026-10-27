// lqcore -- lineage evidence per fragment
//   k-mer   : port of lineage_kmer/snipe_kmer.py (v2: Beta smoothing, per-mate sequences, canonical detection)
//   SNP     : port of lineage_snp/snipe_snp.py (msa / population)
//   fusion  : port of merge/merge_poe.py (entropy-weighted product of experts)
// Probabilities are rounded to 4 decimals wherever the Python pipeline writes and re-reads them.
//
// Difference from the Python pipeline (deliberate, deterministic): a fragment aligned to several transcripts
// gets one SNP row. snipe_snp writes one row per (fragment, transcript) and merge_poe keeps whichever came
// last in an order set by multiprocessing; here the row comes from the transcript with the fewest mismatches,
// then the most SNPs covered, then the first in the reference order (--snp-pick best; also: most, first).
#include "lineage.hpp"

#include <atomic>
#include <cmath>
#include <set>
#include <sys/stat.h>
#include <thread>
#include <unordered_set>

namespace lq {

static inline uint64_t fmix64(uint64_t k) {
    k ^= k >> 33;
    k *= 0xff51afd7ed558ccdULL;
    k ^= k >> 33;
    k *= 0xc4ceb9fe1a85ec53ULL;
    k ^= k >> 33;
    return k;
}
static inline int bcode(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return -1;
    }
}
static inline uint64_t revcomp(uint64_t km, int k) {
    uint64_t r = 0;
    for (int i = 0; i < k; ++i) { r = (r << 2) | (3 - (km & 3)); km >>= 2; }
    return r;
}

inline uint8_t MarkerIndex::get(uint64_t km) const {
    const uint64_t x = fmix64(km);
    const uint64_t bk = x >> (64 - b);
    for (uint32_t i = off[bk], e = off[bk + 1]; i < e; ++i)
        if (h[i] == x) return mask[i];
    return 0;
}

// ================================================================== markers
void MarkerIndex::build(const std::vector<std::pair<std::string, std::string>>& np, int kk, int threads) {
    k = kk;
    const size_t L = np.size();
    if (L == 0 || L > 8) die("between 1 and 8 marker lineages are supported");
    if (k < 1 || k > 32) die("--kmer-len must be 1..32");
    names.clear();
    for (auto& x : np) names.push_back(x.first);
    std::vector<std::vector<uint64_t>> hs(L);
    std::vector<char> canon(L, 0);
    const uint64_t kmask = k == 32 ? ~0ULL : ((1ULL << (2 * k)) - 1);
    auto work = [&](size_t l) {
        LineReader rd(np[l].second);
        std::string_view line;
        std::vector<uint64_t>& v = hs[l];
        std::vector<uint64_t> first;
        std::unordered_set<uint64_t> fs;
        while (rd.next(line)) {
            if (!line.empty() && line[0] == '>') continue;
            size_t b0 = 0, e0 = line.size();
            while (b0 < e0 && isspace((unsigned char)line[b0])) ++b0;
            while (e0 > b0 && isspace((unsigned char)line[e0 - 1])) --e0;
            if (e0 - b0 < (size_t)k) continue;
            uint64_t km = 0;
            bool ok = true;
            for (int i = 0; i < k; ++i) {
                int c = bcode(line[b0 + i]);
                if (c < 0) { ok = false; break; }
                km = ((km << 2) | (uint64_t)c) & kmask;
            }
            if (!ok) continue;
            v.push_back(fmix64(km));
            if (first.size() < 2000 && fs.insert(km).second) first.push_back(km);
        }
        std::sort(v.begin(), v.end());
        v.erase(std::unique(v.begin(), v.end()), v.end());
        size_t both = 0;
        for (uint64_t km : first)
            if (std::binary_search(v.begin(), v.end(), fmix64(revcomp(km, k)))) ++both;
        canon[l] = !first.empty() && (double)both / (double)first.size() < 0.05;
    };
    {
        std::vector<std::thread> th;
        size_t nt = std::max(1, std::min<int>(threads, (int)L));
        for (size_t s = 0; s < nt; ++s)
            th.emplace_back([&, s]() { for (size_t l = s; l < L; l += nt) work(l); });
        for (auto& t : th) t.join();
    }
    count.resize(L);
    size_t tot = 0;
    for (size_t l = 0; l < L; ++l) { count[l] = hs[l].size(); tot += hs[l].size(); }
    canonical = true;
    for (size_t l = 0; l < L; ++l) canonical = canonical && canon[l];
    // k-way merge -> unique hashes with lineage bitmask
    h.clear(); mask.clear();
    h.reserve(tot); mask.reserve(tot);
    std::vector<size_t> pos(L, 0);
    for (;;) {
        uint64_t mn = ~0ULL;
        bool any = false;
        for (size_t l = 0; l < L; ++l)
            if (pos[l] < hs[l].size() && (!any || hs[l][pos[l]] < mn)) { mn = hs[l][pos[l]]; any = true; }
        if (!any) break;
        uint8_t m = 0;
        for (size_t l = 0; l < L; ++l)
            if (pos[l] < hs[l].size() && hs[l][pos[l]] == mn) { m |= (uint8_t)(1u << l); ++pos[l]; }
        h.push_back(mn);
        mask.push_back(m);
    }
    for (auto& v : hs) std::vector<uint64_t>().swap(v);
    if (h.size() >= 0xFFFFFFFFull) die("too many marker k-mers for 32-bit offsets");
    b = 8;
    while (b < 30 && (1ULL << b) < h.size() / 4 + 1) ++b;
    off.assign((1ULL << b) + 1, 0);
    for (uint64_t x : h) off[(x >> (64 - b)) + 1]++;
    for (size_t i = 1; i < off.size(); ++i) off[i] += off[i - 1];
}

void MarkerIndex::clean(const std::string& fasta, const std::unordered_map<std::string, std::string>& anc, int threads) {
    const size_t L = names.size();
    std::unordered_map<std::string, int> li;
    for (size_t l = 0; l < L; ++l) li[upper(names[l])] = (int)l;
    std::vector<std::pair<int, std::string>> seqs;
    size_t n_ref = 0, n_unl = 0;
    {
        LineReader rd(fasta);
        std::string_view l;
        int cur = -2;
        while (rd.next(l)) {
            if (!l.empty() && l[0] == '>') {
                auto ws = split_ws(l.substr(1));
                cur = -1;
                ++n_ref;
                if (!ws.empty()) {
                    std::string nm(ws[0]);
                    auto it = anc.find(nm);
                    if (it == anc.end()) it = anc.find(std::string(rsplit_dot(nm)));
                    if (it != anc.end()) {
                        auto jt = li.find(it->second);
                        if (jt != li.end()) cur = jt->second;
                    }
                }
                if (cur >= 0) seqs.emplace_back(cur, std::string());
                else ++n_unl;
                continue;
            }
            if (cur >= 0) for (char c : l) if (!isspace((unsigned char)c)) seqs.back().second.push_back(c);
        }
    }
    std::vector<uint8_t> allowed(h.size(), 0);
    const uint64_t kmask = k == 32 ? ~0ULL : ((1ULL << (2 * k)) - 1);
    const int shift = 2 * (k - 1);
    auto find = [&](uint64_t km) -> int64_t {
        const uint64_t x = fmix64(km);
        const uint64_t bk = x >> (64 - b);
        for (uint32_t i = off[bk], e = off[bk + 1]; i < e; ++i)
            if (h[i] == x) return (int64_t)i;
        return -1;
    };
    std::atomic<size_t> next{0};
    auto work = [&]() {
        for (;;) {
            size_t q = next.fetch_add(1);
            if (q >= seqs.size()) return;
            const uint8_t bit = (uint8_t)(1u << seqs[q].first);
            uint64_t fw = 0, rc = 0;
            int len = 0;
            for (char ch : seqs[q].second) {
                int c = bcode(ch);
                if (c < 0) { len = 0; fw = rc = 0; continue; }
                fw = ((fw << 2) | (uint64_t)c) & kmask;
                rc = (rc >> 2) | ((uint64_t)(3 - c) << shift);
                if (++len < k) continue;
                if (canonical) {
                    int64_t i = find(fw < rc ? fw : rc);
                    if (i >= 0) __atomic_fetch_or(&allowed[i], bit, __ATOMIC_RELAXED);
                } else {
                    int64_t i = find(fw);
                    if (i >= 0) __atomic_fetch_or(&allowed[i], bit, __ATOMIC_RELAXED);
                    i = find(rc);
                    if (i >= 0) __atomic_fetch_or(&allowed[i], bit, __ATOMIC_RELAXED);
                }
            }
        }
    };
    {
        std::vector<std::thread> th;
        for (int t = 0; t < std::max(1, threads); ++t) th.emplace_back(work);
        for (auto& t : th) t.join();
    }
    std::vector<uint64_t> before(L, 0), removed(L, 0);
    size_t keep = 0;
    for (size_t i = 0; i < h.size(); ++i) {
        uint8_t m = mask[i];
        for (size_t l = 0; l < L; ++l) if (m & (1u << l)) ++before[l];
        if (allowed[i]) {
            uint8_t nm = m & allowed[i];
            for (size_t l = 0; l < L; ++l) if ((m & (1u << l)) && !(nm & (1u << l))) ++removed[l];
            m = nm;
        }
        if (m) { h[keep] = h[i]; mask[keep] = m; ++keep; }
    }
    h.resize(keep);
    mask.resize(keep);
    for (size_t l = 0; l < L; ++l) count[l] = before[l] - removed[l];
    b = 8;
    while (b < 30 && (1ULL << b) < h.size() / 4 + 1) ++b;
    off.assign((1ULL << b) + 1, 0);
    for (uint64_t x : h) off[(x >> (64 - b)) + 1]++;
    for (size_t i = 1; i < off.size(); ++i) off[i] += off[i - 1];
    std::string msg;
    for (size_t l = 0; l < L; ++l)
        msg += sfmt("  %s %s -> %s (-%.2f%%)", names[l].c_str(), commas((long long)before[l]).c_str(), commas((long long)count[l]).c_str(),
                    100.0 * (double)removed[l] / (double)std::max<uint64_t>(before[l], 1));
    logm(sfmt("[INFO] marker cleaning against %s labelled reference sequences (%s without a marker lineage skipped):",
              commas((long long)seqs.size()).c_str(), commas((long long)n_unl).c_str()) + msg);
}

namespace {
const char KM_MAGIC[8] = {'L', 'Q', 'K', 'M', 'E', 'R', '0', '1'};
}

bool MarkerIndex::load_cache(const std::string& path, const std::string& sig) {
    FILE* f = fopen(path.c_str(), "rb");
    if (!f) return false;
    auto rd = [&](void* p, size_t n) { return fread(p, 1, n, f) == n; };
    char m[8];
    uint64_t n = 0;
    bool ok = rd(m, 8) && memcmp(m, KM_MAGIC, 8) == 0;
    std::string s;
    if (ok) { ok = rd(&n, 8); s.resize(n); ok = ok && rd(&s[0], n) && s == sig; }
    int64_t kk = 0, bb = 0;
    uint8_t cn = 0;
    uint64_t L = 0;
    if (ok) ok = rd(&kk, 8) && rd(&bb, 8) && rd(&cn, 1) && rd(&L, 8) && L <= 8;
    if (ok) {
        names.resize(L);
        for (auto& x : names) {
            uint64_t ln = 0;
            ok = ok && rd(&ln, 8);
            if (!ok) break;
            x.resize(ln);
            ok = ok && rd(&x[0], ln);
        }
        count.resize(L);
        ok = ok && rd(count.data(), 8 * L);
    }
    uint64_t N = 0;
    if (ok) ok = rd(&N, 8);
    if (ok) {
        k = (int)kk; b = (int)bb; canonical = cn != 0;
        off.resize((1ULL << b) + 1);
        h.resize(N);
        mask.resize(N);
        ok = rd(off.data(), off.size() * 4) && rd(h.data(), N * 8) && rd(mask.data(), N);
    }
    fclose(f);
    return ok;
}

void MarkerIndex::save_cache(const std::string& path, const std::string& sig) const {
    std::string tmp = path + ".tmp." + std::to_string(getpid());
    FILE* f = fopen(tmp.c_str(), "wb");
    if (!f) { logm("[WARN] cannot write marker cache " + path); return; }
    auto wr = [&](const void* p, size_t n) { if (n && fwrite(p, 1, n, f) != n) die("error writing marker cache"); };
    uint64_t n = sig.size();
    wr(KM_MAGIC, 8); wr(&n, 8); wr(sig.data(), n);
    int64_t kk = k, bb = b;
    uint8_t cn = canonical ? 1 : 0;
    uint64_t L = names.size();
    wr(&kk, 8); wr(&bb, 8); wr(&cn, 1); wr(&L, 8);
    for (auto& x : names) { uint64_t ln = x.size(); wr(&ln, 8); wr(x.data(), ln); }
    wr(count.data(), 8 * L);
    uint64_t N = h.size();
    wr(&N, 8);
    wr(off.data(), off.size() * 4); wr(h.data(), N * 8); wr(mask.data(), N);
    if (fclose(f) != 0) die("error writing marker cache");
    if (rename(tmp.c_str(), path.c_str()) != 0) logm("[WARN] cannot rename marker cache " + tmp);
}

// ================================================================== SNP index
void SnpIndex::load(const std::string& path) {
    LineReader rd(path);
    std::string_view l;
    if (!rd.next(l) || l.rfind("#LQSNP", 0) != 0) die("not an exported SNP DB (run lineage_snp/export_snp_db.py): " + path);
    for (auto f : split_tab(l)) {
        if (f.rfind("mode=", 0) == 0) msa = f.substr(5) == "msa";
        else if (f.rfind("lineages=", 0) == 0) {
            std::string s(f.substr(9));
            size_t b0 = 0;
            for (;;) {
                size_t e = s.find(',', b0);
                lineages.push_back(s.substr(b0, e == std::string::npos ? std::string::npos : e - b0));
                if (e == std::string::npos) break;
                b0 = e + 1;
            }
        }
    }
    const size_t L = lineages.size();
    if (L == 0) die("SNP DB without lineages: " + path);
    std::unordered_map<std::string, uint32_t> seen;
    std::string cur;
    while (rd.next(l)) {
        if (l.empty()) continue;
        auto q = split_tab(l);
        if (q.size() < 3) die("malformed SNP DB line: " + std::string(l));
        if (tids.empty() || q[0] != cur) {
            cur.assign(q[0]);
            if (seen.count(cur)) die("SNP DB not grouped by transcript: " + cur);
            seen.emplace(cur, (uint32_t)tids.size());
            tids.push_back(cur);
            beg.push_back(pos.size());
        }
        pos.push_back((int32_t)parse_ll(q[1]));
        if (msa) {
            std::string_view e = q[2];
            for (size_t i = 0; i < L; ++i) {
                int c = i < e.size() ? bcode(e[i]) : -1;
                exp.push_back((uint8_t)(c < 0 ? 4 : c));
            }
        } else {
            if (q.size() < 5) die("malformed population SNP DB line: " + std::string(l));
            int rc = q[2].size() == 1 ? bcode(q[2][0]) : -1, ac = q[3].size() == 1 ? bcode(q[3][0]) : -1;
            ref.push_back((uint8_t)(rc < 0 ? 4 : rc));
            alt.push_back((uint8_t)(ac < 0 ? 4 : ac));
            std::string s(q[4]);
            size_t b0 = 0;
            for (size_t i = 0; i < L; ++i) {
                size_t e = s.find(',', b0);
                double v = 0.5;
                py_float(std::string_view(s).substr(b0, e == std::string::npos ? std::string::npos : e - b0), v);
                af.push_back(v);
                b0 = e == std::string::npos ? s.size() : e + 1;
            }
        }
    }
    beg.push_back(pos.size());
    for (size_t i = 0; i + 1 < beg.size(); ++i)
        for (uint64_t j = beg[i] + 1; j < beg[i + 1]; ++j)
            if (pos[j] <= pos[j - 1]) die("SNP DB positions not sorted within " + tids[i]);
    for (uint32_t i = 0; i < tids.size(); ++i) resolve[tids[i]] = i;
    for (uint32_t i = 0; i < tids.size(); ++i) {
        std::string base(rsplit_dot(tids[i]));
        if (base != tids[i]) resolve.emplace(base, i);
    }
}

int SnpIndex::lookup(const std::string& r) const {
    auto it = resolve.find(r);
    if (it != resolve.end()) return (int)it->second;
    it = resolve.find(std::string(rsplit_dot(r)));
    return it == resolve.end() ? -1 : (int)it->second;
}

// ================================================================== engine
static int status_class_w(std::string_view s, const Opts& o) {
    auto ist = [&](std::string_view p) {
        if (s.size() < p.size()) return false;
        for (size_t i = 0; i < p.size(); ++i) if (tolower((unsigned char)s[i]) != p[i]) return false;
        return true;
    };
    if (ist("specific")) return 0;
    std::string lc(s);
    for (auto& c : lc) c = (char)tolower((unsigned char)c);
    for (const char* k : {"none", "no_snp", "chimeric", "ambiguous", "weak"})
        if (lc.find(k) != std::string::npos) return 2;
    return 1;
}

void LineageEngine::init(const Opts& o, const Model& M) {
    M_ = &M;
    o_ = o;
    const size_t R = M.refs.size();
    valid_ref_.assign(R, 1);
    ref_gene_dot1_.resize(R);
    ref_gene_rdot_.resize(R);
    for (size_t r = 0; r < R; ++r) {
        const std::string& t = M.refs[r];
        size_t d = t.find('.');
        ref_gene_dot1_[r] = d == std::string::npos ? t : t.substr(0, d);
        ref_gene_rdot_[r] = std::string(rsplit_dot(t));
    }
    if (!o.target_ids.empty()) {
        std::unordered_set<std::string> valid;
        LineReader rd(o.target_ids);
        std::string_view l;
        while (rd.next(l)) {
            size_t b0 = 0, e0 = l.size();
            while (b0 < e0 && isspace((unsigned char)l[b0])) ++b0;
            while (e0 > b0 && isspace((unsigned char)l[e0 - 1])) --e0;
            if (b0 == e0) continue;
            std::string s(l.substr(b0, e0 - b0));
            size_t d = s.find('.');
            valid.insert(d == std::string::npos ? s : s.substr(0, d));
        }
        size_t nv = 0;
        for (size_t r = 0; r < R; ++r) { valid_ref_[r] = valid.count(ref_gene_dot1_[r]) ? 1 : 0; nv += valid_ref_[r]; }
        logm("[INFO] target transcripts for lineage evidence: " + commas((long long)nv) + " / " + commas((long long)R));
    }
    // ---- k-mer markers
    if (!o.markers.empty() && !o.no_kmer) {
        std::vector<std::pair<std::string, std::string>> np;
        std::string sig = "k=" + std::to_string(o.kmer_len);
        for (auto& m : o.markers) {
            size_t c = m.find(':');
            if (c == std::string::npos || c == 0) die("--markers must be NAME:PATH, got " + m);
            std::string nm = m.substr(0, c), p = m.substr(c + 1);
            struct stat sb;
            if (stat(p.c_str(), &sb) != 0) die("marker file not found: " + p);
            np.emplace_back(nm, p);
            sig += "|" + nm + "|" + p + "|" + std::to_string((long long)sb.st_size) + "|" + std::to_string((long long)sb.st_mtime);
        }
        if (o.marker_clean) {   // 1.1.0: the cache records that (and against which reference) the markers were cleaned
            if (o.fasta.empty() || o.ancestry.empty()) die("--marker-clean needs --fasta and -a");
            struct stat sf, sa;
            if (stat(o.fasta.c_str(), &sf) != 0 || stat(o.ancestry.c_str(), &sa) != 0) die("--marker-clean: cannot stat --fasta / -a");
            sig += "|clean|" + o.fasta + "|" + std::to_string((long long)sf.st_size) + "|" + std::to_string((long long)sf.st_mtime) +
                   "|" + o.ancestry + "|" + std::to_string((long long)sa.st_size) + "|" + std::to_string((long long)sa.st_mtime);
        }
        double t0 = wall();
        bool cached = !o.marker_cache.empty() && mk_.load_cache(o.marker_cache, sig);
        if (!cached) {
            mk_.build(np, o.kmer_len, std::max(1, o.threads));
            if (o.marker_clean) mk_.clean(o.fasta, M.anc, std::max(1, o.threads));
            if (!o.marker_cache.empty()) mk_.save_cache(o.marker_cache, sig);
        }
        for (size_t l = 0; l < mk_.names.size(); ++l)
            logm("[INFO] marker " + mk_.names[l] + ": " + commas((long long)mk_.count[l]) + " " + std::to_string(mk_.k) + "-mers");
        logm(std::string("[INFO] marker set is ") + (mk_.canonical ? "canonical -> canonical lookups" : "double-stranded -> both strands scanned") +
             sfmt(" | %s index %.1f GB, %.0fs", cached ? "cached" : "built", (mk_.h.size() * 9.0 + mk_.off.size() * 4.0) / 1e9, wall() - t0));
        logm("[INFO] k-mer: Beta alpha=" + pyfloat(o.kmer_alpha) + "  min_hits=" + std::to_string(o.min_hits) +
             "  specific_thresh=" + pyfloat(o.spec_thresh));
        kmer_on = true;
        kmer_names = mk_.names;
        for (auto& s : kmer_names) kmer_lins.push_back(upper(s));
    }
    // ---- SNP DB
    if (!o.snp_db.empty()) {
        double t0 = wall();
        snp_.load(o.snp_db);
        snp_on = true;
        for (auto& s : snp_.lineages) snp_lins.push_back(upper(s));
        snp_tid_.assign(R, -1);
        size_t nt = 0;
        for (size_t r = 0; r < R; ++r) { snp_tid_[r] = snp_.lookup(M.refs[r]); nt += snp_tid_[r] >= 0; }
        logm("[INFO] SNP DB (" + std::string(snp_.msa ? "msa" : "population") + "): " + commas((long long)snp_.pos.size()) +
             " SNPs on " + commas((long long)snp_.tids.size()) + " transcripts; " + commas((long long)nt) +
             " reference transcripts resolved" + sfmt(" (%.0fs) | error rate %s | pick %s", wall() - t0,
             pyfloat(o.error_rate).c_str(), o.snp_pick.c_str()));
    }
    // ---- lineage axis of the P table / EM
    if (kmer_on && snp_on) {
        if (kmer_lins == snp_lins) em_lins = kmer_lins;
        else {
            std::set<std::string> s(kmer_lins.begin(), kmer_lins.end());
            s.insert(snp_lins.begin(), snp_lins.end());
            em_lins.assign(s.begin(), s.end());
        }
    } else if (kmer_on) em_lins = kmer_lins;
    else if (snp_on) em_lins = snp_lins;
    if (em_lins.size() > 8) die("at most 8 lineages are supported");
    map_k_.assign(em_lins.size(), -1);
    map_s_.assign(em_lins.size(), -1);
    for (size_t i = 0; i < em_lins.size(); ++i) {
        for (size_t j = 0; j < kmer_lins.size(); ++j) if (kmer_lins[j] == em_lins[i]) map_k_[i] = (int)j;
        for (size_t j = 0; j < snp_lins.size(); ++j) if (snp_lins[j] == em_lins[i]) map_s_[i] = (int)j;
    }
    LM_ = std::log(1.0 - o.error_rate);
    LX_ = std::log(o.error_rate / 3.0);
}

void LineageEngine::kmer_eval(const Rec* recs, size_t n, double* p, int& status, long long& hits) const {
    const size_t L = mk_.names.size();
    long long h[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    std::string_view s1, s2;
    for (size_t i = 0; i < n; ++i) {
        const Rec& r = recs[i];
        if ((r.flag & 4) || r.rid == NOCL || !valid_ref_[r.rid] || r.seq == "*" || r.seq.empty()) continue;
        if (r.flag & 0x80) { if (s2.empty()) s2 = r.seq; }
        else if (s1.empty()) s1 = r.seq;
    }
    const int k = mk_.k;
    const uint64_t kmask = k == 32 ? ~0ULL : ((1ULL << (2 * k)) - 1);
    const int shift = 2 * (k - 1);
    for (std::string_view s : {s1, s2}) {
        uint64_t fw = 0, rc = 0;
        int len = 0;
        for (char ch : s) {
            int c = bcode(ch);
            if (c < 0) { len = 0; fw = rc = 0; continue; }
            fw = ((fw << 2) | (uint64_t)c) & kmask;
            rc = (rc >> 2) | ((uint64_t)(3 - c) << shift);
            if (++len < k) continue;
            if (mk_.canonical) {
                uint8_t m = mk_.get(fw < rc ? fw : rc);
                for (size_t l = 0; m; ++l, m >>= 1) h[l] += m & 1;
            } else {
                uint8_t m = mk_.get(fw);
                for (size_t l = 0; m; ++l, m >>= 1) h[l] += m & 1;
                m = mk_.get(rc);
                for (size_t l = 0; m; ++l, m >>= 1) h[l] += m & 1;
            }
        }
    }
    long long total = 0;
    for (size_t l = 0; l < L; ++l) total += h[l];
    hits = total;
    if (total == 0) {
        for (size_t l = 0; l < L; ++l) p[l] = 1.0 / (double)L;
        status = 3;
        return;
    }
    double top = 0;
    for (size_t l = 0; l < L; ++l) {
        p[l] = ((double)h[l] + o_.kmer_alpha) / ((double)total + (double)L * o_.kmer_alpha);
        top = std::max(top, p[l]);
    }
    if (total < o_.min_hits) status = 2;
    else if (top >= o_.spec_thresh) status = 0;
    else status = 1;
}

void LineageEngine::snp_eval(const Rec* recs, const std::vector<uint32_t>& idx, LinScratch& s, double* p, int& status,
                             int& best, long long& nsnp) const {
    const size_t L = snp_.lineages.size();
    best = -1;
    nsnp = 0;
    const int tid = snp_tid_[recs[idx[0]].rid];
    auto uniform = [&](int st) { for (size_t i = 0; i < L; ++i) p[i] = 1.0 / (double)L; status = st; };
    if (tid < 0) { uniform(2); return; }
    s.ll.assign(L, 0.0);
    double* ll = s.ll.data();
    const uint64_t b0 = snp_.beg[tid], b1 = snp_.beg[tid + 1];
    const int32_t* P = snp_.pos.data();
    const double eps = o_.error_rate;
    long long total = 0;
    for (uint32_t ri : idx) {
        const Rec& r = recs[ri];
        if (r.seq == "*") continue;
        long long rpos = r.pos - 1, qpos = 0;
        for_cigar(r.cigar, [&](char op, long long len) {
            if (op == 'M' || op == '=' || op == 'X') {
                const int32_t* it = std::lower_bound(P + b0, P + b1, (int32_t)rpos);
                for (; it != P + b1 && *it < rpos + len; ++it) {
                    long long q = qpos + (*it - rpos);
                    if (q < 0 || q >= (long long)r.seq.size()) continue;
                    int c = bcode(r.seq[q]);
                    if (c < 0) continue;
                    const uint64_t j = (uint64_t)(it - P);
                    if (snp_.msa) {
                        const uint8_t* e = &snp_.exp[j * L];
                        for (size_t i = 0; i < L; ++i) ll[i] += (e[i] == c) ? LM_ : LX_;
                    } else {
                        const int rf = snp_.ref[j], al = snp_.alt[j];
                        const double* a = &snp_.af[j * L];
                        for (size_t i = 0; i < L; ++i) {
                            double af = a[i], lik;
                            if (c == rf) lik = (1 - af) * (1 - eps) + af * (eps / 3);
                            else if (c == al) lik = af * (1 - eps) + (1 - af) * (eps / 3);
                            else lik = eps / 3;
                            ll[i] += std::log(std::max(lik, eps / 3));
                        }
                    }
                    ++total;
                }
                rpos += len;
                qpos += len;
            } else if (op == 'I' || op == 'S') qpos += len;
            else if (op == 'D' || op == 'N') rpos += len;
        });
    }
    nsnp = total;
    if (total == 0) { uniform(2); return; }
    double mx = ll[0];
    for (size_t i = 1; i < L; ++i) mx = std::max(mx, ll[i]);
    if (snp_.msa && total >= 3 && !o_.no_chimeric && mx < std::max(-6.0, -2.0 * (double)total)) { uniform(3); return; }
    std::vector<double> lk(L);
    for (size_t i = 0; i < L; ++i) lk[i] = std::exp(ll[i] - mx);
    double sm = np_pairwise(lk.data(), L);
    for (size_t i = 0; i < L; ++i) p[i] = lk[i] / sm;
    int bi = 0;
    for (size_t i = 1; i < L; ++i) if (p[i] > p[bi]) bi = (int)i;
    best = bi;
    status = p[bi] > (snp_.msa ? 0.99 : 0.9) ? 0 : 1;
}

static double entropy_np(const double* v, size_t n) {
    double s = np_pairwise(v, n);
    if (s <= 0) return 0.0;
    std::vector<double> q;
    for (size_t i = 0; i < n; ++i) { double x = v[i] / s; if (x > 0) q.push_back(x); }
    if (q.size() <= 1) return 0.0;
    for (auto& x : q) x = x * std::log2(x);
    return -np_pairwise(q.data(), q.size());
}

static const char* KST[4] = {"Specific", "Conflict_Weighted", "Weak_Evidence", "Ambiguous_None"};

bool LineageEngine::eval(const Rec* recs, size_t n, std::string_view qname, LinScratch& s, LineageStats& st, double* pv,
                         std::string* dbg_kmer, std::string* dbg_snp, std::string* prow) const {
    if (!kmer_on && !snp_on) return false;
    s.rr.clear();
    for (uint32_t i = 0; i < n; ++i) {
        const Rec& r = recs[i];
        if (!(r.flag & 4) && r.rid != NOCL && valid_ref_[r.rid]) s.rr.emplace_back(r.rid, i);
    }
    if (s.rr.empty()) return false;
    ++st.valid;
    std::sort(s.rr.begin(), s.rr.end());
    char buf[64];
    auto fmt4 = [&](double x) { snprintf(buf, sizeof buf, "%.4f", x); return std::string(buf); };

    // ---- k-mer (one row per fragment: identical for every transcript in snipe_kmer)
    double pk[8];
    int kst = -1;
    long long kh = 0;
    bool has_k = false;
    if (kmer_on) {
        bool any_seq = false;
        for (auto& x : s.rr) if (recs[x.second].seq != "*" && !recs[x.second].seq.empty()) { any_seq = true; break; }
        if (any_seq) {
            kmer_eval(recs, n, pk, kst, kh);
            for (size_t l = 0; l < kmer_lins.size(); ++l) pk[l] = py_round4(pk[l]);
            has_k = true;
            st.kmer_st[kst]++;
            st.hits[kh < 10 ? kh : 10]++;
            if (dbg_kmer) {
                dbg_kmer->append(qname);
                for (size_t l = 0; l < kmer_lins.size(); ++l) { dbg_kmer->push_back('\t'); dbg_kmer->append(fmt4(pk[l])); }
                dbg_kmer->append(sfmt("\t%s\t%lld\n", KST[kst], kh));
            }
        }
    }

    // ---- SNP: one row per (fragment, transcript); the EM gets the row of one deterministic transcript
    double ps[8];
    int sst = -1, sbest = -1;
    bool has_s = false;
    uint32_t s_rid = NOCL;
    if (snp_on) {
        struct Cand { uint32_t rid; int nm; std::vector<uint32_t> idx; };
        std::vector<Cand> cands;
        for (size_t i = 0; i < s.rr.size();) {
            size_t j = i;
            Cand c{s.rr[i].first, 0, {}};
            int m1 = NM_INF, m2 = NM_INF;
            bool seq = false;
            while (j < s.rr.size() && s.rr[j].first == c.rid) {
                const Rec& r = recs[s.rr[j].second];
                c.idx.push_back(s.rr[j].second);
                if (r.seq != "*") seq = true;
                if (r.flag & 0x80) m2 = std::min(m2, r.nm); else m1 = std::min(m1, r.nm);
                ++j;
            }
            c.nm = (m1 < NM_INF ? m1 : 0) + (m2 < NM_INF ? m2 : 0);
            if (seq) cands.push_back(std::move(c));
            i = j;
        }
        if (!cands.empty()) {
            const size_t L = snp_lins.size();
            double tp[8];
            int tst, tbest;
            long long tn;
            long long best_n = -1;
            int best_nm = NM_INF;
            if (dbg_snp) {
                for (auto& c : cands) {
                    snp_eval(recs, c.idx, s, tp, tst, tbest, tn);
                    dbg_snp->append(qname);
                    dbg_snp->push_back('\t');
                    dbg_snp->append(ref_gene_rdot_[c.rid]);
                    for (size_t l = 0; l < L; ++l) { dbg_snp->push_back('\t'); dbg_snp->append(fmt4(tp[l])); }
                    if (tst == 0) dbg_snp->append("\tSpecific_Soft(" + snp_.lineages[tbest] + ")");
                    else dbg_snp->append(tst == 1 ? "\tMulti_Lineage" : (tst == 2 ? "\tNo_SNP_Covered" : "\tChimeric_Drop"));
                    dbg_snp->append(sfmt("\t%lld\n", tn));
                }
            }
            int min_nm = NM_INF;
            for (auto& c : cands) min_nm = std::min(min_nm, c.nm);
            for (auto& c : cands) {
                if (o_.snp_pick == "best" && c.nm != min_nm) continue;
                if (o_.snp_pick == "first" && s_rid != NOCL) break;
                snp_eval(recs, c.idx, s, tp, tst, tbest, tn);
                ++st.snp_rows;
                bool take = false;
                if (s_rid == NOCL) take = true;
                else if (o_.snp_pick == "best") take = tn > best_n;
                else if (o_.snp_pick == "most") take = tn > best_n || (tn == best_n && c.nm < best_nm);
                if (take) {
                    s_rid = c.rid; best_n = tn; best_nm = c.nm;
                    sst = tst; sbest = tbest;
                    for (size_t l = 0; l < L; ++l) ps[l] = tp[l];
                }
            }
            for (size_t l = 0; l < L; ++l) ps[l] = py_round4(ps[l]);
            has_s = true;
            st.snp_st[sst]++;
        }
    }

    // ---- fusion (merge_poe) or the single available route
    const size_t NL = em_lins.size();
    double row[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    std::string status;
    if (has_k && has_s) {
        ++st.both;
        const size_t Lk = kmer_lins.size(), Ls = snp_lins.size();
        double wk = (kst == 0 ? 10.0 : (kst == 1 ? 2.5 : 0.5)) +
                    std::max(0.0, std::log2((double)Lk) - entropy_np(pk, Lk)) * 5.0;
        double ws = (sst == 0 ? 10.0 : (sst == 1 ? 2.5 : 0.5)) +
                    std::max(0.0, std::log2((double)Ls) - entropy_np(ps, Ls)) * 5.0;
        double tot_w = wk + ws;
        double ak = tot_w > 0 ? wk / tot_w : 0.5;
        ak = std::max(o_.poe_min_weight, std::min(1 - o_.poe_min_weight, ak));
        double as = 1 - ak;
        double raw[8], tot = 0;
        for (size_t i = 0; i < NL; ++i) {
            double a1 = map_k_[i] >= 0 ? pk[map_k_[i]] : 1.0 / (double)NL;
            double a2 = map_s_[i] >= 0 ? ps[map_s_[i]] : 1.0 / (double)NL;
            raw[i] = std::pow(a1 + o_.poe_eps, ak) * std::pow(a2 + o_.poe_eps, as);
            tot += raw[i];
        }
        for (size_t i = 0; i < NL; ++i) row[i] = tot > 0 ? raw[i] / tot : 1.0 / (double)NL;
        size_t bi = 0;
        for (size_t i = 1; i < NL; ++i) if (row[i] > row[bi]) bi = i;
        if (row[bi] > 0.9) { status = "Specific_Soft(" + em_lins[bi] + ")"; st.fused_st[0]++; }
        else { status = "Multi_Lineage"; st.fused_st[1]++; }
    } else if (has_k) {
        for (size_t i = 0; i < NL; ++i) row[i] = map_k_[i] >= 0 ? pk[map_k_[i]] : 1.0 / (double)NL;
        if (snp_on) {
            ++st.kmer_only;
            size_t bi = 0;
            for (size_t i = 1; i < NL; ++i) if (row[i] > row[bi]) bi = i;
            if (row[bi] > 0.9) { status = "Specific_Kmer(" + em_lins[bi] + ")"; st.fused_st[2]++; }
            else { status = "Multi_Kmer"; st.fused_st[3]++; }
        } else status = KST[kst];
    } else if (has_s) {
        for (size_t i = 0; i < NL; ++i) row[i] = map_s_[i] >= 0 ? ps[map_s_[i]] : 1.0 / (double)NL;
        if (kmer_on) {
            ++st.snp_only;
            size_t bi = 0;
            for (size_t i = 1; i < NL; ++i) if (row[i] > row[bi]) bi = i;
            if (row[bi] > 0.9) { status = "Specific_SNP(" + em_lins[bi] + ")"; st.fused_st[4]++; }
            else { status = "Multi_SNP"; st.fused_st[5]++; }
        } else if (sst == 0) status = "Specific_Soft(" + snp_.lineages[sbest] + ")";
        else status = sst == 1 ? "Multi_Lineage" : (sst == 2 ? "No_SNP_Covered" : "Chimeric_Drop");
    } else return false;
    for (size_t i = 0; i < NL; ++i) row[i] = py_round4(row[i]);
    if (prow) {
        uint32_t grid = s_rid != NOCL ? s_rid : s.rr[0].first;
        prow->append(qname);
        prow->push_back('\t');
        prow->append(ref_gene_dot1_[grid]);
        for (size_t i = 0; i < NL; ++i) { prow->push_back('\t'); prow->append(fmt4(row[i])); }
        prow->push_back('\t');
        prow->append(status);
        prow->push_back('\n');
    }
    // ---- EM factors (em_global: status weight, p_alpha, p_floor)
    const int wc = status_class_w(status, o_);
    const double w = wc == 0 ? o_.pw_spec : (wc == 1 ? o_.pw_conf : o_.pw_none);
    if (w == 0.0) return false;
    for (size_t i = 0; i < NL; ++i) {
        double x = row[i];
        if (o_.p_alpha > 0) x = (x + o_.p_alpha) / (1.0 + (double)NL * o_.p_alpha);
        x = std::max(x, o_.p_floor);
        if (w != 1.0) x = std::pow(x, w);
        pv[i] = x;
    }
    ++st.pv_used;
    return true;
}

}  // namespace lq
