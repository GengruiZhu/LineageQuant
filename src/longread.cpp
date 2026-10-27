// lqcore lrpick / lrcount -- long-read physical prior (Step 00b) in C++
//
//   lrpick   FASTQ -> FASTQ with one read per molecule (PacBio subreads: every pass of a ZMW is the same molecule).
//            ZMWs come from the read names (movie/zmw/start_end) when they are there. SRA strips them
//            ("SRR37812062.1 1/1"), so otherwise consecutive reads are compared: two passes of one ZMW are reverse
//            complements of each other, so read i is joined to read i-1 when enough of its sampled k-mers are found
//            in the reverse complement of read i-1. A lag-L comparison (reads that can never be one ZMW) measures
//            the chance level; grouping is switched on only when lag-1 matches clearly exceed it. Frequent k-mers
//            (adapters, primers, poly(A), the most abundant transcripts) are masked first.
//            Per group: >= 3 subreads -> the median-length one among the inner passes (first/last are partial);
//            1-2 -> the longest.
//   lrcount  long-read alignments (minimap2 -a with secondaries) -> molecules per allele copy. Same algorithm and
//            byte-identical output as lineagequant/longread/extract_prior.py (LongRead_Counts.tsv, --legacy-out,
//            --summary).
#include "common.hpp"

#include <array>
#include <iterator>
#include <map>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <csignal>
#include <cerrno>
#include <fcntl.h>
#include <sys/wait.h>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace lq {
namespace {

constexpr size_t NPOS = std::string_view::npos;

bool all_digits(std::string_view s) {
    if (s.empty()) return false;
    for (char c : s)
        if (c < '0' || c > '9') return false;
    return true;
}
// Python re.match(r"^(.+/\d+)/(\d+_\d+|ccs.*)$", name): zmw = group 1, sub = group 2 is not "ccs..."
bool pb_name(std::string_view name, std::string_view& zmw, bool& sub) {
    size_t p = name.rfind('/');
    while (p != NPOS && p > 0) {
        std::string_view tail = name.substr(p + 1), head = name.substr(0, p);
        size_t u = tail.find('_');
        bool st = u != NPOS && all_digits(tail.substr(0, u)) && all_digits(tail.substr(u + 1));
        bool cc = tail.size() >= 3 && tail.substr(0, 3) == "ccs";
        if (st || cc) {
            size_t q = head.rfind('/');
            if (q != NPOS && q >= 1 && all_digits(head.substr(q + 1))) {
                zmw = head;
                sub = !cc;
                return true;
            }
        }
        p = name.rfind('/', p - 1);
    }
    return false;
}

std::string sh_quote(const std::string& s) {
    std::string o = "'";
    for (char c : s) o += c == '\'' ? std::string("'\\''") : std::string(1, c);
    return o + "'";
}

// ======================================================================================== lrpick
struct Km {
    uint32_t code, pos;
};
struct PRead {
    std::string name, seq, qual;
    std::string zmw;          // empty: no ZMW in the name
    bool sub = false;
    std::vector<Km> F, R;     // sampled k-mers occurring once in the read / in its reverse complement, by code
};

struct PickOpts {
    std::vector<std::string> in;
    std::string out = "-", stats, label = "-", group = "auto", debug;
    int k = 12, sample_mask = 3, warmup = 20000, lag = 101, min_shared = 4, max_group = 1000, threads = 4;
    double mask_frac = 0.05, min_ovf = 0.8, min_cov = 0.5;
    bool decide_only = false;
};

inline uint32_t mix32(uint32_t x) {
    x ^= x >> 16; x *= 0x7feb352dU; x ^= x >> 15; x *= 0x846ca68bU; x ^= x >> 16;
    return x;
}

class Sketcher {
public:
    int k;
    uint32_t smask, kmask;
    std::vector<uint8_t> masked;          // one byte per k-mer code (k <= 14)
    explicit Sketcher(int k_, int sm) : k(k_), smask((uint32_t)sm), kmask((uint32_t)((1ull << (2 * k_)) - 1)),
                                         masked((size_t)1 << (2 * k_), 0) {}
    static void keep_unique(std::vector<Km>& v) {
        std::sort(v.begin(), v.end(), [](const Km& a, const Km& b) { return a.code < b.code || (a.code == b.code && a.pos < b.pos); });
        size_t o = 0;
        for (size_t i = 0; i < v.size();) {
            size_t j = i + 1;
            while (j < v.size() && v[j].code == v[i].code) ++j;
            if (j == i + 1) v[o++] = v[i];
            i = j;
        }
        v.resize(o);
    }
    // F: k-mers of the read at their start position; R: k-mers of its reverse complement at their position there
    void sketch(std::string_view s, std::vector<Km>& F, std::vector<Km>& R, bool use_mask) const {
        F.clear(); R.clear();
        uint32_t f = 0, r = 0;
        int len = 0;
        const int sh = 2 * (k - 1);
        const uint32_t L = (uint32_t)s.size();
        for (uint32_t i = 0; i < L; ++i) {
            int c;
            switch (s[i]) {
                case 'A': case 'a': c = 0; break;
                case 'C': case 'c': c = 1; break;
                case 'G': case 'g': c = 2; break;
                case 'T': case 't': c = 3; break;
                default: c = -1;
            }
            if (c < 0) { len = 0; f = r = 0; continue; }
            f = ((f << 2) | (uint32_t)c) & kmask;
            r = (r >> 2) | ((uint32_t)(3 - c) << sh);
            if (++len < k) continue;
            const uint32_t st = i + 1 - (uint32_t)k;
            if ((mix32(f) & smask) == 0 && !(use_mask && masked[f])) F.push_back({f, st});
            if ((mix32(r) & smask) == 0 && !(use_mask && masked[r])) R.push_back({r, L - (uint32_t)k - st});
        }
        keep_unique(F);
        keep_unique(R);
    }
};

// Is read A (forward) the reverse complement of read B, i.e. another pass of the same molecule?
// Shared sampled k-mers of A and rc(B) must lie on one diagonal, the shorter read must lie inside the overlap
// that diagonal implies (passes are nested: inner passes are the whole insert, the first / last are partial),
// and they must span most of that overlap (two arrays that only share one transcript match on a short stretch).
struct DiagHit {
    int shared = 0, band = 0;
    double ovf = 0, cov = 0;
    bool ok = false;
};
DiagHit diag_match(const std::vector<Km>& A, size_t La, const std::vector<Km>& B, size_t Lb, int k, const PickOpts& o) {
    DiagHit h;
    std::vector<std::pair<long, long>> dp;          // (diagonal, position in A)
    size_t i = 0, j = 0;
    while (i < A.size() && j < B.size()) {
        if (A[i].code < B[j].code) ++i;
        else if (B[j].code < A[i].code) ++j;
        else { dp.emplace_back((long)A[i].pos - (long)B[j].pos, (long)A[i].pos); ++i; ++j; }
    }
    h.shared = (int)dp.size();
    if (h.shared < o.min_shared) return h;
    std::sort(dp.begin(), dp.end());
    const long mn = (long)std::min(La, Lb);
    const long W = std::max(100L, (long)(0.03 * (double)mn));
    size_t bl = 0, bc = 0;
    for (size_t l = 0, r = 0; l < dp.size(); ++l) {
        if (r < l) r = l;
        while (r + 1 < dp.size() && dp[r + 1].first - dp[l].first <= W) ++r;
        if (r - l + 1 > bc) { bc = r - l + 1; bl = l; }
    }
    h.band = (int)bc;
    if (h.band < o.min_shared) return h;
    long pmin = dp[bl].second, pmax = dp[bl].second;
    for (size_t t = bl; t < bl + bc; ++t) { pmin = std::min(pmin, dp[t].second); pmax = std::max(pmax, dp[t].second); }
    const long dc = dp[bl + bc / 2].first;
    const long lo = std::max(0L, dc), hi = std::min((long)La, (long)Lb + dc);
    if (hi <= lo) return h;
    const double ov = (double)(hi - lo);
    h.ovf = ov / (double)mn;
    h.cov = std::min(1.0, (double)(pmax - pmin + k) / ov);
    h.ok = h.ovf >= o.min_ovf && h.cov >= o.min_cov;
    return h;
}

long long n50(std::vector<uint32_t> v) {
    if (v.empty()) return 0;
    std::sort(v.rbegin(), v.rend());
    long double tot = 0;
    for (auto x : v) tot += x;
    long double acc = 0;
    for (auto x : v) { acc += x; if (acc >= tot / 2) return x; }
    return 0;
}

void pick_usage() {
    std::fprintf(stderr,
        "usage: lqcore lrpick [-i FQ[.gz] ...|-] [-o OUT.fq|-] --stats STATS.tsv [--label NAME]\n"
        "        [--group auto|zmw|adjacent|none] [--k 12] [--warmup 20000] [--lag 101] [--mask-frac 0.05]\n"
        "        [--min-shared 4] [--min-overlap 0.8] [--min-cover 0.5] [-t 4] [--debug-pairs FILE] [--decide-only]\n"
        "  one read per molecule: ZMW from the read name, or (no ZMW in the names, e.g. SRA) consecutive reads\n"
        "  that are reverse complements of each other (passes of one ZMW); STATS also gives the minimap2 preset\n"
        "  --decide-only: read only the warm-up reads, write STATS (kind / preset / grouping) and stop, no FASTQ\n");
}

int lrpick_main(int argc, char** argv) {
    PickOpts o;
    auto need = [&](int& i) -> std::string { if (i + 1 >= argc) { pick_usage(); die(std::string("missing value for ") + argv[i]); } return argv[++i]; };
    for (int i = 0; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "-i") { o.in.push_back(need(i)); while (i + 1 < argc && argv[i + 1][0] != '-') o.in.push_back(argv[++i]); }
        else if (a == "-o") o.out = need(i);
        else if (a == "--stats") o.stats = need(i);
        else if (a == "--label") o.label = need(i);
        else if (a == "--group") o.group = need(i);
        else if (a == "--k") o.k = std::atoi(need(i).c_str());
        else if (a == "--warmup") o.warmup = std::max(100, std::atoi(need(i).c_str()));
        else if (a == "--lag") o.lag = std::max(2, std::atoi(need(i).c_str()));
        else if (a == "--mask-frac") o.mask_frac = std::atof(need(i).c_str());
        else if (a == "--min-overlap") o.min_ovf = std::atof(need(i).c_str());
        else if (a == "--min-cover") o.min_cov = std::atof(need(i).c_str());
        else if (a == "--debug-pairs") o.debug = need(i);
        else if (a == "--decide-only") o.decide_only = true;
        else if (a == "-t" || a == "--threads") o.threads = std::max(1, std::atoi(need(i).c_str()));
        else if (a == "--min-shared") o.min_shared = std::atoi(need(i).c_str());
        else if (a == "-h" || a == "--help") { pick_usage(); return 0; }
        else { pick_usage(); die("unknown option " + a); }
    }
    if (o.in.empty()) o.in.push_back("-");
    if (o.stats.empty()) { pick_usage(); die("--stats is required"); }
    if (o.k < 8 || o.k > 14) die("--k must be 8..14");
    if (o.group != "auto" && o.group != "zmw" && o.group != "adjacent" && o.group != "none") die("--group: auto|zmw|adjacent|none");
    const double T0 = wall();
#ifdef _OPENMP
    omp_set_num_threads(o.threads);
#endif

    // FASTQ to stdout: keep the original stdout for it and send the log lines (logm -> stdout) to stderr
    FILE* out = nullptr;
    if (o.out == "-") {
        std::fflush(stdout);
        int fd = dup(1);
        if (fd < 0 || dup2(2, 1) < 0) die("cannot redirect stdout");
        out = fdopen(fd, "w");
    } else out = fopen(o.out.c_str(), "w");
    if (!out) die("cannot write " + o.out);
    std::vector<char> obuf(1 << 24);
    setvbuf(out, obuf.data(), _IOFBF, obuf.size());

    Sketcher SK(o.k, o.sample_mask);
    long long reads_in = 0, n_sub = 0, n_ccs = 0, n_nozmw = 0, reappeared = 0, groups_out = 0;
    long double bases_in = 0, bases_out = 0, qsum = 0, qn = 0;
    std::vector<uint32_t> lin, lout;
    std::unordered_set<uint64_t> seen;
    std::map<int, long long> hist;

    // ---- reading
    size_t fi = 0;
    std::unique_ptr<LineReader> rd;
    auto next_read = [&](PRead& r) -> bool {
        std::string_view h, s, p, q;
        for (;;) {
            if (!rd) {
                if (fi >= o.in.size()) return false;
                rd.reset(new LineReader(o.in[fi++]));
            }
            if (!rd->next(h)) { rd.reset(); continue; }
            if (h.empty()) continue;
            break;
        }
        if (h[0] != '@') die("not FASTQ: " + std::string(h.substr(0, 80)));
        std::string hs(h);
        if (!rd->next(s)) die("truncated FASTQ");
        std::string ss(s);
        if (!rd->next(p) || !rd->next(q)) die("truncated FASTQ");
        if (q.size() != ss.size()) die("FASTQ: sequence and quality lengths differ near read " + std::to_string(reads_in + 1));
        r.seq = std::move(ss);
        r.qual.assign(q);
        auto toks = split_ws(std::string_view(hs).substr(1));
        if (toks.empty()) die("FASTQ: empty read name");
        r.name.assign(toks[0]);
        r.zmw.clear();
        r.sub = false;
        for (size_t t = 0; t < toks.size() && t < 3; ++t) {
            std::string_view z;
            bool sb;
            if (pb_name(toks[t], z, sb)) { r.name.assign(toks[t]); r.zmw.assign(z); r.sub = sb; break; }
        }
        ++reads_in;
        bases_in += r.seq.size();
        lin.push_back((uint32_t)r.seq.size());
        if (r.zmw.empty()) ++n_nozmw;
        else if (r.sub) ++n_sub;
        else ++n_ccs;
        return true;
    };
    auto emit_group = [&](std::vector<PRead>& g, bool subs) {
        size_t n = g.size();
        hist[(int)std::min<size_t>(n, 11)]++;
        const PRead* p;
        if (n >= 3 && subs) {
            std::vector<const PRead*> inner;
            for (size_t i = 1; i + 1 < n; ++i) inner.push_back(&g[i]);
            std::stable_sort(inner.begin(), inner.end(), [](const PRead* a, const PRead* b) { return a->seq.size() < b->seq.size(); });
            p = inner[inner.size() / 2];
        } else {
            p = &g[0];
            for (size_t i = 1; i < n; ++i) if (g[i].seq.size() > p->seq.size()) p = &g[i];
        }
        fputc('@', out); fwrite(p->name.data(), 1, p->name.size(), out); fputc('\n', out);
        fwrite(p->seq.data(), 1, p->seq.size(), out); fputs("\n+\n", out);
        fwrite(p->qual.data(), 1, p->qual.size(), out); fputc('\n', out);
        lout.push_back((uint32_t)p->seq.size());
        bases_out += p->seq.size();
        ++groups_out;
        g.clear();
    };

    // ---- warm-up: decide how to group
    std::vector<PRead> warm;
    warm.reserve(o.warmup);
    {
        PRead r;
        while ((int)warm.size() < o.warmup && next_read(r)) warm.push_back(std::move(r));
    }
    for (auto& r : warm) for (char c : r.qual) { qsum += (unsigned char)c - 33; qn += 1; }
    const double mean_qv = qn > 0 ? (double)(qsum / qn) : 0.0;
    long long wz = 0;
    for (auto& r : warm) if (!r.zmw.empty()) ++wz;
    std::string mode = o.group;
    double thr = 0, p1w = 0, pLw = 0;
    long long masked_n = 0;
    if (mode == "auto") mode = (!warm.empty() && wz >= 0.9 * (double)warm.size()) ? "zmw" : "test";
    if (mode == "test" || mode == "adjacent") {
        // mask k-mers found in more than mask_frac of the warm-up reads (adapters, primers, poly(A), top transcripts)
        std::vector<uint16_t> cnt((size_t)1 << (2 * o.k), 0);
        std::vector<uint32_t> u;
#pragma omp parallel for schedule(dynamic, 16)
        for (long i = 0; i < (long)warm.size(); ++i) SK.sketch(warm[i].seq, warm[i].F, warm[i].R, false);
        for (auto& r : warm) {
            u.clear();
            for (auto& x : r.F) u.push_back(x.code);
            for (auto& x : r.R) u.push_back(x.code);
            std::sort(u.begin(), u.end());
            u.erase(std::unique(u.begin(), u.end()), u.end());
            for (uint32_t c : u) if (cnt[c] < 65535) ++cnt[c];
        }
        const double lim = std::max(5.0, o.mask_frac * (double)warm.size());
        for (size_t c = 0; c < cnt.size(); ++c) if (cnt[c] > lim) { SK.masked[c] = 1; ++masked_n; }
        long long m1 = 0, mL = 0, n1 = 0, nL = 0;
#pragma omp parallel for schedule(dynamic, 16)
        for (long i = 0; i < (long)warm.size(); ++i) SK.sketch(warm[i].seq, warm[i].F, warm[i].R, true);
#pragma omp parallel for schedule(dynamic, 64) reduction(+ : m1, mL, n1, nL)
        for (long i = 1; i < (long)warm.size(); ++i) {
            ++n1;
            if (diag_match(warm[i].F, warm[i].seq.size(), warm[i - 1].R, warm[i - 1].seq.size(), o.k, o).ok) ++m1;
            if (i >= (long)o.lag) {
                ++nL;
                if (diag_match(warm[i].F, warm[i].seq.size(), warm[i - o.lag].R, warm[i - o.lag].seq.size(), o.k, o).ok) ++mL;
            }
        }
        thr = 1;
        p1w = n1 ? (double)m1 / n1 : 0;
        pLw = nL ? (double)mL / nL : 0;
        if (mode == "test") mode = (p1w >= 0.05 && p1w > 5 * pLw) ? "adjacent" : "none";
    }
    logm(sfmt("[lrpick] %s: warm-up %zu reads, %.0f%% with ZMW names, mean QV %.1f -> grouping: %s", o.label.c_str(),
              warm.size(), warm.empty() ? 0.0 : 100.0 * wz / warm.size(), mean_qv, mode.c_str()));
    if (thr > 0)
        logm(sfmt("[lrpick]   reverse-complement test: consecutive reads = same molecule %.2f%%, "
                  "lag-%d pairs (chance) %.3f%%, masked k-mers %lld", 100 * p1w, o.lag, 100 * pLw, masked_n));
    if (o.decide_only) {
        std::string kind = mode == "zmw" ? (n_sub > std::max(n_ccs, n_nozmw) ? "subreads" : (n_ccs > n_nozmw ? "hifi" : "unknown"))
                         : mode == "adjacent" ? "subreads" : (mean_qv >= 20 ? "hifi" : "unknown");
        FILE* f = fopen(o.stats.c_str(), "w");
        if (!f) die("cannot write " + o.stats);
        std::fprintf(f, "file\t%s\nkind\t%s\npreset\t%s\ngrouping\t%s\nwarmup_reads\t%zu\nmean_qv_warmup\t%.2f\n"
                        "names_with_zmw\t%lld\n", o.label.c_str(), kind.c_str(), kind == "hifi" ? "map-hifi" : "map-pb",
                     mode.c_str(), warm.size(), mean_qv, wz);
        if (thr > 0) std::fprintf(f, "warmup_lag1_match\t%.5f\nwarmup_lagL_match\t%.6f\nmasked_kmers\t%lld\n", p1w, pLw, masked_n);
        fclose(f);
        fclose(out);
        return 0;
    }

    // ---- stream: warm-up reads first, then the rest
    size_t wi = 0;
    const bool sketched_warm = mode == "adjacent";      // warm-up reads already carry masked sketches
    auto fill = [&](std::vector<PRead>& B, size_t n, std::vector<char>& need_sk) -> bool {
        B.clear();
        need_sk.clear();
        while (B.size() < n) {
            PRead r;
            bool from_warm = wi < warm.size();
            if (from_warm) r = std::move(warm[wi++]);
            else {
                if (!warm.empty()) std::vector<PRead>().swap(warm);
                if (!next_read(r)) break;
            }
            need_sk.push_back(!(from_warm && sketched_warm));
            B.push_back(std::move(r));
        }
        return !B.empty();
    };
    FILE* dbg = o.debug.empty() ? nullptr : fopen(o.debug.c_str(), "w");
    if (dbg) std::fputs("prev\tcur\tlen_prev\tlen_cur\tkmers_prev\tkmers_cur\tshared\tband\toverlap_frac\tcover\tjoined\n", dbg);
    std::vector<PRead> grp, B;
    std::vector<char> need_sk;
    long long pairs1 = 0, match1 = 0, pairsL = 0, matchL = 0;
    std::vector<std::vector<Km>> ringR((size_t)o.lag);
    std::vector<size_t> ringL((size_t)o.lag, 0);
    long long idx = 0;                       // global index of the first read of the batch (adjacent mode)
    std::string curz;
    // previous read (last of the previous batch), adjacent mode
    bool have_prev = false;
    std::vector<Km> prevR;
    size_t prevL = 0;
    std::string prevName;
    const size_t BATCH = mode == "adjacent" ? 8192 : 1;
    std::vector<DiagHit> jh;
    std::vector<char> lh;
    while (fill(B, BATCH, need_sk)) {
        if (mode == "adjacent") {
            const long n = (long)B.size();
#pragma omp parallel for schedule(dynamic, 16)
            for (long i = 0; i < n; ++i)
                if (need_sk[i]) SK.sketch(B[i].seq, B[i].F, B[i].R, true);
            jh.assign(n, DiagHit());
            lh.assign(n, 0);
#pragma omp parallel for schedule(dynamic, 32)
            for (long i = 0; i < n; ++i) {
                if (i > 0) jh[i] = diag_match(B[i].F, B[i].seq.size(), B[i - 1].R, B[i - 1].seq.size(), o.k, o);
                else if (have_prev) jh[i] = diag_match(B[i].F, B[i].seq.size(), prevR, prevL, o.k, o);
                const long long g = idx + i;
                if (g >= o.lag) {
                    bool ok = i >= (long)o.lag ? diag_match(B[i].F, B[i].seq.size(), B[i - o.lag].R, B[i - o.lag].seq.size(), o.k, o).ok
                                               : diag_match(B[i].F, B[i].seq.size(), ringR[(size_t)(g % o.lag)], ringL[(size_t)(g % o.lag)], o.k, o).ok;
                    lh[i] = ok ? 2 : 1;
                }
            }
            // ring: the last `lag` reads of this batch (global index g -> slot g % lag)
            for (long i = std::max(0L, n - (long)o.lag); i < n; ++i) {
                const size_t sl = (size_t)((idx + i) % o.lag);
                ringR[sl] = B[i].R;
                ringL[sl] = B[i].seq.size();
            }
            for (long i = 0; i < n; ++i) {
                if (lh[i]) { ++pairsL; if (lh[i] == 2) ++matchL; }
                bool join = false;
                if (i > 0 || have_prev) {
                    ++pairs1;
                    join = jh[i].ok && !grp.empty() && (int)grp.size() < o.max_group;
                    if (join) ++match1;
                    if (dbg) {
                        const std::string& pn = i > 0 ? B[i - 1].name : prevName;
                        size_t pl = i > 0 ? B[i - 1].seq.size() : prevL, pk = i > 0 ? B[i - 1].R.size() : prevR.size();
                        std::fprintf(dbg, "%s\t%s\t%zu\t%zu\t%zu\t%zu\t%d\t%d\t%.3f\t%.3f\t%d\n", pn.c_str(), B[i].name.c_str(), pl,
                                     B[i].seq.size(), pk, B[i].F.size(), jh[i].shared, jh[i].band, jh[i].ovf, jh[i].cov, (int)join);
                    }
                }
                if (!join && !grp.empty()) emit_group(grp, true);
                if (i == n - 1) { prevR = B[i].R; prevL = B[i].seq.size(); prevName = B[i].name; have_prev = true; }
                grp.push_back(std::move(B[i]));
                if (!grp.back().F.empty() || !grp.back().R.empty()) { std::vector<Km>().swap(grp.back().F); }
            }
            idx += n;
            continue;
        }
        for (auto& cur : B) {
            if (mode == "zmw") {
                if (cur.zmw.empty()) {
                    if (!grp.empty()) emit_group(grp, grp[0].sub);
                    curz.clear();
                    grp.push_back(std::move(cur));
                    emit_group(grp, false);
                    continue;
                }
                if (cur.zmw != curz) {
                    if (!grp.empty()) emit_group(grp, grp[0].sub);
                    uint64_t hz = std::hash<std::string>()(cur.zmw);
                    if (!seen.insert(hz).second) ++reappeared;
                    curz = cur.zmw;
                }
                grp.push_back(std::move(cur));
            } else {
                grp.push_back(std::move(cur));
                emit_group(grp, false);
            }
        }
    }
    if (!grp.empty()) emit_group(grp, mode == "adjacent" || grp[0].sub);
    if (fflush(out) != 0 || fclose(out) != 0) die("write error on " + o.out);
    if (dbg) fclose(dbg);

    std::string kind;
    if (mode == "zmw") kind = n_sub > std::max(n_ccs, n_nozmw) ? "subreads" : (n_ccs > n_nozmw ? "hifi" : "unknown");
    else if (mode == "adjacent") kind = "subreads";
    else kind = mean_qv >= 20 ? "hifi" : "unknown";
    const std::string preset = kind == "hifi" ? "map-hifi" : "map-pb";
    FILE* f = fopen(o.stats.c_str(), "w");
    if (!f) die("cannot write " + o.stats);
    auto row = [&](const char* k, const std::string& v) { std::fprintf(f, "%s\t%s\n", k, v.c_str()); };
    auto ll = [](long double x) { return std::to_string((long long)x); };
    row("file", o.label); row("kind", kind); row("preset", preset); row("grouping", mode);
    row("reads_in", ll(reads_in)); row("bases_in", ll(bases_in)); row("n50_in", ll(n50(lin)));
    row("mean_qv_warmup", sfmt("%.2f", mean_qv));
    row("subread_names", ll(n_sub)); row("ccs_names", ll(n_ccs)); row("names_without_zmw", ll(n_nozmw));
    if (mode == "zmw") { row("zmws", ll((long double)seen.size())); row("zmw_reappeared_non_adjacent", ll(reappeared)); }
    if (thr > 0) {
        row("masked_kmers", ll(masked_n));
        row("warmup_lag1_match", sfmt("%.5f", p1w)); row("warmup_lagL_match", sfmt("%.6f", pLw));
    }
    if (mode == "adjacent") {
        row("lag1_match_rate", sfmt("%.5f", pairs1 ? (double)match1 / pairs1 : 0.0));
        row("lagL_match_rate_chance", sfmt("%.6f", pairsL ? (double)matchL / pairsL : 0.0));
    }
    row("reads_out", ll(groups_out)); row("bases_out", ll(bases_out)); row("n50_out", ll(n50(lout)));
    row("reads_per_molecule", sfmt("%.3f", groups_out ? (double)reads_in / groups_out : 0.0));
    for (auto& kv : hist) row(kv.first < 11 ? sfmt("molecules_with_%d_reads", kv.first).c_str() : "molecules_with_11+_reads", ll(kv.second));
    fclose(f);
    logm(sfmt("[lrpick] %s: %s, %s reads -> %s molecules (%.2f reads per molecule), preset %s, %.0fs", o.label.c_str(),
              kind.c_str(), commas(reads_in).c_str(), commas(groups_out).c_str(),
              groups_out ? (double)reads_in / groups_out : 0.0, preset.c_str(), wall() - T0));
    return 0;
}

// ======================================================================================== lrcount
struct CountOpts {
    std::string bam, sam, out = "LongRead_Counts.tsv", legacy, summary, tmpdir, zmw = "auto";
    double as_tol = 0.0, min_id = 0.70;
    int min_len = 150, threads = 4;
};

void count_usage() {
    std::fprintf(stderr,
        "usage: lqcore lrcount (--bam BAM | --sam SAM|-) [-o LongRead_Counts.tsv] [--legacy-out F] [--summary F]\n"
        "        [--as-tol 0] [--min-aln-len 150] [--min-identity 0.70] [--zmw-dedup auto|on|off] [--threads 4]\n"
        "  molecules per allele copy (C++ port of lineagequant/longread/extract_prior.py, same output)\n");
}

// child process "sh -c CMD" whose stdout we read. Unlike popen(): its own process group with SIGPIPE at the
// default action, so stopping early (closing our end, then SIGTERM to the group) always ends the whole pipeline.
// (lqcore 0.5.0 stopped `samtools view` early with pclose(); in the PBS job that left lrcount waiting forever.)
struct Proc {
    pid_t pid = -1;
    int fd = -1;
    std::string cmd;
    void start(const std::string& c) {
        cmd = c;
        int p[2];
        if (pipe(p) != 0) die("pipe() failed");
        pid = fork();
        if (pid < 0) die("fork() failed");
        if (pid == 0) {
            setpgid(0, 0);
            signal(SIGPIPE, SIG_DFL);
            dup2(p[1], 1);
            close(p[0]);
            close(p[1]);
            execl("/bin/sh", "sh", "-c", c.c_str(), (char*)nullptr);
            _exit(127);
        }
        setpgid(pid, pid);
        close(p[1]);
        fd = p[0];
        fcntl(fd, F_SETFD, FD_CLOEXEC);
    }
    // returns the exit status (0 = ok); early = we stopped reading on purpose -> terminate it, status ignored
    int finish(bool early) {
        if (pid <= 0) return 0;
        if (early) kill(-pid, SIGTERM);
        int st = 0;
        while (waitpid(pid, &st, 0) < 0 && errno == EINTR) {}
        pid = -1;
        if (early) return 0;
        return WIFEXITED(st) ? WEXITSTATUS(st) : 128 + (WIFSIGNALED(st) ? WTERMSIG(st) : 0);
    }
};

struct Input {
    Proc proc;
    std::unique_ptr<LineReader> rd;
    ~Input() { stop(); }
    void close() {           // read to the end
        rd.reset();
        if (proc.pid > 0) { int st = proc.finish(false); if (st != 0) die("input pipeline failed (exit " + std::to_string(st) + "): " + proc.cmd); }
    }
    void stop() {            // stopped early on purpose
        rd.reset();
        proc.finish(true);
    }
};

std::string g_collate_prefix;
bool header_grouped(const std::vector<std::string>& hdr) {
    for (auto& l : hdr)
        if (l.rfind("@HD", 0) == 0)
            for (auto f : split_tab(l))
                if (f == "SO:queryname" || f == "GO:query") return true;
    return false;
}
void open_in(const CountOpts& o, Input& in, bool collate) {
    if (!o.sam.empty() && !collate) { in.rd.reset(new LineReader(o.sam)); return; }
    const std::string src = !o.bam.empty() ? o.bam : o.sam;
    std::string th = std::to_string(std::max(1, o.threads));
    std::string cmd = collate ? "samtools collate -O -u -@ " + th + " " + sh_quote(src) + " " + sh_quote(g_collate_prefix) +
                                    " | samtools view -h -@ " + th + " -"
                              : "samtools view -h -@ " + th + " " + sh_quote(src);
    in.proc.start(cmd);
    in.rd.reset(new LineReader(in.proc.fd));     // LineReader owns (and closes) the fd
    in.proc.fd = -1;
}
// reads the header; returns the first record line (or empty)
bool read_header(LineReader& rd, std::vector<std::string>& hdr, std::string& first) {
    std::string_view l;
    while (rd.next(l)) {
        if (!l.empty() && l[0] == '@') { hdr.emplace_back(l); continue; }
        if (l.empty()) continue;
        first.assign(l);
        return true;
    }
    return false;
}

struct LRec {
    bool sec;
    long long qs, qe, as, qa;
    double ident;
    int rid;
};

// SAM line -> fields used here (no allocation)
struct SamL {
    std::string_view qname, rname, cigar;
    int flag = 0;
    bool has_as = false, has_nm = false;
    long long as = 0, nm = 0;
};
bool parse_line(std::string_view l, SamL& r) {
    size_t st[12];
    int nf = 0;
    size_t p = 0;
    st[nf++] = 0;
    while (nf < 12) {
        size_t q = l.find('\t', p);
        if (q == NPOS) break;
        st[nf++] = q + 1;
        p = q + 1;
    }
    if (nf < 11) return false;
    auto fld = [&](int i) { size_t b = st[i]; size_t e = l.find('\t', b); return l.substr(b, (e == NPOS ? l.size() : e) - b); };
    r.qname = fld(0);
    r.flag = (int)parse_ll(fld(1));
    r.rname = fld(2);
    r.cigar = fld(5);
    r.has_as = r.has_nm = false;
    if (nf == 12) {
        size_t b = st[11];
        while (b <= l.size()) {
            size_t e = l.find('\t', b);
            std::string_view tg = l.substr(b, (e == NPOS ? l.size() : e) - b);
            if (tg.size() > 5 && tg[2] == ':' && tg[3] == 'i' && tg[4] == ':') {
                if (tg[0] == 'A' && tg[1] == 'S' && !r.has_as) { r.has_as = true; r.as = parse_ll(tg.substr(5)); }
                else if (tg[0] == 'N' && tg[1] == 'M' && !r.has_nm) { r.has_nm = true; r.nm = parse_ll(tg.substr(5)); }
            }
            if (e == NPOS) break;
            b = e + 1;
        }
    }
    return true;
}
LRec rec_info(const SamL& r, int rid) {
    // cigar ops
    std::vector<std::pair<char, long long>> cg;
    long long n = 0;
    for (char c : r.cigar) {
        if (c >= '0' && c <= '9') n = n * 10 + (c - '0');
        else { cg.emplace_back(c, n); n = 0; }
    }
    long long lead = 0, trail = 0, qa = 0, cols = 0;
    long long i = 0, j = (long long)cg.size() - 1;
    auto clip = [](char c) { return c == 'S' || c == 'H'; };
    while (i <= j && clip(cg[i].first)) { lead += cg[i].second; ++i; }
    while (j >= i && clip(cg[j].first)) { trail += cg[j].second; --j; }
    for (long long t = i; t <= j; ++t) {
        char op = cg[t].first;
        if (op == 'M' || op == 'I' || op == '=' || op == 'X') qa += cg[t].second;
        if (op == 'M' || op == 'I' || op == 'D' || op == '=' || op == 'X') cols += cg[t].second;
    }
    LRec x;
    x.sec = r.flag & 0x100;
    x.qs = (r.flag & 0x10) ? trail : lead;
    x.qe = x.qs + qa;
    x.qa = qa;
    x.ident = (r.has_nm && cols) ? 1.0 - (double)r.nm / (double)cols : 1.0;
    x.as = r.has_as ? r.as : qa - 2 * (r.has_nm ? r.nm : 0);
    x.rid = rid;
    return x;
}

// insertion-ordered string -> T map
template <class T>
struct OMap {
    std::unordered_map<std::string, size_t> ix;
    std::vector<std::string> keys;
    std::vector<T> vals;
    T& at(const std::string& k) {
        auto it = ix.find(k);
        if (it != ix.end()) return vals[it->second];
        ix.emplace(k, keys.size());
        keys.push_back(k);
        vals.emplace_back();
        return vals.back();
    }
};

int lrcount_main(int argc, char** argv) {
    CountOpts o;
    auto need = [&](int& i) -> std::string { if (i + 1 >= argc) { count_usage(); die(std::string("missing value for ") + argv[i]); } return argv[++i]; };
    for (int i = 0; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "-b" || a == "--bam") o.bam = need(i);
        else if (a == "--sam") o.sam = need(i);
        else if (a == "-o" || a == "--out") o.out = need(i);
        else if (a == "--legacy-out") o.legacy = need(i);
        else if (a == "--summary") o.summary = need(i);
        else if (a == "--as-tol") { if (!py_float(need(i), o.as_tol)) die("--as-tol: not a number"); }
        else if (a == "--min-aln-len") o.min_len = std::atoi(need(i).c_str());
        else if (a == "--min-identity") { if (!py_float(need(i), o.min_id)) die("--min-identity: not a number"); }
        else if (a == "--zmw-dedup") o.zmw = need(i);
        else if (a == "--threads") o.threads = std::max(1, std::atoi(need(i).c_str()));
        else if (a == "--tmpdir") o.tmpdir = need(i);
        else if (a == "-h" || a == "--help") { count_usage(); return 0; }
        else { count_usage(); die("unknown option " + a); }
    }
    if (o.bam.empty() == o.sam.empty()) { count_usage(); die("give exactly one of --bam / --sam"); }
    if (o.zmw != "auto" && o.zmw != "on" && o.zmw != "off") die("--zmw-dedup: auto|on|off");
    const double T0 = wall();
    const std::string src = !o.bam.empty() ? o.bam : o.sam;
    if (src != "-" && !file_exists(src)) die("not found: " + src);
    logm("[INFO] long-read alignments: " + src);
    {
        std::string dir = !o.tmpdir.empty() ? o.tmpdir : (o.out.find('/') == std::string::npos ? "." : o.out.substr(0, o.out.rfind('/')));
        const char* job = getenv("PBS_JOBID");
        g_collate_prefix = dir + "/.lqcore_lr_collate_" + (job ? std::string(job) : std::to_string(getpid()));
    }

    // ---- header and grouping
    std::vector<std::string> hdr;
    std::string first;
    bool collate = false;
    std::unique_ptr<LineReader> stdin_rd;       // SAM on stdin: one reader for header and records
    if (src == "-") {
        stdin_rd.reset(new LineReader("-"));
        read_header(*stdin_rd, hdr, first);
    } else if (!o.sam.empty()) {
        LineReader rd(o.sam);
        read_header(rd, hdr, first);
    } else {                                      // BAM: the header alone, read to its end
        Input in;
        in.proc.start("samtools view -H " + sh_quote(src));
        in.rd.reset(new LineReader(in.proc.fd));
        in.proc.fd = -1;
        std::string_view l;
        while (in.rd->next(l)) if (!l.empty()) hdr.emplace_back(l);
        in.close();
    }
    if (!header_grouped(hdr)) {
        if (src == "-") die("SAM on stdin must be grouped by read name (@HD GO:query or SO:queryname)");
        collate = true;
        logm("[INFO] not grouped by read name -> samtools collate (streamed)");
    }
    logm("[INFO] header: " + commas((long long)hdr.size()) + " lines" + (collate ? "" : ", grouped by read name"));
    std::vector<std::string> refs;
    for (auto& l : hdr)
        if (l.rfind("@SQ", 0) == 0)
            for (auto f : split_tab(l))
                if (f.rfind("SN:", 0) == 0) refs.emplace_back(f.substr(3));
    std::unordered_map<std::string, int> r2i;
    r2i.reserve(refs.size() * 2);
    for (size_t i = 0; i < refs.size(); ++i) r2i.emplace(refs[i], (int)i);
    // allele key of each reference (ID before the first dot)
    std::vector<std::string> rkey(refs.size());
    for (size_t i = 0; i < refs.size(); ++i) { size_t d = refs[i].find('.'); rkey[i] = refs[i].substr(0, d); }

    // ---- ZMW de-duplication (subreads of one ZMW = one molecule)
    std::string zmode;
    std::unordered_map<std::string, int> per;   // ZMW (or read name) -> subreads
    bool use_w = false;
    if (o.zmw == "off") zmode = "off";
    else if (src == "-") { zmode = o.zmw + " -> off (stdin)"; }
    else {
        long long n_names = 0, n_pb = 0;
        std::unordered_map<std::string, std::unordered_set<std::string>> smp;
        {
            Input in;
            open_in(o, in, collate);
            std::vector<std::string> h2;
            std::string f2;
            bool have = read_header(*in.rd, h2, f2);
            std::string_view l = f2;
            SamL r;
            while (have && n_names < 200000) {
                if (parse_line(l, r) && !(r.flag & 0x900)) {
                    ++n_names;
                    std::string_view z;
                    bool sb;
                    if (pb_name(r.qname, z, sb) && sb) { ++n_pb; smp[std::string(z)].insert(std::string(r.qname)); }
                }
                if (!in.rd->next(l)) break;
            }
            in.stop();
        }
        bool looks_sub = n_names > 0 && (double)n_pb >= 0.9 * (double)n_names;
        bool multi = false;
        for (auto& kv : smp) if (kv.second.size() > 1) { multi = true; break; }
        if (o.zmw == "auto" && !(looks_sub && multi)) zmode = std::string("auto -> off (") + (!looks_sub ? "not PacBio subread names" : "one subread per ZMW") + ")";
        else if (o.zmw == "on" && n_pb == 0) { logm("[WARN] --zmw-dedup on, but read names carry no ZMW -> off"); zmode = "on -> off (no ZMW in read names)"; }
        else {
            logm("[INFO] subreads of one ZMW = one molecule: counting subreads per ZMW...");
            std::unordered_set<std::string> seen;
            Input in;
            open_in(o, in, collate);
            std::vector<std::string> h2;
            std::string f2;
            bool have = read_header(*in.rd, h2, f2);
            std::string_view l = f2;
            SamL r;
            while (have) {
                if (parse_line(l, r) && !(r.flag & 0x900)) {
                    std::string q(r.qname);
                    if (seen.insert(q).second) {
                        std::string_view z;
                        bool sb;
                        per[pb_name(q, z, sb) ? std::string(z) : q] += 1;
                    }
                }
                if (!in.rd->next(l)) break;
            }
            in.close();
            use_w = true;
            logm(sfmt("[INFO] %s subreads in %s ZMWs (%.2f per ZMW)", commas((long long)seen.size()).c_str(),
                      commas((long long)per.size()).c_str(), (double)seen.size() / std::max<size_t>(per.size(), 1)));
            zmode = o.zmw + " -> on (" + commas((long long)seen.size()) + " subreads / " + commas((long long)per.size()) + " ZMWs)";
        }
    }
    logm("[INFO] ZMW de-duplication: " + zmode);

    // ---- main pass
    std::vector<long long> legacy_ref(refs.size(), 0);
    double st_reads = 0, st_unmapped = 0, st_prim = 0, st_supp = 0, st_sec = 0, st_qcdup = 0;
    double st_seg = 0, st_failed = 0, st_used = 0, st_tied = 0, st_mol = 0;
    OMap<double> ecs;                       // key: allele ids joined (sorted by allele name)
    std::unordered_map<std::string, int> aid;
    std::vector<std::string> anames;
    std::vector<int> rk_id(refs.size(), -1);
    auto allele_of = [&](int rid) -> int {
        int& v = rk_id[rid];
        if (v < 0) {
            auto it = aid.find(rkey[rid]);
            if (it == aid.end()) { it = aid.emplace(rkey[rid], (int)anames.size()).first; anames.push_back(rkey[rid]); }
            v = it->second;
        }
        return v;
    };
    std::vector<std::vector<int>> ec_items;   // parallel to ecs.keys

    std::vector<LRec> recs;
    std::string curq;
    bool have_cur = false, cur_prim = false;
    auto flush = [&]() {
        st_reads += 1;
        if (recs.empty()) { st_unmapped += 1; return; }
        double w = 1.0;
        if (use_w && cur_prim) {
            std::string_view z;
            bool sb;
            auto it = per.find(pb_name(curq, z, sb) ? std::string(z) : curq);
            if (it != per.end()) w = 1.0 / it->second;
        }
        // segments: primary/supplementary sorted by query start (stable); secondaries join the most-overlapping one
        std::vector<std::vector<const LRec*>> segs;
        {
            std::vector<const LRec*> pr;
            for (auto& x : recs) if (!x.sec) pr.push_back(&x);
            std::stable_sort(pr.begin(), pr.end(), [](const LRec* a, const LRec* b) { return a->qs < b->qs; });
            for (auto* x : pr) segs.push_back({x});
        }
        for (auto& x : recs) {
            if (!x.sec) continue;
            long long bo = 0;
            int best = -1;
            for (size_t s = 0; s < segs.size(); ++s) {
                const LRec* p = segs[s][0];
                long long ov = std::min(p->qe, x.qe) - std::max(p->qs, x.qs);
                if (ov > bo) { best = (int)s; bo = ov; }
            }
            if (best >= 0 && (double)bo >= 0.5 * (double)std::max<long long>(1, std::min(x.qe - x.qs, segs[best][0]->qe - segs[best][0]->qs)))
                segs[best].push_back(&x);
            else segs.push_back({&x});
        }
        for (auto& s : segs) {
            st_seg += 1;
            long long top = s[0]->as;
            for (auto* x : s) top = std::max(top, x->as);
            std::vector<const LRec*> cand;
            for (auto* x : s) if ((double)x->as >= (double)top - o.as_tol) cand.push_back(x);
            const LRec* b = cand[0];
            for (auto* x : cand) if (x->as > b->as || (x->as == b->as && x->qa > b->qa)) b = x;
            if (b->qa < o.min_len || b->ident < o.min_id) { st_failed += 1; continue; }
            std::vector<int> ks;
            for (auto* x : cand) ks.push_back(allele_of(x->rid));
            std::sort(ks.begin(), ks.end(), [&](int a, int c) { return anames[a] < anames[c]; });
            ks.erase(std::unique(ks.begin(), ks.end()), ks.end());
            std::string key;
            for (int a : ks) { key.append((const char*)&a, sizeof a); }
            st_used += 1;
            st_mol += w;
            if (ks.size() > 1) st_tied += 1;
            size_t before = ecs.keys.size();
            ecs.at(key) += w;
            if (ecs.keys.size() > before) ec_items.push_back(ks);
        }
    };
    {
        Input in;
        std::vector<std::string> h2;
        std::string f2;
        bool have;
        LineReader* rdp;
        if (stdin_rd) { rdp = stdin_rd.get(); f2 = first; have = !first.empty(); }
        else { open_in(o, in, collate); rdp = in.rd.get(); have = read_header(*rdp, h2, f2); }
        std::string_view l = f2;
        SamL r;
        double ts = wall();
        long long nl = 0;
        while (have) {
            if ((++nl & ((1LL << 25) - 1)) == 0)        // progress every ~33.5 M records
                logm(sfmt("[INFO]   ... %s records, %s reads, %s equivalence classes (%.0fs)", commas(nl).c_str(),
                          commas((long long)st_reads).c_str(), commas((long long)ecs.keys.size()).c_str(), wall() - ts));
            if (!parse_line(l, r)) die("bad SAM line: " + std::string(l.substr(0, 120)));
            if (!have_cur || r.qname != curq) {
                if (have_cur) flush();
                curq.assign(r.qname);
                have_cur = true;
                recs.clear();
                cur_prim = false;
            }
            const int f = r.flag;
            if (!(f & 0x900)) cur_prim = true;
            if (!(f & 0x4)) {
                auto it = r2i.find(std::string(r.rname));
                if (it == r2i.end()) die("reference not in the header: " + std::string(r.rname));
                legacy_ref[it->second] += 1;
                if (f & 0x200 || f & 0x400) st_qcdup += 1;
                else {
                    if (f & 0x100) st_sec += 1;
                    else if (f & 0x800) st_supp += 1;
                    else st_prim += 1;
                    recs.push_back(rec_info(r, it->second));
                }
            }
            if (!rdp->next(l)) break;
        }
        if (have_cur) flush();
        in.close();
        logm(sfmt("[INFO] %s reads, %s alignment records (%.0fs)", commas((long long)st_reads).c_str(),
                  commas((long long)(st_prim + st_supp + st_sec)).c_str(), wall() - ts));
    }

    // ---- EM over equivalence classes (same arithmetic as extract_prior.run_em)
    std::vector<int> all;
    {
        std::vector<char> in_ec(anames.size(), 0);
        for (auto& ks : ec_items) for (int a : ks) in_ec[a] = 1;
        for (size_t a = 0; a < anames.size(); ++a) if (in_ec[a]) all.push_back((int)a);
        std::sort(all.begin(), all.end(), [&](int a, int c) { return anames[a] < anames[c]; });
    }
    const size_t K = all.size();
    std::vector<int> pos(anames.size(), -1);
    for (size_t i = 0; i < K; ++i) pos[all[i]] = (int)i;
    std::vector<double> u(K, 0.0), n(K, 0.0);
    std::vector<int> lens, items, ecid;
    std::vector<double> ew;
    for (size_t e = 0; e < ec_items.size(); ++e) {
        auto& ks = ec_items[e];
        if (ks.size() == 1) { u[pos[ks[0]]] += 0.0 + ecs.vals[e]; continue; }
        int m = (int)lens.size();
        lens.push_back((int)ks.size());
        ew.push_back(ecs.vals[e]);
        for (int a : ks) { items.push_back(pos[a]); ecid.push_back(m); }
    }
    long long iters = 0;
    logm(sfmt("[INFO] EM: %s allele copies, %s tied classes (%s entries), %s single-copy classes", commas((long long)K).c_str(),
              commas((long long)lens.size()).c_str(), commas((long long)items.size()).c_str(),
              commas((long long)(ec_items.size() - lens.size())).c_str()));
    const double tem = wall();
    if (!lens.empty()) {
        const size_t I = items.size(), M = lens.size();
        std::vector<double> acc(K, 0.0), a(I), den(M), nw(K);
        for (size_t i = 0; i < I; ++i) acc[items[i]] += ew[ecid[i]] / (double)lens[ecid[i]];
        for (size_t k2 = 0; k2 < K; ++k2) n[k2] = u[k2] + acc[k2];
        for (long long it = 1; it <= 5000; ++it) {
            iters = it;
            for (size_t i = 0; i < I; ++i) a[i] = n[items[i]];
            std::fill(den.begin(), den.end(), 0.0);
            for (size_t i = 0; i < I; ++i) den[ecid[i]] += a[i];
            std::fill(acc.begin(), acc.end(), 0.0);
            for (size_t i = 0; i < I; ++i) {
                double d = den[ecid[i]];
                double r = d > 0 ? a[i] / d : 1.0 / (double)lens[ecid[i]];
                acc[items[i]] += r * ew[ecid[i]];
            }
            double dmax = 0;
            for (size_t k2 = 0; k2 < K; ++k2) {
                nw[k2] = u[k2] + acc[k2];
                double d = std::fabs(nw[k2] - n[k2]);
                if (d > dmax || std::isnan(d)) dmax = d;
            }
            n.swap(nw);
            if (it % 500 == 0) logm(sfmt("[INFO]   EM iteration %lld: max change %.3g (%.0fs)", it, dmax, wall() - tem));
            if (dmax < 1e-6) break;
        }
        logm(sfmt("[INFO] EM: %lld iterations (%.0fs)", iters, wall() - tem));
    } else n = u;

    // ---- outputs
    long long n_al = 0;
    {
        std::string tmp = o.out + ".tmp";
        FILE* f = fopen(tmp.c_str(), "w");
        if (!f) die("cannot write " + o.out);
        std::fputs("AlleleID\tLongRead_Count\n", f);
        for (size_t i = 0; i < K; ++i) {
            double c = py_round(n[i], 4);
            if (c > 0) { std::fprintf(f, "%s\t%s\n", anames[all[i]].c_str(), py_repr(c).c_str()); ++n_al; }
        }
        if (fclose(f) != 0 || rename(tmp.c_str(), o.out.c_str()) != 0) die("cannot write " + o.out);
    }
    OMap<long long> legacy;
    for (size_t i = 0; i < refs.size(); ++i) legacy.at(rkey[i]) = legacy_ref[i];
    long long leg_total = 0, leg_coll = 0;
    for (auto v : legacy_ref) leg_total += v;
    for (auto v : legacy.vals) leg_coll += v;
    if (!o.legacy.empty()) {
        FILE* f = fopen(o.legacy.c_str(), "w");
        if (!f) die("cannot write " + o.legacy);
        std::fputs("AlleleID\tLongRead_Count\n", f);
        for (size_t i = 0; i < legacy.keys.size(); ++i) std::fprintf(f, "%s\t%lld\n", legacy.keys[i].c_str(), legacy.vals[i]);
        fclose(f);
    }
    std::vector<std::pair<std::string, std::string>> rows = {
        {"reads", std::to_string((long long)st_reads)}, {"reads_without_alignment", std::to_string((long long)st_unmapped)},
        {"records_primary", std::to_string((long long)st_prim)}, {"records_supplementary", std::to_string((long long)st_supp)},
        {"records_secondary", std::to_string((long long)st_sec)}, {"records_qcfail_dup", std::to_string((long long)st_qcdup)},
        {"zmw_dedup", zmode},
        {"segments", std::to_string((long long)st_seg)}, {"segments_failed_qc", std::to_string((long long)st_failed)},
        {"segments_used", std::to_string((long long)st_used)}, {"segments_tied_across_copies", std::to_string((long long)st_tied)},
        {"equivalence_classes", std::to_string((long long)ecs.keys.size())}, {"em_iterations", std::to_string(iters)},
        {"molecules_counted", py_repr(py_round(st_mol, 2))},
        {"alleles_with_reads", std::to_string(n_al)},
        {"legacy_idxstats_total", std::to_string(leg_total)},
        {"legacy_collapsed_total", std::to_string(leg_coll)},
    };
    for (auto& kv : rows) {
        bool num = !kv.second.empty() && all_digits(kv.second);
        std::string v = num ? commas(std::atoll(kv.second.c_str())) : kv.second;
        std::printf("  %-32s%18s\n", kv.first.c_str(), v.c_str());
    }
    if (st_mol > 0) std::printf("  1.0.0 idxstats total / molecules = %.2fx\n", (double)leg_total / st_mol);
    std::fflush(stdout);
    if (!o.summary.empty()) {
        FILE* f = fopen(o.summary.c_str(), "w");
        if (!f) die("cannot write " + o.summary);
        std::fputs("metric\tvalue\n", f);
        for (auto& kv : rows) std::fprintf(f, "%s\t%s\n", kv.first.c_str(), kv.second.c_str());
        fclose(f);
    }
    logm(sfmt("[SUCCESS] %s allele copies with long reads -> %s (%.0fs)", commas(n_al).c_str(), o.out.c_str(), wall() - T0));
    return 0;
}

}  // namespace

int lrpick_entry(int argc, char** argv) { return lrpick_main(argc, argv); }
int lrcount_entry(int argc, char** argv) { return lrcount_main(argc, argv); }

}  // namespace lq
