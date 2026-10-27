// lqcore theta -- port of lineagequant/theta/global_theta.py (Step 01: per-copy expression prior)
//
//   Theta k-mers = k-mers of the target transcripts that belong to exactly one copy (gene id = transcript id
//   up to the first "."). Every read window is looked up in both orientations; each hit adds 1 to its copy,
//   and (paired mode) a copy hit by both mates gets one bonus count. Theta = shrunk hits/probes, normalised
//   within the cluster; clusters with enough long reads use the long-read counts instead.
//
// Output is byte-identical to global_theta.py: same counts, same floating-point operations in the same order,
// same formatting. Python's built-in sum() of floats is compensated (Neumaier) from Python 3.12 on and a plain
// left fold before; --py-sum selects which one to reproduce (the router passes the running interpreter's).
//
// Speed: the index is a sorted array of 2-bit k-mer codes (bucketed by the leading bits), built in parallel
// and cached on disk; the FASTQ pair is decompressed by two reader threads and scanned by all others.
#include "binio.hpp"
#include "common.hpp"

#include <atomic>
#include <cmath>
#include <condition_variable>
#include <deque>
#include <map>
#include <mutex>
#include <parallel/algorithm>
#include <set>
#include <sys/stat.h>
#include <thread>
#include <unordered_map>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace lq {
namespace {

struct TOpts {
    std::string clusters, fasta, r1, r2, lr_prior, outdir = "out", index_dir, index_cache, hits_cache, py_sum = "neumaier";
    int k = 21, threads = 16, min_lr_depth = 5;
    double shrinkage = 10.0;
    bool adaptive = false, no_bonus = false, binary = false;
};

void theta_usage() {
    std::fprintf(stderr,
                 "usage: lqcore theta -c CLUSTERS -f FASTA -r1 R1 -r2 R2 [-k 21] [-t THREADS] [-o OUTDIR]\n"
                 "        [--shrinkage 10.0] [--adaptive_shrinkage] [--no_paired_bonus]\n"
                 "        [--longread_prior LongRead_Counts.tsv] [--min_lr_depth 5]\n"
                 "        [--index_dir DIR | --index-cache FILE] [--hits-cache FILE] [--py-sum neumaier|naive] [--binary]\n"
                 "  same options and output (OUTDIR/Global_Theta_Matrix.tsv) as theta/global_theta.py;\n"
                 "  the k-mer index is cached as DIR/lqcore_theta.k<k>.bin (or FILE);\n"
                 "  --hits-cache keeps the per-copy k-mer hits of this read set, so further runs with another\n"
                 "  --longread_prior (or none) skip the FASTQ scan; the output is the same as with a scan\n"
                 "  --binary writes OUTDIR/Global_Theta_Matrix.lqb instead (same numbers; `lq export` prints the TSV)\n");
}

TOpts parse_theta(int argc, char** argv) {
    TOpts o;
    auto need = [&](int& i) -> std::string {
        if (i + 1 >= argc) { theta_usage(); die(std::string("missing value for ") + argv[i]); }
        return argv[++i];
    };
    for (int i = 0; i < argc; ++i) {
        std::string a = argv[i];
        for (auto& c : a) if (c == '-' && &c != &a[0] && &c != &a[1]) c = '_';   // --index-dir == --index_dir
        if (a == "-c" || a == "--clusters") o.clusters = need(i);
        else if (a == "-f" || a == "--fasta") o.fasta = need(i);
        else if (a == "-r1" || a == "--read1") o.r1 = need(i);
        else if (a == "-r2" || a == "--read2") o.r2 = need(i);
        else if (a == "--longread_prior") o.lr_prior = need(i);
        else if (a == "-k" || a == "--kmer") o.k = std::atoi(need(i).c_str());
        else if (a == "-o" || a == "--outdir") o.outdir = need(i);
        else if (a == "--index_dir") o.index_dir = need(i);
        else if (a == "--index_cache") o.index_cache = need(i);
        else if (a == "--hits_cache") o.hits_cache = need(i);
        else if (a == "-t" || a == "--threads") o.threads = std::max(1, std::atoi(need(i).c_str()));
        else if (a == "--min_lr_depth") o.min_lr_depth = std::atoi(need(i).c_str());
        else if (a == "--shrinkage") { if (!py_float(need(i), o.shrinkage)) die("--shrinkage: not a number"); }
        else if (a == "--adaptive_shrinkage") o.adaptive = true;
        else if (a == "--no_paired_bonus") o.no_bonus = true;
        else if (a == "--py_sum") o.py_sum = need(i);
        else if (a == "--binary") o.binary = true;
        else if (a == "-h" || a == "--help") { theta_usage(); std::exit(0); }
        else { theta_usage(); die("unknown option " + std::string(argv[i])); }
    }
    if (o.clusters.empty() || o.fasta.empty() || o.r1.empty() || o.r2.empty()) { theta_usage(); die("-c, -f, -r1 and -r2 are required"); }
    if (o.k < 1 || o.k > 32) die("-k must be 1..32");
    if (o.py_sum != "neumaier" && o.py_sum != "naive") die("--py-sum must be neumaier or naive");
    return o;
}

// Python built-in sum() over floats (start = int 0; the first float becomes the running value)
struct PySum {
    bool neu;
    double f = 0, c = 0;
    bool any = false;
    explicit PySum(bool neumaier) : neu(neumaier) {}
    void add(double x) {
        if (!any) { f = 0.0 + x; any = true; return; }
        if (!neu) { f += x; return; }
        double t = f + x;
        if (std::fabs(f) >= std::fabs(x)) c += (f - t) + x;
        else c += (x - t) + f;
        f = t;
    }
    double get() const {
        if (!any) return 0.0;
        if (neu && c != 0.0 && std::isfinite(c)) return f + c;
        return f;
    }
};


// Python str.strip() (ASCII whitespace as str.isspace)
inline std::string_view py_strip(std::string_view s) {
    auto ws = [](unsigned char c) { return c == ' ' || (c >= 9 && c <= 13) || (c >= 0x1c && c <= 0x1f); };
    size_t a = 0, b = s.size();
    while (a < b && ws((unsigned char)s[a])) ++a;
    while (b > a && ws((unsigned char)s[b - 1])) --b;
    return s.substr(a, b - a);
}

// 0-3 ACGT (either case), 4 N/n, 5 anything else
struct CodeTab {
    uint8_t t[256];
    CodeTab() {
        for (int i = 0; i < 256; ++i) t[i] = 5;
        t['A'] = t['a'] = 0; t['C'] = t['c'] = 1; t['G'] = t['g'] = 2; t['T'] = t['t'] = 3; t['N'] = t['n'] = 4;
    }
};
const CodeTab CT;

inline char py_upper(char c) { return (c >= 'a' && c <= 'z') ? (char)(c - 32) : c; }
// str.translate(maketrans("ATCGN", "TAGCN")) after upper()
inline char comp_up(char c) {
    c = py_upper(c);
    switch (c) { case 'A': return 'T'; case 'T': return 'A'; case 'C': return 'G'; case 'G': return 'C'; default: return c; }
}

// ------------------------------------------------------------------ clusters (same parsing as global_theta.py)
struct Clusters {
    std::vector<std::string> cids;                       // dict insertion order
    std::vector<std::vector<uint32_t>> alleles;          // per cid, allele (gene) indices with duplicates
    std::vector<std::string> genes;                      // all_targets, first-seen order
    std::unordered_map<std::string, uint32_t> gidx;
    void load(const std::string& path) {
        LineReader rd(path);
        std::string_view l;
        std::unordered_map<std::string, uint32_t> cpos;
        while (rd.next(l)) {
            std::string_view s = py_strip(l);
            if (s.empty() || l.substr(0, 1) == "#") continue;
            auto f = split_tab(s);
            if (f.size() < 2) die("clusters: line without a tab: " + std::string(s.substr(0, 80)));
            std::string cid(f[0]);
            auto it = cpos.find(cid);
            uint32_t ci;
            if (it == cpos.end()) { ci = (uint32_t)cids.size(); cpos.emplace(cid, ci); cids.push_back(cid); alleles.emplace_back(); }
            else { ci = it->second; alleles[ci].clear(); }
            std::string_view al = f[1];
            size_t b = 0;
            for (;;) {
                size_t e = al.find(',', b);
                std::string_view a = py_strip(al.substr(b, e == std::string_view::npos ? std::string_view::npos : e - b));
                size_t d = a.find('.');
                std::string aid(d == std::string_view::npos ? a : a.substr(0, d));
                auto g = gidx.find(aid);
                uint32_t gi;
                if (g == gidx.end()) { gi = (uint32_t)genes.size(); gidx.emplace(aid, gi); genes.push_back(aid); }
                else gi = g->second;
                alleles[ci].push_back(gi);
                if (e == std::string_view::npos) break;
                b = e + 1;
            }
        }
    }
};

// ------------------------------------------------------------------ index
struct ThetaIndex {
    int k = 0, b = 0, shift = 0;
    std::vector<uint64_t> probes;                        // per gene (cluster-unique k-mers)
    std::vector<uint32_t> off;
    std::vector<uint64_t> code;
    std::vector<uint32_t> gid;
    std::unordered_map<std::string, uint32_t> iupac;     // unique k-mers with non-ACGT letters other than N

    inline int64_t get(uint64_t c) const {
        size_t bk = c >> shift;
        const uint64_t* lo = code.data() + off[bk];
        const uint64_t* hi = code.data() + off[bk + 1];
        const uint64_t* p = std::lower_bound(lo, hi, c);
        return (p != hi && *p == c) ? (int64_t)gid[p - code.data()] : -1;
    }
    void finish_offsets() {
        const size_t n = code.size();
        b = 8;
        while (b < 28 && b < 2 * k && (1ULL << b) < n / 8 + 1) ++b;
        if (b > 2 * k) b = 2 * k;
        shift = 2 * k - b;
        off.assign((1ULL << b) + 1, 0);
        for (uint64_t c : code) off[(c >> shift) + 1]++;
        for (size_t i = 1; i < off.size(); ++i) off[i] += off[i - 1];
    }
    void build(const std::string& fasta, const Clusters& C, int threads);
    bool load(const std::string& path, const std::string& sig, size_t ngenes);
    void save(const std::string& path, const std::string& sig) const;
};

void ThetaIndex::build(const std::string& fasta, const Clusters& C, int threads) {
    const double t0 = wall();
    // FASTA as Bio.SeqIO reads it: id = first word of the title; sequence lines right-stripped, spaces and \r removed
    std::string seq;
    std::vector<uint64_t> sbeg{0};
    std::vector<uint32_t> sg;
    {
        LineReader rd(fasta);
        std::string_view l;
        bool keep = false, in = false;
        auto close = [&]() { if (in && keep) { sbeg.push_back(seq.size()); } else if (in) { seq.resize(sbeg.back()); } };
        uint32_t cur = 0;
        while (rd.next(l)) {
            if (!l.empty() && l[0] == '>') {
                close();
                auto w = split_ws(l.substr(1));
                std::string id = w.empty() ? std::string() : std::string(w[0]);
                size_t d = id.find('.');
                std::string g = d == std::string::npos ? id : id.substr(0, d);
                auto it = C.gidx.find(g);
                keep = it != C.gidx.end();
                if (keep) cur = it->second;
                in = true;
                if (keep) sg.push_back(cur);
                continue;
            }
            if (!in || !keep) continue;
            size_t e = l.size();
            while (e > 0 && isspace((unsigned char)l[e - 1])) --e;
            for (size_t i = 0; i < e; ++i)
                if (l[i] != ' ' && l[i] != '\r') seq.push_back(py_upper(l[i]));
        }
        close();
    }
    const size_t R = sg.size();
    const int K = k;
    const uint64_t kmask = K == 32 ? ~0ULL : ((1ULL << (2 * K)) - 1);
    int gb = 1;
    while ((1ULL << gb) < C.genes.size() + 1) ++gb;
    if (2 * K + gb > 64) die(sfmt("lqcore theta: k=%d with %zu genes does not fit the packed index (2k + %d bits > 64); use k <= %d",
                                  K, C.genes.size(), gb, (64 - gb) / 2));
    // pass 1: ACGT windows per record
    std::vector<uint64_t> cnt(R + 1, 0);
#pragma omp parallel for schedule(dynamic, 256) num_threads(threads)
    for (size_t r = 0; r < R; ++r) {
        uint64_t n = 0;
        int run = 0;
        for (uint64_t i = sbeg[r]; i < sbeg[r + 1]; ++i) {
            if (CT.t[(uint8_t)seq[i]] < 4) { if (++run >= K) ++n; }
            else run = 0;
        }
        cnt[r + 1] = n;
    }
    for (size_t r = 0; r < R; ++r) cnt[r + 1] += cnt[r];
    std::vector<uint64_t> all(cnt[R]);
    std::vector<std::vector<std::pair<std::string, uint32_t>>> iu(threads);
    // pass 2: packed (code << gb | gene)
#pragma omp parallel for schedule(dynamic, 256) num_threads(threads)
    for (size_t r = 0; r < R; ++r) {
#ifdef _OPENMP
        auto& myiu = iu[omp_get_thread_num()];
#else
        auto& myiu = iu[0];
#endif
        uint64_t w = cnt[r], fw = 0;
        int run = 0;
        int64_t lastN = -1;
        const uint64_t b0 = sbeg[r], e0 = sbeg[r + 1];
        for (uint64_t i = b0; i < e0; ++i) {
            uint8_t c = CT.t[(uint8_t)seq[i]];
            if (c < 4) { fw = ((fw << 2) | c) & kmask; ++run; }
            else { run = 0; if (c == 4) lastN = (int64_t)i; }
            if (i - b0 + 1 < (uint64_t)K) continue;
            if (run >= K) all[w++] = (fw << gb) | sg[r];
            else if (lastN < (int64_t)(i + 1 - K))      // no N in the window, but another non-ACGT letter
                myiu.emplace_back(seq.substr(i + 1 - K, K), sg[r]);
        }
    }
    std::string().swap(seq);
    __gnu_parallel::sort(all.begin(), all.end());
    all.erase(std::unique(all.begin(), all.end()), all.end());
    probes.assign(C.genes.size(), 0);
    const uint64_t gmask = (1ULL << gb) - 1;
    size_t nu = 0;
    for (size_t i = 0; i < all.size();) {
        size_t j = i + 1;
        while (j < all.size() && (all[j] >> gb) == (all[i] >> gb)) ++j;
        if (j == i + 1) ++nu;
        i = j;
    }
    code.resize(nu);
    gid.resize(nu);
    nu = 0;
    for (size_t i = 0; i < all.size();) {
        size_t j = i + 1;
        while (j < all.size() && (all[j] >> gb) == (all[i] >> gb)) ++j;
        if (j == i + 1) {
            code[nu] = all[i] >> gb;
            gid[nu] = (uint32_t)(all[i] & gmask);
            probes[gid[nu]]++;
            ++nu;
        }
        i = j;
    }
    std::vector<uint64_t>().swap(all);
    std::map<std::string, std::set<uint32_t>> m;
    for (auto& v : iu) for (auto& x : v) m[x.first].insert(x.second);
    for (auto& x : m)
        if (x.second.size() == 1) { iupac.emplace(x.first, *x.second.begin()); probes[*x.second.begin()]++; }
    finish_offsets();
    logm(sfmt("[INFO] Theta index: %s cluster-unique %d-mers%s from %s transcripts (%.0fs)", commas((long long)code.size()).c_str(), K,
              iupac.empty() ? "" : sfmt(" (+%zu with IUPAC letters)", iupac.size()).c_str(), commas((long long)R).c_str(), wall() - t0));
}

const char TH_MAGIC[8] = {'L', 'Q', 'T', 'H', 'E', 'T', 'A', '1'};

void ThetaIndex::save(const std::string& path, const std::string& sig) const {
    std::string tmp = path + sfmt(".tmp%d", (int)getpid());
    FILE* f = fopen(tmp.c_str(), "wb");
    if (!f) { logm("[WARN] cannot write the Theta index cache " + path); return; }
    auto wr = [&](const void* p, size_t n) { return fwrite(p, 1, n, f) == n; };
    uint64_t sl = sig.size(), n = code.size(), np = probes.size(), ni = iupac.size();
    int32_t kk = k;
    bool ok = wr(TH_MAGIC, 8) && wr(&sl, 8) && wr(sig.data(), sl) && wr(&kk, 4) && wr(&np, 8) && wr(probes.data(), np * 8) &&
              wr(&n, 8) && wr(code.data(), n * 8) && wr(gid.data(), n * 4) && wr(&ni, 8);
    for (auto& x : iupac) {
        uint64_t l = x.first.size();
        ok = ok && wr(&l, 8) && wr(x.first.data(), l) && wr(&x.second, 4);
    }
    ok = (fclose(f) == 0) && ok;
    if (ok && rename(tmp.c_str(), path.c_str()) == 0) logm("[INFO] Theta index cached: " + path);
    else { unlink(tmp.c_str()); logm("[WARN] cannot write the Theta index cache " + path); }
}

bool ThetaIndex::load(const std::string& path, const std::string& sig, size_t ngenes) {
    FILE* f = fopen(path.c_str(), "rb");
    if (!f) return false;
    auto rd = [&](void* p, size_t n) { return fread(p, 1, n, f) == n; };
    char m[8];
    uint64_t sl = 0, n = 0, np = 0, ni = 0;
    int32_t kk = 0;
    bool ok = rd(m, 8) && memcmp(m, TH_MAGIC, 8) == 0 && rd(&sl, 8) && sl < (1u << 20);
    std::string s(ok ? sl : 0, '\0');
    ok = ok && rd(&s[0], sl) && s == sig && rd(&kk, 4) && kk == k && rd(&np, 8) && np == ngenes;
    if (ok) { probes.resize(np); ok = rd(probes.data(), np * 8) && rd(&n, 8); }
    if (ok) { code.resize(n); gid.resize(n); ok = rd(code.data(), n * 8) && rd(gid.data(), n * 4) && rd(&ni, 8); }
    for (uint64_t i = 0; ok && i < ni; ++i) {
        uint64_t l = 0;
        ok = rd(&l, 8) && l < 4096;
        std::string km(ok ? l : 0, '\0');
        uint32_t g = 0;
        ok = ok && rd(&km[0], l) && rd(&g, 4);
        if (ok) iupac.emplace(km, g);
    }
    fclose(f);
    if (!ok) { code.clear(); gid.clear(); probes.clear(); iupac.clear(); return false; }
    finish_offsets();
    return true;
}

void mkdir_p(const std::string& d) {
    if (d.empty()) return;
    for (size_t i = 1; i <= d.size(); ++i)
        if (i == d.size() || d[i] == '/') mkdir(d.substr(0, i).c_str(), 0775);
}

std::string file_sig(const std::string& p) {
    struct stat st;
    if (stat(p.c_str(), &st) != 0) die("cannot stat " + p);
    return p + "|" + std::to_string((long long)st.st_size) + "|" + std::to_string((long long)st.st_mtime);
}

// per-copy k-mer hits of one read set: "LQHITS1\n" | u64 siglen | sig | u64 G | i64 units | u64 hits[G]
bool load_hits(const std::string& path, const std::string& sig, size_t G, std::vector<uint64_t>& hits, long long& units) {
    FILE* f = fopen(path.c_str(), "rb");
    if (!f) return false;
    char mg[8];
    uint64_t sl = 0, g = 0;
    bool ok = fread(mg, 1, 8, f) == 8 && memcmp(mg, "LQHITS1\n", 8) == 0 && fread(&sl, 8, 1, f) == 1 && sl < (1u << 20);
    std::string s(ok ? sl : 0, '\0');
    ok = ok && fread(&s[0], 1, sl, f) == sl && s == sig && fread(&g, 8, 1, f) == 1 && g == G && fread(&units, 8, 1, f) == 1;
    if (ok) {
        hits.assign(G, 0);
        ok = fread(hits.data(), 8, G, f) == G;
    }
    fclose(f);
    if (!ok) logm("[INFO] hits cache " + path + " is for other reads/index -> scanning");
    return ok;
}
void save_hits(const std::string& path, const std::string& sig, const std::vector<uint64_t>& hits, long long units) {
    std::string tmp = path + sfmt(".tmp%d", (int)getpid());
    FILE* f = fopen(tmp.c_str(), "wb");
    if (!f) { logm("[WARN] cannot write " + path); return; }
    uint64_t sl = sig.size(), g = hits.size();
    bool ok = fwrite("LQHITS1\n", 1, 8, f) == 8 && fwrite(&sl, 8, 1, f) == 1 && fwrite(sig.data(), 1, sl, f) == sl &&
              fwrite(&g, 8, 1, f) == 1 && fwrite(&units, 8, 1, f) == 1 && fwrite(hits.data(), 8, g, f) == g;
    if (fclose(f) != 0 || !ok || rename(tmp.c_str(), path.c_str()) != 0) { unlink(tmp.c_str()); logm("[WARN] cannot write " + path); }
}

// ------------------------------------------------------------------ FASTQ scan
struct Batch {
    std::string data;
    std::vector<uint32_t> off{0};
    size_t n() const { return off.size() - 1; }
    std::string_view get(size_t i) const { return std::string_view(data.data() + off[i], off[i + 1] - off[i]); }
};

// Reads one FASTQ; `paired`: whole 4-line records only (zip of the two files); otherwise every line i % 4 == 1
class FqReader {
public:
    FqReader(const std::string& path, bool paired, size_t per_batch) : path_(path), paired_(paired), per_(per_batch) {}
    void run() {
        LineReader rd(path_);
        std::string_view l;
        auto b = std::make_unique<Batch>();
        std::string pending;
        long long li = 0;
        bool have = false;
        while (rd.next(l)) {
            int m = (int)(li % 4);
            if (m == 1) { pending.assign(py_strip(l)); have = true; if (!paired_) push_seq(b, pending); }
            if (paired_ && m == 3 && have) { push_seq(b, pending); have = false; }
            ++li;
            if (b->n() >= per_) {
                if (!put(std::move(b))) break;          // consumer stopped (the other file ended)
                b = std::make_unique<Batch>();
            }
        }
        if (b && b->n()) put(std::move(b));
        std::lock_guard<std::mutex> g(mu_);
        done_ = true;
        cv_.notify_all();
    }
    // next batch in order, or nullptr at the end
    std::unique_ptr<Batch> take() {
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
    bool paired_;
    size_t per_;
    std::deque<std::unique_ptr<Batch>> q_;
    std::mutex mu_;
    std::condition_variable cv_;
    bool done_ = false, stop_ = false;
    static void push_seq(std::unique_ptr<Batch>& b, const std::string& s) {
        b->data += s;
        b->off.push_back((uint32_t)b->data.size());
    }
    bool put(std::unique_ptr<Batch> b) {
        std::unique_lock<std::mutex> g(mu_);
        cv_.wait(g, [&] { return q_.size() < 8 || stop_; });
        if (stop_) return false;
        q_.push_back(std::move(b));
        cv_.notify_all();
        return true;
    }
};

struct Counts {
    std::vector<uint64_t> hits, st1, st2;
    uint64_t pid = 0;
    explicit Counts(size_t G) : hits(G, 0), st1(G, 0), st2(G, 0) {}
};

// calls f(gene) for every window hit, forward then reverse-complement orientation per window
template <class F>
inline void scan_read(const ThetaIndex& X, std::string_view s, F&& f) {
    const int K = X.k;
    const size_t L = s.size();
    if (L < (size_t)K) return;
    const uint64_t kmask = K == 32 ? ~0ULL : ((1ULL << (2 * K)) - 1);
    const int top = 2 * K - 2;
    uint64_t fw = 0, rc = 0;
    int run = 0;
    int64_t lastN = -1;
    const bool iu = !X.iupac.empty();
    std::string tmp;
    for (size_t i = 0; i < L; ++i) {
        uint8_t c = CT.t[(uint8_t)s[i]];
        if (c < 4) {
            fw = ((fw << 2) | c) & kmask;
            rc = (rc >> 2) | ((uint64_t)(3 - c) << top);
            ++run;
        } else { run = 0; if (c == 4) lastN = (int64_t)i; }
        if (i + 1 < (size_t)K) continue;
        if (run >= K) {
            int64_t g = X.get(fw);
            if (g >= 0) f((uint32_t)g);
            g = X.get(rc);
            if (g >= 0) f((uint32_t)g);
        } else if (iu && lastN < (int64_t)(i + 1 - K)) {
            tmp.resize(K);
            for (int j = 0; j < K; ++j) tmp[j] = py_upper(s[i + 1 - K + j]);
            auto it = X.iupac.find(tmp);
            if (it != X.iupac.end()) f(it->second);
            for (int j = 0; j < K; ++j) tmp[j] = comp_up(s[i - j]);
            it = X.iupac.find(tmp);
            if (it != X.iupac.end()) f(it->second);
        }
    }
}

// returns total hits per gene; also the number of pairs (or reads) scanned
std::vector<uint64_t> scan_fastq(const TOpts& o, const ThetaIndex& X, size_t G, long long& nunits) {
    const size_t PER = 50000;
    const int W = std::max(1, o.threads - 2);
    std::vector<Counts> cs;
    cs.reserve(W);
    for (int i = 0; i < W; ++i) cs.emplace_back(G);
    std::atomic<long long> units{0};

    auto process = [&](const std::string& r1, const std::string& r2, bool paired) {
        FqReader A(r1, paired, PER);
        std::unique_ptr<FqReader> B(paired ? new FqReader(r2, true, PER) : nullptr);
        std::thread ta([&] { A.run(); });
        std::thread tb;
        if (B) tb = std::thread([&] { B->run(); });
        std::mutex jm;
        bool end = false;
        auto next_job = [&](std::unique_ptr<Batch>& a, std::unique_ptr<Batch>& b) -> bool {
            std::lock_guard<std::mutex> g(jm);
            if (end) return false;
            a = A.take();
            if (B) b = B->take();
            if (!a || (B && !b)) { end = true; return false; }
            if (B && (a->n() < PER || b->n() < PER)) end = true;   // the shorter file ends the pairs
            return true;
        };
        std::vector<std::thread> ws;
        for (int w = 0; w < W; ++w)
            ws.emplace_back([&, w] {
                Counts& C = cs[w];
                std::unique_ptr<Batch> a, b;
                long long nu = 0;
                while (next_job(a, b)) {
                    if (!paired) {
                        for (size_t i = 0; i < a->n(); ++i) { scan_read(X, a->get(i), [&](uint32_t g) { C.hits[g]++; }); ++nu; }
                        continue;
                    }
                    const size_t n = std::min(a->n(), b->n());
                    for (size_t i = 0; i < n; ++i) {
                        std::string_view s1 = a->get(i), s2 = b->get(i);
                        if (s1.empty() || s2.empty()) continue;
                        const uint64_t p = ++C.pid;
                        scan_read(X, s1, [&](uint32_t g) { C.hits[g]++; C.st1[g] = p; });
                        scan_read(X, s2, [&](uint32_t g) {
                            C.hits[g]++;
                            if (C.st2[g] != p) { C.st2[g] = p; if (C.st1[g] == p) C.hits[g]++; }
                        });
                        ++nu;
                    }
                }
                units += nu;
            });
        for (auto& t : ws) t.join();
        A.stop();
        if (B) B->stop();
        // drain so the readers can finish
        while (A.take()) {}
        if (B) while (B->take()) {}
        ta.join();
        if (tb.joinable()) tb.join();
    };
    if (o.no_bonus) { process(o.r1, "", false); process(o.r2, "", false); }
    else process(o.r1, o.r2, true);
    std::vector<uint64_t> tot(G, 0);
    for (auto& C : cs) for (size_t g = 0; g < G; ++g) tot[g] += C.hits[g];
    nunits = units.load();
    return tot;
}

double cluster_shrinkage(const ThetaIndex& X, const std::vector<uint32_t>& al, double base) {
    std::vector<uint64_t> pr;
    for (uint32_t a : al) pr.push_back(X.probes[a]);
    if (pr.empty()) return base;
    std::sort(pr.begin(), pr.end());
    double med = (double)pr[pr.size() / 2];
    if (med <= 0) return base;
    double ad = base * std::sqrt(med / base);
    return std::max(3.0, std::min(30.0, ad));
}

}  // namespace

int theta_main(int argc, char** argv) {
    TOpts o = parse_theta(argc, argv);
    const double T0 = wall();
    const bool neu = o.py_sum == "neumaier";
    mkdir_p(o.outdir);
    std::string cache = !o.index_cache.empty() ? o.index_cache
                        : (!o.index_dir.empty() ? o.index_dir : o.outdir) + "/lqcore_theta.k" + std::to_string(o.k) + ".bin";
    if (!o.index_dir.empty()) mkdir_p(o.index_dir);
    const std::string out = o.outdir + (o.binary ? "/Global_Theta_Matrix.lqb" : "/Global_Theta_Matrix.tsv");

    // long-read priors (optional)
    std::unordered_map<std::string, double> lr;
    const bool hybrid = !o.lr_prior.empty() && file_exists(o.lr_prior);
    if (hybrid) {
        logm("[INFO] Loading long-read physical priors (hybrid mode)...");
        LineReader rd(o.lr_prior);
        std::string_view l;
        bool first = true;
        while (rd.next(l)) {
            if (first) { first = false; continue; }
            std::string_view s = py_strip(l);
            if (s.empty()) continue;
            auto f = split_tab(s);
            if (f.size() != 2) die("long-read prior: expected 2 columns: " + std::string(s.substr(0, 80)));
            double v;
            if (!py_float(f[1], v)) die("long-read prior: not a number: " + std::string(f[1]));
            lr[std::string(f[0])] = v;
        }
    } else logm("[INFO] No long-read prior; pure k-mer Theta mode.");

    Clusters C;
    C.load(o.clusters);
    logm(sfmt("[INFO] clusters: %s | copies: %s | k=%d | threads %d | sum() %s", commas((long long)C.cids.size()).c_str(),
              commas((long long)C.genes.size()).c_str(), o.k, o.threads, neu ? "Python>=3.12 (Neumaier)" : "Python<=3.11 (plain)"));

    ThetaIndex X;
    X.k = o.k;
    const std::string sig = "k=" + std::to_string(o.k) + "|" + file_sig(o.fasta) + "|" + file_sig(o.clusters);
    double ti = wall();
    if (X.load(cache, sig, C.genes.size())) {
        logm(sfmt("[INFO] Theta index loaded from %s: %s k-mers (%.0fs)", cache.c_str(), commas((long long)X.code.size()).c_str(), wall() - ti));
    } else {
        X.build(o.fasta, C, o.threads);
        X.save(cache, sig);
    }

    double ts = wall();
    long long nunits = 0;
    std::vector<uint64_t> hits;
    const std::string hsig = sig + "|" + file_sig(o.r1) + "|" + file_sig(o.r2) + (o.no_bonus ? "|single" : "|paired");
    if (!o.hits_cache.empty() && load_hits(o.hits_cache, hsig, C.genes.size(), hits, nunits)) {
        logm(sfmt("[INFO] k-mer hits loaded from %s (%s %s scanned before)", o.hits_cache.c_str(), commas(nunits).c_str(),
                  o.no_bonus ? "reads" : "pairs"));
    } else {
        hits = scan_fastq(o, X, C.genes.size(), nunits);
        logm(sfmt("[INFO] FASTQ scan (%s): %s %s, %.0fs", o.no_bonus ? "single-end" : "paired-end + bonus",
                  commas(nunits).c_str(), o.no_bonus ? "reads" : "pairs", wall() - ts));
        if (!o.hits_cache.empty()) save_hits(o.hits_cache, hsig, hits, nunits);
    }

    // Theta (same operations and order as global_theta.py)
    ThetaBin TB;
    TB.clusters = C.cids;
    TB.rows.reserve(C.genes.size());
    auto add_row = [&](uint32_t ci, uint32_t a, double theta, bool anchored, double shrink) {
        ThetaRow r;
        r.cl = ci;
        r.allele = C.genes[a];
        r.probes = X.probes[a];
        r.hits = hits[a];
        auto it = lr.find(C.genes[a]);
        if (it != lr.end()) { r.flags |= THF_LR; r.lr = it->second; }
        if (anchored) r.flags |= THF_ANCHORED;
        r.theta = theta;
        r.shrink = shrink;
        TB.rows.push_back(std::move(r));
    };
    long long n_lr = 0, n_km = 0;
    for (size_t ci = 0; ci < C.cids.size(); ++ci) {
        const auto& al = C.alleles[ci];
        double cluster_lr = 0.0;
        bool lr_any = false;
        if (hybrid) {
            PySum s(neu);
            for (uint32_t a : al) { auto it = lr.find(C.genes[a]); if (it != lr.end()) { s.add(it->second); lr_any = true; } }
            cluster_lr = lr_any ? s.get() : 0.0;
        }
        if (hybrid && cluster_lr >= o.min_lr_depth) {
            ++n_lr;
            const double eps = 0.01;
            const double total = cluster_lr + (double)al.size() * eps;
            for (uint32_t a : al) {
                auto it = lr.find(C.genes[a]);
                double v = it == lr.end() ? 0.0 : it->second;
                add_row((uint32_t)ci, a, (v + eps) / total, true, 0.0);
            }
        } else {
            ++n_km;
            // raw = h / p; total_p = sum(p); cmean = sum(raw * p) / total_p
            unsigned long long total_p = 0;
            PySum sm(neu);
            for (uint32_t a : al) {
                const uint64_t p = X.probes[a];
                if (p > 0) {
                    double raw = (double)hits[a] / (double)p;
                    total_p += p;
                    sm.add(raw * (double)p);
                }
            }
            const double cmean = total_p > 0 ? sm.get() / (double)total_p : 0.0;
            const double shrink = o.adaptive ? cluster_shrinkage(X, al, o.shrinkage) : o.shrinkage;
            // shrunk: dict keyed by allele (first-insertion order, last value wins -- identical values anyway)
            std::vector<uint32_t> keys;
            std::unordered_map<uint32_t, double> shr;
            for (uint32_t a : al) {
                const uint64_t p = X.probes[a];
                double v;
                if (p > 0) {
                    double raw = (double)hits[a] / (double)p;
                    double w = (double)p / ((double)p + shrink);
                    double x = w * raw + (1 - w) * cmean;
                    v = (1e-8 > x) ? 1e-8 : x;
                } else {
                    double x = cmean > 0 ? cmean : 1e-6;
                    v = (1e-8 > x) ? 1e-8 : x;
                }
                if (shr.emplace(a, v).second) keys.push_back(a);
                else shr[a] = v;
            }
            PySum st(neu);
            for (uint32_t a : keys) st.add(shr[a]);
            const double tot = st.get();
            for (uint32_t a : al) add_row((uint32_t)ci, a, tot > 0 ? shr[a] / tot : 0.0, false, shrink);
        }
    }
    if (o.binary) write_theta_bin(out, TB);
    else write_text_file(out, [&](FILE* f) {
        theta_tsv_header(f);
        for (auto& r : TB.rows) theta_tsv_row(f, TB, r);
    });
    logm(sfmt("[SUCCESS] %s (%s clusters: %s k-mer inferred%s) %.0fs", out.c_str(), commas((long long)C.cids.size()).c_str(),
              commas(n_km).c_str(), hybrid ? sfmt(", %s long-read anchored", commas(n_lr).c_str()).c_str() : "", wall() - T0));
    return 0;
}

}  // namespace lq
