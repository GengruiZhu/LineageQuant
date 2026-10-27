// lqcore -- model, evidence, parts and the EM solver (port of core/em_global.py)
#include "model.hpp"
#include "binio.hpp"

#include <atomic>
#include <cmath>
#include <memory>
#include <pthread.h>
#include <sched.h>
#include <set>
#include <unordered_set>
#ifdef _OPENMP
#include <omp.h>
#else
static inline int omp_get_num_threads() { return 1; }
static inline int omp_get_thread_num() { return 0; }
#endif

namespace lq {

// ================================================================== options
Opts parse_opts(int argc, char** argv, void (*usage)()) {
    Opts a;
    auto num = [](const std::string& k, const std::string& v) {
        double x;
        if (!py_float(v, x)) die("invalid value for " + k + ": " + v);
        return x;
    };
    for (int i = 0; i < argc; ++i) {
        std::string k = argv[i], v;
        bool has_eq = false;
        if (k.rfind("--", 0) == 0) {
            size_t e = k.find('=');
            if (e != std::string::npos) { v = k.substr(e + 1); k = k.substr(0, e); has_eq = true; }
        }
        auto val = [&]() -> std::string {
            if (has_eq) return v;
            if (i + 1 >= argc) die("missing value for " + k);
            return argv[++i];
        };
        if (k == "-h" || k == "--help") { usage(); std::exit(0); }
        // ---- EM (em_global.py names)
        else if (k == "-t" || k == "--theta") a.theta = val();
        else if (k == "-p" || k == "--pvals") a.pvals = val();
        else if (k == "-a" || k == "--ancestry") a.ancestry = val();
        else if (k == "-o" || k == "--out") a.out = val();
        else if (k == "--bam") a.bam = val();
        else if (k == "--sam" || k == "--sam-in") a.sam = val();
        else if (k == "--fasta") a.fasta = val();
        else if (k == "--no-pool") a.no_pool = true;
        else if (k == "--keep-collate") a.keep_collate = true;
        else if (k == "--fld-weight") a.fld_weight = true;
        else if (k == "--threads") a.threads = (int)num(k, val());
        else if (k == "--workers") a.workers = (int)num(k, val());
        else if (k == "--eff-len") a.eff_len = val();
        else if (k == "--frag-len") a.frag_len = num(k, val());
        else if (k == "--min-eff-len") a.min_eff_len = num(k, val());
        else if (k == "--gc-bias") a.gc_bias = true;
        else if (k == "--gc-hist") a.gc_hist = true;
        else if (k == "--gc-fl-samples") a.gc_fl_samples = (int)num(k, val());
        else if (k == "--eff-out") a.eff_out = val();
        else if (k == "--marker-clean") a.marker_clean = true;
        else if (k == "--max-iter") a.max_iter = (int)num(k, val());
        else if (k == "--min-iter") a.min_iter = (int)num(k, val());
        else if (k == "--tol") a.tol = num(k, val());
        else if (k == "--init") a.init = val();
        else if (k == "--p-weight-specific") a.pw_spec = num(k, val());
        else if (k == "--p-weight-conflict") a.pw_conf = num(k, val());
        else if (k == "--p-weight-none") a.pw_none = num(k, val());
        else if (k == "--p-alpha") a.p_alpha = num(k, val());
        else if (k == "--p-floor") a.p_floor = num(k, val());
        else if (k == "--theta-prior-lr") a.tp_lr = num(k, val());
        else if (k == "--theta-prior-kmer") a.tp_kmer = num(k, val());
        else if (k == "--theta-prior-cap") a.tp_cap = num(k, val());
        else if (k == "--nm-weight") a.nm_weight = val();
        else if (k == "--em-accel") a.em_accel = val();
        else if (k == "--no-pin") a.no_pin = true;
        else if (k == "--aligner") a.aligner = val();
        else if (k == "--map-k") a.map_k = (int)num(k, val());
        else if (k == "--map-stride") a.map_stride = (int)num(k, val());
        else if (k == "--map-max-extra") a.map_max_extra = (int)num(k, val());
        else if (k == "--map-max-occ") a.map_max_occ = (int)num(k, val());
        else if (k == "--map-max-cand") a.map_max_cand = (int)num(k, val());
        else if (k == "--map-max-frag") a.map_max_frag = (int)num(k, val());
        else if (k == "--map-mm-frac") a.map_mm_frac = num(k, val());
        else if (k == "--map-no-rescue") a.map_no_rescue = true;
        else if (k == "--map-cache") a.map_cache = val();
        else if (k == "--salmon") a.salmon = val();
        else if (k == "--salmon-index") a.salmon_index = val();
        else if (k == "--salmon-opts") a.salmon_opts = val();
        else if (k == "--parts") {
            while (i + 1 < argc && std::string(argv[i + 1]).rfind("-", 0) != 0) a.parts.push_back(argv[++i]);
        }
        else if (k == "--idxstats-out") a.idxstats_out = val();
        else if (k == "--write-pvals") a.write_pvals = val();
        // ---- run: alignment
        else if (k == "--r1") a.r1 = val();
        else if (k == "--r2") a.r2 = val();
        else if (k == "-x" || k == "--index") a.index = val();
        else if (k == "--bowtie2") a.bowtie2 = val();
        else if (k == "--bt2-opts") a.bt2_opts = val();
        else if (k == "--bt2-log") a.bt2_log = val();
        else if (k == "--bt2-threads") a.bt2_threads = (int)num(k, val());
        else if (k == "-k" || k == "--max-alignments") a.max_aln = (int)num(k, val());
        else if (k == "--chunk") {
            std::string s = val();
            size_t sl = s.find('/');
            if (sl == std::string::npos) die("--chunk must be i/N (1-based), e.g. 2/4");
            a.chunk_i = atoi(s.substr(0, sl).c_str());
            a.chunk_n = atoi(s.substr(sl + 1).c_str());
            if (a.chunk_n < 1 || a.chunk_i < 1 || a.chunk_i > a.chunk_n) die("--chunk must be i/N with 1 <= i <= N");
        }
        else if (k == "--block-pairs") a.block_pairs = (long long)num(k, val());
        else if (k == "--part-out") a.part_out = val();
        else if (k == "--tmpdir") a.tmpdir = val();
        // ---- run: lineage
        else if (k == "--target-ids") a.target_ids = val();
        else if (k == "--markers") a.markers.push_back(val());
        else if (k == "--kmer-len" || k == "--lineage-kmer") a.kmer_len = (int)num(k, val());
        else if (k == "--kmer-alpha") a.kmer_alpha = num(k, val());
        else if (k == "--min-hits") a.min_hits = (int)num(k, val());
        else if (k == "--specific-thresh") a.spec_thresh = num(k, val());
        else if (k == "--marker-cache") a.marker_cache = val();
        else if (k == "--snp-db") a.snp_db = val();
        else if (k == "--snp-pick") a.snp_pick = val();
        else if (k == "--error-rate") a.error_rate = num(k, val());
        else if (k == "--no-chimeric") a.no_chimeric = true;
        else if (k == "--poe-min-weight") a.poe_min_weight = num(k, val());
        else if (k == "--poe-eps") a.poe_eps = num(k, val());
        else if (k == "--debug-kmer") a.debug_kmer = val();
        else if (k == "--debug-frag") a.debug_frag = val();
        else if (k == "--debug-snp") a.debug_snp = val();
        else a.unknown.push_back(argv[i]);
    }
    if (a.eff_len != "off" && a.eff_len != "len" && a.eff_len != "eff" && a.eff_len != "fld") die("--eff-len must be off|len|eff|fld");
    if (a.gc_bias && a.eff_len == "off") die("--gc-bias corrects effective lengths: not with --eff-len off");
    if (a.gc_fl_samples < 1 || a.gc_fl_samples > 1000) die("--gc-fl-samples must be 1..1000");
    if (a.init != "theta" && a.init != "uniform") die("--init must be theta|uniform");
    if (a.em_accel != "none" && a.em_accel != "squarem") die("--em-accel must be none|squarem");
    if (a.aligner != "bowtie2" && a.aligner != "lqmap" && a.aligner != "sa") die("--aligner must be bowtie2|lqmap|sa");
    if (a.fld_weight) die("--fld-weight (experimental G5) is not implemented in lqcore; use the Python engine");
    if (a.snp_pick != "best" && a.snp_pick != "most" && a.snp_pick != "first") die("--snp-pick must be best|most|first");
    if (a.threads < 1) a.threads = 1;
    return a;
}

void LineageStats::add(const LineageStats& o) {
    valid += o.valid; both += o.both; kmer_only += o.kmer_only; snp_only += o.snp_only; pv_used += o.pv_used;
    snp_rows += o.snp_rows;
    for (int i = 0; i < 4; ++i) { kmer_st[i] += o.kmer_st[i]; snp_st[i] += o.snp_st[i]; }
    for (int i = 0; i < 11; ++i) hits[i] += o.hits[i];
    for (int i = 0; i < 6; ++i) fused_st[i] += o.fused_st[i];
}

// ================================================================== model
void Model::load_theta_ancestry(const Opts& o) {
    if (!o.theta.empty() && file_exists(o.theta) && has_magic(o.theta, THETA_MAGIC)) {
        // binary Theta (theta --binary): the values the text would give after printing and parsing
        ThetaBin TB;
        read_theta_bin(o.theta, TB);
        for (auto& r : TB.rows) {
            const std::string& cl = TB.clusters[r.cl];
            if (!g2c.count(r.allele)) g2c_order.push_back(r.allele);
            g2c[r.allele] = cl;
            theta0[r.allele] = py_round(r.theta, 6);
            lr_count[cl] += (r.flags & THF_LR) ? r.lr : 0.0;
            prior_status[cl] = (r.flags & THF_ANCHORED) ? "LongRead_Anchored" : "Kmer_Inferred";
        }
    } else if (!o.theta.empty() && file_exists(o.theta)) {
        LineReader rd(o.theta);
        std::string_view l;
        bool first = true;
        while (rd.next(l)) {
            if (first) { first = false; continue; }
            auto p = split_tab(l);
            if (p.size() < 6) continue;
            std::string al(p[1]), cl(p[0]);
            if (!g2c.count(al)) g2c_order.push_back(al);
            g2c[al] = cl;
            double v;
            theta0[al] = py_float(p[5], v) ? v : 0.0;
            if (py_float(p[4], v)) lr_count[cl] += v;
            else lr_count.emplace(cl, 0.0);
            if (p.size() > 6) prior_status[cl] = std::string(p[6]);
        }
    }
    if (file_exists(o.ancestry)) {
        LineReader rd(o.ancestry);
        std::string_view l;
        while (rd.next(l)) {
            auto p = split_ws(l);
            if (p.size() >= 2 && l.rfind("GeneID", 0) != 0) {
                std::string u = upper(p[1]);
                anc[std::string(p[0])] = u;
                anc.emplace(std::string(rsplit_dot(p[0])), u);
            }
        }
    }
    std::unordered_set<std::string> cs;
    for (auto& kv : g2c) cs.insert(kv.second);
    logm("[INFO] Theta: " + commas((long long)cs.size()) + " clusters / " + commas((long long)g2c.size()) + " allele copies");
}

void load_ref_gc(RefGC& g, const std::string& fasta, const std::vector<std::string>& refs,
                 const std::unordered_map<std::string, uint32_t>& r2i) {
    const size_t R = refs.size();
    g.off.assign(R, ~0ULL);
    g.len.assign(R, 0);
    std::vector<std::string> seq(R);
    {
        LineReader rd(fasta);
        std::string_view l;
        long long cur = -1;
        while (rd.next(l)) {
            if (!l.empty() && l[0] == '>') {
                auto ws = split_ws(l.substr(1));
                cur = -1;
                if (!ws.empty()) {
                    auto it = r2i.find(std::string(ws[0]));
                    if (it != r2i.end() && seq[it->second].empty()) cur = it->second;
                }
                continue;
            }
            if (cur < 0) continue;
            for (char c : l) if (!isspace((unsigned char)c)) seq[cur].push_back(c);
        }
    }
    uint64_t tot = 0;
    for (size_t r = 0; r < R; ++r) if (!seq[r].empty()) { g.off[r] = tot; g.len[r] = (long long)seq[r].size(); tot += seq[r].size(); }
    g.bits.assign((tot + 63) / 64 + 1, 0);
    size_t nseq = 0;
    for (size_t r = 0; r < R; ++r) {
        if (seq[r].empty()) continue;
        ++nseq;
        const uint64_t o = g.off[r];
        for (size_t i = 0; i < seq[r].size(); ++i) {
            const char c = seq[r][i];
            if (c == 'G' || c == 'C' || c == 'g' || c == 'c' || c == 'S' || c == 's') g.bits[(o + i) >> 6] |= 1ULL << ((o + i) & 63);
        }
        std::string().swap(seq[r]);
    }
    logm(sfmt("[INFO] [GC] reference G/C map: %s of %s sequences, %s bp", commas((long long)nseq).c_str(),
              commas((long long)R).c_str(), commas((long long)tot).c_str()));
}

void Model::set_refs(std::vector<std::string> r, std::vector<long long> l) {
    refs = std::move(r);
    lens = std::move(l);
    r2i.clear();
    r2i.reserve(refs.size() * 2);
    for (uint32_t i = 0; i < refs.size(); ++i) r2i[refs[i]] = i;
}

void Model::build_targets(const Opts& o) {
    std::unordered_map<std::string, uint32_t> tkey;
    ref2t.assign(refs.size(), 0);
    for (size_t r = 0; r < refs.size(); ++r) {
        const std::string& ref = refs[r];
        std::string r1(rsplit_dot(ref));
        const std::string* k = nullptr;
        auto i1 = g2c.find(r1);
        if (i1 != g2c.end()) k = &i1->first;
        else {
            auto i2 = g2c.find(ref);
            if (i2 != g2c.end()) k = &i2->first;
        }
        std::string key = k ? ("A\t" + *k) : ("O\t" + strip_iso(ref));
        auto it = tkey.find(key);
        if (it == tkey.end()) {
            uint32_t id = (uint32_t)names.size();
            tkey.emplace(key, id);
            names.push_back(k ? *k : ref);
            tlen.push_back((double)lens[r]);
            if (k) { tcl.push_back(g2c[*k]); has_cl.push_back(1); }
            else { tcl.push_back(""); has_cl.push_back(0); }
            ref2t[r] = id;
        } else {
            tlen[it->second] = std::max(tlen[it->second], (double)lens[r]);
            ref2t[r] = it->second;
        }
    }
    for (auto& g : g2c_order) {
        std::string key = "A\t" + g;
        if (!tkey.count(key)) {
            uint32_t id = (uint32_t)names.size();
            tkey.emplace(key, id);
            names.push_back(g);
            tlen.push_back(1000.0);
            tcl.push_back(g2c[g]);
            has_cl.push_back(1);
        }
    }
    T = (uint32_t)names.size();
    is_al = has_cl;
    n_al = 0;
    for (uint32_t i = 0; i < T; ++i) n_al += has_cl[i];
    {
        std::set<std::string> s;
        for (uint32_t i = 0; i < T; ++i) if (has_cl[i]) s.insert(tcl[i]);
        cl_sorted.assign(s.begin(), s.end());
    }
    std::unordered_map<std::string, uint32_t> cl_ids;
    for (uint32_t j = 0; j < cl_sorted.size(); ++j) cl_ids[cl_sorted[j]] = j;
    cl_int.assign(T, NOCL);
    for (uint32_t i = 0; i < T; ++i) if (has_cl[i]) cl_int[i] = cl_ids[tcl[i]];
    logm("[INFO] " + commas((long long)refs.size()) + " reference sequences -> " + commas(T) + " genes (" +
         commas(n_al) + " allele copies / " + commas((long long)T - n_al) + " orphans)");

    // ---- [G3] identical-sequence pooling
    POOLG.assign(T, {});
    pool_reps.clear();
    n_pooled = 0;
    pool_on = !o.fasta.empty() && !o.no_pool;
    if (pool_on && file_exists(o.fasta)) {
        std::unordered_map<std::string, Hash128> md5;
        {
            LineReader rd(o.fasta);
            std::string_view l;
            std::string nm;
            bool have = false;
            SeqHasher h;
            while (rd.next(l)) {
                if (!l.empty() && l[0] == '>') {
                    if (have && !nm.empty()) md5[nm] = h.get();
                    auto ws = split_ws(l.substr(1));
                    nm = ws.empty() ? std::string() : std::string(ws[0]);
                    h = SeqHasher();
                    have = true;
                } else {
                    size_t b = 0, e = l.size();
                    while (b < e && isspace((unsigned char)l[b])) ++b;
                    while (e > b && isspace((unsigned char)l[e - 1])) --e;
                    h.add(l.substr(b, e - b));
                }
            }
            if (have && !nm.empty()) md5[nm] = h.get();
        }
        std::vector<std::vector<uint32_t>> t_refs(T);
        std::vector<uint32_t> t_order;
        for (size_t r = 0; r < refs.size(); ++r) {
            uint32_t t = ref2t[r];
            if (t_refs[t].empty()) t_order.push_back(t);
            t_refs[t].push_back((uint32_t)r);
        }
        std::unordered_map<Hash128, uint32_t, Hash128H> by_idx;
        std::vector<std::vector<uint32_t>> by;
        for (uint32_t t : t_order) {
            auto& rs = t_refs[t];
            if (rs.size() != 1) continue;
            auto it = md5.find(refs[rs[0]]);
            if (it == md5.end()) continue;
            auto jt = by_idx.find(it->second);
            if (jt == by_idx.end()) { by_idx.emplace(it->second, (uint32_t)by.size()); by.push_back({t}); }
            else by[jt->second].push_back(t);
        }
        std::vector<uint32_t> POOL(T, NOCL);
        for (auto& ts : by) {
            if (ts.size() < 2) continue;
            uint32_t mn = *std::min_element(ts.begin(), ts.end());
            for (uint32_t t : ts) POOL[t] = mn;
        }
        std::vector<char> seen(T, 0);
        for (auto& ts : by) {
            if (ts.size() < 2) continue;
            for (uint32_t t : ts) {
                uint32_t r = POOL[t];
                if (!seen[r]) { seen[r] = 1; pool_reps.push_back(r); }
                POOLG[r].push_back(t);
                ++n_pooled;
            }
        }
        for (auto& t : ref2t) if (POOL[t] != NOCL) t = POOL[t];
        logm("[INFO] identical-sequence genes pooled: " + commas((long long)n_pooled) + " -> " +
             commas((long long)pool_reps.size()) + " groups (split evenly after convergence)");
    }
    tlin.assign(T, -1);
}

void Model::set_lineages(const std::vector<std::string>& lins) {
    lineages = lins;
    std::unordered_map<std::string, int> lin_rank;
    for (size_t k = 0; k < lins.size(); ++k) lin_rank[lins[k]] = (int)k;
    tlin.assign(T, -1);
    for (uint32_t i = 0; i < T; ++i) {
        if (!is_al[i]) continue;
        const std::string& n = names[i];
        std::string lin = "?";
        auto it = anc.find(n);
        if (it != anc.end()) lin = it->second;
        else {
            auto jt = anc.find(std::string(rsplit_dot(n)));
            if (jt != anc.end()) lin = jt->second;
        }
        auto r = lin_rank.find(lin);
        tlin[i] = r == lin_rank.end() ? -1 : r->second;
    }
}

uint64_t Model::fingerprint() const {
    uint64_t h = 0xcbf29ce484222325ULL;
    auto mix = [&](const void* p, size_t n) {
        const unsigned char* c = (const unsigned char*)p;
        for (size_t i = 0; i < n; ++i) h = (h ^ c[i]) * 0x100000001b3ULL;
    };
    mix(&T, sizeof T);
    for (auto& s : names) { mix(s.data(), s.size()); mix("\n", 1); }
    mix(ref2t.data(), ref2t.size() * 4);
    for (auto& g : POOLG) { size_t n = g.size(); mix(&n, sizeof n); if (n) mix(g.data(), n * 4); }
    mix(tlin.data(), tlin.size() * sizeof(int));
    for (auto& s : lineages) { mix(s.data(), s.size()); mix("\n", 1); }
    return h;
}

// ================================================================== parts
namespace {
const char PART_MAGIC[8] = {'L', 'Q', 'P', 'A', 'R', 'T', '0', '2'};      // 02 (1.1.0): + fragment GC histogram
const char PART_MAGIC1[8] = {'L', 'Q', 'P', 'A', 'R', 'T', '0', '1'};     // 01: 0.x / 1.0 parts, still readable
struct BW {
    FILE* f;
    void raw(const void* p, size_t n) { if (n && fwrite(p, 1, n, f) != n) die("write error (part file)"); }
    void u64(uint64_t v) { raw(&v, 8); }
    void i64(long long v) { raw(&v, 8); }
    void str(const std::string& s) { u64(s.size()); raw(s.data(), s.size()); }
};
struct BR {
    FILE* f;
    std::string path;
    int ver = 2;
    void raw(void* p, size_t n) { if (n && fread(p, 1, n, f) != n) die("truncated part file: " + path); }
    uint64_t u64() { uint64_t v; raw(&v, 8); return v; }
    long long i64() { long long v; raw(&v, 8); return v; }
    std::string str() { std::string s(u64(), '\0'); raw(&s[0], s.size()); return s; }
};
void read_header(BR& r, std::vector<std::string>& refs, std::vector<long long>& lens, std::vector<std::string>& lins,
                 uint64_t& fp) {
    char m[8];
    r.raw(m, 8);
    if (memcmp(m, PART_MAGIC, 8) == 0) r.ver = 2;
    else if (memcmp(m, PART_MAGIC1, 8) == 0) r.ver = 1;
    else die("not an lqcore part file: " + r.path);
    uint64_t R = r.u64();
    refs.resize(R);
    lens.resize(R);
    for (uint64_t i = 0; i < R; ++i) { refs[i] = r.str(); lens[i] = r.i64(); }
    uint64_t L = r.u64();
    lins.resize(L);
    for (auto& s : lins) s = r.str();
    fp = r.u64();
}
}  // namespace

void write_part(const std::string& path, const Model& M, const Evidence& ev, const std::vector<std::string>& lineages) {
    std::string tmp = path + ".tmp";
    FILE* f = fopen(tmp.c_str(), "wb");
    if (!f) die("cannot write " + tmp);
    std::vector<char> buf(1 << 22);
    setvbuf(f, buf.data(), _IOFBF, buf.size());
    BW w{f};
    w.raw(PART_MAGIC, 8);
    w.u64(M.refs.size());
    for (size_t i = 0; i < M.refs.size(); ++i) { w.str(M.refs[i]); w.i64(M.lens[i]); }
    w.u64(lineages.size());
    for (auto& s : lineages) w.str(s);
    w.u64(M.fingerprint());
    w.i64(ev.n_rec); w.i64(ev.n_frag); w.i64(ev.n_unal); w.i64(ev.nm_sum); w.i64(ev.len_sum); w.i64(ev.n_pv_reads);
    w.u64(ev.fl_vals.size());
    w.raw(ev.fl_vals.data(), ev.fl_vals.size() * 4);
    w.u64(ev.single.size());
    w.raw(ev.single.data(), ev.single.size() * 8);
    w.u64(ev.idx_mapped.size());
    w.raw(ev.idx_mapped.data(), ev.idx_mapped.size() * 8);
    w.raw(ev.idx_unmapped.data(), ev.idx_unmapped.size() * 8);
    w.u64(ev.idx_unmapped_star);
    w.raw(&ev.ls, sizeof(LineageStats));
    for (const KeyMap* km : {&ev.eq, &ev.eqx}) {
        w.u64(km->size());
        for (size_t c = 0; c < km->size(); ++c) { w.str(*km->keys[c]); w.u64(km->cnt[c]); }
    }
    w.u64(ev.gc_fg.size());
    w.raw(ev.gc_fg.data(), ev.gc_fg.size() * 8);
    if (fclose(f) != 0) die("error writing " + tmp);
    if (rename(tmp.c_str(), path.c_str()) != 0) die("cannot rename " + tmp);
}

void read_part_refs(const std::string& path, std::vector<std::string>& refs, std::vector<long long>& lens,
                    std::vector<std::string>& lineages) {
    FILE* f = fopen(path.c_str(), "rb");
    if (!f) die("cannot open part " + path);
    BR r{f, path};
    uint64_t fp;
    read_header(r, refs, lens, lineages, fp);
    fclose(f);
}

void merge_part(const std::string& path, const Model& M, Evidence& ev) {
    FILE* f = fopen(path.c_str(), "rb");
    if (!f) die("cannot open part " + path);
    std::vector<char> buf(1 << 22);
    setvbuf(f, buf.data(), _IOFBF, buf.size());
    BR r{f, path};
    std::vector<std::string> refs, lins;
    std::vector<long long> lens;
    uint64_t fp;
    read_header(r, refs, lens, lins, fp);
    if (refs != M.refs) die("part " + path + " was made against different reference sequences");
    if (fp != M.fingerprint()) die("part " + path + " was made with a different model (Theta / ancestry / FASTA / pooling / lineages)");
    ev.n_rec += r.i64(); ev.n_frag += r.i64(); ev.n_unal += r.i64(); ev.nm_sum += r.i64(); ev.len_sum += r.i64();
    ev.n_pv_reads += r.i64();
    std::vector<int32_t> fl(r.u64());
    r.raw(fl.data(), fl.size() * 4);
    for (int32_t v : fl) { if (ev.fl_vals.size() >= 1000000) break; ev.fl_vals.push_back(v); }
    std::vector<double> single(r.u64());
    r.raw(single.data(), single.size() * 8);
    if (single.size() != ev.single.size()) die("part " + path + ": target count differs");
    for (size_t t = 0; t < single.size(); ++t) ev.single[t] += single[t];
    std::vector<uint64_t> im(r.u64()), iu(im.size());
    r.raw(im.data(), im.size() * 8);
    r.raw(iu.data(), iu.size() * 8);
    for (size_t i = 0; i < im.size(); ++i) { ev.idx_mapped[i] += im[i]; ev.idx_unmapped[i] += iu[i]; }
    ev.idx_unmapped_star += r.u64();
    LineageStats ls;
    r.raw(&ls, sizeof ls);
    ev.ls.add(ls);
    for (KeyMap* km : {&ev.eq, &ev.eqx}) {
        uint64_t n = r.u64();
        for (uint64_t c = 0; c < n; ++c) {
            std::string k = r.str();
            km->add(k, r.u64());
        }
    }
    if (r.ver >= 2) {
        std::vector<double> g(r.u64());
        r.raw(g.data(), g.size() * 8);
        if (!g.empty()) {
            if (ev.gc_fg.empty()) ev.gc_fg.assign(g.size(), 0.0);
            if (g.size() != ev.gc_fg.size()) die("part " + path + ": GC histogram size differs");
            for (size_t i = 0; i < g.size(); ++i) ev.gc_fg[i] += g[i];
        }
    }
    fclose(f);
}

void write_idxstats(const std::string& path, const Model& M, const Evidence& ev) {
    ResultBin R;
    R.refs = M.refs;
    R.rlen.assign(M.lens.begin(), M.lens.end());
    R.mapped = ev.idx_mapped;
    R.unmapped = ev.idx_unmapped;
    R.unmapped_star = ev.idx_unmapped_star;
    write_text_file(path, [&](FILE* f) { write_idx_tsv(f, R); });
}

void log_records(const Evidence& ev, const std::string& note) {
    logm("[INFO] " + commas(ev.n_rec) + " records -> " + commas(ev.n_frag) + " fragments (" + commas(ev.n_unal) +
         " unaligned); multi-target classes " + commas((long long)ev.eq.size()) + " + weighted " +
         commas((long long)ev.eqx.size()) + note);
}

void log_lineage_stats(const LineageStats& st, const std::vector<std::string>& lins, long long n_pv) {
    auto pct = [](long long a, long long b) { return sfmt("%.2f%%", 100.0 * (double)a / (double)std::max(b, 1LL)); };
    long long nk = 0, ns = 0;
    for (int i = 0; i < 4; ++i) { nk += st.kmer_st[i]; ns += st.snp_st[i]; }
    logm("[INFO] lineage evidence: " + commas(st.valid) + " fragments on target transcripts");
    if (nk) {
        const char* K[4] = {"Specific", "Conflict_Weighted", "Weak_Evidence", "Ambiguous_None"};
        std::string s = "[INFO]   k-mer:";
        for (int i = 0; i < 4; ++i) s += std::string(" ") + K[i] + " " + commas(st.kmer_st[i]) + " (" + pct(st.kmer_st[i], nk) + ")";
        logm(s);
        s = "[INFO]   k-mer N_Hits:";
        for (int i = 0; i < 11; ++i) s += " " + (i < 10 ? std::to_string(i) : std::string(">=10")) + ":" + pct(st.hits[i], nk);
        logm(s);
    }
    if (ns) {
        const char* S[4] = {"Specific_Soft", "Multi_Lineage", "No_SNP_Covered", "Chimeric_Drop"};
        std::string s = "[INFO]   SNP (one row per fragment):";
        for (int i = 0; i < 4; ++i) s += std::string(" ") + S[i] + " " + commas(st.snp_st[i]) + " (" + pct(st.snp_st[i], ns) + ")";
        logm(s + "; transcripts evaluated " + commas(st.snp_rows));
    }
    if (st.both + st.kmer_only + st.snp_only) {
        const char* F[6] = {"Specific_Soft", "Multi_Lineage", "Specific_Kmer", "Multi_Kmer", "Specific_SNP", "Multi_SNP"};
        long long nf = st.both + st.kmer_only + st.snp_only;
        std::string s = "[INFO]   fused: both " + commas(st.both) + " | k-mer only " + commas(st.kmer_only) + " | SNP only " +
                        commas(st.snp_only) + ";";
        for (int i = 0; i < 6; ++i) if (st.fused_st[i]) s += std::string(" ") + F[i] + " " + pct(st.fused_st[i], nf);
        logm(s);
    }
    std::string ls = "[";
    for (size_t i = 0; i < lins.size(); ++i) ls += (i ? ", '" : "'") + lins[i] + "'";
    logm("[INFO] lineage P: " + ls + "]; reads with a non-zero weight: " + commas(n_pv));
}

// ================================================================== EM
// ================================================================== EM helpers (NUMA layout, deterministic sums)
namespace {
// array whose pages are not touched at allocation (first touch decides the memory node)
template <class T>
struct RawVec {
    std::unique_ptr<T[]> p;
    size_t n = 0;
    explicit RawVec(size_t k) : p(k ? new T[k] : nullptr), n(k) {}
    T& operator[](size_t i) { return p[i]; }
    const T& operator[](size_t i) const { return p[i]; }
};

// [cut[p], cut[p+1]) = block of thread p; blocks balanced by the cost prefix pre(i) (monotone, pre(n) = total)
template <class Fn>
std::vector<size_t> balanced_cuts(size_t n, int P, Fn pre) {
    std::vector<size_t> cut(P + 1, 0);
    cut[P] = n;
    const double tot = (double)pre(n);
    for (int p = 1; p < P; ++p) {
        const double target = tot * p / P;
        size_t lo = cut[p - 1], hi = n;
        while (lo < hi) {
            size_t mid = lo + (hi - lo) / 2;
            if ((double)pre(mid) < target) lo = mid + 1;
            else hi = mid;
        }
        cut[p] = lo;
    }
    return cut;
}

// thread p always gets block p (same thread, same memory, every iteration)
template <class Fn>
void par_blocks(const std::vector<size_t>& cut, Fn&& fn) {
    const int P = (int)cut.size() - 1;
#pragma omp parallel num_threads(P)
    {
        const int nth = omp_get_num_threads(), tid = omp_get_thread_num();
        for (int p = tid; p < P; p += nth) fn(p, cut[p], cut[p + 1]);
    }
}

// sum in fixed blocks of 8192, block sums added in order: independent of the number of threads
template <class Fn>
double det_sum(size_t n, Fn f) {
    const size_t B = 8192, nb = (n + B - 1) / B;
    std::vector<double> bs(nb, 0.0);
#pragma omp parallel for schedule(static)
    for (int64_t b = 0; b < (int64_t)nb; ++b) {
        double s = 0.0;
        const size_t e = std::min(n, (size_t)(b + 1) * B);
        for (size_t i = (size_t)b * B; i < e; ++i) s += f(i);
        bs[b] = s;
    }
    double s = 0.0;
    for (double x : bs) s += x;
    return s;
}

// EM threads pinned one per CPU of the job's CPU set (comput_opt: keep a thread next to its memory). Only inside the
// EM: lqcore never pins the aligner or the scan threads (pinning at start-up put the whole pipeline on one core).
struct EmPin {
    bool on = false;
    cpu_set_t orig;
    void pin(int nth) {
        if (getenv("LQCORE_NO_PIN") || nth < 2) return;
        CPU_ZERO(&orig);
        if (sched_getaffinity(0, sizeof orig, &orig) != 0) return;
        std::vector<int> cpus;
        for (int i = 0; i < CPU_SETSIZE; ++i) if (CPU_ISSET(i, &orig)) cpus.push_back(i);
        if ((int)cpus.size() < nth) {
            logm(sfmt("[INFO] EM: %d threads > %zu CPUs available, threads not pinned", nth, cpus.size()));
            return;
        }
        std::atomic<int> bad{0}, got{0};
#pragma omp parallel num_threads(nth)
        {
            if (omp_get_num_threads() == nth) {
                cpu_set_t s;
                CPU_ZERO(&s);
                CPU_SET(cpus[omp_get_thread_num()], &s);
                if (pthread_setaffinity_np(pthread_self(), sizeof s, &s) != 0) bad++;
                else got++;
            }
        }
        on = got.load() > 0;
        if (on) logm(sfmt("[INFO] EM: %d threads pinned to CPUs %d..%d, memory placed by first touch", got.load(), cpus[0], cpus[nth - 1]));
    }
    void unpin(int nth) {
        if (!on) return;
#pragma omp parallel num_threads(nth)
        { pthread_setaffinity_np(pthread_self(), sizeof orig, &orig); }
        on = false;
    }
};
}  // namespace

void solve_em(const Model& M, Evidence& ev, const Opts& a, double T0) {
    const uint32_t T = M.T;
    const bool NMW = a.nm_weight != "off";
    const long long N_al = ev.n_frag - ev.n_unal;
    double lam = 1.0;
    if (NMW) {
        if (a.nm_weight == "auto") {
            double p_err = std::min(std::max((double)ev.nm_sum / (double)std::max(ev.len_sum, 1LL), 1e-4), 0.1);
            lam = (p_err / 3.0) / (1.0 - p_err);
            logm(sfmt("[INFO] [G4] error rate p = %.5f (best-alignment NM / read length) -> x %.2e per extra mismatch", p_err, lam));
        } else {
            double v;
            if (!py_float(a.nm_weight, v)) die("--nm-weight must be auto|off|number");
            lam = v;
            logm(sfmt("[INFO] [G4] x %.2e per extra mismatch", lam));
        }
        long long nd = 0;
        for (size_t c = 0; c < ev.eqx.size(); ++c) {
            const std::string& k = *ev.eqx.keys[c];
            size_t n;
            memcpy(&n, k.data(), sizeof n);
            size_t off = sizeof n + n * 4;
            bool hp = k[off] != 0;
            off += 1 + (hp ? n * 8 : 0);
            if (k[off] != 0) nd += (long long)ev.eqx.cnt[c];
        }
        logm("[INFO] [G4] multi-mapping fragments whose candidates differ in mismatch count: " + commas(nd));
    }

    // ---- flatten classes: eq first, then eqx
    //   NUMA: the big arrays are allocated untouched and first written by the thread that later computes on that
    //   part (contiguous blocks balanced by entry count), so their pages sit on that thread's memory node.
    const size_t Ceq = ev.eq.size(), C = Ceq + ev.eqx.size();
    auto keyof = [&](size_t c) -> const std::string& { return c < Ceq ? *ev.eq.keys[c] : *ev.eqx.keys[c - Ceq]; };
    std::vector<uint64_t> coff(C + 1, 0);
    for (size_t c = 0; c < C; ++c) {
        size_t n;
        memcpy(&n, keyof(c).data(), sizeof n);
        coff[c + 1] = coff[c] + n;
    }
    const uint64_t E = coff[C];
    if (E >= 0xFFFFFFFFull) die("too many class entries for 32-bit indices");
    const int P = std::max(1, a.threads);
    EmPin pin;
    if (!a.no_pin) pin.pin(P);
    const std::vector<size_t> ccut = balanced_cuts(C, P, [&](size_t c) { return coff[c] + c; });
    RawVec<uint32_t> I(E);
    RawVec<double> F(E), w(E), cnt(C);
    par_blocks(ccut, [&](int, size_t c0, size_t c1) {
        for (size_t c = c0; c < c1; ++c) {
            const std::string& k = keyof(c);
            size_t n;
            memcpy(&n, k.data(), sizeof n);
            const char* p = k.data() + sizeof n;
            uint64_t o = coff[c];
            for (size_t i = 0; i < n; ++i) { memcpy(&I[o + i], p + 4 * i, 4); w[o + i] = 0.0; }
            if (c < Ceq) {
                for (size_t i = 0; i < n; ++i) F[o + i] = 1.0;
                cnt[c] = (double)ev.eq.cnt[c];
            } else {
                size_t off = sizeof n + n * 4;
                bool hp = k[off] != 0;
                off += 1;
                for (size_t i = 0; i < n; ++i) {
                    double wv = 1.0;
                    if (hp) memcpy(&wv, k.data() + off + 8 * i, 8);
                    F[o + i] = wv;
                }
                off += hp ? n * 8 : 0;
                bool hd = k[off] != 0;
                off += 1;
                if (hd)
                    for (size_t i = 0; i < n; ++i) {
                        int d;
                        memcpy(&d, k.data() + off + 4 * i, 4);
                        F[o + i] = F[o + i] * std::pow(lam, (double)d);
                    }
                cnt[c] = (double)ev.eqx.cnt[c - Ceq];
            }
        }
    });
    { KeyMap e1, e2; std::swap(ev.eq, e1); std::swap(ev.eqx, e2); }   // free the class maps

    if (C) {   // aligned-fragment breakdown (log only)
        double s_single = np_sum(ev.single);
        std::vector<double> v_wc, v_xc, v_bo, v_ao;
        for (size_t c = 0; c < C; ++c) {
            bool ha = false, ho = false;
            uint64_t mn = (uint64_t)1 << 40;
            long long mx = -1;
            for (uint64_t e = coff[c]; e < coff[c + 1]; ++e) {
                uint32_t t = I[e];
                if (M.is_al[t]) ha = true; else ho = true;
                long long ci = M.cl_int[t] == NOCL ? -1 : (long long)M.cl_int[t];
                mn = std::min<uint64_t>(mn, ci < 0 ? ((uint64_t)1 << 40) : (uint64_t)ci);
                mx = std::max(mx, ci);
            }
            bool same = (long long)mn == mx;
            if (ha && !ho && same) v_wc.push_back(cnt[c]);
            else if (ha && !ho && !same) v_xc.push_back(cnt[c]);
            else if (!ha && ho) v_bo.push_back(cnt[c]);
            else if (ha && ho) v_ao.push_back(cnt[c]);
        }
        double den = (double)std::max(N_al, 1LL);
        logm(sfmt("[INFO] aligned fragments: unique %.2f%%  within-cluster multi %.2f%%  cross-cluster multi %.2f%%  "
                  "between orphans %.2f%%  allele+orphan %.2f%%",
                  100 * s_single / den, 100 * np_sum(v_wc) / den, 100 * np_sum(v_xc) / den,
                  100 * np_sum(v_bo) / den, 100 * np_sum(v_ao) / den));
    }

    // ---- effective length
    double mu = 0.0;
    // 1.1.0 --eff-len fld: length - mean of the empirical fragment-length distribution truncated at the length + 1
    //   (kallisto / Salmon); S0 / S1 = cumulative count / length-weighted count of the observed template lengths
    const int FLMAX = 2000;
    bool use_fld = false;
    std::vector<double> fl_s0, fl_s1;
    if (a.eff_len == "fld") {
        if (ev.fl_vals.size() >= 1000) {
            use_fld = true;
            std::vector<double> h(FLMAX, 0.0);
            for (int32_t v : ev.fl_vals) if (v > 0 && v < FLMAX) h[v] += 1.0;
            fl_s0.assign(FLMAX, 0.0);
            fl_s1.assign(FLMAX, 0.0);
            for (int x = 1; x < FLMAX; ++x) { fl_s0[x] = fl_s0[x - 1] + h[x]; fl_s1[x] = fl_s1[x - 1] + h[x] * (double)x; }
            std::vector<int32_t> srt(ev.fl_vals.begin(), ev.fl_vals.end());
            std::sort(srt.begin(), srt.end());
            logm(sfmt("[INFO] fragment-length distribution: %s template lengths, mean %.1f, median %d, 2.5-97.5%% %d-%d bp "
                      "-> effective length = length - truncated mean + 1",
                      commas((long long)srt.size()).c_str(), fl_s1[FLMAX - 1] / std::max(fl_s0[FLMAX - 1], 1.0), srt[srt.size() / 2],
                      srt[(size_t)(0.025 * srt.size())], srt[std::min(srt.size() - 1, (size_t)(0.975 * srt.size()))]));
        } else logm("[WARN] --eff-len fld: fewer than 1,000 proper-pair template lengths -> mean-based effective length (eff)");
    }
    if ((a.eff_len == "eff" || a.eff_len == "fld") && !use_fld) {
        if (a.frag_len) mu = a.frag_len;
        else if (ev.fl_vals.size() >= 1000) {
            long long s = 0;
            for (auto v : ev.fl_vals) s += v;
            mu = (double)s / (double)ev.fl_vals.size();
        }
        if (mu <= 0) logm("[WARN] fragment length could not be estimated -> using full transcript length");
        else logm(sfmt("[INFO] mean fragment length %.1f bp", mu));
    }
    std::vector<double> eff(T);
    for (uint32_t t = 0; t < T; ++t) {
        if (a.eff_len == "off") eff[t] = 1.0;
        else if (use_fld) {
            const double L = M.tlen[t];
            const int m = (int)std::min<double>(L, (double)(FLMAX - 1));
            const double mt = (m >= 1 && fl_s0[m] > 0) ? fl_s1[m] / fl_s0[m] : L;
            eff[t] = std::max(L - mt + 1.0, a.min_eff_len);
        } else eff[t] = std::max(mu > 0 ? M.tlen[t] - mu + 1.0 : M.tlen[t], a.min_eff_len);
    }
    std::vector<double> gc_factor;          // 1.1.0 --gc-bias: per-target factor applied to the effective length

    // ---- initial values + prior
    const std::vector<double>& single = ev.single;
    std::vector<double> th0(T), uni(T), alpha(T), prior(T, 0.0);
    for (uint32_t i = 0; i < T; ++i) {
        if (M.is_al[i]) {
            auto it = M.theta0.find(M.names[i]);
            th0[i] = it == M.theta0.end() ? 0.0 : it->second;
        } else th0[i] = 1.0;
    }
    {
        std::vector<double> acc(T, 0.0);
        for (size_t c = 0; c < C; ++c) {
            double wgt = cnt[c] / (double)(coff[c + 1] - coff[c]);
            for (uint64_t e = coff[c]; e < coff[c + 1]; ++e) acc[I[e]] += wgt;
        }
        for (uint32_t t = 0; t < T; ++t) uni[t] = single[t] + acc[t];
    }
    alpha = uni;
    const size_t NC = M.cl_sorted.size();
    if (M.n_al) {
        std::vector<double> Z(NC, 0.0), S(NC, 0.0);
        for (uint32_t t = 0; t < T; ++t) if (M.is_al[t]) { Z[M.cl_int[t]] += th0[t]; S[M.cl_int[t]] += uni[t]; }
        for (uint32_t t = 0; t < T; ++t) {
            if (!M.is_al[t]) continue;
            double z = Z[M.cl_int[t]];
            double tb = z > 0 ? S[M.cl_int[t]] * th0[t] / z : uni[t];
            if (a.init == "theta") alpha[t] = 0.5 * uni[t] + 0.5 * tb;
        }
        std::vector<double> tau_c(NC);
        bool any_tau = false;
        long long n_tau = 0;
        for (size_t j = 0; j < NC; ++j) {
            const std::string& cid = M.cl_sorted[j];
            auto ps = M.prior_status.find(cid);
            double tv;
            if (ps != M.prior_status.end() && ps->second == "LongRead_Anchored") {
                auto lc = M.lr_count.find(cid);
                tv = a.tp_lr * (lc == M.lr_count.end() ? 0.0 : lc->second);
            } else tv = a.tp_kmer;
            tv = std::max(a.tp_cap > 0 ? std::min(tv, a.tp_cap) : tv, 0.0);
            tau_c[j] = tv;
            if (tv != 0.0) { any_tau = true; ++n_tau; }
        }
        if (any_tau) {
            for (uint32_t t = 0; t < T; ++t) {
                if (!M.is_al[t]) continue;
                double z = Z[M.cl_int[t]];
                double thn = z > 0 ? th0[t] / std::max(z, 1e-300) : 0.0;
                prior[t] = tau_c[M.cl_int[t]] * thn;
            }
            logm("[INFO] theta0 prior: " + commas(n_tau) + " clusters with tau > 0");
        }
    }
    if (!M.pool_reps.empty()) {
        std::vector<double> buf;
        for (uint32_t r : M.pool_reps) {
            auto& ms = M.POOLG[r];
            buf.clear(); for (uint32_t m : ms) buf.push_back(alpha[m]);
            double sa = np_sum(buf);
            buf.clear(); for (uint32_t m : ms) buf.push_back(prior[m]);
            double sp = np_sum(buf);
            alpha[r] = sa; prior[r] = sp;
            for (uint32_t m : ms) if (m != r) { alpha[m] = 0.0; prior[m] = 0.0; }
        }
    }

    // ---- target -> entries (increasing entry order)
    std::vector<uint64_t> toff(T + 1, 0);
    for (uint64_t e = 0; e < E; ++e) toff[I[e] + 1]++;
    for (uint32_t t = 0; t < T; ++t) toff[t + 1] += toff[t];
    const std::vector<size_t> tcut = balanced_cuts(T, P, [&](size_t t) { return toff[t] + t; });
    RawVec<uint32_t> tent(E);
    par_blocks(tcut, [&](int, size_t t0, size_t t1) { for (uint64_t k = toff[t0]; k < toff[t1]; ++k) tent[k] = 0; });
    {
        std::vector<uint64_t> fill(toff.begin(), toff.end() - 1);
        for (uint64_t e = 0; e < E; ++e) tent[fill[I[e]]++] = (uint32_t)e;
    }

    // ---- one EM step: phase A per class (E-step weights), phase B per target in entry order (M-step).
    //      Each class / target is computed by one thread with the same operations as em_global.py, so the result
    //      does not depend on the number of threads or on the partition.
    std::vector<double> ae(T), rpart(P);
    std::vector<double> clog;                          // SQUAREM: cnt_c * log(den_c), objective of the step's input
    auto step = [&](const double* ain, double* aout, double* cout, char* bigf, bool want_obj, double& obj) -> double {
        par_blocks(tcut, [&](int, size_t t0, size_t t1) { for (size_t t = t0; t < t1; ++t) ae[t] = ain[t] / eff[t]; });
        if (C) {
            par_blocks(ccut, [&](int, size_t c0, size_t c1) {
                for (size_t c = c0; c < c1; ++c) {
                    const uint64_t b = coff[c], e2 = coff[c + 1];
                    double den = 0.0;
                    for (uint64_t e = b; e < e2; ++e) {
                        double num = ae[I[e]] * F[e];
                        w[e] = num;
                        den += num;
                    }
                    const double d = std::max(den, 1e-300), nk = cnt[c];
                    for (uint64_t e = b; e < e2; ++e) w[e] = w[e] / d * nk;
                    if (want_obj) clog[c] = nk * std::log(d);
                }
            });
        }
        par_blocks(tcut, [&](int p, size_t t0, size_t t1) {
            double r = 0.0;
            for (size_t t = t0; t < t1; ++t) {
                double acc = 0.0;
                for (uint64_t k = toff[t]; k < toff[t + 1]; ++k) acc += w[tent[k]];
                double ct = single[t] + acc;
                cout[t] = ct;
                double nw = ct + prior[t];
                if (nw > 1.0) {
                    double x = std::fabs(nw - ain[t]) / nw;
                    if (x > r) r = x;
                    bigf[t] = 1;
                } else bigf[t] = 0;
                aout[t] = nw;
            }
            rpart[p] = r;
        });
        double r = 0.0;
        for (double x : rpart) r = std::max(r, x);
        if (want_obj) {
            // log-likelihood of the step's input (prior = extra unique fragments): plain EM never decreases it
            double s_cl = det_sum(C, [&](size_t c) { return clog[c]; });
            double s_si = det_sum(T, [&](size_t t) {
                double n = single[t] + prior[t];
                return n > 0 ? n * std::log(ae[t]) : 0.0;
            });
            double z = det_sum(T, [&](size_t t) { return ae[t]; });
            double ntot = det_sum(T, [&](size_t t) { return single[t] + prior[t]; }) + det_sum(C, [&](size_t c) { return cnt[c]; });
            obj = s_cl + s_si - ntot * std::log(z);
            if (!std::isfinite(obj)) obj = -INFINITY;
        }
        return r;
    };

    // ---- 1.1.0 --gc-bias (Salmon-style fragment GC model, one GC model, no sequence-context bias):
    //   foreground = GC percent of the observed fragments (primary proper pair, collected by `run`);
    //   background = GC percent of every possible fragment of every target (FLD quantile lengths), weighted by the target's
    //   expected fragments (alpha) after a first EM pass; w = foreground / background (Laplace-smoothed, clamped to
    //   [0.01, 100], normalised to mean 1 over the background); effective length x mean w over the target's fragments;
    //   then the EM continues from the first-pass solution until it converges again.
    auto gc_correct = [&]() -> bool {
        double nfg = 0.0;
        for (double v : ev.gc_fg) nfg += v;
        if (ev.gc_fg.size() != 101 || nfg < 1000) {
            logm("[WARN] [GC] no fragment GC histogram in the evidence (lqcore run --gc-bias / --gc-hist) -> GC correction skipped");
            return false;
        }
        if (a.fasta.empty() || !file_exists(a.fasta)) die("--gc-bias needs --fasta (target sequences)");
        if (ev.fl_vals.size() < 1000) { logm("[WARN] [GC] fewer than 1,000 template lengths -> GC correction skipped"); return false; }
        double tg = wall();
        RefGC G;
        load_ref_gc(G, a.fasta, M.refs, M.r2i);
        std::vector<uint32_t> rep(T, NOCL);
        for (size_t r = 0; r < M.refs.size(); ++r) {
            if (!G.ok((uint32_t)r)) continue;
            const uint32_t t = M.ref2t[r];
            if (rep[t] == NOCL || G.len[r] > G.len[rep[t]]) rep[t] = (uint32_t)r;
        }
        std::vector<int32_t> fl(ev.fl_vals.begin(), ev.fl_vals.end());
        std::sort(fl.begin(), fl.end());
        const int K = a.gc_fl_samples;
        std::vector<int> ks(K);
        for (int i = 0; i < K; ++i) ks[i] = std::max(1, fl[std::min(fl.size() - 1, (size_t)(((double)i + 0.5) / K * fl.size()))]);
        std::vector<uint32_t> act;
        for (uint32_t t = 0; t < T; ++t) if (alpha[t] > 0 && rep[t] != NOCL) act.push_back(t);
        const size_t NA = act.size();
        std::vector<float> H(NA * 101, 0.0f);
        std::vector<double> Hs(NA, 0.0);
#pragma omp parallel for schedule(dynamic, 64) num_threads(P)
        for (int64_t q = 0; q < (int64_t)NA; ++q) {
            const uint32_t r = rep[act[q]];
            const uint64_t Lr = (uint64_t)G.len[r];
            float* h = &H[(size_t)q * 101];
            double tot = 0.0;
            for (int k : ks) {
                if ((uint64_t)k > Lr) continue;
                uint64_t g = G.count(r, 0, (uint64_t)k);
                const uint64_t kk = (uint64_t)k;
                h[(200 * g + kk) / (2 * kk)] += 1.0f;
                for (uint64_t j = 1; j + kk <= Lr; ++j) {
                    g += (uint64_t)G.bit(r, j + kk - 1);
                    g -= (uint64_t)G.bit(r, j - 1);
                    h[(200 * g + kk) / (2 * kk)] += 1.0f;
                }
                tot += (double)(Lr - kk + 1);
            }
            Hs[q] = tot;
        }
        std::vector<double> bg(101, 0.0);
        for (size_t q = 0; q < NA; ++q) {
            if (Hs[q] <= 0) continue;
            const double wq = alpha[act[q]] / Hs[q];
            const float* h = &H[q * 101];
            for (int b = 0; b <= 100; ++b) bg[b] += wq * (double)h[b];
        }
        double nbg = 0.0;
        for (double v : bg) nbg += v;
        if (nbg <= 0) { logm("[WARN] [GC] empty background -> GC correction skipped"); return false; }
        std::vector<double> w(101);
        for (int b = 0; b <= 100; ++b) {
            const double fp = (ev.gc_fg[b] + 1.0) / (nfg + 101.0), bp = (bg[b] / nbg * nfg + 1.0) / (nfg + 101.0);
            w[b] = std::min(100.0, std::max(0.01, fp / bp));
        }
        double mw = 0.0;
        for (int b = 0; b <= 100; ++b) mw += bg[b] / nbg * w[b];
        for (auto& x : w) x /= mw;
        gc_factor.assign(T, 1.0);
        for (size_t q = 0; q < NA; ++q) {
            if (Hs[q] <= 0) continue;
            const float* h = &H[q * 101];
            double s = 0.0;
            for (int b = 0; b <= 100; ++b) s += (double)h[b] * w[b];
            gc_factor[act[q]] = s / Hs[q];
        }
        for (uint32_t r : M.pool_reps) for (uint32_t m : M.POOLG[r]) gc_factor[m] = gc_factor[r];
        for (uint32_t t = 0; t < T; ++t) eff[t] = std::max(eff[t] * gc_factor[t], a.min_eff_len);
        std::vector<double> fs;
        for (uint32_t t : act) fs.push_back(gc_factor[t]);
        std::sort(fs.begin(), fs.end());
        auto qf = [&](double p) { return fs.empty() ? 1.0 : fs[std::min(fs.size() - 1, (size_t)(p * fs.size()))]; };
        std::string ws;
        for (int b : {20, 30, 40, 50, 60, 70}) ws += sfmt(" %d%%:%.2f", b, w[b]);
        logm(sfmt("[INFO] [GC] %s fragments, %s expressed targets, %d fragment lengths (%d-%d bp) | weight by fragment GC",
                  commas_f0(nfg).c_str(), commas((long long)NA).c_str(), K, ks.front(), ks.back()) + ws);
        logm(sfmt("[INFO] [GC] effective-length factor: 5%% %.3f  median %.3f  95%% %.3f  (%.1fs)", qf(0.05), qf(0.5), qf(0.95),
                  wall() - tg));
        return true;
    };

    double t_em = wall();
    std::vector<double> counts(single), anew(T);
    int it = 0, it_total = 0;
    double rel = 0.0;
    std::vector<char> big(T, 0);
    long long n_cycles = 0, n_accept = 0;
    const int n_rounds = a.gc_bias ? 2 : 1;
    for (int round = 0; round < n_rounds; ++round) {
    if (round == 1) {
        logm(sfmt("[INFO] [GC] first EM pass: %d iterations (last relative change %.2e)", it, rel));
        if (!gc_correct()) break;
        n_cycles = n_accept = 0;
    }
    if (a.em_accel == "none") {
        double dummy = 0;
        for (it = 1; it <= a.max_iter; ++it) {
            rel = step(alpha.data(), anew.data(), counts.data(), big.data(), false, dummy);
            alpha.swap(anew);
            if (it >= a.min_iter && rel < a.tol) break;
        }
        if (it > a.max_iter) it = a.max_iter;
    } else {
        // SQUAREM (Varadhan & Roland 2008, scheme 3); extrapolated values are floored at 1e-6 x the plain-EM value
        // (a target is never set to 0, which EM could not undo).
        // Cycle: p1 = F(p), p2 = F(p1), extrapolate to pe, stabilise pn = F(pe), continue from pn. As in the SQUAREM
        // package, pn is kept only if its log-likelihood is not below that of p (tolerance 1); the check uses the
        // likelihood that the next F(pn) computes anyway, and on failure the cycle restarts from p2 (plain EM, which
        // never lowers the likelihood). All sums use fixed blocks: the path does not depend on the thread count.
        clog.assign(C, 0.0);
        const bool logsp = getenv("LQ_SQUAREM_LOG") != nullptr;   // log-space extrapolation (slower here; kept for tests)
        std::vector<double> p1(T), p2(T), pe(T), pn(T), c1(T), c2(T), cn(T), prev(T), fb(T), fbc(T), fbprev(T);
        std::vector<char> b1(T), b2(T), bn(T), fbb(T);
        double fbrel = 0, l_ref = -INFINITY, st_t = 1;
        bool tentative = false;
        const double stepmax0 = 1.0, mstep = 4.0, inc = 1.0;
        double stepmax = stepmax0;
        int fe = 0;
        auto tr = [&](double x) { return logsp ? std::log(x) : x; };
        while (fe < a.max_iter) {
            double l0 = 0, l1 = 0, le = 0;
            double r1 = step(alpha.data(), p1.data(), c1.data(), b1.data(), true, l0); ++fe;
            if (tentative) {
                tentative = false;
                if (!(l0 >= l_ref - inc)) {          // extrapolation lowered the likelihood: back to plain EM's p2
                    if (st_t == stepmax) stepmax = std::max(stepmax0, stepmax / mstep);
                    alpha.swap(fb); counts.swap(fbc); big.swap(fbb); prev.swap(fbprev); rel = fbrel;
                    --n_accept;
                    if (fe >= a.min_iter && rel < a.tol) break;
                    continue;
                }
                if (st_t == stepmax) stepmax = mstep * stepmax;
            }
            if ((fe >= a.min_iter && r1 < a.tol) || fe >= a.max_iter) { prev = alpha; alpha.swap(p1); counts.swap(c1); big.swap(b1); rel = r1; break; }
            double r2 = step(p1.data(), p2.data(), c2.data(), b2.data(), true, l1); ++fe;
            ++n_cycles;
            auto ok = [&](size_t t) { return alpha[t] > 0 && p1[t] > 0 && p2[t] > 0; };
            double sr2 = det_sum(T, [&](size_t t) { if (!ok(t)) return 0.0; double q = tr(p1[t]) - tr(alpha[t]); return q * q; });
            double sv2 = det_sum(T, [&](size_t t) {
                if (!ok(t)) return 0.0;
                double v = tr(p2[t]) - 2.0 * tr(p1[t]) + tr(alpha[t]);
                return v * v;
            });
            double st = sv2 > 0 ? std::sqrt(sr2 / sv2) : 1.0;
            st = std::max(1.0, std::min(stepmax, st));
            bool ext = false;
            if (std::fabs(st - 1.0) > 0.01 && fe < a.max_iter && !(fe + 1 >= a.min_iter && r2 < a.tol)) {
                bool fin = true;
                par_blocks(tcut, [&](int, size_t t0, size_t t1) {
                    for (size_t t = t0; t < t1; ++t) {
                        if (!ok(t)) { pe[t] = p2[t]; continue; }
                        double l = tr(alpha[t]), q = tr(p1[t]) - l, v = tr(p2[t]) - 2.0 * tr(p1[t]) + l;
                        double x = l + 2.0 * st * q + st * st * v;
                        pe[t] = logsp ? std::exp(x) : std::max(x, 1e-6 * p2[t]);
                    }
                });
                for (uint32_t t = 0; t < T && fin; ++t) fin = std::isfinite(pe[t]);
                if (fin) {
                    double rn = step(pe.data(), pn.data(), cn.data(), bn.data(), false, le); ++fe;
                    // keep p2 as the fallback, continue from pn; confirmed by the next cycle's likelihood
                    fb.swap(p2); fbc.swap(c2); fbb.swap(b2); fbprev = p1; fbrel = r2;
                    prev = pe; alpha.swap(pn); counts.swap(cn); big.swap(bn); rel = rn;
                    l_ref = l0; st_t = st; tentative = true; ext = true; ++n_accept;
                } else if (st == stepmax) stepmax = std::max(stepmax0, stepmax / mstep);
            }
            if (!ext) {
                prev = p1; alpha.swap(p2); counts.swap(c2); big.swap(b2); rel = r2;
                if (st == stepmax && std::fabs(st - 1.0) > 0.01) stepmax = mstep * stepmax;
                else if (std::fabs(st - 1.0) <= 0.01) stepmax = std::max(stepmax, mstep * stepmax0);
            }
            if (!tentative && fe >= a.min_iter && rel < a.tol) break;
        }
        it = fe;
        anew = prev;
        logm(sfmt("[INFO] EM acceleration: SQUAREM (%s, monotone) | %lld cycles, %lld extrapolations kept, %d EM steps",
                  logsp ? "log space" : "linear", n_cycles, n_accept, fe));
    }
    it_total += it;
    }   // rounds
    if (n_rounds > 1) it = it_total;
    pin.unpin(P);
    if (rel < a.tol) {
        logm("[INFO] EM converged: " + std::to_string(it) + " iterations (last relative change " + sfmt("%.2e", rel) +
             sfmt("), %.1fs; assigned ", wall() - t_em) + commas_f0(np_sum(counts)) + " / aligned " + commas(N_al) + " fragments");
    } else {
        long long n_big = 0, n_slow = 0;
        std::vector<double> slow_vals;
        for (uint32_t t = 0; t < T; ++t) {
            if (!big[t]) continue;
            ++n_big;
            double x = std::fabs(alpha[t] - anew[t]) / alpha[t];
            if (x >= a.tol) { ++n_slow; slow_vals.push_back(alpha[t]); }
        }
        logm("[INFO] EM stopped at --max-iter " + std::to_string(a.max_iter) + " before reaching --tol " + pyfloat(a.tol) +
             " (last relative change " + sfmt("%.2e", rel) + sfmt("), %.1fs; assigned ", wall() - t_em) +
             commas_f0(np_sum(counts)) + " / aligned " + commas(N_al) + " fragments");
        logm("[INFO] targets still changing by >= " + pyfloat(a.tol) + " per iteration: " + commas(n_slow) + " of " +
             commas(n_big) + " with > 1 fragment, holding " +
             sfmt("%.3f%%", 100.0 * np_sum(slow_vals) / std::max(np_sum(alpha), 1e-300)) + " of the fragments");
    }

    // ---- unpool, final theta, output
    for (uint32_t r : M.pool_reps) {
        auto& ms = M.POOLG[r];
        double vc = counts[r] / (double)ms.size(), va = alpha[r] / (double)ms.size();
        for (uint32_t m : ms) { counts[m] = vc; alpha[m] = va; }
    }
    std::vector<double> fin(T, 1.0);
    if (M.n_al) {
        std::vector<double> Sx(NC, 0.0);
        for (uint32_t t = 0; t < T; ++t) if (M.is_al[t]) Sx[M.cl_int[t]] += alpha[t];
        for (uint32_t t = 0; t < T; ++t) {
            if (!M.is_al[t]) continue;
            double z = Sx[M.cl_int[t]];
            fin[t] = z > 0 ? alpha[t] / std::max(z, 1e-300) : 0.0;
        }
    }
    ResultBin R;
    R.iter = it;
    R.cluster.resize(T); R.allele.resize(T); R.ancestry.resize(T);
    R.th0.resize(T); R.counts.resize(T); R.fin.resize(T); R.tlen.resize(T); R.eff.resize(T); R.gc.resize(T);
    for (uint32_t i = 0; i < T; ++i) {
        const std::string& n = M.names[i];
        R.cluster[i] = M.has_cl[i] ? M.tcl[i] : std::string("ORPHAN");
        R.allele[i] = n;
        auto ai = M.anc.find(n);
        if (ai != M.anc.end()) R.ancestry[i] = ai->second;
        else {
            auto aj = M.anc.find(std::string(rsplit_dot(n)));
            R.ancestry[i] = aj != M.anc.end() ? aj->second : (M.has_cl[i] ? "UNKNOWN" : "ORPHAN");
        }
        R.th0[i] = th0[i]; R.counts[i] = counts[i]; R.fin[i] = fin[i];
        R.tlen[i] = M.tlen[i]; R.eff[i] = eff[i]; R.gc[i] = gc_factor.empty() ? 1.0 : gc_factor[i];
    }
    if (ends_with(a.out, ".lqres")) {   // 1.1.0: one binary file per sample (EM + effective lengths + per-reference counts)
        R.refs = M.refs;
        R.rlen.assign(M.lens.begin(), M.lens.end());
        R.mapped = ev.idx_mapped;
        R.unmapped = ev.idx_unmapped;
        R.unmapped_star = ev.idx_unmapped_star;
        write_result_bin(a.out, R);
    } else write_text_file(a.out, [&](FILE* f) { write_em_tsv(f, R); });
    if (!a.eff_out.empty()) {   // 1.1.0: effective lengths the EM used (TPM = reads / effective length)
        write_text_file(a.eff_out, [&](FILE* f) { write_eff_tsv(f, R); });
        logm("[INFO] effective lengths: " + a.eff_out);
    }
    logm("[SUCCESS] EM complete: " + a.out + sfmt(" (%.0fs)", wall() - T0));
}

}  // namespace lq
