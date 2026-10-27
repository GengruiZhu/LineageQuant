// lqcore -- model (targets, clusters, pooling), evidence (equivalence classes) and the EM solver
#pragma once
#include "common.hpp"

#include <map>
#include <memory>
#include <unordered_map>
#include <vector>

namespace lq {

static const uint32_t NOCL = 0xFFFFFFFFu;
static const int NM_INF = 1 << 30;

// ------------------------------------------------------------------ options (shared by all subcommands)
struct Opts {
    // EM (same names and defaults as core/em_global.py)
    std::string theta, pvals, ancestry, out = "EM_Final_Expression.tsv", bam, sam, fasta;
    bool no_pool = false, keep_collate = false, fld_weight = false;
    int threads = 16, workers = 0;
    std::string eff_len = "eff";             // off | len | eff (length - mean fragment + 1) | fld (truncated empirical FLD, 1.1.0)
    double frag_len = 0.0, min_eff_len = 25.0;
    bool gc_bias = false;                    // 1.1.0: fragment GC-bias correction of the effective lengths (needs --fasta)
    bool gc_hist = false;                    // run: collect the fragment GC histogram into the evidence (for a later em --gc-bias)
    int gc_fl_samples = 21;                  // fragment lengths (FLD quantiles) used for the expected GC distribution
    std::string eff_out;                     // 1.1.0: per-target effective lengths actually used by the EM (TSV)
    int max_iter = 5000, min_iter = 50;
    double tol = 1e-3;
    std::string init = "theta";
    double pw_spec = 1.0, pw_conf = 0.0, pw_none = 0.0, p_alpha = 0.0, p_floor = 1e-3;
    double tp_lr = 1.0, tp_kmer = 0.0, tp_cap = 0.0;
    std::string nm_weight = "auto";
    std::string em_accel = "none";           // none (= em_global.py) | squarem (log-space SQUAREM, monotone)
    bool no_pin = false;                     // EM: do not pin OpenMP threads to CPUs
    std::vector<std::string> parts;          // em --parts: merge evidence of read chunks
    std::string idxstats_out, write_pvals;   // per-reference counts (samtools idxstats format); fused P table
    // run: alignment (bowtie2 | lqmap: built-in exact-seed + ungapped verification | sa: Salmon selective alignment)
    std::string aligner = "bowtie2";
    std::string r1, r2, index, bowtie2 = "bowtie2", bt2_opts, bt2_log, part_out, tmpdir;
    int map_k = 21, map_stride = 8, map_max_extra = 3, map_max_occ = 4000, map_max_cand = 1000, map_max_frag = 1000;
    double map_mm_frac = 0.1;                // per mate: mismatches (N not counted) <= frac * read length (= --score-min L,0,-0.1)
    bool map_no_rescue = false;
    std::string map_cache;                   // lqmap index cache file
    std::string salmon = "salmon", salmon_index, salmon_opts;
    int max_aln = 200, bt2_threads = 0;
    int chunk_i = 1, chunk_n = 1;
    long long block_pairs = 100000;
    // run: lineage
    std::string target_ids, marker_cache, snp_db, snp_pick = "best";
    std::vector<std::string> markers;        // NAME:PATH
    int kmer_len = 31, min_hits = 3;
    double kmer_alpha = 1.0, spec_thresh = 0.90, error_rate = 0.01, poe_min_weight = 0.05, poe_eps = 1e-6;
    bool no_chimeric = false, no_kmer = false;
    bool marker_clean = false;               // 1.1.0: drop marker k-mers found in reference copies labelled with another lineage
    std::string debug_kmer, debug_snp, debug_frag;
    std::vector<std::string> unknown;
};
Opts parse_opts(int argc, char** argv, void (*usage)());

// ------------------------------------------------------------------ insertion-ordered key -> count map
struct KeyMap {
    std::unordered_map<std::string, uint32_t> idx;
    std::vector<const std::string*> keys;
    std::vector<uint64_t> cnt;
    KeyMap() { idx.reserve(1 << 16); }
    void add(const std::string& k, uint64_t c) {
        auto it = idx.find(k);
        if (it == idx.end()) {
            auto r = idx.emplace(k, (uint32_t)keys.size());
            keys.push_back(&r.first->first);
            cnt.push_back(c);
        } else cnt[it->second] += c;
    }
    size_t size() const { return keys.size(); }
};

// ------------------------------------------------------------------ lineage bookkeeping (logs only)
struct LineageStats {
    long long valid = 0, both = 0, kmer_only = 0, snp_only = 0, pv_used = 0;
    long long kmer_st[4] = {0, 0, 0, 0};      // Specific, Conflict_Weighted, Weak_Evidence, Ambiguous_None
    long long hits[11] = {0};                  // N_Hits 0..9, >=10
    long long snp_st[4] = {0, 0, 0, 0};       // Specific_Soft, Multi_Lineage, No_SNP_Covered, Chimeric_Drop
    long long fused_st[6] = {0, 0, 0, 0, 0, 0};  // Specific_Soft, Multi_Lineage, Specific_Kmer, Multi_Kmer, Specific_SNP, Multi_SNP
    long long snp_rows = 0;                    // (fragment, transcript) SNP evaluations
    void add(const LineageStats& o);
};

// ------------------------------------------------------------------ evidence collected from the reads
struct Evidence {
    std::vector<double> single;            // per target
    KeyMap eq, eqx;                         // multi-target classes (unweighted / weighted), first-appearance order
    long long n_rec = 0, n_frag = 0, n_unal = 0, nm_sum = 0, len_sum = 0;
    std::vector<int32_t> fl_vals;           // first 1,000,000 proper-pair template lengths (0 < tl < 2000)
    std::vector<uint64_t> idx_mapped, idx_unmapped;   // per reference (samtools idxstats semantics)
    uint64_t idx_unmapped_star = 0;
    LineageStats ls;
    long long n_pv_reads = 0;               // reads with lineage P and a non-zero weight (logs)
    std::vector<double> gc_fg;              // 1.1.0: fragments per GC percent (0..100) of the primary alignment; empty = not collected
    void init(size_t T, size_t R) {
        single.assign(T, 0.0);
        idx_mapped.assign(R, 0);
        idx_unmapped.assign(R, 0);
        fl_vals.reserve(1000000);
    }
};

// ------------------------------------------------------------------ GC of the reference sequences (1 bit per base: G/C)
struct RefGC {
    std::vector<uint64_t> bits;             // all references concatenated
    std::vector<uint64_t> off;              // bit offset of reference r (UINT64_MAX = sequence not in the FASTA)
    std::vector<long long> len;
    bool ok(uint32_t r) const { return r < off.size() && off[r] != ~0ULL; }
    // number of G/C in [s, s + n) of reference r
    uint32_t count(uint32_t r, uint64_t s, uint64_t n) const {
        uint64_t a = off[r] + s, b = a + n, c = 0;
        while (a < b) {
            const uint64_t w = a >> 6, sh = a & 63, take = std::min<uint64_t>(64 - sh, b - a);
            uint64_t x = bits[w] >> sh;
            if (take < 64) x &= (1ULL << take) - 1;
            c += (uint64_t)__builtin_popcountll(x);
            a += take;
        }
        return (uint32_t)c;
    }
    bool bit(uint32_t r, uint64_t i) const { uint64_t a = off[r] + i; return (bits[a >> 6] >> (a & 63)) & 1; }
};
// reads FASTA sequences of the given references (by name); references missing from the FASTA get off = UINT64_MAX
void load_ref_gc(RefGC& g, const std::string& fasta, const std::vector<std::string>& refs,
                 const std::unordered_map<std::string, uint32_t>& r2i);

// ------------------------------------------------------------------ model
struct Model {
    std::vector<std::string> refs;
    std::vector<long long> lens;
    std::unordered_map<std::string, uint32_t> r2i;
    // Theta / ancestry
    std::unordered_map<std::string, std::string> g2c;
    std::vector<std::string> g2c_order;
    std::unordered_map<std::string, double> theta0, lr_count;
    std::unordered_map<std::string, std::string> prior_status;
    std::unordered_map<std::string, std::string> anc;
    // targets
    uint32_t T = 0;
    std::vector<std::string> names, tcl;
    std::vector<double> tlen;
    std::vector<char> has_cl, is_al;
    std::vector<uint32_t> ref2t, cl_int;
    std::vector<std::string> cl_sorted;
    long long n_al = 0;
    // [G3] pooling
    std::vector<std::vector<uint32_t>> POOLG;
    std::vector<uint32_t> pool_reps;
    size_t n_pooled = 0;
    bool pool_on = false;
    // lineage axis for the EM factors
    std::vector<std::string> lineages;
    std::vector<int> tlin;

    void load_theta_ancestry(const Opts& o);       // logs "Theta: ..."
    void set_refs(std::vector<std::string> r, std::vector<long long> l);
    void build_targets(const Opts& o);              // logs "N reference sequences -> T genes", pooling
    void set_lineages(const std::vector<std::string>& lins);
    uint64_t fingerprint() const;
};

// ------------------------------------------------------------------ parts (evidence of one read chunk)
void write_part(const std::string& path, const Model& M, const Evidence& ev, const std::vector<std::string>& lineages);
// reads refs (for building the model) without the evidence
void read_part_refs(const std::string& path, std::vector<std::string>& refs, std::vector<long long>& lens,
                    std::vector<std::string>& lineages);
void merge_part(const std::string& path, const Model& M, Evidence& ev);

// ------------------------------------------------------------------ EM + output (identical to em_global.py)
void solve_em(const Model& M, Evidence& ev, const Opts& o, double T0);
void write_idxstats(const std::string& path, const Model& M, const Evidence& ev);
void log_records(const Evidence& ev, const std::string& note);
void log_lineage_stats(const LineageStats& st, const std::vector<std::string>& lins, long long n_pv);

}  // namespace lq
