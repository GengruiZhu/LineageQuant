// lqcore -- per-fragment lineage probabilities (ports of snipe_kmer.py v2, snipe_snp.py, merge_poe.py)
#pragma once
#include "model.hpp"
#include "sam.hpp"

namespace lq {

// ------------------------------------------------------------------ lineage-specific k-mers
// All marker k-mers of all lineages in one array sorted by fmix64(k-mer) (a bijection), bucketed by the
// top bits: a lookup is one offset read plus a short scan of one or two cache lines. ~9 bytes per k-mer.
struct MarkerIndex {
    int k = 0, b = 0;
    bool canonical = false;
    std::vector<std::string> names;          // as given on the command line
    std::vector<uint64_t> count;             // distinct k-mers per lineage
    std::vector<uint32_t> off;
    std::vector<uint64_t> h;
    std::vector<uint8_t> mask;
    void build(const std::vector<std::pair<std::string, std::string>>& name_path, int k, int threads);
    bool load_cache(const std::string& path, const std::string& signature);
    void save_cache(const std::string& path, const std::string& signature) const;
    inline uint8_t get(uint64_t km) const;
    // 1.1.0: drop a k-mer from lineage M's markers if it occurs in a reference sequence labelled with another
    // lineage (k-mers absent from the reference are kept: they are the ones that can point to unassembled copies)
    void clean(const std::string& fasta, const std::unordered_map<std::string, std::string>& anc, int threads);
};

// ------------------------------------------------------------------ lineage-informative SNPs (exported SNP DB)
struct SnpIndex {
    bool msa = true;
    std::vector<std::string> lineages;       // order of the DB (sorted, as snipe_snp._infer_lineages)
    std::vector<std::string> tids;           // DB transcripts, first-appearance order
    std::vector<uint64_t> beg;               // per tid: [beg[i], beg[i+1]) in the arrays below
    std::vector<int32_t> pos;                // 0-based transcript position
    std::vector<uint8_t> exp;                // msa: expected base code per lineage (0-3, 4 = never matches)
    std::vector<uint8_t> ref, alt;           // population
    std::vector<double> af;                  // population: allele frequency per lineage
    std::unordered_map<std::string, uint32_t> resolve;   // name -> tid index (exact, or "base" -> first variant)
    void load(const std::string& path);
    int lookup(const std::string& bam_ref) const;
};

struct LinScratch {
    std::vector<std::pair<uint32_t, uint32_t>> rr;      // (rid, record index)
    std::vector<double> ll;
};

class LineageEngine {
public:
    bool kmer_on = false, snp_on = false;
    std::vector<std::string> kmer_names;                 // as given (debug output)
    std::vector<std::string> kmer_lins, snp_lins, em_lins;   // upper case
    void init(const Opts& o, const Model& M);
    // Evaluates one fragment. Returns true if the EM receives a lineage vector (pv, em_lins order).
    bool eval(const Rec* recs, size_t n, std::string_view qname, LinScratch& s, LineageStats& st, double* pv,
              std::string* dbg_kmer, std::string* dbg_snp, std::string* prow) const;

private:
    const Model* M_ = nullptr;
    Opts o_;
    MarkerIndex mk_;
    SnpIndex snp_;
    std::vector<char> valid_ref_;
    std::vector<int32_t> snp_tid_;                      // per reference: SnpIndex tid or -1
    std::vector<std::string> ref_gene_dot1_, ref_gene_rdot_;   // tid.split(".")[0] / tid.rsplit(".",1)[0]
    std::vector<int> map_k_, map_s_;                     // em_lins -> index in kmer / snp vectors (-1: absent)
    double LM_ = 0, LX_ = 0;
    void kmer_eval(const Rec* recs, size_t n, double* p, int& status, long long& hits) const;
    void snp_eval(const Rec* recs, const std::vector<uint32_t>& idx, LinScratch& s, double* p, int& status, int& best,
                  long long& nsnp) const;
};

}  // namespace lq
