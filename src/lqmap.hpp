// lqcore -- lqmap: the LineageQuant aligner (built in; replaces bowtie2 -k 200 with LineageQuant's settings)
//
// What LineageQuant needs from an aligner is narrower than what bowtie2 -k 200 does:
//   - ungapped, end-to-end alignments (bowtie2 runs with --gbar 99999999, i.e. no gaps);
//   - per mate at most floor(0.1 x read length) mismatches, N not penalised (--score-min L,0,-0.1 --mp 1,1 --np 0);
//   - concordant pairs only, fragment length 1..1000 (--no-mixed --no-discordant -I 1 -X 1000, --fr);
//   - the mismatch count (NM) of every candidate, because the EM weights candidates by it (G4); a candidate 3
//     mismatches worse than the best already carries a weight below ~1e-8, so candidates beyond best + 3 are dropped.
// lqmap finds candidates with exact k-mer seeds taken every --map-stride bases on both strands (index of all k-mer
// positions of the transcriptome), verifies each candidate position by counting mismatches, pairs the mates on the
// same transcript, and rescues a mate by an ungapped scan of the fragment window when seeds missed it.
// Output: the same SAM records LineageQuant reads from bowtie2 (one pair of lines per candidate transcript).
#pragma once
#include "scan.hpp"

#include <array>
#include <atomic>
#include <memory>

namespace lq {

struct MapIndex {
    int k = 21, lowbits = 0, b = 0;
    std::vector<std::string> names;
    std::vector<long long> lens;
    std::vector<uint64_t> tstart;          // offset of each transcript in seq (+ end sentinel)
    std::vector<uint8_t> seq;              // 0-3 ACGT, 4 anything else (upper-cased)
    std::vector<uint32_t> off;             // bucket -> [off[b], off[b+1]) in low / pos
    std::vector<uint32_t> low, pos;        // per k-mer occurrence: low bits of the code (sorted within bucket), position
    void build(const std::string& fasta, int k, int threads);
    bool load(const std::string& path, const std::string& sig);
    void save(const std::string& path, const std::string& sig) const;
    // occurrences of one k-mer code: [*b, *e) into pos
    inline void find(uint64_t code, uint32_t& lo, uint32_t& hi) const {
        const uint64_t bk = code >> lowbits;
        const uint32_t lv = (uint32_t)(code & ((1ULL << lowbits) - 1));
        const uint32_t* L = low.data();
        uint32_t a = off[bk], z = off[bk + 1];
        const uint32_t* p = std::lower_bound(L + a, L + z, lv);
        const uint32_t* q = std::upper_bound(p, L + z, lv);
        lo = (uint32_t)(p - L);
        hi = (uint32_t)(q - L);
    }
    inline uint32_t tid_of(uint64_t g) const {
        return (uint32_t)(std::upper_bound(tstart.begin(), tstart.end(), g) - tstart.begin() - 1);
    }
};

struct MapOpts {
    int stride = 8, max_extra = 3, max_occ = 4000, max_cand = 1000, max_frag = 1000, min_frag = 1;
    double mm_frac = 0.1, n_ceil_frac = 0.5;
    bool rescue = true;
};

struct MapStats {
    std::atomic<long long> pairs{0}, conc0{0}, conc1{0}, concm{0}, rescued{0}, cands{0}, reseeded{0};
};

// one per worker thread (scratch buffers)
class Mapper {
public:
    Mapper(const MapIndex& X, const MapOpts& o, MapStats& st) : X_(X), o_(o), st_(st) {}
    // appends the SAM lines of one read pair (name without /1 /2) to `out`, records line starts in `ls`
    void map_pair(std::string_view name, std::string_view r1, std::string_view r2, std::string& out, std::vector<uint32_t>& ls);
    // same, from candidate positions found by another aligner: given[m] = (transcript, 0-based start, reverse) per mate
    void map_hits(std::string_view name, std::string_view r1, std::string_view r2,
                  const std::vector<std::array<int64_t, 3>>* given, std::string& out, std::vector<uint32_t>& ls);

private:
    struct Hit { uint32_t t; int32_t st; uint8_t rev; int16_t mm, nm; };
    struct PC { uint32_t t, i0, i1; int nm, frag; };
    std::vector<PC> pc_;
    const MapIndex& X_;
    const MapOpts& o_;
    MapStats& st_;
    std::vector<uint8_t> q_[2][2];          // [mate][strand 0 = read, 1 = reverse complement] codes
    std::vector<uint64_t> diag_;
    std::vector<Hit> h_[2];
    void encode(std::string_view s, int m);
    void seeds(int m, int L, int maxmm, bool mm1 = false);
    bool verify(const uint8_t* q, int L, uint64_t g, int maxmm, int& mm, int& nm) const;
    void rescue(int have, const Hit& a, int Lm, int maxmm, std::vector<Hit>& add);
    void finish(std::string_view name, const int* L, const int* maxmm, std::string& out, std::vector<uint32_t>& ls, bool resc);
};

// reads R1/R2 with two decompression threads and emits batches of read pairs (raw = name\0r1\0r2\0, roff = starts)
long long produce_pairs(const std::string& r1, const std::string& r2, int chunk_i, int chunk_n, long long block,
                        size_t per_batch, const ScanEmit& emit);
// aligns every pair of a batch, replacing raw with SAM lines (data / ls / fs)
void map_batch(Mapper& mp, ScanBatch& b);
void log_map_stats(const MapStats& s, const char* who);
MapOpts map_opts_from(const Opts& a);
// builds, or loads from --map-cache, the index (k = 0: sequences only)
void load_map_index(MapIndex& X, const Opts& a, int k, int threads);
// sa: rewrites a batch of Salmon --writeMappings records (one fragment per group) into verified, paired SAM lines
struct SaStats { std::atomic<long long> frags{0}, recs{0}, no_seq{0}, unknown_ref{0}; };
void sa_batch(Mapper& mp, const std::unordered_map<std::string, uint32_t>& name2t, ScanBatch& b, SaStats& st);

}  // namespace lq
