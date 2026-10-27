// lqcore -- multithreaded scan of a name-grouped SAM stream into EM evidence
//   reader thread  -> batches of whole fragments
//   worker threads -> per fragment: candidates, mismatch counts, lineage P, equivalence-class key
//   calling thread -> merges batch results strictly in input order (class numbering = first appearance)
// The result does not depend on the number of workers.
#pragma once
#include "lineage.hpp"
#include "model.hpp"

#include <functional>
#include <memory>

namespace lq {

// em -p: lineage probabilities read from a P_values table (em_global.py semantics)
struct PTable {
    std::vector<std::string> lineages;
    std::unordered_map<std::string, uint32_t> idx;
    std::vector<double> val;
    size_t L = 0;
    void load(const std::string& path, const Opts& o);
};

struct ScanCfg {
    const Model* M = nullptr;
    const Opts* O = nullptr;
    const LineageEngine* LE = nullptr;     // run: lineage evidence computed from the reads
    const PTable* PT = nullptr;            // em -p: lineage evidence from a table
    int workers = 1;
    int chunk_i = 1, chunk_n = 1;          // keep only fragments of blocks b with b % N == i-1 (SAM-level chunking)
    long long block = 100000;
    FILE* dbg_kmer = nullptr;
    FILE* dbg_snp = nullptr;
    FILE* pvals_out = nullptr;
    FILE* dbg_frag = nullptr;              // per fragment: candidates with extra mismatches + lineage P
    const RefGC* gc = nullptr;             // 1.1.0: if set, fragment GC histogram of the primary proper-pair alignment
};

// A batch of whole fragments as SAM text: data = lines ('\n'-terminated), ls = line starts (+ end), fs = first line of
// each fragment (+ sentinel). Producers that start from reads fill `raw` instead and a worker-side `prepare` turns it
// into SAM lines (lqmap); `prepare` may also rewrite SAM lines (sa).
struct ScanBatch {
    size_t id = 0;
    std::string data;
    std::vector<uint32_t> ls, fs;
    std::string raw;                      // producer-defined payload (e.g. read pairs)
    std::vector<uint32_t> roff;
};
using ScanEmit = std::function<void(std::unique_ptr<ScanBatch>)>;
using ScanProducer = std::function<void(const ScanEmit&)>;        // runs on the reader thread, emits batches in order
using ScanPrepare = std::function<void(ScanBatch&, int worker)>;   // runs on a worker before the batch is scanned

// generic ordered scan: the result does not depend on the number of workers
void scan_batches(const ScanProducer& produce, const ScanPrepare& prepare, const ScanCfg& cfg, Evidence& ev);

// reads SAM body lines (header already consumed; `first` is the first body line, may be empty)
void scan_sam(LineReader& rd, const std::string& first, const ScanCfg& cfg, Evidence& ev, const ScanPrepare& prepare = nullptr);

// reads the @SQ header; returns false if the stream has no body line (first is set to the first body line)
bool read_sam_header(LineReader& rd, std::vector<std::string>& refs, std::vector<long long>& lens, std::string& first);

}  // namespace lq
