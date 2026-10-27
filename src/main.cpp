// LineageQuant 1.1.0 -- one executable, C++ only
//
//   lq <command> ...                (bin/lq; bin/lineagequant is the same with the long name)
//   lqcore <command> ...            (bin/lqcore: the engine name of 0.x / 1.0.0 scripts)
//
// The two names differ only in `run`: lq run = whole pipeline, lqcore run = one sample (= quant).
#include "common.hpp"

#include <sched.h>

#ifndef LQ_VERSION
#define LQ_VERSION "1.1.0"
#endif

namespace lq {
int em_main(int argc, char** argv);
int run_main(int argc, char** argv);
int markers_main(int argc, char** argv);
int theta_main(int argc, char** argv);
int map_main(int argc, char** argv);
int lrpick_entry(int argc, char** argv);
int lrcount_entry(int argc, char** argv);
int export_main(int argc, char** argv);
int matrix_main(int argc, char** argv);
int pipeline_main(int argc, char** argv);
int snpdb_msa_main(int argc, char** argv);
int snpdb_vcf_main(int argc, char** argv);
int build_markers_main(int argc, char** argv);
int extract_tx_main(int argc, char** argv);
int normalize_clusters_main(int argc, char** argv);
const char* lq_version() { return LQ_VERSION; }
}

static void usage(const char* prog, bool engine_name) {
    std::fprintf(stderr,
                 "%s %s -- LineageQuant: haplotype-resolved expression quantification for polyploid genomes\n\n"
                 "usage: %s <command> [options]      (<command> -h for its options)\n\n"
                 "pipeline\n"
                 "  %-20s %s\n"
                 "  quant                one sample: alignment -> candidates + lineage evidence -> EM (sample.lqres)\n"
                 "  theta                per-copy prior from cluster-unique k-mers (Step 01)\n"
                 "  markers              lineage-marker index cache (shared by all samples)\n"
                 "  em                   global EM from SAM/BAM or from the --parts of quant --part-out\n"
                 "  map                  lqmap, the built-in aligner: read pairs -> SAM on stdout\n"
                 "  lrpick / lrcount     long reads: one read per molecule / molecules per allele copy (Step 00b)\n\n"
                 "reference preparation (once per reference; run by `lq run` when needed)\n"
                 "  normalize-clusters   allele clusters from other tools -> ClusterID<TAB>g1,g2,... (+ ancestry, orphans)\n"
                 "  extract-transcripts  progenitor genome + GFF -> transcript FASTA (seqid reconciliation, gffread)\n"
                 "  build-markers        lineage-specific k-mer markers (KMC)\n"
                 "  snpdb-msa            SNP database from MAFFT alignments of the copies of each cluster\n"
                 "  snpdb-vcf            SNP database from a population VCF\n\n"
                 "outside the pipeline\n"
                 "  matrix               genome-wide expression matrix (RawReads + TPM per sample) from the .lqres files\n"
                 "  export               .lqres / .lqb binary outputs -> the text files of LineageQuant 1.0.0\n"
                 "  version              print the version\n",
                 prog, LQ_VERSION, prog, "run", engine_name ? "one sample (same as quant; lqcore name)" : "the whole pipeline (references, every sample)");
}

// OMP_PROC_BIND / OMP_PLACES / GOMP_CPU_AFFINITY make libgomp pin the *initial* thread to one CPU when the
// program loads, before main(). Our std::threads and the aligners we launch inherit that mask, so the whole
// pipeline ends up on a single core (POJ root: align+scan 10.8 h instead of ~1 h). Undo it once:
// widen the mask (the kernel intersects it with the job's cpuset), drop the variables, re-exec.
static void undo_omp_binding(char** argv) {
    if (getenv("LQCORE_KEEP_BINDING") || getenv("LQCORE_REEXEC")) return;
    const char* vars[] = {"OMP_PROC_BIND", "OMP_PLACES", "GOMP_CPU_AFFINITY"};
    bool set = false;
    for (const char* v : vars) {
        const char* x = getenv(v);
        if (x && *x && !(std::string(v) == "OMP_PROC_BIND" && (lq::upper(x) == "FALSE"))) set = true;
    }
    if (!set) return;
    cpu_set_t all;
    CPU_ZERO(&all);
    for (int i = 0; i < CPU_SETSIZE; ++i) CPU_SET(i, &all);
    sched_setaffinity(0, sizeof all, &all);
    for (const char* v : vars) unsetenv(v);
    setenv("LQCORE_REEXEC", "1", 1);
    execv("/proc/self/exe", argv);   // only returns on failure; the mask is already widened for this thread
}

int main(int argc, char** argv) {
    undo_omp_binding(argv);
    std::string prog = argv[0];
    size_t sl = prog.find_last_of('/');
    if (sl != std::string::npos) prog = prog.substr(sl + 1);
    const bool engine_name = prog.rfind("lqcore", 0) == 0;
    const char* pname = engine_name ? "lqcore" : prog.rfind("lineagequant", 0) == 0 ? "lineagequant" : "lq";
    if (argc < 2) { usage(pname, engine_name); return 1; }
    std::string cmd = argv[1];
    int n = argc - 2;
    char** v = argv + 2;
    if (cmd == "run") return engine_name ? lq::run_main(n, v) : lq::pipeline_main(n, v);
    if (cmd == "quant") return lq::run_main(n, v);
    if (cmd == "em") return lq::em_main(n, v);
    if (cmd == "markers") return lq::markers_main(n, v);
    if (cmd == "theta") return lq::theta_main(n, v);
    if (cmd == "map") return lq::map_main(n, v);
    if (cmd == "lrpick") return lq::lrpick_entry(n, v);
    if (cmd == "lrcount") return lq::lrcount_entry(n, v);
    if (cmd == "export") return lq::export_main(n, v);
    if (cmd == "matrix") return lq::matrix_main(n, v);
    if (cmd == "snpdb-msa") return lq::snpdb_msa_main(n, v);
    if (cmd == "snpdb-vcf") return lq::snpdb_vcf_main(n, v);
    if (cmd == "build-markers") return lq::build_markers_main(n, v);
    if (cmd == "extract-transcripts") return lq::extract_tx_main(n, v);
    if (cmd == "normalize-clusters") return lq::normalize_clusters_main(n, v);
    if (cmd == "version" || cmd == "--version" || cmd == "-v") {
        std::printf("%s %s\n", pname, LQ_VERSION);
        return 0;
    }
    usage(pname, engine_name);
    return cmd == "-h" || cmd == "--help" ? 0 : 1;
}
