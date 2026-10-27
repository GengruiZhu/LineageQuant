// lqcore em -- global EM engine (Step 03), C++ port of lineagequant/core/em_global.py (1.0.0)
//
// Same command line, same input files, same output columns and the same arithmetic in the same
// order as the Python engine (numpy bincount = sequential accumulation in entry order, numpy sum =
// pairwise summation, Python round(x, 4)). Parallelism never changes a result.
//
// Input (exactly one):
//   --bam BAM        collated on the fly with samtools if not grouped by read name
//   --sam FILE|-     SAM text grouped by read name (e.g. bowtie2 output)
//   --parts P1 P2..  evidence files of read chunks written by `lqcore run --part-out` (summed, then one EM)
#include "scan.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace lq {

static void em_usage() {
    std::fprintf(stderr,
        "usage: lq em -a ANCESTRY (--bam BAM | --sam SAM|- | --parts P1 [P2 ...]) [-t THETA] [-p PVALS]\n"
        "                 [-o OUT.tsv | -o S.lqres (binary)] [--fasta FA] [--no-pool] [--threads N] [--workers N] [--eff-len off|len|eff]\n"
        "                 [--frag-len F] [--min-eff-len L] [--max-iter N] [--min-iter N] [--tol X]\n"
        "                 [--init theta|uniform] [--nm-weight auto|off|LAMBDA]\n"
        "                 [--p-weight-specific W] [--p-weight-conflict W] [--p-weight-none W]\n"
        "                 [--p-alpha A] [--p-floor F] [--theta-prior-lr X] [--theta-prior-kmer X]\n"
        "                 [--theta-prior-cap X] [--idxstats-out FILE]\n"
        "        1.1.0: [--eff-len fld] [--gc-bias (needs --fasta; with --parts the parts need the GC histogram)] [--eff-out FILE]\n");
}

namespace {
struct SamInput {
    FILE* pipe = nullptr;
    std::unique_ptr<LineReader> rd;
    ~SamInput() {
        rd.reset();
        if (pipe) pclose(pipe);
    }
};
std::string sh_quote(const std::string& s) {
    std::string o = "'";
    for (char c : s) o += c == '\'' ? std::string("'\\''") : std::string(1, c);
    return o + "'";
}
bool bam_is_name_grouped(const std::string& bam) {
    std::string cmd = "samtools view -H " + sh_quote(bam);
    FILE* p = popen(cmd.c_str(), "r");
    if (!p) die("cannot run: " + cmd);
    char* line = nullptr;
    size_t cap = 0;
    bool grouped = false;
    while (getline(&line, &cap, p) > 0) {
        std::string_view l(line);
        if (l.rfind("@HD", 0) == 0)
            for (auto f : split_tab(l)) {
                while (!f.empty() && (f.back() == '\n' || f.back() == '\r')) f.remove_suffix(1);
                if (f == "SO:queryname" || f == "GO:query") grouped = true;
            }
    }
    free(line);
    if (pclose(p) != 0) die("samtools view -H failed on " + bam);
    return grouped;
}
void open_input(const Opts& a, SamInput& in) {
    if (!a.sam.empty()) {
        in.rd.reset(new LineReader(a.sam));
        return;
    }
    if (!file_exists(a.bam)) die("BAM not found: " + a.bam);
    std::string th = std::to_string(std::max(1, a.threads));
    std::string cmd;
    if (bam_is_name_grouped(a.bam)) cmd = "samtools view -h -@ " + th + " " + sh_quote(a.bam);
    else {
        std::string dir = a.out.find('/') == std::string::npos ? "." : a.out.substr(0, a.out.rfind('/'));
        const char* job = getenv("PBS_JOBID");
        std::string prefix = dir + "/.lqcore_collate_" + (job ? std::string(job) : std::to_string(getpid()));
        logm("[INFO] collating by read name (samtools collate, streamed)");
        cmd = "samtools collate -O -u -@ " + th + " " + sh_quote(a.bam) + " " + sh_quote(prefix) +
              " | samtools view -h -@ " + th + " -";
    }
    in.pipe = popen(cmd.c_str(), "r");
    if (!in.pipe) die("cannot run: " + cmd);
    in.rd.reset(new LineReader(dup(fileno(in.pipe))));
}
}  // namespace

int em_main(int argc, char** argv) {
    Opts a = parse_opts(argc, argv, em_usage);
    const double T0 = wall();
    if (a.ancestry.empty()) { em_usage(); die("-a/--ancestry is required"); }
    if ((int)!a.bam.empty() + (int)!a.sam.empty() + (int)!a.parts.empty() != 1) {
        em_usage();
        die("give exactly one of --bam / --sam / --parts");
    }
#ifdef _OPENMP
    omp_set_num_threads(a.threads);
#endif
    if (!a.unknown.empty()) {
        std::string s;
        for (auto& u : a.unknown) s += (s.empty() ? "" : " ") + u;
        logm("[INFO] ignoring legacy-engine options: " + s);
    }
    const bool POOL_ON = !a.fasta.empty() && !a.no_pool;
    logm("[INFO] EM (global) | eff-len=" + a.eff_len + (a.gc_bias ? " + GC bias" : "") + " | identical-sequence pooling=" + (POOL_ON ? "on" : "off") +
         " | nm-weight=" + a.nm_weight + " | init=" + a.init + " | tol=" + pyfloat(a.tol) +
         " max_iter=" + std::to_string(a.max_iter) + " | engine=lqcore threads=" + std::to_string(a.threads));

    Model M;
    M.load_theta_ancestry(a);
    Evidence ev;
    if (!a.parts.empty()) {
        std::vector<std::string> refs, lins;
        std::vector<long long> lens;
        read_part_refs(a.parts[0], refs, lens, lins);
        M.set_refs(refs, lens);
        M.build_targets(a);
        M.set_lineages(lins);
        ev.init(M.T, M.refs.size());
        double t0 = wall();
        for (auto& p : a.parts) merge_part(p, M, ev);
        logm("[INFO] merged " + std::to_string(a.parts.size()) + " part(s)" + sfmt(" (%.0fs)", wall() - t0));
        if (!lins.empty()) log_lineage_stats(ev.ls, lins, ev.n_pv_reads);
        log_records(ev, "");
    } else {
        SamInput in;
        open_input(a, in);
        std::vector<std::string> refs;
        std::vector<long long> lens;
        std::string first;
        read_sam_header(*in.rd, refs, lens, first);
        M.set_refs(refs, lens);
        M.build_targets(a);
        PTable pt;
        bool has_pt = !a.pvals.empty() && file_exists(a.pvals);
        if (has_pt) {
            pt.load(a.pvals, a);
            M.set_lineages(pt.lineages);
        }
        ev.init(M.T, M.refs.size());
        ScanCfg cfg;
        cfg.M = &M;
        cfg.O = &a;
        cfg.PT = has_pt ? &pt : nullptr;
        cfg.workers = a.workers > 0 ? a.workers : a.threads;
        RefGC gcmap;
        if (a.gc_bias) {
            if (a.fasta.empty()) die("--gc-bias needs --fasta");
            load_ref_gc(gcmap, a.fasta, M.refs, M.r2i);
            cfg.gc = &gcmap;
        }
        double ts = wall();
        scan_sam(*in.rd, first, cfg, ev);
        in.rd.reset();
        if (in.pipe) {
            int st = pclose(in.pipe);
            in.pipe = nullptr;
            if (st != 0) die("input pipeline (samtools) failed");
        }
        log_records(ev, sfmt("  [scan %.0fs]", wall() - ts));
    }
    if (!a.idxstats_out.empty()) write_idxstats(a.idxstats_out, M, ev);
    solve_em(M, ev, a, T0);
    return 0;
}

}  // namespace lq
