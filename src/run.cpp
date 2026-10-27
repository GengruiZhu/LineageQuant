// lqcore run -- one pass per sample: bowtie2 -> (no sort, no BAM) -> candidates + lineage evidence -> EM
//
//   whole sample on one node:
//     lqcore run --r1 R1 --r2 R2 -x INDEX -a ANC -t THETA --fasta FA [lineage options] -o EM_Final_Expression.tsv
//   split by reads over several jobs/nodes (evidence is additive):
//     lqcore run ... --chunk 2/4 --part-out sample.part2     (x4, any nodes)
//     lqcore em  ... --parts sample.part1 ... sample.part4 -o EM_Final_Expression.tsv
//   Read pairs are dealt to chunks in blocks of --block-pairs (default 100,000), round robin, through two FIFOs.
//
// bowtie2 settings default to LineageQuant 1.0.0 (RSEM's, N not penalised); LQ_BT2_OPTS overrides them like
// core/mapping.sh. --sam-in FILE|- replaces the aligner (testing, or any name-grouped SAM).
//
//   lqcore markers --markers N:F ... --kmer-len 31 --marker-cache FILE   builds the marker index cache once
#include "lqmap.hpp"
#include "scan.hpp"

#include <atomic>
#include <csignal>
#include <fcntl.h>
#include <sched.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <thread>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace lq {

static const char* BT2_DEFAULT =
    "--sensitive --dpad 0 --gbar 99999999 --mp 1,1 --np 0 --n-ceil L,0,0.5 --score-min L,0,-0.1 "
    "-I 1 -X 1000 --no-mixed --no-discordant";

static void run_usage() {
    std::fprintf(stderr,
        "usage: lq quant (= lqcore run) -a ANCESTRY (--r1 R1 --r2 R2 [--aligner bowtie2|lqmap|sa] | --sam-in SAM|-) [-t THETA] [--fasta FA]\n"
        "        [-o EM_OUT | --part-out PART] [--chunk i/N] [--block-pairs N] [--threads N] [--workers N] [--tmpdir DIR]\n"
        "        bowtie2 (default): -x BT2_INDEX [--bt2-threads N] [-k 200] [--bt2-opts \"...\"] [--bt2-log FILE] [--bowtie2 PATH]\n"
        "        lqmap (built in):  --fasta FA [--map-cache FILE] [--map-k 21] [--map-stride 8] [--map-max-extra 3]\n"
        "                           [--map-max-occ 4000] [--map-max-cand 1000] [--map-max-frag 1000] [--map-mm-frac 0.1] [--map-no-rescue]\n"
        "        sa (salmon):       --fasta FA --salmon-index DIR [--salmon PATH] [--salmon-opts \"...\"]\n"
        "        EM acceleration:   [--em-accel none|squarem] [--no-pin]\n"
        "        lineage: [--target-ids FILE] [--markers NAME:FASTA ...] [--kmer-len 31] [--marker-cache FILE]\n"
        "                 [--kmer-alpha 1] [--min-hits 3] [--specific-thresh 0.9] [--snp-db EXPORTED_TSV]\n"
        "                 [--error-rate 0.01] [--no-chimeric] [--snp-pick best|most|first]\n"
        "                 [--poe-min-weight 0.05] [--poe-eps 1e-6]\n"
        "        EM: same options as `lqcore em` (--nm-weight, --eff-len, --p-weight-*, --theta-prior-*, ...)\n"
        "        1.1.0: [--eff-len fld] [--gc-bias] [--gc-hist] [--eff-out FILE] [--marker-clean (with --fasta and -a)]\n"
        "        outputs: -o S.lqres (binary: EM + effective lengths + per-reference counts; lq export)\n"
        "                 or -o EM_Final_Expression.tsv [--idxstats-out FILE] [--eff-out FILE] (text, as 1.0.0)\n"
        "                 [--write-pvals FILE] [--debug-kmer FILE] [--debug-snp FILE] [--debug-frag FILE]\n"
        "        -t THETA: Global_Theta_Matrix.lqb (theta --binary) or .tsv\n");
}

static std::string sh_quote(const std::string& s) {
    std::string o = "'";
    for (char c : s) o += c == '\'' ? std::string("'\\''") : std::string(1, c);
    return o + "'";
}

// writes the selected read records of one FASTQ into a FIFO
static void feed_fifo(const std::string& src, const std::string& fifo, long long block, int ci, int cn,
                      std::atomic<long long>* kept, std::atomic<int>* err) {
    int fd = open(fifo.c_str(), O_WRONLY);
    if (fd < 0) { err->store(1); return; }
    FILE* f = fdopen(fd, "w");
    std::vector<char> buf(1 << 22);
    setvbuf(f, buf.data(), _IOFBF, buf.size());
    LineReader rd(src);
    std::string_view l;
    long long rec = 0, n = 0;
    int li = 0;
    bool keep = false, broken = false;
    while (rd.next(l)) {
        if (li == 0) keep = ((rec / block) % cn) == (ci - 1);
        if (keep && !broken) {
            if (fwrite(l.data(), 1, l.size(), f) != l.size() || fputc('\n', f) == EOF) broken = true;
        }
        if (++li == 4) { li = 0; if (keep) ++n; ++rec; }
        if (broken) break;
    }
    if (li != 0 && !broken) err->store(2);   // truncated FASTQ record
    if (fclose(f) != 0) broken = true;
    if (broken) err->store(3);
    kept->store(n);
}

static void print_tail(const std::string& path, int n) {
    FILE* f = fopen(path.c_str(), "r");
    if (!f) return;
    std::vector<std::string> lines;
    char* line = nullptr;
    size_t cap = 0;
    while (getline(&line, &cap, f) > 0) {
        lines.emplace_back(line);
        if ((int)lines.size() > n) lines.erase(lines.begin());
    }
    free(line);
    fclose(f);
    for (auto& s : lines) { std::printf("  | %s", s.c_str()); }
    std::fflush(stdout);
}

int run_main(int argc, char** argv) {
    Opts a = parse_opts(argc, argv, run_usage);
    const double T0 = wall();
    if (a.ancestry.empty()) { run_usage(); die("-a/--ancestry is required"); }
    const bool from_sam = !a.sam.empty();
    const std::string aln = from_sam ? "sam" : a.aligner;
    if (!from_sam && (a.r1.empty() || a.r2.empty())) { run_usage(); die("give --r1 --r2 (or --sam-in)"); }
    if (aln == "bowtie2" && a.index.empty()) { run_usage(); die("--aligner bowtie2 needs -x BT2_INDEX"); }
    if (aln == "sa" && a.salmon_index.empty()) die("--aligner sa needs --salmon-index DIR (salmon index -k 31 --keepDuplicates)");
    if ((aln == "lqmap" || aln == "sa") && a.fasta.empty()) die("--aligner " + aln + " needs --fasta");
    signal(SIGPIPE, SIG_IGN);
    const int threads = a.threads;
    // threads: bowtie2 / salmon get what the scan does not use; lqmap aligns inside the scan workers
    int workers, ext_threads = 0;
    if (aln == "sam") workers = a.workers > 0 ? a.workers : threads;
    else if (aln == "lqmap") workers = a.workers > 0 ? a.workers : std::max(1, threads - 2);
    else if (aln == "sa") workers = a.workers > 0 ? a.workers : std::max(2, threads / 4);
    else workers = a.workers > 0 ? a.workers : std::max(2, threads / 8);
    if (aln == "bowtie2") ext_threads = a.bt2_threads > 0 ? a.bt2_threads : std::max(1, threads - workers);
    if (aln == "sa") ext_threads = std::max(1, threads - workers);
#ifdef _OPENMP
    omp_set_num_threads(threads);
#endif
    {   // CPUs this process may actually use (job cpuset, taskset, ...); a pinned mask throttles everything
        cpu_set_t cs;
        CPU_ZERO(&cs);
        int ncpu = sched_getaffinity(0, sizeof cs, &cs) == 0 ? CPU_COUNT(&cs) : -1;
        int want = workers + ext_threads + (aln == "lqmap" ? 2 : 0);
        logm("[INFO] CPUs available to this process: " + std::to_string(ncpu) + " (threads requested: " + std::to_string(want) + ")");
        if (ncpu > 0 && ncpu < want)
            logm("[WARN] fewer CPUs than threads: the aligner and the scan will share " + std::to_string(ncpu) +
                 " CPU(s). Check the job's ppn / taskset / OMP_PLACES.");
    }
    std::string outbase = !a.part_out.empty() ? a.part_out : a.out;
    // FIFOs belong on a node-local file system
    std::string tmpdir = !a.tmpdir.empty() ? a.tmpdir : (getenv("TMPDIR") && *getenv("TMPDIR") ? getenv("TMPDIR") : "/tmp");
    std::string what = aln == "sam" ? "input SAM " + a.sam
                       : aln == "bowtie2" ? "bowtie2 -k " + std::to_string(a.max_aln)
                       : aln == "lqmap" ? sfmt("lqmap (k=%d, stride %d, best + %d mismatches)", a.map_k, a.map_stride, a.map_max_extra)
                                        : "salmon selective alignment -> ungapped verification";
    logm("[INFO] lqcore run | " + what + " | chunk " + std::to_string(a.chunk_i) + "/" + std::to_string(a.chunk_n) +
         (a.chunk_n > 1 ? " (blocks of " + commas(a.block_pairs) + (from_sam ? " fragments)" : " read pairs)") : "") +
         " | threads: " + (ext_threads ? aln + " " + std::to_string(ext_threads) + " + " : "") +
         (aln == "lqmap" ? "align+scan " : "scan ") + std::to_string(workers) + " | EM " + std::to_string(threads) +
         (a.em_accel != "none" ? " (" + a.em_accel + ")" : ""));

    Model M;
    M.load_theta_ancestry(a);

    // ---------------------------------------------------------------- input stream / aligner
    FILE* pipe = nullptr;
    std::unique_ptr<LineReader> rd;
    std::vector<std::string> fifos;
    std::vector<std::thread> feeders;
    std::atomic<long long> kept1{0}, kept2{0};
    std::atomic<int> ferr{0};
    std::string ext_log, sa_fifo;
    std::thread sa_watch;
    std::atomic<int> sa_status{-1};
    MapIndex X;
    std::unordered_map<std::string, uint32_t> name2t;
    std::vector<std::string> refs;
    std::vector<long long> lens;
    std::string first;
    auto release_feeders = [&]() {
        // if the aligner died before opening a FIFO, the writer blocks in open(): open + close the read end
        for (auto& f : fifos) { int fd = open(f.c_str(), O_RDONLY | O_NONBLOCK); if (fd >= 0) close(fd); }
    };
    auto make_fifo = [&](const std::string& tag) {
        const char* job = getenv("PBS_JOBID");
        std::string f = tmpdir + "/.lqcore_" + (job ? std::string(job) : std::string("nojob")) + "_" + std::to_string(getpid()) + "_" + tag + ".fifo";
        unlink(f.c_str());
        if (mkfifo(f.c_str(), 0600) != 0) die("cannot create FIFO " + f);
        return f;
    };
    if (aln == "sam") rd.reset(new LineReader(a.sam));
    else if (aln == "lqmap") {
        for (auto& f : {a.r1, a.r2}) if (!file_exists(f)) die("FASTQ not found: " + f);
        load_map_index(X, a, a.map_k, threads);
        refs = X.names;
        lens = X.lens;
    } else {   // bowtie2 or sa: external aligner writing SAM
        for (auto& f : {a.r1, a.r2}) if (!file_exists(f)) die("FASTQ not found: " + f);
        std::string in1 = a.r1, in2 = a.r2;
        if (a.chunk_n > 1) {
            fifos.push_back(make_fifo("R1"));
            fifos.push_back(make_fifo("R2"));
            in1 = fifos[0];
            in2 = fifos[1];
        }
        std::string cmd;
        if (aln == "bowtie2") {
            std::string opts = a.bt2_opts;
            if (opts.empty()) { const char* e = getenv("LQ_BT2_OPTS"); opts = e ? e : BT2_DEFAULT; }
            ext_log = !a.bt2_log.empty() ? a.bt2_log : outbase + ".bowtie2.log";
            cmd = sh_quote(a.bowtie2) + " -p " + std::to_string(ext_threads) + " -k " + std::to_string(a.max_aln) + " " + opts +
                  " -x " + sh_quote(a.index) + " -1 " + sh_quote(in1) + " -2 " + sh_quote(in2) + " 2> " + sh_quote(ext_log);
        } else {
            load_map_index(X, a, 0, threads);           // transcript sequences for the verification
            for (size_t i = 0; i < X.names.size(); ++i) name2t.emplace(X.names[i], (uint32_t)i);
            refs = X.names;
            lens = X.lens;
            sa_fifo = make_fifo("sa");
            ext_log = !a.bt2_log.empty() ? a.bt2_log : outbase + ".salmon.log";
            cmd = sh_quote(a.salmon) + " quant -i " + sh_quote(a.salmon_index) + " -l A -1 " + sh_quote(in1) + " -2 " + sh_quote(in2) +
                  " -p " + std::to_string(ext_threads) + " --writeMappings=" + sh_quote(sa_fifo) +
                  " --writeQualities -o " + sh_quote(outbase + ".salmon") +
                  (a.salmon_opts.empty() ? "" : " " + a.salmon_opts) + " > " + sh_quote(ext_log) + " 2>&1";
        }
        logm("[EXEC] " + cmd);
        pipe = popen(cmd.c_str(), "r");
        if (!pipe) die("cannot start " + aln);
        if (a.chunk_n > 1) {
            feeders.emplace_back(feed_fifo, a.r1, fifos[0], a.block_pairs, a.chunk_i, a.chunk_n, &kept1, &ferr);
            feeders.emplace_back(feed_fifo, a.r2, fifos[1], a.block_pairs, a.chunk_i, a.chunk_n, &kept2, &ferr);
        }
        if (aln == "bowtie2") rd.reset(new LineReader(dup(fileno(pipe))));
        else {
            // salmon writes the mappings into the FIFO; if it dies before opening it, the watcher unblocks our open()
            FILE* pp = pipe;
            pipe = nullptr;
            std::string ff = sa_fifo;
            sa_watch = std::thread([pp, ff, &sa_status]() {
                int st = pclose(pp);
                sa_status = WIFEXITED(st) ? WEXITSTATUS(st) : 255;
                int fd = open(ff.c_str(), O_WRONLY | O_NONBLOCK);
                if (fd >= 0) close(fd);
            });
            rd.reset(new LineReader(sa_fifo));
        }
    }
    auto stop_external = [&]() {
        int st = 0;
        if (pipe) { st = pclose(pipe); pipe = nullptr; }
        if (sa_watch.joinable()) { sa_watch.join(); st = sa_status.load() == 0 ? 0 : (sa_status.load() << 8); }
        release_feeders();
        for (auto& t : feeders) t.join();
        feeders.clear();
        for (auto& f : fifos) unlink(f.c_str());
        if (!sa_fifo.empty()) unlink(sa_fifo.c_str());
        return st;
    };
    if (aln != "lqmap") {
        std::vector<std::string> hr;
        std::vector<long long> hl;
        read_sam_header(*rd, hr, hl, first);
        if (hr.empty()) {
            rd.reset();
            stop_external();
            if (!ext_log.empty()) print_tail(ext_log, 15);
            die("no @SQ header from the aligner (see " + (ext_log.empty() ? a.sam : ext_log) + ")");
        }
        if (aln == "sa") {
            size_t miss = 0;
            for (auto& n : hr) if (!name2t.count(n)) ++miss;
            logm(sfmt("[INFO] salmon index: %s transcripts in the SAM header, %s not in --fasta", commas((long long)hr.size()).c_str(),
                      commas((long long)miss).c_str()));
        } else { refs = hr; lens = hl; }
    }
    M.set_refs(refs, lens);
    M.build_targets(a);

    LineageEngine LE;
    LE.init(a, M);
    const bool lin = LE.kmer_on || LE.snp_on;
    if (lin) M.set_lineages(LE.em_lins);

    // ---------------------------------------------------------------- optional per-fragment tables
    auto open_out = [&](const std::string& p, const std::string& header) -> FILE* {
        if (p.empty()) return nullptr;
        FILE* f = fopen(p.c_str(), "w");
        if (!f) die("cannot write " + p);
        fputs(header.c_str(), f);
        return f;
    };
    std::string hk = "Read_Name", hs = "Read_Name\tTarget_Gene", hp = "Read_Name\tTarget_Gene";
    for (auto& n : LE.kmer_names) hk += "\tP(" + n + ")";
    hk += "\tLineage_Status\tN_Hits\n";
    for (auto& n : LE.snp_lins) hs += "\tP(" + n + ")";
    hs += "\tLineage_Status\tN_SNPs\n";
    for (auto& n : LE.em_lins) hp += "\tP(" + n + ")";
    hp += "\tLineage_Status\n";
    FILE* fk = open_out(a.debug_kmer, hk);
    FILE* fsn = open_out(a.debug_snp, hs);
    FILE* fpv = lin ? open_out(a.write_pvals, hp) : nullptr;
    std::string hf = "Read_Name\tN_Targets\tTargets(extra_mismatches)";
    for (auto& n : LE.em_lins) hf += "\tP(" + n + ")";
    FILE* ffr = open_out(a.debug_frag, hf + "\n");

    // ---------------------------------------------------------------- scan
    Evidence ev;
    ev.init(M.T, M.refs.size());
    ScanCfg cfg;
    cfg.M = &M;
    cfg.O = &a;
    cfg.LE = lin ? &LE : nullptr;
    cfg.workers = workers;
    if (from_sam && a.chunk_n > 1) { cfg.chunk_i = a.chunk_i; cfg.chunk_n = a.chunk_n; cfg.block = a.block_pairs; }
    cfg.dbg_kmer = fk;
    cfg.dbg_snp = fsn;
    cfg.pvals_out = fpv;
    cfg.dbg_frag = ffr;
    RefGC gcmap;
    if (a.gc_bias || a.gc_hist) {   // 1.1.0: fragment GC histogram for the GC-bias model
        if (a.fasta.empty()) die("--gc-bias / --gc-hist need --fasta (transcript sequences)");
        load_ref_gc(gcmap, a.fasta, M.refs, M.r2i);
        cfg.gc = &gcmap;
    }
    double ts = wall();
    MapOpts mo = map_opts_from(a);
    MapStats mst;
    SaStats sst;
    std::vector<std::unique_ptr<Mapper>> mps;
    if (aln == "lqmap" || aln == "sa") for (int i = 0; i < workers; ++i) mps.emplace_back(new Mapper(X, mo, mst));
    long long npairs = -1;
    if (aln == "lqmap") {
        scan_batches([&](const ScanEmit& emit) {
                         npairs = produce_pairs(a.r1, a.r2, a.chunk_i, a.chunk_n, a.block_pairs, 4096, emit);
                     },
                     [&](ScanBatch& b, int w) { map_batch(*mps[w], b); }, cfg, ev);
    } else if (aln == "sa") {
        scan_sam(*rd, first, cfg, ev, [&](ScanBatch& b, int w) { sa_batch(*mps[w], name2t, b, sst); });
    } else {
        scan_sam(*rd, first, cfg, ev);
    }
    rd.reset();
    for (FILE* f : {fk, fsn, fpv, ffr}) if (f && fclose(f) != 0) die("error writing a per-fragment table");
    if (aln == "bowtie2" || aln == "sa") {
        int st = stop_external();
        print_tail(ext_log, 6);
        if (st != 0) die(aln + " failed (exit status " + std::to_string(WIFEXITED(st) ? WEXITSTATUS(st) : (st >> 8)) + "), see " + ext_log);
        if (ferr.load()) die("feeding reads to " + aln + " failed (code " + std::to_string(ferr.load()) + ")");
        if (a.chunk_n > 1) {
            if (kept1.load() != kept2.load()) die("R1 and R2 chunks differ in size: " + commas(kept1.load()) + " vs " + commas(kept2.load()));
            logm("[INFO] chunk " + std::to_string(a.chunk_i) + "/" + std::to_string(a.chunk_n) + ": " + commas(kept1.load()) + " read pairs");
        }
    }
    if (aln == "lqmap") {
        log_map_stats(mst, "lqmap");
        if (a.chunk_n > 1) logm("[INFO] chunk " + std::to_string(a.chunk_i) + "/" + std::to_string(a.chunk_n) + ": " + commas(npairs) + " read pairs");
    }
    if (aln == "sa") {
        log_map_stats(mst, "sa (salmon candidates, ungapped verification)");
        logm(sfmt("[INFO] sa: %s salmon records in %s fragments; %s fragments without a read sequence, %s records on unknown references",
                  commas(sst.recs.load()).c_str(), commas(sst.frags.load()).c_str(), commas(sst.no_seq.load()).c_str(),
                  commas(sst.unknown_ref.load()).c_str()));
    }
    log_records(ev, sfmt("  [align+scan %.0fs]", wall() - ts));
    if (lin) log_lineage_stats(ev.ls, LE.em_lins, ev.n_pv_reads);

    if (!a.part_out.empty()) {
        write_part(a.part_out, M, ev, M.lineages);
        logm("[SUCCESS] evidence of chunk " + std::to_string(a.chunk_i) + "/" + std::to_string(a.chunk_n) + ": " + a.part_out +
             sfmt(" (%.0fs)", wall() - T0));
        return 0;
    }
    if (!a.idxstats_out.empty()) write_idxstats(a.idxstats_out, M, ev);
    solve_em(M, ev, a, T0);
    return 0;
}

int markers_main(int argc, char** argv) {
    Opts a = parse_opts(argc, argv, run_usage);
    if (a.markers.empty() || a.marker_cache.empty())
        die("usage: lqcore markers --markers NAME:FASTA ... --kmer-len K --marker-cache FILE [--marker-clean --fasta REF.fa -a ANCESTRY]");
    Model M;
    if (a.marker_clean) {
        if (a.fasta.empty() || a.ancestry.empty()) die("--marker-clean needs --fasta (reference transcripts) and -a (their lineages)");
        M.load_theta_ancestry(a);
    }
    LineageEngine LE;
    a.snp_db.clear();
    a.target_ids.clear();
    LE.init(a, M);
    logm("[SUCCESS] marker cache: " + a.marker_cache);
    return 0;
}

}  // namespace lq
