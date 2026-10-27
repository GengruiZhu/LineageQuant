// lq run -- the whole pipeline (C++ port of cli.py / config.py / router.py of LineageQuant 1.0.0)
//
// The user declares the ploidy and gives whatever lineage inputs exist; the route follows from them. Every step
// is checkpointed with a .done file and skipped when it exists, so a run resumes where it stopped.
//
//   Step 00b  long-read priors                        lrcount                        (--longread)
//   Step 1a   SNP database (allo; VCF or MSA route)    snpdb-vcf / snpdb-msa          (--vcf / --msa-snp, --phased-fasta)
//   markers   lineage k-mer markers (allo)            extract-transcripts + build-markers (--ancestor) or --markers
//             marker index cache                      markers
//   per sample
//   Step 01   Theta (per-copy prior)                  theta --binary            -> <sample>/Global_Theta_Matrix.lqb
//   Step 00-03 alignment + lineage evidence + EM      quant                     -> <sample>/<sample>.lqres
//
// The genome-wide expression matrix (Step 04 of 1.0.0) is not part of the run: `lq matrix`.
// With the default options every step gives the numbers of LineageQuant 1.0.0 (tests/compare_*.sh), except
// that the MSA SNP database lists every lineage (--legacy-lineages restores 1.0.0; see snpdb-msa -h).
//
// Design notes kept from 1.0.0 (each changed results):
//   [R3] allo without any lineage signal: W = M . theta (as auto).
//   [R4] long-read Theta anchoring for both ploidies.
//   [R5] the Theta k-mer index lives in --data-dir and is shared by all samples.
//   [R6] global EM over read equivalence classes of all genes; effective lengths; identical-sequence pooling.
//   [R7] bowtie2 with RSEM settings is the default aligner; candidates weighted by mismatches. A sample directory
//        records its aligner (aligner.txt); switching aligner on an existing output directory is refused.
#include "prep.hpp"
#include "proc.hpp"

#include <dirent.h>
#include <limits.h>
#include <set>

namespace lq {

const char* lq_version();

namespace {

std::string fmt_dur(double s) {
    long long h = (long long)(s / 3600);
    double rem = s - 3600.0 * (double)h;
    long long m = (long long)(rem / 60);
    double sec = rem - 60.0 * (double)m;
    if (h > 0) return sfmt("%lldh %lldm %.1fs", h, m, sec);
    if (m > 0) return sfmt("%lldm %.1fs", m, sec);
    return sfmt("%.2fs", sec);
}
std::string abspath(const std::string& p) {
    if (p.empty()) return p;
    std::string a = p;
    if (p[0] != '/') {
        char cwd[PATH_MAX];
        if (!getcwd(cwd, sizeof cwd)) die("getcwd failed");
        a = std::string(cwd) + "/" + p;
    }
    // normalise like os.path.abspath: drop ".", "", resolve ".."
    std::vector<std::string> parts;
    size_t b = 0;
    while (b <= a.size()) {
        size_t e = a.find('/', b);
        std::string x = a.substr(b, e == std::string::npos ? std::string::npos : e - b);
        if (x == "..") { if (!parts.empty()) parts.pop_back(); }
        else if (!x.empty() && x != ".") parts.push_back(x);
        if (e == std::string::npos) break;
        b = e + 1;
    }
    std::string o;
    for (auto& x : parts) o += "/" + x;
    return o.empty() ? "/" : o;
}
std::string self_exe() {
    char b[PATH_MAX];
    ssize_t n = readlink("/proc/self/exe", b, sizeof b - 1);
    if (n <= 0) die("cannot find this program (/proc/self/exe)");
    b[n] = 0;
    return b;
}
bool lexists(const std::string& p) { struct stat st; return lstat(p.c_str(), &st) == 0; }
// true when `out` is missing or older than `in` (text export to redo)
bool stale(const std::string& out, const std::string& in) {
    struct stat a, b;
    if (stat(out.c_str(), &a) != 0) return true;
    if (stat(in.c_str(), &b) != 0) return false;
    return a.st_mtim.tv_sec < b.st_mtim.tv_sec || (a.st_mtim.tv_sec == b.st_mtim.tv_sec && a.st_mtim.tv_nsec < b.st_mtim.tv_nsec);
}
std::string base_name(const std::string& p) { size_t s = p.find_last_of('/'); return s == std::string::npos ? p : p.substr(s + 1); }
std::string replace_all(std::string s, const std::string& a, const std::string& b) {
    size_t p = 0;
    while ((p = s.find(a, p)) != std::string::npos) { s.replace(p, a.size(), b); p += b.size(); }
    return s;
}
std::string read_first_line(const std::string& p) {
    FILE* f = fopen(p.c_str(), "r");
    if (!f) return {};
    char* l = nullptr;
    size_t cap = 0;
    std::string s;
    if (getline(&l, &cap, f) > 0) s = l;
    free(l);
    fclose(f);
    while (!s.empty() && (s.back() == '\n' || s.back() == '\r')) s.pop_back();
    return s;
}

// ---------------------------------------------------------------- options
struct Cfg {
    std::string fasta, gff, clusters, ancestry, orphan, target_ids, ploidy = "allo";
    std::vector<std::pair<std::string, std::vector<std::pair<std::string, std::string>>>> ancestors;   // name -> [(fa, gff)]
    std::vector<std::pair<std::string, std::string>> markers;
    bool marker_clean = false, msa_snp = false, legacy_lineages = false, text = true, dry = false;
    std::string vcf, phased_fasta, longread;
    std::vector<std::string> sample_lists;
    std::string data_dir = "./data", out_dir = "./results";
    std::string aligner = "bowtie2", salmon = "salmon";
    int max_aln = 200, threads = 16, kmer = 21, lineage_kmer = 31, min_lr_depth = 5;
    std::string shrinkage = "10.0", error_rate = "0.01";
    bool adaptive = false, no_bonus = false, no_chimeric = false, no_pool = false, gc_bias = false;
    std::string em_accel = "none", nm_weight = "auto", eff_len = "eff", tpm_length;
    std::string pw_spec = "1.0", pw_conf = "0.0", pw_none = "0.0", p_alpha = "0.0", p_floor = "0.001";
    std::string tp_lr = "1.0", tp_kmer = "0.0", tp_cap = "0.0";
    std::string theta_sum = "neumaier";
    std::string engine = "cpp", em_engine = "global";
};

void usage() {
    std::fprintf(stderr,
        "usage: lq run --fasta mRNA.fa --clusters CLUSTERS --ancestry ANCESTRY [--target-ids IDS] [options]\n"
        "\n"
        "required reference files\n"
        "  --fasta FILE            transcriptome (mRNA) FASTA: every allele copy\n"
        "  --clusters FILE         allele clusters ClusterID<TAB>g1,g2,... (lq normalize-clusters)\n"
        "  --ancestry FILE         GeneID<TAB>Lineage\n"
        "  --target-ids FILE       target gene ids, one per line (allo)\n"
        "  --orphan FILE           single-copy genes (accepted; used by `lq matrix`, not by the run)\n"
        "\n"
        "ploidy + lineage inputs\n"
        "  --ploidy auto|allo      auto = copy level only; allo (default) = lineage assignment\n"
        "  --ancestor NAME:fa,gff  progenitor genome for k-mer markers (repeatable; the same NAME twice = several samples)\n"
        "  --markers NAME:fasta    pre-built lineage marker FASTA (repeatable; skips the KMC build)\n"
        "  --marker-clean          drop marker k-mers that occur in a reference copy of another lineage (1.1.0; off)\n"
        "  --vcf FILE --gff FILE --sample-lists NAME:file ...   population-VCF SNP route\n"
        "  --msa-snp               MSA SNP route on the copies in --fasta;  --phased-fasta FILE: on a separate CDS FASTA\n"
        "  --legacy-lineages       MSA SNP database with the lineage list of LineageQuant 1.0.0 (see snpdb-msa -h)\n"
        "  --longread BAM          long-read alignments (minimap2 -a -N 50 -p 0.95 to the transcriptome) anchor Theta\n"
        "\n"
        "input / output\n"
        "  -d, --data-dir DIR      paired FASTQ (<s>_R1.fq.gz / <s>_1.fq.gz / ...) and the shared indexes (./data)\n"
        "  -o, --out-dir DIR       one directory per sample: Global_Theta_Matrix.lqb, <s>.lqres (./results)\n"
        "  --no-text               binary outputs only. Default: after each sample the 1.0.0 text files are written as well\n"
        "                          (EM_Final_Expression.tsv, <s>.idxstats.tsv, <s>.efflen.tsv, Global_Theta_Matrix.tsv;\n"
        "                          lq export), outside the timed quantification\n"
        "\n"
        "alignment\n"
        "  --aligner bowtie2|lqmap|sa   bowtie2 (default; RSEM settings, streamed, no BAM) | lqmap (built in, same rules)\n"
        "                          | sa (Salmon selective alignment proposes, ungapped re-verification)\n"
        "  --salmon PATH           salmon for --aligner sa;  --max-alignments N   bowtie2 -k (200)\n"
        "\n"
        "runtime / tuning\n"
        "  -t, --threads N (16)    --kmer K (21, Theta)    --lineage-kmer K (31)    --shrinkage X (10.0)\n"
        "  --adaptive-shrinkage    --no-paired-bonus       --min-lr-depth N (5)     --error-rate P (0.01)   --no-chimeric\n"
        "  --em-accel none|squarem      none (default) = the EM of 1.0.0; squarem: same maximum, far fewer iterations\n"
        "  --theta-sum compensated|plain   float sums of Theta: compensated (default) = 1.0.0 under Python >= 3.12,\n"
        "                          plain = 1.0.0 under Python <= 3.11 (last-digit differences only)\n"
        "\n"
        "EM\n"
        "  --nm-weight auto|off|X  --no-pool-identical   --eff-len off|len|eff|fld (eff)   --gc-bias (1.1.0; off)\n"
        "  --p-weight-specific 1.0 --p-weight-conflict 0.0 --p-weight-none 0.0 --p-alpha 0.0 --p-floor 1e-3\n"
        "  --theta-prior-lr 1.0    --theta-prior-kmer 0.0  --theta-prior-cap 0.0\n"
        "\n"
        "  --dry-run               print the route and exit\n"
        "\n"
        "Afterwards: lq matrix -r OUT_DIR -o ORPHANS [--tpm-length full|eff]  (expression matrix)\n");
}

Cfg parse(int argc, char** argv) {
    Cfg c;
    auto need = [&](int& i) -> std::string {
        if (i + 1 >= argc) { usage(); die(std::string("missing value for ") + argv[i]); }
        return argv[++i];
    };
    auto many = [&](int& i, std::vector<std::string>& v) {
        int n0 = (int)v.size();
        while (i + 1 < argc && argv[i + 1][0] != '-') v.push_back(argv[++i]);
        if ((int)v.size() == n0) { usage(); die(std::string(argv[i]) + " needs at least one value"); }
    };
    auto num = [](const std::string& k, const std::string& v) {
        double x;
        if (!py_float(v, x)) die(k + ": not a number: " + v);
        return v;
    };
    auto integer = [](const std::string& k, const std::string& v) {
        long long x;
        if (!py_int(v, x)) die(k + ": not an integer: " + v);
        return (int)x;
    };
    for (int i = 0; i < argc; ++i) {
        std::string a = argv[i];
        std::string v;
        size_t eq = a.find('=');
        if (a.rfind("--", 0) == 0 && eq != std::string::npos) { v = a.substr(eq + 1); a = a.substr(0, eq); }
        auto val = [&]() { return !v.empty() ? v : need(i); };
        if (a == "--fasta") c.fasta = val();
        else if (a == "--gff") c.gff = val();
        else if (a == "--clusters") c.clusters = val();
        else if (a == "--ancestry") c.ancestry = val();
        else if (a == "--orphan") c.orphan = val();
        else if (a == "--target-ids") c.target_ids = val();
        else if (a == "--ploidy") { c.ploidy = val(); if (c.ploidy != "auto" && c.ploidy != "allo") die("--ploidy must be auto or allo"); }
        else if (a == "--ancestor") {
            std::string it = val();
            size_t p = it.find(':');
            if (p == std::string::npos || it.find(',') == std::string::npos) die("--ancestor must be NAME:fasta,gff, got: " + it);
            std::string name = it.substr(0, p), paths = it.substr(p + 1);
            size_t q = paths.find(',');
            if (q == std::string::npos) die("--ancestor must be NAME:fasta,gff, got: " + it);
            auto pr = std::make_pair(abspath(paths.substr(0, q)), abspath(paths.substr(q + 1)));
            bool found = false;
            for (auto& x : c.ancestors) if (x.first == name) { x.second.push_back(pr); found = true; }
            if (!found) c.ancestors.push_back({name, {pr}});
        }
        else if (a == "--markers") {
            std::string it = val();
            size_t p = it.find(':');
            if (p == std::string::npos) die("--markers must be NAME:fasta, got: " + it);
            std::string name(py_strip(it.substr(0, p)));
            if (name.empty()) die("--markers needs a lineage name before ':', got: " + it);
            for (auto& x : c.markers) if (x.first == name) die("--markers: lineage '" + name + "' given more than once");
            c.markers.push_back({name, abspath(it.substr(p + 1))});
        }
        else if (a == "--marker-clean") c.marker_clean = true;
        else if (a == "--vcf") c.vcf = val();
        else if (a == "--sample-lists") many(i, c.sample_lists);
        else if (a == "--phased-fasta") c.phased_fasta = val();
        else if (a == "--msa-snp") c.msa_snp = true;
        else if (a == "--legacy-lineages") c.legacy_lineages = true;
        else if (a == "--longread") c.longread = val();
        else if (a == "-d" || a == "--data-dir") c.data_dir = val();
        else if (a == "-o" || a == "--out-dir") c.out_dir = val();
        else if (a == "--text") c.text = true;
        else if (a == "--no-text") c.text = false;
        else if (a == "--aligner") { c.aligner = val(); if (c.aligner != "bowtie2" && c.aligner != "lqmap" && c.aligner != "sa") die("--aligner must be bowtie2, lqmap or sa"); }
        else if (a == "--salmon") c.salmon = val();
        else if (a == "--max-alignments") c.max_aln = integer(a, val());
        else if (a == "-t" || a == "--threads") c.threads = integer(a, val());
        else if (a == "--kmer") c.kmer = integer(a, val());
        else if (a == "--lineage-kmer") c.lineage_kmer = integer(a, val());
        else if (a == "--shrinkage") c.shrinkage = num(a, val());
        else if (a == "--adaptive-shrinkage") c.adaptive = true;
        else if (a == "--no-paired-bonus") c.no_bonus = true;
        else if (a == "--min-lr-depth") c.min_lr_depth = integer(a, val());
        else if (a == "--error-rate") c.error_rate = num(a, val());
        else if (a == "--no-chimeric") c.no_chimeric = true;
        else if (a == "--em-accel") { c.em_accel = val(); if (c.em_accel != "none" && c.em_accel != "squarem") die("--em-accel must be none or squarem"); }
        else if (a == "--theta-sum" || a == "--py-sum") {
            std::string s = val();
            if (s == "compensated" || s == "neumaier") c.theta_sum = "neumaier";
            else if (s == "plain" || s == "naive") c.theta_sum = "naive";
            else die("--theta-sum must be compensated or plain");
        }
        else if (a == "--nm-weight") c.nm_weight = val();
        else if (a == "--no-pool-identical") c.no_pool = true;
        else if (a == "--eff-len") { c.eff_len = val(); if (c.eff_len != "off" && c.eff_len != "len" && c.eff_len != "eff" && c.eff_len != "fld") die("--eff-len must be off, len, eff or fld"); }
        else if (a == "--gc-bias") c.gc_bias = true;
        else if (a == "--tpm-length") c.tpm_length = val();
        else if (a == "--p-weight-specific") c.pw_spec = num(a, val());
        else if (a == "--p-weight-conflict") c.pw_conf = num(a, val());
        else if (a == "--p-weight-none") c.pw_none = num(a, val());
        else if (a == "--p-alpha") c.p_alpha = num(a, val());
        else if (a == "--p-floor") c.p_floor = num(a, val());
        else if (a == "--theta-prior-lr") c.tp_lr = num(a, val());
        else if (a == "--theta-prior-kmer") c.tp_kmer = num(a, val());
        else if (a == "--theta-prior-cap") c.tp_cap = num(a, val());
        else if (a == "--dry-run") c.dry = true;
        // options of 1.0.0 that still parse
        else if (a == "--engine") c.engine = val();
        else if (a == "--em-engine") c.em_engine = val();
        else if (a == "--lqcore" || a == "--theta-threads" || a == "--calib-min-depth" || a == "--calib-det-weight" ||
                 a == "--em-epsilon" || a == "--em-max-iter" || a == "--em-min-iter" || a == "--em-tol" || a == "--em-ll-tol" ||
                 a == "--cross-cluster") (void)val();       // legacy per-cluster EM / SNP calibration: no effect
        else if (a == "--no-calibrate-snp") {}
        else if (a == "--auto-idxstats")
            die("--auto-idxstats (ablation without Theta and EM) is not a route of 1.1.0: run the pipeline, then "
                "lq matrix --mode idxstats (every reported alignment counted in full)");
        else if (a == "-h" || a == "--help") { usage(); std::exit(0); }
        else if (a == "--r1" || a == "--sam-in" || a == "-a" || a == "-x") {
            usage();
            die(a + " is an option of one sample: `lq quant` (= lqcore run). `lq run` is the whole pipeline.");
        }
        else { usage(); die("unrecognized argument: " + std::string(argv[i])); }
    }
    if (c.fasta.empty() || c.clusters.empty() || c.ancestry.empty()) { usage(); die("--fasta, --clusters and --ancestry are required"); }
    if (c.ploidy == "allo" && c.target_ids.empty()) { usage(); die("--target-ids is required (allo)"); }
    // paths
    for (auto* p : {&c.fasta, &c.gff, &c.clusters, &c.ancestry, &c.orphan, &c.target_ids, &c.vcf, &c.phased_fasta, &c.longread,
                    &c.data_dir, &c.out_dir})
        *p = abspath(*p);
    for (auto& s : c.sample_lists) {
        size_t p = s.find(':');
        if (p == std::string::npos) die("--sample-lists must be NAME:path, got: " + s);
        s = s.substr(0, p) + ":" + abspath(s.substr(p + 1));
    }
    if (!c.phased_fasta.empty() && c.msa_snp)
        die("--phased-fasta and --msa-snp are two ways to trigger the same MSA route; give only one. Use --msa-snp if the main "
            "--fasta already holds all allele copies; use --phased-fasta only for a separate CDS file.");
    if (c.msa_snp) c.phased_fasta = c.fasta;
    // validation (config.py)
    if (c.engine == "python") die("--engine python was removed in 1.1.0 (the C++ engine gives the same results); use LineageQuant 1.0.0");
    if (c.em_engine != "global" && c.em_engine != "v6") die("--em-engine legacy (per-cluster EM) was removed in 1.1.0; use LineageQuant 1.0.0");
    if (c.max_aln < 1) die("--max-alignments must be >= 1");
    if (c.nm_weight != "auto" && c.nm_weight != "off") {
        double v;
        if (!py_float(c.nm_weight, v)) die("--nm-weight must be auto, off or a number, got: " + c.nm_weight);
        if (!(v > 0.0 && v <= 1.0)) die("--nm-weight as a number must be in (0, 1]");
    }
    if (c.gc_bias && c.eff_len == "off") die("--gc-bias corrects effective lengths: not with --eff-len off");
    if (c.ploidy == "allo") {
        bool kmer = !c.ancestors.empty() || !c.markers.empty();
        if (!kmer && c.vcf.empty() && c.phased_fasta.empty())
            die("Allopolyploid mode needs at least one lineage route: --ancestor (k-mer), --markers (pre-built k-mer), "
                "--vcf (population SNP), or --phased-fasta / --msa-snp (MSA SNP).");
        if (!c.vcf.empty() && !c.phased_fasta.empty()) die("--vcf and --phased-fasta are mutually exclusive (both provide SNP signal). Pick one.");
        if (!c.vcf.empty() && c.sample_lists.empty()) die("--vcf requires --sample-lists NAME:path [...]");
        if (!c.vcf.empty() && c.gff.empty()) die("--vcf requires --gff");
    }
    if (c.marker_clean && c.ancestors.empty() && c.markers.empty()) die("--marker-clean needs lineage markers (--ancestor or --markers)");
    return c;
}

std::string describe(const Cfg& c, const std::string& exe) {
    std::vector<std::string> L = {"ploidy: " + c.ploidy};
    if (c.aligner == "bowtie2")
        L.push_back("  alignment: bowtie2 (RSEM settings, N not penalised, <= " + std::to_string(c.max_aln) +
                    " alignments/fragment), streamed into quant, no BAM");
    else if (c.aligner == "lqmap")
        L.push_back("  alignment: lqmap (built in: exact k-mer seeds, ungapped mismatch count, concordant pairs, candidates within best + 3 mismatches)");
    else
        L.push_back("  alignment: sa (" + c.salmon + " selective alignment -> ungapped re-verification, candidates within best + 3 mismatches)");
    L.push_back("  engine: " + exe + (c.em_accel != "none" ? " | EM acceleration: " + c.em_accel : std::string()));
    L.push_back("  EM: global over all genes | eff-len=" + c.eff_len + (c.gc_bias ? " + GC bias" : "") + " | mismatch weighting=" +
                c.nm_weight + " | identical-sequence pooling=" + (c.no_pool ? "off" : "on"));
    L.push_back("  outputs: <sample>/Global_Theta_Matrix.lqb, <sample>/<sample>.lqres" +
                std::string(c.text ? " + the 1.0.0 text files (written after each sample, not timed)" : " (binary only; lq export)"));
    if (c.ploidy == "auto") {
        L.push_back(std::string("  copy-level: Theta (") + (c.longread.empty() ? "k-mer" : "long-read anchored") + ") -> EM");
    } else {
        std::vector<std::string> r;
        if (!c.markers.empty()) r.push_back("pre-built k-mer markers (" + std::to_string(c.markers.size()) + " lineages, build skipped)");
        if (!c.ancestors.empty()) {
            size_t n = 0;
            for (auto& x : c.ancestors) n += x.second.size();
            r.push_back("ancestor k-mer (" + std::to_string(c.ancestors.size()) + " lineages, " + std::to_string(n) + " progenitor sample(s))");
        }
        if (c.marker_clean) r.push_back("markers cleaned against the reference copies");
        if (!c.vcf.empty()) r.push_back("population VCF SNP");
        if (!c.phased_fasta.empty()) r.push_back(std::string("MSA SNP") + (c.msa_snp ? " [main fasta]" : " [separate phased CDS]") +
                                                 (c.legacy_lineages ? " [1.0.0 lineage list]" : ""));
        std::string s;
        for (auto& x : r) s += (s.empty() ? "" : ", ") + x;
        L.push_back("  lineage routes: " + s);
        if ((!c.ancestors.empty() || !c.markers.empty()) && (!c.vcf.empty() || !c.phased_fasta.empty()))
            L.push_back("  -> PoE fusion of k-mer + SNP signals");
        L.push_back(std::string("  Theta: ") + (c.longread.empty() ? "k-mer (01)" : "long-read anchored (01_hybrid)"));
    }
    std::string o;
    for (auto& x : L) o += x + "\n";
    return o;
}

// ---------------------------------------------------------------- steps
void run_step(const std::vector<std::string>& cmd, const std::string& name, const std::string& done) {
    if (file_exists(done)) { logm("[SKIP] " + name); return; }
    std::printf("\n");
    logm("[RUN] " + name);
    std::printf("  $ %s\n", join_cmd(cmd).c_str());
    std::fflush(stdout);
    const double t0 = wall();
    RunResult r = run_proc(cmd);
    if (r.status != 0) {
        std::printf("\n");
        logm("[FATAL] " + name + " failed (exit " + std::to_string(r.status) + ")");
        std::exit(1);
    }
    touch_done(done);
    logm("[OK] " + name + " in " + fmt_dur(wall() - t0));
}
// a step done in this process (no command line)
template <class F>
void run_inline(const std::string& name, const std::string& done, F&& f) {
    if (file_exists(done)) { logm("[SKIP] " + name); return; }
    std::printf("\n");
    logm("[RUN] " + name);
    const double t0 = wall();
    f();
    touch_done(done);
    logm("[OK] " + name + " in " + fmt_dur(wall() - t0));
}

struct Sample { std::string name, r1, r2; };
// router.discover_samples(): glob patterns in order (directory order inside a pattern), first suffix that fits
std::vector<Sample> discover_samples(const std::string& dir) {
    const std::vector<std::string> suf = {".R1.fq.gz", "_R1.fq.gz", "_1.fq.gz", ".1.fq.gz", "_R1.fastq.gz", "_1.fastq.gz",
                                          ".1.fastq.gz", "_R1_001.fastq.gz"};
    const std::vector<std::pair<std::string, std::string>> r1r2 = {{".R1.", ".R2."}, {"_R1.", "_R2."}, {"_1.", "_2."},
                                                                   {".1.", ".2."}, {"_R1_", "_R2_"}};
    std::vector<std::string> names;
    if (DIR* d = opendir(dir.c_str())) {
        while (dirent* e = readdir(d)) names.push_back(e->d_name);
        closedir(d);
    }
    std::vector<std::string> r1s;
    std::set<std::string> seen;
    for (auto& s : suf)
        for (auto& n : names) {
            if (n.empty() || n[0] == '.') continue;   // glob skips hidden files
            if (n.size() >= s.size() && n.compare(n.size() - s.size(), s.size(), s) == 0) {
                std::string p = dir + "/" + n;
                if (seen.insert(p).second) r1s.push_back(p);
            }
        }
    std::vector<Sample> out;
    for (auto& r1 : r1s) {
        std::string b = base_name(r1), sname;
        bool got = false;
        for (auto& s : suf) if (b.size() >= s.size() && b.compare(b.size() - s.size(), s.size(), s) == 0) { sname = b.substr(0, b.size() - s.size()); got = true; break; }
        if (!got) sname = b.find('.') != std::string::npos ? b.substr(0, b.rfind('.')) : b;
        std::string r2;
        for (auto& pr : r1r2) {
            std::string c = replace_all(r1, pr.first, pr.second);
            if (c != r1 && file_exists(c)) { r2 = c; break; }
        }
        if (!r2.empty()) out.push_back({sname, r1, r2});
        else logm("[WARN] no R2 for " + r1 + "; skipping");
    }
    return out;
}

// lineages of an MSA SNP database vs the ancestry file (a database built by 1.0.0 may lack some)
void check_snpdb_lineages(const std::string& db, const std::string& ancestry) {
    std::string h = read_first_line(db);
    size_t p = h.find("lineages=");
    if (p == std::string::npos) return;
    std::string ls = h.substr(p + 9);
    ls = ls.substr(0, ls.find('\t'));
    std::set<std::string> have;
    for (auto x : py_split(ls, ',')) have.insert(std::string(x));
    std::set<std::string> want;
    LineReader rd(ancestry);
    std::string_view l;
    while (rd.next(l)) {
        if (py_strip(l).empty() || l.rfind("GeneID", 0) == 0) continue;
        auto f = py_split_ws(py_strip(l));
        if (f.size() >= 2) want.insert(upper(f[1]));
    }
    std::string miss;
    for (auto& w : want) if (!have.count(w)) miss += (miss.empty() ? "" : ", ") + w;
    if (!miss.empty())
        logm("[WARN] " + db + " has no SNP column for lineage(s) " + miss + " (a database built by LineageQuant 1.0.0, which listed "
             "the lineages of the first SNP only). Delete it and its .done file to rebuild with every lineage.");
}

}  // namespace

int pipeline_main(int argc, char** argv) {
    Cfg c = parse(argc, argv);
    const std::string exe = self_exe();
    std::printf("[LineageQuant %s] resolved route:\n%s", lq_version(), describe(c, exe).c_str());
    if (c.dry) { std::printf("[dry-run] exiting without execution.\n"); return 0; }
    if (!c.tpm_length.empty())
        logm("[INFO] --tpm-length belongs to the expression matrix, which is no longer part of the run: lq matrix --tpm-length " + c.tpm_length);
    const double T0 = wall();
    const std::string& data = c.data_dir;
    const std::string& out = c.out_dir;
    mkdirs(data);
    mkdirs(out);
    const bool allo = c.ploidy == "allo";
    const std::string th = std::to_string(c.threads);
    logm(std::string("[INIT] LineageQuant ") + lq_version() + " | ploidy=" + c.ploidy + " (" + (allo ? "allopolyploid" : "autopolyploid") +
         ") | aligner=" + c.aligner + " | threads=" + th);
    if (!allo) logm("[INFO] auto route: mapping -> Theta -> EM; no lineage axis.");

    // ---- external tools
    std::vector<std::string> tools;
    if (c.aligner == "bowtie2") tools = {"bowtie2", "bowtie2-build"};
    else if (c.aligner == "sa") tools = {c.salmon};
    if (allo && !c.phased_fasta.empty()) tools.push_back("mafft");
    if (allo && !c.ancestors.empty()) { tools.push_back("kmc"); tools.push_back("kmc_tools"); tools.push_back("gffread"); }
    if (!c.longread.empty()) tools.push_back("samtools");
    for (auto& t : std::set<std::string>(tools.begin(), tools.end()))
        if (which(t).empty()) die("tool '" + t + "' not on PATH");
    for (auto& f : {c.fasta, c.clusters, c.ancestry}) if (!file_exists(f)) die("not found: " + f);

    // ---- transcriptome index of the aligner (shared by all samples)
    std::string index;
    if (c.aligner == "bowtie2") {
        index = data + "/transcript_index_bt2";
        if (!file_exists(index + ".1.bt2") && !file_exists(index + ".1.bt2l"))
            run_step({"bowtie2-build", "--threads", th, c.fasta, index}, "bowtie2 build index", index + ".done");
    } else if (c.aligner == "sa") {
        index = data + "/salmon_kd_index";
        run_step({c.salmon, "index", "-t", c.fasta, "-i", index, "-k", "31", "-p", th, "--keepDuplicates"},
                 "salmon build index (--keepDuplicates)", index + ".done");
    } else index = data + "/lqmap.k21.bin";   // built and cached by quant

    std::vector<Sample> samples = discover_samples(data);
    if (samples.empty()) die("no FASTQ sample pairs found in " + data);
    {
        std::string s;
        for (auto& x : samples) s += (s.empty() ? "" : ", ") + x.name;
        logm("[INFO] " + std::to_string(samples.size()) + " sample(s): " + s);
    }

    // ---- [R4] long-read prior (both ploidies)
    std::string lr_prior;
    if (!c.longread.empty()) {
        lr_prior = data + "/LongRead_Counts.tsv";
        run_step({exe, "lrcount", "--bam", c.longread, "-o", lr_prior, "--threads", std::to_string(std::min(8, c.threads)),
                  "--legacy-out", data + "/LongRead_Counts.legacy_idxstats.tsv", "--summary", data + "/LongRead_Counts.summary.tsv"},
                 "Step 00b - long-read priors", lr_prior + ".done");
    }

    // ---- allo: SNP database
    std::string snp_db;
    if (allo && !c.vcf.empty()) {
        snp_db = data + "/SNP_DB_population.lq.tsv";
        std::vector<std::string> cmd = {exe, "snpdb-vcf", "--vcf", c.vcf, "--gff3", c.gff, "--out", snp_db, "--sample-lists"};
        cmd.insert(cmd.end(), c.sample_lists.begin(), c.sample_lists.end());
        run_step(cmd, "Step 1a - population SNP DB", snp_db + ".done");
    } else if (allo && !c.phased_fasta.empty()) {
        snp_db = data + "/SNP_DB_msa.lq.tsv";
        std::vector<std::string> cmd = {exe, "snpdb-msa", "--fasta", c.phased_fasta, "--cluster", c.clusters, "--ancestry", c.ancestry,
                                        "--out", snp_db, "--threads", std::to_string(std::max(1, c.threads / 2))};
        if (c.legacy_lineages) cmd.push_back("--legacy-lineages");
        run_step(cmd, "Step 1a - MSA SNP DB", snp_db + ".done");
        if (!c.legacy_lineages) check_snpdb_lineages(snp_db, c.ancestry);
    }

    // ---- allo: lineage k-mer markers
    std::vector<std::pair<std::string, std::string>> marker_paths;
    if (allo && (!c.ancestors.empty() || !c.markers.empty())) {
        if (!c.markers.empty()) {
            for (auto& m : c.markers) {
                if (!nonempty(m.second)) die("--markers " + m.first + ": file not found or empty: " + m.second);
                marker_paths.push_back(m);
            }
            logm("[INFO] using " + std::to_string(marker_paths.size()) + " pre-built lineage marker(s); skipping KMC build");
        } else {
            const std::string mdir = data + "/lineage_markers", mdone = mdir + "/markers.done";
            if (!file_exists(mdone)) {
                const std::string adir = data + "/ancestor_transcripts";
                mkdirs(adir);
                std::vector<std::pair<std::string, std::string>> lf;
                for (auto& a : c.ancestors) {
                    const std::string merged = adir + "/" + a.first + "_transcript.fa";
                    if (!nonempty(merged)) {
                        std::vector<std::string> parts;
                        for (size_t si = 0; si < a.second.size(); ++si) {
                            const std::string tag = a.first + "_s" + std::to_string(si), part = adir + "/_" + tag + "_transcript.fa";
                            if (!nonempty(part)) genome_pair_extract(a.second[si].first, a.second[si].second, adir, tag, part, a.first);
                            parts.push_back(part);
                        }
                        concat_files(parts, merged);
                        long long n = 0;
                        LineReader rd(merged);
                        std::string_view l;
                        while (rd.next(l)) if (!l.empty() && l[0] == '>') ++n;
                        logm("[INFO] lineage " + a.first + ": " + std::to_string(a.second.size()) + " progenitor sample(s) -> " +
                             std::to_string(n) + " transcripts");
                    }
                    lf.push_back({a.first, merged});
                }
                std::vector<std::string> order;
                auto paths = build_markers(lf, mdir, c.lineage_kmer, c.threads, 0.8, {}, "ancestor", false, &order);
                for (auto& n : order) marker_paths.push_back({n, paths[n]});
                touch_done(mdone);
            } else {
                for (auto& a : c.ancestors) {
                    std::string p = mdir + "/" + a.first + ".fasta";
                    if (file_exists(p)) marker_paths.push_back({a.first, p});
                }
            }
        }
    }
    std::string marker_cache;
    if (!marker_paths.empty()) {
        marker_cache = data + "/lqcore_markers.k" + std::to_string(c.lineage_kmer) + (c.marker_clean ? ".clean" : "") + ".bin";
        std::vector<std::string> cmd = {exe, "markers", "--kmer-len", std::to_string(c.lineage_kmer), "--marker-cache", marker_cache, "--threads", th};
        for (auto& m : marker_paths) { cmd.push_back("--markers"); cmd.push_back(m.first + ":" + m.second); }
        if (c.marker_clean) { cmd.push_back("--marker-clean"); cmd.push_back("--fasta"); cmd.push_back(c.fasta); cmd.push_back("-a"); cmd.push_back(c.ancestry); }
        run_step(cmd, std::string("markers - marker index") + (c.marker_clean ? " (cleaned against the reference)" : ""), marker_cache + ".done");
    }

    // ---- per sample
    double t_samples = 0, t_export = 0;
    for (auto& S : samples) {
        const std::string so = out + "/" + S.name;
        mkdirs(so);
        std::printf("\n");
        logm("[SAMPLE] " + S.name);
        for (auto& lr : {std::make_pair(so + "/" + S.name + ".R1.fq.gz", S.r1), std::make_pair(so + "/" + S.name + ".R2.fq.gz", S.r2)})
            if (!lexists(lr.first) && symlink(abspath(lr.second).c_str(), lr.first.c_str()) != 0) die("cannot link " + lr.first);
        // [R7] never mix aligners inside one output directory
        const std::string tag = so + "/aligner.txt";
        std::string prev = file_exists(tag) ? std::string(py_strip(read_first_line(tag))) : std::string();
        if (!prev.empty() && prev != c.aligner)
            die(S.name + ": this output directory was aligned with " + prev + ", but --aligner " + c.aligner +
                " was requested. Use a new --out-dir (or --aligner " + prev + " to resume).");
        if (!file_exists(tag)) write_text_file(tag, [&](FILE* f) { fprintf(f, "%s\n", c.aligner.c_str()); });
        const double ts = wall();

        // Step 01 Theta; [R5] the k-mer index is shared by all samples
        std::string theta = so + "/Global_Theta_Matrix.lqb";
        const std::string theta_txt = so + "/Global_Theta_Matrix.tsv";
        if (!file_exists(theta + ".done") && file_exists(theta_txt + ".done") && file_exists(theta_txt)) {
            theta = theta_txt;
            logm("[SKIP] Step 01 - Theta (" + S.name + "): text Theta of an earlier run");
        } else {
            std::vector<std::string> tc = {exe, "theta", "-c", c.clusters, "-f", c.fasta, "-r1", S.r1, "-r2", S.r2, "-k", std::to_string(c.kmer),
                                           "-t", th, "-o", so, "--shrinkage", c.shrinkage, "--index_dir", data,
                                           "--index-cache", data + "/lqcore_theta.k" + std::to_string(c.kmer) + ".bin", "--py-sum", c.theta_sum, "--binary"};
            if (!lr_prior.empty()) { tc.push_back("--longread_prior"); tc.push_back(lr_prior); tc.push_back("--min_lr_depth"); tc.push_back(std::to_string(c.min_lr_depth)); }
            if (c.adaptive) tc.push_back("--adaptive_shrinkage");
            if (c.no_bonus) tc.push_back("--no_paired_bonus");
            run_step(tc, "Step 01 - Theta (" + S.name + ")", theta + ".done");
        }

        // Steps 00, 02, 02b, 03 in one pass
        const std::string res = so + "/" + S.name + ".lqres";
        const std::string em_txt = so + "/" + (c.longread.empty() ? "EM_Final_Expression.tsv" : "EM_Final_Hybrid.tsv");
        if (!file_exists(res + ".done") && file_exists(em_txt + ".lqcore.done")) {
            logm("[SKIP] Step 00-03 (" + S.name + "): text results of an earlier run");
        } else {
            std::vector<std::string> q = {exe, "quant", "--r1", S.r1, "--r2", S.r2};
            if (c.aligner == "lqmap") q.insert(q.end(), {"--aligner", "lqmap", "--map-cache", index});
            else if (c.aligner == "sa") q.insert(q.end(), {"--aligner", "sa", "--salmon-index", index, "--salmon", c.salmon, "--bt2-log", so + "/" + S.name + ".salmon.log"});
            else q.insert(q.end(), {"-x", index, "-k", std::to_string(c.max_aln), "--bt2-log", so + "/" + S.name + ".bowtie2.log"});
            q.insert(q.end(), {"--threads", th, "-a", c.ancestry, "-t", theta, "-o", res, "--nm-weight", c.nm_weight});
            const bool needs_fasta = c.aligner == "lqmap" || c.aligner == "sa" || c.gc_bias || c.marker_clean;
            if (!c.no_pool) q.insert(q.end(), {"--fasta", c.fasta});
            else if (needs_fasta) q.insert(q.end(), {"--fasta", c.fasta, "--no-pool"});
            if (c.em_accel != "none") q.insert(q.end(), {"--em-accel", c.em_accel});
            if (c.gc_bias) q.push_back("--gc-bias");
            q.insert(q.end(), {"--theta-prior-lr", c.tp_lr, "--theta-prior-kmer", c.tp_kmer, "--theta-prior-cap", c.tp_cap, "--eff-len", c.eff_len,
                               "--p-weight-specific", c.pw_spec, "--p-weight-conflict", c.pw_conf, "--p-weight-none", c.pw_none,
                               "--p-alpha", c.p_alpha, "--p-floor", c.p_floor});
            if (allo) {
                if (!c.target_ids.empty()) q.insert(q.end(), {"--target-ids", c.target_ids});
                for (auto& m : marker_paths) q.insert(q.end(), {"--markers", m.first + ":" + m.second});
                if (!marker_paths.empty()) {
                    q.insert(q.end(), {"--kmer-len", std::to_string(c.lineage_kmer), "--marker-cache", marker_cache});
                    if (c.marker_clean) q.push_back("--marker-clean");
                }
                if (!snp_db.empty()) {
                    q.insert(q.end(), {"--snp-db", snp_db, "--error-rate", c.error_rate});
                    if (c.no_chimeric) q.push_back("--no-chimeric");
                }
                if (marker_paths.empty() && snp_db.empty()) logm("[WARN] allo but no marker/SNP signal -> W = M . theta (same as auto)");
            }
            run_step(q, "Step 00-03 - quant: " + c.aligner + " + lineage + EM (" + S.name + ")", res + ".done");
        }
        t_samples += wall() - ts;
        logm("[SAMPLE DONE] " + S.name + " in " + fmt_dur(wall() - ts));
        if (c.text) {   // the text files of 1.0.0 (human-readable); not counted in the quantification time
            const double tx = wall();
            std::vector<std::string> made;
            const std::string idx_txt = so + "/" + S.name + ".idxstats.tsv", eff_txt = so + "/" + S.name + ".efflen.tsv";
            if (file_exists(res) && (stale(em_txt, res) || stale(idx_txt, res) || stale(eff_txt, res))) {
                RunResult r = run_proc({exe, "export", res, "-o", so, "--name", S.name, "--em-name", base_name(em_txt)}, "-", "");
                if (r.status != 0) die("export of " + res + " failed");
                made.insert(made.end(), {base_name(em_txt), base_name(idx_txt), base_name(eff_txt)});
            }
            if (ends_with(theta, ".lqb") && file_exists(theta) && stale(theta_txt, theta)) {
                RunResult r = run_proc({exe, "export", theta, "-o", theta_txt}, "-", "");
                if (r.status != 0) die("export of " + theta + " failed");
                made.push_back(base_name(theta_txt));
            }
            t_export += wall() - tx;
            if (!made.empty()) {
                std::string m;
                for (auto& x : made) m += (m.empty() ? "" : ", ") + x;
                logm("[TEXT] " + S.name + ": " + m + sfmt(" (%.1fs, not counted)", wall() - tx));
            }
        }
    }
    std::printf("\n");
    const double total = wall() - T0;
    logm("[DONE] quantification " + fmt_dur(total - t_export) + " (samples " + fmt_dur(t_samples) + ")" +
         (c.text ? " | text files " + fmt_dur(t_export) + ", not counted | wall " + fmt_dur(total) : std::string()));
    std::printf("  results: %s/<sample>/<sample>.lqres%s\n", out.c_str(),
                c.text ? " (+ EM_Final_Expression.tsv, <sample>.idxstats.tsv, <sample>.efflen.tsv, Global_Theta_Matrix.tsv)" : "");
    std::printf("  expression matrix (not part of the run): %s matrix -r %s -o %s\n", exe.c_str(), out.c_str(),
                c.orphan.empty() ? "ORPHANS.tsv" : c.orphan.c_str());
    return 0;
}

}  // namespace lq
