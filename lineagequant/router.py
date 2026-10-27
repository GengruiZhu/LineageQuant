#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
router.py -- LineageQuant auto-routing engine (the "brain").

The user declares the ploidy type and supplies whatever inputs are available;
the router decides which modules run. Every external step is checkpointed with a
.done marker and is resumable.

  Step 00   short-read alignment to the transcriptome (bowtie2 by default; hisat2 optional)
  Step 00b  long-read priors (optional)
  Step 01   Theta: per-copy prior (k-mer, or long-read anchored)
  Step 1a/b SNP database (allo, VCF or MSA route) and long-read calibration
  Step 02   per-read lineage probabilities (allo: ancestor k-mer and/or SNP; PoE fusion)
  Step 03   EM: reads -> copies (global engine by default)
  Step 04   aggregation into the genome-wide matrix

Design notes (kept because each changed results):

  [R1] --ploidy auto runs Theta + EM. The pre-v4 auto route was
       mapping -> idxstats -> TPM, i.e. featureCounts -M -O; it is kept only as
       the --auto-idxstats ablation (requires --aligner hisat2).
  [R2] Step 03 always receives the BAM: the compatibility matrix comes from the
       aligner, the lineage step only supplies probabilities.
  [R3] allo without any lineage signal: no -p is passed (W = M . theta), instead of
       an empty P file that zeroed every allocation.
  [R4] Long-read Theta anchoring is used for both ploidies.
  [R5] The Theta k-mer index lives in --data-dir and is shared by all samples.
  [R6] Global EM engine (core/em_global.py): one EM over read equivalence classes of
       all genes, effective-length correction, identical-sequence pooling. Removes the
       double counting of reads shared between clusters / orphan genes.
  [R7] (1.0.0) bowtie2 with RSEM settings is the default aligner. HISAT2 reports only
       a subset of equally good alignments (for 0.01-0.6% of fragments it reports a
       single wrong copy), which the EM cannot repair; bowtie2 -k 200 reports all of
       them plus near-optimal ones, and the EM weights candidates by mismatch count.
       Ambiguous reference bases (N) are not penalised (--np 0 --n-ceil L,0,0.5).
       A sample directory records its aligner (aligner.txt); switching aligner on an
       existing output directory is refused rather than silently mixing BAMs.

Long-read data are an enhancement, not a route: if present they anchor Theta
(01_hybrid) and calibrate MSA SNP AFs (1b).
"""

import os
import re
import sys
import glob
import time
import shutil
import subprocess

_PKG = os.path.dirname(os.path.abspath(__file__))


def _mod(*parts):
    return os.path.join(_PKG, *parts)


def get_now():
    return time.strftime("%Y-%m-%d %H:%M:%S")


def fmt_dur(seconds):
    h, rem = divmod(seconds, 3600)
    m, s = divmod(rem, 60)
    if h > 0:
        return f"{int(h)}h {int(m)}m {s:.1f}s"
    if m > 0:
        return f"{int(m)}m {s:.1f}s"
    return f"{s:.2f}s"


def run_step(cmd, name, done_marker, env=None):
    if os.path.exists(done_marker):
        print(f"[{get_now()}] [SKIP] {name}")
        return
    print(f"\n[{get_now()}] [RUN] {name}")
    print(f"  $ {' '.join(str(x) for x in cmd)}")
    e = os.environ.copy()
    e.setdefault("PYTHONUNBUFFERED", "1")
    if env:
        e.update(env)
    t0 = time.time()
    try:
        subprocess.run(cmd, check=True, env=e)
    except subprocess.CalledProcessError as ex:
        print(f"\n[{get_now()}] [FATAL] {name} failed (exit {ex.returncode})")
        sys.exit(1)
    with open(done_marker, "w") as f:
        f.write("DONE")
    print(f"[{get_now()}] [OK] {name} in {fmt_dur(time.time() - t0)}")


def check_deps(required_tools, required_modules):
    for t in required_tools:
        if shutil.which(t) is None:
            print(f"[FATAL] tool '{t}' not on PATH")
            sys.exit(1)
    for m in required_modules:
        try:
            subprocess.run([sys.executable, "-c", f"import {m}"],
                           check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except subprocess.CalledProcessError:
            print(f"[FATAL] python module '{m}' missing")
            sys.exit(1)


# --------------------------------------------------------------------------- #
def discover_samples(data_dir):
    r1_patterns = ["*.R1.fq.gz", "*_R1.fq.gz", "*_1.fq.gz", "*.1.fq.gz",
                   "*_R1.fastq.gz", "*_1.fastq.gz", "*.1.fastq.gz", "*_R1_001.fastq.gz"]
    r1_files = []
    for pat in r1_patterns:
        for f in glob.glob(os.path.join(data_dir, pat)):
            if f not in r1_files:
                r1_files.append(f)
    suffixes = [".R1.fq.gz", "_R1.fq.gz", "_1.fq.gz", ".1.fq.gz",
                "_R1.fastq.gz", "_1.fastq.gz", ".1.fastq.gz", "_R1_001.fastq.gz"]
    r1_to_r2 = [(r"\.R1\.", ".R2."), (r"_R1\.", "_R2."), (r"_1\.", "_2."),
                (r"\.1\.", ".2."), (r"_R1_", "_R2_")]
    samples = []
    for r1 in r1_files:
        base = os.path.basename(r1)
        sname = next((base[:-len(s)] for s in suffixes if base.endswith(s)),
                     base.rsplit(".", 1)[0])
        r2 = None
        for pat, repl in r1_to_r2:
            cand = re.sub(pat, repl, r1)
            if cand != r1 and os.path.exists(cand):
                r2 = cand
                break
        if r2:
            samples.append((sname, r1, r2))
        else:
            print(f"[WARN] no R2 for {r1}; skipping")
    return samples


# --------------------------------------------------------------------------- #
def route(cfg):
    pipeline_start = time.time()
    data_dir = cfg.data_dir
    out_dir = cfg.out_dir
    os.makedirs(out_dir, exist_ok=True)

    is_allo = (cfg.ploidy == "allo")
    # [R1] pre-v4 auto behaviour (idxstats, no Theta/EM), kept only for ablations
    legacy_auto = bool(getattr(cfg, "auto_idxstats", False))
    if legacy_auto and is_allo:
        print("[FATAL] --auto-idxstats is only meaningful with --ploidy auto")
        sys.exit(1)
    aligner = getattr(cfg, "aligner", "bowtie2")

    from lineagequant import __version__
    print(f"[{get_now()}] [INIT] LineageQuant {__version__} | ploidy={cfg.ploidy} "
          f"({'allopolyploid' if is_allo else 'autopolyploid'}) | aligner={aligner} | "
          f"threads={cfg.threads}")
    if not is_allo:
        if legacy_auto:
            print(f"[{get_now()}] [WARN] --auto-idxstats: mapping -> idxstats -> TPM, "
                  f"no Theta and no EM (featureCounts -M -O equivalent). Ablation only.")
        else:
            print(f"[{get_now()}] [INFO] auto route: mapping -> Theta -> EM -> TPM; "
                  f"no lineage axis, Step 02 skipped.")

    if aligner == "bowtie2":
        tools = ["samtools", "bowtie2", "bowtie2-build"]
    elif aligner == "hisat2":
        tools = ["samtools", "hisat2", "hisat2-build"]
    elif aligner == "sa":
        tools = [getattr(cfg, "salmon", "salmon")]
    else:                                      # lqmap: built into lqcore
        tools = []
    mods = ["pysam", "numpy"]
    if is_allo and cfg.phased_fasta:
        tools.append("mafft")
        mods.append("Bio")
    if is_allo and cfg.ancestors:
        tools += ["kmc", "kmc_tools", "gffread"]
    if not legacy_auto:
        mods.append("Bio")            # global_theta needs biopython
    check_deps(sorted(set(tools)), sorted(set(mods)))
    use_cpp = getattr(cfg, "engine", "python") == "cpp"
    if use_cpp:
        lq = cfg.lqcore
        if not (os.path.isfile(lq) and os.access(lq, os.X_OK)) and shutil.which(lq) is None:
            print(f"[FATAL] --engine cpp: lqcore not found ({lq}); build lq_cpp (bash lq_cpp/hpc/build.sh) "
                  f"and pass --lqcore /path/to/lq_cpp/bin/lqcore or set $LQCORE")
            sys.exit(1)

    # ---- required inputs ----
    if not (cfg.clusters and os.path.exists(cfg.clusters)):
        print("[FATAL] --clusters required (standardize with "
              "lineagequant.cluster.normalize_clusters first)")
        sys.exit(1)
    if not (cfg.ancestry and os.path.exists(cfg.ancestry)):
        print("[FATAL] --ancestry required")
        sys.exit(1)
    if not (cfg.orphan and os.path.exists(cfg.orphan)):
        print("[FATAL] --orphan required")
        sys.exit(1)

    # ---- transcriptome index for the chosen aligner (shared by all samples) ----
    if aligner == "bowtie2":
        index_prefix = os.path.join(data_dir, "transcript_index_bt2")
        if not (os.path.exists(index_prefix + ".1.bt2") or
                os.path.exists(index_prefix + ".1.bt2l")):
            run_step(["bowtie2-build", "--threads", str(cfg.threads), cfg.fasta, index_prefix],
                     "bowtie2 build index", index_prefix + ".done")
    elif aligner == "hisat2":
        index_prefix = os.path.join(data_dir, "transcript_index")
        if not (os.path.exists(index_prefix + ".1.ht2") or
                os.path.exists(index_prefix + ".1.ht2l")):
            run_step(["hisat2-build", "-p", str(cfg.threads), cfg.fasta, index_prefix],
                     "HISAT2 build index", index_prefix + ".done")
    elif aligner == "sa":
        # --engine cpp --aligner sa: Salmon selective alignment proposes candidates, lqcore verifies them ungapped
        index_prefix = os.path.join(data_dir, "salmon_kd_index")
        run_step([getattr(cfg, "salmon", "salmon"), "index", "-t", cfg.fasta, "-i", index_prefix,
                  "-k", "31", "-p", str(cfg.threads), "--keepDuplicates"],
                 "salmon build index (--keepDuplicates)", index_prefix + ".done")
    else:
        # --engine cpp --aligner lqmap: lqcore builds its k-mer position index itself (cached below)
        index_prefix = os.path.join(data_dir, "lqmap.k21.bin")

    samples = discover_samples(data_dir)
    if not samples:
        print(f"[FATAL] no FASTQ sample pairs found in {data_dir}")
        sys.exit(1)
    print(f"[{get_now()}] [INFO] {len(samples)} sample(s): "
          f"{', '.join(s[0] for s in samples)}")

    # ====================================================================== #
    # [R4] long-read physical prior: used for both ploidies (it anchors copy abundance)
    # ====================================================================== #
    lr_prior = None
    if cfg.longread_bam and not legacy_auto:
        # long-read molecules per allele copy (each molecule once: secondary records and extra subread passes are
        # not counted, isoforms are summed); the BAM is read sequentially, any sort order, no index needed
        lr_bam = cfg.longread_bam
        lr_prior = os.path.join(data_dir, "LongRead_Counts.tsv")
        lr_extra = ["--legacy-out", os.path.join(data_dir, "LongRead_Counts.legacy_idxstats.tsv"),
                    "--summary", os.path.join(data_dir, "LongRead_Counts.summary.tsv")]
        if use_cpp:
            lcmd0 = [cfg.lqcore, "lrcount", "--bam", lr_bam, "-o", lr_prior,
                     "--threads", str(min(8, cfg.threads))] + lr_extra
        else:
            lcmd0 = [sys.executable, _mod("longread", "extract_prior.py"), "-b", lr_bam, "-o", lr_prior,
                     "--threads", str(min(8, cfg.threads))] + lr_extra
        run_step(lcmd0, "Step 00b - long-read priors", lr_prior + ".done")

    # ====================================================================== #
    # ALLO ONLY: SNP DB
    # ====================================================================== #
    snp_db_to_use = None
    snp_mode = None
    if is_allo and cfg.vcf:
        snp_mode = "population"
        snp_db = os.path.join(data_dir, "SNP_DB_population.pkl")
        cmd = [sys.executable, _mod("lineage_snp", "build_snp_vcf.py"),
               "--vcf", cfg.vcf, "--gff3", cfg.fasta_gff, "--out", snp_db,
               "--sample-lists"] + cfg.sample_lists
        run_step(cmd, "Step 1a - population SNP DB", snp_db + ".done")
        snp_db_to_use = snp_db
    elif is_allo and cfg.phased_fasta:
        snp_mode = "msa"
        snp_db = os.path.join(data_dir, "SNP_DB_msa.pkl")
        bcmd = [sys.executable, _mod("lineage_snp", "build_snp_msa.py"),
                "--fasta", cfg.phased_fasta, "--cluster", cfg.clusters,
                "--ancestry", cfg.ancestry, "--out", snp_db,
                "--threads", str(max(1, cfg.threads // 2))]
        if use_cpp:      # the text copy for lqcore is written in the same pass (no second load of the pickle)
            bcmd += ["--lq-out", re.sub(r"\.pkl$", "", snp_db) + ".lq.tsv"]
        built_now = not os.path.exists(snp_db + ".done")
        run_step(bcmd, "Step 1a - MSA SNP DB", snp_db + ".done")
        if use_cpp and built_now:
            lq_txt = re.sub(r"\.pkl$", "", snp_db) + ".lq.tsv"
            if os.path.exists(lq_txt):
                with open(lq_txt + ".done", "w") as f:
                    f.write("DONE")
        # Step 1b (long-read calibration of the MSA AFs) is not run: the MSA route scores reads against the
        # expected base of each lineage (entry["bases"]), which calibration never changes, so it had no effect
        # on any result; it only cost a pass over the long reads.
        if cfg.longread_bam and not cfg.no_calibrate_snp:
            print(f"[{get_now()}] [INFO] MSA SNP route: long-read AF calibration skipped "
                  f"(the route uses the expected bases only; calibration never changed a result)")
        snp_db_to_use = snp_db

    # ====================================================================== #
    # ALLO ONLY: lineage k-mer markers
    # ====================================================================== #
    marker_paths = {}
    if is_allo and (cfg.ancestors or cfg.prebuilt_markers):
        if cfg.prebuilt_markers:
            for name, mfa in cfg.prebuilt_markers.items():
                if not (os.path.exists(mfa) and os.path.getsize(mfa)):
                    raise FileNotFoundError(
                        f"--markers {name}: file not found or empty: {mfa}")
                marker_paths[name] = mfa
            print(f"[{get_now()}] [INFO] using {len(marker_paths)} pre-built "
                  f"lineage marker(s); skipping KMC build")
        else:
            marker_dir = os.path.join(data_dir, "lineage_markers")
            marker_done = os.path.join(marker_dir, "markers.done")
            if not os.path.exists(marker_done):
                anc_tx_dir = os.path.join(data_dir, "ancestor_transcripts")
                os.makedirs(anc_tx_dir, exist_ok=True)
                from lineagequant.utils.genome_pair import GenomeAnnotationPair
                lineage_fastas = {}
                for name, anc_samples in cfg.ancestors.items():
                    merged_txfa = os.path.join(anc_tx_dir, f"{name}_transcript.fa")
                    if not (os.path.exists(merged_txfa) and os.path.getsize(merged_txfa)):
                        part_fastas = []
                        for si, (afa, agff) in enumerate(anc_samples):
                            tag = f"{name}_s{si}"
                            part = os.path.join(anc_tx_dir, f"_{tag}_transcript.fa")
                            if not (os.path.exists(part) and os.path.getsize(part)):
                                pair = GenomeAnnotationPair(afa, agff, anc_tx_dir,
                                                            name=tag).prepare()
                                pair.extract_transcripts(part, id_prefix=name)
                            part_fastas.append(part)
                        with open(merged_txfa, "w") as out:
                            for part in part_fastas:
                                with open(part) as fin:
                                    shutil.copyfileobj(fin, out)
                        n = sum(1 for l in open(merged_txfa) if l.startswith(">"))
                        print(f"[{get_now()}] [INFO] lineage {name}: "
                              f"{len(anc_samples)} progenitor sample(s) -> {n} transcripts")
                    lineage_fastas[name] = merged_txfa
                from lineagequant.lineage_kmer.build_markers import build_markers
                marker_paths = build_markers(lineage_fastas, marker_dir,
                                             k=cfg.lineage_kmer, threads=cfg.threads)
                with open(marker_done, "w") as f:
                    f.write("DONE")
            else:
                for name in cfg.ancestors:
                    p = os.path.join(marker_dir, f"{name}.fasta")
                    if os.path.exists(p):
                        marker_paths[name] = p

    # ====================================================================== #
    # C++ engine: SNP DB as text (once, shared by all samples)
    # ====================================================================== #
    snp_db_lq = None
    if use_cpp and snp_db_to_use:
        snp_db_lq = re.sub(r"\.pkl$", "", snp_db_to_use) + ".lq.tsv"
        run_step([sys.executable, _mod("lineage_snp", "export_snp_db.py"),
                  "--snp_db", snp_db_to_use, "--mode", snp_mode, "--out", snp_db_lq],
                 "Step 1c - SNP DB for lqcore", snp_db_lq + ".done")

    # ====================================================================== #
    # PER-SAMPLE
    # ====================================================================== #
    for sname, r1, r2 in samples:
        sout = os.path.join(out_dir, sname)
        os.makedirs(sout, exist_ok=True)
        print(f"\n[{get_now()}] [SAMPLE] {sname}")

        r1_link = os.path.join(sout, f"{sname}.R1.fq.gz")
        r2_link = os.path.join(sout, f"{sname}.R2.fq.gz")
        for lp, rp in [(r1_link, r1), (r2_link, r2)]:
            if not os.path.exists(lp):
                os.symlink(os.path.abspath(rp), lp)

        bam = os.path.join(sout, f"{sname}.sorted.bam")
        # [R7] never mix aligners inside one output directory
        aln_tag = os.path.join(sout, "aligner.txt")
        prev = open(aln_tag).read().strip() if os.path.exists(aln_tag) else (
            "hisat2" if os.path.exists(os.path.join(sout, "00.done")) else None)
        if prev and prev != aligner:
            print(f"[FATAL] {sname}: this output directory was aligned with {prev}, "
                  f"but --aligner {aligner} was requested. Use a new --out-dir "
                  f"(or --aligner {prev} to resume).")
            sys.exit(1)
        if not use_cpp:
            run_step(["bash", _mod("core", "mapping.sh"),
                      r1_link, r2_link, index_prefix, cfg.fasta_gff, sout, str(cfg.threads),
                      aligner, str(cfg.max_alignments)],
                     f"Step 00 - mapping, {aligner} ({sname})", os.path.join(sout, "00.done"))
        if not os.path.exists(aln_tag):
            with open(aln_tag, "w") as f:
                f.write(aligner + "\n")

        # [R1] --auto-idxstats: stop here, Step 04 counts idxstats
        if legacy_auto:
            continue

        # ---- Step 01 Theta (both ploidies) ----
        theta = os.path.join(sout, "Global_Theta_Matrix.tsv")
        tcmd = [sys.executable, _mod("theta", "global_theta.py"),
                "-c", cfg.clusters, "-f", cfg.fasta, "-r1", r1, "-r2", r2,
                "-k", str(cfg.kmer), "-t", str(cfg.theta_threads), "-o", sout,
                "--shrinkage", str(cfg.shrinkage)]
        # [R5] the k-mer index is shared by all samples
        tcmd += ["--index_dir", data_dir]
        if lr_prior:
            tcmd += ["--longread_prior", lr_prior, "--min_lr_depth", str(cfg.min_lr_depth)]
        if cfg.adaptive_shrinkage:
            tcmd.append("--adaptive_shrinkage")
        if cfg.no_paired_bonus:
            tcmd.append("--no_paired_bonus")
        if use_cpp:
            # lqcore theta: same options and output as global_theta.py (byte-identical), all threads, own index cache;
            # --py-sum reproduces this interpreter's built-in sum() (compensated from Python 3.12 on)
            tcmd = [cfg.lqcore, "theta"] + tcmd[2:]
            tcmd[tcmd.index("-t") + 1] = str(cfg.threads)
            tcmd += ["--index-cache", os.path.join(data_dir, f"lqcore_theta.k{cfg.kmer}.bin"),
                     "--py-sum", "neumaier" if sys.version_info >= (3, 12) else "naive"]
        run_step(tcmd, f"Step 01 - Theta ({sname})", theta + ".done")

        # ---- --engine cpp: Steps 00, 02, 02b and 03 in one pass (lqcore run) ----
        if use_cpp:
            em = os.path.join(sout, "EM_Final_Hybrid.tsv" if cfg.longread_bam
                              else "EM_Final_Expression.tsv")
            if aligner == "lqmap":
                aln_args = ["--aligner", "lqmap", "--map-cache", index_prefix]
            elif aligner == "sa":
                aln_args = ["--aligner", "sa", "--salmon-index", index_prefix,
                            "--salmon", getattr(cfg, "salmon", "salmon"),
                            "--bt2-log", os.path.join(sout, f"{sname}.salmon.log")]
            else:
                aln_args = ["-x", index_prefix, "-k", str(cfg.max_alignments),
                            "--bt2-log", os.path.join(sout, f"{sname}.bowtie2.log")]
            lcmd = [cfg.lqcore, "run", "--r1", r1, "--r2", r2] + aln_args + [
                    "--threads", str(cfg.threads),
                    "-a", cfg.ancestry, "-t", theta, "-o", em,
                    "--idxstats-out", os.path.join(sout, f"{sname}.idxstats.tsv"),
                    "--nm-weight", str(getattr(cfg, "nm_weight", "auto"))]
            if getattr(cfg, "em_pool_identical", True):
                lcmd += ["--fasta", cfg.fasta]
            elif aligner in ("lqmap", "sa"):
                lcmd += ["--fasta", cfg.fasta, "--no-pool"]       # the aligner needs the sequences
            if getattr(cfg, "em_accel", "none") != "none":
                lcmd += ["--em-accel", cfg.em_accel]
            for opt, attr in (("--theta-prior-lr", "theta_prior_lr"),
                              ("--theta-prior-kmer", "theta_prior_kmer"),
                              ("--theta-prior-cap", "theta_prior_cap"),
                              ("--eff-len", "eff_len"),
                              ("--p-weight-specific", "p_weight_specific"),
                              ("--p-weight-conflict", "p_weight_conflict"),
                              ("--p-weight-none", "p_weight_none"),
                              ("--p-alpha", "p_alpha"), ("--p-floor", "p_floor")):
                v = getattr(cfg, attr, None)
                if v is not None:
                    lcmd += [opt, str(v)]
            if is_allo:
                if cfg.target_ids:
                    lcmd += ["--target-ids", cfg.target_ids]
                for n, p in marker_paths.items():
                    lcmd += ["--markers", f"{n}:{p}"]
                if marker_paths:
                    lcmd += ["--kmer-len", str(cfg.lineage_kmer), "--marker-cache",
                             os.path.join(data_dir, f"lqcore_markers.k{cfg.lineage_kmer}.bin")]
                if snp_db_lq:
                    lcmd += ["--snp-db", snp_db_lq, "--error-rate", str(cfg.error_rate)]
                    if cfg.no_chimeric:
                        lcmd.append("--no-chimeric")
                if not marker_paths and not snp_db_lq:
                    print(f"[{get_now()}] [WARN] allo but no marker/SNP signal -> "
                          f"W = M . theta (same as auto)")
            run_step(lcmd, f"Step 00-03 - lqcore run: {aligner} + lineage + EM ({sname})",
                     em + ".lqcore.done")
            continue

        # ---- Step 02 lineage assignment (allo only) ----
        pval_for_em = None
        if is_allo:
            kmer_pv = snp_pv = None
            if marker_paths:
                kmer_pv = os.path.join(sout, "P_values_kmer.tsv")
                kcmd = [sys.executable, _mod("lineage_kmer", "snipe_kmer.py"),
                        "-i", cfg.target_ids, "-b", bam, "--kmers"]
                kcmd += [f"{n}:{p}" for n, p in marker_paths.items()]
                kcmd += ["-k", str(cfg.lineage_kmer), "-t", str(cfg.threads), "-o", kmer_pv]
                run_step(kcmd, f"Step 02 - k-mer lineage ({sname})", kmer_pv + ".done")
            if snp_db_to_use:
                snp_pv = os.path.join(sout, "P_values_snp.tsv")
                scmd = [sys.executable, _mod("lineage_snp", "snipe_snp.py"),
                        "-i", cfg.target_ids, "-b", bam, "--snp_db", snp_db_to_use,
                        "--mode", snp_mode, "--error_rate", str(cfg.error_rate),
                        "-t", str(cfg.threads), "-o", snp_pv]
                if cfg.no_chimeric:
                    scmd.append("--no_chimeric")
                run_step(scmd, f"Step 02 - SNP lineage ({sname})", snp_pv + ".done")

            if kmer_pv and snp_pv:
                merged = os.path.join(sout, "P_values_merged.tsv")
                run_step([sys.executable, _mod("merge", "merge_poe.py"),
                          "--kmer", kmer_pv, "--snp", snp_pv, "--out", merged],
                         f"Step 02b - PoE fusion ({sname})", merged + ".done")
                pval_for_em = merged
            else:
                pval_for_em = kmer_pv or snp_pv

            # [R3] no lineage signal -> no -p; EM uses W = M . theta
            if pval_for_em is None:
                print(f"[{get_now()}] [WARN] allo but no marker/SNP signal -> "
                      f"W = M . theta (same as auto)")

        # ---- Step 03 EM (both ploidies; compatibility always from the BAM) ----
        em = os.path.join(sout, "EM_Final_Hybrid.tsv" if cfg.longread_bam
                          else "EM_Final_Expression.tsv")
        emcmd = [sys.executable, _mod("core", "em_dynamic.py"),
                 "-a", cfg.ancestry, "-t", theta, "--bam", bam,
                 "--threads", str(cfg.threads),
                 "--epsilon", str(cfg.em_epsilon),
                 "--max_iter", str(cfg.em_max_iter),
                 "--min_iter", str(cfg.em_min_iter),
                 "--tol", str(cfg.em_tol), "--ll_tol", str(cfg.em_ll_tol),
                 "-o", em]
        for opt, attr in (("--theta-prior-lr", "theta_prior_lr"),
                          ("--theta-prior-kmer", "theta_prior_kmer"),
                          ("--theta-prior-cap", "theta_prior_cap"),
                          ("--eff-len", "eff_len"),
                          ("--cross-cluster", "cross_cluster")):
            v = getattr(cfg, attr, None)
            if v is not None:
                emcmd += [opt, str(v)]
        if pval_for_em:
            emcmd += ["-p", pval_for_em]
            for opt, attr in (("--p-weight-specific", "p_weight_specific"),
                              ("--p-weight-conflict", "p_weight_conflict"),
                              ("--p-weight-none", "p_weight_none"),
                              ("--p-alpha", "p_alpha"),
                              ("--p-floor", "p_floor")):
                v = getattr(cfg, attr, None)
                if v is not None:
                    emcmd += [opt, str(v)]
        em_done = em + ".done"
        if getattr(cfg, "em_engine", "global") == "global":
            # [R6] global engine: one EM over all genes; [R7] mismatch-count weighting
            emcmd = [sys.executable, _mod("core", "em_global.py"),
                     "-a", cfg.ancestry, "-t", theta, "--bam", bam,
                     "--threads", str(cfg.threads), "-o", em,
                     "--nm-weight", str(getattr(cfg, "nm_weight", "auto"))]
            if getattr(cfg, "em_pool_identical", True):
                emcmd += ["--fasta", cfg.fasta]
            for opt, attr in (("--theta-prior-lr", "theta_prior_lr"),
                              ("--theta-prior-kmer", "theta_prior_kmer"),
                              ("--theta-prior-cap", "theta_prior_cap"),
                              ("--eff-len", "eff_len")):
                v = getattr(cfg, attr, None)
                if v is not None:
                    emcmd += [opt, str(v)]
            if pval_for_em:
                emcmd += ["-p", pval_for_em]
                for opt, attr in (("--p-weight-specific", "p_weight_specific"),
                                  ("--p-weight-conflict", "p_weight_conflict"),
                                  ("--p-weight-none", "p_weight_none"),
                                  ("--p-alpha", "p_alpha"), ("--p-floor", "p_floor")):
                    v = getattr(cfg, attr, None)
                    if v is not None:
                        emcmd += [opt, str(v)]
            em_done = em + ".global.done"
        run_step(emcmd, f"Step 03 - EM, {getattr(cfg, 'em_engine', 'global')} engine ({sname})", em_done)

    # ====================================================================== #
    # Step 04: aggregation
    #   allo / auto  -> EM_Allocated_Reads (global engine: every gene; legacy engine:
    #                   multi-copy genes, orphans from idxstats)
    #   idxstats     -> --auto-idxstats ablation (every gene from idxstats)
    # ====================================================================== #
    agg_mode = "idxstats" if legacy_auto else cfg.ploidy
    final = os.path.join(out_dir, "LineageQuant_Final_Merged_TPM.tsv")
    run_step([sys.executable, _mod("core", "aggregator.py"),
              "-r", out_dir, "-o", cfg.orphan, "--out", final,
              "--mode", agg_mode],
             "Step 04 - aggregate TPM", final + ".done")

    print(f"\n[{get_now()}] [DONE] total {fmt_dur(time.time() - pipeline_start)}")
    print(f"  Final TPM: {final}")
    return final
