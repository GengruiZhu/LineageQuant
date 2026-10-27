#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
cli.py -- LineageQuant command-line entry point.

Default invocation runs the auto-routing pipeline. The user declares the ploidy
type and supplies whatever inputs they have; the router decides the rest.

  lineagequant run \
    --fasta mRNA.fasta --gff genome.gff3 \
    --clusters clusters.tsv --ancestry ancestry.tsv --orphan orphans.tsv \
    --target-ids gene_ids.txt \
    --ploidy allo --ancestor A:ancA.fa,ancA.gff D:ancD.fa,ancD.gff \
    --data-dir data --out-dir results -t 32

Helper subcommands:
  lineagequant normalize-clusters ...   # standardize cluster format
  lineagequant build-markers ...        # build lineage k-mer markers only
"""

import sys
import argparse

from lineagequant.config import Config
from lineagequant import router


def _add_run_args(p):
    req = p.add_argument_group("required reference files")
    req.add_argument("--fasta", required=True, help="Transcriptome (mRNA) FASTA")
    req.add_argument("--gff", required=True, help="Reference GFF3 (reserved/mapping)")
    req.add_argument("--clusters", required=True,
                     help="Standardized allele clusters "
                          "(ClusterID<TAB>g1,g2,...). Use normalize-clusters first.")
    req.add_argument("--ancestry", required=True, help="GeneID<TAB>Lineage")
    req.add_argument("--orphan", required=True, help="Orphan/single-copy gene list")
    req.add_argument("--target-ids", dest="target_ids", required=True,
                     help="Target gene ID list (one per line)")

    pl = p.add_argument_group("ploidy + lineage inputs")
    pl.add_argument("--ploidy", choices=["auto", "allo"], default="allo",
                    help="auto=autopolyploid (copy-level only); "
                         "allo=allopolyploid (lineage assignment). Default: allo")
    pl.add_argument("--auto-idxstats", dest="auto_idxstats",
                    action="store_true", default=False,
                    help="[ablation only] with --ploidy auto: mapping -> idxstats -> TPM, "
                         "no Theta and no EM. Every reported alignment is counted in full "
                         "(equivalent to featureCounts -M -O). Requires --aligner hisat2.")
    pl.add_argument("--ancestor", action="append", metavar="NAME:fa,gff",
                    help="Ancestor genome for k-mer markers (repeatable; same NAME "
                         "twice = multi-sample). allo only.")
    pl.add_argument("--markers", action="append", metavar="NAME:fasta",
                    help="Pre-built lineage marker FASTA (repeatable). Skips the "
                         "KMC build step and uses these directly. NAME matches "
                         "--ancestry lineage names. allo only.")
    pl.add_argument("--vcf", default=None, help="Population VCF (allo, SNP route)")
    pl.add_argument("--sample-lists", dest="sample_lists", nargs="+", default=None,
                    metavar="NAME:file", help="Per-lineage sample files (with --vcf)")
    pl.add_argument("--phased-fasta", dest="phased_fasta", default=None,
                    help="Separate phased/haplotype-resolved CDS FASTA for the MSA "
                         "SNP route (allo). Only needed if the alleles' CDS are NOT "
                         "already separated in the main --fasta. If your main "
                         "--fasta already contains all allele copies (the usual "
                         "case), use --msa-snp instead of this.")
    pl.add_argument("--msa-snp", dest="msa_snp", action="store_true",
                    help="Run the MSA SNP route using the main --fasta directly "
                         "(when it already holds all allele-copy CDS). Mutually "
                         "exclusive with --vcf and --phased-fasta.")
    pl.add_argument("--longread", default=None,
                    help="Long-read BAM (minimap2 -a -N 50 -> transcriptome, any sort order); "
                         "long-read molecules per copy anchor Theta")

    io = p.add_argument_group("input / output")
    io.add_argument("-d", "--data-dir", dest="data_dir", default="./data")
    io.add_argument("-o", "--out-dir", dest="out_dir", default="./results")

    al = p.add_argument_group("short-read alignment (Step 00)")
    al.add_argument("--aligner", choices=["bowtie2", "hisat2", "lqmap", "sa"], default="bowtie2",
                    help="bowtie2 (default): RSEM settings -- end-to-end, ungapped, "
                         "--mp 1,1 --score-min L,0,-0.1, ambiguous reference bases not "
                         "penalised (--np 0 --n-ceil L,0,0.5), concordant pairs only, up to "
                         "--max-alignments per fragment. Reports every equally good copy and "
                         "near-optimal ones, which the EM weights by mismatch count. "
                         "hisat2: the pre-1.0 aligner (-k 15 --no-mixed --no-discordant); "
                         "faster, but does not report all equally good copies. "
                         "lqmap (--engine cpp, experimental): lqcore's built-in aligner with the same "
                         "rules as bowtie2 (exact k-mer seeds, ungapped mismatch count, concordant pairs, "
                         "candidates within best + 3 mismatches). "
                         "sa (--engine cpp, experimental): Salmon selective alignment proposes the "
                         "candidates; lqcore re-verifies them with the same ungapped rules.")
    al.add_argument("--salmon", default="salmon", help="salmon binary for --aligner sa")
    al.add_argument("--max-alignments", dest="max_alignments", type=int, default=None,
                    help="Maximum alignments reported per fragment "
                         "(default: 200 for bowtie2, 15 for hisat2).")

    rt = p.add_argument_group("runtime / tuning")
    rt.add_argument("-t", "--threads", type=int, default=16)
    rt.add_argument("--engine", choices=["python", "cpp"], default="python",
                    help="python (default): Steps 00-03 as separate programs with a sorted BAM and "
                         "per-read probability tables in between. cpp: lqcore (C++) runs bowtie2 and "
                         "computes candidates, lineage evidence and the global EM in one pass per "
                         "sample, without BAM or intermediate tables (same results; the SNP row of a "
                         "fragment aligned to several transcripts is chosen deterministically).")
    rt.add_argument("--lqcore", default=None,
                    help="Path to the lqcore binary (default: $LQCORE, then lqcore on PATH)")
    rt.add_argument("--em-accel", dest="em_accel", choices=["none", "squarem"], default="none",
                    help="--engine cpp only. none (default): the EM of 1.0.0, identical results. "
                         "squarem: SQUAREM-accelerated EM (monotone); converges to the same maximum "
                         "in far fewer iterations, results differ at the level of the tolerance.")
    rt.add_argument("--kmer", type=int, default=21, help="Theta k-mer length")
    rt.add_argument("--theta-threads", dest="theta_threads", type=int, default=8,
                    help="Worker processes for the Theta FASTQ scan (default 8, capped at "
                         "--threads). Each worker holds the copy-unique k-mer index, so memory "
                         "grows with this number on large transcriptomes.")
    rt.add_argument("--lineage-kmer", dest="lineage_kmer", type=int, default=31,
                    help="Lineage marker k-mer length")
    rt.add_argument("--shrinkage", type=float, default=10.0)
    rt.add_argument("--adaptive-shrinkage", dest="adaptive_shrinkage",
                    action="store_true", default=False)
    rt.add_argument("--no-paired-bonus", dest="no_paired_bonus",
                    action="store_true", default=False)
    rt.add_argument("--min-lr-depth", dest="min_lr_depth", type=int, default=5)
    rt.add_argument("--error-rate", dest="error_rate", type=float, default=0.01,
                    help="Base error rate for the SNP lineage routes")

    cal = p.add_argument_group("SNP calibration (MSA route, long-read)")
    # kept for old command lines: the MSA route scores expected bases, which calibration never changed,
    # so Step 1b is no longer run and these options have no effect
    cal.add_argument("--no-calibrate-snp", dest="no_calibrate_snp",
                     action="store_true", default=False, help=argparse.SUPPRESS)
    cal.add_argument("--calib-min-depth", dest="calib_min_depth", type=int, default=10, help=argparse.SUPPRESS)
    cal.add_argument("--calib-det-weight", dest="calib_det_weight", type=float, default=0.25, help=argparse.SUPPRESS)
    cal.add_argument("--no-chimeric", dest="no_chimeric",
                     action="store_true", default=False)

    em = p.add_argument_group("EM (Step 03)")
    em.add_argument("--em-engine", dest="em_engine",
                    choices=["global", "legacy", "v6", "v5"], default="global",
                    help="global (default): one EM over read equivalence classes of all "
                         "genes -- every fragment carries total weight 1, single-copy genes "
                         "are included, identical-sequence genes are pooled. "
                         "legacy: the pre-1.0 per-cluster EM with idxstats counts for "
                         "orphan genes. (v6/v5 are accepted as aliases.)")
    em.add_argument("--nm-weight", dest="nm_weight", default="auto",
                    help="global engine: likelihood of a candidate alignment multiplied by "
                         "lambda^(extra mismatches), lambda = (p/3)/(1-p). "
                         "auto (default): p estimated from the best alignments; "
                         "off: all candidates equally likely; or a number = lambda.")
    em.add_argument("--no-pool-identical", dest="no_pool_identical", action="store_true",
                    help="global engine: do not pool genes whose transcripts are identical "
                         "in sequence (default: pooled, split evenly after convergence)")
    em.add_argument("--eff-len", dest="eff_len", choices=["off", "len", "eff"], default="eff",
                    help="Effective-length correction in the E-step: weight = theta / "
                         "(L - mean fragment length + 1). Default eff; off = raw read share.")
    em.add_argument("--p-weight-specific", dest="p_weight_specific",
                    type=float, default=1.0,
                    help="Exponent w on lineage probabilities P^w for reads with "
                         "Lineage_Status=Specific* (allo)")
    em.add_argument("--p-weight-conflict", dest="p_weight_conflict",
                    type=float, default=0.0,
                    help="Exponent w for Conflict/Multi-lineage reads (0 = ignore P)")
    em.add_argument("--p-weight-none", dest="p_weight_none",
                    type=float, default=0.0,
                    help="Exponent w for reads without lineage evidence (0 = ignore P)")
    em.add_argument("--p-alpha", dest="p_alpha", type=float, default=0.0,
                    help="Beta shrinkage of P towards uniform")
    em.add_argument("--p-floor", dest="p_floor", type=float, default=1e-3,
                    help="Lower bound on P, so that P=0 cannot erase a read")
    em.add_argument("--theta-prior-lr", dest="theta_prior_lr",
                    type=float, default=1.0,
                    help="LongRead_Anchored clusters: Dirichlet weight of theta0 in the "
                         "M-step = this value x the cluster's long-read count "
                         "(1.0 = one long read counts as one pseudo-read; 0 = off)")
    em.add_argument("--theta-prior-kmer", dest="theta_prior_kmer",
                    type=float, default=0.0,
                    help="Kmer_Inferred clusters: fixed pseudo-count for theta0 (default 0)")
    em.add_argument("--theta-prior-cap", dest="theta_prior_cap",
                    type=float, default=0.0,
                    help="Upper bound on the Dirichlet weight (0 = no cap)")

    lg = p.add_argument_group("legacy EM engine only (--em-engine legacy)")
    lg.add_argument("--em-epsilon", dest="em_epsilon", type=float, default=0.01)
    lg.add_argument("--em-max-iter", dest="em_max_iter", type=int, default=200)
    lg.add_argument("--em-min-iter", dest="em_min_iter", type=int, default=5)
    lg.add_argument("--em-tol", dest="em_tol", type=float, default=1e-4)
    lg.add_argument("--em-ll-tol", dest="em_ll_tol", type=float, default=0.01)
    lg.add_argument("--cross-cluster", dest="cross_cluster", choices=["keep", "drop", "split"],
                    default=None,
                    help="Reads mapping to several clusters: keep (count in each cluster), "
                         "drop, or split (weight 1 shared across clusters)")

    p.add_argument("--dry-run", action="store_true",
                   help="Print the resolved route and exit without running.")


def main(argv=None):
    from lineagequant import __version__
    parser = argparse.ArgumentParser(
        prog="lineagequant",
        description="Haplotype-resolved expression quantification for polyploid "
                    "genomes (auto-routing).")
    parser.add_argument("--version", action="version",
                        version=f"lineagequant {__version__}")
    sub = parser.add_subparsers(dest="command")

    # run (default)
    p_run = sub.add_parser("run", help="Run the auto-routing pipeline.")
    _add_run_args(p_run)

    # normalize-clusters (delegates to module)
    sub.add_parser("normalize-clusters",
                   help="Standardize allele clusters (see "
                        "lineagequant.cluster.normalize_clusters --help).")
    # build-markers (delegates to module)
    sub.add_parser("build-markers",
                   help="Build lineage k-mer markers only (see "
                        "lineagequant.lineage_kmer.build_markers --help).")

    # `lineagequant --version` / `-h` handled by argparse above.
    # Allow `lineagequant <run-args>` without the `run` subcommand -> treat as `run`.
    args, extra = parser.parse_known_args(argv)

    if args.command == "normalize-clusters":
        from lineagequant.cluster import normalize_clusters
        sys.argv = ["normalize-clusters"] + extra
        return normalize_clusters.main()
    if args.command == "build-markers":
        from lineagequant.lineage_kmer import build_markers
        sys.argv = ["build-markers"] + extra
        return build_markers.main()

    if args.command != "run":
        # No subcommand given: reparse everything as `run`
        p2 = argparse.ArgumentParser(prog="lineagequant")
        _add_run_args(p2)
        args = p2.parse_args(argv)
    elif extra:
        p_run.error(f"unrecognized arguments: {' '.join(extra)}")

    try:
        cfg = Config(args)
    except (ValueError, FileNotFoundError) as ex:
        print(f"[FATAL] {ex}", file=sys.stderr)
        return 2
    print(f"[lineagequant {__version__}] resolved route:")
    print(cfg.describe_route())
    if getattr(args, "dry_run", False):
        print("[dry-run] exiting without execution.")
        return 0
    router.route(cfg)
    return 0


if __name__ == "__main__":
    sys.exit(main())
