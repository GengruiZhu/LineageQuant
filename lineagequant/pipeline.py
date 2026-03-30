#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LineageQuant Master Pipeline
Orchestrates the full lineage-resolved expression quantification workflow.

  Stage 0  HISAT2 alignment + featureCounts
  Stage 1  K-mer Theta initialisation
  Stage 2  BAM-level per-read lineage P-value assignment
  Stage 3  EM deconvolution
  Stage 4  Global TPM matrix aggregation

Supports any number (≥2) of ancestral lineages.
"""

import os
import sys
import glob
import subprocess
import argparse
import time
from collections import OrderedDict

import lineagequant


def get_now():
    return time.strftime("%Y-%m-%d %H:%M:%S")


def run_cmd(cmd_list, step_name):
    """Execute a subprocess with logging and hard-stop on failure."""
    print(f"\n[{get_now()}] [RUNNING] {step_name}")
    print(f"  -> Command: {' '.join(str(x) for x in cmd_list)}")
    try:
        subprocess.run(cmd_list, check=True)
        print(f"[{get_now()}] [SUCCESS] {step_name} completed.")
    except subprocess.CalledProcessError as e:
        print(f"\n[{get_now()}] [FATAL] {step_name} failed (exit {e.returncode}). "
              "Pipeline aborted.")
        sys.exit(1)


def check_files(files_dict):
    """Pre-flight: verify all required files exist before launching."""
    print(f"[{get_now()}] [INFO] Pre-flight file check...")
    missing = [f"{name}: {path}" for name, path in files_dict.items()
               if not os.path.exists(path)]
    if missing:
        print(f"[{get_now()}] [FATAL] Missing required files:")
        for m in missing:
            print(f"  - {m}")
        sys.exit(1)
    print("  -> All required files present.")


def parse_ancestors(ancestor_list):
    """
    Parse '--ancestors NAME:PATH ...' into an OrderedDict.
    Validates that at least 2 lineages are provided.
    """
    result = OrderedDict()
    for item in ancestor_list:
        if ':' not in item:
            print(f"[ERROR] Invalid --ancestors entry '{item}'. Expected NAME:PATH.")
            sys.exit(1)
        name, path = item.split(':', 1)
        result[name.strip()] = path.strip()
    if len(result) < 2:
        print(f"[ERROR] --ancestors requires ≥2 lineages, got {len(result)}.")
        sys.exit(1)
    return result


def main():
    parser = argparse.ArgumentParser(
        description=(
            "LineageQuant v{ver} – Allele-resolved RNA-seq quantification "
            "for polyploid genomes.\n"
            "Supports any number (≥2) of ancestral lineages."
        ).format(ver=lineagequant.__version__),
        formatter_class=argparse.RawTextHelpFormatter
    )

    # ── Reference files (all required) ───────────────────────────────────────
    ref = parser.add_argument_group("Reference files (required)")
    ref.add_argument(
        "--gff", required=True,
        help="Whole-genome GFF3 annotation file"
    )
    ref.add_argument(
        "--fasta", required=True,
        help="Whole-genome reference FASTA file"
    )
    ref.add_argument(
        "--clusters", required=True,
        help="Allele cluster TSV: ClusterID<TAB>Allele1,Allele2,..."
    )
    ref.add_argument(
        "--gene-ids", required=True, dest="gene_ids",
        help="Target gene ID list (one ID per line)"
    )
    ref.add_argument(
        "--ancestry", required=True,
        help="Gene ancestry table: GeneID<TAB>LineageName"
    )
    ref.add_argument(
        "--orphan-genes", required=True, dest="orphan_genes",
        help="Single-copy / orphan gene ID list for featureCounts quantification"
    )
    ref.add_argument(
        "--ancestors", required=True, nargs='+', metavar="NAME:PATH",
        help=(
            "Ancestral lineage k-mer FASTA files (≥2 required).\n"
            "Format : NAME:PATH\n"
            "Example: --ancestors SubA:subA_k15.fa SubB:subB_k15.fa SubC:subC_k15.fa"
        )
    )

    # ── Input / Output ────────────────────────────────────────────────────────
    io = parser.add_argument_group("Input / Output")
    io.add_argument(
        "-d", "--data-dir", dest="data_dir", default="./data",
        help="Directory containing paired FASTQ files (*.R1.fq.gz). Default: ./data"
    )
    io.add_argument(
        "-o", "--out-dir", dest="out_dir", default="./results",
        help="Root output directory. Default: ./results"
    )
    io.add_argument(
        "--index-prefix", dest="index_prefix", default=None,
        help="HISAT2 genome index prefix. Default: <data_dir>/genome_index"
    )

    # ── Runtime parameters ────────────────────────────────────────────────────
    run = parser.add_argument_group("Runtime parameters")
    run.add_argument(
        "-t", "--threads", type=int, default=16,
        help="CPU threads. Default: 16"
    )
    run.add_argument(
        "--kmer", type=int, default=31,
        help="K-mer size for the Theta engine (Module 01). Default: 31"
    )
    run.add_argument(
        "--lineage-kmer", type=int, default=15, dest="lineage_kmer",
        help="K-mer size for lineage k-mer files (Module 02). Default: 15"
    )
    run.add_argument(
        "--epsilon", type=float, default=0.01,
        help="Laplace smoothing pseudo-count for EM. Default: 0.01"
    )
    run.add_argument(
        "--max-iter", type=int, default=200, dest="max_iter",
        help="Maximum EM iterations. Default: 200"
    )
    run.add_argument(
        "--tol", type=float, default=1e-5,
        help="EM convergence tolerance. Default: 1e-5"
    )

    args = parser.parse_args()

    # ── Resolve absolute paths ────────────────────────────────────────────────
    data_dir     = os.path.abspath(args.data_dir)
    out_dir      = os.path.abspath(args.out_dir)
    index_prefix = (
        os.path.abspath(args.index_prefix)
        if args.index_prefix
        else os.path.join(data_dir, "genome_index")
    )
    os.makedirs(out_dir, exist_ok=True)

    # ── Locate module scripts inside the installed package ────────────────────
    pkg_dir    = os.path.dirname(os.path.abspath(__file__))
    modules_dir = os.path.join(pkg_dir, "modules")
    map_sh     = os.path.join(pkg_dir, "scripts", "run_mapping.sh")

    # ── Parse and validate ancestors ──────────────────────────────────────────
    ancestors  = parse_ancestors(args.ancestors)
    kmers_args = [f"{name}:{path}" for name, path in ancestors.items()]

    # ── Pre-flight check ──────────────────────────────────────────────────────
    required = {
        "GFF3"          : args.gff,
        "FASTA"         : args.fasta,
        "Clusters"      : args.clusters,
        "Gene IDs"      : args.gene_ids,
        "Ancestry"      : args.ancestry,
        "Orphan genes"  : args.orphan_genes,
        "Mapping script": map_sh,
    }
    for name, path in ancestors.items():
        required[f"Ancestor k-mer [{name}]"] = path
    check_files(required)

    # ── Build HISAT2 index if absent ──────────────────────────────────────────
    if not (os.path.exists(f"{index_prefix}.1.ht2")
            or os.path.exists(f"{index_prefix}.1.ht2l")):
        print(f"[{get_now()}] [INFO] No HISAT2 index at {index_prefix}. Building...")
        run_cmd(
            ["hisat2-build", "-p", str(args.threads), args.fasta, index_prefix],
            "Build HISAT2 genome index"
        )
    else:
        print(f"[{get_now()}] [SKIP] HISAT2 index detected.")

    # ── Scan FASTQ files ──────────────────────────────────────────────────────
    r1_files = sorted(glob.glob(os.path.join(data_dir, "*.R1.fq.gz")))
    if not r1_files:
        print(f"[{get_now()}] [ERROR] No *.R1.fq.gz files found in {data_dir}")
        sys.exit(1)

    print("\n" + "=" * 75)
    print(f"  LineageQuant  v{lineagequant.__version__}")
    print(f"  Samples   : {len(r1_files)}")
    print(f"  Lineages  : {list(ancestors.keys())}  ({len(ancestors)} total)")
    print(f"  K-mer (θ) : {args.kmer}    K-mer (lineage): {args.lineage_kmer}")
    print(f"  Threads   : {args.threads}")
    print(f"  Output    : {out_dir}")
    print("=" * 75)

    # ── Per-sample loop ───────────────────────────────────────────────────────
    for r1 in r1_files:
        sample_name = os.path.basename(r1).replace(".R1.fq.gz", "")
        r2          = os.path.join(data_dir, f"{sample_name}.R2.fq.gz")

        if not os.path.exists(r2):
            print(f"[{get_now()}] [WARNING] R2 not found for {sample_name}, skipping.")
            continue

        sample_out = os.path.join(out_dir, sample_name)
        os.makedirs(sample_out, exist_ok=True)

        print(f"\n  ── Sample: {sample_name} {'─' * (60 - len(sample_name))}")

        bam_file   = os.path.join(sample_out, f"{sample_name}.sorted.bam")
        fc_file    = os.path.join(sample_out, f"{sample_name}_featureCounts.txt")
        theta_file = os.path.join(sample_out, "Theta_Matrix.tsv")
        pval_file  = os.path.join(sample_out, "P_values.tsv")
        em_file    = os.path.join(sample_out, "EM_Final_Expression.tsv")

        # Stage 0: Alignment + featureCounts
        if not (os.path.exists(bam_file) and os.path.exists(fc_file)):
            run_cmd(
                ["bash", map_sh,
                 r1, r2, index_prefix, args.gff, sample_out, str(args.threads)],
                f"Stage 0: Alignment + featureCounts [{sample_name}]"
            )
        else:
            print(f"[{get_now()}] [SKIP] Stage 0 – BAM and featureCounts exist.")

        # Stage 1: K-mer Theta
        if not os.path.exists(theta_file):
            run_cmd(
                ["python", os.path.join(modules_dir, "theta.py"),
                 "-c", args.clusters,
                 "-g", args.gff,
                 "-f", args.fasta,
                 "-r1", r1, "-r2", r2,
                 "-k", str(args.kmer),
                 "-o", sample_out],
                f"Stage 1: K-mer Theta [{sample_name}]"
            )
        else:
            print(f"[{get_now()}] [SKIP] Stage 1 – Theta matrix exists.")

        # Stage 2: BAM P-value assignment
        if not os.path.exists(pval_file):
            run_cmd(
                ["python", os.path.join(modules_dir, "lineage_pval.py"),
                 "-i", args.gene_ids,
                 "-g", args.gff,
                 "-b", bam_file,
                 "-k", str(args.lineage_kmer),
                 "--kmers"] + kmers_args + ["-o", pval_file],
                f"Stage 2: BAM lineage P-values [{sample_name}]"
            )
        else:
            print(f"[{get_now()}] [SKIP] Stage 2 – P-value table exists.")

        # Stage 3: EM deconvolution
        if not os.path.exists(em_file):
            run_cmd(
                ["python", os.path.join(modules_dir, "em.py"),
                 "-a", args.ancestry,
                 "-t", theta_file,
                 "-p", pval_file,
                 "-o", em_file,
                 "--epsilon", str(args.epsilon),
                 "--max-iter", str(args.max_iter),
                 "--tol",      str(args.tol)],
                f"Stage 3: EM deconvolution [{sample_name}]"
            )
        else:
            print(f"[{get_now()}] [SKIP] Stage 3 – EM output exists.")

    # ── Stage 4: Aggregate TPM matrix ─────────────────────────────────────────
    print(f"\n  ── Stage 4: TPM matrix aggregation {'─' * 38}")
    final_merged = os.path.join(out_dir, "LineageQuant_Final_Merged_TPM.tsv")
    if not os.path.exists(final_merged):
        run_cmd(
            ["python", os.path.join(modules_dir, "aggregator.py"),
             "-g", args.gff,
             "-r", out_dir,
             "-o", args.orphan_genes,
             "--out", "LineageQuant_Final_Merged_TPM.tsv"],
            "Stage 4: TPM matrix aggregation"
        )
    else:
        print(f"[{get_now()}] [SKIP] Stage 4 – Final merged matrix exists.")

    print("\n" + "=" * 75)
    print(f"[{get_now()}] LineageQuant pipeline completed successfully.")
    print(f"  Final matrix: {final_merged}")
    print("=" * 75)


if __name__ == "__main__":
    main()
