# -*- coding: utf-8 -*-
"""
core/aggregator.py -- expression matrix aggregator (Step 04).

Collects per-sample EM allocations into one genome-wide matrix with, per sample,
<sample>_RawReads (EM-assigned fragments) and <sample>_TPM (RawReads per kb of the
longest transcript, scaled to one million).

Mode (--mode):
  - allo / auto : read counts from EM_Allocated_Reads. With the global EM engine
                  (default) every gene is in the EM output; with the legacy engine
                  multi-copy genes come from the EM and single-copy (orphan) genes
                  from samtools idxstats. Genes absent from the EM output fall back
                  to idxstats [A6].
  - idxstats    : ablation only (--auto-idxstats). Every gene is counted from
                  idxstats, i.e. every reported alignment counts in full -- the same
                  numbers as featureCounts -M -O on HISAT2 alignments. It runs none of
                  LineageQuant's algorithms.

Type column: Polyploid(EM) (allele copy, from the EM), Orphan(EM) (single-copy gene,
from the global EM), Orphan(BAM_Base) / Polyploid(idxstats) (idxstats counts).

EM file search order: EM_Final_Hybrid.tsv -> EM_Final_Expression_v2.tsv ->
EM_Final_Expression.tsv (covers hybrid/long-read, variant, and k-mer runs).
"""
import os
import glob
import datetime
import argparse
import pandas as pd
import pysam


def get_now():
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def parse_orphan_genes(orphan_file):
    print(f"[INFO] Parsing orphan database: {os.path.basename(orphan_file)}")
    orphans = set()
    with open(orphan_file) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                for g in parts[1].split(","):
                    orphans.add(g.rsplit(".",1)[0])
            elif len(parts) == 1:
                orphans.add(parts[0].rsplit(".",1)[0])
    return orphans


def find_em_file(sample_dir):
    for name in ("EM_Final_Hybrid.tsv", "EM_Final_Expression_v2.tsv",
                 "EM_Final_Expression.tsv"):
        p = os.path.join(sample_dir, name)
        if os.path.exists(p):
            return p
    return None


def main():
    parser = argparse.ArgumentParser(description="TPM aggregator (Step 04)")
    parser.add_argument("-r", "--results_dir", required=True)
    parser.add_argument("-o", "--orphan", required=True)
    parser.add_argument("--out", default="LineageQuant_Final_Merged_TPM.tsv")
    parser.add_argument("--mode", choices=["auto", "allo", "idxstats"], default="allo",
                        help="allo/auto = EM for multi-copy + idxstats for orphans; "
                             "idxstats = ablation: every gene from idxstats (featureCounts -M -O equivalent)")
    args = parser.parse_args()

    use_em = args.mode in ("allo", "auto")
    orphan_set = parse_orphan_genes(args.orphan)

    bam_files = glob.glob(os.path.join(args.results_dir, "**", "*.sorted.bam"),
                          recursive=True)
    # --engine cpp writes no BAM; lqcore writes <sample>.idxstats.tsv (samtools idxstats format)
    idx_files = glob.glob(os.path.join(args.results_dir, "**", "*.idxstats.tsv"),
                          recursive=True)
    if not bam_files and not idx_files:
        print("[ERROR] No sorted.bam (or lqcore .idxstats.tsv) files found in results directories.")
        return
    bam_by_sample = {}
    for bf in bam_files:
        sname = os.path.basename(bf).replace(".sorted.bam", "")
        if sname not in bam_by_sample or bf.count(os.sep) < bam_by_sample[sname].count(os.sep):
            bam_by_sample[sname] = bf
    for bf in idx_files:
        sname = os.path.basename(bf)[:-len(".idxstats.tsv")]
        if sname not in bam_by_sample:
            bam_by_sample[sname] = bf
    sample_names = sorted(bam_by_sample.keys())
    print(f"[INFO] Found {len(sample_names)} samples. Mode={args.mode} "
          f"(EM={'yes' if use_em else 'no'}). Aggregating...")
    if not use_em:
        print("[WARN] mode=idxstats: no EM. Output equals featureCounts -M -O and does "
              "not reflect LineageQuant's quantification.")

    master_counts = {}
    gene_lengths = {}
    gene_src = {}          # [A6] source of each gene's count: em / em_orphan / idxstats
    n_missing_em = 0
    for sample in sample_names:
        bam_file = bam_by_sample[sample]
        idx_counts = {}
        idx_text = (pysam.idxstats(bam_file) if bam_file.endswith(".bam")
                    else open(bam_file).read())
        for line in idx_text.splitlines():
            parts = line.split("\t")
            if len(parts) >= 4 and parts[0] != "*":
                gene_id = parts[0].rsplit(".",1)[0]
                length_bp = int(parts[1])
                if gene_id not in gene_lengths or length_bp > gene_lengths[gene_id]:
                    gene_lengths[gene_id] = length_bp
                idx_counts[gene_id] = idx_counts.get(gene_id, 0) + int(parts[2]) / 2.0
        if not use_em:
            for g, c in idx_counts.items():
                master_counts.setdefault(g, {})[sample] = c; gene_src.setdefault(g, "idxstats")
            continue

        if use_em:
            sample_dir = os.path.dirname(bam_file)
            em_file = find_em_file(sample_dir)
            if em_file is None:
                n_missing_em += 1
                print(f"[WARN] {sample}: no EM output found; multi-copy genes will be "
                      f"missing for this sample (did Step 03 run?)")
                for g, c in idx_counts.items():
                    if g in orphan_set:
                        master_counts.setdefault(g, {})[sample] = c; gene_src.setdefault(g, "idxstats")
                continue
            df_em = pd.read_csv(em_file, sep="\t")
            df_em["Clean_Gene"] = df_em["AlleleID"].astype(str).str.rsplit(".",n=1).str[0]
            # transcripts of one gene (.t1/.t2/...) are summed per gene, not overwritten
            em_sum = df_em.groupby("Clean_Gene")["EM_Allocated_Reads"].sum()
            orph_em = set(df_em.loc[df_em["ClusterID"].astype(str) == "ORPHAN", "Clean_Gene"])
            for g, reads in em_sum.items():
                master_counts.setdefault(g, {})
                master_counts[g][sample] = float(reads)
                gene_src[g] = "em_orphan" if g in orph_em else "em"
            # [A6] fallback: genes not covered by the EM are counted from idxstats
            #      (previously only genes in orphans.tsv were, so genes in neither
            #      clusters.tsv nor orphans.tsv vanished). The global engine covers every gene.
            for g, c in idx_counts.items():
                if g not in em_sum.index:
                    master_counts.setdefault(g, {})[sample] = c
                    gene_src.setdefault(g, "idxstats")

    if use_em and n_missing_em == len(sample_names):
        print("[ERROR] no sample has EM output. Run Step 03, or use --mode idxstats.")
        return

    print("[INFO] Computing TPM normalization...")
    final_rows = []
    for gene, counts_dict in master_counts.items():
        length_bp = gene_lengths.get(gene, 1000)
        length_kb = length_bp / 1000.0
        src = gene_src.get(gene, "idxstats")
        if src == "em_orphan":
            gene_type = "Orphan(EM)"
        elif src == "em":
            gene_type = "Polyploid(EM)"
        elif use_em or gene in orphan_set:
            gene_type = "Orphan(BAM_Base)"
        else:
            gene_type = "Polyploid(idxstats)"
        row = {"GeneID": gene, "Type": gene_type, "Length_bp": length_bp}
        for sample in sample_names:
            reads = counts_dict.get(sample, 0.0)
            row[f"{sample}_RawReads"] = reads
            row[f"{sample}_RPK"] = reads / length_kb if length_kb > 0 else 0.0
        final_rows.append(row)

    merged_df = pd.DataFrame(final_rows)
    for sample in sample_names:
        rpk_col = f"{sample}_RPK"
        tpm_col = f"{sample}_TPM"
        total_rpk = merged_df[rpk_col].sum()
        scaling = total_rpk / 1_000_000.0 if total_rpk > 0 else 1.0
        merged_df[tpm_col] = (merged_df[rpk_col] / scaling).round(2)
        merged_df.drop(columns=[rpk_col], inplace=True)

    merged_df.sort_values(by=["Type", "GeneID"], ascending=[False, True], inplace=True)
    out_path = os.path.join(args.results_dir, args.out)
    merged_df.to_csv(out_path, sep="\t", index=False)
    n_poly = int((merged_df["Type"] == "Polyploid(EM)").sum())
    print(f"[INFO] {len(merged_df):,} genes ({n_poly:,} allele copies / "
          f"{len(merged_df)-n_poly:,} single-copy)")
    print(f"[{get_now()}] [SUCCESS] TPM matrix: {out_path}")


if __name__ == "__main__":
    main()
