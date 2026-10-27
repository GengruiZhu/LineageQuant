# -*- coding: utf-8 -*-
# TEST REFERENCE ONLY: the Python aggregator of the 1.1.0 development router (adds --tpm-length eff to the 1.0.0
# aggregator); tests/compare_matrix.sh compares `lineagequant matrix` with it.
"""
core/aggregator.py -- expression matrix aggregator (Step 04), LineageQuant 1.1.0.

Collects per-sample EM allocations (lqcore run: EM_Final_*.tsv, <sample>.idxstats.tsv, <sample>.efflen.tsv) into one
genome-wide matrix with, per sample, <sample>_RawReads (EM-assigned fragments) and <sample>_TPM.

  --tpm-length full (default, = 1.0.0): TPM from RawReads per kb of the longest transcript.
  --tpm-length eff  (1.1.0): TPM from RawReads / the effective length the EM used (lqcore --eff-out; includes the
                    --gc-bias factor). Genes without an effective length (not in the EM) fall back to full length.

Genes absent from the EM output are counted from idxstats [A6]. Type column: Polyploid(EM) (allele copy),
Orphan(EM) (single-copy gene, from the global EM), Orphan(BAM_Base) (idxstats fallback).

EM file search order: EM_Final_Hybrid.tsv -> EM_Final_Expression_v2.tsv -> EM_Final_Expression.tsv.
"""
import os
import glob
import datetime
import argparse
import pandas as pd


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
                    orphans.add(g.rsplit(".", 1)[0])
            elif len(parts) == 1:
                orphans.add(parts[0].rsplit(".", 1)[0])
    return orphans


def find_em_file(sample_dir):
    for name in ("EM_Final_Hybrid.tsv", "EM_Final_Expression_v2.tsv",
                 "EM_Final_Expression.tsv"):
        p = os.path.join(sample_dir, name)
        if os.path.exists(p):
            return p
    return None


def read_efflen(path):
    """lqcore --eff-out: Target, Length, EffectiveLength, GC_Factor -> {gene: [(target, eff)]}"""
    eff = {}
    if not (path and os.path.exists(path)):
        return eff
    with open(path) as f:
        f.readline()
        for l in f:
            x = l.rstrip("\n").split("\t")
            if len(x) >= 3:
                eff[x[0]] = float(x[2])
    return eff


def main():
    parser = argparse.ArgumentParser(description="TPM aggregator (Step 04)")
    parser.add_argument("-r", "--results_dir", required=True)
    parser.add_argument("-o", "--orphan", required=True)
    parser.add_argument("--out", default="LineageQuant_Final_Merged_TPM.tsv")
    parser.add_argument("--mode", choices=["auto", "allo"], default="allo",
                        help="allo / auto: EM_Allocated_Reads for every gene in the EM, idxstats for the rest")
    parser.add_argument("--tpm-length", dest="tpm_length", choices=["full", "eff"], default="full",
                        help="full (default, 1.0.0): reads per kb of the longest transcript; "
                             "eff: reads per effective length used by the EM (<sample>.efflen.tsv)")
    args = parser.parse_args()

    orphan_set = parse_orphan_genes(args.orphan)
    idx_files = glob.glob(os.path.join(args.results_dir, "**", "*.idxstats.tsv"), recursive=True)
    if not idx_files:
        print("[ERROR] No <sample>.idxstats.tsv (lqcore run --idxstats-out) found in the results directories.")
        return
    idx_by_sample = {}
    for bf in idx_files:
        sname = os.path.basename(bf)[:-len(".idxstats.tsv")]
        if sname not in idx_by_sample or bf.count(os.sep) < idx_by_sample[sname].count(os.sep):
            idx_by_sample[sname] = bf
    sample_names = sorted(idx_by_sample.keys())
    print(f"[INFO] Found {len(sample_names)} samples. Mode={args.mode} (EM=yes). Aggregating...")

    master_counts = {}
    master_rate = {}       # --tpm-length eff: sum over the gene's EM targets of reads / effective length
    gene_lengths = {}
    gene_src = {}          # [A6] source of each gene's count: em / em_orphan / idxstats
    n_missing_em = n_no_eff = 0
    for sample in sample_names:
        idx_file = idx_by_sample[sample]
        idx_counts = {}
        with open(idx_file) as f:
            idx_text = f.read()
        for line in idx_text.splitlines():
            parts = line.split("\t")
            if len(parts) >= 4 and parts[0] != "*":
                gene_id = parts[0].rsplit(".", 1)[0]
                length_bp = int(parts[1])
                if gene_id not in gene_lengths or length_bp > gene_lengths[gene_id]:
                    gene_lengths[gene_id] = length_bp
                idx_counts[gene_id] = idx_counts.get(gene_id, 0) + int(parts[2]) / 2.0

        sample_dir = os.path.dirname(idx_file)
        em_file = find_em_file(sample_dir)
        if em_file is None:
            n_missing_em += 1
            print(f"[WARN] {sample}: no EM output found; multi-copy genes will be "
                  f"missing for this sample (did lqcore run finish?)")
            for g, c in idx_counts.items():
                if g in orphan_set:
                    master_counts.setdefault(g, {})[sample] = c
                    gene_src.setdefault(g, "idxstats")
            continue
        df_em = pd.read_csv(em_file, sep="\t")
        df_em["Clean_Gene"] = df_em["AlleleID"].astype(str).str.rsplit(".", n=1).str[0]
        # transcripts of one gene (.t1/.t2/...) are summed per gene, not overwritten
        em_sum = df_em.groupby("Clean_Gene")["EM_Allocated_Reads"].sum()
        orph_em = set(df_em.loc[df_em["ClusterID"].astype(str) == "ORPHAN", "Clean_Gene"])
        for g, reads in em_sum.items():
            master_counts.setdefault(g, {})
            master_counts[g][sample] = float(reads)
            gene_src[g] = "em_orphan" if g in orph_em else "em"
        if args.tpm_length == "eff":
            eff = read_efflen(os.path.join(sample_dir, f"{sample}.efflen.tsv"))
            if not eff:
                n_no_eff += 1
                print(f"[WARN] {sample}: no {sample}.efflen.tsv -> TPM by full length for this sample")
            else:
                for t, g, r in zip(df_em["AlleleID"].astype(str), df_em["Clean_Gene"], df_em["EM_Allocated_Reads"]):
                    e = eff.get(t)
                    if e and e > 0:
                        master_rate.setdefault(g, {})
                        master_rate[g][sample] = master_rate[g].get(sample, 0.0) + float(r) / e
        # [A6] genes not covered by the EM are counted from idxstats
        for g, c in idx_counts.items():
            if g not in em_sum.index:
                master_counts.setdefault(g, {})[sample] = c
                gene_src.setdefault(g, "idxstats")

    if n_missing_em == len(sample_names):
        print("[ERROR] no sample has EM output (lqcore run).")
        return

    print("[INFO] Computing TPM normalization" + (" (effective lengths)" if args.tpm_length == "eff" else "") + "...")
    final_rows = []
    for gene, counts_dict in master_counts.items():
        length_bp = gene_lengths.get(gene, 1000)
        length_kb = length_bp / 1000.0
        src = gene_src.get(gene, "idxstats")
        if src == "em_orphan":
            gene_type = "Orphan(EM)"
        elif src == "em":
            gene_type = "Polyploid(EM)"
        else:
            gene_type = "Orphan(BAM_Base)"
        row = {"GeneID": gene, "Type": gene_type, "Length_bp": length_bp}
        for sample in sample_names:
            reads = counts_dict.get(sample, 0.0)
            row[f"{sample}_RawReads"] = reads
            rate = master_rate.get(gene, {}).get(sample) if args.tpm_length == "eff" else None
            if rate is not None:
                row[f"{sample}_RPK"] = rate * 1000.0          # per kb of effective length
            else:
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
