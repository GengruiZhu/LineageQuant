# -*- coding: utf-8 -*-
"""
LineageQuant Module 04 – TPM Matrix Aggregator
Merges EM-deconvolved polyploid gene counts with featureCounts-based
orphan (single-copy) gene counts, then computes global TPM values
across all samples.
"""

import os
import glob
import pandas as pd


def parse_orphan_genes(orphan_file):
    """
    Read the orphan (single-copy) gene list.
    Expected format: ClusterID<TAB>Gene1[,Gene2,...] (one gene per cluster).
    """
    print(f"[INFO] Loading orphan gene list: {os.path.basename(orphan_file)}...")
    orphans = set()
    with open(orphan_file, 'r') as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            for g in parts[1].split(','):
                orphans.add(g.split('.')[0])
    print(f"  -> {len(orphans)} orphan genes loaded.")
    return orphans


def parse_gff_lengths(gff_file):
    """Extract gene physical lengths (bp) from GFF3."""
    print("[INFO] Parsing gene lengths from GFF3...")
    lengths = {}
    with open(gff_file, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 9 or parts[2] != 'gene':
                continue
            attr    = parts[8]
            gene_id = None
            for item in attr.split(';'):
                if item.startswith('ID='):
                    gene_id = item.replace('ID=', '').strip()
                    break
            if gene_id:
                lengths[gene_id] = int(parts[4]) - int(parts[3]) + 1
    return lengths


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="LineageQuant Module 04: TPM Matrix Aggregator"
    )
    parser.add_argument("-g", "--gff",         required=True,
                        help="Whole-genome GFF3 annotation file")
    parser.add_argument("-r", "--results-dir", required=True, dest="results_dir",
                        help="Root results directory containing per-sample sub-directories")
    parser.add_argument("-o", "--orphan",      required=True,
                        help="Orphan gene TSV list (single-copy genes for featureCounts fallback)")
    parser.add_argument("--out",               default="LineageQuant_Final_Merged_TPM.tsv",
                        help=(
                            "Output filename. If relative, written inside --results-dir. "
                            "If absolute, written directly to that path. "
                            "(default: LineageQuant_Final_Merged_TPM.tsv)"
                        ))
    args = parser.parse_args()

    # ---- Step 1: Prepare reference data ----
    gene_lengths = parse_gff_lengths(args.gff)
    orphan_set   = parse_orphan_genes(args.orphan)

    # ---- Step 2: Discover sample directories via featureCounts files ----
    fc_files = glob.glob(os.path.join(args.results_dir, "*", "*_featureCounts.txt"))
    if not fc_files:
        print("[ERROR] No featureCounts output files found in results directory.")
        return

    sample_names = [os.path.basename(os.path.dirname(f)) for f in fc_files]
    print(f"[INFO] Found {len(sample_names)} sample(s). Building merged matrix...")

    # master_counts[gene][sample] = raw read count
    master_counts = {}

    # ---- Step 3: Per-sample data extraction ----
    for sample in sample_names:
        sample_dir = os.path.join(args.results_dir, sample)

        # A) Polyploid genes: EM-deconvolved counts
        em_file = os.path.join(sample_dir, "EM_Final_Expression.tsv")
        if os.path.exists(em_file):
            df_em = pd.read_csv(em_file, sep='\t')
            for _, row in df_em.iterrows():
                gene  = row['AlleleID'].split('.')[0]
                reads = row['EM_Allocated_Reads']
                if gene not in master_counts:
                    master_counts[gene] = {}
                master_counts[gene][sample] = reads

        # B) Orphan (single-copy) genes: featureCounts raw counts
        fc_file = os.path.join(sample_dir, f"{sample}_featureCounts.txt")
        if os.path.exists(fc_file):
            df_fc = pd.read_csv(fc_file, sep='\t', comment='#')
            df_fc.rename(
                columns={df_fc.columns[0]: 'Geneid', df_fc.columns[-1]: 'Counts'},
                inplace=True
            )
            for _, row in df_fc.iterrows():
                gene = str(row['Geneid']).split('.')[0]
                # Only include genes on the orphan list; never overwrite EM counts
                if gene in orphan_set:
                    reads = float(row['Counts'])
                    if gene not in master_counts:
                        master_counts[gene] = {}
                    master_counts[gene][sample] = reads

    # ---- Step 4: Compute global TPM ----
    print("[INFO] Computing global TPM values...")
    final_rows = []

    for gene, counts_dict in master_counts.items():
        length_bp  = gene_lengths.get(gene, 1000)
        length_kb  = length_bp / 1000.0
        gene_type  = "Orphan(FC)" if gene in orphan_set else "Polyploid(EM)"
        row        = {'GeneID': gene, 'Type': gene_type, 'Length_bp': length_bp}

        for sample in sample_names:
            reads            = counts_dict.get(sample, 0.0)
            row[f"{sample}_RawReads"] = reads
            row[f"{sample}_RPK"]      = reads / length_kb if length_kb > 0 else 0.0

        final_rows.append(row)

    merged_df = pd.DataFrame(final_rows)

    for sample in sample_names:
        rpk_col        = f"{sample}_RPK"
        tpm_col        = f"{sample}_TPM"
        total_rpk      = merged_df[rpk_col].sum()
        scaling_factor = total_rpk / 1_000_000.0 if total_rpk > 0 else 1.0
        merged_df[tpm_col] = (merged_df[rpk_col] / scaling_factor).round(2)
        merged_df.drop(columns=[rpk_col], inplace=True)

    # ---- Step 5: Sort and write output ----
    merged_df.sort_values(
        by=['Type', 'GeneID'], ascending=[False, True], inplace=True
    )

    out_path = os.path.join(args.results_dir, args.out)
    merged_df.to_csv(out_path, sep='\t', index=False)

    print("\n" + "=" * 60)
    print(f"[SUCCESS] Global TPM matrix complete.")
    print(f"  Polyploid (EM) genes : {len([g for g in master_counts if g not in orphan_set])}")
    print(f"  Orphan (FC) genes    : {len([g for g in master_counts if g in orphan_set])}")
    print(f"  Output               : {out_path}")
    print("=" * 60)


if __name__ == "__main__":
    main()
