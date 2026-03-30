# -*- coding: utf-8 -*-
"""
LineageQuant Module 03 – EM Deconvolution Engine
Iteratively allocates ambiguously mapped reads to alleles using an
Expectation-Maximization algorithm, producing final allele-resolved
expression values.

Supports any number (≥2) of ancestral lineages; lineage columns are
auto-detected from the P-value table header (columns named P(NAME)).
"""

import sys
import os
from collections import defaultdict


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="LineageQuant Module 03: EM Deconvolution Engine"
    )
    parser.add_argument("-t", "--theta",    required=True,
                        help="Theta matrix TSV produced by Module 01")
    parser.add_argument("-p", "--pvals",    required=True,
                        help="P-value matrix TSV produced by Module 02")
    parser.add_argument("-a", "--ancestry", required=True,
                        help="Gene ancestry table: GeneID<TAB>LineageName")
    parser.add_argument("-o", "--out",      default="EM_Final_Expression.tsv",
                        help="Output file path (default: EM_Final_Expression.tsv)")
    parser.add_argument("--epsilon",   type=float, default=0.01,
                        help="Laplace smoothing pseudo-count (default: 0.01)")
    parser.add_argument("--max-iter",  dest="max_iter", type=int, default=200,
                        help="Maximum EM iterations (default: 200)")
    parser.add_argument("--tol",       type=float, default=1e-5,
                        help="Convergence tolerance (default: 1e-5)")
    args = parser.parse_args()

    print("\n" + "=" * 60)
    print("[INFO] LineageQuant EM Deconvolution Engine starting...")

    # ---- Step 1: Load ancestry table ----
    ancestry_dict = {}
    with open(args.ancestry, 'r', encoding='utf-8') as f:
        for line in f:
            if not line.strip() or line.startswith('GeneID'):
                continue
            parts = line.strip().split()
            if len(parts) >= 2:
                gene_id = parts[0].split('.')[0]
                ancestry_dict[gene_id] = parts[1].upper()
    print(f"  -> Loaded ancestry labels for {len(ancestry_dict)} genes.")

    # ---- Step 2: Load initial Theta values ----
    cluster_alleles = defaultdict(list)
    initial_hits    = {}
    gene_to_cluster = {}

    with open(args.theta, 'r', encoding='utf-8') as f:
        f.readline()  # skip header
        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split('\t')
            cid, aid = parts[0], parts[1]
            cluster_alleles[cid].append(aid)
            initial_hits[aid]    = float(parts[3])
            gene_to_cluster[aid] = cid

    # ---- Step 3: Auto-detect lineage columns from P-value header ----
    cluster_reads = defaultdict(dict)
    p_cols        = {}   # {lineage_name: column_index}

    with open(args.pvals, 'r', encoding='utf-8') as f:
        header = f.readline().strip().split('\t')

        for i, col in enumerate(header):
            if col.startswith("P(") and col.endswith(")"):
                lin_name       = col[2:-1].upper()
                p_cols[lin_name] = i

        print(f"  -> Detected lineage dimensions: {list(p_cols.keys())}")

        for line in f:
            if not line.strip():
                continue
            parts      = line.strip().split('\t')
            read_name  = parts[0]
            target_gene = parts[1]

            if target_gene in gene_to_cluster:
                cid = gene_to_cluster[target_gene]
                if read_name not in cluster_reads[cid]:
                    p_dict = {lin: float(parts[idx]) for lin, idx in p_cols.items()}
                    cluster_reads[cid][read_name] = p_dict

    total_reads = sum(len(reads) for reads in cluster_reads.values())
    print(f"  -> Matched {total_reads} reads to cluster loci.")

    # ---- Step 4: EM iterations per cluster ----
    print("\n[INFO] Running EM iterations per cluster...")
    final_results  = []
    num_lineages   = len(p_cols) if p_cols else 2

    for cid, alleles in cluster_alleles.items():
        smoothed_hits = {a: (initial_hits[a] + args.epsilon) for a in alleles}
        total_smoothed = sum(smoothed_hits.values())
        theta          = {a: smoothed_hits[a] / total_smoothed for a in alleles}

        reads_info = cluster_reads.get(cid, {})
        if not reads_info:
            for a in alleles:
                final_results.append(
                    (cid, a, ancestry_dict.get(a, "UNKNOWN"),
                     initial_hits[a], 0.0, theta[a])
                )
            continue

        # EM loop
        for _ in range(args.max_iter):
            old_theta      = theta.copy()
            expected_counts = {a: 0.0 for a in alleles}

            # E-step: compute expected counts
            for read_name, p_dict in reads_info.items():
                weights    = {}
                weight_sum = 0.0
                for a in alleles:
                    lin    = ancestry_dict.get(a, "UNKNOWN")
                    # Fall back to uniform prior for unrecognised lineage labels
                    prior_p = p_dict.get(lin, 1.0 / num_lineages)
                    w       = theta[a] * prior_p
                    weights[a] = w
                    weight_sum += w

                if weight_sum > 0:
                    for a in alleles:
                        expected_counts[a] += weights[a] / weight_sum

            # M-step: update Theta
            total_expected = sum(expected_counts.values()) + len(alleles) * args.epsilon
            for a in alleles:
                theta[a] = (expected_counts[a] + args.epsilon) / total_expected

            # Convergence check
            diff = sum(abs(theta[a] - old_theta[a]) for a in alleles)
            if diff < args.tol:
                break

        for a in alleles:
            final_results.append(
                (cid, a, ancestry_dict.get(a, "UNKNOWN"),
                 initial_hits[a], expected_counts[a], theta[a])
            )

    # ---- Step 5: Write output ----
    with open(args.out, 'w', encoding='utf-8') as fout:
        fout.write("ClusterID\tAlleleID\tAncestry\t"
                   "Initial_Kmer_Hits\tEM_Allocated_Reads\tFinal_Theta\n")
        for res in final_results:
            fout.write(
                f"{res[0]}\t{res[1]}\t{res[2]}\t"
                f"{res[3]:.1f}\t{res[4]:.2f}\t{res[5]:.6f}\n"
            )

    print(f"\n[SUCCESS] EM convergence complete. Output written to: {args.out}")


if __name__ == "__main__":
    main()
