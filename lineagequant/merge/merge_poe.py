# -*- coding: utf-8 -*-
"""
merge_poe.py -- Product-of-Experts fusion of k-mer + SNP lineage signals (Step 02b).

Aligns two independent Stage-02 outputs by read_name and fuses their lineage
probability vectors via entropy-weighted PoE. The fuser is AGNOSTIC to the SNP
source (MSA vs VCF) -- it only consumes probability vectors. Lineage names are
read from the P(...) column headers, so any number of lineages is supported.

  P_merged(L) propto (P_kmer(L)+eps)^alpha_kmer * (P_snp(L)+eps)^alpha_snp

alpha derives from per-read information weight = discrete base + entropy bonus.
"""
import argparse
import numpy as np
from collections import defaultdict


def entropy(vec):
    vec = np.asarray(vec, dtype=float)
    s = vec.sum()
    if s <= 0:
        return 0.0
    vec = vec / s
    vec = vec[vec > 0]
    if len(vec) <= 1:
        return 0.0
    return float(-np.sum(vec * np.log2(vec)))


def load_pvalues(path, source_label):
    data = {}
    with open(path) as f:
        header = f.readline().strip().split("\t")
        p_cols = {}
        for i, col in enumerate(header):
            if col.startswith("P(") and col.endswith(")"):
                p_cols[col[2:-1].upper()] = i
        lineage_names = list(p_cols.keys())
        n_lin = len(lineage_names)
        # status column by name: snipe_kmer writes Lineage_Status, N_Hits, Best_AS and
        # snipe_snp writes Lineage_Status, N_SNPs, so the status is not the last column
        st_idx = next((i for i, c in enumerate(header) if c == "Lineage_Status"),
                      next((i for i, c in enumerate(header) if "status" in c.lower()), len(header) - 1))
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            qname, gene_id = parts[0], parts[1]
            status = parts[st_idx] if st_idx < len(parts) else parts[-1]
            probs = [float(parts[p_cols[lin]]) for lin in lineage_names]
            if source_label == "kmer":
                base = 10.0 if status.startswith("Specific") else \
                       (2.5 if status.startswith("Conflict") else 0.5)
            else:
                base = 10.0 if status.startswith("Specific_Soft") else \
                       (2.5 if status.startswith("Multi") else 0.5)
            h = entropy(probs)
            continuous = max(0.0, (np.log2(n_lin) - h)) * 5.0
            data[qname] = (gene_id, probs, base + continuous)
    print(f"  [{source_label}] {len(data):,} fragments from {path}")
    return data, lineage_names


def main():
    parser = argparse.ArgumentParser(description="PoE merge of k-mer + SNP signals (Step 02b)")
    parser.add_argument("--kmer", required=True)
    parser.add_argument("--snp", required=True)
    parser.add_argument("--out", default="P_values_merged.tsv")
    parser.add_argument("--min_weight", type=float, default=0.05)
    parser.add_argument("--eps", type=float, default=1e-6)
    args = parser.parse_args()

    kmer_data, kmer_lins = load_pvalues(args.kmer, "kmer")
    snp_data, snp_lins = load_pvalues(args.snp, "snp")
    all_lins = kmer_lins if kmer_lins == snp_lins else sorted(set(kmer_lins + snp_lins))
    n_lin = len(all_lins)
    eps = args.eps

    all_reads = set(kmer_data) | set(snp_data)
    both = kmer_only = snp_only = 0
    with open(args.out, "w") as out_f:
        out_f.write("Read_Name\tTarget_Gene\t" +
                    "\t".join(f"P({l})" for l in all_lins) + "\tLineage_Status\n")
        for qname in all_reads:
            in_k, in_s = qname in kmer_data, qname in snp_data
            if in_k and in_s:
                both += 1
                gene_k, p_k, w_k = kmer_data[qname]
                gene_s, p_s, w_s = snp_data[qname]
                gene_id = gene_k or gene_s
                tot_w = w_k + w_s
                alpha_k = (w_k / tot_w) if tot_w > 0 else 0.5
                alpha_k = max(args.min_weight, min(1 - args.min_weight, alpha_k))
                alpha_s = 1 - alpha_k
                merged_raw = []
                for lin in all_lins:
                    pk = p_k[kmer_lins.index(lin)] if lin in kmer_lins else 1.0 / n_lin
                    ps = p_s[snp_lins.index(lin)] if lin in snp_lins else 1.0 / n_lin
                    merged_raw.append((pk + eps) ** alpha_k * (ps + eps) ** alpha_s)
                tot = sum(merged_raw)
                merged = [v / tot for v in merged_raw] if tot > 0 else [1.0 / n_lin] * n_lin
                mp_ = max(merged)
                status = (f"Specific_Soft({all_lins[merged.index(mp_)]})"
                          if mp_ > 0.9 else "Multi_Lineage")
            elif in_k:
                kmer_only += 1
                gene_id, p_k, _ = kmer_data[qname]
                merged = [p_k[kmer_lins.index(l)] if l in kmer_lins else 1.0 / n_lin
                          for l in all_lins]
                mp_ = max(merged)
                status = (f"Specific_Kmer({all_lins[merged.index(mp_)]})"
                          if mp_ > 0.9 else "Multi_Kmer")
            else:
                snp_only += 1
                gene_id, p_s, _ = snp_data[qname]
                merged = [p_s[snp_lins.index(l)] if l in snp_lins else 1.0 / n_lin
                          for l in all_lins]
                mp_ = max(merged)
                status = (f"Specific_SNP({all_lins[merged.index(mp_)]})"
                          if mp_ > 0.9 else "Multi_SNP")
            out_f.write(f"{qname}\t{gene_id}\t" +
                        "\t".join(f"{merged[i]:.4f}" for i in range(n_lin)) +
                        f"\t{status}\n")

    total = both + kmer_only + snp_only
    print(f"\n[SUMMARY] {total:,} fragments | both:{both:,} kmer-only:{kmer_only:,} "
          f"snp-only:{snp_only:,}")
    print(f"[SUCCESS] {args.out}")


if __name__ == "__main__":
    main()
