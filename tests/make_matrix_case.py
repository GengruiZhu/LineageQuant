#!/usr/bin/env python3
"""Synthetic lqcore text outputs for many genes (the matrix sums over > 8192 genes run in numpy's blocks).
usage: make_matrix_case.py OUTDIR N_GENES N_SAMPLES SEED"""
import os, random, sys
out, G, NS, seed = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
random.seed(seed)
genes = []
for i in range(G):
    g = f"G{i:06d}"
    k = random.choice([1, 1, 1, 2, 3])
    tx = [g] if k == 1 else [f"{g}.t{j + 1}" for j in range(k)]
    genes.append((g, tx, random.random() < 0.3))
os.makedirs(out, exist_ok=True)
with open(os.path.join(out, "orphans.tsv"), "w") as f:
    for i, (g, tx, orph) in enumerate(genes):
        if orph:
            f.write(f"OrphanGrp_{i:06d}\t{tx[0]}\n")
for s in range(NS):
    d = os.path.join(out, "res", f"R{s}")
    os.makedirs(d, exist_ok=True)
    em = open(os.path.join(d, "EM_Final_Expression.tsv"), "w")
    em.write("ClusterID\tAlleleID\tAncestry\tInitial_Theta\tEM_Allocated_Reads\tFinal_Theta\tEM_Iter\tEM_LogLik\n")
    ef = open(os.path.join(d, f"R{s}.efflen.tsv"), "w")
    ef.write("Target\tLength\tEffectiveLength\tGC_Factor\n")
    ix = open(os.path.join(d, f"R{s}.idxstats.tsv"), "w")
    for g, tx, orph in genes:
        in_em = random.random() < 0.9
        for t in tx:
            L = random.randint(200, 8000)
            ix.write(f"{t}\t{L}\t{random.randint(0, 5000)}\t0\n")
            if in_em:
                r = random.choice([0.0, random.random() * 10 ** random.randint(0, 5)])
                em.write(f"{'ORPHAN' if orph else 'Global_Cluster_' + g}\t{t}\t{'ORPHAN' if orph else 'A'}\t{random.random():.6f}\t{r:.2f}\t{random.random():.6f}\t57\t0\n")
                e = max(0.0, L - random.random() * 400 + 1) if random.random() > 0.01 else 0.0
                ef.write(f"{t}\t{L}\t{e:.4f}\t{random.random() + 0.5:.5f}\n")
    ix.write("*\t0\t0\t12345\n")
