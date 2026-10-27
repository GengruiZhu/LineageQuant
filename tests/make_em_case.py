#!/usr/bin/env python3
"""Synthetic input for comparing `lqcore em` with lineagequant/core/em_global.py.

Covers: allele clusters + orphans, genes with several transcripts (.t1/.t2 and .1/.2),
identical-sequence genes (pooling), Theta-only alleles without any transcript, long-read
anchored clusters (tau prior), lineage P rows with mixed statuses and repeated reads,
mismatch differences between candidates, secondary records with SEQ '*', unmapped pairs,
proper / improper pairs.

usage: make_em_case.py OUTDIR SEED [N_FRAGMENTS]
"""
import os, random, sys

out, seed = sys.argv[1], int(sys.argv[2])
NF = int(sys.argv[3]) if len(sys.argv) > 3 else 20000
random.seed(seed)
os.makedirs(out, exist_ok=True)
B = "ACGT"
LIN = ["SOF", "SRO", "SSP"]


def rseq(n):
    return "".join(random.choice(B) for _ in range(n))


def mutate(s, rate):
    return "".join(random.choice(B) if random.random() < rate else c for c in s)


genes = {}          # gene -> list of (transcript name, seq)
clusters = {}       # cluster -> genes
lineage = {}
gi = 0
for c in range(40):
    base = rseq(random.randint(300, 2500))
    members = []
    for k in range(random.randint(2, 7)):
        gi += 1
        g = f"Sc{gi:05d}"
        if random.random() < 0.15 and members:                       # identical copy -> pooling
            seq = genes[members[0]][0][1]
        else:
            seq = mutate(base, random.choice([0.002, 0.01, 0.03]))
        style = random.random()
        if style < 0.6:
            tx = [(g, seq)]
        elif style < 0.85:
            tx = [(f"{g}.1", seq)]
        else:
            tx = [(f"{g}.t1", seq), (f"{g}.t2", seq[: max(200, len(seq) // 2)])]
        genes[g] = tx
        members.append(g)
        lineage[g] = random.choice(LIN)
    clusters[f"Global_Cluster_{c:06d}"] = members
for o in range(60):                                                  # orphans
    gi += 1
    g = f"Or{gi:05d}"
    seq = rseq(random.randint(300, 3000))
    if random.random() < 0.3:
        genes[g] = [(f"{g}.t1", seq), (f"{g}.t2", seq[100:])]
    elif random.random() < 0.1 and o > 0:                            # identical orphans
        prev = [x for x in genes if x.startswith("Or")][-1]
        genes[g] = [(g, genes[prev][0][1])]
    else:
        genes[g] = [(g, seq)]
# Theta-only alleles (in clusters, no transcript in the reference)
clusters["Global_Cluster_000000"].append("ScGhost1")
lineage["ScGhost1"] = "SRO"

refs = [(t, s) for g in genes for t, s in genes[g]]
with open(f"{out}/ref.fa", "w") as f:
    for t, s in refs:
        f.write(f">{t} some description\n")
        for i in range(0, len(s), 70):
            f.write(s[i:i + 70].lower() if random.random() < 0.1 else s[i:i + 70])
            f.write("\n")

with open(f"{out}/theta.tsv", "w") as f:
    f.write("Cluster_ID\tAllele_ID\tUnique_Probes\tRaw_Hits\tLR_Count\tTheta\tPrior_Status\tShrinkage_Used\n")
    for c, mem in clusters.items():
        anchored = random.random() < 0.3
        for g in mem:
            lr = random.randint(0, 40) if anchored else 0
            th = random.random()
            f.write(f"{c}\t{g}\t{random.randint(0, 50)}\t{random.randint(0, 500)}\t{lr}\t{th:.6f}\t"
                    f"{'LongRead_Anchored' if anchored else 'Kmer_Inferred'}\tNA\n")

with open(f"{out}/ancestry.tsv", "w") as f:
    f.write("GeneID\tAncestry\n")
    for g, l in lineage.items():
        name = g if random.random() < 0.8 else g + ".1"
        f.write(f"{name}\t{l.lower() if random.random() < 0.2 else l}\n")

# ---------------------------------------------------------------- fragments
names = [t for t, _ in refs]
tx_of = {t: g for g in genes for t, _ in genes[g]}
gene_list = list(genes)
cl_of = {g: c for c, m in clusters.items() for g in m}
recs = []
pv_rows = []
for i in range(NF):
    q = f"frag{i:07d}"
    if random.random() < 0.05:                                       # unmapped pair
        recs.append(f"{q}\t77\t*\t0\t0\t*\t*\t0\t0\t{rseq(100)}\t{'I'*100}\tYT:Z:UP")
        recs.append(f"{q}\t141\t*\t0\t0\t*\t*\t0\t0\t{rseq(100)}\t{'I'*100}\tYT:Z:UP")
        continue
    g = random.choice(gene_list)
    tgt = [random.choice(genes[g])[0]]
    if g in cl_of and random.random() < 0.8:
        for o in clusters[cl_of[g]]:
            if o in genes and o != g and random.random() < 0.6:
                tgt.append(random.choice(genes[o])[0])
    if random.random() < 0.08:                                       # cross-cluster / orphan hit
        tgt.append(random.choice(names))
    tgt = list(dict.fromkeys(tgt))
    tl = random.randint(150, 600) if random.random() < 0.97 else random.randint(2000, 5000)
    b1, b2 = random.choice([0, 0, 0, 1, 2]), random.choice([0, 0, 0, 1])
    for j, t in enumerate(tgt):
        nm1 = b1 if j == 0 else b1 + random.choice([0, 0, 0, 1, 2, 3])
        nm2 = b2 if j == 0 else b2 + random.choice([0, 0, 0, 1, 2])
        sec = 0 if j == 0 else 256
        prop = 2 if random.random() < 0.95 else 0
        f1, f2 = 1 | prop | 64 | 32 | sec, 1 | prop | 128 | 16 | sec
        s1 = rseq(100) if (sec == 0 or random.random() < 0.7) else "*"
        s2 = rseq(100) if (sec == 0 or random.random() < 0.7) else "*"
        recs.append(f"{q}\t{f1}\t{t}\t1\t255\t100M\t=\t300\t{tl}\t{s1}\t{'I'*len(s1) if s1!='*' else '*'}\tAS:i:0\tNM:i:{nm1}\tYT:Z:CP")
        recs.append(f"{q}\t{f2}\t{t}\t300\t255\t100M\t=\t1\t{-tl}\t{s2}\t{'I'*len(s2) if s2!='*' else '*'}\tAS:i:0\tXS:i:0\tNM:i:{nm2}\tYT:Z:CP")
    if random.random() < 0.35:
        for rep in range(random.randint(1, 3)):
            p = [random.random() for _ in LIN]
            s = sum(p)
            p = [x / s for x in p]
            st = random.choice(["Specific", "Specific(SOF)", "Conflict", "Weak", "Ambiguous", "Multi_Lineage",
                                "Specific_Soft(SRO)", "No_SNP_Covered", "Chimeric_Drop"])
            gene_col = tx_of[tgt[0]]
            pv_rows.append(f"{q}\t{gene_col}\t" + "\t".join(f"{x:.4f}" for x in p) + f"\t{st}\t{random.randint(0, 9)}")
            if random.random() < 0.02:
                pv_rows.append(f"{q}\t{gene_col}\tNA\t0.1\t0.2\tSpecific\t1")      # unparsable -> skipped

with open(f"{out}/case.sam", "w") as f:
    f.write("@HD\tVN:1.0\tSO:unsorted\tGO:query\n")
    for t, s in refs:
        f.write(f"@SQ\tSN:{t}\tLN:{len(s)}\n")
    f.write("@PG\tID:bowtie2\tPN:bowtie2\n")
    f.write("\n".join(recs) + "\n")
with open(f"{out}/pvals.tsv", "w") as f:
    f.write("Read_Name\tTarget_Gene\tP(SOF)\tP(SRO)\tP(SSP)\tLineage_Status\tN_Hits\n")
    f.write("\n".join(pv_rows) + "\n")
print(f"{out}: {len(refs)} transcripts, {len(genes)} genes, {len(clusters)} clusters, {NF} fragments, {len(pv_rows)} P rows")
