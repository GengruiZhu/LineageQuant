#!/usr/bin/env python3
"""Synthetic allopolyploid for the end-to-end comparison (tests/compare_pipeline.sh).

  ref/tx.fa clusters.tsv ancestry.tsv orphans.tsv target_ids.txt   copies of lineages A and B (+ a second A copy),
                                                                  some with two isoforms, and single-copy genes
  ref/progA.fa/.gff  ref/progB1.fa/.gff  ref/progB2.fa/.gff         progenitor genomes (>= 1 Mb chromosome + a
                                                                  scaffold); B has two samples
  ref/genome.fa/.gff3 ref/pop.vcf.gz ref/A.samples ref/B.samples    the polyploid genome with every copy as a gene
                                                                  and a population VCF (lineage-informative SNPs)
  fq/S1_R1.fq.gz fq/S1_R2.fq.gz fq/S2_1.fq.gz fq/S2_2.fq.gz         read pairs (100 bp, fragments ~250 bp, 0.3% errors)

usage: make_pipeline_case.py OUTDIR SEED [N_PAIRS]
"""
import gzip, os, random, sys

out, seed = sys.argv[1], int(sys.argv[2])
NP = int(sys.argv[3]) if len(sys.argv) > 3 else 15000
random.seed(seed)
R = os.path.join(out, "ref"); F = os.path.join(out, "fq")
os.makedirs(R, exist_ok=True); os.makedirs(F, exist_ok=True)
B = "ACGT"
rseq = lambda n: "".join(random.choice(B) for _ in range(n))
mutate = lambda s, r: "".join(random.choice(B) if random.random() < r else c for c in s)
wrap = lambda s, w=60: "\n".join(s[i:i + w] for i in range(0, len(s), w))
rc = lambda s: s.translate(str.maketrans("ACGT", "TGCA"))[::-1]

tx, clusters, anc, orph, genes_of = [], [], [], [], {}
for c in range(36):
    base = rseq(random.randint(500, 1600))
    members = []
    for lin, k in (("A", 1), ("B", 1), ("A", 2)):
        if k == 2 and random.random() < 0.5:
            continue
        g = f"Sc{c + 1:03d}{lin}{k}"
        s = mutate(base, 0.03)
        if random.random() < 0.2:
            tx.append((f"{g}.t1", s)); tx.append((f"{g}.t2", s[: len(s) * 2 // 3]))
        else:
            tx.append((f"{g}.t1", s))
        members.append(g); anc.append((g, lin)); genes_of[g] = s
    clusters.append((f"Global_Cluster_{c:06d}", members))
for o in range(20):
    g = f"Or{o + 1:03d}"
    s = rseq(random.randint(400, 1500))
    tx.append((f"{g}.t1", s)); orph.append(g)
with open(os.path.join(R, "tx.fa"), "w") as f:
    for n, s in tx: f.write(f">{n}\n{wrap(s)}\n")
with open(os.path.join(R, "clusters.tsv"), "w") as f:
    for cid, m in clusters: f.write(f"{cid}\t{','.join(sorted(m))}\n")
with open(os.path.join(R, "ancestry.tsv"), "w") as f:
    f.write("GeneID\tAncestry\n")
    for g, l in anc: f.write(f"{g}\t{l}\n")
with open(os.path.join(R, "orphans.tsv"), "w") as f:
    for i, g in enumerate(orph, 1): f.write(f"OrphanGrp_{i:06d}\t{g}\n")
with open(os.path.join(R, "target_ids.txt"), "w") as f:
    for g, _ in anc: f.write(g + "\n")

# progenitor genomes: the copies of their lineage (diverged 2%) as single-exon genes
def genome(path, lin, rate, tag):
    chrom = list(rseq(1_050_000))
    gff = ["##gff-version 3"]
    pos, n = 2000, 0
    for g, l in anc:
        if l != lin or g.endswith("2"):
            continue
        n += 1
        s = mutate(genes_of[g], rate)
        strand = "+" if n % 3 else "-"
        seg = s if strand == "+" else rc(s)
        chrom[pos:pos + len(seg)] = list(seg)
        gid = f"{tag}g{n:03d}"
        gff += [f"Chr01\tsim\tgene\t{pos + 1}\t{pos + len(seg)}\t.\t{strand}\t.\tID={gid}",
                f"Chr01\tsim\tmRNA\t{pos + 1}\t{pos + len(seg)}\t.\t{strand}\t.\tID={gid}.1;Parent={gid}",
                f"Chr01\tsim\texon\t{pos + 1}\t{pos + len(seg)}\t.\t{strand}\t.\tParent={gid}.1"]
        pos += len(seg) + random.randint(2000, 9000)
    sc = rseq(30000)
    gff += [f"scf1\tsim\tgene\t10\t500\t.\t+\t.\tID={tag}scf", f"scf1\tsim\tmRNA\t10\t500\t.\t+\t.\tID={tag}scf.1;Parent={tag}scf",
            f"scf1\tsim\texon\t10\t500\t.\t+\t.\tParent={tag}scf.1"]
    with open(path + ".fa", "w") as f:
        f.write(f">Chr01 {tag} chromosome 1\n{wrap(''.join(chrom), 80)}\n>scf1\n{wrap(sc, 80)}\n")
    open(path + ".gff", "w").write("\n".join(gff) + "\n")
genome(os.path.join(R, "progA"), "A", 0.02, "pA")
genome(os.path.join(R, "progB1"), "B", 0.02, "pB")
genome(os.path.join(R, "progB2"), "B", 0.025, "pC")

# polyploid genome (every copy a gene; transcript ids as in tx.fa) + population VCF
chroms, gff, sites = {}, ["##gff-version 3"], []
for ci, part in enumerate((tx[: len(tx) // 2], tx[len(tx) // 2:])):
    name = f"chr{ci + 1}"
    seq, pos = [rseq(300)], 300
    for t, s in part:
        if not t.endswith(".t1"):
            continue
        g = t[:-3]
        strand = "+" if len(seq) % 2 else "-"
        seg = s if strand == "+" else rc(s)
        gff += [f"{name}\tsim\tgene\t{pos + 1}\t{pos + len(seg)}\t.\t{strand}\t.\tID={g}",
                f"{name}\tsim\tmRNA\t{pos + 1}\t{pos + len(seg)}\t.\t{strand}\t.\tID={t};Parent={g}",
                f"{name}\tsim\texon\t{pos + 1}\t{pos + len(seg)}\t.\t{strand}\t.\tParent={t}",
                f"{name}\tsim\tCDS\t{pos + 1}\t{pos + len(seg)}\t.\t{strand}\t0\tParent={t}"]
        for _ in range(random.randint(0, 6)):
            sites.append((name, pos + random.randint(1, len(seg)), g))
        seq.append(seg); pos += len(seg)
        sp = rseq(random.randint(100, 400)); seq.append(sp); pos += len(sp)
    chroms[name] = "".join(seq)
with open(os.path.join(R, "genome.fa"), "w") as f:
    for n, s in chroms.items(): f.write(f">{n}\n{wrap(s)}\n")
open(os.path.join(R, "genome.gff3"), "w").write("\n".join(gff) + "\n")
lin_of = dict(anc)
samples = {"A": ["a1", "a2", "a3", "a4"], "B": ["b1", "b2", "b3", "b4"]}
for l, ss in samples.items(): open(os.path.join(R, f"{l}.samples"), "w").write("\n".join(ss) + "\n")
allss = samples["A"] + samples["B"]
vl = ["##fileformat=VCFv4.2", '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="a">',
      '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="d">', "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(allss)]
for name, p, g in sorted(set(sites), key=lambda x: (x[0], x[1])):
    ref = chroms[name][p - 1]
    alt = random.choice([b for b in B if b != ref])
    hot = lin_of.get(g, random.choice("AB"))
    cols = []
    for s in allss:
        dp = random.randint(2, 60)
        f_alt = 0.9 if (s[0].upper() == hot) else 0.05
        a = sum(random.random() < f_alt for _ in range(dp))
        cols.append(f"{dp - a},{a}:{dp}")
    vl.append(f"{name}\t{p}\t.\t{ref}\t{alt}\t50\tPASS\t.\tAD:DP\t" + "\t".join(cols))
with gzip.open(os.path.join(R, "pop.vcf.gz"), "wt") as f: f.write("\n".join(vl) + "\n")

# reads
def qual(n): return "I" * n
def err(s): return "".join(random.choice(B) if random.random() < 0.003 else c for c in s)
for si, (n1, n2) in enumerate((("S1_R1.fq.gz", "S1_R2.fq.gz"), ("S2_1.fq.gz", "S2_2.fq.gz"))):
    w = [random.lognormvariate(0, 1.2) for _ in tx]
    tot = sum(w)
    f1, f2 = gzip.open(os.path.join(F, n1), "wt"), gzip.open(os.path.join(F, n2), "wt")
    for i in range(NP):
        r = random.random() * tot
        k = 0
        while r > w[k] and k < len(tx) - 1: r -= w[k]; k += 1
        t, s = tx[k]
        fl = max(120, min(len(s), int(random.gauss(250, 40))))
        if len(s) < 120: continue
        st = random.randint(0, len(s) - fl)
        frag = s[st:st + fl]
        a, b = err(frag[:100]), err(rc(frag)[:100])
        if random.random() < 0.5: a, b = b, a
        f1.write(f"@p{si}_{i}/1\n{a}\n+\n{qual(len(a))}\n"); f2.write(f"@p{si}_{i}/2\n{b}\n+\n{qual(len(b))}\n")
    f1.close(); f2.close()
print(f"{out}: {len(tx)} transcripts, {len(clusters)} clusters, {len(sites)} VCF sites, {NP} pairs x 2 samples")
