#!/usr/bin/env python3
"""Synthetic inputs for the reference-preparation steps (tests/compare_prep.sh).

  copies.fa / clusters.tsv / ancestry.tsv   MSA SNP route: 3 lineages, isoforms (.t1/.t2 of different length),
                                            multi-dot ids (style=dots), identical copies, members without sequence or
                                            lineage, a first cluster without lineage L3 (the 1.0.0 lineage bug), N bases
  genome.fa / genome.gff3 / pop.vcf.gz / L*.samples   VCF route: genes on both strands with exon + CDS features,
                                            multi-allelic / indel / missing DP and AD records, records outside genes
  ancX.fa / ancX.gff (X = A, B)             progenitor genomes: descriptions in headers, seqids that differ between
                                            FASTA and GFF (by name, by length), a short scaffold, strand '?', CRLF (B)
  pairs.txt / std.tsv / table.tsv / all_genes.txt   normalize-clusters inputs

usage: make_prep_case.py OUTDIR SEED [dots|plain]
"""
import gzip, os, random, sys

out, seed = sys.argv[1], int(sys.argv[2])
style = sys.argv[3] if len(sys.argv) > 3 else "dots"
random.seed(seed)
os.makedirs(out, exist_ok=True)
B = "ACGT"
LIN = ["L1", "L2", "L3"]


def rseq(n):
    return "".join(random.choice(B) for _ in range(n))


def mutate(s, rate):
    return "".join(random.choice(B) if random.random() < rate else c for c in s)


def wrap(s, w=60):
    return "\n".join(s[i:i + w] for i in range(0, len(s), w))


def rc(s):
    return s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


# ------------------------------------------------------------------ MSA route
fa, cl, an = [], [], ["GeneID\tAncestry"]
gi = 0
for c in range(70):
    base = rseq(random.randint(150, 900))
    members = []
    lins = LIN[:2] if c == 0 else LIN                  # first cluster: no L3
    ident = random.random() < 0.08
    for k in range(random.randint(2, 5)):
        gi += 1
        lin = lins[k % len(lins)] if k < len(lins) else random.choice(lins)
        if style == "dots" and random.random() < 0.25:
            g = f"Sp.{lin}{gi:05d}"                      # two dots in the transcript id below
        else:
            g = f"{lin}g{gi:05d}"
        s = base if ident else mutate(base, random.choice([0.005, 0.02, 0.05]))
        if random.random() < 0.05:
            s = s[:50] + "N" + s[51:]
        r = random.random()
        if r < 0.5:
            fa.append((f"{g}.t1", s))
        elif r < 0.75:
            fa.append((f"{g}.t1", s[: len(s) // 2]))
            fa.append((f"{g}.t2", s))
        elif r < 0.9:
            fa.append((g, s))
        else:
            fa.append((f"{g}.1", s))
        members.append(g)
        if random.random() < 0.95:
            an.append(f"{g}\t{lin.lower() if random.random() < 0.2 else lin}")
    if random.random() < 0.1:
        members.append(f"NoSeq{c:04d}")
    cl.append(f"Global_Cluster_{c:06d}\t" + ",".join(m + (".t1" if random.random() < 0.3 else "") for m in members))
cl.insert(5, "# comment line")
with open(os.path.join(out, "copies.fa"), "w") as f:
    for n, s in fa:
        f.write(f">{n} some description\n{wrap(s)}\n")
open(os.path.join(out, "clusters.tsv"), "w").write("\n".join(cl) + "\n\n")
open(os.path.join(out, "ancestry.tsv"), "w").write("\n".join(an) + "\n")

# ------------------------------------------------------------------ VCF route (reference genome + GFF3 + VCF)
chroms = {f"chr{i}": rseq(random.randint(30000, 60000)) for i in (1, 2)}
gff, genes_pos = ["##gff-version 3"], []
for ch, s in chroms.items():
    pos = 500
    n = 0
    while pos < len(s) - 6000:
        n += 1
        strand = random.choice("+-")
        gid = f"{ch}g{n:03d}"
        exons, p = [], pos
        for e in range(random.randint(1, 4)):
            L = random.randint(80, 600)
            exons.append((p, p + L - 1))
            p += L + random.randint(50, 700)
        gend = exons[-1][1]
        gff.append(f"{ch}\tsim\tgene\t{pos}\t{gend}\t.\t{strand}\t.\tID={gid}")
        gff.append(f"{ch}\tsim\tmRNA\t{pos}\t{gend}\t.\t{strand}\t.\tID={gid}.t1;Parent={gid}")
        with_exon = random.random() < 0.5
        for (b, e) in exons:
            if with_exon:
                gff.append(f"{ch}\tsim\texon\t{b}\t{e}\t.\t{strand}\t.\tParent={gid}.t1")
            gff.append(f"{ch}\tsim\tCDS\t{b}\t{e}\t.\t{strand}\t0\tID=cds.{gid};Parent={gid}.t1")
        genes_pos.append((ch, exons))
        pos = gend + random.randint(200, 2000)
with open(os.path.join(out, "genome.fa"), "w") as f:
    for ch, s in chroms.items():
        f.write(f">{ch}\n{wrap(s)}\n")
open(os.path.join(out, "genome.gff3"), "w").write("\n".join(gff) + "\n")
samples = [f"S{i}" for i in range(1, 13)]
groups = {"L1": samples[0:4], "L2": samples[4:8], "L3": samples[8:11]}     # S12 in no group
for l, ss in groups.items():
    open(os.path.join(out, f"{l}.samples"), "w").write("\n".join(ss + ([" "] if l == "L2" else [])) + "\n")
recs = []
for ch, exons in genes_pos:
    for (b, e) in exons:
        for _ in range(random.randint(0, 4)):
            recs.append((ch, random.randint(b, e)))
    recs.append((ch, exons[0][0] - 20))                 # outside the gene
recs.sort(key=lambda x: (x[0], x[1]))
lines = ["##fileformat=VCFv4.2", '##FORMAT=<ID=GT,Number=1,Type=String,Description="g">',
         '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="a">', '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="d">',
         "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples)]
seen = set()
for ch, p in recs:
    if (ch, p) in seen:
        continue
    seen.add((ch, p))
    ref = chroms[ch][p - 1]
    alt = random.choice([b for b in B if b != ref])
    kind = random.random()
    if kind < 0.05:
        alt = alt + "," + random.choice([b for b in B if b not in (ref, alt)])
    elif kind < 0.08:
        alt = alt + "T"                                  # insertion
    elif kind < 0.1:
        ref = ref + "A"                                  # deletion
    elif kind < 0.12:
        alt = "."
    elif kind < 0.13:
        ref = ref.lower()
    fmt = random.choice(["GT:AD:DP"] * 8 + ["GT:DP", "GT:AD"])
    cols = []
    for s in samples:
        dp = random.choice([random.randint(0, 250), random.randint(3, 60), random.randint(3, 60)])
        a = random.randint(0, dp)
        adv = f"{dp - a},{a}" if random.random() > 0.05 else "."
        dps = str(dp) if random.random() > 0.05 else "."
        if random.random() < 0.03:
            cols.append("0/1")                           # trailing fields missing
            continue
        v = {"GT": "0/1", "AD": adv, "DP": dps}
        cols.append(":".join(v[k] for k in fmt.split(":")))
    lines.append(f"{ch}\t{p}\t.\t{ref}\t{alt}\t50\tPASS\t.\t{fmt}\t" + "\t".join(cols))
with gzip.open(os.path.join(out, "pop.vcf.gz"), "wt") as f:
    f.write("\n".join(lines) + "\n")

# ------------------------------------------------------------------ progenitor genomes for markers / extract-transcripts
def progenitor(tag, crlf, base_genes):
    seqs, gl = {}, ["##gff-version 3"]
    names = [("Chr01", "CHR_01"), ("Chr02", "CM00%s2.1" % tag)]       # by token / by length
    for ci, (fname, gname) in enumerate(names):
        s = rseq(1_000_000 + 20_000 * ci + random.randint(1000, 4000))   # lengths > 0.5% apart
        pos = 1000
        n = 0
        for g in base_genes:
            if pos > len(s) - 5000:
                break
            n += 1
            gid = f"{tag}{fname}g{n:03d}"
            strand = random.choice("+-?") if random.random() < 0.05 else random.choice("+-")
            seg = mutate(g, 0.03)
            s = s[:pos] + seg + s[pos + len(seg):]
            gl.append(f"{gname}\tsim\tgene\t{pos + 1}\t{pos + len(seg)}\t.\t{strand}\t.\tID={gid}")
            gl.append(f"{gname}\tsim\tmRNA\t{pos + 1}\t{pos + len(seg)}\t.\t{strand}\t.\tID={gid}.1;Parent={gid}")
            gl.append(f"{gname}\tsim\texon\t{pos + 1}\t{pos + len(seg)}\t.\t{strand}\t.\tParent={gid}.1")
            pos += len(seg) + random.randint(500, 3000)
        # a feature near the end: the GFF max coordinate fixes the length match for CM00x2.1
        L = len(s)
        gl.append(f"{gname}\tsim\tgene\t{L - 400}\t{L - 100}\t.\t+\t.\tID={tag}{fname}gEnd")
        gl.append(f"{gname}\tsim\tmRNA\t{L - 400}\t{L - 100}\t.\t+\t.\tID={tag}{fname}gEnd.1;Parent={tag}{fname}gEnd")
        gl.append(f"{gname}\tsim\texon\t{L - 400}\t{L - 100}\t.\t+\t.\tParent={tag}{fname}gEnd.1")
        seqs[fname] = s
    seqs["scaffold_9"] = rseq(20000)
    gl.append(f"scaffold_9\tsim\tgene\t100\t900\t.\t+\t.\tID={tag}scg1")
    gl.append(f"scaffold_9\tsim\tmRNA\t100\t900\t.\t+\t.\tID={tag}scg1.1;Parent={tag}scg1")
    gl.append(f"scaffold_9\tsim\texon\t100\t900\t.\t+\t.\tParent={tag}scg1.1")
    gl.append("unknown_seq\tsim\tgene\t1\t50\t.\t+\t.\tID=x")
    nl = "\r\n" if crlf else "\n"
    with open(os.path.join(out, f"anc{tag}.fa"), "w", newline="") as f:
        for n, s in seqs.items():
            f.write(f">{n} {tag} progenitor\tlen={len(s)}{nl}" + nl.join(s[i:i + 80] for i in range(0, len(s), 80)) + nl)
    open(os.path.join(out, f"anc{tag}.gff"), "w").write("\n".join(gl) + "\n")


core = [rseq(random.randint(300, 1200)) for _ in range(60)]
progenitor("A", False, core)
progenitor("B", True, [mutate(g, 0.05) for g in core])

# ------------------------------------------------------------------ normalize-clusters inputs
gn = [f"Gh_{random.choice('AD')}{i:04d}" for i in range(400)]
with open(os.path.join(out, "pairs.txt"), "w") as f:
    f.write("# jcvi anchors\n")
    for i in range(300):
        if i % 50 == 0:
            f.write("###\n")
        a, b = random.choice(gn), random.choice(gn)
        sep = "\t" if random.random() < 0.8 else "  "
        f.write(f"{a}.1{sep}{b}.t2{sep}{random.randint(1, 900)}\n")
    f.write("\n")
with open(os.path.join(out, "std.tsv"), "w") as f:
    for i in range(80):
        f.write(f"c{i}\t" + ",".join(random.sample(gn, random.randint(1, 4))) + ("," if i % 9 == 0 else "") + "\n")
    f.write("lonely\n")
with open(os.path.join(out, "table.tsv"), "w") as f:
    f.write("id\tA\tD\tE\n")
    for i in range(60):
        row = [f"row{i}"] + [random.choice(gn + ["."]) if random.random() > 0.15 else "" for _ in range(3)]
        f.write("\t".join(row) + "\n")
with open(os.path.join(out, "all_genes.txt"), "w") as f:
    for g in gn + [f"Gh_X{i:03d}" for i in range(20)]:
        f.write(g + ".1\n")
print(f"{out}: {len(fa)} transcripts, {len(cl)} cluster lines, {len(recs)} VCF positions")
