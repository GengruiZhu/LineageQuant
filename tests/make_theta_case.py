#!/usr/bin/env python3
"""Synthetic case for `lqcore theta` vs theta/global_theta.py, with the awkward inputs on purpose:
isoforms (.1/.2) of one copy, identical copies, a repeated cluster id, a trailing comma, an allele missing
from the FASTA, IUPAC letters and N in transcripts, lowercase and multi-line FASTA with spaces, reads with N,
lowercase and IUPAC letters, empty sequence lines, R2 longer than R1, a truncated last record."""
import gzip, os, random, sys

out = sys.argv[1]
seed = int(sys.argv[2]) if len(sys.argv) > 2 else 1
random.seed(seed)
os.makedirs(out, exist_ok=True)
B = "ACGT"


def rnd(n):
    return "".join(random.choice(B) for _ in range(n))


def mutate(s, rate):
    s = list(s)
    for i in range(len(s)):
        if random.random() < rate:
            s[i] = random.choice(B.replace(s[i], "") if s[i] in B else B)
    return "".join(s)


def rc(s):
    t = {"A": "T", "C": "G", "G": "C", "T": "A"}
    return "".join(t.get(c, c) for c in reversed(s))


tx = {}            # transcript id -> sequence
clusters = []      # (cid, [allele ids as written])
gi = 0
for c in range(60):
    base = rnd(random.randint(300, 1500))
    n = random.randint(2, 6)
    al = []
    for a in range(n):
        gi += 1
        g = f"Gene{gi:05d}"
        s = mutate(base, random.choice([0.0, 0.002, 0.01, 0.03]))
        if a == 1 and c % 7 == 0:
            s = tx[al[0] + ".1"]                       # identical copy
        if c % 5 == 0 and a == 0:
            p = random.randrange(50, len(s) - 50)
            s = s[:p] + random.choice("RYKMSW") + s[p + 1:]   # IUPAC letter
        if c % 6 == 0 and a == 0:
            p = random.randrange(50, len(s) - 50)
            s = s[:p] + "NNN" + s[p + 3:]
        tx[g + ".1"] = s
        if a == 0 and c % 4 == 0:
            tx[g + ".2"] = s[: len(s) // 2] + rnd(100)   # isoform of the same copy
        al.append(g)
    clusters.append((f"CL{c:04d}", al))
# FASTA: lowercase chunks, line wrapping, a space, a description
with open(f"{out}/tx.fa", "w") as f:
    for i, (t, s) in enumerate(tx.items()):
        w = s.lower() if i % 9 == 0 else s
        f.write(f">{t} desc={i}\n")
        for j in range(0, len(w), 61):
            line = w[j:j + 61]
            if i % 13 == 0 and j == 0:
                line = line[:10] + " " + line[10:]
            f.write(line + "\n")
with open(f"{out}/clusters.txt", "w") as f:
    f.write("# comment line\n\n")
    for k, (cid, al) in enumerate(clusters):
        names = [a + ".1" if k % 2 else a for a in al]
        if k == 3:
            names.append("Missing00001")
        if k == 4:
            names.append(names[0])                     # duplicate allele within a cluster
        line = cid + "\t" + ",".join(names) + ("," if k == 5 else "")
        f.write(line + "\n")
    f.write(f"{clusters[2][0]}\t{','.join(clusters[2][1][:2])}\n")   # repeated cluster id (Python dict reset)
# reads
ids = list(tx)
w = [random.random() ** 3 for _ in ids]


def frag():
    t = random.choices(ids, weights=w)[0]
    s = tx[t]
    L = random.randint(150, min(400, len(s)))
    p = random.randint(0, len(s) - L)
    fr = s[p:p + L]
    r1, r2 = fr[:100], rc(fr)[:100]
    if random.random() < 0.5:
        r1, r2 = r2, r1
    return mutate(r1, 0.005), mutate(r2, 0.005)


def rec(i, s):
    return f"@r{i}\n{s}\n+\n{'I' * len(s)}\n"


N = 6000
with gzip.open(f"{out}/R1.fq.gz", "wt") as f1, gzip.open(f"{out}/R2.fq.gz", "wt") as f2:
    for i in range(N):
        a, b = frag()
        if i % 50 == 0:
            a = a[:30] + "N" + a[31:]
        if i % 77 == 0:
            b = b.lower()
        if i % 91 == 0:
            a = a[:40] + "R" + a[41:]
        if i % 333 == 0:
            a = ""                                     # empty sequence line
        if i % 401 == 0:
            b = "  " + b + " "                        # whitespace around the sequence
        f1.write(rec(i, a))
        f2.write(rec(i, b))
    for i in range(N, N + 7):                          # R2 longer
        f2.write(rec(i, frag()[1]))
    f1.write("@trunc\nACGTACGTACGTACGTACGTACGTACGT\n")    # truncated last record of R1
# long-read prior for some copies
with open(f"{out}/lr.tsv", "w") as f:
    f.write("Allele_ID\tCount\n")
    for k, (cid, al) in enumerate(clusters):
        if k % 3 == 0:
            for a in al:
                f.write(f"{a}\t{random.choice([0.0, 1.0, 2.5, 7.0, 12.0, 0.333, 150000.0])}\n")
print(f"{len(tx)} transcripts, {len(clusters)} clusters, {N} pairs -> {out}")
