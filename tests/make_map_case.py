#!/usr/bin/env python3
"""Synthetic polyploid-like transcriptome + read pairs for testing lqmap against an exhaustive aligner.
Clusters of 6-10 copies diverged 0-4 %, identical copies, isoforms, N runs in transcripts; read pairs with
sequencing errors, from 'unseen alleles' (5 % divergent), with N, near transcript ends, and from nowhere."""
import gzip, os, random, sys

out, seed = sys.argv[1], int(sys.argv[2])
npairs = int(sys.argv[3]) if len(sys.argv) > 3 else 600
rl = int(sys.argv[4]) if len(sys.argv) > 4 else 100
random.seed(seed)
os.makedirs(out, exist_ok=True)
B = "ACGT"
rnd = lambda n: "".join(random.choice(B) for _ in range(n))


def mutate(s, rate):
    return "".join(random.choice(B.replace(c, "")) if (c in B and random.random() < rate) else c for c in s)


def rc(s):
    t = {"A": "T", "C": "G", "G": "C", "T": "A"}
    return "".join(t.get(c, "N") for c in reversed(s))


tx = {}
for c in range(18):
    base = rnd(random.randint(500, 1600))
    for a in range(random.randint(6, 10)):
        s = mutate(base, random.choice([0.0, 0.005, 0.01, 0.02, 0.04]))
        if a == 1 and c % 4 == 0:
            s = tx[f"C{c:02d}_0"]
        if a == 2 and c % 5 == 0:
            p = random.randrange(50, len(s) - 60)
            s = s[:p] + "N" * 5 + s[p + 5:]
        tx[f"C{c:02d}_{a}"] = s
        if a == 0 and c % 3 == 0:
            tx[f"C{c:02d}_{a}.iso"] = s[: len(s) * 2 // 3]
with open(f"{out}/tx.fa", "w") as f:
    for i, (n, s) in enumerate(tx.items()):
        f.write(f">{n} d{i}\n")
        w = s.lower() if i % 7 == 0 else s
        for j in range(0, len(w), 70):
            f.write(w[j:j + 70] + "\n")
names = list(tx)
q = "I" * rl
with gzip.open(f"{out}/R1.fq.gz", "wt") as f1, gzip.open(f"{out}/R2.fq.gz", "wt") as f2:
    for i in range(npairs):
        kind = random.random()
        if kind < 0.04:
            r1, r2 = rnd(rl), rnd(rl)
        else:
            src = tx[random.choice(names)]
            if kind < 0.20:
                src = mutate(src, 0.05)                 # unseen allele
            L = min(len(src), random.randint(rl + 20, 450))
            p = 0 if random.random() < 0.05 else random.randint(0, len(src) - L)
            fr = src[p:p + L]
            r1, r2 = fr[:rl], rc(fr)[:rl]
            if len(r1) < rl or len(r2) < rl:
                r1, r2 = rnd(rl), rnd(rl)
            r1, r2 = mutate(r1, 0.006), mutate(r2, 0.006)
            if random.random() < 0.5:
                r1, r2 = r2, r1
        if i % 37 == 0:
            r1 = r1[:20] + "N" + r1[21:]
        if i % 53 == 0:
            r2 = r2.lower()
        f1.write(f"@p{i}/1\n{r1}\n+\n{q}\n")
        f2.write(f"@p{i}/2\n{r2}\n+\n{q}\n")
print(f"{len(tx)} transcripts, {npairs} pairs of {rl} bp -> {out}")
