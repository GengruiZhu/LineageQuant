#!/usr/bin/env python3
"""Synthetic paired-end SAM with fragment GC bias, for lqcore 1.1.0 --eff-len fld / --gc-bias.

  make_gc_case.py OUTDIR SEED [N_FRAG] [BIAS]

Transcripts: 300 'genes', lengths 300-4000 bp, each with GC-rich and AT-rich stretches (so fragment GC varies
within and between transcripts); a quarter have a second isoform (shared prefix + own tail) -> multi-mapping.
Fragments: transcript ~ abundance x (length - fragment + 1), fragment length ~ N(260, 40) truncated to [80, 600],
start uniform, kept with probability w(GC) (BIAS = strong: exp(-((gc - 0.45) / 0.10)^2); weak; none).
Each fragment is written for every transcript that contains its sequence exactly (best first = primary).
Outputs: ref.fa, ancestry.tsv (empty), case.sam (name-grouped), truth.tsv (transcript, true TPM).
"""
import math
import os
import random
import sys

out, seed = sys.argv[1], int(sys.argv[2])
nf = int(sys.argv[3]) if len(sys.argv) > 3 else 60000
bias = sys.argv[4] if len(sys.argv) > 4 else "strong"
R = random.Random(seed)
os.makedirs(out, exist_ok=True)


def chunk(n, gc):
    return "".join(R.choice("GC") if R.random() < gc else R.choice("AT") for _ in range(n))


def rc(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


tx = {}
for g in range(300):
    L = R.randint(300, 4000)
    s, n = [], 0
    while n < L:
        m = min(L - n, R.randint(80, 400))
        s.append(chunk(m, R.choice([0.25, 0.4, 0.5, 0.6, 0.75])))
        n += m
    seq = "".join(s)
    tx[f"g{g}"] = seq
    if g % 4 == 0 and L > 800:
        cut = R.randint(L // 3, 2 * L // 3)
        tail = chunk(R.randint(200, 1200), R.choice([0.3, 0.5, 0.7]))
        tx[f"g{g}iso"] = seq[:cut] + tail
names = list(tx)
ab = {t: math.exp(R.gauss(0, 1.3)) for t in names}
for t in R.sample(names, 30):
    ab[t] = 0.0


def w(gc):
    if bias == "none":
        return 1.0
    if bias == "weak":
        return math.exp(-((gc - 0.47) / 0.25) ** 2)
    return math.exp(-((gc - 0.45) / 0.10) ** 2)


def eff_true(L):   # E[#starts] under the truncated N(260, 40) [80, 600]
    s = z = 0.0
    for k in range(80, 601):
        p = math.exp(-0.5 * ((k - 260) / 40) ** 2)
        z += p
        if k <= L:
            s += p * (L - k + 1)
    return s / z


eff = {t: eff_true(len(tx[t])) for t in names}
wt = [ab[t] * eff[t] for t in names]
with open(f"{out}/ref.fa", "w") as f:
    for t in names:
        f.write(f">{t}\n{tx[t]}\n")
with open(f"{out}/ancestry.tsv", "w") as f:
    f.write("GeneID\tAncestry\n")
tot = sum(ab[t] for t in names)
with open(f"{out}/truth.tsv", "w") as f:
    f.write("transcript\ttpm\n")
    for t in names:
        f.write(f"{t}\t{ab[t] / tot * 1e6:.6f}\n")
lines, made, tries = [], 0, 0
RL = 100
while made < nf and tries < nf * 50:
    tries += 1
    t = R.choices(names, weights=wt)[0]
    s = tx[t]
    k = int(round(R.gauss(260, 40)))
    if k < 80 or k > 600 or k > len(s):
        continue
    st = R.randint(0, len(s) - k)
    frag = s[st:st + k]
    gc = (frag.count("G") + frag.count("C")) / k
    if R.random() > w(gc):
        continue
    made += 1
    q = f"f{made}"
    m1, m2 = frag[:RL], rc(frag[-RL:])
    hits = [(t, st)]
    u = t[:-3] if t.endswith("iso") else t + "iso"     # only an isoform pair shares sequence
    if u in tx:
        p = tx[u].find(frag)
        if p >= 0:
            hits.append((u, p))
    for j, (u, p) in enumerate(hits):
        sec = 256 if j else 0
        p1, p2 = p + 1, p + k - RL + 1
        lines.append(f"{q}\t{99 | sec}\t{u}\t{p1}\t255\t{RL}M\t=\t{p2}\t{k}\t{m1}\t*\tNM:i:0")
        lines.append(f"{q}\t{147 | sec}\t{u}\t{p2}\t255\t{RL}M\t=\t{p1}\t{-k}\t{rc(m2)}\t*\tNM:i:0")
with open(f"{out}/case.sam", "w") as f:
    f.write("@HD\tVN:1.6\tGO:query\n")
    for t in names:
        f.write(f"@SQ\tSN:{t}\tLN:{len(tx[t])}\n")
    f.write("\n".join(lines) + "\n")
print(f"{out}: {len(names)} transcripts, {made} fragments ({bias} GC bias)")
