#!/usr/bin/env python3
"""Simulated PacBio subreads of MAS-like arrays, for lqcore lrpick.

  make_subreads.py SEED OUT.fq TRUTH.tsv [sra|pb|hifi] [N_ZMW] [ERR] [N_GENES]

Transcripts: NG (default 3000) 'genes' with a skewed abundance (a few very abundant ones), each with 2-4 homoeolog copies
(~2 % divergent). A molecule is an array of 1-8 transcripts joined by shared adapters, with a shared
5' primer and a poly(A) tail per segment. A ZMW gives 1-6 passes alternating strand; the first and the
last pass are partial. Errors: ERR (default 0.12) split 70 % indels / 30 % substitutions.
TRUTH: read name -> molecule (ZMW) id.
"""
import random
import sys

seed, out, truth = int(sys.argv[1]), sys.argv[2], sys.argv[3]
mode = sys.argv[4] if len(sys.argv) > 4 else "sra"
nz = int(sys.argv[5]) if len(sys.argv) > 5 else 3000
err = float(sys.argv[6]) if len(sys.argv) > 6 else 0.12
R = random.Random(seed)
B = "ACGT"
rs = lambda n: "".join(R.choice(B) for _ in range(n))
comp = str.maketrans("ACGT", "TGCA")
rc = lambda s: s.translate(comp)[::-1]
adapter = [rs(20) for _ in range(8)]
primer = rs(25)
tx = []
NG = int(sys.argv[7]) if len(sys.argv) > 7 else 3000
for g in range(NG):
    base = rs(R.randint(600, 2500))
    for c in range(R.randint(2, 4)):
        s = list(base)
        for i in range(len(s)):
            if R.random() < 0.02:
                s[i] = R.choice(B)
        tx.append("".join(s))
w = [1.0 / (i + 1) ** 1.1 for i in range(len(tx))]
R.shuffle(w)


def noisy(s, e):
    o = []
    for ch in s:
        r = R.random()
        if r < e * 0.35:            # deletion
            continue
        if r < e * 0.7:             # insertion
            o.append(R.choice(B))
            o.append(ch)
            continue
        if r < e:
            o.append(R.choice(B))
            continue
        o.append(ch)
    return "".join(o)


with open(out, "w") as f, open(truth, "w") as t:
    n = 0
    for z in range(1, nz + 1):
        segs = R.choices(range(len(tx)), weights=w, k=R.randint(1, 8))
        mol = ""
        for i, k in enumerate(segs):
            seg = primer + tx[k] + "A" * 30
            if R.random() < 0.5:
                seg = rc(seg)
            mol += adapter[i % 8] + seg
        if mode == "hifi":
            passes = [mol]
            e = 0.002
        else:
            np_ = R.choice([1, 1, 2, 2, 3, 4, 6])
            passes = []
            for p in range(np_):
                s = mol if p % 2 == 0 else rc(mol)
                if p == 0 and np_ > 1:
                    s = s[R.randint(len(s) // 3, len(s) - 50):]
                if p == np_ - 1 and np_ > 1:
                    s = s[: R.randint(50, max(51, 2 * len(s) // 3))]
                passes.append(s)
            e = err
        s0 = 0
        for p, s in enumerate(passes):
            n += 1
            seq = noisy(s, e)
            nm = f"SRR9.{n} {n}/1" if mode in ("sra", "hifi") else f"mv/{z}/{s0}_{s0 + len(seq)}"
            s0 += len(seq) + 45
            q = "~" * len(seq) if mode == "hifi" else "!" * len(seq)
            f.write(f"@{nm}\n{seq}\n+\n{q}\n")
            t.write(f"{nm.split()[0]}\t{z}\n")
