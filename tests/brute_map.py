#!/usr/bin/env python3
"""Exhaustive ungapped aligner with lqmap's rules (independent implementation, numpy), for testing.
Per mate: every position of every transcript on both strands; mismatches between two ACGT bases are scored,
an ambiguous base on either side is not scored but counts in NM; per mate scored mismatches <= floor(0.1 L) and
ambiguous positions <= floor(0.5 L). Pairs: same transcript, opposite strands, forward mate starts at or before the
reverse mate, fragment = max(end) - forward start in [1, 1000]. Per transcript the best pair (lowest NM sum);
fragment keeps transcripts within best + 3. Output: name, then transcript:nm ... (sorted)."""
import gzip, sys
import numpy as np

fa, r1p, r2p = sys.argv[1:4]
EXTRA, MAXF = 3, 1000
code = np.full(256, 4, dtype=np.int8)
for i, c in enumerate("ACGT"):
    code[ord(c)] = i
    code[ord(c.lower())] = i


def enc(s):
    return code[np.frombuffer(s.encode(), dtype=np.uint8)]


names, seqs = [], []
cur = []
for l in open(fa):
    l = l.strip()
    if l.startswith(">"):
        if names: seqs.append("".join(cur))
        names.append(l[1:].split()[0]); cur = []
    else:
        cur.append(l.replace(" ", ""))
seqs.append("".join(cur))
ref = [enc(s) for s in seqs]
R = np.concatenate(ref)
starts = np.cumsum([0] + [len(r) for r in ref])
_valid = {}


def valid_mask(L):
    """positions p where the window [p, p+L) lies inside one transcript"""
    if L not in _valid:
        n = len(R) - L + 1
        p = np.arange(n)
        t0 = np.searchsorted(starts, p, side="right") - 1
        t1 = np.searchsorted(starts, p + L - 1, side="right") - 1
        _valid[L] = (t0 == t1, t0)
    return _valid[L]


def hits(q):
    """all (t, start, mm, nm) with scored mm <= floor(0.1 L) and ambiguous <= floor(0.5 L), whole transcriptome at once"""
    L = len(q)
    maxmm, nceil = int(0.1 * L + 1e-9), int(0.5 * L)
    n = len(R) - L + 1
    mm = np.zeros(n, dtype=np.int16); nm = np.zeros(n, dtype=np.int16); amb = np.zeros(n, dtype=np.int16)
    for i in range(L):
        w = R[i:i + n]
        ab = (w == 4) | (q[i] == 4)
        ne = w != q[i]
        amb += ab; nm += (ne | ab); mm += (ne & ~ab)
    ok, tid = valid_mask(L)
    sel = np.nonzero(ok & (mm <= maxmm) & (amb <= nceil))[0]
    return [(int(tid[p]), int(p - starts[tid[p]]), int(mm[p]), int(nm[p])) for p in sel]


def rd(p):
    f = gzip.open(p, "rt")
    while True:
        h = f.readline()
        if not h: return
        s = f.readline().strip(); f.readline(); f.readline()
        n = h[1:].split()[0]
        if n.endswith("/1") or n.endswith("/2"): n = n[:-2]
        yield n, s


rcmap = np.array([3, 2, 1, 0, 4], dtype=np.int8)
for (n, s1), (_, s2) in zip(rd(r1p), rd(r2p)):
    q1, q2 = enc(s1), enc(s2)
    H = []
    for q in (q1, q2):
        f = hits(q); r = hits(rcmap[q[::-1]])
        H.append([(t, s, mm, nm, 0) for t, s, mm, nm in f] + [(t, s, mm, nm, 1) for t, s, mm, nm in r])
    L = (len(q1), len(q2))
    best_t = {}
    for a in H[0]:
        for b in H[1]:
            if a[0] != b[0] or a[4] == b[4]: continue
            if a[4] == 0: f, fl, r, rl = a, L[0], b, L[1]
            else: f, fl, r, rl = b, L[1], a, L[0]
            if f[1] > r[1]: continue
            fr = max(f[1] + fl, r[1] + rl) - f[1]
            if fr < 1 or fr > MAXF: continue
            v = a[3] + b[3]
            if a[0] not in best_t or v < best_t[a[0]]: best_t[a[0]] = v
    if best_t:
        b = min(best_t.values())
        keep = sorted(f"{names[t]}:{v}" for t, v in best_t.items() if v <= b + EXTRA)
    else:
        keep = []
    print(n + "\t" + " ".join(keep))
