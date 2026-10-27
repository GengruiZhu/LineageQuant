#!/usr/bin/env python3
"""Stand-in for bowtie2 in tests: reads -1/-2 (gzip, plain or FIFO) pair by pair and prints, for each pair,
the records of that read name from a library SAM (env FAKE_BT2_SAM), in input order. Writes a bowtie2-like
summary to stderr. Options it does not need are ignored."""
import gzip, os, sys

args = sys.argv[1:]
r1 = args[args.index("-1") + 1]
r2 = args[args.index("-2") + 1]
lib = os.environ["FAKE_BT2_SAM"]
head, recs = [], {}
with open(lib) as f:
    for l in f:
        if l.startswith("@"):
            head.append(l)
        else:
            recs.setdefault(l.split("\t", 1)[0], []).append(l)


def op(p):
    f = open(p, "rb")
    m = f.peek(2)[:2] if hasattr(f, "peek") else b""
    if m == b"\x1f\x8b":
        return gzip.open(f, "rt")
    return open(p, "rt") if False else __import__("io").TextIOWrapper(f)


out = sys.stdout
out.write("".join(head))
n = al = 0
with op(r1) as a, op(r2) as b:
    while True:
        h1 = a.readline()
        if not h1:
            break
        a.readline(); a.readline(); a.readline()
        h2 = b.readline(); b.readline(); b.readline(); b.readline()
        q1, q2 = h1[1:].split()[0], h2[1:].split()[0]
        if q1.endswith("/1"): q1, q2 = q1[:-2], q2[:-2]
        assert q1 == q2, (q1, q2)
        n += 1
        rs = recs.get(q1, [])
        if rs and not int(rs[0].split("\t")[1]) & 4: al += 1
        out.write("".join(rs))
sys.stderr.write(f"{n} reads; of these:\n  {n} (100.00%) were paired; of these:\n"
                 f"    {n - al} aligned concordantly 0 times\n{100.0 * al / max(n, 1):.2f}% overall alignment rate\n")
