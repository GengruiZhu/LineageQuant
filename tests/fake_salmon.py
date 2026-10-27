#!/usr/bin/env python3
"""Stand-in for `salmon quant ... --writeMappings=FILE` in tests: takes the candidates of `lqcore map` (env
FAKE_SALMON_LQ = lqcore binary, FAKE_SALMON_FASTA = transcriptome) and writes them the way salmon does: no NM tag,
SEQ '*' on secondary records, plus noise that LineageQuant must discard (positions shifted by 2, orphan mates,
records on an unknown reference). Writes quant.sf-like output into -o. Other options are ignored."""
import os, random, subprocess, sys

a = sys.argv[1:]
r1, r2 = a[a.index("-1") + 1], a[a.index("-2") + 1]
wm = [x for x in a if x.startswith("--writeMappings=")][0].split("=", 1)[1]
od = a[a.index("-o") + 1]
os.makedirs(od, exist_ok=True)
lq, fa = os.environ["FAKE_SALMON_LQ"], os.environ["FAKE_SALMON_FASTA"]
sam = subprocess.run([lq, "map", "--fasta", fa, "--r1", r1, "--r2", r2, "--threads", "2"],
                     stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, check=True, text=True).stdout.splitlines()
random.seed(3)
with open(wm, "w") as f:
    cur, grp = None, []

    def flush():
        if not grp: return
        out = []
        for l in grp:
            q = l.split("\t")
            if int(q[1]) & 4: continue                       # salmon writes no unmapped records here
            q = q[:11]                                       # drop NM
            if int(q[1]) & 256: q[9] = "*"
            out.append("\t".join(q))
        if out and random.random() < 0.3:                    # a shifted (gapped-looking) candidate
            q = out[0].split("\t"); q[3] = str(int(q[3]) + 2); q[1] = str(int(q[1]) | 256); q[9] = "*"
            out.append("\t".join(q))
        if out and random.random() < 0.2:                    # an orphan mate
            q = out[0].split("\t"); q[1] = str((int(q[1]) | 256) & ~2); q[2] = "UNKNOWN_TX" if random.random() < 0.5 else q[2]
            q[3] = str(max(1, int(q[3]) + 300)); q[9] = "*"
            out.append("\t".join(q))
        for l in out: f.write(l + "\n")

    for l in sam:
        if l.startswith("@"):
            if l.startswith("@SQ") or l.startswith("@HD"): f.write(l + "\n")
            continue
        n = l.split("\t", 1)[0]
        if n != cur:
            flush(); cur, grp = n, []
        grp.append(l)
    flush()
open(os.path.join(od, "quant.sf"), "w").write("Name\tLength\tEffectiveLength\tTPM\tNumReads\n")
print("[fake salmon] done", file=sys.stderr)
