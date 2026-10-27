#!/usr/bin/env python3
"""Random long-read alignment SAM for lqcore lrcount vs extract_prior.py.

  make_lr_case.py SEED OUT.sam [pb|sra] [grouped|coord]

Reads have 1-4 segments (primary + supplementary, hard or soft clipped, either strand), secondaries on other
copies with tied / lower AS (some overlapping no segment), missing AS/NM tags, QC-fail / duplicate records,
unmapped reads; pb names put 1-5 subreads in a ZMW, sra names carry no ZMW.
"""
import random
import sys

seed, out = int(sys.argv[1]), sys.argv[2]
names = sys.argv[3] if len(sys.argv) > 3 else "pb"
order = sys.argv[4] if len(sys.argv) > 4 else "grouped"
R = random.Random(seed)
genes = [f"g{i:03d}" for i in range(40)]
refs = []
for g in genes:
    for iso in range(R.choice([1, 1, 2, 3])):
        refs.append(f"{g}.{iso + 1}")
lines = []
recs = []
zmw = 0
nread = 0
while nread < 600:
    zmw += 1
    nsub = R.choice([1, 1, 2, 3, 5]) if names == "pb" else 1
    s0 = 0
    for k in range(nsub):
        nread += 1
        L = R.randint(300, 6000)
        q = f"mv1/{zmw}/{s0}_{s0 + L}" if names == "pb" else f"SRR1.{nread}"
        s0 += L + 40
        if R.random() < 0.05:
            recs.append([q, "4", "*", "0", "0", "*", "*", "0", "0", "*", "*"])
            continue
        nseg = R.choice([1, 1, 1, 2, 3, 4])
        cuts = sorted(R.sample(range(50, L - 50), nseg - 1)) if nseg > 1 else []
        bounds = [0] + cuts + [L]
        first = True
        for si in range(nseg):
            a, b = bounds[si], bounds[si + 1]
            ln = b - a
            if ln < 30:
                continue
            ref = R.choice(refs)
            rev = R.random() < 0.5
            lead, trail = (L - b, a) if rev else (a, L - b)
            clip = "S" if first else R.choice("HS")
            cig = (f"{lead}{clip}" if lead else "") + f"{ln}M" + (f"{trail}{clip}" if trail else "")
            if R.random() < 0.2 and ln > 60:
                d = R.randint(1, 5)
                cig = (f"{lead}{clip}" if lead else "") + f"{ln // 2}M{d}D{ln - ln // 2 - 2}M2I" + (f"{trail}{clip}" if trail else "")
            nm = R.randint(0, max(1, ln // 4))
            AS = ln * 2 - nm * 6
            flag = (16 if rev else 0) | (0 if first else 2048)
            if R.random() < 0.03:
                flag |= R.choice([512, 1024])
            tags = []
            if R.random() > 0.05:
                tags.append(f"AS:i:{AS}")
            if R.random() > 0.05:
                tags.append(f"NM:i:{nm}")
            recs.append([q, str(flag), ref, "1", "60", cig, "*", "0", "0", "*", "*"] + tags)
            first = False
            # secondaries: ties (same AS) on other copies / isoforms, lower ones, some shifted
            for _ in range(R.choice([0, 0, 1, 2, 4])):
                r2 = R.choice(refs)
                as2 = AS if R.random() < 0.5 else AS - R.randint(1, 40)
                sh = 0 if R.random() < 0.8 else R.randint(0, L)
                l2, t2 = (min(lead + sh, L - 1), max(trail - sh, 0)) if not rev else (lead, trail)
                m2 = max(1, L - l2 - t2)
                fl2 = 256 | (16 if R.random() < 0.5 else 0)
                tg = [f"AS:i:{as2}"] if R.random() > 0.05 else []
                tg += [f"NM:i:{R.randint(0, max(1, m2 // 4))}"]
                recs.append([q, str(fl2), r2, "1", "0", (f"{l2}S" if l2 else "") + f"{m2}M" + (f"{t2}S" if t2 else ""),
                             "*", "0", "0", "*", "*"] + tg)
        if first:   # every segment was too short: unmapped read
            recs.append([q, "4", "*", "0", "0", "*", "*", "0", "0", "*", "*"])
hd = "@HD\tVN:1.6\tSO:unsorted\tGO:query" if order == "grouped" else "@HD\tVN:1.6\tSO:coordinate"
if order != "grouped":
    R.shuffle(recs)
with open(out, "w") as f:
    f.write(hd + "\n")
    for r in refs:
        f.write(f"@SQ\tSN:{r}\tLN:10000\n")
    for r in recs:
        f.write("\t".join(r) + "\n")
