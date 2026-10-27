#!/usr/bin/env python3
"""Synthetic allopolyploid case for comparing `lqcore run` with snipe_kmer / snipe_snp / merge_poe / em_global.

Clusters of copies from three lineages (SOF/SRO/SSP) share a base sequence; each lineage has its own alleles
at diagnostic sites (the MSA SNP DB) and each copy has a few private mutations. Lineage-specific 31-mers are
written as markers (canonical by default). Read pairs (2x100, fragment 180-400) are sampled from copies with
sequencing errors and "aligned" to every copy of their cluster within 8 mismatches (no indels; the copies are
colinear), plus unmapped pairs. Output: ref.fa, case.sam (name-grouped), theta.tsv, ancestry.tsv,
snp_db.pkl, markers/{SOF,SRO,SSP}.fa, ids_all.txt, ids_one.txt (one copy per cluster + orphans).

usage: make_lineage_case.py OUTDIR SEED [N_FRAGMENTS] [--noncanonical]
"""
import os, pickle, random, sys

out, seed = sys.argv[1], int(sys.argv[2])
NF = int(sys.argv[3]) if len(sys.argv) > 3 and not sys.argv[3].startswith("-") else 6000
CANON = "--noncanonical" not in sys.argv
random.seed(seed)
os.makedirs(f"{out}/markers", exist_ok=True)
B = "ACGT"
LIN = ["SOF", "SRO", "SSP"]
K = 31
COMP = str.maketrans("ACGTN", "TGCAN")


def rc(s): return s.translate(COMP)[::-1]
def rseq(n): return "".join(random.choice(B) for _ in range(n))


tx = []            # (name, seq, cluster or None, lineage)
db = {}
clusters = {}
gi = 0
for c in range(25):
    L = random.randint(400, 1500)
    base = list(rseq(L))
    sites = sorted(random.sample(range(L), max(3, L // 40)))
    allele = {s: {l: random.choice(B) for l in LIN} for s in sites}
    cid = f"Global_Cluster_{c:05d}"
    mem = []
    for k in range(random.randint(2, 6)):
        gi += 1
        lin = LIN[k % 3] if k < 3 else random.choice(LIN)
        s = base[:]
        for p in sites: s[p] = allele[p][lin]
        for p in random.sample(range(L), random.randint(0, 4)): s[p] = random.choice(B)
        name = f"POJ{gi:06d}" if random.random() < 0.8 else f"POJ{gi:06d}.1"
        tx.append((name, "".join(s), cid, lin))
        mem.append(name.split(".")[0])
        for p in sites:
            if random.random() < 0.9:
                db[(name, p)] = {"bases": dict(allele[p])}
    clusters[cid] = mem
for o in range(30):
    gi += 1
    tx.append((f"POJ{gi:06d}", rseq(random.randint(300, 1500)), None, random.choice(LIN)))

with open(f"{out}/ref.fa", "w") as f:
    for n, s, _, _ in tx:
        f.write(f">{n}\n{s}\n")
with open(f"{out}/theta.tsv", "w") as f:
    f.write("Cluster_ID\tAllele_ID\tUnique_Probes\tRaw_Hits\tLR_Count\tTheta\tPrior_Status\tShrinkage_Used\n")
    for cid, mem in clusters.items():
        for g in mem:
            f.write(f"{cid}\t{g}\t3\t30\t0\t{random.random():.6f}\tKmer_Inferred\tNA\n")
with open(f"{out}/ancestry.tsv", "w") as f:
    f.write("GeneID\tAncestry\n")
    for n, s, c, l in tx:
        f.write(f"{n.split('.')[0]}\t{l}\n")
with open(f"{out}/snp_db.pkl", "wb") as f:
    pickle.dump(db, f)

# lineage-specific k-mers (present in copies of one lineage only)
occ = {}
for n, s, c, l in tx:
    for i in range(len(s) - K + 1):
        km = s[i:i + K]
        key = min(km, rc(km))
        occ.setdefault(key, set()).add(l)
for l in LIN:
    with open(f"{out}/markers/{l}.fa", "w") as f:
        j = 0
        for key, ls in occ.items():
            if ls == {l} and random.random() < 0.7:
                km = key if CANON else (key if random.random() < 0.5 else rc(key))
                f.write(f">{j}\n{km}\n"); j += 1
                if not CANON and random.random() < 0.5:
                    f.write(f">{j}\n{rc(km)}\n"); j += 1

with open(f"{out}/ids_all.txt", "w") as f:
    for n, *_ in tx: f.write(n.split(".")[0] + "\n")
with open(f"{out}/ids_one.txt", "w") as f:
    seen = set()
    for n, s, c, l in tx:
        if c is None or c not in seen:
            f.write(n.split(".")[0] + "\n"); seen.add(c)

# reads
by_c = {}
for t in tx: by_c.setdefault(t[2] or t[0], []).append(t)
recs = []
def mm(a, b): return sum(x != y for x, y in zip(a, b))
for i in range(NF):
    q = f"frag{i:07d}"
    if random.random() < 0.04:
        recs.append(f"{q}\t77\t*\t0\t0\t*\t*\t0\t0\t{rseq(100)}\t{'I'*100}\tYT:Z:UP")
        recs.append(f"{q}\t141\t*\t0\t0\t*\t*\t0\t0\t{rseq(100)}\t{'I'*100}\tYT:Z:UP")
        continue
    src = random.choice(tx)
    s = src[1]
    fl = random.randint(180, min(400, len(s)))
    p = random.randint(0, len(s) - fl)
    frag = list(s[p:p + fl])
    for j in range(fl):
        if random.random() < 0.006: frag[j] = random.choice(B)
    if random.random() < 0.02: frag[random.randrange(fl)] = "N"
    frag = "".join(frag)
    m1, m2 = frag[:100], frag[-100:]          # mate2 as it appears on the forward strand
    cands = [t for t in by_c[src[2] or src[0]] if len(t[1]) >= p + fl]
    hits = []
    for t in cands:
        n1, n2 = mm(m1, t[1][p:p + 100]), mm(m2, t[1][p + fl - 100:p + fl])
        if n1 <= 8 and n2 <= 8: hits.append((n1 + n2, t[0], n1, n2))
    if not hits:
        recs.append(f"{q}\t77\t*\t0\t0\t*\t*\t0\t0\t{m1}\t{'I'*100}\tYT:Z:UP")
        recs.append(f"{q}\t141\t*\t0\t0\t*\t*\t0\t0\t{rc(m2)}\t{'I'*100}\tYT:Z:UP")
        continue
    hits.sort()
    rev = random.random() < 0.5
    for j, (tot, t, n1, n2) in enumerate(hits):
        sec = 256 if j else 0
        if not rev:   # mate1 forward at p, mate2 reverse at p+fl-100
            recs.append(f"{q}\t{99 | sec}\t{t}\t{p + 1}\t255\t100M\t=\t{p + fl - 99}\t{fl}\t{m1}\t{'I'*100}\tAS:i:{-n1}\tNM:i:{n1}\tYT:Z:CP")
            recs.append(f"{q}\t{147 | sec}\t{t}\t{p + fl - 99}\t255\t100M\t=\t{p + 1}\t{-fl}\t{m2}\t{'I'*100}\tAS:i:{-n2}\tNM:i:{n2}\tYT:Z:CP")
        else:         # mate1 is the reverse read
            recs.append(f"{q}\t{83 | sec}\t{t}\t{p + fl - 99}\t255\t100M\t=\t{p + 1}\t{-fl}\t{m2}\t{'I'*100}\tAS:i:{-n2}\tNM:i:{n2}\tYT:Z:CP")
            recs.append(f"{q}\t{163 | sec}\t{t}\t{p + 1}\t255\t100M\t=\t{p + fl - 99}\t{fl}\t{m1}\t{'I'*100}\tAS:i:{-n1}\tNM:i:{n1}\tYT:Z:CP")
with open(f"{out}/case.sam", "w") as f:
    f.write("@HD\tVN:1.0\tSO:unsorted\tGO:query\n")
    for n, s, *_ in tx: f.write(f"@SQ\tSN:{n}\tLN:{len(s)}\n")
    f.write("\n".join(recs) + "\n")
print(f"{out}: {len(tx)} transcripts, {len(clusters)} clusters, {len(db)} SNPs, {NF} fragments, markers {'canonical' if CANON else 'both strands'}")
