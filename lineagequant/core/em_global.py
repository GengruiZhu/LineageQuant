#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
core/em_global.py -- global EM engine (Step 03), default since LineageQuant 1.0.0

  [G1] All genes (allele copies in clusters and single-copy "orphan" genes) enter one
       EM over read equivalence classes. EM factorises over the connected components
       of the read-compatibility graph, so this equals running EM per component.
       Every fragment carries total weight 1: no double counting of reads shared
       between clusters, no idxstats counting of multi-mapping orphan reads, and no
       genes missing from the output because they are absent from orphans.tsv.
  [G2] Effective-length correction (default): E-step weight = theta / (L - mean
       fragment length + 1), with L the longest transcript of the gene.
  [G3] Genes whose (single) transcripts are identical in sequence are pooled into one
       target (needs --fasta) and split evenly after convergence. Their read evidence
       is identical, so without pooling a handful of asymmetrically reported reads
       would be amplified by the EM into 0/1 splits.
  [G4] Mismatch-count weighting (--nm-weight, default auto): when the candidate
       alignments of a fragment differ in mismatch count (sum of NM over both mates),
       each extra mismatch multiplies the likelihood by lambda = (p/3)/(1-p), with p
       the sequencing error rate estimated from the best alignments. This is the
       likelihood model behind Salmon's --scoreExp and RSEM's error model, and it is
       what makes bowtie2 -k 200 (which also reports near-optimal alignments) usable.
       With HISAT2 input it has essentially no effect (HISAT2 reports equal-score
       candidates only).
  [G5] --fld-weight (experimental, off): additionally weight candidates by the
       empirical fragment-length distribution. No measurable gain in benchmarks.

  Kept from the legacy engine: theta0 from Global_Theta_Matrix as initial values,
  the Dirichlet prior tau for LongRead_Anchored clusters, and the lineage
  probabilities (-p) with the same --p-weight-* / --p-alpha / --p-floor options.

  Output has the same columns as the legacy EM_Final_Expression.tsv; single-copy
  genes have ClusterID = "ORPHAN". The BAM may be coordinate-sorted (it is collated
  by read name into a temporary file first) or already grouped by read name.
"""
import os, re, sys, time, argparse, itertools, subprocess, tempfile, hashlib
from collections import Counter, defaultdict
import numpy as np

ap = argparse.ArgumentParser(description="LineageQuant global EM engine (Step 03)")
ap.add_argument("-t", "--theta", default=None)
ap.add_argument("-p", "--pvals", default=None)
ap.add_argument("-a", "--ancestry", required=True)
ap.add_argument("-o", "--out", default="EM_Final_Expression.tsv")
ap.add_argument("--bam", required=True)
ap.add_argument("--fasta", default=None, help="transcriptome FASTA; if given, identical-sequence genes are pooled (G3)")
ap.add_argument("--no-pool", action="store_true")
ap.add_argument("--threads", type=int, default=16)
ap.add_argument("--eff-len", choices=["off", "len", "eff"], default="eff")
ap.add_argument("--frag-len", type=float, default=0.0)
ap.add_argument("--min-eff-len", type=float, default=25.0)
ap.add_argument("--max-iter", type=int, default=5000)
ap.add_argument("--min-iter", type=int, default=50)
ap.add_argument("--tol", type=float, default=1e-3)
ap.add_argument("--init", choices=["theta", "uniform"], default="theta",
                help="initial abundances: theta (default; Theta-weighted within clusters) or uniform "
                     "(multi-mapping fragments split equally) -- for checking that results do not depend on the start")
ap.add_argument("--p-weight-specific", type=float, default=1.0)
ap.add_argument("--p-weight-conflict", type=float, default=0.0)
ap.add_argument("--p-weight-none", type=float, default=0.0)
ap.add_argument("--p-alpha", type=float, default=0.0)
ap.add_argument("--p-floor", type=float, default=1e-3)
ap.add_argument("--theta-prior-lr", type=float, default=1.0)
ap.add_argument("--theta-prior-kmer", type=float, default=0.0)
ap.add_argument("--theta-prior-cap", type=float, default=0.0)
ap.add_argument("--keep-collate", action="store_true")
ap.add_argument("--nm-weight", default="auto",
                help="auto (default) | off | lambda: likelihood ratio per extra mismatch (G4)")
ap.add_argument("--fld-weight", action="store_true", help="experimental (G5)")
a, unknown = ap.parse_known_args()
T0 = time.time()
def log(m): print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {m}", flush=True)
if unknown: log(f"[INFO] ignoring legacy-engine options: {' '.join(unknown)}")
def strip_iso(g): return re.sub(r'\.t?\d+$', '', g.split()[0])
log(f"[INFO] EM (global) | eff-len={a.eff_len} | identical-sequence pooling={'on' if (a.fasta and not a.no_pool) else 'off'} | "
    f"nm-weight={a.nm_weight} | init={a.init} | "
    f"tol={a.tol} max_iter={a.max_iter}")

# ---------------------------------------------------------------- Theta / ancestry
g2c, theta0, lr_count, prior_status = {}, {}, defaultdict(float), {}
if a.theta and os.path.exists(a.theta):
    with open(a.theta) as f:
        f.readline()
        for l in f:
            p = l.rstrip("\n").split("\t")
            if len(p) < 6: continue
            g2c[p[1]] = p[0]
            try: theta0[p[1]] = float(p[5])
            except ValueError: theta0[p[1]] = 0.0
            try: lr_count[p[0]] += float(p[4])
            except (ValueError, IndexError): pass
            if len(p) > 6: prior_status[p[0]] = p[6]
anc = {}
if os.path.exists(a.ancestry):
    for l in open(a.ancestry):
        p = l.split()
        if len(p) >= 2 and not l.startswith("GeneID"):
            anc[p[0]] = p[1].upper(); anc.setdefault(p[0].rsplit(".", 1)[0], p[1].upper())
log(f"[INFO] Theta: {len(set(g2c.values())):,} clusters / {len(g2c):,} allele copies")

# ---------------------------------------------------------------- input BAM (collated by read name if needed)
def reader(path):
    if path.endswith(".sam"):                         # plain SAM (tests)
        refs, lens, body = [], [], []
        fh = open(path)
        for l in fh:
            if l.startswith("@SQ"):
                d = dict(x.split(":", 1) for x in l.rstrip("\n").split("\t")[1:])
                refs.append(d["SN"]); lens.append(int(d["LN"]))
            elif not l.startswith("@"):
                body.append(l); break
        yield refs, lens
        r2i = {r: i for i, r in enumerate(refs)}
        for l in itertools.chain(body, fh):
            q = l.rstrip("\n").split("\t"); fl = int(q[1])
            nm = next((int(x[5:]) for x in q[11:] if x.startswith("NM:i:")), 0)
            yield (q[0], (-1 if fl & 4 else r2i[q[2]]), abs(int(q[8])), bool(fl & 2) and not fl & 0x900,
                   nm, 2 if fl & 0x80 else 1, (len(q[9]) if q[9] != "*" and not fl & 0x900 else 0))
        return
    import pysam
    bam = pysam.AlignmentFile(path, "rb", check_sq=False)
    yield list(bam.references), list(bam.lengths)
    for r in bam.fetch(until_eof=True):
        prim = not r.is_secondary and not r.is_supplementary
        yield (r.query_name, -1 if r.is_unmapped else r.reference_id, abs(r.template_length),
               r.is_proper_pair and prim, (r.get_tag("NM") if r.has_tag("NM") else 0),
               2 if r.is_read2 else 1, (r.query_length if prim else 0))
    bam.close()

bam_in, tmp_collate = a.bam, None
if not a.bam.endswith(".sam"):
    import pysam
    hd = pysam.AlignmentFile(a.bam, "rb", check_sq=False).header.to_dict().get("HD", {})
    if hd.get("SO") != "queryname" and hd.get("GO") != "query":
        tmpdir = os.path.dirname(os.path.abspath(a.out)) or "."
        tmp_collate = os.path.join(tmpdir, ".v6_collate.bam")
        t = time.time()
        subprocess.run(["samtools", "collate", "-@", str(a.threads), "-o", tmp_collate, a.bam,
                        os.path.join(tmpdir, ".v6_collate_tmp")], check=True)
        bam_in = tmp_collate
        log(f"[INFO] samtools collate done ({time.time()-t:.0f}s)")

rd = reader(bam_in)
refs, lens = next(rd)
tkey, names, tlen, tcl = {}, [], [], []
ref2t = []
for ref, ln in zip(refs, lens):
    r1 = ref.rsplit(".", 1)[0]
    k = r1 if r1 in g2c else (ref if ref in g2c else None)
    key = ("A", k) if k is not None else ("O", strip_iso(ref))
    if key not in tkey:
        tkey[key] = len(names); names.append(k if k is not None else ref)
        tlen.append(ln); tcl.append(g2c[k] if k is not None else None)
    else:
        i = tkey[key]; tlen[i] = max(tlen[i], ln)
    ref2t.append(tkey[key])
for g in g2c:
    if ("A", g) not in tkey:
        tkey[("A", g)] = len(names); names.append(g); tlen.append(1000); tcl.append(g2c[g])
T = len(names)
is_al = np.array([c is not None for c in tcl])
cl_ids = {c: i for i, c in enumerate(sorted({c for c in tcl if c is not None}))}
cl_int = np.array([cl_ids[c] if c is not None else -1 for c in tcl], dtype=np.int64)
log(f"[INFO] {len(refs):,} reference sequences -> {T:,} genes ({int(is_al.sum()):,} allele copies / {int((~is_al).sum()):,} orphans)")

# ---------------------------------------------------------------- [G3] identical-sequence pooling
POOL = {}
if a.fasta and not a.no_pool and os.path.exists(a.fasta):
    md5 = {}; nm = None; buf = []
    for l in open(a.fasta):
        if l.startswith(">"):
            if nm: md5[nm] = hashlib.md5("".join(buf).upper().encode()).hexdigest()
            nm = l[1:].split()[0]; buf = []
        else: buf.append(l.strip())
    if nm: md5[nm] = hashlib.md5("".join(buf).upper().encode()).hexdigest()
    t_refs = defaultdict(list)
    for ref, t in zip(refs, ref2t): t_refs[t].append(ref)
    by = defaultdict(list)
    for t, rs in t_refs.items():
        if len(rs) == 1 and rs[0] in md5: by[md5[rs[0]]].append(t)
    for ts_ in by.values():
        if len(ts_) > 1:
            for t in ts_: POOL[t] = min(ts_)
    ref2t = [POOL.get(t, t) for t in ref2t]
    log(f"[INFO] identical-sequence genes pooled: {len(POOL):,} -> {len(set(POOL.values())):,} groups (split evenly after convergence)")
POOLG = defaultdict(list)
for m_, r_ in POOL.items(): POOLG[r_].append(m_)

# ---------------------------------------------------------------- lineage probabilities (optional)
lineages, Pread, WC = [], {}, {"specific": a.p_weight_specific, "conflict": a.p_weight_conflict, "none": a.p_weight_none}
_SPEC = re.compile(r"^specific", re.I); _NONE = re.compile(r"none|no_snp|chimeric|ambiguous|weak", re.I)
def status_class(s):
    if _SPEC.match(s): return "specific"
    if _NONE.search(s): return "none"
    return "conflict"
if a.pvals and os.path.exists(a.pvals):
    with open(a.pvals) as fh:
        cols = fh.readline().rstrip("\n").split("\t")
        pc = [(m.group(1).upper(), i) for i, c in enumerate(cols) for m in [re.fullmatch(r"P\((.+)\)", c.strip())] if m]
        si = max([i for i, c in enumerate(cols) if "status" in c.lower()], default=None)
        lineages = [l for l, _ in pc]
        for l in fh:
            q = l.rstrip("\n").split("\t")
            if q[0] in Pread: continue
            wc = status_class(q[si]) if si is not None and len(q) > si else "conflict"
            w = WC[wc]
            if w == 0.0: continue
            try: pv = np.array([float(q[i]) for _, i in pc])
            except (ValueError, IndexError): continue
            if a.p_alpha > 0: pv = (pv + a.p_alpha) / (1.0 + len(pc) * a.p_alpha)
            pv = np.maximum(pv, a.p_floor)
            if w != 1.0: pv = pv ** w
            Pread[q[0]] = pv
    log(f"[INFO] lineage P: {lineages}; reads with a non-zero weight: {len(Pread):,}")
L = len(lineages); lin_rank = {l: k for k, l in enumerate(lineages)}
tlin = np.array([lin_rank.get(anc.get(n, anc.get(n.rsplit(".", 1)[0], "?")), -1) if is_al[i] else -1
                 for i, n in enumerate(names)])
def factors(ts, pv):
    out = []
    for t in ts:
        mem = POOLG.get(t, [t])
        out.append(float(np.mean([pv[tlin[m]] if tlin[m] >= 0 else 1.0 / L for m in mem])))
    return tuple(round(x, 4) for x in out)

# ---------------------------------------------------------------- scan -> equivalence classes
single = [0] * T
eq = Counter(); eqx = Counter()          # eqx: weighted classes, key = (targets, lineage factors|None, extra mismatches|None, fragment-length bins|None)
fl_vals = []
n_rec = n_frag = n_unal = 0
NMW = a.nm_weight != "off"; FLW = a.fld_weight
INF = 1 << 30
nm_sum = 0; len_sum = 0; n_wcls = 0
cur, ts, td, qlen = None, set(), {}, 0
def flush():
    global n_frag, n_unal, nm_sum, len_sum
    n_frag += 1
    if not ts: n_unal += 1; return
    dn = tb = None
    if (NMW or FLW) and td:
        tot = {t: (e[0] if e[0] < INF else 0) + (e[1] if e[1] < INF else 0) for t, e in td.items()}
        best = min(tot.values())
        if qlen: nm_sum += best; len_sum += qlen
    if len(ts) == 1: single[next(iter(ts))] += 1; return
    k = tuple(sorted(ts))
    if NMW:
        d_ = tuple(tot[t] - best for t in k)
        if any(d_): dn = d_
    if FLW:
        b_ = tuple(td[t][2] // 10 for t in k)
        if len(set(b_)) > 1: tb = b_
    pf = factors(k, Pread[cur]) if (Pread and cur in Pread) else None
    if pf is None and dn is None and tb is None: eq[k] += 1
    else: eqx[(k, pf, dn, tb)] += 1
for qn, rid, tl, pp, nm, mate, ql in rd:
    n_rec += 1
    if qn != cur:
        if cur is not None: flush()
        cur, ts, td, qlen = qn, set(), {}, 0
    qlen += ql
    if rid >= 0:
        t = ref2t[rid]; ts.add(t)
        if NMW or FLW:
            e = td.get(t)
            if e is None: e = td[t] = [INF, INF, 0]
            if mate == 1: e[0] = min(e[0], nm)
            else: e[1] = min(e[1], nm)
            if tl and not e[2]: e[2] = tl
    if pp and 0 < tl < 2000 and len(fl_vals) < 1_000_000: fl_vals.append(tl)
if cur is not None: flush()
if tmp_collate and not a.keep_collate:
    try: os.remove(tmp_collate)
    except OSError: pass
single = np.array(single, dtype=np.float64)
N_al = n_frag - n_unal
log(f"[INFO] {n_rec:,} records -> {n_frag:,} fragments ({n_unal:,} unaligned); multi-target classes {len(eq):,} + weighted {len(eqx):,}")
lam = 1.0
if NMW:
    if a.nm_weight == "auto":
        p_err = min(max(nm_sum / max(len_sum, 1), 1e-4), 0.1)
        lam = (p_err / 3.0) / (1.0 - p_err)
        log(f"[INFO] [G4] error rate p = {p_err:.5f} (best-alignment NM / read length) -> x {lam:.2e} per extra mismatch")
    else:
        lam = float(a.nm_weight); log(f"[INFO] [G4] x {lam:.2e} per extra mismatch")
    nd = sum(c for (k, pf, dn, tb), c in eqx.items() if dn is not None)
    log(f"[INFO] [G4] multi-mapping fragments whose candidates differ in mismatch count: {nd:,}")
fld = None
if FLW and len(fl_vals) >= 1000:
    h = np.bincount(np.array(fl_vals) // 10, minlength=200).astype(float) + 1.0
    fld = h / h.max()
    log(f"[INFO] [G5] fragment-length distribution: {len(fl_vals):,} proper-pair records, mode {10*int(np.argmax(h))}-{10*int(np.argmax(h))+9} bp")
def wvec(k, pf, dn, tb):
    w = np.ones(len(k)) if pf is None else np.array(pf, dtype=float)
    if dn is not None: w = w * lam ** np.array(dn, dtype=float)
    if tb is not None and fld is not None:
        f_ = np.array([fld[min(b, len(fld) - 1)] for b in tb]); w = w * f_ / f_.max()
    return tuple(w)
keys = list(eq.keys()) + [x[0] for x in eqx.keys()]
facs = [None] * len(eq) + [wvec(*x) for x in eqx.keys()]
cnt = np.array(list(eq.values()) + list(eqx.values()), dtype=np.float64)
C = len(keys)
klen = np.array([len(k) for k in keys], dtype=np.int64)
I = np.fromiter(itertools.chain.from_iterable(keys), dtype=np.int64, count=int(klen.sum()))
K = np.repeat(np.arange(C), klen)
F = np.fromiter(itertools.chain.from_iterable(f if f is not None else (1.0,) * len(k) for k, f in zip(keys, facs)),
                dtype=np.float64, count=int(klen.sum()))
del eq, eqx, keys, facs
if C:
    ha = np.bincount(K, weights=is_al[I].astype(float), minlength=C) > 0
    ho = np.bincount(K, weights=(~is_al[I]).astype(float), minlength=C) > 0
    starts = np.concatenate([[0], np.cumsum(klen)[:-1]])
    ci = cl_int[I]
    same = np.minimum.reduceat(np.where(ci < 0, 1 << 40, ci), starts) == np.maximum.reduceat(ci, starts)
    rows = [("unique", single.sum()), ("within-cluster multi", cnt[ha & ~ho & same].sum()), ("cross-cluster multi", cnt[ha & ~ho & ~same].sum()),
            ("between orphans", cnt[~ha & ho].sum()), ("allele+orphan", cnt[ha & ho].sum())]
    log("[INFO] aligned fragments: " + "  ".join(f"{k} {100*v/max(N_al,1):.2f}%" for k, v in rows))

# ---------------------------------------------------------------- effective length
tlen = np.array(tlen, dtype=np.float64)
mu = 0.0
if a.eff_len == "eff":
    mu = a.frag_len or (float(np.mean(fl_vals)) if len(fl_vals) >= 1000 else 0.0)
    if mu <= 0: log("[WARN] fragment length could not be estimated -> using full transcript length")
    else: log(f"[INFO] mean fragment length {mu:.1f} bp")
eff = (np.maximum(tlen - mu + 1.0 if mu > 0 else tlen, a.min_eff_len) if a.eff_len != "off" else np.ones(T))

# ---------------------------------------------------------------- initial values + prior
th0 = np.array([theta0.get(n, 0.0) if is_al[i] else 1.0 for i, n in enumerate(names)])
uni = single + (np.bincount(I, weights=(cnt / klen)[K], minlength=T) if C else 0.0)
alpha = uni.copy()
prior = np.zeros(T)
if is_al.any():
    ca = cl_int[is_al]; ta = th0[is_al]
    Z = np.bincount(ca, weights=ta); S = np.bincount(ca, weights=uni[is_al])
    ok = Z[ca] > 0
    tb = uni[is_al].copy(); tb[ok] = S[ca][ok] * ta[ok] / Z[ca][ok]
    if a.init == "theta":
        alpha[is_al] = 0.5 * uni[is_al] + 0.5 * tb
    cl_name = {v: k for k, v in cl_ids.items()}
    def tau_of(cid):
        t = (a.theta_prior_lr * lr_count.get(cid, 0.0) if prior_status.get(cid, "") == "LongRead_Anchored"
             else a.theta_prior_kmer)
        return max(min(t, a.theta_prior_cap) if a.theta_prior_cap > 0 else t, 0.0)
    tau_c = np.array([tau_of(cl_name[j]) for j in range(len(cl_ids))])
    if tau_c.any():
        thn = np.where(Z[ca] > 0, ta / np.maximum(Z[ca], 1e-300), 0.0)
        prior[is_al] = tau_c[ca] * thn
        log(f"[INFO] theta0 prior: {int((tau_c > 0).sum()):,} clusters with tau > 0")
if POOL:                                                  # non-representative pool members stay out of the EM
    for r_, ms in POOLG.items():
        others = [m for m in ms if m != r_]
        alpha[r_] = alpha[ms].sum(); prior[r_] = prior[ms].sum(); alpha[others] = 0.0; prior[others] = 0.0

# ---------------------------------------------------------------- EM
t = time.time(); it = 0; rel = 0.0
nK = cnt[K] if C else np.zeros(0)
counts = single.copy()
for it in range(1, a.max_iter + 1):
    if C:
        num = (alpha / eff)[I] * F
        den = np.bincount(K, weights=num, minlength=C)
        w = num / np.maximum(den[K], 1e-300) * nK
        counts = single + np.bincount(I, weights=w, minlength=T)
    new = counts + prior
    big = new > 1.0
    dchg = np.abs(new[big] - alpha[big]) / new[big] if big.any() else np.zeros(0)
    rel = float(dchg.max()) if dchg.size else 0.0
    alpha = new
    if it >= a.min_iter and rel < a.tol: break
if rel < a.tol:
    log(f"[INFO] EM converged: {it} iterations (last relative change {rel:.2e}), {time.time()-t:.1f}s; assigned {counts.sum():,.0f} / aligned {N_al:,} fragments")
else:
    slow = dchg >= a.tol
    log(f"[INFO] EM stopped at --max-iter {a.max_iter} before reaching --tol {a.tol:g} (last relative change {rel:.2e}), "
        f"{time.time()-t:.1f}s; assigned {counts.sum():,.0f} / aligned {N_al:,} fragments")
    log(f"[INFO] targets still changing by >= {a.tol:g} per iteration: {int(slow.sum()):,} of {int(big.sum()):,} with > 1 fragment, "
        f"holding {alpha[big][slow].sum() / max(alpha.sum(), 1e-300):.3%} of the fragments")

def unpool(v):
    v = v.copy()
    for r_, ms in POOLG.items(): v[ms] = v[r_] / len(ms)
    return v
counts = unpool(counts); alpha = unpool(alpha)
fin = np.ones(T)
if is_al.any():
    Sx = np.bincount(cl_int[is_al], weights=alpha[is_al]); z = Sx[cl_int[is_al]]
    fin[is_al] = np.where(z > 0, alpha[is_al] / np.maximum(z, 1e-300), 0.0)
with open(a.out, "w") as fo:
    fo.write("ClusterID\tAlleleID\tAncestry\tInitial_Theta\tEM_Allocated_Reads\tFinal_Theta\tEM_Iter\tEM_LogLik\n")
    for i in range(T):
        n = names[i]
        c = tcl[i] if tcl[i] is not None else "ORPHAN"
        lin = anc.get(n, anc.get(n.rsplit(".", 1)[0], "UNKNOWN" if tcl[i] is not None else "ORPHAN"))
        fo.write(f"{c}\t{n}\t{lin}\t{th0[i]:.6f}\t{counts[i]:.2f}\t{fin[i]:.6f}\t{it}\t0\n")
log(f"[SUCCESS] EM complete: {a.out} ({time.time()-T0:.0f}s)")
