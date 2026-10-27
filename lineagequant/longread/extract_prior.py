# -*- coding: utf-8 -*-
"""
extract_prior.py -- long-read physical prior (Step 00b): long-read molecules per allele copy.

1.0.0 took `samtools idxstats` (every mapped record: secondary and supplementary alignments
included, and every subread pass of a PacBio ZMW) and collapsed transcripts to alleles with
`=` (the last isoform in the BAM header overwrote the others). That inflates the counts (and
with them the Dirichlet weight tau = theta_prior_lr x cluster count in the EM) and flattens the
copy proportions, because a read from one homoeolog also leaves secondary records on the others.

Now every molecule is counted once:
  * records are grouped per read (the BAM is collated first unless it is already grouped by
    name, e.g. straight minimap2 output);
  * primary + supplementary alignments of a read are its segments (a MAS/Kinnex array or a
    chimera has several); each secondary alignment joins the segment it overlaps on the read;
  * inside a segment only the alignments with the best AS (within --as-tol) are candidates;
    a segment whose best alignment is shorter than --min-aln-len or below --min-identity is dropped;
  * candidates are allele copies (ID before the first dot, the same key as the Theta matrix),
    so isoforms of one copy add up; segments tied across copies are shared by a small EM over
    equivalence classes (sequence-identical copies split evenly);
  * raw PacBio subreads: all subreads of a ZMW are one molecule, so each subread gets weight
    1 / (subreads of its ZMW) (--zmw-dedup auto|on|off; auto turns it on when read names are
    PacBio subread names and some ZMW has more than one subread).

Output (same format as before): AlleleID<TAB>LongRead_Count (molecules, may be fractional).
--legacy-out writes the 1.0.0 numbers for the same BAM (idxstats + `=`), for comparison.
"""
import os
import re
import sys
import time
import shutil
import argparse
import tempfile
import subprocess
from collections import defaultdict

import numpy as np
import pysam

ZMW_RE = re.compile(r"^(.+/\d+)/(\d+_\d+|ccs.*)$")
Q_OPS = (0, 1, 7, 8)          # M I = X: consume the query inside the alignment
COL_OPS = (0, 1, 2, 7, 8)     # M I D = X: alignment columns
CLIP_OPS = (4, 5)             # S H


def now():
    return time.strftime("%Y-%m-%d %H:%M:%S")


def akey(ref):
    """allele copy of a transcript: ID before the first dot (as in global_theta.py)"""
    return ref.split(".")[0]


def zmw_of(name):
    m = ZMW_RE.match(name)
    return (m.group(1), not m.group(2).startswith("ccs")) if m else (None, False)


def grouped_by_name(path):
    with pysam.AlignmentFile(path, "rb", check_sq=False) as b:
        hd = b.header.to_dict().get("HD", {})
    return hd.get("SO") == "queryname" or hd.get("GO") == "query"


def open_grouped(path, threads, tmpdir):
    """BAM grouped by read name; collate into a temporary file if needed"""
    if grouped_by_name(path):
        return path, None
    td = tempfile.mkdtemp(prefix=".lr_collate_", dir=tmpdir)
    out = os.path.join(td, "collated.bam")
    print(f"[{now()}] [INFO] BAM is not grouped by read name -> samtools collate ({td})", flush=True)
    subprocess.run(["samtools", "collate", "-@", str(threads), "-u", "-o", out, path,
                    os.path.join(td, "tmp")], check=True)
    return out, td


def zmw_weights(path, sample=200000, mode="auto"):
    """read name -> weight 1/(subreads of the ZMW), or None when no de-duplication is needed"""
    if mode == "off":
        return None, "off"
    names, n_pb = [], 0
    with pysam.AlignmentFile(path, "rb", check_sq=False) as b:
        for r in b.fetch(until_eof=True):
            if r.flag & 0x900:
                continue
            names.append(r.query_name)
            if len(names) >= sample:
                break
    per = defaultdict(set)
    for q in names:
        z, sub = zmw_of(q)
        if z is not None and sub:
            n_pb += 1
            per[z].add(q)
    looks_sub = names and n_pb >= 0.9 * len(names)
    multi = any(len(v) > 1 for v in per.values())
    if mode == "auto" and not (looks_sub and multi):
        why = ("not PacBio subread names" if not looks_sub else "one subread per ZMW")
        return None, f"auto -> off ({why})"
    if mode == "on" and not n_pb:
        print(f"[{now()}] [WARN] --zmw-dedup on, but read names carry no ZMW -> off", flush=True)
        return None, "on -> off (no ZMW in read names)"
    print(f"[{now()}] [INFO] subreads of one ZMW = one molecule: counting subreads per ZMW...", flush=True)
    per = defaultdict(int)
    seen = set()
    with pysam.AlignmentFile(path, "rb", check_sq=False) as b:
        for r in b.fetch(until_eof=True):
            if r.flag & 0x900:
                continue
            q = r.query_name
            if q in seen:
                continue
            seen.add(q)
            z, _ = zmw_of(q)
            per[z if z is not None else q] += 1
    w = {}
    for q in seen:
        z, _ = zmw_of(q)
        w[q] = 1.0 / per[z if z is not None else q]
    n_z = len(per)
    print(f"[{now()}] [INFO] {len(seen):,} subreads in {n_z:,} ZMWs ({len(seen) / max(n_z, 1):.2f} per ZMW)", flush=True)
    return w, f"{mode} -> on ({len(seen):,} subreads / {n_z:,} ZMWs)"


def rec_info(r):
    """(is_secondary, q_start, q_end, AS, aligned query length, identity, reference id) on the read's own strand"""
    cig = r.cigartuples
    i, j = 0, len(cig) - 1
    lead = trail = 0
    while i <= j and cig[i][0] in CLIP_OPS:
        lead += cig[i][1]; i += 1
    while j >= i and cig[j][0] in CLIP_OPS:
        trail += cig[j][1]; j -= 1
    qa = cols = 0
    for op, ln in cig[i:j + 1]:
        if op in Q_OPS:
            qa += ln
        if op in COL_OPS:
            cols += ln
    qs = trail if r.is_reverse else lead
    a_s = r.get_tag("AS") if r.has_tag("AS") else None
    nm = r.get_tag("NM") if r.has_tag("NM") else None
    ident = 1.0 - nm / cols if (nm is not None and cols) else 1.0
    if a_s is None:
        a_s = qa - 2 * (nm or 0)
    return (bool(r.flag & 0x100), qs, qs + qa, a_s, qa, ident, r.reference_id)


def segments_of(recs):
    """primary/supplementary alignments are segments; each secondary joins the one it overlaps most"""
    segs = [[x] for x in sorted((x for x in recs if not x[0]), key=lambda x: x[1])]
    for x in recs:
        if not x[0]:
            continue
        best, bo = None, 0
        for s in segs:
            p = s[0]
            ov = min(p[2], x[2]) - max(p[1], x[1])
            if ov > bo:
                best, bo = s, ov
        if best is not None and bo >= 0.5 * max(1, min(x[2] - x[1], best[0][2] - best[0][1])):
            best.append(x)
        else:
            segs.append([x])
    return segs


def run_em(ecs, max_iter=5000, tol=1e-6):
    """ecs: {tuple(sorted allele keys): weight} -> {allele: molecules}"""
    uniq = defaultdict(float)
    multi = []
    for ks, w in ecs.items():
        if len(ks) == 1:
            uniq[ks[0]] += w
        else:
            multi.append((ks, w))
    keys = sorted(set(uniq) | {k for ks, _ in multi for k in ks})
    idx = {k: i for i, k in enumerate(keys)}
    K = len(keys)
    u = np.zeros(K)
    for k, w in uniq.items():
        u[idx[k]] += w
    it = 0
    if multi:
        lens = np.array([len(ks) for ks, _ in multi], dtype=np.int64)
        items = np.fromiter((idx[k] for ks, _ in multi for k in ks), dtype=np.int64, count=int(lens.sum()))
        ec = np.repeat(np.arange(len(multi)), lens)
        w = np.array([w for _, w in multi])
        n = u + np.bincount(items, weights=(w / lens)[ec], minlength=K)          # start: even split
        for it in range(1, max_iter + 1):
            a = n[items]
            den = np.bincount(ec, weights=a, minlength=len(multi))[ec]
            r = np.where(den > 0, a / np.where(den > 0, den, 1.0), 1.0 / lens[ec])
            new = u + np.bincount(items, weights=r * w[ec], minlength=K)
            d = np.abs(new - n).max()
            n = new
            if d < tol:
                break
    else:
        n = u
    return {k: float(n[i]) for i, k in enumerate(keys)}, it


def main():
    ap = argparse.ArgumentParser(description="Long-read physical prior (Step 00b): molecules per allele copy")
    ap.add_argument("-b", "--bam", required=True, help="long-read BAM against the transcriptome (minimap2 -a, "
                    "with secondary alignments, e.g. -N 50); any sort order")
    ap.add_argument("-o", "--out", default="LongRead_Counts.tsv")
    ap.add_argument("--as-tol", type=float, default=0.0,
                    help="candidates: AS >= best AS of the segment - this (default 0 = exact ties only)")
    ap.add_argument("--min-aln-len", type=int, default=150, help="drop segments whose best alignment is shorter")
    ap.add_argument("--min-identity", type=float, default=0.70,
                    help="drop segments whose best alignment has lower identity (1 - NM/columns)")
    ap.add_argument("--zmw-dedup", choices=["auto", "on", "off"], default="auto",
                    help="PacBio subreads: weight 1/(subreads of the ZMW) per subread")
    ap.add_argument("--legacy-out", default=None, help="also write the 1.0.0 counts (idxstats + '=') of this BAM")
    ap.add_argument("--summary", default=None, help="summary table (TSV)")
    ap.add_argument("--threads", type=int, default=4, help="samtools collate threads")
    ap.add_argument("--tmpdir", default=None)
    args = ap.parse_args()

    t0 = time.time()
    print(f"[{now()}] [INFO] long-read BAM: {args.bam}", flush=True)
    path, tmp = open_grouped(args.bam, args.threads, args.tmpdir or os.path.dirname(os.path.abspath(args.out)))
    try:
        weights, zmode = zmw_weights(path, mode=args.zmw_dedup)
        print(f"[{now()}] [INFO] ZMW de-duplication: {zmode}", flush=True)
        bam = pysam.AlignmentFile(path, "rb", check_sq=False)
        refs = list(bam.references)
        rkey = [akey(r) for r in refs]
        legacy_ref = np.zeros(len(refs), dtype=np.int64)
        st = defaultdict(float)
        ecs = defaultdict(float)

        def flush(name, recs):
            st["reads"] += 1
            if not recs:
                st["reads_unmapped"] += 1
                return
            w = 1.0 if weights is None else weights.get(name, 1.0)
            for s in segments_of(recs):
                st["segments"] += 1
                top = max(x[3] for x in s)
                cand = [x for x in s if x[3] >= top - args.as_tol]
                b = max(cand, key=lambda x: (x[3], x[4]))
                if b[4] < args.min_aln_len or b[5] < args.min_identity:
                    st["segments_failed_qc"] += 1
                    continue
                ks = tuple(sorted({rkey[x[6]] for x in cand}))
                st["segments_used"] += 1
                st["molecules"] += w
                if len(ks) > 1:
                    st["segments_tied"] += 1
                ecs[ks] += w

        cur, recs = None, []
        for r in bam.fetch(until_eof=True):
            q = r.query_name
            if q != cur:
                if cur is not None:
                    flush(cur, recs)
                cur, recs = q, []
            f = r.flag
            if f & 0x4:
                continue
            legacy_ref[r.reference_id] += 1          # idxstats: every mapped record
            if f & 0x200 or f & 0x400:
                st["records_qcfail_dup"] += 1
                continue
            st["rec_secondary" if f & 0x100 else ("rec_supplementary" if f & 0x800 else "rec_primary")] += 1
            recs.append(rec_info(r))
        if cur is not None:
            flush(cur, recs)
        bam.close()
    finally:
        if tmp:
            shutil.rmtree(tmp, ignore_errors=True)

    counts, iters = run_em(ecs)
    with open(args.out, "w") as f:
        f.write("AlleleID\tLongRead_Count\n")
        for k in sorted(counts):
            c = round(counts[k], 4)
            if c > 0:
                f.write(f"{k}\t{c}\n")

    legacy = {}
    for i, k in enumerate(rkey):           # 1.0.0: idxstats order, '=' (last isoform wins)
        legacy[k] = int(legacy_ref[i])
    if args.legacy_out:
        with open(args.legacy_out, "w") as f:
            f.write("AlleleID\tLongRead_Count\n")
            for k, c in legacy.items():
                f.write(f"{k}\t{c}\n")

    n_al = sum(1 for v in counts.values() if round(v, 4) > 0)
    rows = [
        ("reads", int(st["reads"])), ("reads_without_alignment", int(st["reads_unmapped"])),
        ("records_primary", int(st["rec_primary"])), ("records_supplementary", int(st["rec_supplementary"])),
        ("records_secondary", int(st["rec_secondary"])), ("records_qcfail_dup", int(st["records_qcfail_dup"])),
        ("zmw_dedup", zmode),
        ("segments", int(st["segments"])), ("segments_failed_qc", int(st["segments_failed_qc"])),
        ("segments_used", int(st["segments_used"])), ("segments_tied_across_copies", int(st["segments_tied"])),
        ("equivalence_classes", len(ecs)), ("em_iterations", iters),
        ("molecules_counted", round(st["molecules"], 2)),
        ("alleles_with_reads", n_al),
        ("legacy_idxstats_total", int(legacy_ref.sum())),
        ("legacy_collapsed_total", int(sum(legacy.values()))),
    ]
    for k, v in rows:
        print(f"  {k:<32}{v:>18,}" if isinstance(v, int) else f"  {k:<32}{v!s:>18}", flush=True)
    if st["molecules"] > 0:
        print(f"  1.0.0 idxstats total / molecules = {legacy_ref.sum() / st['molecules']:.2f}x", flush=True)
    if args.summary:
        with open(args.summary, "w") as f:
            f.write("metric\tvalue\n")
            for k, v in rows:
                f.write(f"{k}\t{v}\n")
    print(f"[{now()}] [SUCCESS] {n_al:,} allele copies with long reads -> {args.out} ({time.time() - t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
