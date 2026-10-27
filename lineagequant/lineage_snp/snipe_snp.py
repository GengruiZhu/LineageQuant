# -*- coding: utf-8 -*-
"""
snipe_snp.py -- per-read lineage probability from SNPs (Step 02_v2 / 02_v2_hybrid).

Merged sniper covering both SNP routes via --mode:
  --mode msa  (default): base-matching against MSA-derived fixed differences.
      Each SNP entry has bases:{lineage:base}; read_base==expected -> P=1-eps.
      Supports calibrated AFs (--calibrated) and chimeric detection.
  --mode population: genotype-likelihood against population VCF allele freqs.
      P(base|lineage) = af*(1-eps) + (1-af)*(eps/3) per SNP.

Lineage names are read from the SNP DB itself (the per-entry AF / bases keys),
so any number N of lineages is supported -- no hard-coded SOF/SRO/SSP.

Memory (1.0.0 fix, 2026-10-14): the SNP DB is converted once, in the parent, into
a compact per-transcript index (sorted positions + small numeric arrays) and the
Python dict is freed before the worker pool is forked. Workers share the compact
index read-only; previously every worker rebuilt a full dict copy of the DB, so
memory grew with the number of threads (~800 GB at 32 threads on hybrid
sugarcane). Per fragment only the SNP observations are kept, not the read
sequence and the full reference-position list. Output is identical to the
previous version.
"""
import gc
import sys
import pickle
import argparse
import bisect
import multiprocessing as mp
from collections import defaultdict
import pysam
import numpy as np

BASES = {"A", "C", "G", "T"}
_CODE = {"A": 0, "C": 1, "G": 2, "T": 3}
_OTHER = 4                      # expected base that can never match a read base

_index = _resolve = _lineages = _bam = _eps = _mode = None
_calibrated = _no_chimeric = None
_LM = _LX = None                # log(1-eps), log(eps/3)


def _infer_lineages(snp_db):
    """Infer lineage names from SNP DB entries (bases keys or AF keys)."""
    for entry in snp_db.values():
        if "bases" in entry and entry["bases"]:
            return sorted(entry["bases"].keys())
        keys = [k for k in entry if k not in ("bases", "ref", "alt")]
        if keys:
            return sorted(keys)
    return []


def _code(b):
    return _CODE.get(str(b).upper(), _OTHER)


def build_compact_index(snp_db, lineages, mode, eps):
    """{tid: record} with positions sorted ascending.

    msa:        (pos int64[n], exp uint8[n, n_lin])        exp = expected base code
    population: (pos int64[n], ref uint8[n], alt uint8[n],
                 logref float64[n, n_lin], logalt float64[n, n_lin])
    The per-SNP log-likelihood terms are computed with exactly the expressions the
    previous version evaluated per read, so the results are bit-identical.
    """
    n_lin = len(lineages)
    tmp = {}                                    # keeps first-appearance order of tids
    for (tid, pos), entry in snp_db.items():
        if mode == "msa":
            expected = entry.get("bases", {})
            if not expected:
                ref_a = entry["ref"].upper()
                alt_a = entry["alt"].upper()
                expected = {l: (alt_a if entry.get(l, 0.0) > 0.5 else ref_a)
                            for l in lineages}
            row = bytes(_code(expected.get(l, "N")) for l in lineages)
            rec = (int(pos), row)
        else:
            ref_a = entry["ref"].upper()
            alt_a = entry["alt"].upper()
            lr, la = [], []
            for l in lineages:
                af = entry.get(l, 0.5)
                lik_r = (1 - af) * (1 - eps) + af * (eps / 3)
                lik_a = af * (1 - eps) + (1 - af) * (eps / 3)
                lr.append(np.log(max(lik_r, eps / 3)))
                la.append(np.log(max(lik_a, eps / 3)))
            rec = (int(pos), _code(ref_a) if len(ref_a) == 1 else _OTHER,
                   _code(alt_a) if len(alt_a) == 1 else _OTHER, lr, la)
        tmp.setdefault(tid, []).append(rec)

    index = {}
    for tid, recs in tmp.items():
        recs.sort(key=lambda r: r[0])
        pos = np.fromiter((r[0] for r in recs), dtype=np.int64, count=len(recs))
        if mode == "msa":
            exp = np.frombuffer(b"".join(r[1] for r in recs), dtype=np.uint8)
            index[tid] = (pos, exp.reshape(len(recs), n_lin).copy())
        else:
            index[tid] = (pos,
                          np.array([r[1] for r in recs], dtype=np.uint8),
                          np.array([r[2] for r in recs], dtype=np.uint8),
                          np.array([r[3] for r in recs], dtype=np.float64).reshape(len(recs), n_lin),
                          np.array([r[4] for r in recs], dtype=np.float64).reshape(len(recs), n_lin))
        tmp[tid] = None
    return index


def _build_resolve(idx):
    resolve = {}
    base_to_var = defaultdict(list)
    for tid in idx:
        resolve[tid] = tid
        base = tid.rsplit(".", 1)[0]
        if base != tid:
            base_to_var[base].append(tid)
    for base, variants in base_to_var.items():
        resolve.setdefault(base, variants[0])
    return resolve


def init_worker(index, resolve, lineages, bam, eps, mode, calibrated, no_chimeric):
    global _index, _resolve, _lineages, _bam, _eps, _mode, _calibrated, _no_chimeric, _LM, _LX
    _index, _resolve = index, resolve          # shared read-only, no per-worker copy
    _lineages, _bam, _eps, _mode = lineages, bam, eps, mode
    _calibrated, _no_chimeric = calibrated, no_chimeric
    _LM, _LX = np.log(1.0 - eps), np.log(eps / 3.0)


def _resolve_tid(bam_tid):
    r = _resolve.get(bam_tid)
    if r is None:
        r = _resolve.get(bam_tid.rsplit(".", 1)[0])
    return r


class _TxTable:
    """Per-transcript lookup built in the worker (small, transient)."""
    __slots__ = ("pos", "vec", "nlin")

    def __init__(self, rec, n_lin):
        self.nlin = n_lin
        self.pos = rec[0].tolist()
        n = len(self.pos)
        b4 = np.arange(4, dtype=np.uint8)
        # vec[j, b] = log-likelihood vector added when a read shows base b at SNP j
        if _mode == "msa":
            exp = rec[1]
            self.vec = np.where(exp[:, None, :] == b4[None, :, None], _LM, _LX)
        else:
            ref, alt, lr, la = rec[1], rec[2], rec[3], rec[4]
            vec = np.full((n, 4, n_lin), _LX, dtype=np.float64)
            for b in range(4):
                m_ref = ref == b
                m_alt = (alt == b) & ~m_ref         # ref is tested first, as before
                vec[m_ref, b, :] = lr[m_ref]
                vec[m_alt, b, :] = la[m_alt]
            self.vec = vec


def _collect_hits(read, seq, tab, hits):
    """Append (snp_j, base_code) for every SNP the read covers, in query order."""
    pos = tab.pos
    rs = read.reference_start
    re_ = read.reference_end
    if rs is None or re_ is None:
        return
    lo = bisect.bisect_left(pos, rs)
    hi = bisect.bisect_left(pos, re_)
    if lo >= hi:
        return
    cig = read.cigartuples
    if cig is not None and len(cig) == 1 and cig[0][0] in (0, 7, 8) and cig[0][1] == len(seq):
        # ungapped, unclipped alignment (the bowtie2 setting LineageQuant uses)
        for j in range(lo, hi):
            c = _CODE.get(seq[pos[j] - rs].upper())
            if c is not None:
                hits.append((j, c))
        return
    for qpos, rpos in enumerate(read.get_reference_positions(full_length=True)):
        if rpos is None or rpos < rs or rpos >= re_:
            continue
        j = bisect.bisect_left(pos, rpos, lo, hi)
        if j < hi and pos[j] == rpos:
            c = _CODE.get(seq[qpos].upper())
            if c is not None:
                hits.append((j, c))


def _fragment_probs(hits, tab):
    n_lin = len(_lineages)
    if tab is None or not hits:
        return {l: 1.0 / n_lin for l in _lineages}, "No_SNP_Covered", 0
    log_lik = np.zeros(n_lin)
    vec = tab.vec
    for j, c in hits:
        log_lik += vec[j, c]
    total = len(hits)
    max_ll = log_lik.max()
    if _mode == "msa" and total >= 3 and not _no_chimeric:
        if max_ll < max(-6.0, -2.0 * total):
            return {l: 1.0 / n_lin for l in _lineages}, "Chimeric_Drop", total
    lik = np.exp(log_lik - max_ll)
    lik /= lik.sum()
    probs = {l: float(lik[i]) for i, l in enumerate(_lineages)}
    best = max(probs, key=probs.get)
    thresh = 0.99 if _mode == "msa" else 0.9
    status = f"Specific_Soft({best})" if probs[best] > thresh else "Multi_Lineage"
    return probs, status, total


def process_transcript_chunk(id_chunk):
    bam = pysam.AlignmentFile(_bam, "rb")
    out = []
    counts = defaultdict(int)
    n_lin = len(_lineages)
    for tid in id_chunk:
        gene_id = tid.rsplit(".", 1)[0]
        rid = _resolve_tid(tid)
        rec = _index.get(rid) if rid is not None else None
        tab = _TxTable(rec, n_lin) if rec is not None else None
        frag = {}
        try:
            for read in bam.fetch(tid):
                if read.is_unmapped:
                    continue
                seq = read.query_sequence
                if seq is None:
                    continue
                hits = frag.get(read.query_name)
                if hits is None:
                    hits = frag[read.query_name] = []
                if tab is not None:
                    _collect_hits(read, seq, tab, hits)
        except ValueError:
            continue
        for qname, hits in frag.items():
            probs, status, n = _fragment_probs(hits, tab)
            p_str = "\t".join(f"{probs[l]:.4f}" for l in _lineages)
            out.append(f"{qname}\t{gene_id}\t{p_str}\t{status}\t{n}\n")
            counts[status] += 1
        frag = None
    bam.close()
    return "".join(out), dict(counts), len(out)


def main():
    parser = argparse.ArgumentParser(description="SNP lineage sniper (Step 02_v2)")
    parser.add_argument("-i", "--ids", required=True)
    parser.add_argument("-b", "--bam", required=True)
    parser.add_argument("--snp_db", required=True)
    parser.add_argument("-o", "--out", default="P_values_snp.tsv")
    parser.add_argument("-t", "--threads", type=int, default=16)
    parser.add_argument("--error_rate", type=float, default=0.01)
    parser.add_argument("--mode", choices=["msa", "population"], default="msa",
                        help="msa: MSA fixed-difference base matching; "
                             "population: VCF allele-frequency genotype likelihood")
    parser.add_argument("--calibrated", action="store_true", default=False,
                        help="(msa) SNP DB AFs were calibrated by long reads")
    parser.add_argument("--no_chimeric", action="store_true", default=False,
                        help="(msa) disable chimeric detection")
    parser.add_argument("--chunk", type=int, default=100,
                        help="transcripts per worker task (default 100)")
    args = parser.parse_args()

    print(f"[INFO] Loading SNP DB: {args.snp_db} (mode={args.mode})")
    with open(args.snp_db, "rb") as f:
        snp_db = pickle.load(f)
    print(f"  {len(snp_db):,} SNP entries")
    lineages = _infer_lineages(snp_db)
    if not lineages:
        print("[FATAL] Could not infer lineage names from SNP DB.")
        sys.exit(1)
    print(f"  {len(lineages)} lineages: {lineages}")

    index = build_compact_index(snp_db, lineages, args.mode, args.error_rate)
    del snp_db
    gc.collect()
    resolve = _build_resolve(index)
    nbytes = sum(a.nbytes for rec in index.values() for a in rec)
    print(f"[INFO] compact SNP index: {len(index):,} transcripts, "
          f"{nbytes / 1e9:.2f} GB arrays (shared by all workers)")

    with open(args.ids) as f:
        valid = {line.strip().split(".")[0] for line in f if line.strip()}
    bam_info = pysam.AlignmentFile(args.bam, "rb")
    all_tids = bam_info.references
    bam_info.close()
    target_tids = [t for t in all_tids if t.split(".")[0] in valid]
    print(f"[INFO] {len(target_tids)} target transcripts in BAM.")

    chunk = max(1, args.chunk)
    id_chunks = [target_tids[i:i + chunk] for i in range(0, len(target_tids), chunk)]
    try:
        ctx = mp.get_context("fork")            # workers inherit the index copy-on-write
    except ValueError:
        ctx = mp
    if hasattr(gc, "freeze"):
        gc.freeze()                             # keep the GC from touching shared pages
    pool = ctx.Pool(args.threads, initializer=init_worker,
                    initargs=(index, resolve, lineages, args.bam, args.error_rate,
                              args.mode, args.calibrated, args.no_chimeric))
    total = 0
    status_counts = defaultdict(int)
    with open(args.out, "w") as out_f:
        out_f.write("Read_Name\tTarget_Gene\t" +
                    "\t".join(f"P({l})" for l in lineages) + "\tLineage_Status\tN_SNPs\n")
        for text, counts, n in pool.imap_unordered(process_transcript_chunk, id_chunks):
            out_f.write(text)
            total += n
            for s, c in counts.items():
                status_counts[s] += c
    pool.close()
    pool.join()
    print(f"\n[SUMMARY] {total:,} fragments")
    for s, c in sorted(status_counts.items(), key=lambda x: -x[1]):
        print(f"  {s}: {c:,} ({c/max(total,1)*100:.1f}%)")
    print(f"[SUCCESS] {args.out}")


if __name__ == "__main__":
    main()
