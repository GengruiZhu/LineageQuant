# -*- coding: utf-8 -*-
"""
calibrate_snp.py -- long-read empirical calibration of MSA SNP AFs (Step 02c).

Replaces deterministic 0/1 allele frequencies (from MSA) with empirical
estimates from long-read pileup at each SNP position. This eliminates the
information-theoretic mismatch that causes systematic Chimeric_Drop in the MSA
SNP route. Only meaningful for the MSA route (VCF AFs are already empirical).

Lineage names are read from the ancestry file (any number N).
Output has identical structure to the input SNP DB -> drop-in for snipe_snp.
"""
import os
import sys
import pickle
import argparse
import time
from collections import defaultdict
import numpy as np

try:
    import pysam
except ImportError:
    print("[FATAL] pysam not found.")
    sys.exit(1)


def get_now():
    return time.strftime("%Y-%m-%d %H:%M:%S")


def load_ancestry(path):
    amap = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("GeneID") or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 2:
                amap[parts[0].split(".")[0]] = parts[1].upper()
    return amap


def collect_pileup(bam_path, snp_db, ancestry):
    tx_pos = defaultdict(dict)
    tx_lineage = {}
    for (tid, pos), entry in snp_db.items():
        tx_pos[tid][pos] = entry
        tx_lineage[tid] = ancestry.get(tid.rsplit(".", 1)[0], "UNKNOWN")
    alt_c, ref_c, entries_out = defaultdict(int), defaultdict(int), {}
    resolve = {}
    bam = pysam.AlignmentFile(bam_path, "rb")
    n_hit = 0
    # whole file in its own order (no index needed, any sort order); a secondary alignment is the
    # same bases again (minimap2 writes it without SEQ anyway), QC-fail / duplicate reads are skipped
    for read in bam.fetch(until_eof=True):
        if (read.is_unmapped or read.is_secondary or read.is_qcfail or read.is_duplicate
                or read.query_sequence is None):
            continue
        bam_tid = bam.get_reference_name(read.reference_id)
        rid = resolve.get(bam_tid)
        if rid is None:
            if bam_tid in tx_pos:
                rid = bam_tid
            else:
                base = bam_tid.rsplit(".", 1)[0]
                rid = base if base in tx_pos else False
            resolve[bam_tid] = rid
        if rid is False:
            continue
        snp_map = tx_pos[rid]
        ref_pos = read.get_reference_positions(full_length=True)
        seq = read.query_sequence
        for qpos, rpos in enumerate(ref_pos):
            if rpos is None:
                continue
            entry = snp_map.get(rpos)
            if entry is None:
                continue
            b = seq[qpos].upper()
            if b not in {"A", "C", "G", "T"}:
                continue
            key = (rid, rpos)
            entries_out[key] = entry
            if b == entry["ref"].upper():
                ref_c[key] += 1
            elif b == entry["alt"].upper():
                alt_c[key] += 1
            n_hit += 1
    bam.close()
    print(f"[{get_now()}] [INFO] pileup: {n_hit:,} SNP-hits, "
          f"{len(set(alt_c)|set(ref_c)):,} positions")
    return alt_c, ref_c, tx_lineage, entries_out


def calibrate(snp_db, alt_c, ref_c, tx_lineage, entries_out, lineages,
              min_depth=10, det_weight=0.25, eps=0.5):
    calibrated, per_pos_acc, covered = {}, [], set()
    for key, entry in entries_out.items():
        tid, pos = key
        lin = tx_lineage.get(tid, "UNKNOWN")
        det_af = entry.get(lin, 0.0)
        depth = alt_c.get(key, 0) + ref_c.get(key, 0)
        new = dict(entry)
        if depth >= min_depth:
            emp_alt = (alt_c.get(key, 0) + eps) / (depth + 2 * eps)
            cov_factor = min(1.0, depth / (min_depth * 3))
            eff_w = det_weight + (1 - det_weight) * (1 - cov_factor) * 0.33
            caf = max(0.005, min(0.995, eff_w * det_af + (1 - eff_w) * emp_alt))
            if lin in lineages:
                new[lin] = caf
            per_pos_acc.append(emp_alt if det_af > 0.5 else 1 - emp_alt)
        calibrated[key] = new
        covered.add(key)
    for key, entry in snp_db.items():
        if key not in covered:
            calibrated[key] = dict(entry)
    global_acc = float(np.median(per_pos_acc)) if per_pos_acc else 0.98
    # global softening of remaining deterministic AFs
    for key, entry in calibrated.items():
        tid, pos = key
        lin = tx_lineage.get(tid, "UNKNOWN")
        for l in lineages:
            depth = alt_c.get(key, 0) + ref_c.get(key, 0)
            if l == lin and depth >= min_depth:
                continue
            det = snp_db.get(key, {}).get(l, 0.0)
            entry[l] = global_acc if det > 0.5 else max(0.001, (1 - global_acc) / 3)
    return calibrated, global_acc


def main():
    parser = argparse.ArgumentParser(description="Long-read SNP AF calibration (Step 02c)")
    parser.add_argument("--snp_db", required=True)
    parser.add_argument("--bam", required=True, help="Long-read BAM -> phased ref")
    parser.add_argument("--ancestry", required=True)
    parser.add_argument("--out", default="Calibrated_SNP_DB.pkl")
    parser.add_argument("--stats_out", default=None)
    parser.add_argument("--min_depth", type=int, default=10)
    parser.add_argument("--det_weight", type=float, default=0.25)
    parser.add_argument("--eps", type=float, default=0.5)
    args = parser.parse_args()

    print(f"[{get_now()}] [LOAD] {args.snp_db}")
    with open(args.snp_db, "rb") as f:
        snp_db = pickle.load(f)
    ancestry = load_ancestry(args.ancestry)
    lineages = sorted(set(ancestry.values()))
    print(f"  {len(snp_db):,} SNPs | {len(lineages)} lineages: {lineages}")

    bam = args.bam     # read sequentially: sorted or not, no index needed

    alt_c, ref_c, tx_lineage, entries_out = collect_pileup(bam, snp_db, ancestry)
    calibrated, global_acc = calibrate(snp_db, alt_c, ref_c, tx_lineage,
                                       entries_out, lineages,
                                       min_depth=args.min_depth,
                                       det_weight=args.det_weight, eps=args.eps)
    with open(args.out, "wb") as f:
        pickle.dump(calibrated, f, protocol=pickle.HIGHEST_PROTOCOL)
    print(f"[{get_now()}] [SUCCESS] {args.out} | global_acc={global_acc:.4f}")
    if args.stats_out:
        with open(args.stats_out, "w") as f:
            f.write(f"metric\tvalue\nglobal_accuracy\t{global_acc}\n"
                    f"total_snps\t{len(snp_db)}\n")


if __name__ == "__main__":
    main()
