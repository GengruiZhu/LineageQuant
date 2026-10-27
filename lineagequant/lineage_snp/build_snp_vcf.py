# -*- coding: utf-8 -*-
"""
build_snp_vcf.py -- population-VCF SNP allele-frequency DB (Step 1a, VCF route).

Builds a lineage-specific SNP database from a population VCF, keyed by
(transcript_id, transcript_0based_pos). Genomic coordinates are projected onto
transcript CDS coordinates via GFF3 exon structure to match what
pysam.get_reference_positions() returns on transcriptome-aligned BAMs.

Lineage sample groups are passed as --sample-lists NAME:file [NAME:file ...]
(any number N), replacing the hard-coded SRO/SSP/SOF arguments.
"""
import os
import sys
import pickle
import argparse
from collections import defaultdict
import numpy as np


def build_genomic_to_transcript_map(gff3_path):
    tx_exons = defaultdict(list)
    with open(gff3_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            p = line.strip().split("\t")
            if len(p) < 9 or p[2].lower() not in ("cds", "exon"):
                continue
            chrom, start, end, strand = p[0], int(p[3]), int(p[4]), p[6]
            tx_id = None
            for kv in p[8].split(";"):
                if kv.strip().startswith("Parent="):
                    tx_id = kv.strip().split("=")[1].split(".")[0]
                    break
            if tx_id:
                tx_exons[(chrom, tx_id, strand)].append((start, end))
    g2t = defaultdict(list)
    for (chrom, tx_id, strand), exons in tx_exons.items():
        exons_sorted = sorted(exons, key=lambda x: x[0], reverse=(strand == "-"))
        offset = 0
        for start, end in exons_sorted:
            for gpos in range(start, end + 1):
                tpos = offset + (end - gpos if strand == "-" else gpos - start)
                g2t[(chrom, gpos)].append((tx_id, tpos))
            offset += end - start + 1
    return dict(g2t)


def build_snp_db(vcf_path, g2t, sample_groups, lineages,
                 min_dp=3, min_called=3, max_dp=200, min_af_diff=0.0):
    from cyvcf2 import VCF
    vcf = VCF(vcf_path)
    all_samples = vcf.samples
    idx = {l: np.array([i for i, s in enumerate(all_samples) if s in sample_groups[l]],
                       dtype=int) for l in lineages}
    for l in lineages:
        if idx[l].size == 0:
            print(f"[WARNING] lineage '{l}' has no matching samples in VCF header.")

    snp_db = {}
    stats = {"total": 0, "cds": 0, "entries": 0}

    def calc_af(arr, dp_all, ad_all):
        if arr.size == 0:
            return 0.5
        gdp, gad = dp_all[arr], ad_all[arr]
        mask = (gdp >= min_dp) & (gdp <= max_dp)
        if np.sum(mask) < min_called:
            return 0.5
        v = gad[mask]
        alt, ref = np.sum(v[:, 1]), np.sum(v[:, 0])
        tot = alt + ref
        return float(alt) / tot if tot > 0 else 0.5

    for rec in vcf:
        if len(rec.ALT) != 1 or len(rec.REF) > 1 or len(rec.ALT[0]) > 1:
            continue
        stats["total"] += 1
        entries = g2t.get((rec.CHROM, rec.POS))
        if not entries:
            continue
        ad_all = rec.format("AD")
        dp_all = rec.format("DP")
        if ad_all is None or dp_all is None:
            continue
        dp_all = dp_all.flatten()
        afs = {l: calc_af(idx[l], dp_all, ad_all) for l in lineages}
        if min_af_diff > 0 and (max(afs.values()) - min(afs.values())) < min_af_diff:
            continue
        stats["cds"] += 1
        for tx_id, tpos in entries:
            entry = {"ref": rec.REF, "alt": rec.ALT[0]}
            entry.update(afs)
            snp_db[(tx_id, tpos)] = entry
            stats["entries"] += 1
    return snp_db, stats


def main():
    parser = argparse.ArgumentParser(description="Population VCF SNP DB (Step 1a)")
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--gff3", required=True)
    parser.add_argument("--sample-lists", nargs="+", required=True, metavar="NAME:FILE",
                        help="Per-lineage sample-name files (>=2). Any number N.")
    parser.add_argument("--out", required=True)
    parser.add_argument("--min_dp", type=int, default=3)
    parser.add_argument("--min_called", type=int, default=3)
    parser.add_argument("--max_dp", type=int, default=200)
    parser.add_argument("--min_af_diff", type=float, default=0.0)
    args = parser.parse_args()

    def load_list(path):
        with open(path) as f:
            return set(line.strip() for line in f if line.strip())

    sample_groups = {}
    for item in args.sample_lists:
        if ":" not in item:
            parser.error(f"--sample-lists must be NAME:FILE, got {item}")
        name, path = item.split(":", 1)
        sample_groups[name] = load_list(path)
    lineages = sorted(sample_groups.keys())
    print(f"[INFO] {len(lineages)} lineages: {lineages}")

    print("[STEP 1] Projecting genomic->transcript coordinates...")
    g2t = build_genomic_to_transcript_map(args.gff3)
    print(f"  {len(g2t):,} genomic CDS positions")

    print("[STEP 2] Building SNP DB from VCF...")
    snp_db, stats = build_snp_db(args.vcf, g2t, sample_groups, lineages,
                                 min_dp=args.min_dp, min_called=args.min_called,
                                 max_dp=args.max_dp, min_af_diff=args.min_af_diff)
    print(f"\n[SUMMARY] biallelic={stats['total']:,} in_CDS={stats['cds']:,} "
          f"entries={stats['entries']:,}")

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "wb") as f:
        pickle.dump(snp_db, f, protocol=pickle.HIGHEST_PROTOCOL)
    print(f"[SUCCESS] {len(snp_db):,} entries -> {args.out} "
          f"({os.path.getsize(args.out)/1e6:.1f} MB)")


if __name__ == "__main__":
    main()
