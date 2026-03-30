# -*- coding: utf-8 -*-
"""
LineageQuant Module 01 – Global K-mer Theta Engine
Builds allele-specific k-mer probe dictionaries and computes
per-allele normalised abundance scores (Theta) from raw FASTQ reads.

Features:
  - Checkpoint / resume: skips k-mer index rebuild if master_index.pkl exists
  - Double-strand FASTQ scanning (forward + reverse-complement)
  - Case-insensitive (soft-masked sequences handled via forced uppercase)
"""

import os
import sys
import gzip
import pickle
import subprocess
from collections import defaultdict

try:
    from Bio import SeqIO
except ImportError:
    print("[FATAL] biopython not found. Install with: pip install biopython")
    sys.exit(1)


def get_now():
    import datetime
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="LineageQuant Module 01: Global K-mer Theta Engine"
    )
    parser.add_argument("-c", "--clusters",  required=True,
                        help="Allele cluster TSV: ClusterID<TAB>Allele1,Allele2,...")
    parser.add_argument("-g", "--gff",       required=True, help="Whole-genome GFF3")
    parser.add_argument("-f", "--fasta",     required=True, help="Whole-genome FASTA")
    parser.add_argument("-r1", "--read1",    required=True, help="RNA-seq R1 FASTQ (.fq.gz)")
    parser.add_argument("-r2", "--read2",    required=True, help="RNA-seq R2 FASTQ (.fq.gz)")
    parser.add_argument("-k",  "--kmer",     type=int, default=31,
                        help="K-mer size (default: 31)")
    parser.add_argument("-o", "--outdir",    default="theta_out",
                        help="Output directory (default: theta_out)")
    parser.add_argument("-n", "--out-name",  dest="out_name",
                        default="Theta_Matrix.tsv",
                        help="Output filename inside --outdir (default: Theta_Matrix.tsv)")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    k = args.kmer

    # Checkpoint file paths
    temp_gff    = os.path.join(args.outdir, "global_temp.gff3")
    raw_mrna_fa = os.path.join(args.outdir, "global_mrna.fa")
    index_pkl   = os.path.join(args.outdir, "master_index.pkl")
    out_tsv     = os.path.join(args.outdir, args.out_name)

    # ---- Step 1: Parse cluster table ----
    print(f"[{get_now()}] [STEP 1] Parsing cluster table...")
    cluster_dict   = {}
    all_target_ids = set()
    gene_to_cid    = {}

    with open(args.clusters, 'r') as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            parts  = line.strip().split('\t')
            cid    = parts[0]
            alleles = parts[1].split(',')
            cluster_dict[cid] = []
            for a in alleles:
                aid = a.strip().split('.')[0]
                cluster_dict[cid].append(aid)
                all_target_ids.add(aid)
                gene_to_cid[aid] = cid

    print(f"  -> Loaded {len(cluster_dict)} clusters, {len(all_target_ids)} alleles.")

    # ---- Steps 2-4: Build k-mer index (or load from checkpoint) ----
    if os.path.exists(index_pkl):
        print(f"[{get_now()}] [INFO] Checkpoint detected ({index_pkl}). Loading k-mer index...")
        with open(index_pkl, 'rb') as f:
            global_master_dict, probe_counts = pickle.load(f)
        print(f"  -> Loaded {len(global_master_dict)} unique probes.")
    else:
        print(f"[{get_now()}] [STEP 2] Extracting mRNA sequences (gffread)...")
        # Write a reduced GFF containing only target genes to speed up gffread
        with open(args.gff, 'r') as fin, open(temp_gff, 'w') as fout:
            for line in fin:
                if line.startswith('#'):
                    continue
                for tid in all_target_ids:
                    if tid in line:
                        fout.write(line)
                        break

        cmd = ["gffread", temp_gff, "-g", args.fasta, "-w", raw_mrna_fa]
        print(f"  -> Executing: {' '.join(cmd)}")
        subprocess.run(cmd, check=True)

        print(f"[{get_now()}] [STEP 3] Building k-mer dictionary (strictly unique probes)...")
        kmer_to_alleles = defaultdict(set)

        for record in SeqIO.parse(raw_mrna_fa, "fasta"):
            header_id = record.id
            matched_aid = None
            for tid in all_target_ids:
                if tid in header_id:
                    matched_aid = tid
                    break
            if not matched_aid:
                continue

            # Force uppercase to handle soft-masked (lowercase) genome sequences
            seq_str = str(record.seq).upper()
            for i in range(len(seq_str) - k + 1):
                kmer = seq_str[i:i + k]
                if 'N' not in kmer:
                    kmer_to_alleles[kmer].add(matched_aid)

        print(f"[{get_now()}] [STEP 4] Filtering to uniquely mapping probes and saving checkpoint...")
        global_master_dict = {}
        probe_counts       = defaultdict(int)

        for kmer, aids in kmer_to_alleles.items():
            if len(aids) == 1:
                aid  = list(aids)[0]
                cid  = gene_to_cid[aid]
                global_master_dict[kmer] = (cid, aid)
                probe_counts[aid] += 1

        with open(index_pkl, 'wb') as f:
            pickle.dump((global_master_dict, probe_counts), f)
        print(f"  -> Retained {len(global_master_dict)} unique probes. Checkpoint saved.")

    if not global_master_dict:
        print("[FATAL] No unique probes found. Check that cluster IDs match the FASTA headers.")
        sys.exit(1)

    # ---- Step 5: Single-pass double-strand FASTQ scan ----
    print(f"[{get_now()}] [STEP 5] Scanning FASTQ reads (double-strand, uppercase)...")
    hit_counts = defaultdict(int)

    trans = str.maketrans('ATCGN', 'TAGCN')

    def rev_comp(seq):
        return seq.translate(trans)[::-1]

    def scan_fastq(fq_file):
        print(f"  -> Scanning {os.path.basename(fq_file)}...")
        open_fn = gzip.open if fq_file.endswith('.gz') else open
        with open_fn(fq_file, 'rt') as f:
            for i, line in enumerate(f):
                if i % 4 != 1:
                    continue
                # Force uppercase to match the probe dictionary
                read_seq = line.strip().upper()
                rc_seq   = rev_comp(read_seq)

                for j in range(len(read_seq) - k + 1):
                    kmer = read_seq[j:j + k]
                    if kmer in global_master_dict:
                        _, aid = global_master_dict[kmer]
                        hit_counts[aid] += 1

                for j in range(len(rc_seq) - k + 1):
                    kmer = rc_seq[j:j + k]
                    if kmer in global_master_dict:
                        _, aid = global_master_dict[kmer]
                        hit_counts[aid] += 1

    scan_fastq(args.read1)
    scan_fastq(args.read2)

    # ---- Step 6: Compute Theta and write output ----
    print(f"[{get_now()}] [STEP 6] Computing Theta matrix...")
    with open(out_tsv, "w") as out_f:
        out_f.write("Cluster_ID\tAllele_ID\tUnique_Probes\tRaw_Hits\t"
                    "Normalized_Score\tTheta\n")

        for cid, alleles in cluster_dict.items():
            raw_scores  = {}
            total_score = 0.0

            for aid in alleles:
                probes = probe_counts.get(aid, 0)
                hits   = hit_counts.get(aid, 0)
                score  = (hits / probes) if probes > 0 else 0.0
                raw_scores[aid] = score
                total_score    += score

            for aid in alleles:
                theta = (raw_scores[aid] / total_score) if total_score > 0 else 0.0
                out_f.write(
                    f"{cid}\t{aid}\t{probe_counts.get(aid, 0)}\t"
                    f"{hit_counts.get(aid, 0)}\t{raw_scores[aid]:.4f}\t{theta:.6f}\n"
                )

    # Clean up large temporary GFF
    if os.path.exists(temp_gff):
        os.remove(temp_gff)

    print(f"[{get_now()}] [SUCCESS] Theta matrix written to: {out_tsv}")
    print("=" * 60)


if __name__ == "__main__":
    main()
