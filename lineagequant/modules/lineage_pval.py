# -*- coding: utf-8 -*-
"""
LineageQuant Module 02 – BAM Lineage P-value Assignment
Assigns a per-lineage probability vector to every read overlapping
target gene loci, based on ancestral k-mer signature matching.

Supports any number (≥2) of ancestral lineages via --kmers NAME:PATH.
"""

import sys
import os
import pysam
import argparse


def load_kmer_set(fasta_file, k=15):
    """Load lineage k-mers from a FASTA file into a set (forced uppercase)."""
    kmer_set = set()
    with open(fasta_file, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                kmer_set.add(line.strip().upper()[:k])
    return kmer_set


def calculate_dynamic_P(read_seq, lineage_dicts, k=15):
    """
    Assign a lineage probability vector to a single read.

    Strategy:
      - Count k-mer hits (forward + reverse-complement) for each lineage.
      - Unique hit  → probability 1.0 to that lineage.
      - Multi-hit   → probability split equally among hit lineages.
      - No hit      → probability split equally among all lineages (uniform prior).

    Parameters
    ----------
    read_seq      : str  – raw read sequence
    lineage_dicts : dict – {lineage_name: set_of_kmers}
    k             : int  – k-mer size (must match the loaded k-mer sets)

    Returns
    -------
    probs  : dict  – {lineage_name: float}
    status : str   – human-readable assignment category
    """
    hits = {name: 0 for name in lineage_dicts}

    trans   = str.maketrans('ATCGN', 'TAGCN')
    fwd_seq = read_seq.upper()
    rev_seq = fwd_seq.translate(trans)[::-1]
    seq_len = len(fwd_seq)

    for i in range(seq_len - k + 1):
        kmer_fwd = fwd_seq[i:i + k]
        kmer_rev = rev_seq[i:i + k]
        for name, kmer_set in lineage_dicts.items():
            if kmer_fwd in kmer_set:
                hits[name] += 1
            if kmer_rev in kmer_set:
                hits[name] += 1

    hit_lineages      = [name for name, count in hits.items() if count > 0]
    num_total_lineages = len(lineage_dicts)
    probs              = {name: 0.0 for name in lineage_dicts}

    if len(hit_lineages) == 1:
        probs[hit_lineages[0]] = 1.0
        status = f"{hit_lineages[0]}_Specific"
    elif len(hit_lineages) > 1:
        val = 1.0 / len(hit_lineages)
        for name in hit_lineages:
            probs[name] = val
        status = f"Conflict_{len(hit_lineages)}"
    else:
        val = 1.0 / num_total_lineages
        for name in lineage_dicts:
            probs[name] = val
        status = "Ambiguous_None"

    return probs, status


def main():
    parser = argparse.ArgumentParser(
        description="LineageQuant Module 02: BAM Lineage P-value Assignment"
    )
    parser.add_argument("-i", "--ids",   required=True,
                        help="Target gene ID list (one ID per line)")
    parser.add_argument("-g", "--gff",   required=True, help="Whole-genome GFF3")
    parser.add_argument("-b", "--bam",   required=True, help="Sorted and indexed BAM file")
    parser.add_argument("--kmers",       nargs='+', required=True, metavar="NAME:PATH",
                        help=(
                            "Ancestral lineage k-mer FASTA files. "
                            "Format: NAME:PATH  (≥2 required)  "
                            "Example: --kmers SubA:subA_k15.fa SubB:subB_k15.fa"
                        ))
    parser.add_argument("-k", "--kmer",  type=int, default=15,
                        help="K-mer size used when loading kmer files (default: 15)")
    parser.add_argument("-o", "--out",   default="P_values.tsv",
                        help="Output file path (default: P_values.tsv)")
    args = parser.parse_args()

    # ---- Step 1: Load ancestral k-mer dictionaries ----
    lineage_dicts = {}
    lineage_names = []
    print("[INFO] Loading ancestral lineage k-mer dictionaries...")
    for item in args.kmers:
        if ':' not in item:
            print(f"[ERROR] Invalid --kmers entry: '{item}'. Expected NAME:PATH.")
            sys.exit(1)
        name, path = item.split(':', 1)
        lineage_names.append(name)
        lineage_dicts[name] = load_kmer_set(path, k=args.kmer)
        print(f"  -> Loaded {name}: {len(lineage_dicts[name])} k-mers from {path}")

    if len(lineage_dicts) < 2:
        print("[ERROR] --kmers requires at least 2 lineages.")
        sys.exit(1)

    # ---- Step 2: Load target gene IDs ----
    with open(args.ids) as f:
        target_ids = {line.strip().split('.')[0] for line in f if line.strip()}
    print(f"[INFO] Loaded {len(target_ids)} target gene IDs.")

    # ---- Step 3: Parse gene coordinates from GFF ----
    gene_coords = {}
    with open(args.gff, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 9 or parts[2] != 'gene':
                continue
            attr = parts[8]
            for tid in target_ids:
                if tid in attr:
                    gene_coords[tid] = (parts[0], int(parts[3]), int(parts[4]))
                    break
    print(f"[INFO] Resolved coordinates for {len(gene_coords)} genes.")

    # ---- Step 4: Scan BAM and assign P-values ----
    bam_file         = pysam.AlignmentFile(args.bam, "rb")
    processed_reads  = set()
    total_fetched    = 0
    total_saved      = 0

    with open(args.out, "w") as out_f:
        p_headers  = [f"P({name})" for name in lineage_names]
        header_str = "Read_Name\tTarget_Gene\t" + "\t".join(p_headers) + "\tLineage_Status\n"
        out_f.write(header_str)

        for gene_id, (chrom, start, end) in gene_coords.items():
            gene_fetched = 0
            for read in bam_file.fetch(chrom, start, end):
                total_fetched += 1
                gene_fetched  += 1

                if read.query_name in processed_reads:
                    continue
                processed_reads.add(read.query_name)

                seq = read.query_sequence
                if not seq:
                    continue

                total_saved += 1
                probs, status = calculate_dynamic_P(seq, lineage_dicts, k=args.kmer)
                p_values_str  = "\t".join([str(probs[name]) for name in lineage_names])
                out_f.write(f"{read.query_name}\t{gene_id}\t{p_values_str}\t{status}\n")

            print(f"  -> {gene_id}: fetched {gene_fetched} reads.")

    bam_file.close()
    print(f"\n[REPORT] Total reads fetched : {total_fetched}")
    print(f"[REPORT] Unique reads written: {total_saved}")
    print(f"[SUCCESS] P-value table written to: {args.out}")


if __name__ == "__main__":
    main()
