#!/bin/bash
# ==============================================================================
# LineageQuant Module 00: Alignment & Baseline Counting
# Input : R1/R2 FASTQ, HISAT2 index, GFF3, output directory, thread count
# Output: Sorted BAM + featureCounts count table
# ==============================================================================
set -e

if [ "$#" -ne 6 ]; then
    echo "Usage: $0 <R1.fq.gz> <R2.fq.gz> <Index_Prefix> <GFF3> <Out_Dir> <Threads>"
    exit 1
fi

R1=$1
R2=$2
IDX=$3
GFF=$4
OUT_DIR=$5
THREADS=$6

# Derive sample name by stripping the .R1.fq.gz suffix
sample_name=$(basename "$R1" ".R1.fq.gz")
bam_file="$OUT_DIR/${sample_name}.sorted.bam"
fc_out="$OUT_DIR/${sample_name}_featureCounts.txt"

echo "[MODULE 00] Starting alignment pipeline for: $sample_name"

# ---- Step 1: HISAT2 alignment and BAM sorting ----
if [ ! -f "${bam_file}.bai" ]; then
    echo "  -> Running HISAT2 alignment and samtools sort..."
    hisat2 -p "$THREADS" -x "$IDX" -1 "$R1" -2 "$R2" | \
    samtools view -bS - | \
    samtools sort -@ "$THREADS" -o "$bam_file" -
    samtools index "$bam_file"
    echo "  -> BAM file ready."
else
    echo "  -> BAM file already exists, skipping alignment."
fi

# ---- Step 2: featureCounts global baseline count ----
# Note: uses -t gene -g ID for GFF3 files where genes have an ID= attribute.
# Adjust -t (feature type) and -g (meta-feature attribute) if your GFF3
# uses a different structure (e.g. -t exon -g Parent).
if [ ! -f "$fc_out" ]; then
    echo "  -> Running featureCounts..."
    featureCounts -T "$THREADS" -p -t gene -g ID -a "$GFF" -o "$fc_out" "$bam_file"
    echo "  -> featureCounts complete."
else
    echo "  -> featureCounts output already exists, skipping."
fi

echo "[MODULE 00] Alignment pipeline complete for: $sample_name"
