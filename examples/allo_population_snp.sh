#!/bin/bash
# ==============================================================================
# LineageQuant example. All paths are placeholders -- edit the variables below.
# No account-specific or absolute paths are hard-coded.
# ==============================================================================
set -euo pipefail

THREADS=${THREADS:-16}
DATA_DIR=${DATA_DIR:-./data}        # directory holding paired FASTQ files
OUT_DIR=${OUT_DIR:-./results}

# --- shared reference inputs (edit these) ---
FASTA=${FASTA:-transcriptome.fasta} # transcriptome (mRNA) FASTA
GFF=${GFF:-genome.gff3}             # reference GFF3
CLUSTERS=${CLUSTERS:-clusters.tsv}  # standardized allele clusters
ANCESTRY=${ANCESTRY:-ancestry.tsv}  # GeneID<TAB>Lineage
ORPHAN=${ORPHAN:-orphans.tsv}       # single-copy gene list
TARGET_IDS=${TARGET_IDS:-target_gene_ids.txt}

VCF=${VCF:-population.vcf.gz}
SAMP_A=${SAMP_A:-samples_A.txt}   # sample names for lineage A (one per line)
SAMP_D=${SAMP_D:-samples_D.txt}

# Route B2: allopolyploid, population VCF SNP lineage assignment.
lineagequant run \
  --fasta "$FASTA" --gff "$GFF" \
  --clusters "$CLUSTERS" --ancestry "$ANCESTRY" --orphan "$ORPHAN" \
  --target-ids "$TARGET_IDS" \
  --ploidy allo \
  --vcf "$VCF" --sample-lists "A:$SAMP_A" "D:$SAMP_D" \
  --data-dir "$DATA_DIR" --out-dir "$OUT_DIR" --threads "$THREADS"
