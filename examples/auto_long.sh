#!/bin/bash
# ==============================================================================
# LineageQuant 1.1.0 example (bin/ of the build on PATH). All paths are placeholders -- edit the variables below.
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

LONGREAD=${LONGREAD:-longreads.bam}  # long-read BAM (minimap2 -a -N 50 -p 0.95 -> transcriptome; any sort order)

# Route A2: autopolyploid, long-read anchored Theta.
lq run \
  --fasta "$FASTA" --gff "$GFF" \
  --clusters "$CLUSTERS" --ancestry "$ANCESTRY" --orphan "$ORPHAN" \
  --target-ids "$TARGET_IDS" \
  --ploidy auto --longread "$LONGREAD" \
  --data-dir "$DATA_DIR" --out-dir "$OUT_DIR" --threads "$THREADS"

# expression matrix (outside the run): <OUT_DIR>/LineageQuant_Final_Merged_TPM.tsv
lq matrix -r "$OUT_DIR" -o "$ORPHAN"
