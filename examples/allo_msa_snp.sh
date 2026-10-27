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

PHASED=${PHASED:-phased_transcriptome.fasta}  # phased allele copies

# Route: allopolyploid, MSA SNP lineage assignment (no long-read calibration).
lineagequant run \
  --fasta "$FASTA" --gff "$GFF" \
  --clusters "$CLUSTERS" --ancestry "$ANCESTRY" --orphan "$ORPHAN" \
  --target-ids "$TARGET_IDS" \
  --ploidy allo --phased-fasta "$PHASED" \
  --data-dir "$DATA_DIR" --out-dir "$OUT_DIR" --threads "$THREADS"
