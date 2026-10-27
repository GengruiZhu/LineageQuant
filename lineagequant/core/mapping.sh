#!/bin/bash
# ==============================================================================
# LineageQuant core/mapping.sh -- Step 00: short-read alignment to the transcriptome
#
# Args: <R1.fq.gz> <R2.fq.gz> <Index_prefix> <GFF3(reserved)> <Out_Dir> <Threads>
#       [Aligner: bowtie2 (default) | hisat2] [Max alignments per fragment]
#
# bowtie2 (default since 1.0.0) -- RSEM's settings, plus ambiguous reference
#   bases (N) not penalised:
#     --sensitive --dpad 0 --gbar 99999999 --mp 1,1 --np 0 --n-ceil L,0,0.5
#     --score-min L,0,-0.1 -I 1 -X 1000 --no-mixed --no-discordant -k 200
#   End-to-end and ungapped; reports every equally good copy and near-optimal
#   ones. Step 03 weights the candidates of a fragment by mismatch count.
#
# hisat2 (pre-1.0 behaviour):
#     -k 15 --no-mixed --no-discordant
#
# Expert overrides (environment):
#   LQ_BT2_OPTS     replaces the bowtie2 option string above (except -p/-k/-x/-1/-2)
#   LQ_HISAT_OPTS   replaces "--no-mixed --no-discordant" for hisat2
#   LQ_HISAT_K / LQ_BT2_K   override the alignment cap (the CLI's --max-alignments wins)
# ==============================================================================
set -e
set -o pipefail

if [ "$#" -lt 6 ]; then
    echo "Usage: $0 <R1.fq.gz> <R2.fq.gz> <Index_prefix> <GFF3> <Out_Dir> <Threads> [bowtie2|hisat2] [max_alignments]"
    exit 1
fi

R1=$1; R2=$2; IDX=$3; GFF=$4; OUT_DIR=$5; THREADS=$6   # GFF reserved for compatibility
ALIGNER=${7:-bowtie2}
MAXK=$8
sample_name=$(basename "$R1" ".R1.fq.gz")
bam_file="$OUT_DIR/${sample_name}.sorted.bam"

echo "[INFO] Transcriptome mapping for $sample_name ($ALIGNER)"
A_THREADS=$(( THREADS * 2 / 3 )); S_THREADS=$(( THREADS - A_THREADS ))
[ "$A_THREADS" -lt 2 ] && A_THREADS=2
[ "$S_THREADS" -lt 2 ] && S_THREADS=2
SORT_MEM="1G"

if [ ! -f "${bam_file}.bai" ]; then
    rm -f "${bam_file}.tmp."*.bam
    case "$ALIGNER" in
        bowtie2)
            K=${MAXK:-${LQ_BT2_K:-200}}
            OPTS=${LQ_BT2_OPTS:---sensitive --dpad 0 --gbar 99999999 --mp 1,1 --np 0 --n-ceil L,0,0.5 --score-min L,0,-0.1 -I 1 -X 1000 --no-mixed --no-discordant}
            echo "[EXEC] bowtie2 -k $K $OPTS"
            bowtie2 -p "$A_THREADS" -k "$K" $OPTS -x "$IDX" -1 "$R1" -2 "$R2" 2> "$OUT_DIR/${sample_name}.bowtie2.log" | \
                samtools sort -@ "$S_THREADS" -m "$SORT_MEM" -o "$bam_file" -
            tail -n 6 "$OUT_DIR/${sample_name}.bowtie2.log" | sed 's/^/  | /'
            ;;
        hisat2)
            K=${MAXK:-${LQ_HISAT_K:-15}}
            OPTS=${LQ_HISAT_OPTS:---no-mixed --no-discordant}
            echo "[EXEC] HISAT2 -k $K $OPTS"
            hisat2 -p "$A_THREADS" -k "$K" -x "$IDX" -1 "$R1" -2 "$R2" $OPTS 2> "$OUT_DIR/${sample_name}.hisat2.log" | \
                samtools sort -@ "$S_THREADS" -m "$SORT_MEM" -o "$bam_file" -
            tail -n 6 "$OUT_DIR/${sample_name}.hisat2.log" | sed 's/^/  | /'
            ;;
        *)
            echo "[FATAL] unknown aligner: $ALIGNER (bowtie2 | hisat2)"; exit 1 ;;
    esac
    echo "[EXEC] Indexing BAM"
    samtools index -@ "$THREADS" "$bam_file"
fi
echo "[SUCCESS] Mapping complete: $bam_file"
