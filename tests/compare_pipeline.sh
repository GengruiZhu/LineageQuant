#!/bin/bash
# End to end: `lq run` (C++) vs the pipeline of LineageQuant 1.0.0 (cli.py / router.py, --engine cpp:
# Python for the SNP database, markers, genome/annotation handling and the matrix; lqcore for theta / run).
# Same data (tests/make_pipeline_case.py), four routes; compares the shared reference files, per sample the Theta,
# EM and idxstats tables (the .lqb / .lqres binaries exported to text) and the expression matrix
# (`lq matrix` vs aggregator.py). External tools: tests/fake_tools (+ fake bowtie2 for route 4).
#   bash tests/compare_pipeline.sh bin/lineagequant P100_PKG_DIR [WORKDIR]     (P100_PKG_DIR = .../LineageQuant_v1.0.0/lineagequant)
set -u
LQ=$(readlink -f "${1:?lineagequant binary}"); P=$(readlink -f "${2:?LineageQuant 1.0.0 package dir}"); W=${3:-/tmp/lq_pipeline_test}
LQC=$(dirname "$LQ")/lqcore
[ -x "$LQC" ] || { echo "needs $LQC (the engine name, built next to lineagequant)"; exit 2; }
HERE=$(cd "$(dirname "$0")" && pwd)
rm -rf "$W"; mkdir -p "$W/bin"; fail=0
cp "$HERE"/fake_tools/* "$W/bin/" && chmod +x "$W/bin"/* || exit 1   # copies: a checkout may have lost the exec bits
printf '#!/bin/bash\nexec python3 %s "$@"\n' "$HERE/fake_bowtie2.py" > "$W/bin/bowtie2"
printf '#!/bin/bash\nfor a; do p=$a; done\necho fake > "$p.1.bt2"\n' > "$W/bin/bowtie2-build"
chmod +x "$W/bin/bowtie2" "$W/bin/bowtie2-build"
export PATH=$W/bin:$PATH
export PYTHONPATH=$HERE/fakepysam:$HERE/fakebio:$HERE/fakecyvcf2:$(dirname "$P")${PYTHONPATH:+:$PYTHONPATH}
export PYTHONDONTWRITEBYTECODE=1
ok () { printf "  %-62s %s\n" "$1" "$2"; [ "${2%% *}" = identical ] || fail=1; }
same () {
  if [ ! -e "$1" ] && [ ! -e "$2" ]; then echo "identical (absent on both sides)"
  elif cmp -s "$1" "$2"; then echo "identical ($(wc -l < "$1") lines)"; else echo "DIFF: diff $1 $2"; fi; }

python3 "$HERE/make_pipeline_case.py" "$W/case" 5 15000 > /dev/null || exit 1
R=$W/case/ref
# library SAM for the fake bowtie2 (route 4): lqmap alignments of both samples
"$LQ" map --fasta "$R/tx.fa" --r1 "$W/case/fq/S1_R1.fq.gz" --r2 "$W/case/fq/S1_R2.fq.gz" --threads 4 > "$W/lib1.sam" 2> /dev/null
"$LQ" map --fasta "$R/tx.fa" --r1 "$W/case/fq/S2_1.fq.gz" --r2 "$W/case/fq/S2_2.fq.gz" --threads 4 2> /dev/null | grep -v '^@' > "$W/lib2.sam"
cat "$W/lib1.sam" "$W/lib2.sam" > "$W/lib.sam"; export FAKE_BT2_SAM=$W/lib.sam

BASE=(--fasta "$R/tx.fa" --gff "$R/genome.gff3" --clusters "$R/clusters.tsv" --ancestry "$R/ancestry.tsv" --orphan "$R/orphans.tsv" --target-ids "$R/target_ids.txt")
route () {   # name "python-only args" -- common args...
  local r=$1; shift
  echo "== $r"
  mkdir -p "$W/$r/dpy" "$W/$r/dcpp"
  for f in "$W"/case/fq/*; do ln -s "$f" "$W/$r/dpy/"; ln -s "$f" "$W/$r/dcpp/"; done
  local t0=$(date +%s.%N)
  python3 -m lineagequant.cli run "${BASE[@]}" "$@" --engine cpp --lqcore "$LQC" -d "$W/$r/dpy" -o "$W/$r/opy" -t 4 > "$W/$r/py.log" 2>&1 \
    || { echo "  1.0.0 pipeline failed"; tail -5 "$W/$r/py.log"; fail=1; return; }
  local t1=$(date +%s.%N)
  "$LQ" run "${BASE[@]}" "$@" --theta-sum plain -d "$W/$r/dcpp" -o "$W/$r/ocpp" -t 4 > "$W/$r/cpp.log" 2>&1 \
    || { echo "  lq run failed"; tail -5 "$W/$r/cpp.log"; fail=1; return; }
  local t2=$(date +%s.%N)
  for f in SNP_DB_msa.lq.tsv SNP_DB_population.lq.tsv lineage_markers/A.fasta lineage_markers/B.fasta \
           ancestor_transcripts/A_transcript.fa ancestor_transcripts/B_transcript.fa; do
    [ -e "$W/$r/dpy/$f" ] || [ -e "$W/$r/dcpp/$f" ] && ok "data/$f" "$(same "$W/$r/dpy/$f" "$W/$r/dcpp/$f")"
  done
  for s in S1 S2; do
    for f in Global_Theta_Matrix.tsv EM_Final_Expression.tsv $s.idxstats.tsv; do
      ok "$s/$f" "$(same "$W/$r/opy/$s/$f" "$W/$r/ocpp/$s/$f")"
    done
  done
  "$LQ" matrix -r "$W/$r/ocpp" -o "$R/orphans.tsv" --out "$W/$r/cpp_matrix.tsv" > /dev/null 2>&1 || fail=1
  ok "matrix (from .lqres) vs LineageQuant_Final_Merged_TPM.tsv" "$(same "$W/$r/opy/LineageQuant_Final_Merged_TPM.tsv" "$W/$r/cpp_matrix.tsv")"
  printf "  wall: 1.0.0 %.1fs (incl. Step 04) | lq run %.1fs\n" "$(echo "$t1 - $t0" | bc)" "$(echo "$t2 - $t1" | bc)"
}
route allo_ancestor_msa --ploidy allo --ancestor "A:$R/progA.fa,$R/progA.gff" --ancestor "B:$R/progB1.fa,$R/progB1.gff" \
      --ancestor "B:$R/progB2.fa,$R/progB2.gff" --msa-snp --aligner lqmap
route auto --ploidy auto --aligner lqmap
M=$W/allo_ancestor_msa/dpy/lineage_markers
route allo_markers_vcf --ploidy allo --markers "A:$M/A.fasta" --markers "B:$M/B.fasta" --vcf "$R/pop.vcf.gz" \
      --sample-lists "A:$R/A.samples" "B:$R/B.samples" --aligner lqmap --em-accel squarem --nm-weight 0.02 --p-weight-conflict 0.5
route allo_markers_bowtie2 --ploidy allo --markers "A:$M/A.fasta" --markers "B:$M/B.fasta" --aligner bowtie2 --max-alignments 50 \
      --no-pool-identical --eff-len len --shrinkage 5 --theta-prior-kmer 2

# resume: a second run skips every step and changes nothing
cp "$W/auto/ocpp/S1/S1.lqres" "$W/s1.before"
"$LQ" run "${BASE[@]}" --ploidy auto --aligner lqmap --theta-sum plain -d "$W/auto/dcpp" -o "$W/auto/ocpp" -t 4 > "$W/auto/cpp2.log" 2>&1 || fail=1
n_run=$(grep -c "\[RUN\]" "$W/auto/cpp2.log"); n_skip=$(grep -c "\[SKIP\]" "$W/auto/cpp2.log")
ok "resume: second run ($n_skip skipped, $n_run run)" "$([ "$n_run" = 0 ] && cmp -s "$W/s1.before" "$W/auto/ocpp/S1/S1.lqres" && echo identical || echo DIFF)"
n_txt=$(grep -c "\[TEXT\]" "$W/auto/cpp2.log")
ok "resume: text files up to date, not rewritten ($n_txt exports)" "$([ "$n_txt" = 0 ] && echo identical || echo DIFF)"
# --no-text: binary only, same .lqres
"$LQ" run "${BASE[@]}" --ploidy auto --aligner lqmap --theta-sum plain --no-text -d "$W/auto/dcpp" -o "$W/auto/ocpp_bin" -t 4 > "$W/auto/cpp4.log" 2>&1 || fail=1
ok "--no-text: binary only, same results" "$([ ! -e "$W/auto/ocpp_bin/S1/EM_Final_Expression.tsv" ] && cmp -s "$W/auto/ocpp_bin/S1/S1.lqres" "$W/auto/ocpp/S1/S1.lqres" && echo identical || echo DIFF)"
# switching the aligner on an existing output directory is refused
"$LQ" run "${BASE[@]}" --ploidy auto --aligner bowtie2 -d "$W/auto/dcpp" -o "$W/auto/ocpp" -t 2 > "$W/auto/cpp3.log" 2>&1
ok "aligner switch refused" "$(grep -q "was aligned with lqmap" "$W/auto/cpp3.log" && echo identical || echo DIFF)"

echo; [ $fail = 0 ] && echo "ALL PIPELINE CHECKS PASSED" || echo "SOME PIPELINE CHECKS FAILED"
exit $fail
