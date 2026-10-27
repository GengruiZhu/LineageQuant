#!/bin/bash
# lqcore lrcount vs lineagequant/longread/extract_prior.py on random SAMs (fake pysam / samtools, no HTS libs needed)
#   bash tests/compare_lrcount.sh [N_SEEDS]
set -u
H=$(cd "$(dirname "$0")/.." && pwd); T=$(mktemp -d); N=${1:-8}
LQ=$H/bin/lqcore; EP=${EP:-$H/../package/LineageQuant_v1.0.0/lineagequant/longread/extract_prior.py}   # EP=<1.0.0>/lineagequant/longread/extract_prior.py
mkdir -p "$T/fakebin" && cp "$H"/tests/fakesamtools/* "$T/fakebin/" && chmod +x "$T/fakebin"/* || exit 1   # exec bits may be lost
export PATH=$T/fakebin:$PATH
pass=0; fail=0
for s in $(seq 1 "$N"); do
  for nm in pb sra; do for od in grouped coord; do
    for opt in "" "--as-tol 15" "--min-identity 0.9 --min-aln-len 400" "--zmw-dedup off" "--zmw-dedup on"; do
      c=$T/c$s.$nm.$od.sam; [ -s "$c" ] || python3 "$H/tests/make_lr_case.py" "$s" "$c" $nm $od
      rm -f $T/p.tsv $T/p.leg $T/p.sum $T/c.tsv $T/c.leg $T/c.sum   # never compare the outputs of the previous case
      ran=1
      PYTHONPATH=$H/tests/fakepysam_lr python3 "$EP" -b "$c" -o $T/p.tsv --legacy-out $T/p.leg --summary $T/p.sum --tmpdir $T $opt > $T/p.log 2>&1 || { ran=0; echo "python failed: $c $opt"; tail -3 $T/p.log; }
      "$LQ" lrcount --sam "$c" -o $T/c.tsv --legacy-out $T/c.leg --summary $T/c.sum --tmpdir $T $opt > $T/c.log 2>&1 || { ran=0; echo "lqcore failed: $c $opt"; tail -3 $T/c.log; }
      if [ $ran = 1 ] && [ -s $T/p.tsv ] && cmp -s $T/p.tsv $T/c.tsv && cmp -s $T/p.leg $T/c.leg && cmp -s $T/p.sum $T/c.sum; then pass=$((pass+1));
      else fail=$((fail+1)); echo "DIFF seed=$s $nm $od '$opt'"; diff $T/p.tsv $T/c.tsv | head -4; diff $T/p.sum $T/c.sum | head -6; fi
    done
  done; done
done
echo "lrcount vs extract_prior.py: $pass identical, $fail different"; rm -rf "$T"; [ $fail = 0 ]
