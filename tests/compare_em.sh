#!/bin/bash
# Compare `lqcore em` with lineagequant/core/em_global.py on synthetic cases.
#   bash tests/compare_em.sh PY_EM_GLOBAL LQCORE [WORKDIR]
# Every configuration must give byte-identical EM_Final_Expression tables, and lqcore must give
# byte-identical output for 1, 2 and 4 threads.
set -u
PY=${1:?path to em_global.py}; LQ=${2:?path to lqcore}; W=${3:-/tmp/lqcore_em_test}
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$W"; fail=0; n=0
run_case () {   # tag seed nfrag extra-args...
  local tag=$1 seed=$2 nf=$3; shift 3
  local d=$W/case_$seed
  [ -s "$d/case.sam" ] || python3 "$HERE/make_em_case.py" "$d" "$seed" "$nf" > /dev/null
  local base=(-a "$d/ancestry.tsv" -t "$d/theta.tsv" -p "$d/pvals.tsv" --fasta "$d/ref.fa")
  python3 "$PY" "${base[@]}" --bam "$d/case.sam" -o "$d/py_$tag.tsv" "$@" > "$d/py_$tag.log" 2>&1 || { echo "  python failed: $tag"; tail -3 "$d/py_$tag.log"; fail=1; return; }
  for th in 1 2 4; do
    "$LQ" em "${base[@]}" --sam "$d/case.sam" -o "$d/cpp_${tag}_t$th.tsv" --threads $th "$@" > "$d/cpp_${tag}_t$th.log" 2>&1 \
      || { echo "  lqcore failed: $tag"; tail -3 "$d/cpp_${tag}_t$th.log"; fail=1; return; }
  done
  n=$((n+1))
  local st="identical"
  cmp -s "$d/py_$tag.tsv" "$d/cpp_${tag}_t1.tsv" || st="DIFF ($(diff <(cut -f1-6 "$d/py_$tag.tsv") <(cut -f1-6 "$d/cpp_${tag}_t1.tsv") | grep -c '^<') rows)"
  local thr="threads 1=2=4"
  { cmp -s "$d/cpp_${tag}_t1.tsv" "$d/cpp_${tag}_t2.tsv" && cmp -s "$d/cpp_${tag}_t1.tsv" "$d/cpp_${tag}_t4.tsv"; } || thr="THREADS DIFFER"
  [ "$st" = identical ] && [ "$thr" = "threads 1=2=4" ] || fail=1
  local it=$(grep -oE "EM (converged|stopped)[^,]*" "$d/py_$tag.log" | head -1)
  printf "  %-34s seed %-3s %-12s %-15s | %s\n" "$tag" "$seed" "$st" "$thr" "$it"
}
for seed in ${SEEDS:-1 2 3}; do
  run_case default            $seed 20000
  run_case uniform            $seed 20000 --init uniform
  run_case nm_off             $seed 20000 --nm-weight off
  run_case nm_fixed           $seed 20000 --nm-weight 0.01
  run_case no_pool            $seed 20000 --no-pool
  run_case efflen_len         $seed 20000 --eff-len len
  run_case efflen_off         $seed 20000 --eff-len off
  run_case lr_prior           $seed 20000 --theta-prior-lr 1.0 --theta-prior-cap 25
  run_case kmer_prior         $seed 20000 --theta-prior-kmer 2.5
  run_case p_weights          $seed 20000 --p-weight-conflict 0.5 --p-weight-none 0.25 --p-alpha 0.01 --p-floor 0.01
  run_case maxiter_cap        $seed 20000 --max-iter 7 --min-iter 3
  run_case fraglen            $seed 20000 --frag-len 250 --min-eff-len 50
done
echo; [ $fail = 0 ] && echo "ALL $n CONFIGURATIONS IDENTICAL" || echo "SOME CONFIGURATIONS DIFFER"
exit $fail
