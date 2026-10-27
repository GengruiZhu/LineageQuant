#!/bin/bash
# Binary outputs and the expression matrix:
#   1. quant/em -o S.lqres, then `lq export` == the text files written directly (EM, idxstats, efflen)
#   2. `lq matrix` == core/aggregator.py of LineageQuant 1.0.0 (--mode allo, idxstats) and of the
#      Python router of 1.1.0 (--tpm-length full / eff), on the text outputs and on the .lqres files
#   bash tests/compare_matrix.sh LINEAGEQUANT P100_PKG_DIR [AGGREGATOR_110] [WORKDIR]
#     P100_PKG_DIR = .../LineageQuant_v1.0.0/lineagequant ; AGGREGATOR_110 (default tests/ref_py/aggregator_110.py) =
#     the 1.0.0 aggregator plus --tpm-length eff
set -u
LQ=${1:?lineagequant binary}; P100=${2:?LineageQuant 1.0.0 package dir}
HERE=$(cd "$(dirname "$0")" && pwd)
AGG110=${3:-$HERE/ref_py/aggregator_110.py}; W=${4:-/tmp/lq_matrix_test}
export PYTHONPATH=$HERE/fakepysam:$HERE/fakebio${PYTHONPATH:+:$PYTHONPATH}
export PYTHONDONTWRITEBYTECODE=1
rm -rf "$W"; mkdir -p "$W"; fail=0
ok () { printf "  %-60s %s\n" "$1" "$2"; [ "$2" = identical ] || fail=1; }
same () { cmp -s "$1" "$2" && echo identical || echo "DIFF ($1 vs $2)"; }

# three references / read sets (different gene sets: the matrix takes the union)
for s in 1 2 3; do python3 "$HERE/make_em_case.py" "$W/case_$s" $s 20000 > /dev/null; done
declare -A OPT=([1]="" [2]="--eff-len fld" [3]="--eff-len fld --gc-bias")
mkdir -p "$W/txt" "$W/bin"
for s in 1 2 3; do
  d=$W/case_$s; nm=S$s
  base=(-a "$d/ancestry.tsv" -t "$d/theta.tsv" -p "$d/pvals.tsv" --fasta "$d/ref.fa" --sam "$d/case.sam" --threads 2)
  mkdir -p "$W/txt/$nm" "$W/bin/$nm"
  "$LQ" em "${base[@]}" ${OPT[$s]} -o "$W/txt/$nm/EM_Final_Expression.tsv" --idxstats-out "$W/txt/$nm/$nm.idxstats.tsv" \
       --eff-out "$W/txt/$nm/$nm.efflen.tsv" > "$W/txt/$nm.log" 2>&1 || { echo "em (text) failed: $nm"; tail -3 "$W/txt/$nm.log"; exit 1; }
  "$LQ" em "${base[@]}" ${OPT[$s]} -o "$W/bin/$nm/$nm.lqres" > "$W/bin/$nm.log" 2>&1 || { echo "em (bin) failed: $nm"; tail -3 "$W/bin/$nm.log"; exit 1; }
  "$LQ" export "$W/bin/$nm/$nm.lqres" -o "$W/exp/$nm" > /dev/null || exit 1
  for f in EM_Final_Expression.tsv $nm.idxstats.tsv $nm.efflen.tsv; do ok "export $nm/$f" "$(same "$W/txt/$nm/$f" "$W/exp/$nm/$f")"; done
done
# a sample without EM output (text only): idxstats of S1
mkdir -p "$W/txt/sub/S4"; cp "$W/txt/S1/S1.idxstats.tsv" "$W/txt/sub/S4/S4.idxstats.tsv"
# orphans: OrphanGrp<TAB>gene (the ORPHAN rows of S1) + one-column lines
awk -F'\t' 'NR>1 && $1=="ORPHAN"{n++; if (n%3) print "OrphanGrp_" n "\t" $2; else print $2}' "$W/txt/S1/EM_Final_Expression.tsv" > "$W/orphans.tsv"
echo "# comment" >> "$W/orphans.tsv"

for mode in allo idxstats; do
  python3 "$P100/core/aggregator.py" -r "$W/txt" -o "$W/orphans.tsv" --out "$W/py100_$mode.tsv" --mode $mode > "$W/py100_$mode.log" 2>&1 \
    || { echo "aggregator 1.0.0 failed"; tail -3 "$W/py100_$mode.log"; fail=1; }
  "$LQ" matrix -r "$W/txt" -o "$W/orphans.tsv" --out "$W/cpp_txt_$mode.tsv" --mode $mode > /dev/null 2>&1 || fail=1
  ok "matrix (text inputs) --mode $mode vs aggregator 1.0.0" "$(same "$W/py100_$mode.tsv" "$W/cpp_txt_$mode.tsv")"
done
for tl in full eff; do
  python3 "$AGG110" -r "$W/txt" -o "$W/orphans.tsv" --out "$W/py110_$tl.tsv" --tpm-length $tl > "$W/py110_$tl.log" 2>&1 \
    || { echo "aggregator (1.1.0 Python) failed"; tail -3 "$W/py110_$tl.log"; fail=1; }
  "$LQ" matrix -r "$W/txt" -o "$W/orphans.tsv" --out "$W/cpp_txt_$tl.tsv" --tpm-length $tl > /dev/null 2>&1 || fail=1
  ok "matrix (text inputs) --tpm-length $tl vs aggregator" "$(same "$W/py110_$tl.tsv" "$W/cpp_txt_$tl.tsv")"
done
# .lqres inputs (+ the EM-less S4 as text): same matrix as the text inputs
mkdir -p "$W/bin/sub/S4"; cp "$W/txt/sub/S4/S4.idxstats.tsv" "$W/bin/sub/S4/"
for tl in full eff; do
  "$LQ" matrix -r "$W/bin" -o "$W/orphans.tsv" --out "$W/cpp_bin_$tl.tsv" --tpm-length $tl > /dev/null 2>&1 || fail=1
  ok "matrix (.lqres inputs) --tpm-length $tl vs aggregator" "$(same "$W/py110_$tl.tsv" "$W/cpp_bin_$tl.tsv")"
done
"$LQ" matrix -r "$W/bin" -o "$W/orphans.tsv" --out "$W/cpp_bin_idxstats.tsv" --mode idxstats > /dev/null 2>&1 || fail=1
ok "matrix (.lqres inputs) --mode idxstats vs aggregator 1.0.0" "$(same "$W/py100_idxstats.tsv" "$W/cpp_bin_idxstats.tsv")"
# many genes (numpy sums over > 8192 values in blocks)
python3 "$HERE/make_matrix_case.py" "$W/big" 30000 3 7
for tl in full eff; do
  python3 "$AGG110" -r "$W/big/res" -o "$W/big/orphans.tsv" --out "$W/big/py_$tl.tsv" --tpm-length $tl > /dev/null 2>&1 || fail=1
  "$LQ" matrix -r "$W/big/res" -o "$W/big/orphans.tsv" --out "$W/big/cpp_$tl.tsv" --tpm-length $tl > /dev/null 2>&1 || fail=1
  ok "matrix 30,000 genes --tpm-length $tl vs aggregator" "$(same "$W/big/py_$tl.tsv" "$W/big/cpp_$tl.tsv")"
done
printf "  (matrix: %s genes x %s samples)\n" "$(($(wc -l < "$W/cpp_txt_allo.tsv") - 1))" 4

echo; [ $fail = 0 ] && echo "ALL IDENTICAL" || echo "SOME DIFFER"
exit $fail
