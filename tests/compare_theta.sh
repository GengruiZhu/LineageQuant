#!/bin/bash
# `lqcore theta` vs theta/global_theta.py on synthetic cases (tests/make_theta_case.py), byte for byte;
# also theta --binary followed by `export` (1.1.0).
#   bash tests/compare_theta.sh /path/to/lineagequant bin/lqcore
# Biopython is replaced by tests/fakebio when it is not installed (same FASTA rules).
# Python >= 3.12 sums floats with Neumaier compensation; older interpreters are run through py312sum_run.py
# for the neumaier configurations, natively for the naive ones.
set -o pipefail
PKG=$(readlink -f "$1"); LQ=$(readlink -f "$2"); T=$(cd "$(dirname "$0")" && pwd)
GT=$PKG/theta/global_theta.py
[ -s "$GT" ] && [ -x "$LQ" ] || { echo "usage: $0 LINEAGEQUANT_PKG_DIR LQCORE"; exit 2; }
python3 -c "import Bio" 2>/dev/null || export PYTHONPATH=$T/fakebio${PYTHONPATH:+:$PYTHONPATH}
W=$(mktemp -d); trap 'rm -rf "$W"' EXIT
fail=0
run_cfg () {   # name py_sum extra-args...
  local name=$1 ps=$2; shift 2
  local d=$W/$name; mkdir -p "$d/py" "$d/cpp" "$d/idx"
  local A=(-c "$C/clusters.txt" -f "$C/tx.fa" -r1 "$C/R1.fq.gz" -r2 "$C/R2.fq.gz" "$@")
  if [ "$ps" = neumaier ]; then python3 "$T/py312sum_run.py" "$GT" "${A[@]}" -t 2 -o "$d/py" --index_dir "$d/py" > "$d/py.log" 2>&1
  else python3 "$GT" "${A[@]}" -t 2 -o "$d/py" --index_dir "$d/py" > "$d/py.log" 2>&1; fi || { echo "  python failed ($name)"; tail "$d/py.log"; fail=1; return; }
  for th in 1 4; do
    rm -rf "$d/cpp"; mkdir -p "$d/cpp"
    HC=(); [ $th = 4 ] && HC=(--hits-cache "$C/hits.bin")    # threads 4: k-mer hits cache shared by the configurations of a case
    "$LQ" theta "${A[@]}" -t $th -o "$d/cpp" --index-cache "$d/idx/i.bin" --py-sum "$ps" "${HC[@]}" > "$d/cpp$th.log" 2>&1 || { echo "  lqcore failed ($name)"; tail "$d/cpp$th.log"; fail=1; return; }
    if cmp -s "$d/py/Global_Theta_Matrix.tsv" "$d/cpp/Global_Theta_Matrix.tsv"; then
      printf "  %-34s threads %d  %s  identical  (%s rows)\n" "$name" $th "$(grep -q 'index loaded from' "$d/cpp$th.log" && echo 'cache' || echo 'build')$(grep -q 'hits loaded' "$d/cpp$th.log" && echo '+hits' || true)" \
        "$(($(wc -l < "$d/py/Global_Theta_Matrix.tsv") - 1))"
    else
      echo "  $name threads $th: DIFFERENT"; diff "$d/py/Global_Theta_Matrix.tsv" "$d/cpp/Global_Theta_Matrix.tsv" | head -6; fail=1
    fi
  done
  # 1.1.0: binary Theta (--binary) -> `export` gives the same TSV
  rm -rf "$d/bin"; mkdir -p "$d/bin"
  "$LQ" theta "${A[@]}" -t 2 -o "$d/bin" --index-cache "$d/idx/i.bin" --py-sum "$ps" --binary > "$d/bin.log" 2>&1 \
    && "$LQ" export "$d/bin/Global_Theta_Matrix.lqb" -o "$d/bin/exported.tsv" > /dev/null 2>&1 \
    && cmp -s "$d/py/Global_Theta_Matrix.tsv" "$d/bin/exported.tsv" \
    && printf "  %-34s --binary + export identical\n" "$name" || { echo "  $name: --binary / export DIFFERENT"; fail=1; }
}
for seed in 1 2 3; do
  C=$W/case$seed; python3 "$T/make_theta_case.py" "$C" $seed > /dev/null || exit 1
  echo "== case seed $seed"
  run_cfg "s$seed-paired-naive" naive
  run_cfg "s$seed-paired-neumaier" neumaier
  run_cfg "s$seed-adaptive-neumaier" neumaier --adaptive_shrinkage --shrinkage 7.5
  run_cfg "s$seed-single-end" neumaier --no_paired_bonus
  run_cfg "s$seed-hybrid-lr" neumaier --longread_prior "$C/lr.tsv" --min_lr_depth 3
  run_cfg "s$seed-k15-naive" naive -k 15
done
echo
[ $fail = 0 ] && echo "ALL THETA CHECKS PASSED" || { echo "THETA CHECKS FAILED"; exit 1; }
