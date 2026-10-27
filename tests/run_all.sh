#!/bin/bash
# All equivalence tests of LineageQuant 1.1.0 against the Python steps of LineageQuant 1.0.0 (synthetic data, minutes).
#   bash tests/run_all.sh /path/to/LineageQuant_v1.0.0/lineagequant      (after bash build.sh)
# Needs python3 with numpy and pandas; Biopython / pysam / cyvcf2 and the external tools are replaced by the stand-ins
# in tests/ (both sides use the same ones). Each test prints its own table; the summary is at the end.
set -u
P100=$(readlink -f "${1:?path to the lineagequant/ directory of LineageQuant 1.0.0}")
H=$(cd "$(dirname "$0")/.." && pwd); T=$H/tests; LQ=$H/bin/lq; LQC=$H/bin/lqcore
[ -x "$LQ" ] && [ -x "$LQC" ] || { echo "build first: bash build.sh"; exit 2; }
[ "$("$LQ" version)" = "lq 1.1.0" ] && [ "$("$H/bin/lineagequant" version)" = "lineagequant 1.1.0" ] && [ "$("$LQC" version)" = "lqcore 1.1.0" ] \
  || { echo "bin/lq, bin/lineagequant, bin/lqcore must all be 1.1.0: bash build.sh"; exit 2; }
[ -s "$P100/core/em_global.py" ] || { echo "$P100 is not the lineagequant/ directory of LineageQuant 1.0.0"; exit 2; }
python3 -c "import numpy, pandas" || { echo "python3 with numpy and pandas needed"; exit 2; }
export PYTHONDONTWRITEBYTECODE=1
W=$(mktemp -d); trap 'rm -rf "$W"' EXIT
declare -a NAME=() RES=()   # : must be empty: an exported NAME=... would otherwise become element 0
t () {   # name command...
  local n=$1; shift
  echo; echo "################ $n"
  local t0=$(date +%s)
  if "$@"; then RES+=("ok"); else RES+=("FAILED"); fi
  NAME+=("$n ($(( $(date +%s) - t0 ))s)")
}
PP=$T/fakepysam:$T/fakebio:$(dirname "$P100")
t "EM (36 configurations)"            env PYTHONPATH=$PP bash "$T/compare_em.sh" "$P100/core/em_global.py" "$LQC" "$W/em"
t "lineage chain (k-mer, SNP, PoE, EM)" env PYTHONPATH=$PP bash "$T/compare_lineage.sh" "$P100" "$LQC"
t "Theta (36 configurations + binary)" bash "$T/compare_theta.sh" "$P100" "$LQC"
t "long-read counts (80)"             env EP=$P100/longread/extract_prior.py PYTHONPATH=$T/fakepysam_lr bash "$T/compare_lrcount.sh" 4
t "lqmap vs exhaustive alignment"     bash "$T/compare_map.sh" "$LQC"
t "1.1.0 options (fld, GC, marker clean)" bash "$T/compare_new.sh" "$LQC"
t "binary outputs + matrix"           bash "$T/compare_matrix.sh" "$LQ" "$P100" "$T/ref_py/aggregator_110.py" "$W/matrix"
t "reference preparation"             bash "$T/compare_prep.sh" "$LQ" "$P100" "$W/prep"
t "whole pipeline (4 routes)"         bash "$T/compare_pipeline.sh" "$LQ" "$P100" "$W/pipeline"
echo; echo "================ summary"
fail=0
for i in "${!NAME[@]}"; do printf "  %-50s %s\n" "${NAME[$i]}" "${RES[$i]}"; [ "${RES[$i]}" = ok ] || fail=1; done
[ $fail = 0 ] && echo "ALL TESTS PASSED" || echo "SOME TESTS FAILED"
exit $fail
