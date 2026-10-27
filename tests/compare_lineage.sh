#!/bin/bash
# Compare `lqcore run` (lineage evidence + EM in one pass) with the Python pipeline
#   snipe_kmer.py -> snipe_snp.py -> merge_poe.py -> em_global.py   (pysam replaced by tests/fakepysam)
#   bash tests/compare_lineage.sh PKG_DIR(lineagequant) LQCORE [WORKDIR]
# A  every transcript is a lineage target: k-mer rows (per fragment) and SNP rows (per fragment x transcript)
#    must be identical to the Python rows
# B  one target transcript per cluster (one SNP row per fragment, so the Python choice is unambiguous):
#    fused P table and EM table must be identical to the Python chain
# C  read chunks: 3 parts merged with `lqcore em --parts` vs the single run; threads 1 vs 4
set -u
PKG=${1:?lineagequant package dir}; LQ=${2:?lqcore}; W=${3:-/tmp/lq_lineage_test}
HERE=$(cd "$(dirname "$0")" && pwd); FP=$HERE/fakepysam
mkdir -p "$W"; fail=0
cmp_rows () {   # kind python_file cpp_file
python3 - "$@" <<'PY'
import sys
from collections import Counter
kind, pf, cf = sys.argv[1:4]
def rows(p):
    with open(p) as f:
        f.readline()
        return [l.rstrip("\n").split("\t") for l in f if l.strip()]
P, C = rows(pf), rows(cf)
if kind == "kmer":          # python: q gene P.. status hits AS  ->  q P.. status hits
    d = {}
    bad = 0
    for r in P:
        v = tuple(r[2:-1])
        if r[0] in d and d[r[0]] != v: bad += 1
        d[r[0]] = v
    c = {r[0]: tuple(r[1:]) for r in C}
    diff = [q for q in set(d) | set(c) if d.get(q) != c.get(q)]
    print(f"    k-mer: python {len(d):,} fragments (rows disagreeing within a fragment: {bad}), lqcore {len(c):,}, different {len(diff):,}")
    for q in diff[:3]: print(f"      {q}: py {d.get(q)} | cpp {c.get(q)}")
    sys.exit(1 if diff or bad else 0)
if kind == "snp":           # q gene P.. status n
    a, b = Counter(tuple(r) for r in P), Counter(tuple(r) for r in C)
    d1, d2 = a - b, b - a
    print(f"    SNP: python {sum(a.values()):,} rows, lqcore {sum(b.values()):,}, only-python {sum(d1.values()):,}, only-lqcore {sum(d2.values()):,}")
    for r in list(d1)[:3]: print("      py  ", r)
    for r in list(d2)[:3]: print("      cpp ", r)
    sys.exit(1 if d1 or d2 else 0)
if kind == "merged":        # q gene P.. status
    d = {r[0]: tuple(r[2:]) for r in P}
    c = {r[0]: tuple(r[2:]) for r in C}
    diff = [q for q in set(d) | set(c) if d.get(q) != c.get(q)]
    print(f"    fused P: python {len(d):,} fragments, lqcore {len(c):,}, different {len(diff):,}")
    for q in diff[:3]: print(f"      {q}: py {d.get(q)} | cpp {c.get(q)}")
    sys.exit(1 if diff else 0)
PY
}
maxdiff () {
python3 - "$1" "$2" <<'PY'
import sys
a = [l.split("\t") for l in open(sys.argv[1])][1:]; b = [l.split("\t") for l in open(sys.argv[2])][1:]
assert len(a) == len(b) and all(x[:3] == y[:3] for x, y in zip(a, b))
d = max(abs(float(x[4]) - float(y[4])) for x, y in zip(a, b))
n = sum(x != y for x, y in zip(a, b))
print(f"    rows differing {n} of {len(a)}, max |EM_Allocated_Reads diff| {d:.4f}")
sys.exit(0 if d <= 0.02 else 1)
PY
}
for mode in canonical noncanonical; do
  d=$W/$mode; mkdir -p "$d"
  extra=""; [ $mode = noncanonical ] && extra="--noncanonical"
  [ -s "$d/case.sam" ] || python3 "$HERE/make_lineage_case.py" "$d" 11 6000 $extra
  python3 "$PKG/lineage_snp/export_snp_db.py" --snp_db "$d/snp_db.pkl" --mode msa --out "$d/snp.lq.tsv" > /dev/null || { echo "export failed"; exit 1; }
  MK=(--markers "SOF:$d/markers/SOF.fa" --markers "SRO:$d/markers/SRO.fa" --markers "SSP:$d/markers/SSP.fa")
  PMK=("SOF:$d/markers/SOF.fa" "SRO:$d/markers/SRO.fa" "SSP:$d/markers/SSP.fa")
  BASE=(-a "$d/ancestry.tsv" -t "$d/theta.tsv" --fasta "$d/ref.fa")
  echo "== $mode markers"
  # ---- A
  PYTHONPATH=$FP python3 "$PKG/lineage_kmer/snipe_kmer.py" -i "$d/ids_all.txt" -b "$d/case.sam" --kmers "${PMK[@]}" -k 31 -o "$d/py_kmer_all.tsv" -t 2 > "$d/py_kmer_all.log" 2>&1 || { echo "  snipe_kmer failed"; tail -5 "$d/py_kmer_all.log"; exit 1; }
  PYTHONPATH=$FP python3 "$PKG/lineage_snp/snipe_snp.py" -i "$d/ids_all.txt" -b "$d/case.sam" --snp_db "$d/snp_db.pkl" --mode msa -o "$d/py_snp_all.tsv" -t 2 > "$d/py_snp_all.log" 2>&1 || { echo "  snipe_snp failed"; tail -5 "$d/py_snp_all.log"; exit 1; }
  "$LQ" run --sam-in "$d/case.sam" "${BASE[@]}" --target-ids "$d/ids_all.txt" "${MK[@]}" --snp-db "$d/snp.lq.tsv" \
      --debug-kmer "$d/cpp_kmer_all.tsv" --debug-snp "$d/cpp_snp_all.tsv" -o "$d/cpp_em_all.tsv" --threads 2 > "$d/cpp_all.log" 2>&1 || { echo "  lqcore failed"; tail -5 "$d/cpp_all.log"; exit 1; }
  echo "  A (all transcripts are targets)"
  cmp_rows kmer "$d/py_kmer_all.tsv" "$d/cpp_kmer_all.tsv" || fail=1
  cmp_rows snp "$d/py_snp_all.tsv" "$d/cpp_snp_all.tsv" || fail=1
  # ---- B
  PYTHONPATH=$FP python3 "$PKG/lineage_kmer/snipe_kmer.py" -i "$d/ids_one.txt" -b "$d/case.sam" --kmers "${PMK[@]}" -k 31 -o "$d/py_kmer_one.tsv" -t 2 > /dev/null 2>&1
  PYTHONPATH=$FP python3 "$PKG/lineage_snp/snipe_snp.py" -i "$d/ids_one.txt" -b "$d/case.sam" --snp_db "$d/snp_db.pkl" --mode msa -o "$d/py_snp_one.tsv" -t 2 > /dev/null 2>&1
  python3 "$PKG/merge/merge_poe.py" --kmer "$d/py_kmer_one.tsv" --snp "$d/py_snp_one.tsv" --out "$d/py_merged_one.tsv" > /dev/null
  python3 "$PKG/core/em_global.py" "${BASE[@]}" -p "$d/py_merged_one.tsv" --bam "$d/case.sam" -o "$d/py_em_one.tsv" > "$d/py_em_one.log" 2>&1
  "$LQ" run --sam-in "$d/case.sam" "${BASE[@]}" --target-ids "$d/ids_one.txt" "${MK[@]}" --snp-db "$d/snp.lq.tsv" \
      --write-pvals "$d/cpp_pv_one.tsv" -o "$d/cpp_em_one.tsv" --threads 4 > "$d/cpp_one.log" 2>&1
  echo "  B (one target transcript per cluster: whole chain)"
  cmp_rows merged "$d/py_merged_one.tsv" "$d/cpp_pv_one.tsv" || fail=1
  if cmp -s "$d/py_em_one.tsv" "$d/cpp_em_one.tsv"; then echo "    EM table: identical to python"; else echo "    EM table: DIFFERENT"; maxdiff "$d/py_em_one.tsv" "$d/cpp_em_one.tsv"; fail=1; fi
  # ---- C
  "$LQ" run --sam-in "$d/case.sam" "${BASE[@]}" --target-ids "$d/ids_one.txt" "${MK[@]}" --snp-db "$d/snp.lq.tsv" -o "$d/cpp_em_one_t1.tsv" --threads 1 > /dev/null 2>&1
  for i in 1 2 3; do
    "$LQ" run --sam-in "$d/case.sam" "${BASE[@]}" --target-ids "$d/ids_one.txt" "${MK[@]}" --snp-db "$d/snp.lq.tsv" \
        --chunk $i/3 --block-pairs 500 --part-out "$d/part$i" --threads 2 > "$d/part$i.log" 2>&1 || { echo "  part $i failed"; tail -5 "$d/part$i.log"; fail=1; }
  done
  "$LQ" em "${BASE[@]}" --parts "$d/part1" "$d/part2" "$d/part3" -o "$d/cpp_em_parts.tsv" --idxstats-out "$d/parts.idxstats" --threads 2 > "$d/parts_em.log" 2>&1 || { echo "  em --parts failed"; tail -5 "$d/parts_em.log"; fail=1; }
  echo "  C (threads, chunks)"
  cmp -s "$d/cpp_em_one.tsv" "$d/cpp_em_one_t1.tsv" && echo "    threads 4 = threads 1" || { echo "    THREADS DIFFER"; fail=1; }
  echo "    3 chunks merged vs single run:"; maxdiff "$d/cpp_em_one.tsv" "$d/cpp_em_parts.tsv" || fail=1
  grep -h "records ->" "$d/cpp_one.log" "$d/parts_em.log" | cut -c22- | sed 's/^/    /'
done
echo; [ $fail = 0 ] && echo "ALL LINEAGE CHECKS PASSED" || echo "SOME LINEAGE CHECKS FAILED"
exit $fail
