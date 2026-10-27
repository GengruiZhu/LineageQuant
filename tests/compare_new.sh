#!/bin/bash
# lqcore 1.1.0 additions (all opt-in; the defaults stay byte-identical with LineageQuant 1.0.0 -> compare_em.sh etc.)
#   bash tests/compare_new.sh LQCORE [WORKDIR]
#   1  --eff-len fld   effective length = length - mean of the observed template lengths truncated at the length + 1;
#                      checked against a direct computation from the SAM (every target, 4 decimals)
#   2  --gc-bias       synthetic fragments with no / weak / strong fragment-GC bias (tests/make_gc_case.py):
#                      TPM = reads / effective length vs the true TPM (MARD, Spearman); GC correction must help with
#                      bias and must not hurt without it; part files (run --gc-hist -> em --parts --gc-bias) = one pass
#   3  --marker-clean  marker k-mers found in a reference copy of another lineage are dropped, the rest kept
set -u
LQ=$(readlink -f "${1:?lqcore}"); W=${2:-/tmp/lq_new_test}; H=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$W"; fail=0

echo "== 1  --eff-len fld"
d=$W/gc_strong; [ -s "$d/case.sam" ] || python3 "$H/make_gc_case.py" "$d" 1 60000 strong
"$LQ" em -a "$d/ancestry.tsv" --sam "$d/case.sam" --fasta "$d/ref.fa" --no-pool --eff-len fld --eff-out "$d/eff_fld.tsv" \
   -o "$d/em_fld.tsv" --threads 2 > "$d/fld.log" 2>&1 || { echo "  lqcore failed"; tail -3 "$d/fld.log"; exit 1; }
python3 - "$d/case.sam" "$d/eff_fld.tsv" <<'PY' || fail=1
import sys
fl = []
for l in open(sys.argv[1]):
    if l.startswith("@") or not l.strip(): continue
    x = l.split("\t"); f = int(x[1]); t = abs(int(x[8]))
    if (f & 2) and not (f & 0x900) and 0 < t < 2000 and len(fl) < 1000000: fl.append(t)
h = [0] * 2000
for t in fl: h[t] += 1
s0 = [0.0] * 2000; s1 = [0.0] * 2000
for x in range(1, 2000): s0[x] = s0[x - 1] + h[x]; s1[x] = s1[x - 1] + h[x] * x
bad = n = 0
for l in list(open(sys.argv[2]))[1:]:
    nm, L, e, g = l.split("\t"); L = float(L); m = int(min(L, 1999))
    mt = s1[m] / s0[m] if m >= 1 and s0[m] > 0 else L
    exp = max(L - mt + 1, 25.0)
    n += 1
    if abs(exp - float(e)) > 1e-3: bad += 1
print(f"  {n} targets, {bad} differ from the direct computation")
sys.exit(1 if bad or n == 0 else 0)
PY

echo "== 2  --gc-bias (TPM = reads / effective length vs truth)"
score () {   # dir em eff -> "MARD Spearman"
  python3 - "$1/truth.tsv" "$2" "$3" <<'PY'
import sys, math
tr = {l.split()[0]: float(l.split()[1]) for l in list(open(sys.argv[1]))[1:]}
rd = {}
for l in list(open(sys.argv[2]))[1:]:
    x = l.split("\t"); rd[x[1]] = float(x[4])
ef = {}
for l in list(open(sys.argv[3]))[1:]:
    x = l.split("\t"); ef[x[0]] = float(x[2])
v = {t: rd.get(t, 0.0) / ef[t] for t in tr if t in ef}
z = sum(v.values()); est = {t: x / z * 1e6 for t, x in v.items()}
ts = sorted(tr)
a = [tr[t] for t in ts]; b = [est.get(t, 0.0) for t in ts]
mard = sum(0 if x == y == 0 else abs(x - y) / (x + y) for x, y in zip(a, b)) / len(a)
def rank(u):
    o = sorted(range(len(u)), key=lambda i: u[i]); r = [0.0] * len(u); i = 0
    while i < len(o):
        j = i
        while j + 1 < len(o) and u[o[j + 1]] == u[o[i]]: j += 1
        for k in range(i, j + 1): r[o[k]] = (i + j) / 2
        i = j + 1
    return r
ra, rb = rank(a), rank(b); ma, mb = sum(ra) / len(ra), sum(rb) / len(rb)
c = sum((x - ma) * (y - mb) for x, y in zip(ra, rb)) / math.sqrt(sum((x - ma) ** 2 for x in ra) * sum((y - mb) ** 2 for y in rb))
print(f"{mard:.4f} {c:.4f}")
PY
}
for b in none weak strong; do
  d=$W/gc_$b; [ -s "$d/case.sam" ] || python3 "$H/make_gc_case.py" "$d" 1 60000 $b > /dev/null
  base=(-a "$d/ancestry.tsv" --sam "$d/case.sam" --fasta "$d/ref.fa" --no-pool --threads 2)
  "$LQ" em "${base[@]}" --eff-out "$d/e_eff.tsv" -o "$d/m_eff.tsv" > "$d/l_eff.log" 2>&1 || { echo "  failed eff $b"; fail=1; continue; }
  "$LQ" em "${base[@]}" --eff-len fld --eff-out "$d/e_fld.tsv" -o "$d/m_fld.tsv" > "$d/l_fld.log" 2>&1 || { echo "  failed fld $b"; fail=1; continue; }
  "$LQ" em "${base[@]}" --eff-len fld --gc-bias --eff-out "$d/e_gc.tsv" -o "$d/m_gc.tsv" > "$d/l_gc.log" 2>&1 || { echo "  failed gc $b"; tail -3 "$d/l_gc.log"; fail=1; continue; }
  s1=$(score "$d" "$d/m_eff.tsv" "$d/e_eff.tsv"); s2=$(score "$d" "$d/m_fld.tsv" "$d/e_fld.tsv"); s3=$(score "$d" "$d/m_gc.tsv" "$d/e_gc.tsv")
  printf "  bias %-7s  MARD/Spearman  eff %s   fld %s   fld+gc %s\n" "$b" "$s1" "$s2" "$s3"
  grep -h "\[GC\] effective-length factor" "$d/l_gc.log" | sed 's/^.*\[GC\]/            [GC]/'
  m2=${s2% *}; m3=${s3% *}
  case $b in
    strong) python3 -c "import sys; sys.exit(0 if $m3 < 0.8 * $m2 else 1)" || { echo "  ✗ GC correction should clearly help under strong bias"; fail=1; } ;;
    none)   python3 -c "import sys; sys.exit(0 if $m3 < 1.05 * $m2 + 0.002 else 1)" || { echo "  ✗ GC correction hurts without bias"; fail=1; } ;;
  esac
done
# the same through run --gc-hist -> part -> em --parts --gc-bias
d=$W/gc_strong
"$LQ" run -a "$d/ancestry.tsv" --sam-in "$d/case.sam" --fasta "$d/ref.fa" --no-pool --threads 2 --gc-hist --part-out "$d/c.part" > "$d/l_part.log" 2>&1 \
  && "$LQ" em -a "$d/ancestry.tsv" --parts "$d/c.part" --fasta "$d/ref.fa" --no-pool --threads 2 --eff-len fld --gc-bias \
       --eff-out "$d/e_gcp.tsv" -o "$d/m_gcp.tsv" > "$d/l_gcp.log" 2>&1 || { echo "  ✗ run --gc-hist / em --parts failed"; tail -3 "$d/l_gcp.log"; fail=1; }
if cmp -s <(cut -f1-6 "$d/m_gc.tsv") <(cut -f1-6 "$d/m_gcp.tsv") && cmp -s "$d/e_gc.tsv" "$d/e_gcp.tsv"; then echo "  run --gc-hist -> em --parts --gc-bias = em --sam --gc-bias"
else echo "  ✗ part route differs from the single pass"; fail=1; fi
for th in 1 4; do
  "$LQ" em -a "$d/ancestry.tsv" --sam "$d/case.sam" --fasta "$d/ref.fa" --no-pool --threads $th --eff-len fld --gc-bias \
     --eff-out "$d/e_t$th.tsv" -o "$d/m_t$th.tsv" > /dev/null 2>&1
done
cmp -s "$d/m_t1.tsv" "$d/m_t4.tsv" && cmp -s "$d/e_t1.tsv" "$d/e_t4.tsv" && echo "  threads 1 = 4" || { echo "  ✗ threads change the GC result"; fail=1; }

echo "== 3  --marker-clean"
d=$W/mk; mkdir -p "$d"
python3 - "$d" <<'PY'
import random, sys
R = random.Random(5); d = sys.argv[1]
rs = lambda n: "".join(R.choice("ACGT") for _ in range(n))
A, B, U = rs(400), rs(400), rs(400)          # A-labelled copy, B-labelled copy, sequence absent from the reference
open(f"{d}/ref.fa", "w").write(f">a1\n{A}\n>b1\n{B}\n>o1\n{rs(300)}\n")
open(f"{d}/anc.tsv", "w").write("GeneID\tAncestry\na1\tA\nb1\tB\n")
K = 25
rc = lambda s: s.translate(str.maketrans("ACGT", "TGCA"))[::-1]
kms = lambda s: [min(s[i:i + K], rc(s[i:i + K])) for i in range(len(s) - K + 1)]   # canonical, as KMC writes them
mA = kms(A)[:100] + kms(B)[:30] + kms(U)[:40]     # 30 of A's markers are in the B copy -> must go
mB = kms(B)[100:200] + kms(A)[200:210]            # 10 of B's markers are in the A copy -> must go
with open(f"{d}/A.fa", "w") as f:
    for i, k in enumerate(mA): f.write(f">{i}\n{k}\n")
with open(f"{d}/B.fa", "w") as f:
    for i, k in enumerate(mB): f.write(f">{i}\n{k}\n")
PY
"$LQ" markers --markers "A:$d/A.fa" --markers "B:$d/B.fa" --kmer-len 25 --marker-cache "$d/raw.bin" > "$d/raw.log" 2>&1
"$LQ" markers --markers "A:$d/A.fa" --markers "B:$d/B.fa" --kmer-len 25 --marker-cache "$d/clean.bin" --marker-clean --fasta "$d/ref.fa" -a "$d/anc.tsv" > "$d/clean.log" 2>&1
r=$(grep -oE "A [0-9,]+ -> [0-9,]+ \(-[0-9.]+%\)  B [0-9,]+ -> [0-9,]+" "$d/clean.log")
echo "  $r"
[ "$r" = "A 170 -> 140 (-17.65%)  B 110 -> 100" ] || { echo "  ✗ expected A 170 -> 140, B 110 -> 100"; cat "$d/clean.log"; fail=1; }
cmp -s "$d/raw.bin" "$d/clean.bin" && { echo "  ✗ cleaned cache identical to the raw one"; fail=1; }

echo; [ $fail = 0 ] && echo "ALL 1.1.0 CHECKS PASSED" || echo "SOME 1.1.0 CHECKS FAILED"
exit $fail
