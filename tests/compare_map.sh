#!/bin/bash
# lqmap vs an exhaustive ungapped aligner (tests/brute_map.py) on synthetic polyploid-like data.
#   bash tests/compare_map.sh bin/lqcore
# Per fragment, the set of candidate transcripts and their pair NM must be identical; threads 1 = 4.
set -o pipefail
LQ=$(readlink -f "$1"); T=$(cd "$(dirname "$0")" && pwd)
[ -x "$LQ" ] || { echo "usage: $0 LQCORE"; exit 2; }
W=$(mktemp -d); trap 'rm -rf "$W"' EXIT
fail=0
sam2sets () {   # SAM -> name \t transcript:pairNM ... (sorted), one line per fragment
  awk -F'\t' '!/^@/ {
      n=$1; if (n!=cur) { if (cur!="") emit(); cur=n; delete nm; delete seen; k=0 }
      if (and($2,4)) next
      m = and($2,64) ? 1 : 2; nmv=0
      for (i=12;i<=NF;i++) if ($i ~ /^NM:i:/) nmv=substr($i,6)+0
      key=$3; if (!(key in seen)) { seen[key]=1; order[++k]=key; nm[key]=0 }
      nm[key]+=nmv
    }
    function emit(   i,s,a,n2) { s=""; n2=0; for (i=1;i<=k;i++) a[++n2]=order[i] ":" nm[order[i]]; asort(a); for (i=1;i<=n2;i++) s = s (i>1?" ":"") a[i]; print cur "\t" s; delete a }
    END { if (cur!="") emit() }' "$1"
}
for cfg in "1 400 100" "2 300 150" "3 400 100"; do
  set -- $cfg; seed=$1 np=$2 rl=$3
  d=$W/c$seed; python3 "$T/make_map_case.py" "$d" "$seed" "$np" "$rl" > /dev/null || exit 1
  python3 "$T/brute_map.py" "$d/tx.fa" "$d/R1.fq.gz" "$d/R2.fq.gz" > "$d/brute.txt" || { echo "brute failed"; exit 1; }
  for th in 1 4; do
    "$LQ" map --fasta "$d/tx.fa" --r1 "$d/R1.fq.gz" --r2 "$d/R2.fq.gz" --threads $th > "$d/lq$th.sam" 2> "$d/lq$th.log" || { echo "lqcore map failed"; cat "$d/lq$th.log"; exit 1; }
  done
  sam2sets "$d/lq1.sam" > "$d/lq.txt"
  cmp -s "$d/lq1.sam" "$d/lq4.sam" && thr="threads 1=4" || { thr="THREADS DIFFER"; fail=1; }
  tot=$(wc -l < "$d/brute.txt")
  same=$(paste "$d/brute.txt" "$d/lq.txt" | awk -F'\t' '$2==$4' | wc -l)
  aln=$(awk -F'\t' '$2!=""' "$d/brute.txt" | wc -l)
  printf "  seed %s  %3s bp  %4s pairs (%s aligned)  identical candidate sets %s/%s  %s  | %s\n" $seed $rl $tot $aln $same $tot "$thr" \
    "$(grep -oE 'overall [0-9.]+%; [0-9.]+ candidates per aligned pair; re-seeded [0-9,]+, rescued [0-9,]+' "$d/lq1.log")"
  if [ "$same" != "$tot" ]; then
    fail=1
    paste "$d/brute.txt" "$d/lq.txt" | awk -F'\t' '$2!=$4 {print "    brute: " $1 " " $2 "\n    lqmap: " $3 " " $4}' | head -8
  fi
done
echo
[ $fail = 0 ] && echo "ALL LQMAP CHECKS PASSED" || { echo "LQMAP CHECKS FAILED"; exit 1; }
