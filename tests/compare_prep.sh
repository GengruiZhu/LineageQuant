#!/bin/bash
# Reference preparation: C++ ports vs the Python modules of LineageQuant 1.0.0, byte for byte, on synthetic data
# (tests/make_prep_case.py). External tools are replaced by test stand-ins (tests/fake_tools: mafft, kmc, kmc_tools,
# gffread, samtools faidx; tests/fakecyvcf2, tests/fakebio) -- both sides call the same ones.
#   bash tests/compare_prep.sh LINEAGEQUANT P100_PKG_DIR [WORKDIR]      (P100_PKG_DIR = .../LineageQuant_v1.0.0/lineagequant)
set -u
LQ=$(readlink -f "${1:?lineagequant binary}"); P=$(readlink -f "${2:?LineageQuant 1.0.0 package dir}"); W=${3:-/tmp/lq_prep_test}
HERE=$(cd "$(dirname "$0")" && pwd)
export PYTHONPATH=$HERE/fakecyvcf2:$HERE/fakebio:$(dirname "$P")${PYTHONPATH:+:$PYTHONPATH}
export PYTHONDONTWRITEBYTECODE=1
rm -rf "$W"; mkdir -p "$W/fakebin"; fail=0
# stand-ins copied and made executable here: a checkout may have lost the exec bits (git on Windows, zip files)
cp "$HERE"/fake_tools/* "$W/fakebin/" && chmod +x "$W/fakebin"/* || exit 1
export PATH=$W/fakebin:$PATH
ok () { printf "  %-66s %s\n" "$1" "$2"; [ "${2%% *}" = identical ] || fail=1; }
same () {
  if [ ! -e "$1" ] && [ ! -e "$2" ]; then echo "identical (absent on both sides)"
  elif cmp -s "$1" "$2"; then echo "identical ($(wc -l < "$1") lines)"; else echo "DIFF: diff $1 $2"; fi; }

for case in dots plain; do
  C=$W/$case; python3 "$HERE/make_prep_case.py" "$C" $([ $case = dots ] && echo 1 || echo 2) $case > /dev/null || exit 1
  echo "== case $case"

  # ---------------------------------------------------------------- normalize-clusters
  N=$C/norm; mkdir -p "$N"
  nc () {   # tag args...
    local tag=$1; shift
    python3 "$P/cluster/normalize_clusters.py" "$@" --clusters-out "$N/py_$tag.c" --ancestry-out "$N/py_$tag.a" --orphan-out "$N/py_$tag.o" \
       --all-genes "$C/all_genes.txt" --subgenome-tag 'A:^Gh_A' --subgenome-tag 'D:_D\d+$' > /dev/null 2>&1 || { echo "  python normalize failed ($tag)"; fail=1; return; }
    "$LQ" normalize-clusters "$@" --clusters-out "$N/cpp_$tag.c" --ancestry-out "$N/cpp_$tag.a" --orphan-out "$N/cpp_$tag.o" \
       --all-genes "$C/all_genes.txt" --subgenome-tag 'A:^Gh_A' --subgenome-tag 'D:_D\d+$' > /dev/null 2>&1 || { echo "  lineagequant normalize failed ($tag)"; fail=1; return; }
    for x in c a o; do ok "normalize-clusters $tag .$x" "$(same "$N/py_$tag.$x" "$N/cpp_$tag.$x")"; done
  }
  nc standard --input "$C/std.tsv" --format standard
  nc pairs --input "$C/pairs.txt" --format pairs
  nc jcvi --input "$C/pairs.txt" --format jcvi --min-size 3
  nc tsv_all --input "$C/table.tsv" --format tsv
  nc tsv_cols --input "$C/table.tsv" --format tsv --tsv-cols A 2 E --tsv-id-col D

  # ---------------------------------------------------------------- extract-transcripts (genome_pair)
  X=$C/extract; mkdir -p "$X/py" "$X/cpp"
  for t in A B; do
    python3 - "$C/anc$t.fa" "$C/anc$t.gff" "$X/py" "${t}_s0" "$X/py/_${t}_s0_transcript.fa" "$t" <<'PY' > /dev/null 2>&1 || { echo "  python genome_pair failed ($t)"; fail=1; }
import sys
from lineagequant.utils.genome_pair import GenomeAnnotationPair
fa, gff, wd, name, out, pre = sys.argv[1:]
GenomeAnnotationPair(fa, gff, wd, name=name).prepare().extract_transcripts(out, id_prefix=pre)
PY
    "$LQ" extract-transcripts --fasta "$C/anc$t.fa" --gff "$C/anc$t.gff" --workdir "$X/cpp" --name "${t}_s0" \
         --out "$X/cpp/_${t}_s0_transcript.fa" --id-prefix "$t" > /dev/null 2>&1 || { echo "  lq extract-transcripts failed ($t)"; fail=1; }
    for f in "_${t}_s0.genome.fa" "_${t}_s0.genome.fa.fai" "_${t}_s0.clean.gff" "_${t}_s0_transcript.fa"; do
      ok "extract-transcripts $f" "$(same "$X/py/$f" "$X/cpp/$f")"
    done
  done

  # ---------------------------------------------------------------- build-markers (KMC)
  M=$C/markers; mkdir -p "$M"
  A=$X/py/_A_s0_transcript.fa; Bf=$X/py/_B_s0_transcript.fa
  python3 - "$A" "$Bf" "$M/py_anc" <<'PY' > /dev/null 2>&1 || { echo "  python build_markers failed"; fail=1; }
import sys
from lineagequant.lineage_kmer.build_markers import build_markers
a, b, out = sys.argv[1:]
print(build_markers({"A": a, "B": b}, out, k=31, threads=2))        # the router's call (mode ancestor)
PY
  "$LQ" build-markers --lineage A:"$A" --lineage B:"$Bf" -o "$M/cpp_anc" -k 31 -t 2 > /dev/null 2>&1 || { echo "  lq build-markers failed"; fail=1; }
  for l in A B; do ok "build-markers (ancestor) $l.fasta" "$(same "$M/py_anc/$l.fasta" "$M/cpp_anc/$l.fasta")"; done
  # ancestor_x_poly with two inputs for one lineage, k 21, kept intermediates; poly
  awk '/^>/{n++} n%3' "$A" > "$M/polyA.fa"; awk '/^>/{n++} n%4' "$Bf" > "$M/polyB.fa"    # subsets of the progenitor transcripts
  python3 "$P/lineage_kmer/build_markers.py" --lineage A:"$A" --lineage A:"$M/polyB.fa" --lineage B:"$Bf" --polyploid A:"$M/polyA.fa" \
      --polyploid B:"$M/polyB.fa" --mode ancestor_x_poly -o "$M/py_axp" -k 21 -t 2 --keep-intermediate > /dev/null 2>&1 || { echo "  python axp failed"; fail=1; }
  "$LQ" build-markers --lineage A:"$A" --lineage A:"$M/polyB.fa" --lineage B:"$Bf" --polyploid A:"$M/polyA.fa" \
      --polyploid B:"$M/polyB.fa" --mode ancestor_x_poly -o "$M/cpp_axp" -k 21 -t 2 --keep-intermediate > /dev/null 2>&1 || { echo "  lineagequant axp failed"; fail=1; }
  for f in A.fasta B.fasta merged.txt merged_lc.txt _anc_A_merged.fa; do ok "build-markers (ancestor_x_poly, k 21) $f" "$(same "$M/py_axp/$f" "$M/cpp_axp/$f")"; done
  sed "s|$M/py_axp|D|g" "$M/py_axp/A.cfg" > "$M/py_A.cfg"; sed "s|$M/cpp_axp|D|g" "$M/cpp_axp/A.cfg" > "$M/cpp_A.cfg"
  ok "build-markers (ancestor_x_poly) A.cfg (directory names replaced)" "$(same "$M/py_A.cfg" "$M/cpp_A.cfg")"
  python3 "$P/lineage_kmer/build_markers.py" --polyploid A:"$M/polyA.fa" --polyploid B:"$M/polyB.fa" --mode poly -o "$M/py_poly" -k 25 > /dev/null 2>&1 || { echo "  python poly failed"; fail=1; }
  "$LQ" build-markers --polyploid A:"$M/polyA.fa" --polyploid B:"$M/polyB.fa" --mode poly -o "$M/cpp_poly" -k 25 > /dev/null 2>&1 || { echo "  lineagequant poly failed"; fail=1; }
  for l in A B; do ok "build-markers (poly, k 25) $l.fasta" "$(same "$M/py_poly/$l.fasta" "$M/cpp_poly/$l.fasta")"; done

  # ---------------------------------------------------------------- snpdb-msa (MAFFT)
  S=$C/snp; mkdir -p "$S"
  python3 "$P/lineage_snp/build_snp_msa.py" --fasta "$C/copies.fa" --cluster "$C/clusters.tsv" --ancestry "$C/ancestry.tsv" \
      --out "$S/msa.pkl" --threads 1 --lq-out "$S/py_msa.lq.tsv" > "$S/py_msa.log" 2>&1 || { echo "  python build_snp_msa failed"; tail -3 "$S/py_msa.log"; fail=1; }
  "$LQ" snpdb-msa --fasta "$C/copies.fa" --cluster "$C/clusters.tsv" --ancestry "$C/ancestry.tsv" --out "$S/cpp_msa_legacy.lq.tsv" \
      --threads 4 --legacy-lineages > "$S/cpp_msa_legacy.log" 2>&1 || { echo "  lq snpdb-msa failed"; tail -3 "$S/cpp_msa_legacy.log"; fail=1; }
  ok "snpdb-msa --legacy-lineages (4 threads) vs build_snp_msa.py" "$(same "$S/py_msa.lq.tsv" "$S/cpp_msa_legacy.lq.tsv")"
  "$LQ" snpdb-msa --fasta "$C/copies.fa" --cluster "$C/clusters.tsv" --ancestry "$C/ancestry.tsv" --out "$S/cpp_msa.lq.tsv" \
      --threads 3 > "$S/cpp_msa.log" 2>&1 || { echo "  lq snpdb-msa failed"; fail=1; }
  echo "    lineages 1.0.0: $(head -1 "$S/py_msa.lq.tsv" | cut -f3)   1.1.0: $(head -1 "$S/cpp_msa.lq.tsv" | cut -f3)"
  if [ $case = plain ]; then
    # without multi-dot ids the copies' lineages agree; only the lineage list differs -> the 1.0.0 DB written with all lineages
    python3 - "$S/msa.pkl" "$C/ancestry.tsv" "$S/py_msa_all.lq.tsv" <<'PY' || fail=1
import pickle, sys
from lineagequant.lineage_snp import export_snp_db as E
from lineagequant.utils.lineage import lineage_names_from_ancestry
db = pickle.load(open(sys.argv[1], "rb"))
lins = lineage_names_from_ancestry(sys.argv[2])
E.infer_lineages = lambda d: lins
E.write_lq(db, "msa", sys.argv[3])
PY
    ok "snpdb-msa (all lineages) vs 1.0.0 DB written with all lineages" "$(same "$S/py_msa_all.lq.tsv" "$S/cpp_msa.lq.tsv")"
  fi

  # ---------------------------------------------------------------- snpdb-vcf
  python3 "$P/lineage_snp/build_snp_vcf.py" --vcf "$C/pop.vcf.gz" --gff3 "$C/genome.gff3" --sample-lists L2:"$C/L2.samples" \
      L1:"$C/L1.samples" L3:"$C/L3.samples" --out "$S/pop.pkl" > "$S/py_vcf.log" 2>&1 \
    && python3 "$P/lineage_snp/export_snp_db.py" --snp_db "$S/pop.pkl" --mode population --out "$S/py_pop.lq.tsv" >> "$S/py_vcf.log" 2>&1 \
    || { echo "  python VCF route failed"; tail -3 "$S/py_vcf.log"; fail=1; }
  "$LQ" snpdb-vcf --vcf "$C/pop.vcf.gz" --gff3 "$C/genome.gff3" --sample-lists L2:"$C/L2.samples" L1:"$C/L1.samples" L3:"$C/L3.samples" \
      --out "$S/cpp_pop.lq.tsv" > "$S/cpp_vcf.log" 2>&1 || { echo "  lq snpdb-vcf failed"; tail -3 "$S/cpp_vcf.log"; fail=1; }
  ok "snpdb-vcf vs build_snp_vcf.py + export_snp_db.py" "$(same "$S/py_pop.lq.tsv" "$S/cpp_pop.lq.tsv")"
  python3 "$P/lineage_snp/build_snp_vcf.py" --vcf "$C/pop.vcf.gz" --gff3 "$C/genome.gff3" --sample-lists L1:"$C/L1.samples" L3:"$C/L3.samples" \
      --out "$S/pop2.pkl" --min_dp 5 --max_dp 40 --min_called 2 --min_af_diff 0.2 > /dev/null 2>&1 \
    && python3 "$P/lineage_snp/export_snp_db.py" --snp_db "$S/pop2.pkl" --mode population --out "$S/py_pop2.lq.tsv" > /dev/null 2>&1 || fail=1
  "$LQ" snpdb-vcf --vcf "$C/pop.vcf.gz" --gff3 "$C/genome.gff3" --sample-lists L1:"$C/L1.samples" L3:"$C/L3.samples" \
      --out "$S/cpp_pop2.lq.tsv" --min_dp 5 --max_dp 40 --min_called 2 --min_af_diff 0.2 > /dev/null 2>&1 || fail=1
  ok "snpdb-vcf (min_dp 5, max_dp 40, min_called 2, min_af_diff 0.2)" "$(same "$S/py_pop2.lq.tsv" "$S/cpp_pop2.lq.tsv")"
done
echo; [ $fail = 0 ] && echo "ALL PREPARATION STEPS IDENTICAL" || echo "SOME PREPARATION STEPS DIFFER"
exit $fail
