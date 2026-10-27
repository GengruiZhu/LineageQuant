# LineageQuant 1.1.0 file formats

## Binary outputs

Little-endian. `u8/u32/u64` unsigned integers, `i64` signed, `f64` IEEE-754 double, `str` = `u32` length +
bytes (no terminator). Every file starts with an 8-byte magic. Writers go through a temporary file and a
rename, so a reader never sees a half-written file. The source of truth is `src/binio.hpp`.

The numbers are the ones the text writers of 1.0.0 print, before printing: `lq export` turns
each file into the 1.0.0 text files byte for byte, and every reader that takes the text also takes the binary
(`quant -t`, `em -t`, `matrix`). Readers apply the same rounding the text would impose (Theta `%.6f`,
EM-assigned fragments `%.2f`, effective length `%.4f`), so results do not depend on which form was read.

### `Global_Theta_Matrix.lqb` — `theta --binary`

```
"LQTHETA1"
u32 n_clusters
  str cluster_id                 x n_clusters
u64 n_rows                       (one per allele copy, in the order of the text file)
  u32 cluster                    index into the cluster list
  str allele_id
  u64 unique_probes              Unique_Probes
  u64 raw_hits                   Raw_Hits
  u8  flags                      1 = has a long-read count, 2 = LongRead_Anchored (else Kmer_Inferred)
  f64 lr_count                   LR_Count (text: Python repr, "0" without a count)
  f64 theta                      Theta (text: %.6f)
  f64 shrinkage                  Shrinkage_Used (text: %.2f; "NA" when anchored)
```

### `<sample>.lqres` — `quant -o <sample>.lqres` (also `em -o`)

```
"LQRES001"
i64 em_iterations                EM_Iter
u64 n_targets                    (EM order)
  str cluster_id                 ClusterID ("ORPHAN" for single-copy genes)
  str allele_id                  AlleleID
  str ancestry                   Ancestry
  f64 initial_theta              Initial_Theta (%.6f)
  f64 em_reads                   EM_Allocated_Reads (%.2f)
  f64 final_theta                Final_Theta (%.6f)
  f64 length                     efflen Length (%.0f)
  f64 effective_length           efflen EffectiveLength (%.4f)
  f64 gc_factor                  efflen GC_Factor (%.5f)
u64 n_references                 (samtools idxstats semantics)
  str name
  i64 length
  u64 mapped
  u64 unmapped
u64 unmapped_without_reference   the "*" line
```

`lq export S.lqres -o DIR` writes `DIR/EM_Final_Expression.tsv` (`--em-name` to rename, e.g.
`EM_Final_Hybrid.tsv`), `DIR/S.idxstats.tsv`, `DIR/S.efflen.tsv`;
`lq export Global_Theta_Matrix.lqb` writes `Global_Theta_Matrix.tsv` next to it.

## Text files

| File | Columns |
|---|---|
| `EM_Final_Expression.tsv` | `ClusterID AlleleID Ancestry Initial_Theta EM_Allocated_Reads Final_Theta EM_Iter EM_LogLik` (header line) |
| `<sample>.idxstats.tsv` | `name length mapped unmapped` (no header; last line `*`) |
| `<sample>.efflen.tsv` | `Target Length EffectiveLength GC_Factor` (header line) |
| `Global_Theta_Matrix.tsv` | `Cluster_ID Allele_ID Unique_Probes Raw_Hits LR_Count Theta Prior_Status Shrinkage_Used` |
| `LineageQuant_Final_Merged_TPM.tsv` | `GeneID Type Length_bp <s>_RawReads... <s>_TPM...` (floats as Python repr, rows sorted by Type descending, GeneID) |

### SNP database (`snpdb-msa`, `snpdb-vcf`; read by `quant --snp-db`)

```
#LQSNP  mode=msa         lineages=L1,L2,...  n=N
tid  pos  EXPECTED                 msa: one character per lineage, A/C/G/T or '.' (never matches)
tid  pos  REF  ALT  AF,AF,...      population: allele frequency per lineage (Python repr)
```

Grouped by transcript in the order the transcripts first appear in the database, positions (0-based, on the
transcript) ascending.

### Marker FASTA (`build-markers`)

`OUTDIR/<lineage>.fasta`, two lines per k-mer: `><lineage>-1` and the k-mer (canonical, as KMC dumps it).

### Indexes and caches in `--data-dir` (internal)

`transcript_index_bt2.*` (bowtie2), `salmon_kd_index/`, `lqmap.k21.bin`, `lqcore_theta.k<k>.bin`,
`lqcore_markers.k<k>[.clean].bin`, `lineage_markers/`, `ancestor_transcripts/`, `SNP_DB_*.lq.tsv`,
`LongRead_Counts.tsv`. The `.bin` caches record the size and modification time of their inputs and are
rebuilt when those change; the other entries are rebuilt only after their `.done` file is deleted.
