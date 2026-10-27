# LineageQuant

**Haplotype-resolved expression quantification for polyploid genomes**

[![Version 1.0.0](https://img.shields.io/badge/version-1.0.0-orange)](CHANGELOG.md)
[![License: Non-Commercial Research](https://img.shields.io/badge/License-Non--Commercial-blue.svg)](LICENSE)
[![Python 3.8+](https://img.shields.io/badge/Python-3.8%2B-brightgreen)](https://www.python.org/)

LineageQuant quantifies **subgenome / allele-copy expression** in polyploids. It is
an expression *quantifier*, not a read classifier: every homoeologous or allelic
copy is a separate target, reads that are compatible with several copies are
apportioned by expectation–maximisation (EM), and the output is a read-count and
TPM matrix at copy level that can be summed to subgenome or gene level.

It integrates short-read allopolyploid quantification, long-read/hybrid
quantification and lineage-specific k-mer marker construction into a single
**auto-routing** package: you declare the ploidy type and supply whatever inputs
you have, and the router runs the matching modules.

> **1.0.0** is the first stable release. Short reads are now aligned with bowtie2
> (RSEM settings), and a global EM with mismatch-count weighting assigns them.
> See [CHANGELOG.md](CHANGELOG.md) for everything that changed since 0.1.0 and
> how to reproduce earlier results.

---

## Contents

- [Core concepts](#core-concepts)
- [Pipeline overview](#pipeline-overview)
- [How reads are assigned](#how-reads-are-assigned)
- [Installation](#installation)
- [Quick start](#quick-start)
- [Routing table](#routing-table)
- [Input file formats](#input-file-formats)
- [Output](#output)
- [Key parameters](#key-parameters)
- [Run time and resources](#run-time-and-resources)
- [Reproducing pre-1.0 results](#reproducing-pre-10-results)
- [Benchmark summary](#benchmark-summary)
- [Limitations](#limitations)
- [Troubleshooting](#troubleshooting)

---

## Core concepts

### Autopolyploid vs allopolyploid (the first-level fork)

- **Autopolyploid (`--ploidy auto`)** — subgenomes are not differentiated, so
  there is no "lineage" concept. The pipeline performs **copy-level
  quantification only**: a Theta engine estimates per-copy priors, then the EM
  allocates reads among copies. The lineage-assignment stage is skipped.
- **Allopolyploid (`--ploidy allo`)** — subgenomes carry distinct ancestral
  lineages. The pipeline adds **lineage assignment** (which lineage did a read
  come from?) and passes per-read lineage probabilities into the same EM.

### Two distinct k-mer systems (don't conflate them)

| | Theta k-mer (system 1) | Lineage k-mer (system 2) |
|---|---|---|
| Extraction | cluster-internal unique k-mers | KMC N-way subtraction of ancestor genomes |
| Serves | Step 01 (per-copy Theta prior) | Step 02 (lineage probability) |
| Used by | both ploidies | allopolyploid only |
| Replaceable by | long-read counts (01_hybrid) | MSA SNP or population VCF SNP |

### Three lineage-assignment routes (allopolyploid only)

1. **Ancestor k-mer** (`--ancestor NAME:fa,gff ...`, or pre-built `--markers NAME:fasta`)
   — KMC complex subtraction of ancestral genomes yields lineage-specific marker
   k-mers; reads are scanned against them.
2. **Population VCF SNP** (`--vcf` + `--sample-lists NAME:file ...`) — per-lineage
   allele frequencies from a population VCF; genotype likelihood per read.
3. **MSA SNP** (`--msa-snp`, or `--phased-fasta` for a separate CDS file) — MAFFT
   alignment of the allele copies yields fixed differences; reads are base-matched
   against them.

MSA and VCF are mutually exclusive (both are SNP signals). Either SNP route can
stack with the ancestor k-mer route, triggering Product-of-Experts (PoE) fusion.

### Long reads are an enhancement, not a route

A long-read BAM (`--longread`) anchors the Theta prior with physical counts
(01_hybrid) and calibrates MSA SNP allele frequencies (1b). Without long reads
the k-mer Theta engine and deterministic SNP frequencies are used.

---

## Pipeline overview

```
FASTQ (R1+R2) ──Step 00 bowtie2 (RSEM settings, -k 200)──→ transcriptome BAM
                                                              │
Theta engine (Step 01 / 01_hybrid) ───────────────────────────┤ per-copy prior (EM start values)
                                                              │
[allopolyploid only] lineage assignment                       │
   ancestor k-mer (02)    ┐                                   │
   SNP: MSA / VCF (02)    ┼── PoE (02b) ─────────────────────┤ per-read lineage probabilities
                          ┘                                   │
                                                              ▼
                         global EM over all genes (Step 03): effective length,
                         mismatch-count weighting, identical-sequence pooling
                                                              │
                                    aggregation (Step 04)     ▼
                                    LineageQuant_Final_Merged_TPM.tsv
```

All steps are checkpointed (`.done` markers); re-running resumes from the last
completed step.

---

## How reads are assigned

**Alignment (Step 00).** By default reads are aligned to the transcriptome with
bowtie2 using RSEM's settings — end-to-end, ungapped, `--mp 1,1 --score-min
L,0,-0.1 -I 1 -X 1000 --no-mixed --no-discordant -k 200` — with ambiguous
reference bases not penalised (`--np 0 --n-ceil L,0,0.5`). This reports every
equally good copy of a fragment and also near-optimal ones. HISAT2 (the pre-1.0
aligner, still available with `--aligner hisat2`) returns only a subset of
equally good alignments: in our benchmarks it reported a single, wrong copy for
0.01–0.6% of fragments, and those reads accumulated on lowly expressed copies
next to highly expressed homoeologues or paralogues — an error no EM can repair.

**EM (Step 03, global engine).**

- One EM over read equivalence classes of **all** genes — allele copies and
  single-copy genes alike. Each fragment carries total weight 1, so reads shared
  between clusters or between a cluster and a single-copy gene are never counted
  twice.
- **Mismatch-count weighting** (`--nm-weight auto`): when the candidates of a
  fragment differ in mismatch count, each extra mismatch multiplies the
  likelihood by λ = (p/3)/(1−p), with the error rate p estimated from the best
  alignments. This is what lets the near-optimal alignments reported by bowtie2
  be used correctly.
- **Effective length** (`--eff-len eff`): E-step weight = θ / (L − mean fragment
  length + 1).
- **Identical-sequence pooling**: genes whose transcripts are identical in
  sequence cannot be distinguished by any read; they are estimated as one target
  and split evenly.
- Theta priors from Step 01 give the starting values; long-read anchored clusters
  also contribute a Dirichlet prior; in allopolyploids, per-read lineage
  probabilities multiply the likelihood (`--p-weight-*`).

---

## Installation

```bash
git clone https://github.com/GengruiZhu/LineageQuant.git
cd LineageQuant
conda env create -f environment.yml
conda activate lineagequant
pip install -e .

lineagequant --version        # lineagequant 1.0.0
lineagequant run --help
```

Requirements: Python ≥ 3.8 (numpy, pandas, pysam, biopython), **bowtie2 ≥ 2.4**,
samtools ≥ 1.10. Optional by route: hisat2 (`--aligner hisat2`), kmc +
kmc_tools + gffread (ancestor k-mer route), mafft (MSA SNP route), cyvcf2 +
bcftools (VCF route), minimap2 (long reads). `lineagequant run` checks the tools
it needs for the selected route and names any that are missing.

---

## Quick start

### Step 1 — standardize your allele clusters

LineageQuant does not build clusters. Produce them with your tool of choice
(GraphAllele, AlleleFinder, jcvi synteny, custom) and standardize the output:

```bash
lineagequant normalize-clusters \
  --input  jcvi_pairs.tsv \
  --format pairs \
  --clusters-out clusters.tsv \
  --ancestry-out ancestry.tsv \
  --subgenome-tag A:_A --subgenome-tag D:_D \
  --orphan-out orphans.tsv --all-genes all_gene_ids.txt
```

Supported `--format`: `standard`, `jcvi`, `tsv`, `pairs`. Pair-based inputs are
merged into connected components, so clusters spanning more than two subgenomes
are handled correctly.

### Step 2 — run the pipeline

Allopolyploid with ancestor k-mer markers:

```bash
lineagequant run \
  --fasta      transcriptome.fasta \
  --gff        genome.gff3 \
  --clusters   clusters.tsv \
  --ancestry   ancestry.tsv \
  --orphan     orphans.tsv \
  --target-ids target_gene_ids.txt \
  --ploidy     allo \
  --ancestor   A:ancestorA.fa,ancestorA.gff \
  --ancestor   D:ancestorD.fa,ancestorD.gff \
  --data-dir   data \
  --out-dir    results \
  --threads    32
```

Autopolyploid (copy-level only):

```bash
lineagequant run \
  --fasta transcriptome.fasta --gff genome.gff3 \
  --clusters clusters.tsv --ancestry ancestry.tsv --orphan orphans.tsv \
  --target-ids target_gene_ids.txt \
  --ploidy auto --data-dir data --out-dir results --threads 32
```

`--data-dir` holds the paired FASTQ files (one pair per sample) and receives the
shared indices. Add `--dry-run` to print the resolved route without executing.
More examples are in [`examples/`](examples/); a Chinese usage guide is in
[`docs/USAGE_zh.md`](docs/USAGE_zh.md).

---

## Routing table

The router selects a path from the inputs you provide. Lineage names are read
from the ancestry file (any number N).

### Autopolyploid (`--ploidy auto`)

| # | short | long | route |
|---|---|---|---|
| A1 | ✓ | | Theta k-mer (01) → EM |
| A2 | ✓ | ✓ | long-read Theta (01_hybrid) → EM |

### Allopolyploid (`--ploidy allo`)

| # | short | long | ancestor | VCF | phased | route |
|---|---|---|---|---|---|---|
| B1 | ✓ | | ✓ | | | Theta k-mer + ancestor-k-mer lineage → EM |
| B2 | ✓ | | | ✓ | | Theta k-mer + VCF SNP lineage → EM |
| B3 | ✓ | | ✓ | ✓ | | Theta k-mer + ancestor-k-mer + VCF SNP → **PoE** → EM |
| B4 | ✓ | ✓ | ✓ | | | long-read Theta + ancestor-k-mer lineage → EM |
| B5 | ✓ | ✓ | | | ✓ | long-read Theta + MSA SNP (calibrated) → EM |
| B6 | ✓ | ✓ | ✓ | | ✓ | long-read Theta + ancestor-k-mer + MSA SNP (calibrated) → **PoE** → EM |
| B7 | ✓ | ✓ | | ✓ | | long-read Theta + VCF SNP lineage → EM |
| B8 | ✓ | ✓ | ✓ | ✓ | | long-read Theta + ancestor-k-mer + VCF SNP → **PoE** → EM |

With `--ploidy allo` but no lineage signal, the EM runs without lineage
probabilities (equivalent to the auto route).

---

## Input file formats

| Argument | Format |
|---|---|
| `--fasta` | Transcriptome (mRNA) FASTA — every allele copy a separate record |
| `--gff` | Reference GFF3 |
| `--clusters` | `ClusterID<TAB>Allele1,Allele2,...` (one cluster per line) |
| `--ancestry` | `GeneID<TAB>LineageName` |
| `--orphan` | `OrphanID<TAB>GeneID` (single-copy genes) |
| `--target-ids` | one gene ID per line |
| `--ancestor` | `NAME:fasta,gff` per ancestral lineage (≥2; allopolyploid) |
| `--markers` | `NAME:fasta` pre-built lineage markers (skips the KMC build) |
| `--vcf` | population VCF (allopolyploid SNP route) |
| `--sample-lists` | `NAME:file` per lineage, sample names one per line (with `--vcf`) |
| `--phased-fasta` / `--msa-snp` | MSA SNP route: separate phased CDS FASTA / use the main `--fasta` |
| `--longread` | sorted long-read BAM (minimap2 → transcriptome) |
| FASTQ (in `--data-dir`) | paired, gzipped; `.R1.fq.gz`, `_R1.fastq.gz`, `_1.fq.gz`, `_R1_001.fastq.gz`, … auto-detected |

Gene IDs are derived from transcript IDs by removing the last `.N` suffix, so
`Gene.1` and `Gene.2` are isoforms of `Gene`.

---

## Output

Per sample (`<out-dir>/<sample>/`):

| File | Description |
|---|---|
| `<sample>.sorted.bam` (+ `.bai`) | transcriptome alignments (all reported alignments) |
| `<sample>.bowtie2.log` / `.hisat2.log` | aligner summary |
| `aligner.txt` | aligner used for this sample directory |
| `Global_Theta_Matrix.tsv` | per-copy Theta priors |
| `P_values_*.tsv` | per-read lineage probabilities (allopolyploid) |
| `EM_Final_Expression.tsv` (`EM_Final_Hybrid.tsv` with long reads) | `ClusterID, AlleleID, Ancestry, Initial_Theta, EM_Allocated_Reads, Final_Theta, EM_Iter, EM_LogLik`; single-copy genes have `ClusterID = ORPHAN` |

Final (`<out-dir>/LineageQuant_Final_Merged_TPM.tsv`), one row per gene/copy:

| Column | Description |
|---|---|
| `GeneID` | gene / allele copy |
| `Type` | `Polyploid(EM)` allele copy · `Orphan(EM)` single-copy gene · `Orphan(BAM_Base)` / `Polyploid(idxstats)` counted from idxstats (fallback / ablation) |
| `Length_bp` | longest transcript length |
| `<sample>_RawReads` | EM-assigned fragments (use these for count-based analyses such as DESeq2/edgeR) |
| `<sample>_TPM` | RawReads per kb of `Length_bp`, scaled to one million |

Subgenome-level expression is the sum of `RawReads` over the copies of a
lineage within a cluster (lineages from `--ancestry`).

---

## Key parameters

```
--ploidy {auto,allo}        autopolyploid (copy-level) vs allopolyploid (lineage)
-t, --threads N             CPU threads (default 16)
--theta-threads N           Theta FASTQ-scan workers (default 8; each holds the k-mer index)

Alignment
--aligner {bowtie2,hisat2}  default bowtie2 (RSEM settings); hisat2 = pre-1.0, faster
--max-alignments N          alignments per fragment (default 200 bowtie2 / 15 hisat2)

EM
--em-engine {global,legacy} default global; legacy = pre-1.0 per-cluster EM
--nm-weight {auto,off,λ}    mismatch-count weighting (default auto)
--eff-len {eff,len,off}     effective-length correction (default eff)
--no-pool-identical         do not pool identical-sequence genes
--p-weight-specific/-conflict/-none, --p-alpha, --p-floor   lineage-probability weighting
--theta-prior-lr/-kmer/-cap Dirichlet prior from Theta (long-read anchored clusters)

Theta / lineage
--kmer N                    Theta k-mer length (default 21)
--lineage-kmer N            lineage marker k-mer length (default 31)
--shrinkage F               Theta shrinkage (default 10.0); --adaptive-shrinkage
--min-lr-depth N            min long reads to anchor a cluster's Theta (default 5)
--error-rate F              SNP base error rate (default 0.01)
--no-calibrate-snp, --no-chimeric   MSA route options

--dry-run                   print the resolved route and exit
```

Expert overrides for Step 00 (environment variables): `LQ_BT2_OPTS` replaces the
bowtie2 option string, `LQ_HISAT_OPTS` the HISAT2 one. See `lineagequant run
--help` for the full list.

---

## Run time and resources

**C++ engine (`--engine cpp`).** On large data, add `--engine cpp --lqcore /path/to/lq_cpp/bin/lqcore`.

- Build it once with `bash lq_cpp/hpc/build.sh` (needs GCC >= 7 and zlib).
- bowtie2's output is consumed as it is produced, with no sorting, no BAM and no
  per-read probability tables.
- Results match the Python steps; see CHANGELOG for the one deliberate difference, the
  SNP row of a multi-mapping fragment.
- To split one sample over several nodes, use `lqcore run --chunk i/N --part-out` on each
  node, then `lqcore em --parts`.

bowtie2 with `-k 200` is slower and writes larger BAMs than HISAT2: in our
benchmarks (2–17 million read pairs per sample, 16 threads) alignment took 3–10×
the HISAT2 time and the BAMs were 1.5–3× larger. Building the bowtie2 index is a
one-off cost per transcriptome. The global EM runs in minutes with modest memory,
even for transcriptomes of ~200,000 targets. Use `--aligner
hisat2` when speed matters more than accuracy among near-identical copies.

---

## Reproducing pre-1.0 results

1.0.0 refuses to reuse an output directory aligned with another aligner (earlier
versions always used HISAT2): use a new `--out-dir`, or pass `--aligner hisat2`
to resume. The behaviour of the development versions that preceded 1.0.0
(HISAT2 alignment, per-cluster EM, orphan genes counted from idxstats) is
available as:

```bash
lineagequant run ... --aligner hisat2 --em-engine legacy --eff-len off
```

The EM of the 0.1.0 release (lineage probabilities without an alignment
compatibility matrix) is not retained; see [CHANGELOG.md](CHANGELOG.md).

The `--auto-idxstats` ablation (featureCounts -M -O equivalent counting) requires
`--aligner hisat2`.

---

## Benchmark summary

Twelve species spanning three ploidy classes (three diploids plus a reference
diploid, two autopolyploids, six allopolyploids including the AAB *Musa*
triploids and pentaploid *Rosa canina*), reads simulated from each species'
transcriptome with known read-level truth (three replicates), compared with
Kallisto, Salmon (`--keepDuplicates`) and featureCounts (`-M -O`):

- **Gene level** (fraction of genes within ±10% of truth): best or tied in 9 of
  12 species; the three remaining differences are ≤ 0.13 percentage points.
  Mean absolute percentage error: best or tied in 11 of 12.
- **Subgenome level** (six allopolyploids): lowest homoeolog-expression-bias
  call error in 6 of 6; mean absolute error of subgenome shares lower than
  Kallisto in 6 of 6 and lower than Salmon in 3 of 6 (others within 0.016
  percentage points).

Details are in the accompanying manuscript.

---

## Limitations

- Alignment is end-to-end and pairs must be concordant: fragments with a mate
  outside the annotated transcript (unannotated isoforms or UTR extensions) are
  not counted. A complete transcriptome, with every allele copy as a separate
  record, is assumed.
- Copies identical in sequence cannot be resolved by any read; they are reported
  with equal shares.
- Effective length uses the longest transcript of a gene.

---

## Troubleshooting

- **`tool 'bowtie2' not on PATH`** — activate the conda environment, or install
  with `conda install -c bioconda "bowtie2>=2.4"`.
- **`this output directory was aligned with hisat2`** — the directory comes from
  an earlier run; use a new `--out-dir` (or `--aligner hisat2` to resume it).
- **`--auto-idxstats ... Add --aligner hisat2`** — the ablation is defined on
  HISAT2 alignments only.
- **Resuming** — delete the relevant `.done` marker in the sample directory to
  re-run a single step.

---

## Helper subcommands

```
lineagequant normalize-clusters ...   # standardize allele clusters
lineagequant build-markers ...        # build lineage k-mer markers only
```

---

## Citation

Manuscript in preparation. Until then, please cite this repository and the
version used (`lineagequant --version`).

## Authors

Yi Chen, Gengrui Zhu

## License

This project is licensed under a **Non-Commercial Research License**. Free for
academic and research use only; commercial use is strictly prohibited. See
[LICENSE](LICENSE).
