# LineageQuant

**Haplotype-resolved expression quantification for polyploid genomes**

[![Version 1.1.0](https://img.shields.io/badge/version-1.1.0-orange)](CHANGELOG.md)
[![License: Non-Commercial Research](https://img.shields.io/badge/License-Non--Commercial-blue.svg)](LICENSE)
[![C++17](https://img.shields.io/badge/C%2B%2B17-one%20executable-blue)](docs/ENGINE.md)

LineageQuant quantifies **subgenome / allele-copy expression** in polyploids. It is
an expression *quantifier*, not a read classifier: every homoeologous or allelic
copy is a separate target, reads that are compatible with several copies are
apportioned by expectation–maximisation (EM), and the output is a read-count and
TPM matrix at copy level that can be summed to subgenome or gene level.

It integrates short-read allopolyploid quantification, long-read/hybrid
quantification and lineage-specific k-mer marker construction into a single
**auto-routing** program: you declare the ploidy type and supply whatever inputs
you have, and the pipeline runs the matching steps.

> **1.1.0 is one C++ executable** (`lq`, long name `lineagequant`; no Python, no pip). With the default
> options it gives the results of 1.0.0 — every step is compared byte for byte with the
> 1.0.0 Python steps by `tests/run_all.sh` — at a fraction of the time and memory
> (POJ2878 root: 6 h → 1 h). The quantification writes binary per-sample results (`.lqres`, `.lqb`);
> the 1.0.0 text files are written from them after each sample, outside the timed quantification
> (`--no-text`: binary only). The expression matrix is built afterwards by `lq matrix`.
> New, opt-in: `--eff-len fld`, `--gc-bias`, `--marker-clean`, `matrix --tpm-length eff`.
> LineageQuant 1.0.0 (Python) stays available unchanged. See [CHANGELOG.md](CHANGELOG.md).

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
- [Commands](#commands)
- [Key parameters](#key-parameters)
- [Run time and resources](#run-time-and-resources)
- [Equivalence with 1.0.0](#equivalence-with-100)
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

Lineage markers add information only where the alignment cannot separate the copies (copies missing from
the reference, collapsed assemblies). On diploids, autopolyploids and allopolyploids with a complete
reference they do not improve the result, so no marker route is switched on by default.

### Long reads are an enhancement, not a route

A long-read BAM (`--longread`) anchors the Theta prior with physical counts
(01_hybrid). Without long reads the k-mer Theta engine is used.

---

## Pipeline overview

```
reference preparation (once, shared by all samples; skipped when done)
   ├── bowtie2 / Salmon / lqmap index of the transcriptome
   ├── [allo] SNP database        snpdb-msa (MAFFT) | snpdb-vcf
   ├── [allo] lineage markers     extract-transcripts (gffread) + build-markers (KMC) | --markers
   └── [long reads] Step 00b      lrcount: molecules per allele copy

per sample
   ├── Step 01     theta    per-copy prior (EM start values)            -> Global_Theta_Matrix.lqb
   └── Step 00-03  quant    one pass, no BAM:                            -> <sample>.lqres
         bowtie2 -k 200 (RSEM settings) streamed in | lqmap (built in) | sa (Salmon)
           → candidates + mismatch counts per fragment
           → [allo] ancestor k-mer + SNP (MSA / VCF) lineage probabilities → PoE
           → global EM over all genes: effective length, mismatch-count weighting,
             identical-sequence pooling (opt-in: --eff-len fld, --gc-bias)

afterwards (not part of the run)
   └── matrix       genome-wide RawReads + TPM matrix                   -> LineageQuant_Final_Merged_TPM.tsv
```

All steps are checkpointed (`.done` markers); re-running the same command resumes from the last
completed step.

---

## How reads are assigned

**Alignment (Step 00).** By default reads are aligned to the transcriptome with
bowtie2 using RSEM's settings — end-to-end, ungapped, `--mp 1,1 --score-min
L,0,-0.1 -I 1 -X 1000 --no-mixed --no-discordant -k 200` — with ambiguous
reference bases not penalised (`--np 0 --n-ceil L,0,0.5`). This reports every
equally good copy of a fragment and also near-optimal ones. `quant` reads bowtie2's
output as it is produced (nothing is sorted or written). HISAT2 (the pre-1.0
aligner, available in 1.0.0) returns only a subset of equally good alignments: in
our benchmarks it reported a single, wrong copy for 0.01–0.6% of fragments, and
those reads accumulated on lowly expressed copies next to highly expressed
homoeologues or paralogues — an error no EM can repair. `--aligner lqmap` (the
built-in aligner, same rules) and `--aligner sa` (Salmon proposes candidates,
LineageQuant re-verifies them) are faster alternatives.

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
  length + 1). `--eff-len fld` (1.1.0) uses the mean of the observed fragment-length
  distribution truncated at L (kallisto / Salmon). `--gc-bias` (1.1.0) rescales each
  effective length by the expected sampling weight of the target's fragments under a
  fragment-GC model learned from the sample (Salmon's foreground / background model).
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
conda env create -f environment.yml      # external tools + compiler (or install them yourself)
conda activate lineagequant
bash build.sh                            # seconds; GCC >= 7 (or clang >= 6) and zlib
export PATH=$PWD/bin:$PATH

lq version                     # lq 1.1.0
lq run --help
```

`build.sh` installs `bin/lq.<time>`, `bin/lineagequant.<time>` and `bin/lqcore.<time>` (three copies of one program)
and switches the links `bin/lq` / `bin/lineagequant` / `bin/lqcore` to them atomically; older builds are kept (overwriting a binary that runs on
another node of a network file system crashes it there). `make && make install PREFIX=$HOME/.local` also works.
On a cluster, build on a login node of the cluster itself: a binary built on a newer system may need a newer
glibc than the compute nodes have.

**Prebuilt binary** (nothing to compile): the [GitHub release](https://github.com/GengruiZhu/LineageQuant/releases/latest)
has a fully static build (`bash build.sh --static`: no glibc, libstdc++ or zlib needed at run time; any x86-64 Linux
with kernel ≥ 3.2 and a CPU with SSE4.2/POPCNT, i.e. any x86-64 server CPU since about 2009). The external tools
still come from `environment.yml` or your own installation.

```bash
tar xzf lq-1.1.0-linux-x86_64-static.tar.gz       # lq, lineagequant, lqcore (hard links of one file), LICENSE, CHANGELOG.md
export PATH=$PWD/lq-1.1.0-linux-x86_64-static:$PATH
lq version                                         # lq 1.1.0
```

`lq` is the command; `lineagequant` is the same with the long name. `lqcore` is the same program under the engine
name of 0.x / 1.0.0 scripts: `lqcore run` runs one sample (= `lq quant`), `lq run` runs the pipeline. Everything
else is the same under all three names.

External tools, by route (`lq run` checks the ones it needs and names any that are missing):

| Tool | Needed for |
|---|---|
| bowtie2 ≥ 2.4 | `--aligner bowtie2` (default) |
| salmon | `--aligner sa` |
| kmc + kmc_tools, gffread | `--ancestor` (marker build) |
| mafft | `--msa-snp` / `--phased-fasta` |
| samtools | `--longread` (reading the BAM) |
| minimap2 | making the long-read BAM |

`--aligner lqmap` needs no external tool.

---

## Quick start

### Step 1 — standardize your allele clusters

LineageQuant does not build clusters. Produce them with your tool of choice
(GraphAllele, AlleleFinder, jcvi synteny, custom) and standardize the output:

```bash
lq normalize-clusters \
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
lq run \
  --fasta      transcriptome.fasta \
  --clusters   clusters.tsv \
  --ancestry   ancestry.tsv \
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
lq run \
  --fasta transcriptome.fasta --clusters clusters.tsv --ancestry ancestry.tsv \
  --ploidy auto --data-dir data --out-dir results --threads 32
```

`--data-dir` holds the paired FASTQ files (one pair per sample) and receives the
shared indices. Add `--dry-run` to print the resolved route without executing.

### Step 3 — expression matrix

```bash
lq matrix -r results -o orphans.tsv            # -> results/LineageQuant_Final_Merged_TPM.tsv
```

More examples are in [`examples/`](examples/); a Chinese usage guide is in
[`docs/USAGE_zh.md`](docs/USAGE_zh.md).

---

## Routing table

The pipeline selects a path from the inputs you provide. Lineage names are read
from the ancestry file (any number N).

### Autopolyploid (`--ploidy auto`)

| # | short | long | route |
|---|---|---|---|
| A1 | ✓ | | Theta k-mer (01) → EM |
| A2 | ✓ | ✓ | long-read Theta (01_hybrid) → EM |

### Allopolyploid (`--ploidy allo`)

| # | short | long | ancestor | VCF | MSA | route |
|---|---|---|---|---|---|---|
| B1 | ✓ | | ✓ | | | Theta k-mer + ancestor-k-mer lineage → EM |
| B2 | ✓ | | | ✓ | | Theta k-mer + VCF SNP lineage → EM |
| B3 | ✓ | | ✓ | ✓ | | Theta k-mer + ancestor-k-mer + VCF SNP → **PoE** → EM |
| B4 | ✓ | ✓ | ✓ | | | long-read Theta + ancestor-k-mer lineage → EM |
| B5 | ✓ | ✓ | | | ✓ | long-read Theta + MSA SNP → EM |
| B6 | ✓ | ✓ | ✓ | | ✓ | long-read Theta + ancestor-k-mer + MSA SNP → **PoE** → EM |
| B7 | ✓ | ✓ | | ✓ | | long-read Theta + VCF SNP lineage → EM |
| B8 | ✓ | ✓ | ✓ | ✓ | | long-read Theta + ancestor-k-mer + VCF SNP → **PoE** → EM |

The MSA route works without long reads as well. With `--ploidy allo` but no lineage
signal, the EM runs without lineage probabilities (equivalent to the auto route).

---

## Input file formats

| Argument | Format |
|---|---|
| `--fasta` | Transcriptome (mRNA) FASTA — every allele copy a separate record |
| `--clusters` | `ClusterID<TAB>Allele1,Allele2,...` (one cluster per line) |
| `--ancestry` | `GeneID<TAB>LineageName` |
| `--target-ids` | one gene ID per line (allopolyploid) |
| `--orphan` | `OrphanID<TAB>GeneID` (single-copy genes; used by `matrix`) |
| `--ancestor` | `NAME:fasta,gff` per ancestral lineage (≥2; allopolyploid) |
| `--markers` | `NAME:fasta` pre-built lineage markers (skips the KMC build) |
| `--vcf` + `--gff` | population VCF (plain or bgzip, not BCF) and the reference GFF3 (VCF route) |
| `--sample-lists` | `NAME:file` per lineage, sample names one per line (with `--vcf`) |
| `--phased-fasta` / `--msa-snp` | MSA SNP route: separate phased CDS FASTA / use the main `--fasta` |
| `--longread` | long-read BAM (`minimap2 -a -N 50 -p 0.95` → transcriptome), any sort order |
| FASTQ (in `--data-dir`) | paired, gzipped; `.R1.fq.gz`, `_R1.fastq.gz`, `_1.fq.gz`, `_R1_001.fastq.gz`, … auto-detected |

Gene IDs are derived from transcript IDs by removing the last `.N` suffix, so
`Gene.1` and `Gene.2` are isoforms of `Gene`.

---

## Output

Per sample (`<out-dir>/<sample>/`):

| File | Description |
|---|---|
| `<sample>.lqres` | binary result: EM per target (Theta, EM-assigned fragments, final Theta), effective lengths the EM used, per-reference alignment counts |
| `Global_Theta_Matrix.lqb` | binary per-copy Theta priors |
| `EM_Final_Expression.tsv` (`EM_Final_Hybrid.tsv` with long reads) | text of the EM: `ClusterID, AlleleID, Ancestry, Initial_Theta, EM_Allocated_Reads, Final_Theta, EM_Iter, EM_LogLik`; single-copy genes have `ClusterID = ORPHAN` |
| `<sample>.idxstats.tsv`, `<sample>.efflen.tsv`, `Global_Theta_Matrix.tsv` | text of the per-reference counts (samtools idxstats format), effective lengths (`Target, Length, EffectiveLength, GC_Factor`) and Theta |
| `<sample>.bowtie2.log` / `.salmon.log` | aligner summary |
| `aligner.txt` | aligner used for this sample directory |

The quantification writes only the binary files; the text files (the files of 1.0.0, byte for byte) are
written from them after each sample and are not counted in the quantification time — the log reports both
(`[DONE] quantification ... | text files ..., not counted`). `--no-text` skips them; `lq export`
writes them from a `.lqres` / `.lqb` at any time. Layouts: [docs/FORMATS.md](docs/FORMATS.md).

Expression matrix (`lq matrix`, `<out-dir>/LineageQuant_Final_Merged_TPM.tsv`), one row per gene/copy,
the same file as Step 04 of 1.0.0:

| Column | Description |
|---|---|
| `GeneID` | gene / allele copy |
| `Type` | `Polyploid(EM)` allele copy · `Orphan(EM)` single-copy gene · `Orphan(BAM_Base)` counted from the per-reference counts (fallback) |
| `Length_bp` | longest transcript length |
| `<sample>_RawReads` | EM-assigned fragments (use these for count-based analyses such as DESeq2/edgeR) |
| `<sample>_TPM` | RawReads per kb of `Length_bp` (`--tpm-length eff`: per effective length), scaled to one million |

Subgenome-level expression is the sum of `RawReads` over the copies of a
lineage within a cluster (lineages from `--ancestry`). When the reference lacks copies (collapsed or
unannotated), copy-level shares are biased while lineage-level shares stay robust; read results at the
lineage level first.

---

## Commands

```
pipeline
  lq run                  the whole pipeline (reference preparation + every sample)
  lq quant                one sample (= lqcore run): alignment -> lineage evidence -> EM
  lq theta | markers | em | map | lrpick | lrcount     the engine steps
reference preparation (run by `run` when needed; usable on their own)
  lq normalize-clusters   allele clusters from other tools -> standard format (+ ancestry, orphans)
  lq extract-transcripts  progenitor genome + GFF -> transcript FASTA (seqid reconciliation, gffread)
  lq build-markers        lineage-specific k-mer markers (KMC; modes ancestor / ancestor_x_poly / poly)
  lq snpdb-msa            SNP database from MAFFT alignments of the copies of each cluster
  lq snpdb-vcf            SNP database from a population VCF
outside the run
  lq matrix               expression matrix from the .lqres files (or 1.0.0 text outputs)
  lq export               .lqres / .lqb -> the text files of 1.0.0
```

`<command> -h` prints the options of each. Splitting one sample over several nodes:
`lq quant ... --chunk i/N --part-out s.part<i>` per node, then
`lq em ... --parts s.part1 ... -o s.lqres` (equivalence-class counts are additive).

---

## Key parameters

```
--ploidy {auto,allo}        autopolyploid (copy-level) vs allopolyploid (lineage)
-t, --threads N             CPU threads (default 16)
--em-accel {none,squarem}   SQUAREM-accelerated EM (default none = the 1.0.0 EM)

Alignment
--aligner {bowtie2,lqmap,sa} default bowtie2 (RSEM settings); lqmap = built in; sa = Salmon candidates
--max-alignments N          bowtie2 alignments per fragment (default 200)

EM
--nm-weight {auto,off,λ}    mismatch-count weighting (default auto)
--eff-len {eff,fld,len,off} effective-length correction (default eff; fld = truncated FLD mean)
--gc-bias                   fragment GC-bias correction of the effective lengths (1.1.0, off)
--no-pool-identical         do not pool identical-sequence genes
--p-weight-specific/-conflict/-none, --p-alpha, --p-floor   lineage-probability weighting
--theta-prior-lr/-kmer/-cap Dirichlet prior from Theta (long-read anchored clusters)

Theta / lineage
--kmer N                    Theta k-mer length (default 21)
--lineage-kmer N            lineage marker k-mer length (default 31)
--shrinkage F               Theta shrinkage (default 10.0); --adaptive-shrinkage
--min-lr-depth N            min long reads to anchor a cluster's Theta (default 5)
--theta-sum {compensated,plain}  float sums of Theta (default compensated = 1.0.0 under Python >= 3.12;
                            plain = 1.0.0 under Python <= 3.11; last-digit differences)
--error-rate F              SNP base error rate (default 0.01)
--no-chimeric               SNP routes: keep lineage-chimeric reads
--marker-clean              drop marker k-mers present in a reference copy of another lineage (1.1.0, off)
--legacy-lineages           MSA SNP database with the lineage list of 1.0.0 (see Equivalence)

--no-text                   binary outputs only (default: the 1.0.0 text files too, written after each sample, not timed)
--dry-run                   print the resolved route and exit
```

Expert override for Step 00 (environment variable): `LQ_BT2_OPTS` replaces the
bowtie2 option string. `lq run --help` lists every option; options of 1.0.0 without
an effect in 1.1.0 (`--engine cpp`, `--lqcore`, `--theta-threads`, the calibration and legacy-EM options)
still parse.

---

## Run time and resources

- bowtie2's output is consumed as it is produced, with no sorting, no BAM and no
  per-read probability tables. POJ2878 root (~38 M read pairs, 32 threads): 6 h 05 m with
  the 1.0.0 Python engine, 1 h 04 m with the C++ engine (bowtie2); a POJ tissue with
  `--aligner lqmap --em-accel squarem` takes 15–30 min and ~11 GB.
- Human transcriptome (~200,000 transcripts, 40 M read pairs, `--aligner lqmap --em-accel squarem`): about
  2 min per sample, 3.5 GB.
- bowtie2 with `-k 200` dominates the run time with the default aligner; `--aligner lqmap`
  applies the same alignment rules several times faster.
- The per-sample outputs are a few MB of binary; the expression matrix is not part of the run
  (as kallisto and Salmon write per-sample estimates and leave the matrix to the user).

---

## Equivalence with 1.0.0

`bash tests/run_all.sh /path/to/LineageQuant_v1.0.0/lineagequant` (minutes, synthetic data; needs python3 with
numpy and pandas; MAFFT, KMC, gffread, Biopython, pysam and cyvcf2 are replaced by stand-ins used by both
sides) compares every step with the Python of 1.0.0, byte for byte:

| Test | Compared with (1.0.0) |
|---|---|
| EM, 36 configurations, 1/2/4 threads | `core/em_global.py` |
| k-mer rows, SNP rows, fused table, EM, read chunks | `snipe_kmer.py`, `snipe_snp.py`, `merge_poe.py` + `em_global.py` |
| Theta, 36 configurations (+ binary + export) | `theta/global_theta.py` |
| long-read counts, 80 comparisons | `longread/extract_prior.py` |
| lqmap | an exhaustive ungapped aligner |
| binary outputs → export; matrix (incl. 30,000 genes) | the text writers; `core/aggregator.py` |
| normalize-clusters, extract-transcripts, build-markers (3 modes), snpdb-msa, snpdb-vcf | the Python modules |
| the whole pipeline, 4 routes (ancestor + MSA, auto, markers + VCF + SQUAREM, bowtie2) | `cli.py` / `router.py` |

Two deliberate differences:

- The SNP row of a fragment aligned to several transcripts comes from one fixed transcript (fewest
  mismatches → most SNPs covered → lowest index); 1.0.0 took the last row of an intermediate file, which
  depended on process scheduling.
- The MSA SNP database lists **every lineage of the ancestry file**. 1.0.0 took the lineages of the
  first SNP in the database only, and looked the copies up by their id without its last "." field; a lineage
  missing from the first cluster (or ids with several dots) silently lost its SNP evidence.
  `--legacy-lineages` reproduces 1.0.0 exactly (tested); when the first cluster holds every lineage
  the two are identical (tested). The pipeline warns when it finds a database built the 1.0.0 way.

Pre-1.0 behaviour (HISAT2, per-cluster EM) and the `--auto-idxstats` ablation need LineageQuant 1.0.0;
`lq matrix --mode idxstats` gives the same kind of counts (every reported alignment counted in
full, featureCounts `-M -O`) from 1.1.0 outputs.

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
- Effective length uses the longest transcript of a gene (allele copies with several
  isoforms are one EM target).
- `--gc-bias` uses one fragment-GC model (Salmon additionally conditions it on the GC of
  the fragment ends and models 5′/3′ sequence bias).
- `--subgenome-tag` regexes use ECMAScript syntax (the usual `^PREFIX` / `_A$` patterns are the same as
  in Python).

---

## Troubleshooting

- **`tool 'bowtie2' not on PATH`** — activate the conda environment, or install
  with `conda install -c bioconda "bowtie2>=2.4"` (or use `--aligner lqmap`).
- **`GLIBC_2.xx not found`** — build on the machine (cluster) where it will run: `bash build.sh`, or use the
  static binary of the GitHub release.
- **`this output directory was aligned with ...`** — the directory comes from a run
  with another aligner; use a new `--out-dir`.
- **`--r1 is an option of one sample`** — `lq run` is the pipeline; one sample is
  `lq quant` (or `lqcore run`).
- **Resuming** — rerun the same command; delete a `.done` marker to re-run that step. Sample directories
  written by 1.0.0 with `--engine cpp` are resumed too (their text Theta / EM files are used).

---

## Citation

Manuscript in preparation. Until then, please cite this repository and the
version used (`lq version`).

## Authors

Yi Chen, Gengrui Zhu

## License

This project is licensed under a **Non-Commercial Research License**. Free for
academic and research use only; commercial use is strictly prohibited. See
[LICENSE](LICENSE).
