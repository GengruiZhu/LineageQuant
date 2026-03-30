# LineageQuant

**Allele-resolved RNA-seq quantification for polyploid and hybrid genomes**

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python 3.7+](https://img.shields.io/badge/Python-3.7%2B-brightgreen)](https://www.python.org/)

---

## Description

LineageQuant is a bioinformatics pipeline for allele-resolved expression
quantification in polyploid and hybrid genomes. It deconvolves RNA-seq
signals from multi-copy gene families across ancestral subgenomes using a
three-stage framework:

1. **K-mer Theta engine** – Subgenome-specific k-mer probes are constructed
   from allele sequences to compute an initial per-allele read abundance
   score (Theta).
2. **BAM lineage P-value assignment** – Each read overlapping a target locus
   is scanned against ancestral k-mer dictionaries (forward and
   reverse-complement) and assigned a per-lineage probability vector.
3. **Expectation-Maximization (EM) deconvolution** – Ambiguously mapping
   reads are iteratively reallocated using the Theta priors and per-read
   lineage probabilities, converging on final allele-resolved expression
   values.

Orphan (single-copy) genes are quantified in parallel by featureCounts and
merged into a single unified TPM output matrix.

LineageQuant supports **any number (≥ 2) of ancestral lineages** and is
applicable to polyploid crops and hybrid organisms such as sugarcane,
wheat, and cotton.

---

## Pipeline overview

```
FASTQ (R1 + R2)
     │
     ▼  Stage 0 – HISAT2 + featureCounts
  Sorted BAM  +  raw count table
     │
     ▼  Stage 1 – K-mer Theta engine
  Theta_Matrix.tsv          (per-allele initial abundance)
     │
     ▼  Stage 2 – BAM lineage P-value assignment
  P_values.tsv              (per-read lineage probability vectors)
     │
     ▼  Stage 3 – EM deconvolution
  EM_Final_Expression.tsv   (allele-resolved allocated read counts + Theta)
     │
     ▼  Stage 4 – TPM matrix aggregation
  LineageQuant_Final_Merged_TPM.tsv
```

All stages support **checkpoint / resume**: if an output file already
exists it is skipped automatically.

---

## Installation

### Option A – Conda (recommended)

```bash
git clone https://github.com/<your-username>/LineageQuant.git
cd LineageQuant
conda env create -f environment.yml
conda activate lineagequant
pip install -e .
```

### Option B – pip (tools must be installed separately)

```bash
pip install biopython pysam pandas
# also install: hisat2, samtools, subread (featureCounts), gffread
pip install -e .
```

---

## Quick start

```bash
lineagequant \
  --gff        /path/to/genome.gff3 \
  --fasta      /path/to/genome.fasta \
  --clusters   /path/to/allele_clusters.tsv \
  --gene-ids   /path/to/target_gene_ids.txt \
  --ancestry   /path/to/gene_ancestry.tsv \
  --orphan-genes /path/to/single_copy_genes.tsv \
  --ancestors  SubA:subA_k15.fa  SubB:subB_k15.fa  SubC:subC_k15.fa \
  --data-dir   /path/to/fastq_directory \
  --out-dir    /path/to/results \
  --threads    16
```

---

## Input file formats

| Argument | Format description |
|---|---|
| `--gff` | Standard GFF3; genes annotated with `type=gene` and `ID=` attribute |
| `--fasta` | Indexed whole-genome FASTA (`.fai` index optional, `gffread` will use it) |
| `--clusters` | Tab-separated: `ClusterID<TAB>Allele1,Allele2,...` (one cluster per line) |
| `--gene-ids` | Plain text, one gene ID per line |
| `--ancestry` | Tab-separated: `GeneID<TAB>LineageName` |
| `--orphan-genes` | Same format as `--clusters`; lists single-copy genes to quantify with featureCounts |
| `--ancestors` | `NAME:PATH` pairs, one per ancestral lineage (≥2 required); PATH points to a FASTA file where each record is a single k-mer sequence |
| FASTQ | Paired-end, gzipped, named `<sample>.R1.fq.gz` / `<sample>.R2.fq.gz` |

---

## Output files

Per-sample (inside `<out_dir>/<sample>/`):

| File | Description |
|---|---|
| `<sample>.sorted.bam` | HISAT2-aligned, sorted BAM |
| `<sample>_featureCounts.txt` | Raw featureCounts table |
| `master_index.pkl` | K-mer probe index checkpoint (Stage 1) |
| `Theta_Matrix.tsv` | Per-allele initial Theta values |
| `P_values.tsv` | Per-read lineage probability vectors |
| `EM_Final_Expression.tsv` | EM-allocated read counts and final Theta |

Final merged output (inside `<out_dir>/`):

| File | Description |
|---|---|
| `LineageQuant_Final_Merged_TPM.tsv` | Whole-genome TPM matrix across all samples |

---

## All parameters

```
Reference files (required):
  --gff            Whole-genome GFF3 annotation
  --fasta          Whole-genome reference FASTA
  --clusters       Allele cluster TSV
  --gene-ids       Target gene ID list
  --ancestry       Gene ancestry table
  --orphan-genes   Orphan gene list
  --ancestors      Lineage k-mer FASTA files: NAME:PATH [NAME:PATH ...]

Input / Output:
  -d, --data-dir   FASTQ directory  (default: ./data)
  -o, --out-dir    Results directory (default: ./results)
  --index-prefix   HISAT2 index prefix (default: <data-dir>/genome_index)

Runtime:
  -t, --threads    CPU threads                (default: 16)
  --kmer           K-mer size for Theta engine (default: 31)
  --lineage-kmer   K-mer size for lineage files (default: 15)
  --epsilon        EM Laplace pseudo-count     (default: 0.01)
  --max-iter       Maximum EM iterations       (default: 200)
  --tol            EM convergence tolerance    (default: 1e-5)
```

---

## Authors

Yi Chen, Gengrui Zhu

---

## License

[MIT](LICENSE)
