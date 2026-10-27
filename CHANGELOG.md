# Changelog

All notable changes to LineageQuant are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[Semantic Versioning](https://semver.org/).

## [1.1.0] — 2026-10-26

### Highlights

- **One C++ executable, no Python.** `lq` (long name `lineagequant`, the same program; C++17, OpenMP, zlib;
  `bash build.sh`) runs the whole pipeline
  and every step of it. `pip`, `setup.py`, numpy, pandas, Biopython, pysam and cyvcf2 are gone; the external tools
  (bowtie2 / salmon, KMC, gffread, MAFFT, samtools for long reads) are called as before. `bin/lqcore` is the same
  program under the engine name of earlier scripts (`lqcore run` = one sample).
- **Same results as 1.0.0 with the default options.** `tests/run_all.sh` compares every step with the Python of
  1.0.0, byte for byte: EM (36 configurations), the lineage chain and read chunks, Theta (36 configurations), long-read
  counts (80), the reference-preparation steps (normalize-clusters, extract-transcripts, build-markers in three modes,
  snpdb-msa, snpdb-vcf), the expression matrix, and the whole pipeline on four routes. Two deliberate differences:
  the SNP row of a multi-mapping fragment (below) and the lineage list of the MSA SNP database (Fixed).
- **Binary per-sample outputs; the matrix is outside the run.** The quantification of each sample writes
  `Global_Theta_Matrix.lqb` and `<sample>.lqres` (EM + effective lengths + per-reference counts) instead of four text
  tables; the 1.0.0 text files are then written from them byte for byte, outside the timed quantification
  (`--no-text` to skip; `lq export` at any time). The genome-wide matrix (Step 04 of
  1.0.0) is `lq matrix`, run afterwards (same file as `core/aggregator.py`); like kallisto and Salmon,
  the run ends with the per-sample estimates.
- Speed (POJ2878 root, ~38 M read pairs, 32 threads): 6 h 05 m (1.0.0 Python engine) -> 1 h 04 m (bowtie2);
  with `--aligner lqmap --em-accel squarem` a POJ tissue takes 15–30 min.

### Commands

- `lq run` (pipeline), `quant` (one sample; = `lqcore run`), `theta`, `markers`, `em`, `map`, `lrpick`,
  `lrcount`.
- Reference preparation, ported from the Python modules with the same outputs: `normalize-clusters`
  (`cluster/normalize_clusters.py`), `extract-transcripts` (`utils/genome_pair.py`; writes the FASTA index
  itself, so samtools is no longer needed for it), `build-markers` (`lineage_kmer/build_markers.py` +
  `filter_low_complexity.py` + `txt2fasta.py`), `snpdb-msa` (`lineage_snp/build_snp_msa.py` + `export_snp_db.py`,
  MAFFT in parallel), `snpdb-vcf` (`build_snp_vcf.py` + `export_snp_db.py`; plain or bgzip VCF). The SNP
  databases are written directly in the text form the engine reads (no pickle).
- `matrix` (`core/aggregator.py`, plus `--mode idxstats` and `--tpm-length eff`), `export`.

### Removed

- The Python package (`lineagequant/`, `setup.py`, `lineagequant build-core`): the 1.0.0 Python engine and router
  live on in LineageQuant 1.0.0, unchanged.
- From `lq run`: Step 04 (now `lq matrix`), `--tpm-length` (a `matrix` option; accepted with a
  notice), `--engine python`, `--em-engine legacy`, `--aligner hisat2`, `--auto-idxstats`. `--engine cpp`, `--lqcore`,
  `--theta-threads`, `--em-epsilon / -max-iter / -min-iter / -tol / -ll-tol`, `--cross-cluster` and the calibration
  options still parse and do nothing. `--orphan` is accepted (the matrix uses it); `--gff` is needed only by the VCF
  route; `--target-ids` only for allo.
- Step 1b (long-read calibration of the MSA SNP frequencies): it never changed a result (see below).

### Added (the new options are off by default)

- `--eff-len fld`: effective length = length − mean of the observed fragment-length distribution truncated at the
  length + 1 (kallisto / Salmon), instead of length − overall mean + 1. Matters for transcripts shorter than a few
  fragment lengths.
- `--gc-bias`: fragment GC-bias correction, Salmon's foreground / background model with one GC model:
  - foreground = GC percent of the observed fragments (primary proper pair, collected during the scan);
  - background = GC percent of every possible fragment of every expressed target at FLD quantile lengths, weighted by
    the target's expected fragments after a first EM pass;
  - weight = foreground / background (Laplace-smoothed, clamped to [0.01, 100], mean 1 over the background); each
    target's effective length is multiplied by the mean weight of its possible fragments and the EM continues to
    convergence.
  - Synthetic test (`tests/compare_new.sh`): TPM MARD 0.2055 -> 0.1100 under strong GC bias, 0.0986 -> 0.0756 under
    weak bias, unchanged without bias; identical through `quant --gc-hist` + `em --parts`; independent of the thread
    count.
- `--marker-clean`: marker k-mers that occur in a reference copy labelled with another lineage are dropped (k-mers
  absent from the reference are kept). Markers built from progenitor genomes contain such k-mers wherever a progenitor
  diverged at a site the subgenomes share, and they give false lineage calls (rosa4 hold-out with ancestor markers from
  a 0.5–2 % diverged source: wrong-lineage marker calls on 1–2 % of the reads whose source copy is in the reference).
- `matrix --tpm-length eff`: TPM from reads per effective length used by the EM (includes the GC factor) instead of
  reads per kb of the longest transcript.
- `run --theta-sum compensated|plain`: the float sums of Theta as Python >= 3.12 (default) or <= 3.11 computed them
  in 1.0.0 (the router of 1.0.0 followed the running interpreter; last-digit differences).
- `run` writes the 1.0.0 text files after each sample by default, from the binaries and outside the timed
  quantification (the log reports quantification time and text time separately); `--no-text` keeps the binaries
  only. `run --legacy-lineages`.
- `quant / em --eff-out FILE`, `quant --gc-hist`, `markers --marker-clean --fasta FA -a ANCESTRY`, `theta --binary`,
  `quant / em -o S.lqres`; `quant -t` / `em -t` read the binary or the text Theta.
- Part files (`quant --part-out`) are format 02 (+ fragment GC histogram); format 01 parts are still read.
- `bash build.sh --static`: fully static binary for releases (x86-64 Linux, kernel ≥ 3.2, SSE4.2/POPCNT).
- `build.sh`: versioned binaries `bin/lineagequant.<time>`, `bin/lqcore.<time>` with atomic links, older builds kept;
  `make install PREFIX=...`.

### Fixed

- **MSA SNP database: every lineage.** `export_snp_db.py` of 1.0.0 took the lineage list from the first SNP of the
  database, and `build_snp_msa.py` looked the copies up by their id without its last "." field while it selected
  them by the id up to the first ".". A lineage absent from the first aligned cluster, or copies with ids holding
  several dots, lost all SNP evidence without a warning (seen on POJ: a database with SOF and SSP only).
  `snpdb-msa` lists every lineage of the ancestry file and uses one id rule; `--legacy-lineages` (and
  `run --legacy-lineages`) reproduces 1.0.0 exactly. With the first cluster holding every lineage and ids with one
  dot, the two databases are identical. `lq run` warns when it finds a database built the 1.0.0 way.
- `snpdb-msa` merges the MAFFT results in cluster order: the database no longer depends on which worker finished
  first (1.0.0 merged in completion order).
- numpy-compatible sums over more than 8192 values (logs, matrix) follow numpy's blocked summation.

### Added (development builds 0.2–0.5, then `--engine cpp`; now the only engine)

- `--engine cpp`: the C++ engine **lqcore** (separate source tree `lq_cpp/`) runs bowtie2
  and computes, in one pass over its output, the candidates and mismatch counts of every
  fragment, the k-mer and SNP lineage probabilities, their fusion, and the global EM.
  No BAM is sorted or written, and no per-read probability table is written between the
  steps. Same model and arithmetic as the Python steps, tested on synthetic data:
  - The EM is byte-identical to `core/em_global.py`.
  - The k-mer rows are byte-identical to `snipe_kmer.py`; the SNP rows to `snipe_snp.py`.
  - The fused table and the EM output are byte-identical to the
    `merge_poe.py` + `em_global.py` chain.
  - Results do not depend on the number of threads.
- Reads can be split over several jobs. `lqcore run --chunk i/N --part-out` writes the
  evidence of one chunk, and `lqcore em --parts` sums the chunks and runs one EM, because
  equivalence-class counts are additive.
- `lineage_snp/export_snp_db.py`: writes the SNP DB pickle as text for lqcore.
- With `--engine cpp`, Step 01 (Theta) runs as `lqcore theta` instead of `theta/global_theta.py`:
  - output byte-identical (tested on synthetic data with isoforms, IUPAC letters, N, lowercase,
    unequal R1/R2, long-read anchoring, adaptive shrinkage, single-end mode);
  - uses all threads; its own k-mer index cache (`lqcore_theta.k<k>.bin` in the data directory);
  - `--py-sum` reproduces the running interpreter's built-in `sum()` (compensated from Python 3.12 on).
- lqcore 0.4.0, EM memory layout: arrays are first written by the thread that later computes on them (NUMA
  first-touch), work is split into fixed balanced blocks, and EM threads are pinned. Results are unchanged
  (byte-identical to `em_global.py`, independent of the thread count).
- `--em-accel squarem` (`--engine cpp`, default `none`): SQUAREM-accelerated EM. Extrapolations are kept only when
  the likelihood does not fall (monotone), and the result does not depend on the thread count. It converges to the
  same maximum as the plain EM in several-fold fewer iterations; results differ at the level of the tolerance.
- `--aligner lqmap` (`--engine cpp`, experimental): lqcore's built-in aligner with the bowtie2 rules used by
  LineageQuant (ungapped, end-to-end, at most 0.1 x read length mismatches, concordant pairs, fragments up to 1,000 bp).
  - Candidate search: exact 21-mer seeds every 8 bases on both strands; each candidate is verified by its mismatch
    count, then mate rescue, then 1-substitution seeds for divergent fragments.
  - Candidates more than 3 mismatches worse than the best are dropped.
  - Equal to exhaustive search on synthetic polyploid data. Needs the 12-species benchmark before it can become a
    default.
- `--aligner sa` (`--engine cpp`, experimental): Salmon selective alignment proposes candidates; lqcore re-verifies
  them with the same ungapped rules.
- The aggregator reads `<sample>.idxstats.tsv` (written by lqcore) when a sample has no BAM.

### Changed (lqcore)

- **SNP row choice**: a fragment aligned to several transcripts gets the SNP row of one
  transcript, chosen deterministically: fewest mismatches, then most SNPs covered, then
  first in reference order (`--snp-pick`). The Python steps write one row per
  (fragment, transcript), and `merge_poe.py` keeps whichever row came last, an order set by
  `multiprocessing`.

### Fixed (long-read prior)

- `longread/extract_prior.py` counted `samtools idxstats`, i.e. every mapped record:
  secondary and supplementary alignments, and every subread pass of a PacBio ZMW. It then
  collapsed isoforms to allele copies with `=`, so the last isoform overwrote the others.
  Both inflated the Dirichlet weight tau (= `--theta-prior-lr` x cluster long-read count),
  and the secondary records flattened the copy proportions. It now counts molecules:
  - Records are grouped per read; primary and supplementary alignments are the read's
    segments (MAS/Kinnex arrays give several), and secondary alignments join the segment
    they overlap.
  - Within a segment only the best-AS alignments (`--as-tol`) are candidates. Segments
    tied across copies are shared by an EM over equivalence classes; sequence-identical
    copies split evenly.
  - Isoforms of a copy are summed.
  - Subreads of one ZMW count as one molecule (`--zmw-dedup`).
  - Segments with a short or low-identity best alignment are dropped
    (`--min-aln-len`, `--min-identity`).
  - `--legacy-out` writes the old numbers for comparison. The output format is unchanged
    (counts may be fractional).
- `lineage_snp/calibrate_snp.py` reads the BAM sequentially (any sort order, no index needed)
  and skips secondary, QC-fail and duplicate records explicitly. Results are unchanged for
  minimap2 BAMs, whose secondary records carry no sequence.
- **Step 1b (long-read calibration of MSA SNP AFs) is no longer run.**
  - The MSA route scores reads only against each lineage's expected base (`entry["bases"]`),
    and calibration never changes that. `snipe_snp.py --calibrated` was accepted but never
    used. So calibration had no effect on any result and only cost a pass over the long reads.
  - `--no-calibrate-snp`, `--calib-min-depth` and `--calib-det-weight` are still accepted
    for old command lines; they now do nothing.
- The long-read BAM is no longer indexed by the router (nothing needs the index).

### Added (lqcore 0.5.0)

- `lqcore lrcount`: a C++ port of the new `extract_prior.py` with byte-identical output.
  Tested in 200 comparisons on random SAMs covering ties, MAS segments, hard/soft clips, missing AS/NM,
  QC-fail/duplicate records, ZMW names and collation. With `--engine cpp` the router uses it
  for Step 00b.
- `lqcore lrcount` no longer uses `popen`/`pclose` for samtools.
  - It reads the BAM header with `samtools view -H` to its end.
  - It runs samtools in its own process group with SIGPIPE at the default action, and ends that group with
    SIGTERM when it stops reading early.
  - The first version stopped `samtools view` early with `pclose()`; in a PBS job samtools never exited and
    lrcount waited at 0 % CPU (POJ Iso-Seq, 36 h).
  - It now reports progress every ~33.5 M records and every 500 EM iterations.
- `lqcore lrpick`: keeps one read per molecule from PacBio subreads (FASTQ in, FASTQ out).
  - When read names carry `movie/zmw/start_end`, it groups by ZMW.
  - SRA strips those names. In that case consecutive reads count as one molecule when the
    sampled k-mers of one read and of the reverse complement of the other share one
    diagonal, the shorter read lies inside the implied overlap, and the shared k-mers span
    most of it.
  - Pairs 101 reads apart measure the chance level; grouping turns on only when consecutive
    matches clearly exceed it.
  - On simulated MAS-like subreads (8–15 % error): read-level inflation 2.7x → 1.3–1.4x;
    base-level inflation 2.1x → 1.05x; 0.5–0.7 % of molecules merged by mistake.
  - HiFi input is left as it is, and the stats file names the minimap2 preset.
  - `--decide-only` reads only the warm-up reads and writes the decision, so HiFi files can go to
    minimap2 unchanged.
- `lqcore run --debug-frag FILE`: writes one row per fragment.
  - Columns: the candidate copies, how many extra mismatches each has over the best, and the lineage factors
    that enter the EM (NA when the fragment is not Specific).
  - Diagnostic only; results are unchanged.
- `lqcore theta --hits-cache FILE`: stores the per-copy k-mer hits of a read set. Further Θ
  runs on the same reads, with another or no long-read prior, skip the FASTQ scan; output is
  identical.
- `build_snp_msa.py --lq-out`: writes the lqcore text copy of the MSA SNP DB in the same
  pass. With `--engine cpp` the router no longer loads the pickle a second time to export it.

## [1.0.0] — 2026-10-12

First stable release. It supersedes the 0.1.0 release on GitHub and all
intermediate development builds. (The 0.1.0 package metadata carried the
internal string "2.0.0"; from this release on, the package version, the
`--version` output and the Git tag are the same.)

Every change below was driven by a benchmark with read-level ground truth: twelve
species spanning three ploidy classes, reads simulated from each species'
transcriptome, scored against Kallisto, Salmon and featureCounts, plus a second
simulation with sequencing models learned from real RNA-seq of eleven of the
species.

### Highlights

- **bowtie2 is the default short-read aligner** (RSEM settings, `-k 200`,
  ambiguous reference bases not penalised). HISAT2 does not report all equally
  good alignments: for 0.01–0.6% of fragments it reported a single, wrong copy,
  and these reads accumulated on lowly expressed copies next to highly expressed
  homoeologues or paralogues. This was the entire remaining gap to Salmon in the
  benchmark.
- **Global EM engine** (`core/em_global.py`, default): one EM over read
  equivalence classes of all genes, with effective-length correction,
  mismatch-count weighting and identical-sequence pooling.
- Gene-level accuracy (fraction of genes within ±10% of truth) is now the best or
  tied among the four tools in 9 of 12 species; homoeolog-expression-bias calls
  have the lowest error in all six allopolyploids.

### Added

- `--theta-threads` (default 8): worker processes for the Theta FASTQ scan.
- `core/em_global.py --init {theta,uniform}` (default `theta`): starting values of the
  EM. `uniform` splits every multi-mapping fragment equally among its candidates; it
  exists to check that copy-level results do not depend on the starting point.
- `--aligner {bowtie2,hisat2}` (default `bowtie2`) and `--max-alignments N`
  (default 200 for bowtie2, 15 for hisat2).
  bowtie2 options: `--sensitive --dpad 0 --gbar 99999999 --mp 1,1 --np 0
  --n-ceil L,0,0.5 --score-min L,0,-0.1 -I 1 -X 1000 --no-mixed
  --no-discordant -k 200`. The bowtie2 index (`transcript_index_bt2.*`) is built
  once in `--data-dir` and shared by all samples.
- `--em-engine {global,legacy}` (default `global`; `v6`/`v5` accepted as
  aliases). The global engine:
  - runs one EM over all genes, so every fragment carries total weight 1;
  - corrects for effective length (`--eff-len`, default `eff`);
  - weights candidate alignments by mismatch count (`--nm-weight`, default
    `auto`: each extra mismatch multiplies the likelihood by λ = (p/3)/(1−p),
    with p estimated from the best alignments; `off` or a fixed λ are accepted);
  - pools genes whose transcripts are identical in sequence and splits them
    evenly after convergence (`--no-pool-identical` to disable);
  - keeps the Theta start values, the long-read Dirichlet prior and the
    lineage-probability weighting of the legacy engine.
- `--eff-len {off,len,eff}` and `--cross-cluster {keep,drop,split}` (legacy engine).
- `--p-weight-specific`, `--p-weight-conflict`, `--p-weight-none`, `--p-alpha`,
  `--p-floor`: weighting of per-read lineage probabilities by their evidence
  status.
- `--theta-prior-lr`, `--theta-prior-kmer`, `--theta-prior-cap`: Theta as a
  Dirichlet prior in the M-step (long-read anchored clusters).
- `--auto-idxstats`: ablation that reproduces featureCounts-style counting
  (requires `--aligner hisat2`).
- Each sample directory records its aligner in `aligner.txt`; re-running an
  output directory with a different aligner stops with a clear message instead
  of silently mixing BAMs.
- Expert environment overrides for Step 00: `LQ_BT2_OPTS`, `LQ_HISAT_OPTS`,
  `LQ_BT2_K`, `LQ_HISAT_K`.
- `CHANGELOG.md`; Chinese usage guide `docs/USAGE_zh.md`.

### Changed

- `--ploidy auto` now runs Theta and EM. Before, the auto route was mapping →
  idxstats → TPM, numerically identical to featureCounts `-M -O`; that path is
  kept only as the `--auto-idxstats` ablation.
- The EM always builds its compatibility matrix from the BAM; the lineage step
  supplies probabilities only (`W = M · P^w · θ`). In 0.1.0 lineage
  probabilities replaced the alignment evidence altogether.
- The Theta k-mer index is written to `--data-dir` and shared by all samples
  (`global_theta.py --index_dir`).
- `snipe_kmer.py` (Step 02, k-mer lineage): Beta-smoothed probabilities with the
  status decided by the amount of evidence (`Specific`, `Conflict_Weighted`,
  `Weak_Evidence`, `Ambiguous_None`); mates de-duplicated; canonical marker sets
  detected automatically; new `N_Hits` and `Best_AS` columns. `--legacy`
  reproduces the 0.1.0 behaviour.
- `build_markers.py`: optional `--polyploid NAME:FASTA` and `--mode
  {ancestor,ancestor_x_poly,poly}`, so that markers must also be specific among
  the polyploid's own subgenome transcripts.
- Aggregator (Step 04): EM counts are summed per gene across transcripts; with
  the global engine single-copy genes are labelled `Orphan(EM)`.
- Command-line help is in English; configuration errors are reported without a
  traceback; unknown options to `lineagequant run` are rejected instead of being
  ignored.

### Fixed

- Reads compatible with copies in several clusters were counted once per cluster
  (legacy `--cross-cluster keep`), and multi-mapping reads of single-copy genes
  were counted once per target (idxstats). The global EM counts every fragment
  exactly once.
- Genes listed in neither `clusters.tsv` nor `orphans.tsv` were missing from the
  output matrix; they are now quantified (and any gene not covered by the EM
  falls back to idxstats).
- Gene IDs containing dots were truncated at the first dot in the aggregator
  (`split(".")[0]`); only the last `.N` suffix is removed now.
- Allopolyploid runs without any lineage signal produced an empty
  `P_values_auto.tsv`, which zeroed every EM allocation; the EM now runs without
  lineage probabilities in that case.
- Copies with identical sequence were split 0/1 by the unconstrained EM when a
  few reads happened to be reported for only one of them; they are pooled now.
- Fragments overlapping runs of N in the reference transcriptome were discarded
  by the aligner; N is no longer penalised.
- The Theta step (`global_theta.py`) queued the whole FASTQ in memory before its
  workers caught up and started one worker per `--threads`, each holding the
  copy-unique k-mer index; large samples (e.g. 45 million read pairs against a
  333,000-transcript reference) were killed by the out-of-memory killer. At most
  2 × workers read chunks are now in flight, and the number of workers is set by
  the new `--theta-threads` (default 8). Results are unchanged.
- The SNP lineage step (`snipe_snp.py`) rebuilt a full Python copy of the SNP
  database in every worker process, so its memory grew with `--threads`; on hybrid
  sugarcane at 32 threads it reached about 800 GB and stalled in swap. The database
  is now converted once into a compact per-transcript index (sorted positions and
  small numeric arrays) that all workers share, and only the SNP observations of
  each fragment are kept instead of the read sequence and full position list.
  Output is byte-identical to the previous version (tested on both `--mode msa`
  and `--mode population`, with clipped and gapped alignments).
- The k-mer + SNP fusion (`merge_poe.py`) read the per-fragment status from the last
  column of each input, which since the per-read evidence columns were added
  (`N_Hits`, `Best_AS`, `N_SNPs`) is no longer the status. Every fragment therefore
  got the lowest reliability tier and the experts were weighted by entropy alone.
  The status column is now located by name.
- The global EM reported "EM converged" also when it stopped at `--max-iter` before
  meeting `--tol`. It now says which of the two happened and, in the second case, how
  many targets were still changing and what share of the fragments they hold.
- `--markers` accepted an empty lineage name (`--markers :file.fa`) and silently
  kept only the last file given for a repeated name; both are now rejected.

### Compatibility

- Output directories created by earlier versions were aligned with HISAT2. Use a
  new `--out-dir` with 1.0.0, or pass `--aligner hisat2` to resume them.
- Behaviour of the development versions before 1.0.0:
  `--aligner hisat2 --em-engine legacy --eff-len off`.
- bowtie2 alignment is slower than HISAT2 (3–10× in the benchmark) and the BAMs
  are larger (1.5–3×).
- New dependency: bowtie2 ≥ 2.4 (in `environment.yml`); hisat2 is optional.

## [0.1.0] — 2026-06-25

Initial public release.

- Auto-routing pipeline that integrates short-read allopolyploid quantification,
  long-read/hybrid quantification and lineage-specific k-mer marker
  construction (KMC N-way subtraction).
- Ploidy fork (`--ploidy auto|allo`), three lineage routes (ancestor k-mer,
  population VCF SNP, MSA SNP) with Product-of-Experts fusion, long-read Theta
  anchoring and MSA SNP calibration.
- HISAT2 alignment (`-k 15 --no-mixed --no-discordant`), per-cluster EM,
  TPM aggregation; `normalize-clusters` and `build-markers` helper commands.

[1.0.0]: https://github.com/GengruiZhu/LineageQuant/releases/tag/v1.0.0
[0.1.0]: https://github.com/GengruiZhu/LineageQuant/releases/tag/v0.1.0
