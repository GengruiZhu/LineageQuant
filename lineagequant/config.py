#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
config.py -- configuration container and input resolution for LineageQuant.

Holds all resolved parameters and the parsed input groups that the router uses
to decide which modules to run. Path resolution (relative -> absolute, NAME:PATH
parsing) happens here so the router stays focused on dispatch logic.
"""

import os


def _abs(p, base=None):
    if p is None:
        return None
    if os.path.isabs(p):
        return p
    return os.path.abspath(os.path.join(base, p)) if base else os.path.abspath(p)


class Config:
    """Resolved configuration for a LineageQuant run."""

    def __init__(self, args):
        # ---- core paths ----
        self.data_dir = _abs(args.data_dir)
        self.out_dir = _abs(args.out_dir)
        self.fasta = _abs(args.fasta)            # transcriptome FASTA (mRNA)
        self.fasta_gff = _abs(args.gff)          # reference GFF (reserved/mapping)
        self.clusters = _abs(args.clusters)
        self.ancestry = _abs(args.ancestry)
        self.orphan = _abs(args.orphan)
        self.target_ids = _abs(args.target_ids)

        # ---- ploidy ----
        self.ploidy = args.ploidy                # 'auto' | 'allo'

        # ---- lineage assignment inputs (allopolyploid) ----
        # ancestors: {name: [(ancestor_fasta, ancestor_gff), ...]}  (>=1 sample per lineage)
        # Repeating --ancestor with the same NAME adds another progenitor sample
        # to that lineage (widens k-mer evidence). All samples of one lineage are
        # merged before KMC subtraction.
        self.ancestors = {}
        for item in (args.ancestor or []):
            # format NAME:fa,gff
            if ":" not in item or "," not in item:
                raise ValueError(f"--ancestor must be NAME:fasta,gff, got: {item}")
            name, paths = item.split(":", 1)
            fa, gff = paths.split(",", 1)
            self.ancestors.setdefault(name, []).append((_abs(fa), _abs(gff)))

        # Pre-built lineage markers (skip building from ancestors).
        # markers: {name: marker_fasta}.  --markers NAME:path  (repeatable)
        # If given, the KMC build step is skipped and these are used directly.
        # NAMES should match the lineage names in --ancestry.
        self.prebuilt_markers = {}
        for item in (getattr(args, "markers", None) or []):
            if ":" not in item:
                raise ValueError(f"--markers must be NAME:fasta, got: {item}")
            name, path = item.split(":", 1)
            name = name.strip()
            if not name:
                raise ValueError(f"--markers needs a lineage name before ':', got: {item}")
            if name in self.prebuilt_markers:
                raise ValueError(f"--markers: lineage '{name}' given more than once")
            self.prebuilt_markers[name] = _abs(path)

        # population VCF route
        self.vcf = _abs(args.vcf)
        self.sample_lists = []                   # list of 'NAME:path'
        for item in (args.sample_lists or []):
            if ":" not in item:
                raise ValueError(f"--sample-lists must be NAME:path, got: {item}")
            name, path = item.split(":", 1)
            self.sample_lists.append(f"{name}:{_abs(path)}")

        # MSA route.
        # Two ways to trigger it:
        #   --phased-fasta <fa> : alleles' CDS live in a SEPARATE file
        #   --msa-snp           : alleles' CDS already in the main --fasta (usual case);
        #                         avoids re-supplying the same sequences.
        self.msa_snp = bool(getattr(args, "msa_snp", False))
        if getattr(args, "phased_fasta", None) and self.msa_snp:
            raise ValueError(
                "--phased-fasta and --msa-snp are two ways to trigger the same MSA "
                "route; give only one. Use --msa-snp if the main --fasta already "
                "holds all allele copies; use --phased-fasta only for a separate "
                "CDS file.")
        if self.msa_snp:
            # MSA route runs on the main transcriptome directly.
            self.phased_fasta = self.fasta
        else:
            self.phased_fasta = _abs(args.phased_fasta)

        # long-read
        self.longread_bam = _abs(args.longread)

        # ---- runtime / tuning ----
        self.threads = args.threads
        self.kmer = args.kmer
        self.theta_threads = max(1, min(getattr(args, "theta_threads", 8) or 8, args.threads))
        self.lineage_kmer = args.lineage_kmer
        self.shrinkage = args.shrinkage
        self.adaptive_shrinkage = args.adaptive_shrinkage
        self.no_paired_bonus = args.no_paired_bonus
        self.min_lr_depth = args.min_lr_depth
        self.error_rate = args.error_rate

        # SNP calibration
        self.no_calibrate_snp = args.no_calibrate_snp
        self.calib_min_depth = args.calib_min_depth
        self.calib_det_weight = args.calib_det_weight
        self.no_chimeric = args.no_chimeric

        # ---- short-read alignment (Step 00) ----
        self.aligner = getattr(args, "aligner", "bowtie2") or "bowtie2"
        self.salmon = getattr(args, "salmon", None) or "salmon"
        default_k = 15 if self.aligner == "hisat2" else 200
        self.max_alignments = getattr(args, "max_alignments", None) or default_k

        # ---- ablation: pre-1.0 auto behaviour (mapping -> idxstats -> TPM) ----
        # Equivalent to featureCounts -M -O; no Theta, no EM. For comparisons only.
        self.auto_idxstats = bool(getattr(args, 'auto_idxstats', False))

        # ---- exponent on lineage probabilities, by Lineage_Status ----
        self.p_weight_specific = getattr(args, 'p_weight_specific', 1.0)
        self.p_weight_conflict = getattr(args, 'p_weight_conflict', 0.0)
        self.p_weight_none = getattr(args, 'p_weight_none', 0.0)
        self.p_alpha = getattr(args, 'p_alpha', 0.0)
        self.p_floor = getattr(args, 'p_floor', 1e-3)

        # ---- theta0 as a Dirichlet prior in the M-step ----
        # tau = theta_prior_lr x cluster long-read count (LongRead_Anchored)
        #     = theta_prior_kmer                         (Kmer_Inferred)
        self.theta_prior_lr = getattr(args, 'theta_prior_lr', 1.0)
        self.theta_prior_kmer = getattr(args, 'theta_prior_kmer', 0.0)
        self.theta_prior_cap = getattr(args, 'theta_prior_cap', 0.0)

        # EM
        self.em_epsilon = args.em_epsilon
        self.em_max_iter = args.em_max_iter
        self.em_min_iter = args.em_min_iter
        self.em_tol = args.em_tol
        self.em_ll_tol = args.em_ll_tol
        eng = getattr(args, "em_engine", "global") or "global"
        self.em_engine = {"v6": "global", "v5": "legacy"}.get(eng, eng)
        self.em_pool_identical = not getattr(args, "no_pool_identical", False)
        self.nm_weight = str(getattr(args, "nm_weight", "auto") or "auto")
        self.eff_len = getattr(args, "eff_len", None)
        self.cross_cluster = getattr(args, "cross_cluster", None)
        # C++ engine (lqcore): one pass per sample, no BAM
        self.engine = getattr(args, "engine", "python") or "python"
        import shutil as _sh
        self.lqcore = (getattr(args, "lqcore", None) or os.environ.get("LQCORE") or _sh.which("lqcore") or "lqcore")
        self.em_accel = getattr(args, "em_accel", "none") or "none"

        self._validate()

    def _validate(self):
        if self.auto_idxstats and self.ploidy != "auto":
            raise ValueError("--auto-idxstats is only meaningful with --ploidy auto")
        if self.auto_idxstats and self.aligner != "hisat2":
            raise ValueError(
                "--auto-idxstats counts every reported alignment in full; with bowtie2 "
                "(which also reports near-optimal alignments) that is not featureCounts-"
                "equivalent. Add --aligner hisat2 to reproduce the pre-1.0 ablation.")
        if self.aligner in ("lqmap", "sa") and self.engine != "cpp":
            raise ValueError(f"--aligner {self.aligner} needs --engine cpp")
        if self.em_accel != "none" and self.engine != "cpp":
            raise ValueError("--em-accel needs --engine cpp")
        if self.engine == "cpp":
            if self.aligner not in ("bowtie2", "lqmap", "sa"):
                raise ValueError("--engine cpp supports --aligner bowtie2 | lqmap | sa (not hisat2)")
            if self.em_engine != "global":
                raise ValueError("--engine cpp implements the global EM engine only")
            if self.auto_idxstats:
                raise ValueError("--engine cpp does not support the --auto-idxstats ablation")
        if self.max_alignments < 1:
            raise ValueError("--max-alignments must be >= 1")
        if self.nm_weight not in ("auto", "off"):
            try:
                v = float(self.nm_weight)
            except ValueError:
                raise ValueError(f"--nm-weight must be auto, off or a number, got: {self.nm_weight}")
            if not 0.0 < v <= 1.0:
                raise ValueError("--nm-weight as a number must be in (0, 1]")
        if self.ploidy == "allo":
            has_kmer = bool(self.ancestors or self.prebuilt_markers)
            has_route = bool(has_kmer or self.vcf or self.phased_fasta)
            if not has_route:
                raise ValueError(
                    "Allopolyploid mode needs at least one lineage route: "
                    "--ancestor (k-mer), --markers (pre-built k-mer), "
                    "--vcf (population SNP), or --phased-fasta (MSA SNP).")
            if self.vcf and self.phased_fasta:
                raise ValueError(
                    "--vcf and --phased-fasta are mutually exclusive "
                    "(both provide SNP signal). Pick one.")
            if self.vcf and not self.sample_lists:
                raise ValueError("--vcf requires --sample-lists NAME:path [...]")

    def describe_route(self):
        """Human-readable summary of the route the router will take."""
        lines = [f"ploidy: {self.ploidy}"]
        if self.aligner == "bowtie2":
            lines.append(f"  alignment: bowtie2 (RSEM settings, N not penalised, "
                         f"<= {self.max_alignments} alignments/fragment)")
        elif self.aligner == "lqmap":
            lines.append("  alignment: lqmap (built in: exact k-mer seeds, ungapped mismatch count, "
                         "concordant pairs, candidates within best + 3 mismatches)")
        elif self.aligner == "sa":
            lines.append(f"  alignment: sa ({self.salmon} selective alignment -> ungapped re-verification, "
                         f"candidates within best + 3 mismatches)")
        else:
            lines.append(f"  alignment: hisat2 (-k {self.max_alignments} --no-mixed --no-discordant)")
        if self.engine == "cpp":
            lines.append(f"  engine: lqcore ({self.lqcore}) -- {self.aligner} -> candidates + lineage evidence -> "
                         f"EM in one pass, no BAM" + (f" | EM acceleration: {self.em_accel}" if self.em_accel != "none" else ""))
        if not self.auto_idxstats:
            if self.em_engine == "global":
                lines.append(f"  EM: global over all genes | eff-len={self.eff_len} | "
                             f"mismatch weighting={self.nm_weight} | "
                             f"identical-sequence pooling={'on' if self.em_pool_identical else 'off'}")
            else:
                lines.append("  EM: legacy per-cluster engine (orphans by idxstats)")
        if self.ploidy == "auto":
            if self.auto_idxstats:
                lines.append("  [ablation] mapping -> idxstats -> TPM "
                             "(no Theta, no EM; == featureCounts -M -O)")
            else:
                mode = "long-read anchored" if self.longread_bam else "k-mer"
                lines.append(f"  copy-level: Theta ({mode}) -> EM (M from BAM) -> TPM")
            return "\n".join(lines)
        routes = []
        if self.prebuilt_markers:
            routes.append(f"pre-built k-mer markers "
                          f"({len(self.prebuilt_markers)} lineages, build skipped)")
        if self.ancestors:
            n_samp = sum(len(v) for v in self.ancestors.values())
            routes.append(f"ancestor k-mer ({len(self.ancestors)} lineages, "
                          f"{n_samp} progenitor sample(s))")
        if self.vcf:
            routes.append("population VCF SNP")
        if self.phased_fasta:
            src = " [main fasta]" if self.msa_snp else " [separate phased CDS]"
            routes.append(f"MSA SNP{src}")
        lines.append(f"  lineage routes: {', '.join(routes)}")
        if (self.ancestors or self.prebuilt_markers) and (self.vcf or self.phased_fasta):
            lines.append("  -> PoE fusion of k-mer + SNP signals")
        theta = "long-read anchored (01_hybrid)" if self.longread_bam else "k-mer (01)"
        lines.append(f"  Theta: {theta}")
        return "\n".join(lines)
