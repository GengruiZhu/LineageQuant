#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
genome_pair.py -- Genome/Annotation compatibility layer.

Binds a (FASTA, GFF) pair and makes them interoperable regardless of how the
sequence names differ between the two files. Solves the common pitfall where a
FASTA uses chromosome names like "Chr01" while its GFF uses accession names like
"CM001740.1" for the same physical sequence.

Responsibilities:
  1. Clean FASTA headers (strip whitespace/TAB suffixes, e.g. GWH format).
  2. Build a .fai index and select chromosome-scale sequences (>= min_chr_len),
     excluding scaffolds/contigs/organelles.
  3. Auto-build a seqid mapping between FASTA names and GFF column-1 names:
       - exact-name match first
       - length-based + fuzzy-name fallback (Chr01 <-> CM001740.1)
  4. Filter the GFF to chromosome seqids and drop strand '?' (organelle
     trans-splicing) features.
  5. Provide a unified extract_transcripts() that runs gffread on the cleaned
     pair, independent of the original naming.

This module shells out to `samtools faidx` and `gffread`, which must be on PATH
(both are declared in environment.yml).
"""

import os
import re
import sys
import shutil
import subprocess
from collections import OrderedDict


def _run(cmd, **kw):
    return subprocess.run(cmd, check=True,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                          text=True, **kw)


class GenomeAnnotationPair:
    """A FASTA + GFF bound together with automatic seqid reconciliation."""

    def __init__(self, fasta, gff, workdir, name="genome",
                 min_chr_len=1_000_000, drop_unknown_strand=True):
        if not os.path.exists(fasta):
            raise FileNotFoundError(f"FASTA not found: {fasta}")
        if not os.path.exists(gff):
            raise FileNotFoundError(f"GFF not found: {gff}")
        for tool in ("samtools", "gffread"):
            if shutil.which(tool) is None:
                raise RuntimeError(f"Required tool '{tool}' not found on PATH.")

        self.src_fasta = os.path.abspath(fasta)
        self.src_gff = os.path.abspath(gff)
        self.name = name
        self.min_chr_len = min_chr_len
        self.drop_unknown_strand = drop_unknown_strand
        self.workdir = os.path.abspath(workdir)
        os.makedirs(self.workdir, exist_ok=True)

        # Derived artifacts (filled by prepare())
        self.clean_fasta = os.path.join(self.workdir, f"_{name}.genome.fa")
        self.clean_gff = os.path.join(self.workdir, f"_{name}.clean.gff")
        self.fasta_lengths = OrderedDict()   # fasta_seqid -> length
        self.chrom_ids = []                  # fasta seqids >= min_chr_len
        self.seqid_map = {}                  # gff_seqid -> fasta_seqid
        self._prepared = False

    # ------------------------------------------------------------------ #
    # FASTA header cleaning + indexing
    # ------------------------------------------------------------------ #
    def _clean_fasta_headers(self):
        """Strip whitespace/TAB suffix from FASTA headers (keep first token)."""
        if os.path.exists(self.clean_fasta) and os.path.getsize(self.clean_fasta):
            return
        with open(self.src_fasta) as fin, open(self.clean_fasta, "w") as fout:
            for line in fin:
                if line.startswith(">"):
                    fout.write(re.sub(r"[ \t].*", "", line.rstrip("\n")) + "\n")
                else:
                    fout.write(line)

    def _index_and_select_chroms(self):
        """Build .fai and select sequences >= min_chr_len."""
        fai = self.clean_fasta + ".fai"
        if not (os.path.exists(fai) and os.path.getsize(fai)):
            _run(["samtools", "faidx", self.clean_fasta])
        self.fasta_lengths.clear()
        with open(fai) as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 2:
                    self.fasta_lengths[parts[0]] = int(parts[1])
        self.chrom_ids = [sid for sid, ln in self.fasta_lengths.items()
                          if ln >= self.min_chr_len]

    # ------------------------------------------------------------------ #
    # GFF seqid collection
    # ------------------------------------------------------------------ #
    def _collect_gff_seqids(self):
        """Return {gff_seqid: max_coordinate_seen} to estimate lengths."""
        gff_span = {}
        with open(self.src_gff) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                p = line.rstrip("\n").split("\t")
                if len(p) < 5:
                    continue
                sid = p[0]
                try:
                    end = int(p[4])
                except ValueError:
                    continue
                if sid not in gff_span or end > gff_span[sid]:
                    gff_span[sid] = end
        return gff_span

    # ------------------------------------------------------------------ #
    # seqid mapping: fasta <-> gff
    # ------------------------------------------------------------------ #
    @staticmethod
    def _norm_token(s):
        """Normalize a seqid for fuzzy matching: lowercase, drop non-alnum,
        strip common chromosome prefixes."""
        s = s.lower()
        s = re.sub(r"[^a-z0-9]", "", s)
        s = re.sub(r"^(chromosome|chrom|chr|scaffold|contig|seq)", "", s)
        return s

    def _build_seqid_map(self, gff_span):
        """
        Map every GFF seqid to a FASTA seqid.

        Strategy (in order):
          1. Exact match (gff seqid == fasta seqid).
          2. Normalized-token match (Chr01 ~ chr01 ~ 01).
          3. Length match: GFF max-coordinate close to a FASTA length
             (within 0.5%), unique candidate only.
        Unmapped GFF seqids are simply not added (their features get dropped).
        """
        self.seqid_map = {}
        fasta_ids = list(self.fasta_lengths.keys())
        fasta_norm = {self._norm_token(s): s for s in fasta_ids}

        # length lookup with small tolerance
        def length_match(gff_end):
            best = None
            for fid, flen in self.fasta_lengths.items():
                # GFF max coord should be <= sequence length, and close to it
                if gff_end <= flen and (flen - gff_end) / flen <= 0.005:
                    if best is not None:
                        return None  # ambiguous
                    best = fid
            return best

        for gid, gend in gff_span.items():
            # 1. exact
            if gid in self.fasta_lengths:
                self.seqid_map[gid] = gid
                continue
            # 2. normalized token
            nt = self._norm_token(gid)
            if nt in fasta_norm:
                self.seqid_map[gid] = fasta_norm[nt]
                continue
            # 3. length match
            lm = length_match(gend)
            if lm is not None:
                self.seqid_map[gid] = lm
                continue
            # else: unmapped -> dropped

    # ------------------------------------------------------------------ #
    # GFF rewriting (rename seqids to FASTA names + filter)
    # ------------------------------------------------------------------ #
    def _rewrite_gff(self):
        """Write clean GFF: rename seqids to FASTA names, keep only chromosome
        seqids, drop strand '?' features."""
        chrom_set = set(self.chrom_ids)
        with open(self.src_gff) as fin, open(self.clean_gff, "w") as fout:
            for line in fin:
                if line.startswith("#"):
                    continue
                p = line.rstrip("\n").split("\t")
                if len(p) < 8:
                    continue
                gid = p[0]
                fid = self.seqid_map.get(gid)
                if fid is None or fid not in chrom_set:
                    continue
                if self.drop_unknown_strand and p[6] == "?":
                    continue
                p[0] = fid  # rename to FASTA seqid
                fout.write("\t".join(p) + "\n")

    # ------------------------------------------------------------------ #
    # Public API
    # ------------------------------------------------------------------ #
    def prepare(self):
        """Run the full reconciliation pipeline. Idempotent."""
        if self._prepared:
            return self
        self._clean_fasta_headers()
        self._index_and_select_chroms()
        gff_span = self._collect_gff_seqids()
        self._build_seqid_map(gff_span)
        self._rewrite_gff()
        self._prepared = True

        n_map = len(self.seqid_map)
        n_chr = len(self.chrom_ids)
        n_all = len(self.fasta_lengths)
        sys.stderr.write(
            f"[genome_pair:{self.name}] chromosomes {n_chr}/{n_all} "
            f"(>= {self.min_chr_len:,}bp), seqid map {n_map} gff->fasta\n")
        return self

    def extract_transcripts(self, out_fa, feature="mRNA", id_prefix=None):
        """
        Run gffread on the cleaned (fasta, gff) to produce a transcript FASTA.

        feature : kept for API symmetry; gffread -w extracts spliced exons
                  (mRNA/transcript). CDS-only extraction would use -x.
        id_prefix : if set, prepend "<prefix>__" to every transcript ID
                    (used to tag ancestor origin for lineage tracing).
        """
        if not self._prepared:
            self.prepare()
        log = os.path.join(self.workdir, f"_gffread_{self.name}.log")
        with open(log, "w") as lg:
            subprocess.run(
                ["gffread", "-w", out_fa, "-g", self.clean_fasta, self.clean_gff],
                check=True, stdout=lg, stderr=lg, text=True)
        if id_prefix:
            tmp = out_fa + ".tmp"
            with open(out_fa) as fin, open(tmp, "w") as fout:
                for line in fin:
                    if line.startswith(">"):
                        sid = line[1:].split()[0]
                        fout.write(f">{id_prefix}__{sid}\n")
                    else:
                        fout.write(line)
            os.replace(tmp, out_fa)
        return out_fa

    def n_transcripts(self, fa):
        return sum(1 for l in open(fa) if l.startswith(">"))


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(
        description="Reconcile a FASTA+GFF pair and extract transcripts.")
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--gff", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--name", default="genome")
    ap.add_argument("--out", help="If set, extract transcripts to this FASTA.")
    ap.add_argument("--id-prefix", default=None)
    ap.add_argument("--min-chr-len", type=int, default=1_000_000)
    args = ap.parse_args()

    pair = GenomeAnnotationPair(
        args.fasta, args.gff, args.workdir, name=args.name,
        min_chr_len=args.min_chr_len).prepare()

    if args.out:
        pair.extract_transcripts(args.out, id_prefix=args.id_prefix)
        print(f"[done] {pair.n_transcripts(args.out)} transcripts -> {args.out}")
