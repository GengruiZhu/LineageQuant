#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
txt2fasta.py -- convert a merged k-mer table to per-lineage marker FASTAs.

Input : merged TSV WITH header, columns = Kmer  FinalScore  Species
Output: <outdir>/<Species>.fasta, one per lineage; each record is 2 lines:
          >Species-FinalScore
          <kmer>

Exposes write_fastas() for programmatic use; also runnable as a CLI.
Replaces the original txt2fasta.sh (awk) with a pure-Python equivalent so the
package has no awk dependency.
"""

import os
import sys
import argparse


def write_fastas(input_txt, outdir):
    """Write per-lineage FASTAs. Return {lineage_name: fasta_path}."""
    os.makedirs(outdir, exist_ok=True)
    with open(input_txt) as f:
        header = f.readline().rstrip('\n').split('\t')
    if 'Kmer' not in header or 'Species' not in header:
        raise ValueError(f"header needs Kmer & Species: {header}")
    ck = header.index('Kmer')
    csp = header.index('Species')
    cs = header.index('FinalScore') if 'FinalScore' in header else None
    need = max(ck, csp, cs if cs is not None else 0)

    handles = {}
    paths = {}
    try:
        with open(input_txt) as f:
            f.readline()
            for line in f:
                p = line.rstrip('\n').split('\t')
                if len(p) <= need:
                    continue
                sp = p[csp]
                score = p[cs] if cs is not None else "1"
                if sp not in handles:
                    path = os.path.join(outdir, f"{sp}.fasta")
                    handles[sp] = open(path, 'w')
                    paths[sp] = path
                handles[sp].write(f">{sp}-{score}\n{p[ck]}\n")
    finally:
        for h in handles.values():
            h.close()

    for sp, path in paths.items():
        n = sum(1 for l in open(path) if l.startswith('>'))
        sys.stderr.write(f"  {os.path.basename(path)}: {n} kmers\n")
    return paths


def main():
    ap = argparse.ArgumentParser(description="merged k-mer table -> per-lineage FASTAs")
    ap.add_argument("input_txt")
    ap.add_argument("outdir")
    args = ap.parse_args()
    write_fastas(args.input_txt, args.outdir)


if __name__ == "__main__":
    main()
