#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
filter_low_complexity.py -- drop low-complexity k-mers from a merged
(Kmer<TAB>FinalScore<TAB>Species) file.

Criteria (cheap-first for speed):
  - any mononucleotide 6+ run (AAAAAA / TTTTTT / CCCCCC / GGGGGG)
  - AT or GC content > threshold (default 80%)
  - single base (AAAA...)
  - full dinucleotide repeat (ATATAT...)
  - any 2-10 bp motif repeated >= 3x (e.g. CCAATA CCAATA CCAATA)

Exposes filter_file() for programmatic use; also runnable as a CLI.
"""

import re
import sys
import argparse
from collections import defaultdict

_REPEAT_RE = re.compile(r'([ACGT]{2,10})\1{2,}')


def is_low_complexity(kmer, k, thr):
    if 'AAAAAA' in kmer or 'TTTTTT' in kmer or 'CCCCCC' in kmer or 'GGGGGG' in kmer:
        return True
    a = kmer.count('A'); t = kmer.count('T'); c = kmer.count('C'); g = kmer.count('G')
    if (a + t) > thr or (c + g) > thr:
        return True
    if len(set(kmer)) == 1:
        return True
    di = kmer[:2]
    if di * (k // 2) == kmer[:(k // 2) * 2]:
        return True
    if _REPEAT_RE.search(kmer):
        return True
    return False


def filter_file(input_file, output_file, k=31, at_gc_frac=0.8, verbose=True):
    """Filter a merged Kmer/FinalScore/Species TSV; return (kept, removed) totals."""
    thr = int(at_gc_frac * k)
    with open(input_file) as f:
        header = f.readline().rstrip('\n').split('\t')
    if 'Kmer' not in header or 'Species' not in header:
        raise ValueError(f"header needs 'Kmer' and 'Species': {header}")
    ck = header.index('Kmer'); csp = header.index('Species')
    need = max(ck, csp)

    kept = defaultdict(int)
    removed = defaultdict(int)
    with open(input_file) as f, open(output_file, 'w') as out:
        out.write('\t'.join(header) + '\n')
        f.readline()
        for line in f:
            p = line.rstrip('\n').split('\t')
            if len(p) <= need:
                continue
            sp = p[csp]
            if is_low_complexity(p[ck], k, thr):
                removed[sp] += 1
            else:
                out.write(line)
                kept[sp] += 1

    tot_k = sum(kept.values()); tot_r = sum(removed.values())
    if verbose:
        sys.stderr.write("[filter_low_complexity] per-lineage kept / removed:\n")
        for sp in sorted(set(list(kept) + list(removed))):
            k_ = kept[sp]; r_ = removed[sp]; tot = k_ + r_
            pct = 100.0 * r_ / tot if tot else 0.0
            sys.stderr.write(
                f"  {sp:16s} kept={k_:>10d}  removed={r_:>10d}  ({pct:.1f}% LC)\n")
        g = tot_k + tot_r
        sys.stderr.write(
            f"  {'TOTAL':16s} kept={tot_k:>10d}  removed={tot_r:>10d}  "
            f"({100.0*tot_r/g if g else 0:.1f}%)\n")
    return tot_k, tot_r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--input_file', required=True, help="merged Kmer FinalScore Species")
    ap.add_argument('--output_file', required=True)
    ap.add_argument('--k', type=int, default=31)
    ap.add_argument('--at_gc_frac', type=float, default=0.8)
    args = ap.parse_args()
    filter_file(args.input_file, args.output_file, k=args.k, at_gc_frac=args.at_gc_frac)


if __name__ == '__main__':
    main()
