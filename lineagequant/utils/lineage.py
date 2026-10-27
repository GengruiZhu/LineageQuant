#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lineage.py -- lineage-name handling, fully parameterized.

Replaces the hard-coded LINEAGE_NAMES = ["SOF","SRO","SSP"] scattered across the
original scripts. Lineage names are read dynamically from the Ancestry file
(column 2), de-duplicated and sorted, so the software supports any number N of
lineages with arbitrary names.

Ancestry file format:  GeneID<TAB>LineageName   (one per line, optional header
line starting with 'GeneID').
"""

import sys


def load_ancestry(path):
    """Return {clean_gene_id: lineage_name}."""
    amap = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("GeneID") or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 2:
                amap[parts[0].split(".")[0]] = parts[1].upper()
    return amap


def lineage_names_from_ancestry(path):
    """Return sorted unique lineage names found in an Ancestry file."""
    names = sorted(set(load_ancestry(path).values()))
    if not names:
        raise ValueError(f"No lineage names found in ancestry file: {path}")
    return names


def parse_name_path_pairs(items):
    """
    Parse a list of 'NAME:PATH' strings into an ordered dict {name: path}.
    Used for --ancestor, --lineages, --sample-lists.
    """
    out = {}
    for item in items:
        if ":" not in item:
            raise ValueError(f"Expected NAME:PATH, got: {item}")
        name, path = item.split(":", 1)
        out[name] = path
    return out


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="Inspect lineage names in an ancestry file.")
    ap.add_argument("ancestry")
    args = ap.parse_args()
    names = lineage_names_from_ancestry(args.ancestry)
    sys.stderr.write(f"[lineage] {len(names)} lineages: {names}\n")
    print(",".join(names))
