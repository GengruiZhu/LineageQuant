# -*- coding: utf-8 -*-
"""
normalize_clusters.py -- standardize allele clusters from various tools.

LineageQuant does NOT build clusters. Users produce allele clusters with their
preferred tool (AlleleFinder, jcvi synteny, custom pipelines) and this module
normalizes those outputs into the standard LineageQuant cluster format.

Standard cluster format (one cluster per line):
    ClusterID<TAB>Allele1,Allele2,...

Also derives, when --ancestry-out / --orphan-out are given:
    Ancestry  : GeneID<TAB>LineageName   (lineage = subgenome inferred per input)
    Orphan    : OrphanGrp_N<TAB>GeneID   (single-copy genes, standard cluster shape)

Supported input formats (--format):
  - standard : already ClusterID<TAB>g1,g2,...  (validate + renumber only)
  - jcvi     : jcvi.compara.catalog ortholog .anchors / .lifted.anchors or the
               two-column pair TSV (geneA<TAB>geneB[<TAB>score]); pairs are
               merged into connected components (handles >2 subgenomes).
  - tsv      : generic TSV where each row lists the alleles of one cluster across
               columns; --tsv-cols / --tsv-id-col control parsing.
  - pairs    : generic two-column geneA<TAB>geneB list -> connected components.

Lineage assignment for Ancestry (optional): provide --lineage-of as a Python-
style mapping is avoided; instead --subgenome-tag NAME:REGEX pairs assign a gene
to a lineage by matching its ID against a regex (e.g. A:^GH_A  D:^GH_D). Genes
not matching any tag get lineage "UNKNOWN".
"""
import os
import re
import sys
import argparse
from collections import defaultdict


# --------------------------------------------------------------------------- #
# Union-Find for connected components (pair-based inputs)
# --------------------------------------------------------------------------- #
class UnionFind:
    def __init__(self):
        self.parent = {}

    def find(self, x):
        self.parent.setdefault(x, x)
        root = x
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[x] != root:
            self.parent[x], x = root, self.parent[x]
        return root

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[rb] = ra

    def components(self):
        comp = defaultdict(list)
        for node in list(self.parent):
            comp[self.find(node)].append(node)
        return list(comp.values())


def _clean(g):
    return g.strip().split(".")[0]


# --------------------------------------------------------------------------- #
# Parsers -> list of clusters (each a list of gene IDs)
# --------------------------------------------------------------------------- #
def parse_standard(path):
    clusters = []
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            genes = [_clean(g) for g in parts[1].split(",") if g.strip()]
            if genes:
                clusters.append(genes)
    return clusters


def parse_pairs(path, comment_prefix="#"):
    uf = UnionFind()
    n_pairs = 0
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith(comment_prefix):
                continue
            # jcvi .anchors blocks separated by lines starting with '###'
            if line.startswith("###"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                parts = line.split()
            if len(parts) < 2:
                continue
            a, b = _clean(parts[0]), _clean(parts[1])
            if a and b:
                uf.union(a, b)
                n_pairs += 1
    sys.stderr.write(f"[normalize] parsed {n_pairs} pairs -> components\n")
    return [c for c in uf.components() if len(c) >= 1]


def parse_tsv(path, cols=None, id_col=None):
    clusters = []
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        col_idx = list(range(len(header)))
        if cols:
            col_idx = [header.index(c) if c in header else int(c) for c in cols]
        for line in f:
            if not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            genes = []
            for ci in col_idx:
                if ci < len(p) and p[ci].strip() and p[ci] != ".":
                    genes.append(_clean(p[ci]))
            if id_col is not None:
                idc = header.index(id_col) if id_col in header else int(id_col)
                genes = [g for j, g in enumerate(genes) if col_idx[j] != idc]
            if len(genes) >= 1:
                clusters.append(genes)
    return clusters


# --------------------------------------------------------------------------- #
# Output writers
# --------------------------------------------------------------------------- #
def write_clusters(clusters, out_path, min_size=2):
    n = 0
    with open(out_path, "w") as out:
        for genes in clusters:
            uniq = sorted(set(genes))
            if len(uniq) < min_size:
                continue
            out.write(f"Global_Cluster_{n:06d}\t{','.join(uniq)}\n")
            n += 1
    sys.stderr.write(f"[normalize] wrote {n} clusters (>= {min_size}) -> {out_path}\n")
    return n


def write_orphans(clusters, all_genes, out_path):
    """Genes not in any multi-member cluster -> orphan groups."""
    in_cluster = set()
    for genes in clusters:
        if len(set(genes)) >= 2:
            in_cluster.update(_clean(g) for g in genes)
    orphans = sorted(set(all_genes) - in_cluster) if all_genes else []
    with open(out_path, "w") as out:
        for i, g in enumerate(orphans, start=1):
            out.write(f"OrphanGrp_{i:06d}\t{g}\n")
    sys.stderr.write(f"[normalize] wrote {len(orphans)} orphan genes -> {out_path}\n")
    return len(orphans)


def write_ancestry(clusters, subgenome_tags, out_path):
    """Assign each clustered gene a lineage by regex tag; write Ancestry."""
    compiled = [(name, re.compile(rx)) for name, rx in subgenome_tags]
    genes = sorted({_clean(g) for c in clusters for g in c if len(set(c)) >= 2})
    n_unknown = 0
    with open(out_path, "w") as out:
        out.write("GeneID\tAncestry\n")
        for g in genes:
            lin = "UNKNOWN"
            for name, rx in compiled:
                if rx.search(g):
                    lin = name
                    break
            if lin == "UNKNOWN":
                n_unknown += 1
            out.write(f"{g}\t{lin}\n")
    if n_unknown:
        sys.stderr.write(f"[normalize] WARNING: {n_unknown} genes unmatched by any "
                         f"--subgenome-tag (lineage=UNKNOWN)\n")
    sys.stderr.write(f"[normalize] wrote ancestry for {len(genes)} genes -> {out_path}\n")
    return len(genes)


def main():
    ap = argparse.ArgumentParser(
        description="Standardize allele clusters into LineageQuant format.")
    ap.add_argument("--input", required=True, help="raw cluster file")
    ap.add_argument("--format", required=True,
                    choices=["standard", "jcvi", "tsv", "pairs"])
    ap.add_argument("--clusters-out", required=True)
    ap.add_argument("--ancestry-out", default=None)
    ap.add_argument("--orphan-out", default=None)
    ap.add_argument("--all-genes", default=None,
                    help="file of all gene IDs (one per line) for orphan derivation")
    ap.add_argument("--subgenome-tag", action="append", default=[], metavar="NAME:REGEX",
                    help="assign lineage by regex on gene ID, e.g. A:^GH_A (repeatable)")
    ap.add_argument("--tsv-cols", nargs="+", default=None,
                    help="(tsv) columns holding alleles (names or 0-based indices)")
    ap.add_argument("--tsv-id-col", default=None, help="(tsv) cluster-ID column to exclude")
    ap.add_argument("--min-size", type=int, default=2)
    args = ap.parse_args()

    if args.format == "standard":
        clusters = parse_standard(args.input)
    elif args.format in ("jcvi", "pairs"):
        clusters = parse_pairs(args.input)
    else:
        clusters = parse_tsv(args.input, cols=args.tsv_cols, id_col=args.tsv_id_col)

    write_clusters(clusters, args.clusters_out, min_size=args.min_size)

    if args.ancestry_out:
        tags = []
        for item in args.subgenome_tag:
            if ":" not in item:
                ap.error(f"--subgenome-tag must be NAME:REGEX, got {item}")
            name, rx = item.split(":", 1)
            tags.append((name, rx))
        if not tags:
            sys.stderr.write("[normalize] no --subgenome-tag; all lineages UNKNOWN\n")
        write_ancestry(clusters, tags, args.ancestry_out)

    if args.orphan_out:
        all_genes = None
        if args.all_genes:
            with open(args.all_genes) as f:
                all_genes = [_clean(l) for l in f if l.strip()]
        write_orphans(clusters, all_genes, args.orphan_out)


if __name__ == "__main__":
    main()
