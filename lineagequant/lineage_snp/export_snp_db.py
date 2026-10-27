# -*- coding: utf-8 -*-
"""
export_snp_db.py -- write a SNP DB pickle as plain text for the C++ engine (lqcore).

    python export_snp_db.py --snp_db SNP_DB_msa.pkl --mode msa --out SNP_DB_msa.lq.tsv

Output (grouped by transcript in the DB's own order, positions ascending):
    #LQSNP  mode=msa         lineages=L1,L2,...  n=N
    tid  pos  EXPECTED          (msa: one character per lineage, A/C/G/T or '.' = never matches)
    tid  pos  REF  ALT  AF,AF   (population: allele frequency per lineage, exact repr)

The expected bases follow snipe_snp.py exactly: entry["bases"], or, without bases, alt if the
lineage AF > 0.5 else ref; lineage order = snipe_snp._infer_lineages (sorted).
"""
import argparse
import os
import pickle
import sys


def infer_lineages(snp_db):
    for entry in snp_db.values():
        if "bases" in entry and entry["bases"]:
            return sorted(entry["bases"].keys())
        keys = [k for k in entry if k not in ("bases", "ref", "alt")]
        if keys:
            return sorted(keys)
    return []


def one(b):
    b = str(b).upper()
    return b if b in ("A", "C", "G", "T") else "."


def write_lq(db, mode, out):
    """write the SNP DB dict as lqcore text; returns (n_snps, n_transcripts, lineages)"""
    lins = infer_lineages(db)
    if not lins:
        sys.exit("[FATAL] could not infer lineage names from the SNP DB")
    by = {}
    for (tid, pos), e in db.items():
        by.setdefault(tid, []).append((int(pos), e))
    tmp = out + ".tmp"
    n = 0
    with open(tmp, "w") as f:
        f.write(f"#LQSNP\tmode={mode}\tlineages={','.join(lins)}\tn={sum(len(v) for v in by.values())}\n")
        for tid, rows in by.items():
            rows.sort(key=lambda x: x[0])
            for pos, e in rows:
                if mode == "msa":
                    exp = e.get("bases", {})
                    if not exp:
                        ref_a, alt_a = e["ref"].upper(), e["alt"].upper()
                        exp = {l: (alt_a if e.get(l, 0.0) > 0.5 else ref_a) for l in lins}
                    f.write(f"{tid}\t{pos}\t{''.join(one(exp.get(l, 'N')) for l in lins)}\n")
                else:
                    ref_a, alt_a = e["ref"].upper(), e["alt"].upper()
                    afs = ",".join(repr(float(e.get(l, 0.5))) for l in lins)
                    f.write(f"{tid}\t{pos}\t{ref_a if len(ref_a) == 1 else '.'}\t{alt_a if len(alt_a) == 1 else '.'}\t{afs}\n")
                n += 1
    os.replace(tmp, out)
    return n, len(by), lins


def main():
    ap = argparse.ArgumentParser(description="SNP DB pickle -> lqcore text")
    ap.add_argument("--snp_db", required=True)
    ap.add_argument("--mode", choices=["msa", "population"], default="msa")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    with open(a.snp_db, "rb") as f:
        db = pickle.load(f)
    n, nt, lins = write_lq(db, a.mode, a.out)
    print(f"[SUCCESS] {n:,} SNPs on {nt:,} transcripts, lineages {lins} -> {a.out}")


if __name__ == "__main__":
    main()
