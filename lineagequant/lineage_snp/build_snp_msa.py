# -*- coding: utf-8 -*-
"""
build_snp_msa.py -- phased-genome MSA SNP extractor (Step 1a, MSA route).

For each allele cluster, runs MAFFT MSA on the phased CDS of its alleles and
records columns where lineages differ as SNPs. Lineage names are read from the
ancestry file (any number N), replacing the hard-coded SOF/SRO/SSP.

SNP entry: {"bases": {lineage: base}, "ref": b0, "alt": b1, <lineage>: AF...}
Keyed by (full_transcript_id, transcript_0based_pos) for BAM coordinate match.

Multiprocess: N parallel MAFFT workers. Shells out to `mafft`.
"""
import os
import sys
import pickle
import argparse
import subprocess
import tempfile
import shutil
import multiprocessing as mp
from collections import defaultdict
from Bio import SeqIO, AlignIO


def get_now():
    import datetime
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def check_mafft():
    try:
        r = subprocess.run(["mafft", "--version"], stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE, text=True, timeout=10)
        return r.returncode == 0
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False


_seq_dict = _ancestry = _temp_root = _mafft_opts = _lineages = None


def _init_worker(seq_dict, ancestry, temp_root, mafft_opts, lineages):
    global _seq_dict, _ancestry, _temp_root, _mafft_opts, _lineages
    _seq_dict, _ancestry, _temp_root, _mafft_opts, _lineages = \
        seq_dict, ancestry, temp_root, mafft_opts, lineages


def _process_cluster_batch(cluster_batch):
    local_db = {}
    stats = {"clusters": 0, "snps_found": 0, "skipped_gaps": 0, "skipped_identical": 0}
    worker_temp = tempfile.mkdtemp(dir=_temp_root, prefix="mafft_")
    for cluster in cluster_batch:
        allele_map = cluster["allele_map"]
        if len(allele_map) < 2:
            continue
        full_ids = list(allele_map.values())
        seqs = [_seq_dict[fid] for fid in full_ids]
        if all(s == seqs[0] for s in seqs[1:]):
            stats["skipped_identical"] += 1
            continue
        stats["clusters"] += 1
        tmp_fa = os.path.join(worker_temp, f"{cluster['cid']}.fa")
        with open(tmp_fa, "w") as fa:
            for fid in full_ids:
                fa.write(f">{fid}\n{_seq_dict[fid]}\n")
        try:
            r = subprocess.run(["mafft"] + _mafft_opts + [tmp_fa],
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               text=True, timeout=120)
        except subprocess.TimeoutExpired:
            os.remove(tmp_fa)
            continue
        os.remove(tmp_fa)
        if r.returncode != 0 or not r.stdout.strip():
            continue
        aln_path = os.path.join(worker_temp, f"{cluster['cid']}.aln")
        with open(aln_path, "w") as af:
            af.write(r.stdout)
        try:
            alignment = AlignIO.read(aln_path, "fasta")
        except Exception:
            os.remove(aln_path)
            continue
        recs = list(alignment)
        col_to_pos = defaultdict(dict)
        for rec in recs:
            up = 0
            for ci, base in enumerate(str(rec.seq)):
                if base != "-":
                    col_to_pos[rec.id][ci] = up
                    up += 1
        for ci in range(alignment.get_alignment_length()):
            bases = [str(rec.seq)[ci] for rec in recs]
            if "-" in bases or "N" in bases:
                stats["skipped_gaps"] += 1
                continue
            if len(set(bases)) > 1:
                ref_base = bases[0]
                lineage_base = {}
                for b, rec in zip(bases, recs):
                    lin = _ancestry.get(rec.id.rsplit(".", 1)[0], "UNKNOWN")
                    lineage_base.setdefault(lin, b)
                alt_base = next(iter(set(bases) - {ref_base}))
                afs = {l: 0.0 for l in _lineages}
                for lin, b in lineage_base.items():
                    if lin in afs and b == alt_base:
                        afs[lin] = 1.0
                for rec in recs:
                    pos = col_to_pos[rec.id][ci]
                    entry = {"bases": dict(lineage_base), "ref": ref_base, "alt": alt_base}
                    entry.update(afs)
                    local_db[(rec.id, pos)] = entry
                stats["snps_found"] += 1
        os.remove(aln_path)
    shutil.rmtree(worker_temp, ignore_errors=True)
    return local_db, stats


def main():
    parser = argparse.ArgumentParser(description="MSA SNP extractor (Step 1a)")
    parser.add_argument("--fasta", required=True, help="phased transcriptome FASTA")
    parser.add_argument("--cluster", required=True)
    parser.add_argument("--ancestry", required=True)
    parser.add_argument("--out", default="Perfect_SNP_DB.pkl")
    parser.add_argument("--threads", type=int, default=0)
    parser.add_argument("--lq-out", dest="lq_out", default=None,
                        help="also write the DB as lqcore text (same as export_snp_db.py --mode msa), "
                             "so the C++ engine needs no second pass over the pickle")
    parser.add_argument("--mafft_mode", default="retree1",
                        choices=["auto", "retree1", "fast"])
    args = parser.parse_args()

    n_threads = args.threads if args.threads > 0 else max(1, mp.cpu_count() // 2)
    if not check_mafft():
        print(f"[{get_now()}] [FATAL] MAFFT not found on PATH.")
        sys.exit(1)
    mafft_opts = {"auto": ["--auto", "--quiet"],
                  "retree1": ["--retree", "1", "--quiet"],
                  "fast": ["--retree", "1", "--maxiterate", "0", "--quiet"]}[args.mafft_mode]

    # Ancestry + dynamic lineage names
    ancestry = {}
    with open(args.ancestry) as f:
        for line in f:
            if not line.strip() or line.startswith("GeneID"):
                continue
            parts = line.strip().split()
            if len(parts) >= 2:
                ancestry[parts[0].split(".")[0]] = parts[1].upper()
    lineages = sorted(set(ancestry.values()))
    print(f"[{get_now()}] [INIT] MSA SNP | {n_threads} workers | "
          f"{len(lineages)} lineages: {lineages}")

    # Load FASTA, pick best transcript per gene
    seq_dict, clean_to_full = {}, defaultdict(list)
    for rec in SeqIO.parse(args.fasta, "fasta"):
        seq_dict[rec.id] = str(rec.seq).upper()
        clean_to_full[rec.id.split(".")[0]].append(rec.id)
    best = {cid: max(fids, key=lambda fid: len(seq_dict[fid]))
            for cid, fids in clean_to_full.items()}

    # Parse clusters
    batches, dropped_seq, dropped_anc = [], 0, 0
    with open(args.cluster) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            cid = parts[0]
            alleles = [a.split(".")[0] for a in parts[1].split(",")]
            allele_map = {}
            for a in alleles:
                fid = best.get(a)
                if fid is None:
                    dropped_seq += 1
                elif a not in ancestry:
                    dropped_anc += 1
                else:
                    allele_map[a] = fid
            if len(allele_map) >= 2:
                batches.append({"cid": cid, "allele_map": allele_map})
    print(f"[{get_now()}] [INFO] {len(batches):,} clusters to align")

    chunk = 50
    chunks = [batches[i:i + chunk] for i in range(0, len(batches), chunk)]
    temp_root = tempfile.mkdtemp(prefix="snpdb_")
    pool = mp.Pool(n_threads, initializer=_init_worker,
                   initargs=(seq_dict, ancestry, temp_root, mafft_opts, lineages))
    snp_db = {}
    agg = {"clusters": 0, "snps_found": 0, "skipped_gaps": 0, "skipped_identical": 0}
    n_done = 0
    try:
        for local_db, st in pool.imap_unordered(_process_cluster_batch, chunks):
            snp_db.update(local_db)
            for kk in agg:
                agg[kk] += st.get(kk, 0)
            n_done += 1
            if n_done % 50 == 0:
                print(f"  ... {n_done}/{len(chunks)} chunks ({agg['snps_found']:,} SNPs)")
    finally:
        pool.close()
        pool.join()
        shutil.rmtree(temp_root, ignore_errors=True)

    print(f"\n[SUMMARY] clusters={agg['clusters']:,} identical={agg['skipped_identical']:,} "
          f"SNPs={agg['snps_found']:,} DB_entries={len(snp_db):,}")
    with open(args.out, "wb") as f:
        pickle.dump(snp_db, f, protocol=pickle.HIGHEST_PROTOCOL)
    print(f"[{get_now()}] [SUCCESS] {args.out} ({os.path.getsize(args.out)/1e6:.1f} MB)")
    if args.lq_out:
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from export_snp_db import write_lq
        n, nt, _ = write_lq(snp_db, "msa", args.lq_out)
        print(f"[{get_now()}] [SUCCESS] {args.lq_out} ({n:,} SNPs on {nt:,} transcripts, lqcore text)")


if __name__ == "__main__":
    main()
