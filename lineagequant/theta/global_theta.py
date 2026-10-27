# -*- coding: utf-8 -*-
"""
global_theta.py -- per-copy expression prior (Theta) engine (Step 01 / 01_hybrid).

Unified Theta engine covering both modes:
  - Hybrid (long-read available): clusters with >= min_lr_depth long reads use
    the physical long-read counts as Theta (strong physical anchor); clusters
    below threshold fall back to k-mer statistical inference.
  - K-mer only (no long-read): all clusters use paired-end Theta-k-mer harvesting
    with cluster-adaptive shrinkage (the original Step 01 behaviour).

Theta k-mers (system 1) = cluster-internal UNIQUE k-mers (owned by exactly one
copy). Distinct from lineage k-mers (system 2, KMC subtraction). Output column 5
is Theta, for compatibility with core/em_dynamic.py.
"""
import os
import sys
import gzip
import math
import pickle
import argparse
import multiprocessing as mp
from collections import defaultdict, Counter, deque

try:
    from Bio import SeqIO
except ImportError:
    print("[FATAL] biopython not found. pip install biopython")
    sys.exit(1)

global_master_dict = None
global_k = None
trans_table = str.maketrans("ATCGN", "TAGCN")


def get_now():
    import datetime
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def init_worker(shared_dict, kmer_len):
    global global_master_dict, global_k
    global_master_dict, global_k = shared_dict, kmer_len


def process_read_pair_chunk(pairs):
    local_hits = defaultdict(int)
    n_pairs = len(pairs) // 2
    for p in range(n_pairs):
        r1 = pairs[p * 2].upper()
        r2 = pairs[p * 2 + 1].upper()
        r1_rc = r1.translate(trans_table)[::-1]
        r2_rc = r2.translate(trans_table)[::-1]
        hits_r1, hits_r2 = set(), set()
        for seq in (r1, r1_rc):
            for j in range(len(seq) - global_k + 1):
                km = seq[j:j + global_k]
                if km in global_master_dict:
                    aid = global_master_dict[km][1]
                    local_hits[aid] += 1
                    hits_r1.add(aid)
        for seq in (r2, r2_rc):
            for j in range(len(seq) - global_k + 1):
                km = seq[j:j + global_k]
                if km in global_master_dict:
                    aid = global_master_dict[km][1]
                    local_hits[aid] += 1
                    hits_r2.add(aid)
        for aid in (hits_r1 & hits_r2):
            local_hits[aid] += 1
    return local_hits


def process_read_chunk(reads):
    local_hits = defaultdict(int)
    for read_seq in reads:
        read_seq = read_seq.upper()
        rc = read_seq.translate(trans_table)[::-1]
        for j in range(len(read_seq) - global_k + 1):
            kf = read_seq[j:j + global_k]
            if kf in global_master_dict:
                local_hits[global_master_dict[kf][1]] += 1
            kr = rc[j:j + global_k]
            if kr in global_master_dict:
                local_hits[global_master_dict[kr][1]] += 1
    return local_hits


def fastq_pair_generator(r1_path, r2_path):
    o1 = gzip.open if r1_path.endswith(".gz") else open
    o2 = gzip.open if r2_path.endswith(".gz") else open
    with o1(r1_path, "rt") as f1, o2(r2_path, "rt") as f2:
        while True:
            try:
                r1l = [next(f1) for _ in range(4)]
                r2l = [next(f2) for _ in range(4)]
            except StopIteration:
                break
            s1, s2 = r1l[1].strip(), r2l[1].strip()
            if s1 and s2:
                yield s1, s2


def _drain(pending, final, keep):
    """collect finished chunks until at most `keep` are outstanding (bounded memory)"""
    while len(pending) > keep:
        i, r = pending.popleft()
        try:
            final.update(r.get(timeout=3600))
        except Exception as e:
            print(f"[{get_now()}] [WARN] chunk {i} failed: {e}")


def scan_fastq_paired(r1, r2, pool, chunk_size=25000, max_pending=None):
    # 1.0.0: the reader used to queue the whole FASTQ in memory before the workers caught up
    # (killed by the OOM killer on large samples); at most 2 x workers chunks are now in flight.
    final = Counter()
    max_pending = max_pending or 2 * pool._processes
    chunk, pending, n = [], deque(), 0
    for s1, s2 in fastq_pair_generator(r1, r2):
        chunk.append(s1)
        chunk.append(s2)
        if len(chunk) >= chunk_size * 2:
            pending.append((n, pool.apply_async(process_read_pair_chunk, (chunk,)))); n += 1
            chunk = []
            _drain(pending, final, max_pending)
    if chunk:
        pending.append((n, pool.apply_async(process_read_pair_chunk, (chunk,))))
    _drain(pending, final, 0)
    return final


def estimate_cluster_shrinkage(probe_counts, alleles, base_shrinkage=10.0):
    probes = [probe_counts.get(a, 0) for a in alleles]
    median_p = float(sorted(probes)[len(probes) // 2]) if probes else 0
    if median_p <= 0:
        return base_shrinkage
    adaptive = base_shrinkage * math.sqrt(median_p / base_shrinkage)
    return max(3.0, min(30.0, adaptive))


def main():
    parser = argparse.ArgumentParser(description="Theta engine (Step 01 / hybrid)")
    parser.add_argument("-c", "--clusters", required=True)
    parser.add_argument("-f", "--fasta", required=True, help="transcriptome FASTA")
    parser.add_argument("-r1", "--read1", required=True)
    parser.add_argument("-r2", "--read2", required=True)
    parser.add_argument("--longread_prior", default=None,
                        help="LongRead_Counts.tsv from Step 00b. If provided, enables "
                             "hybrid physical-anchor mode; otherwise pure k-mer mode.")
    parser.add_argument("-k", "--kmer", type=int, default=21)
    parser.add_argument("-o", "--outdir", default="out")
    parser.add_argument("--index_dir", default=None,
                        help="k-mer 索引 (transcript_index_k*.pkl) 的存放目录。\n不给则用 --outdir。索引只依赖转录组和 k, 与样本无关, 多样本应共用。")
    parser.add_argument("-t", "--threads", type=int, default=16)
    parser.add_argument("--min_lr_depth", type=int, default=5)
    parser.add_argument("--shrinkage", type=float, default=10.0)
    parser.add_argument("--adaptive_shrinkage", action="store_true", default=False)
    parser.add_argument("--no_paired_bonus", action="store_true", default=False)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    k = args.kmer
    idx_dir = args.index_dir or args.outdir
    os.makedirs(idx_dir, exist_ok=True)
    index_pkl = os.path.join(idx_dir, f"transcript_index_k{k}.pkl")
    out_tsv = os.path.join(args.outdir, "Global_Theta_Matrix.tsv")

    # Long-read priors (optional)
    lr_counts = defaultdict(float)
    use_hybrid = bool(args.longread_prior and os.path.exists(args.longread_prior))
    if use_hybrid:
        print(f"[{get_now()}] [INFO] Loading long-read physical priors (hybrid mode)...")
        with open(args.longread_prior) as f:
            f.readline()
            for line in f:
                if line.strip():
                    aid, c = line.strip().split("\t")
                    lr_counts[aid] = float(c)
    else:
        print(f"[{get_now()}] [INFO] No long-read prior; pure k-mer Theta mode.")

    # Parse clusters
    print(f"[{get_now()}] [STEP 1] Parsing cluster map...")
    cluster_dict, gene_to_cid, all_targets = {}, {}, set()
    with open(args.clusters) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            cid, alleles = parts[0], parts[1].split(",")
            cluster_dict[cid] = []
            for a in alleles:
                aid = a.strip().split(".")[0]
                cluster_dict[cid].append(aid)
                gene_to_cid[aid] = cid
                all_targets.add(aid)

    # Harvest unique k-mers (Theta k-mers, system 1)
    if os.path.exists(index_pkl):
        print(f"[{get_now()}] [INFO] Loading cached index {index_pkl}")
        with open(index_pkl, "rb") as f:
            master_dict, probe_counts = pickle.load(f)
    else:
        print(f"[{get_now()}] [STEP 2] Harvesting cluster-unique k-mers...")
        kmer_to_alleles = defaultdict(set)
        for rec in SeqIO.parse(args.fasta, "fasta"):
            gid = rec.id.split(".")[0]
            if gid not in all_targets:
                continue
            s = str(rec.seq).upper()
            for i in range(len(s) - k + 1):
                km = s[i:i + k]
                if "N" not in km:
                    kmer_to_alleles[km].add(gid)
        master_dict, probe_counts = {}, defaultdict(int)
        for km, ids in kmer_to_alleles.items():
            if len(ids) == 1:
                gid = next(iter(ids))
                master_dict[km] = (gene_to_cid[gid], gid)
                probe_counts[gid] += 1
        with open(index_pkl, "wb") as f:
            pickle.dump((master_dict, probe_counts), f)

    # Scan FASTQ
    pool = mp.Pool(args.threads, initializer=init_worker, initargs=(master_dict, k))
    if args.no_paired_bonus:
        print(f"[{get_now()}] [STEP 3] FASTQ scan (single-end)...")
        final_hits = Counter()
        def scan_single(fq):
            chunk, pending, n = [], deque(), 0
            of = gzip.open if fq.endswith(".gz") else open
            with of(fq, "rt") as f:
                for i, line in enumerate(f):
                    if i % 4 == 1:
                        chunk.append(line.strip())
                        if len(chunk) >= 50000:
                            pending.append((n, pool.apply_async(process_read_chunk, (chunk,)))); n += 1
                            chunk = []
                            _drain(pending, final_hits, 2 * args.threads)
                if chunk:
                    pending.append((n, pool.apply_async(process_read_chunk, (chunk,))))
            _drain(pending, final_hits, 0)
        scan_single(args.read1)
        scan_single(args.read2)
    else:
        print(f"[{get_now()}] [STEP 3] FASTQ scan (paired-end + bonus)...")
        final_hits = scan_fastq_paired(args.read1, args.read2, pool)
    pool.close()
    pool.join()

    # Compute Theta
    print(f"[{get_now()}] [STEP 4] Computing Theta "
          f"({'hybrid' if use_hybrid else 'k-mer'} mode)...")
    with open(out_tsv, "w") as out_f:
        out_f.write("Cluster_ID\tAllele_ID\tUnique_Probes\tRaw_Hits\tLR_Count\t"
                    "Theta\tPrior_Status\tShrinkage_Used\n")
        for cid, alleles in cluster_dict.items():
            cluster_lr = sum(lr_counts.get(a, 0) for a in alleles) if use_hybrid else 0
            if use_hybrid and cluster_lr >= args.min_lr_depth:
                status = "LongRead_Anchored"
                eps = 0.01
                total = cluster_lr + len(alleles) * eps
                for aid in alleles:
                    lr = lr_counts.get(aid, 0)
                    p = probe_counts.get(aid, 0)
                    h = final_hits.get(aid, 0)
                    theta = (lr + eps) / total
                    out_f.write(f"{cid}\t{aid}\t{p}\t{h}\t{lr}\t{theta:.6f}\t{status}\tNA\n")
            else:
                status = "Kmer_Inferred"
                allele_data = []
                for aid in alleles:
                    p = probe_counts.get(aid, 0)
                    h = final_hits.get(aid, 0)
                    raw = (h / p) if p > 0 else 0.0
                    allele_data.append((aid, raw, p, h))
                total_p = sum(p for _, _, p, _ in allele_data if p > 0)
                cmean = (sum(r * p for _, r, p, _ in allele_data if p > 0) / total_p
                         if total_p > 0 else 0.0)
                shrink = (estimate_cluster_shrinkage(probe_counts, alleles, args.shrinkage)
                          if args.adaptive_shrinkage else args.shrinkage)
                shrunk = {}
                for aid, raw, p, h in allele_data:
                    if p > 0:
                        w = p / (p + shrink)
                        shrunk[aid] = max(w * raw + (1 - w) * cmean, 1e-8)
                    else:
                        shrunk[aid] = max(cmean if cmean > 0 else 1e-6, 1e-8)
                tot = sum(shrunk.values())
                for aid, raw, p, h in allele_data:
                    theta = (shrunk[aid] / tot) if tot > 0 else 0.0
                    lr = lr_counts.get(aid, 0)
                    out_f.write(f"{cid}\t{aid}\t{p}\t{h}\t{lr}\t{theta:.6f}\t{status}\t{shrink:.2f}\n")
    print(f"[{get_now()}] [SUCCESS] {out_tsv}")


if __name__ == "__main__":
    main()
