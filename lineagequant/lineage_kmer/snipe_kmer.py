#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
snipe_kmer.py (v2) -- per-read lineage probability via lineage-specific k-mers (Step 02)

相对原版的改动 (每一条都有实测支撑, cotton 模拟):

  [S1] ★ 原版: probs[n] = hits[n] / total_hits, 且
              status = "Specific" if len(hit_counts) == 1
       -> 一条 read 只要命中 1 个 A marker、0 个 D marker, 就得到 P(A)=1.0000 且判 Specific。
          命中 1 次和命中 20 次输出完全一样, 没有任何置信度。
       实测后果: cluster 内多映射的 read 里 47.93% 被标成 Specific, 而这批 read 的
                 AUC = 0.5218 (抛硬币)。它们撞上的是 A2/D5 祖先分歧的残留 k-mer,
                 那个 k-mer 在多倍体的 At/Dt 转录本上两边都存在。
       修复: Beta 平滑  P(l) = (h_l + alpha) / (sum_h + L*alpha)
             status 由证据量决定:
               sum_h == 0                         -> Ambiguous_None
               sum_h <  --min-hits                -> Weak_Evidence   (P 已平滑, 接近均匀)
               max P  >= --specific-thresh        -> Specific
               否则                                -> Conflict_Weighted
             并输出 N_Hits 列, 让下游 (EM / merge_poe) 能按证据量给权重。

  [S2] 原版 pe_cache[qname] 是 set(query_sequence)。同一 tid 上若存在反向链的次要记录,
       其 RC 序列会作为第 3、4 条进入集合 -> k-mer 命中数被重复计数。
       修复: 按 mate (is_read1) 去重, 每个 fragment 最多两条序列。

  [S3] 原版正反链各扫一遍。若 marker 集是 KMC 的 canonical k-mer, 一次真实命中会被
       算两次 / 或系统性漏掉。
       修复: 启动时自动检测 marker 集是否 canonical, 是则只查 canonical form。

  [S4] 新增 Best_AS 列 (该 read 在该基因上的最优 alignment score)。
       让 EM 可以把 0/1 兼容矩阵换成 exp((AS - AS_best)/T) 的软兼容矩阵。

输出列:
    Read_Name  Target_Gene  P(<lin>)...  Lineage_Status  N_Hits  Best_AS
(尾部两列。下游 em_dynamic 必须用正则抓 P(...) 列, 不能再用 header[2:-1]。)

--legacy 可完全复现原版行为 (用于对照)。
"""
import sys
import argparse
import multiprocessing as mp
from collections import defaultdict

import pysam

_dicts = _names = _bam = _k = None
_alpha = _min_hits = _spec_th = _canon = _legacy = None
_TT = str.maketrans("ATCGN", "TAGCN")


def _rc(s):
    return s.translate(_TT)[::-1]


def load_kmer_set(fa, k):
    ks = set()
    with open(fa) as f:
        for line in f:
            if line.startswith(">"):
                continue
            s = line.strip().upper()
            if len(s) >= k:
                ks.add(s[:k])
    return ks


def detect_canonical(kset, k, n_probe=2000):
    """marker 集是不是 canonical(每个 k-mer 只存 min(kmer, rc) 一种形式)?"""
    both = 0
    seen = 0
    for km in kset:
        if seen >= n_probe:
            break
        seen += 1
        if _rc(km) in kset:
            both += 1
    if seen == 0:
        return False
    frac = both / seen
    return frac < 0.05          # 几乎没有互补对同时存在 -> canonical


def init_worker(d, n, b, k, alpha, min_hits, spec_th, canon, legacy):
    global _dicts, _names, _bam, _k, _alpha, _min_hits, _spec_th, _canon, _legacy
    _dicts, _names, _bam, _k = d, n, b, k
    _alpha, _min_hits, _spec_th, _canon, _legacy = alpha, min_hits, spec_th, canon, legacy


def count_hits(seq):
    hits = {n: 0 for n in _names}
    f = seq.upper()
    L = len(f)
    if L < _k:
        return hits
    if _canon and not _legacy:
        for i in range(L - _k + 1):
            km = f[i:i + _k]
            if "N" in km:
                continue
            c = min(km, _rc(km))
            for n in _names:
                if c in _dicts[n]:
                    hits[n] += 1
    else:
        r = _rc(f)
        for i in range(L - _k + 1):
            kf, kr = f[i:i + _k], r[i:i + _k]
            for n in _names:
                if kf in _dicts[n]:
                    hits[n] += 1
                if kr in _dicts[n]:
                    hits[n] += 1
    return hits


def calc_P(seq):
    hits = count_hits(seq)
    total = sum(hits.values())
    L = len(_names)

    if _legacy:                       # 原版行为, 仅用于对照
        nz = {n: c for n, c in hits.items() if c > 0}
        if total > 0:
            probs = {n: hits[n] / total for n in _names}
            status = "Specific" if len(nz) == 1 else "Conflict_Weighted"
        else:
            probs = {n: 1.0 / L for n in _names}
            status = "Ambiguous_None"
        return probs, status, total

    if total == 0:
        return {n: 1.0 / L for n in _names}, "Ambiguous_None", 0

    # Beta 平滑: 命中数少 -> 自动靠近均匀先验
    probs = {n: (hits[n] + _alpha) / (total + L * _alpha) for n in _names}
    top = max(probs.values())
    if total < _min_hits:
        status = "Weak_Evidence"
    elif top >= _spec_th:
        status = "Specific"
    else:
        status = "Conflict_Weighted"
    return probs, status, total


def process_chunk(id_chunk):
    bam = pysam.AlignmentFile(_bam, "rb")
    out = []
    for tid in id_chunk:
        gene_id = tid.split(".")[0]
        seqs = defaultdict(dict)      # qname -> {is_read1: seq}
        best_as = {}                  # qname -> best AS
        try:
            for read in bam.fetch(tid):
                if read.is_unmapped or not read.query_sequence:
                    continue
                q = read.query_name
                # 按 mate 去重, 避免反向链的次要记录把序列重复塞进来
                seqs[q].setdefault(read.is_read1, read.query_sequence)
                try:
                    a = read.get_tag("AS")
                except KeyError:
                    a = 0
                if q not in best_as or a > best_as[q]:
                    best_as[q] = a
        except ValueError:
            continue
        for q, mates in seqs.items():
            probs, status, nh = calc_P("N".join(mates.values()))
            p_str = "\t".join(f"{probs[n]:.4f}" for n in _names)
            out.append(f"{q}\t{gene_id}\t{p_str}\t{status}\t{nh}\t{best_as.get(q, 0)}\n")
    bam.close()
    return out


def main():
    ap = argparse.ArgumentParser(description="K-mer lineage sniper (Step 02, v2)")
    ap.add_argument("-i", "--ids", required=True)
    ap.add_argument("-b", "--bam", required=True)
    ap.add_argument("--kmers", nargs="+", required=True, metavar="NAME:PATH")
    ap.add_argument("-k", "--kmer_len", type=int, required=True)
    ap.add_argument("-o", "--out", default="P_values_kmer.tsv")
    ap.add_argument("-t", "--threads", type=int, default=16)
    ap.add_argument("--alpha", type=float, default=1.0,
                    help="Beta 平滑强度; 越大越保守。0 = 不平滑")
    ap.add_argument("--min-hits", type=int, default=3,
                    help="低于这个命中数一律判 Weak_Evidence")
    ap.add_argument("--specific-thresh", type=float, default=0.90,
                    help="平滑后最大 P 超过它才判 Specific")
    ap.add_argument("--legacy", action="store_true",
                    help="完全复现原版 (hits/total, 1 个命中就 Specific)")
    args = ap.parse_args()

    dicts, names = {}, []
    for item in args.kmers:
        name, path = item.split(":", 1)
        names.append(name)
        dicts[name] = load_kmer_set(path, args.kmer_len)
        print(f"[INFO] marker {name}: {len(dicts[name]):,} 个 {args.kmer_len}-mer")

    canon = all(detect_canonical(dicts[n], args.kmer_len) for n in names)
    print(f"[INFO] marker 集看起来是 {'canonical' if canon else '双链全存'}"
          f" -> {'只查 canonical form' if canon and not args.legacy else '正反链各扫一遍'}")
    if args.legacy:
        print("[WARN] --legacy: 复现原版行为 (P=hits/total, 单一血缘命中即 Specific)")
    else:
        print(f"[INFO] Beta alpha={args.alpha}  min_hits={args.min_hits}  "
              f"specific_thresh={args.specific_thresh}")

    with open(args.ids) as f:
        valid = {l.strip().split(".")[0] for l in f if l.strip()}
    bi = pysam.AlignmentFile(args.bam, "rb")
    tids = [t for t in bi.references if t.split(".")[0] in valid]
    bi.close()
    print(f"[INFO] {len(tids):,} target transcripts in BAM index.")

    chunks = [tids[i:i + 100] for i in range(0, len(tids), 100)]
    pool = mp.Pool(args.threads, initializer=init_worker,
                   initargs=(dicts, names, args.bam, args.kmer_len, args.alpha,
                             args.min_hits, args.specific_thresh, canon, args.legacy))
    total = 0
    st_ct = defaultdict(int)
    hit_hist = defaultdict(int)
    with open(args.out, "w") as fo:
        fo.write("Read_Name\tTarget_Gene\t" + "\t".join(f"P({n})" for n in names) +
                 "\tLineage_Status\tN_Hits\tBest_AS\n")
        for res in pool.imap_unordered(process_chunk, chunks):
            for line in res:
                fo.write(line)
                total += 1
                p = line.rstrip("\n").split("\t")
                st_ct[p[-3]] += 1
                h = int(p[-2])
                hit_hist[min(h, 10)] += 1
    pool.close()
    pool.join()

    print(f"\n[SUMMARY] {total:,} (fragment, gene) 行")
    for s, c in sorted(st_ct.items(), key=lambda x: -x[1]):
        print(f"  {s:<20} {c:>12,}  ({c/max(total,1)*100:5.2f}%)")
    print("\n  N_Hits 分布 (>=10 合并):")
    for h in sorted(hit_hist):
        lab = f"{h}" if h < 10 else ">=10"
        print(f"    {lab:>4}: {hit_hist[h]:>12,}  ({hit_hist[h]/max(total,1)*100:5.2f}%)")
    print("\n  ↑ 如果 N_Hits==1 的比例很高, 那就是原版把它们全判成 Specific 且 P=1.0 的地方")
    print(f"[SUCCESS] {args.out}")


if __name__ == "__main__":
    main()
