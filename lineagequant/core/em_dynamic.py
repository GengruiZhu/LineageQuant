#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
core/em_dynamic.py (v5) -- EM convergence engine (Step 03)

v5 新增两条, 都是 cotton 实测逼出来的:

  [B10] ★ 有效长度校正 (--eff-len)

    实测 (cotton_rep1, 同一份 HISAT2 BAM):
        HISAT2比对 + salmon的EM     0.8585   slope 0.9674
        HISAT2比对 + LQ的EM         0.8575   slope 0.9657
        salmon 全流程                0.8586   slope 0.9678
    换比对器只值 +0.0001, 换 EM 值 +0.0010。**salmon 领先的是模型, 不是比对。**

    差在哪: salmon 的 E-step 权重是 alpha_i / L_eff_i (每碱基覆盖率),
    LQ 用的是裸 theta_i (read 比例)。簇内两个拷贝长度不同时, 模糊 read 的分配就偏了。

        E: w_ri ∝ M_ri · P_ri · (alpha_i / L_eff_i)
        M: alpha_i = Σ_r w_ri + tau·theta0_i + eps
        L_eff_i = max(L_i − meanFragLen + 1, min_eff)

    --eff-len eff 为默认 (2026-09-29 起, 差距分解: 簇内基因 +1.6~+4.5 pp);
    --eff-len off 时 L_eff ≡ 1, **精确退回 v4**。
    输出的 Final_Theta 仍是 read 比例 (alpha 归一化), 口径不变。

  [B11] 跨 cluster 多映射 read (--cross-cluster)

    实测 cotton: 28,812 条 (0.26%) 的 fragment 比对到 >1 个 cluster 的基因。
    v3/v4 会在每个 cluster 里各存一份, 各贡献一条完整 read -> 重复计数。
    cotton 上影响很小; POJ (簇均 10.04 等位, 31.6% 多映射) 不会小。
      keep(默认)  = 旧行为, 每个 cluster 各记一条 (向后兼容)
      drop        = 丢掉这些 read
      split       = 按候选基因数在各 cluster 间分权重 (一条 read 总权重仍为 1)

---- v4 / v3 的修复 ----
  [B9] theta0 作为 M-step 的 Dirichlet 先验, 权重 tau = 长读深度。
       POJ cluster 0: 313 条长读 vs 1.98 条短读, 旧代码采信后者。
  [B1] W = M · P^w · theta。原版没有比对兼容矩阵 M, 血缘 P 替换了比对证据。
       cotton r 0.7765 -> 0.9445, slope 0.9114 -> 0.9956。
  [B6] 按 Read_Name 去重  [B3/B3b] --p-weight-* / --p-floor  [B5] 正则抓 P(...) 列
  [NEW v3] --bam: M 直接从 BAM 建
"""
import re
import sys
import argparse
import datetime
from collections import defaultdict

import numpy as np


def now():
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def log(msg):
    print(f"[{now()}] {msg}", flush=True)


_SPECIFIC = re.compile(r"^specific", re.I)
_CONFLICT = re.compile(r"conflict|multi", re.I)
_NONE = re.compile(r"none|no_snp|chimeric|ambiguous|weak", re.I)


def status_class(s):
    if _SPECIFIC.match(s):
        return "specific"
    if _NONE.search(s):
        return "none"
    if _CONFLICT.search(s):
        return "conflict"
    return "conflict"


def load_ancestry(path):
    anc = {}
    with open(path, encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith("GeneID"):
                continue
            p = line.split()
            if len(p) >= 2:
                anc[p[0]] = p[1].upper()
                anc.setdefault(p[0].rsplit(".",1)[0], p[1].upper())
    return anc


def load_theta(path):
    cluster_alleles = defaultdict(list)
    initial_thetas, gene_to_cluster, allele_pos = {}, {}, {}
    lr_count = defaultdict(float)
    prior_status = {}
    with open(path, encoding="utf-8") as f:
        f.readline()
        for line in f:
            if not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            cid, aid = p[0], p[1]
            allele_pos[aid] = len(cluster_alleles[cid])
            cluster_alleles[cid].append(aid)
            initial_thetas[aid] = float(p[5])
            gene_to_cluster[aid] = cid
            try:
                lr_count[cid] += float(p[4])
            except (IndexError, ValueError):
                pass
            if len(p) > 6:
                prior_status[cid] = p[6]
    return (cluster_alleles, initial_thetas, gene_to_cluster, allele_pos,
            lr_count, prior_status)


def parse_pval_header(line):
    cols = line.rstrip("\n").split("\t")
    p_cols = []
    for i, c in enumerate(cols):
        m = re.fullmatch(r"P\((.+)\)", c.strip())
        if m:
            p_cols.append((m.group(1).upper(), i))
    status_i = None
    for i, c in enumerate(cols):
        if "status" in c.strip().lower():
            status_i = i
    return 0, 1, status_i, p_cols


def estimate_frag_len(bam_path, n_sample=1000000):
    """从 proper pair 的 TLEN 估平均片段长度。拿不到就返回 None。"""
    import pysam
    bam = pysam.AlignmentFile(bam_path, "rb")
    vals = []
    for r in bam.fetch(until_eof=True):
        if r.is_proper_pair and not r.is_secondary and not r.is_supplementary \
                and r.template_length and 0 < r.template_length < 2000:
            vals.append(r.template_length)
            if len(vals) >= n_sample:
                break
    bam.close()
    if len(vals) < 1000:
        return None
    return float(np.mean(vals))


def load_bam(bam_path, gene_to_cluster, allele_pos, min_mapq=0, primary_only=False):
    """
    返回 reads: {cid: {qname: [mask, pv, wc]}}
         gene_len: {gene: 最长转录本长度}
         n_cid_of: {qname: 该 read 落到几个 cluster}
    """
    import pysam
    reads = defaultdict(dict)
    n_cid = defaultdict(int)
    gene_len = {}
    n_rec = n_used = n_skip = 0
    bam = pysam.AlignmentFile(bam_path, "rb")
    for name, ln in zip(bam.references, bam.lengths):
        g = name.rsplit(".",1)[0]
        if g not in gene_len or ln > gene_len[g]:
            gene_len[g] = ln
    for r in bam.fetch(until_eof=True):
        n_rec += 1
        if r.is_unmapped:
            continue
        if primary_only and (r.is_secondary or r.is_supplementary):
            continue
        if min_mapq and r.mapping_quality < min_mapq:
            continue
        g = r.reference_name.rsplit(".",1)[0]
        cid = gene_to_cluster.get(g)
        if cid is None:
            n_skip += 1
            continue
        n_used += 1
        j = allele_pos[g]
        d = reads[cid]
        rec = d.get(r.query_name)
        if rec is None:
            d[r.query_name] = [1 << j, None, "none"]
            n_cid[r.query_name] += 1
        else:
            rec[0] |= (1 << j)
    bam.close()
    n_reads = sum(len(v) for v in reads.values())
    log(f"[INFO] BAM {n_rec:,} 条记录, {n_used:,} 条落在 cluster 基因上 "
        f"({n_skip:,} 条落在孤儿/非目标基因) -> {n_reads:,} 条 read")
    return reads, gene_len, n_cid


def load_pvals(path, gene_to_cluster, allele_pos, reads, use_dedup, bam_mode, n_cid):
    fh = open(path, encoding="utf-8")
    name_i, gene_i, status_i, p_cols = parse_pval_header(fh.readline())
    lineages = [l for l, _ in p_cols]
    p_idx = [i for _, i in p_cols]
    L = len(lineages)
    log(f"[INFO] lineages={lineages}" if L else "[INFO] P 文件没有 P(...) 列 -> 血缘均匀")

    n_rows = n_skip = n_orphan = 0
    row_id = 0
    stat_ct = defaultdict(int)
    for line in fh:
        if not line.strip():
            continue
        parts = line.rstrip("\n").split("\t")
        cid = gene_to_cluster.get(parts[gene_i])
        if cid is None:
            n_skip += 1
            continue
        n_rows += 1
        j = allele_pos[parts[gene_i]]
        pv, wc = None, "none"
        if L:
            try:
                pv = [float(parts[i]) for i in p_idx]
            except (ValueError, IndexError):
                n_skip += 1
                continue
            st = parts[status_i] if (status_i is not None and len(parts) > status_i) else ""
            wc = status_class(st)
            stat_ct[wc] += 1
        d = reads[cid]
        if bam_mode:
            rec = d.get(parts[name_i])
            if rec is None:
                n_orphan += 1
                continue
            if rec[1] is None:
                rec[1], rec[2] = pv, wc
        else:
            key = parts[name_i] if use_dedup else row_id
            row_id += 1
            rec = d.get(key)
            if rec is None:
                d[key] = [1 << j, pv, wc]
                n_cid[key] += 1
            else:
                rec[0] |= (1 << j)
    fh.close()
    log(f"[INFO] P 文件 {n_rows:,} 行有效, {n_skip:,} 跳过"
        + (f", {n_orphan:,} 行的 read 不在 BAM 里" if n_orphan else ""))
    return L, lineages, stat_ct


# --------------------------------------------------------------- EM
def em_one(M, P, wgt, theta0, tau, eff, eps, max_iter, min_iter, tol, ll_tol):
    """
    E: w_ri ∝ M_ri · P_ri · (alpha_i / eff_i)
    M: alpha_i = Σ_r w_ri·wgt_r + tau·theta0_i + eps
    eff ≡ 1 时精确等于 v4。
    """
    A = theta0.shape[0]
    n = M.shape[0]
    alpha = theta0 * max(float(wgt.sum()), 1.0) + eps
    base = M * P
    prior = tau * theta0
    prev_ll = -np.inf
    counts = np.zeros(A)
    wcol = wgt.reshape(-1, 1)
    n_iter = 0
    for n_iter in range(max_iter):
        rate = alpha / eff
        W = base * rate
        Wsum = W.sum(axis=1, keepdims=True)
        ll = float(np.sum(wcol[Wsum.ravel() > 0] *
                          np.log(np.maximum(Wsum[Wsum > 0], 1e-300)).reshape(-1, 1)))
        ok = (Wsum > 0).ravel()
        counts = np.zeros(A)
        if ok.any():
            counts = ((W[ok] / Wsum[ok]) * wcol[ok]).sum(axis=0)
        new_alpha = counts + prior + eps
        if n_iter >= min_iter:
            f_old = alpha / alpha.sum()
            f_new = new_alpha / new_alpha.sum()
            if float(np.max(np.abs(f_new - f_old))) < tol:
                alpha = new_alpha
                break
            imp = ll - prev_ll
            if 0 <= imp < ll_tol:
                alpha = new_alpha
                break
        alpha = new_alpha
        prev_ll = ll
    theta = alpha / alpha.sum()
    return counts, theta, n_iter + 1, ll


def main():
    ap = argparse.ArgumentParser(description="EM engine (Step 03, v5)")
    ap.add_argument("-t", "--theta", required=True)
    ap.add_argument("-p", "--pvals", default=None)
    ap.add_argument("-a", "--ancestry", required=True)
    ap.add_argument("-o", "--out", default="EM_Final_Expression.tsv")
    ap.add_argument("--bam", default=None)
    ap.add_argument("--min-mapq", type=int, default=0)
    ap.add_argument("--primary-only", action="store_true")
    ap.add_argument("--epsilon", type=float, default=0.01)
    ap.add_argument("--max_iter", type=int, default=200)
    ap.add_argument("--min_iter", type=int, default=5)
    ap.add_argument("--tol", type=float, default=1e-4)
    ap.add_argument("--ll_tol", type=float, default=0.01)
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--p-weight-specific", type=float, default=1.0)
    ap.add_argument("--p-weight-conflict", type=float, default=0.0)
    ap.add_argument("--p-weight-none", type=float, default=0.0)
    ap.add_argument("--p-alpha", type=float, default=0.0)
    ap.add_argument("--p-floor", type=float, default=1e-3)
    ap.add_argument("--theta-prior-lr", type=float, default=1.0)
    ap.add_argument("--theta-prior-kmer", type=float, default=0.0)
    ap.add_argument("--theta-prior-cap", type=float, default=0.0)
    # --- v5 ---
    ap.add_argument("--eff-len", choices=["off", "len", "eff"], default="eff",
                    help="E-step 权重除以长度。off=v4 行为; len=用转录本长度; "
                         "eff=L - meanFragLen + 1 (从 BAM 的 TLEN 估)。需要 --bam")
    ap.add_argument("--frag-len", type=float, default=0.0,
                    help="手动指定平均片段长度 (配 --eff-len eff)。0=从 BAM 估")
    ap.add_argument("--min-eff-len", type=float, default=25.0)
    ap.add_argument("--cross-cluster", choices=["keep", "drop", "split"], default="keep",
                    help="跨 cluster 多映射 read 的处理。keep=旧行为(每簇各记一条, 重复计数); "
                         "drop=丢弃; split=按簇均分权重")
    ap.add_argument("--no-align", action="store_true")
    ap.add_argument("--no-dedup", action="store_true")
    args = ap.parse_args()

    if not args.pvals and not args.bam:
        ap.error("至少要给 -p/--pvals 或 --bam 之一")
    if args.eff_len != "off" and not args.bam:
        log("[WARN] 没有 --bam, 拿不到转录本长度 -> --eff-len off")
        args.eff_len = "off"
    if args.no_align and args.bam:
        ap.error("--no-align 和 --bam 互斥")

    use_align = not args.no_align
    use_dedup = not args.no_dedup
    bam_mode = args.bam is not None
    WC = {"specific": args.p_weight_specific,
          "conflict": args.p_weight_conflict,
          "none": args.p_weight_none}

    log(f"[INFO] EM v5 | M源={'BAM' if bam_mode else ('P文件' if use_align else '无(原版)')} "
        f"dedup={use_dedup} w={WC['specific']}/{WC['conflict']}/{WC['none']} "
        f"floor={args.p_floor} eps={args.epsilon}")
    log(f"[INFO] theta 先验 lr={args.theta_prior_lr} kmer={args.theta_prior_kmer} | "
        f"eff-len={args.eff_len} | cross-cluster={args.cross_cluster}")

    ancestry = load_ancestry(args.ancestry)
    (cluster_alleles, initial_thetas, gene_to_cluster, allele_pos,
     lr_count, prior_status) = load_theta(args.theta)
    log(f"[INFO] {len(cluster_alleles):,} clusters, {len(gene_to_cluster):,} alleles")

    gene_len, n_cid = {}, defaultdict(int)
    if bam_mode:
        reads, gene_len, n_cid = load_bam(args.bam, gene_to_cluster, allele_pos,
                                          args.min_mapq, args.primary_only)
    else:
        reads = defaultdict(dict)

    L, lineages, stat_ct = (0, [], {})
    if args.pvals:
        L, lineages, stat_ct = load_pvals(args.pvals, gene_to_cluster, allele_pos,
                                          reads, use_dedup, bam_mode, n_cid)
    else:
        log("[INFO] 没有 -p -> 血缘均匀, W = M · (alpha/eff)")

    # ---- 有效长度 ----
    mean_fl = 0.0
    if args.eff_len == "eff":
        mean_fl = args.frag_len or (estimate_frag_len(args.bam) or 0.0)
        if mean_fl <= 0:
            log("[WARN] 估不出片段长度, 退回 --eff-len len")
            args.eff_len = "len"
        else:
            log(f"[INFO] 平均片段长度 {mean_fl:.1f} bp (从 proper pair 的 TLEN 估)")

    # ---- 跨 cluster read ----
    n_cross = sum(1 for q, c in n_cid.items() if c > 1)
    if n_cross:
        tot = len(n_cid)
        log(f"[INFO] 跨 cluster 多映射 read: {n_cross:,} / {tot:,} "
            f"({n_cross/tot*100:.2f}%) -> {args.cross_cluster}")

    lin_rank = {l: k for k, l in enumerate(lineages)}
    n_reads = sum(len(v) for v in reads.values())
    log(f"[INFO] 共 {n_reads:,} 条 read 进入 EM")
    if stat_ct:
        tot = sum(stat_ct.values())
        for k in ("specific", "conflict", "none"):
            log(f"[INFO]   status={k:<9} {stat_ct[k]:>12,} ({stat_ct[k]/tot:6.2%}) -> w={WC[k]}")

    def tau_of(cid):
        t = (args.theta_prior_lr * lr_count.get(cid, 0.0)
             if prior_status.get(cid, "") == "LongRead_Anchored" else args.theta_prior_kmer)
        return max(min(t, args.theta_prior_cap) if args.theta_prior_cap > 0 else t, 0.0)

    results, iters, multi = [], [], 0
    for cid, alleles in cluster_alleles.items():
        A = len(alleles)
        theta0 = np.array([initial_thetas[a] for a in alleles], dtype=np.float64)
        eff = np.ones(A)
        if args.eff_len != "off":
            Ls = np.array([float(gene_len.get(a, 1000)) for a in alleles])
            eff = np.maximum(Ls - mean_fl + 1.0 if args.eff_len == "eff" else Ls,
                             args.min_eff_len)
        tau = tau_of(cid)
        lin_of = [lin_rank.get(ancestry.get(a, ancestry.get(a.rsplit(".",1)[0], "?")), -1)
                  for a in alleles] if L else [-1] * A

        d = reads.get(cid)
        if not d:
            th = (theta0 + args.epsilon) / (theta0 + args.epsilon).sum()
            results.append((cid, alleles, theta0, np.zeros(A), th, 0, 0.0))
            continue

        items = list(d.items())
        keep = []
        for q, rec in items:
            c = n_cid.get(q, 1)
            if c > 1 and args.cross_cluster == "drop":
                continue
            keep.append((rec, 1.0 / c if (c > 1 and args.cross_cluster == "split") else 1.0))
        if not keep:
            th = (theta0 + args.epsilon) / (theta0 + args.epsilon).sum()
            results.append((cid, alleles, theta0, np.zeros(A), th, 0, 0.0))
            continue

        n = len(keep)
        Mm = np.ones((n, A))
        P = np.ones((n, A))
        wgt = np.ones(n)
        for r, ((mask, pv, wc), w) in enumerate(keep):
            wgt[r] = w
            if use_align:
                for j in range(A):
                    Mm[r, j] = 1.0 if (mask >> j) & 1 else 0.0
                if bin(mask).count("1") > 1:
                    multi += 1
            if not L or pv is None:
                continue
            pw = WC[wc]
            if pw == 0.0:
                continue
            pvec = np.asarray(pv, dtype=np.float64)
            if args.p_alpha > 0:
                pvec = (pvec + args.p_alpha) / (1.0 + L * args.p_alpha)
            pvec = np.maximum(pvec, args.p_floor)
            if pw != 1.0:
                pvec = pvec ** pw
            for j in range(A):
                P[r, j] = pvec[lin_of[j]] if lin_of[j] >= 0 else (1.0 / L)

        counts, theta, it, ll = em_one(Mm, P, wgt, theta0, tau, eff, args.epsilon,
                                       args.max_iter, args.min_iter, args.tol, args.ll_tol)
        iters.append(it)
        results.append((cid, alleles, theta0, counts, theta, it, ll))

    if iters:
        arr = np.array(iters)
        log(f"[INFO] EM convergence: min={arr.min()} median={np.median(arr):.0f} "
            f"mean={arr.mean():.1f} max={arr.max()} at_max={int((arr >= args.max_iter).sum())}")
    if use_align:
        log(f"[INFO] cluster 内多映射 read: {multi:,} / {n_reads:,} "
            f"({multi/max(n_reads,1)*100:.2f}%)")

    results.sort(key=lambda x: x[0])
    with open(args.out, "w", encoding="utf-8") as fo:
        fo.write("ClusterID\tAlleleID\tAncestry\tInitial_Theta\tEM_Allocated_Reads\t"
                 "Final_Theta\tEM_Iter\tEM_LogLik\n")
        for cid, alleles, t0, cnt, th, it, ll in results:
            for i, a in enumerate(alleles):
                lin = ancestry.get(a, ancestry.get(a.rsplit(".",1)[0], "UNKNOWN"))
                fo.write(f"{cid}\t{a}\t{lin}\t{t0[i]:.6f}\t{cnt[i]:.2f}\t"
                         f"{th[i]:.6f}\t{it}\t{ll:.2f}\n")
    log(f"[SUCCESS] EM complete: {args.out}")


if __name__ == "__main__":
    main()
