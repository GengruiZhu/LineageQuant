#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_markers.py (v2) -- Lineage-specific k-mer marker construction (KMC).

原版的问题 (cotton 实测):
    marker 只做了 "祖先 L 的转录本 − 其他祖先的转录本" 的减法。
    可是我们要判的是 read 来自**多倍体内部**哪个亚基因组拷贝, 不是来自哪个祖先。
    结果:
        mkA ∩ At       = 9,099,738   -> 只有 50.1% 的 A marker 出现在 At 转录本里
        mkD ∩ Dt       = 8,891,405   -> 只有 27.1%
        mkA ∩ At ∩ Dt  = 1,176,002   -> 声称"A2 特异", 却在 At 和 Dt 上同时存在
        mkD ∩ Dt ∩ At  = 1,184,140

    那 118 万个 k-mer 就是全部问题所在: 一段在 At/Dt 上完全相同的序列(read 会多映射),
    因为 D5 祖先基因组里那个基因丢失/分化了, 于是整段被标成"A2 特异"。
    多映射 read 在这段上密集命中十几个假 marker -> P=1.0 + status=Specific,
    而它们的真实判别力 AUC = 0.52 (抛硬币)。

v2 的改动:
    新增 --polyploid NAME:FASTA (多倍体自己的亚基因组转录本, 每条血缘一个)
    marker 定义变成:
        L_spec = (anc_L − 其他 anc) ∩ poly_L − 其他 poly
    即一个 marker 必须同时:
        1. 在祖先 L 的转录本里有、其他祖先里没有   (原有条件)
        2. 出现在多倍体的 L 亚基因组转录本里        (新增)
        3. 不出现在多倍体的其他亚基因组转录本里    (新增)

    --mode poly  则完全不用祖先: L_spec = poly_L − 其他 poly
                 (最干净, 覆盖率最高; 但无法检测同源交换/基因转换)

cotton 实测效果 (mode=ancestor_x_poly):
    A: 18,162,395 -> 7,923,736 (43.63%)     D: 32,847,193 -> 7,707,265 (23.46%)
    A:D 比例 1:1.80 -> 1:0.97  (那个不对称本来就是祖先基因组重复序列含量的差异)
    多映射 read: Specific 47.93% -> 3.36%,  Ambiguous_None 41.16% -> 95.93%
    残留的 3.36% Specific 的 AUC = 0.96 (真信息)
"""

import os
import sys
import shutil
import argparse
import subprocess
from collections import OrderedDict

try:
    from . import filter_low_complexity as lcmod
    from . import txt2fasta as t2fmod
except ImportError:
    _here = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, _here)
    import filter_low_complexity as lcmod
    import txt2fasta as t2fmod


def _run(cmd, quiet=True):
    out = subprocess.DEVNULL if quiet else None
    subprocess.run(cmd, check=True, stdout=out, stderr=out)


def _check_tools():
    for t in ("kmc", "kmc_tools"):
        if shutil.which(t) is None:
            raise RuntimeError(f"Required tool '{t}' not found on PATH.")


def _db_size(db):
    """k-mer 计数, 用于日志"""
    hist = db + ".hist"
    try:
        _run(["kmc_tools", "transform", db, "histogram", hist])
        with open(hist) as f:
            return sum(int(l.split()[1]) for l in f if l.strip())
    except Exception:
        return -1


def build_kmc_db(fa, db_prefix, k, threads, max_count=1_000_000):
    if os.path.exists(db_prefix + ".kmc_pre"):
        return db_prefix
    tmp = db_prefix + "_kmctmp"
    os.makedirs(tmp, exist_ok=True)
    try:
        # 注意: 不要用 -cs1, 会关掉 counter, kmc_tools 直接拒绝后续运算
        _run(["kmc", f"-k{k}", "-ci1", f"-cs{max_count}", f"-t{threads}",
              "-fm", fa, db_prefix, tmp])
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    return db_prefix


def complex_expr(name, expr_lines, out_db, outdir, threads):
    cfg = os.path.join(outdir, f"{name}.cfg")
    with open(cfg, "w") as f:
        f.write("\n".join(expr_lines) + "\n")
    _run(["kmc_tools", f"-t{threads}", "complex", cfg])
    return out_db


def build_markers(anc_fastas, outdir, k=31, threads=8, at_gc_frac=0.8,
                  poly_fastas=None, mode="ancestor", keep_intermediate=False):
    """
    anc_fastas  : {lineage: 祖先转录本 fasta}
    poly_fastas : {lineage: 多倍体对应亚基因组转录本 fasta}  (mode != 'ancestor' 时必需)
    mode        : ancestor         原版行为
                  ancestor_x_poly  (anc_L − 其他anc) ∩ poly_L − 其他poly   <- 推荐
                  poly             poly_L − 其他poly                        <- 最干净
    """
    _check_tools()
    os.makedirs(outdir, exist_ok=True)

    if mode == "poly":
        if not poly_fastas or len(poly_fastas) < 2:
            raise ValueError("mode=poly 需要 >=2 个 --polyploid")
        names = list(poly_fastas.keys())
    else:
        if len(anc_fastas) < 2:
            raise ValueError("需要 >=2 条血缘")
        names = list(anc_fastas.keys())
        if mode == "ancestor_x_poly":
            missing = [n for n in names if n not in (poly_fastas or {})]
            if missing:
                raise ValueError(f"mode=ancestor_x_poly 缺 --polyploid: {missing}")

    log = sys.stderr.write

    anc_db, poly_db = {}, {}
    if mode != "poly":
        for n in names:
            log(f"[build_markers] KMC count ancestor: {n}\n")
            anc_db[n] = build_kmc_db(anc_fastas[n], os.path.join(outdir, f"anc_{n}_k{k}"),
                                     k, threads)
            log(f"    {n}: {_db_size(anc_db[n]):,} k-mer\n")
    if mode != "ancestor":
        for n in names:
            log(f"[build_markers] KMC count polyploid subgenome: {n}\n")
            poly_db[n] = build_kmc_db(poly_fastas[n], os.path.join(outdir, f"poly_{n}_k{k}"),
                                      k, threads)
            log(f"    {n}: {_db_size(poly_db[n]):,} k-mer\n")

    txts = {}
    for n in names:
        others = [o for o in names if o != n]
        spec = os.path.join(outdir, f"{n}_spec")
        txt = os.path.join(outdir, f"{n}.txt")
        lines = ["INPUT:"]
        terms = []

        if mode == "poly":
            lines.append(f"s = {poly_db[n]} -ci1")
            for i, o in enumerate(others, 1):
                lines.append(f"q{i} = {poly_db[o]} -ci1")
            expr = "s" + "".join(f" - q{i}" for i in range(1, len(others) + 1))
        else:
            lines.append(f"s = {anc_db[n]} -ci1")
            for i, o in enumerate(others, 1):
                lines.append(f"o{i} = {anc_db[o]} -ci1")
            expr = "s" + "".join(f" - o{i}" for i in range(1, len(others) + 1))
            if mode == "ancestor_x_poly":
                lines.append(f"p = {poly_db[n]} -ci1")
                for i, o in enumerate(others, 1):
                    lines.append(f"q{i} = {poly_db[o]} -ci1")
                expr = f"({expr}) * p" + "".join(f" - q{i}" for i in range(1, len(others) + 1))

        lines += ["OUTPUT:", f"{spec} = {expr}"]
        log(f"[build_markers] {n}: {expr}\n")
        complex_expr(n, lines, spec, outdir, threads)
        log(f"    -> {_db_size(spec):,} k-mer\n")
        _run(["kmc_tools", "transform", spec, "dump", txt])
        txts[n] = txt

    merged = os.path.join(outdir, "merged.txt")
    merged_lc = os.path.join(outdir, "merged_lc.txt")
    with open(merged, "w") as out:
        out.write("Kmer\tFinalScore\tSpecies\n")
        for n in names:
            with open(txts[n]) as fin:
                for line in fin:
                    kmer = line.split("\t", 1)[0].split()[0]
                    if kmer:
                        out.write(f"{kmer}\t1\t{n}\n")

    lcmod.filter_file(merged, merged_lc, k=k, at_gc_frac=at_gc_frac)
    marker_paths = t2fmod.write_fastas(merged_lc, outdir)

    log("[build_markers] final marker counts:\n")
    for n in names:
        mf = os.path.join(outdir, f"{n}.fasta")
        c = sum(1 for l in open(mf) if l.startswith(">")) if os.path.exists(mf) else 0
        log(f"    {n}: {c:,}\n")

    if not keep_intermediate:
        for p in [merged, merged_lc] + list(txts.values()):
            try:
                os.remove(p)
            except OSError:
                pass
    return marker_paths


def _group(items, flag):
    g = OrderedDict()
    for it in items or []:
        if ":" not in it:
            raise SystemExit(f"{flag} 必须是 NAME:PATH, 收到: {it}")
        n, p = it.split(":", 1)
        g.setdefault(n, []).append(p)
    return g


def _merge_multi(grouped, outdir, tag):
    res = {}
    for n, paths in grouped.items():
        if len(paths) == 1:
            res[n] = paths[0]
        else:
            os.makedirs(outdir, exist_ok=True)
            m = os.path.join(outdir, f"_{tag}_{n}_merged.fa")
            with open(m, "w") as out:
                for p in paths:
                    with open(p) as fin:
                        shutil.copyfileobj(fin, out)
            sys.stderr.write(f"[build_markers] {tag} {n}: merged {len(paths)} inputs\n")
            res[n] = m
    return res


def main():
    ap = argparse.ArgumentParser(
        description="Build lineage-specific marker k-mers via KMC (v2).")
    ap.add_argument("--lineage", action="append", metavar="NAME:FASTA",
                    help="祖先转录本 FASTA, NAME:PATH (可重复)")
    ap.add_argument("--polyploid", action="append", metavar="NAME:FASTA",
                    help="多倍体自身该亚基因组的转录本 FASTA, NAME:PATH (可重复)")
    ap.add_argument("--mode", choices=["ancestor", "ancestor_x_poly", "poly"],
                    default="ancestor_x_poly",
                    help="ancestor=原版(不推荐); ancestor_x_poly=推荐; poly=只用多倍体自身")
    ap.add_argument("-o", "--outdir", required=True)
    ap.add_argument("-k", "--kmer", type=int, default=31)
    ap.add_argument("-t", "--threads", type=int, default=8)
    ap.add_argument("--at-gc-frac", type=float, default=0.8)
    ap.add_argument("--keep-intermediate", action="store_true")
    args = ap.parse_args()

    anc = _merge_multi(_group(args.lineage, "--lineage"), args.outdir, "anc")
    poly = _merge_multi(_group(args.polyploid, "--polyploid"), args.outdir, "poly")

    if args.mode != "poly" and len(anc) < 2:
        ap.error("--lineage 至少要两条血缘")
    if args.mode != "ancestor" and len(poly) < 2:
        ap.error(f"mode={args.mode} 需要至少两条 --polyploid")

    paths = build_markers(anc, args.outdir, k=args.kmer, threads=args.threads,
                          at_gc_frac=args.at_gc_frac, poly_fastas=poly,
                          mode=args.mode, keep_intermediate=args.keep_intermediate)
    print("[done] marker FASTAs:")
    for n, p in paths.items():
        print(f"  {n}: {p}")


if __name__ == "__main__":
    main()
