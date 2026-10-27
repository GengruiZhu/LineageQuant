# LineageQuant 使用说明（v1.0.0）

多倍体亚基因组 / 等位拷贝表达定量软件。作者：Yi Chen, Gengrui Zhu

本文档讲：安装与检查、输入准备、按数据组合怎么跑、全部参数、输出文件、常见问题。
英文版与设计说明见仓库根目录 `README.md`、`CHANGELOG.md` 与 `docs/DESIGN.md`。

---

## 0. 这个软件做什么

LineageQuant 定量多倍体里"每个亚基因组 / 等位拷贝"的表达量。它是**定量器**，不是 read 分类器：
参考转录组里每个同源 / 等位拷贝都是独立目标，兼容多个拷贝的 read 由 EM 按概率分配，
最终输出拷贝级的 read 数与 TPM，可以再按血统加和得到亚基因组表达。

自动路由：你只声明物种类型（同源 / 异源），再把手头有的输入交给它，软件自己决定走哪条流程。
缺什么输入就跳过对应模块，所有路线最后都汇到 EM → 表达矩阵。

- **同源多倍体（`--ploidy auto`）**：亚基因组无分化，没有"血统"概念，只做拷贝级定量：
  Theta 给每个拷贝初值 → EM 把 read 分到拷贝。
- **异源多倍体（`--ploidy allo`）**：亚基因组有血统分化。先算每条 read 来自各血统的概率，
  再连同比对证据一起进 EM。

两套独立的 k-mer 系统：

| | Theta k-mer（`--kmer`，默认 21） | 血统 k-mer（`--lineage-kmer`，默认 31） |
|---|---|---|
| 作用 | 簇内区分拷贝，给 EM 初值 | 区分血统，给每条 read 的血统概率 |
| 来源 | 转录组里每个拷贝独有的 k-mer | KMC 祖先基因组相减 |
| 用于 | 两种倍性 | 只异源 |

### 1.0.0 的定量方式（与以前版本的主要区别）

- **Step 00 比对**：默认用 bowtie2，参数与 RSEM 默认相同（端到端、不开空位、`--mp 1,1`、
  `--score-min L,0,-0.1`、只保留成对一致的比对、每个片段最多 200 个比对位置），并且参考序列中的
  N 不罚分（`--np 0 --n-ceil L,0,0.5`）。它会把同样好的拷贝全部报出来，也报出稍差的位置。
  以前的 HISAT2 只报其中一部分：基准测试里有 0.01–0.6% 的片段只被报到一个错误的拷贝上，
  这些 read 堆到高表达拷贝旁边的低表达拷贝上，EM 无法纠正。HISAT2 仍可用 `--aligner hisat2`。
- **Step 03 EM（global 引擎）**：全部基因（等位拷贝与单拷贝基因）一起做 EM，每个片段总权重为 1，
  不会重复计数；有效长度校正；候选比对按错配数加权（每多一个错配，似然乘以 λ = (p/3)/(1−p)，
  p 由最佳比对的错配率估计）；序列完全相同的基因合并估计、收敛后平分。

---

## 1. 安装

```bash
cd /path/to/LineageQuant_v1.0.0      # 解压后的目录（里面有 setup.py）
conda env create -f environment.yml  # 新建环境，所有工具一起装好（推荐）
conda activate lineagequant
pip install -e .                     # 注册 lineagequant 命令

lineagequant --version               # 应输出: lineagequant 1.0.0
lineagequant run --help
```

依赖（`environment.yml` 已包含）：

- Python ≥ 3.8：numpy、pandas、pysam、biopython（VCF 路线另需 cyvcf2）
- 必需工具：**bowtie2 ≥ 2.4**、samtools ≥ 1.10
- 按路线需要：hisat2（`--aligner hisat2`）、kmc + kmc_tools + gffread（祖先 k-mer 路线）、
  mafft（MSA SNP 路线）、bcftools / tabix / htslib（VCF 路线）、minimap2（长读长）

注意：`pip install -e .` 只装 Python 包，不装 conda 工具。已有环境里请自己 `conda install`
缺的工具；运行时 `lineagequant run` 会检查本次路线需要的工具，缺哪个就报哪个，例如
`[FATAL] tool 'bowtie2' not on PATH`。

建等位簇的工具（GraphAllele / AlleleFinder / jcvi 等）不在依赖里。自己产出等位簇后，用
`lineagequant normalize-clusters` 标准化（见第 6 节）。

---

## 2. 安装后怎么检查

```bash
# [1] 所有模块能否编译
for f in $(find lineagequant -name "*.py"); do python3 -m py_compile "$f"; done && echo 编译OK

# [2] 命令是否可用
lineagequant --version
lineagequant run --help

# [3] 路由是否正确（--dry-run 只打印路由，不真跑）
mkdir -p /tmp/lq_check && cd /tmp/lq_check
touch mRNA.fa genome.gff clusters.tsv ancestry.tsv orphans.tsv ids.txt ancA.fa ancA.gff ancD.fa ancD.gff
mkdir -p data results
lineagequant run --fasta mRNA.fa --gff genome.gff --clusters clusters.tsv \
  --ancestry ancestry.tsv --orphan orphans.tsv --target-ids ids.txt \
  --ploidy allo --ancestor A:ancA.fa,ancA.gff --ancestor D:ancD.fa,ancD.gff \
  --data-dir data --out-dir results --dry-run
```

应输出类似：

```
[lineagequant 1.0.0] resolved route:
ploidy: allo
  alignment: bowtie2 (RSEM settings, N not penalised, <= 200 alignments/fragment)
  EM: global over all genes | eff-len=eff | mismatch weighting=auto | identical-sequence pooling=on
  lineage routes: ancestor k-mer (2 lineages, 2 progenitor sample(s))
  Theta: k-mer (01)
[dry-run] exiting without execution.
```

---

## 3. 输入文件准备

必需（所有模式）：

| 参数 | 内容 |
|---|---|
| `--fasta` | 转录组（mRNA）FASTA；每个等位拷贝是一条独立记录。可从基因组 + GFF 用 gffread 提取 |
| `--gff` | 参考 GFF3 注释 |
| `--clusters` | 标准化等位簇，每行一簇：`ClusterID<TAB>Allele1,Allele2,...` |
| `--ancestry` | `GeneID<TAB>LineageName`；血统名从第 2 列动态读取，数量与名字不限 |
| `--orphan` | 单拷贝基因：`OrphanID<TAB>GeneID` |
| `--target-ids` | 目标基因 ID，每行一个 |

基因 ID 由转录本 ID 去掉最后一个 `.N` 后缀得到（`Gene.1`、`Gene.2` 是 `Gene` 的同工型）。

血统路线输入（异源，至少给一种）：

- `--ancestor NAME:fasta,gff`：祖先基因组，KMC 相减得血统特异 k-mer；可重复，≥ 2 个。
  同一个 NAME 重复给多次 = 该血统有多个祖先样本，合并后再建 marker。
- `--markers NAME:fasta`：预建血统 marker，跳过 KMC 建库。NAME 要与 `--ancestry` 的血统名一致；
  marker 的 k-mer 长度须与 `--lineage-kmer` 一致。
- `--vcf FILE` + `--sample-lists NAME:file ...`：群体 VCF，每个血统一个样本名文件。
- `--msa-snp`：MSA SNP 路线，直接用主 `--fasta`（主转录组已含所有等位拷贝时用这个）；
  `--phased-fasta FILE` 只在等位拷贝 CDS 在单独文件里时用。两者互斥，也都与 `--vcf` 互斥。

长读长（可选）：`--longread` 只接受已比对到转录组的 sorted BAM，作用是锚定 Theta 和校准 MSA SNP 频率：

```bash
minimap2 -ax map-pb -N 50 -p 0.95 transcriptome.fa reads.fq.gz | samtools sort -o longreads.sorted.bam
samtools index longreads.sorted.bam
```

HiFi / FLNC 用 `-x map-hifi`。`-N 50 -p 0.95` 让多拷贝间打平的比对作为 secondary 写出来，长读计数靠它判断平局。每条分子只计一次：secondary 不计，同一 ZMW 的多条 subread 合计 1。

FASTQ：放在 `--data-dir` 目录里，双端 gzip，自动识别 `.R1.fq.gz`、`_R1.fq.gz`、`_1.fq.gz`、
`_R1.fastq.gz`、`_R1_001.fastq.gz` 等命名。共用的索引（bowtie2 / Theta k-mer）也写在这个目录里。

---

## 4. 怎么跑（按数据组合选，10 条路线）

先 `--dry-run` 看路由，再去掉它正式跑。所有步骤可断点续跑（`.done` 标记），中断后重跑同样命令即可。

### 同源多倍体（auto）

**A1 同源 + 仅短读长（拷贝级定量）**

```bash
lineagequant run \
  --fasta transcriptome.fasta --gff genome.gff3 \
  --clusters clusters.tsv --ancestry ancestry.tsv --orphan orphans.tsv \
  --target-ids target_gene_ids.txt \
  --ploidy auto --data-dir data --out-dir results --threads 32
```

**A2 同源 + 长读长**：在 A1 上加 `--longread longreads.sorted.bam`。

### 异源多倍体（allo）

**B1 异源 + 祖先 k-mer**

```bash
lineagequant run \
  --fasta transcriptome.fasta --gff genome.gff3 \
  --clusters clusters.tsv --ancestry ancestry.tsv --orphan orphans.tsv \
  --target-ids target_gene_ids.txt \
  --ploidy allo \
  --ancestor A:ancestorA.fa,ancestorA.gff \
  --ancestor D:ancestorD.fa,ancestorD.gff \
  --data-dir data --out-dir results --threads 32
```

| 路线 | 在基本命令上的输入 |
|---|---|
| B2 群体 VCF | `--vcf population.vcf.gz --sample-lists A:samples_A.txt D:samples_D.txt` |
| B3 祖先 k-mer + VCF（PoE 融合） | B1 + B2 的输入 |
| B4 祖先 k-mer + 长读长 | B1 + `--longread longreads.sorted.bam` |
| B5 MSA SNP + 长读长 | `--msa-snp --longread longreads.sorted.bam` |
| B6 祖先 k-mer + MSA SNP + 长读长（PoE） | B1 + `--msa-snp --longread ...` |
| B7 VCF + 长读长 | B2 + `--longread ...` |
| B8 祖先 k-mer + VCF + 长读长（PoE） | B3 + `--longread ...` |

任意 B 路线都可以把 `--ancestor ...` 换成预建 marker（`--markers SOF:SOF_31.fa --markers SRO:SRO_31.fa ...`），
跳过 KMC 建库。`examples/` 里有这些路线的现成脚本。

---

## 5. 全部参数说明

必需参考文件：`--fasta`、`--gff`、`--clusters`、`--ancestry`、`--orphan`、`--target-ids`（见第 3 节）。

倍性 + 血统输入：

| 参数 | 说明 |
|---|---|
| `--ploidy {auto,allo}` | auto = 同源（只拷贝级）；allo = 异源（血统分配）。默认 allo |
| `--ancestor NAME:fa,gff` | 祖先基因组（可重复） |
| `--markers NAME:fasta` | 预建血统 marker（可重复） |
| `--vcf` / `--sample-lists NAME:file` | 群体 VCF 路线 |
| `--phased-fasta FILE` / `--msa-snp` | MSA SNP 路线 |
| `--longread FILE` | 长读长 BAM |
| `--auto-idxstats` | 仅对照用：同源模式下 比对 → idxstats → TPM，不跑 Theta / EM，等价于 featureCounts -M -O；必须配 `--aligner hisat2` |

输入输出：`-d/--data-dir`（FASTQ 与共用索引，默认 `./data`）、`-o/--out-dir`（默认 `./results`）。

比对（Step 00）：

| 参数 | 说明 |
|---|---|
| `--aligner {bowtie2,hisat2}` | 默认 bowtie2（RSEM 参数，N 不罚分）；hisat2 = 1.0 之前的比对方式，更快 |
| `--max-alignments N` | 每个片段最多报几个比对（默认 bowtie2 200、hisat2 15） |

高级：环境变量 `LQ_BT2_OPTS` 整体替换 bowtie2 参数串，`LQ_HISAT_OPTS` 替换 HISAT2 的。

EM（Step 03）：

| 参数 | 说明 |
|---|---|
| `--em-engine {global,legacy}` | 默认 global（全部基因一起 EM）；legacy = 1.0 之前的按簇 EM + 单拷贝基因用 idxstats。`v6` / `v5` 是别名 |
| `--nm-weight {auto,off,λ}` | 按错配数加权，默认 auto；off = 候选一视同仁；也可直接给 λ（0–1） |
| `--eff-len {eff,len,off}` | 有效长度校正，默认 eff |
| `--no-pool-identical` | 不合并序列完全相同的基因 |
| `--p-weight-specific / -conflict / -none` | 按血统证据状态给血统概率加指数权重（默认 1 / 0 / 0） |
| `--p-alpha`、`--p-floor` | 血统概率向均匀收缩的强度（默认 0）与下限（默认 1e-3） |
| `--theta-prior-lr / -kmer / -cap` | Theta 作为 Dirichlet 先验的权重（长读长锚定的簇）与上限 |
| `--em-epsilon`、`--em-max-iter`、`--em-min-iter`、`--em-tol`、`--em-ll-tol`、`--cross-cluster` | 只对 legacy 引擎有效 |

运行 / Theta / 血统：

| 参数 | 说明 |
|---|---|
| `-t, --threads N` | 线程数（默认 16） |
| `--kmer N` | Theta k-mer 长度（默认 21） |
| `--theta-threads N` | Theta 扫描 FASTQ 的进程数（默认 8，不超过 `--threads`）；每个进程都持有一份 k-mer 索引，大转录组上调大会很占内存 |
| `--lineage-kmer N` | 血统 marker k-mer 长度（默认 31） |
| `--shrinkage F` / `--adaptive-shrinkage` | Theta 收缩强度（默认 10.0）/ 簇自适应收缩 |
| `--no-paired-bonus` | 关闭双端配对加成 |
| `--min-lr-depth N` | 用长读长锚定一个簇的最少长读数（默认 5） |
| `--error-rate F` | SNP 路线的碱基错误率（默认 0.01） |
| `--no-calibrate-snp`、`--calib-min-depth`、`--calib-det-weight`、`--no-chimeric` | MSA 路线 + 长读长的校准选项 |
| `--dry-run` | 只打印路由，不执行 |

---

## 6. 辅助命令

**normalize-clusters**：把建簇工具的输出标准化。

```bash
lineagequant normalize-clusters \
  --input jcvi_pairs.tsv --format pairs \
  --clusters-out clusters.tsv --ancestry-out ancestry.tsv \
  --subgenome-tag A:_A --subgenome-tag D:_D \
  --orphan-out orphans.tsv --all-genes all_gene_ids.txt
```

`--format`：`standard`（已是标准格式，只校验和重编号）、`jcvi`、`pairs`（两列 pair，合并连通分量，
支持 > 2 个亚基因组）、`tsv`。`--subgenome-tag NAME:REGEX` 按基因 ID 正则分血统（可重复）。

**build-markers**：只构建血统 k-mer，不跑全流程。

```bash
lineagequant build-markers \
  --lineage A:ancestorA_transcript.fa --lineage D:ancestorD_transcript.fa \
  -o markers_out -k 31 -t 16
```

输出 `markers_out/A.fasta`、`markers_out/D.fasta`，之后 `run` 用 `--markers` 复用。
`--polyploid NAME:FASTA` 与 `--mode {ancestor,ancestor_x_poly,poly}` 可以要求 marker
在多倍体自己的各亚基因组转录本之间也是特异的。

---

## 7. 输出文件

每个样本（`<out-dir>/<sample>/`）：

| 文件 | 内容 |
|---|---|
| `<sample>.sorted.bam`（+ `.bai`） | 比对到转录组的排序 BAM（包含全部报出的比对） |
| `<sample>.bowtie2.log` / `.hisat2.log` | 比对器摘要 |
| `aligner.txt` | 这个样本目录用的比对器 |
| `Global_Theta_Matrix.tsv` | 每拷贝 Theta 先验 |
| `P_values_kmer.tsv` / `P_values_snp.tsv` / `P_values_merged.tsv` | 每条 read 的血统概率（异源） |
| `EM_Final_Expression.tsv`（有长读长时为 `EM_Final_Hybrid.tsv`） | `ClusterID, AlleleID, Ancestry, Initial_Theta, EM_Allocated_Reads, Final_Theta, EM_Iter, EM_LogLik`；单拷贝基因的 ClusterID 为 `ORPHAN` |

数据目录（`<data-dir>/`）里可复用的中间产物：`transcript_index_bt2.*`（bowtie2 索引）、
`transcript_index.*.ht2`（HISAT2 索引，仅 `--aligner hisat2`）、`transcript_index_k*.pkl`（Theta k-mer 索引）、
`lineage_markers/{NAME}.fasta`（建好的 marker）、`SNP_DB_*.pkl`、`LongRead_Counts.tsv`。

最终（`<out-dir>/LineageQuant_Final_Merged_TPM.tsv`）：

| 列 | 内容 |
|---|---|
| `GeneID` | 基因 / 等位拷贝 |
| `Type` | `Polyploid(EM)` 等位拷贝；`Orphan(EM)` 单拷贝基因；`Orphan(BAM_Base)` / `Polyploid(idxstats)` 由 idxstats 计数（兜底或对照模式） |
| `Length_bp` | 最长转录本长度 |
| `<sample>_RawReads` | EM 分配到的片段数（做差异表达等计数分析用这一列） |
| `<sample>_TPM` | RawReads 除以 Length_bp（kb）后缩放到一百万 |

亚基因组表达 = 同一簇内同一血统各拷贝的 RawReads 之和（血统见 `--ancestry`）。

---

## 8. 常见问题

**Q：路由走错了？**
先 `--dry-run` 看 resolved route，检查 `--ploidy` 和对应输入。异源没给任何血统路线会报错；
`--vcf` 与 `--phased-fasta` / `--msa-snp` 不能同时给。

**Q：中断了要重跑？**
直接重跑同样命令，已完成的步骤会跳过。要强制重跑某一步，删掉对应的 `.done` 文件。

**Q：报 `this output directory was aligned with hisat2`？**
这个输出目录是以前版本（HISAT2）跑的。1.0.0 不允许在同一个目录里混用两种比对器：
换一个新的 `--out-dir`，或者加 `--aligner hisat2` 接着用旧目录。

**Q：怎么复现 1.0 之前的结果？**
`--aligner hisat2 --em-engine legacy --eff-len off`。`--auto-idxstats` 对照也必须配 `--aligner hisat2`。

**Q：bowtie2 太慢？**
bowtie2 `-k 200` 的比对时间约为 HISAT2 的 3–10 倍，BAM 约大 1.5–3 倍。更看重速度时用
`--aligner hisat2`，代价是近乎相同的拷贝之间分配不如 bowtie2 准。

**Q：有一部分片段没被计数？**
比对是端到端的且要求两端成对：一端落在注释转录本之外（未注释的同工型、UTR 延伸）的片段不计数。
请尽量提供完整的转录组（每个等位拷贝一条记录）。

**Q：Theta k-mer 和血统 k-mer 有什么区别？**
`--kmer`（默认 21）区分簇内拷贝、给 EM 初值，两种倍性都用；`--lineage-kmer`（默认 31）区分血统，只异源用。

**Q：同源多倍体能用祖先 k-mer 吗？**
不能。同源亚基因组间无分化，祖先 k-mer 相减会减没。同源只做拷贝级定量（`--ploidy auto`）。

**Q：工具报 not on PATH？**
`pip` 不装 conda 工具。用 `conda env create -f environment.yml` 装全，或单独 `conda install -c bioconda <工具>`。
