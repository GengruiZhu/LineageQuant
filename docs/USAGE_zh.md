# LineageQuant 使用说明（v1.1.0）

多倍体亚基因组 / 等位拷贝表达定量软件。作者：Yi Chen, Gengrui Zhu

本文档讲：安装与检查、输入准备、按数据组合怎么跑、全部参数、输出文件、常见问题。
英文版见 `README.md`；改动见 `CHANGELOG.md`；设计见 `docs/DESIGN.md`（路由）与 `docs/ENGINE.md`（引擎）；
二进制格式见 `docs/FORMATS.md`。

---

## 0. 这个软件做什么

LineageQuant 定量多倍体里"每个亚基因组 / 等位拷贝"的表达量。它是**定量器**，不是 read 分类器：
参考转录组里每个同源 / 等位拷贝都是独立目标，兼容多个拷贝的 read 由 EM 按概率分配，
最终输出拷贝级的 read 数与 TPM，可以再按血统加和得到亚基因组表达。

自动路由：你只声明物种类型（同源 / 异源），再把手头有的输入交给它，软件自己决定走哪条流程。
缺什么输入就跳过对应模块，所有路线最后都汇到 EM。

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

### 1.1.0 与 1.0.0 的区别

- **只有一个 C++ 程序**：命令是 `lq`（全名 `lineagequant` 也能用，两者完全相同；`bash build.sh` 编译，十几秒）。没有 Python、没有 pip，
  不再需要 numpy / pandas / Biopython / pysam / cyvcf2。外部工具照旧调用（bowtie2 或 salmon、KMC、gffread、
  MAFFT，长读长时 samtools）。`bin/lqcore` 也是同一个程序，旧脚本里的 `lqcore run`（单个样本）照常可用。
- **默认参数下结果与 1.0.0 相同**：`tests/run_all.sh` 把每一步都和 1.0.0 的 Python 逐字节对比——EM 36 组配置、
  谱系整条链和分块、Theta 36 组配置、长读计数 80 组、建库各步（normalize-clusters、extract-transcripts、
  build-markers 三种模式、snpdb-msa、snpdb-vcf）、表达矩阵、整条流程 4 条路线。刻意的不同只有两处（见第 8 节）。
- **定量只写二进制**：`Global_Theta_Matrix.lqb` 与 `<样本>.lqres`（EM + 有效长度 + 每条参考的比对数）。
  人类可读的文本（与 1.0.0 的文件逐字节相同）默认在每个样本定量结束后由二进制转写出来，**不计入定量耗时**
  （日志 `[DONE] quantification … | text files …, not counted` 分开报）；`--no-text` 只留二进制，
  事后也可以用 `lq export` 转。
- **表达矩阵不在流程里**：`lq run` 跑到每个样本的 EM 为止（kallisto、Salmon 也是这样），
  全基因组矩阵之后用 `lq matrix` 生成，与 1.0.0 的 Step 04 输出同一个文件。
- **快得多**：POJ2878 root（约 3,800 万对 reads，32 线程）1.0.0 Python 引擎 6 小时，1.1.0 约 1 小时；
  用 `--aligner lqmap --em-accel squarem` 一个 POJ 组织 15–30 分钟。
- **新增、默认关闭的选项**：`--eff-len fld`、`--gc-bias`、`--marker-clean`、`matrix --tpm-length eff`（第 5 节）。

### 定量方式（沿用 1.0.0）

- **Step 00 比对**：默认用 bowtie2，参数与 RSEM 默认相同（端到端、不开空位、`--mp 1,1`、
  `--score-min L,0,-0.1`、只保留成对一致的比对、每个片段最多 200 个比对位置），参考序列中的
  N 不罚分（`--np 0 --n-ceil L,0,0.5`）。它会把同样好的拷贝全部报出来，也报出稍差的位置。
  bowtie2 的输出边产生边读，不排序、不写 BAM、不写逐 read 的概率表。
  `--aligner lqmap`（自带比对器，规则相同）和 `--aligner sa`（Salmon 提候选、再按同样规则验证）更快。
- **Step 03 EM（global 引擎）**：全部基因（等位拷贝与单拷贝基因）一起做 EM，每个片段总权重为 1，
  不会重复计数；有效长度校正；候选比对按错配数加权（每多一个错配，似然乘以 λ = (p/3)/(1−p)，
  p 由最佳比对的错配率估计）；序列完全相同的基因合并估计、收敛后平分。

---

## 1. 安装

```bash
cd /path/to/LineageQuant_v1.1.0       # 解压后的目录（里面有 build.sh、src/、tests/）
conda env create -f environment.yml   # 外部工具 + 编译器（推荐；也可以自己装）
conda activate lineagequant
bash build.sh                         # 编译 -> bin/lq、bin/lineagequant、bin/lqcore（同一程序；需要 GCC >= 7 和 zlib）
export PATH=$PWD/bin:$PATH            # 写进 ~/.bashrc 或作业脚本

lq version                  # 应输出: lq 1.1.0
lq run --help
```

**集群上一定要在集群自己的登录节点上编译**：在别的（更新的）系统上编出来的程序，计算节点会报
`GLIBC_2.xx not found`。`build.sh` 把新程序装成 `bin/lq.<时间>`、`bin/lineagequant.<时间>`、`bin/lqcore.<时间>`（同一程序的三份拷贝），再原子地把
`bin/lq`、`bin/lineagequant`、`bin/lqcore` 三个链接指过去，旧版本保留（网络文件系统上删掉或覆盖正在别的节点运行的程序，
那边会直接崩）。也可以 `make && make install PREFIX=$HOME/.local`。

别人怎么装（两种方式，结果相同）：

- **从源码编译**（上面这样）：需要 GCC ≥ 7 和 zlib，十几秒。
- **直接下载编好的程序**：GitHub release 里的 `lq-1.1.0-linux-x86_64-static.tar.gz`（`bash build.sh --static`
  编出来的完全静态程序，不依赖 glibc / libstdc++ / zlib；x86-64 Linux、内核 ≥ 3.2、CPU 支持 SSE4.2/POPCNT
  ——2009 年以后的 x86-64 服务器 CPU 都满足）。解压后里面是 `lq`、`lineagequant`、`lqcore`（同一文件的硬链接）：
  `tar xzf lq-1.1.0-linux-x86_64-static.tar.gz && export PATH=$PWD/lq-1.1.0-linux-x86_64-static:$PATH`。
  外部工具仍按下表自己装（或 `conda env create -f environment.yml`）。

外部工具（按路线需要；`lq run` 会检查本次路线要用的工具，缺哪个就报哪个，例如
`[FATAL] tool 'bowtie2' not on PATH`）：

| 工具 | 什么时候需要 |
|---|---|
| bowtie2 ≥ 2.4 | `--aligner bowtie2`（默认） |
| salmon | `--aligner sa` |
| kmc + kmc_tools、gffread | `--ancestor`（建 marker） |
| mafft | `--msa-snp` / `--phased-fasta` |
| samtools | `--longread`（读长读 BAM） |
| minimap2 | 生成长读 BAM |

`--aligner lqmap` 不需要任何外部工具。

建等位簇的工具（GraphAllele / AlleleFinder / jcvi 等）不在依赖里。自己产出等位簇后，用
`lq normalize-clusters` 标准化（见第 6 节）。

---

## 2. 安装后怎么检查

```bash
# [1] 与 1.0.0 的等价性（需要一份 LineageQuant 1.0.0 的目录和带 numpy、pandas 的 python3；合成数据，约 3 分钟）
bash tests/run_all.sh /path/to/LineageQuant_v1.0.0/lineagequant
#   ...
#   ================ summary
#     EM (36 configurations) (27s)                       ok
#     ...
#     whole pipeline (4 routes) (13s)                    ok
#   ALL TESTS PASSED

# [2] 路由是否正确（--dry-run 只打印路由，不真跑）
lq run --fasta mRNA.fa --clusters clusters.tsv --ancestry ancestry.tsv --target-ids ids.txt \
  --ploidy allo --ancestor A:ancA.fa,ancA.gff --ancestor D:ancD.fa,ancD.gff \
  --data-dir data --out-dir results --dry-run
```

应输出类似：

```
[LineageQuant 1.1.0] resolved route:
ploidy: allo
  alignment: bowtie2 (RSEM settings, N not penalised, <= 200 alignments/fragment), streamed into quant, no BAM
  engine: /path/to/LineageQuant_v1.1.0/bin/lq.20261027_090000
  EM: global over all genes | eff-len=eff | mismatch weighting=auto | identical-sequence pooling=on
  outputs: <sample>/Global_Theta_Matrix.lqb, <sample>/<sample>.lqres + the 1.0.0 text files (written after each sample, not timed)
  lineage routes: ancestor k-mer (2 lineages, 2 progenitor sample(s))
  Theta: k-mer (01)
[dry-run] exiting without execution.
```

测试里 MAFFT、KMC、gffread、Biopython、pysam、cyvcf2 都换成了 `tests/` 里的替身程序（两边用同一套），
所以不需要装这些工具也能跑。

---

## 3. 输入文件准备

必需：

| 参数 | 内容 |
|---|---|
| `--fasta` | 转录组（mRNA）FASTA；每个等位拷贝是一条独立记录。可从基因组 + GFF 用 gffread 提取 |
| `--clusters` | 标准化等位簇，每行一簇：`ClusterID<TAB>Allele1,Allele2,...` |
| `--ancestry` | `GeneID<TAB>LineageName`；血统名从第 2 列动态读取，数量与名字不限 |
| `--target-ids` | 目标基因 ID，每行一个（异源必需） |

其他：`--orphan`（单拷贝基因 `OrphanID<TAB>GeneID`；流程本身不用，`lq matrix` 要用，给了会原样提示）、
`--gff`（参考 GFF3，只有群体 VCF 路线需要）。

基因 ID 由转录本 ID 去掉最后一个 `.N` 后缀得到（`Gene.1`、`Gene.2` 是 `Gene` 的同工型）。

血统路线输入（异源，至少给一种）：

- `--ancestor NAME:fasta,gff`：祖先基因组，KMC 相减得血统特异 k-mer；可重复，≥ 2 个。
  同一个 NAME 重复给多次 = 该血统有多个祖先样本，合并后再建 marker。
  基因组里 ≥ 1 Mb 的序列才算染色体；GFF 的序列名与 FASTA 不一致时自动对应（同名 → 规范化名字如
  `Chr01 ~ chr1 ~ CHR_01` → 长度相差 ≤ 0.5% 且唯一）；链为 `?` 的特征丢掉。
- `--markers NAME:fasta`：预建血统 marker，跳过 KMC 建库。NAME 要与 `--ancestry` 的血统名一致；
  marker 的 k-mer 长度须与 `--lineage-kmer` 一致。
- `--vcf FILE` + `--sample-lists NAME:file ...` + `--gff`：群体 VCF（纯文本或 bgzip，不支持 BCF），每个血统一个样本名文件。
- `--msa-snp`：MSA SNP 路线，直接用主 `--fasta`（主转录组已含所有等位拷贝时用这个）；
  `--phased-fasta FILE` 只在等位拷贝 CDS 在单独文件里时用。两者互斥，也都与 `--vcf` 互斥。

**marker 的来源**：血统 marker 必须来自独立于参考组装的序列（祖先基因组 / 转录组）。不要用参考转录组本身
"按亚基因组相减"来造 marker 再拿去救参考里缺失的拷贝：犬蔷薇留出实验里，从缺拷贝的参考重建的 marker
在需要它的 read 上反而指错祖先（缺失拷贝与另一祖先拷贝共有的 k-mer 变成了"另一祖先特异"）。

**marker 不是越多越好**：二倍体、同源多倍体和参考完整的简单异源多倍体上，比对本身已经能区分拷贝，
血统 marker 不提供额外信息，开了反而可能变差；所以任何 marker 路线都不默认打开，按物种情况选择。

长读长（可选）：`--longread` 接受比对到转录组的 BAM，任意排序，不需要索引；作用是锚定 Theta：

```bash
minimap2 -ax map-hifi -N 50 -p 0.95 transcriptome.fa reads.fq.gz | samtools view -b -o longreads.bam -
```

subreads 用 `-x map-pb`。`-N 50 -p 0.95` 让多拷贝间打平的比对作为 secondary 写出来，长读计数靠它判断平局。
每条分子只计一次：secondary 不计，同一 ZMW 的多条 subread 合计 1，同工型加到拷贝上。
多个 FASTQ 文件请先 `cat` 成一条流再给 minimap2（minimap2 把多个 query 文件当成同一片段的几段）。

FASTQ：放在 `--data-dir` 目录里，双端 gzip，自动识别 `.R1.fq.gz`、`_R1.fq.gz`、`_1.fq.gz`、
`_R1.fastq.gz`、`_R1_001.fastq.gz` 等命名。共用的索引（bowtie2 / lqmap / Theta / marker）也写在这个目录里。

---

## 4. 怎么跑（按数据组合选，10 条路线）

先 `--dry-run` 看路由，再去掉它正式跑。所有步骤可断点续跑（`.done` 标记），中断后重跑同样命令即可。
跑完之后用 `lq matrix` 出表达矩阵（第 7 节）。

### 同源多倍体（auto）

**A1 同源 + 仅短读长（拷贝级定量）**

```bash
lq run \
  --fasta transcriptome.fasta --clusters clusters.tsv --ancestry ancestry.tsv \
  --ploidy auto --data-dir data --out-dir results --threads 32
lq matrix -r results -o orphans.tsv
```

**A2 同源 + 长读长**：在 A1 上加 `--longread longreads.bam`。

### 异源多倍体（allo）

**B1 异源 + 祖先 k-mer**

```bash
lq run \
  --fasta transcriptome.fasta \
  --clusters clusters.tsv --ancestry ancestry.tsv --target-ids target_gene_ids.txt \
  --ploidy allo \
  --ancestor A:ancestorA.fa,ancestorA.gff \
  --ancestor D:ancestorD.fa,ancestorD.gff \
  --data-dir data --out-dir results --threads 32
lq matrix -r results -o orphans.tsv
```

| 路线 | 在基本命令上的输入 |
|---|---|
| B2 群体 VCF | `--vcf population.vcf.gz --sample-lists A:samples_A.txt D:samples_D.txt --gff genome.gff3` |
| B3 祖先 k-mer + VCF（PoE 融合） | B1 + B2 的输入 |
| B4 祖先 k-mer + 长读长 | B1 + `--longread longreads.bam` |
| B5 MSA SNP + 长读长 | `--msa-snp --longread longreads.bam`（不给长读长也能跑） |
| B6 祖先 k-mer + MSA SNP + 长读长（PoE） | B1 + `--msa-snp --longread ...` |
| B7 VCF + 长读长 | B2 + `--longread ...` |
| B8 祖先 k-mer + VCF + 长读长（PoE） | B3 + `--longread ...` |

任意 B 路线都可以把 `--ancestor ...` 换成预建 marker（`--markers SOF:SOF_31.fa --markers SRO:SRO_31.fa ...`），
跳过 KMC 建库。`examples/` 里有这些路线的现成脚本。

大数据提速：加 `--aligner lqmap --em-accel squarem`（lqmap 与 bowtie2 规则相同；SQUAREM 收敛到同一个极大值，
差别在收敛容差量级）。

**一个样本拆到多个节点**（read 证据可加和）：

```bash
# 每块一个作业（i = 1..N）
lq quant --r1 R1.fq.gz --r2 R2.fq.gz --aligner lqmap --map-cache data/lqmap.k21.bin --fasta mRNA.fa \
  -a ancestry.tsv -t results/S/Global_Theta_Matrix.lqb [血统选项] --chunk i/N --part-out S.part$i --threads 32
# 全部块完成后
lq em -a ancestry.tsv -t results/S/Global_Theta_Matrix.lqb --fasta mRNA.fa --parts S.part1 ... S.partN -o results/S/S.lqres
```

---

## 5. 全部参数说明

`lq run -h` 列出全部；下面按用途分组。

倍性 + 血统输入：

| 参数 | 说明 |
|---|---|
| `--ploidy {auto,allo}` | auto = 同源（只拷贝级）；allo = 异源（血统分配）。默认 allo |
| `--ancestor NAME:fa,gff` | 祖先基因组（可重复） |
| `--markers NAME:fasta` | 预建血统 marker（可重复） |
| `--marker-clean` | （1.1.0，默认关）去掉在"另一血统"的参考拷贝里出现过的 marker k-mer；参考里没有的 k-mer 保留。见下 |
| `--vcf` / `--sample-lists NAME:file` / `--gff` | 群体 VCF 路线 |
| `--phased-fasta FILE` / `--msa-snp` | MSA SNP 路线 |
| `--legacy-lineages` | MSA SNP 库用 1.0.0 的血统列表（只为复现 1.0.0；见第 8 节） |
| `--longread FILE` | 长读长 BAM |

输入输出：`-d/--data-dir`（FASTQ 与共用索引，默认 `./data`）、`-o/--out-dir`（默认 `./results`）、
`--no-text`（只写二进制；默认每个样本定量结束后另外写出 1.0.0 的文本文件，不计入耗时）。

比对（Step 00，在 `quant` 里）：

| 参数 | 说明 |
|---|---|
| `--aligner {bowtie2,lqmap,sa}` | 默认 bowtie2（RSEM 参数，N 不罚分，输出流式读入）；lqmap = 自带（同样规则，快）；sa = Salmon 提候选、再验证 |
| `--max-alignments N` | bowtie2 每个片段最多报几个比对（默认 200） |
| `--salmon PATH` | `--aligner sa` 用的 salmon |

高级：环境变量 `LQ_BT2_OPTS` 整体替换 bowtie2 参数串。

EM（Step 03）：

| 参数 | 说明 |
|---|---|
| `--nm-weight {auto,off,λ}` | 按错配数加权，默认 auto；off = 候选一视同仁；也可直接给 λ（0–1） |
| `--eff-len {eff,fld,len,off}` | 有效长度。eff（默认，1.0.0）= 长度 − 平均片段长 + 1；fld（1.1.0）= 长度 − 观测片段长分布截断到该长度后的均值 + 1（kallisto / Salmon 的做法） |
| `--gc-bias` | （1.1.0，默认关）片段 GC 偏好校正，见下 |
| `--em-accel {none,squarem}` | none（默认）= 1.0.0 的 EM；squarem = 加速，结果差在收敛容差量级 |
| `--no-pool-identical` | 不合并序列完全相同的基因 |
| `--p-weight-specific / -conflict / -none` | 按血统证据状态给血统概率加指数权重（默认 1 / 0 / 0） |
| `--p-alpha`、`--p-floor` | 血统概率向均匀收缩的强度（默认 0）与下限（默认 1e-3） |
| `--theta-prior-lr / -kmer / -cap` | Theta 作为 Dirichlet 先验的权重（长读长锚定的簇）与上限 |

运行 / Theta / 血统：

| 参数 | 说明 |
|---|---|
| `-t, --threads N` | 线程数（默认 16） |
| `--kmer N` | Theta k-mer 长度（默认 21） |
| `--lineage-kmer N` | 血统 marker k-mer 长度（默认 31） |
| `--shrinkage F` / `--adaptive-shrinkage` | Theta 收缩强度（默认 10.0）/ 簇自适应收缩 |
| `--no-paired-bonus` | 关闭双端配对加成 |
| `--min-lr-depth N` | 用长读长锚定一个簇的最少长读分子数（默认 5） |
| `--theta-sum {compensated,plain}` | Theta 里浮点求和的方式。1.0.0 跟着运行它的 Python：≥ 3.12 是补偿求和（= compensated，默认），≤ 3.11 是普通求和（= plain）。只影响最后一位 |
| `--error-rate F` | SNP 路线的碱基错误率（默认 0.01） |
| `--no-chimeric` | SNP 路线：不丢弃谱系嵌合的 read |
| `--dry-run` | 只打印路由，不执行 |

1.0.0 里已经不起作用的选项（`--engine cpp`、`--lqcore`、`--theta-threads`、长读校准选项、旧 EM 的
`--em-epsilon` 等）照样能写，不报错；`--tpm-length` 现在属于 `matrix`，写在 `run` 上会提示。

**`--gc-bias` 做什么**（Salmon 的前景 / 背景 GC 模型，单一 GC 模型）：

1. 扫描 read 时记下每个片段（主比对的一致配对区间）的 GC 百分比 = 前景分布；
2. 先跑一遍 EM，按各目标的期望片段数，统计"每个表达目标所有可能片段"（片段长取观测分布的分位数）的 GC 分布 = 背景；
3. 权重 = 前景 / 背景（拉普拉斯平滑，限制在 0.01–100，背景上均值为 1）；每个目标的有效长度乘以它所有可能片段的平均权重；
4. EM 从第一遍的解继续，直到再次收敛。

合成数据（`tests/compare_new.sh`）：强 GC 偏好下 TPM 的 MARD 0.2055 → 0.1100，弱偏好 0.0986 → 0.0756，
无偏好不变。只校正片段整体 GC，不做 5′/3′ 序列偏好。默认关闭，等人类数据（Polyester 强 / 弱 GC）验证后再决定默认值。

**`--marker-clean` 做什么**：从祖先基因组建的 marker 里，凡是在"另一血统"的参考拷贝中出现过的 k-mer 都去掉。
祖先与多倍体里的亚基因组已经分化：亚基因组共有的位点，若某个祖先在那里突变过，这段 k-mer 就会变成另一个祖先的
"特异" marker，于是给共有区域的 read 打上错误的血统。犬蔷薇留出实验（祖先 marker 来自加了 0.5–2% 随机替换的序列）里，
来源拷贝在参考里的 read 有 1–2% 被这样判错。参考里根本没有的 k-mer 保留——它们正是能指向参考缺失拷贝的那部分。

---

## 6. 其他命令

都在同一个程序里，`lq <命令> -h` 看各自的参数。

**normalize-clusters**：把建簇工具的输出标准化。

```bash
lq normalize-clusters \
  --input jcvi_pairs.tsv --format pairs \
  --clusters-out clusters.tsv --ancestry-out ancestry.tsv \
  --subgenome-tag A:_A --subgenome-tag D:_D \
  --orphan-out orphans.tsv --all-genes all_gene_ids.txt
```

`--format`：`standard`（已是标准格式，只校验和重编号）、`jcvi`、`pairs`（两列 pair，合并连通分量，
支持 > 2 个亚基因组）、`tsv`。`--subgenome-tag NAME:REGEX` 按基因 ID 正则分血统（可重复；ECMAScript 语法，
常用的 `^前缀`、`_A$` 与 Python 相同）。

**extract-transcripts**：祖先基因组 + GFF → 转录本 FASTA（`run` 建 marker 时自动调用）。

```bash
lq extract-transcripts --fasta ancA.fa --gff ancA.gff --workdir work --out A_tx.fa --id-prefix A
```

**build-markers**：只构建血统 k-mer，不跑全流程。

```bash
lq build-markers \
  --lineage A:ancestorA_transcript.fa --lineage D:ancestorD_transcript.fa \
  -o markers_out -k 31 -t 16
```

输出 `markers_out/A.fasta`、`markers_out/D.fasta`，之后 `run` 用 `--markers` 复用。
`--polyploid NAME:FASTA` 与 `--mode {ancestor,ancestor_x_poly,poly}` 可以要求 marker
在多倍体自己的各亚基因组转录本之间也是特异的（给了 `--polyploid` 时默认 ancestor_x_poly，否则 ancestor）。

**snpdb-msa / snpdb-vcf**：只建 SNP 库（`run` 会自动调用）。

```bash
lq snpdb-msa --fasta mRNA.fa --cluster clusters.tsv --ancestry ancestry.tsv --out SNP_DB_msa.lq.tsv --threads 16
lq snpdb-vcf --vcf pop.vcf.gz --gff3 genome.gff3 --sample-lists A:A.txt D:D.txt --out SNP_DB_population.lq.tsv
```

**export**：二进制 → 1.0.0 的文本文件（逐字节相同）。

```bash
lq export results/S1/S1.lqres                     # -> EM_Final_Expression.tsv, S1.idxstats.tsv, S1.efflen.tsv
lq export results/S1/Global_Theta_Matrix.lqb      # -> Global_Theta_Matrix.tsv
```

---

## 7. 输出文件

每个样本（`<out-dir>/<sample>/`）：

| 文件 | 内容 |
|---|---|
| `<sample>.lqres` | 二进制结果：每个目标的 EM（初始 Theta、EM 分配的片段数、最终 Theta）、EM 用的有效长度、每条参考序列的比对数 |
| `Global_Theta_Matrix.lqb` | 二进制的每拷贝 Theta 先验 |
| `EM_Final_Expression.tsv` 等文本 | 默认由上面两个二进制转写（见下），`--no-text` 时没有 |
| `<sample>.bowtie2.log` / `.salmon.log` | 比对器摘要 |
| `aligner.txt` | 这个样本目录用的比对器 |

默认写出（或 `lq export` 还原出）的文本文件：`EM_Final_Expression.tsv`（有长读长时为
`EM_Final_Hybrid.tsv`；`ClusterID, AlleleID, Ancestry, Initial_Theta, EM_Allocated_Reads, Final_Theta, EM_Iter, EM_LogLik`；
单拷贝基因的 ClusterID 为 `ORPHAN`）、`<sample>.idxstats.tsv`（samtools idxstats 格式）、`<sample>.efflen.tsv`
（`Target, Length, EffectiveLength, GC_Factor`）、`Global_Theta_Matrix.tsv`。

不写 BAM，也不写逐 read 的血统概率表；需要逐片段诊断时用 `lq quant --debug-frag FILE`。

数据目录（`<data-dir>/`）里可复用的中间产物：`transcript_index_bt2.*`（bowtie2 索引）、`lqmap.k21.bin`（lqmap 索引）、
`lqcore_theta.k*.bin`（Theta k-mer 索引）、`ancestor_transcripts/`、`lineage_markers/{NAME}.fasta`（建好的 marker）、
`lqcore_markers.k*[.clean].bin`（marker 索引）、`SNP_DB_*.lq.tsv`、`LongRead_Counts.tsv`。

**表达矩阵**：`lq matrix -r <out-dir> -o orphans.tsv`（几秒；不算在定量耗时里），
写 `<out-dir>/LineageQuant_Final_Merged_TPM.tsv`，与 1.0.0 Step 04 的文件相同：

| 列 | 内容 |
|---|---|
| `GeneID` | 基因 / 等位拷贝 |
| `Type` | `Polyploid(EM)` 等位拷贝；`Orphan(EM)` 单拷贝基因；`Orphan(BAM_Base)` 由每条参考的比对数计数（兜底） |
| `Length_bp` | 最长转录本长度 |
| `<sample>_RawReads` | EM 分配到的片段数（做差异表达等计数分析用这一列） |
| `<sample>_TPM` | RawReads 除以 Length_bp（kb）后缩放到一百万；`--tpm-length eff` 时除以有效长度 |

`matrix` 也能读 1.0.0 的文本输出目录；`--mode idxstats` 是 1.0.0 `--auto-idxstats` 那种对照（每个报出的比对都算满，
等于 featureCounts `-M -O`，不经过 LineageQuant 的任何算法）。

亚基因组表达 = 同一簇内同一血统各拷贝的 RawReads 之和（血统见 `--ancestry`）。
参考里缺拷贝时（组装塌缩 / 漏注释），拷贝层面的分配会被带偏，而血统（祖先）层面的比例稳得多
（犬蔷薇留出实验、甘蔗 Iso-Seq 对照），解读时以血统层面为主。

---

## 8. 与 1.0.0 的两处刻意不同

1. **多重比对片段的 SNP 证据**：一条 read 比到多个转录本时，只取一个转录本的 SNP 证据，规则固定为错配最少 →
   覆盖 SNP 最多 → 参考序号最小。1.0.0 取中间文件里最后一行，由多进程的先后决定。
2. **MSA SNP 库的血统列表**：1.0.0 的 `export_snp_db.py` 只取库里第一个 SNP 的血统作为整个库的血统列表，
   `build_snp_msa.py` 又用"去掉最后一个 . 之后的部分"去查拷贝的血统（而筛选拷贝时用的是"第一个 . 之前的部分"）。
   第一个比对上的簇里缺某个血统，或者转录本 ID 里有多个点时，那个血统的 SNP 证据会悄悄全部丢掉
   （POJ 上出现过：库里只有 SOF 和 SSP）。1.1.0 用 ancestry 里的全部血统、一种 ID 规则；
   `--legacy-lineages` 可以完全复现 1.0.0。第一个簇含全部血统、ID 只有一个点时两者完全相同。
   `run` 发现数据目录里是 1.0.0 方式建的库时会给出警告（删掉库和它的 `.done` 就会按新方式重建）。

---

## 9. 常见问题

**Q：报 `GLIBC_2.xx not found`？**
程序是在别的机器上编的。在要运行的集群上执行 `bash build.sh`，或者改用 GitHub release 里的静态程序。

**Q：报 `--r1 is an option of one sample`？**
`lq run` 是整条流程；单个样本用 `lq quant`（或 `lqcore run`）。

**Q：路由走错了？**
先 `--dry-run` 看 resolved route，检查 `--ploidy` 和对应输入。异源没给任何血统路线会报错；
`--vcf` 与 `--phased-fasta` / `--msa-snp` 不能同时给。

**Q：中断了要重跑？**
直接重跑同样命令，已完成的步骤会跳过。要强制重跑某一步，删掉对应的 `.done` 文件。
1.0.0 用 `--engine cpp` 跑出来的样本目录也能接着跑（沿用里面的文本 Theta / EM）。

**Q：报 `this output directory was aligned with ...`？**
这个输出目录是用另一种比对器跑的，同一个目录里不允许混用：换一个新的 `--out-dir`。

**Q：想看文本结果？**
默认就有（每个样本目录里的 `EM_Final_Expression.tsv` 等）；用了 `--no-text` 的话，`lq export <样本>.lqres`。

**Q：怎么复现 1.0.0 / 1.0 之前的结果？**
1.0.0：1.1.0 的默认参数即可（MSA SNP 路线加 `--legacy-lineages`，Theta 用 `--theta-sum` 对上当时的 Python 版本）。
1.0 之前（HISAT2、按簇 EM）需要 LineageQuant 1.0.0：`--aligner hisat2 --em-engine legacy --eff-len off`。

**Q：bowtie2 太慢？**
用 `--aligner lqmap`（比对规则相同，快好几倍），再加 `--em-accel squarem`。

**Q：有一部分片段没被计数？**
比对是端到端的且要求两端成对：一端落在注释转录本之外（未注释的同工型、UTR 延伸）的片段不计数。
请尽量提供完整的转录组（每个等位拷贝一条记录）。

**Q：Theta k-mer 和血统 k-mer 有什么区别？**
`--kmer`（默认 21）区分簇内拷贝、给 EM 初值，两种倍性都用；`--lineage-kmer`（默认 31）区分血统，只异源用。

**Q：同源多倍体能用祖先 k-mer 吗？**
不能。同源亚基因组间无分化，祖先 k-mer 相减会减没。同源只做拷贝级定量（`--ploidy auto`）。

**Q：工具报 not on PATH？**
用 `conda env create -f environment.yml` 装全，或单独 `conda install -c bioconda <工具>`。
