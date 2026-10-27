# LineageQuant 整合架构设计 (v3 — Auto-Routing Edition)

> **1.1.0 说明（2026-10-26）**：LineageQuant 1.1.0 是一个 C++ 程序（`lineagequant`，源码 `src/`），路由逻辑与本文一致，
> 每一步换成同名子命令：Theta = `theta`，Step 00/02/02b/03 = `quant`（一趟、不写 BAM 与逐 read 概率表），
> Step 00b = `lrcount`，SNP 库 = `snpdb-msa` / `snpdb-vcf`，marker = `extract-transcripts` + `build-markers`，
> 等位簇标准化 = `normalize-clusters`；Step 04 聚合改为流程外的 `matrix`。下文提到的 Python 模块
> （`router.py`、`em_global.py`、`snipe_kmer.py`、`genome_pair.py`、`build_markers.py` 等）只在 LineageQuant 1.0.0 里，
> 1.1.0 的对应子命令在默认参数下与它们逐字节一致（`tests/run_all.sh`）；引擎设计见 `docs/ENGINE.md`，
> 文件格式见 `docs/FORMATS.md`。HISAT2、按簇 EM（legacy）与 `--auto-idxstats` 只在 1.0.0 里有。
>
> **1.0.0 说明（2026-10）**：本文是 v3 自动路由架构的设计文档，路由逻辑沿用至今。1.0.0 在此基础上改了两处定量核心，
> 详见仓库根目录 `CHANGELOG.md`：
> (1) Step 00 默认改用 bowtie2（RSEM 参数，`-k 200`，参考序列中的 N 不罚分），HISAT2 保留为 `--aligner hisat2`；
> (2) Step 03 默认改用 global EM 引擎（`core/em_global.py`）：全部基因一起做 EM、有效长度校正、按错配数加权、
> 序列相同的基因合并。原按簇 EM 保留为 `--em-engine legacy`。下文涉及 HISAT2 与按簇 EM 的描述按 1.0.0 理解即可。


> 把早期的短读长/异源流程 (v7) + 长读长/hybrid 流程 (v8) + 血缘特异
> k-mer 构建 (KMC 复合相减) 整合成一个参数化、自动路由的完整软件。
>
> 用词约定: homolog(同源基因) 与 allele(等位基因) 在多倍体中边界不清晰,
> 本文档及代码保守统一用中性的"等位簇 / allele cluster / cluster",
> 不武断断言"同源"或"等位"。亚基因组间对应拷贝统称"等位拷贝(allele copy)"。

---

## 一、核心理念：输入驱动的自动路由

用户**不手动选模块**。只声明物种类型 (同源/异源)，软件检查**给了哪些输入**，
自动决定走哪条流程。允许空参数，缺啥跳过对应模块，所有路最后汇到 EM → TPM。

---

## 二、同源 vs 异源：根本性分叉

这是整个路由系统的第一级决策，决定了管线是否需要"血统分配"。

### 同源多倍体 (auto)

亚基因组间无分化，不存在"血统 (lineage)"概念。管线只做一件事：**拷贝级表达定量**——
把 reads 分配给等位簇内的各个拷贝，输出每个拷贝的表达量。不涉及 KMC 祖先相减、
不涉及 MSA、不涉及群体 VCF、不涉及双专家融合。本质上是"多倍体版 Kallisto"。

核心思路：利用等位拷贝间的 CDS 序列差异（虽然是微小的）做 reads 分配。

### 异源多倍体 (allo)

亚基因组间有明确血统分化（如甘蔗 SOF/SRO/SSP）。管线做两件事：
1. **血统级分配**：每条 read 来自哪个血统？
2. **拷贝级定量**：在同一血统内，read 来自哪个具体拷贝？

血统分配是异源独有的维度。拷贝级定量在血统分配之后进行。

---

## 三、两套独立的 K-mer 系统

这是最容易被混淆的概念。LQ 中使用了两套完全不同的 K-mer，
分别服务于不同层级，提取方式、作用对象、使用阶段都不同。

| | 系统 1: Theta K-mer | 系统 2: 血统 K-mer |
|---|---|---|
| **提取方式** | 簇内收割唯一 K-mer（仅此拷贝拥有） | KMC N 元复合相减（祖先种间差异） |
| **输入** | mRNA.fasta（所有等位拷贝的 CDS） | 祖先基因组 FASTA + GFF |
| **作用对象** | 单个等位拷贝 | 整个血统 |
| **服务阶段** | Step 01（Theta 先验） | Step 02 K-mer 路线（血统概率） |
| **输出** | `Initial_Theta`（每个拷贝的初始表达量猜测） | `P(SOF), P(SRO), P(SSP)`（每条 read 的血统概率） |
| **同源/异源** | 两者都用 | 仅异源多倍体 |
| **可被替代** | 长读长 idxstats（01_hybrid） | MSA SNP 或群体 VCF SNP |

**关键认知**：即使异源多倍体走了双专家融合，"Theta K-mer"和"血统 K-mer"
也仍然是两套独立的系统。融合发生在血统概率层面（系统 2 + SNP），
Theta 先验（系统 1 或其长读长替代）以独立维度进入 EM。

---

## 四、血统界定路线（仅异源多倍体）

异源多倍体有三种方式获得 read 的血统归属概率。

### 路线 A：祖先 K-mer（KMC 复合相减）

**触发条件**：`--ancestor NAME:fa,gff [...]` 存在

**流程**：
1. `build_markers.py`：对祖先种基因组做 KMC N 元复合相减
   - K_A ∩ K_B ∩ ... → 共享（无区分力，丢弃）
   - K_A \ (K_B ∪ K_C ...) → A 血统特异 K-mer
   - 对每个祖先重复，输出各血统特异 K-mer FASTA
2. `snipe_kmer.py`：对每条 read 扫描血统特异 K-mer
   - 命中血统 A 的 K-mer 数 / 总命中数 → P(A | read)
   - 无命中 → status = NoSignal

**适用性**：异源多倍体的祖先种间序列分化足够大（通常 >2%），KMC 相减能产生足够的血统特异 K-mer。
同源多倍体无法使用该路线（祖先不存在或序列分化不够）。

### 路线 B：群体 VCF SNP

**触发条件**：`--vcf <vcf> --sample-lists NAME:list [...]` 存在

**流程**：
1. `build_snp_vcf.py`：从群体 VCF 构建 SNP 数据库
   - 对每个 SNP 位点，按血统分组计算等位基因频率 (AF)
   - 通过 GFF3 外显子结构，将基因组坐标投射到转录本 CDS 坐标
   - 输出 `{(transcript_id, CDS_pos): {SOF_af, SRO_af, SSP_af}}`
2. `snipe_snp.py`：对每条 read 覆盖的 SNP，计算基因型似然
   - `P(read_base | lineage) = af × (1-ε) + (1-af) × (ε/3)`
   - 联合所有 SNP 的对数似然，归一化得到三血统概率

**特点**：AF 是经验值（0-1 连续），产出的概率向量有梯度，而非一刀切。

**适用性**：有群体重测序数据 + 明确的血统分组样本。不需要祖先基因组。

### 路线 C：MSA SNP（Phased 转录组固定差异）

**触发条件**：`--phased-fasta <fa>` 存在

**流程**：
1. `build_snp_msa.py`：对每个等位簇的 phased CDS 做 MAFFT 多序列比对
   - 扫描每个比对列：若不同血统间存在碱基差异，记录为 SNP
   - 输出 `bases: {SOF: 'A', SRO: 'C', SSP: 'G'}` + 确定性 AF (0/1)
2. `snipe_snp_hybrid.py`：对每条 read 覆盖的 SNP，逐碱基匹配
   - `read_base == expected_base → P = 1-ε; 否则 → P = ε/3`
   - 处理大小写匹配和三等位位点
   - 联合似然，归一化得到三血统概率

**特点**：基于固定差异，产出的概率向量更尖锐（接近 0/1）。确定性 AF
会导致 Chimeric_Drop 问题——见第五节长读长校准。

**适用性**：有 phased 转录组数据。同源多倍体无法使用（血统间无固定差异），
异源多倍体（包括甘蔗 POJ 这种三血缘分化足够大的情况）可以有效使用。

### 三条路线的互斥与叠加

- **MSA 和 VCF 互斥**：二者都提供 SNP 信号，由输入决定走哪条
  （`--phased-fasta` vs `--vcf`），不会同时运行
- **SNP 和祖先 K-mer 可叠加**：任何一条 SNP 路线都可以与祖先 K-mer 路线共存，
  共存时触发双专家 PoE 融合（见第六节）

---

## 五、长读长 = 锚定增强（非独立路线）

长读长不产生新的血统信号路线，而是在两个位置增强现有路线。

### 位置 1：Theta 先验增强（01_hybrid 替代 01）

- 用 PacBio BAM 的 `idxstats` 物理计数作为每个等位拷贝的初始 Theta
- 替代短读长路线中基于 Theta K-mer 探针收割的统计推断
- 物理先验比统计先验更可靠，尤其在低探针基因上
- 无长读长 → 退回 Theta K-mer 统计推断（01）

### 位置 2：MSA SNP 校准（02c，仅 MSA 路线）

- MSA 路线的 AF 是确定性的（0 或 1）
- 问题：如果 SSP 拷贝在某个位点也偶尔出现 SOF 的预期碱基，
  read 会被错误地完全判归 SOF → Chimeric_Drop
- `02c_calibrate_snp_db.py` 用长读长 pileup 重新估计 AF：
  在 SNP 位点上统计实际碱基频率，用经验值替换确定性的 0/1
- 校准后 SNP 数据库的 AF 变成连续值，减少 Chimeric_Drop
- 仅对 MSA 路线有效：VCF 路线的 AF 已经是经验值，不需要校准
- 无长读长 → 跳过校准，直接用确定性的 0/1 bases

---

## 六、双专家 PoE 融合（02b，条件触发）

**触发条件**：祖先 K-mer 信号（路线 A）和 SNP 信号（路线 B 或 C）**同时产出成功**。

**融合器 `02b_merge_lineage_probs.py` 的核心设计**：

1. **按 read_name 对齐**：两路输出都以 read_name 为 key，各含三血统概率向量
2. **计算每条 read 在每路中的信息权重**：

   ```
   base_weight:
     K-mer 路线: Specific=10, Conflict=2.5, NoSignal=0.5
     SNP 路线:   Specific_Soft=10, Multi=2.5, 其他=0.5

   continuous = max(0, log₂(N_lineages) - H(P)) × 5.0
   n_info = base + continuous
   ```

   熵越低（概率越集中）→ continuous 越高 → 该路线的权重越大。
   例如 `(0.9, 0.05, 0.05)` 的连续权重远高于 `(0.4, 0.3, 0.3)`。

3. **动态归一化 + 最小权重约束**：

   ```
   alpha_kmer = w_kmer / (w_kmer + w_snp)
   alpha_snp = 1 - alpha_kmer
   alpha = clamp(alpha, 0.05, 0.95)
   ```

   min_weight=0.05 确保弱信号不会被完全忽略。

4. **Product of Experts 融合**：

   ```
   P_merged(lineage) ∝ (P_kmer(lineage) + ε)^α_kmer × (P_snp(lineage) + ε)^α_snp
   然后归一化
   ```

   PoE 的三种行为模式：
   - 一方强一方弱 → 弱信号指数 ≈ 0.05，`x^0.05 ≈ 1`（近似常数），不稀释强信号
   - 双方一致指向同一血统 → 乘积互相放大，合并结果更尖锐
   - 双方分歧 → 乘积同时压制双方高值，产生保守的中间分布

5. **仅一路有该 read**：直接使用该路概率，不做融合（实践中大量
   reads 只覆盖其中一路）。

**关键设计性质**：融合器**不感知 SNP 信号的来源**（MSA 还是 VCF）。
它只消费概率向量 `(P_SOF, P_SRO, P_SSP)`。上游是 MSA 还是 VCF
只影响概率向量的尖锐度（MSA 更 one-hot，VCF 更梯度化），
不影响融合数学。

---

## 七、汇合点（所有路线共用）

```
Theta先验 (01 或 01_hybrid) → 同源: 拷贝级先验
                              异源: 拷贝级先验 × 血统归属
                                      │
                                      ▼  03_EM_dynamic
                                EM 迭代分配 reads 到等位拷贝
                                      │
                               EM_Final_Expression.tsv
                                      │
                                      ▼  04_matrix_aggregator
                                孤儿基因: idxstats 直接计数
                                多倍体基因: EM 分配 reads
                                      │
                           LineageQuant_Final_Merged_TPM.tsv
```

### 同源多倍体的 EM

输入：仅 Theta 先验（拷贝级），无血统概率。

EM 直接将 reads 分配给簇内的各个等位拷贝。分配依据是 reads 和拷贝 CDS 序列的
匹配程度（类似 Kallisto 的 likelihood 计算）。输出每个拷贝的 Final_Theta。

### 异源多倍体的 EM

输入：Theta 先验（拷贝级）+ 血统概率矩阵（read × lineage）。

EM 首先将每条 read 按血统概率分配到对应血统的拷贝子集，
再在子集内按拷贝级先验做分配。输出每个拷贝的 Final_Theta
和 EM_Allocated_Reads。

---

## 八、完整路由表（10 条有效路径）

### 同源多倍体 (auto)

| # | 短读长 | 长读长 | 路线 |
|---|--------|--------|------|
| A1 | ✓ | | Theta K-mer (01) → EM |
| A2 | ✓ | ✓ | 长读长 Theta (01_hybrid) → EM |

### 异源多倍体 (allo)

| # | 短读长 | 长读长 | 祖先 | VCF | Phased | 路线 |
|---|--------|--------|------|-----|--------|------|
| B1 | ✓ | | ✓ | | | Theta K-mer (01) + 祖先 K-mer 血统 → EM |
| B2 | ✓ | | | ✓ | | Theta K-mer (01) + VCF SNP 血统 → EM |
| B3 | ✓ | | ✓ | ✓ | | Theta K-mer (01) + 祖先 K-mer + VCF SNP → **PoE 融合** → EM |
| B4 | ✓ | ✓ | ✓ | | | 长读长 Theta (01_hybrid) + 祖先 K-mer 血统 → EM |
| B5 | ✓ | ✓ | | | ✓ | 长读长 Theta (01_hybrid) + MSA SNP (02c 校准) → EM |
| B6 | ✓ | ✓ | ✓ | | ✓ | 长读长 Theta (01_hybrid) + 祖先 K-mer + MSA SNP (02c 校准) → **PoE 融合** → EM |
| B7 | ✓ | ✓ | | ✓ | | 长读长 Theta (01_hybrid) + VCF SNP 血统 → EM |
| B8 | ✓ | ✓ | ✓ | ✓ | | 长读长 Theta (01_hybrid) + 祖先 K-mer + VCF SNP → **PoE 融合** → EM |

### 路径注释

- **B1**：最基本的异源路线。仅靠祖先 K-mer 做血统分配，无长读长增强。
- **B3, B6, B8**：三条双专家融合路线。B3 是 VCF+K-mer（短读长），B6 是 MSA+K-mer（长读长增强），B8 是 VCF+K-mer（长读长增强）。融合器对三者使用完全相同的 PoE 算法。
- **B5**：甘蔗 POJ2878 实际跑通的路径。长读长提供 Theta 锚定 + MSA SNP 校准，单 SNP 路线不融合 K-mer。
- **B6**：POJ 的潜在优化方向。在 B5 基础上启用祖先 K-mer，实现 MSA SNP + 祖先 K-mer 双专家融合。
- **B4**：有长读长 Theta 锚定和祖先 K-mer，但没有 MSA 或 VCF。适合只有祖先基因组而无 phased/群体数据的场景。

---

## 九、自动路由决策树

```
输入检查:
  必需: --genome <fa> --gff <gff> --reads R1,R2
  可选: --ancestor NAME:fa,gff [...]   (祖先基因组, KMC 相减)
        --vcf <vcf> --sample-lists NAME:list [...]   (群体变异)
        --phased-fasta <fa>   (phased 转录组, MSA)
        --longread <bam>   (长读长)
  声明: --ploidy {auto|allo}

═══ 决策 ═══

STEP 1 [等位簇]:
  if 已提供 --clusters: 直接用
  elif 提供参考基因组 + GFF: jcvi 共线性建簇 → clusters/ancestry/orphan
  else: 报错 (无法定义等位簇)

STEP 2 [血统界定] — 仅异源多倍体; 同源跳过此步:
  kmer_sig = snp_sig = None

  # 路线 A: 祖先 K-mer
  if --ploidy == allo and --ancestor 存在:
      build_markers(KMC N 元相减) → marker fastas
      kmer_sig = snipe_kmer(bam, markers)

  # 路线 B 或 C: SNP (MSA 与 VCF 互斥)
  if --ploidy == allo and --vcf 存在:
      build_snp_vcf(vcf, gff, sample-lists) → snp_db
      snp_sig = snipe_snp(bam, snp_db)
  elif --ploidy == allo and --phased-fasta 存在:
      build_snp_msa(phased, clusters, ancestry) → snp_db
      if --longread: calibrate_snp(snp_db, longread) → snp_db   # 02c
      snp_sig = snipe_snp_hybrid(bam, snp_db)

  # 双专家融合 (条件触发)
  if kmer_sig and snp_sig:
      lineage_prob = merge_poe(kmer_sig, snp_sig)    # 02b
  else:
      lineage_prob = kmer_sig or snp_sig

STEP 3 [Theta 先验]:
  if --longread:
      lr_prior = extract_prior(longread)             # 00b
      theta = global_theta_hybrid(reads, clusters, lr_prior)   # 01_hybrid
  else:
      theta = global_theta(reads, clusters)          # 01 (Theta K-mer 统计推断)

STEP 4 [汇合]:
  if --ploidy == auto:
      # 同源: 仅拷贝级先验, 无血统概率
      em = em_dynamic(theta, clusters)               # 03 (Kallisto 式)
  else:
      # 异源: 拷贝级先验 + 血统概率
      em = em_dynamic(theta, lineage_prob, ancestry) # 03
  tpm = aggregate(em, orphan, bams)                  # 04
```

---

## 十、工程化层

### 10.1 格式兼容层 (utils/genome_pair.py；1.1.0：`lineagequant extract-transcripts`，`src/prep.cpp`)

解决"同一基因组和注释要能互通读写"（如参考基因组 fasta 用 Chr01，GFF 用 CM001740）。

```python
class GenomeAnnotationPair:
    def __init__(fasta, gff, min_chr_len=1_000_000):
        # 1. fasta 头清洗: 切 TAB/空格后缀 (GWH 格式)
        # 2. 建 fai, 选 ≥min_chr_len 的序列名 (排除 scaffold/contig/organelle)
        # 3. 自动建 seqid 映射: fasta 序列名 ↔ gff 第 1 列
        #    - 直接名匹配 优先
        #    - 长度匹配 + 模糊名 兜底 (Chr01 ↔ CM001740)
        # 4. gff 过滤: 只留染色体 seqid + strand≠? (排除 organelle trans-splicing)
    def extract_transcripts(out_fa, feature='mRNA'):
        # gffread 从清洗后的 (fa,gff) 提转录本, 统一接口
    def seqid_map():  # 返回 fasta↔gff 名映射表
```

不管原始命名，绑定 (fa, gff) 一组，内部统一读写。

### 10.2 血统名彻底参数化

现有脚本中 `LINEAGE_NAMES=["SOF","SRO","SSP"]` 是硬编码的。整合后：
- 血统名从 **Ancestry 文件第 2 列**动态读取（去重排序），不限数量
- `build_snp_vcf` / `snipe_snp` / `merge_poe` / `calibrate` 全部改成接受 `lineage_names` 列表参数
- SNP entry 的 `bases` / AF 字段用血统名做 key，通用化
  （已部分支持，如 `build_perfect_snp_db.py` 的 `bases` dict）

---

## 十一、目录结构

1.0.0（Python 包）的目录结构见 LineageQuant 1.0.0。1.1.0：

```
LineageQuant_v1.1.0/
├── build.sh                 # 编译 -> bin/lineagequant(.<时间>) + bin/lqcore(.<时间>)，原子切换、保留旧版本
├── Makefile                 # make / make install PREFIX=...
├── src/
│   ├── main.cpp             # 子命令分派（程序名为 lqcore 时 run = 单样本）
│   ├── pipeline.cpp         # ★ run：自动路由（1.0.0 的 cli.py + config.py + router.py）
│   ├── prep.cpp/.hpp        # 建库：normalize-clusters / extract-transcripts / build-markers / snpdb-msa / snpdb-vcf
│   ├── proc.hpp             # 调外部工具（fork/exec、超时、PATH 查找）
│   ├── theta.cpp            # 01 / 01_hybrid Theta（global_theta.py）
│   ├── run.cpp              # quant：00 比对流 + 02/02b 血统证据 + 03 EM（一趟）；markers
│   ├── scan.cpp/.hpp        # 比对记录 -> 候选、错配数、等价类、血统因子
│   ├── lineage.cpp/.hpp     # k-mer 血统（snipe_kmer.py）、SNP 血统（snipe_snp.py）、PoE（merge_poe.py）、marker 清洗
│   ├── model.cpp/.hpp       # 目标 / 簇 / 合并；global EM（em_global.py）；有效长度、GC 偏好；part 文件
│   ├── em.cpp               # em：SAM/BAM 或 parts -> EM
│   ├── lqmap.cpp/.hpp       # 自带比对器
│   ├── longread.cpp         # lrpick / lrcount（extract_prior.py）
│   ├── tools.cpp            # matrix（aggregator.py）、export
│   ├── binio.hpp            # 二进制输出 .lqb / .lqres 与文本写出
│   ├── common.hpp sam.hpp   # 日志、格式、与 Python/numpy 一致的数值细节
├── tests/                   # 与 1.0.0 Python 逐字节对拍（run_all.sh）；外部工具与 Python 依赖的替身
├── docs/
│   ├── DESIGN.md            # 本文档（路由设计）
│   ├── ENGINE.md            # 引擎设计（C++ 重写的取舍、确定性、并行）
│   ├── FORMATS.md           # 二进制与文本格式
│   └── USAGE_zh.md          # 中文使用说明
├── examples/                # 各路线的示例脚本
├── environment.yml          # 外部工具 + 编译器
├── README.md  CHANGELOG.md  LICENSE
```

---

## 十二、整合原则

1. **公共脚本以 v8 为准**（v8 bin 覆盖 v7 bin：2 个完全相同，3 个 v8 是超集）。
2. **所有路径参数化**，绝不留绝对路径/账号信息。
3. **断点续跑**沿用（`.done` marker）。
4. **依赖检查**沿用 v8 的 preflight（按路线检查 bowtie2 / salmon、samtools、mafft、kmc、gffread；1.1.0 不再有 Python 模块）。
5. 保留现有仓库的入口名 `lineagequant` 和参数风格 `NAME:PATH`。
6. **同源和异源的 EM 共用同一套代码**，差异仅在于是否传入血统概率参数。
7. **融合器不感知 SNP 来源**——MSA 和 VCF 的融合算法完全一致。

---

## 附录：v2 → v3 修订列表

| 修订项 | v2 (旧) | v3 (新) |
|--------|---------|---------|
| 同源多倍体路线 | MSA SNP → EM（方向性错误） | Theta K-mer → EM（Kallisto 式拷贝级定量） |
| MSA 适用场景 | 标注为"同源多倍体"专属 | 改为"有 phased 转录组"，同源/异源均可，异源需血统分化足够 |
| 决策树中 MSA 位置 | `elif --phased-fasta`（与 VCF 共用 elif，暗示与 KMC 互斥） | 保持与 VCF 互斥，但从 VCF 的 `elif` 中独立出来，允许与 KMC 叠加 |
| 两套 K-mer | 混称为"k-mer 信号" | 显式区分 Theta K-mer（系统 1）和血统 K-mer（系统 2） |
| 路由表 | 6 行，缺失多条路径 | 10 行，覆盖所有有效输入组合 |
| 同源 EM | 未描述 | 明确为拷贝级定量，无血统概率输入 |
| 融合路线数量 | 仅列出 VCF + K-mer | 列出 VCF+K-mer、MSA+K-mer 两条，明确融合器不感知 SNP 来源 |
| 示例脚本 | 3 个 | 7 个（补充同源、MSA+K-mer 融合、完整 POJ 版） |


---

## 附录：v3 → 1.0.0 修订列表

| 修订项 | v3 | 1.0.0 |
|--------|----|-------|
| 短读长比对 | HISAT2 -k 15 | bowtie2（RSEM 参数，-k 200，N 不罚分）；HISAT2 可选 |
| EM | 按簇 EM，孤儿基因用 idxstats 计数 | global EM：全部基因一起，每个片段总权重 1 |
| 有效长度 | 无 | 默认开（θ / (L − 平均片段长 + 1)） |
| 候选比对权重 | 一视同仁 | 按错配数加权（λ = (p/3)/(1−p)，p 由数据估计） |
| 序列完全相同的拷贝 | 各自估计（易出现 0/1 分法） | 合并估计、组内平分 |
| auto 路线 | 比对 → idxstats → TPM | 比对 → Theta → EM → TPM（旧行为保留为 --auto-idxstats 对照） |
| 比对器混用保护 | 无 | 每个样本目录记录 aligner.txt，换比对器时拒绝续跑 |

---

## 附录：1.0.0 → 1.1.0 修订列表

| 修订项 | 1.0.0 | 1.1.0 |
|--------|-------|-------|
| 实现 | Python 包 + 可选 C++ 引擎（`--engine cpp`） | 一个 C++ 程序；`lqcore` 是同一程序（`lqcore run` = 单样本） |
| 每样本输出 | 文本：Theta、EM、idxstats（BAM 或 lqcore 文本） | 二进制 `.lqb` / `.lqres`；`export` 逐字节还原文本 |
| Step 04 聚合 | 流程最后一步 | 流程外的 `matrix`（同一文件），不计入定量耗时 |
| SNP 库 | pickle + 导出文本 | 直接写引擎读的文本；MSA 并行结果按簇顺序合并 |
| MSA SNP 血统列表 | 库中第一个 SNP 的血统（可能漏血统） | ancestry 的全部血统；`--legacy-lineages` 复现 1.0.0 |
| 基因组 fai | samtools faidx | 程序自己写（与 htslib 相同的格式与规则） |
| Step 1b 长读校准 | 运行（对结果无影响） | 不运行 |
| 新选项（默认关） | — | `--eff-len fld`、`--gc-bias`、`--marker-clean`、`matrix --tpm-length eff` |
