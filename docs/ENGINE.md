# LineageQuant 引擎设计：C++ 重写（lqcore → lq）

2026-10-19 起。目标：资源和用时降一个数量级，结果与 LineageQuant 1.0.0 逐项一致。

## 1. 为什么重写，先改哪里

POJ2878 root（约 3,800 万对 reads，32 线程）1.0.0 实测 329 分钟：

| 步骤 | 用时 | 现在的实现 | 主要浪费 |
|---|---|---|---|
| bowtie2 比对 + 排序建索引 | 48 min | C++ | -k 200 的大 BAM 还要排序、建索引 |
| Θ | 8.5 min | Python | （不改，影响小） |
| k-mer 谱系 | 33 min | Python | 按比对记录算，同一条 read 最多重复 200 次 |
| SNP 谱系 | 15–24 min | Python | 每个转录本单独 fetch |
| 谱系融合 PoE | 78 min | Python | 把两张 1 亿多行的文本表对齐 |
| EM | 137 min | Python + numpy | 读 BAM 建等价类约 58 min；5,000 轮全体陪跑约 78 min |

另外，内存峰值约 500 GB（PBS 统计），要交 fat 队列。

**实测（2026-10-20，POJ root，node43，32 线程，从头跑）**：

| | Python 引擎 | lqcore 0.2.1 |
|---|---|---|
| 总用时 | 6 h 04 m 50 s | 1 h 03 m 49 s |
| 分步 | 比对 98 min、Θ 24.5 min、k-mer 26.7 min、SNP 5 min、PoE 75.6 min、EM 134 min | SNP 导出 1.6 min、Θ 23.9 min、比对+谱系+EM 37.9 min（其中比对+扫描 23.5 min，EM 13.7 min） |
| 单进程峰值内存 | 24.4 GB | 18.7 GB（Θ，Python）；lqcore 本身 5.3 GB |

- EM 对拍：`lqcore em` 吃 Python 那次的 BAM、P 表和 Θ，结果与 Python 的 `EM_Final_Expression.tsv` 逐字节相同。
- 两个引擎之间的差异（基因 r = 0.99984，0.1% 的基因差 >10%）全部来自两处刻意的输入差别：SNP 证据取哪一行；估片段长度用的 reads 顺序。
- Θ（Python，没改）占 C++ 总用时的 37%。0.3.0 起改为 `lqcore theta`：
  - 输出与 `global_theta.py` 逐字节相同（`tests/compare_theta.sh`，36 组配置）；
  - 索引是排序的 2-bit k-mer 数组，并行建好后缓存到磁盘；两个线程解压 R1/R2，其余线程扫描；
  - Python 3.12 起内置 `sum()` 对浮点做 Neumaier 补偿求和，`--py-sum` 跟随路由所用的解释器；
  - 本地实测：单个扫描线程每秒约 3.6 万对；预计 POJ root 上 Θ 从 24 分钟降到 2–3 分钟（瓶颈是 gzip 解压）。

**先后顺序**（每一步都能单独替换进现有流程、单独验收）：

1. `lqcore em`：全局 EM。✓ 本地 36 组配置与 Python 逐字节相同。
2. 谱系证据：k-mer + SNP + PoE 在一趟里算完，不写中间概率表。✓ 对拍结果：
   - k-mer 行、SNP 行与 snipe_kmer、snipe_snp 逐字节相同；
   - 融合表和 EM 结果与 merge_poe + em_global 整条链逐字节相同。
3. `lqcore run`：直接接 bowtie2 的输出流，不排序、不写 BAM。✓
   - `--chunk i/N` 按 reads 分块（两个 FIFO 喂数），part 文件可以相加，`lqcore em --parts` 合并后做一次 EM；
   - 3 块合并与整块跑结果相同，线程数也不影响结果。
4. 自写比对器（精确种子 + 汉明距离的无空位比对，见第 6 节）：未做。它会改变方法本身，要在 12 个物种上重新验证。
   眼下比对阶段的提速靠按 reads 拆分到多个节点（第 3 项）。

## 1b. 0.4.0（2026-10-20）：EM 的 NUMA 布局、EM 加速、自带比对器、Salmon 前端

**EM 的 NUMA 布局：结果逐字节不变。** 同一个 EM 在 node43 上测过 803 秒、823 秒和 1,991 秒。原因：
- 大数组由主线程一次性写满，内存页全落在一个 NUMA 节点上；
- E 步是动态调度，线程每轮读到的都是别处的内存。

现在的做法：
- 数组分配时不初始化，由之后负责计算那一段的线程第一次写入（first-touch），内存页就落在该线程所在的节点；
- 类、目标按条目数切成均衡的连续块，第 p 个线程每轮固定算第 p 块；
- EM 线程逐个钉到作业 CPU 集合里的核上。只在 EM 里绑，比对和扫描阶段不绑：启动时就绑核会把整条流水线压到 1 个核上。

验证：`compare_em.sh` 36 组配置仍与 Python 逐字节相同，线程数不影响结果。

**`--em-accel squarem`（默认关闭）：**
- 方法：SQUAREM 方案 3，按单调性保留外推（同 SQUAREM 软件包的做法）；外推值不低于普通 EM 值的 1e-6 倍，目标永远不会被置零；所有求和按固定分块进行，线程数不影响结果。
- 合成的难例（3,000 个拷贝，10 拷贝簇，reads 大多同时比到全部拷贝）：

| 收敛容差 | 普通 EM | SQUAREM |
|---|---|---|
| 1e-3 | 377 步 | 88 步 |
| 1e-5 | 1,489 步 | 160 步 |
| 1e-7 / 1e-8 | 1e-7：6,610 步、97 秒 | 1e-8：742 步、4.2 秒 |

- 两者收敛到同一个解：最大差 0.01 条 reads。
- POJ 上普通 EM 跑满 5,000 步也没收敛，加速的空间在这里。

**lqmap（`--aligner lqmap`，`lqcore map`）：只做 LQ 需要的事。**
- 比对规则与 bowtie2 的 RSEM 参数相同：
  - 无空位、端到端；
  - 每个 mate 打分错配 ≤ floor(0.1 L)，N 不罚分但计入 NM；
  - 一致配对，片段长度 1–1000，不允许 dovetail。
- 只保留比最好比对多 3 个错配以内的候选：再差的 G4 权重小于 1e-11。
- 索引：转录组全部 21-mer 的位置，按前若干位分桶，桶内按低位排序；建好后可缓存。
- 找候选：
  1. 每个 mate 两条链，每 8 bp 取一个精确种子，得到对角线；
  2. 按汉明距离验证；
  3. 同一转录本上配对；
  4. mate rescue：一个 mate 已落在 best + 3 以内、另一个没找到时，在片段窗口里无空位扫描；
  5. 高分歧片段（没有配对，或最好的配对已用掉一半错配预算）再用带 1 个错配的种子补找。
- 输出：与 bowtie2 相同格式的 SAM 记录，直接送进现有扫描，谱系和 EM 代码不变。
- 验证：`tests/compare_map.sh`，与独立的穷举比对器 `tests/brute_map.py` 对拍。合成的多倍体式转录组：
  - 6–10 拷贝的簇，分歧 0–4%；含相同拷贝、异构体、含 N 的转录本；
  - reads 含测序错误、N、小写，16% 来自 5% 分歧的"未见等位基因"。
  
  结果：1,100 个片段的候选集和 NM 全部相同，1 线程与 4 线程输出相同。
- 敏感度：比 bowtie2 `--sensitive` 更密（bowtie2 是 22-mer、每 12–15 bp 一个种子），所以不会更少。真实数据上与 bowtie2 的差别要在 12 个物种的基准上确认。

**sa（`--aligner sa`）：Salmon 的 selective alignment 只负责提出候选。**
- lqcore 从 SAM 里取 (转录本, 位置, 链)，用 lqmap 同一套规则无空位重新验证、重新配对、做 best + 3 过滤；
- 需要空位才成立的候选会被丢弃，模型与 bowtie2 版保持一致；
- 本地只编译通过，没有实跑：沙盒磁盘满，假 salmon 的对拍没跑成；
- 日志里逐项计数：salmon 记录数、缺序列的片段、未知参考。

**拆分到多个节点：** `--chunk i/N` + `lqcore em --parts`，0.2 起就有。开发时的集群脚本（不在发布包里）一次按依赖提交 prep → N 块 → merge，并测墙钟。


## 1c. 1.1.0（2026-10-26）：只剩一个 C++ 程序；二进制输出；有效长度、GC 偏好、marker 清洗

- **整个软件是一个程序** `lq`（= `lineagequant`）：1.0.0 的 Python 路由（cli / config / router）移植成 `run`（src/pipeline.cpp），
  建库的 Python 模块移植成 `normalize-clusters`、`extract-transcripts`、`build-markers`、`snpdb-msa`、`snpdb-vcf`
  （src/prep.cpp），聚合器移植成 `matrix`（src/tools.cpp）。外部工具（bowtie2 / salmon、KMC、gffread、MAFFT、samtools）照旧以子进程调用
  （src/proc.hpp：fork/exec，不经 shell；MAFFT 120 s 超时）。`lqcore` 是同一个程序：程序名以 lqcore 开头时 `run` = 单样本
  （= `quant`），其余子命令两种名字下相同。build.sh 装两份拷贝而不是链接——很多脚本会 `readlink -f`，链接会把名字解析掉。
- **与 Python 一致的细节**（移植时逐条照搬，tests/compare_prep.sh、compare_pipeline.sh 逐字节验证）：
  - Python 的 `str.strip()` 空白集合（含 \x1c–\x1f）、`int()`（空白、符号、数字间下划线）、文本模式的通用换行（\r\n、\r → \n）；
  - Biopython FASTA 规则（id = 标题第一个词，序列行右去空白再删空格与 \r）；
  - `sorted()` 的稳定性与 `reverse=True` 时相等元素保持原序；dict 插入顺序（覆盖不改位置）；
  - pandas：`groupby().sum()` 是逐组 Kahan 求和、组按键排序；`Series.sum()` 是 numpy 的分块成对求和（每 8192 个一块，块内成对、块间左折叠）；
    `round(2)` = rint(x·100)/100；`to_csv` 的浮点 = Python repr、NaN 写空串、csv 最小引用；
  - htslib fai_build_core 的建索引规则（行长一致性、名字取到第一个空白、偏移）；cyvcf2 的整型 FORMAT 哨兵值（缺失 INT32_MIN、补位 INT32_MIN+1）；
  - normalize_clusters.py 的并查集分量顺序（按分量里最早插入的节点）、`--tsv-id-col` 按"第 j 个请求列"过滤的原始行为。
- **二进制输出**：Theta（LQTHETA1）与每样本结果（LQRES001：EM + 有效长度 + idxstats）存的是文本打印前的数；读二进制时按文本的精度
  取整（%.6f / %.2f / %.4f），所以读文本还是读二进制结果相同。`export` 用与直接写文本同一段代码打印（binio.hpp），逐字节一致。
- **MSA SNP 库**：结果按簇顺序合并（1.0.0 按进程完成顺序），血统列表取 ancestry 的全部血统（1.0.0 取库里第一个 SNP 的血统，见 CHANGELOG），
  `--legacy-lineages` 复现 1.0.0。
- 1.0.0 的回归测试在新格式下照常通过（compare_em 36 组、compare_lineage 整条链与分块、compare_theta 36 组 + 二进制、compare_lrcount 80 组）。
- `--eff-len fld`：S0[x] / S1[x] = 观测模板长度（前 100 万个一致配对主比对，0 < TLEN < 2000）的累计数 / 长度加权累计数；
  长度 L 的有效长度 = L − S1[min(L,1999)] / S0[min(L,1999)] + 1（下限 --min-eff-len）。`compare_new.sh` 用 SAM 直接算，逐目标一致。
- `--gc-bias`（Salmon 的片段 GC 模型，单一模型）：
  - 扫描时对每个比上的片段取 mate1 的主一致配对记录，区间 [min(POS, PNEXT) − 1, + |TLEN|) 在参考上的 G/C 数（1 bit / 碱基、popcount），
    按 GC 百分比计 101 格 -> Evidence.gc_fg（part 格式 02 末尾）；
  - 第一遍 EM 收敛后：每个 alpha > 0 的目标取最长参考序列，片段长取模板长度的 21 个分位数，滑窗统计全部可能片段的 GC 分布 h_t；
    背景 = Σ alpha_t · h_t / |h_t|；w = (前景 + 1) / (背景缩放到前景总数 + 1)，限制在 [0.01, 100]，背景上均值归一；
  - 有效长度 × Σ h_t·w / Σ h_t，EM 从第一遍的解继续到收敛；所有求和顺序固定，线程数不影响结果。
  - 合成数据（tests/make_gc_case.py，360 条转录本、6 万片段）：强偏好 MARD 0.2055 -> 0.1100，弱 0.0986 -> 0.0756，无偏好 0.0752 -> 0.0751。
- `--marker-clean`：建好 marker 索引后，扫描带谱系标签的参考序列的全部 k-mer（规范 / 双链与查询时相同），记下每个 marker k-mer 出现在哪些
  谱系的参考拷贝里（allowed 位掩码，原子 OR）；allowed ≠ 0 时 mask &= allowed，mask 变 0 的删掉，重建分桶。动机：犬蔷薇留出 v3，
  祖先 marker 来自加了 0.5–2% 替换的序列时，亚基因组共有位点上的祖先突变变成另一祖先的"特异" k-mer，来源在参考里的 read 有 1–2% 被判错谱系。

## 2. 语言：C++17 + OpenMP + zlib

- 速度：C++ 和 Rust 都编译成本地机器码，同样的数据结构速度相同，差别在写法，不在语言。
- 选 C++ 的理由：
  - bowtie2 是 C++，以后要嵌入或改写它的代码，同语言最省事；
  - 集群 GCC + OpenMP 现成，不用额外装工具链；
  - 依赖只有 zlib。libstdc++ 静态链接，产物是单个可执行文件。
- 不依赖 htslib：输入是 SAM 文本，生产环境直接来自 bowtie2 的输出流；对拍或读已有 BAM 时用 `samtools view` 管道。

## 3. 确定性与等价性（硬性要求）

- **与 Python 结果一致**：Python 引擎里每一步浮点运算的顺序都照搬：
  - numpy `bincount` 是按条目顺序依次累加；
  - numpy `sum` 是 8 路成对求和（`np_pairwise`），超过 8192 个元素时按 8192 一块分块（`np_add_reduce`）；
  - Python 的 `round(x, 4)`；
  - 等价类按首次出现的顺序编号，先无权重类、后带权重类。
  - 目前唯一可能不一致的地方是 `pow`：numpy 在 AVX-512 机器上可能用 SVML 的实现，最后一位可能不同。
- **与线程数无关**：EM 每轮分两步，每一步的写入对象都只归一个线程：
  - A 步：每个等价类由一个线程算 E 步权重；
  - B 步：每个目标按条目顺序累加。
  
  所以 1 个线程和 32 个线程逐字节相同。
- **对拍**：
  - `tests/compare_em.sh`：合成数据，12 种配置 × 3 个种子，再加强制 3,000–4,000 轮的长跑，全部逐字节相同；
  - POJ 实数据：开发时用集群脚本和已有的正式结果比（不在发布包里）。

## 4. 数据结构

- **等价类**：用 CSR 存——每类一个起始偏移；每条条目存目标（u32）和权重（f64）；另有每类计数（f64）。再建一份"目标 → 条目"的转置索引（同样是 CSR，条目号升序），B 步就不用加锁、也不用原子操作。
- **扫描阶段**：每个片段只保留一个小向量（目标，两个 mate 的最小 NM），等价类键是紧凑的字节串。
- **相同序列合并**：判断序列是否完全相同，用 128 位序列哈希代替 MD5，只用于分组，不改变分组结果。
- **以后（里程碑 2–3）**：
  - k-mer 用 2-bit 编码存成 u64，查表用开放寻址哈希或排序数组；每个片段只查一次，不按比对记录查。
  - SNP 每个转录本存一个排序的位置数组和 uint8 期望碱基，就是 snipe_snp 里已经做过的那套紧凑索引。

## 5. 并行与 NUMA（借鉴 KmerGenoPhaser 的 comput_opt）

- **first-touch 分区**：所有大数组用和计算时相同的静态分区来初始化，内存页自然分散到各个 NUMA 域。
- **线程绑定：不用**。设了 `OMP_PROC_BIND` / `OMP_PLACES`，libgomp 在程序启动时就把主线程钉在 1 个核上，扫描线程（std::thread）和 popen 出来的 bowtie2 全部继承。POJ root 实测：align+scan 10.8 小时，约 1.4 万条记录/秒，本该约 1 小时。
  - comput_opt 的绑核适用于纯 OpenMP 程序；lqcore 混用了 OpenMP、std::thread 和子进程，所以不适用。
  - 0.2.1 起，lqcore 检测到这些变量会先放开 CPU 掩码、删掉变量、重新 exec 自己；设 `LQCORE_KEEP_BINDING=1` 可以保留绑核。
  - `run` 开头打印本进程可用的 CPU 数，比请求的线程数少时报警。
- **数据并行**：
  - 等价类计数可以直接相加。所以拆 reads 并行时，各块各自建类，最后合并、做一次 EM。这相当于只做一次 all-reduce，不像 SGD 每一步都要同步。
  - EM 按连通分量天然可拆，可以选择让已收敛的分量先停。这个选项会改变结果（差异在 tol 量级），默认关闭。
- **参数**：线程数等全部走 CLI，不写死。
- **临时文件**：用 `PBS_JOBID` 加进程号保证唯一，即 comput_opt 文件 rendezvous 的思路。
- **只读大索引**（marker k-mer、比对索引）以后放 `/dev/shm` 用 mmap 共享：同一节点上的多个组织作业只占一份内存。

## 6. 比对阶段的思路（里程碑 4）

- 按参考坐标或按簇切分参考，不减少总工作量：多重比对每条 read 仍要和全部分片比。还会丢跨簇的候选，POJ 里有 6.6–7.1% 的片段跨簇。
- **该拆的是 reads**：FASTQ 分块，多节点各跑 bowtie2 + lqcore 的扫描，产出等价类计数，然后合并。
- **自己建索引**：我们的 bowtie2 参数本来就不开空位（`--gbar 99999999`），比对就是固定偏移上的错配计数。
  - G4 下只需要找全"比最好比对多 2 个错配以内"的候选，因为多 3 个错配时权重已经乘到约 3×10⁻⁹。
  - 鸽巢原理：150 bp 切成 5 段 × 30 bp，错配 ≤ 4 的位置至少有一段完全一致。所以用 30-mer 精确查表找候选不会漏，再用 2-bit 序列做异或 + popcount 算错配数。
  - 错配多的少数 reads 用短种子兜底。
  - 必须先在 12 个物种的基准上证明候选不比 bowtie2 少，才能作为默认。

## 7. 已知要在流式版本里定下来的问题

- 现在 EM 读 `P_values_merged.tsv` 时，每个 read 用的是"文件里第一条非零权重的行"。
- 文件顺序由 `imap_unordered` 决定，跑两次可能不同；SNP 证据又是按基因算的，所以选中的行不确定。
- 流式版本要改成确定的规则，例如取最佳比对目标那一行，或者把各候选目标的 SNP 证据合起来。
- 改了规则会改变结果，要单独评估，并写进 CHANGELOG。

## 8. 目录

```
src/main.cpp        子命令入口（程序名决定 run 的含义）
src/pipeline.cpp    run：自动路由
src/prep.cpp/.hpp   建库各步；proc.hpp 外部工具
src/theta.cpp       theta
src/run.cpp         quant（单样本一趟）、markers
src/scan.cpp/.hpp   比对记录 -> 候选与等价类
src/lineage.cpp/.hpp  k-mer / SNP 血统、PoE、marker 清洗
src/model.cpp/.hpp  模型、global EM、有效长度、GC、part 文件
src/em.cpp          em
src/lqmap.cpp/.hpp  自带比对器
src/longread.cpp    lrpick / lrcount
src/tools.cpp       matrix、export
src/binio.hpp       .lqb / .lqres 与文本写出
src/common.hpp      日志、格式化、读行（plain/gz/管道）、numpy / Python 兼容的数值细节
tests/              合成数据 + 与 1.0.0 Python 对拍（run_all.sh）
build.sh / Makefile 编译
```
