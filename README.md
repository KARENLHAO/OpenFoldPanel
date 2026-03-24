# OpenFoldPanel

`OpenFoldPanel` 用来比较多组 AI 生成的蛋白结构模型，并把分散在结构文件里的关键信息整理成统一的二维 flat panel 报告。面对同一目标蛋白的多组候选模型，单看 top-ranked model 往往不够，尤其在多聚体界面、局部低置信度区段和不同预测器结果不一致时，更需要横向比较。

`OpenFoldPanel` 会围绕参考链，把二级结构、可及性、疏水性、分子接触、可用时的置信度，以及可选的 MSA / conservation 信息并排展示，帮助你更快看出哪些区域稳定、哪些区域需要复核。输出既包括浏览友好的 `report.html`，也包括适合归档和汇报的链级 PDF、`tracks.json`、`summary.txt` 和 `logs.txt`。

这个项目受 FoldScript 的比较与展示思路启发，但实现形态是开源、可本地运行、可接入批处理流程的 CLI，而不是 web server。

## 核心能力

- 支持单个结构文件输入，也支持 `.zip`、`.tar.gz` 等压缩包批量输入。
- 支持多模型、多链和多 job 处理，适合把不同预测结果放在同一套视图里比较。
- 生成 `report.html`、`reference-chain-<CHAIN>.pdf`、`tracks.json`、`summary.txt` 和 `logs.txt`。
- 提供二级结构、可及性、疏水性、接触、可选 MSA / conservation 轨道，以及可用时的置信度轨道。
- 可覆盖一级、二级，以及三级/四级相关结构线索的汇总与比较：一级结构以参考链序列、查询序列和可选 MSA 行展示，二级结构以 DSSP 或几何回退得到的 secondary structure 轨道展示，三级/四级则通过可及性、疏水性、二硫键、置信度和跨链/配体/核酸/离子接触提供相互作用摘要。
- 支持 `PDB`、`CIF`、`mmCIF` 以及常见归档格式，便于接入现有预测和整理流程。

## 安装与依赖

建议使用单独的 `conda` 或 `mamba` 环境。这样更容易把 Python 依赖、DSSP 和 MSA 工具放在一套可复现的环境里。

### 1. 创建环境

```bash
mamba create -n openfoldpanel -c conda-forge -c bioconda \
  python=3.10 pip \
  cmake ninja cxx-compiler pkg-config git

conda activate openfoldpanel
```

### 2. 安装项目

```bash
python -m pip install -r requirements.txt
python -m pip install -e .
```

### 3. 可选增强依赖

`OpenFoldPanel` 的核心报告流程只依赖 Python 包就能跑起来，但下面这些工具会明显提升结果完整性：

- `mkdssp` / `dssp`
  用于更可靠的 secondary structure 和 accessibility 计算。缺失时程序会根据原子几何关系粗略估计二级结构和可及性，保证报告还能生成，但精度不如 DSSP。
- `blastp` 或 `mmseqs`
  用于同源序列检索。
- `clustalo`
  用于多序列比对和 conservation 计算。

如果你希望把这些增强能力放进同一个环境，可以继续安装：

```bash
mamba install -n openfoldpanel -c conda-forge -c bioconda \
  blast mmseqs2 clustalo
```

### 4. 安装 DSSP（推荐）

如果只是想先把流程跑通，也可以先用社区包临时安装 DSSP：

```bash
mamba install -n openfoldpanel -c sbl dssp
```

如果你把 DSSP 安装在项目本地目录，例如 `./.local/dssp/bin`，记得把它加入 `PATH`：

```bash
export PATH="$(pwd)/.local/dssp/bin:$PATH"
```

正式环境更推荐使用官方仓库源码编译，路径和版本都更容易统一管理：

- 官方仓库：`https://github.com/PDB-REDO/dssp`

本项目默认建议将 DSSP 安装在项目目录中，便于环境自包含和团队复现。推荐目录结构如下：

- 源码目录：`./vendor/dssp`
- 构建目录：`./vendor/dssp/build`
- 安装目录：`./.local/dssp`


先在项目根目录下准备 DSSP 源码：

```bash
mkdir -p ./vendor
git clone https://github.com/PDB-REDO/dssp.git ./vendor/dssp
cd ./vendor/dssp
```

然后将 DSSP 编译并安装到项目目录下的 `./.local/dssp`：

```bash
cmake -S . -B build -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="$(cd ../.. && pwd)/.local/dssp"

cmake --build build -j"$(nproc)"
cmake --install build
```

安装完成后，项目内常见的可执行文件路径是：

- `./.local/dssp/bin/mkdssp`
- `./.local/dssp/bin/dssp`

### 5. 验证外部工具

先把项目内 DSSP 的 `bin` 目录加入当前 shell 的 `PATH`：

```bash
export PATH="$(pwd)/.local/dssp/bin:$PATH"
```

然后确认当前环境内可以直接找到这些命令：

```bash
which mkdssp || which dssp
blastp -version
mmseqs version
clustalo --version
```

如果这些命令都能正常输出版本信息，说明本地增强依赖已经基本就绪。


## 快速开始

### 单个结构文件

```bash
python -m openfoldpanel \
  --input ./model.pdb \
  --outdir ./out
```

### 压缩包输入，根目录多个结构文件视为一个 job

```bash
python -m openfoldpanel \
  --input ./models.tar.gz \
  --outdir ./out
```

### 压缩包输入，根目录多个一级子目录视为多个 job

```bash
python -m openfoldpanel \
  --input ./batch_jobs.zip \
  --outdir ./out
```

### 带 MSA 的增强示例

```bash
python -m openfoldpanel \
  --input ./models.zip \
  --outdir ./out \
  --msa-db ./blastdb/swissprot_fasta/uniprot_sprot.fasta \
  --max-homologs-displayed 5 \
  --evalue 1e-6
```

## 输入与输出

### 推荐输入前提

`OpenFoldPanel` 很适合整理和比较由 AlphaFold 2 或 AlphaFold 3 生成的结构结果：

- AlphaFold 2 常见输入是 `PDB`，AlphaFold 3 常见输入是 `CIF / mmCIF`，这几类格式当前都支持。
- 如果你要在同一个 job 里比较多份模型，最好让这些文件包含相同数量的蛋白链，并保持对应链的序列一致。当前实现会以首个成功解析的模型作为参考轴；当残基编号对不上时，会尝试做保守对齐，但输入越规整，结果越稳定。
- 单体、同源多聚体和杂多聚体都可以处理；结构里可以带修饰残基、核酸、配体或离子。
- 从报告可读性和人工复核体验出发，建议把单个比较批次控制在 `25` 个模型文件以内。当前 CLI 的输入入口是单个结构文件或压缩包；如果你要比较多份结果，建议把它们整理到同一个 job 目录或同一个压缩包里再运行。

### 支持的输入

- 结构文件：`.pdb`、`.cif`、`.mmcif`
- 压缩格式：`.zip`、`.tar`、`.tar.gz`、`.tgz`、`.tar.bz2`、`.tbz2`、`.tar.xz`、`.txz`

### job 识别规则

- 单个结构文件本身就是一个 job。
- 如果压缩包根目录下是多个一级子目录，每个一级子目录会被视为一个 job。
- 如果压缩包根目录下直接是多个结构文件，整个压缩包会被视为一个 job。
- 同一个 job 里的结构文件会按自然排序处理，非结构文件会被忽略并写入日志。

### 每个 job 的典型输出

- `report.html`
- `reference-chain-<CHAIN>.pdf`
- `tracks.json`
- `summary.txt`
- `logs.txt`

如果 PDF 导出依赖不可用，程序会跳过 PDF，但仍然生成 HTML、JSON 和文本摘要，并把该 job 记为 `partial_success`。

## 关键参数速览

| 参数 | 默认值 | 说明 |
| --- | --- | --- |
| `--input PATH` | 无 | 输入结构文件或压缩包。 |
| `--outdir OUTDIR` | 无 | 输出目录。每个 job 会生成独立结果子目录。 |
| `--chain ALL\|CHAIN_ID` | `ALL` | 为所有蛋白链出图，或只分析指定参考链。 |
| `--columns INT` | `80` | 每个渲染 block 显示的残基列数。 |
| `--font-size INT` | `12` | 报告和图板的基础字号。 |
| `--hyd-window INT` | `3` | 疏水性轨道的平滑窗口大小。 |
| `--msa-db PATH` | 无 | 本地 BLAST / MMseqs 数据库前缀，或蛋白 FASTA 文件。 |
| `--max-homologs-displayed INT` | `5` | 同源序列检索和展示上限，允许范围 `0-25`。 |
| `--evalue VALUE` | `1e-6` | BLAST / MMseqs 显著性阈值，只接受固定枚举值。 |
| `--disable-msa` | 关闭 | 显式关闭同源搜索、比对和 conservation 计算。 |
| `--keep-temp` | 关闭 | 保留解压、检索和比对过程中生成的临时目录。 |
| `--contact-cutoff FLOAT` | `3.7` | 弱接触判定阈值，单位为 Angstrom。 |
| `--strong-contact-cutoff FLOAT` | `3.2` | 强接触判定阈值，单位为 Angstrom。 |
| `--verbose` | 关闭 | 输出详细日志。 |

几个最常用的参数值得单独留意：

- `--chain`
  默认值是 `ALL`。程序会收集首个成功解析结构中的蛋白链，并为每条链分别出图；如果你只关心某一条链，可以显式指定，例如 `--chain B`。
- `--msa-db`
  可以传 BLAST / MMseqs 数据库前缀，也可以直接传蛋白 FASTA 文件。数据库下载、建库和传参示例见 [blastdb/README.md](/workspaces/re-foldscript/blastdb/README.md)。
- `--max-homologs-displayed`
  控制最多检索多少条、最多显示多少条 homolog。设为 `0` 时会跳过 homolog 行，只保留 query 行。
- `--evalue`
  当前允许值为 `1e-4`、`1e-5`、`1e-6`、`1e-7`、`1e-8`、`1e-9`、`1e-10`、`1e-11`、`1e-12`。它控制命中是否通过显著性筛选，不直接控制结果条数。
- `--disable-msa`
  如果你只想快速生成结构图板，或当前环境没有 BLAST / MMseqs / Clustal Omega，可以直接关闭 MSA 流程。

## 致谢与参考文献

`OpenFoldPanel` 在产品思路上受 FoldScript 很大启发。感谢 FoldScript 团队把“比较多组 AI 蛋白模型并快速读出差异”这件事做得足够清楚，也给了这个项目很明确的参考方向。

- Robert, X., Guillon, C., Gouet, P. (2025). *FoldScript: a web server for the efficient analysis of AI-generated 3D protein models*. *Nucleic Acids Research*, 53(W1), W277-W282. DOI: [10.1093/nar/gkaf326](https://doi.org/10.1093/nar/gkaf326)
