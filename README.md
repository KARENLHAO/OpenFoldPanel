# OpenFoldPanel

`OpenFoldPanel` 是一个本地可运行的开源命令行工具，用于从单个结构文件或结构压缩包生成接近 FoldScript 左侧 flat representation panel 的结果图。

它不调用 FoldScript 在线服务，不依赖任何私有组件。项目重点是：

- 正确处理单文件和批量压缩包输入
- 输出统一 `tracks.json`
- 原生 SVG 排版与 block wrapping
- 导出 SVG / PDF / PNG
- DSSP / MSA 缺失时优雅降级

## 原 FoldScript 组件到当前实现的映射

| FoldScript 侧能力 | 当前实现 |
| --- | --- |
| 结构解析 / 链与残基处理 | Gemmi 为主，纯 Python 解析回退 |
| PDB / mmCIF | Gemmi 为主，纯 Python PDB/mmCIF `_atom_site` 回退 |
| DSSP 特征 | 外部 `mkdssp` / `dssp` 调用；缺失时用几何近似回退 |
| 接触分析 | 本地重原子距离搜索，兼容 Gemmi/Biopython 思路 |
| 同源搜索 | `blastp` 或 `mmseqs easy-search` |
| 多序列比对 | `clustalo` |
| flat panel 绘制 | 本地原生 SVG 渲染 |
| PDF / PNG 导出 | CairoSVG |

## 已实现能力

- 单个 `.pdb` / `.cif` / `.mmcif` 输入
- `.zip` / `.tar` / `.tar.gz` / `.tgz` / `.tar.bz2` / `.tbz2` / `.tar.xz` / `.txz` 输入
- 安全解压，防 Zip Slip / Tar 路径穿越
- 批量 job 发现和自然排序
- 参考链自动选择
- 多模型统一到同一序列轴
- `tracks.json` 输出
- 二级结构轨道
- query + homolog MSA 轨道
- accessibility 轨道
- hydropathy 轨道
- contacts 轨道
- SVG block wrapping
- PDF / PNG 导出接口

## 安装

推荐使用单独的 `conda/mamba` 环境统一安装 OpenFoldPanel、BLAST+、MMseqs2、Clustal Omega，以及编译 DSSP 所需工具链。这样可以避免 `.venv`、系统 PATH 和外部命令混用。

### 1. 创建统一 conda 环境

建议 Python 3.10 及以上。

```bash
mamba create -n openfoldpanel -c conda-forge -c bioconda \
  python=3.10 pip \
  blast mmseqs2 clustalo \
  cmake ninja cxx-compiler pkg-config git

conda activate openfoldpanel
```

上面这一步会安装：

- OpenFoldPanel 运行所需 Python 版本
- BLAST+
  - `blastp`
- MMseqs2
  - `mmseqs`
- Clustal Omega
  - `clustalo`
- DSSP 源码编译所需工具链
  - `cmake`
  - `ninja`
  - C/C++ 编译器
  - `pkg-config`
  - `git`

### 2. 安装 OpenFoldPanel 本身

在同一个 conda 环境中安装项目：

```bash
python -m pip install -r requirements.txt
python -m pip install -e .
```

### 3. 安装 DSSP

本项目推荐 DSSP 采用官方仓库源码编译路线：

- 官方仓库：`https://github.com/PDB-REDO/dssp`

本项目默认建议将 DSSP 安装在项目目录中，便于环境自包含和团队复现。推荐目录结构如下：

- 源码目录：`./vendor/dssp`
- 构建目录：`./vendor/dssp/build`
- 安装目录：`./.local/dssp`

建议将 `vendor/dssp` 和 `.local/dssp` 加入 `.gitignore`，避免误提交第三方源码和构建产物。

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

安装完成后，项目内通常会生成：

- `./.local/dssp/bin/mkdssp`
或
- `./.local/dssp/bin/dssp`

### 4. 验证外部工具

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

如果上面都能正常输出版本信息，OpenFoldPanel 的完整增强路径就可用了。

### 5. 运行 OpenFoldPanel

后续请始终先激活同一个环境，并确保项目内 DSSP 在 `PATH` 中：

```bash
conda activate openfoldpanel
export PATH="$(pwd)/.local/dssp/bin:$PATH"
python -m openfoldpanel --input ./model.pdb --outdir ./out
```

### 6. DSSP 兜底方案

如果 DSSP 官方源码编译失败，可以临时采用社区 conda 包作为兜底：

```bash
mamba install -n openfoldpanel -c sbl dssp
```

说明：

- 这是兜底方案，不是默认方案
- 默认仍推荐使用 `PDB-REDO/dssp` 官方仓库源码编译版

### 7. 可选的 Ubuntu 系统包方案

如果你更偏向系统包管理，也可以用 Ubuntu 自带包安装部分工具：

```bash
sudo apt-get update
sudo apt-get install -y dssp ncbi-blast+ mmseqs2 clustalo
```

但需要注意：

- Ubuntu 20.04 自带的 `dssp` 版本通常较旧
- 因此系统包方案不作为本项目默认推荐路线
- 本 README 的主流程仍是统一 `conda/mamba` 环境 + DSSP 官方源码编译

## 运行方式

### 单文件

```bash
python -m openfoldpanel --input ./model.pdb --outdir ./out
```

### 压缩包，根目录多个结构文件，视为单 job

```bash
python -m openfoldpanel --input ./models.tar.gz --outdir ./out
```

### 压缩包，根目录多个一级子目录，视为多个 job

```bash
python -m openfoldpanel --input ./batch_jobs.zip --outdir ./out
```

### 带更多参数

```bash
python -m openfoldpanel \
  --input ./models.zip \
  --outdir ./out \
  --chain B \
  --columns 80 \
  --font-size 12 \
  --hyd-window 3 \
  --msa-db /path/to/local/db \
  --max-hits 10 \
  --contact-cutoff 3.7 \
  --strong-contact-cutoff 3.2 \
  --verbose
```

## CLI 参数

- `--input PATH`
- `--outdir OUTDIR`
- `--chain ALL|A`
- `--columns 80`
- `--font-size 12`
- `--hyd-window 3`
- `--msa-db /path/to/db`
- `--max-hits 10`
- `--disable-msa`
- `--keep-temp`
- `--contact-cutoff 3.7`
- `--strong-contact-cutoff 3.2`
- `--verbose`

## 输入规则

### 支持的结构输入

- `.pdb`
- `.cif`
- `.mmcif`

### 支持的压缩格式

- `.zip`
- `.tar`
- `.tar.gz`
- `.tgz`
- `.tar.bz2`
- `.tbz2`
- `.tar.xz`
- `.txz`

### 批量 job 规则

1. 单个结构文件本身就是一个 job。
2. 如果输入是压缩包：
   - 根目录存在多个一级子目录时，每个一级子目录视为一个 job。
   - 根目录直接是多个结构文件时，整个压缩包视为一个 job。
   - 某个一级子目录下有多个结构文件时，这些文件属于该 job 的多个模型。
3. 每个 job 中的结构文件按自然排序处理，例如 `ranked_1` 在 `ranked_10` 前面。
4. 非结构文件会被忽略，并记录到日志。

## 参考链规则

默认 `--chain ALL`。

`ALL` 规则：

1. 收集首个解析成功结构中的所有蛋白链
2. 每条蛋白链都会生成独立 PDF
3. HTML 报告中可通过按钮切换不同参考链

若显式指定 `--chain A`，则要求该链存在且是蛋白链。

`AUTO` 已废弃，不再接受作为用户输入参数。HTML 默认激活链仍会在内部按旧规则自动选择：

1. 优先选择蛋白链 `A`
2. 否则选择第一个蛋白链

## 轨道说明

### 1. Secondary structure

- 优先使用 DSSP
- 若 DSSP 不可用，使用几何近似回退
- 语义：
  - `H/G/I` -> helix
  - `E/B` -> strand
  - `T` -> turn
  - 其他 -> coil

### 2. Accessibility

- 优先使用 DSSP ASA 换算相对可及性
- 若 DSSP 不可用，使用局部原子拥挤度近似
- 分类：
  - `buried < 0.1`
  - `intermediate 0.1 - 0.4`
  - `accessible 0.4 - 1.0`
  - `highly_exposed > 1.0`

### 3. Hydropathy

- Kyte-Doolittle
- 默认窗口 `3`
- 分类：
  - `hydrophilic < -1.5`
  - `intermediate -1.5 ~ 1.5`
  - `hydrophobic > 1.5`

### 4. Contacts

- 默认仅重原子
- 计算参考链每个残基到其他链 / 核酸 / 配体 / 离子的最短原子距离
- 阈值：
  - 强接触：`< 3.2 A`
  - 弱接触：`3.2 - 3.7 A`
- 图上主符号按最短距离选择
- JSON 保存所有命中

符号规则：

- 蛋白-蛋白：对方链 ID
- 跨链同名同号残基：`#`
- 核酸：`*`
- 卟啉类：`:`
- 糖：`"`
- 离子：`+`
- 其他小分子：`^`

### 5. MSA

完整模式需要：

- 未指定 `--disable-msa`
- 提供 `--msa-db`
- 本地可用 `blastp` 或 `mmseqs`
- 本地可用 `clustalo`

否则自动跳过 MSA，但 query 行仍会显示。

保守性着色规则：

- `100% identity` -> 白字红底
- `70%-99%` 同类保守 -> 红字黄底
- 其他 -> 黑字透明底

相似性分组：

- `HKR`
- `DE`
- `STNQ`
- `AVLIM`
- `FYW`

## 输出内容

每个 job 生成独立目录，至少包含：

- `report.html`
- `reference-chain-A.pdf`
- `reference-chain-B.pdf` 等
- `tracks.json`
- `summary.txt`
- `logs.txt`

说明：

- 不再默认落盘 `panel.svg` 或 `panel.png`
- SVG 仅作为内存中的中间渲染格式，用于 PDF 导出和 HTML 内联
- 若 CairoSVG 未安装，则 PDF 会跳过，但 `report.html`、`tracks.json`、`summary.txt`、`logs.txt` 仍会生成，并将 job 标记为 `partial_success`

## tracks.json 结构

顶层字段包括：

- `job_name`
- `default_reference_chain`
- `chain_panels`
- `warnings`
- `status`

`chain_panels[]` 内部保留单参考链视图，字段至少包括：

- `reference_chain`
- `sequence_axis`
- `models`
- `msa`
- `hydropathy`
- `render_config`
- `warnings`
- `status`

模型字段至少包括：

- `name`
- `display_name`
- `chain`
- `secondary_structure`
- `plddt`
- `accessibility`
- `contacts`

## 渲染布局

- 单个 job 输出一个静态 HTML 报告，顶部提供参考链切换按钮
- 每条参考链输出一个科学报告风 PDF
- 核心 panel 仍为固定左标签列 + 右侧残基网格
- 标签列宽度会根据可见标签动态扩展
- 支持 `--columns` 控制每个 block 列数
- 超出宽度自动 block wrapping，不做横向压缩
- 默认行顺序：
  - residue positions
  - 每个模型的 secondary structure
  - query sequence
  - homolog alignment rows
  - residue accessibility
  - hydropathy profile
  - 每个模型的 inter-chain contacts

## 测试

```bash
conda activate openfoldpanel
pytest
```

当前测试覆盖：

- 安全解压
- job 发现与自然排序
- hydropathy
- contact 逻辑
- layout wrapping
- CLI smoke

## 已知限制

- 纯 Python mmCIF 回退解析只覆盖常见 `_atom_site` 场景，不是完整 mmCIF 解释器。
- DSSP 缺失时 secondary structure 和 accessibility 使用近似算法，适合离线预览，不应替代高质量结构注释。
- PDF 导出依赖 CairoSVG；若环境缺失，仅保证 HTML、JSON 与文本摘要。
- MSA 搜索依赖本地数据库与外部程序，不会自动下载数据库。
- 当前接触分类表覆盖了常见核酸、糖、金属离子和卟啉类，复杂配体家族仍可能落入 `other_ligand`。
