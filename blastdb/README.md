## blastdb 目录说明

这个目录用来存放 `OpenFoldPanel` 使用的序列数据库、下载脚本、建库脚本，以及 `PDB cluster` 的原始数据文件。

当前实现里，`--msa-db` 既支持：

- 现成的 BLAST 数据库前缀，例如 `./blastdb/pdbaa/pdbaa`
- 原始蛋白 FASTA，例如 `./blastdb/swissprot_fasta/uniprot_sprot.fasta`

不过这里推荐的长期使用方式是：先下载原始序列文件，再建成本地 BLAST 数据库，最后把数据库前缀传给 `--msa-db`。

## 目录里的资源

| 资源 | 路径 | 用途 | 推荐用法 |
| --- | --- | --- | --- |
| UniProt Swiss-Prot FASTA | `swissprot_fasta/` | 通用蛋白原始序列输入 | 建库后传生成的前缀 |
| NCBI swissprot BLAST 库 | `swissprot/` | 已建好的 Swiss-Prot BLAST 数据库 | 直接传 `./blastdb/swissprot/swissprot` |
| NCBI pdbaa BLAST 库 | `pdbaa/` | 已建好的 PDBAA BLAST 数据库 | 直接传 `./blastdb/pdbaa/pdbaa` |
| PDB cluster 原始文件 | `pdb_cluster_src/` | `pdb_seqres.txt` 和 `clusters-by-entity-{50,70,90,95}.txt` | 当前仅保留原始数据，后续建库脚本 |
| 下载脚本 | `Download_db.sh` | 下载原始文件和现成 BLAST 库 | 先下载资源 |
| 通用建库脚本 | `Build_blastdb.sh` | 把 FASTA 或序列文件建成 BLAST 数据库 | FASTA -> BLAST 前缀 |

## 依赖要求

运行本目录下的脚本前，建议先确认这些命令可用：

- `update_blastdb.pl`
  来自 NCBI BLAST+，用于下载 `pdbaa` 等现成 BLAST 数据库
- `curl` 或 `wget`
  用于下载 FASTA 和 PDB cluster 原始文件
- `gunzip`
  用于解压 `.gz` 文件
- `makeblastdb`
  用于把本地序列文件建成可直接传给 `--msa-db` 的 BLAST 数据库前缀
如果你使用 `conda` / `mamba` 环境，通常安装 `blast` 或 `ncbi-blast+` 后就会带上 `update_blastdb.pl` 和 `makeblastdb`。

## 推荐工作流

这里有两条常见路径：

1. 直接使用现成的 `pdbaa` BLAST 数据库
2. 下载原始 FASTA 后，用 `Build_blastdb.sh` 建库

如果你只想尽快跑通流程，最简单的是直接使用 `pdbaa`。
`PDB cluster` 相关原始数据会继续保留在 `pdb_cluster_src/`，具体处理方法将在后续版本更新。

## 下载脚本用法

默认会下载推荐组合：`swissprot-fasta` 和 `pdbaa`。

```bash
bash ./blastdb/Download_db.sh
```

如果你只想下载特定资源，可以显式指定目标：

```bash
bash ./blastdb/Download_db.sh swissprot-fasta
bash ./blastdb/Download_db.sh pdbaa
bash ./blastdb/Download_db.sh pdb-clusters
```

如果你想把脚本支持的资源一次性都拉下来：

```bash
bash ./blastdb/Download_db.sh all
```

如果某些普通文件已经存在，但你想强制重新下载：

```bash
bash ./blastdb/Download_db.sh --force pdb-clusters
```

脚本支持的目标如下：

- `recommended`
  下载推荐组合：`swissprot-fasta` 和 `pdbaa`
- `swissprot-fasta`
  下载 `uniprot_sprot.fasta`
- `swissprot-blast`
  下载 NCBI 提供的 `swissprot` BLAST 数据库
- `pdbaa`
  下载 NCBI 提供的 `pdbaa` BLAST 数据库
- `pdb-clusters`
  下载 `pdb_seqres.txt` 和 `clusters-by-entity-{50,70,90,95}.txt`
- `all`
  下载上面所有资源

## 通用建库脚本用法

如果你已经有本地序列文件，想直接建一个 BLAST 数据库前缀，使用：

```bash
bash ./blastdb/Build_blastdb.sh --input <source_file> --out-prefix <prefix>
```

例如，把 Swiss-Prot FASTA 建成 BLAST 数据库：

```bash
bash ./blastdb/Build_blastdb.sh \
  --input ./blastdb/swissprot_fasta/uniprot_sprot.fasta \
  --out-prefix ./blastdb/swissprot/swissprot \
  --title SWISSPROT
```

建库完成后，可直接这样传给 `OpenFoldPanel`：

```bash
python -m openfoldpanel \
  --input ./model.pdb \
  --outdir ./out \
  --msa-db ./blastdb/build/swissprot/swissprot \
  --max-homologs-displayed 5
```

## PDB cluster 数据说明

`pdb_cluster_src/` 中保留了 `pdb_seqres.txt` 和 `clusters-by-entity-{50,70,90,95}.txt` 等原始文件，方便后续继续整理 `PDB cluster` 相关数据库流程。

当前版本暂不提供 `PDBAA50/70/90/95` 的正式构建方法，相关处理方式将在后续更新中补充。

## 直接使用现成 BLAST 数据库

如果你不需要自己建库，也可以直接使用已经下载好的 BLAST 数据库前缀。

使用 `swissprot`：

```bash
python -m openfoldpanel \
  --input ./model.pdb \
  --outdir ./out \
  --msa-db ./blastdb/swissprot/swissprot \
  --max-homologs-displayed 5
```

使用 `pdbaa`：

```bash
python -m openfoldpanel \
  --input ./model.pdb \
  --outdir ./out \
  --msa-db ./blastdb/pdbaa/pdbaa \
  --max-homologs-displayed 5
```

## 补充说明

- 传 BLAST 数据库时，`--msa-db` 应该指向数据库前缀，而不是某个具体后缀文件；例如传 `./blastdb/pdbaa/pdbaa`，不是 `./blastdb/pdbaa/pdbaa.psq`
- 代码层面仍然兼容直接传 FASTA；如果本地有 `blastp` 和 `makeblastdb`，程序会临时建库；如果没有 `makeblastdb` 但有 `mmseqs`，程序会直接把 FASTA 交给 `mmseqs easy-search`
- 作为长期可复现的使用方式，更推荐把原始序列先建成你自己的本地 BLAST 数据库，再传数据库前缀
