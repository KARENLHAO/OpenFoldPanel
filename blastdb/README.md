## blastdb 目录说明

这个目录用来存放 `OpenFoldPanel` 使用的序列数据库、下载脚本、建库脚本和配套说明。

当前实现里，`--msa-db` 既支持：

- 现成的 BLAST 数据库前缀，例如 `./blastdb/pdbaa/pdbaa`
- 原始蛋白 FASTA，例如 `./blastdb/swissprot_fasta/uniprot_sprot.fasta`

不过这里推荐的流程是：先下载原始序列文件，再用 `makeblastdb` 建成本地 BLAST 数据库，最后把数据库前缀传给 `--msa-db`。

## 目录里的资源

| 资源 | 下载目标 | 用途 | 推荐用法 |
| --- | --- | --- | --- |
| UniProt Swiss-Prot FASTA | `swissprot_fasta/` | 推荐的通用原始序列输入 | 建库后传生成的前缀 |
| NCBI swissprot BLAST 库 | `swissprot/` | 已建好的 Swiss-Prot BLAST 数据库 | 直接传 `./blastdb/swissprot/swissprot` |
| NCBI pdbaa BLAST 库 | `pdbaa/` | PDB 派生蛋白序列数据库 | 直接传 `./blastdb/pdbaa/pdbaa` |
| PDB cluster 原始文件 | `pdb_cluster_src/` | 可选的聚类源数据和 `pdb_seqres` 原始序列文件 | `pdb_seqres.txt` 建库后再传前缀 |

## 依赖要求

运行 [Download_db.sh](/workspaces/re-foldscript/blastdb/Download_db.sh) 和 [Build_blastdb.sh](/workspaces/re-foldscript/blastdb/Build_blastdb.sh) 前，建议先确认这些命令可用：

- `update_blastdb.pl`
  来自 NCBI BLAST+，用于下载 `swissprot` 和 `pdbaa` 这类 BLAST 数据库
- `curl` 或 `wget`
  用于下载 FASTA 和 PDB cluster 原始文件
- `gunzip`
  用于解压 `.gz` 文件
- `makeblastdb`
  用于把本地序列文件建成可直接传给 `--msa-db` 的 BLAST 数据库前缀

如果你使用 `conda` / `mamba` 环境，通常安装 `blast` 或 `ncbi-blast+` 后就会带上 `update_blastdb.pl`。

## 推荐工作流

推荐按下面三步走：

1. 用 [Download_db.sh](/workspaces/re-foldscript/blastdb/Download_db.sh) 下载原始序列文件或现成数据库
2. 如果拿到的是原始序列文件，用 [Build_blastdb.sh](/workspaces/re-foldscript/blastdb/Build_blastdb.sh) 建成 BLAST 数据库
3. 把数据库前缀传给 `--msa-db`

## 下载脚本用法

默认会下载推荐组合：`swissprot-fasta` 和 `pdbaa`。

```bash
bash ./blastdb/Download_db.sh
```

如果你只想下载特定资源，可以显式指定目标：

```bash
bash ./blastdb/Download_db.sh swissprot-fasta
bash ./blastdb/Download_db.sh pdbaa
bash ./blastdb/Download_db.sh swissprot-blast
```

如果你想把脚本支持的资源一次性都拉下来：

```bash
bash ./blastdb/Download_db.sh all
```

如果某些普通文件已经存在，但你想强制重新下载：

```bash
bash ./blastdb/Download_db.sh --force swissprot-fasta
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
  下载 RCSB 提供的 `pdb_seqres` 和 cluster 源文件
- `all`
  下载上面所有资源

## 建库脚本用法

如果你决定统一先建 BLAST 数据库，再把数据库前缀传给 `--msa-db`，使用：

```bash
bash ./blastdb/Build_blastdb.sh --input <source_file> --out-prefix <prefix>
```

### 1. 由 Swiss-Prot FASTA 建库

```bash
bash ./blastdb/Build_blastdb.sh \
  --input ./blastdb/swissprot_fasta/uniprot_sprot.fasta \
  --out-prefix ./blastdb/build/swissprot/swissprot \
  --title SWISSPROT
```

建库完成后，传参方式是：

```bash
python -m openfoldpanel \
  --input ./model.pdb \
  --outdir ./out \
  --msa-db ./blastdb/build/swissprot/swissprot \
  --max-homologs-displayed 5
```

### 2. 由 pdb_seqres.txt 建库

```bash
bash ./blastdb/Build_blastdb.sh \
  --input ./blastdb/pdb_cluster_src/pdb_seqres.txt \
  --out-prefix ./blastdb/build/pdb_seqres/pdb_seqres \
  --title PDB_SEQRES
```

建库完成后，传参方式是：

```bash
python -m openfoldpanel \
  --input ./model.pdb \
  --outdir ./out \
  --msa-db ./blastdb/build/pdb_seqres/pdb_seqres \
  --max-homologs-displayed 5
```

### 3. 使用现成 BLAST 数据库时不需要建库

如果你下载的是现成 BLAST 数据库，例如 `swissprot` 或 `pdbaa`，就不需要再运行 `Build_blastdb.sh`，可以直接传前缀。

## 给 OpenFoldPanel 的传参方式

### 1. 使用已建好的 Swiss-Prot BLAST 数据库

```bash
python -m openfoldpanel \
  --input ./model.pdb \
  --outdir ./out \
  --msa-db ./blastdb/swissprot/swissprot \
  --max-homologs-displayed 5
```

### 2. 使用已建好的 PDBAA BLAST 数据库

```bash
python -m openfoldpanel \
  --input ./model.pdb \
  --outdir ./out \
  --msa-db ./blastdb/pdbaa/pdbaa \
  --max-homologs-displayed 5
```

## 补充说明

- 传 BLAST 数据库时，`--msa-db` 应该指向数据库前缀，而不是某个具体后缀文件；例如传 `./blastdb/pdbaa/pdbaa`，不是 `./blastdb/pdbaa/pdbaa.psq`
- `pdb_cluster_src/pdb_seqres.txt` 不能直接传给 `--msa-db`。它是原始序列文本，推荐先用 `Build_blastdb.sh` 建库，再传生成的前缀
- 代码层面仍然兼容直接传 FASTA；如果本地有 `blastp` 和 `makeblastdb`，程序会临时建库；如果没有 `makeblastdb` 但有 `mmseqs`，程序会直接把 FASTA 交给 `mmseqs easy-search`
- 作为长期可复现的使用方式，更推荐把原始序列先建成你自己的本地 BLAST 数据库，再传数据库前缀
