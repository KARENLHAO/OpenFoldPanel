## Swiss-Prot / PDBA A 数据库说明

本目录里现在同时可能存在两类 MSA 输入：

- 传统 BLAST 数据库前缀，例如 `./swissprot/swissprot`
- 原始蛋白 FASTA，例如 `./swissprot_fasta/uniprot_sprot.fasta`

对 OpenFoldPanel 来说，`--msa-db` 两种都可以传。

### Swiss-Prot

如果你希望直接以 UniProt/Swiss-Prot 的 FASTA 作为同源搜索输入，可以这样传：

```bash
python -m openfoldpanel \
  --input ./model.pdb \
  --outdir ./out \
  --msa-db ./blastdb/swissprot_fasta/uniprot_sprot.fasta \
  --max-homologs-displayed 5
```

行为说明：

- 若本地有 `blastp` 和 `makeblastdb`，程序会先把 `uniprot_sprot.fasta` 临时建成 BLAST 蛋白库
- 若没有 `makeblastdb` 但有 `mmseqs`，则会直接用 FASTA 跑 `mmseqs easy-search`
- homolog 标签会优先从 FASTA 标识头恢复，例如 `sp|P01674|KV3AM_MOUSE`

如果你已经手工准备好了 BLAST 库，也可以继续传数据库前缀：

```bash
--msa-db ./blastdb/swissprot/swissprot
```
