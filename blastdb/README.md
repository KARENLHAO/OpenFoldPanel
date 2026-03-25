## `blastdb` Directory Guide

This directory stores the sequence databases, download scripts, build scripts, and raw `PDB cluster` files used by `OpenFoldPanel`.

In the current implementation, `--msa-db` supports both:

- A ready-to-use BLAST database prefix, such as `./blastdb/pdbaa/pdbaa`
- A raw protein FASTA file, such as `./blastdb/swissprot_fasta/uniprot_sprot.fasta`

For long-term reproducible use, the recommended workflow is to download the raw sequence files first, build a local BLAST database from them, and then pass the resulting database prefix to `--msa-db`.

## Resources in This Directory

| Resource | Path | Purpose | Recommended Use |
| --- | --- | --- | --- |
| UniProt Swiss-Prot FASTA | `swissprot_fasta/` | General protein sequence source input | Build a database, then pass the generated prefix |
| NCBI `swissprot` BLAST DB | `swissprot/` | Prebuilt Swiss-Prot BLAST database | Pass `./blastdb/swissprot/swissprot` directly |
| NCBI `pdbaa` BLAST DB | `pdbaa/` | Prebuilt PDBAA BLAST database | Pass `./blastdb/pdbaa/pdbaa` directly |
| Raw PDB cluster files | `pdb_cluster_src/` | `pdb_seqres.txt` and `clusters-by-entity-{50,70,90,95}.txt` | Raw inputs retained for future database tooling |
| Download script | `Download_db.sh` | Downloads raw files and prebuilt BLAST databases | Fetch resources first |
| Generic build script | `Build_blastdb.sh` | Builds a BLAST database from FASTA or sequence files | FASTA -> BLAST prefix |

## Dependency Requirements

Before running scripts in this directory, verify that the following commands are available:

- `update_blastdb.pl`
  Provided by NCBI BLAST+, used to download prebuilt BLAST databases such as `pdbaa`
- `curl` or `wget`
  Used to download FASTA files and raw PDB cluster files
- `gunzip`
  Used to extract `.gz` files
- `makeblastdb`
  Used to turn local sequence files into BLAST database prefixes that can be passed directly to `--msa-db`

If you use a `conda` / `mamba` environment, installing `blast` or `ncbi-blast+` usually provides both `update_blastdb.pl` and `makeblastdb`.

## Recommended Workflow

Two common paths are supported:

1. Use the prebuilt `pdbaa` BLAST database directly
2. Download a raw FASTA file and build your own database with `Build_blastdb.sh`

If you only want the fastest way to get the pipeline working, using `pdbaa` directly is the simplest option.
Raw `PDB cluster` files remain in `pdb_cluster_src/`, and their downstream handling will be expanded in a future version.

## Download Script Usage

By default, the download script fetches the recommended combination: `swissprot-fasta` and `pdbaa`.

```bash
bash ./blastdb/Download_db.sh
```

If you only want a specific resource, pass the target explicitly:

```bash
bash ./blastdb/Download_db.sh swissprot-fasta
bash ./blastdb/Download_db.sh pdbaa
bash ./blastdb/Download_db.sh pdb-clusters
```

If you want to fetch every supported resource in one pass:

```bash
bash ./blastdb/Download_db.sh all
```

If some ordinary files already exist but you want to redownload them:

```bash
bash ./blastdb/Download_db.sh --force pdb-clusters
```

Supported targets are:

- `recommended`
  Downloads the recommended combination: `swissprot-fasta` and `pdbaa`
- `swissprot-fasta`
  Downloads `uniprot_sprot.fasta`
- `swissprot-blast`
  Downloads the NCBI `swissprot` BLAST database
- `pdbaa`
  Downloads the NCBI `pdbaa` BLAST database
- `pdb-clusters`
  Downloads `pdb_seqres.txt` and `clusters-by-entity-{50,70,90,95}.txt`
- `all`
  Downloads everything above

## Generic Database Build Script Usage

If you already have a local sequence file and want to create a BLAST database prefix directly, use:

```bash
bash ./blastdb/Build_blastdb.sh --input <source_file> --out-prefix <prefix>
```

For example, to build a BLAST database from the Swiss-Prot FASTA:

```bash
bash ./blastdb/Build_blastdb.sh \
  --input ./blastdb/swissprot_fasta/uniprot_sprot.fasta \
  --out-prefix ./blastdb/swissprot/swissprot \
  --title SWISSPROT
```

After the database is built, pass it to `OpenFoldPanel` like this:

```bash
python -m openfoldpanel \
  --input ./model.pdb \
  --outdir ./out \
  --msa-db ./blastdb/build/swissprot/swissprot \
  --max-homologs-displayed 5
```

## PDB Cluster Notes

`pdb_cluster_src/` keeps raw files such as `pdb_seqres.txt` and `clusters-by-entity-{50,70,90,95}.txt` so that future `PDB cluster` database workflows can be built on top of them.

The current version does not yet provide an official build path for `PDBAA50/70/90/95`. That workflow will be documented in a later update.

## Use a Prebuilt BLAST Database Directly

If you do not need to build your own database, you can point `--msa-db` directly at an existing BLAST database prefix.

Using `swissprot`:

```bash
python -m openfoldpanel \
  --input ./model.pdb \
  --outdir ./out \
  --msa-db ./blastdb/swissprot/swissprot \
  --max-homologs-displayed 5
```

Using `pdbaa`:

```bash
python -m openfoldpanel \
  --input ./model.pdb \
  --outdir ./out \
  --msa-db ./blastdb/pdbaa/pdbaa \
  --max-homologs-displayed 5
```

## Additional Notes

- When passing a BLAST database, `--msa-db` should point to the database prefix, not to one of the suffix files. For example, use `./blastdb/pdbaa/pdbaa`, not `./blastdb/pdbaa/pdbaa.psq`.
- The code still supports passing FASTA directly. If `blastp` and `makeblastdb` are available locally, the program builds a temporary database. If `makeblastdb` is unavailable but `mmseqs` is present, the FASTA is passed to `mmseqs easy-search`.
- For long-term reproducibility, building your own local BLAST database from raw sequence data and then passing that prefix remains the recommended workflow.
