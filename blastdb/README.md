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
| Raw PDB cluster files | `pdb_cluster_src/` | `pdb_seqres.txt` and `clusters-by-entity-{50,70,90,95}.txt` | Inputs for building chain-level `PDBAA50/70/90/95` FASTA and BLAST databases |
| Download script | `Download_db.sh` | Downloads raw files and prebuilt BLAST databases | Fetch resources first |
| Generic build script | `Build_blastdb.sh` | Builds a BLAST database from FASTA or sequence files | FASTA -> BLAST prefix |
| PDBAA cluster build script | `Build_pdbaa_clusters.py` | Builds `PDBAA50/70/90/95` chain-level FASTA and optional BLAST prefixes from raw PDB cluster data | Use after `pdb-clusters` download |

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

Three common paths are supported:

1. Use the prebuilt `pdbaa` BLAST database directly
2. Download a raw FASTA file and build your own database with `Build_blastdb.sh`
3. Download `pdb-clusters` and build chain-level `PDBAA50/70/90/95` with `Build_pdbaa_clusters.py`

If you only want the fastest way to get the pipeline working, using `pdbaa` directly is the simplest option.
If you need representative-only PDB sequence databases that align with the RCSB entity clusters at `50 / 70 / 90 / 95%` identity, use the `PDBAA50/70/90/95` workflow below.

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

## Build PDBAA50 / 70 / 90 / 95

`Build_pdbaa_clusters.py` converts RCSB cluster representatives into chain-level FASTA headers like `pdb|5B8C|A`, then optionally builds BLAST database prefixes from those FASTA files.

The workflow is:

1. Treat `clusters-by-entity-{50,70,90,95}.txt` as the authoritative representative definition
2. Keep only the first token from each cluster line
3. Keep only PDB representative tokens such as `5B8C_2`; skip `AF_...` and `MA_...`
4. Use `polymer_entity` metadata to resolve `entry + entity -> canonical auth_asym_id`
5. Pull the sequence itself from local `pdb_seqres.txt` protein-chain records only
6. Write FASTA headers as `pdb|ENTRY|CHAIN`

Required input directory contents:

- `pdb_seqres.txt`
- `clusters-by-entity-50.txt`
- `clusters-by-entity-70.txt`
- `clusters-by-entity-90.txt`
- `clusters-by-entity-95.txt`

Build a single identity FASTA:

```bash
python ./blastdb/Build_pdbaa_clusters.py \
  --identity 95 \
  --pdb-cluster-src ./blastdb/pdb_cluster_src
```

Build the FASTA and BLAST database prefix together:

```bash
python ./blastdb/Build_pdbaa_clusters.py \
  --identity 95 \
  --pdb-cluster-src ./blastdb/pdb_cluster_src \
  --build-blastdb
```

Build every supported identity in one pass:

```bash
python ./blastdb/Build_pdbaa_clusters.py \
  --all \
  --pdb-cluster-src ./blastdb/pdb_cluster_src \
  --build-blastdb
```

Outputs are written under `./blastdb/build/`:

- `./blastdb/build/pdbaa50/pdbaa50.fasta`
- `./blastdb/build/pdbaa50/pdbaa50`
- `./blastdb/build/pdbaa70/pdbaa70.fasta`
- `./blastdb/build/pdbaa70/pdbaa70`
- `./blastdb/build/pdbaa90/pdbaa90.fasta`
- `./blastdb/build/pdbaa90/pdbaa90`
- `./blastdb/build/pdbaa95/pdbaa95.fasta`
- `./blastdb/build/pdbaa95/pdbaa95`

The build cache is shared across identities:

- `./blastdb/build/pdbaa_cache/entity_chain_map.json`

Recommended acceptance checks:

- FASTA headers match `^pdb\|[0-9A-Z]{4}\|[^|[:space:]]+$`
- FASTA does not contain `AF_`, `MA_`, or `mol:na`
- Multi-chain entries such as `5B8C` resolve to chain-level headers rather than entity tokens
- `blastdbcmd` can retrieve a record directly by chain-level identifier

Example:

```bash
blastdbcmd \
  -db ./blastdb/build/pdbaa95/pdbaa95 \
  -entry 'pdb|5B8C|A'
```

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
