# OpenFoldPanel

`OpenFoldPanel` compares multiple AI-generated protein structure models and organizes the key information scattered across structure files into a unified 2D flat-panel report. When you have several candidate models for the same target, looking only at the top-ranked structure is often not enough, especially around multimer interfaces, locally low-confidence regions, or disagreements between prediction systems.

Centered on a reference chain, `OpenFoldPanel` lays out secondary structure, accessibility, hydropathy, molecular contacts, confidence when available, and optional MSA / conservation information side by side. The goal is to make it easier to spot stable regions, regions that need manual review, and differences across model sets. Outputs include a browser-friendly `report.html`, chain-level PDFs, `tracks.json`, `summary.txt`, and `logs.txt`.

The project is inspired by FoldScript's comparison and presentation ideas, but it is implemented as an open-source, local-first CLI that can plug into batch workflows rather than a web server.

## Core Capabilities

- Accepts a single structure file or batch inputs packaged as archives such as `.zip` and `.tar.gz`.
- Supports multi-model, multi-chain, and multi-job processing so results from different predictors can be compared in one view.
- Generates `report.html`, `reference-chain-<CHAIN>.pdf`, `tracks.json`, `summary.txt`, and `logs.txt`.
- Provides tracks for secondary structure, accessibility, hydropathy, contacts, optional MSA / conservation, and confidence when present.
- Summarizes sequence-level, secondary-structure-level, and tertiary/quaternary interaction clues in one report: reference-chain and MSA rows for sequence context, DSSP or geometry-based secondary-structure tracks, and interaction summaries from accessibility, hydropathy, disulfides, confidence, and cross-chain / ligand / nucleic-acid / ion contacts.
- Supports `PDB`, `CIF`, `mmCIF`, and common archive formats so it can be inserted into existing prediction and curation workflows.

## Installation and Dependencies

Using an isolated `conda` or `mamba` environment is recommended. It makes it easier to keep Python packages, DSSP, and MSA tools in one reproducible setup.

### 1. Create an environment

```bash
mamba create -n openfoldpanel -c conda-forge -c bioconda \
  python=3.10 pip \
  cmake ninja cxx-compiler pkg-config git

conda activate openfoldpanel
```

### 2. Install the project

```bash
python -m pip install -r requirements.txt
python -m pip install -e .
```

### 3. Optional enhancement dependencies

The core `OpenFoldPanel` reporting flow only needs Python packages, but the tools below significantly improve feature completeness:

- `mkdssp` / `dssp`
  Used for more reliable secondary-structure and accessibility calculation. If unavailable, the program falls back to rough geometry-based estimates so reports can still be generated, but DSSP is more accurate.
- `blastp` or `mmseqs`
  Used for homolog search.
- `clustalo`
  Used for multiple-sequence alignment and conservation scoring.

If you want these tools in the same environment, continue with:

```bash
mamba install -n openfoldpanel -c conda-forge -c bioconda \
  blast mmseqs2 clustalo
```

### 4. Install DSSP (recommended)

If you only want to get the pipeline running quickly, you can temporarily install DSSP from a community package:

```bash
mamba install -n openfoldpanel -c sbl dssp
```

If DSSP is installed into a project-local directory such as `./.local/dssp/bin`, remember to add it to `PATH`:

```bash
export PATH="$(pwd)/.local/dssp/bin:$PATH"
```

For a managed production setup, building from the official repository is recommended because paths and versions are easier to keep consistent:

- Official repository: `https://github.com/PDB-REDO/dssp`

This project recommends installing DSSP under the project directory so the environment remains self-contained and easy to reproduce across a team. A suggested layout is:

- Source directory: `./vendor/dssp`
- Build directory: `./vendor/dssp/build`
- Install directory: `./.local/dssp`

Prepare the DSSP source tree from the project root:

```bash
mkdir -p ./vendor
git clone https://github.com/PDB-REDO/dssp.git ./vendor/dssp
cd ./vendor/dssp
```

Then build and install DSSP into `./.local/dssp`:

```bash
cmake -S . -B build -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="$(cd ../.. && pwd)/.local/dssp"

cmake --build build -j"$(nproc)"
cmake --install build
```

After installation, the most common project-local executables are:

- `./.local/dssp/bin/mkdssp`
- `./.local/dssp/bin/dssp`

### 5. Verify external tools

First, add the project-local DSSP `bin` directory to the current shell `PATH`:

```bash
export PATH="$(pwd)/.local/dssp/bin:$PATH"
```

Then verify that the current environment can resolve these commands directly:

```bash
which mkdssp || which dssp
blastp -version
mmseqs version
clustalo --version
```

If they all print version information successfully, the local enhancement toolchain is mostly ready.

For database download, database building, and `--msa-db` configuration, see [blastdb/README.md](./blastdb/README.md).

## Quick Start

### Single structure file

```bash
python -m openfoldpanel \
  --input ./model.pdb \
  --outdir ./out
```

### Archive input: multiple structure files at the archive root are treated as one job

```bash
python -m openfoldpanel \
  --input ./models.tar.gz \
  --outdir ./out
```

### Archive input: multiple first-level subdirectories at the archive root are treated as multiple jobs

```bash
python -m openfoldpanel \
  --input ./batch_jobs.zip \
  --outdir ./out
```

### Enhanced example with MSA

```bash
python -m openfoldpanel \
  --input ./models.zip \
  --outdir ./out \
  --msa-db ./blastdb/swissprot_fasta/uniprot_sprot.fasta \
  --max-homologs-displayed 5 \
  --evalue 1e-6
```

## Inputs and Outputs

### Recommended input conditions

`OpenFoldPanel` works especially well for organizing and comparing structures produced by AlphaFold 2 or AlphaFold 3:

- AlphaFold 2 commonly outputs `PDB`, while AlphaFold 3 commonly outputs `CIF / mmCIF`; all of these formats are currently supported.
- If you want to compare multiple models inside one job, it is best if they contain the same number of protein chains and matching chain sequences. The current implementation uses the first successfully parsed model as the reference axis; it will attempt conservative alignment if residue numbering differs, but cleaner inputs give more stable outputs.
- Monomers, homo-oligomers, and hetero-oligomers are all supported. Structures may include modified residues, nucleic acids, ligands, or ions.
- For report readability and manual review, keeping a comparison batch to `25` model files or fewer is recommended. The CLI currently accepts a single structure file or archive as input; if you need to compare many results, place them inside one job directory or one archive before running the tool.

### Supported inputs

- Structure files: `.pdb`, `.cif`, `.mmcif`
- Archive formats: `.zip`, `.tar`, `.tar.gz`, `.tgz`, `.tar.bz2`, `.tbz2`, `.tar.xz`, `.txz`

### Job discovery rules

- A single structure file is itself a job.
- If the archive root contains multiple first-level subdirectories, each first-level subdirectory is treated as a separate job.
- If the archive root directly contains multiple structure files, the whole archive is treated as one job.
- Structure files inside the same job are processed in natural sort order. Non-structure files are ignored and written to the log.

### Typical outputs for each job

- `report.html`
- `reference-chain-<CHAIN>.pdf`
- `tracks.json`
- `summary.txt`
- `logs.txt`

If PDF export dependencies are unavailable, the program skips PDF generation but still writes HTML, JSON, and text summaries, and marks the job as `partial_success`.

## Key CLI Parameters

| Parameter | Default | Description |
| --- | --- | --- |
| `--input PATH` | None | Input structure file or archive. |
| `--outdir OUTDIR` | None | Output directory. Each job gets its own result subdirectory. |
| `--chain ALL\|CHAIN_ID` | `ALL` | Render all protein chains, or analyze only a specific reference chain. |
| `--columns INT` | `80` | Residue columns shown in each rendered block. |
| `--font-size INT` | `12` | Base font size for reports and panels. |
| `--hyd-window INT` | `3` | Smoothing window size for the hydropathy track. |
| `--msa-db PATH` | None | Local BLAST / MMseqs database prefix, or a protein FASTA file. |
| `--max-homologs-displayed INT` | `5` | Maximum homolog sequences to retrieve and render. Allowed range: `0-25`. |
| `--evalue VALUE` | `1e-6` | BLAST / MMseqs significance threshold. Only fixed enumerated values are accepted. |
| `--disable-msa` | Off | Explicitly disable homolog search, alignment, and conservation calculation. |
| `--keep-temp` | Off | Keep temporary directories generated during extraction, search, and alignment. |
| `--contact-cutoff FLOAT` | `3.7` | Weak-contact distance cutoff in angstroms. |
| `--strong-contact-cutoff FLOAT` | `3.2` | Strong-contact distance cutoff in angstroms. |
| `--verbose` | Off | Enable verbose logging. |

A few commonly used parameters deserve extra attention:

- `--chain`
  The default is `ALL`. The program collects protein chains from the first successfully parsed structure and renders each chain separately. If you only care about one chain, specify it explicitly, for example `--chain B`.
- `--msa-db`
  You can pass either a BLAST / MMseqs database prefix or a protein FASTA file directly. For database download, build steps, and examples, see [blastdb/README.md](/workspaces/re-foldscript/blastdb/README.md).
- `--max-homologs-displayed`
  Controls how many homologs can be retrieved and how many rows are ultimately rendered. If set to `0`, homolog rows are skipped and only the query row is kept.
- `--evalue`
  The currently accepted values are `1e-4`, `1e-5`, `1e-6`, `1e-7`, `1e-8`, `1e-9`, `1e-10`, `1e-11`, and `1e-12`. This threshold controls significance filtering, not the number of returned results.
- `--disable-msa`
  If you only want a quick structural panel, or your environment does not include BLAST / MMseqs / Clustal Omega, you can disable the MSA stage directly.

## Acknowledgements and Reference

`OpenFoldPanel` takes strong inspiration from FoldScript at the product-design level. Thanks to the FoldScript team for showing how clear and useful AI protein-model comparison can be, and for providing a concrete reference point for this project.

- Robert, X., Guillon, C., Gouet, P. (2025). *FoldScript: a web server for the efficient analysis of AI-generated 3D protein models*. *Nucleic Acids Research*, 53(W1), W277-W282. DOI: [10.1093/nar/gkaf326](https://doi.org/10.1093/nar/gkaf326)
