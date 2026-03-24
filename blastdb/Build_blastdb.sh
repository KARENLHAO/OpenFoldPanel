#!/usr/bin/env bash

set -euo pipefail

INPUT=""
OUT_PREFIX=""
TITLE=""
FORCE=0

usage() {
  cat <<'EOF'
Usage:
  bash blastdb/Build_blastdb.sh --input <source_file> --out-prefix <prefix> [options]

Required:
  --input <source_file>   Local sequence file to build into a protein BLAST database.
  --out-prefix <prefix>   Output BLAST database prefix, for example:
                          ./blastdb/build/swissprot/swissprot

Options:
  --title <db_title>      Optional BLAST database title.
  --force                 Remove existing files under the same prefix before rebuilding.
  -h, --help              Show this help message.

Examples:
  bash blastdb/Build_blastdb.sh \
    --input ./blastdb/swissprot_fasta/uniprot_sprot.fasta \
    --out-prefix ./blastdb/build/swissprot/swissprot \
    --title SWISSPROT

  bash blastdb/Build_blastdb.sh \
    --input ./blastdb/pdb_cluster_src/pdb_seqres.txt \
    --out-prefix ./blastdb/build/pdb_seqres/pdb_seqres \
    --title PDB_SEQRES
EOF
}

log() {
  printf '[build-blastdb] %s\n' "$*"
}

die() {
  printf '[build-blastdb] ERROR: %s\n' "$*" >&2
  exit 1
}

has_command() {
  command -v "$1" >/dev/null 2>&1
}

require_command() {
  has_command "$1" || die "Required command not found: $1"
}

while (($# > 0)); do
  case "$1" in
    --input)
      shift
      (($# > 0)) || die "Missing value for --input"
      INPUT="$1"
      ;;
    --out-prefix)
      shift
      (($# > 0)) || die "Missing value for --out-prefix"
      OUT_PREFIX="$1"
      ;;
    --title)
      shift
      (($# > 0)) || die "Missing value for --title"
      TITLE="$1"
      ;;
    --force)
      FORCE=1
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      usage >&2
      die "Unknown argument: $1"
      ;;
  esac
  shift
done

[[ -n "${INPUT}" ]] || {
  usage >&2
  die "--input is required"
}

[[ -n "${OUT_PREFIX}" ]] || {
  usage >&2
  die "--out-prefix is required"
}

require_command makeblastdb

if [[ ! -f "${INPUT}" ]]; then
  die "Input file does not exist: ${INPUT}"
fi

OUT_DIR="$(dirname -- "${OUT_PREFIX}")"
mkdir -p "${OUT_DIR}"

existing_outputs=(
  "${OUT_PREFIX}.pdb"
  "${OUT_PREFIX}.phr"
  "${OUT_PREFIX}.pin"
  "${OUT_PREFIX}.pjs"
  "${OUT_PREFIX}.pog"
  "${OUT_PREFIX}.pos"
  "${OUT_PREFIX}.pot"
  "${OUT_PREFIX}.ppd"
  "${OUT_PREFIX}.ppi"
  "${OUT_PREFIX}.psq"
  "${OUT_PREFIX}.ptf"
  "${OUT_PREFIX}.pto"
)

has_existing_output=0
for path in "${existing_outputs[@]}"; do
  if [[ -e "${path}" ]]; then
    has_existing_output=1
    break
  fi
done

if [[ "${has_existing_output}" -eq 1 && "${FORCE}" -ne 1 ]]; then
  die "Output prefix already exists. Re-run with --force to rebuild: ${OUT_PREFIX}"
fi

if [[ "${has_existing_output}" -eq 1 && "${FORCE}" -eq 1 ]]; then
  log "Removing existing BLAST database files under prefix ${OUT_PREFIX}"
  rm -f "${existing_outputs[@]}"
fi

command=(
  makeblastdb
  -in "${INPUT}"
  -dbtype prot
  -parse_seqids
  -out "${OUT_PREFIX}"
)

if [[ -n "${TITLE}" ]]; then
  command+=(-title "${TITLE}")
fi

log "Building protein BLAST database"
log "Input: ${INPUT}"
log "Output prefix: ${OUT_PREFIX}"

"${command[@]}"

for required in "${OUT_PREFIX}.pin" "${OUT_PREFIX}.phr" "${OUT_PREFIX}.psq"; do
  [[ -f "${required}" ]] || die "BLAST database build appears incomplete; missing file: ${required}"
done

log "Build completed successfully."
log "Use this prefix with OpenFoldPanel:"
log "  --msa-db ${OUT_PREFIX}"
