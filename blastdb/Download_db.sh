#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
BLASTDB_DIR="${SCRIPT_DIR}"
FORCE=0

usage() {
  cat <<'EOF'
Usage:
  bash blastdb/Download_db.sh [--force] [target ...]

Note:
  This script only downloads raw files or prebuilt BLAST archives.
  If you want to build your own BLAST database prefix from a local sequence file,
  use: bash blastdb/Build_blastdb.sh --input <source_file> --out-prefix <prefix>
  PDB cluster raw files are downloaded into blastdb/pdb_cluster_src/.
  The downstream handling method for PDB cluster data will be updated later.

Targets:
  recommended      Download the recommended set: swissprot-fasta and pdbaa.
  all              Download every supported resource.
  swissprot-fasta  Download UniProt Swiss-Prot FASTA to blastdb/swissprot_fasta/.
  swissprot-blast  Download the NCBI swissprot BLAST database to blastdb/swissprot/.
  pdbaa            Download the NCBI pdbaa BLAST database to blastdb/pdbaa/.
  pdb-clusters     Download optional PDB cluster source files to blastdb/pdb_cluster_src/.

Options:
  --force          Re-download plain files even if they already exist.
  -h, --help       Show this help message.

Examples:
  bash blastdb/Download_db.sh
  bash blastdb/Download_db.sh swissprot-fasta pdbaa
  bash blastdb/Download_db.sh all
EOF
}

log() {
  printf '[blastdb] %s\n' "$*"
}

die() {
  printf '[blastdb] ERROR: %s\n' "$*" >&2
  exit 1
}

has_command() {
  command -v "$1" >/dev/null 2>&1
}

require_command() {
  has_command "$1" || die "Required command not found: $1"
}

download_file() {
  local url="$1"
  local output_path="$2"

  if [[ -f "${output_path}" && "${FORCE}" -eq 0 ]]; then
    log "Skipping existing file: ${output_path}"
    return
  fi

  mkdir -p "$(dirname -- "${output_path}")"

  if has_command curl; then
    log "Downloading ${url}"
    curl -L --fail --retry 3 --output "${output_path}" "${url}"
    return
  fi

  if has_command wget; then
    log "Downloading ${url}"
    wget -O "${output_path}" "${url}"
    return
  fi

  die "Neither curl nor wget is available. Please install one of them first."
}

download_gzip_as_plain() {
  local url="$1"
  local output_path="$2"
  local gzip_path="${output_path}.gz"

  require_command gunzip

  if [[ -f "${output_path}" && "${FORCE}" -eq 0 ]]; then
    log "Skipping existing file: ${output_path}"
    return
  fi

  download_file "${url}" "${gzip_path}"
  log "Extracting ${gzip_path}"
  gunzip -f "${gzip_path}"
}

download_blast_database() {
  local database_name="$1"
  local target_dir="$2"

  require_command update_blastdb.pl
  mkdir -p "${target_dir}"

  log "Downloading BLAST database ${database_name} into ${target_dir}"
  (
    cd "${target_dir}"
    update_blastdb.pl --decompress "${database_name}"
  )
}

download_swissprot_fasta() {
  local target_dir="${BLASTDB_DIR}/swissprot_fasta"
  local output_path="${target_dir}/uniprot_sprot.fasta"

  download_gzip_as_plain \
    "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz" \
    "${output_path}"
}

download_swissprot_blast() {
  download_blast_database "swissprot" "${BLASTDB_DIR}/swissprot"
}

download_pdbaa() {
  download_blast_database "pdbaa" "${BLASTDB_DIR}/pdbaa"
}

download_pdb_clusters() {
  local target_dir="${BLASTDB_DIR}/pdb_cluster_src"
  mkdir -p "${target_dir}"

  download_gzip_as_plain \
    "https://files.rcsb.org/pub/pdb/derived_data/pdb_seqres.txt.gz" \
    "${target_dir}/pdb_seqres.txt"
  download_file \
    "https://cdn.rcsb.org/resources/sequence/clusters/clusters-by-entity-95.txt" \
    "${target_dir}/clusters-by-entity-95.txt"
  download_file \
    "https://cdn.rcsb.org/resources/sequence/clusters/clusters-by-entity-90.txt" \
    "${target_dir}/clusters-by-entity-90.txt"
  download_file \
    "https://cdn.rcsb.org/resources/sequence/clusters/clusters-by-entity-70.txt" \
    "${target_dir}/clusters-by-entity-70.txt"
  download_file \
    "https://cdn.rcsb.org/resources/sequence/clusters/clusters-by-entity-50.txt" \
    "${target_dir}/clusters-by-entity-50.txt"
}

TARGETS=()
while (($# > 0)); do
  case "$1" in
    -h|--help)
      usage
      exit 0
      ;;
    --force)
      FORCE=1
      ;;
    *)
      TARGETS+=("$1")
      ;;
  esac
  shift
done

if ((${#TARGETS[@]} == 0)); then
  TARGETS=("recommended")
fi

declare -a EXPANDED_TARGETS=()
for target in "${TARGETS[@]}"; do
  case "${target}" in
    recommended)
      EXPANDED_TARGETS+=("swissprot-fasta" "pdbaa")
      ;;
    all)
      EXPANDED_TARGETS+=("swissprot-fasta" "swissprot-blast" "pdbaa" "pdb-clusters")
      ;;
    swissprot-fasta|swissprot-blast|pdbaa|pdb-clusters)
      EXPANDED_TARGETS+=("${target}")
      ;;
    *)
      usage >&2
      die "Unknown target: ${target}"
      ;;
  esac
done

declare -A SEEN=()
declare -a ORDERED_TARGETS=()
for target in "${EXPANDED_TARGETS[@]}"; do
  if [[ -z "${SEEN[${target}]:-}" ]]; then
    ORDERED_TARGETS+=("${target}")
    SEEN["${target}"]=1
  fi
done

for target in "${ORDERED_TARGETS[@]}"; do
  case "${target}" in
    swissprot-fasta)
      download_swissprot_fasta
      ;;
    swissprot-blast)
      download_swissprot_blast
      ;;
    pdbaa)
      download_pdbaa
      ;;
    pdb-clusters)
      download_pdb_clusters
      ;;
  esac
done

log "Finished requested downloads."
