#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import time
from pathlib import Path
from typing import NoReturn
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen


SCRIPT_DIR = Path(__file__).resolve().parent
BUILD_DIR = SCRIPT_DIR / "build"
BUILD_BLASTDB_SCRIPT = SCRIPT_DIR / "Build_blastdb.sh"
SUPPORTED_IDENTITIES = (50, 70, 90, 95)
API_URL_TEMPLATE = "https://data.rcsb.org/rest/v1/core/polymer_entity/{entry_id}/{entity_id}"
PDB_SEQRES_FILENAME = "pdb_seqres.txt"
PDB_ENTITY_TOKEN_RE = re.compile(r"^[0-9A-Z]{4}_[0-9]+$")


def log(message: str) -> None:
    print(f"[pdbaa-clusters] {message}")


def warn(message: str) -> None:
    print(f"[pdbaa-clusters] WARNING: {message}", file=sys.stderr)


def die(message: str) -> NoReturn:
    print(f"[pdbaa-clusters] ERROR: {message}", file=sys.stderr)
    raise SystemExit(1)


def _parse_sleep_seconds(value: str) -> float:
    try:
        parsed = float(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("--sleep-seconds must be a number.") from exc
    if parsed < 0:
        raise argparse.ArgumentTypeError("--sleep-seconds must be non-negative.")
    return parsed


def _parse_cluster_dir(value: str) -> Path:
    path = Path(value).expanduser().resolve()
    if not path.exists():
        raise argparse.ArgumentTypeError(f"--pdb-cluster-src does not exist: {path}")
    if not path.is_dir():
        raise argparse.ArgumentTypeError(f"--pdb-cluster-src must be a directory: {path}")
    return path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build PDBAA50/70/90/95 FASTA files and optional BLAST databases from RCSB cluster files."
    )
    selection = parser.add_mutually_exclusive_group(required=True)
    selection.add_argument("--identity", type=int, choices=SUPPORTED_IDENTITIES, help="Cluster identity to build.")
    selection.add_argument("--all", action="store_true", help="Build all supported identities: 50, 70, 90, 95.")
    parser.add_argument(
        "--pdb-cluster-src",
        type=_parse_cluster_dir,
        required=True,
        help="Path to the directory containing pdb_seqres.txt and clusters-by-entity-{50,70,90,95}.txt.",
    )
    parser.add_argument(
        "--build-blastdb",
        action="store_true",
        help="After generating FASTA, call Build_blastdb.sh to build a BLAST database prefix.",
    )
    parser.add_argument("--force", action="store_true", help="Overwrite existing FASTA and BLAST outputs.")
    parser.add_argument(
        "--sleep-seconds",
        type=_parse_sleep_seconds,
        default=0.02,
        help="Delay between RCSB Data API requests.",
    )
    return parser


def cluster_file_path(cluster_dir: Path, identity: int) -> Path:
    return cluster_dir / f"clusters-by-entity-{identity}.txt"


def pdb_seqres_path(cluster_dir: Path) -> Path:
    return cluster_dir / PDB_SEQRES_FILENAME


def fasta_output_path(identity: int) -> Path:
    return BUILD_DIR / f"pdbaa{identity}" / f"pdbaa{identity}.fasta"


def blast_prefix_path(identity: int) -> Path:
    return BUILD_DIR / f"pdbaa{identity}" / f"pdbaa{identity}"


def entity_chain_cache_path() -> Path:
    return BUILD_DIR / "pdbaa_cache" / "entity_chain_map.json"


def iter_representative_entities(path: Path) -> list[str]:
    representatives: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            token = line.split()[0]
            if not PDB_ENTITY_TOKEN_RE.fullmatch(token):
                warn(f"{path.name}:{line_number} representative token was skipped: {token}")
                continue
            representatives.append(token)
    return representatives


def load_pdb_protein_chain_sequences(path: Path) -> dict[tuple[str, str], str]:
    sequences: dict[tuple[str, str], str] = {}
    header: str | None = None
    sequence_chunks: list[str] = []

    def commit_record(raw_header: str | None, chunks: list[str]) -> None:
        if raw_header is None:
            return
        header_text = raw_header.strip()
        if "mol:protein" not in header_text.lower():
            return

        token = header_text.split(None, 1)[0]
        if "_" not in token:
            warn(f"{path.name} has an unexpected FASTA header and it was skipped: {header_text}")
            return

        entry_id, chain_id = token.split("_", 1)
        sequence = re.sub(r"\s+", "", "".join(chunks))
        if not entry_id or not chain_id or not sequence:
            warn(f"{path.name} has an incomplete protein FASTA record and it was skipped: {header_text}")
            return

        sequences[(entry_id.upper(), chain_id)] = sequence

    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith(">"):
                commit_record(header, sequence_chunks)
                header = line[1:]
                sequence_chunks = []
                continue
            if header is not None:
                sequence_chunks.append(line.strip())

    commit_record(header, sequence_chunks)
    return sequences


def load_entity_chain_cache(path: Path) -> dict[str, dict[str, str]]:
    if not path.exists():
        return {}

    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        die(f"Entity-chain cache is not valid JSON: {path} ({exc})")

    if not isinstance(payload, dict):
        die(f"Entity-chain cache must contain a JSON object: {path}")

    normalized: dict[str, dict[str, str]] = {}
    for raw_key, raw_value in payload.items():
        if not isinstance(raw_key, str) or not isinstance(raw_value, dict):
            warn(f"Malformed entity-chain cache entry was skipped: {raw_key}")
            continue

        entry_id = raw_value.get("entry_id")
        entity_id = raw_value.get("entity_id")
        chain_id = raw_value.get("chain_id")
        if not all(isinstance(item, str) and item.strip() for item in (entry_id, entity_id, chain_id)):
            warn(f"Incomplete entity-chain cache entry was skipped: {raw_key}")
            continue

        cache_key = raw_key.upper()
        normalized[cache_key] = {
            "entry_id": entry_id.upper(),
            "entity_id": entity_id,
            "chain_id": chain_id.strip(),
        }

    return normalized


def save_entity_chain_cache(path: Path, cache: dict[str, dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp_path = path.with_suffix(".json.tmp")
    serialized = json.dumps(cache, indent=2, sort_keys=True)
    temp_path.write_text(serialized + "\n", encoding="utf-8")
    temp_path.replace(path)


def _first_chain_id(payload: dict, entity_token: str) -> str:
    container = payload.get("rcsb_polymer_entity_container_identifiers") or {}
    for key in ("auth_asym_ids", "asym_ids"):
        value = container.get(key)
        if isinstance(value, list):
            for item in value:
                if isinstance(item, str) and item.strip():
                    return item.strip()
        elif isinstance(value, str) and value.strip():
            return value.strip().split(",")[0].strip()
    raise ValueError(f"no chain identifier was found for {entity_token}")


def fetch_polymer_entity_chain_mapping(entity_token: str) -> dict[str, str]:
    entry_id, entity_id = entity_token.rsplit("_", 1)
    request = Request(
        API_URL_TEMPLATE.format(entry_id=entry_id, entity_id=entity_id),
        headers={"User-Agent": "OpenFoldPanel-blastdb/1.0"},
    )
    with urlopen(request, timeout=30) as response:
        payload = json.load(response)

    return {
        "entry_id": entry_id.upper(),
        "entity_id": entity_id,
        "chain_id": _first_chain_id(payload, entity_token),
    }


def resolve_entity_chain_mapping(
    entity_token: str,
    cache: dict[str, dict[str, str]],
    cache_path: Path,
) -> tuple[dict[str, str], bool]:
    cache_key = entity_token.upper()
    cached = cache.get(cache_key)
    if cached is not None:
        return cached, False

    mapping = fetch_polymer_entity_chain_mapping(cache_key)
    cache[cache_key] = mapping
    save_entity_chain_cache(cache_path, cache)
    return mapping, True


def write_fasta_record(handle, header: str, sequence: str) -> None:
    handle.write(f">{header}\n")
    for start in range(0, len(sequence), 80):
        handle.write(sequence[start : start + 80] + "\n")


def load_existing_fasta_headers(path: Path) -> set[str]:
    if not path.exists():
        return set()

    headers: set[str] = set()
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if line.startswith(">") and len(line) > 1:
                headers.add(line[1:])
    return headers


def build_identity_fasta(identity: int, cluster_dir: Path, *, sleep_seconds: float, force: bool) -> tuple[Path, Path]:
    cluster_path = cluster_file_path(cluster_dir, identity)
    seqres_path = pdb_seqres_path(cluster_dir)
    fasta_path = fasta_output_path(identity)
    blast_prefix = blast_prefix_path(identity)
    cache_path = entity_chain_cache_path()

    if not cluster_path.exists():
        die(f"Cluster file does not exist: {cluster_path}")
    if not seqres_path.exists():
        die(f"PDB SEQRES file does not exist: {seqres_path}")

    if fasta_path.exists() and not force:
        die(f"FASTA output already exists. Re-run with --force to rebuild: {fasta_path}")

    fasta_path.parent.mkdir(parents=True, exist_ok=True)
    temp_path = fasta_path.with_suffix(".fasta.tmp")
    if force and temp_path.exists():
        log(f"PDBAA{identity}: removing existing temporary FASTA to rebuild from scratch")
        temp_path.unlink()

    representatives = iter_representative_entities(cluster_path)
    if not representatives:
        die(f"No representative entities were found in: {cluster_path}")

    protein_sequences = load_pdb_protein_chain_sequences(seqres_path)
    cache = load_entity_chain_cache(cache_path)
    existing_headers = load_existing_fasta_headers(temp_path)
    resumed_records = len(existing_headers)
    if resumed_records:
        log(f"PDBAA{identity}: resuming from existing temporary FASTA with {resumed_records} record(s)")

    written = 0
    skipped = 0
    open_mode = "a" if temp_path.exists() else "w"
    with temp_path.open(open_mode, encoding="utf-8") as handle:
        for index, entity_token in enumerate(representatives, start=1):
            try:
                mapping, fetched = resolve_entity_chain_mapping(entity_token, cache, cache_path)
            except (HTTPError, URLError, TimeoutError, ValueError, json.JSONDecodeError) as exc:
                die(f"{entity_token} chain mapping could not be fetched: {exc}")

            header = f"pdb|{mapping['entry_id']}|{mapping['chain_id']}"
            if header in existing_headers:
                continue

            sequence = protein_sequences.get((mapping["entry_id"], mapping["chain_id"]))
            if sequence is None:
                warn(
                    f"{entity_token} mapped to {mapping['entry_id']} chain {mapping['chain_id']}, "
                    "but no local protein chain sequence was found; it was skipped."
                )
                skipped += 1
                continue

            write_fasta_record(handle, header, sequence)
            existing_headers.add(header)
            written += 1

            if fetched and sleep_seconds > 0 and index != len(representatives):
                time.sleep(sleep_seconds)

    total_written = len(existing_headers)
    if total_written == 0:
        temp_path.unlink(missing_ok=True)
        die(f"No FASTA records were written for PDBAA{identity}.")

    temp_path.replace(fasta_path)
    message = f"PDBAA{identity}: wrote {written} new representative sequence(s) to {fasta_path}; total {total_written}"
    if resumed_records:
        message += f"; resumed {resumed_records}"
    if skipped:
        message += f"; skipped {skipped}"
    log(message)
    return fasta_path, blast_prefix


def build_blast_database(identity: int, fasta_path: Path, blast_prefix: Path, *, force: bool) -> None:
    command = [
        "bash",
        str(BUILD_BLASTDB_SCRIPT),
        "--input",
        str(fasta_path),
        "--out-prefix",
        str(blast_prefix),
        "--title",
        f"PDBAA{identity}",
    ]
    if force:
        command.append("--force")

    log(f"PDBAA{identity}: building BLAST database prefix {blast_prefix}")
    subprocess.run(command, check=True)


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    identities = list(SUPPORTED_IDENTITIES) if args.all else [args.identity]
    built_prefixes: list[Path] = []

    for identity in identities:
        fasta_path, blast_prefix = build_identity_fasta(
            identity,
            args.pdb_cluster_src,
            sleep_seconds=args.sleep_seconds,
            force=args.force,
        )
        if args.build_blastdb:
            build_blast_database(identity, fasta_path, blast_prefix, force=args.force)
            built_prefixes.append(blast_prefix)

    if args.build_blastdb and built_prefixes:
        log("Generated BLAST database prefixes:")
        for prefix in built_prefixes:
            log(f"  --msa-db {prefix}")
    else:
        log("FASTA generation completed. Use Build_blastdb.sh if you want BLAST database prefixes.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
