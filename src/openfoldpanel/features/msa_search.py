"""Optional homolog search wrappers for BLAST+ or MMseqs2."""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

from openfoldpanel.constants import validate_evalue
from openfoldpanel.utils.subprocess import ExternalToolError, MissingExecutableError, run_command, which

_UNIPROT_PATTERN = re.compile(r"((?:sp|tr)\|[^|;\s]+\|[^|;\s]+)", re.IGNORECASE)
_UNIPROT_ACCESSION_PATTERN = re.compile(r"^(?:[A-Z0-9]{6}|[A-Z0-9]{10})(?:\.[0-9]+)?$", re.IGNORECASE)
_PDB_PATTERN = re.compile(r"(pdb\|[^|;\s]+\|[^|;\s]+)", re.IGNORECASE)
_FASTA_SUFFIXES = {".fasta", ".fa", ".faa", ".fsa", ".fas"}


@dataclass(slots=True)
class PreparedDatabase:
    search_path: Path
    fasta_sequence_index: "FastaSequenceIndex | None" = None
    uniprot_header_map: dict[str, str] | None = None
    preparation_warnings: list[str] | None = None


@dataclass(slots=True)
class FastaSequenceIndex:
    fasta_path: Path
    token_offsets: dict[str, int]
    supported_identifier_tokens: dict[str, str]
    accession_tokens: dict[str, str]


def search_homologs(
    query_fasta: Path,
    database: Path,
    *,
    max_homologs_displayed: int,
    evalue: str,
    workdir: Path,
    logger: logging.Logger,
) -> tuple[list[tuple[str, str]], list[str]]:
    """Search a local sequence database for homologs."""

    validated_evalue = validate_evalue(evalue)
    preparation_warnings: list[str] = []
    prepared_blast_database: PreparedDatabase | None = None
    fasta_sequence_index, uniprot_header_map = _load_fasta_lookup_support(database)
    if which("blastp"):
        prepared_blast_database = _prepare_blast_database(
            database,
            workdir,
            logger,
            fasta_sequence_index=fasta_sequence_index,
            uniprot_header_map=uniprot_header_map,
        )
        preparation_warnings.extend(prepared_blast_database.preparation_warnings or [])
        if _can_use_blastp(database, prepared_blast_database):
            logger.info("MSA search: using blastp")
            logger.info("BLAST evalue threshold: %s", validated_evalue)
            rows, warnings = _search_with_blastp(
                query_fasta,
                prepared_blast_database.search_path,
                max_homologs_displayed=max_homologs_displayed,
                evalue=validated_evalue,
                logger=logger,
                fasta_sequence_index=prepared_blast_database.fasta_sequence_index,
                uniprot_header_map=prepared_blast_database.uniprot_header_map,
            )
            return rows, [*preparation_warnings, *warnings]

    if which("mmseqs"):
        logger.info("MSA search: using mmseqs easy-search")
        logger.info("MMseqs evalue threshold: %s", validated_evalue)
        rows, warnings = _search_with_mmseqs(
            query_fasta,
            database,
            max_homologs_displayed=max_homologs_displayed,
            evalue=validated_evalue,
            workdir=workdir,
            logger=logger,
            fasta_sequence_index=fasta_sequence_index,
            uniprot_header_map=uniprot_header_map,
        )
        return rows, [*preparation_warnings, *warnings]

    if preparation_warnings:
        return [], preparation_warnings
    return [], ["Neither blastp nor mmseqs was found; MSA search was skipped."]


def _is_fasta_database(database: Path) -> bool:
    return database.suffix.casefold() in _FASTA_SUFFIXES


def _load_fasta_lookup_support(
    database: Path,
) -> tuple[FastaSequenceIndex | None, dict[str, str] | None]:
    if not _is_fasta_database(database):
        return None, None
    return _load_fasta_sequence_index(database), _load_uniprot_header_map_from_fasta(database)


def _can_use_blastp(database: Path, prepared_database: PreparedDatabase) -> bool:
    if not _is_fasta_database(database):
        return True
    return prepared_database.search_path != database


def _prepare_blast_database(
    database: Path,
    workdir: Path,
    logger: logging.Logger,
    *,
    fasta_sequence_index: FastaSequenceIndex | None,
    uniprot_header_map: dict[str, str] | None,
) -> PreparedDatabase:
    if not _is_fasta_database(database):
        return PreparedDatabase(search_path=database, preparation_warnings=[])

    if which("makeblastdb") is None:
        warning = "makeblastdb was not found; FASTA MSA database input could not be prepared for blastp."
        logger.warning(warning)
        return PreparedDatabase(
            search_path=database,
            fasta_sequence_index=fasta_sequence_index,
            uniprot_header_map=uniprot_header_map,
            preparation_warnings=[warning],
        )

    blastdb_dir = workdir.parent / "_msa_blastdb"
    blastdb_dir.mkdir(parents=True, exist_ok=True)
    prefix = blastdb_dir / database.stem
    required_files = [Path(f"{prefix}.{suffix}") for suffix in ("pin", "phr", "psq")]
    if not all(path.exists() for path in required_files):
        logger.info("MSA search: building temporary BLAST database from FASTA %s", database)
        try:
            run_command(
                [
                    "makeblastdb",
                    "-in",
                    str(database),
                    "-dbtype",
                    "prot",
                    "-parse_seqids",
                    "-out",
                    str(prefix),
                ]
            )
        except (MissingExecutableError, ExternalToolError) as exc:
            warning = f"makeblastdb failed for FASTA MSA database {database}: {exc}"
            logger.warning(warning)
            return PreparedDatabase(
                search_path=database,
                fasta_sequence_index=fasta_sequence_index,
                uniprot_header_map=uniprot_header_map,
                preparation_warnings=[warning],
            )
    return PreparedDatabase(
        search_path=prefix,
        fasta_sequence_index=fasta_sequence_index,
        uniprot_header_map=uniprot_header_map,
        preparation_warnings=[],
    )


@lru_cache(maxsize=8)
def _load_uniprot_header_map_from_fasta(fasta_path: Path) -> dict[str, str]:
    fasta_index = _load_fasta_sequence_index(fasta_path)
    if fasta_index is None:
        return {}
    return dict(fasta_index.accession_tokens)


@lru_cache(maxsize=8)
def _load_fasta_sequence_index(fasta_path: Path) -> FastaSequenceIndex | None:
    if not fasta_path.is_file():
        return None

    token_offsets: dict[str, int] = {}
    supported_identifier_tokens: dict[str, str] = {}
    accession_tokens: dict[str, str] = {}
    with fasta_path.open("rb") as handle:
        while True:
            offset = handle.tell()
            line = handle.readline()
            if not line:
                break
            if not line.startswith(b">"):
                continue
            header = line[1:].decode("utf-8", errors="replace").strip()
            token = header.split(maxsplit=1)[0]
            if not token:
                continue
            token_offsets.setdefault(token, offset)
            supported = _extract_supported_identifier(token) or _extract_supported_identifier(header)
            if supported:
                supported_identifier_tokens.setdefault(supported, token)
            accession = _normalize_uniprot_accession(token) or _normalize_uniprot_accession(header)
            if accession:
                accession_tokens.setdefault(accession, token)

    return FastaSequenceIndex(
        fasta_path=fasta_path,
        token_offsets=token_offsets,
        supported_identifier_tokens=supported_identifier_tokens,
        accession_tokens=accession_tokens,
    )


def _extract_supported_identifier(header: str) -> str | None:
    """Extract a supported structured identifier from a header-like string."""

    normalized = header.strip().lstrip(">")
    if not normalized:
        return None
    for pattern in (_UNIPROT_PATTERN, _PDB_PATTERN):
        match = pattern.search(normalized)
        if match:
            return match.group(1)
    return None


def _normalize_uniprot_accession(text: str) -> str | None:
    normalized = text.strip().lstrip(">")
    if not normalized:
        return None

    token = _extract_supported_identifier(normalized)
    accession: str | None = None
    if token is not None and token.casefold().startswith(("sp|", "tr|")):
        accession = token.split("|", 2)[1]
    else:
        compact = normalized.split()[0].split(";", 1)[0]
        parts = compact.split("|")
        if len(parts) >= 2 and parts[0].casefold() in {"sp", "tr"}:
            accession = parts[1]
        elif "|" not in compact:
            accession = compact

    if not accession or _UNIPROT_ACCESSION_PATTERN.fullmatch(accession) is None:
        return None
    return accession.split(".", 1)[0] or None


def _blast_lookup_candidates(raw_identifier: str) -> list[str]:
    """Build a small set of lookup candidates for blastdbcmd recovery."""

    normalized = raw_identifier.strip().lstrip(">")
    if not normalized:
        return []

    candidates: list[str] = []

    def register(candidate: str) -> None:
        if candidate and candidate not in candidates:
            candidates.append(candidate)

    register(normalized)
    register(normalized.rstrip("|"))

    prefix, separator, remainder = normalized.partition("|")
    if separator and prefix.lower() in {"sp", "tr"}:
        accession = remainder.partition("|")[0]
        if accession:
            register(accession)
            accession_base = accession.split(".", 1)[0]
            register(f"{prefix}|{accession_base}|")
            register(f"{prefix}|{accession_base}")
            register(accession_base)

    return candidates


def _lookup_text_candidates(*texts: str) -> list[str]:
    candidates: list[str] = []

    def register(candidate: str) -> None:
        normalized = candidate.strip().lstrip(">")
        if normalized and normalized not in candidates:
            candidates.append(normalized)

    for text in texts:
        normalized = text.strip().lstrip(">")
        if not normalized:
            continue
        register(normalized)
        token = normalized.split(maxsplit=1)[0]
        register(token)
        for fragment in token.split(";"):
            register(fragment)
    return candidates


def _read_fasta_payload(text: str) -> tuple[str | None, str]:
    header_token: str | None = None
    sequence_parts: list[str] = []
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith(">"):
            if header_token is not None:
                break
            header_token = stripped[1:].strip().split(maxsplit=1)[0] or None
            continue
        if header_token is not None:
            sequence_parts.append(stripped)
    return header_token, "".join(sequence_parts)


def _read_indexed_fasta_sequence(index: FastaSequenceIndex, token: str) -> str | None:
    offset = index.token_offsets.get(token)
    if offset is None:
        return None

    sequence_parts: list[str] = []
    with index.fasta_path.open("rb") as handle:
        handle.seek(offset)
        header_line = handle.readline()
        if not header_line.startswith(b">"):
            return None
        while True:
            line = handle.readline()
            if not line or line.startswith(b">"):
                break
            sequence_parts.append(line.decode("utf-8", errors="replace").strip())
    sequence = "".join(sequence_parts)
    return sequence or None


def _recover_full_sequence_from_fasta_index(
    fasta_sequence_index: FastaSequenceIndex,
    *,
    raw_identifier: str,
    resolved_identifier: str,
    candidates: tuple[str, ...],
) -> tuple[str, str] | None:
    token: str | None = None
    for candidate in _lookup_text_candidates(resolved_identifier, raw_identifier, *candidates):
        if candidate in fasta_sequence_index.token_offsets:
            token = candidate
            break
        supported = _extract_supported_identifier(candidate)
        if supported and supported in fasta_sequence_index.supported_identifier_tokens:
            token = fasta_sequence_index.supported_identifier_tokens[supported]
            break
        accession = _normalize_uniprot_accession(candidate)
        if accession and accession in fasta_sequence_index.accession_tokens:
            token = fasta_sequence_index.accession_tokens[accession]
            break

    if token is None:
        return None

    sequence = _read_indexed_fasta_sequence(fasta_sequence_index, token)
    if sequence is None:
        return None
    return token, sequence


def _recover_full_sequence_from_blastdb(
    database: Path,
    *,
    raw_identifier: str,
    resolved_identifier: str,
    candidates: tuple[str, ...],
    logger: logging.Logger,
) -> str | None:
    if which("blastdbcmd") is None:
        return None

    lookup_candidates: list[str] = []
    for candidate in _lookup_text_candidates(resolved_identifier, raw_identifier, *candidates):
        for lookup in (candidate, * _blast_lookup_candidates(candidate)):
            if lookup not in lookup_candidates:
                lookup_candidates.append(lookup)

    for candidate in lookup_candidates:
        try:
            result = run_command(
                [
                    "blastdbcmd",
                    "-db",
                    str(database),
                    "-entry",
                    candidate,
                    "-outfmt",
                    "%f",
                ]
            )
        except (MissingExecutableError, ExternalToolError) as exc:
            logger.debug("blastdbcmd sequence recovery failed for %s via %s: %s", raw_identifier, candidate, exc)
            continue

        _, sequence = _read_fasta_payload(result.stdout)
        if sequence:
            return sequence
    return None


def _recover_blast_identifier(database: Path, raw_identifier: str, logger: logging.Logger) -> str | None:
    """Recover a supported structured identifier from the BLAST database header."""

    if which("blastdbcmd") is None:
        return None

    for candidate in _blast_lookup_candidates(raw_identifier):
        try:
            result = run_command(
                [
                    "blastdbcmd",
                    "-db",
                    str(database),
                    "-entry",
                    candidate,
                    "-outfmt",
                    "%f",
                ]
            )
        except (MissingExecutableError, ExternalToolError) as exc:
            logger.debug("blastdbcmd header recovery failed for %s via %s: %s", raw_identifier, candidate, exc)
            continue

        first_line = next((line for line in result.stdout.splitlines() if line.startswith(">")), "")
        recovered = _extract_supported_identifier(first_line)
        if recovered:
            return recovered
    return None


def _resolve_full_sequence(
    *,
    raw_identifier: str,
    resolved_identifier: str,
    candidates: tuple[str, ...],
    aligned_fragment: str,
    fasta_sequence_index: FastaSequenceIndex | None,
    database: Path | None,
    logger: logging.Logger,
) -> tuple[str, str]:
    if fasta_sequence_index is not None:
        recovered = _recover_full_sequence_from_fasta_index(
            fasta_sequence_index,
            raw_identifier=raw_identifier,
            resolved_identifier=resolved_identifier,
            candidates=candidates,
        )
        if recovered is not None:
            _, sequence = recovered
            return sequence, "fasta_index"

    if database is not None:
        recovered = _recover_full_sequence_from_blastdb(
            database,
            raw_identifier=raw_identifier,
            resolved_identifier=resolved_identifier,
            candidates=candidates,
            logger=logger,
        )
        if recovered is not None:
            return recovered, "blastdbcmd"

    return aligned_fragment.replace("-", ""), "aligned_fragment"


def _resolve_identifier(
    *,
    raw_identifier: str,
    candidates: tuple[str, ...],
    database: Path | None,
    logger: logging.Logger,
    uniprot_header_map: dict[str, str] | None = None,
) -> tuple[str, bool]:
    """Resolve a final display identifier from supported database header formats."""

    if uniprot_header_map:
        for candidate in (*candidates, raw_identifier):
            accession = _normalize_uniprot_accession(candidate)
            if accession and accession in uniprot_header_map:
                return uniprot_header_map[accession], True

    for candidate in (*candidates, raw_identifier):
        resolved = _extract_supported_identifier(candidate)
        if resolved:
            return resolved, True

    if database is not None:
        recovered = _recover_blast_identifier(database, raw_identifier, logger)
        if recovered:
            return recovered, True

    return raw_identifier, False


def _search_with_blastp(
    query_fasta: Path,
    database: Path,
    *,
    max_homologs_displayed: int,
    evalue: str,
    logger: logging.Logger,
    fasta_sequence_index: FastaSequenceIndex | None = None,
    uniprot_header_map: dict[str, str] | None = None,
) -> tuple[list[tuple[str, str]], list[str]]:
    try:
        result = run_command(
            [
                "blastp",
                "-query",
                str(query_fasta),
                "-db",
                str(database),
                "-max_target_seqs",
                str(max_homologs_displayed),
                "-evalue",
                evalue,
                "-outfmt",
                "6 sseqid sallseqid stitle sseq",
            ]
        )
    except (MissingExecutableError, ExternalToolError) as exc:
        logger.warning("blastp search failed: %s", exc)
        return [], ["blastp search failed; MSA track was skipped."]

    rows: list[tuple[str, str]] = []
    identifier_fallback_count = 0
    full_length_from_fasta_count = 0
    full_length_from_blastdb_count = 0
    fragment_fallback_count = 0
    for line in result.stdout.splitlines():
        if not line.strip():
            continue
        parts = line.strip().split("\t", 3)
        if len(parts) == 4:
            raw_identifier, all_identifiers, header, sequence = parts
        elif len(parts) == 3:
            raw_identifier, header, sequence = parts
            all_identifiers = ""
        else:
            continue
        identifier, recovered = _resolve_identifier(
            raw_identifier=raw_identifier,
            candidates=(all_identifiers, header),
            database=database,
            logger=logger,
            uniprot_header_map=uniprot_header_map,
        )
        if not recovered:
            identifier_fallback_count += 1
        full_sequence, recovery_source = _resolve_full_sequence(
            raw_identifier=raw_identifier,
            resolved_identifier=identifier,
            candidates=(all_identifiers, header),
            aligned_fragment=sequence,
            fasta_sequence_index=fasta_sequence_index,
            database=database,
            logger=logger,
        )
        if recovery_source == "fasta_index":
            full_length_from_fasta_count += 1
        elif recovery_source == "blastdbcmd":
            full_length_from_blastdb_count += 1
        else:
            fragment_fallback_count += 1
        rows.append((identifier, full_sequence))

    warnings: list[str] = []
    if full_length_from_fasta_count:
        logger.info(
            "blastp full-length sequence recovery: resolved %s homolog hit(s) from FASTA index",
            full_length_from_fasta_count,
        )
    if full_length_from_blastdb_count:
        logger.info(
            "blastp full-length sequence recovery: resolved %s homolog hit(s) via blastdbcmd",
            full_length_from_blastdb_count,
        )
    if identifier_fallback_count:
        message = (
            f"blastp header recovery failed for {identifier_fallback_count} homolog hit(s); raw identifiers were used as fallback."
        )
        logger.warning(message)
        warnings.append(message)
    if fragment_fallback_count:
        message = (
            f"blastp full-length sequence recovery fell back to aligned fragments for {fragment_fallback_count} homolog hit(s)."
        )
        logger.warning(message)
        warnings.append(message)
    return rows, warnings


def _search_with_mmseqs(
    query_fasta: Path,
    database: Path,
    *,
    max_homologs_displayed: int,
    evalue: str,
    workdir: Path,
    logger: logging.Logger,
    fasta_sequence_index: FastaSequenceIndex | None = None,
    uniprot_header_map: dict[str, str] | None = None,
) -> tuple[list[tuple[str, str]], list[str]]:
    output_tsv = workdir / "mmseqs_hits.tsv"
    tmp_dir = workdir / "mmseqs_tmp"
    try:
        run_command(
            [
                "mmseqs",
                "easy-search",
                str(query_fasta),
                str(database),
                str(output_tsv),
                str(tmp_dir),
                "--format-output",
                "target,theader,alntseq",
                "--max-seqs",
                str(max_homologs_displayed),
                "-e",
                evalue,
            ]
        )
    except (MissingExecutableError, ExternalToolError) as exc:
        logger.warning("mmseqs search failed: %s", exc)
        return [], ["mmseqs search failed; MSA track was skipped."]

    rows: list[tuple[str, str]] = []
    identifier_fallback_count = 0
    full_length_from_fasta_count = 0
    fragment_fallback_count = 0
    if output_tsv.exists():
        for line in output_tsv.read_text(encoding="utf-8").splitlines():
            if not line.strip():
                continue
            parts = line.split("\t", 2)
            if len(parts) < 3:
                continue
            raw_identifier, header, sequence = parts
            identifier, recovered = _resolve_identifier(
                raw_identifier=raw_identifier,
                candidates=(header,),
                database=None,
                logger=logger,
                uniprot_header_map=uniprot_header_map,
            )
            if not recovered:
                identifier_fallback_count += 1
            full_sequence, recovery_source = _resolve_full_sequence(
                raw_identifier=raw_identifier,
                resolved_identifier=identifier,
                candidates=(header,),
                aligned_fragment=sequence,
                fasta_sequence_index=fasta_sequence_index,
                database=None,
                logger=logger,
            )
            if recovery_source == "fasta_index":
                full_length_from_fasta_count += 1
            else:
                fragment_fallback_count += 1
            rows.append((identifier, full_sequence))

    warnings: list[str] = []
    if full_length_from_fasta_count:
        logger.info(
            "mmseqs full-length sequence recovery: resolved %s homolog hit(s) from FASTA index",
            full_length_from_fasta_count,
        )
    if identifier_fallback_count:
        message = (
            f"mmseqs header recovery failed for {identifier_fallback_count} homolog hit(s); raw identifiers were used as fallback."
        )
        logger.warning(message)
        warnings.append(message)
    if fragment_fallback_count:
        message = (
            f"mmseqs full-length sequence recovery fell back to aligned fragments for {fragment_fallback_count} homolog hit(s)."
        )
        logger.warning(message)
        warnings.append(message)
    return rows, warnings
