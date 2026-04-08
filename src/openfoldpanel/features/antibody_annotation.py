"""Antibody-chain detection and CDR annotation helpers."""

from __future__ import annotations

from collections.abc import Sequence
from importlib import import_module

from openfoldpanel.constants import SUPPORTED_ANTIBODY_NUMBERING
from openfoldpanel.models import AntibodyAnnotation, RegionAnnotation, SequenceAxisPosition

REGION_ORDER = ("FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4")
CDR_REGION_NAMES = ("CDR1", "CDR2", "CDR3")
NON_ANTIBODY_MARKERS = (
    "variable chain sequence not recognized",
    "not recognized",
)


def annotate_antibody_chain(
    sequence: str,
    sequence_axis: Sequence[SequenceAxisPosition],
    *,
    chain_id: str,
) -> tuple[dict[str, AntibodyAnnotation], list[str]]:
    """Build multi-scheme CDR annotations for an antibody-like reference chain."""

    try:
        chain_class, chain_parse_error = _load_abnumber_types()
    except Exception as exc:
        return {}, [_warning_message(chain_id, scheme, f"abnumber is unavailable ({exc})") for scheme in SUPPORTED_ANTIBODY_NUMBERING]

    annotations: dict[str, AntibodyAnnotation] = {}
    warnings: list[str] = []
    non_antibody_schemes: set[str] = set()
    for scheme in SUPPORTED_ANTIBODY_NUMBERING:
        normalized_scheme = normalize_antibody_scheme(scheme)
        try:
            numbered_chain = chain_class(sequence, scheme=normalized_scheme)
        except chain_parse_error as exc:
            if _looks_like_non_antibody_sequence(exc):
                non_antibody_schemes.add(normalized_scheme)
                continue
            warnings.append(_warning_message(chain_id, normalized_scheme, str(exc)))
            continue
        except Exception as exc:
            warnings.append(_warning_message(chain_id, normalized_scheme, str(exc)))
            continue

        chain_type = _normalize_chain_type(numbered_chain)
        if chain_type is None:
            non_antibody_schemes.add(normalized_scheme)
            continue

        try:
            regions = _extract_cdr_regions(numbered_chain, sequence_axis, normalized_scheme)
        except Exception as exc:
            warnings.append(_warning_message(chain_id, normalized_scheme, str(exc)))
            continue

        if not regions:
            warnings.append(_warning_message(chain_id, normalized_scheme, "no CDR regions were found after numbering"))
            continue

        annotations[normalized_scheme] = AntibodyAnnotation(
            scheme=normalized_scheme,
            chain_type=chain_type,
            regions=regions,
        )

    if annotations:
        return annotations, warnings
    if non_antibody_schemes and len(non_antibody_schemes) == len(SUPPORTED_ANTIBODY_NUMBERING):
        return {}, []
    return {}, warnings


def normalize_antibody_scheme(scheme: str) -> str:
    normalized = scheme.strip().lower()
    if normalized not in SUPPORTED_ANTIBODY_NUMBERING:
        allowed = ", ".join(SUPPORTED_ANTIBODY_NUMBERING)
        raise ValueError(f'Unsupported antibody numbering scheme "{scheme}". Expected one of: {allowed}.')
    return normalized


def antibody_scheme_label(scheme: str) -> str:
    normalized = normalize_antibody_scheme(scheme)
    if normalized == "imgt":
        return "IMGT"
    return normalized.capitalize()


def _warning_message(chain_id: str, scheme: str, detail: str) -> str:
    return (
        f"Antibody numbering for reference chain {chain_id} using {antibody_scheme_label(scheme)} "
        f"was skipped: {detail}."
    )


def _load_abnumber_types():
    abnumber_module = import_module("abnumber")
    exceptions_module = import_module("abnumber.exceptions")
    return abnumber_module.Chain, exceptions_module.ChainParseError


def _looks_like_non_antibody_sequence(exc: BaseException) -> bool:
    message = str(exc).lower()
    return any(marker in message for marker in NON_ANTIBODY_MARKERS)


def _normalize_chain_type(numbered_chain) -> str | None:
    chain_type = str(getattr(numbered_chain, "chain_type", "")).upper()
    if chain_type == "H":
        return "heavy"
    if chain_type in {"K", "L"}:
        return "light"
    return None


def _extract_cdr_regions(
    numbered_chain,
    sequence_axis: Sequence[SequenceAxisPosition],
    scheme: str,
) -> list[RegionAnnotation]:
    regions_mapping = getattr(numbered_chain, "regions", None)
    if not isinstance(regions_mapping, dict):
        raise ValueError("abnumber Chain.regions is unavailable")

    variable_length = 0
    annotations: list[RegionAnnotation] = []
    for region_name in REGION_ORDER:
        region = regions_mapping.get(region_name)
        region_length = len(region) if region is not None else 0
        if region_name in CDR_REGION_NAMES and region_length > 0:
            annotations.append(
                RegionAnnotation(
                    name=region_name,
                    start=variable_length,
                    end=variable_length + region_length,
                    display_label=f"{region_name} - {antibody_scheme_label(scheme)}",
                )
            )
        variable_length += region_length

    if variable_length <= 0:
        raise ValueError("numbered antibody variable region was empty")
    if variable_length > len(sequence_axis):
        raise ValueError(
            f"numbered variable region length {variable_length} exceeds sequence axis length {len(sequence_axis)}"
        )
    return annotations
