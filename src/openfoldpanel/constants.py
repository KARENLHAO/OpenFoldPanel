"""Project-wide constants for openfoldpanel."""

from __future__ import annotations

SUPPORTED_STRUCTURE_SUFFIXES = {".pdb", ".cif", ".mmcif"}
SUPPORTED_ARCHIVE_SUFFIXES = {
    ".zip",
    ".tar",
    ".tgz",
    ".tbz2",
    ".txz",
    ".tar.gz",
    ".tar.bz2",
    ".tar.xz",
}

DEFAULT_FONT_FAMILY = '"WenQuanYi Zen Hei", "Noto Sans CJK SC", "Microsoft YaHei", "Liberation Sans", "Nimbus Sans", sans-serif'
DEFAULT_HEADING_FONT_FAMILY = '"Noto Serif CJK SC", "Source Han Serif SC", "Songti SC", "STSong", "WenQuanYi Zen Hei", "Liberation Serif", "Nimbus Roman", serif'
DEFAULT_COLUMNS = 80
DEFAULT_FONT_SIZE = 12
DEFAULT_HYDROPATHY_WINDOW = 3
DEFAULT_MAX_HOMOLOGS_DISPLAYED = 5
DEFAULT_EVALUE = "1e-6"  # BLAST hit significance threshold
MAX_HOMOLOGS_DISPLAYED_LIMIT = 25
DEFAULT_CONTACT_CUTOFF = 3.7
DEFAULT_STRONG_CONTACT_CUTOFF = 3.2

ALLOWED_EVALUES = ["1e-4", "1e-5", "1e-6", "1e-7", "1e-8", "1e-9", "1e-10", "1e-11", "1e-12"]


def validate_evalue(value: str, *, parameter_name: str = "evalue") -> str:
    """Validate the fixed BLAST/MMseqs hit significance threshold enum."""

    if value not in ALLOWED_EVALUES:
        allowed = ", ".join(ALLOWED_EVALUES)
        raise ValueError(f"{parameter_name} must be one of: {allowed}.")
    return value

DEFAULT_CELL_WIDTH_RATIO = 0.68
DEFAULT_ROW_HEIGHT_RATIO = 1.6
DEFAULT_LABEL_COLUMNS = 18
DEFAULT_MARGIN = 18
DEFAULT_BLOCK_GAP = 24
DEFAULT_TICK_HEIGHT = 18

ACCESSIBILITY_THRESHOLDS = {
    "buried": (None, 0.1),
    "intermediate": (0.1, 0.4),
    "accessible": (0.4, 1.0),
    "highly_exposed": (1.0, None),
}

HYDROPATHY_THRESHOLDS = {
    "hydrophilic": (None, -1.5),
    "intermediate": (-1.5, 1.5),
    "hydrophobic": (1.5, None),
}

PLDDT_THRESHOLDS = {
    "very_high": 90.0,
    "confident": 70.0,
    "low": 50.0,
}

COLORS = {
    "background": "#eef3f4",
    "surface": "#fbfcfb",
    "border": "#c4d1d5",
    "grid": "#d8e2e6",
    "text": "#132531",
    "muted_text": "#5a6d77",
    "accent": "#1f6a63",
    "accent_soft": "#deece9",
    "accent_border": "#93bab4",
    "warning": "#9d621f",
    "warning_bg": "#f7eedc",
    "warning_border": "#d4b98d",
    "contact_bg": "#f7eedc",
    "helix_fill": "#1f6a63",
    "helix_stroke": "#184f4b",
    "strand_fill": "#2e5fa0",
    "strand_stroke": "#234a7c",
    "turn_text": "#9a5a18",
    "coil": "#e5ecef",
    "msa_query_bg": "#214a5b",
    "msa_query_text": "#f3f8fa",
    "msa_identity_bg": "#1f6a63",
    "msa_identity_text": "#f3fbfa",
    "msa_similar_bg": "#dfece7",
    "msa_similar_text": "#214b48",
    "msa_default_bg": "#fbfcfb",
    "msa_default_text": "#1a2c38",
    "accessibility_buried": "#223746",
    "accessibility_intermediate": "#7496ab",
    "accessibility_accessible": "#bfd6df",
    "accessibility_highly_exposed": "#e1c18a",
    "hydropathy_hydrophilic": "#6d98ad",
    "hydropathy_intermediate": "#d6dde0",
    "hydropathy_hydrophobic": "#b87921",
    "contact_strong": "#9a5b1c",
    "contact_weak": "#2f68a6",
    "contact_multi_outline": "#2d7068",
    "disulfide_symbol": "#5FA79A",
    "plddt_very_high": "#174d86",
    "plddt_confident": "#5387aa",
    "plddt_low": "#c59133",
    "plddt_very_low": "#b56a39",
}

MSA_SIMILARITY_GROUPS = [
    frozenset("HKR"),
    frozenset("DE"),
    frozenset("STNQ"),
    frozenset("AVLIM"),
    frozenset("FYW"),
]

COMMON_IONS = {
    "NA",
    "K",
    "CA",
    "MG",
    "MN",
    "ZN",
    "FE",
    "CU",
    "CO",
    "CL",
    "CD",
    "NI",
}

COMMON_SUGARS = {
    "NAG",
    "BMA",
    "MAN",
    "FUC",
    "GAL",
    "GLC",
    "NDG",
    "SIA",
}

PORPHYRIN_LIKE = {
    "HEM",
    "HEC",
    "HEA",
    "HEB",
    "CLA",
    "BCL",
    "PQN",
}

NUCLEIC_ACID_RESNAMES = {
    "A",
    "C",
    "G",
    "U",
    "I",
    "DA",
    "DC",
    "DG",
    "DT",
    "DU",
}
