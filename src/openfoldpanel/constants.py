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
DEFAULT_MAX_HITS = 10
DEFAULT_MSA_DISPLAY_ROWS = 1
DEFAULT_CONTACT_CUTOFF = 3.7
DEFAULT_STRONG_CONTACT_CUTOFF = 3.2

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
    "background": "#f6f4ee",
    "surface": "#fffdf8",
    "border": "#cad5dc",
    "grid": "#dbe3e7",
    "text": "#152433",
    "muted_text": "#5b6d79",
    "accent": "#0f766e",
    "accent_soft": "#d7efea",
    "accent_border": "#96c9be",
    "warning": "#8a5a12",
    "warning_bg": "#fff4de",
    "warning_border": "#e0bc7f",
    "helix_fill": "#0f766e",
    "helix_stroke": "#115e59",
    "strand_fill": "#2563eb",
    "strand_stroke": "#1d4ed8",
    "turn_text": "#b45309",
    "coil": "#edf2f5",
    "msa_identity_bg": "#0f766e",
    "msa_identity_text": "#ffffff",
    "msa_similar_bg": "#d7efea",
    "msa_similar_text": "#0f3f45",
    "msa_default_bg": "#fffdf8",
    "msa_default_text": "#152433",
    "accessibility_buried": "#243746",
    "accessibility_intermediate": "#7fb6d6",
    "accessibility_accessible": "#b7d8e9",
    "accessibility_highly_exposed": "#f2d39b",
    "hydropathy_hydrophilic": "#6baed6",
    "hydropathy_intermediate": "#d8dde0",
    "hydropathy_hydrophobic": "#d97706",
    "contact_strong": "#b45309",
    "contact_weak": "#1d4ed8",
    "contact_multi_outline": "#0f766e",
    "plddt_very_high": "#0b4da2",
    "plddt_confident": "#4f9dd5",
    "plddt_low": "#d9a441",
    "plddt_very_low": "#c96b32",
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
