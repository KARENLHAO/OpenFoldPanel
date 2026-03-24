"""Typed data models used across the pipeline."""

from __future__ import annotations

from dataclasses import dataclass, field, fields as dataclass_fields, is_dataclass
from pathlib import Path
from typing import Any


@dataclass(slots=True)
class PipelineConfig:
    input_path: Path
    outdir: Path
    chain: str
    columns: int
    max_homologs_displayed: int
    evalue: str
    font_size: int
    hyd_window: int
    msa_db: Path | None
    disable_msa: bool
    keep_temp: bool
    contact_cutoff: float
    strong_contact_cutoff: float
    verbose: bool = False


@dataclass(slots=True)
class JobDefinition:
    name: str
    root_dir: Path
    structure_files: list[Path]
    ignored_files: list[str] = field(default_factory=list)


@dataclass(slots=True)
class ResidueId:
    chain_id: str
    seq_id: int
    insertion_code: str = ""

    @property
    def label(self) -> str:
        return f"{self.seq_id}{self.insertion_code}".strip()


@dataclass(slots=True)
class AtomRecord:
    atom_name: str
    element: str
    x: float
    y: float
    z: float
    bfactor: float | None = None
    occupancy: float | None = None
    is_hetatm: bool = False


@dataclass(slots=True)
class ResidueRecord:
    residue_id: ResidueId
    resname: str
    atoms: list[AtomRecord]
    residue_type: str
    one_letter: str
    auth_seq_id: str


@dataclass(slots=True)
class ChainRecord:
    chain_id: str
    residues: list[ResidueRecord]
    entity_type: str

    @property
    def sequence(self) -> str:
        return "".join(res.one_letter for res in self.residues)


@dataclass(slots=True)
class ParsedStructure:
    name: str
    source_path: Path
    chains: dict[str, ChainRecord]
    format: str


@dataclass(slots=True)
class SequenceAxisPosition:
    residue_index: int
    chain: str
    seq_id: int
    insertion_code: str
    resname: str
    one_letter: str
    label: str


@dataclass(slots=True)
class SecondaryStructureEntry:
    residue_index: int
    dssp_code: str | None
    category: str


@dataclass(slots=True)
class AccessibilityEntry:
    residue_index: int
    absolute: float | None
    relative: float | None
    category: str | None


@dataclass(slots=True)
class HydropathyEntry:
    residue_index: int
    residue: str
    value: float | None
    category: str | None


@dataclass(slots=True)
class ContactHit:
    partner_type: str
    partner_chain: str | None
    partner_resname: str
    partner_resid: str
    min_distance: float
    symbol: str
    strength_category: str


@dataclass(slots=True)
class ContactEntry:
    residue_index: int
    partner_type: str | None
    partner_chain: str | None
    partner_resname: str | None
    partner_resid: str | None
    min_distance: float | None
    symbol: str | None
    strength_category: str | None
    is_multi_contact: bool = False
    all_contacts: list[ContactHit] = field(default_factory=list)


@dataclass(slots=True)
class DisulfideBond:
    residue_index_a: int
    residue_index_b: int
    chain_a: str
    chain_b: str


@dataclass(slots=True)
class ModelTracks:
    name: str
    source_path: str
    chain: str
    secondary_structure: list[SecondaryStructureEntry]
    plddt: list[float | None]
    accessibility: list[AccessibilityEntry]
    contacts: list[ContactEntry]
    disulfides: list[DisulfideBond] = field(default_factory=list)
    display_name: str | None = None


@dataclass(slots=True)
class MSARow:
    identifier: str
    sequence: str
    is_query: bool = False


@dataclass(slots=True)
class ConservationEntry:
    residue_index: int
    identity_fraction: float
    similarity_fraction: float
    style: str


@dataclass(slots=True)
class MSAData:
    enabled: bool
    query: str
    rows: list[MSARow] = field(default_factory=list)
    conservation: list[ConservationEntry] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    leading_display_overrides: list[str | None] = field(
        default_factory=list,
        repr=False,
        compare=False,
        metadata={"serialize": False},
    )


@dataclass(slots=True)
class RenderConfig:
    columns: int
    max_homologs_displayed: int
    font_size: int
    cell_width: float
    row_height: float
    label_width: float
    margin: float
    colors: dict[str, str]
    font_family: str
    heading_font_family: str


@dataclass(slots=True)
class JobPanelData:
    job_name: str
    reference_chain: str
    sequence_axis: list[SequenceAxisPosition]
    models: list[ModelTracks]
    msa: MSAData
    hydropathy: list[HydropathyEntry]
    render_config: RenderConfig
    warnings: list[str] = field(default_factory=list)
    status: str = "success"


@dataclass(slots=True)
class JobReportData:
    job_name: str
    default_reference_chain: str
    chain_panels: list[JobPanelData]
    warnings: list[str] = field(default_factory=list)
    status: str = "success"


@dataclass(slots=True)
class JobRunResult:
    job_name: str
    status: str
    output_dir: str
    warnings: list[str] = field(default_factory=list)
    error: str | None = None
    artifacts: list[str] = field(default_factory=list)


def dataclass_to_dict(value: Any) -> Any:
    """Convert nested dataclasses into plain serializable dictionaries."""

    if is_dataclass(value):
        return {
            field_info.name: dataclass_to_dict(getattr(value, field_info.name))
            for field_info in dataclass_fields(value)
            if field_info.metadata.get("serialize", True)
        }
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, list):
        return [dataclass_to_dict(item) for item in value]
    if isinstance(value, dict):
        return {key: dataclass_to_dict(item) for key, item in value.items()}
    return value
