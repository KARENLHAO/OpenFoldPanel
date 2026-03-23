"""Reference-axis and residue-alignment helpers."""

from __future__ import annotations

from dataclasses import dataclass

from openfoldpanel.models import ChainRecord, ResidueRecord, SequenceAxisPosition


@dataclass(slots=True)
class ChainAlignment:
    chain: ChainRecord
    residue_by_axis_index: dict[int, ResidueRecord]
    warnings: list[str]


def build_sequence_axis(chain: ChainRecord) -> list[SequenceAxisPosition]:
    """Create a display axis directly from a reference chain."""

    axis: list[SequenceAxisPosition] = []
    for index, residue in enumerate(chain.residues):
        axis.append(
            SequenceAxisPosition(
                residue_index=index,
                chain=chain.chain_id,
                seq_id=residue.residue_id.seq_id,
                insertion_code=residue.residue_id.insertion_code,
                resname=residue.resname,
                one_letter=residue.one_letter,
                label=residue.residue_id.label,
            )
        )
    return axis


def align_chain_to_axis(axis: list[SequenceAxisPosition], chain: ChainRecord) -> ChainAlignment:
    """Map a chain onto the reference axis by residue keys or global alignment."""

    axis_lookup = {
        (position.seq_id, position.insertion_code, position.resname): index
        for index, position in enumerate(axis)
    }
    residue_by_axis_index: dict[int, ResidueRecord] = {}
    warnings: list[str] = []

    direct_match = True
    for residue in chain.residues:
        key = (residue.residue_id.seq_id, residue.residue_id.insertion_code, residue.resname)
        if key not in axis_lookup:
            direct_match = False
            break

    if direct_match and len(chain.residues) == len(axis):
        for residue in chain.residues:
            key = (residue.residue_id.seq_id, residue.residue_id.insertion_code, residue.resname)
            residue_by_axis_index[axis_lookup[key]] = residue
        return ChainAlignment(chain=chain, residue_by_axis_index=residue_by_axis_index, warnings=warnings)

    alignment = global_align("".join(pos.one_letter for pos in axis), chain.sequence)
    axis_index = -1
    chain_index = -1
    for ref_char, model_char in zip(alignment.reference, alignment.query):
        if ref_char != "-":
            axis_index += 1
        if model_char != "-":
            chain_index += 1
        if ref_char != "-" and model_char != "-" and axis_index >= 0 and chain_index >= 0:
            residue_by_axis_index[axis_index] = chain.residues[chain_index]
    warnings.append(
        f"Chain {chain.chain_id} was aligned conservatively to the reference axis because residue identifiers did not match directly."
    )
    return ChainAlignment(chain=chain, residue_by_axis_index=residue_by_axis_index, warnings=warnings)


@dataclass(slots=True)
class AlignmentResult:
    reference: str
    query: str


def global_align(reference: str, query: str) -> AlignmentResult:
    """A small Needleman-Wunsch implementation for residue-axis alignment."""

    match_score = 2
    mismatch_score = -1
    gap_score = -2

    rows = len(reference) + 1
    cols = len(query) + 1
    scores = [[0] * cols for _ in range(rows)]
    trace = [[""] * cols for _ in range(rows)]

    for i in range(1, rows):
        scores[i][0] = i * gap_score
        trace[i][0] = "U"
    for j in range(1, cols):
        scores[0][j] = j * gap_score
        trace[0][j] = "L"

    for i in range(1, rows):
        for j in range(1, cols):
            diag = scores[i - 1][j - 1] + (match_score if reference[i - 1] == query[j - 1] else mismatch_score)
            up = scores[i - 1][j] + gap_score
            left = scores[i][j - 1] + gap_score
            best = max(diag, up, left)
            scores[i][j] = best
            trace[i][j] = "D" if best == diag else ("U" if best == up else "L")

    aligned_ref: list[str] = []
    aligned_query: list[str] = []
    i = len(reference)
    j = len(query)
    while i > 0 or j > 0:
        move = trace[i][j] if i >= 0 and j >= 0 else ""
        if i > 0 and j > 0 and move == "D":
            aligned_ref.append(reference[i - 1])
            aligned_query.append(query[j - 1])
            i -= 1
            j -= 1
        elif i > 0 and (j == 0 or move == "U"):
            aligned_ref.append(reference[i - 1])
            aligned_query.append("-")
            i -= 1
        else:
            aligned_ref.append("-")
            aligned_query.append(query[j - 1])
            j -= 1
    return AlignmentResult(reference="".join(reversed(aligned_ref)), query="".join(reversed(aligned_query)))
