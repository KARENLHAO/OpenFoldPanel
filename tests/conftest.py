from __future__ import annotations

from pathlib import Path


def format_atom_line(
    serial: int,
    record: str,
    atom_name: str,
    resname: str,
    chain: str,
    resseq: int,
    x: float,
    y: float,
    z: float,
    *,
    element: str,
    bfactor: float = 90.0,
) -> str:
    return (
        f"{record:<6}{serial:>5} {atom_name:^4}{' ':1}{resname:>3} {chain:1}"
        f"{resseq:>4}{' ':1}   "
        f"{x:>8.3f}{y:>8.3f}{z:>8.3f}"
        f"{1.00:>6.2f}{bfactor:>6.2f}          "
        f"{element:>2}"
    )


def build_test_pdb() -> str:
    lines = [
        format_atom_line(1, "ATOM", "N", "ALA", "A", 1, 0.0, 0.0, 0.0, element="N"),
        format_atom_line(2, "ATOM", "CA", "ALA", "A", 1, 1.5, 0.0, 0.0, element="C"),
        format_atom_line(3, "ATOM", "C", "ALA", "A", 1, 2.5, 0.0, 0.0, element="C"),
        format_atom_line(4, "ATOM", "O", "ALA", "A", 1, 3.5, 0.0, 0.0, element="O"),
        format_atom_line(5, "ATOM", "N", "GLY", "A", 2, 4.5, 0.0, 0.0, element="N"),
        format_atom_line(6, "ATOM", "CA", "GLY", "A", 2, 5.5, 0.0, 0.0, element="C"),
        format_atom_line(7, "ATOM", "C", "GLY", "A", 2, 6.5, 0.0, 0.0, element="C"),
        format_atom_line(8, "ATOM", "O", "GLY", "A", 2, 7.5, 0.0, 0.0, element="O"),
        format_atom_line(9, "ATOM", "N", "SER", "A", 3, 8.5, 0.0, 0.0, element="N"),
        format_atom_line(10, "ATOM", "CA", "SER", "A", 3, 9.5, 0.0, 0.0, element="C"),
        format_atom_line(11, "ATOM", "C", "SER", "A", 3, 10.5, 0.0, 0.0, element="C"),
        format_atom_line(12, "ATOM", "O", "SER", "A", 3, 11.5, 0.0, 0.0, element="O"),
        format_atom_line(13, "ATOM", "CA", "TYR", "B", 1, 1.8, 2.8, 0.0, element="C"),
        format_atom_line(14, "ATOM", "CB", "TYR", "B", 1, 1.8, 3.4, 0.0, element="C"),
        format_atom_line(15, "HETATM", "ZN", "ZN", "Z", 1, 1.0, 2.0, 0.0, element="ZN"),
    ]
    return "\n".join(lines) + "\nEND\n"


def write_test_pdb(path: Path) -> Path:
    path.write_text(build_test_pdb(), encoding="utf-8")
    return path


def build_test_mmcif(*, chain_id: str = "A", seq_id_start: int = 1) -> str:
    atoms = [
        ("ATOM", 1, "N", "N", "ALA", chain_id, seq_id_start, 0.0, 0.0, 0.0, 90.0),
        ("ATOM", 2, "C", "CA", "ALA", chain_id, seq_id_start, 1.5, 0.0, 0.0, 90.0),
        ("ATOM", 3, "C", "C", "ALA", chain_id, seq_id_start, 2.5, 0.0, 0.0, 90.0),
        ("ATOM", 4, "O", "O", "ALA", chain_id, seq_id_start, 3.5, 0.0, 0.0, 90.0),
        ("ATOM", 5, "N", "N", "GLY", chain_id, seq_id_start + 1, 4.5, 0.0, 0.0, 90.0),
        ("ATOM", 6, "C", "CA", "GLY", chain_id, seq_id_start + 1, 5.5, 0.0, 0.0, 90.0),
        ("ATOM", 7, "C", "C", "GLY", chain_id, seq_id_start + 1, 6.5, 0.0, 0.0, 90.0),
        ("ATOM", 8, "O", "O", "GLY", chain_id, seq_id_start + 1, 7.5, 0.0, 0.0, 90.0),
        ("ATOM", 9, "N", "N", "SER", chain_id, seq_id_start + 2, 8.5, 0.0, 0.0, 90.0),
        ("ATOM", 10, "C", "CA", "SER", chain_id, seq_id_start + 2, 9.5, 0.0, 0.0, 90.0),
        ("ATOM", 11, "C", "C", "SER", chain_id, seq_id_start + 2, 10.5, 0.0, 0.0, 90.0),
        ("ATOM", 12, "O", "O", "SER", chain_id, seq_id_start + 2, 11.5, 0.0, 0.0, 90.0),
    ]
    lines = [
        "data_demo",
        "#",
        "loop_",
        "_atom_site.group_PDB",
        "_atom_site.id",
        "_atom_site.type_symbol",
        "_atom_site.label_atom_id",
        "_atom_site.label_alt_id",
        "_atom_site.label_comp_id",
        "_atom_site.label_asym_id",
        "_atom_site.label_entity_id",
        "_atom_site.label_seq_id",
        "_atom_site.pdbx_PDB_ins_code",
        "_atom_site.Cartn_x",
        "_atom_site.Cartn_y",
        "_atom_site.Cartn_z",
        "_atom_site.occupancy",
        "_atom_site.B_iso_or_equiv",
        "_atom_site.auth_seq_id",
        "_atom_site.auth_comp_id",
        "_atom_site.auth_asym_id",
        "_atom_site.auth_atom_id",
        "_atom_site.pdbx_PDB_model_num",
    ]
    for record, serial, element, atom_name, resname, current_chain, seq_id, x, y, z, bfactor in atoms:
        label_seq = seq_id - seq_id_start + 1
        lines.append(
            f"{record} {serial} {element} {atom_name} . {resname} {current_chain} 1 {label_seq} ? "
            f"{x:.3f} {y:.3f} {z:.3f} 1.00 {bfactor:.2f} {seq_id} {resname} {current_chain} {atom_name} 1"
        )
    lines.append("#")
    return "\n".join(lines) + "\n"


def write_test_mmcif(path: Path, *, chain_id: str = "A", seq_id_start: int = 1) -> Path:
    path.write_text(build_test_mmcif(chain_id=chain_id, seq_id_start=seq_id_start), encoding="utf-8")
    return path
