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
