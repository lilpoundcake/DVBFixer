from pathlib import Path

from dvbfixer.pdbutils.inference import _apply_filter, _load_atom_table


def _atom(serial: int, name: str, residue: str, chain: str, resid: int, x: float) -> str:
    element = name[0]
    return (f"ATOM  {serial:5d} {name:>4s} {residue:>3s} {chain}{resid:4d}    "
            f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00          {element:>2s}\n")


def test_filter_drops_cn_contact_between_nonadjacent_residues(tmp_path: Path) -> None:
    pdb = tmp_path / "close-fold.pdb"
    pdb.write_text(
        _atom(1, "C", "ALA", "A", 1, 0.0)
        + _atom(2, "N", "GLY", "A", 2, 1.3)
        + _atom(3, "C", "GLY", "A", 2, 2.6)
        + _atom(4, "N", "ASN", "A", 3, 3.9)
        + _atom(5, "C", "ASN", "A", 3, 5.2)
        + "END\n"
    )
    atoms, by_serial = _load_atom_table(pdb)
    assert _apply_filter({(1, 2), (3, 4), (1, 4)}, atoms, by_serial) == {
        (1, 2), (3, 4)
    }
