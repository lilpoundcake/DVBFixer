from pathlib import Path

from dvbfixer.diagnose.structural import check_duplicate_chain_coordinates
from dvbfixer.pdbutils.duplicates import duplicate_protein_chain_coordinates

FIXTURES = Path(__file__).parent / "fixtures"


def _atom(serial: int, atom: str, chain: str, x: float) -> str:
    element = atom[0]
    return (f"ATOM  {serial:5d} {atom:>4s} ALA {chain}   1    "
            f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00          {element:>2s}\n")


def test_coordinate_identical_protein_chains_are_reported(tmp_path: Path) -> None:
    pdb = tmp_path / "merged_frames.pdb"
    pdb.write_text(
        _atom(1, "N", "A", 0.0) + _atom(2, "CA", "A", 1.4)
        + _atom(3, "N", "B", 0.0) + _atom(4, "CA", "B", 1.4) + "END\n"
    )
    assert duplicate_protein_chain_coordinates(pdb) == [("A", "B", 2)]
    findings = check_duplicate_chain_coordinates(pdb)
    assert len(findings) == 1
    assert findings[0].severity.value == "WARNING"
    assert "MODEL/ENDMDL" in findings[0].message


def test_same_sequence_at_different_coordinates_is_not_reported(tmp_path: Path) -> None:
    pdb = tmp_path / "real_dimer.pdb"
    pdb.write_text(
        _atom(1, "N", "A", 0.0) + _atom(2, "CA", "A", 1.4)
        + _atom(3, "N", "B", 10.0) + _atom(4, "CA", "B", 11.4) + "END\n"
    )
    assert duplicate_protein_chain_coordinates(pdb) == []


def test_real_8dis_missing_model_separators_is_reported() -> None:
    pdb = FIXTURES / "overlap" / "8dis_t_u.pdb"
    assert duplicate_protein_chain_coordinates(pdb) == [("d", "D", 3695)]
