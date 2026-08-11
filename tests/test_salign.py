import pytest

from dvbfixer.salign import extract_chain, parse_template_spec, run_biopython_superposition


def _atom(serial: int, chain: str, residue: int) -> str:
    return (
        f"ATOM  {serial:5d}  CA  ALA {chain}{residue:4d}    "
        "   0.000   0.000   0.000  1.00 20.00           C  \n"
    )


def test_template_spec_accepts_path_and_chain(tmp_path):
    pdb = tmp_path / "template.pdb"
    pdb.write_text(_atom(1, "A", 1))
    assert parse_template_spec(str(pdb)) == (pdb.resolve(), None)
    assert parse_template_spec(f"{pdb}:A") == (pdb.resolve(), "A")


def test_extract_chain_keeps_only_requested_atoms(tmp_path):
    source = tmp_path / "source.pdb"
    source.write_text("HEADER test\n" + _atom(1, "A", 1) + _atom(2, "B", 1) + "END\n")
    output = tmp_path / "chain.pdb"
    extract_chain(source, "B", output)
    text = output.read_text()
    assert "HEADER test" in text
    assert "ALA B" in text
    assert "ALA A" not in text


def test_extract_chain_rejects_unknown_chain(tmp_path):
    source = tmp_path / "source.pdb"
    source.write_text(_atom(1, "A", 1))
    with pytest.raises(ValueError, match="no atoms"):
        extract_chain(source, "Z", tmp_path / "missing.pdb")


def test_biopython_superposition_writes_fitted_structures(tmp_path, monkeypatch):
    reference = tmp_path / "reference.pdb"
    mobile = tmp_path / "mobile.pdb"
    reference.write_text("".join(_atom(index, "A", index) for index in range(1, 4)) + "END\n")
    mobile.write_text("".join(_atom(index, "B", index) for index in range(1, 4)) + "END\n")

    def fake_alignment(input_path, output_path, engine, output_format, verbose):
        output_path.write_text(input_path.read_text())
        return "test"

    monkeypatch.setattr("dvbfixer.msa.run_alignment", fake_alignment)
    output = tmp_path / "alignment.pir"
    fit_dir = tmp_path / "fitted"
    fitted = run_biopython_superposition(
        [f"{reference}:A", f"{mobile}:B"], output, fit_dir,
    )

    assert len(fitted) == 2
    assert all(path.exists() for path in fitted)
    assert output.read_text().count(">P1;") == 2
