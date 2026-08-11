import pytest

from dvbfixer import msa


def test_fasta_validation_rejects_duplicate_ids(tmp_path):
    fasta = tmp_path / "input.fasta"
    fasta.write_text(">same\nAAAA\n>same\nAAAT\n")
    with pytest.raises(ValueError, match="unique"):
        msa._read_fasta(fasta)


def test_mafft_adapter_validates_and_writes_alignment(tmp_path, monkeypatch):
    fasta = tmp_path / "input.fasta"
    fasta.write_text(">target\nACDE\n>template\nACE\n")
    executable = tmp_path / "mafft"
    executable.write_text("#!/bin/sh\nprintf '>target\\nACDE\\n>template\\nAC-E\\n'\n")
    executable.chmod(0o755)
    monkeypatch.setattr(msa, "available_engines", lambda: {
        "mafft": str(executable), "muscle": None, "clustalo": None,
    })
    output = tmp_path / "aligned.fasta"
    assert msa.run_alignment(fasta, output) == "mafft"
    assert output.read_text() == ">target\nACDE\n>template\nAC-E\n"


def test_pir_output(tmp_path, monkeypatch):
    fasta = tmp_path / "input.fasta"
    fasta.write_text(">a\nAAA\n>b\nABA\n")
    executable = tmp_path / "mafft"
    executable.write_text("#!/bin/sh\nprintf '>a\\nAAA\\n>b\\nABA\\n'\n")
    executable.chmod(0o755)
    monkeypatch.setattr(msa, "available_engines", lambda: {
        "mafft": str(executable), "muscle": None, "clustalo": None,
    })
    output = tmp_path / "aligned.pir"
    msa.run_alignment(fasta, output, output_format="pir")
    assert ">P1;a" in output.read_text()
    assert "AAA*" in output.read_text()


def test_missing_selected_engine_is_actionable(tmp_path, monkeypatch):
    fasta = tmp_path / "input.fasta"
    fasta.write_text(">a\nAAA\n>b\nABA\n")
    monkeypatch.setattr(msa, "available_engines", lambda: {
        "mafft": None, "muscle": None, "clustalo": None,
    })
    with pytest.raises(RuntimeError, match="not available"):
        msa.run_alignment(fasta, tmp_path / "out.fasta", engine="muscle")


def test_template_sequence_can_supply_second_record(tmp_path, monkeypatch):
    fasta = tmp_path / "input.fasta"
    fasta.write_text(">target\nAAA\n")
    pdb = tmp_path / "template.pdb"
    pdb.write_text(
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C  \n"
        "ATOM      2  CA  ALA A   2       1.000   0.000   0.000  1.00 20.00           C  \n"
        "ATOM      3  CA  ALA A   3       2.000   0.000   0.000  1.00 20.00           C  \n"
    )
    executable = tmp_path / "mafft"
    executable.write_text("#!/bin/sh\ncat \"$2\"\n")
    executable.chmod(0o755)
    monkeypatch.setattr(msa, "available_engines", lambda: {
        "mafft": str(executable), "muscle": None, "clustalo": None,
    })
    output = tmp_path / "aligned.fasta"
    msa.run_alignment(fasta, output, templates=[f"{pdb}:A"])
    assert ">template_1_template_A" in output.read_text()
