import json

from dvbfixer.ffutils.artificial_terminals import (
    normalize_fasta_truncated_terminal_variants,
)


def _atom(serial, atom, resname, chain, resid):
    element = atom[0]
    return (
        f"ATOM  {serial:5d} {atom:^4s} {resname:>3s} {chain}{resid:4d}    "
        f"   0.000   0.000   0.000  1.00  0.00          {element:>2s}\n"
    )


def test_normalizes_only_fasta_truncated_c_terminal_variant(tmp_path):
    pdb = tmp_path / "prepared.pdb"
    pdb.write_text(
        _atom(1, "N", "ALA", "B", 10)
        + _atom(2, "CA", "ALA", "B", 10)
        + _atom(3, "N", "GLH", "B", 11)
        + _atom(4, "CA", "GLH", "B", 11)
        + _atom(5, "HE2", "GLH", "B", 11)
        + _atom(6, "OXT", "GLH", "B", 11)
        + "TER\nEND\n"
    )
    fasta = tmp_path / "full.fasta"
    fasta.write_text(">chain_B\nAEG\n")
    dat = tmp_path / "prepared.dat"
    dat.write_text(json.dumps({
        "description": "test", "total_added": 2,
        "residue_summary": {"B/GLH11": {"heavy": 1, "hydrogen": 1}},
        "added_atoms": [
            {"chain": "B", "resid": "11", "icode": "", "resname": "GLH", "atom": "HE2", "element": "H"},
            {"chain": "B", "resid": "11", "icode": "", "resname": "GLH", "atom": "OXT", "element": "O"},
        ],
        "variant_overrides": {"B:11:": "GLH"},
    }))

    changed = normalize_fasta_truncated_terminal_variants(pdb, fasta)

    assert changed == [("B", "11", "GLH", "GLU")]
    text = pdb.read_text()
    assert " GLU B  11" in text
    assert "GLH" not in text
    assert "HE2" not in text
    data = json.loads(dat.read_text())
    assert "variant_overrides" not in data
    assert [atom["atom"] for atom in data["added_atoms"]] == ["OXT"]
    assert data["added_atoms"][0]["resname"] == "GLU"


def test_keeps_variant_at_true_fasta_terminus(tmp_path):
    pdb = tmp_path / "prepared.pdb"
    pdb.write_text(
        _atom(1, "N", "ALA", "A", 1)
        + _atom(2, "CA", "ALA", "A", 1)
        + _atom(3, "N", "GLH", "A", 2)
        + _atom(4, "CA", "GLH", "A", 2)
        + _atom(5, "HE2", "GLH", "A", 2)
        + "END\n"
    )
    fasta = tmp_path / "full.fasta"
    fasta.write_text(">chain_A\nAE\n")

    assert normalize_fasta_truncated_terminal_variants(pdb, fasta) == []
    assert "GLH" in pdb.read_text()


def test_normalizes_uncapped_physical_terminal_without_fasta(tmp_path):
    pdb = tmp_path / "prepared.pdb"
    pdb.write_text(_atom(1, "N", "ALA", "A", 1) + _atom(2, "N", "GLH", "A", 2))

    assert normalize_fasta_truncated_terminal_variants(pdb) == [
        ("A", "2", "GLH", "GLU")
    ]


def test_preserves_variant_when_capped(tmp_path):
    pdb = tmp_path / "prepared.pdb"
    pdb.write_text(
        _atom(1, "N", "ACE", "A", 0)
        + _atom(2, "N", "ASH", "A", 1)
        + _atom(3, "N", "GLH", "A", 2)
        + _atom(4, "N", "NME", "A", 3)
    )

    assert normalize_fasta_truncated_terminal_variants(pdb) == []
    assert "ASH" in pdb.read_text()
    assert "GLH" in pdb.read_text()


def test_charmm_terminal_protonation_is_not_changed(tmp_path):
    pdb = tmp_path / "prepared.pdb"
    pdb.write_text(_atom(1, "N", "ALA", "A", 1) + _atom(2, "N", "GLH", "A", 2))

    assert normalize_fasta_truncated_terminal_variants(
        pdb, force_field="charmm"
    ) == []
    assert "GLH" in pdb.read_text()
