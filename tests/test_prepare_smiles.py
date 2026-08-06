"""Optional SMILES-guided small-molecule preparation."""

from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("rdkit", reason="SMILES guidance requires RDKit")
pytest.importorskip("openmm", reason="SMILES guidance requires OpenMM")

from openmm import Vec3
from openmm.app import Topology, element
from openmm.unit import Quantity, nanometer

from dvbfixer.prepare.cli import parse_args
from dvbfixer.prepare.pipeline import main as prepare_main
from dvbfixer.prepare.smiles import (
    SmilesPreparationError,
    add_hydrogens_from_smiles,
    parse_smiles_mappings,
)
from dvbfixer.zbs import _run_pipeline
from dvbfixer.zbs import parse_args as parse_zbs_args


def _topology(residues, bonds):
    top = Topology()
    chain = top.addChain("A")
    atoms = {}
    xyz = []
    for resname, resid, atom_specs in residues:
        residue = top.addResidue(resname, chain, str(resid))
        for name, elem, position in atom_specs:
            atoms[(resid, name)] = top.addAtom(name, elem, residue)
            xyz.append(Vec3(*position))
    for left, right in bonds:
        top.addBond(atoms[left], atoms[right])
    return top, Quantity(xyz, nanometer)


def _atoms_by_residue(topology):
    return {
        (res.name, res.id): [(atom.name, atom.element.symbol) for atom in res.atoms()]
        for res in topology.residues()
    }


def test_smiles_is_optional_and_repeatable_in_prepare_and_zbs():
    assert parse_args(["input.pdb"]).smiles == []
    assert parse_zbs_args(["input.pdb"]).smiles == []
    values = ["LIG=C=C", "DRG=[NH4+]"]
    assert parse_args(["input.pdb", "--smiles", values[0], "--smiles", values[1]]).smiles == values
    assert parse_zbs_args(["input.pdb", "--smiles", values[0]]).smiles == values[:1]
    assert parse_smiles_mappings(values) == {"LIG": "C=C", "DRG": "[NH4+]"}


@pytest.mark.parametrize("value, message", [
    ("CC", "RESNAME=SMILES"),
    ("LIG=not-a-smiles", "invalid SMILES"),
    ("LIG=CC.O", "one connected molecule"),
])
def test_invalid_smiles_mapping_fails(value, message):
    with pytest.raises(SmilesPreparationError, match=message):
        parse_smiles_mappings([value])


def test_duplicate_smiles_mapping_fails():
    with pytest.raises(SmilesPreparationError, match="duplicate"):
        parse_smiles_mappings(["LIG=CC", "LIG=C=C"])


def test_smiles_adds_alkene_h_and_preserves_heavy_coordinates():
    top, positions = _topology(
        [("LIG", 1, [
            ("C1", element.carbon, (0.0, 0.0, 0.0)),
            ("C2", element.carbon, (0.134, 0.0, 0.0)),
        ])],
        [((1, "C1"), (1, "C2"))],
    )
    new_top, new_positions = add_hydrogens_from_smiles(top, positions, {"LIG": "C=C"})
    atoms = list(new_top.atoms())
    assert [(a.name, a.element.symbol) for a in atoms[:3]] == [
        ("C1", "C"), ("H1", "H"), ("H12", "H"),
    ]
    assert sum(a.element.symbol == "H" for a in atoms) == 4
    heavy_positions = [new_positions[a.index] for a in atoms if a.element.symbol != "H"]
    assert tuple(heavy_positions[0].value_in_unit(nanometer)) == pytest.approx((0.0, 0.0, 0.0))
    assert tuple(heavy_positions[1].value_in_unit(nanometer)) == pytest.approx((0.134, 0.0, 0.0))
    cc_bond = next(b for b in new_top.bonds()
                   if b[0].element.symbol == "C" and b[1].element.symbol == "C")
    assert cc_bond.order == 2.0


def test_smiles_formal_charge_controls_hydrogen_count():
    top, positions = _topology(
        [("AMM", 1, [("N1", element.nitrogen, (0, 0, 0))])], [],
    )
    new_top, _ = add_hydrogens_from_smiles(top, positions, {"AMM": "[NH4+]"})
    assert sum(a.element.symbol == "H" for a in new_top.atoms()) == 4


def test_smiles_replaces_existing_h_but_leaves_unmapped_ligand_unchanged():
    top, positions = _topology(
        [
            ("LIG", 1, [
                ("C1", element.carbon, (0.0, 0.0, 0.0)),
                ("C2", element.carbon, (0.134, 0.0, 0.0)),
                ("OLD", element.hydrogen, (-0.1, 0.0, 0.0)),
            ]),
            ("UNK", 2, [
                ("C1", element.carbon, (0.5, 0.0, 0.0)),
                ("KEEP", element.hydrogen, (0.6, 0.0, 0.0)),
            ]),
        ],
        [
            ((1, "C1"), (1, "C2")),
            ((1, "C1"), (1, "OLD")),
            ((2, "C1"), (2, "KEEP")),
        ],
    )
    new_top, _ = add_hydrogens_from_smiles(top, positions, {"LIG": "C=C"})
    residues = _atoms_by_residue(new_top)
    assert "OLD" not in {name for name, _ in residues[("LIG", "1")]}
    assert ("KEEP", "H") in residues[("UNK", "2")]


def test_mapping_applies_to_every_matching_instance():
    top, positions = _topology(
        [
            ("LIG", 1, [("C1", element.carbon, (0, 0, 0)),
                         ("C2", element.carbon, (0.134, 0, 0))]),
            ("LIG", 2, [("C1", element.carbon, (0.5, 0, 0)),
                         ("C2", element.carbon, (0.634, 0, 0))]),
        ],
        [((1, "C1"), (1, "C2")), ((2, "C1"), (2, "C2"))],
    )
    new_top, _ = add_hydrogens_from_smiles(top, positions, {"LIG": "C=C"})
    for residue in new_top.residues():
        assert sum(a.element.symbol == "H" for a in residue.atoms()) == 4


def test_symmetric_aromatic_mapping_is_accepted():
    specs = []
    for i in range(6):
        specs.append((f"C{i + 1}", element.carbon, (i * 0.14, 0, 0)))
    bonds = [((1, f"C{i + 1}"), (1, f"C{(i + 1) % 6 + 1}")) for i in range(6)]
    top, positions = _topology([("BEN", 1, specs)], bonds)
    new_top, _ = add_hydrogens_from_smiles(top, positions, {"BEN": "c1ccccc1"})
    assert sum(a.element.symbol == "H" for a in new_top.atoms()) == 6
    assert all(b.order == 1.5 for b in new_top.bonds()
               if b[0].element.symbol != "H" and b[1].element.symbol != "H")


def test_chemically_distinct_symmetric_mapping_fails():
    top, positions = _topology(
        [("ACT", 1, [
            ("C1", element.carbon, (0, 0, 0)),
            ("C2", element.carbon, (0.15, 0, 0)),
            ("O1", element.oxygen, (0.27, 0.08, 0)),
            ("O2", element.oxygen, (0.27, -0.08, 0)),
        ])],
        [
            ((1, "C1"), (1, "C2")), ((1, "C2"), (1, "O1")),
            ((1, "C2"), (1, "O2")),
        ],
    )
    with pytest.raises(SmilesPreparationError, match="chemically distinct"):
        add_hydrogens_from_smiles(top, positions, {"ACT": "CC(=O)O"})


def test_ionized_resonance_oxygens_map_without_false_ambiguity():
    top, positions = _topology(
        [("ACT", 1, [
            ("C1", element.carbon, (0, 0, 0)),
            ("C2", element.carbon, (0.15, 0, 0)),
            ("O1", element.oxygen, (0.27, 0.08, 0)),
            ("O2", element.oxygen, (0.27, -0.08, 0)),
        ])],
        [
            ((1, "C1"), (1, "C2")), ((1, "C2"), (1, "O1")),
            ((1, "C2"), (1, "O2")),
        ],
    )
    new_top, _ = add_hydrogens_from_smiles(top, positions, {"ACT": "CC(=O)[O-]"})
    oxygen_names = {a.name for a in new_top.atoms() if a.element.symbol == "O"}
    oxygen_h_parents = {
        b[0].name for b in new_top.bonds()
        if b[0].element.symbol == "O" and b[1].element.symbol == "H"
    } | {
        b[1].name for b in new_top.bonds()
        if b[1].element.symbol == "O" and b[0].element.symbol == "H"
    }
    assert oxygen_names == {"O1", "O2"}
    assert oxygen_h_parents == set()


def test_missing_mismatched_and_covalent_targets_fail_clearly():
    top, positions = _topology(
        [("LIG", 1, [("C1", element.carbon, (0, 0, 0)),
                      ("C2", element.carbon, (0.14, 0, 0))])],
        [((1, "C1"), (1, "C2"))],
    )
    with pytest.raises(SmilesPreparationError, match="no ligand residue"):
        add_hydrogens_from_smiles(top, positions, {"DRG": "CC"})
    with pytest.raises(SmilesPreparationError, match="does not match"):
        add_hydrogens_from_smiles(top, positions, {"LIG": "CO"})

    covalent_top, covalent_positions = _topology(
        [
            ("LIG", 1, [("C1", element.carbon, (0, 0, 0)),
                         ("C2", element.carbon, (0.14, 0, 0))]),
            ("ALA", 2, [("CA", element.carbon, (0.28, 0, 0))]),
        ],
        [((1, "C1"), (1, "C2")), ((1, "C2"), (2, "CA"))],
    )
    with pytest.raises(SmilesPreparationError, match="external heavy-atom bond"):
        add_hydrogens_from_smiles(covalent_top, covalent_positions, {"LIG": "CC"})


@pytest.mark.parametrize("flag", [
    "--strip-heterogens", "--no-heterogen-h",
])
def test_prepare_rejects_incompatible_smiles_flags(tmp_path, flag):
    source = tmp_path / "input.pdb"
    source.write_text("END\n")
    with pytest.raises(SystemExit) as exc:
        prepare_main([str(source), "--smiles", "LIG=CC", flag])
    assert exc.value.code == 2


def test_zbs_forwards_smiles_to_prepare(monkeypatch, tmp_path):
    source = tmp_path / "input.pdb"
    source.write_text("END\n")
    output = tmp_path / "output.pdb"
    captured = []

    def fake_prepare(argv):
        captured.extend(argv)
        destination = Path(argv[argv.index("-o") + 1])
        destination.write_text("END\n")

    monkeypatch.setattr("dvbfixer.prepare.main", fake_prepare)
    args = parse_zbs_args([
        str(source), "-o", str(output), "--skip-renumber", "--skip-model",
        "--skip-minimize", "--no-postflight", "--no-align-to-input",
        "--smiles", "LIG=C=C",
    ])
    _run_pipeline(args, source)
    index = captured.index("--smiles")
    assert captured[index + 1] == "LIG=C=C"
    assert output.exists()
