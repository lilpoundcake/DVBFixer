"""Pure B-factor policy: no OpenMM required in this test module."""

from pathlib import Path

import pytest

from dvbfixer.model.diffusion.anchor_policy import (
    AnchorMovementLimits,
    AnchorPositionRestraint,
    build_anchor_restraints,
    check_anchor_movements,
)
from dvbfixer.model.diffusion.contract import AtomIdentity, ResidueIdentity

SOURCE = Path(__file__).parent / "fixtures/protein_context/9cqd_split.pdb"


def observed(text):
    return tuple(AtomIdentity(l[21], l[22:26].strip(), l[26].strip(), l[12:16].strip())
                 for l in text.splitlines() if l.startswith("ATOM  "))


def test_chain_reference_is_independent_of_anchor_mask_and_uses_source_only():
    text = SOURCE.read_text()
    fixed = observed(text)
    left, right = ResidueIdentity("W", "25"), ResidueIdentity("X", "139")
    single = build_anchor_restraints(text, (left,), observed_atoms=fixed)
    both = build_anchor_restraints(text, (left, right), observed_atoms=fixed)
    assert len(single) == 6 and len(both) == 13
    assert single == both[:6]
    assert single[0].backbone_b == pytest.approx(88.725)
    assert single[0].chain_reference_b == pytest.approx(74.7925)
    assert single[0].k_kj_mol_nm2 == pytest.approx(10000 * 74.7925 / 88.725)
    assert single[0].k_kj_mol_nm2 != 10000  # single-anchor weighting must not cancel
    assert all(x.fallback_reason is None for x in both)


@pytest.mark.parametrize("header,reason", [
    ("", "provenance"), ("EXPDTA    ELECTRON MICROSCOPY\n", "provenance"),
    ("EXPDTA    X-RAY DIFFRACTION\nREMARK pLDDT stored in B\n", "provenance"),
])
def test_unknown_and_confidence_provenance_falls_back(header, reason):
    text = header + "".join(l for l in SOURCE.read_text().splitlines(keepends=True) if l.startswith("ATOM  "))
    result = build_anchor_restraints(text, (ResidueIdentity("W", "25"),), observed_atoms=observed(text))
    assert all(x.k_kj_mol_nm2 == 10000 and reason in x.fallback_reason for x in result)


def test_explicit_disable_and_override():
    text = SOURCE.read_text()
    args = {"anchors": (ResidueIdentity("W", "25"),), "observed_atoms": observed(text)}
    disabled = build_anchor_restraints(text, **args, bfactor_source="disabled")
    assert all(x.k_kj_mol_nm2 == 10000 and x.fallback_reason == "B weighting disabled" for x in disabled)
    no_header = "".join(l for l in text.splitlines(keepends=True) if l.startswith("ATOM  "))
    override = build_anchor_restraints(no_header, **args, bfactor_source="experimental")
    assert override[0].k_kj_mol_nm2 < 10000 and override[0].fallback_reason is None


@pytest.mark.parametrize("field", ["      ", "   nan", " -1.00", "  0.00"])
def test_invalid_anchor_b_values_fall_back_without_invalid_metadata(field):
    text = "".join(l[:60] + field + l[66:] if l.startswith("ATOM  ") and l[21] == "W" and
                   l[22:26].strip() == "25" and l[12:16].strip() == "CA" else l
                   for l in SOURCE.read_text().splitlines(keepends=True))
    result = build_anchor_restraints(text, (ResidueIdentity("W", "25"),), observed_atoms=observed(text))
    assert all(x.k_kj_mol_nm2 == 10000 and x.backbone_b is None and "invalid" in x.fallback_reason for x in result)


def test_incomplete_backbone_and_small_population_fall_back():
    text = SOURCE.read_text()
    fixed = observed(text)
    result = build_anchor_restraints(text, (ResidueIdentity("W", "25"),),
                                    observed_atoms=tuple(a for a in fixed if not (a.chain == "W" and a.residue_number == "25" and a.atom_name == "O")))
    assert all("incomplete" in x.fallback_reason for x in result)
    small = build_anchor_restraints(text, (ResidueIdentity("W", "25"),),
                                   observed_atoms=tuple(a for a in fixed if a.chain == "W" and a.residue_number == "25"))
    assert all("insufficient" in x.fallback_reason for x in small)


def test_constant_distribution_and_clipping():
    original = SOURCE.read_text()
    constant = "".join(l[:60] + " 50.00" + l[66:] if l.startswith("ATOM  ") else l
                       for l in original.splitlines(keepends=True))
    args = {"anchors": (ResidueIdentity("W", "25"),), "observed_atoms": observed(original)}
    assert all("constant" in r.fallback_reason for r in build_anchor_restraints(constant, **args))
    for value, expected in (("  1.00", 20000), ("999.00", 5000)):
        text = "".join(l[:60] + value + l[66:] if l.startswith("ATOM  ") and l[21] == "W"
                       and l[22:26].strip() == "25" else l for l in original.splitlines(keepends=True))
        assert all(r.k_kj_mol_nm2 == expected for r in build_anchor_restraints(text, **args))


def test_identity_case_and_insertion_codes_preserved():
    text = SOURCE.read_text()
    selected = "".join(l for l in text.splitlines(keepends=True) if l.startswith("ATOM  ")
                       and l[21] == "W" and l[22:26].strip() == "25")
    text += "".join(l[:21] + "w" + l[22:26] + "A" + l[27:] for l in selected.splitlines(keepends=True))
    result = build_anchor_restraints(text, (ResidueIdentity("w", "25", "A"),), observed_atoms=observed(text))
    assert len(result) == 6 and all(r.atom.chain == "w" and r.atom.insertion_code == "A" for r in result)


@pytest.mark.parametrize("kwargs", [{"k_kj_mol_nm2": 0}, {"k_kj_mol_nm2": float("nan")},
                                    {"reference_angstrom": (0, float("inf"), 0)}])
def test_restraint_rejects_invalid_parameters(kwargs):
    values = {"atom": AtomIdentity("W", "25", "", "CA"), "reference_angstrom": (0, 0, 0), "k_kj_mol_nm2": 10000}
    with pytest.raises(ValueError):
        AnchorPositionRestraint(**(values | kwargs))


def test_per_anchor_movement_not_diluted_over_structure_or_aligned_away():
    text = SOURCE.read_text()
    restraints = build_anchor_restraints(text, (ResidueIdentity("W", "25"), ResidueIdentity("X", "139")),
                                         observed_atoms=observed(text))
    xyz = {r.atom: r.reference_angstrom for r in restraints}
    for r in restraints:
        if r.atom.chain == "X":
            x, y, z = xyz[r.atom]
            xyz[r.atom] = (x + 1, y, z)
    # Arbitrary test bounds, not scientifically calibrated production defaults.
    results = check_anchor_movements(restraints, xyz, AnchorMovementLimits(0.5, 0.8))
    assert results[0].within_limits and results[0].backbone_rmsd_angstrom == 0
    assert not results[1].within_limits and results[1].backbone_rmsd_angstrom == pytest.approx(1)
    assert results[1].heavy_displacement_angstrom_max == pytest.approx(1)


def test_movement_limits_require_complete_finite_measurements():
    atom = AtomIdentity("W", "25", "", "CA")
    restraint = AnchorPositionRestraint(atom, (0, 0, 0), 1000)
    limits = AnchorMovementLimits(0.5, 1)
    with pytest.raises(ValueError, match="missing or nonfinite"):
        check_anchor_movements((restraint,), {}, limits)
    with pytest.raises(ValueError, match="complete backbone"):
        check_anchor_movements((restraint,), {atom: (0, 0, 0)}, limits)
    with pytest.raises(ValueError):
        AnchorMovementLimits(float("inf"), 1)


def test_approved_default_limits_and_exact_boundary():
    limits = AnchorMovementLimits()
    assert limits.backbone_rmsd_angstrom_max == 0.5
    assert limits.ca_displacement_angstrom_max == 0.5
    assert limits.heavy_displacement_angstrom_max == 1.25
    atoms = tuple(AtomIdentity("W", "29", "", name) for name in ("N", "CA", "C", "O", "CB"))
    restraints = tuple(AnchorPositionRestraint(atom, (0, 0, 0), 10000) for atom in atoms)
    coordinates = {atom: (0.5, 0, 0) if atom.atom_name in {"N", "CA", "C", "O"}
                   else (1.25, 0, 0) for atom in atoms}
    result, = check_anchor_movements(restraints, coordinates, limits)
    assert result.within_limits
    assert result.ca_displacement_angstrom == 0.5


@pytest.mark.parametrize("failure", ["ca", "backbone", "heavy"])
def test_each_default_limit_independently_blocks_excess_motion(failure):
    atoms = tuple(AtomIdentity("W", "29", "", name) for name in ("N", "CA", "C", "O", "CB"))
    restraints = tuple(AnchorPositionRestraint(atom, (0, 0, 0), 10000) for atom in atoms)
    coordinates = {atom: (0, 0, 0) for atom in atoms}
    for atom in atoms:
        if failure == "ca" and atom.atom_name == "CA":
            coordinates[atom] = (0.500001, 0, 0)
        elif failure == "backbone" and atom.atom_name in {"N", "C", "O"}:
            coordinates[atom] = (0.6, 0, 0)
        elif failure == "heavy" and atom.atom_name == "CB":
            coordinates[atom] = (1.250001, 0, 0)
    result, = check_anchor_movements(restraints, coordinates, AnchorMovementLimits())
    assert not result.within_limits
    if failure == "ca":
        assert result.backbone_rmsd_angstrom < 0.5 and result.heavy_displacement_angstrom_max < 1.25
    if failure == "backbone":
        assert result.ca_displacement_angstrom == 0 and result.heavy_displacement_angstrom_max < 1.25
