"""Narrow tests for research weighting, identity and OpenMM unit conversion."""
from pathlib import Path

import numpy as np
import pytest
from scripts.cc89_anchor_restraints import add_anchor_force, bfactor_weights

from dvbfixer.model.diffusion.contract import AtomIdentity, ResidueIdentity

SOURCE = Path(__file__).parent / "fixtures/protein_context/9cqd_split.pdb"


def anchor_atoms():
    wanted = {("W", "25"), ("W", "29"), ("X", "139"), ("X", "142")}
    return tuple(AtomIdentity(l[21], l[22:26].strip(), l[26].strip(), l[12:16].strip())
                 for l in SOURCE.read_text().splitlines() if l.startswith("ATOM  ")
                 and (l[21], l[22:26].strip()) in wanted)


def test_weights_use_source_backbone_and_have_equal_atom_weighted_mean():
    atoms = anchor_atoms()
    weights, means = bfactor_weights(SOURCE.read_text(), atoms)
    assert len(atoms) == 24
    assert np.mean(list(weights.values())) == pytest.approx(1)
    assert means[ResidueIdentity("W", "25")] == pytest.approx(88.725)
    assert means[ResidueIdentity("X", "139")] == pytest.approx(97.325)
    assert weights[atoms[0]] > weights[next(a for a in atoms if a.chain == "X")]
    for r in means:
        assert len({weights[a] for a in atoms if (a.chain, a.residue_number, a.insertion_code)
                    == (r.chain, r.residue_number, r.insertion_code)}) == 1


@pytest.mark.parametrize("value", ["  0.00", " -1.00", "   nan"])
def test_reject_unusable_bfactor(value):
    lines = SOURCE.read_text().splitlines(keepends=True)
    text = "".join(l[:60] + value + l[66:] if l.startswith("ATOM  ") and l[21] == "W"
                   and l[22:26].strip() == "25" and l[12:16].strip() == "CA" else l for l in lines)
    with pytest.raises(ValueError, match="finite and positive"):
        bfactor_weights(text, anchor_atoms())


def test_incomplete_backbone_rejected():
    text = "".join(l for l in SOURCE.read_text().splitlines(keepends=True)
                   if not (l.startswith("ATOM  ") and l[21] == "W" and l[22:26].strip() == "25"
                           and l[12:16].strip() == "O"))
    with pytest.raises(ValueError, match="complete"):
        bfactor_weights(text, anchor_atoms())


def test_case_and_insertion_code_are_not_collapsed():
    atoms = anchor_atoms()[:6]
    text = "".join(l for l in SOURCE.read_text().splitlines(keepends=True)
                   if l.startswith("ATOM  ") and l[21] == "W" and l[22:26].strip() == "25")
    other = text.replace(" W  25 ", " w  25A")
    other_atoms = tuple(AtomIdentity("w", a.residue_number, "A", a.atom_name) for a in atoms)
    weights, means = bfactor_weights(text + other, (*atoms, *other_atoms))
    assert len(means) == 2 and len(weights) == 12


def test_restraint_energy_and_zero_at_source():
    openmm = pytest.importorskip("openmm")
    from openmm.app import Topology, element
    from openmm.unit import kilojoule_per_mole, nanometer

    topology = Topology()
    chain = topology.addChain("W")
    res = topology.addResidue("GLY", chain, "25")
    topology.addAtom("CA", element.carbon, res)
    atom = AtomIdentity("W", "25", "", "CA")
    system = openmm.System()
    system.addParticle(12)
    force = add_anchor_force(system, topology, {atom: np.array([10., 20., 30.])}, {atom: 1000.})
    assert force.getNumParticles() == 1
    context = openmm.Context(system, openmm.VerletIntegrator(0.001), openmm.Platform.getPlatformByName("Reference"))
    context.setPositions([[1., 2., 3.]] * nanometer)
    assert context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kilojoule_per_mole) == pytest.approx(0)
    context.setPositions([[1.1, 2., 3.]] * nanometer)
    # 1 Å = 0.1 nm: 0.5 * 1000 * 0.1² = 5 kJ/mol, not 500.
    assert context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kilojoule_per_mole) == pytest.approx(5)
    with pytest.raises(ValueError, match="missing"):
        add_anchor_force(system, topology, {atom: np.zeros(3)}, {AtomIdentity("w", "25", "", "CA"): 1000})
