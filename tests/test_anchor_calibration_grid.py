from types import SimpleNamespace

import numpy as np
import pytest
from openmm import NonbondedForce, System, XmlSerializer
from openmm.app import Topology, element
from openmm.unit import dalton
from scripts import anchor_calibration_grid as grid
from scripts.anchor_calibration_preflight import variants

from dvbfixer.model.diffusion.contract import AtomIdentity, ResidueIdentity


def line(serial, name, number, point):
    return (f"ATOM  {serial:5d} {name:^4s} GLY a{number:4d}    "
            + "".join(f"{x:8.3f}" for x in point)
            + f"{1.:6.2f}{10.:6.2f}          {name[0]:>2s}  \n")


def tiny_case():
    topology, system = Topology(), System()
    chain = topology.addChain("a")
    nonbonded = NonbondedForce()
    names, positions, lines, identities = [], [], [], []
    for number in (1, 2, 3):
        residue = topology.addResidue("GLY", chain, str(number))
        atoms = {}
        for name, elem, delta in (("N", element.nitrogen, (0, 0, 0)), ("CA", element.carbon, (1, 1, 0)),
                                  ("C", element.carbon, (2, 0, 0)), ("O", element.oxygen, (2, -1, 0))):
            atoms[name] = topology.addAtom(name, elem, residue)
            point = np.array(delta)+np.array([3*(number-1), 0, 0])
            positions.append(point/10)
            lines.append(line(len(lines)+1, name, number, point))
            identities.append(AtomIdentity("a", str(number), "", name))
            system.addParticle(elem.mass)
            nonbonded.addParticle(0., .3, .1)
        for first, second in (("N", "CA"), ("CA", "C"), ("C", "O")):
            topology.addBond(atoms[first], atoms[second])
        if names:
            topology.addBond(names[-1]["C"], atoms["N"])
        names.append(atoms)
    system.addForce(nonbonded)
    residues = [ResidueIdentity("a", str(i)) for i in (1, 2, 3)]
    request = SimpleNamespace(
        target_sequences=(SimpleNamespace(chain="a", sequence="GGG"),),
        sequence_placements=(SimpleNamespace(chain="a", observed_target_indices=(0, 2), observed_residues=(residues[0], residues[2])),),
        gaps=(SimpleNamespace(chain="a", target_interval=SimpleNamespace(start=1, stop=2),
                              generated_residues=(residues[1],), left_anchor=residues[0], right_anchor=residues[2]),),
        fixed_atoms=tuple(a for a in identities if a.residue_number != "2"),
        generated_atoms=tuple(a for a in identities if a.residue_number == "2"))
    source = "".join(l for l in lines if l[22:26].strip() != "2")+"END\n"
    return request, source, "".join(lines)+"END\n", topology, system, np.array(positions)


def test_all25_cells_reuse_start_and_freeze_only_correct_observations():
    from dvbfixer.model.diffusion.geometry_conditioning import FULL_NONBONDED_GROUP

    request, source, proposal, topology, base, initial = tiny_case()
    raw = XmlSerializer.serialize(base)
    before = initial.copy()
    for cell in variants():
        system = XmlSerializer.deserialize(raw)
        restraints, limits, mobile = grid.build_cell(request, source, proposal, topology, system, initial, cell)
        assert set(request.generated_atoms) <= mobile
        assert len(restraints) == (0 if cell["window"] == 0 else 8)
        fixed = {grid.identity(a) for a in topology.atoms() if system.getParticleMass(a.index).value_in_unit(dalton) == 0}
        assert fixed == set(request.fixed_atoms)-mobile
        assert any(isinstance(f, NonbondedForce) and f.getForceGroup() == FULL_NONBONDED_GROUP for f in system.getForces())
        assert np.array_equal(initial, before)
        if cell["envelope"]:
            assert limits.backbone_rmsd_angstrom_max == cell["envelope"][0]
            assert limits.ca_displacement_angstrom_max == cell["envelope"][1]
            assert limits.heavy_displacement_angstrom_max == cell["envelope"][2]
    assert XmlSerializer.serialize(base) == raw


def test_cell_parameters_cannot_be_silently_changed():
    req, source, proposal, topology, system, positions = tiny_case()
    with pytest.raises(ValueError, match="parameters"):
        grid.build_cell(req, source, proposal, topology, system, positions, {**variants()[1], "strength": 123.})
    with pytest.raises(ValueError, match="undeclared"):
        grid.cell_by_id("more-free")


def test_native_and_serialized_motion_are_audited_separately(tmp_path):
    from dvbfixer.model.diffusion.anchor_policy import AnchorMovementLimits, build_anchor_restraints
    from dvbfixer.model.diffusion.runner_refinement import rewrite_generated_coordinates

    req, source, proposal, _, _, _ = tiny_case()
    allowed = {AtomIdentity("a", "1", "", "CA")}
    changed = rewrite_generated_coordinates(proposal, {next(iter(allowed)): np.array([1.25, 1, 0])})
    grid.independent_coordinate_audit(proposal, changed, allowed)
    with pytest.raises(ValueError, match="outside"):
        grid.independent_coordinate_audit(proposal, changed, set())
    source_path, final_path = tmp_path / "source.pdb", tmp_path / "final.pdb"
    source_path.write_text(source)
    final_path.write_text(changed)
    restraints = build_anchor_restraints(source, (ResidueIdentity("a", "1"),), observed_atoms=req.fixed_atoms,
                                         bfactor_source="disabled")
    assert grid.independent_movements(source_path, final_path, restraints, AnchorMovementLimits())[0]["within_limits"]
    assert not grid.independent_movements(source_path, final_path, restraints,
        AnchorMovementLimits(ca_displacement_angstrom_max=.1))[0]["within_limits"]


def test_descriptive_environment_metrics_do_not_invent_secondary_assignments_for_short_fragments():
    req, _, proposal, topology, _, _ = tiny_case()
    coordinates = {a: np.array(row.xyz) for a, row in grid.read_atoms(proposal).items()}
    metrics = grid.observed_environment_metrics(topology, coordinates, coordinates, set(req.fixed_atoms))
    assert metrics["lost_contacts"] == [] and metrics["new_contacts"] == []
    assert metrics["secondary_changes"] == []
    assert "no intersegment/interchain" in metrics["secondary_method"]


def test_shared_preparation_can_only_shorten_existing_executor_budget():
    from scripts.run_mac_modeller_comparison import job_timeout_seconds

    assert job_timeout_seconds({}) == 1800.
    assert job_timeout_seconds({"timeout_seconds": 1740.5}) == 1740.5
    for invalid in (0, -1, 1801, float("nan"), float("inf"), True, "1800"):
        with pytest.raises(ValueError, match="timeout"):
            job_timeout_seconds({"timeout_seconds": invalid})


def test_same_named_case_cannot_substitute_unbound_snapshot_or_output(tmp_path):
    case, snapshot, output = tmp_path / "case", tmp_path / "snapshot", tmp_path / "output"
    cell = grid.cell_by_id("fixed")
    job = {"case_id": "case-fixed", "case": str(case),
           "command": ["worker", "--snapshot", str(snapshot), "--output", str(output)]}
    plan = {"jobs": [job]}
    assert grid.declared_worker_job(plan, case, snapshot, cell, output) == job
    with pytest.raises(ValueError, match="paths"):
        grid.declared_worker_job(plan, tmp_path / "other" / "case", snapshot, cell, output)
    with pytest.raises(ValueError, match="paths"):
        grid.declared_worker_job(plan, case, tmp_path / "different-snapshot", cell, output)
