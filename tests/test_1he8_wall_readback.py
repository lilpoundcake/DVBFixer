from types import SimpleNamespace

import pytest
from openmm import System, XmlSerializer
from scripts.audit_1he8_traced_mobility import verify_wall_system

from dvbfixer.model.diffusion.anchor_policy import AnchorMovementLimits, AnchorPositionRestraint
from dvbfixer.model.diffusion.contract import AtomIdentity
from dvbfixer.model.diffusion.geometry_conditioning import add_anchor_boundary_walls


def test_actual_wall_system_not_only_json_is_checked(tmp_path):
    residue = SimpleNamespace(chain=SimpleNamespace(id="a"), id="82", insertionCode="A")
    atoms = [SimpleNamespace(index=i, name=name, residue=residue) for i, name in enumerate(("N", "CA", "C", "O"))]
    topology = SimpleNamespace(atoms=lambda: iter(atoms))
    restraints = tuple(AnchorPositionRestraint(AtomIdentity("a", "82", "A", atom.name), (1., 2., 3.), 10000.) for atom in atoms)
    particles = [{"index": a.index, "chain": "a", "residue_number": "82", "insertion_code": "A", "atom_name": a.name} for a in atoms]
    system = System()
    for _ in atoms:
        system.addParticle(12.)
    inner = AnchorMovementLimits(.45, 1.20, .45)
    add_anchor_boundary_walls(system, topology, restraints, inner)
    path = tmp_path / "actual-system.xml"
    path.write_text(XmlSerializer.serialize(system))
    report = verify_wall_system(path, particles, restraints, inner)
    assert report["particle_wall_terms"] == 4 and report["backbone_wall_terms"] == 1
    force = system.getForce(1)
    ids, parameters = force.getBondParameters(0)
    changed = list(parameters)
    changed[1] = .05  # Actual wall no longer agrees with the declared .045nm.
    force.setBondParameters(0, ids, changed)
    path.write_text(XmlSerializer.serialize(system))
    with pytest.raises(ValueError, match="backbone wall/reference differs"):
        verify_wall_system(path, particles, restraints, inner)
