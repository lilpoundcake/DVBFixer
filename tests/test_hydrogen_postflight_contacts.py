"""O–H nonparent contacts from the unchanged SER regression's second predicate."""
import math

import pytest

OH_LENGTHS = (0.97, 1.2458133889150464)
CONTACTS = (0.719, 0.999, 1.0, 1.01, 1.1, 1.5)


def contact_graph(oh_length=0.97, contact=0.719, parent_symbol='O'):
    from openmm import Vec3, unit
    from openmm.app import Topology, element

    topology = Topology()
    residue = topology.addResidue('SER', topology.addChain('B'), '126', 'A')
    parent = topology.addAtom('OG', element.Element.getBySymbol(parent_symbol), residue)
    hydrogen = topology.addAtom('HG', element.hydrogen, residue)
    carbon = topology.addAtom('CB', element.carbon, residue)
    topology.addAtom('contactO', element.oxygen, residue)
    topology.addBond(parent, hydrogen)
    topology.addBond(parent, carbon)
    theta = math.radians(109.5)
    positions = [Vec3(0., 0., 0.), Vec3(oh_length/10, 0., 0.),
        Vec3(0.14*math.cos(theta), 0., 0.14*math.sin(theta)), Vec3(oh_length/10, contact/10, 0.)]*unit.nanometer
    return topology, positions


def test_valid_length_does_not_hide_close_nonparent_oxygen():
    pytest.importorskip('openmm')
    from dvbfixer.ffutils import geometry

    top, positions = contact_graph()
    findings = geometry.detect_misplaced_hydrogens(top, positions)
    assert len(findings) == 1
    assert findings[0].reason == 'oh_nonparent_contact'
    assert findings[0].coincident_with.name == 'contactO'
    assert findings[0].contact_distance_nm == pytest.approx(0.0719)
    with pytest.raises(geometry.HydrogenGeometryError):
        geometry.repair_misplaced_hydrogens(top, positions)


@pytest.mark.parametrize('oh_length', OH_LENGTHS)
@pytest.mark.parametrize('contact', CONTACTS)
def test_complete_contact_grid_is_safe_or_explicitly_refused(oh_length, contact):
    pytest.importorskip('openmm')
    from openmm import unit

    from dvbfixer.ffutils import geometry

    top, positions = contact_graph(oh_length, contact)
    before = [list(p) for p in positions.value_in_unit(unit.angstrom)]
    invalid_input = oh_length >= 1.2 or contact <= 1.0
    try:
        geometry.repair_misplaced_hydrogens(top, positions)
    except geometry.HydrogenGeometryError:
        assert invalid_input
    else:
        after = [list(p) for p in positions.value_in_unit(unit.angstrom)]
        assert math.dist(after[0], after[1]) < 1.2
        assert math.dist(after[1], after[3]) > 1.0
        if not invalid_input:
            assert after == before
    after = [list(p) for p in positions.value_in_unit(unit.angstrom)]
    assert after[0] == before[0] and after[2:] == before[2:]


def test_contact_policy_does_not_change_nonoxygen_parent_threshold():
    pytest.importorskip('openmm')
    from openmm import unit

    from dvbfixer.ffutils import geometry

    top, positions = contact_graph(1.09, 0.719, 'C')
    before = [list(p) for p in positions.value_in_unit(unit.nanometer)]
    assert geometry.detect_misplaced_hydrogens(top, positions) == []
    assert geometry.repair_misplaced_hydrogens(top, positions) == 0
    assert [list(p) for p in positions.value_in_unit(unit.nanometer)] == before


def test_diagnose_reports_nonparent_contact_not_a_false_bond_length_error():
    pytest.importorskip('openmm')
    from dvbfixer.diagnose.structural import check_misplaced_hydrogens

    top, positions = contact_graph()
    findings = check_misplaced_hydrogens(top, positions)
    assert len(findings) == 1
    assert findings[0].chain == 'B' and findings[0].resid == '126A'
    assert 'contactO' in findings[0].message
    assert 'non-parent' in findings[0].message
