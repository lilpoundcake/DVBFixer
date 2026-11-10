"""H-angle readback is descriptive and separate from existing heavy gates."""
import math

import pytest
from scripts.audit_methyl_whole_control import hydrogen_angles


def data():
    openmm = pytest.importorskip('openmm')
    rows = [dict(key=['A', '16', '', 'GLN', name], element=element)
            for name, element in [('HA', 'H'), ('CA', 'C'), ('CB', 'C'), ('C', 'C')]]
    points = [[1., 0., 0.], [0., 0., 0.], [.5, math.sqrt(3)/2, 0.], [1., 1., 1.]]
    system = openmm.System()
    for _ in rows:
        system.addParticle(1.)
    force = openmm.HarmonicAngleForce()
    force.addAngle(0, 1, 2, math.pi/2, 100.)
    force.addAngle(1, 2, 3, math.pi/2, 100.)
    system.addForce(force)
    return system, rows, points


def test_h_angle_inventory_counts_only_h_containing_terms_and_reports_strain():
    system, rows, points = data()
    result = hydrogen_angles(system, rows, points)
    assert result['count'] == 1
    assert result['maximum_absolute_deviation_degrees'] == pytest.approx(30.)
    assert result['rms_deviation_degrees'] == pytest.approx(30.)
    assert result['largest_deviations'][0]['observed_degrees'] == pytest.approx(60.)
    assert result['largest_deviations'][0]['reference_degrees'] == pytest.approx(90.)
    assert result['largest_deviations'][0]['deviation_degrees'] == pytest.approx(-30.)


def test_undefined_h_angle_is_not_reported_as_clean_geometry():
    system, rows, points = data()
    points[0] = points[1][:]
    with pytest.raises(ValueError, match='undefined'):
        hydrogen_angles(system, rows, points)
