"""The proposed HA control changes one generated H, not observed heavy geometry."""
import math

import pytest
from scripts.research_ha_start_control import ARMS, initial_points, outward_ha, require_generated_ha


def source():
    rows = [dict(key=['A', '16', '', 'GLN', name], element=element)
            for name, element in [('CA', 'C'), ('N', 'N'), ('C', 'C'), ('CB', 'C'), ('HA', 'H')]]
    radius = math.sqrt(8/9)
    points = [[0., 0., 0.], [radius, 0., 1/3], [-radius/2, radius*math.sqrt(3)/2, 1/3],
        [-radius/2, -radius*math.sqrt(3)/2, 1/3], [0., 0., .1]]
    return rows, points


def test_outward_construction_preserves_radius_and_all_other_coordinates():
    rows, points = source()
    original = [p[:] for p in points]
    corrected, index = outward_ha(rows, points)
    assert index == 4 and corrected[index] == pytest.approx([0., 0., -.1])
    assert math.dist(corrected[0], corrected[4]) == pytest.approx(math.dist(points[0], points[4]))
    assert corrected[:4] == points[:4] == original[:4]
    assert points == original
    assert ARMS == ('unchanged', 'outward_HA')


@pytest.mark.parametrize('fault', ['missing', 'duplicate', 'chain_case', 'element', 'zero_radius', 'nonfinite'])
def test_unsupported_ha_identity_or_geometry_is_not_guessed(fault):
    rows, points = source()
    if fault == 'missing':
        rows, points = rows[:-1], points[:-1]
    elif fault == 'duplicate':
        rows, points = rows+rows[-1:], points+points[-1:]
    elif fault == 'chain_case':
        rows[0]['key'][0] = 'a'
    elif fault == 'element':
        rows[-1]['element'] = 'C'
    elif fault == 'zero_radius':
        points[-1] = points[0][:]
    else:
        points[-1][0] = math.nan
    with pytest.raises(ValueError):
        outward_ha(rows, points)


def test_original_ha_is_never_an_allowed_generated_intervention(tmp_path):
    source_file = tmp_path/'source.pdb'
    source_file.write_text('ATOM      1  HA  GLN A  16       1.000   2.000   3.000  1.00  0.00           H\n')
    with pytest.raises(ValueError, match='original observation'):
        require_generated_ha(source_file)
    source_file.write_text('ATOM      1  CA  GLN A  16       1.000   2.000   3.000  1.00  0.00           C\n')
    require_generated_ha(source_file)


def test_saved_state_must_match_captured_coordinates_exactly():
    openmm = pytest.importorskip('openmm')
    from openmm.unit import nanometer

    system = openmm.System()
    system.addParticle(12)
    integrator = openmm.VerletIntegrator(.001)
    context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName('Reference'))
    try:
        context.setPositions([[.1, .2, .3]]*nanometer)
        state = context.getState(positions=True)
        rows = [dict(xyz=[1., 2., 3.])]
        assert initial_points(rows, state)[0] == pytest.approx([.1, .2, .3])
        rows[0]['xyz'][0] += 1e-10
        with pytest.raises(ValueError, match='starting coordinates'):
            initial_points(rows, state)
    finally:
        del context
        del integrator
