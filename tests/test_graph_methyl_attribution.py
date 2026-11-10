"""Generic focal frame and saved-only pair closure against native readback."""
import gzip
import json
import math

import pytest
from scripts.attribute_graph_methyl_ha import frame, run


def source():
    rows = [dict(key=['a', '82', 'A', 'ASP', n], element=e) for n, e in
        [('CA', 'C'), ('N', 'N'), ('C', 'C'), ('CB', 'C'), ('HA', 'H'), ('OTHER', 'H')]]
    points = [[0., 0., 0.], [.14, 0., 0.], [0., .14, 0.], [0., 0., .14], [-.06, -.06, -.06], [-.19, -.1, -.13]]
    return rows, points


def test_focal_frame_uses_exact_case_and_insertion_code():
    rows, points = source()
    ca, ha, outward = frame(rows, points, ['a', '82', 'A', 'ASP', 'HA'])
    assert (ca, ha) == (0, 4)
    assert outward == pytest.approx([-1/math.sqrt(3)]*3)
    with pytest.raises(ValueError):
        frame(rows, points, ['A', '82', 'A', 'ASP', 'HA'])
    with pytest.raises(ValueError):
        frame(rows, points, ['a', '82', '', 'ASP', 'HA'])


def test_saved_only_native_pair_closure_preserves_source_and_never_minimizes(tmp_path, monkeypatch):
    openmm = pytest.importorskip('openmm')
    rows, points = source()
    system = openmm.System()
    force = openmm.NonbondedForce()
    force.setNonbondedMethod(openmm.NonbondedForce.CutoffNonPeriodic)
    force.setCutoffDistance(1.)
    for i in range(len(rows)):
        system.addParticle(1.)
        force.addParticle(.1 if i%2 else -.1, .1, .1)
    force.addException(0, 4, 0., .1, 0.)
    system.addForce(force)
    folder = tmp_path/'source'
    folder.mkdir()
    xml = openmm.XmlSerializer.serialize(system)
    (folder/'H-system.xml').write_text(xml)
    (folder/'H-ownership.json').write_text(json.dumps(dict(rows=rows, graph=dict(groups=[]))))
    path = folder/'path.json.gz'
    path.write_bytes(gzip.compress(json.dumps(dict(initial_nm=points, final_nm=points, callbacks=[])).encode()))
    plan = tmp_path/'plan.json'
    plan.write_text(json.dumps(dict(inputs={}, source_hashes={}, openmm_version=openmm.__version__,
        target=['a', '82', 'A', 'ASP', 'HA'], policy={},
        cases=[dict(id='synthetic', seed=1, arm='unchanged', source=str(folder), path=str(path))])))
    def forbidden(*args, **kwargs):
        raise AssertionError('saved-only attribution must not minimize')
    monkeypatch.setattr(openmm.LocalEnergyMinimizer, 'minimize', forbidden)
    output = tmp_path/'attribution'
    run(plan, output)
    result = json.loads((output/'evidence.json').read_bytes())
    assert result['summary']['nodes'] == 2 and result['summary']['HA_pair_observations'] == 10
    assert result['summary']['maximum_native_pair_closure_error_kj_mol_nm'] < 1e-7
    assert (folder/'H-system.xml').read_text() == xml
    assert not result['minimization_executed']
