"""Native graph, exact identity and generated/mobile ownership fail closed."""
import copy
import gzip
import json

import pytest
from scripts.carbon_H_capture_adapter import adapt, labels, validate_native


def fixture():
    openmm = pytest.importorskip('openmm')
    rows = [dict(key=['d', '82', 'A', 'ALA', n], element=e)
        for n, e in [('CA', 'C'), ('HA', 'H')]]
    points = [[.1, .2, .3], [.2, .2, .3]]
    pre = [dict(rows[0], xyz=[10*x for x in points[0]])]
    system = openmm.System()
    system.addParticle(0.)
    system.addParticle(1.)
    force = openmm.HarmonicBondForce()
    force.addBond(0, 1, .1, 1000.)
    system.addForce(force)
    return rows, [[0, 1]], pre, points, system


def test_exact_native_graph_and_pre_H_ownership():
    rows, bonds, pre, points, system = fixture()
    assert validate_native(rows, bonds, pre, points, points, system) == ([1], [1])


@pytest.mark.parametrize('fault', ['chain_case', 'icode', 'coordinate', 'mobile_heavy', 'fixed_new_H', 'graph', 'endpoint'])
def test_aliases_graph_mass_or_observation_drift_refuse(fault):
    rows, bonds, pre, points, system = fixture()
    final = copy.deepcopy(points)
    if fault == 'chain_case':
        pre[0]['key'] = ['D', *pre[0]['key'][1:]]
    elif fault == 'icode':
        pre[0]['key'] = ['d', '82', '', 'ALA', 'CA']
    elif fault == 'coordinate':
        pre[0]['xyz'][0] += 1e-10
    elif fault == 'mobile_heavy':
        system.setParticleMass(0, 12.)
    elif fault == 'fixed_new_H':
        system.setParticleMass(1, 0.)
    elif fault == 'graph':
        bonds = []
    else:
        final[0][0] += 1e-10
    with pytest.raises(ValueError):
        validate_native(rows, bonds, pre, points, final, system)


def test_ambiguous_complete_atom_identity_refuses():
    with pytest.raises(ValueError, match='ambiguous'):
        labels([dict(key=['d', '82', 'A', 'ALA', 'CA'], element='C')]*2)


@pytest.mark.parametrize('kind', ['geometry-reconstruction-stress', 'prepared-protein-view'])
def test_both_saved_capture_adapters_keep_exact_nm_and_ownership(tmp_path, kind):
    openmm = pytest.importorskip('openmm')

    rows, bonds, pre, points, system = fixture()
    atoms = [dict(r, xyz=[10*x for x in p]) for r, p in zip(rows, points, strict=True)]
    def save(name, data):
        path = tmp_path/name
        raw = json.dumps(data).encode()
        path.write_bytes(gzip.compress(raw) if path.suffix == '.gz' else raw)
        return str(path)
    case = dict(status='captured', context_kind=kind,
        atoms=save('atoms.json.gz', dict(atoms=atoms)),
        topology=save('topology.json', dict(atoms=rows, bonds=bonds)),
        pre_H=save('input.json.gz', dict(atoms=pre)),
        baseline_path=save('path.json.gz', dict(initial_nm=points, final_nm=points)),
        system=str(tmp_path/'system.xml'))
    (tmp_path/'system.xml').write_text(openmm.XmlSerializer.serialize(system))
    save('phase-1-before.json.gz', dict(atoms=atoms))
    save('h-internal-after.json.gz', dict(atoms=atoms))
    case['ownership'] = save('ownership.json', dict(rows=rows, bonds=bonds, new_h=[1], mobile=[1]))
    data = adapt(case)
    assert data['initial_nm'] == data['final_nm'] == points
    assert data['rows'] == rows and data['new_h'] == [1]
    if kind == 'geometry-reconstruction-stress':
        atoms[1]['key'] = ['D', '82', 'A', 'ALA', 'HA']
        save('phase-1-before.json.gz', dict(atoms=atoms))
        with pytest.raises(ValueError, match='label bijection'):
            adapt(case)
    else:
        save('ownership.json', dict(rows=rows, bonds=bonds, new_h=[], mobile=[1]))
        with pytest.raises(ValueError, match='ownership'):
            adapt(case)


def test_inherited_refusal_does_not_open_any_coordinate_file():
    with pytest.raises(ValueError, match='inherited source refusal'):
        adapt(dict(status='refused_source'))


def test_kept_H_remains_fixed_and_excluded_from_new_ownership():
    openmm = pytest.importorskip('openmm')
    rows, bonds, pre, points, system = fixture()
    rows.append(dict(key=['d', '82', 'A', 'ALA', 'HB1'], element='H'))
    points.append([.1, .3, .3])
    pre.append(dict(rows[2], xyz=[10*x for x in points[2]]))
    system.addParticle(0.)
    force = next(f for f in system.getForces() if isinstance(f, openmm.HarmonicBondForce))
    force.addBond(0, 2, .1, 1000.)
    bonds.append([0, 2])
    assert validate_native(rows, bonds, pre, points, points, system) == ([1], [1])
    changed = copy.deepcopy(points)
    changed[2][0] += 1e-10
    with pytest.raises(ValueError, match='fixed observation'):
        validate_native(rows, bonds, pre, points, changed, system)
