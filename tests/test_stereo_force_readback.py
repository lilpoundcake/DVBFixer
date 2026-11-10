"""Independent readback must reject bad arithmetic and changed switch coordinates."""
import copy
import json
import math

import pytest
import scripts.audit_stereo_force_readback as reader
from scripts.research_dunbrack_prior import sha


@pytest.mark.parametrize('bad', [math.nan, math.inf, 5.0])
def test_readback_does_not_accept_nonfinite_or_inconsistent_arithmetic(bad):
    with pytest.raises(ValueError):
        reader.close(1.0, bad)


@pytest.mark.parametrize('corruption', [None, 'coordinates', 'energy', 'phase_scale'])
def test_saved_only_phase_checks(tmp_path, monkeypatch, corruption):
    # Synthetic arithmetic fixture only; source binding is separately exercised
    # by the real retained readback. No scientific files or OpenMM needed here.
    monkeypatch.setattr(reader, 'verify', lambda record: None)
    monkeypatch.setattr(reader, 'historical_verify', lambda record, root: None)
    parent = tmp_path/'parent.json'
    parent.write_text('{}')
    source = tmp_path/'source'
    folder = source/'case'/'seed-1'/'raw'
    folder.mkdir(parents=True)
    (folder/'phase-1-after.json.gz').write_bytes(b'identical saved coordinates')
    (folder/'phase-2-before.json.gz').write_bytes(b'changed' if corruption == 'coordinates' else b'identical saved coordinates')
    state = dict(phase=1, side='before', total_energy_kj_mol=5.,
        force_components=[dict(index=0, type='CustomExternalForce', name='restraint', energy_kj_mol=5.)],
        centres=[dict(label='target_beta', residue=['d', '82', 'A', 'THR'], names=['CB', 'CA', 'OG1', 'CG2'],
            volume_angstrom3=2., angle_degrees=55., component_projections_kj_mol_nm=[2.],
            total_projection_kj_mol_nm=2., guard=None)])
    states = []
    for phase, side in ((1, 'before'), (1, 'after'), (2, 'before'), (2, 'after')):
        current = copy.deepcopy(state)
        current.update(phase=phase, side=side)
        if phase == 2:
            current['total_energy_kj_mol'] = .5
            current['force_components'][0]['energy_kj_mol'] = .5
            current['centres'][0]['component_projections_kj_mol_nm'] = [.2]
            current['centres'][0]['total_projection_kj_mol_nm'] = .2
        states.append(current)
    if corruption == 'energy':
        states[-1]['total_energy_kj_mol'] = .6
    elif corruption == 'phase_scale':
        states[2]['centres'][0]['component_projections_kj_mol_nm'] = [.3]
        states[2]['centres'][0]['total_projection_kj_mol_nm'] = .3
    cell = dict(id='case', seed=1, arm='raw', status='evaluated', states=states)
    plan = dict(inputs={}, source_hashes={}, executions=[dict(id='case', seed=1, arm='raw')],
                run_root=str(source), parent_plan=str(parent), historical_root=str(tmp_path))
    plan_path = tmp_path/'plan.json'
    plan_raw = json.dumps(plan).encode()
    plan_path.write_bytes(plan_raw)
    run = tmp_path/'run'
    run.mkdir()
    (run/'executed-plan.json').write_bytes(plan_raw)
    (run/'cell-01.json').write_text(json.dumps(cell))
    (run/'evidence.json').write_text(json.dumps(dict(plan_sha256=sha(plan_raw), cases=[cell],
        summary=dict(enrolled=1, evaluated=1, readback_blockers=0, states=4, centres=4))))
    output = tmp_path/'audited.json'
    if corruption:
        with pytest.raises(ValueError):
            reader.audit(plan_path, run, output)
        assert not output.exists()
    else:
        reader.audit(plan_path, run, output)
        result = json.loads(output.read_bytes())
        assert result['summary']['phase_switches_checked'] == 1
        assert result['phase_switches'][0]['centres'][0]['net_change_kj_mol_nm'] == pytest.approx(-1.8)
