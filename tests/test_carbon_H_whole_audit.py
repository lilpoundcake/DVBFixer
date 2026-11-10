"""Saved whole readback refuses malformed frames and unsafe trace metadata."""
import gzip
import json
import math

import pytest
from scripts.audit_carbon_H_whole import read_frame, trace_contract
from scripts.research_dunbrack_prior import sha


@pytest.mark.parametrize('fault', ['count', 'coordinate', 'energy', 'hash'])
def test_native_frame_inventory_energy_and_hash_are_not_trusted(tmp_path, fault):
    data = dict(positions_nm=[[0., 0., 0.]], energy_kj_mol=0.)
    if fault == 'count':
        data['positions_nm'].append([0., 0., 0.])
    elif fault == 'coordinate':
        data['positions_nm'][0][0] = math.inf
    elif fault == 'energy':
        data['energy_kj_mol'] = math.nan
    raw = gzip.compress(json.dumps(data).encode())
    xml = b'<System/>'
    (tmp_path/'initial.json.gz').write_bytes(raw)
    (tmp_path/'initial-system.xml').write_bytes(xml)
    frame = dict(stage='initial', coordinate_sha256=sha(raw), system_sha256=sha(xml))
    if fault == 'hash':
        frame['coordinate_sha256'] = '0'*64
    with pytest.raises(ValueError):
        read_frame(tmp_path, frame, 1)


@pytest.mark.parametrize('fault', ['reflection', 'restore', 'stage', 'identity'])
def test_whole_trace_never_silently_accepts_recovery_or_wrong_identity(fault):
    case = dict(id='source', seed=1, conditioning='raw')
    trace = dict(case_id=['source', 1, 'raw'], arm='unchanged', plan_sha256='frozen',
        emergency_reflection=False, observations_restored=False, production_changed=False,
        frames=[dict(stage=s) for s in ('initial', 'phase1', 'phase2')])
    trace_contract(trace, case, 'unchanged', 'frozen')
    if fault == 'reflection':
        trace['emergency_reflection'] = True
    elif fault == 'restore':
        trace['observations_restored'] = True
    elif fault == 'stage':
        trace['frames'].pop()
    else:
        trace['case_id'] = ['other', 1, 'raw']
    with pytest.raises(ValueError):
        trace_contract(trace, case, 'unchanged', 'frozen')
