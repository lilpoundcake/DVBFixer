"""Preparation controls refuse silent platform/seed/observer drift."""
import copy

import pytest
from scripts.audit_methyl_diversity_admission import check_controls


def controls():
    return dict(id='case', plan_sha256='digest', hooks_restored=True, native_coordinates_read=False,
        fixer_instances=[dict(platform='Reference')],
        H_calls=[dict(platform='Reference', seed=1, state_restored=True)])


def test_exact_preparation_controls_are_accepted_without_mutation():
    record = controls()
    before = copy.deepcopy(record)
    check_controls(record, 'case', 'digest')
    assert record == before


@pytest.mark.parametrize('fault', ['seed', 'platform', 'state', 'hooks', 'empty', 'binding'])
def test_preparation_control_drift_is_not_silently_accepted(fault):
    record = controls()
    if fault == 'seed':
        record['H_calls'][0]['seed'] = 2
    elif fault == 'platform':
        record['fixer_instances'][0]['platform'] = 'CPU'
    elif fault == 'state':
        record['H_calls'][0]['state_restored'] = False
    elif fault == 'hooks':
        record['hooks_restored'] = False
    elif fault == 'empty':
        record['H_calls'] = []
    else:
        record['plan_sha256'] = 'other'
    with pytest.raises(ValueError):
        check_controls(record, 'case', 'digest')
