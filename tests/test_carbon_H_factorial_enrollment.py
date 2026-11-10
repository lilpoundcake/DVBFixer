"""Factorial cohort stays fixed without private coordinates or scientific calls."""
import json

import pytest
from scripts.enroll_carbon_H_factorial import ARMS, SOURCES, STRESS, enroll
from scripts.research_dunbrack_prior import sha


def fixture(tmp_path):
    legacy = tmp_path/'legacy.json'
    legacy.write_text(json.dumps(dict(cases=[dict(id=label, seed=seed, source_arm=arm,
        source=f'old/{label}/{seed}/{arm}', path=f'old/{label}/{seed}/{arm}/path.json.gz')
        for label in STRESS for seed in (1,2) for arm in ('raw','beta_only','full_generated')])))
    capture = tmp_path/'capture.json'
    capture.write_text(json.dumps(dict(schema='dvbfixer.methyl-diversity-capture-plan.v2',cases=[dict(id=label) for label,_ in SOURCES])))
    root = tmp_path/'capture'
    root.mkdir()
    (root/'evidence.json').write_text(json.dumps(dict(plan_sha256=sha(capture.read_bytes()),cases=[
        dict(id=label,seed=seed,status='refused_source' if label in {'1UBQ','2VLN','2VLQ'} else 'captured')
        for label,_ in SOURCES for seed in (1,2)])))
    return legacy,capture,root


def test_all_stress_conditions_and_source_refusals_are_enrolled_without_cleanup(tmp_path):
    legacy,capture,root = fixture(tmp_path)
    output = tmp_path/'factorial.json'
    enroll(legacy,capture,root,output)
    record = json.loads(output.read_bytes())
    assert len(record['cases']) == 40
    assert record['summary']['potential_replays'] == 136 and record['summary']['inherited_refused_arm_cells'] == 24
    assert ARMS == ('unchanged','all_methyl_zero','all_methylene','combined_zero')
    assert sum(c['status']=='refused_source' for c in record['cases']) == 6
    assert sum(c['context_kind']=='geometry-reconstruction-stress' for c in record['cases']) == 24
    assert not record['coordinates_read'] and not record['groups_discovered'] and not record['arms_evaluated']
    assert record['policy']['H_seeds'] == [1,2] and not record['production_changed']


def test_omitted_stress_context_is_not_silently_replaced(tmp_path):
    legacy,capture,root = fixture(tmp_path)
    record = json.loads(legacy.read_bytes())
    record['cases'].pop()
    legacy.write_text(json.dumps(record))
    with pytest.raises(ValueError,match='ALL24'):
        enroll(legacy,capture,root,tmp_path/'factorial.json')
