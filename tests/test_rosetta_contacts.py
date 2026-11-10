"""Backbone-only selection and complete pre-reconstruction enrollment."""
import copy
import json
import sys
from dataclasses import replace

import pytest
from scripts import research_rosetta_contacts as module


def test_ranking_ignores_sidechain_coordinates():
    original = module.project(module.Path('tests/fixtures/protein_context/scoring_followup/1AKI.pdb').read_bytes())
    expected = module.ranked_pairs(original)
    assert expected
    changed = []
    for line in original.splitlines(keepends=True):
        if line.startswith(b'ATOM  ') and line[12:16].strip() != b'CA':
            line = line[:30]+b' 999.000 999.000 999.000'+line[54:]
        changed.append(line)
    assert module.ranked_pairs(b''.join(changed)) == expected
    assert all(p['ca_distance_angstrom'] <= 8 and p['file_order_indices'][1]-p['file_order_indices'][0] >= 3 for p in expected)


def test_enrollment_never_previews_reconstruction_or_scores(tmp_path, monkeypatch):
    from dvbfixer.protein_completion import mixed_cohort

    def never(*args, **kwargs):
        raise AssertionError('enrollment must not preview outcomes')

    monkeypatch.setattr(mixed_cohort, 'plan_mixed_cohort', never)
    monkeypatch.setattr(mixed_cohort, 'complete_mixed_cohort', never)
    before = 'pyrosetta' in sys.modules
    module.enroll(tmp_path/'enrollment')
    record = json.loads((tmp_path/'enrollment/enrollment.json').read_bytes())
    assert ('pyrosetta' in sys.modules) == before
    assert len(record['cases']) == sum(min(4, len(rows))*3 for rows in record['input_geometry_rankings'].values())
    for system, rows in record['input_geometry_rankings'].items():
        assert rows == sorted(rows, key=lambda p: (p['ca_distance_angstrom'], p['file_order_indices']))
        for rank, pair in enumerate(rows[:4]):
            cases = [c for c in record['cases'] if c['system'] == system and c['selection_rank'] == rank]
            assert len(cases) == 3
            assert all(c['selection_geometry'] == pair for c in cases)
    with pytest.raises(FileExistsError):
        module.enroll(tmp_path/'enrollment')


def test_independent_selector_replay_rejects_missing_case_even_with_valid_payload_hashes(tmp_path):
    from scripts.audit_rosetta_contacts import audit_enrollment

    root = tmp_path/'enrollment'
    module.enroll(root)
    report = audit_enrollment(root)
    assert report['eligible_pairs'] == {'1CRN': 0, '1AKI': 4}
    path = root/'enrollment.json'
    record = json.loads(path.read_bytes())
    record['cases'].pop()
    path.write_text(json.dumps(record))
    with pytest.raises(ValueError, match='omitted or duplicated'):
        audit_enrollment(root)


def test_native_measured_geometry_is_not_a_staggered_generator_candidate(tmp_path):
    from scripts.audit_rosetta_contacts import measured_native_audits

    from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef
    from dvbfixer.pdb_component_reconstruction import _component_ref, _parse_atom_lines
    from dvbfixer.protein_completion.contracts import Addition, Candidate
    from dvbfixer.protein_completion.mixed_cohort import _joint, plan_mixed_cohort

    root = tmp_path/'enrollment'
    module.enroll(root)
    record = json.loads((root/'enrollment.json').read_bytes())
    case = record['cases'][0]
    plan = plan_mixed_cohort((root/case['source']).read_bytes(), tuple(ComponentInstanceRef(**t['residue']) for t in case['targets']))
    atoms = _parse_atom_lines((root/case['reference']).read_text().splitlines(keepends=True))
    native = []
    for target in plan.targets:
        request = target.request
        points = {a.name: a.position for a in atoms if _component_ref(a) == request.residue}
        native.append(Candidate(request.digest, tuple(Addition(ExactAtomRef(request.residue, n), n[0], points[n]) for n in request.missing_atoms), None))
    before = copy.deepcopy(native)
    generator = _joint(plan, tuple(native))
    measured = measured_native_audits(plan, tuple(native), rounded=False)
    for old, new in zip(generator, measured, strict=True):
        assert 'torsion_metadata_or_grid_mismatch' in old.findings
        assert set(new['findings']) == set(old.findings)-{'torsion_metadata_or_grid_mismatch'}
    assert any('generated_clash' in finding for a in measured for finding in a['findings'])
    assert native == before
    changed = tuple([replace(native[0], additions=tuple(replace(a, xyz=(999., 999., 999.)) for a in native[0].additions)), *native[1:]])
    assert not all(a['accepted'] for a in measured_native_audits(plan, changed, rounded=False))
