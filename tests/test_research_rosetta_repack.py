"""Mechanism harness contracts; tests do not need or initialize PyRosetta."""
import hashlib
import json
import sys

import pytest
from scripts import research_rosetta_repack as module


@pytest.fixture
def frozen(tmp_path):
    root = tmp_path/'frozen'
    root.mkdir()
    data = b'ATOM      1  N   VAL A  39       1.000   2.000   3.000  1.00 20.00           N  \nTER\nEND\n'
    cases = [dict(id=name, source=name+'.pdb', reference='ref.pdb') for name in module.CASES]
    artifacts = {name: hashlib.sha256(data).hexdigest() for c in cases for name in (c['source'], c['reference'])}
    for name in artifacts:
        (root/name).write_bytes(data)
    (root/'plan.json').write_text(json.dumps(dict(cases=cases, artifacts=artifacts)))
    return root


def test_plan_is_outcome_independent_and_does_not_import_pyrosetta(tmp_path, frozen):
    imported = 'pyrosetta' in sys.modules
    module.plan(frozen, tmp_path/'plan')
    plan = json.loads((tmp_path/'plan/plan.json').read_text())
    assert [c['id'] for c in plan['cases']] == list(module.CASES)
    assert plan['seeds'] == [7] and plan['scientific_accuracy_benchmark'] is False
    assert ('pyrosetta' in sys.modules) == imported
    for name, digest in plan['artifacts'].items():
        assert hashlib.sha256((tmp_path/'plan'/name).read_bytes()).hexdigest() == digest


def test_movement_keeps_case_insertion_code_residue_and_atom_identity():
    atom = dict(chain='A', resid=82, icode='A', resname='VAL', atom='CG1', xyz=[1., 2., 3.])
    assert module.movement([atom], [atom])['observations_preserved']
    for changed in (dict(atom, chain='a'), dict(atom, icode=' '), dict(atom, atom='CG2'), dict(atom, resname='THR')):
        audit = module.movement([atom], [changed])
        assert not audit['observations_preserved'] and audit['missing_observations'] == [atom]
    moved = module.movement([atom], [dict(atom, xyz=[2., 2., 3.])])
    assert moved['maximum_angstrom'] == 1 and len(moved['moved_above_1e_7_angstrom']) == 1


def test_invalid_seed_and_changed_inputs_refuse_before_pyrosetta_or_output(tmp_path, frozen):
    module.plan(frozen, tmp_path/'plan')
    with pytest.raises(ValueError, match='seed not predeclared'):
        module.probe(tmp_path/'plan', module.CASES[0], tmp_path/'out', auto=False, preserve_cb=False, seed=19)
    (tmp_path/'plan'/f'{module.CASES[0]}.pdb').write_bytes(b'changed')
    with pytest.raises(ValueError, match='planned input changed'):
        module.probe(tmp_path/'plan', module.CASES[0], tmp_path/'out', auto=False, preserve_cb=False, seed=7)
    assert not (tmp_path/'out').exists()
