"""Preparation wrapper preserves borrowed inputs and canonical DatRecord loading."""
import json

import pytest
from scripts import prepare_methyl_diversity as module

from dvbfixer.ffutils.dat import DatRecord


def fixture(tmp_path, monkeypatch, status='admitted'):
    projected = tmp_path/'projection.pdb'
    projected.write_bytes(b'END\n')
    plan = tmp_path/'source-plan.json'
    runtime = {'synthetic': True}
    plan.write_text(json.dumps(dict(inputs={}, source_hashes={}, runtime=runtime,
        cases=[dict(id='synthetic', status=status, projection=str(projected))])))
    output = tmp_path/'prepared'
    output.mkdir()
    monkeypatch.setattr(module, 'preflight', lambda: dict(runtime))
    return plan, output, projected


def test_wrapper_uses_existing_prepare_and_loads_its_dat_without_rewriting_input(tmp_path, monkeypatch):
    plan, output, projected = fixture(tmp_path, monkeypatch)
    calls = []
    def prepared(plan_path, label, folder):
        calls.append((plan_path, label, folder))
        (folder/'prepared.pdb').write_bytes((folder/'input.pdb').read_bytes())
        DatRecord().save(folder/'prepared.dat')
    monkeypatch.setattr(module, 'prepare_worker', prepared)
    module.worker(plan, 'synthetic', output)
    assert calls == [(plan, 'synthetic', output)]
    assert (output/'input.pdb').read_bytes() == projected.read_bytes() == b'END\n'


def test_refused_source_never_enters_prepare(tmp_path, monkeypatch):
    plan, output, _ = fixture(tmp_path, monkeypatch, 'refused_source')
    def forbidden(*args):
        raise AssertionError('must not prepare a source refusal')
    monkeypatch.setattr(module, 'prepare_worker', forbidden)
    with pytest.raises(ValueError, match='refused source'):
        module.worker(plan, 'synthetic', output)
    assert not (output/'input.pdb').exists()


def test_prepare_input_mutation_is_a_blocker_not_an_admitted_result(tmp_path, monkeypatch):
    plan, output, projected = fixture(tmp_path, monkeypatch)
    def mutated(plan_path, label, folder):
        (folder/'input.pdb').write_bytes(b'CHANGED\n')
        (folder/'prepared.pdb').write_bytes(b'END\n')
        DatRecord().save(folder/'prepared.dat')
    monkeypatch.setattr(module, 'prepare_worker', mutated)
    with pytest.raises(ValueError, match='borrowed input'):
        module.worker(plan, 'synthetic', output)
    assert projected.read_bytes() == b'END\n'


def test_resource_drift_fails_before_starting_preparation(tmp_path, monkeypatch):
    plan, output, _ = fixture(tmp_path, monkeypatch)
    monkeypatch.setattr(module, 'preflight', lambda: {'changed': True})
    with pytest.raises(ValueError, match='runtime/resources'):
        module.worker(plan, 'synthetic', output)
    assert not (output/'input.pdb').exists()
