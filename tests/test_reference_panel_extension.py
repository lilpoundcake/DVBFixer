"""Orchestration declaration only; no checkpoints, subprocesses or model claims."""
import json
from pathlib import Path

from scripts import extend_reference_reconstruction_panel as module
from scripts.prepare_reference_reconstruction_case import prepare

from dvbfixer.model.diffusion.contract import DiffusionRequest

FIXTURE = Path(__file__).parent / 'fixtures/protein_context/1UBQ.pdb'


def test_reference_extension_keeps_seeds_requests_and_output_paths_disjoint(tmp_path, monkeypatch):
    root = tmp_path / 'evidence'
    case = root / 'case'
    prepare(FIXTURE, case, chain='A', first=1, last=72, hidden_first=10, hidden_last=13, case_id='unit-control')
    pilot = root / 'reference-pilot'
    pilot.mkdir()
    jobs = []
    for arm in ('rfd3', 'boltz', 'protenix'):
        work = pilot / arm
        arguments = {
            'rfd3': ['design', 'seed=7', f'out_dir={work}/run/predictions'],
            'boltz': ['predict', 'input.yaml', '--seed', '7', '--out_dir', str(work / 'run/predictions')],
            'protenix': ['--request', str(case / 'evaluation-mask.json'), '--output', str(work / 'run/predictions')],
        }[arm]
        jobs.append({'case': str(case), 'arm': arm, 'work': str(work), 'seed': 7,
                     'expected_fixed': 546, 'expected_generated': 28,
                     'command': ['python', 'guard.py', '--output', str(work / 'run'), '--', *arguments]})
    (pilot / 'plan.json').write_text(json.dumps({'jobs': jobs}))
    original = (pilot / 'plan.json').read_bytes()
    monkeypatch.setattr(module.subprocess, 'check_output', lambda *_args, **_kwargs: b'fake-source-diff')
    output = root / 'extension'
    module.declare(root, output)
    plan = json.loads((output / 'plan.json').read_text())
    assert len(plan['jobs']) == 15
    assert {j['seed'] for j in plan['jobs']} == {42, 101, 102, 103, 104}
    assert len({j['work'] for j in plan['jobs']}) == 15
    assert (pilot / 'plan.json').read_bytes() == original
    for job in plan['jobs']:
        command = job['command']
        assert str(Path(job['work']) / 'run/predictions') in ' '.join(command)
        if job['arm'] == 'protenix':
            request_path = Path(command[command.index('--request')+1])
            request = DiffusionRequest.from_json(request_path.read_text())
            assert request.seeds == (job['seed'],)
            assert (request_path.parent / request.normalized_pdb.path).read_bytes() == (case / 'source.pdb').read_bytes()
        elif job['arm'] == 'rfd3':
            assert f'seed={job["seed"]}' in command
        else:
            assert command[command.index('--seed')+1] == str(job['seed'])
