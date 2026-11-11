#!/usr/bin/env python3
"""Five additional coupled-group seeds for the retained 1HE8 context projection.

Original seed7 and its explicit metadata recovery remain unchanged. New joint
groups are unreviewed controls, not covered by old cc94 visual approvals. Other
eight donors enter only post-sampling assessment, never model conditioning.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
import traceback
from dataclasses import asdict
from pathlib import Path

from scripts.materialize_rfd3_local_controls import native_group

from dvbfixer.ffutils.geometry import ChiralityError, assert_all_l
from dvbfixer.model.diffusion.anchor_acceptance import _structure
from dvbfixer.model.diffusion.contract import DiffusionRequest
from dvbfixer.model.diffusion.runner_refinement import rewrite_generated_coordinates
from dvbfixer.model.diffusion.staged_runtime import assessment_adapter


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sampling_arguments(arguments, *, inputs, work, seed):
    """Rewrite semantic arguments, never replace shared path-name prefixes."""
    if sum(p.startswith('inputs=') for p in arguments) != 1 or sum(p.startswith('out_dir=') for p in arguments) != 1:
        raise ValueError('original context command has ambiguous input/output arguments')
    return [f'inputs={inputs / "rfd3-input.json"}' if p.startswith('inputs=') else
            f'out_dir={work / "run/predictions"}' if p.startswith('out_dir=') else
            f'seed={seed}' if p == 'seed=7' else p for p in arguments]


def run(root, output):
    inputs = root / '1he8-rfd3-pair-context-input'
    old = root / '1he8-rfd3-pair-context'
    original = json.loads((old / 'command.json').read_text())
    manifest = json.loads((inputs / 'manifest.json').read_text())
    baseline = root / '1he8-visual-compatible-selection'
    request = DiffusionRequest.from_json((baseline / 'evaluation-mask.json').read_text())
    source, base_text = (baseline / 'source.pdb').read_bytes().decode(), (baseline / 'selected.pdb').read_bytes().decode()
    if manifest['source_sha256'] != request.normalized_pdb.sha256 or hashlib.sha256(source.encode()).hexdigest() != request.normalized_pdb.sha256:
        raise ValueError('context/full-source custody mismatch')
    if manifest['selected_gap_indices'] != [2, 3] or manifest['generated_atoms'] != 201:
        raise ValueError('declared coupled group changed')
    output.mkdir(exist_ok=False)
    seeds = (42, 101, 102, 103, 104)
    engine = Path('.artifacts/mac-reconstruction-2026-11-11/engines/rfd3').absolute()
    files = [inputs / p for p in ('manifest.json', 'conditioning.pdb', 'source-projection.pdb', 'rfd3-input.json')]
    files += [Path(__file__), Path('scripts/trace_rfd3_seed_transport.py'), Path('scripts/rfd3_mps_seeded_research.py'),
              Path('scripts/materialize_rfd3_local_controls.py'), Path('deploy/rfd3/source-component-width.patch'),
              engine / 'models/rfd3/src/rfd3/utils/inference.py']
    declaration = {'purpose': 'complete coupled 1HE8 context-group development panel; no public qualification',
        'new_seeds': seeds, 'seed7': 'retain old context run plus explicit mapping recovery, no resampling',
        'reserved_seeds_unused': True, 'joint_gap_indices': [2, 3], 'group_atomic': True,
        'source_only_conditioning': True, 'model_residues': 168, 'fixed_atoms': 1136, 'generated_atoms': 201,
        'context_entity_cuts': 'ten explicit noncontiguous model entities, not original whole-chain topology',
        'seed_contract': 'private-mps-seed-at-model-v1', 'schedule': '200 configured timesteps, two recycles, sample1',
        'source_width_amendment': 'new metadata uses maintained width patch; original seed7 recovery unchanged; built-input equality was checked, not every downstream feature',
        'visual_approval': 'new RFD3 groups unreviewed; old user cc94 allow-lists do not apply',
        'baseline_other_eight_donors_sha256': digest(baseline / 'selected.pdb'),
        'files': {str(p): digest(p) for p in files}, 'serial_ml': True, 'resource_policy': 'pressure-v2',
        'automatic_retries': False, 'new_refinement': False, 'qualification': False}
    (output / 'declaration.json').write_text(json.dumps(declaration, indent=2)+'\n')
    (output / 'source.diff').write_bytes(subprocess.check_output(['git', 'diff'], cwd=engine))
    (output / 'entry-script.py').write_bytes(Path(__file__).read_bytes())
    results = []
    for seed in seeds:
        work = output / f'seed-{seed}'
        args = sampling_arguments(original[4:], inputs=inputs, work=work, seed=seed)
        command = [sys.executable, 'scripts/run_mac_reconstruction_guarded.py', '--python', original[0],
                   '--module', original[3], '--output', str(work / 'run'), '--resource-policy', 'pressure-v2',
                   '--timeout', '1800', '--', *args]
        env = os.environ.copy()
        env.update(PYTHONPATH=os.pathsep.join(str(Path(p).absolute()) for p in ['.', 'src', 'scripts']),
            DVBFIXER_RFD3_MODEL_SEED=str(seed), DVBFIXER_RFD3_EXPECTED_FIXED='1136',
            DVBFIXER_RFD3_EXPECTED_GENERATED='201', DVBFIXER_RFD3_TRANSPORT_OUTPUT=str(work / 'run/transport'))
        code = subprocess.run(command, env=env).returncode
        row = {'seed': seed, 'returncode': code}
        if code == 0:
            try:
                xyz, custody = native_group(inputs, work / 'run')
                if len(xyz) != 201 or not set(xyz) <= set(request.generated_atoms):
                    raise ValueError('coupled group ownership differs')
                text = rewrite_generated_coordinates(base_text, xyz)
                assessment = assessment_adapter(request, source)(text)
                structure = _structure(text, 'context-group.pdb')
                try:
                    assert_all_l(structure.pdb.topology, structure.pdb.positions)
                    row['all_l'] = True
                except ChiralityError as exc:
                    row['all_l'] = False
                    row['chirality_error'] = str(exc)
                (work / 'raw-assembly.pdb').write_text(text)
                row.update(custody=custody, assessment=asdict(assessment),
                           actionable_findings=sum(i.actionable for i in assessment.issues),
                           group_atomic=True, unreviewed_control=True)
            except Exception:
                row['audit_failure'] = traceback.format_exc()
        results.append(row)
        (output / 'results.json').write_text(json.dumps(results, indent=2)+'\n')
        execution_path = work / 'run/execution.json'
        execution = json.loads(execution_path.read_text()) if execution_path.exists() else {}
        if code != 0 or not execution or execution.get('guard_stop_reason'):
            raise RuntimeError('coupled panel stopped on technical/guard failure; no automatic retries')
    (output / 'completed.json').write_text(json.dumps({'executed': len(results), 'qualification': False})+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    run(args.root.absolute(), args.output.absolute())
