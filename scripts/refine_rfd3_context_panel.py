#!/usr/bin/env python3
"""Bounded independent admission/ranking and shared repair of new coupled groups.

At most two of the five new unreviewed RFD3 groups are refined. Rank independently
reconstructed whole-assembly geometry counts, then seed; no confidence or hidden
reference reward. Segment/surface diagnostics are not qualified vetoes. Original
seed7 remains a separately retained comparator, not a new trial or retry.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import traceback
from pathlib import Path

from scripts.materialize_rfd3_local_controls import native_group
from scripts.refine_mac_reconstruction_smoke import run as refine

from dvbfixer.ffutils.geometry import ChiralityError, assert_all_l
from dvbfixer.model.diffusion.anchor_acceptance import _structure
from dvbfixer.model.diffusion.contract import DiffusionRequest
from dvbfixer.model.diffusion.runner_refinement import rewrite_generated_coordinates
from dvbfixer.model.diffusion.staged_runtime import assessment_adapter


def run(root, output):
    panel = root / '1he8-rfd3-context-development-panel-attempt02'
    completion = json.loads((panel / 'completed.json').read_text())
    if completion['executed'] != 5:
        raise ValueError('all five newly declared coupled calls must complete first')
    declaration = json.loads((panel / 'declaration.json').read_text())
    for filename, expected in declaration['files'].items():
        if hashlib.sha256(Path(filename).read_bytes()).hexdigest() != expected:
            raise ValueError(f'declared context input/runtime changed: {filename}')
    baseline = root / '1he8-visual-compatible-selection'
    request = DiffusionRequest.from_json((baseline / 'evaluation-mask.json').read_text())
    source = (baseline / 'source.pdb').read_bytes().decode()
    text = (baseline / 'selected.pdb').read_bytes().decode()
    if (hashlib.sha256(text.encode()).hexdigest() != declaration['baseline_other_eight_donors_sha256']
            or hashlib.sha256(source.encode()).hexdigest() != request.normalized_pdb.sha256):
        raise ValueError('full assembly/source differs from declared conditioning frame')
    assess = assessment_adapter(request, source)
    rows = json.loads((panel / 'results.json').read_text())
    output.mkdir(exist_ok=False)
    (output / 'entry-script.py').write_bytes(Path(__file__).read_bytes())
    protocol = {'max_new_refinements': 2, 'seeds_considered': [42,101,102,103,104],
        'seed7': 'retained separate comparator, not rerun', 'selection': 'independent source custody + global L + native fixed drift <=0.001 A; rank complete-assembly actionable count, then seed',
        'generated_atoms_mobile_per_trial': len(request.generated_atoms), 'observed_atoms_mobile': False,
        'coupled_gaps': [2,3], 'group_atomic': True, 'unreviewed_controls': True,
        'platform': 'CPU', 'threads': 1, 'refinement_seed': 7, 'phase_caps': [100,75,75,1000], 'restarts': 0,
        'source_results_sha256': hashlib.sha256((panel / 'results.json').read_bytes()).hexdigest(),
        'hidden_reference_or_confidence_selection': False, 'surface_alarms_are_not_qualified_vetoes': True,
        'qualification': False}
    (output / 'protocol.json').write_text(json.dumps(protocol, indent=2)+'\n')
    eligible, admission = [], []
    for row in rows:
        record = {'seed': row['seed'], 'eligible_for_repair': False}
        if row['returncode'] != 0:
            record['reason'] = 'technical failure'
            admission.append(record)
            continue
        try:
            work = panel / f'seed-{row["seed"]}'
            xyz, custody = native_group(root / '1he8-rfd3-pair-context-input', work / 'run')
            if 'custody' not in row or any(custody[k] != row['custody'][k] for k in
                    ('native_sha256', 'input_manifest_sha256', 'metadata_sha256')):
                raise ValueError('native group differs from retained sampling custody')
            if custody['native_observed_drift_max_angstrom'] > .001:
                raise ValueError('native fixed observations exceed declared precision')
            rebuilt = rewrite_generated_coordinates(text, xyz)
            if rebuilt.encode() != (work / 'raw-assembly.pdb').read_bytes():
                raise ValueError('retained assembly differs from independent group materialization')
            assessment = assess(rebuilt)
            if assessment.admission_failures:
                raise ValueError('whole-assembly source admission failure')
            structure = _structure(rebuilt, 'coupled-admission.pdb')
            assert_all_l(structure.pdb.topology, structure.pdb.positions)
            score = sum(i.actionable for i in assessment.issues)
            record.update(eligible_for_repair=True, actionable_findings=score, native_custody=custody)
            eligible.append((score, row['seed'], rebuilt))
        except (ValueError, ChiralityError) as exc:
            record['reason'] = str(exc)
        admission.append(record)
    selected = sorted(eligible, key=lambda r: (r[0], r[1]))[:2]
    (output / 'independent-admission-and-selection.json').write_text(json.dumps({
        'all_five_retained': admission, 'selected_seeds': [r[1] for r in selected],
        'bounded_not_optimal_over_all_donor_combinations': True}, indent=2)+'\n')
    results = []
    for score, seed, candidate in selected:
        work = output / f'seed-{seed}'
        work.mkdir()
        (work / 'source-restored.pdb').write_text(candidate)
        (work / 'evaluation-mask.json').write_text(request.to_json())
        (work / request.normalized_pdb.path).write_text(source)
        (work / 'source.pdb').write_text(source)
        result = {'sampling_seed': seed, 'raw_actionable_findings': score}
        try:
            refine(work, work / 'shared-refinement', platform='CPU', geometry_preconditioning=True,
                   force_shared_refinement=True)
            result['status'] = json.loads((work / 'shared-refinement/status.json').read_text())
        except Exception:
            result['failure'] = traceback.format_exc()
        results.append(result)
        (output / 'results.json').write_text(json.dumps(results, indent=2)+'\n')
    (output / 'completed.json').write_text(json.dumps({'considered':len(rows), 'refined':len(results), 'qualified':False})+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    run(args.root.absolute(), args.output.absolute())
