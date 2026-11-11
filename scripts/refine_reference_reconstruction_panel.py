#!/usr/bin/env python3
"""Common generated-only derivatives of completed reference development panels.

No resampling, held-out scoring during refinement, observed motion or threshold
changes. Raw-clean proposals are not rebuilt. Original seed-7 derivatives are
reused; original raw/stock-native results remain authoritative and unchanged.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import traceback
from pathlib import Path

from scripts.audit_mac_reconstruction_smoke import audit
from scripts.refine_mac_reconstruction_smoke import run as refine


def run(root, output):
    panel = root / 'reference-development-extension'
    completion = json.loads((panel / 'completed.json').read_text())
    if completion['executed'] != 30:
        raise ValueError('all declared reference sampling calls must finish first')
    rows = json.loads((panel / 'results.json').read_text())
    output.mkdir(exist_ok=False)
    repair_cases = [r for r in rows if 'audit' in r and
                    not r['audit']['assessments']['source-restored']['geometry_passed']]
    declaration = {'raw_results_sha256': hashlib.sha256((panel / 'results.json').read_bytes()).hexdigest(),
        'purpose': 'post-result common repair arm; original raw failures unchanged',
        'cases_to_repair': [(r['case'], r['arm'], r['seed']) for r in repair_cases],
        'platform': 'Reference', 'refinement_seed': 7, 'max_iterations': 1000, 'restarts': 0,
        'geometry_preconditioning': False, 'all_generated_atoms_mobile': True,
        'observed_atoms_mobile': False, 'reference_coordinates_supplied_to_refiner': False,
        'raw_clean_cases_not_rebuilt': True, 'seed7_derivatives': 'reuse reference-pilot-generated-only-repair, no rerun',
        'qualification': False}
    (output / 'declaration.json').write_text(json.dumps(declaration, indent=2)+'\n')
    (output / 'entry-script.py').write_bytes(Path(__file__).read_bytes())
    results = []
    for row in rows:
        name = f'{row["case"]}-{row["arm"]}-seed{row["seed"]}'
        result = {'case': row['case'], 'arm': row['arm'], 'sampling_seed': row['seed'], 'refinement_seed': 7}
        if 'audit' not in row:
            result['status'] = 'raw-not-auditable-no-coordinate-repair-authorized'
        elif row['audit']['assessments']['source-restored']['geometry_passed']:
            result['status'] = 'raw-clean-not-rebuilt'
        else:
            try:
                work = output / name
                refine(panel / name / 'audit', work)
                result['refinement'] = json.loads((work / 'status.json').read_text())
                case = root / 'reference-panel-preparation' / row['case']
                result['audit'] = audit(case, work / 'refined.pdb', work / 'audit', seed=row['seed'])
            except Exception:
                result['failure'] = traceback.format_exc()
        results.append(result)
        (output / 'results.json').write_text(json.dumps(results, indent=2)+'\n')
    (output / 'completed.json').write_text(json.dumps({'processed': len(results), 'qualified': False})+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    run(args.root.absolute(), args.output.absolute())
