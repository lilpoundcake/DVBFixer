#!/usr/bin/env python3
"""Independent complete-panel accounting; never choose donors or tune thresholds.

All declared seeds stay in execution and geometry denominators. Held-out errors
are reported for every auditable output, not just accepted or best candidates.
Missing/failed outputs are counted separately, never imputed as zero RMSD.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import statistics
from pathlib import Path

import numpy as np
from scripts.audit_mac_reconstruction_smoke import inventory


def summarize_rows(rows):
    """Pure accounting, including technical failures with no coordinate audit."""
    audited = [row for row in rows if 'audit' in row]
    errors = [r['audit']['heldout_generated_rmsd_angstrom'] for r in audited]
    accepted = [r for r in audited if r['audit']['assessments']['source-restored']['geometry_passed']
                and not r['audit']['assessments']['source-restored']['admission_failures']]
    native_audited = [r for r in audited if 'aligned-native' in r['audit']['assessments']]
    native_accepted = [r for r in native_audited if r['audit']['assessments']['aligned-native']['geometry_passed']
                       and not r['audit']['assessments']['aligned-native']['admission_failures']]
    return {'declared': len(rows), 'executed_successfully': sum(r['returncode'] == 0 for r in rows),
        'independently_audited': len(audited), 'source_restored_geometry_passes': len(accepted),
        'native_source_and_geometry_audited': len(native_audited),
        'native_source_and_geometry_passes': len(native_accepted),
        'native_semantics': 'post-reinjection output for Protenix; untouched proposal before restoration for stock Boltz',
        'geometry_denominator_includes_all_declared_seeds': True,
        'technical_or_audit_failures': [r['seed'] for r in rows if r['returncode'] != 0 or 'audit' not in r],
        'all_audited_heldout_heavy_rmsd_median': statistics.median(errors) if errors else None,
        'all_audited_heldout_heavy_rmsd_range': [min(errors), max(errors)] if errors else None,
        'per_seed': [{'seed': r['seed'], 'returncode': r['returncode'],
                      'heldout_heavy_rmsd': r.get('audit', {}).get('heldout_generated_rmsd_angstrom'),
                      'source_restored_geometry_pass': r in accepted,
                      'native_source_geometry_pass': r in native_accepted if r in native_audited else None}
                     for r in sorted(rows, key=lambda r: r['seed'])],
        'qualified': False}


def run(root, output):
    extension = root / 'reference-development-extension'
    completion = json.loads((extension / 'completed.json').read_text())
    if completion['executed'] != 30:
        raise ValueError('the declared extension has not completed all 30 calls')
    original = json.loads((root / 'reference-pilot/results.json').read_text())
    rows = [{**r, 'seed': 7, 'case': Path(r['case']).name} for r in original]
    rows += json.loads((extension / 'results.json').read_text())
    case_names = sorted({r['case'] for r in rows})
    summary = {'scope': 'two cropped development structures, no independent family validation',
               'reserved_seeds_used': False, 'reference_only_for_posthoc_scoring': True,
               'no_best_seed_selection_or_threshold_tuning': True, 'groups': [], 'qualified': False}
    for name in case_names:
        case = root / 'reference-panel-preparation' / name
        manifest = json.loads((case / 'manifest.json').read_text())
        for filename in ('source.pdb', 'reference-projection-NOT-CONDITIONING.pdb'):
            if hashlib.sha256((case / filename).read_bytes()).hexdigest() != manifest['files'][filename]:
                raise ValueError('declared reference/source projection bytes changed')
        reference = inventory(case / 'reference-projection-NOT-CONDITIONING.pdb')
        source = inventory(case / 'source.pdb')
        axis = tuple(sorted(reference))
        generated = tuple(a for a in axis if a not in source)
        backbone = tuple(a for a in generated if a.atom_name in ('N', 'CA', 'C', 'O'))
        index = {a: i for i, a in enumerate(axis)}
        for arm in ('rfd3', 'boltz', 'protenix'):
            selected = [r for r in rows if r['case'] == name and r['arm'] == arm]
            if sorted(r['seed'] for r in selected) != [7, 42, 101, 102, 103, 104]:
                raise ValueError('panel has duplicate or missing declared seeds')
            group = {'case': name, 'backend': arm, **summarize_rows(selected)}
            backbone_errors = []
            for row in selected:
                if 'audit' not in row:
                    continue
                work = (root / 'reference-pilot' / f'{name}-{arm}' if row['seed'] == 7 else
                        extension / f'{name}-{arm}-seed{row["seed"]}')
                with np.load(work / 'audit/coordinates.npz') as arrays:
                    aligned = arrays['aligned']
                if aligned.shape != (len(axis), 3) or not np.isfinite(aligned).all():
                    raise ValueError('independent coordinate readback has a different/nonfinite axis')
                def rmsd(atoms, coordinates=aligned, lookup=index, expected=reference):
                    return float(np.sqrt(np.mean([np.sum((coordinates[lookup[a]]-expected[a][2])**2) for a in atoms])))
                heavy, bb = rmsd(generated), rmsd(backbone)
                if not np.isclose(heavy, row['audit']['heldout_generated_rmsd_angstrom'], atol=1e-10, rtol=0):
                    raise ValueError('independent held-out error differs from original audit')
                next(r for r in group['per_seed'] if r['seed'] == row['seed'])['heldout_backbone_rmsd'] = bb
                backbone_errors.append(bb)
            group['all_audited_heldout_backbone_rmsd_median'] = statistics.median(backbone_errors) if backbone_errors else None
            summary['groups'].append(group)
    with output.open('x') as stream:
        stream.write(json.dumps(summary, indent=2)+'\n')
    print(json.dumps([{k: g[k] for k in ('case', 'backend', 'source_restored_geometry_passes',
          'native_source_and_geometry_passes',
          'all_audited_heldout_heavy_rmsd_median')} for g in summary['groups']], indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    run(args.root.absolute(), args.output.absolute())
