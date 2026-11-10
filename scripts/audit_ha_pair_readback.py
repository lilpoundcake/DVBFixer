"""Independent complete private-pair inventory/sum/publication audit, saved-only."""
from __future__ import annotations

import argparse
import gzip
import json
import math
from pathlib import Path

from scripts.audit_ha_minimizer_path import outward_frame, projection
from scripts.audit_legacy_start_controls import load_atoms
from scripts.research_dunbrack_prior import sha, write
from scripts.research_ha_minimizer_path import verify


def audit(plan_path, root, output):
    raw = plan_path.read_bytes()
    plan = json.loads(raw)
    verify(plan)
    run_raw = (root/'evidence.json').read_bytes()
    run = json.loads(run_raw)
    prior = json.loads(Path(plan['prior']).read_bytes())
    parent = json.loads(Path(plan['path_plan']).read_bytes())
    if run['plan_sha256'] != sha(raw) or (root/'executed-plan.json').read_bytes() != raw:
        raise ValueError('pair decomposition plan differs')
    expected = [(c['id'], c['seed'], c['source_arm']) for c in parent['cases']]
    if [(r['id'], r['seed'], r['source_arm']) for r in run['cases']] != expected:
        raise ValueError('pair decomposition omitted/reordered cells')
    errors, total_pairs = [], 0
    for case, cell, native in zip(parent['cases'], run['cases'], prior['cases'], strict=True):
        if any(tuple(r[k] for k in ('id', 'seed', 'source_arm')) != tuple(case[k] for k in ('id', 'seed', 'source_arm')) for r in (cell, native)):
            raise ValueError('native/pair cell identity differs')
        if cell['status'] != 'evaluated':
            raise ValueError('pair decomposition has an unresolved blocker')
        artifact = root/cell['raw_artifact']
        if sha(artifact.read_bytes()) != cell['raw_sha256']:
            raise ValueError('complete private pair artifact changed')
        rows = load_atoms(Path(case['source'])/'h-internal-before.json.gz')
        focal = [i for i, r in enumerate(rows) if r['key'] == ['A', '16', '', 'GLN', 'HA']]
        if len(focal) != 1:
            raise ValueError('unique HA missing')
        _, outward = outward_frame(rows)
        with gzip.open(artifact, 'rt') as stream:
            for summary, bound in zip(cell['nodes'], native['nodes'], strict=True):
                if (summary['label'], summary['ordinal']) != (bound['label'], bound['ordinal']):
                    raise ValueError('native/pair node identity differs')
                private = json.loads(next(stream))
                terms = private.pop('terms')
                if private != summary:
                    raise ValueError('private pair summary differs from publication')
                expected_indices = set(range(len(rows)))-set(focal)
                if len(terms) != len(expected_indices) or {t['partner_index'] for t in terms} != expected_indices:
                    raise ValueError('complete HA-touched pair inventory differs')
                for term in terms:
                    atom = rows[term['partner_index']]
                    if term['partner'] != atom['key'] or term['element'] != atom['element']:
                        raise ValueError('pair atom identity differs')
                    if not math.isclose(projection(term['force_kj_mol_nm'], outward), term['outward_projection_kj_mol_nm'], rel_tol=1e-10, abs_tol=1e-8):
                        raise ValueError('pair vector and published axis projection differ')
                    if not math.isclose(term['LJ_outward_projection_kj_mol_nm']+term['Coulomb_outward_projection_kj_mol_nm'], term['outward_projection_kj_mol_nm'], rel_tol=1e-10, abs_tol=1e-8):
                        raise ValueError('pair LJ/Coulomb projections do not close')
                force = [sum(t['force_kj_mol_nm'][j] for t in terms) for j in range(3)]
                physical = next(f['HA_force_kj_mol_nm'] for f in bound['force_components'] if f['type'] == 'NonbondedForce')
                for a, b in zip(force, physical, strict=True):
                    if not math.isclose(a, b, rel_tol=1e-8, abs_tol=1e-4):
                        raise ValueError('private pair sum differs from independent unmodified native force')
                    errors.append(abs(a-b))
                if force != summary['total_force_kj_mol_nm'] or len(terms) != summary['pair_count']:
                    raise ValueError('published complete pair totals differ')
                best = sorted((t for t in terms if t['outward_projection_kj_mol_nm'] < 0), key=lambda t: (t['outward_projection_kj_mol_nm'], t['partner_index']))[:3]
                if best != summary['largest_inward_axis_pairs']:
                    raise ValueError('published inspection view differs from frozen three-pair rule')
                for name, subset in (
                        ('H_partner', [t for t in terms if t['element'] in {'H', 'D'}]),
                        ('heavy_partner', [t for t in terms if t['element'] not in {'H', 'D'}]),
                        ('same_residue', [t for t in terms if t['partner'][:4] == rows[focal[0]]['key'][:4]]),
                        ('other_residue', [t for t in terms if t['partner'][:4] != rows[focal[0]]['key'][:4]])):
                    group = summary['categories'][name]
                    if group['count'] != len(subset) or group['outward_projection_kj_mol_nm'] != sum(t['outward_projection_kj_mol_nm'] for t in subset):
                        raise ValueError('complete pair category totals differ')
                total_pairs += len(terms)
            if stream.read():
                raise ValueError('unreported extra private pair nodes')
    if total_pairs != run['summary']['HA_pair_observations']:
        raise ValueError('cohort pair observation count differs')
    summary = dict(run['summary'], independently_checked_pairs=total_pairs,
        independent_maximum_force_closure_error_kj_mol_nm=max(errors, default=None))
    write(output, dict(schema='dvbfixer.ha-nonbonded-pairs-audit.v1', cases=run['cases'], summary=summary,
        plan_sha256=sha(raw), run_sha256=sha(run_raw), audit_script_sha256=sha(Path(__file__).read_bytes()),
        decomposition_script_sha256=plan['source_hashes'][str(Path('scripts/audit_ha_nonbonded_pairs.py').resolve())],
        openmm_version=plan['openmm_version'], protocol=plan['policy'],
        production_changed=False, native_coordinates_read=False, scientific_qualification=False,
        original_thr2_causality_resolved=False))
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--readback', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    audit(args.plan, args.readback, args.output)
