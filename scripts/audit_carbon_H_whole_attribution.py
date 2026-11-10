"""Independent saved-force sums and finite-difference volume projections.

No Context, minimizer or coordinate mutation. Complete native arrays stay private.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from scripts.audit_legacy_stereo_trace import determinant
from scripts.carbon_H_capture_adapter import read_json
from scripts.research_dunbrack_prior import sha, write
from scripts.research_ha_minimizer_path import verify


def scalar_volume(points):
    return determinant(dict(zip(('origin', 'a', 'b', 'c'), points, strict=True)), 'origin', 'a', 'b', 'c')


def independent_gradient(points):
    gradient = []
    for i in range(4):
        row = []
        for j in range(3):
            plus, minus = [list(p) for p in points], [list(p) for p in points]
            plus[i][j] += 1e-7
            minus[i][j] -= 1e-7
            row.append((scalar_volume(plus)-scalar_volume(minus))/2e-7)
        gradient.append(row)
    return gradient


def require_close(actual, expected, message, relative=1e-8, absolute=1e-5):
    if not math.isfinite(actual) or not math.isfinite(expected) or not math.isclose(actual, expected, rel_tol=relative, abs_tol=absolute):
        raise ValueError(message)
    return abs(actual-expected)


def projection(array, indices, gradient, norm):
    return sum(array[i][j]*g[j] for i, g in zip(indices, gradient) for j in range(3))/norm


def audit(plan_path, root, output):
    raw = plan_path.read_bytes()
    plan = json.loads(raw)
    verify(plan)
    run_raw = (root/'evidence.json').read_bytes()
    run = json.loads(run_raw)
    expected = [[c['id'], c['seed'], c['conditioning'], c['arm']] for c in plan['cases']]
    if (run['plan_sha256'] != sha(raw) or (root/'executed-plan.json').read_bytes() != raw
            or [c['case_id'] for c in run['cases']] != expected):
        raise ValueError('complete saved force-attribution enrollment differs')
    native_errors, projection_errors, gradient_errors = [], [], []
    for case, cell in zip(plan['cases'], run['cases'], strict=True):
        if [f['stage'] for f in cell['frames']] != ['initial', 'phase1', 'phase2']:
            raise ValueError('saved force-state inventory differs')
        rows = case['rows']
        by_identity = {tuple(r['key']):i for i, r in enumerate(rows)}
        if len(by_identity) != len(rows):
            raise ValueError('ambiguous attribution identity')
        for frame in cell['frames']:
            artifact = root/frame['raw_artifact']
            if sha(artifact.read_bytes()) != frame['raw_sha256']:
                raise ValueError('retained native force arrays changed')
            source = read_json(Path(case['source'])/f"{frame['stage']}.json.gz")
            points = source['positions_nm']
            forces = read_json(artifact)
            original, components = forces['original_forces_kj_mol_nm'], forces['component_forces_kj_mol_nm']
            if (len(components) != len(frame['force_components']) or len(original) != len(rows)
                    or any(len(f) != len(rows) for f in components)
                    or any(len(vector) != 3 or not all(math.isfinite(x) for x in vector)
                        for array in (original, *components) for vector in array)):
                raise ValueError('all-particle native force inventory malformed')
            for i in range(len(rows)):
                for j in range(3):
                    native_errors.append(require_close(sum(f[i][j] for f in components), original[i][j], 'native force-vector closure differs'))
            require_close(frame['original_energy_kj_mol'], source['energy_kj_mol'], 'native energy snapshot differs', 1e-10, 1e-6)
            require_close(sum(c['energy_kj_mol'] for c in frame['force_components']), frame['original_energy_kj_mol'], 'native component energy closure differs', 1e-10, 1e-6)
            if [c['label'] for c in frame['centres']] != [c['label'] for c in case['centres']]:
                raise ValueError('complete focal centre inventory differs')
            for declared, centre in zip(case['centres'], frame['centres'], strict=True):
                indices = [by_identity[tuple(declared['residue']+[n])] for n in declared['names']]
                if indices != centre['indices'] or centre['residue'] != declared['residue']:
                    raise ValueError('full focal identity/index differs')
                focal = [points[i] for i in indices]
                require_close(centre['volume_angstrom3'], 1000*scalar_volume(focal), 'independent signed volume differs', 1e-10, 1e-10)
                numerical = independent_gradient(focal)
                analytic = centre['gradient_nm2']
                for observed, exact in zip(numerical, analytic, strict=True):
                    for a, b in zip(observed, exact, strict=True):
                        gradient_errors.append(require_close(a, b, 'independent Cartesian derivative differs', 1e-6, 1e-9))
                norm = math.sqrt(sum(x*x for g in analytic for x in g))
                require_close(norm, centre['gradient_norm_nm2'], 'gradient normalization differs', 1e-10, 1e-12)
                total = projection(original, indices, analytic, norm)
                projection_errors.append(require_close(total, centre['original_native_projection_kj_mol_nm'], 'native signed-volume projection differs'))
                values = [projection(f, indices, analytic, norm) for f in components]
                if len(values) != len(centre['component_projections_kj_mol_nm']):
                    raise ValueError('native force component projection omitted')
                for a, b in zip(values, centre['component_projections_kj_mol_nm'], strict=True):
                    require_close(a, b, 'force component projection differs')
                projection_errors.append(require_close(sum(values), total, 'projected native force sum differs'))
                guard = centre['guard']
                if guard is not None:
                    q = scalar_volume(focal)/guard['vref_nm3']
                    expected_force = guard['k']*max(guard['qmin']-q, 0)/guard['vref_nm3']*norm
                    require_close(values[guard['group']], expected_force, 'analytic guard projection differs')
    summary = dict(paths=len(run['cases']), states=sum(len(c['frames']) for c in run['cases']),
        focal_measurements=sum(len(f['centres']) for c in run['cases'] for f in c['frames']),
        maximum_native_force_closure_error_kj_mol_nm=max(native_errors),
        maximum_projection_closure_error_kj_mol_nm=max(projection_errors),
        maximum_independent_gradient_error_nm2=max(gradient_errors))
    verify(plan)
    write(output, dict(schema='dvbfixer.carbon-H-whole-attribution-audit.v1', cases=run['cases'], summary=summary,
        plan_sha256=sha(raw), run_sha256=sha(run_raw), policy=plan['policy'],
        audit_script_sha256=sha(Path(__file__).read_bytes()), production_changed=False,
        minimization_executed=False, scientific_qualification=False, original_thr2_causality_resolved=False))
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    audit(args.plan, args.run, args.output)
