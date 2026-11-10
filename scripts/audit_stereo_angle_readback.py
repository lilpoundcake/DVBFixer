"""Independent angle-record closure and finite-difference force audit.

The finite differences are numerical derivative controls on temporary coordinate
copies, not minimized/rebuilt/published structural candidates. No Context is used.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from scripts.audit_legacy_start_controls import load_atoms
from scripts.audit_stereo_force_projection import resolve, verify, volume_gradient
from scripts.research_dunbrack_prior import sha, write

STEP_NM = 1e-7


def angle_energy(points, reference, stiffness):
    if len(points) != 3 or any(len(p) != 3 or not all(math.isfinite(x) for x in p) for p in points):
        raise ValueError('invalid finite-difference points')
    if not math.isfinite(reference) or not math.isfinite(stiffness) or stiffness < 0:
        raise ValueError('invalid finite-difference angle parameters')
    a, b = [tuple(x-y for x, y in zip(p, points[1])) for p in (points[0], points[2])]
    length = math.hypot(*a)*math.hypot(*b)
    if length <= 1e-24:
        raise ValueError('undefined finite-difference angle')
    cosine = max(-1., min(1., sum(x*y for x, y in zip(a, b))/length))
    return .5*stiffness*(math.acos(cosine)-reference)**2


def directional_force(points, directions, reference, stiffness):
    first = [[x+STEP_NM*d for x, d in zip(p, direction)] for p, direction in zip(points, directions, strict=True)]
    second = [[x-STEP_NM*d for x, d in zip(p, direction)] for p, direction in zip(points, directions, strict=True)]
    value = -(angle_energy(first, reference, stiffness)-angle_energy(second, reference, stiffness))/(2*STEP_NM)
    if not math.isfinite(value):
        raise ValueError('nonfinite directional derivative')
    return value


def audit(plan_path, root, output):
    import openmm
    from openmm.unit import kilojoules_per_mole, radian

    plan_raw = plan_path.read_bytes()
    plan = json.loads(plan_raw)
    verify(plan)
    report_raw = (root/'evidence.json').read_bytes()
    report = json.loads(report_raw)
    if report['plan_sha256'] != sha(plan_raw) or (root/'executed-plan.json').read_bytes() != plan_raw:
        raise ValueError('angle readback plan differs')
    prior = json.loads(Path(plan['prior_readback']).read_bytes())
    expected = [(r['id'], r['seed'], r['arm']) for r in prior['cases']]
    if [(r['id'], r['seed'], r['arm']) for r in report['cases']] != expected:
        raise ValueError('angle readback omitted an enrolled cell')
    differences = []
    term_count = 0
    for cell in report['cases']:
        raw_bytes = (root/cell['raw_artifact']).read_bytes()
        if sha(raw_bytes) != cell['raw_sha256']:
            raise ValueError('raw angle cell changed')
        raw = json.loads(raw_bytes)
        if raw['status'] != cell['status']:
            raise ValueError('raw angle cell status differs')
        if cell['status'] != 'evaluated':
            continue
        source = Path(plan['run_root'])/cell['id']/f"seed-{cell['seed']}"/cell['arm']
        if [(s['phase'], s['side']) for s in cell['states']] != [(1, 'before'), (1, 'after'), (2, 'before'), (2, 'after')]:
            raise ValueError('partial or reordered saved states')
        for public_state, raw_state in zip(cell['states'], raw['states'], strict=True):
            phase, side = public_state['phase'], public_state['side']
            rows = load_atoms(source/f'phase-{phase}-{side}.json.gz')
            positions = [tuple(x/10 for x in row['xyz']) for row in rows]
            system = openmm.XmlSerializer.deserialize((source/f'phase-{phase}-{side}-system.xml').read_text())
            for centre, raw_centre in zip(public_state['centres'], raw_state['centres'], strict=True):
                terms = raw_centre['terms']
                if {k: v for k, v in raw_centre.items() if k != 'terms'} != centre:
                    raise ValueError('complete raw terms do not match the summary')
                indices = resolve(rows, centre)
                focal = set(indices)
                expected_terms = set()
                for group, force in enumerate(system.getForces()):
                    if isinstance(force, openmm.HarmonicAngleForce):
                        for index in range(force.getNumAngles()):
                            a, b, c, *_ = force.getAngleParameters(index)
                            if focal.intersection((a, b, c)):
                                expected_terms.add((group, index))
                if len(terms) != len(expected_terms) or {(t['group'], t['term']) for t in terms} != expected_terms:
                    raise ValueError('touched-angle term inventory differs')
                _, gradient, norm = volume_gradient([positions[i] for i in indices])
                directions = {i: tuple(x/norm for x in g) for i, g in zip(indices, gradient, strict=True)}
                for term in terms:
                    force = system.getForce(term['group'])
                    a, b, c, reference, stiffness = force.getAngleParameters(term['term'])
                    particles = [a, b, c]
                    elements = [rows[i]['element'] for i in particles]
                    if term['particles'] != particles or term['elements'] != elements or term['atoms'] != [rows[i]['key'] for i in particles]:
                        raise ValueError('angle term identities differ')
                    kind = 'H_containing' if any(e in {'H', 'D'} for e in elements) else 'heavy_only'
                    if term['category'] != kind:
                        raise ValueError('angle H/heavy classification differs')
                    numeric = directional_force([positions[i] for i in particles],
                        [directions.get(i, (0., 0., 0.)) for i in particles],
                        reference.value_in_unit(radian), stiffness.value_in_unit(kilojoules_per_mole/radian**2))
                    value = term['projection_kj_mol_nm']
                    if not math.isclose(value, numeric, rel_tol=1e-5, abs_tol=1e-3):
                        raise ValueError(f'angle projection differs from numerical derivative: {value} vs {numeric}')
                    differences.append(abs(value-numeric))
                    term_count += 1
                for kind in ('H_containing', 'heavy_only'):
                    selected = [t for t in terms if t['category'] == kind]
                    totals = centre['categories'][kind]
                    if totals['count'] != len(selected) or not math.isclose(totals['projection_kj_mol_nm'], sum(t['projection_kj_mol_nm'] for t in selected), abs_tol=1e-8):
                        raise ValueError('angle category totals differ')
                total = sum(t['projection_kj_mol_nm'] for t in terms)
                if not math.isclose(total, centre['prior_projection_kj_mol_nm'], rel_tol=1e-8, abs_tol=1e-6):
                    raise ValueError('angle sum differs from the recorded OpenMM reference')
                selected = sorted((t for t in terms if t['projection_kj_mol_nm'] < 0),
                    key=lambda t: (t['projection_kj_mol_nm'], t['group'], t['term']))[:3]
                if selected != centre['largest_flattening_terms']:
                    raise ValueError('published inspection terms differ from declared rule')
    summary = dict(report['summary'], finite_difference_terms_checked=term_count,
        maximum_finite_difference_error_kj_mol_nm=max(differences, default=None))
    write(output, dict(schema='dvbfixer.stereo-angle-terms-audit.v1', cases=report['cases'], summary=summary,
        plan_sha256=sha(plan_raw), run_sha256=sha(report_raw), audit_script_sha256=sha(Path(__file__).read_bytes()),
        finite_difference_step_nm=STEP_NM, analytic_saved_only=True, native_coordinates_read=False,
        production_changed=False, scientific_qualification=False, original_thr2_causality_resolved=False))
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--readback', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    audit(args.plan, args.readback, args.output)
