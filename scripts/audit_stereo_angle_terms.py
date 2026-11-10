"""Retrospective harmonic-angle decomposition at every retained force snapshot.

No Context, integration, minimization, H placement or native-reference scoring.
Analytic three-body derivatives must reproduce the previously measured OpenMM
angle-force projection. Full touched-angle records remain in private raw cells;
the summary retains ALL states/centres and both H-containing/heavy-only totals.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from scripts.audit_legacy_start_controls import load_atoms
from scripts.audit_stereo_force_projection import (
    historical_verify,
    projected_force,
    resolve,
    verify,
    volume_gradient,
)
from scripts.research_dunbrack_prior import sha, write


def angle_energy_forces(points, reference_radians, stiffness):
    """Harmonic E/forces for nm coordinates and k in kJ/mol/rad²."""
    if len(points) != 3 or any(len(p) != 3 or not all(math.isfinite(x) for x in p) for p in points):
        raise ValueError('three finite three-dimensional points required')
    if not math.isfinite(reference_radians) or not math.isfinite(stiffness) or stiffness < 0:
        raise ValueError('invalid harmonic-angle parameters')
    a, b = [tuple(x-y for x, y in zip(p, points[1])) for p in (points[0], points[2])]
    ra, rb = math.hypot(*a), math.hypot(*b)
    if min(ra, rb) <= 1e-12:
        raise ValueError('zero angle vector')
    u, v = [tuple(x/r for x in p) for p, r in ((a, ra), (b, rb))]
    cosine = max(-1., min(1., sum(x*y for x, y in zip(u, v))))
    sine = math.sqrt(max(0., 1-cosine*cosine))
    if sine <= 1e-8:
        raise ValueError('undefined near-collinear angle derivative')
    angle = math.acos(cosine)
    delta = angle-reference_radians
    scale = stiffness*delta/sine
    first = tuple(scale*(y-cosine*x)/ra for x, y in zip(u, v))
    third = tuple(scale*(x-cosine*y)/rb for x, y in zip(u, v))
    middle = tuple(-x-y for x, y in zip(first, third))
    energy = .5*stiffness*delta*delta
    if not math.isfinite(energy) or not all(math.isfinite(x) for f in (first, middle, third) for x in f):
        raise ValueError('nonfinite harmonic-angle energy/force')
    return angle, energy, [first, middle, third]


def term_projection(particles, forces, focal_indices, gradient, norm):
    sparse = {i: [0., 0., 0.] for i in focal_indices}
    for i, force in zip(particles, forces, strict=True):
        if i in sparse:
            sparse[i] = [x+y for x, y in zip(sparse[i], force)]
    return projected_force(sparse, focal_indices, gradient, norm)


def hydrogen_class(elements):
    if any(element is None or not element for element in elements):
        raise ValueError('unknown angle element')
    return 'H_containing' if any(e in {'H', 'D'} for e in elements) else 'heavy_only'


def decompose_state(folder, prior):
    import openmm
    from openmm.unit import kilojoules_per_mole, radian

    phase, side = prior['phase'], prior['side']
    rows = load_atoms(folder/f'phase-{phase}-{side}.json.gz')
    system = openmm.XmlSerializer.deserialize((folder/f'phase-{phase}-{side}-system.xml').read_text())
    if len(rows) != system.getNumParticles() or any(system.isVirtualSite(i) for i in range(len(rows))):
        raise ValueError('particle inventory differs or virtual-site forces need separate treatment')
    positions = [tuple(x/10 for x in row['xyz']) for row in rows]
    angles = []
    angle_groups = []
    for group, force in enumerate(system.getForces()):
        if not isinstance(force, openmm.HarmonicAngleForce):
            continue
        if force.usesPeriodicBoundaryConditions():
            raise ValueError('periodic angle derivatives are outside this readback scope')
        angle_groups.append(group)
        for i in range(force.getNumAngles()):
            a, b, c, reference, stiffness = force.getAngleParameters(i)
            angles.append((group, i, [a, b, c], reference.value_in_unit(radian),
                           stiffness.value_in_unit(kilojoules_per_mole/radian**2)))
    if not angle_groups:
        raise ValueError('no captured harmonic-angle group')
    for group in angle_groups:
        if prior['force_components'][group]['type'] != 'HarmonicAngleForce':
            raise ValueError('prior angle-group mapping differs')
    public_centres, raw_centres = [], []
    for centre in prior['centres']:
        indices = resolve(rows, centre)
        if indices != centre['indices']:
            raise ValueError('focal atom indices differ from bound prior readback')
        _, gradient, norm = volume_gradient([positions[i] for i in indices])
        terms = []
        focal = set(indices)
        for group, term, particles, reference, stiffness in angles:
            if not focal.intersection(particles):
                continue
            observed, energy, forces = angle_energy_forces([positions[i] for i in particles], reference, stiffness)
            category = hydrogen_class([rows[i]['element'] for i in particles])
            projection = term_projection(particles, forces, indices, gradient, norm)
            terms.append(dict(group=group, term=term, particles=particles,
                atoms=[rows[i]['key'] for i in particles], elements=[rows[i]['element'] for i in particles],
                category=category, observed_degrees=math.degrees(observed), reference_degrees=math.degrees(reference),
                energy_kj_mol=energy, projection_kj_mol_nm=projection))
        expected = sum(centre['component_projections_kj_mol_nm'][i] for i in angle_groups)
        total = sum(t['projection_kj_mol_nm'] for t in terms)
        if not math.isclose(total, expected, rel_tol=1e-8, abs_tol=1e-6):
            raise ValueError(f'angle terms do not close to saved OpenMM projection: {total} vs {expected}')
        categories = {kind: dict(count=sum(t['category']==kind for t in terms),
            projection_kj_mol_nm=sum(t['projection_kj_mol_nm'] for t in terms if t['category']==kind))
            for kind in ('H_containing', 'heavy_only')}
        summary = dict(label=centre['label'], residue=centre['residue'], names=centre['names'],
            angle_degrees=centre['angle_degrees'], volume_angstrom3=centre['volume_angstrom3'],
            term_count=len(terms), categories=categories, total_projection_kj_mol_nm=total,
            prior_projection_kj_mol_nm=expected, closure_error_kj_mol_nm=total-expected,
            largest_flattening_terms=sorted((t for t in terms if t['projection_kj_mol_nm'] < 0),
                key=lambda t: (t['projection_kj_mol_nm'], t['group'], t['term']))[:3])
        public_centres.append(summary)
        raw_centres.append(dict(summary, terms=terms))
    return dict(phase=phase, side=side, centres=public_centres), dict(phase=phase, side=side, centres=raw_centres)


def plan(parent_plan_path, prior_path, output):
    import openmm

    parent_raw = parent_plan_path.read_bytes()
    parent = json.loads(parent_raw)
    verify(parent)
    historical_verify(json.loads(Path(parent['parent_plan']).read_bytes()), Path(parent['historical_root']))
    prior_raw = prior_path.read_bytes()
    prior = json.loads(prior_raw)
    if prior['plan_sha256'] != sha(parent_raw) or prior['summary']['evaluated'] != 24:
        raise ValueError('complete previous force readback not bound')
    expected = [(r['id'], r['seed'], r['arm']) for r in parent['executions']]
    if [(r['id'], r['seed'], r['arm']) for r in prior['cases']] != expected:
        raise ValueError('previous readback omitted or reordered a cohort cell')
    if any(r['status'] != 'evaluated' or [(s['phase'], s['side']) for s in r['states']] !=
           [(1, 'before'), (1, 'after'), (2, 'before'), (2, 'after')] for r in prior['cases']):
        raise ValueError('previous readback has incomplete cells/states')
    inputs = {**parent['inputs'], str(parent_plan_path): sha(parent_raw), str(prior_path): sha(prior_raw)}
    sources = {**parent['source_hashes'], str(Path(__file__)): sha(Path(__file__).read_bytes())}
    write(output, dict(schema='dvbfixer.stereo-angle-terms-plan.v1', inputs=inputs, source_hashes=sources,
        prior_readback=str(prior_path), parent_plan=str(parent_plan_path), run_root=parent['run_root'],
        openmm_version=openmm.__version__,
        policy=dict(timing='retrospective to known force-group outcomes; frozen before term decomposition',
            enrollment='ALL24 cells,96 states,240 existing focal centres; every touched harmonic-angle term retained privately',
            derivative='analytic harmonic angle Cartesian force; sum must match bound OpenMM Reference angle projection',
            groups='H-containing includes any explicit H/D; otherwise heavy-only, no atom-name or chain normalization',
            limits='no Context/integration/minimization/H placement/repair/native scoring/parameter-grid or causal replay',
            refusal='unknown elements, virtual sites, periodic angles, undefined derivatives or force-sum mismatch retain blocker',
            publication='all centre totals and class counts; exactly three most-negative terms per centre for inspection, complete raw terms retained'),
        production_changed=False, native_coordinates_read=False))


def run(plan_path, output):
    import openmm

    plan_raw = plan_path.read_bytes()
    record = json.loads(plan_raw)
    verify(record)
    parent = json.loads(Path(record['parent_plan']).read_bytes())
    historical_verify(json.loads(Path(parent['parent_plan']).read_bytes()), Path(parent['historical_root']))
    if openmm.__version__ != record['openmm_version']:
        raise ValueError('OpenMM version differs from enrollment')
    prior = json.loads(Path(record['prior_readback']).read_bytes())
    output.mkdir(parents=True, exist_ok=False)
    (output/'executed-plan.json').write_bytes(plan_raw)
    results = []
    for i, case in enumerate(prior['cases'], 1):
        row = {k: case[k] for k in ('id', 'seed', 'arm')}
        row.update(status='evaluated', states=[])
        raw_row = dict(row, states=[])
        folder = Path(record['run_root'])/case['id']/f"seed-{case['seed']}"/case['arm']
        try:
            for state in case['states']:
                summary, raw = decompose_state(folder, state)
                row['states'].append(summary)
                raw_row['states'].append(raw)
        except Exception as error:
            row.update(status='readback_blocker', error_type=type(error).__name__, error=str(error))
            raw_row.update({k: row[k] for k in ('status', 'error_type', 'error')})
        raw_name = f'cell-{i:02d}-angles.json'
        write(output/raw_name, raw_row)
        row.update(raw_artifact=raw_name, raw_sha256=sha((output/raw_name).read_bytes()))
        results.append(row)
    verify(record)
    centres = [c for row in results for state in row['states'] for c in state['centres']]
    summary = dict(enrolled=len(results), evaluated=sum(r['status']=='evaluated' for r in results),
        readback_blockers=sum(r['status']=='readback_blocker' for r in results),
        states=sum(len(r['states']) for r in results), centres=len(centres),
        touched_angle_observations=sum(c['term_count'] for c in centres),
        maximum_projection_closure_error_kj_mol_nm=max((abs(c['closure_error_kj_mol_nm']) for c in centres), default=None))
    write(output/'evidence.json', dict(schema='dvbfixer.stereo-angle-terms.v1', cases=results,
        plan_sha256=sha(plan_raw), summary=summary, analytic_saved_only=True,
        production_changed=False, native_coordinates_read=False, scientific_qualification=False,
        original_thr2_causality_resolved=False))
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('plan')
    p.add_argument('--parent', type=Path, required=True)
    p.add_argument('--prior', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p = sub.add_parser('run')
    p.add_argument('--plan', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.command == 'plan':
        plan(args.parent, args.prior, args.output)
    else:
        run(args.plan, args.output)
