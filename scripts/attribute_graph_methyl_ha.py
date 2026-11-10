"""Retrospective saved-only ASP49 force/pair attribution on ALL six 1TIM paths.

Fresh readback Contexts only; no minimization, integration, ablation or new H.
Original pair formulas and native-force closure are reused without tuning.
"""
from __future__ import annotations

import argparse
import gzip
import json
import math
from pathlib import Path

from scripts.audit_ha_nonbonded_pairs import pair_force, parameters
from scripts.audit_stereo_force_projection import grouped_system
from scripts.research_dunbrack_prior import sha, write
from scripts.research_ha_minimizer_path import verify


def frame(rows, points, target):
    indices = {tuple(r['key']): i for i, r in enumerate(rows)}
    keys = [(*target[:4], name) for name in ('CA', 'N', 'C', 'CB', target[-1])]
    if len(indices) != len(rows) or any(key not in indices for key in keys):
        raise ValueError('exact focal HA/heavy frame missing')
    ca, n, c, cb, ha = [indices[key] for key in keys]
    summed = [0., 0., 0.]
    for i in (n, c, cb):
        vector = [x-y for x, y in zip(points[i], points[ca])]
        length = math.hypot(*vector)
        if length <= 1e-12:
            raise ValueError('undefined focal heavy bond direction')
        summed = [x+y/length for x, y in zip(summed, vector)]
    norm = math.hypot(*summed)
    if norm <= 1e-12:
        raise ValueError('undefined focal outward frame')
    return ca, ha, [-x/norm for x in summed]


def projection(vector, axis):
    return sum(x*y for x, y in zip(vector, axis, strict=True))


def plan(h_plan, h_root, h_audit, whole_audit, output):
    import openmm

    raw = h_plan.read_bytes()
    parent = json.loads(raw)
    verify(parent)
    run_raw = (h_root/'evidence.json').read_bytes()
    run = json.loads(run_raw)
    audit_raw = h_audit.read_bytes()
    prior = json.loads(audit_raw)
    whole_raw = whole_audit.read_bytes()
    whole = json.loads(whole_raw)
    if prior['plan_sha256'] != sha(raw) or prior['run_sha256'] != sha(run_raw):
        raise ValueError('upstream H audit differs')
    if prior['summary']['statuses'] != {'completed': 30, 'refused_source': 18} or whole['summary']['statuses'] != {'completed': 30, 'refused_source': 18}:
        raise ValueError('complete H/downstream cohort not bound')
    inputs = {**parent['inputs'], str(h_plan): sha(raw), str(h_root/'evidence.json'): sha(run_raw),
        str(h_audit): sha(audit_raw), str(whole_audit): sha(whole_raw)}
    cases = []
    for source in parent['cases']:
        if source['id'] != '1TIM':
            continue
        if source['seed'] not in (1, 2) or source['status'] != 'captured':
            raise ValueError('declared matched 1TIM contexts unavailable')
        for arm in parent['policy']['arms']:
            result = next(r for r in run['cases'] if (r['id'], r['seed'], r['arm']) == ('1TIM', source['seed'], arm))
            if result['status'] != 'completed':
                raise ValueError('matched 1TIM H path has blocker')
            folder = h_root/'1TIM'/f"seed-{source['seed']}"/arm
            for name, digest in result['artifacts'].items():
                path = folder/name
                if sha(path.read_bytes()) != digest:
                    raise ValueError('bound H path changed')
                inputs[str(path)] = digest
            system = openmm.XmlSerializer.deserialize((Path(source['source'])/'H-system.xml').read_text())
            parameters(system)
            cases.append(dict(id='1TIM', seed=source['seed'], arm=arm, source=source['source'], path=str(folder/'path.json.gz')))
    if len(cases) != 6:
        raise ValueError('all three arms and both original seeds required')
    sources = dict(parent['source_hashes'])
    for name in ('scripts/attribute_graph_methyl_ha.py', 'scripts/audit_ha_nonbonded_pairs.py',
                 'scripts/audit_stereo_force_projection.py'):
        sources[name] = sha(Path(name).read_bytes())
    write(output, dict(schema='dvbfixer.graph-methyl-ASP49-attribution-plan.v1', cases=cases, inputs=inputs,
        source_hashes=sources, openmm_version=openmm.__version__, target=['A', '49', '', 'ASP', 'HA'],
        policy=dict(timing='retrospective to known ASP49 H/whole outcomes, frozen before pair-contributor inspection',
            enrollment='all three 1TIM H arms, seeds1/2, every saved initial/callback/final node and every HA-other-particle pair; no favorable endpoint/seed filter',
            measurement='fresh Reference grouped-force readback Context, exact serialized FF restoration after measurement grouping; analytic pair sum closes to original native NonbondedForce',
            physics='shared tested unswitched CutoffNonPeriodic reaction-field Coulomb/Lorentz-Berthelot LJ with exact native exceptions; no force/exclusion/parameter edits',
            reporting='all node group forces/class totals and exactly three most-negative outward pair projections; complete pair terms private; original observed heavy frame fixed',
            limits='no minimization, integration, new coordinates/H/seeds, ablation, native score or production repair; do not assume methyl contributor or original THR2 cause'),
        production_changed=False, scientific_qualification=False))


def run(plan_path, output):
    import openmm
    from openmm.unit import kilojoules_per_mole, nanometer

    raw = plan_path.read_bytes()
    record = json.loads(raw)
    verify(record)
    if openmm.__version__ != record['openmm_version']:
        raise ValueError('readback runtime differs')
    output.mkdir(parents=True, exist_ok=False)
    (output/'executed-plan.json').write_bytes(raw)
    results, errors, pair_count = [], [], 0
    for ordinal, case in enumerate(record['cases']):
        source = Path(case['source'])
        ownership = json.loads((source/'H-ownership.json').read_bytes())
        rows = ownership['rows']
        path = json.loads(gzip.decompress(Path(case['path']).read_bytes()))
        ca, ha, outward = frame(rows, path['initial_nm'], record['target'])
        methyl_H = {i for group in ownership['graph']['groups'] for i in group['hydrogens']}
        system = openmm.XmlSerializer.deserialize((source/'H-system.xml').read_text())
        particles, exceptions, cutoff, dielectric = parameters(system)
        grouped, groups = grouped_system(system)
        integrator = openmm.VerletIntegrator(0.)
        context = openmm.Context(grouped, integrator, openmm.Platform.getPlatformByName('Reference'))
        nodes = [('initial', path['initial_nm'])]
        nodes.extend((f"callback-{cb['ordinal']}", [cb['positions_nm'][3*i:3*i+3] for i in range(len(rows))]) for cb in path['callbacks'])
        nodes.append(('final', path['final_nm']))
        published = []
        raw_path = output/f'cell-{ordinal+1:02d}-pairs.jsonl.gz'
        try:
            with raw_path.open('xb') as stream, gzip.GzipFile(fileobj=stream, mode='wb', mtime=0) as compressed:
                for label, points in nodes:
                    if points[ca] != path['initial_nm'][ca]:
                        raise ValueError('focal fixed heavy frame moved')
                    context.setPositions(points*nanometer)
                    components = []
                    nb_vector = None
                    for group in groups:
                        vector = list(context.getState(forces=True, groups=1 << group['index']).getForces().value_in_unit(kilojoules_per_mole/nanometer)[ha])
                        components.append(dict(type=group['type'], index=group['index'], force_kj_mol_nm=vector,
                            outward_projection_kj_mol_nm=projection(vector, outward)))
                        if group['type'] == 'NonbondedForce':
                            nb_vector = vector
                    if nb_vector is None:
                        raise ValueError('native nonbonded measurement unavailable')
                    terms = []
                    q1, s1, e1 = particles[ha]
                    for i, partner in enumerate(rows):
                        if i == ha:
                            continue
                        q2, s2, e2 = particles[i]
                        exception = exceptions.get(tuple(sorted((ha, i))))
                        charge, sigma, epsilon = exception if exception is not None else (q1*q2, (s1+s2)/2, math.sqrt(e1*e2))
                        vector = [a-b for a, b in zip(points[ha], points[i])]
                        lj, coulomb = pair_force(vector, charge, sigma, epsilon, cutoff, dielectric, exception is not None)
                        total = [x+y for x, y in zip(lj, coulomb)]
                        terms.append(dict(index=i, partner=partner['key'], element=partner['element'],
                            methyl_H=i in methyl_H, distance_angstrom=10*math.hypot(*vector), exception=exception is not None,
                            force_kj_mol_nm=total, outward_projection_kj_mol_nm=projection(total, outward),
                            LJ_outward_projection_kj_mol_nm=projection(lj, outward), Coulomb_outward_projection_kj_mol_nm=projection(coulomb, outward)))
                    summed = [sum(t['force_kj_mol_nm'][j] for t in terms) for j in range(3)]
                    if any(not math.isclose(a, b, rel_tol=1e-8, abs_tol=1e-4) for a, b in zip(summed, nb_vector)):
                        raise ValueError('all-pair 3D sum does not close to native NonbondedForce')
                    errors.extend(abs(a-b) for a, b in zip(summed, nb_vector))
                    categories = {}
                    for name, chosen in (
                            ('methyl_H', [t for t in terms if t['methyl_H']]),
                            ('other_H', [t for t in terms if t['element'] in {'H', 'D'} and not t['methyl_H']]),
                            ('heavy', [t for t in terms if t['element'] not in {'H', 'D'}]),
                            ('same_residue', [t for t in terms if t['partner'][:4] == record['target'][:4]]),
                            ('other_residue', [t for t in terms if t['partner'][:4] != record['target'][:4]])):
                        categories[name] = dict(count=len(chosen), outward_projection_kj_mol_nm=sum(t['outward_projection_kj_mol_nm'] for t in chosen))
                    h_vector = [x-y for x, y in zip(points[ha], points[ca])]
                    radius = math.hypot(*h_vector)
                    summary = dict(stage=label, force_components=components, categories=categories,
                        CA_H_angstrom=10*radius, outward_cosine=projection(h_vector, outward)/radius,
                        total_nonbonded_force_kj_mol_nm=summed,
                        largest_inward_pairs=sorted((t for t in terms if t['outward_projection_kj_mol_nm'] < 0),
                            key=lambda t: (t['outward_projection_kj_mol_nm'], t['index']))[:3])
                    compressed.write((json.dumps(dict(summary=summary, pairs=terms), sort_keys=True, allow_nan=False)+'\n').encode())
                    published.append(summary)
                    pair_count += len(terms)
        finally:
            del context
            del integrator
        results.append(dict(case=case, nodes=published, raw_artifact=raw_path.name, raw_sha256=sha(raw_path.read_bytes())))
    verify(record)
    summary = dict(contexts=len(results), nodes=sum(len(c['nodes']) for c in results), HA_pair_observations=pair_count,
        maximum_native_pair_closure_error_kj_mol_nm=max(errors, default=None))
    write(output/'evidence.json', dict(schema='dvbfixer.graph-methyl-ASP49-attribution.v1', cases=results,
        summary=summary, plan_sha256=sha(raw), policy=record['policy'], production_changed=False,
        minimization_executed=False, scientific_qualification=False, original_thr2_causality_resolved=False))
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('plan')
    p.add_argument('--h-plan', type=Path, required=True)
    p.add_argument('--h-root', type=Path, required=True)
    p.add_argument('--h-audit', type=Path, required=True)
    p.add_argument('--whole-audit', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p = sub.add_parser('run')
    p.add_argument('--plan', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.command == 'plan':
        plan(args.h_plan, args.h_root, args.h_audit, args.whole_audit, args.output)
    else:
        run(args.plan, args.output)
