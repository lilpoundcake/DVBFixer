"""Retrospective exact HA-touched nonperiodic nonbonded pair decomposition.

All 24 cells/1,248 bound path nodes. No Context, minimization, coordinate change,
force ablation or native scoring. Complete pair terms stay private, no FF tables.
"""
from __future__ import annotations

import argparse
import gzip
import json
import math
from pathlib import Path

from scripts.audit_ha_minimizer_path import outward_frame, projection, sphere_tangent
from scripts.audit_legacy_start_controls import load_atoms
from scripts.research_dunbrack_prior import sha, write
from scripts.research_ha_minimizer_path import verify

COULOMB = 138.93545764438198  # OpenMM MD units, checked against synthetic native force.


def pair_force(vector_nm, charge_product, sigma_nm, epsilon_kj_mol, cutoff_nm, dielectric, exception=False):
    if len(vector_nm) != 3 or not all(math.isfinite(x) for x in (*vector_nm, charge_product, sigma_nm, epsilon_kj_mol, cutoff_nm, dielectric)):
        raise ValueError('nonfinite or malformed pair')
    if sigma_nm < 0 or epsilon_kj_mol < 0 or cutoff_nm <= 0 or dielectric <= 0:
        raise ValueError('unsupported pair parameters')
    radius = math.hypot(*vector_nm)
    zero = [0., 0., 0.]
    if (not exception and radius >= cutoff_nm) or (charge_product == 0 and epsilon_kj_mol == 0):
        return zero[:], zero[:]
    if radius <= 1e-12:
        raise ValueError('undefined interacting coincident pair')
    sixth = (sigma_nm/radius)**6
    lj_scale = 24*epsilon_kj_mol*(2*sixth*sixth-sixth)/radius**2
    reaction = 0. if exception else (dielectric-1)/(2*dielectric+1)/cutoff_nm**3
    coulomb_scale = COULOMB*charge_product*(1/radius**3-2*reaction)
    return [lj_scale*x for x in vector_nm], [coulomb_scale*x for x in vector_nm]


def parameters(system):
    import openmm
    from openmm.unit import elementary_charge, kilojoules_per_mole, nanometer

    forces = [f for f in system.getForces() if isinstance(f, openmm.NonbondedForce)]
    if len(forces) != 1:
        raise ValueError('original nonbonded force not unique')
    force = forces[0]
    if force.getNonbondedMethod() != openmm.NonbondedForce.CutoffNonPeriodic or force.getUseSwitchingFunction() or force.getNumParticleParameterOffsets() or force.getNumExceptionParameterOffsets():
        raise ValueError('only unswitched nonperiodic reaction-field/no-offset pair force is supported')
    particles = []
    for i in range(force.getNumParticles()):
        q, sigma, epsilon = force.getParticleParameters(i)
        particles.append((q.value_in_unit(elementary_charge), sigma.value_in_unit(nanometer), epsilon.value_in_unit(kilojoules_per_mole)))
    exceptions = {}
    for i in range(force.getNumExceptions()):
        a, b, charge, sigma, epsilon = force.getExceptionParameters(i)
        key = tuple(sorted((a, b)))
        if key in exceptions:
            raise ValueError('duplicate nonbonded exception')
        exceptions[key] = (charge.value_in_unit(elementary_charge**2), sigma.value_in_unit(nanometer), epsilon.value_in_unit(kilojoules_per_mole))
    return particles, exceptions, force.getCutoffDistance().value_in_unit(nanometer), force.getReactionFieldDielectric()


def plan(path_plan_path, path_run, prior_path, output):
    import openmm

    parent_raw = path_plan_path.read_bytes()
    parent = json.loads(parent_raw)
    verify(parent)
    prior_raw = prior_path.read_bytes()
    prior = json.loads(prior_raw)
    run_raw = (path_run/'evidence.json').read_bytes()
    run = json.loads(run_raw)
    if prior['plan_sha256'] != sha(parent_raw) or prior['run_sha256'] != sha(run_raw) or prior['summary']['measured_nodes'] != 1248:
        raise ValueError('complete bound native-path force readback differs')
    if len(prior['cases']) != 24 or any(r['status'] != 'completed' for r in prior['cases']):
        raise ValueError('pair attribution requires all original path cells')
    inputs = {**parent['inputs'], str(path_plan_path): sha(parent_raw), str(prior_path): sha(prior_raw), str(path_run/'evidence.json'): sha(run_raw)}
    for case, row in zip(parent['cases'], run['cases'], strict=True):
        folder = path_run/case['id']/f"seed-{case['seed']}"/case['source_arm']
        for name, digest in row['artifacts'].items():
            path = folder/name
            if sha(path.read_bytes()) != digest:
                raise ValueError('path artifact differs before pair enrollment')
            inputs[str(path)] = digest
        system = openmm.XmlSerializer.deserialize((Path(case['source'])/'h-internal-system.xml').read_text())
        parameters(system)
    sources = {**parent['source_hashes']}
    for name in ('scripts/audit_ha_nonbonded_pairs.py', 'scripts/audit_ha_minimizer_path.py',
                 'scripts/audit_legacy_start_controls.py', 'scripts/research_dunbrack_prior.py'):
        sources[str(Path(name).resolve())] = sha(Path(name).read_bytes())
    write(output, dict(schema='dvbfixer.ha-nonbonded-pairs-plan.v1', inputs=inputs, source_hashes=sources,
        path_plan=str(path_plan_path), path_run=str(path_run), prior=str(prior_path), openmm_version=openmm.__version__,
        policy=dict(timing='retrospective to known nonbonded-group dominance, frozen before identifying individual pair contributors',
            enrollment='ALL24 cells/1248 original nodes; every HA-other-particle pair, including exceptions and beyond-cutoff zeros',
            physics='native Lorentz-Berthelot LJ and reaction-field Coulomb; exception pairs use their own parameters and unmodified Coulomb, without cutoff',
            checks='analytic full three-dimensional pair sum must match every bound unmodified NonbondedForce measurement',
            publication='every node/class total and exactly three most-negative outward-axis pair projections per node; complete per-pair records private',
            retention='retain all blockers/partial files, no new coordinates, no exclusion edits or physical/causal pair ablation'),
        production_changed=False, native_coordinates_read=False))


def run(plan_path, output):
    import openmm

    raw = plan_path.read_bytes()
    record = json.loads(raw)
    verify(record)
    if openmm.__version__ != record['openmm_version']:
        raise ValueError('OpenMM version differs')
    parent = json.loads(Path(record['path_plan']).read_bytes())
    prior = json.loads(Path(record['prior']).read_bytes())
    output.mkdir(parents=True, exist_ok=False)
    (output/'executed-plan.json').write_bytes(raw)
    results, errors = [], []
    total_pairs = 0
    for cell_index, (case, prior_case) in enumerate(zip(parent['cases'], prior['cases'], strict=True)):
        identity = (case['id'], case['seed'], case['source_arm'])
        if tuple(prior_case[k] for k in ('id', 'seed', 'source_arm')) != identity:
            raise ValueError('prior native force cohort omitted/reordered cells')
        source = Path(case['source'])
        folder = Path(record['path_run'])/case['id']/f"seed-{case['seed']}"/case['source_arm']
        path = json.loads(gzip.decompress((folder/'path.json.gz').read_bytes()))
        rows = load_atoms(source/'h-internal-before.json.gz')
        index = [i for i, r in enumerate(rows) if r['key'] == ['A', '16', '', 'GLN', 'HA']]
        if len(index) != 1:
            raise ValueError('HA pair identity differs')
        focal = index[0]
        system = openmm.XmlSerializer.deserialize((source/'h-internal-system.xml').read_text())
        particles, exceptions, cutoff, dielectric = parameters(system)
        if len(particles) != len(rows):
            raise ValueError('particle/atom inventory differs')
        nodes = [('initial', None, path['initial_nm'])]
        nodes.extend(('callback', cb['ordinal'], [cb['positions_nm'][3*i:3*i+3] for i in range(len(rows))]) for cb in path['callbacks'])
        nodes.append(('final', None, path['final_nm']))
        centre, outward = outward_frame(rows)
        published = []
        raw_path = output/f'cell-{cell_index+1:02d}-pairs.jsonl.gz'
        with raw_path.open('xb') as stream, gzip.GzipFile(fileobj=stream, mode='wb', mtime=0, compresslevel=6) as compressed:
            for (label, ordinal, points), bound in zip(nodes, prior_case['nodes'], strict=True):
                if (label, ordinal) != (bound['label'], bound['ordinal']):
                    raise ValueError('prior path node inventory differs')
                _, tangent = sphere_tangent(centre, [10*x for x in points[focal]], outward)
                terms = []
                for i, partner in enumerate(rows):
                    if i == focal:
                        continue
                    exception = exceptions.get(tuple(sorted((focal, i))))
                    q1, s1, e1 = particles[focal]
                    q2, s2, e2 = particles[i]
                    charge, sigma, epsilon = exception if exception is not None else (q1*q2, (s1+s2)/2, math.sqrt(e1*e2))
                    vector = [a-b for a, b in zip(points[focal], points[i])]
                    lj, coulomb = pair_force(vector, charge, sigma, epsilon, cutoff, dielectric, exception is not None)
                    force = [a+b for a, b in zip(lj, coulomb)]
                    terms.append(dict(partner_index=i, partner=partner['key'], element=partner['element'],
                        distance_angstrom=10*math.hypot(*vector), exception=exception is not None,
                        force_kj_mol_nm=force, outward_projection_kj_mol_nm=projection(force, outward),
                        LJ_outward_projection_kj_mol_nm=projection(lj, outward),
                        Coulomb_outward_projection_kj_mol_nm=projection(coulomb, outward),
                        sphere_tangent_projection_kj_mol_nm=projection(force, tangent)))
                summed = [sum(t['force_kj_mol_nm'][j] for t in terms) for j in range(3)]
                original = [p for p in bound['force_components'] if p['type'] == 'NonbondedForce']
                if len(original) != 1 or any(not math.isclose(a, b, rel_tol=1e-8, abs_tol=1e-4) for a, b in zip(summed, original[0]['HA_force_kj_mol_nm'])):
                    raise ValueError('complete HA pair force sum differs from unmodified native NonbondedForce')
                errors.extend(abs(a-b) for a, b in zip(summed, original[0]['HA_force_kj_mol_nm']))
                groups = {}
                for name, selected in (
                        ('H_partner', [t for t in terms if t['element'] in {'H', 'D'}]),
                        ('heavy_partner', [t for t in terms if t['element'] not in {'H', 'D'}]),
                        ('same_residue', [t for t in terms if t['partner'][:4] == rows[focal]['key'][:4]]),
                        ('other_residue', [t for t in terms if t['partner'][:4] != rows[focal]['key'][:4]])):
                    groups[name] = dict(count=len(selected), outward_projection_kj_mol_nm=sum(t['outward_projection_kj_mol_nm'] for t in selected))
                best = sorted((t for t in terms if t['outward_projection_kj_mol_nm'] < 0), key=lambda t: (t['outward_projection_kj_mol_nm'], t['partner_index']))[:3]
                summary = dict(label=label, ordinal=ordinal, pair_count=len(terms), categories=groups,
                    total_force_kj_mol_nm=summed, total_outward_projection_kj_mol_nm=projection(summed, outward),
                    total_sphere_tangent_projection_kj_mol_nm=projection(summed, tangent),
                    largest_inward_axis_pairs=best)
                compressed.write((json.dumps(dict(summary, terms=terms), sort_keys=True, allow_nan=False)+'\n').encode())
                published.append(summary)
                total_pairs += len(terms)
        results.append(dict(id=case['id'], seed=case['seed'], source_arm=case['source_arm'], status='evaluated',
            nodes=published, raw_artifact=raw_path.name, raw_sha256=sha(raw_path.read_bytes())))
    verify(record)
    summary = dict(enrolled=len(results), evaluated=len(results), nodes=sum(len(r['nodes']) for r in results),
        HA_pair_observations=total_pairs, maximum_force_component_closure_error_kj_mol_nm=max(errors, default=None))
    write(output/'evidence.json', dict(schema='dvbfixer.ha-nonbonded-pairs.v1', cases=results, summary=summary,
        plan_sha256=sha(raw), production_changed=False, native_coordinates_read=False,
        scientific_qualification=False, original_thr2_causality_resolved=False))
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('plan')
    p.add_argument('--path-plan', type=Path, required=True)
    p.add_argument('--path-run', type=Path, required=True)
    p.add_argument('--prior', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p = sub.add_parser('run')
    p.add_argument('--plan', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.command == 'plan':
        plan(args.path_plan, args.path_run, args.prior, args.output)
    else:
        run(args.plan, args.output)
