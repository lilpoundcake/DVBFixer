"""Saved native H Systems: atom/bond identity, singular pair and force readback.

ALL six captures and BOTH before/returned states; no minimizer or force edits
beyond verified measurement group metadata. Never infer native termination enums.
"""
from __future__ import annotations

import argparse
import json
import math
from collections import deque
from pathlib import Path

from scripts.audit_stereo_force_projection import grouped_system
from scripts.audit_thr2_native_capture import read, verify_files
from scripts.research_dunbrack_prior import sha, write
from scripts.research_ser_hg_capture import check


def graph_distance(bonds, start, stop):
    adjacent = {}
    for a, b in bonds:
        adjacent.setdefault(a, []).append(b)
        adjacent.setdefault(b, []).append(a)
    queue, seen = deque([(start, 0)]), {start}
    while queue:
        atom, distance = queue.popleft()
        if atom == stop:
            return distance
        for other in adjacent.get(atom, []):
            if other not in seen:
                seen.add(other)
                queue.append((other, distance+1))
    return None


def audit(plan_path, root, output):
    import openmm
    from openmm import unit

    plan = read(plan_path)
    check(plan)
    execution = read(root/'evidence.json')
    if execution['plan_sha256'] != sha(plan_path.read_bytes()) or [(c['label'], c['seed']) for c in execution['cases']] != [(c['label'], c['seed']) for c in plan['cases']]:
        raise ValueError('complete frozen six-case execution differs')
    reports = []
    for result in execution['cases']:
        source = root/result['label']
        verify_files(source, result)
        capture = read(source/'capture.json')
        report = dict(label=result['label'], status=result['status'], native=[])
        for call in capture['native']:
            prefix = call['before']['prefix']
            rows = read(source/f'{prefix}-atoms.json.gz')['atoms']
            target = {r['key'][-1]:i for i, r in enumerate(rows) if r['key'][:3] == ['B', '126', '']}
            hg, og, oxt = [target[name] for name in ('HG', 'OG', 'OXT')]
            xml = (source/f'{prefix}-system.xml').read_text()
            system = openmm.XmlSerializer.deserialize(xml)
            bonds = read(source/f"native-{call['number']}-topology.json")['bonds']
            nonbonded = [f for f in system.getForces() if isinstance(f, openmm.NonbondedForce)]
            if len(nonbonded) != 1:
                raise ValueError('exactly one native nonbonded force required')
            force = nonbonded[0]
            pair = {name:force.getParticleParameters(index) for name, index in (('HG', hg), ('OXT', oxt))}
            parameters = {name:dict(charge_e=float(p[0].value_in_unit(unit.elementary_charge)),
                sigma_nm=float(p[1].value_in_unit(unit.nanometer)), epsilon_kj_mol=float(p[2].value_in_unit(unit.kilojoules_per_mole))) for name, p in pair.items()}
            exceptions = [i for i in range(force.getNumExceptions()) if set(force.getExceptionParameters(i)[:2]) == {hg, oxt}]
            bonded = []
            for f in system.getForces():
                if isinstance(f, openmm.HarmonicBondForce):
                    for i in range(f.getNumBonds()):
                        a, b, length, strength = f.getBondParameters(i)
                        if hg in (a, b):
                            bonded.append(dict(atoms=[rows[a]['key'], rows[b]['key']],
                                equilibrium_angstrom=length.value_in_unit(unit.angstrom),
                                strength_kj_mol_nm2=strength.value_in_unit(unit.kilojoules_per_mole/unit.nanometer**2)))
            mobile = [i for i in range(system.getNumParticles()) if system.getParticleMass(i).value_in_unit(unit.dalton) > 0]
            item = dict(number=call['number'], system_sha256=sha(xml.encode()), atoms=len(rows), mobile=len(mobile),
                hg_bonds=bonded, hg_oxt_graph_distance=graph_distance(bonds, hg, oxt),
                hg_oxt_nonbonded_exceptions=exceptions, nonbonded_method=int(force.getNonbondedMethod()),
                cutoff_nm=force.getCutoffDistance().value_in_unit(unit.nanometer), particle_parameters=parameters,
                masses_dalton={name:system.getParticleMass(i).value_in_unit(unit.dalton) for name, i in (('HG', hg), ('OG', og), ('OXT', oxt))}, states=[])
            product = parameters['HG']['charge_e']*parameters['OXT']['charge_e']
            mixed_epsilon = math.sqrt(parameters['HG']['epsilon_kj_mol']*parameters['OXT']['epsilon_kj_mol'])
            # The finite reaction-field correction cannot cancel an attractive
            # 1/r singularity; no pair-energy estimate is called a whole objective.
            item['pair_limiting_behavior'] = dict(charge_product_e2=product, mixed_epsilon_kj_mol=mixed_epsilon,
                attractive_coulomb_without_pair_lj=(product < 0 and mixed_epsilon == 0 and not exceptions),
                qualification='local pair limit as r approaches zero with other distances finite, not a blanket force-field defect claim')
            grouped, groups = grouped_system(system)
            initial_state = openmm.XmlSerializer.deserialize((source/f'{prefix}-state.xml').read_text())
            initial_xyz = initial_state.getPositions(asNumpy=True).value_in_unit(unit.nanometer).tolist()
            epsilon = call['tolerance']/max(1., math.sqrt(sum(x*x for p in initial_xyz for x in p)/len(rows)))
            for side in ('before', 'after'):
                if side not in call:
                    continue
                state = openmm.XmlSerializer.deserialize((source/f"{call[side]['prefix']}-state.xml").read_text())
                xyz = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer).tolist()
                direction = [xyz[oxt][j]-xyz[hg][j] for j in range(3)]
                distance = math.sqrt(sum(x*x for x in direction))
                if distance <= 0:
                    raise ValueError('exact collision: force projection undefined; retain failure')
                unit_vector = [x/distance for x in direction]
                for mode in ('captured', 'Reference'):
                    declared = call['before']['platform']
                    name = 'Reference' if mode == 'Reference' else declared['name']
                    properties = {} if mode == 'Reference' else declared['constructor_properties']
                    integrator = openmm.VerletIntegrator(0.)
                    context = openmm.Context(grouped, integrator, openmm.Platform.getPlatformByName(name), properties)
                    try:
                        context.setState(state)
                        terms = []
                        for group in groups:
                            measured = context.getState(energy=True, forces=True, groups=1 << group['index'])
                            vector = measured.getForces(asNumpy=True).value_in_unit(unit.kilojoules_per_mole/unit.nanometer)[hg]
                            terms.append(dict(type=group['type'], energy_kj_mol=measured.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole),
                                hg_force_toward_oxt_kj_mol_nm=sum(x*y for x, y in zip(vector, unit_vector, strict=True))))
                        measured = context.getState(energy=True, forces=True)
                        vectors = measured.getForces(asNumpy=True).value_in_unit(unit.kilojoules_per_mole/unit.nanometer)
                        projection = sum(x*y for x, y in zip(vectors[hg], unit_vector, strict=True))
                        energy = measured.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
                        force_norm = math.sqrt(sum(vectors[i][j]**2 for i in mobile for j in range(3)))
                        ratio = force_norm/max(1., math.sqrt(sum(x*x for p in xyz for x in p)))
                        if not all(math.isfinite(v) for v in (projection, energy, force_norm, ratio)):
                            raise ValueError('nonfinite physical readback; do not report a pass')
                        item['states'].append(dict(side=side, mode=mode, hg_oxt_angstrom=10*distance,
                            og_hg_angstrom=10*math.dist(xyz[og], xyz[hg]),
                            energy_kj_mol=energy, hg_force_toward_oxt_kj_mol_nm=float(projection),
                            mobile_force_l2_kj_mol_nm=force_norm, native_scaled_gradient_ratio=ratio,
                            native_epsilon=epsilon, satisfies_native_gradient_test=ratio <= epsilon, terms=terms,
                            force_group_projection_sum_error_kj_mol_nm=float(abs(projection-sum(t['hg_force_toward_oxt_kj_mol_nm'] for t in terms))),
                            force_group_energy_sum_error_kj_mol=abs(energy-sum(t['energy_kj_mol'] for t in terms))))
                    finally:
                        del context
            report['native'].append(item)
        reports.append(report)
    nodes = [n for r in reports for c in r['native'] for n in c['states']]
    calls = [c for r in reports for c in r['native']]
    summary = dict(captures=len(reports), native_calls=len(calls), readback_states=len(nodes),
        distinct_native_systems=len({c['system_sha256'] for c in calls}),
        pair_limit_cases=sum(c['pair_limiting_behavior']['attractive_coulomb_without_pair_lj'] for c in calls),
        returned_nodes_satisfying_native_gradient_test=sum(n['satisfies_native_gradient_test'] for n in nodes if n['side'] == 'after'))
    write(output, dict(schema='dvbfixer.ser-hg-forces.v1', plan_sha256=sha(plan_path.read_bytes()),
        execution_sha256=sha((root/'evidence.json').read_bytes()), audit_script_sha256=sha(Path(__file__).read_bytes()),
        summary=summary, cases=reports, additional_native_minimization=False, production_changed=False,
        limits='saved-state force/identity evidence for fresh native collisions, not original lost-run attribution or observed native return/error enums'))
    check(plan)
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('plan', 'root', 'output'):
        parser.add_argument('--'+name, type=Path, required=True)
    args = parser.parse_args()
    audit(args.plan.resolve(), args.root.resolve(), args.output.resolve())
