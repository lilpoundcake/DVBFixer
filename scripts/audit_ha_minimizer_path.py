"""Saved-only independent HA path/force readback at ALL native callbacks.

No minimization/integration/coordinate perturbation. Components are Cartesian
forces, not L-BFGS search directions, integrated work or pairwise causality.
"""
from __future__ import annotations

import argparse
import gzip
import json
import math
from collections import Counter
from pathlib import Path

from scripts.audit_ha_start_control import ha_geometry, inventory
from scripts.audit_legacy_start_controls import coordinate_delta, load_atoms
from scripts.audit_stereo_force_projection import grouped_system
from scripts.research_dunbrack_prior import sha, write
from scripts.research_ha_minimizer_path import verify


def outward_frame(rows):
    mapped = inventory(rows)
    centre = mapped[('A', '16', '', 'GLN', 'CA')]['xyz']
    unit_sum = [0., 0., 0.]
    for name in ('N', 'C', 'CB'):
        vector = [x-y for x, y in zip(mapped[('A', '16', '', 'GLN', name)]['xyz'], centre)]
        norm = math.hypot(*vector)
        if norm <= 1e-12:
            raise ValueError('undefined fixed heavy outward frame')
        unit_sum = [a+b/norm for a, b in zip(unit_sum, vector)]
    norm = math.hypot(*unit_sum)
    if norm <= 1e-12:
        raise ValueError('undefined fixed heavy outward sum')
    return centre, [-x/norm for x in unit_sum]


def sphere_tangent(centre, point, outward):
    vector = [x-y for x, y in zip(point, centre)]
    radius = math.hypot(*vector)
    if radius <= 1e-12 or not all(math.isfinite(x) for x in vector+list(outward)):
        raise ValueError('undefined HA sphere tangent')
    radial = [x/radius for x in vector]
    cosine = sum(x*y for x, y in zip(outward, radial))
    tangent = [x-cosine*y for x, y in zip(outward, radial)]
    norm = math.hypot(*tangent)
    return cosine, [x/norm for x in tangent] if norm > 1e-10 else None


def projection(force, direction):
    return sum(x*y for x, y in zip(force, direction)) if direction is not None else None


def force_readback(source_system, template, nodes):
    import openmm
    from openmm.unit import kilojoules_per_mole, nanometer

    system = openmm.XmlSerializer.deserialize(source_system)
    if system.getNumConstraints() or system.usesPeriodicBoundaryConditions() or any(system.isVirtualSite(i) for i in range(system.getNumParticles())):
        raise ValueError('path readback outside zero-constraint/nonperiodic scope')
    grouped, groups = grouped_system(system)
    indices = [i for i, row in enumerate(template) if row['key'] == ['A', '16', '', 'GLN', 'HA']]
    if len(indices) != 1 or grouped.getNumParticles() != len(template):
        raise ValueError('exact HA particle identity differs')
    index = indices[0]
    centre, outward = outward_frame(template)
    integrator = openmm.VerletIntegrator(0.)
    context = openmm.Context(grouped, integrator, openmm.Platform.getPlatformByName('Reference'))
    summaries = []
    try:
        for node in nodes:
            points = node['points_nm']
            rows = [dict(row, xyz=[10*x for x in p]) for row, p in zip(template, points, strict=True)]
            cosine, tangent = sphere_tangent(centre, rows[index]['xyz'], outward)
            context.setPositions(points*nanometer)
            parts = []
            for group in groups:
                state = context.getState(forces=True, energy=True, groups=1 << group['index'])
                force = list(state.getForces().value_in_unit(kilojoules_per_mole/nanometer)[index])
                parts.append(dict(index=group['index'], type=group['type'], name=group['name'],
                    energy_kj_mol=state.getPotentialEnergy().value_in_unit(kilojoules_per_mole),
                    HA_force_kj_mol_nm=force, outward_projection_kj_mol_nm=projection(force, outward),
                    sphere_tangent_projection_kj_mol_nm=projection(force, tangent)))
            combined = context.getState(forces=True, energy=True)
            force = list(combined.getForces().value_in_unit(kilojoules_per_mole/nanometer)[index])
            energy = combined.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
            summed = [sum(p['HA_force_kj_mol_nm'][j] for p in parts) for j in range(3)]
            if any(not math.isclose(a, b, rel_tol=1e-8, abs_tol=1e-5) for a, b in zip(summed, force)):
                raise ValueError('component HA forces do not sum to combined force')
            if not math.isclose(sum(p['energy_kj_mol'] for p in parts), energy, rel_tol=1e-9, abs_tol=1e-5):
                raise ValueError('component energies do not close')
            callback = node.get('callback')
            gradient_error = None
            if callback is not None:
                gradient = callback['gradient_kj_mol_nm'][3*index:3*index+3]
                if any(not math.isclose(a, -b, rel_tol=1e-7, abs_tol=1e-4) for a, b in zip(force, gradient)):
                    raise ValueError('native HA objective gradient differs from physical force')
                stats = callback['statistics']
                if stats['restraint energy'] != 0 or stats['max constraint error'] != 0:
                    raise ValueError('unexpected constraint contribution to objective')
                if not math.isclose(stats['system energy'], energy, rel_tol=1e-8, abs_tol=1e-4):
                    raise ValueError('native callback energy differs from saved-coordinate readback')
                gradient_error = max(abs(a+b) for a, b in zip(force, gradient))
            summaries.append(dict(label=node['label'], ordinal=node.get('ordinal'), iteration=node.get('iteration'),
                geometry=ha_geometry(rows), outward_cosine=cosine,
                signed_outward_distance_angstrom=projection([x-y for x, y in zip(rows[index]['xyz'], centre)], outward),
                sphere_tangent_defined=tangent is not None, force_components=parts,
                total_HA_force_kj_mol_nm=force, total_energy_kj_mol=energy,
                outward_projection_kj_mol_nm=projection(force, outward),
                sphere_tangent_projection_kj_mol_nm=projection(force, tangent),
                native_gradient_maximum_error_kj_mol_nm=gradient_error))
    finally:
        del context
        del integrator
    return summaries


def audit(plan_path, root, output):
    raw = plan_path.read_bytes()
    plan = json.loads(raw)
    verify(plan)
    run_raw = (root/'evidence.json').read_bytes()
    run = json.loads(run_raw)
    expected = [(c['id'], c['seed'], c['source_arm']) for c in plan['cases']]
    if run['plan_sha256'] != sha(raw) or (root/'executed-plan.json').read_bytes() != raw:
        raise ValueError('reported H-path enrollment differs')
    if [(c['id'], c['seed'], c['source_arm']) for c in run['cases']] != expected:
        raise ValueError('reported H path omitted/reordered cells')
    results = []
    for case, row in zip(plan['cases'], run['cases'], strict=True):
        folder = root/case['id']/f"seed-{case['seed']}"/case['source_arm']
        source = Path(case['source'])
        for name, digest in row['artifacts'].items():
            if sha((folder/name).read_bytes()) != digest:
                raise ValueError('retained H-path artifact changed')
        result = {k: row[k] for k in ('id', 'seed', 'source_arm', 'status', 'returncode')}
        if row['status'] != 'completed':
            results.append(result)
            continue
        path = json.loads(gzip.decompress((folder/'path.json.gz').read_bytes()))
        trace = json.loads((folder/'path-trace.json').read_bytes())
        if trace['case'] != case or trace['plan_sha256'] != sha(raw) or trace['reporter_early_stop'] or trace['integration']:
            raise ValueError('reported path contract differs')
        template = load_atoms(source/'h-internal-before.json.gz')
        for label in ('initial', 'final'):
            if len(path[label+'_nm']) != len(template) or any(row['xyz'] != [10*x for x in point]
                    for row, point in zip(path[label+'_atoms'], path[label+'_nm'], strict=True)):
                raise ValueError('native nm endpoint differs from audited Angstrom endpoint')
        for expected_path, observed in [('h-internal-before.json.gz', path['initial_atoms']), ('h-internal-after.json.gz', path['final_atoms'])]:
            delta = coordinate_delta(load_atoms(source/expected_path), observed)
            if any(not delta[k]['exact'] for k in ('heavy', 'hydrogen')):
                raise ValueError('reported path start/endpoint not reproduced exactly')
        callbacks = path['callbacks']
        if not callbacks or len(callbacks) != trace['callbacks'] or len(callbacks) > plan['policy']['maximum_callbacks_per_cell']:
            raise ValueError('native callback count outside declared inventory')
        if callbacks[-1]['positions_nm'] != [x for p in path['final_nm'] for x in p]:
            raise ValueError('final native callback differs from returned endpoint')
        nodes = [dict(label='initial', points_nm=path['initial_nm'])]
        for ordinal, cb in enumerate(callbacks):
            if cb['ordinal'] != ordinal or len(cb['positions_nm']) != 3*len(template) or len(cb['gradient_kj_mol_nm']) != 3*len(template):
                raise ValueError('native callback inventory differs')
            points = [cb['positions_nm'][3*i:3*i+3] for i in range(len(template))]
            if any(p != path['initial_nm'][i] for i, p in enumerate(points) if template[i]['element'] not in {'H', 'D'}):
                raise ValueError('massless heavy coordinate moved at a native callback')
            nodes.append(dict(label='callback', ordinal=ordinal, iteration=cb['iteration'], points_nm=points, callback=cb))
        nodes.append(dict(label='final', points_nm=path['final_nm']))
        source_system = (source/'h-internal-system.xml').read_text()
        if (folder/'h-internal-system.xml').read_text() != source_system:
            raise ValueError('reported H System differs from original')
        states = force_readback(source_system, template, nodes)
        result.update(nodes=states, callbacks=len(callbacks),
            first_negative_outward_callback=next((s['ordinal'] for s in states if s['label'] == 'callback' and s['outward_cosine'] < 0), None),
            first_all_acute_callback=next((s['ordinal'] for s in states if s['label'] == 'callback' and max(s['geometry']['ha_angles_degrees'].values()) < 90), None),
            original_start_endpoint_exact=True, all_heavy_callback_coordinates_exact=True)
        results.append(result)
    states = [s for r in results if r['status'] == 'completed' for s in r['nodes']]
    errors = [s['native_gradient_maximum_error_kj_mol_nm'] for s in states if s['native_gradient_maximum_error_kj_mol_nm'] is not None]
    summary = dict(enrolled=len(results), statuses=dict(Counter(r['status'] for r in results)),
        callbacks=sum(r.get('callbacks', 0) for r in results), measured_nodes=len(states),
        maximum_native_HA_gradient_error_kj_mol_nm=max(errors, default=None))
    write(output, dict(schema='dvbfixer.ha-minimizer-path-audit.v1', cases=results, summary=summary,
        plan_sha256=sha(raw), run_sha256=sha(run_raw), audit_script_sha256=sha(Path(__file__).read_bytes()),
        driver_sha256=plan['source_hashes'][str(Path('scripts/research_ha_minimizer_path.py').resolve())],
        openmm_version=plan['openmm_version'], protocol=plan['policy'],
        production_changed=False, native_coordinates_read=False, scientific_qualification=False,
        original_thr2_causality_resolved=False))
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    audit(args.plan, args.run, args.output)
