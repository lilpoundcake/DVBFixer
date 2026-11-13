"""Independent paired whole-System/ownership audit; no optimization or H draws."""
from __future__ import annotations

import argparse
import io
import json
import math
from pathlib import Path

from scripts.audit_thr2_native_capture import read, verify_files
from scripts.qualify_glycam_fallback import restore
from scripts.qualify_hydrogen_initialization_guard import atoms_from_rows
from scripts.research_dunbrack_prior import sha, write
from scripts.research_post_h_1vcu import check, patch_prepared
from scripts.research_thr2_hb_only import transport

from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef
from dvbfixer.ffutils.variants import VARIANT_TO_PARENT
from dvbfixer.model.diffusion.protein_atoms import CANONICAL_HEAVY_ATOMS


def digest(value):
    return sha(json.dumps(value, sort_keys=True, separators=(',', ':')).encode())


def identities(atoms):
    return [ExactAtomRef(ComponentInstanceRef(**a['identity']['component']), a['identity']['name']) for a in atoms]


def penalty(points, constraints, strength):
    """Analytic distance-restraint energy/gradient, independent of native code."""
    import numpy as np

    gradient = np.zeros_like(np.asarray(points, dtype=float))
    energy, max_error = 0., 0.
    allowance = np.zeros_like(gradient)
    for a, b, distance in constraints:
        delta = np.asarray(points[a])-np.asarray(points[b])
        radius = float(np.linalg.norm(delta))
        if radius <= 1e-15 or distance <= 0:
            raise ValueError('degenerate distance constraint')
        dr = radius-distance
        force_gradient = strength*dr*delta/radius
        gradient[a] += force_gradient
        gradient[b] -= force_gradient
        energy += 0.5*strength*dr*dr
        max_error = max(max_error, abs(dr)/distance)
        rounding = 32*np.finfo(np.float32).eps*abs(strength)*(radius+distance)
        allowance[a] += rounding
        allowance[b] += rounding
    return energy, gradient, max_error, allowance


def target_geometry(atoms, points):
    import numpy as np

    target = {a['identity']['name']:i for i, a in enumerate(atoms)
              if a['identity']['component']['chain_id'] == 'A'
              and a['identity']['component']['sequence_number'] == '0'
              and a['identity']['component']['insertion_code'] == ''
              and VARIANT_TO_PARENT.get(a['residue_name'], a['residue_name']) == 'SER'}
    if not all(n in target for n in ('N', 'CA', 'C', 'CB', 'HA')):
        return None
    p = np.asarray(points)
    heavy = [target[n] for n in ('N', 'C', 'CB')]
    hv = float(1000*np.linalg.det(p[heavy]-p[target['CA']]))
    fv = float(1000*np.linalg.det(p[heavy]-p[target['HA']]))
    return dict(heavy_ca_volume_angstrom3=hv, ha_four_volume_angstrom3=fv,
                sign_conflict=hv*fv < 0, ha_distance_angstrom=10*float(np.linalg.norm(p[target['HA']]-p[target['CA']])))


def audit_h(payload, committed, enabled):
    import numpy as np

    if payload['stage'] == 'api-error':
        return dict(stage='api-error', error_type=payload['error_type'], error=payload['error'])
    if payload['raw'] is None or payload.get('admission_error'):
        return dict(stage='guard-refusal', error=payload.get('admission_error'), guard_admitted=False)
    before, raw, proposal = payload['before'], payload['raw'], payload['result']
    old_ids, new_ids = identities(before['atoms']), identities(raw['atoms'])
    lookup = {key:i for i, key in enumerate(new_ids)}
    if len(lookup) != len(new_ids) or not set(old_ids) <= set(new_ids):
        raise ValueError('API deleted/merged input identity')
    new = [i for i, key in enumerate(new_ids) if key not in set(old_ids)]
    if any(raw['atoms'][i]['element'] != 'H' for i in new):
        raise ValueError('new heavy atom in public H API')
    for i, identity in enumerate(old_ids):
        a, b = before['atoms'][i], raw['atoms'][lookup[identity]]
        if (a['element'] != b['element'] or VARIANT_TO_PARENT.get(a['residue_name'], a['residue_name'])
                != VARIANT_TO_PARENT.get(b['residue_name'], b['residue_name'])):
            raise ValueError('existing atom chemistry changed')
    transport(before['positions_nm'], [raw['positions_nm'][lookup[k]] for k in old_ids],
              'Reference' if payload['reference_requested'] else 'captured')
    a, b = np.asarray(raw['positions_nm']), np.asarray(proposal['positions_nm'])
    moved = np.where(np.any(a != b, axis=1))[0].tolist()
    if not set(moved) <= set(new):
        raise ValueError('guard moved a pre-call observation')
    adjacent = [[] for _ in new_ids]
    for x, y in raw['bonds']:
        adjacent[x].append(y)
        adjacent[y].append(x)
    expected, supported = [], []
    for h in new:
        parent, = adjacent[h]
        center = raw['atoms'][parent]
        canonical = VARIANT_TO_PARENT.get(center['residue_name'], center['residue_name'])
        inventory = CANONICAL_HEAVY_ATOMS.get(canonical, ())
        heavy = [i for i in adjacent[parent] if raw['atoms'][i]['element'] not in ('H', 'D')]
        hs = [i for i in adjacent[parent] if raw['atoms'][i]['element'] in ('H', 'D')]
        if (center['element'] != 'C' or len(heavy) != 3 or hs != [h]
                or new_ids[parent].name not in inventory
                or any(new_ids[i].component != new_ids[parent].component for i in [h]+heavy)
                or any(new_ids[i].name not in inventory for i in heavy)):
            continue
        supported.append(h)
        vectors = a[heavy]-a[parent]
        four = a[heavy]-a[h]
        hv, fv = float(np.linalg.det(vectors)), float(np.linalg.det(four))
        lengths, radius = np.linalg.norm(vectors, axis=1), float(np.linalg.norm(a[h]-a[parent]))
        if np.any(lengths < 0.10) or np.any(lengths > 0.20) or not 0.05 <= radius <= 0.15:
            raise ValueError('accepted unsupported distance')
        if abs(hv/np.prod(lengths)) <= 1e-4 or abs(fv/np.prod(np.linalg.norm(four, axis=1))) <= 1e-4:
            raise ValueError('accepted degenerate supported center')
        if hv*fv < 0:
            expected.append(h)
            direction = np.sum(a[parent]-a[heavy], axis=0)
            wanted = a[parent]+radius*direction/np.linalg.norm(direction)
            if np.linalg.norm(b[h]-wanted) > 1e-14 or hv*float(np.linalg.det(b[heavy]-b[h])) <= 0:
                raise ValueError('independent post-H formula/sign mismatch')
            if abs(float(np.linalg.norm(b[h]-b[parent]))-radius) > 1e-14:
                raise ValueError('C-H radius changed')
            other = [i for i in range(len(a)) if i not in (h, parent)]
            if np.min(np.linalg.norm(b[other]-b[h], axis=1)) < 0.05:
                raise ValueError('gross proposal overlap')
    if moved != expected or [raw['atoms'][i]['identity'] for i in moved] != proposal['moved']:
        raise ValueError('independent selectivity mismatch')
    old_set = set(old_ids)
    old_graph = {tuple(sorted((old_ids[i], old_ids[j]))) for i, j in before['bonds']}
    after_graph = {tuple(sorted((new_ids[i], new_ids[j]))) for i, j in raw['bonds'] if new_ids[i] in old_set and new_ids[j] in old_set}
    if old_graph != after_graph:
        raise ValueError('API changed pre-call bond graph')
    old_types = {tuple(sorted((old_ids[i], old_ids[j]))):(kind, order)
                 for i, j, kind, order in payload['before_bond_metadata']}
    raw_types = {tuple(sorted((new_ids[i], new_ids[j]))):(kind, order)
                 for i, j, kind, order in payload['raw_bond_metadata'] if new_ids[i] in old_set and new_ids[j] in old_set}
    if old_types != raw_types:
        raise ValueError('API changed observed bond metadata')
    if committed is None:
        raise ValueError('successful public H call lacks committed snapshot')
    committed_ids = [x.identity for x in atoms_from_rows(committed['atoms'])]
    if committed_ids != new_ids or committed['bonds'] != raw['bonds']:
        raise ValueError('commit identity/graph mismatch')
    if committed['positions_nm'] != (proposal['positions_nm'] if enabled else raw['positions_nm']):
        raise ValueError('raw/candidate commit mismatch')
    target_id = ExactAtomRef(ComponentInstanceRef(1, 'A', '0'), 'HA')
    return dict(stage='guard-proposal', guard_admitted=True, new_h=len(new), supported=len(supported),
        moved=[raw['atoms'][i]['identity'] for i in moved], observation_coordinates_exact_after_guard=True,
        target_ha_present_before=target_id in old_set, target_ha_new_in_call=target_id in {new_ids[i] for i in new},
        target_before_guard=target_geometry(raw['atoms'], a), target_after_proposal=target_geometry(raw['atoms'], b),
        target_committed=target_geometry(raw['atoms'], committed['positions_nm']))


def output_summary(path):
    from openmm.app import PDBFile

    from dvbfixer.ffutils.geometry import find_d_residues

    text = path.read_text()
    keys = [(line[21], line[22:26].strip(), line[26].strip(), line[17:20].strip(), line[12:16].strip())
            for line in text.splitlines() if line.startswith(('ATOM  ', 'HETATM'))]
    if len(keys) != len(set(keys)):
        raise ValueError('duplicate final PDB atom identity')
    pdb = PDBFile(io.StringIO(text))
    names = sorted({r.name for r in pdb.topology.residues()})
    return dict(atoms=len(keys), negative_alpha=[list(x) for x in find_d_residues(pdb.topology, pdb.positions)],
                ligand_names_present=[n for n in ('DAN', 'EPE') if n in names])


def isolated_anchors(path):
    """Audit clone only; retain/report anchors, never alter native Systems."""
    import openmm

    system = openmm.XmlSerializer.deserialize(path.read_text())
    anchors = []
    for force in system.getForces():
        if not isinstance(force, openmm.CustomExternalForce):
            continue
        names = [force.getPerParticleParameterName(i) for i in range(force.getNumPerParticleParameters())]
        if names != ['k', 'x0', 'y0', 'z0']:
            continue
        for i in range(force.getNumParticles()):
            particle, params = force.getParticleParameters(i)
            anchors.append([particle, *params[1:]])
            force.setParticleParameters(i, particle, [params[0], 0., 0., 0.])
    return sha(openmm.XmlSerializer.serialize(system).encode()), anchors


def audit_whole(plan, folder, call):
    import numpy as np
    import openmm
    from openmm.unit import dalton, kilojoules_per_mole, nanometer

    from dvbfixer.ffutils.geometry import find_d_residues

    number, prefix = call['number'], call['before']['prefix']
    snap = read(folder/f'{prefix}-atoms.json.gz')
    topology, _ = restore(snap)
    atoms = [dict(identity=as_identity(a), residue_name=r['key'][3], element=r['element'])
             for a, r in zip(atoms_from_rows(snap['atoms']), snap['atoms'], strict=True)]
    system = openmm.XmlSerializer.deserialize((folder/f'{prefix}-system.xml').read_text())
    integrator = openmm.XmlSerializer.deserialize((folder/f'{prefix}-integrator.xml').read_text())
    state = openmm.XmlSerializer.deserialize((folder/f'{prefix}-state.xml').read_text())
    p = call['before']['platform']
    context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName(p['name']), p['constructor_properties'])
    context.setState(state)
    path = read(folder/f'whole-{number}-path.json.gz')
    mobile = [i for i in range(system.getNumParticles()) if system.getParticleMass(i).value_in_unit(dalton) > 0]
    constraints = []
    for i in range(system.getNumConstraints()):
        x, y, distance = system.getConstraintParameters(i)
        constraints.append((x, y, float(distance.value_in_unit(nanometer))))
    if any(system.isVirtualSite(i) for i in range(system.getNumParticles())):
        raise ValueError('virtual-site audit is outside this experiment')
    single = p['all_properties'].get('Precision') == 'single'
    entries = [('initial', path['initial_nm'], None)]
    for index, cb in enumerate(path['callbacks']):
        if cb['ordinal'] != index:
            raise ValueError('callback ordinal mismatch')
        entries.append((f'callback-{index}', np.asarray(cb['positions_nm']).reshape((-1, 3)).tolist(), cb))
    entries.append(('returned', path['final_nm'], None))
    nodes, max_gradient, max_energy, max_penalty, physics_passed = [], 0., 0., 0., True
    epsilon = call['tolerance']/max(1., math.sqrt(sum(x*x for p0 in path['initial_nm'] for x in p0)/len(path['initial_nm'])))
    try:
        for label, points, cb in entries:
            points = np.asarray(points)
            if points.shape != (system.getNumParticles(), 3) or not np.all(np.isfinite(points)):
                raise ValueError('malformed native coordinates')
            context.setPositions(points.tolist()*nanometer)
            measured = context.getState(energy=True, forces=True)
            energy = float(measured.getPotentialEnergy().value_in_unit(kilojoules_per_mole))
            physical = -measured.getForces(asNumpy=True).value_in_unit(kilojoules_per_mole/nanometer)
            if not math.isfinite(energy) or not np.all(np.isfinite(physical[mobile])):
                raise ValueError('nonfinite independent physical readback')
            constraint_error = penalty(points, constraints, 0.)[2]
            negative = find_d_residues(topology, points.tolist()*nanometer)
            node = dict(stage=label, target=target_geometry(atoms, points), energy_kj_mol=energy,
                        negative_alpha=[list(x) for x in negative], maximum_relative_constraint_error=constraint_error,
                        physical_mobile_gradient_l2=float(np.linalg.norm(physical[mobile])))
            if not constraints:
                scaled_gradient = float(np.linalg.norm(physical[mobile]))/max(1., float(np.linalg.norm(points)))
                node.update(native_scaled_gradient=scaled_gradient, native_epsilon=epsilon,
                            satisfies_unconstrained_native_gradient_test=scaled_gradient <= epsilon)
            if cb is not None:
                stats = cb['statistics']
                ep, gp, _, allowance = penalty(points, constraints, stats['restraint strength'])
                expected = physical+gp
                native = np.asarray(cb['gradient_kj_mol_nm']).reshape((-1, 3))
                error = np.abs(expected-native)
                bound = (plan['policy']['gradient_absolute_bound']+plan['policy']['gradient_relative_bound']
                         *np.maximum(1., np.maximum(np.abs(expected), np.abs(native))))
                if single:
                    bound += allowance
                energy_error = abs(energy-stats['system energy'])
                penalty_error = abs(ep-stats['restraint energy'])
                passed = (bool(np.all(error[mobile] <= bound[mobile]))
                    and energy_error <= plan['policy']['energy_absolute_bound']+plan['policy']['energy_relative_bound']*max(1., abs(energy))
                    and penalty_error <= plan['policy']['penalty_energy_absolute_bound']+plan['policy']['penalty_energy_relative_bound']*max(1., abs(stats['restraint energy'])))
                physics_passed &= passed
                max_gradient = max(max_gradient, float(np.max(error[mobile])))
                max_energy = max(max_energy, energy_error)
                max_penalty = max(max_penalty, penalty_error)
                node.update(gradient_maximum_error=float(np.max(error[mobile])), energy_readback_error=energy_error,
                    penalty_energy_readback_error=penalty_error, numerical_envelopes_passed=passed,
                    native_statistics=stats, augmented_mobile_gradient_l2=float(np.linalg.norm(expected[mobile])))
            nodes.append(node)
    finally:
        del context
    returned = openmm.XmlSerializer.deserialize((folder/f"{call['after']['prefix']}-state.xml").read_text()).getPositions(asNumpy=True).value_in_unit(nanometer).tolist()
    if returned != path['final_nm']:
        raise ValueError('returned State/path mismatch')
    last = entries[-2][1]
    delta = float(np.max(np.linalg.norm(np.asarray(last)-np.asarray(path['final_nm']), axis=1)))
    memberships = sorted({r['key'][3] for r in snap['atoms']})
    return dict(number=number, returned=path['returned'], callbacks=len(path['callbacks']), constraints=len(constraints),
        max_iterations=call['max_iterations'], last_reported_iteration=path['callbacks'][-1]['iteration'] if path['callbacks'] else None,
        final_relative_constraint_error=nodes[-1]['maximum_relative_constraint_error'],
        native_integrator_constraint_tolerance=integrator.getConstraintTolerance(),
        physics_envelopes_passed=physics_passed, maximum_mobile_gradient_error=max_gradient,
        maximum_system_energy_error=max_energy, maximum_penalty_energy_error=max_penalty,
        last_callback_to_return_maximum_nm=delta, constrained_stationarity_certified=False if constraints else None,
        ligand_names_present=[n for n in ('DAN', 'EPE') if n in memberships], nodes=nodes,
        system_sha256=sha((folder/f'{prefix}-system.xml').read_bytes()),
        topology_identity_sha256=sha(str([(r['key'], r['element']) for r in snap['atoms']]).encode()))


def as_identity(atom):
    from dataclasses import asdict
    return asdict(atom.identity)


def audit(plan_path, root, output):
    plan, execution = read(plan_path), read(root/'evidence.json')
    check(plan)
    if execution['plan_sha256'] != sha(plan_path.read_bytes()) or (root/'executed-plan.json').read_bytes() != plan_path.read_bytes():
        raise ValueError('execution binding mismatch')
    if len(execution['cases']) != len(plan['cases']):
        raise ValueError('incomplete two-arm accounting')
    reports = []
    for index, (case, result) in enumerate(zip(plan['cases'], execution['cases'], strict=True)):
        if result['arm'] != case['arm'] or result['index'] != index:
            raise ValueError('case order mismatch')
        folder = root/f'case-{index:02d}'
        verify_files(folder, result)
        report = {k:result[k] for k in ('index', 'arm', 'status', 'returncode', 'inherited_baseline_cache')}
        if not (folder/'capture.json').exists():
            report['missing_capture'] = True
            reports.append(report)
            continue
        capture = read(folder/'capture.json')
        if capture['plan_sha256'] != sha(plan_path.read_bytes()):
            raise ValueError('worker capture mismatch')
        expected = patch_prepared(Path(plan['prepared']).read_bytes(), read(Path(plan['proposal'])), case['enabled'])
        if (folder/'input.pdb').read_bytes() != expected or (folder/'input.dat').read_bytes() != Path(plan['dat']).read_bytes():
            raise ValueError('prepared/DatRecord handoff changed unexpectedly')
        calls = []
        raw_returns = []
        for h in capture['h_calls']:
            journal = folder/f"h-{h['number']}-journal.json.gz"
            if not journal.exists():
                calls.append(dict(number=h['number'], stage='no-journal', api_returned=h['returned'],
                    error_type=h.get('error_type'), error=h.get('error')))
                continue
            payload = read(journal)
            committed_path = folder/f"h-{h['number']}-committed.json.gz"
            committed = read(committed_path) if committed_path.exists() else None
            observed = audit_h(payload, committed, case['enabled'])
            observed.update(number=h['number'], forcefield_provided=h['forcefield_provided'], xmls=h['xmls'],
                            api_returned=h['returned'], caller_rng_restored=h['caller_rng_restored'])
            calls.append(observed)
            if h['returned'] and payload.get('raw') is not None:
                raw_returns.append(dict(number=h['number'], forcefield_provided=h['forcefield_provided'], xmls=h['xmls'],
                    atoms_sha256=digest(payload['raw']['atoms']), positions_sha256=digest(payload['raw']['positions_nm'])))
        whole = [audit_whole(plan, folder, c) for c in capture['whole_calls'] if c.get('after')]
        report.update(h_calls=calls, whole_calls=whole, output_exists=(folder/'minimized.pdb').is_file(),
            final_output=output_summary(folder/'minimized.pdb') if (folder/'minimized.pdb').is_file() else None,
            raw_api_returns=raw_returns,
            native_create_errors=[c for c in capture['creates'] if not c['returned']],
            whole_complex_present_in_all_phases=bool(whole) and all(len(c['ligand_names_present']) == 2 for c in whole))
        reports.append(report)
    system_comparison, api_comparison = [], []
    if len(reports) == 2:
        for a, b in zip(reports[0].get('whole_calls', []), reports[1].get('whole_calls', [])):
            number = a['number']
            left, la = isolated_anchors(root/f'case-00/whole-{number}-before-system.xml')
            right, ra = isolated_anchors(root/f'case-01/whole-{number}-before-system.xml')
            system_comparison.append(dict(number=a['number'], exact_system=a['system_sha256'] == b['system_sha256'],
                exact_atom_identity=a['topology_identity_sha256'] == b['topology_identity_sha256'],
                exact_parameters_with_only_restraint_anchors_isolated=left == right, exact_restraint_anchors=la == ra,
                maximum_anchor_component_difference_nm=max((abs(x-y) for p, q in zip(la, ra, strict=True)
                    for x, y in zip(p[1:], q[1:], strict=True)), default=0.) if len(la) == len(ra) else None))
        import numpy as np
        for a, b in zip(reports[0].get('raw_api_returns', []), reports[1].get('raw_api_returns', [])):
            left = read(root/f"case-00/h-{a['number']}-journal.json.gz")
            right = read(root/f"case-01/h-{b['number']}-journal.json.gz")
            lr = read(root/f"case-00/h-{a['number']}-rng.json.gz")
            rr = read(root/f"case-01/h-{b['number']}-rng.json.gz")
            same_atoms = a['atoms_sha256'] == b['atoms_sha256']
            delta = np.asarray(left['raw']['positions_nm'])-np.asarray(right['raw']['positions_nm']) if same_atoms else None
            api_comparison.append(dict(baseline_call=a['number'], guarded_call=b['number'], exact_atom_identity=same_atoms,
                same_ff_route=(a['forcefield_provided'], a['xmls']) == (b['forcefield_provided'], b['xmls']),
                same_initialized_rng=lr['initialized'] == rr['initialized'], same_final_rng=lr['after'] == rr['after'],
                exact_raw_coordinates=a['positions_sha256'] == b['positions_sha256'],
                maximum_raw_coordinate_difference_nm=float(np.max(np.linalg.norm(delta, axis=1))) if delta is not None else None))
    write(output, dict(schema='dvbfixer.private-post-h-1vcu-audit.v1', plan_sha256=sha(plan_path.read_bytes()),
        execution_sha256=sha((root/'evidence.json').read_bytes()), script_sha256=sha(Path(__file__).read_bytes()),
        cases=reports, paired_systems=system_comparison, paired_api_returns=api_comparison, production_changed=False,
        limits='One newly constructed paired 1VCU experiment, not a historical replay or accuracy/rate estimate. Constraint-aware numerical envelopes are not bit equality or global convergence.'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('plan', 'root', 'output'):
        parser.add_argument('--'+name, type=Path, required=True)
    args = parser.parse_args()
    audit(args.plan.resolve(), args.root.resolve(), args.output.resolve())
