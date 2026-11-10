"""Private factorial stage tracing on explicitly retained template failures.

No installed files or production functions are edited. All process-local hooks
are restored. Rebuild calls/seeds belong to the shared helper. Observational
stock runs must reproduce every saved untraced seed state before interpretation.
"""
from __future__ import annotations

import argparse
import gzip
import io
import json
from contextlib import ExitStack, redirect_stderr, redirect_stdout
from pathlib import Path
from unittest.mock import patch

import numpy as np
import openmm as mm
import pdbfixer.pdbfixer as pf
from openmm import unit
from scripts.branched_stereo_oracle import rdkit_labels
from scripts.branched_template_baseline import preflight
from scripts.research_rosetta_joint import sha
from scripts.research_rosetta_repack import write_json

from dvbfixer.ffutils.geometry import rebuild_missing_atoms_with_retry

PARENT = Path('tmp/pdbfixer-research/independent-comparison-01')
CASES = {
    '2CI2-VAL28-CB-CG1-CG2': 'previously confirmed D-CA endpoint',
    '2CI2-THR55-CB-OG1-CG2': 'previously confirmed D-CA endpoint',
    '1R69-THR27-CB-OG1-CG2': 'previously small-volume CA-unassigned endpoint, not confirmed D',
    '2CI2-THR22-CB-OG1-CG2': 'first whole-sidechain case by original enrollment order, not success-selected',
}


def atom_key(atom):
    residue = atom.residue
    return (residue.chain.id, str(residue.id), residue.insertionCode.strip(), atom.name)


def signed_volume(points):
    n, ca, c, cb = (np.asarray(points[name], dtype=float) for name in ('N', 'CA', 'C', 'CB'))
    return float(np.dot(np.cross(n-ca, c-ca), cb-ca))


def proper_overlay(points1, points2):
    if len(points1) < 3:
        raise ValueError('proper rotation control requires at least three anchors')
    p, q = np.asarray(points1, dtype=float), np.asarray(points2, dtype=float)
    u, _, vt = np.linalg.svd((q-q.mean(axis=0)).T @ (p-p.mean(axis=0)))
    correction = np.diag([1., 1., np.linalg.det(u @ vt)])
    return -q.mean(axis=0), (u @ correction @ vt).T, p.mean(axis=0)


def plan(output):
    raw = (PARENT/'evidence.json').read_bytes()
    parent = json.loads(raw)
    environment = preflight()
    if environment != parent['environment']:
        raise ValueError('scientific environment differs from retained untraced comparator')
    cases = []
    for identifier, scope in CASES.items():
        root = PARENT/identifier
        request = json.loads((root/'source-plan.json').read_bytes())['request']
        proposal = json.loads((root/'template-proposal.json').read_bytes())
        cases.append(dict(id=identifier, scope=scope, target=request['residue'], residue_name=request['residue_name'],
                          input_path=str(root/'input.pdb'), input_sha256=sha((root/'input.pdb').read_bytes()),
                          prior_proposal_path=str(root/'template-proposal.json'),
                          prior_proposal_sha256=sha((root/'template-proposal.json').read_bytes()),
                          prior_seed_count=len(proposal['snapshots'])))
    record = dict(schema='dvbfixer.template-causality-plan.v1', cases=cases, environment=environment,
                            parent_comparison_sha256=sha(raw), parent_comparison_path=str(PARENT/'evidence.json'),
                            arms=[dict(overlay=overlay, relaxation=relaxation) for overlay in ('stock', 'proper')
                                  for relaxation in ('normal', 'none')],
                            seeds='unchanged shared helper 1..5; max five attempts per arm, no independent extra rebuild',
                            interventions='proper determinant-corrected fit; relaxation-none bypasses BOTH minimization and escape dynamics, not a geometry remedy',
                            tracing='exact heavy/all construction boundaries, before/after minimization, accepted minimizer reporter iterates, before/after escape step calls, returned helper seed states',
                            admissibility='original source observations must survive; stock-normal must exactly reproduce all saved untraced seed atoms/positions; mismatch is an instrumentation blocker',
                            scope='retrospectively chosen endpoint failures plus identity-chosen first whole control; causal mechanisms only, not new accuracy qualification',
                            source_hashes={filename: sha(Path(filename).read_bytes()) for filename in (
                                str(Path(__file__).resolve().relative_to(Path.cwd())),
                                'scripts/branched_stereo_oracle.py', 'scripts/branched_template_baseline.py',
                                'src/dvbfixer/ffutils/geometry.py')},
                            installed_pdbfixer_sha256=sha(Path(pf.__file__).read_bytes()),
                            public_capability=False)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open('xb') as stream:
        stream.write((json.dumps(record, sort_keys=True, indent=2, allow_nan=False)+'\n').encode())


def run_arm(declaration, arm, destination):
    destination.mkdir(parents=True, exist_ok=False)
    source = Path(declaration['input_path']).read_bytes()
    (destination/'input.pdb').write_bytes(source)
    target = declaration['target']
    residue_key = target['chain_id'], target['sequence_number'], target['insertion_code']
    events, overlays, iterations, returns = [], [], [], []
    active_topology = None
    active_seed = None
    original_overlay = pf._overlayPoints
    original_minimize = mm.LocalEnergyMinimizer.minimize
    original_integrator = mm.LangevinIntegrator

    class Fixer(pf.PDBFixer):
        def _addAtomsToTopology(self, heavyAtomsOnly, omitUnknownMolecules):
            nonlocal active_topology
            result = super()._addAtomsToTopology(heavyAtomsOnly, omitUnknownMolecules)
            active_topology = result[0]
            capture('construction-heavy' if heavyAtomsOnly else 'construction-all', result[0], result[1])
            return result

    fixer = Fixer(pdbfile=io.StringIO(source.decode('ascii')))
    fixer.platform = mm.Platform.getPlatformByName('Reference')
    fixer.findMissingResidues()
    fixer.missingResidues = {}
    fixer.findMissingAtoms()
    observed_positions = np.asarray(fixer.positions.value_in_unit(unit.angstrom), dtype=float)
    observed = {atom_key(a): observed_positions[a.index] for a in fixer.topology.atoms()}

    def measured(topology, xyz):
        keys = [atom_key(a) for a in topology.atoms()]
        if len(keys) != len(set(keys)) or len(xyz) != len(keys):
            raise ValueError('ambiguous/mismatched stage atom inventory')
        by_key = dict(zip(keys, xyz, strict=True))
        if not observed.keys() <= by_key.keys():
            raise ValueError('stage lost original observations')
        movement = max(float(np.linalg.norm(by_key[k]-p)) for k, p in observed.items())
        points = {k[3]: list(map(float, p)) for k, p in by_key.items() if k[:3] == residue_key}
        return keys, points, movement

    def capture(stage, topology, positions, **extra):
        xyz = np.asarray(positions.value_in_unit(unit.angstrom), dtype=float)
        keys, points, movement = measured(topology, xyz)
        filename = f'event-{len(events):04d}.json.gz'
        payload = dict(atoms=[dict(key=k, element=a.element.symbol, xyz=list(map(float, xyz[a.index])))
                             for k, a in zip(keys, topology.atoms(), strict=True)])
        data = gzip.compress((json.dumps(payload, sort_keys=True, allow_nan=False)+'\n').encode(), mtime=0)
        (destination/filename).write_bytes(data)
        events.append(dict(stage=stage, seed=active_seed, coordinates=filename, coordinate_sha256=sha(data),
                           target_points_angstrom=points, ca_volume_angstrom3=signed_volume(points),
                           target_cip=rdkit_labels(points, declaration['residue_name']),
                           observed_max_movement_angstrom=movement, **extra))
        if movement > 1e-6:
            raise ValueError('stage moved original observed atoms')

    def overlay(points1, points2):
        # The low-anchor fallback is explicitly the ORIGINAL, never recursion
        # through the temporarily installed hook. These cases have >=3 anchors.
        if arm['overlay'] == 'proper' and len(points1) >= 3:
            result = proper_overlay(points1, points2)
        else:
            result = original_overlay(points1, points2)
        overlays.append(dict(seed=active_seed, anchors=len(points1), determinant=float(np.linalg.det(result[1])),
                             reference_points_nm=np.asarray(points1).tolist(), template_points_nm=np.asarray(points2).tolist(),
                             translate_template_nm=np.asarray(result[0]).tolist(), rotation=np.asarray(result[1]).tolist(),
                             translate_reference_nm=np.asarray(result[2]).tolist()))
        return result

    class Integrator(original_integrator):
        def step(self, steps):
            context = self._research_context
            capture('before-dynamics', active_topology, context.getState(getPositions=True).getPositions(), steps=steps)
            if arm['relaxation'] == 'normal':
                super().step(steps)
            capture('after-dynamics', active_topology, context.getState(getPositions=True).getPositions(), steps=steps,
                    bypassed=arm['relaxation'] == 'none')

    def minimize(context, *args, **kwargs):
        context.getIntegrator()._research_context = context
        call_index = sum(e['stage'] == 'before-minimize' for e in events)
        capture('before-minimize', active_topology, context.getState(getPositions=True).getPositions(), call_index=call_index)

        class Reporter(mm.MinimizationReporter):
            def report(self, iteration, x, grad, values):
                xyz = np.asarray(x, dtype=float).reshape(-1, 3)*10.
                _, points, movement = measured(active_topology, xyz)
                iterations.append(dict(seed=active_seed, call_index=call_index, iteration=iteration,
                                       target_points_angstrom=points, ca_volume_angstrom3=signed_volume(points),
                                       observed_max_movement_angstrom=movement))
                if movement > 1e-6:
                    raise ValueError('minimizer iterate moved an observed atom')
                return False

        if arm['relaxation'] == 'normal':
            original_minimize(context, *args, reporter=Reporter(), **kwargs)
        capture('after-minimize', active_topology, context.getState(getPositions=True).getPositions(), call_index=call_index,
                bypassed=arm['relaxation'] == 'none')

    original_add = fixer.addMissingAtoms

    def traced_add(seed=None):
        nonlocal active_seed
        active_seed = seed
        original_add(seed=seed)  # Only invoked through the shared helper below.
        capture('returned-seed', fixer.topology, fixer.positions)
        returns.append(events[-1])

    fixer.addMissingAtoms = traced_add
    try:
        with ExitStack() as stack:
            stack.enter_context(patch.object(pf, '_overlayPoints', overlay))
            stack.enter_context(patch.object(mm, 'LangevinIntegrator', Integrator))
            stack.enter_context(patch.object(mm.LocalEnergyMinimizer, 'minimize', minimize))
            with (destination/'helper.log').open('x') as log, redirect_stdout(log), redirect_stderr(log):
                rebuild_missing_atoms_with_retry(fixer, verbose=True, max_attempts=5, log_prefix='causal probe: ')
    finally:
        fixer.addMissingAtoms = original_add
    if (pf._overlayPoints is not original_overlay or mm.LocalEnergyMinimizer.minimize is not original_minimize
            or mm.LangevinIntegrator is not original_integrator):
        raise ValueError('process-local hooks were not restored')
    reproduced = None
    if arm == dict(overlay='stock', relaxation='normal'):
        prior = json.loads(Path(declaration['prior_proposal_path']).read_bytes())
        if len(returns) != len(prior['snapshots']):
            raise ValueError('instrumentation changed the untraced retry count')
        for event, snapshot in zip(returns, prior['snapshots'], strict=True):
            current = json.loads(gzip.decompress((destination/event['coordinates']).read_bytes()))['atoms']
            old = {(a['component']['chain_id'], a['component']['sequence_number'], a['component']['insertion_code'], a['name']): (element, xyz)
                   for a, element, xyz in snapshot['atoms']}
            if event['seed'] != snapshot['seed'] or {tuple(a['key']): (a['element'], a['xyz']) for a in current} != old:
                raise ValueError('instrumentation changed original seed coordinates/inventory')
        reproduced = True
    result = dict(schema='dvbfixer.template-causality-arm.v1', case_id=declaration['id'], arm=arm,
                  input_sha256=sha(source), events=events, overlays=overlays, minimizer_iterations=iterations,
                  stock_normal_exact_seed_reproduction=reproduced, hooks_restored=True, public_capability=False)
    write_json(destination/'evidence.json', result)
    return result


def run(plan_path, output):
    raw = plan_path.read_bytes()
    declaration = json.loads(raw)
    for filename, digest in declaration['source_hashes'].items():
        if sha(Path(filename).read_bytes()) != digest:
            raise ValueError('frozen experiment implementation changed')
    if preflight() != declaration['environment'] or sha(Path(pf.__file__).read_bytes()) != declaration['installed_pdbfixer_sha256']:
        raise ValueError('scientific environment changed')
    if sha(Path(declaration['parent_comparison_path']).read_bytes()) != declaration['parent_comparison_sha256']:
        raise ValueError('original comparison changed')
    output.mkdir(parents=True, exist_ok=False)
    (output/'executed-plan.json').write_bytes(raw)
    results = {}
    for case in declaration['cases']:
        if sha(Path(case['input_path']).read_bytes()) != case['input_sha256'] or sha(Path(case['prior_proposal_path']).read_bytes()) != case['prior_proposal_sha256']:
            raise ValueError('original causal-control input/proposal changed')
        for arm in declaration['arms']:
            name = case['id']+'-'+arm['overlay']+'-'+arm['relaxation']
            try:
                result = run_arm(case, arm, output/name)
            except Exception as error:
                write_json(output/'execution-error.json', dict(case=name, type=type(error).__name__, error=str(error)))
                raise
            results[name] = dict(case_id=case['id'], arm=arm, evidence_sha256=sha((output/name/'evidence.json').read_bytes()),
                                 seed_count=sum(e['stage'] == 'returned-seed' for e in result['events']))
    write_json(output/'evidence.json', dict(schema='dvbfixer.template-causality-run.v1', plan_sha256=sha(raw),
                                           results=results, installed_source_unchanged=sha(Path(pf.__file__).read_bytes()) == declaration['installed_pdbfixer_sha256']))
    print(json.dumps(dict(cases=len(declaration['cases']), arms=len(results), returned_seeds=sum(r['seed_count'] for r in results.values())), indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('plan')
    p.add_argument('--output', type=Path, required=True)
    p = sub.add_parser('run')
    p.add_argument('--plan', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.command == 'plan':
        plan(args.output)
    else:
        run(args.plan, args.output)
