"""Frozen real-input preparation cohort; matched postflight replay is separate.

No production patch, new force, minimizer option or native reporter. Each case
runs once; baseline/candidate geometry consume the same saved H coordinates.
"""
from __future__ import annotations

import argparse
import functools
import inspect
import json
import os
import random
import signal
import subprocess
import sys
import time
from pathlib import Path
from unittest.mock import patch

from scripts.research_dunbrack_prior import sha, write
from scripts.research_ha_minimizer_path import verify
from scripts.research_legacy_start_controls import compressed, snapshot_atoms
from scripts.research_thr2_native_capture import capture_context, update_owned

ROOT = Path(__file__).resolve().parents[1]
BASELINE = ROOT/'tmp/hydrogen-postflight-contract-01/baseline_geometry.py'
SOURCES = (
    ('ubiquitin-water', 'protein_context/1UBQ.pdb', 'amber', 7.),
    ('antibody-disulfide', 'regressions/1DQJ.pdb', 'amber', 7.),
    ('acidic-context-ph3', 'c_glh/8cde_t_u.pdb', 'amber', 3.),
    ('pdb-glycoprotein', '3ry6.pdb', 'amber+glycam', 7.),
    ('protein-ligands', '1VCU.pdb', 'amber', 7.),
    ('glycam-glycoprotein', 'hinge_CH3_glycosylated.pdb', 'amber+glycam', 7.),
)


def source_inventory(path):
    atoms = [line for line in path.read_text().splitlines() if line.startswith(('ATOM  ', 'HETATM'))]
    return dict(atoms=len(atoms), alternate_location_records=sum(line[16] != ' ' for line in atoms),
        hydrogens=sum(line[76:78].strip() in ('H', 'D') for line in atoms),
        chains=sorted({line[21] for line in atoms}),
        residue_names=sorted({line[17:21].strip() for line in atoms}))


def snapshot(output, prefix, topology, positions):
    from openmm.unit import nanometer

    compressed(output/f'{prefix}.json.gz', dict(atoms=snapshot_atoms(topology, positions),
        positions_nm=[[float(x) for x in p] for p in positions.value_in_unit(nanometer)],
        bonds=[[a.index, b.index] for a, b in topology.bonds()]))


def plan(output):
    import openmm.app.forcefield
    import openmm.app.modeller
    import openmm.version

    manifest = ROOT/'tests/fixtures/MANIFEST.sha256'
    checksums = {}
    for line in manifest.read_text().splitlines():
        if line.strip():
            digest, name = line.split(maxsplit=1)
            checksums[name.lstrip('*')] = digest
    inputs = {str(manifest):sha(manifest.read_bytes()), str(ROOT/'tests/fixtures/README.md'):sha((ROOT/'tests/fixtures/README.md').read_bytes()),
        str(BASELINE):sha(BASELINE.read_bytes())}
    protocol = ROOT/'docs/plans/oh-native-qualification.md'
    inputs[str(protocol)] = sha(protocol.read_bytes())
    if inputs[str(BASELINE)] != '8004da55b3a31bd5b4ad92c5ab3981ef0715af5f13ab09d85f45a13d8df17a1c':
        raise ValueError('pre-implementation baseline bytes differ')
    cases = []
    for label, relative, amber, ph in SOURCES:
        source = ROOT/'tests/fixtures'/relative
        digest = sha(source.read_bytes())
        if checksums.get(str(source.relative_to(ROOT))) != digest:
            raise ValueError(f'tracked input manifest differs: {source}')
        inputs[str(source)] = digest
        for request in (amber, 'charmm'):
            cases.append(dict(label=label, source=str(source), source_sha256=digest,
                requested_ff=request, ph=ph, seed=1, source_inventory=source_inventory(source)))
    for obj in (openmm.app.forcefield, openmm.app.modeller):
        p = Path(obj.__file__).resolve()
        inputs[str(p)] = sha(p.read_bytes())
    for p in sorted((Path(openmm.app.forcefield.__file__).parent/'data').rglob('*.xml')):
        inputs[str(p.resolve())] = sha(p.read_bytes())
    sources = {str(p):sha(p.read_bytes()) for p in sorted((ROOT/'src/dvbfixer').rglob('*.py'))}
    test = ROOT/'tests/test_qualify_oh_native.py'
    sources[str(test)] = sha(test.read_bytes())
    for name in ('qualify_oh_native.py', 'audit_oh_native.py', 'qualify_hydrogen_postflight.py',
                 'research_dunbrack_prior.py', 'research_ha_minimizer_path.py',
                 'research_legacy_start_controls.py', 'research_thr2_native_capture.py'):
        p = ROOT/'scripts'/name
        sources[str(p)] = sha(p.read_bytes())
    write(output, dict(schema='dvbfixer.oh-native-plan.v1', python=sys.executable,
        inputs=inputs, source_hashes=sources, baseline=str(BASELINE), cases=cases,
        runtime=dict(full_version=openmm.version.version, revision=openmm.version.git_revision),
        protocol=dict(cases=12, unique_structural_sources=6, per_case_seconds=1200, total_seconds=7200,
            maximum_h_calls=8, maximum_native_h_calls=8,
            enrollment='ALL six named real tracked inputs x two requested FF routes, H seed1, exactly one prep each; no outcome-dependent replacement or retry',
            cli='prepare copied input --ff REQUEST --ph P --keep-water --atom-naming standard; other legacy defaults unchanged',
            observation='native H boundaries and actual production repair before/after, exact nm positions and RNG state; no native reporter or energy/force query',
            comparison='offline old/new repair on identical saved pre-repair nm coordinates; require exact candidate replay of the actual in-process result',
            ff_coverage='record actual XMLs and forcefield-free fallback, not just requested FF; CHARMM naming via GLYCAM bridge is NOT native CHARMM glycan sampling',
            chemical_coverage='measure actual OH residue/component classes; late RDKit/OpenBabel ligand H is outside matched Modeller postflight and must not be claimed covered',
            validation='same-residue OH contract, immutable heavy/valid H/bonds/identity, independent short-OH and cross-residue contacts as warnings, DatRecord load and actual output status',
            readiness='any unexpected matched replay, valid-H/heavy/identity movement or unjustified new refusal blocks qualification; every refusal and upstream blocker remains visible',
            limits='six source structures and one H draw label per route, no rate/accuracy claim, independent native repeats, minimizer qualification or automatic merge'),
        production_changed=False, main_merge_authorized=False))


def check(record):
    import openmm.version

    verify(record)
    if openmm.version.git_revision != record['runtime']['revision'] or openmm.version.version != record['runtime']['full_version']:
        raise ValueError('native runtime provenance changed')


def worker(plan_path, index, output):
    import openmm
    from openmm.app import ForceField, Modeller
    from openmm.unit import kilojoules_per_mole, nanometer

    import dvbfixer.ffutils
    from dvbfixer.ffutils import geometry
    from dvbfixer.prepare import main

    record = json.loads(plan_path.read_bytes())
    check(record)
    case = record['cases'][index]
    old_h, old_create, old_min, old_repair = Modeller.addHydrogens, ForceField.createSystem, openmm.LocalEnergyMinimizer.minimize, geometry.repair_misplaced_hydrogens
    old_builder = dvbfixer.ffutils.build_glycam_system
    origins, registry, h_calls, native, repairs = {}, {}, [], [], []
    active = dict(hydrogen=False, h_call=None)

    def progress():
        update_owned(output/'capture.json', dict(case=case, h_calls=h_calls, native=native, repairs=repairs,
            plan_sha256=sha(plan_path.read_bytes())))

    @functools.wraps(old_builder)
    def builder(xmls, *args, **kwargs):
        ff = old_builder(xmls, *args, **kwargs)
        origins[id(ff)] = [str(p) for p in xmls]
        return ff

    @functools.wraps(old_create)
    def create(ff, topology, *args, **kwargs):
        system = old_create(ff, topology, *args, **kwargs)
        registry[int(system.this)] = (system, topology, origins.get(id(ff)))
        return system

    @functools.wraps(old_h)
    def hydrogen(model, *args, **kwargs):
        if len(h_calls) >= record['protocol']['maximum_h_calls']:
            raise SystemExit('research H-call bound exceeded; no automatic retry')
        bound = inspect.signature(old_h).bind(model, *args, **kwargs)
        bound.apply_defaults()
        ff = bound.arguments['forcefield']
        number = len(h_calls)+1
        row = dict(number=number, returned=False, forcefield_provided=ff is not None,
            forcefield_xmls=origins.get(id(ff)) if ff is not None else None,
            variants=bound.arguments['variants'], before_rng=random.getstate())
        snapshot(output, f'h-{number}-before', model.topology, model.positions)
        random.seed(case['seed'])
        row['initialized_rng'] = random.getstate()
        h_calls.append(row)
        active.update(hydrogen=True, h_call=number)
        progress()
        try:
            result = old_h(model, *args, **kwargs)
            row['returned'] = True
            return result
        except BaseException as error:
            row.update(error_type=type(error).__name__, error=str(error))
            raise
        finally:
            row['after_rng'] = random.getstate()
            random.setstate(row['before_rng'])
            row['restored_rng'] = random.getstate()
            active.update(hydrogen=False, h_call=None)
            snapshot(output, f'h-{number}-after', model.topology, model.positions)
            progress()

    @functools.wraps(old_min)
    def minimize(context, *args, **kwargs):
        if not active['hydrogen']:
            return old_min(context, *args, **kwargs)
        if len(native) >= record['protocol']['maximum_native_h_calls']:
            raise SystemExit('research native H-call bound exceeded; no automatic retry')
        bound = inspect.signature(old_min).bind(context, *args, **kwargs)
        bound.apply_defaults()
        if bound.arguments['reporter'] is not None:
            raise SystemExit('unexpected native reporter')
        system = context.getSystem()
        if int(system.this) in registry:
            topology = registry[int(system.this)][1]
        else:
            frame = inspect.currentframe()
            try:
                topology = frame.f_back.f_locals['newTopology']
            finally:
                del frame
        if topology.getNumAtoms() != system.getNumParticles():
            raise SystemExit('native atom/System inventory differs')
        number = len(native)+1
        tolerance = bound.arguments['tolerance']
        if hasattr(tolerance, 'value_in_unit'):
            tolerance = tolerance.value_in_unit(kilojoules_per_mole/nanometer)
        row = dict(number=number, h_call=active['h_call'], tolerance=float(tolerance),
            maximum_iterations=bound.arguments['maxIterations'], returned=False)
        native.append(row)
        write(output/f'native-{number}-bonds.json', dict(bonds=[[a.index, b.index] for a, b in topology.bonds()]))
        row['before'] = capture_context(context, output, f'native-{number}-before', topology)
        progress()
        try:
            result = old_min(context, *args, **kwargs)
            row['returned'] = True
            return result
        finally:
            row['after'] = capture_context(context, output, f'native-{number}-after', topology)
            progress()

    @functools.wraps(old_repair)
    def repair(topology, positions, *args, **kwargs):
        number = len(repairs)+1
        snapshot(output, f'repair-{number}-before', topology, positions)
        row = dict(number=number, returned=False)
        repairs.append(row)
        progress()
        started = time.perf_counter()
        try:
            result = old_repair(topology, positions, *args, **kwargs)
            row.update(returned=True, repairs=result)
            return result
        except BaseException as error:
            row.update(error_type=type(error).__name__, error=str(error))
            raise
        finally:
            row['repair_call_seconds'] = time.perf_counter()-started
            snapshot(output, f'repair-{number}-after', topology, positions)
            progress()

    args = [str(output/'input.pdb'), '-o', str(output/'prepared.pdb'), '--ff', case['requested_ff'],
        '--ph', str(case['ph']), '--keep-water', '--atom-naming', 'standard']
    write(output/'command.json', dict(args=args, seed=case['seed']))
    with patch.object(dvbfixer.ffutils, 'build_glycam_system', builder), patch.object(ForceField, 'createSystem', create), \
            patch.object(Modeller, 'addHydrogens', hydrogen), patch.object(openmm.LocalEnergyMinimizer, 'minimize', minimize), \
            patch.object(geometry, 'repair_misplaced_hydrogens', repair):
        main(args)
    check(record)
    progress()


def run(plan_path, output):
    record = json.loads(plan_path.read_bytes())
    check(record)
    output.mkdir(parents=True, exist_ok=False)
    (output/'executed-plan.json').write_bytes(plan_path.read_bytes())
    env = dict(os.environ, PATH=str(Path(record['python']).parent)+os.pathsep+os.environ.get('PATH', ''),
        PYTHONPATH=os.pathsep.join((str(ROOT/'src'), str(ROOT))))
    started, results = time.monotonic(), []
    for index, case in enumerate(record['cases']):
        folder = output/f'case-{index:02d}'
        folder.mkdir()
        (folder/'input.pdb').write_bytes(Path(case['source']).read_bytes())
        row = dict(index=index, label=case['label'], requested_ff=case['requested_ff'], status='budget_not_executed', returncode=None)
        remaining = record['protocol']['total_seconds']-(time.monotonic()-started)
        if remaining > 0:
            command = [record['python'], str(Path(__file__).resolve()), 'worker', '--plan', str(plan_path), '--index', str(index), '--output', str(folder)]
            with (folder/'prepare.log').open('xb') as stream:
                try:
                    process = subprocess.Popen(command, cwd=ROOT, env=env, stdout=stream,
                        stderr=subprocess.STDOUT, start_new_session=True)
                    code = process.wait(timeout=min(remaining, record['protocol']['per_case_seconds']))
                    row.update(status='completed' if code == 0 else 'execution_blocker', returncode=code)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGKILL)
                    process.wait()
                    row['status'] = 'timeout'
        row['artifacts'] = {str(p.relative_to(folder)):sha(p.read_bytes()) for p in sorted(folder.rglob('*')) if p.is_file()}
        write(folder/'result.json', row)
        results.append(row)
        update_owned(output/'progress.json', dict(plan_sha256=sha(plan_path.read_bytes()), cases=results))
        print(json.dumps({k:v for k, v in row.items() if k != 'artifacts'}), flush=True)
    check(record)
    write(output/'evidence.json', dict(plan_sha256=sha(plan_path.read_bytes()), cases=results, production_changed=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('plan')
    p.add_argument('--output', type=Path, required=True)
    for name in ('run', 'worker'):
        p = sub.add_parser(name)
        p.add_argument('--plan', type=Path, required=True)
        p.add_argument('--output', type=Path, required=True)
        if name == 'worker':
            p.add_argument('--index', type=int, required=True)
    args = parser.parse_args()
    if args.command == 'plan':
        plan(args.output.resolve())
    elif args.command == 'run':
        run(args.plan.resolve(), args.output.resolve())
    else:
        worker(args.plan.resolve(), args.index, args.output.resolve())
