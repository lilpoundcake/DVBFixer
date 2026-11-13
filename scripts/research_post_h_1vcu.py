"""Private paired 1VCU public-API continuation, not a historical-System replay."""
from __future__ import annotations

import argparse
import functools
import inspect
import os
import random
import shutil
import signal
import subprocess
import time
from dataclasses import asdict
from pathlib import Path
from unittest.mock import patch

from scripts.audit_thr2_native_capture import read
from scripts.post_hydrogen_public_api import apply
from scripts.qualify_glycam_fallback import snapshot as topology_snapshot
from scripts.research_dunbrack_prior import sha, write
from scripts.research_ha_minimizer_path import reporter_type, verify
from scripts.research_legacy_start_controls import compressed, platform_properties
from scripts.research_post_hydrogen_guard import check as check_previous
from scripts.research_thr2_native_capture import update_owned

from dvbfixer.ffutils.dat import DatRecord

ROOT = Path(__file__).resolve().parents[1]


def patch_prepared(raw, saved, enabled):
    """Bind the retained post-H HA/neighbors, alter exactly one coordinate span."""
    moved = saved['moved_indices']
    if len(moved) != 1:
        raise ValueError('one retained 1VCU HA proposal required')
    h = moved[0]
    atom = saved['atoms'][h]
    component = atom['identity']['component']
    if (component['chain_id'], component['sequence_number'], component['insertion_code'],
            atom['residue_name'], atom['identity']['name']) != ('A', '0', '', 'SER', 'HA'):
        raise ValueError('wrong retained target identity')
    parents = [b if a == h else a for a, b in saved['bonds'] if h in (a, b)]
    if len(parents) != 1:
        raise ValueError('unique HA parent required')
    parent = parents[0]
    heavy = [b if a == parent else a for a, b in saved['bonds'] if parent in (a, b)
             and saved['atoms'][b if a == parent else a]['element'] not in ('H', 'D')]
    required = [h, parent]+heavy
    if len(heavy) != 3:
        raise ValueError('three heavy neighbors required')
    lines = raw.decode('utf-8').splitlines(keepends=True)
    found = {}
    for index, line in enumerate(lines):
        if line.startswith(('ATOM  ', 'HETATM')):
            key = line[21], line[22:26].strip(), line[26].strip(), line[17:20].strip(), line[12:16].strip()
            if key in found:
                raise ValueError('prepared PDB has ambiguous atom identity')
            found[key] = index
    for i in required:
        a = saved['atoms'][i]
        c = a['identity']['component']
        key = c['chain_id'], c['sequence_number'], c['insertion_code'], a['residue_name'], a['identity']['name']
        line = lines[found[key]]
        observed = [float(line[j:j+8]) for j in (30, 38, 46)]
        if any(abs(x-10*y) > 0.0005+1e-12 for x, y in zip(observed, saved['initial_nm'][i], strict=True)):
            raise ValueError('prepared target geometry is not the retained post-H handoff')
    target = found[('A', '0', '', 'SER', 'HA')]
    if enabled:
        xyz = ''.join(f'{10*x:8.3f}' for x in saved['proposed_nm'][h])
        if len(xyz) != 24:
            raise ValueError('proposal exceeds PDB coordinate fields')
        line = lines[target]
        lines[target] = line[:30]+xyz+line[54:]
    result = ''.join(lines).encode('utf-8')
    original = raw.decode('utf-8').splitlines(keepends=True)
    if any(line != original[i] for i, line in enumerate(lines) if i != target):
        raise ValueError('non-target PDB bytes changed')
    if lines[target][:30]+lines[target][54:] != original[target][:30]+original[target][54:]:
        raise ValueError('target non-coordinate PDB bytes changed')
    return result


def cache_hashes(folder):
    if not folder.exists():
        return {}
    result = {}
    for p in sorted(folder.rglob('*')):
        if p.is_symlink():
            raise ValueError('owned ligand cache may not contain symlinks')
        if p.is_file():
            result[str(p.relative_to(folder))] = sha(p.read_bytes())
    return result


def check(plan):
    import openmm.version

    verify(plan)
    check_previous(read(Path(plan['parent_plan'])))
    if (openmm.version.version, openmm.version.git_revision) != (plan['runtime']['version'], plan['runtime']['revision']):
        raise ValueError('OpenMM runtime changed')


def freeze(output):
    import openmm.version

    parent_path = ROOT/'docs/research/post-hydrogen-guard-static-plan.json'
    parent = read(parent_path)
    check_previous(parent)
    case = parent['cases'][14]
    if case['label'] != 'protein-ligands' or case['requested_ff'] != 'amber':
        raise ValueError('declared 1VCU fallback case missing')
    prepared_root = Path(case['source'])
    proposal = ROOT/'tmp/pdbfixer-research/post-hydrogen-guard-static-01/case-14.json.gz'
    original = ROOT/'tests/fixtures/1VCU.pdb'
    raw_source = original.read_text().splitlines()
    if (sum(line.startswith('MODEL ') for line in raw_source) > 1
            or any(line.startswith(('ATOM  ', 'HETATM')) and line[16] != ' ' for line in raw_source)):
        raise ValueError('private 1VCU experiment requires one altloc-unambiguous model')
    if any(line.startswith(('ATOM  ', 'HETATM')) and line[76:78].strip() in ('H', 'D') for line in raw_source):
        raise ValueError('this qualification requires a source with zero H/D observations')
    DatRecord.load(prepared_root/'prepared.dat')
    patch_prepared((prepared_root/'prepared.pdb').read_bytes(), read(proposal), False)
    inputs = {**parent['inputs'], **parent['source_hashes']}
    for path in (parent_path, original, proposal, prepared_root/'prepared.pdb', prepared_root/'prepared.dat',
                 ROOT/'docs/plans/post-hydrogen-public-api.md'):
        inputs[str(path)] = sha(path.read_bytes())
    sources = {str(p):sha(p.read_bytes()) for p in sorted((ROOT/'src/dvbfixer').rglob('*.py'))}
    for name in ('post_hydrogen_public_api.py', 'research_post_h_1vcu.py', 'audit_post_h_1vcu.py',
                 'qualify_glycam_fallback.py', 'research_thr2_native_capture.py', 'research_legacy_start_controls.py',
                 'research_ha_minimizer_path.py', 'research_dunbrack_prior.py'):
        p = ROOT/'scripts'/name
        sources[str(p)] = sha(p.read_bytes())
    for name in ('test_post_hydrogen_public_api.py', 'test_post_h_1vcu_handoff.py'):
        p = ROOT/'tests'/name
        sources[str(p)] = sha(p.read_bytes())
    import sys
    write(output, dict(schema='dvbfixer.private-post-h-1vcu-plan.v1', parent_plan=str(parent_path),
        inputs=inputs, source_hashes=sources, python=sys.executable,
        runtime=dict(version=openmm.version.version, revision=openmm.version.git_revision),
        prepared=str(prepared_root/'prepared.pdb'), dat=str(prepared_root/'prepared.dat'),
        original_source=str(original), proposal=str(proposal), source_observed_h=0,
        cases=[dict(arm=arm, enabled=enabled, h_seed=1) for arm, enabled in (('baseline', False), ('guarded', True))],
        cli=['--ff', 'amber', '--ph', '7', '--no-solvent', '--keep-water', '--atom-naming', 'standard'],
        policy=dict(per_case_seconds=1800, total_seconds=3600, maximum_h_calls=8,
            maximum_whole_calls=8, maximum_callbacks=4096, native_retries=0,
            energy_absolute_bound=1e-3, energy_relative_bound=2e-6,
            gradient_absolute_bound=1e-3, gradient_relative_bound=2e-5,
            penalty_energy_absolute_bound=1e-3, penalty_energy_relative_bound=2e-5,
            fp32_constraint_arithmetic_factor=32,
            label='one new paired public-minimize experiment, not an exact historical replay'),
        production_changed=False))


def table(value):
    if value is None:
        return None
    return dict(atoms=[asdict(a) for a in value.atoms], positions_nm=value.positions_nm, bonds=value.bonds)


def capture(context, topology, output, prefix):
    import openmm

    state = context.getState(positions=True, velocities=True, parameters=True, integratorParameters=True)
    for name, value in (('system', context.getSystem()), ('integrator', context.getIntegrator()), ('state', state)):
        with (output/f'{prefix}-{name}.xml').open('x') as stream:
            stream.write(openmm.XmlSerializer.serialize(value))
    with (output/f'{prefix}.chk').open('xb') as stream:
        stream.write(context.createCheckpoint())
    topology_snapshot(output, prefix+'-atoms', topology, state.getPositions())
    return dict(prefix=prefix, platform=platform_properties(context.getPlatform(), context))


def worker(plan_path, index, output):
    from openmm.app import ForceField, Modeller, Simulation
    from openmm.unit import nanometer

    from dvbfixer.minimize.pipeline import main

    plan = read(plan_path)
    check(plan)
    case = plan['cases'][index]
    (output/'input.pdb').write_bytes(patch_prepared(Path(plan['prepared']).read_bytes(), read(Path(plan['proposal'])), case['enabled']))
    (output/'input.dat').write_bytes(Path(plan['dat']).read_bytes())
    DatRecord.load(output/'input.dat')
    original_h, original_sim = Modeller.addHydrogens, Simulation.minimizeEnergy
    original_init, original_create = ForceField.__init__, ForceField.createSystem
    ff_origins, h_calls, whole_calls, creates = {}, [], [], []
    active = dict(h_call=None)

    def progress():
        update_owned(output/'capture.json', dict(plan_sha256=sha(plan_path.read_bytes()),
            case=case, h_calls=h_calls, whole_calls=whole_calls, creates=creates))

    @functools.wraps(original_init)
    def init(ff, *files, **kwargs):
        result = original_init(ff, *files, **kwargs)
        ff_origins[id(ff)] = (ff, [str(f) for f in files])
        return result

    @functools.wraps(original_create)
    def create(ff, topology, *args, **kwargs):
        row = dict(number=len(creates)+1, h_call=active['h_call'],
                   xmls=ff_origins.get(id(ff), (None, None))[1], atoms=topology.getNumAtoms(), returned=False)
        creates.append(row)
        try:
            system = original_create(ff, topology, *args, **kwargs)
            row.update(returned=True, particles=system.getNumParticles())
            return system
        except BaseException as error:
            row.update(error_type=type(error).__name__, error=str(error))
            raise
        finally:
            progress()

    @functools.wraps(original_h)
    def hydrogen(model, *args, **kwargs):
        if len(h_calls) >= plan['policy']['maximum_h_calls']:
            raise SystemExit('public H-call budget exceeded')
        bound = inspect.signature(original_h).bind(model, *args, **kwargs)
        bound.apply_defaults()
        ff = bound.arguments['forcefield']
        number = len(h_calls)+1
        row = dict(number=number, returned=False, forcefield_provided=ff is not None,
                   xmls=ff_origins.get(id(ff), (None, None))[1] if ff is not None else None)
        h_calls.append(row)
        previous_active = active['h_call']
        active['h_call'] = number
        before_rng = random.getstate()
        random.seed(case['h_seed'])
        initialized_rng = random.getstate()

        def journal(event):
            raw = dict(event)
            raw['before'], raw['raw'] = table(event['before']), table(event['raw'])
            if 'result' in raw:
                raw['result'] = asdict(raw['result']) if raw['result'] is not None else None
            compressed(output/f'h-{number}-journal.json.gz', raw)
            row.update(stage=event['stage'], admission_error=event.get('admission_error'))
            if event.get('result') is not None:
                row.update(moved=[asdict(a) for a in event['result'].moved], supported=len(event['result'].supported),
                           unsupported=len(event['result'].unsupported))
            row['api_observation_rounding_nm'] = event.get('api_observation_rounding_nm')
            progress()

        progress()
        try:
            outcome = apply(model, original_h, *args, protected_hydrogens=frozenset(),
                            enabled=case['enabled'], journal=journal, **kwargs)
            topology_snapshot(output, f'h-{number}-committed', model.topology, model.positions)
            row.update(returned=True, selected_variants=outcome.returned)
            return outcome.returned
        except BaseException as error:
            row.update(error_type=type(error).__name__, error=str(error))
            raise
        finally:
            compressed(output/f'h-{number}-rng.json.gz', dict(before=before_rng, initialized=initialized_rng, after=random.getstate()))
            random.setstate(before_rng)
            row['caller_rng_restored'] = random.getstate() == before_rng
            active['h_call'] = previous_active
            progress()

    @functools.wraps(original_sim)
    def simulate(simulation, *args, **kwargs):
        if len(whole_calls) >= plan['policy']['maximum_whole_calls']:
            raise SystemExit('whole-call budget exceeded')
        bound = inspect.signature(original_sim).bind(simulation, *args, **kwargs)
        bound.apply_defaults()
        if bound.arguments.get('reporter') is not None:
            raise SystemExit('unexpected existing native reporter')
        number = len(whole_calls)+1
        import openmm
        from openmm.unit import kilojoules_per_mole
        tolerance = bound.arguments['tolerance']
        row = dict(number=number, returned=False, max_iterations=bound.arguments['maxIterations'],
                   tolerance=float(tolerance.value_in_unit(kilojoules_per_mole/nanometer)))
        whole_calls.append(row)
        row['before'] = capture(simulation.context, simulation.topology, output, f'whole-{number}-before')
        initial = simulation.context.getState(positions=True).getPositions(asNumpy=True).value_in_unit(nanometer).tolist()
        reporter = reporter_type(openmm)(len(initial), plan['policy']['maximum_callbacks'])
        bound.arguments['reporter'] = reporter
        progress()
        try:
            result = original_sim(*bound.args, **bound.kwargs)
            row['returned'] = True
            return result
        except BaseException as error:
            row.update(error_type=type(error).__name__, error=str(error))
            raise
        finally:
            final = simulation.context.getState(positions=True).getPositions(asNumpy=True).value_in_unit(nanometer).tolist()
            compressed(output/f'whole-{number}-path.json.gz', dict(initial_nm=initial, final_nm=final,
                       callbacks=reporter.records, returned=row['returned']))
            row['after'] = capture(simulation.context, simulation.topology, output, f'whole-{number}-after')
            row['callbacks'] = len(reporter.records)
            progress()

    args = [str(output/'input.pdb'), '-o', str(output/'minimized.pdb'), '--dat', str(output/'input.dat'), *plan['cli']]
    write(output/'command.json', dict(args=args, case=case))
    with patch.object(ForceField, '__init__', init), patch.object(ForceField, 'createSystem', create), \
            patch.object(Modeller, 'addHydrogens', hydrogen), patch.object(Simulation, 'minimizeEnergy', simulate):
        main(args)
    progress()


def run(plan_path, root):
    plan = read(plan_path)
    check(plan)
    root.mkdir(parents=True, exist_ok=False)
    (root/'executed-plan.json').write_bytes(plan_path.read_bytes())
    started, rows = time.monotonic(), []
    for index, case in enumerate(plan['cases']):
        folder = root/f'case-{index:02d}'
        folder.mkdir()
        cache = folder/'ligand-cache'
        inherited = False
        if index:
            previous = root/'case-00/ligand-cache'
            cache_hashes(previous)
            if previous.exists():
                shutil.copytree(previous, cache)
                inherited = True
        if not cache.exists():
            cache.mkdir()
        result = dict(index=index, arm=case['arm'], status='not-executed-budget', returncode=None,
                      inherited_baseline_cache=inherited, cache_before=cache_hashes(cache))
        remaining = plan['policy']['total_seconds']-(time.monotonic()-started)
        if remaining > 0:
            env = dict(os.environ, PYTHONPATH=os.pathsep.join((str(ROOT/'src'), str(ROOT))),
                       PATH=str(Path(plan['python']).parent)+os.pathsep+os.environ.get('PATH', ''),
                       DVBFIXER_LIG_CACHE=str(cache))
            command = [plan['python'], str(Path(__file__).resolve()), 'worker', '--plan', str(plan_path),
                       '--index', str(index), '--output', str(folder)]
            with (folder/'run.log').open('xb') as log:
                process = subprocess.Popen(command, cwd=folder, env=env, stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
                try:
                    code = process.wait(timeout=min(remaining, plan['policy']['per_case_seconds']))
                    result.update(status='completed' if code == 0 else 'execution-error', returncode=code)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGKILL)
                    result.update(status='timeout', returncode=process.wait())
        result['cache_after'] = cache_hashes(cache)
        result['artifacts'] = {str(p.relative_to(folder)):sha(p.read_bytes()) for p in sorted(folder.rglob('*')) if p.is_file()}
        write(root/f'case-{index:02d}-result.json', result)
        rows.append(result)
        update_owned(root/'progress.json', dict(cases=rows))
    write(root/'evidence.json', dict(plan_sha256=sha(plan_path.read_bytes()), cases=rows))
    check(plan)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=('plan', 'run', 'worker'))
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--plan', type=Path)
    parser.add_argument('--index', type=int)
    args = parser.parse_args()
    if args.command == 'plan':
        freeze(args.output.resolve())
    elif args.command == 'worker':
        worker(args.plan.resolve(), args.index, args.output.resolve())
    else:
        run(args.plan.resolve(), args.output.resolve())
