"""Changed-candidate cohort with explicit topology occurrences and recovery trace."""
from __future__ import annotations

import argparse
import functools
import json
import os
import signal
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path
from unittest.mock import patch

from scripts import audit_oh_native, qualify_oh_native
from scripts.research_dunbrack_prior import sha, write
from scripts.research_ha_minimizer_path import verify
from scripts.research_legacy_start_controls import compressed, platform_properties, snapshot_atoms
from scripts.research_thr2_native_capture import update_owned

ROOT = Path(__file__).resolve().parents[1]
PRIVATE = ROOT/'tmp/hydrogen-postflight-contract-01/glycam-fallback-01'


def snapshot(output, prefix, topology, positions):
    from openmm.unit import nanometer

    atoms = list(topology.atoms())
    rows = snapshot_atoms(topology, positions)
    for row, atom in zip(rows, atoms, strict=True):
        row.update(chain_index=atom.residue.chain.index, residue_index=atom.residue.index)
    compressed(output/f'{prefix}.json.gz', dict(atoms=rows,
        positions_nm=[[float(x) for x in p] for p in positions.value_in_unit(nanometer)],
        bonds=[[a.index, b.index] for a, b in topology.bonds()],
        bond_metadata=[dict(type=str(b.type) if b.type is not None else None, order=b.order) for b in topology.bonds()]))


def restore(saved):
    from openmm import Vec3, unit
    from openmm.app import Amide, Aromatic, Double, Single, Topology, Triple, element

    topology = Topology()
    chains, residues, atoms = {}, {}, []
    for row in saved['atoms']:
        chain_id, resid, icode, resname, name = row['key']
        ci, ri = row['chain_index'], row['residue_index']
        if ci not in chains:
            chains[ci] = topology.addChain(chain_id)
        if ri not in residues:
            residues[ri] = topology.addResidue(resname, chains[ci], resid, icode)
        res = residues[ri]
        if (res.chain is not chains[ci] or res.chain.id != chain_id
                or (res.id, res.insertionCode, res.name) != (resid, icode, resname)):
            raise ValueError('Inconsistent serialized topology occurrence')
        atoms.append(topology.addAtom(name, element.Element.getBySymbol(row['element']), res))
    types = {str(t):t for t in (Amide, Aromatic, Double, Single, Triple)}
    for (a, b), metadata in zip(saved['bonds'], saved['bond_metadata'], strict=True):
        label = metadata['type']
        if label is not None and label not in types:
            raise ValueError('Unsupported explicit bond type')
        topology.addBond(atoms[a], atoms[b], types[label] if label is not None else None, metadata['order'])
    if len(atoms) != len(saved['positions_nm']):
        raise ValueError('Atom/position transport mismatch')
    return topology, [Vec3(*p) for p in saved['positions_nm']]*unit.nanometer


def freeze(output):
    import openmm.version

    old_path = ROOT/'docs/plans/oh-native-qualification-plan.json'
    old = json.loads(old_path.read_bytes())
    inputs = dict(old['inputs'])
    for p in (old_path, ROOT/'docs/plans/glycam-hydrogen-fallback.md', PRIVATE/'baseline.json',
              ROOT/'docs/research/glycam-fallback-saved-state.json',
              ROOT/'docs/research/oh-native-qualification-summary.json'):
        inputs[str(p)] = sha(p.read_bytes())
    baseline = json.loads((PRIVATE/'baseline.json').read_bytes())
    saved_control = json.loads((ROOT/'docs/research/glycam-fallback-saved-state.json').read_bytes())
    if saved_control['helper_sha256'] != sha((ROOT/'src/dvbfixer/prepare/fallback.py').read_bytes()):
        raise ValueError('Saved-state control does not bind the current helper')
    for name, digest in baseline['historical_inputs'].items():
        if sha(Path(name).read_bytes()) != digest:
            raise ValueError('Original failed input changed')
        inputs[name] = digest
    sources = {str(p):sha(p.read_bytes()) for p in sorted((ROOT/'src/dvbfixer').rglob('*.py'))}
    sources.update({str(p):sha(p.read_bytes()) for p in sorted((ROOT/'scripts').glob('*.py'))})
    for name in ('test_prepare_glycam_fallback.py', 'test_pdbfixer_usage_inventory.py', 'test_qualify_glycam_fallback.py'):
        p = ROOT/'tests'/name
        sources[str(p)] = sha(p.read_bytes())
    record = dict(schema='dvbfixer.glycam-fallback-plan.v1', python=sys.executable,
        inputs=inputs, source_hashes=sources, baseline=old['baseline'], cases=old['cases'],
        runtime=dict(full_version=openmm.version.version, revision=openmm.version.git_revision),
        protocol=dict(old['protocol'],
            enrollment='Changed-candidate prospective cohort: same twelve declared preparations once each, seed1 per H call; original 11/12 cohort remains immutable and incomplete',
            observation='Explicit chain/residue indices, exact nm positions, bond types/orders, recovery before/after; no native reporter or altered native FF/minimizer',
            readiness='Require 12/12 completed preparations, unique raw PDB identities and valid DatRecord, exact postflight replay and immutable recovery atoms/coordinates/inherited external bonds',
            limits='One draw label per six sources x two FF requests; forcefield-free remains forcefield-free, not native CHARMM glycan or direct ligand-OH qualification'),
        main_merge_authorized=False)
    verify(record)
    write(output, record)


def worker(plan, index, output):
    import openmm

    from dvbfixer.prepare import fallback

    record = json.loads(plan.read_bytes())
    qualify_oh_native.check(record)
    original = fallback.restore_glycam_fallback_bonds
    recoveries = []

    @functools.wraps(original)
    def recovery(topology, positions, *a, **kw):
        number = len(recoveries)+1
        snapshot(output, f'recovery-{number}-before', topology, positions)
        row = dict(number=number, returned=False)
        recoveries.append(row)
        try:
            result = original(topology, positions, *a, **kw)
            row.update(returned=True, added_bonds=result)
            return result
        except BaseException as error:
            row.update(error_type=type(error).__name__, error=str(error))
            raise
        finally:
            snapshot(output, f'recovery-{number}-after', topology, positions)
            update_owned(output/'recovery.json', dict(recoveries=recoveries))

    def native_capture(context, folder, prefix, topology):
        state = context.getState(positions=True, velocities=True, parameters=True, integratorParameters=True)
        for name, value in (('system', context.getSystem()), ('integrator', context.getIntegrator()), ('state', state)):
            with (folder/f'{prefix}-{name}.xml').open('x') as stream:
                stream.write(openmm.XmlSerializer.serialize(value))
        with (folder/f'{prefix}.chk').open('xb') as stream:
            stream.write(context.createCheckpoint())
        snapshot(folder, prefix+'-atoms', topology, state.getPositions())
        return dict(prefix=prefix, platform=platform_properties(context.getPlatform(), context))

    with patch.object(qualify_oh_native, 'snapshot', snapshot), \
            patch.object(qualify_oh_native, 'capture_context', native_capture), \
            patch.object(fallback, 'restore_glycam_fallback_bonds', recovery):
        qualify_oh_native.worker(plan, index, output)


def run(plan, output):
    record = json.loads(plan.read_bytes())
    qualify_oh_native.check(record)
    output.mkdir(parents=True, exist_ok=False)
    (output/'executed-plan.json').write_bytes(plan.read_bytes())
    env = dict(os.environ, PATH=str(Path(record['python']).parent)+os.pathsep+os.environ.get('PATH', ''),
        PYTHONPATH=os.pathsep.join((str(ROOT/'src'), str(ROOT))))
    started, rows = time.monotonic(), []
    for index, case in enumerate(record['cases']):
        folder = output/f'case-{index:02d}'
        folder.mkdir()
        (folder/'input.pdb').write_bytes(Path(case['source']).read_bytes())
        row = dict(index=index, label=case['label'], requested_ff=case['requested_ff'],
                   status='budget_not_executed', returncode=None)
        remaining = record['protocol']['total_seconds']-(time.monotonic()-started)
        if remaining > 0:
            command = [record['python'], str(Path(__file__).resolve()), 'worker', '--plan', str(plan),
                       '--index', str(index), '--output', str(folder)]
            with (folder/'prepare.log').open('xb') as log:
                process = subprocess.Popen(command, cwd=ROOT, env=env, stdout=log,
                                           stderr=subprocess.STDOUT, start_new_session=True)
                try:
                    code = process.wait(timeout=min(remaining, record['protocol']['per_case_seconds']))
                    row.update(status='completed' if code == 0 else 'execution_blocker', returncode=code)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGKILL)
                    process.wait()
                    row['status'] = 'timeout'
        row['artifacts'] = {str(p.relative_to(folder)):sha(p.read_bytes()) for p in sorted(folder.rglob('*')) if p.is_file()}
        write(folder/'result.json', row)
        rows.append(row)
        update_owned(output/'progress.json', dict(plan_sha256=sha(plan.read_bytes()), cases=rows))
        print(json.dumps({k:v for k, v in row.items() if k != 'artifacts'}), flush=True)
    qualify_oh_native.check(record)
    write(output/'evidence.json', dict(plan_sha256=sha(plan.read_bytes()), cases=rows))


def audit(plan, run_path, raw_output, output):
    record = json.loads(plan.read_bytes())
    verify(record)
    with patch.object(audit_oh_native, 'restore', restore):
        audit_oh_native.audit(plan, run_path, raw_output)
    result = json.loads(raw_output.read_bytes())
    extras = []
    for case in result['cases']:
        folder = run_path/f"case-{case['index']:02d}"
        recoveries = json.loads((folder/'recovery.json').read_bytes())['recoveries'] if (folder/'recovery.json').exists() else []
        rows = []
        for row in recoveries:
            before = audit_oh_native.read(folder/f"recovery-{row['number']}-before.json.gz")
            after = audit_oh_native.read(folder/f"recovery-{row['number']}-after.json.gz")
            old = {tuple(sorted(pair)) for pair in before['bonds']}
            new = {tuple(sorted(pair)) for pair in after['bonds']}
            same = before['atoms'] == after['atoms'] and before['positions_nm'] == after['positions_nm']
            internal = all(before['atoms'][a]['residue_index'] == before['atoms'][b]['residue_index'] for a, b in new-old)
            old_metadata = dict(zip((tuple(sorted(p)) for p in before['bonds']), before['bond_metadata'], strict=True))
            new_metadata = dict(zip((tuple(sorted(p)) for p in after['bonds']), after['bond_metadata'], strict=True))
            inherited = old <= new and all(new_metadata[k] == value for k, value in old_metadata.items())
            rows.append(dict(number=row['number'], returned=row['returned'], added_bonds=row.get('added_bonds'),
                added_bond_identities=[[before['atoms'][a]['key'], before['atoms'][b]['key']] for a, b in sorted(new-old)],
                atoms_and_exact_nm_unchanged=same, inherited_bonds_preserved=inherited,
                external_bonds_unchanged=internal and old <= new,
                passed=same and inherited and internal and row['returned']
                    and row.get('added_bonds') == len(new-old)))
        keys, model = [], 0
        if case['output_present']:
            for line in (folder/'prepared.pdb').read_text().splitlines():
                if line.startswith('MODEL'):
                    model += 1
                if line.startswith(('ATOM  ', 'HETATM')):
                    keys.append((model, line[21:22], line[22:26], line[26:27], line[17:21].strip(), line[12:16].strip(), line[16:17]))
        duplicate = [dict(identity=list(key), count=count) for key, count in Counter(keys).items() if count > 1]
        extras.append(dict(index=case['index'], recoveries=rows, raw_output_duplicate_identities=duplicate,
                           raw_output_atom_records=len(keys)))
    valid = all(not c['raw_output_duplicate_identities'] and c['raw_output_atom_records'] > 0
                and all(r['passed'] for r in c['recoveries']) for c in extras)
    result.update(schema='dvbfixer.glycam-fallback-audit.v1', recovery_and_output_audits=extras,
        frozen_scalar_audit_sha256=sha(raw_output.read_bytes()), adapter_sha256=sha(Path(__file__).read_bytes()),
        historical_cohort_status_changed=False, native_charmm_glycan_qualification=False,
        direct_ligand_oh_qualification=False)
    result['summary']['recovery_and_unique_output_passed'] = valid
    result['summary']['qualification_passed'] = result['summary']['qualification_passed'] and valid
    write(output, result)
    verify(record)
    print(json.dumps(result['summary'], indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('plan')
    p.add_argument('--output', type=Path, required=True)
    for name in ('worker', 'run', 'audit'):
        p = sub.add_parser(name)
        p.add_argument('--plan', type=Path, required=True)
        p.add_argument('--output', type=Path, required=True)
        if name == 'worker':
            p.add_argument('--index', type=int, required=True)
        if name == 'audit':
            p.add_argument('--run', type=Path, required=True)
            p.add_argument('--raw-output', type=Path, required=True)
    args = parser.parse_args()
    if args.command == 'plan':
        freeze(args.output.resolve())
    elif args.command == 'worker':
        worker(args.plan.resolve(), args.index, args.output.resolve())
    elif args.command == 'run':
        run(args.plan.resolve(), args.output.resolve())
    else:
        audit(args.plan.resolve(), args.run.resolve(), args.raw_output.resolve(), args.output.resolve())
