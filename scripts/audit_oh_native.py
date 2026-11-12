"""Independent bounded OH measurements and exact saved-state matched replay."""
from __future__ import annotations

import argparse
import gzip
import json
import math
from collections import Counter
from pathlib import Path

from scripts.qualify_hydrogen_postflight import load, reconstruct
from scripts.research_dunbrack_prior import sha, write
from scripts.research_ha_minimizer_path import verify


def read(path):
    return json.loads(gzip.decompress(path.read_bytes()))


def restore(saved):
    from openmm import Vec3, unit

    top, _ = reconstruct(saved['atoms'], saved['bonds'])
    return top, [Vec3(*p) for p in saved['positions_nm']]*unit.nanometer


def measure(saved):
    from dvbfixer.ffutils import PROTEIN_RESIDUES
    from dvbfixer.residue_registry import is_sugar_resname

    atoms, points = saved['atoms'], saved['positions_nm']
    groups = {}
    for i, row in enumerate(atoms):
        groups.setdefault(tuple(row['key'][:4]), []).append(i)
    counts, violations, short, cross = Counter(), [], [], []
    for a, b in saved['bonds']:
        if atoms[a]['element'] == 'H' and atoms[b]['element'] == 'O':
            a, b = b, a
        if atoms[a]['element'] != 'O' or atoms[b]['element'] != 'H':
            continue
        key = atoms[b]['key']
        name = key[3]
        category = ('water' if name in ('HOH', 'WAT', 'SOL') else 'sugar' if is_sugar_resname(name)
                    else 'protein' if name in PROTEIN_RESIDUES else 'other')
        counts[category+':'+name] += 1
        distance = 10*math.dist(points[a], points[b])
        contacts = [(10*math.dist(points[i], points[b]), atoms[i]['key'])
                    for i in groups[tuple(key[:4])] if i != a and atoms[i]['element'] != 'H']
        nearest = min(contacts) if contacts else None
        if not math.isfinite(distance) or distance >= 1.2 or (nearest and nearest[0] <= 1.0):
            violations.append(dict(hydrogen=key, oh_angstrom=distance,
                nearest_nonparent_angstrom=nearest[0] if nearest else None))
        if distance < 0.8:
            short.append(dict(hydrogen=key, oh_angstrom=distance))
        for i, atom in enumerate(atoms):
            if atom['element'] == 'H' or atom['key'][:4] == key[:4]:
                continue
            d = 10*math.dist(points[i], points[b])
            if d <= 1.0:
                cross.append(dict(hydrogen=key, heavy=atom['key'], distance_angstrom=d))
    return dict(oh_counts=dict(counts), violations=violations, short_oh_warnings=short,
        cross_residue_warnings=cross, passed=not violations)


def replay(saved, module):
    from openmm.unit import nanometer
    from scripts.research_legacy_start_controls import snapshot_atoms

    top, positions = restore(saved)
    flagged = {f.h_atom.index for f in module.detect_misplaced_hydrogens(top, positions)}
    result = dict(status='completed')
    try:
        result['repairs'] = module.repair_misplaced_hydrogens(top, positions)
    except Exception as error:
        result.update(status='error', error_type=type(error).__name__)
    after = [[float(x) for x in p] for p in positions.value_in_unit(nanometer)]
    changed = [i for i, (a, b) in enumerate(zip(saved['positions_nm'], after, strict=True)) if a != b]
    result.update(measurement=measure(dict(saved, positions_nm=after)),
        heavy_preserved=all(saved['atoms'][i]['element'] == 'H' for i in changed),
        valid_h_preserved=all(i in flagged for i in changed),
        bonds_preserved=[[a.index, b.index] for a, b in top.bonds()] == saved['bonds'],
        identity_preserved=[r['key'] for r in snapshot_atoms(top, positions)] == [r['key'] for r in saved['atoms']],
        changed_atom_keys=[saved['atoms'][i]['key'] for i in changed])
    return result, after


def audit(plan_path, run_path, output):
    from openmm.app import PDBFile
    from openmm.unit import nanometer
    from scripts.research_legacy_start_controls import snapshot_atoms

    from dvbfixer.ffutils import geometry
    from dvbfixer.ffutils.dat import DatRecord

    plan = json.loads(plan_path.read_bytes())
    verify(plan)
    evidence = json.loads((run_path/'evidence.json').read_bytes())
    if evidence['plan_sha256'] != sha(plan_path.read_bytes()):
        raise ValueError('execution plan mismatch')
    baseline = load('oh_native_baseline', Path(plan['baseline']))
    reports = []
    for case in evidence['cases']:
        folder = run_path/f"case-{case['index']:02d}"
        for name, digest in case['artifacts'].items():
            if sha((folder/name).read_bytes()) != digest:
                raise ValueError(f'changed retained artifact: {folder/name}')
        report = {k:v for k, v in case.items() if k != 'artifacts'}
        capture = json.loads((folder/'capture.json').read_bytes()) if (folder/'capture.json').exists() else {}
        report['h_routes'] = [{k:r.get(k) for k in ('number', 'returned', 'forcefield_provided', 'forcefield_xmls', 'error_type')}
                              for r in capture.get('h_calls', [])]
        report['native_platforms'] = [r['before']['platform'] for r in capture.get('native', [])]
        report['repairs'] = []
        for row in capture.get('repairs', []):
            before = read(folder/f"repair-{row['number']}-before.json.gz")
            actual = read(folder/f"repair-{row['number']}-after.json.gz")
            old, _ = replay(before, baseline)
            new, positions = replay(before, geometry)
            # Atom snapshots include coordinates; compare identity separately, not their xyz fields.
            exact = (positions == actual['positions_nm'] and before['bonds'] == actual['bonds']
                     and [r['key'] for r in before['atoms']] == [r['key'] for r in actual['atoms']]
                     and (new['status'] == 'completed') == row['returned']
                     and new.get('repairs') == row.get('repairs') and new.get('error_type') == row.get('error_type'))
            report['repairs'].append(dict(number=row['number'], before=measure(before), baseline=old,
                candidate=new, exact_actual_replay=exact,
                new_refusal=new['status'] != 'completed' and old['status'] == 'completed'))
        dat = folder/'prepared.dat'
        report['dat'] = dict(present=dat.exists(), valid=False)
        if dat.exists():
            try:
                record = DatRecord.load(dat)
                report['dat'].update(valid=True, variant_overrides=record.variant_overrides)
            except Exception as error:
                report['dat']['error_type'] = type(error).__name__
        report['output_present'] = (folder/'prepared.pdb').exists()
        if report['output_present']:
            try:
                pdb = PDBFile(str(folder/'prepared.pdb'))
                report['output_measurement'] = measure(dict(atoms=snapshot_atoms(pdb.topology, pdb.positions),
                    positions_nm=[list(p) for p in pdb.positions.value_in_unit(nanometer)],
                    bonds=[[a.index, b.index] for a, b in pdb.topology.bonds()]))
            except Exception as error:
                report['output_measurement'] = dict(passed=False, error_type=type(error).__name__)
        reports.append(report)
    rows = [r for c in reports for r in c['repairs']]
    blockers = [dict(case=c['index'], repair=r['number']) for c in reports for r in c['repairs']
        if not r['exact_actual_replay'] or r['new_refusal'] or not all(r['candidate'][k]
            for k in ('heavy_preserved', 'valid_h_preserved', 'bonds_preserved', 'identity_preserved'))
        or not r['candidate']['measurement']['passed']]
    summary = dict(declared_cases=len(plan['cases']), retained_cases=len(reports),
        completed_cases=sum(c['status'] == 'completed' for c in reports), matched_repair_states=len(rows),
        exact_replays=sum(r['exact_actual_replay'] for r in rows), blockers=blockers,
        qualification_passed=len(reports) == len(plan['cases']) and not blockers and bool(rows)
            and all(c['status'] == 'completed' and c['dat']['valid'] and c['output_present'] and c['repairs']
                and c['output_measurement']['passed'] for c in reports))
    write(output, dict(schema='dvbfixer.oh-native-audit.v1', plan_sha256=sha(plan_path.read_bytes()),
        evidence_sha256=sha((run_path/'evidence.json').read_bytes()), summary=summary, cases=reports,
        limits=plan['protocol']['limits'], main_merged=False))
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('plan', 'run', 'output'):
        parser.add_argument('--'+name, type=Path, required=True)
    args = parser.parse_args()
    audit(args.plan.resolve(), args.run.resolve(), args.output.resolve())
