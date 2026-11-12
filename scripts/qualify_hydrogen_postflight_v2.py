"""Amended qualification: retained v1 failure plus distance and contact controls."""
from __future__ import annotations

import argparse
import gzip
import io
import json
import math
from pathlib import Path

from scripts.qualify_hydrogen_postflight import EVIDENCE, ROOT, load, reconstruct
from scripts.qualify_hydrogen_postflight import run as distance_run
from scripts.research_dunbrack_prior import sha, write
from scripts.research_ha_minimizer_path import verify
from scripts.research_legacy_start_controls import compressed, snapshot_atoms


def contract(rows, bonds):
    """Independent explicit identity/distance checks, not production findings."""
    adjacent = {i:[] for i in range(len(rows))}
    for a, b in bonds:
        adjacent[a].append(b)
        adjacent[b].append(a)
    values = []
    for i, row in enumerate(rows):
        if row['element'] != 'H':
            continue
        parents = [j for j in adjacent[i] if rows[j]['element'] != 'H']
        if len(parents) != 1:
            raise ValueError('qualification H must have one explicit heavy parent')
        parent = parents[0]
        if rows[parent]['element'] != 'O':
            continue
        distance = math.dist(row['xyz'], rows[parent]['xyz'])
        contacts = [(math.dist(row['xyz'], other['xyz']), other['key']) for j, other in enumerate(rows)
            if j != parent and other['element'] != 'H' and other['key'][:4] == row['key'][:4]]
        nearest = min(contacts) if contacts else None
        values.append(dict(hydrogen=row['key'], oh_angstrom=distance,
            nearest_nonparent_heavy_angstrom=nearest[0] if nearest else None,
            nearest_nonparent_heavy=nearest[1] if nearest else None,
            passed=math.isfinite(distance) and distance < 1.2 and (nearest is None or nearest[0] > 1.0)))
    return dict(passed=all(v['passed'] for v in values), hydrogens=values)


def freeze(output):
    from openmm.app import PDBFile

    from dvbfixer.force_field_naming import NamingConversionRequest, convert_force_field_naming

    old_plan = EVIDENCE/'qualification-plan.json'
    old = json.loads(old_plan.read_bytes())
    verify(old)
    cases_path = ROOT/'tests/test_hydrogen_postflight_contacts.py'
    cases = load('hydrogen_postflight_contact_cases', cases_path)
    folder = EVIDENCE/'failed-candidate-v1'
    original = folder/'repaired.pdb'
    normalized = convert_force_field_naming(NamingConversionRequest(original.read_text(), 'amber19sb', 'standard')).pdb_text
    old_lines = [line for line in original.read_text().splitlines() if line.startswith(('ATOM  ', 'HETATM'))]
    new_lines = [line for line in normalized.splitlines() if line.startswith(('ATOM  ', 'HETATM'))]
    if len(old_lines) != len(new_lines) or any(a[30:54] != b[30:54] for a, b in zip(old_lines, new_lines, strict=True)):
        raise ValueError('name normalization changed observed coordinates')
    canonical = folder/'canonical.pdb'
    with canonical.open('x') as stream:
        stream.write(normalized)
    pdb = PDBFile(io.StringIO(normalized))
    snapshot = folder/'canonical-topology.json.gz'
    compressed(snapshot, dict(atoms=snapshot_atoms(pdb.topology, pdb.positions), bonds=[[a.index, b.index] for a, b in pdb.topology.bonds()]))
    protocol = EVIDENCE/'protocol-v2-before-implementation.md'
    with protocol.open('xb') as stream:
        stream.write((ROOT/'docs/plans/hydrogen-postflight-contract.md').read_bytes())
    paths = [old_plan, cases_path, protocol, ROOT/'src/dvbfixer/force_field_naming.py', ROOT/'src/dvbfixer/domain/force_field_naming.py']
    paths.extend(p for p in folder.iterdir() if p.is_file())
    inputs = {**old['inputs'], **old['source_hashes'], **{str(p):sha(p.read_bytes()) for p in paths}}
    write(output, dict(schema='dvbfixer.hydrogen-postflight-v2-plan.v1', old_plan=str(old_plan),
        baseline=old['baseline'], failed_candidate=str(folder/'geometry.py'), retained_failure=str(snapshot),
        cases_source=str(cases_path), inputs=inputs, source_hashes={str(Path(__file__).resolve()):sha(Path(__file__).read_bytes())},
        contact_cases=[dict(oh_length=length, contact=contact) for length in cases.OH_LENGTHS for contact in cases.CONTACTS],
        policy=dict(distance_and_native_records=54, contact_records=12, real_failure_records=1,
            native_sampling=False, oh_upper_angstrom=1.2, nonparent_lower_angstrom=1.0,
            radial_admission_clearance_angstrom=1.01,
            allowed_refusal='typed HydrogenGeometryError only for initially invalid blocked contact toy graphs; real failure and original 54 must repair/retain successfully',
            limits='deterministic local graph and saved-artifact qualification, not broad hydroxyl force-field accuracy or stereochemistry qualification')))


def run(plan_path, output, public_output):
    from dvbfixer.ffutils import geometry as candidate

    plan = json.loads(plan_path.read_bytes())
    verify(plan)
    candidate_hash = sha(Path(candidate.__file__).read_bytes())
    output.mkdir(parents=True, exist_ok=False)
    (output/'executed-plan.json').write_bytes(plan_path.read_bytes())
    distance_file = output/'distance-evidence.json'
    distance_run(Path(plan['old_plan']), output/'distance-stage', distance_file)
    distance = json.loads(distance_file.read_bytes())
    distance_checks = []
    for row in distance['cases']:
        folder = output/'distance-stage'/f"case-{row['index']:02d}"
        raw = json.loads(gzip.decompress((folder/'input.json.gz').read_bytes()))
        repaired = json.loads(gzip.decompress((folder/'candidate.json.gz').read_bytes()))
        distance_checks.append(dict(index=row['index'], **contract(repaired['atoms'], raw['bonds'])))
    baseline = load('hydrogen_postflight_v2_baseline', Path(plan['baseline']))
    failed = load('hydrogen_postflight_v2_failed', Path(plan['failed_candidate']))
    cases = load('hydrogen_postflight_contact_cases', Path(plan['cases_source']))
    enrollment = []
    for case in plan['contact_cases']:
        top, positions = cases.contact_graph(**case)
        enrollment.append(dict(kind='contact-graph', **case, rows=snapshot_atoms(top, positions), bonds=[[a.index, b.index] for a, b in top.bonds()]))
    original = json.loads(gzip.decompress(Path(plan['retained_failure']).read_bytes()))
    enrollment.append(dict(kind='retained-v1-SER-failure', rows=original['atoms'], bonds=original['bonds']))
    reports = []
    for index, case in enumerate(enrollment):
        folder = output/f'extra-{index:02d}'
        folder.mkdir()
        report = {k:v for k, v in case.items() if k not in ('rows', 'bonds')}
        report['before'] = contract(case['rows'], case['bonds'])
        report['arms'] = {}
        for name, module in (('baseline', baseline), ('failed-v1', failed), ('candidate', candidate)):
            top, positions = reconstruct(case['rows'], case['bonds'])
            before = snapshot_atoms(top, positions)
            outcome = dict(status='completed')
            try:
                outcome['repairs'] = module.repair_misplaced_hydrogens(top, positions)
            except Exception as error:
                outcome.update(status='error', error_type=type(error).__name__, error=str(error))
            after = snapshot_atoms(top, positions)
            compressed(folder/f'{name}.json.gz', dict(atoms=after))
            outcome.update(after=contract(after, case['bonds']), identity_preserved=[r['key'] for r in before] == [r['key'] for r in after],
                heavy_preserved=all(a['xyz'] == b['xyz'] for a, b in zip(before, after, strict=True) if a['element'] != 'H'),
                valid_input_unchanged=before == after if report['before']['passed'] else None)
            outcome['accepted_refusal'] = (case['kind'] == 'contact-graph' and not report['before']['passed']
                and outcome.get('error_type') == 'HydrogenGeometryError')
            report['arms'][name] = outcome
        reports.append(report)
    outcomes = [r['arms']['candidate'] for r in reports]
    real = reports[-1]['arms']['candidate']
    valid = all(a['heavy_preserved'] and a['identity_preserved'] and a['valid_input_unchanged'] is not False
        and (a['accepted_refusal'] or (a['status'] == 'completed' and a['after']['passed'])) for a in outcomes)
    summary = dict(records=67, original_distance_stage_passed=distance['summary']['qualification_passed'],
        original_54_contact_checks_passed=sum(r['passed'] for r in distance_checks),
        extra_completed_safe=sum(a['status'] == 'completed' and a['after']['passed'] for a in outcomes),
        extra_typed_refusals=sum(a['accepted_refusal'] for a in outcomes),
        real_retained_failure_repaired=real['status'] == 'completed' and real['after']['passed'],
        qualification_passed=valid and distance['summary']['qualification_passed'] and all(r['passed'] for r in distance_checks)
            and real['status'] == 'completed' and real['after']['passed'])
    write(public_output, dict(schema='dvbfixer.hydrogen-postflight-v2-qualification.v1', plan_sha256=sha(plan_path.read_bytes()),
        candidate_source_sha256=candidate_hash, script_sha256=sha(Path(__file__).read_bytes()), summary=summary,
        distance_stage=distance, distance_stage_contact_checks=distance_checks, extra_cases=reports,
        native_sampling=False, force_field_changed=False, limits=plan['policy']['limits']))
    verify(plan)
    if sha(Path(candidate.__file__).read_bytes()) != candidate_hash:
        raise ValueError('candidate source changed during qualification')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('plan')
    p.add_argument('--output', type=Path, required=True)
    p = sub.add_parser('run')
    for name in ('plan', 'output', 'public-output'):
        p.add_argument('--'+name, type=Path, required=True)
    args = parser.parse_args()
    if args.command == 'plan':
        freeze(args.output.resolve())
    else:
        run(args.plan.resolve(), args.output.resolve(), args.public_output.resolve())
