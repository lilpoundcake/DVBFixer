"""Saved-only SER-HG boundary audit plus deterministic detector-contract controls."""
from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from pathlib import Path

from scripts.audit_thr2_native_capture import read, verify_files
from scripts.research_dunbrack_prior import sha, write
from scripts.research_ser_hg_capture import check, findings, measurements


def static_detector_sweep():
    """Toy local graphs, not reconstructed originals or extra preparation samples."""
    from openmm import Vec3, unit
    from openmm.app import Topology, element

    from dvbfixer.ffutils import geometry

    lengths = (0.97, 1.199, 1.2, 1.2458133889150464, 1.551, 1.552, 1.553, 1.7)
    result = []
    for length in lengths:
        topology = Topology()
        residue = topology.addResidue('SER', topology.addChain('B'), '126')
        cb = topology.addAtom('CB', element.carbon, residue)
        og = topology.addAtom('OG', element.oxygen, residue)
        hg = topology.addAtom('HG', element.hydrogen, residue)
        topology.addAtom('OXT', element.oxygen, residue)
        topology.addBond(cb, og)
        topology.addBond(og, hg)
        positions = [Vec3(-0.14, 0., 0.), Vec3(0., 0., 0.), Vec3(length/10, 0., 0.), Vec3(0., 0.4, 0.)]*unit.nanometer
        before = [list(p) for p in positions.value_in_unit(unit.nanometer)]
        detected = findings(topology, positions)
        repaired = geometry.repair_misplaced_hydrogens(topology, positions)
        after = [list(p) for p in positions.value_in_unit(unit.nanometer)]
        result.append(dict(input_og_hg_angstrom=length, test_distance_predicate=length < 1.2,
            findings=detected, repairs=repaired, final_og_hg_angstrom=10*math.dist(after[1], after[2]),
            coordinates_unchanged=before == after))
    return dict(enrollment='all eight predeclared O-H distances, one-H/one-heavy-neighbor toy graphs with no close nonparents',
        canonical_oh_angstrom=10*geometry._BOND_LEN_NM['O'], parent_distance_factor=geometry._PARENT_DIST_FACTOR,
        coincidence_tolerance_angstrom=10*geometry._COINCIDENT_TOL_NM,
        distance_only_cutoff_angstrom=10*geometry._BOND_LEN_NM['O']*geometry._PARENT_DIST_FACTOR,
        rows=result, limits='detector/test contract gap, not proof of the lost failure geometry or its stochastic cause')


def pdb_rows(path):
    rows = []
    for line in path.read_text().splitlines():
        if line.startswith(('ATOM  ', 'HETATM')):
            rows.append(dict(key=[line[21], line[22:26].strip(), line[26].strip(), line[17:20].strip(), line[12:16].strip()],
                element=line[76:78].strip(), xyz=[float(line[a:b]) for a, b in ((30, 38), (38, 46), (46, 54))]))
    return rows


def audit(plan_path, root, output):
    from openmm import XmlSerializer
    from openmm.unit import angstrom

    plan = read(plan_path)
    check(plan)
    execution = read(root/'evidence.json')
    if execution['plan_sha256'] != sha(plan_path.read_bytes()) or (root/'executed-plan.json').read_bytes() != plan_path.read_bytes():
        raise ValueError('executed SER-HG plan differs')
    if [(c['label'], c['seed']) for c in execution['cases']] != [(c['label'], c['seed']) for c in plan['cases']]:
        raise ValueError('complete six-case enrollment differs')
    reports = []
    for result in execution['cases']:
        source = root/result['label']
        verify_files(source, result)
        if sha((source/'input.pdb').read_bytes()) != plan['synthetic_input_sha256']:
            raise ValueError('fresh input differs from frozen tracked literal')
        report = {k:v for k, v in result.items() if k != 'artifacts'}
        capture_path = source/'capture.json'
        if capture_path.exists():
            capture = read(capture_path)
            if capture['plan_sha256'] != sha(plan_path.read_bytes()) or capture['case'] != dict(label=result['label'], seed=result['seed']):
                raise ValueError('capture plan/label differs')
            h_calls = []
            for h in capture['h_calls']:
                expected_rng = h['after_rng'] if result['seed'] is None else h['before_rng']
                if h['restored_rng'] != expected_rng:
                    raise ValueError('Python H RNG scope was not preserved')
                before = read(source/f"h-{h['number']}-before-atoms.json.gz")['atoms']
                after = read(source/f"h-{h['number']}-after-atoms.json.gz")['atoms']
                if measurements(before) != h['before'] or measurements(after) != h['after']:
                    raise ValueError('H boundary measurements disagree with coordinate payload')
                h_calls.append(dict(number=h['number'], returned=h['returned'], before=measurements(before), after=measurements(after)))
            native = []
            for call in capture['native']:
                sides = {}
                for side in ('before', 'after'):
                    if side not in call:
                        continue
                    prefix = call[side]['prefix']
                    rows = read(source/f'{prefix}-atoms.json.gz')['atoms']
                    state = XmlSerializer.deserialize((source/f'{prefix}-state.xml').read_text())
                    points = state.getPositions(asNumpy=True).value_in_unit(angstrom).tolist()
                    if points != [r['xyz'] for r in rows]:
                        raise ValueError('serialized native State differs from atom payload')
                    sides[side] = measurements(rows)
                native.append(dict(number=call['number'], returned=call['returned'], tolerance=call['tolerance'],
                    max_iterations=call['max_iterations'], platform=call.get('before', {}).get('platform'), **sides))
            repairs = []
            for repair in capture['repairs']:
                sides = {side:measurements(read(source/f"repair-{repair['number']}-{side}-atoms.json.gz")['atoms']) for side in ('before', 'after')}
                if any(sides[side] != repair[side] for side in sides):
                    raise ValueError('repair measurements disagree with coordinate payload')
                repairs.append({k:v for k, v in repair.items() if k not in ('before', 'after')} | sides)
            report.update(h_calls=h_calls, native=native, repairs=repairs)
        prepared = source/'prepared.pdb'
        if prepared.exists():
            final = measurements(pdb_rows(prepared))
            report['final_pdb'] = final
            report['satisfies_both_original_geometry_assertions'] = final.get('satisfies_original_og_hg_assertion', False) and final.get('satisfies_original_hg_oxt_assertion', False)
        reports.append(report)
    summary = dict(enrolled=len(reports), statuses=dict(Counter(r['status'] for r in reports)),
        native_h_calls=sum(len(r.get('native', [])) for r in reports),
        final_both_assertions_pass=sum(r.get('satisfies_both_original_geometry_assertions', False) for r in reports),
        final_og_hg_assertion_fail=sum('final_pdb' in r and not r['final_pdb'].get('satisfies_original_og_hg_assertion', False) for r in reports))
    write(output, dict(schema='dvbfixer.ser-hg-audit.v1', plan_sha256=sha(plan_path.read_bytes()),
        execution_sha256=sha((root/'evidence.json').read_bytes()), audit_script_sha256=sha(Path(__file__).read_bytes()),
        summary=summary, cases=reports, static_detector=static_detector_sweep(),
        production_changed=False, scientific_thresholds_changed=False, additional_native_minimization=False,
        limits=plan['protocol']['limits'], original_failed_artifacts=plan['original_failed_artifacts']))
    check(plan)
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('plan', 'root', 'output'):
        parser.add_argument('--'+name, type=Path, required=True)
    args = parser.parse_args()
    audit(args.plan.resolve(), args.root.resolve(), args.output.resolve())
