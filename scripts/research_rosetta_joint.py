"""Bounded exact joint controls; no production change or interaction approximation.

Enrollment precedes new joint reconstruction. Two already studied systems and
residues are deliberately reused: this is a joint-mechanism control, not a new
independent-family accuracy cohort. All failures and every enumerated combination
are retained. Frozen single-target scoring is invoked twice to audit both targets.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path

from scripts.research_rosetta_fixed_scoring import score_cell
from scripts.research_rosetta_repack import identity, write_json
from scripts.rosetta_fixed_contract import choose_scored_own

POLICY = Path('tmp/pdbfixer-research/rosetta-fixed-plan-03/plan.json')
SYSTEMS = ('1CRN', '1AKI')


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def source_hashes() -> dict:
    paths = [Path(__file__).resolve().relative_to(Path('.').resolve()), *sorted(Path('src/dvbfixer/protein_completion').glob('*.py')),
             Path('src/dvbfixer/ffutils/dat.py'),
             Path('src/dvbfixer/diagnose/geometry.py'), Path('src/dvbfixer/domain/structure_identity.py'),
             Path('src/dvbfixer/pdb_component_reconstruction.py'),
             *[Path('scripts')/p for p in ('prospective_branched_cohort.py', 'unfiltered_branched_comparison.py',
                                          'research_rosetta_fixed_scoring.py', 'research_rosetta_repack.py',
                                          'rosetta_fixed_contract.py', 'branched_stereo_oracle.py')]]
    return {str(p): sha(p.read_bytes()) for p in paths}


def mask_classes(names: tuple[str, str]) -> tuple:
    from dvbfixer.protein_completion.branched import BRANCHES

    a, b = (BRANCHES[n] for n in names)
    pair_a, pair_b = (a.first, a.second), (b.first, b.second)
    return (('distal-distal', (pair_a, pair_b)),
            ('whole-whole', (('CB', *pair_a), ('CB', *pair_b))),
            ('anchored-distal', ((a.first,), pair_b)))


def enroll(output: Path) -> None:
    from dataclasses import asdict

    from scripts.prospective_branched_cohort import project
    from scripts.unfiltered_branched_comparison import select_targets

    policy_raw = POLICY.read_bytes()
    policy = json.loads(policy_raw)
    for filename, digest in policy['source_hashes'].items():
        if sha(Path(filename).read_bytes()) != digest:
            raise ValueError('chosen single-target score implementation changed')
    payloads, cases = {}, []
    for system in SYSTEMS:
        original = Path(f'tests/fixtures/protein_context/scoring_followup/{system}.pdb').read_bytes()
        reference = project(original)
        payloads[system+'-reference.pdb'] = reference
        targets = select_targets(reference.decode().splitlines(keepends=True))
        for group in (0, 1):
            chosen = targets[group*2:group*2+2]
            names = tuple(n for n, _ in chosen)
            for label, masks in mask_classes(names):
                identifier = f'{system}-group{group}-{label}'
                removed = {(t.chain_id, t.sequence_number, t.insertion_code): set(mask)
                           for (_, t), mask in zip(chosen, masks, strict=True)}
                source = b''.join(l for l in reference.splitlines(keepends=True) if not (
                    l.startswith(b'ATOM  ') and l[12:16].strip().decode() in removed.get(
                        (l[21:22].decode(), l[22:26].strip().decode(), l[26:27].strip().decode()), set())))
                filename = identifier+'/source.pdb'
                payloads[filename] = source
                cases.append(dict(id=identifier, system=system, mask_class=label,
                                  targets=[dict(residue_name=n, residue=asdict(t), mask=mask)
                                           for (n, t), mask in zip(chosen, masks, strict=True)],
                                  source=filename, reference=system+'-reference.pdb'))
    output.mkdir(parents=True, exist_ok=False)
    for name, data in payloads.items():
        p = output/name
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(data)
    write_json(output/'enrollment.json', dict(
        schema='dvbfixer.exact-joint-enrollment.v1', cases=cases,
        artifacts={n: sha(data) for n, data in payloads.items()}, source_hashes=source_hashes(),
        original_source_hashes={s: sha(Path(f'tests/fixtures/protein_context/scoring_followup/{s}.pdb').read_bytes()) for s in SYSTEMS},
        single_score_policy_sha256=sha(policy_raw), score_policy=policy,
        selection='first two VAL/first two THR in file order, consecutive two-target groups; three mask classes; no outcome filtering',
        max_targets=2, max_combinations=9,
        candidate_gate='unchanged mixed individual/raw-joint/rounded/global stereo guards; no scoring rescue',
        selection_policy='minimum full-pose ref2015 over every old-gate-feasible joint combination; stable enumeration tie; any missing feasible score blocks selection',
        baselines=['independent per-target clearance choices, which may fail joint guards', 'existing raw-joint-clearance choice with canonical rounded materialization'],
        coupling_diagnostic='E(i,j)-E(i,j0)-E(i0,j)+E(i0,j0); reference is lexicographically first source-feasible indices, never native; full Cartesian score matrix required; a single-option axis is uninformative',
        scope='previously studied systems/constituents; new joint mechanisms only, not independent prospective accuracy validation',
        public_capability=False))


def verify(record: dict, root: Path) -> None:
    for filename, digest in record['source_hashes'].items():
        if sha(Path(filename).read_bytes()) != digest:
            raise ValueError('joint implementation changed after enrollment')
    for filename, digest in record['artifacts'].items():
        if sha((root/filename).read_bytes()) != digest:
            raise ValueError('joint input artifact changed')


def prepare(enrollment_root: Path, output: Path) -> None:
    from dataclasses import asdict

    from scripts.branched_stereo_oracle import rdkit_labels
    from scripts.research_rosetta_repack import observed_atoms

    from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef
    from dvbfixer.pdb_component_reconstruction import _component_ref, _parse_atom_lines
    from dvbfixer.protein_completion.contracts import Addition, Candidate
    from dvbfixer.protein_completion.mixed_cohort import (
        complete_mixed_cohort,
        materialize_mixed_cohort,
        plan_mixed_cohort,
    )
    from dvbfixer.protein_completion.source import _insert_candidates
    from dvbfixer.protein_completion.validation import sidechain_bonds

    raw = (enrollment_root/'enrollment.json').read_bytes()
    enrollment = json.loads(raw)
    verify(enrollment, enrollment_root)
    output.mkdir(parents=True, exist_ok=False)
    cases, artifacts = [], {}

    def save(name, data):
        p = output/name
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(data)
        artifacts[name] = sha(data)

    for declaration in enrollment['cases']:
        name = declaration['id']
        source = (enrollment_root/declaration['source']).read_bytes()
        save(name+'/source.pdb', source)
        save(name+'/reference.pdb', (enrollment_root/declaration['reference']).read_bytes())
        targets = tuple(ComponentInstanceRef(**t['residue']) for t in declaration['targets'])
        try:
            plan = plan_mixed_cohort(source, targets)
        except ValueError as error:
            cases.append(dict(declaration, status='admission_refused', reason=str(error)))
            continue
        result = complete_mixed_cohort(source, plan)
        if len(plan.targets) != 2 or len(result.attempts) > enrollment['max_combinations']:
            raise ValueError('enrolled exact two-target combination budget exceeded')
        save(name+'/source-plan.json', (json.dumps(asdict(plan), indent=2)+'\n').encode())
        save(name+'/reconstruction.json', (json.dumps(asdict(result), indent=2)+'\n').encode())
        requests = tuple(t.request for t in plan.targets)
        cells = []
        for i, attempt in enumerate(result.attempts):
            candidates = tuple(t.attempts[j].candidate for t, j in zip(result.targets, attempt.indices, strict=True))
            pdb = _insert_candidates(source, requests, candidates)
            reason, eligible = None, False
            if attempt.accepted:
                try:
                    materialized = materialize_mixed_cohort(source, plan, candidates)
                except ValueError as error:
                    reason = str(error)
                else:
                    if materialized.pdb != pdb:
                        raise ValueError('canonical materialization differs from diagnostic insertion')
                    path = output/f'{name}/joint-{i}.dat'
                    materialized.dat_record().save(path, verbose=False)
                    artifacts[f'{name}/joint-{i}.dat'] = sha(path.read_bytes())
                    eligible = True
            else:
                reason = 'raw_joint_geometry_rejected'
            # Independent target CIP is an additional non-negotiable check.
            atoms = _parse_atom_lines(pdb.decode().splitlines(keepends=True))
            cips = []
            for target in plan.targets:
                points = {a.name: a.position for a in atoms if _component_ref(a) == target.request.residue}
                labels = rdkit_labels(points, target.request.residue_name)
                expected = {'CA': 'S', 'CB': 'R'} if target.request.residue_name == 'THR' else {'CA': 'S'}
                cips.append(labels)
                if labels != expected:
                    eligible = False
                    reason = 'independent_target_CIP_mismatch'
            filename = f'{name}/joint-{i}.pdb'
            save(filename, pdb)
            cells.append(dict(id=f'joint-{i}', arm='own', index=i, indices=attempt.indices, pdb=filename,
                              old_source_accepted=attempt.accepted, common_accepted=eligible,
                              raw_joint_audits=[asdict(a) for a in attempt.audits],
                              rounded_refusal=reason, target_cip=cips))
        reference_atoms = _parse_atom_lines((enrollment_root/declaration['reference']).read_text().splitlines(keepends=True))
        native = []
        for request in requests:
            points = {a.name: a.position for a in reference_atoms if _component_ref(a) == request.residue}
            native.append(Candidate(request.digest, tuple(Addition(ExactAtomRef(request.residue, n), n[0], points[n])
                                                          for n in request.missing_atoms), None))
        native_pdb = _insert_candidates(source, requests, tuple(native))
        save(name+'/native-control.pdb', native_pdb)
        cells.append(dict(id='native-control', arm='native_control', index=-1, indices=None,
                          pdb=name+'/native-control.pdb', old_source_accepted=False, common_accepted=False))
        # Use the complete native control solely to index graph identities here;
        # its coordinates never enter candidates or score selection.
        lookup = {(a['chain'], str(a['resid']), a['icode'].strip(), a['atom']): identity(a)
                  for a in observed_atoms(native_pdb)}

        def key(atom, lookup=lookup):
            c = atom.component
            return lookup[(c.chain_id, c.sequence_number, c.insertion_code, atom.name)]

        views = []
        for request in requests:
            internal = [(key(ExactAtomRef(request.residue, a)), key(ExactAtomRef(request.residue, b)))
                        for a, b in sidechain_bonds(request.residue_name)
                        if (request.residue.chain_id, request.residue.sequence_number, request.residue.insertion_code, a) in lookup
                        and (request.residue.chain_id, request.residue.sequence_number, request.residue.insertion_code, b) in lookup]
            views.append(dict(target_key=key(ExactAtomRef(request.residue, 'CA'))[:4], target_internal=internal,
                              expected_external=[(key(a), key(b)) for a, b in request.links]))
        independent = [t.selected_index for t in result.targets]
        independent_joint_index = next((i for i, a in enumerate(result.attempts) if list(a.indices) == independent), None)
        cases.append(dict(declaration, status='planned', views=views, cells=cells,
                          reconstruction_status=result.status, reconstruction_reason=result.reason,
                          individual_attempt_counts=[len(r.attempts) for r in result.targets],
                          individual_feasible_counts=[sum(a.audit.accepted for a in r.attempts) for r in result.targets],
                          independent_indices=independent, independent_joint_index=independent_joint_index,
                          raw_clearance_joint_index=result.selected_index,
                          generated=[key(ExactAtomRef(r.residue, n)) for r in requests for n in r.missing_atoms]))
    (output/'executed-enrollment.json').write_bytes(raw)
    write_json(output/'plan.json', dict(schema='dvbfixer.exact-joint-score-plan.v1', cases=cases,
                                      artifacts=artifacts, source_hashes=enrollment['source_hashes'],
                                      enrollment_sha256=sha(raw), enrollment=enrollment))


def merge_views(records: list, snapshots: list) -> dict:
    """No target may be silently dropped, scored in a different pose or averaged."""
    if len(records) != 2 or len(snapshots) != 2:
        raise ValueError('exactly two independently audited target views required')
    if any(r['score_status'] != 'scored' or r['import_audit']['accepted'] is not True for r in records):
        return dict(score_status='target_view_import_blocker', total_score=None)
    for key in ('heavy_atoms', 'hydrogens', 'virtual_atoms', 'residues', 'external', 'allowed_extra'):
        if snapshots[0][key] != snapshots[1][key]:
            return dict(score_status='inconsistent_target_view_pose', total_score=None)
    a, b = (r['total_score'] for r in records)
    if any(not isinstance(v, (int, float)) or isinstance(v, bool) or not math.isfinite(v) for v in (a, b)) or a != b:
        return dict(score_status='inconsistent_target_view_score', total_score=None)
    return dict(score_status='scored', total_score=a)


def context_signature(snapshot: dict, views: list) -> str:
    targets = {tuple(v['target_key'][:3]) for v in views}
    allowed = {tuple(k) for k in snapshot['allowed_extra']}
    value = dict(residues=snapshot['residues'], external=snapshot['external'],
                 extra_heavy=[a for a in snapshot['heavy_atoms'] if identity(a) in allowed],
                 virtual_atoms=snapshot['virtual_atoms'],
                 fixed_context_hydrogens=[h for h in snapshot['hydrogens']
                                          if (h['chain'], h['resid'], h['icode']) not in targets])
    return sha(json.dumps(value, sort_keys=True, allow_nan=False).encode())


def coupling_diagnostic(rows: list[dict]) -> dict:
    """Conditional full-score nonadditivity, not an atom-pair or kcal/mol energy."""
    own = [r for r in rows if r['arm'] == 'own']
    if not own:
        return dict(status='empty_source_feasible_matrix')
    matrix = {tuple(r['indices']): r for r in own}
    if len(matrix) != len(own):
        raise ValueError('duplicate joint indices')
    left = sorted({i for i, _ in matrix})
    right = sorted({j for _, j in matrix})
    if len(matrix) != len(left)*len(right) or any(r['score_status'] != 'scored' for r in own):
        return dict(status='incomplete_score_matrix')
    if any(isinstance(r['total_score'], bool) or not isinstance(r['total_score'], (float, int))
           or not math.isfinite(r['total_score']) for r in own):
        return dict(status='incomplete_score_matrix')
    a, b = left[0], right[0]
    values = [dict(indices=[i, j], delta=matrix[i, j]['total_score']-matrix[i, b]['total_score']
                   -matrix[a, j]['total_score']+matrix[a, b]['total_score']) for i in left for j in right]
    return dict(status='informative_two_axis_matrix' if min(len(left), len(right)) > 1 else 'uninformative_single_option_axis',
                reference_indices=[a, b], dimensions=[len(left), len(right)],
                max_absolute_delta=max(abs(v['delta']) for v in values), values=values)


def score(plan_root: Path, output: Path) -> None:
    raw = (plan_root/'plan.json').read_bytes()
    plan = json.loads(raw)
    verify(plan, plan_root)
    policy = plan['enrollment']['score_policy']
    import pyrosetta

    if policy['expected_rosetta_revision'] not in pyrosetta.version():
        raise ValueError('unexpected Rosetta build')
    pyrosetta.init(policy['init_flags'])
    sf = pyrosetta.create_score_function(policy['score_function'])
    output.mkdir(parents=True, exist_ok=False)
    (output/'executed-plan.json').write_bytes(raw)
    write_json(output/'environment.json', dict(version=pyrosetta.version(),
                                              weights={str(t): sf.get_weight(t) for t in sf.get_nonzero_weighted_scoretypes()}))
    results = {}
    for case in plan['cases']:
        folder = output/case['id']
        folder.mkdir()
        if case['status'] != 'planned':
            results[case['id']] = case
            continue
        rows, signature = [], None

        def execute(cell, case=case, folder=folder, rows=rows):
            nonlocal signature
            cell_root = folder/cell['id']
            cell_root.mkdir()
            records, snapshots = [], []
            try:
                for i, view in enumerate(case['views']):
                    destination = cell_root/f'view-{i}'
                    destination.mkdir()
                    record = score_cell(pyrosetta, sf, plan_root, destination, view, cell)
                    write_json(destination/'score.json', record)
                    records.append(record)
                    snapshots.append(json.loads(gzip.decompress((destination/'pose-snapshot.json.gz').read_bytes())))
                merged = merge_views(records, snapshots)
                if merged['score_status'] == 'scored':
                    current = context_signature(snapshots[0], case['views'])
                    merged['context_signature'] = current
                    if signature is None:
                        signature = current
                    elif current != signature:
                        merged.update(score_status='inconsistent_context_blocker', total_score=None)
            except Exception as error:
                merged = dict(score_status='execution_blocker', total_score=None, error=str(error), error_type=type(error).__name__)
            row = dict(cell, **merged, views=records)
            write_json(cell_root/'score.json', row)
            rows.append(row)

        for cell in case['cells']:
            if cell['arm'] == 'own':
                execute(cell)
        selection = choose_scored_own(rows)
        write_json(folder/'selection-before-native.json', selection)
        if selection['selected_index'] is not None:
            chosen = next(c for c in case['cells'] if c['arm'] == 'own' and c['index'] == selection['selected_index'])
            for suffix in ('.pdb', '.dat'):
                with (folder/('selected-candidate'+suffix)).open('xb') as stream:
                    stream.write((plan_root/Path(chosen['pdb']).with_suffix(suffix)).read_bytes())
        for cell in case['cells']:
            if cell['arm'] == 'native_control':
                execute(cell)
        results[case['id']] = dict(status='evaluated', cells=rows, selection=selection)
        write_json(folder/'result.json', results[case['id']])
    write_json(output/'evidence.json', dict(schema='dvbfixer.exact-joint-raw.v1', plan_sha256=sha(raw),
                                          results=results, public_capability=False))
    print(json.dumps(dict(cases=len(results), cells=sum(len(r.get('cells', [])) for r in results.values()),
                         scored=sum(c['score_status'] == 'scored' for r in results.values() for c in r.get('cells', []))), indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('enroll')
    p.add_argument('--output', type=Path, required=True)
    p = sub.add_parser('prepare')
    p.add_argument('--enrollment', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p = sub.add_parser('score')
    p.add_argument('--plan', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.command == 'enroll':
        enroll(args.output)
    elif args.command == 'prepare':
        prepare(args.enrollment, args.output)
    else:
        score(args.plan, args.output)


if __name__ == '__main__':
    main()
