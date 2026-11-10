"""Independent all-whole-mask stage/geometry audit, then held-out reporting.

Commit the geometry-only ledger before opening a native reference. No rebuilding,
optimizer, reflection, score evaluation or native-based candidate selection.
"""
from __future__ import annotations

import argparse
import gzip
import json
import math
import statistics
from dataclasses import asdict
from pathlib import Path

import pdbfixer.pdbfixer as pf
from scripts.audit_template_anchors import fit_check
from scripts.audit_template_causality import (
    audit_points,
    determinant,
    points_equal,
    read_pdb,
    volume,
)
from scripts.branched_template_baseline import preflight
from scripts.compare_branched_completion import measure
from scripts.research_rosetta_joint import sha
from scripts.research_rosetta_repack import write_json
from scripts.rosetta_fixed_contract import choose_scored_own
from scripts.summarize_rosetta_fixed_scoring import source_records_preserved
from scripts.unfiltered_branched_comparison import reporting_scores

from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef
from dvbfixer.pdb_component_reconstruction import _component_ref, _parse_atom_lines
from dvbfixer.protein_completion.branched_source import plan_branched_source
from dvbfixer.protein_completion.contracts import Addition


def median(values):
    values = list(values)
    return statistics.median(values) if values else None


def accuracy(request, additions, pdb, truth):
    points = {a.name: a.position for a in _parse_atom_lines(pdb.decode().splitlines(keepends=True)) if _component_ref(a) == request.residue}
    rounded = tuple(Addition(ExactAtomRef(request.residue, atom), atom[0], points[atom]) for atom in request.missing_atoms)
    return dict(raw=reporting_scores(request, additions, truth), rounded=reporting_scores(request, rounded, truth),
                native_reporting_only=True)


def audit(plan_path, run_root, output):
    plan_raw = plan_path.read_bytes()
    plan = json.loads(plan_raw)
    for filename, digest in {**plan['source_hashes'], **plan['inputs']}.items():
        if sha(Path(filename).read_bytes()) != digest:
            raise ValueError('frozen expansion input/implementation changed')
    if preflight() != plan['environment']:
        raise ValueError('scientific environment changed')
    if (run_root/'executed-plan.json').read_bytes() != plan_raw:
        raise ValueError('executed expansion plan differs')
    run_raw = (run_root/'evidence.json').read_bytes()
    run = json.loads(run_raw)
    if run['plan_sha256'] != sha(plan_raw) or not run['installed_source_unchanged']:
        raise ValueError('expansion run/source binding differs')
    expected_ids = {c['id']+'-'+method for c in plan['cases'] for method in plan['methods']}
    if set(run['results']) != expected_ids:
        raise ValueError('method/case omitted or added')
    mask_root = Path(plan['mask_root'])
    masks = json.loads((mask_root/'plan.json').read_bytes())
    for filename, digest in masks['source_hashes'].items():
        if sha(Path(filename).read_bytes()) != digest:
            raise ValueError('original mask implementation changed')
    expected_whole = [c['id'] for c in masks['cases'] if 'CB' in c['mask'] and len(c['mask']) == 3]
    if expected_whole != [c['id'] for c in plan['cases']]:
        raise ValueError('complete original whole-mask inventory differs')
    score_plan_raw = Path(plan['prior_score_plan']).read_bytes()
    score_plan = json.loads(score_plan_raw)
    score_root = Path(plan['prior_score_plan']).parent
    for filename, digest in score_plan['source_hashes'].items():
        if sha(Path(filename).read_bytes()) != digest:
            raise ValueError('prior scoring implementation changed')
    score_run = json.loads(Path(plan['prior_score_run']).read_bytes())
    if score_run['plan_sha256'] != sha(score_plan_raw):
        raise ValueError('prior own-scored baseline plan differs')
    score_declarations = {c['id']: c for c in score_plan['cases']}
    ledger, replay, stage_count, iterate_count, seed_count, baseline_seeds, seen_repeats, fits = {}, {}, 0, 0, 0, 0, 0, 0
    for case in plan['cases']:
        identifier, target = case['id'], ComponentInstanceRef(**case['target'])
        key = target.chain_id, target.sequence_number, target.insertion_code
        source = Path(case['input_path']).read_bytes()
        if sha(source) != case['input_sha256'] or source != (mask_root/case['source']).read_bytes():
            raise ValueError('source differs from original masked bytes')
        request = plan_branched_source(source, target).request
        if set(request.missing_atoms) != set(case['mask']) or len(request.missing_atoms) != 3:
            raise ValueError('actual missing-atom mask differs from whole-sidechain declaration')
        observed = {tuple(a['key']): a['xyz'] for a in read_pdb(source)}
        target_observed = {k[3]: p for k, p in observed.items() if k[:3] == key}
        data = case['anchor_data']
        template_raw = (Path(pf.__file__).parent/data['resource']).read_bytes()
        if sha(template_raw) != data['resource_sha256']:
            raise ValueError('canonical template resource changed')
        template_atoms = read_pdb(template_raw)
        template = {a['key'][3]: a['xyz'] for a in template_atoms}
        names = [a['key'][3] for a in template_atoms if a['key'][3] in target_observed]
        p = [[x/10. for x in target_observed[n]] for n in names]
        q = [[x/10. for x in template[n]] for n in names]
        if data['names'] != names or data['reference_points_nm'] != p or data['template_points_nm'] != q:
            raise ValueError('source/resource anchor mapping differs')
        case_ledger = dict(fixture=case['fixture'], ncac_previously_executed=case['ncac_previously_executed'], methods={})
        case_replay = dict(request=request, source=source, methods={})
        for method in plan['methods']:
            name = identifier+'-'+method
            root = run_root/name
            trace_raw = (root/'evidence.json').read_bytes()
            trace = json.loads(trace_raw)
            if sha(trace_raw) != run['results'][name]['trace_sha256'] or not trace['hooks_restored']:
                raise ValueError('trace digest/hook restoration differs')
            if (root/'input.pdb').read_bytes() != source:
                raise ValueError('method input bytes differ')
            by_seed, payloads = {}, {}
            for event in trace['events']:
                raw = (root/event['coordinates']).read_bytes()
                if sha(raw) != event['coordinate_sha256']:
                    raise ValueError('stage coordinate hash differs')
                atoms = json.loads(gzip.decompress(raw))['atoms']
                coords = {tuple(a['key']): a['xyz'] for a in atoms}
                if len(coords) != len(atoms) or not observed.keys() <= coords.keys():
                    raise ValueError('stage identities lost or duplicated')
                movement = max(math.dist(coords[k], xyz) for k, xyz in observed.items())
                if movement > 1e-6 or not math.isclose(movement, event['observed_max_movement_angstrom'], abs_tol=1e-12):
                    raise ValueError('stage moved observed coordinates')
                points = {k[3]: xyz for k, xyz in coords.items() if k[:3] == key}
                if points != event['target_points_angstrom']:
                    raise ValueError('stage target differs from full coordinates')
                audit_points(points, event['ca_volume_angstrom3'], case['residue_name'], event['target_cip'])
                by_seed.setdefault(event['seed'], []).append((event, atoms))
                payloads[event['coordinates']] = raw
                stage_count += 1
            if list(by_seed) != list(range(1, len(by_seed)+1)) or len(by_seed) > 5 or len(by_seed) != run['results'][name]['seed_count']:
                raise ValueError('helper seed inventory differs')
            seed_count += len(by_seed)
            for iteration in trace['minimizer_iterations']:
                audit_points(iteration['target_points_angstrom'], iteration['ca_volume_angstrom3'])
                if iteration['observed_max_movement_angstrom'] > 1e-6 or any(math.dist(iteration['target_points_angstrom'][n], xyz) > 1e-6 for n, xyz in target_observed.items()):
                    raise ValueError('minimizer iterate moved original observations')
                iterate_count += 1
            constructions = [e for e in trace['events'] if e['stage'].startswith('construction-')]
            if len(constructions) != 2*len(by_seed) or len(trace['overlays']) != len(constructions):
                raise ValueError('construction/overlay count differs')
            calls = None
            if method == 'ncac':
                transport_raw = (root/'anchor-policy.json').read_bytes()
                if sha(transport_raw) != run['results'][name]['transport_sha256']:
                    raise ValueError('used-anchor transport digest differs')
                transport = json.loads(transport_raw)
                calls = transport['calls']
                if not transport['research_hook_restored'] or len(calls) != len(constructions):
                    raise ValueError('used-anchor fit/restoration inventory differs')
                if case['ncac_previously_executed']:
                    previous = Path(plan['seen_anchor_run'])/(identifier+'-ncac-normal')
                    if not run['results'][name]['seen_ncac_trace_exact'] or trace_raw != (previous/'evidence.json').read_bytes():
                        raise ValueError('previously seen NCAC trace does not reproduce exactly')
                    for filename, raw in payloads.items():
                        if raw != (previous/filename).read_bytes():
                            raise ValueError('previously seen NCAC stage payload differs')
                    seen_repeats += 1
            for i, (overlay, stage) in enumerate(zip(trace['overlays'], constructions, strict=True)):
                if method == 'ncac':
                    fit_check(calls[i], overlay, names, p, q, 'ncac')
                    fits += 1
                if overlay['seed'] != stage['seed']:
                    raise ValueError('fit/construction seed order differs')
                matrix = overlay['rotation']
                det = determinant(matrix)
                if not math.isclose(abs(det), 1., abs_tol=1e-10):
                    raise ValueError('overlay transform determinant invalid')
                transformed = {}
                for atom, xyz in template.items():
                    centered = [xyz[k]/10.+overlay['translate_template_nm'][k] for k in range(3)]
                    transformed[atom] = [10.*(sum(matrix[k][j]*centered[j] for j in range(3))+overlay['translate_reference_nm'][k]) for k in range(3)]
                hybrid = dict(transformed)
                hybrid.update(target_observed)
                if not points_equal({atom: hybrid[atom] for atom in stage['target_points_angstrom']}, stage['target_points_angstrom']):
                    raise ValueError('saved fit does not reproduce actual fixed-backbone hybrid')
                if not math.isclose(volume(transformed), det*volume(template), rel_tol=1e-11, abs_tol=1e-11):
                    raise ValueError('full-template volume law differs')
            returns = []
            for seed, events in by_seed.items():
                initial = next(e for e, _ in events if e['stage'] == 'construction-heavy')
                event, atoms = next((e, atoms) for e, atoms in events if e['stage'] == 'returned-seed')
                if method == 'stock':
                    prior = json.loads(Path(case['prior_proposal_path']).read_bytes())
                    snapshot = next(s for s in prior['snapshots'] if s['seed'] == seed)
                    old = {(a['component']['chain_id'], a['component']['sequence_number'], a['component']['insertion_code'], a['name']): (element, xyz)
                           for a, element, xyz in snapshot['atoms']}
                    if {tuple(a['key']): (a['element'], a['xyz']) for a in atoms} != old or not trace['stock_normal_exact_seed_reproduction'] or len(prior['snapshots']) != len(by_seed):
                        raise ValueError('stock seed state differs from untraced original')
                    baseline_seeds += 1
                added = tuple(Addition(ExactAtomRef(request.residue, atom), atom[0], tuple(event['target_points_angstrom'][atom])) for atom in request.missing_atoms)
                measured, pdb = measure(source, request, added, anchored=False)
                returns.append(dict(seed=seed, initial_volume=initial['ca_volume_angstrom3'], initial_cip=initial['target_cip'],
                                    returned_volume=event['ca_volume_angstrom3'], returned_cip=event['target_cip'],
                                    common_audit=asdict(measured), dynamics_calls=sum(e['stage'] == 'before-dynamics' for e, _ in events)))
                final_added, final_pdb = added, pdb
            case_ledger['methods'][method] = dict(seed_count=len(returns), seeds=returns,
                                                  selected_returned_seed=returns[-1]['seed'], common_accepted=returns[-1]['common_audit']['accepted'],
                                                  trace_sha256=sha(trace_raw))
            case_replay['methods'][method] = dict(additions=final_added, pdb=final_pdb)
        declaration = score_declarations[identifier]
        own_rows = [r for r in score_run['results'][identifier]['cells'] if r['arm'] == 'own']
        own_selection = choose_scored_own(own_rows)
        if own_selection != score_run['results'][identifier]['selection']:
            raise ValueError('prior own score selection replay differs')
        selected = own_selection['selected_index']
        if selected is None:
            case_ledger['methods']['own_scored'] = dict(common_accepted=False, reason='prior own pool empty')
        else:
            cell = next(c for c in declaration['cells'] if c['arm'] == 'own' and c['index'] == selected)
            pdb = (score_root/cell['pdb']).read_bytes()
            if sha(pdb) != score_plan['artifacts'][cell['pdb']]:
                raise ValueError('prior own-scored artifact digest differs')
            if not source_records_preserved(source, pdb, declaration['generated']):
                raise ValueError('prior own-scored candidate changed source observations')
            points = {a.name: a.position for a in _parse_atom_lines(pdb.decode().splitlines(keepends=True)) if _component_ref(a) == request.residue}
            added = tuple(Addition(ExactAtomRef(request.residue, atom), atom[0], points[atom]) for atom in request.missing_atoms)
            measured, checked = measure(source, request, added, anchored=False)
            if checked != pdb:
                raise ValueError('prior own-scored rounded materialization differs')
            case_ledger['methods']['own_scored'] = dict(common_accepted=measured.accepted, common_audit=asdict(measured), selected_index=selected)
            case_replay['methods']['own_scored'] = dict(additions=added, pdb=pdb)
        ledger[identifier] = case_ledger
        replay[identifier] = case_replay
    # All proposals and all input-only/common geometry are now committed. Native
    # coordinate files are opened ONLY below this publication boundary.
    geometry_path = output.with_name(output.stem+'-geometry.json')
    if output.exists() or geometry_path.exists():
        raise FileExistsError('exclusive new geometry/audit destinations required')
    write_json(geometry_path, dict(schema='dvbfixer.ncac-expansion-geometry.v1', cases=ledger, native_coordinates_read=False,
                                  experiment_plan_sha256=sha(plan_raw), run_sha256=sha(run_raw), audit_script_sha256=sha(Path(__file__).read_bytes())))
    reports = {}
    for case in plan['cases']:
        identifier = case['id']
        request, source = replay[identifier]['request'], replay[identifier]['source']
        reference_raw = (mask_root/case['reference']).read_bytes()
        if sha(reference_raw) != masks['artifacts'][case['reference']]:
            raise ValueError('native reporting-only reference digest differs')
        truth = {a.name: a.position for a in _parse_atom_lines(reference_raw.decode().splitlines(keepends=True)) if _component_ref(a) == request.residue}
        native = tuple(Addition(ExactAtomRef(request.residue, atom), atom[0], truth[atom]) for atom in request.missing_atoms)
        native_audit, _ = measure(source, request, native, anchored=False)
        methods = {}
        for method, saved in replay[identifier]['methods'].items():
            methods[method] = dict(common_accepted=ledger[identifier]['methods'][method]['common_accepted'],
                                   accuracy=accuracy(request, saved['additions'], saved['pdb'], truth))
        reports[identifier] = dict(fixture=case['fixture'], ncac_previously_executed=case['ncac_previously_executed'],
                                   native_common_audit=asdict(native_audit), methods=methods)
    strata = {}
    for label, subset in [('all', list(reports.values())), ('seen_four', [r for r in reports.values() if r['ncac_previously_executed']]),
                          ('new_to_variant_seven', [r for r in reports.values() if not r['ncac_previously_executed']])]:
        coverage = {method: sum(r['methods'].get(method, {}).get('common_accepted', False) for r in subset) for method in ('stock', 'ncac', 'own_scored')}
        mutual = [r for r in subset if all(r['methods'].get(method, {}).get('common_accepted', False) for method in ('ncac', 'own_scored'))]
        stock_mutual = [r for r in subset if all(r['methods'].get(method, {}).get('common_accepted', False) for method in ('stock', 'ncac'))]
        strata[label] = dict(cases=len(subset), common_passes=coverage, native_common_passes=sum(r['native_common_audit']['accepted'] for r in subset),
                             ncac_vs_own_mutual_pass_count=len(mutual), stock_vs_ncac_mutual_pass_count=len(stock_mutual),
                             ncac_vs_own_mutual_rounded_exact_medians={m: median(r['methods'][m]['accuracy']['rounded']['exact_name_rmsd_angstrom'] for r in mutual) for m in ('ncac', 'own_scored')},
                             ncac_vs_own_mutual_generated_val_symmetry_medians={m: median(r['methods'][m]['accuracy']['rounded']['generated_only_val_symmetry_rmsd_angstrom'] for r in mutual) for m in ('ncac', 'own_scored')},
                             stock_vs_ncac_mutual_rounded_exact_medians={m: median(r['methods'][m]['accuracy']['rounded']['exact_name_rmsd_angstrom'] for r in stock_mutual) for m in ('stock', 'ncac')},
                             all_raw_diagnostic_rounded_exact_medians={m: median(r['methods'][m]['accuracy']['rounded']['exact_name_rmsd_angstrom'] for r in subset if m in r['methods']) for m in ('stock', 'ncac', 'own_scored')})
    aggregate = dict(cases=len(plan['cases']), method_arms=len(run['results']), returned_seeds=seed_count, full_stage_snapshots=stage_count,
                     target_minimizer_iterates=iterate_count, used_ncac_fits=fits, exact_stock_seed_repeats=baseline_seeds, exact_seen_ncac_trace_repeats=seen_repeats,
                     dynamics_step_calls=sum(s['dynamics_calls'] for c in ledger.values() for m in ('stock', 'ncac') for s in c['methods'][m]['seeds']),
                     strata=strata)
    write_json(output, dict(schema='dvbfixer.ncac-expansion-audit.v1', aggregate=aggregate, cases=reports,
                            geometry_ledger_name=geometry_path.name, geometry_ledger_sha256=sha(geometry_path.read_bytes()),
                            experiment_plan=plan, experiment_plan_sha256=sha(plan_raw), run_path=str(run_root), run_evidence_sha256=sha(run_raw),
                            audit_script_sha256=sha(Path(__file__).read_bytes()), production_changed=False, scientific_qualification=False,
                            scope=plan['scope']))
    print(json.dumps(aggregate, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    audit(args.plan, args.run, args.output)
