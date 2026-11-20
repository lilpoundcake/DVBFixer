"""Frozen nine-cell retained-context heavy comparison; H is auxiliary, not scored."""
from __future__ import annotations

import argparse
import hashlib
import json
import statistics
from dataclasses import replace
from pathlib import Path

from scripts.chemical_reconstruction_all20_controls import (
    FFBONDED,
    HEAVY_POLICY,
    authority_for,
    heldout_after_proposal,
    load_parent,
)
from scripts.chemical_reconstruction_ccd_controls import ROOT, _json, _write
from scripts.chemical_reconstruction_incoming_type_controls import check_pins as check_type_pins

from dvbfixer.chemical_reconstruction.amide_authority import FFBONDED_SHA256, load_amide_authority
from dvbfixer.chemical_reconstruction.candidate_audit import audit_candidate
from dvbfixer.chemical_reconstruction.comparison_metrics import measure_common_heavy
from dvbfixer.chemical_reconstruction.construction import proposal_admission, propose
from dvbfixer.chemical_reconstruction.incoming_source import (
    audit_incoming_source,
    bind_incoming_source,
)
from dvbfixer.chemical_reconstruction.incoming_type import candidate_family, prove_incoming_type
from dvbfixer.chemical_reconstruction.metric_projection import propose_projection
from dvbfixer.chemical_reconstruction.observed_peptide import (
    _parse_source,
    compile_peptide_request_plan,
    select_interior_target,
)
from dvbfixer.chemical_reconstruction.planar_port import (
    load_planar_port_authority,
    planar_port_constraints,
)
from dvbfixer.chemical_reconstruction.projected_hydrogen_audit import (
    audit_projected_hydrogen_result,
)
from dvbfixer.chemical_reconstruction.projected_hydrogen_contracts import (
    ProjectedHeavyParent,
    ProjectedHydrogenResult,
)
from dvbfixer.chemical_reconstruction.projected_hydrogens import complete_projected_hydrogens
from dvbfixer.chemical_reconstruction.projection_audit import audit_projection_result
from dvbfixer.chemical_reconstruction.projection_contracts import (
    ProjectionResult,
    append_projection_constraints,
)
from dvbfixer.chemical_reconstruction.proposals import ProposalResult
from dvbfixer.chemical_reconstruction.research_decode import (
    decode_research_record,
    parse_research_json,
)
from dvbfixer.chemical_reconstruction.source_preparation import (
    audit_source_transport,
    prepare_source,
    prepare_source_probe,
)
from dvbfixer.chemical_reconstruction.template_comparison import (
    SeededHeavyBaseline,
    audit_seeded_baseline,
    masked_probe_source,
    seeded_heavy_baseline,
)
from dvbfixer.chemical_reconstruction.type_proven_port import TypeProvenPlanarPortDeclaration
from dvbfixer.domain.structure_identity import ExactAtomRef

PROTOCOL = ROOT/'docs/plans/chemical-reconstruction-comparison-v1.md'
COHORT = ROOT/'docs/plans/chemical-reconstruction-comparison-cohort-v1.json'
PROTOCOL_SHA256 = '212b611e2936cacc31c685cfc83c6c4159c0b4f2006a50bb01d6bd88c1360826'
COHORT_SHA256 = 'a59f537e181e68ad92e24da37ac75b437e9d8d5ccf98b3a26314cc180a26f24d'


def check_pins():
    catalog = check_type_pins()
    for path, pin in ((PROTOCOL, PROTOCOL_SHA256), (COHORT, COHORT_SHA256)):
        if hashlib.sha256(path.read_bytes()).hexdigest() != pin:
            raise ValueError('comparison-frozen-protocol-or-cohort')
    cohort = json.loads(COHORT.read_bytes())
    for row in cohort['sources']:
        if hashlib.sha256((ROOT/row['path']).read_bytes()).hexdigest() != row['sha256']:
            raise ValueError('comparison-frozen-source')
    return cohort, catalog


def prepare_cell(row, mask, catalog):
    data = (ROOT/row['path']).read_bytes()
    preparation = prepare_source(data)
    _, atoms = _parse_source(data)
    chain = next(a.chain_id for a in atoms if a.record == 'ATOM  ')
    identity = select_interior_target(data, 'PRO', chain=chain)
    if (identity.chain_id, identity.sequence_number, identity.insertion_code) != (row['chain'], row['residue_number'], row['insertion_code']):
        raise ValueError('comparison-target-enrollment-binding')
    authority = authority_for('PRO', 'PRO', identity, {})
    probe = prepare_source_probe(preparation, authority, withheld=frozenset(mask))
    if {a.name for a in probe.request.ownership.generated} != set(mask):
        raise ValueError('comparison-natural-incomplete-target-is-not-frozen-mask')
    if any(probe.request.ownership.generated & {c.link.component_atom, c.link.other_atom} for c in probe.request.connections):
        raise ValueError('comparison-unsupported-target-sidechain-link')
    source = bind_incoming_source(preparation, probe, authority)
    proof = prove_incoming_type(probe.request, source, catalog, load_parent(candidate_family(source.residue_name, catalog)[0]))
    if audit_incoming_source(preparation, probe, authority, source):
        raise ValueError('comparison-incoming-byte-admission')
    parameters = load_planar_port_authority(FFBONDED)
    types = {r[0]: r[1] for s in catalog.states if s.name == 'PRO' for r in s.atoms}
    n, ca, cd = (ExactAtomRef(identity, a) for a in ('N', 'CA', 'CD'))
    declaration = TypeProvenPlanarPortDeclaration(source.incoming, n, ca, cd,
        (proof.atom_type, types['N'], types['CA'], types['CD']), authority.digest, parameters.digest, incoming_proof=proof)
    request = append_projection_constraints(probe.request, planar_port_constraints(probe.request, authority, declaration, parameters))
    if not audit_source_transport(preparation, replace(probe, request=request)).passed:
        raise ValueError('comparison-port-source-transport')
    masked, removed = masked_probe_source(preparation.prepared_text.encode('ascii'), probe.request)
    return preparation, authority, probe, declaration, parameters, request, masked, removed


def measure_arms(request, baseline, projection, template, reference):
    measurements = {}
    for arm, result in (('original', baseline), ('projection', projection), ('seeded-PDBFixer', template)):
        for phase in ('raw', 'rounded'):
            supplied = getattr(result, phase) if result is not None else None
            additions = supplied if arm == 'seeded-PDBFixer' else supplied.additions if supplied is not None else ()
            complete = len(additions or ()) == len(request.ownership.generated) and {o.identity for o in additions or ()} == request.ownership.generated
            measurements[f'{arm}-{phase}'] = (measure_common_heavy(request, HEAVY_POLICY, additions, reference)
                if complete else {'status': 'unscored-incomplete-or-unavailable'})
    measurements['native-reference-diagnostic-only'] = measure_common_heavy(request, HEAVY_POLICY, reference, reference)
    return measurements


def run_controls(output: Path):
    cohort, catalog = check_pins()
    H_parameters = load_amide_authority(FFBONDED, expected_sha256=FFBONDED_SHA256)
    output.mkdir(parents=True, exist_ok=False)
    records = []
    for row in cohort['sources']:
        for mask in cohort['masks']:
            name = row['id']+'-'+ '-'.join(mask)
            directory = output/name
            directory.mkdir()
            item = {'id': name, 'source': row['id'], 'mask': mask, 'baseline_called': False,
                    'projection_called': False, 'H_called': False, 'template_called': False}
            evidence = {'control': item}
            (directory/'original.pdb').write_bytes((ROOT/row['path']).read_bytes())
            try:
                try:
                    preparation, authority, probe, declaration, parameters, request, masked, removed = prepare_cell(row, mask, catalog)
                except ValueError as exc:
                    item.update(status='source-refused', reason=str(exc))
                    continue
                (directory/'prepared.pdb').write_bytes(preparation.prepared_text.encode('ascii'))
                (directory/'masked-template-input.pdb').write_bytes(masked)
                evidence.update(authority=authority, probe=probe, declaration=declaration, port_parameters=parameters,
                                request=request, removed_mask_records=removed)
                original = probe.request
                plan = compile_peptide_request_plan(original, HEAVY_POLICY)
                baseline = propose(original, plan)
                evidence.update(baseline_plan=plan, baseline=baseline)
                item.update(baseline_called=True, baseline_status=baseline.status)
                projection = None
                if baseline.status == 'declared-geometry-passed':
                    anchors = tuple(ExactAtomRef(authority.identity, atom) for atom in ('CA', 'N', 'C'))
                    projection = propose_projection(original, request, plan, baseline, anchors)
                    evidence['projection'] = projection
                    item.update(projection_called=True, projection_status=projection.status)
                    if audit_projection_result(original, request, plan, baseline, anchors, projection):
                        raise RuntimeError('comparison-projection-readback-failed')
                    if projection.status == 'projected-geometry-passed':
                        parent = ProjectedHeavyParent(original, request, projection, authority, declaration, parameters)
                        result = complete_projected_hydrogens(parent, H_parameters, backbone_types=None)
                        evidence.update(parent=parent, H_result=result, H_parameters=H_parameters)
                        item.update(H_called=True, H_status=result.status)
                        if audit_projected_hydrogen_result(parent, H_parameters, None, result):
                            raise RuntimeError('comparison-H-readback-failed')
                # Independent heavy arm runs even if our baseline/projection rejects.
                template = seeded_heavy_baseline(original, masked)
                evidence['template'] = template
                item.update(template_called=True, template_status=template.status, template_reason=template.reason)
                # First withheld coordinate extraction, AFTER all arm/H phases.
                reference = tuple(o for o in heldout_after_proposal('observed', authority, directory) if o.identity in original.ownership.generated)
                measurements = measure_arms(request, baseline, projection, template, reference)
                evidence.update(postproposal_reference=reference, measurements=measurements)
                item.update(status='measured', reason='all-available-arms-retained-not-qualification')
            finally:
                _write(directory/'evidence.json', evidence)
                records.append(item)
    check_pins()
    summary = {'schema': 'dvbfixer.chemical-comparison-controls.v1', 'protocol_sha256': PROTOCOL_SHA256,
               'cohort_sha256': COHORT_SHA256, 'declared_cells': 9, 'source_systems': 3, 'controls': records,
               'H_comparison_or_native_quality_qualified': False, 'production_enabled': False}
    _write(output/'summary.json', summary)
    return summary


def readback_controls(output: Path):
    cohort, catalog = check_pins()
    H_parameters = load_amide_authority(FFBONDED, expected_sha256=FFBONDED_SHA256)
    records = []
    for row in cohort['sources']:
        for mask in cohort['masks']:
            name = row['id']+'-'+ '-'.join(mask)
            directory = output/name
            saved = parse_research_json((directory/'evidence.json').read_bytes())
            def dto(key, cls, saved=saved):
                return decode_research_record(json.dumps(saved[key], allow_nan=False).encode(), cls)
            item = saved['control']
            if (item['id'], item['source'], item['mask']) != (name, row['id'], mask):
                raise ValueError('comparison-readback-enrollment')
            if (directory/'original.pdb').read_bytes() != (ROOT/row['path']).read_bytes():
                raise ValueError('comparison-readback-original-bytes')
            records.append(item)
            try:
                preparation, authority, probe, declaration, parameters, request, masked, removed = prepare_cell(row, mask, catalog)
            except ValueError as exc:
                if (set(saved) != {'control'} or item['status'] != 'source-refused' or item['reason'] != str(exc)
                        or any(item[k] for k in ('baseline_called', 'projection_called', 'H_called', 'template_called'))):
                    raise ValueError('comparison-readback-source-refusal') from exc
                continue
            if ((directory/'prepared.pdb').read_bytes() != preparation.prepared_text.encode('ascii')
                    or (directory/'masked-template-input.pdb').read_bytes() != masked
                    or saved['removed_mask_records'] != list(removed)
                    or any(saved[key] != _json(value) for key, value in
                           [('authority', authority), ('probe', probe), ('declaration', declaration), ('port_parameters', parameters), ('request', request)])):
                raise ValueError('comparison-readback-source-masked-authority-binding')
            baseline = dto('baseline', ProposalResult)
            plan = compile_peptide_request_plan(probe.request, HEAVY_POLICY)
            if (saved['baseline_plan'] != _json(plan) or baseline.plan_digest != plan.digest
                    or baseline.admission != proposal_admission(probe.request, plan)
                    or not item['baseline_called'] or item['baseline_status'] != baseline.status):
                raise ValueError('comparison-readback-baseline-plan')
            for phase in ('raw', 'rounded'):
                candidate = getattr(baseline, phase)
                if candidate is not None and audit_candidate(probe.request, plan, candidate) != getattr(baseline, phase+'_audit'):
                    raise ValueError('comparison-readback-baseline-audit')
            projection = dto('projection', ProjectionResult) if item['projection_called'] else None
            if (projection is not None) != (baseline.status == 'declared-geometry-passed'):
                raise ValueError('comparison-readback-projection-stage')
            if projection is not None:
                anchors = tuple(ExactAtomRef(authority.identity, a) for a in ('CA', 'N', 'C'))
                if audit_projection_result(probe.request, request, plan, baseline, anchors, projection) or item['projection_status'] != projection.status:
                    raise ValueError('comparison-readback-projection')
            if item['H_called'] != bool(projection and projection.status == 'projected-geometry-passed'):
                raise ValueError('comparison-readback-H-stage')
            if item['H_called']:
                parent = ProjectedHeavyParent(probe.request, request, projection, authority, declaration, parameters)
                result = dto('H_result', ProjectedHydrogenResult)
                if (saved['parent'] != _json(parent) or saved['H_parameters'] != _json(H_parameters)
                        or audit_projected_hydrogen_result(parent, H_parameters, None, result) or item['H_status'] != result.status):
                    raise ValueError('comparison-readback-H')
            template = dto('template', SeededHeavyBaseline)
            if (template.request_digest != probe.request.digest or template.masked_source_sha256 != hashlib.sha256(masked).hexdigest()
                    or not item['template_called'] or item['template_status'] != template.status or item['template_reason'] != template.reason):
                raise ValueError('comparison-readback-template-provenance')
            audit_seeded_baseline(probe.request, masked, template)
            reference = tuple(o for o in heldout_after_proposal('observed', authority, directory) if o.identity in probe.request.ownership.generated)
            if (saved['postproposal_reference'] != _json(reference)
                    or saved['measurements'] != _json(measure_arms(request, baseline, projection, template, reference))
                    or item['status'] != 'measured' or item['reason'] != 'all-available-arms-retained-not-qualification'):
                raise ValueError('comparison-readback-measurements')
    summary = parse_research_json((output/'summary.json').read_bytes())
    expected = {'schema': 'dvbfixer.chemical-comparison-controls.v1', 'protocol_sha256': PROTOCOL_SHA256,
                'cohort_sha256': COHORT_SHA256, 'declared_cells': 9, 'source_systems': 3, 'controls': records,
                'H_comparison_or_native_quality_qualified': False, 'production_enabled': False}
    if summary != expected:
        raise ValueError('comparison-readback-summary')
    return expected


def summarize_controls(output: Path):
    """Descriptive complete-pair accounting; nested masks are NOT independent."""
    summary = readback_controls(output)
    rows = []
    for item in summary['controls']:
        saved = parse_research_json((output/item['id']/'evidence.json').read_bytes())
        measurements = saved.get('measurements', {})
        rows.append({'id': item['id'], 'source': item['source'], 'mask': item['mask'],
                     'control': item, 'measurements': measurements})
    aggregate = {}
    for phase in ('raw', 'rounded'):
        counts, comparisons = {}, {}
        for arm in ('original', 'projection', 'seeded-PDBFixer'):
            available = [r for r in rows if 'rmsd_angstrom' in r['measurements'].get(arm+'-'+phase, {})]
            counts[arm] = {'complete_scored_cells': [r['id'] for r in available],
                'common_geometry_passes': sum(r['measurements'][arm+'-'+phase]['geometry']['passed'] for r in available),
                'severe_overlap_pair_count': sum(len(r['measurements'][arm+'-'+phase]['severe_overlap_pairs']) for r in available)}
        for left, right in (('projection', 'original'), ('projection', 'seeded-PDBFixer'), ('seeded-PDBFixer', 'original')):
            pairs = [r for r in rows if all('rmsd_angstrom' in r['measurements'].get(arm+'-'+phase, {}) for arm in (left, right))]
            deltas = [r['measurements'][left+'-'+phase]['rmsd_angstrom']-r['measurements'][right+'-'+phase]['rmsd_angstrom'] for r in pairs]
            comparisons[left+'-minus-'+right] = {'paired_cells': [r['id'] for r in pairs],
                'excluded_cells': [r['id'] for r in rows if r not in pairs],
                'RMSD_deltas_angstrom': deltas, 'worsened_cells': sum(d > 0 for d in deltas),
                'median_delta_angstrom': statistics.median(deltas) if deltas else None}
        aggregate[phase] = {'arms': counts, 'complete_pair_comparisons': comparisons}
    return {'schema': 'dvbfixer.chemical-comparison-descriptive.v1', 'declared_cells': 9, 'source_systems': 3,
            'rows': rows, 'aggregate': aggregate, 'independent_trial_count_is_not_nine': True,
            'H_native_or_general_superiority_qualified': False}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    parser.add_argument('--readback', action='store_true')
    args = parser.parse_args()
    print(json.dumps(readback_controls(args.output) if args.readback else run_controls(args.output), indent=2))
