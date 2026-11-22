"""Frozen saved-coordinate diagnosis, not a reconstruction or counterfactual search."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from scripts import chemical_reconstruction_comparison_controls as comparison
from scripts.chemical_reconstruction_ccd_controls import ROOT, _json, _write

from dvbfixer.chemical_reconstruction.closure_diagnostics import (
    fixed_neighbour_support,
    measure_terms,
    trace_frames,
)
from dvbfixer.chemical_reconstruction.contracts import (
    Observation,
    ReconstructionRequest,
    digest_record,
)
from dvbfixer.chemical_reconstruction.observed_peptide import (
    ObservedPeptideProbe,
    compile_peptide_request_plan,
)
from dvbfixer.chemical_reconstruction.projection_contracts import ProjectionResult
from dvbfixer.chemical_reconstruction.proposals import GeometryPlan, ProposalResult
from dvbfixer.chemical_reconstruction.research_decode import (
    decode_research_record,
    parse_research_json,
)
from dvbfixer.chemical_reconstruction.template_comparison import SeededHeavyBaseline

PROTOCOL = ROOT/'docs/plans/chemical-reconstruction-closure-diagnosis-v1.md'
LEDGER = ROOT/'docs/plans/chemical-reconstruction-closure-diagnosis-inputs-v1.json'
PROTOCOL_SHA256 = '5c68d2e8a380fce7bbed3aef78594b0b3dbefc48046be5f5105306b2a4cfea8a'
LEDGER_SHA256 = '0bfa42be42d2d8ab472f2c715d6d2cd73e8714270f3e52ac9044ef92a9c8964c'


def _pinned(path: Path, pin: str) -> bytes:
    with path.open('rb') as handle:
        data = handle.read(64*1024*1024+1)
    if len(data) > 64*1024*1024 or hashlib.sha256(data).hexdigest() != pin:
        raise ValueError('closure-diagnosis-frozen-artifact-pin:'+str(path))
    return data


def qualified_inputs():
    _pinned(PROTOCOL, PROTOCOL_SHA256)
    ledger = json.loads(_pinned(LEDGER, LEDGER_SHA256))
    root = ROOT/ledger['root']
    _pinned(root/'summary.json', ledger['summary_sha256'])
    evidence = [(row, _pinned(root/row['id']/'evidence.json', row['sha256'])) for row in ledger['entries']]
    # The independent existing reader rebinds pinned sources, all chemistry,
    # requests/masks/coordinates, template snapshots/global L and native scores.
    summary = comparison.readback_controls(root)
    if [r['id'] for r in summary['controls']] != [r['id'] for r in ledger['entries']]:
        raise ValueError('closure-diagnosis-ledger-denominator')
    return evidence, summary


def diagnose_cell(row, data):
    saved = parse_research_json(data)
    def dto(key, cls):
        return decode_research_record(json.dumps(saved[key], allow_nan=False).encode(), cls)
    probe = dto('probe', ObservedPeptideProbe)
    original = probe.request
    common = dto('request', ReconstructionRequest)
    plan = dto('baseline_plan', GeometryPlan)
    if plan != compile_peptide_request_plan(original, plan.policy):
        raise ValueError('closure-diagnosis-recompiled-original-plan')
    common_plan = GeometryPlan(common.digest, plan.strategy, plan.policy)
    baseline = dto('baseline', ProposalResult)
    template = dto('template', SeededHeavyBaseline)
    projection = dto('projection', ProjectionResult) if 'projection' in saved else None
    reference = tuple(decode_research_record(json.dumps(o, allow_nan=False).encode(), Observation) for o in saved['postproposal_reference'])
    phases = {}
    for phase in ('raw', 'rounded'):
        candidate = getattr(baseline, phase)
        if candidate is None:
            raise ValueError('closure-diagnosis-needs-frozen-complete-original')
        own = measure_terms(original, plan, candidate.additions)
        if own.audit_findings != getattr(baseline, phase+'_audit').findings:
            raise ValueError('closure-diagnosis-original-audit-binding')
        frames, omitted = trace_frames(original, plan, candidate.additions)
        tables = {'original-own': own, 'original-common': measure_terms(common, common_plan, candidate.additions)}
        supplied = getattr(template, phase)
        if supplied:
            tables.update({'template-own': measure_terms(original, plan, supplied),
                           'template-common': measure_terms(common, common_plan, supplied)})
        if projection is not None and getattr(projection, phase) is not None:
            tables['projection-common'] = measure_terms(common, common_plan, getattr(projection, phase).additions)
        phases[phase] = {'original_candidate_digest': digest_record(candidate), 'tables': tables,
                         'original_frames': frames, 'original_unenforced_length_edges': omitted}
    return {'id': row['id'], 'input_sha256': row['sha256'], 'original_request_digest': original.digest,
            'common_request_digest': common.digest, 'original_plan_digest': plan.digest,
            'old_control': saved['control'], 'phases': phases,
            'native_diagnostic_only': {'own': measure_terms(original, plan, reference), 'common': measure_terms(common, common_plan, reference)},
            'fixed_neighbour_necessary_support': fixed_neighbour_support(original, plan.policy)}


def build_report():
    inputs, summary = qualified_inputs()
    records = [diagnose_cell(row, data) for row, data in inputs]
    return {'schema': 'dvbfixer.chemical-closure-diagnosis.v1', 'protocol_sha256': PROTOCOL_SHA256,
        'ledger_sha256': LEDGER_SHA256, 'declared_saved_cells': 9,
        'historical_projection_calls': sum(r['projection_called'] for r in summary['controls']),
        'historical_original_rejections': [r['id'] for r in summary['controls'] if r['baseline_status'] == 'rejected'],
        'new_coordinate_proposals': 0, 'old_results_changed': False, 'records': records,
        'necessary_support_is_not_complete_feasibility': True, 'production_enabled': False}


def run_diagnosis(output: Path):
    report = build_report()  # Qualify ALL inputs before reserving new output.
    output.mkdir(parents=True, exist_ok=False)
    _write(output/'diagnosis.json', report)
    return report


def readback_diagnosis(output: Path):
    report = build_report()
    with (output/'diagnosis.json').open('rb') as handle:
        saved = parse_research_json(handle.read(64*1024*1024+1))
    if saved != _json(report):
        raise ValueError('closure-diagnosis-saved-report-mismatch')
    return {'schema': 'dvbfixer.chemical-closure-diagnosis-readback.v1', 'saved_cells': 9,
            'new_coordinate_proposals': 0, 'report_sha256': hashlib.sha256((output/'diagnosis.json').read_bytes()).hexdigest()}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    parser.add_argument('--readback', action='store_true')
    args = parser.parse_args()
    if args.readback:
        print(json.dumps(readback_diagnosis(args.output), indent=2))
    else:
        report = run_diagnosis(args.output)
        print('diagnosed', len(report['records']), 'saved cells; no new proposals')
