"""Frozen graph-grid/steric comparison; native coordinates are postselection only."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from scripts.chemical_reconstruction_all20_controls import (
    HEAVY_POLICY,
    PAIRS,
    heldout_after_proposal,
    observed_input,
)
from scripts.chemical_reconstruction_all20_controls import check_pins as check_all20_pins
from scripts.chemical_reconstruction_ccd_controls import ROOT, _write

from dvbfixer.chemical_reconstruction.conformer_search import audit_search_result, search_conformers
from dvbfixer.chemical_reconstruction.contracts import digest_record
from dvbfixer.chemical_reconstruction.heldout_scoring import (
    HoldoutScoringError,
    score_heavy_holdout,
)
from dvbfixer.chemical_reconstruction.observed_peptide import compile_peptide_request_plan
from dvbfixer.chemical_reconstruction.steric_ranking import radii_digest

PROTOCOL = ROOT / 'docs/plans/chemical-reconstruction-conformer-grid-v1.md'
PROTOCOL_SHA256 = '205fdfdbb1a21b2db3ae734d83e3cd760d6ed16fc8ecf6eb2b1c86ea22142536'
LEDGER = ROOT / 'docs/plans/chemical-reconstruction-conformer-baselines-v1.json'
LEDGER_SHA256 = '8941e430a8f9c44076e13a6c97a730fcbc59282b562824572bdd9a93361a5737'


def check_pins():
    check_all20_pins()
    for path, pin in ((PROTOCOL, PROTOCOL_SHA256), (LEDGER, LEDGER_SHA256)):
        if hashlib.sha256(path.read_bytes()).hexdigest() != pin:
            raise ValueError(f'frozen conformer protocol/ledger changed: {path}')
    ledger = json.loads(LEDGER.read_text())
    if hashlib.sha256((ROOT/ledger['radii_owner']).read_bytes()).hexdigest() != ledger['radii_owner_sha256']:
        raise ValueError('frozen VDW radii owner changed')
    if radii_digest() != ledger['radii_values_digest']:
        raise ValueError('frozen VDW radii values changed')


def compare_after_selection(request, base_plan, authority, result, directory, evidence):
    # No reference extraction or heldout scoring until all search/selection ends.
    reference = heldout_after_proposal('observed', authority, directory)
    evidence['heldout_reference_after_selection'] = reference
    metrics = {'baseline': {}, 'selected': {}, 'attempts': []}
    for kind in ('raw', 'rounded'):
        candidate = getattr(result.baseline, kind)
        if candidate is not None:
            metrics['baseline'][kind] = score_heavy_holdout(request, base_plan, candidate, reference)
    if result.plan is not None:
        for attempt in result.attempts:
            record = {'index': attempt.index}
            for kind in ('raw', 'rounded'):
                candidate = getattr(attempt, kind)
                if candidate is None:
                    continue
                try:
                    score = score_heavy_holdout(request, result.plan, candidate, reference)
                except HoldoutScoringError as exc:
                    record[kind] = {'status': 'not-scored', 'reason': str(exc)}
                else:
                    record[kind] = {'status': 'measured-not-a-ranking-input', 'score': score}
                    if attempt.index == result.selected_index:
                        metrics['selected'][kind] = score
            metrics['attempts'].append(record)
    evidence['postselection_native_metrics'] = metrics


def run_controls(output: Path):
    check_pins()
    output.mkdir(parents=True, exist_ok=False)
    baselines = {r['source_component']: r for r in json.loads(LEDGER.read_text())['entries']}
    records = []
    for source, state in PAIRS:
        directory = output/f'observed-{source}-{state}'
        directory.mkdir()
        item = {'id': directory.name, 'source_component': source, 'requested_state': state, 'search_called': False}
        evidence = {'control': item, 'frozen_baseline': baselines[source]}
        try:
            request, authority, _incoming = observed_input(source, state, 'heavy', directory, evidence)
            evidence['request'] = request
            if evidence['original_sha256'] != baselines[source]['original_source_sha256']:
                raise ValueError('frozen-source-byte-binding')
            if not request.ownership.generated:
                if baselines[source]['status'] != 'not-applicable-no-sidechain-heavy-atoms':
                    raise ValueError('unexpected-empty-heavy-mask')
                item.update(status='not-applicable-no-sidechain-heavy-atoms', reason='GLY-empty-sidechain-mask')
            else:
                if (authority.identity.sequence_number != baselines[source]['target']
                        or sorted(a.name for a in request.ownership.generated) != baselines[source]['heavy_mask']):
                    raise ValueError('frozen-target-or-heavy-mask-binding')
                base_plan = compile_peptide_request_plan(request, HEAVY_POLICY)
                evidence['base_plan'] = base_plan
                roots = tuple(sorted(o.identity for o in request.observations if o.identity.component == authority.identity))
                item['search_called'] = True
                result = search_conformers(request, base_plan, authority, roots,
                                           expected_baseline_signature=baselines[source]['rounded_additions_sha256'])
                evidence['result'] = result
                item.update(status=result.status, reason=result.reason, target=authority.identity.sequence_number,
                            attempts=len(result.attempts), selected_index=result.selected_index,
                            selectable_attempts=sum(a.selectable for a in result.attempts))
                if result.plan is not None:
                    audit = audit_search_result(request, authority, result)
                    evidence['search_measurement_audit'] = audit
                    if audit:
                        raise ValueError('saved-search-measurement-audit-failed:'+','.join(audit))
                    item['eligible_axes'] = len(result.plan.topology.rotors)
                    assert digest_record(result.baseline.rounded.additions) == baselines[source]['rounded_additions_sha256']
                compare_after_selection(request, base_plan, authority, result, directory, evidence)
        except ValueError as exc:
            item.update(status='refused', reason=str(exc))
        _write(directory/'evidence.json', evidence)
        records.append(item)
    check_pins()
    summary = {'schema': 'dvbfixer.fragment-torsion-grid-controls.v1', 'protocol_sha256': PROTOCOL_SHA256,
               'baseline_ledger_sha256': LEDGER_SHA256, 'declared': 20, 'executed_entries': len(records), 'controls': records,
               'empirical_rotamer_native_H_or_production_qualified': False, 'production_enabled': False}
    _write(output/'summary.json', summary)
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    print(json.dumps(run_controls(parser.parse_args().output), indent=2))
