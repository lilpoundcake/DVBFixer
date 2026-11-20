"""Frozen WP2L saved-parent continuation; no heavy search or native-score input."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from scripts.chemical_reconstruction_all20_controls import FFBONDED, backbone_types
from scripts.chemical_reconstruction_all20_controls import check_pins as check_all20_pins
from scripts.chemical_reconstruction_ccd_controls import ROOT, _json, _write

from dvbfixer.chemical_reconstruction.amide_authority import FFBONDED_SHA256, load_amide_authority
from dvbfixer.chemical_reconstruction.conformer_search import ConformerSearchResult
from dvbfixer.chemical_reconstruction.contracts import ReconstructionRequest, digest_record
from dvbfixer.chemical_reconstruction.declared_peptide import DeclaredPeptideAuthority
from dvbfixer.chemical_reconstruction.observed_peptide import ObservedPeptideProbe
from dvbfixer.chemical_reconstruction.research_decode import (
    MAX_BYTES,
    decode_research_record,
    parse_research_json,
)
from dvbfixer.chemical_reconstruction.source_preparation import (
    audit_source_transport,
    prepare_source,
)
from dvbfixer.chemical_reconstruction.staged_audit import audit_staged_result
from dvbfixer.chemical_reconstruction.staged_completion import complete_staged
from dvbfixer.chemical_reconstruction.staged_contracts import StagedResult

PROTOCOL = ROOT/'docs/plans/chemical-reconstruction-heavy-h-v1.md'
PROTOCOL_SHA256 = 'ebcc9f8f87703b6c131f82b761a21dccb1438c537b08ef23434d2d68a3eda3af'
LEDGER = ROOT/'docs/plans/chemical-reconstruction-heavy-h-parents-v1.json'
LEDGER_SHA256 = '834756ca0c0d3477872a20979071e25bf1fcef926b961c6d09977a2388aa4af7'


def _read_bounded(path: Path) -> bytes:
    if path.stat().st_size > MAX_BYTES:
        raise ValueError('staged-parent-file-byte-budget')
    with path.open('rb') as handle:
        data = handle.read(MAX_BYTES + 1)
    if len(data) > MAX_BYTES:
        raise ValueError('staged-parent-file-byte-budget')
    return data


def _read_pinned(path: Path, pin: str) -> bytes:
    data = _read_bounded(path)
    if hashlib.sha256(data).hexdigest() != pin:
        raise ValueError('staged-parent-artifact-hash-mismatch')
    return data


def check_pins():
    check_all20_pins()
    _read_pinned(PROTOCOL, PROTOCOL_SHA256)
    return json.loads(_read_pinned(LEDGER, LEDGER_SHA256))


def decode_parent(data: bytes):
    # Top-level retrospective/native fields deliberately never leave this I/O adapter.
    saved = parse_research_json(data)
    if not isinstance(saved, dict):
        raise ValueError('staged-parent-needs-object')
    def dto(key, expected):
        return decode_research_record(json.dumps(saved[key], allow_nan=False).encode(), expected)
    request = dto('request', ReconstructionRequest)
    authority = dto('authority', DeclaredPeptideAuthority)
    probe = dto('probe', ObservedPeptideProbe)
    parent = dto('result', ConformerSearchResult)
    if digest_record(request) != digest_record(probe.request) or not probe.real_peptide_partners:
        raise ValueError('staged-original-probe-request-binding')
    return request, authority, probe, parent


def verify_source(parent_path, saved, probe):
    original = (parent_path.parent/'original.pdb').read_bytes()
    if hashlib.sha256(original).hexdigest() != saved['original_sha256']:
        raise ValueError('staged-original-source-bytes-mismatch')
    source = ROOT/saved['source']
    if source.read_bytes() != original:
        raise ValueError('staged-tracked-source-byte-binding')
    prepared = prepare_source(original)
    if (prepared.prepared_text.encode('ascii') != (parent_path.parent/'prepared.pdb').read_bytes()
            or not audit_source_transport(prepared, probe).passed):
        raise ValueError('staged-original-source-transport-refused')
    return original, prepared.prepared_text.encode('ascii')


def run_controls(output: Path):
    ledger = check_pins()
    # All twenty saved parents must exist and match BEFORE any stage or reservation.
    parents = [(row, ROOT/row['parent_artifact'], _read_pinned(ROOT/row['parent_artifact'], row['parent_sha256']))
               for row in ledger['entries']]
    output.mkdir(parents=True, exist_ok=False)
    parameters = load_amide_authority(FFBONDED, expected_sha256=FFBONDED_SHA256)
    records = []
    for row, parent_path, data in parents:
        directory = output/row['id']
        directory.mkdir()
        (directory/'parent-evidence.json').write_bytes(data)
        item = {'id': row['id'], 'parent_status': row['status'], 'parent_sha256': row['parent_sha256'],
                'selected_index': row['selected_index'], 'stage_called': False}
        evidence = {'control': item, 'frozen_parent': row}
        if row['selected_index'] is None:
            item.update(status='not-applicable' if 'GLY' in row['id'] else 'parent-unavailable',
                        reason='unchanged-frozen-parent-disposition-no-H-only-fallback')
        else:
            request, authority, probe, parent = decode_parent(data)
            if parent.selected_index != row['selected_index']:
                raise ValueError('staged-frozen-parent-selection-mismatch')
            checkpoint = parent.attempts[parent.selected_index].rounded
            if checkpoint is None or digest_record(checkpoint.additions) != row['rounded_additions_sha256']:
                raise ValueError('staged-frozen-rounded-heavy-signature-mismatch')
            saved = parse_research_json(data)
            original, prepared = verify_source(parent_path, saved, probe)
            (directory/'original.pdb').write_bytes(original)
            (directory/'prepared.pdb').write_bytes(prepared)
            types = None if authority.inventory.state_name == 'PRO' else backbone_types(authority.inventory.state_name)
            incoming = probe.real_peptide_partners[0]
            result = complete_staged(request, parent, authority, incoming, parameters, backbone_types=types)
            findings = audit_staged_result(request, parent, authority, incoming, parameters, types, result)
            if findings:
                raise ValueError('staged-saved-result-integrity:'+','.join(findings))
            item.update(stage_called=True, status=result.status, reason=result.reason,
                        combined_passed=bool(result.audit and result.audit.passed))
            evidence.update(result=result, explicit_backbone_types=types, parameters=parameters,
                            independent_readback_findings=findings)
        _write(directory/'evidence.json', evidence)
        records.append(item)
    summary = {'schema': 'dvbfixer.staged-heavy-H-controls.v1', 'protocol_sha256': PROTOCOL_SHA256,
               'parent_ledger_sha256': LEDGER_SHA256, 'declared': 20, 'entries': records,
               'new_heavy_searches': 0, 'heavy_reselection_or_movement': False,
               'native_H_or_production_qualified': False}
    check_pins()
    _write(output/'summary.json', summary)
    return summary


def readback_controls(output: Path):
    ledger = check_pins()
    parameters = load_amide_authority(FFBONDED, expected_sha256=FFBONDED_SHA256)
    records = []
    for row in ledger['entries']:
        directory = output/row['id']
        data = _read_pinned(directory/'parent-evidence.json', row['parent_sha256'])
        saved = parse_research_json(_read_bounded(directory/'evidence.json'))
        if not isinstance(saved, dict):
            raise ValueError('staged-saved-evidence-needs-object')
        if saved['frozen_parent'] != row:
            raise ValueError('staged-saved-parent-ledger-binding')
        if row['selected_index'] is None:
            if 'result' in saved or saved['control']['stage_called']:
                raise ValueError('staged-unexecuted-parent-disposition')
            records.append({'id': row['id'], 'status': saved['control']['status'], 'stage_called': False})
            continue
        request, authority, probe, parent = decode_parent(data)
        original, prepared = verify_source(directory/'parent-evidence.json', parse_research_json(data), probe)
        if (directory/'original.pdb').read_bytes() != original or (directory/'prepared.pdb').read_bytes() != prepared:
            raise ValueError('staged-saved-source-transport')
        types = None if authority.inventory.state_name == 'PRO' else backbone_types(authority.inventory.state_name)
        result = decode_research_record(json.dumps(saved['result'], allow_nan=False).encode(), StagedResult)
        findings = audit_staged_result(request, parent, authority, probe.real_peptide_partners[0], parameters, types, result)
        if findings:
            raise ValueError('staged-independent-result-integrity:'+','.join(findings))
        checkpoint = parent.attempts[row['selected_index']].rounded
        if checkpoint is None or digest_record(checkpoint.additions) != row['rounded_additions_sha256']:
            raise ValueError('staged-independent-checkpoint-signature')
        if (saved['control']['status'] != result.status or not saved['control']['stage_called']
                or saved['control']['selected_index'] != row['selected_index']
                or saved['explicit_backbone_types'] != (list(types) if types is not None else None)
                or saved['parameters'] != _json(parameters)):
            raise ValueError('staged-control-parameter-or-disposition-mismatch')
        records.append({'id': row['id'], 'status': result.status, 'stage_called': True,
                        'raw_findings': result.audit.findings if result.audit else (result.reason,),
                        'heavy_checkpoint_signature': row['rounded_additions_sha256']})
    return {'schema': 'dvbfixer.staged-heavy-H-readback.v1', 'entries': records}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    parser.add_argument('--readback', action='store_true')
    arguments = parser.parse_args()
    print(json.dumps(readback_controls(arguments.output) if arguments.readback else run_controls(arguments.output), indent=2))
