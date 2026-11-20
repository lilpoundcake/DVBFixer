"""Frozen WP2N projected-parent H continuation; never solve/project/search heavy."""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from dataclasses import replace
from pathlib import Path

from scripts.chemical_reconstruction_all20_controls import FFBONDED
from scripts.chemical_reconstruction_ccd_controls import ROOT, _json, _write
from scripts.chemical_reconstruction_projection_controls import check_pins as check_projection_pins
from scripts.chemical_reconstruction_staged_controls import _read_bounded, _read_pinned

from dvbfixer.chemical_reconstruction.amide_authority import FFBONDED_SHA256, load_amide_authority
from dvbfixer.chemical_reconstruction.contracts import ReconstructionRequest, digest_record
from dvbfixer.chemical_reconstruction.declared_peptide import DeclaredPeptideAuthority
from dvbfixer.chemical_reconstruction.observed_peptide import ObservedPeptideProbe
from dvbfixer.chemical_reconstruction.planar_port import PlanarPortAuthority, PlanarPortDeclaration
from dvbfixer.chemical_reconstruction.projected_hydrogen_audit import (
    audit_projected_hydrogen_result,
)
from dvbfixer.chemical_reconstruction.projected_hydrogen_contracts import (
    ProjectedHeavyParent,
    ProjectedHydrogenResult,
)
from dvbfixer.chemical_reconstruction.projected_hydrogens import complete_projected_hydrogens
from dvbfixer.chemical_reconstruction.projection_contracts import ProjectionResult
from dvbfixer.chemical_reconstruction.research_decode import (
    decode_research_record,
    parse_research_json,
)
from dvbfixer.chemical_reconstruction.source_preparation import (
    audit_source_transport,
    prepare_source,
)

PROTOCOL = ROOT/'docs/plans/chemical-reconstruction-projected-h-v1.md'
PROTOCOL_SHA256 = '9e65a418c8142c521e1688c276a00084bf24b591d6f3ef984115ddb26eda9944'
LEDGER = ROOT/'docs/plans/chemical-reconstruction-projected-h-parents-v1.json'
LEDGER_SHA256 = '9f6c6690912c1bfd5e074b9005ec3627d24b908cccbed0e10b8e36cd0fa5c69b'


def check_pins():
    check_projection_pins()
    _read_pinned(PROTOCOL, PROTOCOL_SHA256)
    return json.loads(_read_pinned(LEDGER, LEDGER_SHA256))


def decode_parent(data: bytes, row):
    saved = parse_research_json(data)
    if not isinstance(saved, dict):
        raise ValueError('projected-H-parent-needs-object')
    def dto(key, cls):
        return decode_research_record(json.dumps(saved[key], allow_nan=False).encode(), cls)
    # These are the ONLY fields entering the scientific parent; native metrics stay here.
    parent = ProjectedHeavyParent(dto('original_request', ReconstructionRequest), dto('request', ReconstructionRequest),
        dto('result', ProjectionResult), dto('authority', DeclaredPeptideAuthority),
        dto('port_declaration', PlanarPortDeclaration), dto('parameters', PlanarPortAuthority))
    probe = dto('probe', ObservedPeptideProbe)
    if (parent.original.digest != probe.request.digest or digest_record(parent.projection) != row['projection_digest']
            or parent.projection.rounded is None
            or digest_record(parent.projection.rounded.additions) != row['rounded_additions_sha256']):
        raise ValueError('projected-H-frozen-source-projection-checkpoint-binding')
    return parent, probe


def verify_source(directory: Path, parent_saved, parent: ProjectedHeavyParent, probe: ObservedPeptideProbe):
    original = (directory/'original.pdb').read_bytes()
    control = parent_saved['control']
    if (hashlib.sha256(original).hexdigest() != control['original_sha256']
            or original != (ROOT/control['source']).read_bytes()):
        raise ValueError('projected-H-original-source-byte-binding')
    preparation = prepare_source(original)
    if (preparation.prepared_text.encode('ascii') != (directory/'prepared.pdb').read_bytes()
            or not audit_source_transport(preparation, probe).passed
            or not audit_source_transport(preparation, replace(probe, request=parent.request)).passed):
        raise ValueError('projected-H-source-transport-refused')


def run_controls(output: Path):
    ledger = check_pins()
    # Hash ALL three dispositions before reserving output or making an H proposal.
    parents = [(row, ROOT/row['parent_artifact'], _read_pinned(ROOT/row['parent_artifact'], row['parent_sha256']))
               for row in ledger['entries']]
    output.mkdir(parents=True, exist_ok=False)
    parameters = load_amide_authority(FFBONDED, expected_sha256=FFBONDED_SHA256)
    records = []
    for row, parent_path, data in parents:
        directory = output/row['id']
        directory.mkdir()
        (directory/'parent-evidence.json').write_bytes(data)
        for name in ('original.pdb', 'prepared.pdb'):
            if (parent_path.parent/name).is_file():
                shutil.copyfile(parent_path.parent/name, directory/name)
        item = {'id': row['id'], 'parent_status': row['status'], 'parent_sha256': row['parent_sha256'], 'H_stage_called': False}
        evidence = {'control': item, 'frozen_parent': row}
        if row['rounded_additions_sha256'] is None:
            item.update(status='parent-unavailable', reason=row['reason'])
        else:
            parent, probe = decode_parent(data, row)
            saved = parse_research_json(data)
            verify_source(directory, saved, parent, probe)
            result = complete_projected_hydrogens(parent, parameters, backbone_types=None)
            findings = audit_projected_hydrogen_result(parent, parameters, None, result)
            if findings:
                raise ValueError('projected-H-result-readback:'+','.join(findings))
            item.update(H_stage_called=True, status=result.status, reason=result.reason,
                        combined_passed=bool(result.audit and result.audit.passed))
            evidence.update(parent=parent, H_parameters=parameters, explicit_backbone_types=None, result=result,
                            independent_readback_findings=findings)
        _write(directory/'evidence.json', evidence)
        records.append(item)
    check_pins()
    summary = {'schema': 'dvbfixer.projected-parent-H-controls.v1', 'protocol_sha256': PROTOCOL_SHA256,
               'parent_ledger_sha256': LEDGER_SHA256, 'declared_entries': 3, 'controls': records,
               'new_heavy_solves_searches_or_samples': 0, 'checkpoint_movement_or_reselection': False,
               'native_H_or_production_qualified': False}
    _write(output/'summary.json', summary)
    return summary


def readback_controls(output: Path):
    ledger = check_pins()
    parameters = load_amide_authority(FFBONDED, expected_sha256=FFBONDED_SHA256)
    records = []
    saved_controls = []
    for row in ledger['entries']:
        directory = output/row['id']
        data = _read_pinned(directory/'parent-evidence.json', row['parent_sha256'])
        evidence = parse_research_json(_read_bounded(directory/'evidence.json'))
        if not isinstance(evidence, dict) or evidence['frozen_parent'] != row:
            raise ValueError('projected-H-ledger-evidence-binding')
        control = evidence['control']
        saved_controls.append(control)
        if (control['id'] != row['id'] or control['parent_status'] != row['status']
                or control['parent_sha256'] != row['parent_sha256']):
            raise ValueError('projected-H-parent-control-binding')
        saved_parent = parse_research_json(data)
        original = (directory/'original.pdb').read_bytes()
        if (hashlib.sha256(original).hexdigest() != saved_parent['control']['original_sha256']
                or original != (ROOT/saved_parent['control']['source']).read_bytes()):
            raise ValueError('projected-H-saved-original-source-binding')
        if row['rounded_additions_sha256'] is None:
            if ('result' in evidence or control['H_stage_called'] is not False or control['status'] != 'parent-unavailable'
                    or control['reason'] != row['reason']):
                raise ValueError('projected-H-unavailable-parent-disposition')
            try:
                preparation = prepare_source(original)
            except ValueError as exc:
                if str(exc) != row['reason'] or (directory/'prepared.pdb').exists():
                    raise ValueError('projected-H-refused-source-preparation-mismatch') from exc
            else:
                if (directory/'prepared.pdb').read_bytes() != preparation.prepared_text.encode('ascii'):
                    raise ValueError('projected-H-unavailable-parent-preparation-transport')
            records.append({'id': row['id'], 'status': control['status'], 'reason': control['reason'], 'H_stage_called': False})
            continue
        parent, probe = decode_parent(data, row)
        verify_source(directory, saved_parent, parent, probe)
        typed_parent = decode_research_record(json.dumps(evidence['parent'], allow_nan=False).encode(), ProjectedHeavyParent)
        if typed_parent.digest != parent.digest:
            raise ValueError('projected-H-parent-kind-or-binding')
        result = decode_research_record(json.dumps(evidence['result'], allow_nan=False).encode(), ProjectedHydrogenResult)
        findings = audit_projected_hydrogen_result(parent, parameters, None, result)
        if findings:
            raise ValueError('projected-H-independent-readback:'+','.join(findings))
        if (evidence['H_parameters'] != _json(parameters) or evidence['explicit_backbone_types'] is not None
                or control['status'] != result.status or control['reason'] != result.reason or control['H_stage_called'] is not True
                or control['combined_passed'] is not bool(result.audit and result.audit.passed)):
            raise ValueError('projected-H-parameter-or-control-mismatch')
        records.append({'id': row['id'], 'status': result.status, 'H_stage_called': True, 'parent_digest': parent.digest,
                        'rounded_heavy_signature': row['rounded_additions_sha256'],
                        'raw_findings': result.audit.findings if result.audit else (result.reason,)})
    summary = parse_research_json(_read_bounded(output/'summary.json'))
    expected_summary = {'schema': 'dvbfixer.projected-parent-H-controls.v1', 'protocol_sha256': PROTOCOL_SHA256,
        'parent_ledger_sha256': LEDGER_SHA256, 'declared_entries': 3, 'controls': saved_controls,
        'new_heavy_solves_searches_or_samples': 0, 'checkpoint_movement_or_reselection': False,
        'native_H_or_production_qualified': False}
    if json.dumps(summary, sort_keys=True, allow_nan=False) != json.dumps(expected_summary, sort_keys=True, allow_nan=False):
        raise ValueError('projected-H-summary-binding')
    return {'schema': 'dvbfixer.projected-parent-H-readback.v1', 'controls': records}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    parser.add_argument('--readback', action='store_true')
    args = parser.parse_args()
    print(json.dumps(readback_controls(args.output) if args.readback else run_controls(args.output), indent=2))
