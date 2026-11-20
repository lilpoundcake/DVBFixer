"""Frozen WP2O type-proven port continuation; native references only after heavy/H."""
from __future__ import annotations

import argparse
import hashlib
import json
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
from scripts.chemical_reconstruction_hydrogen_controls import RTP
from scripts.chemical_reconstruction_projection_controls import SOURCES
from scripts.chemical_reconstruction_projection_controls import check_pins as check_old_pins

from dvbfixer.chemical_reconstruction.amide_authority import FFBONDED_SHA256, load_amide_authority
from dvbfixer.chemical_reconstruction.construction import propose
from dvbfixer.chemical_reconstruction.heldout_scoring import score_heavy_holdout
from dvbfixer.chemical_reconstruction.incoming_source import (
    audit_incoming_source,
    bind_incoming_source,
)
from dvbfixer.chemical_reconstruction.incoming_type import (
    IncomingCarbonylTypeProof,
    candidate_family,
    load_incoming_type_catalog,
    prove_incoming_type,
)
from dvbfixer.chemical_reconstruction.metric_projection import propose_projection
from dvbfixer.chemical_reconstruction.observed_peptide import (
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
from dvbfixer.chemical_reconstruction.projection_contracts import append_projection_constraints
from dvbfixer.chemical_reconstruction.research_decode import (
    decode_research_record,
    parse_research_json,
)
from dvbfixer.chemical_reconstruction.source_preparation import (
    audit_source_transport,
    prepare_source,
    prepare_source_probe,
)
from dvbfixer.chemical_reconstruction.type_proven_port import TypeProvenPlanarPortDeclaration
from dvbfixer.domain.structure_identity import ExactAtomRef

PROTOCOL = ROOT/'docs/plans/chemical-reconstruction-incoming-type-v1.md'
PROTOCOL_SHA256 = 'fb2326ff975b631a569679d23b39944a26c2a64f217f9f03cd39221b655a92b7'


def check_pins():
    check_old_pins()
    if hashlib.sha256(PROTOCOL.read_bytes()).hexdigest() != PROTOCOL_SHA256:
        raise ValueError('incoming-type-protocol-pin')
    return load_incoming_type_catalog(RTP)


def scores_after_proposals(directory, parent):
    reference = heldout_after_proposal('observed', parent.peptide_authority, directory)
    projection = parent.projection
    scores = {}
    for arm, request, plan, result in (
        ('baseline', parent.original, projection.plan.baseline_plan, projection.plan.baseline),
        ('projection', parent.request, projection.plan, projection),
    ):
        for phase in ('raw', 'rounded'):
            candidate = getattr(result, phase)
            if candidate is not None:
                scores[f'{arm}-{phase}'] = score_heavy_holdout(request, plan, candidate, reference)
    return reference, scores


def run_controls(output: Path):
    catalog = check_pins()
    port_parameters = load_planar_port_authority(FFBONDED)
    H_parameters = load_amide_authority(FFBONDED, expected_sha256=FFBONDED_SHA256)
    output.mkdir(parents=True, exist_ok=False)
    records = []
    for name, relative, pin in SOURCES:
        directory = output/name
        directory.mkdir()
        item = {'id': name, 'source': relative, 'original_sha256': pin,
                'baseline_called': False, 'projection_called': False, 'H_called': False}
        evidence = {'control': item}
        try:
            data = (ROOT/relative).read_bytes()
            (directory/'original.pdb').write_bytes(data)
            preparation = prepare_source(data)
            (directory/'prepared.pdb').write_bytes(preparation.prepared_text.encode('ascii'))
            identity = select_interior_target(data, 'PRO', chain='A')
            authority = authority_for('PRO', 'PRO', identity, evidence)
            probe = prepare_source_probe(preparation, authority, withheld=frozenset({'CB', 'CG', 'CD'}))
            source = bind_incoming_source(preparation, probe, authority)
            proof = prove_incoming_type(probe.request, source, catalog, load_parent(candidate_family(source.residue_name, catalog)[0]))
            if audit_incoming_source(preparation, probe, authority, source):
                raise ValueError('incoming-source-independent-admission-failed')
            evidence.update(probe=probe, proof=proof)
            original = probe.request
            baseline_plan = compile_peptide_request_plan(original, HEAVY_POLICY)
            baseline = propose(original, baseline_plan)
            item.update(target=identity.sequence_number, baseline_called=True, baseline_status=baseline.status,
                        incoming_label=source.residue_name, candidate_states=proof.candidate_states)
            evidence.update(baseline_plan=baseline_plan, baseline=baseline)
            if baseline.status != 'declared-geometry-passed':
                raise ValueError('unchanged-baseline-not-raw-rounded-qualified')
            pro_types = {r[0]: r[1] for s in catalog.states if s.name == 'PRO' for r in s.atoms}
            n, ca, cd = (ExactAtomRef(identity, atom) for atom in ('N', 'CA', 'CD'))
            declaration = TypeProvenPlanarPortDeclaration(source.incoming, n, ca, cd,
                (proof.atom_type, pro_types['N'], pro_types['CA'], pro_types['CD']),
                authority.digest, port_parameters.digest, incoming_proof=proof)
            request = append_projection_constraints(original, planar_port_constraints(original, authority, declaration, port_parameters))
            if not audit_source_transport(preparation, replace(probe, request=request)).passed:
                raise ValueError('type-proven-port-source-transport-failed')
            anchors = ca, n, ExactAtomRef(identity, 'C')
            projection = propose_projection(original, request, baseline_plan, baseline, anchors)
            parent = ProjectedHeavyParent(original, request, projection, authority, declaration, port_parameters)
            evidence['parent'] = parent
            item.update(projection_called=True, projection_status=projection.status)
            if audit_projection_result(original, request, baseline_plan, baseline, anchors, projection):
                raise ValueError('type-proven-projection-readback-failed')
            if projection.status == 'projected-geometry-passed':
                result = complete_projected_hydrogens(parent, H_parameters, backbone_types=None)
                evidence.update(H_result=result, H_parameters=H_parameters)
                item.update(H_called=True, status=result.status, reason=result.reason)
                if audit_projected_hydrogen_result(parent, H_parameters, None, result):
                    raise ValueError('type-proven-H-readback-failed')
            else:
                item.update(status=projection.status, reason=projection.reason)
            # Never expose withheld coordinates to either proposal phase.
            reference, scores = scores_after_proposals(directory, parent)
            evidence.update(postproposal_reference=reference, postproposal_native_scores=scores)
        except ValueError as exc:
            item.update(status='refused', reason=str(exc))
        _write(directory/'evidence.json', evidence)
        records.append(item)
    check_pins()
    summary = {'schema': 'dvbfixer.incoming-type-controls.v1', 'protocol_sha256': PROTOCOL_SHA256,
               'declared': 3, 'controls': records, 'charge_or_state_assigned': False,
               'public_or_production_qualified': False}
    _write(output/'summary.json', summary)
    return summary


def readback_controls(output: Path):
    catalog = check_pins()
    port_parameters = load_planar_port_authority(FFBONDED)
    H_parameters = load_amide_authority(FFBONDED, expected_sha256=FFBONDED_SHA256)
    controls = []
    for name, relative, pin in SOURCES:
        directory = output/name
        data = (directory/'original.pdb').read_bytes()
        if data != (ROOT/relative).read_bytes() or hashlib.sha256(data).hexdigest() != pin:
            raise ValueError('incoming-readback-source-pin')
        saved = parse_research_json((directory/'evidence.json').read_bytes())
        def dto(key, cls, saved=saved):
            return decode_research_record(json.dumps(saved[key], allow_nan=False).encode(), cls)
        item = saved['control']
        if (item['id'], item['source'], item['original_sha256']) != (name, relative, pin):
            raise ValueError('incoming-readback-control-source')
        controls.append(item)
        try:
            preparation = prepare_source(data)
        except ValueError as exc:
            if (item != {'id': name, 'source': relative, 'original_sha256': pin, 'baseline_called': False,
                         'projection_called': False, 'H_called': False, 'status': 'refused', 'reason': str(exc)}
                    or set(saved) != {'control'} or (directory/'prepared.pdb').exists()):
                raise ValueError('incoming-readback-source-refusal') from exc
            continue
        if (directory/'prepared.pdb').read_bytes() != preparation.prepared_text.encode('ascii'):
            raise ValueError('incoming-readback-preparation')
        parent = dto('parent', ProjectedHeavyParent)
        proof = dto('proof', IncomingCarbonylTypeProof)
        identity = select_interior_target(data, 'PRO', chain='A')
        authority = authority_for('PRO', 'PRO', identity, {})
        probe = prepare_source_probe(preparation, authority, withheld=frozenset({'CB', 'CG', 'CD'}))
        source = bind_incoming_source(preparation, probe, authority)
        expected = prove_incoming_type(probe.request, source, catalog, load_parent(candidate_family(source.residue_name, catalog)[0]))
        if (proof.digest != expected.digest or audit_incoming_source(preparation, probe, authority, proof.source)
                or _json(probe) != saved['probe'] or parent.original.digest != probe.request.digest
                or parent.peptide_authority.digest != authority.digest or parent.port_parameters != port_parameters
                or not isinstance(parent.port_declaration, TypeProvenPlanarPortDeclaration)
                or parent.port_declaration.incoming_proof.digest != proof.digest):
            raise ValueError('incoming-readback-pinned-authority-source-proof')
        request = append_projection_constraints(probe.request, planar_port_constraints(probe.request, authority, parent.port_declaration, port_parameters))
        projection = parent.projection
        plan = compile_peptide_request_plan(probe.request, HEAVY_POLICY)
        anchors = tuple(ExactAtomRef(identity, atom) for atom in ('CA', 'N', 'C'))
        if (request.digest != parent.request.digest or projection.plan.baseline_plan != plan
                or saved['baseline_plan'] != _json(plan) or saved['baseline'] != _json(projection.plan.baseline)
                or audit_projection_result(probe.request, request, plan, projection.plan.baseline, anchors, projection)
                or not audit_source_transport(preparation, replace(probe, request=request)).passed):
            raise ValueError('incoming-readback-parent-projection')
        if projection.status == 'projected-geometry-passed':
            result = dto('H_result', ProjectedHydrogenResult)
            if (audit_projected_hydrogen_result(parent, H_parameters, None, result)
                    or saved['H_parameters'] != _json(H_parameters) or not item['H_called']
                    or item['status'] != result.status or item['reason'] != result.reason):
                raise ValueError('incoming-readback-H')
        elif item['H_called'] or 'H_result' in saved or item['status'] != projection.status or item['reason'] != projection.reason:
            raise ValueError('incoming-readback-projection-negative')
        if (not item['baseline_called'] or not item['projection_called'] or item['target'] != identity.sequence_number
                or item['incoming_label'] != source.residue_name or item['candidate_states'] != list(proof.candidate_states)
                or item['baseline_status'] != projection.plan.baseline.status or item['projection_status'] != projection.status):
            raise ValueError('incoming-readback-stage-accounting')
        reference, scores = scores_after_proposals(directory, parent)
        if saved['postproposal_reference'] != _json(reference) or saved['postproposal_native_scores'] != _json(scores):
            raise ValueError('incoming-readback-native-scores')
    summary = parse_research_json((output/'summary.json').read_bytes())
    if summary != {'schema': 'dvbfixer.incoming-type-controls.v1', 'protocol_sha256': PROTOCOL_SHA256,
                   'declared': 3, 'controls': controls, 'charge_or_state_assigned': False, 'public_or_production_qualified': False}:
        raise ValueError('incoming-readback-summary')
    return {'schema': 'dvbfixer.incoming-type-readback.v1', 'controls': controls}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    parser.add_argument('--readback', action='store_true')
    args = parser.parse_args()
    print(json.dumps(readback_controls(args.output) if args.readback else run_controls(args.output), indent=2))
