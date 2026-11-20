"""Prospective three-source planar-port heavy projection; no H/native selection."""
from __future__ import annotations

import argparse
import hashlib
from dataclasses import replace
from pathlib import Path

from scripts.chemical_reconstruction_all20_controls import (
    FFBONDED,
    HEAVY_POLICY,
    authority_for,
    heldout_after_proposal,
)
from scripts.chemical_reconstruction_all20_controls import check_pins as check_all20_pins
from scripts.chemical_reconstruction_ccd_controls import ROOT, _write
from scripts.chemical_reconstruction_hydrogen_controls import RTP

from dvbfixer.chemical_reconstruction.construction import propose
from dvbfixer.chemical_reconstruction.contracts import digest_record
from dvbfixer.chemical_reconstruction.heldout_scoring import score_heavy_holdout
from dvbfixer.chemical_reconstruction.metric_projection import propose_projection
from dvbfixer.chemical_reconstruction.observed_peptide import (
    _parse_source,
    compile_peptide_request_plan,
    select_interior_target,
)
from dvbfixer.chemical_reconstruction.planar_port import (
    PlanarPortDeclaration,
    load_planar_port_authority,
    planar_port_constraints,
)
from dvbfixer.chemical_reconstruction.projection_audit import audit_projection_result
from dvbfixer.chemical_reconstruction.projection_contracts import append_projection_constraints
from dvbfixer.chemical_reconstruction.source_preparation import (
    audit_source_transport,
    prepare_source,
    prepare_source_probe,
)
from dvbfixer.domain.structure_identity import ExactAtomRef
from dvbfixer.pdb_component_reconstruction import _component_ref
from dvbfixer.rtp_parser import parse_rtp

PROTOCOL = ROOT/'docs/plans/chemical-reconstruction-port-projection-v1.md'
PROTOCOL_SHA256 = '7ec9dc0154a53f54bd7788f94ec0e94e431711025b52954cab35b40a5bbe1863'
SOURCES = (
    ('1UBQ', 'tests/fixtures/protein_context/1UBQ.pdb', 'd4a6812d8951cf6594e6a0763f089e35f5a80b62acb3c117b2c5565228a7b161'),
    ('1FR2', 'tests/fixtures/regressions/1FR2.pdb', 'ba93115f20b5e16eddbef8d70f1de7c25971b3c1c505acbaec5fb3670acf889a'),
    ('1A3N', 'tests/fixtures/protein_context/1A3N.pdb', 'a4db006269e9b6f9db0eb5447a5d40e59d67857484a9294b9ca71d642ffe9861'),
)
RTP_SHA256 = 'ece03322c6df812c1467bfdb67240e694230ea28ba675d86cccb8a4823e9a419'


def check_pins():
    check_all20_pins()
    for path, pin in ((PROTOCOL, PROTOCOL_SHA256), (RTP, RTP_SHA256), *((ROOT/p, pin) for _, p, pin in SOURCES)):
        if hashlib.sha256(path.read_bytes()).hexdigest() != pin:
            raise ValueError(f'frozen-projection-source-pin:{path}')


def run_controls(output: Path):
    import scipy

    check_pins()
    parameters = load_planar_port_authority(FFBONDED)
    _, templates = parse_rtp(RTP)
    output.mkdir(parents=True, exist_ok=False)
    records = []
    for name, relative, pin in SOURCES:
        directory = output/name
        directory.mkdir()
        item = {'id': name, 'source': relative, 'original_sha256': pin, 'projection_called': False}
        evidence = {'control': item, 'parameters': parameters}
        try:
            original_bytes = (ROOT/relative).read_bytes()
            (directory/'original.pdb').write_bytes(original_bytes)
            preparation = prepare_source(original_bytes)
            (directory/'prepared.pdb').write_bytes(preparation.prepared_text.encode('ascii'))
            identity = select_interior_target(original_bytes, 'PRO', chain='A')
            authority = authority_for('PRO', 'PRO', identity, evidence)
            probe = prepare_source_probe(preparation, authority, withheld=frozenset({'CB', 'CG', 'CD'}))
            source_audit = audit_source_transport(preparation, probe)
            evidence.update(probe=probe, source_transport_audit=source_audit)
            if not source_audit.passed:
                raise ValueError('projection-source-transport-failed')
            original = probe.request
            baseline_plan = compile_peptide_request_plan(original, HEAVY_POLICY)
            baseline = propose(original, baseline_plan)
            evidence.update(original_request=original, baseline_plan=baseline_plan, baseline=baseline)
            item.update(target=identity.sequence_number, baseline_status=baseline.status)
            if baseline.status != 'declared-geometry-passed':
                raise ValueError('unchanged-baseline-not-raw-rounded-qualified')
            incoming = probe.real_peptide_partners[0]
            _, atoms = _parse_source(original_bytes)
            incoming_resname = next(a.residue_name for a in atoms if _component_ref(a) == incoming.component and a.name == incoming.name)
            if incoming_resname not in templates:
                raise ValueError('actual-incoming-template-unresolved:'+incoming_resname)
            target_types = {row[0]: row[1] for row in templates['PRO'].atoms}
            incoming_types = {row[0]: row[1] for row in templates[incoming_resname].atoms}
            types = (incoming_types['C'], target_types['N'], target_types['CA'], target_types['CD'])
            n, ca, cd = (ExactAtomRef(identity, atom) for atom in ('N', 'CA', 'CD'))
            declaration = PlanarPortDeclaration(incoming, n, ca, cd, types, authority.digest, parameters.digest)
            extra = planar_port_constraints(original, authority, declaration, parameters)
            request = append_projection_constraints(original, extra)
            if not audit_source_transport(preparation, replace(probe, request=request)).passed:
                raise ValueError('new-port-request-changed-source-transport')
            anchors = ca, n, ExactAtomRef(identity, 'C')
            result = propose_projection(original, request, baseline_plan, baseline, anchors)
            item['projection_called'] = True
            evidence.update(port_declaration=declaration, request=request, anchors=anchors, result=result,
                            incoming_residue_name=incoming_resname, original_observations_digest=digest_record(original.observations))
            findings = audit_projection_result(original, request, baseline_plan, baseline, anchors, result)
            if findings:
                raise ValueError('saved-projection-audit:'+','.join(findings))
            item.update(status=result.status, reason=result.reason)
            # Native heavy reference extraction begins ONLY after proposal completion.
            reference = heldout_after_proposal('observed', authority, directory)
            evidence['postproposal_reference'] = reference
            scores = {}
            for arm, req, plan, proposal in (('baseline', original, baseline_plan, baseline), ('projection', request, result.plan, result)):
                for phase in ('raw', 'rounded'):
                    candidate = getattr(proposal, phase)
                    if candidate is not None and plan is not None:
                        scores[f'{arm}-{phase}'] = score_heavy_holdout(req, plan, candidate, reference)
            evidence['postproposal_native_scores'] = scores
        except ValueError as exc:
            item.update(status='refused', reason=str(exc))
        _write(directory/'evidence.json', evidence)
        records.append(item)
    check_pins()
    summary = {'schema': 'dvbfixer.generated-heavy-port-projection-controls.v1', 'protocol_sha256': PROTOCOL_SHA256,
               'scipy_version': scipy.__version__, 'declared': 3, 'controls': records, 'H_proposals': 0,
               'historical_WP2L_parent_replaced': False, 'native_H_or_production_qualified': False}
    _write(output/'summary.json', summary)
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    print(run_controls(parser.parse_args().output))
