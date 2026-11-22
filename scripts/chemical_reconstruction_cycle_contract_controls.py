"""Nine source/mask planning controls; no candidates, native scores or solvers."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from scripts.chemical_reconstruction_all20_controls import FFBONDED, authority_for, load_parent
from scripts.chemical_reconstruction_ccd_controls import ROOT, _json, _write
from scripts.chemical_reconstruction_comparison_controls import check_pins as check_comparison_pins

from dvbfixer.chemical_reconstruction.contracts import digest_record
from dvbfixer.chemical_reconstruction.cycle_contract import (
    CyclePlanningAdmission,
    plan_cycle_construction,
)
from dvbfixer.chemical_reconstruction.incoming_source import (
    audit_incoming_source,
    bind_incoming_source,
)
from dvbfixer.chemical_reconstruction.incoming_type import candidate_family, prove_incoming_type
from dvbfixer.chemical_reconstruction.observed_peptide import _parse_source, select_interior_target
from dvbfixer.chemical_reconstruction.planar_port import (
    load_planar_port_authority,
    planar_port_constraints,
)
from dvbfixer.chemical_reconstruction.research_decode import (
    decode_research_record,
    parse_research_json,
)
from dvbfixer.chemical_reconstruction.source_preparation import prepare_source, prepare_source_probe
from dvbfixer.chemical_reconstruction.type_proven_port import TypeProvenPlanarPortDeclaration
from dvbfixer.domain.structure_identity import ExactAtomRef

PROTOCOL = ROOT/'docs/plans/chemical-reconstruction-cycle-contract-v1.md'
PROTOCOL_SHA256 = '18f005bc4f4ec17ec52b8edf1650dda507bda907dd1b46afff96ef8aef5d59c4'


def _canonical_json(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), ensure_ascii=False, allow_nan=False)


def check_pins():
    if hashlib.sha256(PROTOCOL.read_bytes()).hexdigest() != PROTOCOL_SHA256:
        raise ValueError('cycle-contract-protocol-pin')
    return check_comparison_pins()


def prepare_cell(row, mask, catalog):
    """Only source/chemical adapters, never an old result/candidate/reference."""
    data = (ROOT/row['path']).read_bytes()
    preparation = prepare_source(data)
    _, atoms = _parse_source(data)
    chain = next(a.chain_id for a in atoms if a.record == 'ATOM  ')
    identity = select_interior_target(data, 'PRO', chain=chain)
    if (identity.chain_id, identity.sequence_number, identity.insertion_code) != (row['chain'], row['residue_number'], row['insertion_code']):
        raise ValueError('cycle-contract-enrolled-target-binding')
    authority = authority_for('PRO', 'PRO', identity, {})
    probe = prepare_source_probe(preparation, authority, withheld=frozenset(mask))
    if {a.name for a in probe.request.ownership.generated} != set(mask):
        raise ValueError('cycle-contract-natural-incompleteness-not-frozen-mask')
    source = bind_incoming_source(preparation, probe, authority)
    proof = prove_incoming_type(probe.request, source, catalog, load_parent(candidate_family(source.residue_name, catalog)[0]))
    if audit_incoming_source(preparation, probe, authority, source):
        raise ValueError('cycle-contract-independent-incoming-byte-admission')
    parameters = load_planar_port_authority(FFBONDED)
    types = {r[0]: r[1] for s in catalog.states if s.name == 'PRO' for r in s.atoms}
    n, ca, cd = (ExactAtomRef(identity, a) for a in ('N', 'CA', 'CD'))
    declaration = TypeProvenPlanarPortDeclaration(source.incoming, n, ca, cd,
        (proof.atom_type, types['N'], types['CA'], types['CD']), authority.digest, parameters.digest, incoming_proof=proof)
    extra = planar_port_constraints(probe.request, authority, declaration, parameters)
    anchors = ca, n, ExactAtomRef(identity, 'C')
    evidence = {'original_source_sha256': preparation.original_sha256, 'prepared_source_sha256': preparation.prepared_sha256,
                'preparation_digest': preparation.digest, 'probe_digest': digest_record(probe),
                'source_request_digest': probe.request.digest, 'port_declaration': declaration, 'port_parameters': parameters}
    return probe.request, anchors, extra, evidence


def run_controls(output: Path):
    cohort, catalog = check_pins()  # All tracked sources/authorities before reservation.
    output.mkdir(parents=True, exist_ok=False)
    records = []
    for row in cohort['sources']:
        for mask in cohort['masks']:
            name = row['id']+'-'+ '-'.join(mask)
            directory = output/name
            directory.mkdir()
            item = {'id': name, 'source': row['id'], 'mask': mask, 'problem_compilation_called': False}
            evidence = {'control': item}
            try:
                source, anchors, extra, source_evidence = prepare_cell(row, mask, catalog)
            except ValueError as exc:
                item.update(status='source-refused', reason=str(exc))
            else:
                result = plan_cycle_construction(source, anchors, extra_constraints=extra)
                evidence.update(source_evidence=source_evidence, result=result)
                item.update(problem_compilation_called=True, status=result.status, reason=result.reason)
                if result.problem is not None:
                    coverage = result.problem.coverage
                    item.update(problem_digest=result.problem.digest, generated_atoms=len(coverage.variables),
                        all_cycle_edges=len(coverage.cycle_edges), generated_cycle_edges=len(coverage.generated_cycle_edges),
                        generated_blocks=len(coverage.blocks), total_scalar_terms=len(coverage.metrics)+coverage.constraint_scalar_count)
            _write(directory/'evidence.json', evidence)
            records.append(item)
    check_pins()
    summary = {'schema': 'dvbfixer.cycle-contract-controls.v1', 'protocol_sha256': PROTOCOL_SHA256,
               'declared_cells': 9, 'source_systems': 3, 'controls': records, 'new_coordinate_proposals': 0,
               'initializer_or_solver_defined': False, 'geometry_or_source_coverage_qualified': False, 'production_enabled': False}
    _write(output/'summary.json', summary)
    return summary


def readback_controls(output: Path):
    cohort, catalog = check_pins()
    records = []
    for row in cohort['sources']:
        for mask in cohort['masks']:
            name = row['id']+'-'+ '-'.join(mask)
            with (output/name/'evidence.json').open('rb') as handle:
                saved = parse_research_json(handle.read(64*1024*1024+1))
            item = saved['control']
            if (item['id'], item['source'], item['mask']) != (name, row['id'], mask):
                raise ValueError('cycle-contract-saved-enrollment')
            records.append(item)
            try:
                source, anchors, extra, source_evidence = prepare_cell(row, mask, catalog)
            except ValueError as exc:
                expected_refusal = {'id': name, 'source': row['id'], 'mask': mask,
                    'problem_compilation_called': False, 'status': 'source-refused', 'reason': str(exc)}
                if set(saved) != {'control'} or _canonical_json(item) != _canonical_json(expected_refusal):
                    raise ValueError('cycle-contract-saved-source-refusal') from exc
                continue
            result = decode_research_record(json.dumps(saved['result'], allow_nan=False).encode(), CyclePlanningAdmission)
            expected = plan_cycle_construction(source, anchors, extra_constraints=extra)
            expected_item = {'id': name, 'source': row['id'], 'mask': mask, 'problem_compilation_called': True,
                             'status': expected.status, 'reason': expected.reason}
            if expected.problem is not None:
                coverage = expected.problem.coverage
                expected_item.update(problem_digest=expected.problem.digest, generated_atoms=len(coverage.variables),
                    all_cycle_edges=len(coverage.cycle_edges), generated_cycle_edges=len(coverage.generated_cycle_edges),
                    generated_blocks=len(coverage.blocks), total_scalar_terms=len(coverage.metrics)+coverage.constraint_scalar_count)
            if (digest_record(result) != digest_record(expected)
                    or _canonical_json(saved['source_evidence']) != _canonical_json(_json(source_evidence))
                    or _canonical_json(item) != _canonical_json(expected_item)
                    or set(saved) != {'control', 'source_evidence', 'result'}):
                raise ValueError('cycle-contract-saved-source-or-problem-binding')
    summary = {'schema': 'dvbfixer.cycle-contract-controls.v1', 'protocol_sha256': PROTOCOL_SHA256,
               'declared_cells': 9, 'source_systems': 3, 'controls': records, 'new_coordinate_proposals': 0,
               'initializer_or_solver_defined': False, 'geometry_or_source_coverage_qualified': False, 'production_enabled': False}
    with (output/'summary.json').open('rb') as handle:
        if _canonical_json(parse_research_json(handle.read(64*1024*1024+1))) != _canonical_json(summary):
            raise ValueError('cycle-contract-saved-summary')
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    parser.add_argument('--readback', action='store_true')
    args = parser.parse_args()
    print(json.dumps(readback_controls(args.output) if args.readback else run_controls(args.output), indent=2))
