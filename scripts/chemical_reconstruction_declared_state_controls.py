"""Frozen explicit PRO/HIP declarations; no neutral histidine graph transformation."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from dataclasses import replace
from pathlib import Path

from scripts.chemical_reconstruction_ccd_controls import ROOT, _write
from scripts.chemical_reconstruction_extended_h_controls import FIXTURES
from scripts.chemical_reconstruction_extended_h_controls import check_pins as check_extended_pins
from scripts.chemical_reconstruction_hydrogen_controls import (
    FFBONDED,
    FIXED_POLICY,
    GENERATED_POLICY,
    _execute,
    backbone_types,
)
from scripts.chemical_reconstruction_naming_controls import raw_inventory
from scripts.chemical_reconstruction_peptide_controls import load_ccd
from scripts.chemical_reconstruction_source_controls import _bytes

from dvbfixer.chemical_reconstruction.amide_authority import FFBONDED_SHA256, load_amide_authority
from dvbfixer.chemical_reconstruction.ccd_authority import load_geometry_authority
from dvbfixer.chemical_reconstruction.contracts import (
    ChemistryKind,
    Component,
    Connection,
    LinkKind,
    Observation,
    Operation,
    Origin,
    Ownership,
    ReconstructionRequest,
    Unit,
    digest_record,
)
from dvbfixer.chemical_reconstruction.declared_peptide import (
    PeptideStateDeclaration,
    compile_declared_hydrogens,
    derive_declared_peptide,
)
from dvbfixer.chemical_reconstruction.inventory_naming import transport_inventory
from dvbfixer.chemical_reconstruction.measurement import _dot, _sub, _unit
from dvbfixer.chemical_reconstruction.observed_peptide import select_interior_target
from dvbfixer.chemical_reconstruction.peptide_state import PeptideSubstitution, derive_local_peptide
from dvbfixer.chemical_reconstruction.source_preparation import (
    audit_source_transport,
    prepare_source,
    prepare_source_probe,
)
from dvbfixer.domain.force_field_naming import ForceFieldTarget, NamingProfile
from dvbfixer.domain.nonprotein_reconstruction import ExternalLink
from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef

PROTOCOL = ROOT / "docs/plans/chemical-reconstruction-declared-state-v1.md"
PROTOCOL_SHA256 = "eca02e9ee7d6e5efcf366345e95a97b4e149bba7ef9596ca889ed8dcad4d8e24"
ADDENDUM = ROOT / 'docs/plans/chemical-reconstruction-declared-state-metrics-v1.md'
ADDENDUM_SHA256 = '89d2bcac72afb3c5436789d390700aa8f2b04a80b20edc6135c9fc18d6cbc08b'
SOURCE_METRICS = ROOT / 'docs/plans/chemical-reconstruction-source-metrics-v1.md'
SOURCE_METRICS_SHA256 = '226c9a95716fa9b8b84a12c35175c82ac47469c95d6c720d8a72d14fe253efe7'
VARIANTS = ROOT / "src/dvbfixer/ffutils/variants.py"
VARIANTS_SHA256 = "5eccaf1aa6d8d056de85250f74ba8ddfa4128832d01c000ac65d964e41c53039"
SOURCE = "tests/fixtures/protein_context/1UBQ.pdb"
CASES = (("PRO", "PRO", "H"), ("HIS", "HIP", "H2"))


def check_pins():
    check_extended_pins()
    for path, pin in ((PROTOCOL, PROTOCOL_SHA256), (ADDENDUM, ADDENDUM_SHA256),
                      (SOURCE_METRICS, SOURCE_METRICS_SHA256), (VARIANTS, VARIANTS_SHA256)):
        if hashlib.sha256(path.read_bytes()).hexdigest() != pin:
            raise ValueError(f"frozen declared-state pin mismatch: {path}")


def declaration_inputs(source_name, state, identity, leaving):
    pins = json.loads((FIXTURES / "manifest.json").read_text())["components"]
    parent = load_geometry_authority(FIXTURES / f"{source_name}.cif", source_name, expected_sha256=pins[source_name]["sha256"])
    alanine = load_ccd("ALA")
    inventory = transport_inventory(raw_inventory(state), target=ForceFieldTarget.AMBER19SB, profile=NamingProfile.STANDARD)
    left = ComponentInstanceRef(identity.model, "__left_chemical_witness__", "1")
    right = ComponentInstanceRef(identity.model, "__right_chemical_witness__", "1")
    substitutions = tuple(PeptideSubstitution(ExactAtomRef(a, "C"), ExactAtomRef(b, "N"),
                          frozenset((ExactAtomRef(a, "OXT"), ExactAtomRef(a, "HXT"), ExactAtomRef(b, h))))
                          for a, b, h in ((left, identity, leaving), (identity, right, "H2")))
    declaration = PeptideStateDeclaration(identity, state, parent.digest, alanine.digest, digest_record(inventory.canonical), substitutions)
    return parent, alanine, inventory, declaration


def authority_for(source_name, state, identity, leaving, evidence):
    parent, alanine, inventory, declaration = declaration_inputs(source_name, state, identity, leaving)
    evidence.update(parent_authority=parent, alanine_authority=alanine, inventory_transport=inventory, declaration=declaration)
    try:
        derive_local_peptide(parent, alanine, identity, inventory.canonical)
    except ValueError as exc:
        evidence['old_derivation'] = {'status': 'refused', 'reason': str(exc)}
    else:
        raise AssertionError('old PRO/HIP refusal unexpectedly changed')
    authority = derive_declared_peptide(parent, alanine, identity, inventory.canonical, declaration)
    evidence['authority'] = authority
    return authority


def synthetic_base(source_name, state, leaving, evidence):
    identity = ComponentInstanceRef(1, "a", "82", "B", occurrence=2)
    authority = authority_for(source_name, state, identity, leaving, evidence)
    graph = authority.state.graph
    points = {a.name: a.ideal_position for a in graph.atoms if a.element not in {"H", "D"}}
    n = points['N']
    if state == 'PRO':
        neighbours = sorted(b.atom2 if b.atom1 == 'N' else b.atom1 for b in graph.bonds if 'N' in {b.atom1, b.atom2})
        a, b = (_unit(_sub(points[name], n)) for name in neighbours)
        direction = _unit(tuple(-a[i]-b[i] for i in range(3)))
    else:
        u, c = _unit(_sub(points['CA'], n)), _unit(_sub(points['C'], n))
        projection = _dot(c, u)
        v = _unit(tuple(c[i]-projection*u[i] for i in range(3)))
        angle = math.radians(121.96)
        direction = tuple(math.cos(angle)*u[i]+math.sin(angle)*v[i] for i in range(3))
    incoming = ExactAtomRef(ComponentInstanceRef(1, 'A', '81'), 'C')
    position = tuple(n[i]+1.335*direction[i] for i in range(3))
    observations = tuple(Observation(ExactAtomRef(identity, a.name), a.element, points[a.name], Unit.ANGSTROM)
                         for a in graph.atoms if a.element not in {'H', 'D'}) + (Observation(incoming, 'C', position, Unit.ANGSTROM),)
    pin = digest_record(('synthetic-declared-peptide-heavy-only-v1', state, observations, authority.declaration.digest))
    base = ReconstructionRequest(pin, Operation.HEAVY,
            (Component(identity, ChemistryKind.PROTEIN, (authority.state,)), Component(incoming.component, ChemistryKind.UNKNOWN)),
            observations, frozenset({identity}), Ownership(frozenset()),
            (Connection(ExternalLink(incoming, ExactAtomRef(identity, 'N'), 1.), LinkKind.COVALENT),),
            authority.constraints(minimum_volume=.05, planar_tolerance_angstrom=.05),
            origin=Origin('synthetic-declared-peptide-heavy-only-v1', pin))
    return base, authority, incoming


def observed_base(source_name, state, leaving, directory, evidence):
    original = (ROOT / SOURCE).read_bytes()
    _bytes(directory / 'original.pdb', original)
    preparation = prepare_source(original)
    _bytes(directory / 'prepared.pdb', preparation.prepared_text.encode('ascii'))
    evidence.update(source=SOURCE, original_sha256=preparation.original_sha256,
                    preparation={'digest': preparation.digest, 'prepared_sha256': preparation.prepared_sha256,
                                 'removed_metadata': preparation.removed_metadata})
    identity = select_interior_target(original, source_name, chain='A', number=None)
    authority = authority_for(source_name, state, identity, leaving, evidence)
    probe = prepare_source_probe(preparation, authority, withheld=frozenset())
    audit = audit_source_transport(preparation, probe)
    evidence['source_transport_audit'] = audit
    if not audit.passed:
        raise ValueError('declared-state-source-transport-failed')
    return probe.request, authority, probe.real_peptide_partners[0]


def compile_view(base, authority, incoming, evidence):
    parameters = load_amide_authority(FFBONDED, expected_sha256=FFBONDED_SHA256)
    types = None if authority.inventory.state_name == 'PRO' else backbone_types(authority.inventory.state_name)
    evidence.update(base_request=base, parameters=parameters, explicit_backbone_types=types)
    return compile_declared_hydrogens(base, authority, incoming, parameters, backbone_types=types,
                                     fixed_policy=FIXED_POLICY, generated_policy=GENERATED_POLICY)


def run_controls(output: Path):
    check_pins()
    output.mkdir(parents=True, exist_ok=False)
    records = []
    for arm in ('synthetic', 'observed'):
        for source_name, state, leaving in CASES:
            directory = output / f'{arm}-{state}'
            directory.mkdir()
            item = {'id': directory.name, 'arm': arm, 'source_component': source_name, 'requested_state': state, 'proposal_entry_called': False}
            evidence = {'control': item}
            try:
                base, authority, incoming = (synthetic_base(source_name, state, leaving, evidence) if arm == 'synthetic' else
                                             observed_base(source_name, state, leaving, directory, evidence))
                view = compile_view(base, authority, incoming, evidence)
                item.update(target=authority.identity.sequence_number, declaration_digest=authority.declaration.digest,
                            observations_unchanged=view.request.observations == base.observations,
                            expected_charge=authority.state.graph.net_charge)
                _execute(view, authority, item, evidence)
            except ValueError as exc:
                item.update(status='refused', reason=str(exc))
            _write(directory / 'evidence.json', evidence)
            records.append(item)
    for control, source_name, state, leaving in (('unchanged-HID', 'HIS', 'HID', 'H2'), ('unchanged-HIE', 'HIS', 'HIE', 'H2'),
                                                ('wrong-PRO-H2', 'PRO', 'PRO', 'H2'), ('wrong-source-pin', 'HIS', 'HIP', 'H2')):
        item = {'id': control, 'arm': 'chemical-proof-only', 'proposal_entry_called': False}
        evidence = {'control': item}
        parent, alanine, inventory, declaration = declaration_inputs(source_name, state, ComponentInstanceRef(1, 'A', '82', 'B'), leaving)
        if control == 'wrong-source-pin':
            declaration = replace(declaration, parent_digest='0'*64)
        evidence.update(parent_authority=parent, alanine_authority=alanine, inventory_transport=inventory, declaration=declaration)
        try:
            evidence['authority'] = derive_declared_peptide(parent, alanine, declaration.target, inventory.canonical, declaration)
        except ValueError as exc:
            item.update(status='refused', reason=str(exc))
        else:
            item.update(status='chemical-proof-passed')
        _write(output / f'{control}.json', evidence)
        records.append(item)
    check_pins()
    summary = {'schema': 'dvbfixer.declared-peptide-state-controls.v1', 'protocol_sha256': PROTOCOL_SHA256,
               'integration_addendum_sha256': ADDENDUM_SHA256,
               'source_metric_protocol_sha256': SOURCE_METRICS_SHA256,
               'declared': 8, 'executed_entries': len(records), 'controls': records, 'production_enabled': False,
               'neutral_histidine_transformations_or_full_polymer_qualified': False}
    _write(output / 'summary.json', summary)
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    print(json.dumps(run_controls(parser.parse_args().output), indent=2))
