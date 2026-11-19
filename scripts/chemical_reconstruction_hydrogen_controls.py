"""Frozen complete-H controls with fixed heavy observations; no public route."""
from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from dataclasses import replace
from pathlib import Path

from scripts.chemical_reconstruction_amide_controls import FFBONDED, synthetic_case
from scripts.chemical_reconstruction_ccd_controls import ROOT, _write
from scripts.chemical_reconstruction_naming_controls import raw_inventory
from scripts.chemical_reconstruction_peptide_controls import RTP, load_ccd
from scripts.chemical_reconstruction_source_controls import _bytes
from scripts.chemical_reconstruction_source_controls import check_pins as check_source_pins

from dvbfixer.chemical_reconstruction.amide_authority import FFBONDED_SHA256, load_amide_authority
from dvbfixer.chemical_reconstruction.construction import propose
from dvbfixer.chemical_reconstruction.contracts import Operation, Origin, Ownership, digest_record
from dvbfixer.chemical_reconstruction.hydrogen_completion import (
    audit_hydrogen_cip,
    compile_complete_hydrogens,
)
from dvbfixer.chemical_reconstruction.inventory_naming import transport_inventory
from dvbfixer.chemical_reconstruction.observed_peptide import select_interior_target
from dvbfixer.chemical_reconstruction.peptide_state import derive_local_peptide
from dvbfixer.chemical_reconstruction.proposals import AuditPolicy
from dvbfixer.chemical_reconstruction.source_preparation import (
    audit_source_transport,
    prepare_source,
    prepare_source_probe,
)
from dvbfixer.domain.force_field_naming import ForceFieldTarget, NamingProfile
from dvbfixer.rtp_parser import parse_rtp

PROTOCOL = ROOT / "docs/plans/chemical-reconstruction-hydrogens-v1.md"
PROTOCOL_SHA256 = "758504f53d172a581d25b64ac3670095a6e1b8fd393f63786a0eadba3562b35c"
FIXED_POLICY = AuditPolicy(.1, 12., .7, .02)
GENERATED_POLICY = AuditPolicy(.02, 5., .7, .02)
OBSERVED = (
    ("1UBQ-ALA", "tests/fixtures/protein_context/1UBQ.pdb", "ALA", None),
    ("1UBQ-ILE", "tests/fixtures/protein_context/1UBQ.pdb", "ILE", None),
    ("1UBQ-LYS", "tests/fixtures/protein_context/1UBQ.pdb", "LYS", None),
    ("1UBQ-PHE45", "tests/fixtures/protein_context/1UBQ.pdb", "PHE", "45"),
    ("1FR2-THR", "tests/fixtures/regressions/1FR2.pdb", "THR", None),
    ("1A3N-ALA", "tests/fixtures/protein_context/1A3N.pdb", "ALA", None),
)


def backbone_types(name):
    _, templates = parse_rtp(RTP)
    types = {row[0]: row[1] for row in templates[name].atoms}
    return types["N"], types["H"], types["CA"]


def synthetic_view(name: str, *, incoming_angle: float = 121.96):
    old, _, evidence = synthetic_case(name, incoming_angle=incoming_angle)
    observations = tuple(o for o in old.observations if o.identity.component not in old.targets or o.element not in {"H", "D"})
    source_digest = digest_record(("synthetic-complete-H-heavy-only-v1", name, incoming_angle, observations))
    base = replace(old, source_sha256=source_digest, operation=Operation.HEAVY, ownership=Ownership(frozenset()),
                   observations=observations, origin=Origin("synthetic-complete-H-heavy-only-v1", source_digest))
    authority = evidence["authority"]
    incoming = old.connections[0].link.component_atom
    view = compile_complete_hydrogens(base, authority, incoming, evidence["parameters"], backbone_types=backbone_types(name),
                                       fixed_policy=FIXED_POLICY, generated_policy=GENERATED_POLICY)
    return view, authority, evidence


def _execute(view, authority, item, evidence):
    evidence["view"] = view
    item.update(expected_hydrogens=len(view.request.ownership.generated), groups=dict(Counter(g.kind for g in view.groups)),
                proposal_entry_called=True)
    result = propose(view.request, view.plan)
    evidence["result"] = result
    raw_cip = audit_hydrogen_cip(view, authority, result.raw) if result.raw is not None else None
    rounded_cip = audit_hydrogen_cip(view, authority, result.rounded) if result.rounded is not None else None
    evidence.update(raw_cip=raw_cip, rounded_cip=rounded_cip)
    status, reason = result.status, result.reason
    if status == "declared-geometry-passed" and not (raw_cip and raw_cip.passed and rounded_cip and rounded_cip.passed):
        status, reason = "conditional-CIP-rejected", "raw-or-rounded-conditional-CIP-failed"
    item.update(status=status, reason=reason, geometry_status=result.status, coordinates_produced=result.coordinates_produced,
                coordinate_constructor_executed=result.raw is not None,
                raw_cip_passed=bool(raw_cip and raw_cip.passed), rounded_cip_passed=bool(rounded_cip and rounded_cip.passed),
                proposed_hydrogens=len(result.raw.additions) if result.raw is not None else 0)


def check_pins():
    check_source_pins()
    if hashlib.sha256(PROTOCOL.read_bytes()).hexdigest() != PROTOCOL_SHA256:
        raise ValueError("frozen complete-H protocol changed")
    load_amide_authority(FFBONDED, expected_sha256=FFBONDED_SHA256)


def run_controls(output: Path) -> dict:
    check_pins()
    output.mkdir(parents=True, exist_ok=False)
    records = []
    for name, angle in (*((n, 121.96) for n in ("ALA", "THR", "ILE", "LYS", "PHE")), ("ALA", 100.)):
        item = {"id": f"synthetic-{name}-{angle:g}", "arm": "synthetic", "proposal_entry_called": False}
        evidence = {"control": item}
        try:
            view, authority, source = synthetic_view(name, incoming_angle=angle)
            evidence.update(source)
            _execute(view, authority, item, evidence)
        except ValueError as exc:
            item.update(status="refused", reason=str(exc))
        _write(output / f"{item['id']}.json", evidence)
        records.append(item)
    for name, path, residue, number in OBSERVED:
        directory = output / name
        directory.mkdir()
        original = (ROOT / path).read_bytes()
        source_sha = hashlib.sha256(original).hexdigest()
        _bytes(directory / "original.pdb", original)
        item = {"id": name, "arm": "observed", "source": path, "original_sha256": source_sha, "proposal_entry_called": False}
        evidence = {"control": item}
        try:
            prepared = prepare_source(original)
            _bytes(directory / "prepared.pdb", prepared.prepared_text.encode("ascii"))
            evidence["preparation"] = {"digest": prepared.digest, "original_sha256": prepared.original_sha256,
                                       "prepared_sha256": prepared.prepared_sha256, "removed_metadata": prepared.removed_metadata}
            target = select_interior_target(original, residue, chain="A", number=number)
            inventory = transport_inventory(raw_inventory(residue), target=ForceFieldTarget.AMBER19SB, profile=NamingProfile.STANDARD)
            authority = derive_local_peptide(load_ccd(residue), load_ccd("ALA"), target, inventory.canonical)
            evidence.update(inventory_transport=inventory, authority=authority)
            base_probe = prepare_source_probe(prepared, authority, withheld=frozenset())
            parameters = load_amide_authority(FFBONDED, expected_sha256=FFBONDED_SHA256)
            evidence["parameters"] = parameters
            view = compile_complete_hydrogens(base_probe.request, authority, base_probe.real_peptide_partners[0], parameters,
                                               backbone_types=backbone_types(residue), fixed_policy=FIXED_POLICY,
                                               generated_policy=GENERATED_POLICY)
            source_audit = audit_source_transport(prepared, replace(base_probe, request=view.request))
            evidence["source_transport_audit"] = source_audit
            if not source_audit.passed:
                raise ValueError("complete-H-source-transport-audit-failed")
            item.update(target=target.sequence_number, source_observations=len(view.request.observations),
                        fractional_context=len(view.request.context_policy.fractional_atoms))
            _execute(view, authority, item, evidence)
        except ValueError as exc:
            item.update(status="refused", reason=str(exc))
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == source_sha
        _write(directory / "evidence.json", evidence)
        records.append(item)
    check_pins()
    summary = {"schema": "dvbfixer.complete-selected-H-controls.v1", "protocol_sha256": PROTOCOL_SHA256,
               "declared": 12, "executed_entries": len(records), "controls": records,
               "full_polymer_or_general_H_qualification": False, "production_enabled": False}
    _write(output / "summary.json", summary)
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    print(json.dumps(run_controls(parser.parse_args().output), indent=2))
