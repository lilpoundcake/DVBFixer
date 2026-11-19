"""Frozen joint sp3 comparison; old views are inputs, old generated H are not."""
from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import replace
from pathlib import Path

from scripts.chemical_reconstruction_ccd_controls import ROOT, _write
from scripts.chemical_reconstruction_hydrogen_controls import (
    FFBONDED,
    FIXED_POLICY,
    GENERATED_POLICY,
    OBSERVED,
    _execute,
    backbone_types,
    synthetic_view,
)
from scripts.chemical_reconstruction_hydrogen_controls import check_pins as check_H_pins
from scripts.chemical_reconstruction_naming_controls import raw_inventory
from scripts.chemical_reconstruction_peptide_controls import load_ccd
from scripts.chemical_reconstruction_source_controls import _bytes

from dvbfixer.chemical_reconstruction.amide_authority import FFBONDED_SHA256, load_amide_authority
from dvbfixer.chemical_reconstruction.hydrogen_completion import (
    compile_complete_hydrogens,
    compile_joint_hydrogens,
)
from dvbfixer.chemical_reconstruction.inventory_naming import transport_inventory
from dvbfixer.chemical_reconstruction.observed_peptide import select_interior_target
from dvbfixer.chemical_reconstruction.peptide_state import derive_local_peptide
from dvbfixer.chemical_reconstruction.source_preparation import (
    audit_source_transport,
    prepare_source,
    prepare_source_probe,
)
from dvbfixer.domain.force_field_naming import ForceFieldTarget, NamingProfile

PROTOCOL = ROOT / "docs/plans/chemical-reconstruction-joint-h-v1.md"
PROTOCOL_SHA256 = "b9450d98fe02afad1d87ae9cd3d7674eae2c112aa33b43cd94d6598b8adf42ef"


def check_pins():
    check_H_pins()
    if hashlib.sha256(PROTOCOL.read_bytes()).hexdigest() != PROTOCOL_SHA256:
        raise ValueError("frozen joint-H protocol changed")


def source_view(original: bytes, residue: str, number: str | None):
    prepared = prepare_source(original)
    target = select_interior_target(original, residue, chain="A", number=number)
    inventory = transport_inventory(raw_inventory(residue), target=ForceFieldTarget.AMBER19SB, profile=NamingProfile.STANDARD)
    authority = derive_local_peptide(load_ccd(residue), load_ccd("ALA"), target, inventory.canonical)
    probe = prepare_source_probe(prepared, authority, withheld=frozenset())
    parameters = load_amide_authority(FFBONDED, expected_sha256=FFBONDED_SHA256)
    old = compile_complete_hydrogens(probe.request, authority, probe.real_peptide_partners[0], parameters,
                                     backbone_types=backbone_types(residue), fixed_policy=FIXED_POLICY, generated_policy=GENERATED_POLICY)
    return old, authority, prepared, probe, inventory, parameters


def _joint(old, authority, item, evidence):
    view = compile_joint_hydrogens(old)
    evidence["old_view"] = old
    item.update(joint_groups=len(view.plan.joint_groups), request_and_constraints_unchanged=view.request == old.request,
                old_plan_digest=old.plan.digest)
    _execute(view, authority, item, evidence)


def run_controls(output: Path) -> dict:
    check_pins()
    output.mkdir(parents=True, exist_ok=False)
    records = []
    for name, angle in (*((n, 121.96) for n in ("ALA", "THR", "ILE", "LYS", "PHE")), ("ALA", 100.)):
        item = {"id": f"synthetic-{name}-{angle:g}", "arm": "synthetic", "proposal_entry_called": False}
        evidence = {"control": item}
        try:
            old, authority, source = synthetic_view(name, incoming_angle=angle)
            evidence.update(source)
            _joint(old, authority, item, evidence)
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
            old, authority, prepared, probe, inventory, parameters = source_view(original, residue, number)
            _bytes(directory / "prepared.pdb", prepared.prepared_text.encode("ascii"))
            evidence["preparation"] = {"digest": prepared.digest, "original_sha256": prepared.original_sha256,
                                       "prepared_sha256": prepared.prepared_sha256, "removed_metadata": prepared.removed_metadata}
            evidence.update(inventory_transport=inventory, authority=authority, parameters=parameters)
            source_audit = audit_source_transport(prepared, replace(probe, request=old.request))
            evidence["source_transport_audit"] = source_audit
            if not source_audit.passed:
                raise ValueError("joint-H-source-transport-audit-failed")
            item.update(target=authority.identity.sequence_number, source_observations=len(old.request.observations),
                        fractional_context=len(old.request.context_policy.fractional_atoms))
            _joint(old, authority, item, evidence)
        except ValueError as exc:
            item.update(status="refused", reason=str(exc))
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == source_sha
        _write(directory / "evidence.json", evidence)
        records.append(item)
    check_pins()
    summary = {"schema": "dvbfixer.joint-sp3-H-controls.v1", "protocol_sha256": PROTOCOL_SHA256,
               "declared": 12, "executed_entries": len(records), "controls": records,
               "retrospective_cases_with_prospective_algorithm": True,
               "full_polymer_or_general_H_qualification": False, "production_enabled": False}
    _write(output / "summary.json", summary)
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    print(json.dumps(run_controls(parser.parse_args().output), indent=2))
