"""Prospective source-bound heavy probes; no atom editing or public output route."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from scripts.chemical_reconstruction_ccd_controls import ROOT, _write
from scripts.chemical_reconstruction_naming_controls import check_pins as check_naming_pins
from scripts.chemical_reconstruction_naming_controls import raw_inventory
from scripts.chemical_reconstruction_peptide_controls import OBSERVED, POLICY, load_ccd

from dvbfixer.chemical_reconstruction.admission import admit
from dvbfixer.chemical_reconstruction.construction import propose
from dvbfixer.chemical_reconstruction.inventory_naming import transport_inventory
from dvbfixer.chemical_reconstruction.observed_peptide import (
    compile_peptide_plan,
    select_interior_target,
)
from dvbfixer.chemical_reconstruction.peptide_state import derive_local_peptide
from dvbfixer.chemical_reconstruction.source_preparation import (
    audit_source_transport,
    prepare_source,
    prepare_source_probe,
)
from dvbfixer.domain.force_field_naming import ForceFieldTarget, NamingProfile

PROTOCOL = ROOT / "docs/plans/chemical-reconstruction-source-v1.md"
PROTOCOL_SHA256 = "5584293d5b08ff138155e6142670d7cb3242f9f6ef2f78de0ae421a4787325cb"


def check_pins() -> None:
    check_naming_pins()
    if hashlib.sha256(PROTOCOL.read_bytes()).hexdigest() != PROTOCOL_SHA256:
        raise ValueError("frozen source protocol changed")
    pins = {path: digest for digest, path in (line.split("  ", 1) for line in
            (ROOT / "tests/fixtures/MANIFEST.sha256").read_text().splitlines())}
    for _, path, _, _, _ in OBSERVED:
        if hashlib.sha256((ROOT / path).read_bytes()).hexdigest() != pins[path]:
            raise ValueError("observed source pin mismatch")


def _bytes(path: Path, value: bytes) -> None:
    with path.open("xb") as handle:
        handle.write(value)


def run_controls(output: Path) -> dict:
    check_pins()
    output.mkdir(parents=True, exist_ok=False)
    records = []
    for name, path, residue, number, mask in OBSERVED:
        directory = output / name
        directory.mkdir()
        original = (ROOT / path).read_bytes()
        source_sha = hashlib.sha256(original).hexdigest()
        _bytes(directory / "original.pdb", original)
        item = {"id": name, "source": path, "original_sha256": source_sha,
                "preparation_completed": False, "request_created": False, "constructor_executed": False,
                "coordinates_produced": False}
        evidence = {"control": item}
        try:
            preparation = prepare_source(original)
            _bytes(directory / "prepared.pdb", preparation.prepared_text.encode("ascii"))
            evidence["preparation"] = {
                "digest": preparation.digest, "original_sha256": preparation.original_sha256,
                "prepared_sha256": preparation.prepared_sha256, "removed_metadata": preparation.removed_metadata}
            item.update(preparation_completed=True, removed_metadata_records=len(preparation.removed_metadata),
                        prepared_sha256=preparation.prepared_sha256)
            target = select_interior_target(original, residue, chain="A", number=number)
            transport = transport_inventory(raw_inventory(residue), target=ForceFieldTarget.AMBER19SB, profile=NamingProfile.STANDARD)
            authority = derive_local_peptide(load_ccd(residue), load_ccd("ALA"), target, transport.canonical)
            evidence.update(inventory_transport=transport, authority=authority)
            probe = prepare_source_probe(preparation, authority, withheld=mask)
            evidence.update(probe=probe, planning=admit(probe.request))
            transport_audit = audit_source_transport(preparation, probe)
            evidence["source_transport_audit"] = transport_audit
            if not transport_audit.passed:
                raise ValueError("source-transport-audit-failed")
            item.update(request_created=True, retained_atom_count=len(probe.request.observations),
                        source_atom_count=probe.source_atom_count, withheld_count=len(probe.withheld),
                        fractional_context_count=len(probe.request.context_policy.fractional_atoms))
            plan = compile_peptide_plan(probe, POLICY)
            evidence["plan"] = plan
            item["constructor_executed"] = True
            result = propose(probe.request, plan)
            evidence["result"] = result
            item.update(status=result.status, reason=result.reason, coordinates_produced=result.coordinates_produced)
        except ValueError as exc:
            item.update(status="refused", reason=str(exc))
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == source_sha
        _write(directory / "evidence.json", evidence)
        records.append(item)
    check_pins()
    summary = {"schema": "dvbfixer.prepared-observed-source-controls.v1", "protocol_sha256": PROTOCOL_SHA256,
               "declared": 3, "executed_entries": len(records), "controls": records,
               "full_complex_or_h_qualification": False, "production_enabled": False}
    _write(output / "summary.json", summary)
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    print(json.dumps(run_controls(parser.parse_args().output), indent=2))
