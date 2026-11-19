"""Frozen private naming controls; no coordinate constructor or source preparation."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from scripts.chemical_reconstruction_ccd_controls import FIXTURES, ROOT, _write
from scripts.chemical_reconstruction_peptide_controls import RTP, RTP_SHA256, load_ccd

from dvbfixer.chemical_reconstruction.inventory_naming import transport_inventory
from dvbfixer.chemical_reconstruction.peptide_state import TemplateInventory, derive_local_peptide
from dvbfixer.domain.force_field_naming import ForceFieldTarget, NamingProfile
from dvbfixer.domain.structure_identity import ComponentInstanceRef
from dvbfixer.rtp_parser import parse_rtp

PROTOCOL = ROOT / "docs/plans/chemical-reconstruction-naming-v1.md"
PROTOCOL_SHA256 = "4944476be83f67cc440894326dd066c4b59d8d10f720a70a180d2cbca36f7230"
OWNER = ROOT / "src/dvbfixer/domain/force_field_naming.py"
OWNER_SHA256 = "bc35adb9c9513356fbeefdf7bce88dd78883144697df9bc6fbed585bda5be566"
CCD_MANIFEST_SHA256 = "29c4f83b28435a34b666b5c1489c4f3423c53beb276e75d16cd0d8c3c055dece"
PAIRS = (("ALA", "ALA"), ("THR", "THR"), ("ILE", "ILE"),
         ("LYS", "LYS"), ("PHE", "PHE"), ("LYS", "LYN"))


def raw_inventory(name: str) -> TemplateInventory:
    if hashlib.sha256(RTP.read_bytes()).hexdigest() != RTP_SHA256:
        raise ValueError("pinned RTP source mismatch")
    _, templates = parse_rtp(RTP)
    template = templates[name]
    atoms = frozenset(a[0] for a in template.atoms)
    if len(atoms) != len(template.atoms):
        raise ValueError("duplicate raw RTP atom")
    bonds, external = [], []
    for row in (*template.bonds, *template.impropers, *template.cmap):
        if any(n.startswith(("-", "+")) for n in row):
            external.append(tuple(row))
        if any(n not in atoms and not n.startswith(("-", "+")) for n in row):
            raise ValueError("unknown local RTP token")
    for a, b in template.bonds:
        if a in atoms and b in atoms:
            bonds.append(frozenset((a, b)))
    if len(set(bonds)) != len(bonds):
        raise ValueError("duplicate raw RTP bond")
    return TemplateInventory(str(RTP), RTP_SHA256, name, atoms, frozenset(bonds), external_tokens=tuple(external))


def check_pins() -> None:
    for path, digest in ((PROTOCOL, PROTOCOL_SHA256), (OWNER, OWNER_SHA256), (RTP, RTP_SHA256),
                         (FIXTURES / "manifest.json", CCD_MANIFEST_SHA256)):
        if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError(f"frozen naming source pin mismatch: {path.name}")
    pins = json.loads((FIXTURES / "manifest.json").read_bytes())["components"]
    for name in {source for source, _ in PAIRS}:
        if hashlib.sha256((FIXTURES / f"{name}.cif").read_bytes()).hexdigest() != pins[name]["sha256"]:
            raise ValueError("pinned CCD source mismatch")


def run_controls(output: Path) -> dict:
    check_pins()
    output.mkdir(parents=True, exist_ok=False)
    witness = load_ccd("ALA")
    records = []
    for source, state in PAIRS:
        item = {"id": f"{source}-{state}", "source": source, "state": state}
        evidence = {"control": item}
        try:
            raw = raw_inventory(state)
            transport = transport_inventory(raw, target=ForceFieldTarget.AMBER19SB, profile=NamingProfile.STANDARD)
            evidence["transport"] = transport
            authority = derive_local_peptide(load_ccd(source), witness,
                                              ComponentInstanceRef(1, "a", "82", "B", occurrence=2), transport.canonical)
        except ValueError as exc:
            item.update(status="refused", reason=str(exc))
        else:
            evidence["authority"] = authority
            item.update(status="conditional-inventory-proof-passed", transport_digest=transport.digest,
                        authority_digest=authority.digest, renamed_atoms=len(transport.canonical.aliases))
        _write(output / f"{item['id']}.json", evidence)
        records.append(item)
    check_pins()
    summary = {"schema": "dvbfixer.inventory-naming-controls.v1", "protocol_sha256": PROTOCOL_SHA256,
               "owner_sha256": OWNER_SHA256, "rtp_sha256": RTP_SHA256,
               "declared": 6, "executed_entries": len(records), "controls": records,
               "coordinate_constructor_executions": 0, "h_geometry_qualified": False, "production_enabled": False}
    _write(output / "summary.json", summary)
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    print(json.dumps(run_controls(parser.parse_args().output), indent=2))
