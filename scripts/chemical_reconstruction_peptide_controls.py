"""Conditional local chemistry and whole-observed-source probes; no public workflow."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from scripts.chemical_reconstruction_ccd_controls import FIXTURES, ROOT, _write

from dvbfixer.chemical_reconstruction.ccd_authority import load_geometry_authority
from dvbfixer.chemical_reconstruction.construction import propose
from dvbfixer.chemical_reconstruction.observed_peptide import (
    compile_peptide_plan,
    prepare_observed_probe,
    select_interior_target,
)
from dvbfixer.chemical_reconstruction.peptide_state import TemplateInventory, derive_local_peptide
from dvbfixer.chemical_reconstruction.proposals import AuditPolicy
from dvbfixer.domain.force_field_naming import GROMACS_AMBER19SB_ATOM_RENAMES
from dvbfixer.domain.structure_identity import ComponentInstanceRef
from dvbfixer.rtp_parser import parse_rtp

RTP = ROOT / "FF/amber19sb.ff/aminoacids.rtp"
RTP_SHA256 = "ece03322c6df812c1467bfdb67240e694230ea28ba675d86cccb8a4823e9a419"
POLICY = AuditPolicy(.10, 12., .7, .02)
PROTOCOL_SHA256 = "65578a5517e53b5d7da71e9cc65ccae8b87e468cf23dd1fe502ffdd8a511e9ad"
OBSERVED = (
    ("1UBQ-PHE45", "tests/fixtures/protein_context/1UBQ.pdb", "PHE", "45", frozenset({"CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"})),
    ("1FR2-first-THR", "tests/fixtures/regressions/1FR2.pdb", "THR", None, frozenset({"CB", "OG1", "CG2"})),
    ("1A3N-first-ALA", "tests/fixtures/protein_context/1A3N.pdb", "ALA", None, frozenset({"CB"})),
)


def inventory(name: str) -> TemplateInventory:
    content = RTP.read_bytes()
    if hashlib.sha256(content).hexdigest() != RTP_SHA256:
        raise ValueError("pinned RTP source mismatch")
    _, templates = parse_rtp(RTP)
    template = templates[name]
    # Canonical owner maps PDB -> GROMACS; this inventory transport uses its
    # explicit inverse. Do not try alternate maps after a mismatch.
    aliases = {value: key for key, value in GROMACS_AMBER19SB_ATOM_RENAMES.get(name, {}).items()}
    raw_names = {a[0] for a in template.atoms}
    names = {aliases.get(n, n) for n in raw_names}
    if len(names) != len(raw_names):
        raise ValueError("non-bijective canonical inventory name map")
    bonds = set()
    external = []
    for row in (*template.bonds, *template.impropers, *template.cmap):
        if any(n.startswith(("-", "+")) for n in row):
            external.append(tuple(row))
    for first, second in template.bonds:
        if first in raw_names and second in raw_names:
            bonds.add(frozenset((aliases.get(first, first), aliases.get(second, second))))
    return TemplateInventory(str(RTP), RTP_SHA256, name, frozenset(names), frozenset(bonds),
                             tuple(sorted(aliases.items())), tuple(external))


def load_ccd(name: str):
    pins = json.loads((FIXTURES / "manifest.json").read_text())["components"]
    return load_geometry_authority(FIXTURES / f"{name}.cif", name, expected_sha256=pins[name]["sha256"])


def run_probes(output: Path) -> dict:
    protocol = ROOT / "docs/plans/chemical-reconstruction-peptide-v1.md"
    if hashlib.sha256(protocol.read_bytes()).hexdigest() != PROTOCOL_SHA256:
        raise ValueError("frozen peptide protocol changed")
    source_pins = {path: digest for digest, path in (line.split("  ", 1) for line in
                   (ROOT / "tests/fixtures/MANIFEST.sha256").read_text().splitlines())}
    for _, path, _, _, _ in OBSERVED:
        if hashlib.sha256((ROOT / path).read_bytes()).hexdigest() != source_pins[path]:
            raise ValueError("observed source pin mismatch")
    output.mkdir(parents=True, exist_ok=False)
    records = []
    witness = load_ccd("ALA")
    for source, state in (("ALA", "ALA"), ("THR", "THR"), ("ILE", "ILE"), ("LYS", "LYS"), ("PHE", "PHE"), ("LYS", "LYN")):
        item = {"id": f"chemistry-{source}-{state}", "kind": "conditional-chemistry", "source": source, "state": state}
        evidence = {"control": item}
        try:
            parent, selected_inventory = load_ccd(source), inventory(state)
            evidence.update(parent_authority=parent, inventory=selected_inventory)
            authority = derive_local_peptide(parent, witness, ComponentInstanceRef(1, "A", "82", "B", occurrence=2), selected_inventory)
        except ValueError as exc:
            item.update(status="refused", reason=str(exc))
        else:
            evidence["authority"] = authority
            item.update(status="conditional-chemistry-passed", authority_digest=authority.digest,
                        removed=sorted(authority.removed_atoms), witness_geometry_qualified=False,
                        actual_neighbour_chemistry_verified=False)
        _write(output / f"{item['id']}.json", evidence)
        records.append(item)
    for name, path, residue, number, mask in OBSERVED:
        original = (ROOT / path).read_bytes()
        source_sha = hashlib.sha256(original).hexdigest()
        for arm in ("heavy", "hydrogens"):
            item = {"id": f"{name}-{arm}", "kind": "observed-probe", "source": path, "source_sha256": source_sha}
            evidence = {"control": item}
            try:
                target = select_interior_target(original, residue, chain="A", number=number)
                authority = derive_local_peptide(load_ccd(residue), witness, target, inventory(residue))
                evidence["authority"] = authority
                if arm == "hydrogens":
                    from dvbfixer.chemical_reconstruction.contracts import Operation
                    prepare_observed_probe(original, authority, withheld=frozenset(), operation=Operation.HYDROGENS)
                    raise AssertionError("H operation unexpectedly admitted")
                probe = prepare_observed_probe(original, authority, withheld=mask)
                evidence["probe"] = probe
                from dvbfixer.chemical_reconstruction.admission import admit
                evidence["planning"] = admit(probe.request)
                item.update(source_atom_count=probe.source_atom_count, retained_atom_count=len(probe.request.observations),
                            withheld_count=len(probe.withheld), unassessed_internal_connections=len(probe.unassessed_internal_connections))
                plan = compile_peptide_plan(probe, POLICY)
                evidence["plan"] = plan
                result = propose(probe.request, plan)
                evidence["result"] = result
                item.update(status=result.status, reason=result.reason)
            except ValueError as exc:
                item.update(status="refused", reason=str(exc))
            assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == source_sha
            _write(output / f"{item['id']}.json", evidence)
            records.append(item)
    summary = {"schema": "dvbfixer.conditional-peptide-probes.v1", "protocol_sha256": hashlib.sha256(protocol.read_bytes()).hexdigest(),
               "declared": 12, "executed_entries": len(records), "controls": records,
               "full_polymer_or_h_or_complex_qualification": False, "production_enabled": False}
    _write(output / "summary.json", summary)
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    summary = run_probes(parser.parse_args().output)
    print(json.dumps({k: v for k, v in summary.items() if k != "controls"}, indent=2))
