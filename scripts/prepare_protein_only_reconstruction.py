#!/usr/bin/env python3
"""Explicit private protein-only preparation BEFORE freezing a new gap request.

The source and historical full-context requests are never modified. No generated
proposal is accepted here, and removed context must not be spliced back after
sampling/refinement. Production CLI/defaults are not changed.
"""
from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

from scripts.mac_reconstruction_comparison import digest, save

from dvbfixer.ffutils import PROTEIN_RESIDUES

POLICY = "private-protein-only-before-gap-reconstruction-v1"


def residue_key(line):
    return line[21], line[22:26].strip(), line[26].strip()


def protein_only(text):
    """Filter records without reserializing any retained observation.

    Canonical/recognized protein HETATMs are kept, not converted. Unknown ATOM
    residues are refused: a modified protein cannot be silently classified as a
    removable ligand. Explicit protein/nonprotein covalent links are refused.
    """
    lines = text.splitlines(keepends=True)
    if any(line.startswith(("MODEL ", "ENDMDL")) for line in lines):
        raise ValueError("select one model explicitly before protein-only preparation")
    atoms, kept, removed = {}, set(), []
    residue_classes = {}
    for line in lines:
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        serial = int(line[6:11])
        if serial in atoms:
            raise ValueError("ambiguous duplicate atom serial; do not rewrite observations implicitly")
        name = line[17:20].strip()
        is_protein = name in PROTEIN_RESIDUES
        if line.startswith("ATOM  ") and not is_protein:
            raise ValueError(f"unclassified ATOM residue {residue_key(line)} {name}; review modified protein/chemistry")
        key = residue_key(line)
        if key in residue_classes and residue_classes[key] != is_protein:
            raise ValueError("protein/nonprotein residue identity collision")
        residue_classes[key] = is_protein
        atoms[serial] = line
        if is_protein:
            kept.add(serial)
        else:
            removed.append({"serial": serial, "chain": key[0], "resid": key[1], "icode": key[2],
                            "resname": name, "atom": line[12:16].strip(), "record": line})
    if not kept:
        raise ValueError("no recognized protein atoms")
    # LINK is explicit source evidence: do not silently detach a PTM, glycan or
    # metal connection just to produce a canonical-only neural-network input.
    for line in lines:
        if line.startswith("LINK  "):
            first = (line[21], line[22:26].strip(), line[26].strip())
            second = (line[51], line[52:56].strip(), line[56].strip())
            classes = [residue_classes.get(key) for key in (first, second)]
            if None in classes:
                raise ValueError("unresolved LINK endpoints; review before cleaning")
            if classes[0] != classes[1]:
                raise ValueError("protein/nonprotein LINK requires an explicit chemistry decision")
        if line.startswith("CONECT"):
            endpoints = [int(line[i:i+5]) for i in range(6, len(line.rstrip()), 5) if line[i:i+5].strip()]
            # Follow the existing dangling-CONECT rule: never partially rewrite
            # a record with a serial absent from the atom inventory.
            if endpoints and all(i in atoms for i in endpoints):
                if any(i in kept for i in endpoints) and not all(i in kept for i in endpoints):
                    raise ValueError("protein/nonprotein CONECT requires an explicit chemistry decision")
    retained_names = {atoms[i][17:20].strip() for i in kept}
    output, discarded, last_serial = [], Counter(), None
    for line in lines:
        record = line[:6].strip()
        keep = True
        if record in {"ATOM", "HETATM"}:
            last_serial = int(line[6:11])
            keep = last_serial in kept
        elif record in {"ANISOU", "SIGATM", "SIGUIJ"}:
            serial = int(line[6:11])
            if serial not in atoms or line[6:27] != atoms[serial][6:27]:
                raise ValueError("unresolved atom-associated metadata")
            keep = serial in kept
        elif record == "TER":
            keep = last_serial in kept
        elif record == "CONECT":
            endpoints = [int(line[i:i+5]) for i in range(6, len(line.rstrip()), 5) if line[i:i+5].strip()]
            keep = bool(endpoints) and all(i in kept for i in endpoints)
        elif record == "LINK":
            keep = residue_classes[(line[21], line[22:26].strip(), line[26].strip())]
        elif record == "SSBOND":
            ends = [(line[15], line[17:21].strip(), line[21].strip()),
                    (line[29], line[31:35].strip(), line[35].strip())]
            if not all(residue_classes.get(key) for key in ends):
                raise ValueError("SSBOND does not resolve to retained protein residues")
        elif record in {"HET", "HETNAM", "HETSYN", "FORMUL"}:
            start = {"HET": 7, "HETNAM": 11, "HETSYN": 11, "FORMUL": 12}[record]
            keep = line[start:start+3].strip() in retained_names
        elif record == "MASTER":
            keep = False  # Source record counts no longer describe the projection.
        if keep:
            output.append(line)
        else:
            discarded[record] += 1
    projected = "".join(output)
    actual = [line for line in projected.splitlines(keepends=True) if line.startswith(("ATOM  ", "HETATM"))]
    expected = [line for line in lines if line.startswith(("ATOM  ", "HETATM")) and int(line[6:11]) in kept]
    if actual != expected:
        raise ValueError("protein observation preservation invariant failed")
    return projected, {"policy": POLICY, "retained_protein_atoms": len(kept),
        "retained_protein_hetatm_atoms": sum(line.startswith("HETATM") for line in actual),
        "removed_atoms": removed, "discarded_record_counts": dict(discarded),
        "retained_chains": sorted({line[21] for line in actual}),
        "retained_atom_records_exact": True, "reattach_removed_context": False,
        "limitation": "protein-only reconstruction, not full-complex or ligand-compatible qualification"}


def prepare(source, output):
    source = source.resolve()
    if output.exists():
        raise FileExistsError("protein-only preparation requires a fresh output directory")
    original = source.read_bytes()
    projected, report = protein_only(original.decode("ascii"))
    output.mkdir(parents=True, exist_ok=False)
    (output / "original.pdb").write_bytes(original)
    (output / "protein-only.pdb").write_bytes(projected.encode("ascii"))
    save(output / "preparation.json", {**report, "original_source": str(source),
        "original_sha256": digest(source), "retained_original_sha256": digest(output / "original.pdb"),
        "protein_only_sha256": digest(output / "protein-only.pdb"), "code_sha256": digest(Path(__file__)),
        "fresh_request_required": True, "historical_requests_modified": False,
        "coordinate_generation_calls": 0, "minimization_calls": 0})
    if source.read_bytes() != original:
        raise ValueError("original source changed during preparation")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    prepare(args.source, args.output)
