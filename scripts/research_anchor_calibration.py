#!/usr/bin/env python3
"""Execute the predeclared experimental-coordinate pilot without changing policy."""
from __future__ import annotations

import argparse
import hashlib
import json
import urllib.request
from collections import Counter
from dataclasses import asdict
from pathlib import Path


def parse_model(text, model, chain):
    from dvbfixer.model.diffusion.contract import AtomIdentity, ResidueIdentity
    from dvbfixer.model.diffusion.protein_atoms import canonical_heavy_atom_identities

    active = 1
    coordinates, residue_names, exclusions, order, rows = {}, {}, {}, [], []
    for line in text.splitlines(keepends=True):
        if line.startswith("MODEL "):
            active = int(line[10:14])
        elif line.startswith("ENDMDL"):
            active = None
        elif line.startswith("ATOM  ") and active == model and line[21] == chain:
            if (line[76:78].strip() or line[12:16].strip()[:1]) in {"H", "D"}:
                continue
            rows.append(line)
            residue = ResidueIdentity(chain, line[22:26].strip(), line[26].strip())
            if residue not in residue_names:
                order.append(residue)
                residue_names[residue] = line[17:20]
            if residue_names[residue] != line[17:20]:
                raise ValueError("ambiguous experimental residue identity")
            atom = AtomIdentity(chain, residue.residue_number, residue.insertion_code, line[12:16].strip())
            if atom in coordinates or line[16:17].strip() or float(line[54:60]) != 1.0:
                exclusions[residue] = "altloc, duplicate atom or occupancy not exactly one"
            coordinates[atom] = tuple(float(line[i:i + 8]) for i in (30, 38, 46))
    if not coordinates:
        raise ValueError("declared model/chain is absent")
    for residue, name in residue_names.items():
        try:
            expected = set(canonical_heavy_atom_identities(residue, name))
        except ValueError:
            exclusions[residue] = "noncanonical residue"
            continue
        if not expected <= coordinates.keys():
            exclusions[residue] = "incomplete canonical heavy atoms"
    return coordinates, residue_names, exclusions, order, "".join(rows) + "END\n"


def main():
    import numpy as np
    from openmm.app import PDBFile

    from dvbfixer.diagnose.chemistry import check_bond_lengths, check_ca_chirality
    from dvbfixer.model.diffusion.anchor_calibration import PROPOSED_LIMITS, compare_anchor
    from dvbfixer.model.diffusion.contract import ResidueIdentity

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--declaration", type=Path, default=Path("docs/research/anchor-calibration-pilot-v1.json"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    declaration_bytes = args.declaration.read_bytes()
    declaration = json.loads(declaration_bytes)
    grid = [{"name": name, "backbone_rmsd": limit.backbone_rmsd_angstrom_max,
        "ca_displacement": limit.ca_displacement_angstrom_max,
        "max_heavy": limit.heavy_displacement_angstrom_max} for name, limit in PROPOSED_LIMITS]
    if grid != declaration["candidate_limits_angstrom"]:
        raise ValueError("calibration grid differs from predeclaration")
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / "declaration.json").write_bytes(declaration_bytes)
    (args.output / "run-script.py").write_bytes(Path(__file__).read_bytes())
    import dvbfixer.model.diffusion.anchor_calibration as module

    (args.output / "anchor_calibration.py").write_bytes(Path(module.__file__).read_bytes())
    texts, metadata, custody = {}, {}, {}
    for pdb in sorted({c[role] for c in declaration["cases"] for role in ("reference", "mobile")}):
        for kind, url in (("pdb", f"https://files.rcsb.org/download/{pdb}.pdb"),
                          ("metadata.json", f"https://data.rcsb.org/rest/v1/core/entry/{pdb}")):
            with urllib.request.urlopen(url, timeout=60) as response:
                data = response.read(4_000_001)
            if len(data) > 4_000_000:
                raise ValueError("calibration download exceeds declared input bound")
            (args.output / f"{pdb}.{kind}").write_bytes(data)
            custody[f"{pdb}.{kind}"] = {"url": url, "sha256": hashlib.sha256(data).hexdigest(), "bytes": len(data)}
            if kind == "pdb":
                texts[pdb] = data.decode("ascii")
            else:
                metadata[pdb] = json.loads(data)
    (args.output / "input-custody.json").write_text(json.dumps(custody, indent=2) + "\n")
    cases = []
    for index, case in enumerate(declaration["cases"]):
        ref, refnames, refexclude, order, reftext = parse_model(texts[case["reference"]], case["reference_model"], declaration["chain"])
        mobile, mobnames, mobexclude, moborder, mobtext = parse_model(texts[case["mobile"]], case["mobile_model"], declaration["chain"])
        if order != moborder or refnames != mobnames:
            raise ValueError("pilot requires identical sequence/identity ordering; no guessed mapping")
        excluded = {**refexclude, **mobexclude}
        valid = set(order) - excluded.keys()

        def residue_of(atom):
            return ResidueIdentity(atom.chain, atom.residue_number, atom.insertion_code)

        shared = tuple(sorted(a for a in ref.keys() & mobile.keys() if residue_of(a) in valid))
        # A terminal OXT mismatch also excludes the anchor, rather than hiding
        # an unmatched experimental atom in the maximum-heavy statistic.
        for residue in tuple(valid):
            a = {atom for atom in ref if residue_of(atom) == residue}
            b = {atom for atom in mobile if residue_of(atom) == residue}
            if a != b:
                valid.remove(residue)
                excluded[residue] = "heavy-atom axes differ (including terminal atoms)"
        shared = tuple(a for a in shared if residue_of(a) in valid)
        measurements, failures = [], []
        for i, residue in enumerate(order):
            if residue not in valid:
                continue
            contacts = tuple(a for a in shared if residue_of(a) not in set(order[max(0, i - 1):i + 2]))
            for frame in ("global-leave-anchor-out", "local-window-4-leave-anchor-out"):
                fit_residues = valid - {residue}
                if frame.startswith("local"):
                    fit_residues &= set(order[max(0, i - 4):i + 5])
                fit = tuple(a for a in shared if residue_of(a) in fit_residues and a.atom_name in {"N", "CA", "C", "O"})
                try:
                    measurement = compare_anchor(ref, mobile, residue, fit, contact_atoms=contacts,
                        contact_cutoff_angstrom=declaration["contact_cutoff_angstrom"])
                    measurements.append({"frame": frame, **asdict(measurement)})
                except ValueError as exc:
                    failures.append({"frame": frame, "residue": asdict(residue), "error": str(exc)})
        from io import StringIO

        geometry = []
        for label, text in (("reference", reftext), ("mobile", mobtext)):
            structure = PDBFile(StringIO(text))
            findings = [*check_ca_chirality(structure.topology, structure.positions),
                        *check_bond_lengths(structure.topology, structure.positions)]
            geometry.append({"role": label, "findings": [asdict(f) for f in findings]})
        summaries = []
        for frame in ("global-leave-anchor-out", "local-window-4-leave-anchor-out"):
            rows = [m for m in measurements if m["frame"] == frame]
            if not rows:
                raise ValueError("calibration produced no measurements in a declared frame")
            summaries.append({"frame": frame, "residues": len(rows),
                "limits_counts": dict(Counter(name for m in rows for name in m["limits_passed"])),
                "quantiles": {metric: {str(q): float(np.quantile([m[metric] for m in rows], q))
                    for q in (0.5, 0.9, 0.95, 1.0)} for metric in
                    ("backbone_rmsd_angstrom", "ca_displacement_angstrom", "maximum_heavy_displacement_angstrom")}})
        report = {"case": case, "measurements": measurements, "fit_failures": failures,
            "excluded": [{"residue": asdict(r), "reason": reason} for r, reason in excluded.items()],
            "geometry": geometry, "summaries": summaries}
        (args.output / f"case-{index}.json").write_text(json.dumps(report, indent=2) + "\n")
        cases.append({"case": case, "summaries": summaries, "exclusions": len(excluded)})
        print(json.dumps(cases[-1]), flush=True)
    report = {"scope": "one-protein methodological pilot; no recommended or adopted movement ceiling",
        "declaration_sha256": hashlib.sha256(declaration_bytes).hexdigest(), "cases": cases,
        "sources": {pdb: {"method": value["exptl"],
            "resolution_angstrom": value["rcsb_entry_info"].get("resolution_combined"),
            "citation": value.get("rcsb_primary_citation")} for pdb, value in metadata.items()},
        "production_limits_unchanged": True}
    (args.output / "summary.json").write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
