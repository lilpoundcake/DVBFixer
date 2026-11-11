#!/usr/bin/env python3
"""Admission of all three deposited 1BRS dimers, preserving every refusal.

This is not a search over masks to rescue a model. All three author-declared
assemblies and the same source span/mask policy are frozen before preparation.
They are copies of one complex, never independent family-level validation cases.
"""
from __future__ import annotations

import argparse
import json
import traceback
from pathlib import Path

import numpy as np
from scripts.audit_mac_reconstruction_smoke import inventory
from scripts.mac_reconstruction_comparison import digest, save
from scripts.prepare_mac_comparison_expansion import prepare
from scripts.prepare_mac_interface_control import EXPECTED, SOURCE

from dvbfixer.biological_assembly import parse_biological_assemblies


def run(output: Path):
    project = Path(__file__).resolve().parents[1]
    source = project / SOURCE
    if digest(source) != EXPECTED:
        raise ValueError("retained source differs from pinned experimental download")
    declared = (("1", "A", "D"), ("2", "B", "E"), ("3", "C", "F"))
    assemblies = parse_biological_assemblies(source.read_text().splitlines(keepends=True))
    cases = []
    for identifier, target, partner in declared:
        groups = [group for group in assemblies[identifier].groups if group.chains]
        if len(groups) != 1 or groups[0].chains != [target, partner]:
            raise ValueError("author-declared assembly chain identities changed")
        operators = groups[0].operators()
        if len(operators) != 1 or not np.array_equal(operators[0][1], np.eye(3)) or np.any(operators[0][2]):
            raise ValueError("this explicit projection requires the deposited identity operator")
        cases.append({"id": f"1BRS-assembly{identifier}-{target}3-110-{partner}3-89-gap83-90",
            "stratum": "deposited barnase/barstar biological dimer", "assembly": identifier,
            "fixture": SOURCE, "spans": [[target, 3, 110], [partner, 3, 89]], "hidden": [[target, 83, 90]]})
    output.mkdir(parents=True, exist_ok=False)
    save(output / "protocol.json", {"cases": cases, "source_sha256": EXPECTED,
        "all_three_author_deposited_assemblies_in_denominator": True, "mask_or_span_tuning": False,
        "existing_assembly1_refusal_unchanged": True, "family_independence_groups": 1, "sampling_calls": 0,
        "script_sha256": digest(Path(__file__)), "qualification": False})
    (output / "entry-script.py").write_bytes(Path(__file__).read_bytes())
    rows = []
    for spec in cases:
        try:
            manifest = prepare(spec, output / spec["id"])
            target, partner = spec["spans"][0][0], spec["spans"][1][0]
            source_coords = inventory(output / spec["id"] / "source.pdb")
            anchors = np.array([v[2] for a, v in source_coords.items() if a.chain == target and a.residue_number in {"82", "91"}])
            partners = np.array([v[2] for a, v in source_coords.items() if a.chain == partner])
            distances = np.linalg.norm(anchors[:, None]-partners[None, :], axis=2)
            rows.append({"case": spec["id"], "mask_admitted": True, "manifest": manifest,
                "observed_anchor_partner_min_angstrom": float(distances.min()), "source_only_contact_pairs_6A": int((distances <= 6).sum())})
        except Exception:
            error = traceback.format_exc()
            save(output / spec["id"] / "refusal.json", {"error": error, "sampling_executed": False})
            rows.append({"case": spec["id"], "mask_admitted": False, "error": error})
    save(output / "results.json", rows)
    print(json.dumps([{k: v for k, v in row.items() if k != "manifest"} for row in rows], indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    run(args.output.resolve())
