#!/usr/bin/env python3
"""Compare preserved native cc94 arrays with raw PDBs, without resampling.

Checks whether gross generated-backbone bond defects precede coordinate
materialization and measures changes to junctions after fixed-atom restoration.
No coordinate files or acceptance decisions are rewritten.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import re
import tarfile
from pathlib import Path

import numpy as np


def key(atom):
    return atom["chain"], atom["residue_number"], atom["insertion_code"], atom["atom_name"]


def pdb_coordinates(data):
    result = {}
    for line in data.decode().splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            identity = line[21], line[22:26].strip(), line[26].strip(), line[12:16].strip()
            if identity in result:
                raise ValueError("duplicate raw PDB atom identity")
            result[identity] = np.array([float(line[i:i + 8]) for i in (30, 38, 46)])
    return result


def intra_residue_bonds(data, gap):
    from io import StringIO

    from openmm.app import PDBFile

    pdb = PDBFile(StringIO(data.decode()))
    generated = {(r["chain"], r["residue_number"], r["insertion_code"]) for r in gap["generated_residues"]}
    result = []
    for a, b in pdb.topology.bonds():
        residue = a.residue
        prefix = residue.chain.id, str(residue.id), residue.insertionCode.strip()
        if residue is b.residue and prefix in generated:
            result.append(((*prefix, a.name), (*prefix, b.name)))
    return result


def analyze_gap(gap, data, bond_pairs):
    axis = json.loads(data["native-axis-0.json"])
    with np.load(io.BytesIO(data["native-0.npz"]), allow_pickle=False) as payload:
        coordinates = payload["coordinates"]
    if coordinates.shape != (len(axis), 3) or not np.isfinite(coordinates).all():
        raise ValueError("invalid native coordinate shape")
    native = {key(atom): xyz for atom, xyz in zip(axis, coordinates, strict=True)}
    if len(native) != len(axis):
        raise ValueError("duplicate native atom identity")
    raw = pdb_coordinates(data["candidate.pdb"])
    differences, bonds = [], []
    for first, second in bond_pairs:
        if first not in native or second not in native:
            raise ValueError("native proposal omitted a bonded generated atom")
        n = float(np.linalg.norm(native[first] - native[second]))
        p = float(np.linalg.norm(raw[first] - raw[second]))
        differences.append(abs(n - p))
        bonds.append({"residue": {"chain": first[0], "residue_number": first[1], "insertion_code": first[2]},
                      "bond": [first[3], second[3]], "native_angstrom": n, "raw_pdb_angstrom": p})
    junctions = []
    for side, anchor, residue, first, second in (
            ("left", gap["left_anchor"], gap["generated_residues"][0], "C", "N"),
            ("right", gap["right_anchor"], gap["generated_residues"][-1], "N", "C")):
        if anchor is None:
            continue
        a = anchor["chain"], anchor["residue_number"], anchor["insertion_code"], first
        b = residue["chain"], residue["residue_number"], residue["insertion_code"], second
        n = float(np.linalg.norm(native[a] - native[b]))
        p = float(np.linalg.norm(raw[a] - raw[b]))
        junctions.append({"side": side, "native_angstrom": n, "raw_pdb_angstrom": p,
                          "native_in_1_20_1_45": 1.2 <= n <= 1.45, "raw_in_1_20_1_45": 1.2 <= p <= 1.45})
    summary = json.loads(data["summary.json"])
    return {"max_intra_residue_distance_change_angstrom": max(differences),
            "shortest_native_bond": min(bonds, key=lambda b: b["native_angstrom"]),
            "longest_native_bond": max(bonds, key=lambda b: b["native_angstrom"]),
            "junctions": junctions, "native_conditioning_fit": summary["native_conditioning_fit"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    cohort = json.loads((args.root / "cohort-final-outcomes.json").read_text())
    rows = []
    for case in cohort["rows"]:
        job = json.loads((args.root / case["job_record"]).read_text())
        execution = args.root / ("1N8Z-managed-followup" if case["pdb"] == "1N8Z" else "managed-execution")
        run = execution / "data/projects" / job["workspaceId"] / job["outputDir"]
        manifest = json.loads((run / "run-summary.json").read_text())
        expected = {r["file"]: r for r in manifest["diagnostics"]["files"]}
        request = json.loads((args.root / "cases" / case["pdb"] / "request-before-preparation.json").read_text())
        groups = {}
        with tarfile.open(run / "diagnostics.tar.gz", "r:gz") as archive:
            for member in archive:
                match = re.search(r"/seed-(\d+)/raw/gap-(\d+)/(native-0\.npz|native-axis-0\.json|candidate\.pdb|summary\.json)$", member.name)
                if not match:
                    continue
                if not member.isfile() or member.size > 16_000_000:
                    raise ValueError("native audit member is not bounded/regular")
                data = archive.extractfile(member).read()
                if hashlib.sha256(data).hexdigest() != expected[member.name]["sha256"]:
                    raise ValueError("native geometry artifact digest mismatch")
                group = groups.setdefault((int(match[1]), int(match[2])), {})
                if match[3] in group:
                    raise ValueError("duplicate native artifact")
                group[match[3]] = data
        if len(groups) != len(request["gaps"]) * 5:
            raise ValueError("native artifact inventory differs from declared five-seed gap plan")
        bond_cache = {}
        for (seed, gap_index), data in sorted(groups.items()):
            if set(data) != {"native-0.npz", "native-axis-0.json", "candidate.pdb", "summary.json"}:
                raise ValueError("missing native artifacts")
            gap = request["gaps"][gap_index]
            if gap_index not in bond_cache:
                bond_cache[gap_index] = intra_residue_bonds(data["candidate.pdb"], gap)
            rows.append({"pdb": case["pdb"], "seed": seed, "gap_index": gap_index,
                "gap_kind": gap["gap_kind"], "gap_length": len(gap["generated_residues"]),
                **analyze_gap(gap, data, bond_cache[gap_index])})
        print(case["pdb"], "native gap samples audited:", len(groups), flush=True)
    junctions = [j for row in rows for j in row["junctions"]]
    result = {"scope": "all generated intra-residue bonds from OpenMM topology, digest-bound native arrays versus raw PDBs; no sampling or minimization",
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "gap_samples": len(rows), "max_intra_residue_distance_change_angstrom": max(
            row["max_intra_residue_distance_change_angstrom"] for row in rows),
        "junctions": len(junctions), "native_junctions_outside_1_20_1_45": sum(
            not j["native_in_1_20_1_45"] for j in junctions),
        "raw_junctions_outside_1_20_1_45": sum(not j["raw_in_1_20_1_45"] for j in junctions),
        "junctions_good_native_bad_after_restoration": sum(j["native_in_1_20_1_45"] and not j["raw_in_1_20_1_45"] for j in junctions),
        "rows": rows}
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(json.dumps({k: v for k, v in result.items() if k != "rows"}, indent=2), flush=True)


if __name__ == "__main__":
    main()
