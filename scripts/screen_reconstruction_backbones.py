#!/usr/bin/env python3
"""Retain experimental 3D segment diagnostics; never claim threading acceptance."""
from __future__ import annotations

import argparse
import hashlib
import json
from collections import defaultdict
from dataclasses import asdict
from pathlib import Path

import numpy as np

from dvbfixer.model.diffusion.backbone_intersections import BackboneSegment, screen_segments
from dvbfixer.model.diffusion.contract import AtomIdentity, DiffusionRequest
from dvbfixer.model.diffusion.loop_surface_probe import probe_loop_surface
from dvbfixer.model.diffusion.protein_atoms import CANONICAL_HEAVY_ATOMS
from dvbfixer.pdbutils.sequence_evidence import observed_peptide_continuity


def project(text):
    """Canonical heavy-backbone bonds; actual peptide evidence, never number gaps."""
    lines = text.splitlines()
    if any(line.startswith(("MODEL ", "ENDMDL")) for line in lines):
        raise ValueError("select one explicit model before screening")
    residues = defaultdict(dict)
    chains = defaultdict(list)
    serials = {}
    skipped = set()
    for line in lines:
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        if len(line) < 54:
            raise ValueError("short atom record")
        chain, number, code = line[21], line[22:26].strip(), line[26].strip()
        if line[17:20].strip() not in CANONICAL_HEAVY_ATOMS:
            skipped.add((chain, number, code, line[17:20].strip()))
            continue
        key = chain, number, code
        identity = AtomIdentity(*key, line[12:16].strip())
        serial = int(line[6:11])
        if serial in serials or identity.atom_name in residues[key]:
            raise ValueError("ambiguous atom serial/identity")
        serials[serial] = identity
        if not residues[key]:
            chains[chain].append((number, code))
        xyz = tuple(float(line[i:i+8]) for i in (30, 38, 46))
        if identity.atom_name in ("N", "CA", "C") and (line[16].strip() or float(line[54:60]) != 1):
            raise ValueError("select complete-occupancy, unambiguous backbone before screening")
        residues[key][identity.atom_name] = identity, xyz
    segments = []
    breaks = []
    bonds = set()
    for key, atoms in residues.items():
        if not {"N", "CA", "C"} <= atoms.keys():
            raise ValueError(f"incomplete backbone: {key}")
        for left, right in (("N", "CA"), ("CA", "C")):
            a, b = atoms[left], atoms[right]
            segments.append(BackboneSegment(a[0], b[0], a[1], b[1]))
    for chain, keys in chains.items():
        evidence = observed_peptide_continuity(lines, chain, keys)
        for left, right, connected in zip(keys, keys[1:], evidence):
            if connected is True:
                a, b = residues[(chain, *left)]["C"], residues[(chain, *right)]["N"]
                segments.append(BackboneSegment(a[0], b[0], a[1], b[1]))
            else:
                breaks.append({"chain": chain, "left": left, "right": right, "connectivity": connected})
    # Explicit source links are retained, including sulfur-sulfur CONECT edges.
    # No distance-only disulfide inference and no closure/topology claim.
    for line in lines:
        if line.startswith("CONECT"):
            values = [int(line[i:i+5]) for i in range(6, len(line), 5) if line[i:i+5].strip()]
            if values and all(n in serials for n in values):
                bonds.update(frozenset((serials[values[0]], serials[n])) for n in values[1:] if n != values[0])
    return tuple(segments), frozenset(bonds), {"unconnected_boundaries": breaks, "noncanonical_records_not_screened": sorted(skipped),
        "disulfide_policy": "explicit CONECT preserved; SSBOND/sidechain-surface topology not tested", "topology_qualified": False}


def run(request_path, source_path, candidate_path, output):
    request = DiffusionRequest.from_json(request_path.read_text())
    source_bytes = source_path.read_bytes()
    if hashlib.sha256(source_bytes).hexdigest() != request.normalized_pdb.sha256:
        raise ValueError("source digest differs from declared ownership")
    source = project(source_bytes.decode())
    candidate_bytes = candidate_path.read_bytes()
    candidate = project(candidate_bytes.decode())
    coordinates = {atom: xyz for segment in candidate[0] for atom, xyz in
                   ((segment.first, segment.start), (segment.second, segment.stop))}
    surfaces = []
    generated = set(request.generated_atoms)
    for gap_index, gap in enumerate(request.gaps):
        if gap.left_anchor is None or gap.right_anchor is None:
            surfaces.append({"gap_index": gap_index, "unsupported": "terminal gap has no two-anchor closure"})
            continue
        residues = (gap.left_anchor, *gap.generated_residues, gap.right_anchor)
        vertices = [coordinates[AtomIdentity(r.chain, r.residue_number, r.insertion_code, "CA")] for r in residues]
        own = {(r.chain, r.residue_number, r.insertion_code) for r in residues}
        polygon = np.asarray(vertices)
        selected = []
        for index, segment in enumerate(candidate[0]):
            if any((a.chain, a.residue_number, a.insertion_code) in own for a in (segment.first, segment.second)):
                continue
            endpoints = np.asarray((segment.start, segment.stop))
            if np.all(endpoints.max(axis=0) >= polygon.min(axis=0)-1e-5) and np.all(endpoints.min(axis=0) <= polygon.max(axis=0)+1e-5):
                selected.append(index)
        paths = tuple(candidate[0][i] for i in selected)
        result = asdict(probe_loop_surface(vertices, paths))
        for hit in result["piercings"]:
            index = selected[hit["segment_index"]]
            segment = candidate[0][index]
            hit["candidate_segment_index"] = index
            hit["path_ownership"] = "generated" if {segment.first, segment.second} & generated else "observed"
        surfaces.append({"gap_index": gap_index, "result": result, "own_residues_excluded": True,
                         "source_gap_has_no_surface_control": True})
    report = {"source_sha256": request.normalized_pdb.sha256, "candidate_sha256": hashlib.sha256(candidate_bytes).hexdigest(),
        "request_sha256": hashlib.sha256(request_path.read_bytes()).hexdigest(), "diagnostic_only": True,
        "limits": "not calibrated; open-path topology, loop closure ambiguity and trajectory crossing remain untested",
        "source": {"projection": source[2], "screen": asdict(screen_segments(source[0], frozenset(), bonded_atoms=source[1]))},
        "candidate": {"projection": candidate[2], "screen": asdict(screen_segments(candidate[0], frozenset(request.generated_atoms), bonded_atoms=candidate[1]))},
        "segments": [asdict(segment) for segment in candidate[0]], "loop_surface_probes": surfaces}
    serialized = json.dumps(report, indent=2) + "\n"
    with output.open("x") as stream:
        stream.write(serialized)
    print(json.dumps({"candidate_segment_alarms": len(report["candidate"]["screen"]["alarms"]),
                      "source_segment_alarms": len(report["source"]["screen"]["alarms"])}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("request", "source", "candidate", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    run(args.request, args.source, args.candidate, args.output)
