"""Execute synthetic single-chi authority transport; NOT real Dunbrack sampling.

No scientific database rows, held-out experimental coordinates or downloads are
used. The arbitrary smooth numbers below exercise parser, phi/psi lookup,
generation, immutable observations and independent evidence/geometry readback.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
from dataclasses import asdict
from pathlib import Path

from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef
from dvbfixer.protein_completion.backbone import backbone_angles
from dvbfixer.protein_completion.contracts import CompletionRequest, Observation
from dvbfixer.protein_completion.geometry import l_cb
from dvbfixer.protein_completion.rotamers import LibraryAuthority, load_single_chi
from dvbfixer.protein_completion.service import complete_sidechain
from dvbfixer.protein_completion.validation import audit_candidate


def synthetic_table(category):
    lines = ["# Authored synthetic transport control. NOT Dunbrack statistics.\n"]
    for phi in range(-180, 180, 10):
        for psi in range(-180, 180, 10):
            p1 = .4+.1*math.cos(math.radians(phi))
            p2 = .35+.05*math.sin(math.radians(psi))
            means = (60+3*math.sin(math.radians(phi)), 175+3*math.cos(math.radians(psi)),
                     -65+3*math.cos(math.radians(phi+psi)))
            for label, (probability, mean, sigma) in enumerate(zip((p1, p2, 1-p1-p2), means, (8, 10, 12)), 1):
                lines.append(f"{category} {phi} {psi} 0 {label} 0 0 0 {probability:.9f} {mean:.9f} 0 0 0 {sigma} 0 0 0\n")
    return "".join(lines).encode()


def request(name, with_cb):
    component = ComponentInstanceRef(1, "a", "82", "A", occurrence=2)
    positions = {"N": (-1.46, 0., 0.), "CA": (0., 0., 0.),
                 "C": (0.55, 1.43, 0.), "O": (1.76, 1.66, 0.)}
    if with_cb:
        positions["CB"] = l_cb(positions["N"], positions["CA"], positions["C"])
    observed = tuple(Observation(ExactAtomRef(component, n), n[0], p) for n, p in positions.items())
    prev = Observation(ExactAtomRef(ComponentInstanceRef(1, "a", "81"), "C"), "C", (-2.1, -0.7, 0.8))
    after = Observation(ExactAtomRef(ComponentInstanceRef(1, "a", "83"), "N"), "N", (0., 2.15, -0.9))
    observed += (prev, after)
    source_digest = hashlib.sha256(json.dumps([asdict(o) for o in observed], sort_keys=True).encode()).hexdigest()
    last = "OG" if name == "SER" else "SG"
    return CompletionRequest(source_digest, component, name, observed,
                             (last,) if with_cb else ("CB", last),
                             ((prev.atom, ExactAtomRef(component, "N")),
                              (ExactAtomRef(component, "C"), after.atom)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path)
    parser.add_argument("--authority-probe", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    notice = b"Authored numerical integration controls, not externally licensed Dunbrack data.\n"
    (args.output/"synthetic-notice.txt").write_bytes(notice)
    root = Path(__file__).resolve().parents[1]
    report = {"schema": "dvbfixer.synthetic-rotamer-transport.v1", "python": platform.python_version(),
              "real_dunbrack_statistics_loaded": False, "public_backend_enabled": False,
              "qualifies_empirical_library": False,
              "source_hashes": {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
                                for p in [Path(__file__).resolve(), *sorted((root/"src/dvbfixer/protein_completion").glob("*.py"))]},
              "libraries": {}, "cases": {}}
    if args.authority_probe:
        report["authority_probe"] = json.loads(args.authority_probe.read_text())
        report["authority_probe_sha256"] = hashlib.sha256(args.authority_probe.read_bytes()).hexdigest()
    for category, residue in (("SER", "SER"), ("CYH", "CYS")):
        data = synthetic_table(category)
        path = args.output/f"synthetic-{category}.lib"
        path.write_bytes(data)
        authority = LibraryAuthority("authored-integration-control-v1", category,
                                     "scripts/probe_rotamer_transport.py, no experimental statistics",
                                     "synthetic-control", "synthetic-test-not-dunbrack",
                                     hashlib.sha256(notice).hexdigest())
        library = load_single_chi(path, expected_sha256=hashlib.sha256(data).hexdigest(), authority=authority)
        report["libraries"][category] = {"authority": asdict(authority), "data_sha256": library.data_sha256,
                                          "view_sha256": library.view_sha256, "model_digest": library.model_digest,
                                          "cells": len(library.cells)}
        for with_cb in (False, True):
            req = request(residue, with_cb)
            before = req.digest
            result = complete_sidechain(req, library=library)
            if result.status != "completed" or req.digest != before:
                raise RuntimeError("synthetic completion or observation binding failed")
            for attempt in result.attempts:
                if audit_candidate(req, attempt.candidate, library=library) != attempt.audit:
                    raise RuntimeError("independent synthetic readback disagrees")
            report["cases"][f"{residue}-{'terminal' if with_cb else 'whole-sidechain'}"] = {
                "request": asdict(req), "backbone_angles": backbone_angles(req),
                "result": asdict(result), "observations_digest_unchanged": req.digest == before}
    evidence = args.output/"evidence.json"
    evidence.write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
    if args.summary:
        portable = report | {"raw_record_sha256": hashlib.sha256(evidence.read_bytes()).hexdigest(),
                             "retained_directory": str(args.output)}
        with args.summary.open("x") as handle:
            json.dump(portable, handle, indent=2, allow_nan=False)
            handle.write("\n")
    print(json.dumps({k: v["result"]["status"] for k, v in report["cases"].items()}))


if __name__ == "__main__":
    main()
