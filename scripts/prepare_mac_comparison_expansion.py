#!/usr/bin/env python3
"""Freeze new terminal/pair/interface reference masks without neural sampling.

No source repair, occupancy edits, automatic case replacement or public admission
bypass. Scientific readiness and backend transport support are separate outputs.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import traceback
from dataclasses import asdict
from pathlib import Path

import numpy as np
from scripts.mac_reconstruction_comparison import DEVELOPMENT_SEEDS, digest, save

from dvbfixer.diagnose.report import Severity
from dvbfixer.ffutils.geometry import assert_all_l
from dvbfixer.model.diffusion.anchor_acceptance import _structure
from dvbfixer.model.diffusion.contract import (
    DIFFUSION_SCHEMA_VERSION,
    ArtifactReference,
    AtomIdentity,
    DiffusionRequest,
    GapKind,
    GapRegion,
    ResidueIdentity,
    SequencePlacement,
    TargetInterval,
    TargetSequence,
)
from dvbfixer.model.diffusion.input_geometry import assess_input_geometry
from dvbfixer.model.diffusion.protein_atoms import _CANONICAL_ONE_TO_THREE, CANONICAL_HEAVY_ATOMS
from dvbfixer.model.diffusion.staged_runtime import assessment_adapter

PROJECT = Path(__file__).resolve().parents[1]
CASES = (
    {"id": "1UBQ-A1-72-Nterm1-4", "stratum": "N-terminal", "fixture": "tests/fixtures/protein_context/1UBQ.pdb",
     "spans": [["A", 1, 72]], "hidden": [["A", 1, 4]]},
    {"id": "1UBQ-A1-72-pair10-13_46-49", "stratum": "two-gap candidate; interaction not assumed",
     "fixture": "tests/fixtures/protein_context/1UBQ.pdb", "spans": [["A", 1, 72]], "hidden": [["A", 10, 13], ["A", 46, 49]]},
    {"id": "1A3N-A1-59-B1-146-gap30-37", "stratum": "protein-partner candidate; interface proximity must be reported",
     "fixture": "tests/fixtures/protein_context/1A3N.pdb", "spans": [["A", 1, 59], ["B", 1, 146]], "hidden": [["A", 30, 37]]},
)


def prepare(spec: dict, output: Path, *, allow_observed_terminal_oxt=False):
    output.mkdir(parents=True, exist_ok=False)
    reference = PROJECT / spec["fixture"]
    save(output / "prospective-case.json", {**spec, "fixture_sha256": digest(reference), "sampling": False})
    spans = {chain: (first, last) for chain, first, last in spec["spans"]}
    if len(spans) != len(spec["spans"]) or sum(last-first+1 for first, last in spans.values()) > 256:
        raise ValueError("declare unique case-sensitive chains and <=256 total residues")
    hidden = {(chain, str(number), "") for chain, first, last in spec["hidden"] for number in range(first, last+1)}
    records, names, atoms = {}, {}, {}
    for line in reference.read_text().splitlines(keepends=True):
        if not line.startswith("ATOM  ") or line[21] not in spans:
            continue
        chain = line[21]
        number = int(line[22:26])
        if not spans[chain][0] <= number <= spans[chain][1]:
            continue
        residue = chain, str(number), line[26].strip()
        atom_name, residue_name = line[12:16].strip(), line[17:20].strip()
        if line[16].strip() or residue[2] or float(line[54:60]) != 1:
            raise ValueError(f"ambiguous/partial occupancy at {residue}; source not rewritten")
        if residue_name not in CANONICAL_HEAVY_ATOMS:
            raise ValueError("noncanonical source")
        identity = AtomIdentity(*residue, atom_name)
        if identity in records:
            raise ValueError("duplicate source atom")
        if not np.isfinite([float(line[i:i+8]) for i in (30, 38, 46)]).all():
            raise ValueError("nonfinite source coordinate")
        if residue in names and names[residue] != residue_name:
            raise ValueError("ambiguous residue chemistry")
        names[residue] = residue_name
        records[identity] = line
        atoms.setdefault(residue, set()).add(atom_name)
    expected_residues = {(chain, str(i), "") for chain, (first, last) in spans.items() for i in range(first, last+1)}
    if set(names) != expected_residues or not hidden < expected_residues:
        raise ValueError("incomplete reference interval or invalid mask; no native coordinates invented")
    for residue, name in names.items():
        expected_atoms = set(CANONICAL_HEAVY_ATOMS[name])
        if (allow_observed_terminal_oxt and residue[1] == str(spans[residue[0]][1])
                and residue not in hidden and "OXT" in atoms[residue]):
            expected_atoms.add("OXT")
        if atoms[residue] != expected_atoms:
            raise ValueError(f"incomplete/extra canonical heavy atoms at {residue}; no source completion")
    symbols = {name: symbol for symbol, name in _CANONICAL_ONE_TO_THREE.items()}
    targets, placements, gaps = [], [], []
    for chain, (first, last) in spans.items():
        ordered = tuple(ResidueIdentity(chain, str(i), "") for i in range(first, last+1))
        sequence = "".join(symbols[names[chain, r.residue_number, ""]] for r in ordered)
        targets.append(TargetSequence(chain, sequence))
        observed = tuple(i for i, r in enumerate(ordered) if (r.chain, r.residue_number, "") not in hidden)
        placements.append(SequencePlacement(chain, len(sequence), observed, tuple(ordered[i] for i in observed)))
        for mask_chain, start, stop in spec["hidden"]:
            if mask_chain != chain:
                continue
            left, right = start-first, stop-first+1
            if not 0 <= left < right <= len(ordered) or right-left == len(ordered):
                raise ValueError("invalid hidden interval")
            left_anchor = ordered[left-1] if left else None
            right_anchor = ordered[right] if right < len(ordered) else None
            generated = ordered[left:right]
            kind = GapKind.N_TERMINAL if left_anchor is None else GapKind.C_TERMINAL if right_anchor is None else GapKind.INTERNAL
            gaps.append(GapRegion(chain, TargetInterval(left, right), left_anchor, right_anchor, generated,
                tuple(r for r in (left_anchor, *generated, right_anchor) if r is not None), kind))
    header = f"REMARK 999 EXPLICIT PRIVATE {spec['id']} PROJECTION; NOT FULL ENTRY\n"
    def materialize(include_hidden):
        lines = [header]
        for chain in spans:
            lines.extend(line for atom, line in records.items() if atom.chain == chain and
                         (include_hidden or (atom.chain, atom.residue_number, atom.insertion_code) not in hidden))
            lines.append("TER\n")
        return "".join(lines)+"END\n"
    native, source = materialize(True), materialize(False)
    (output / "source.pdb").write_text(source)
    (output / "reference-projection-NOT-CONDITIONING.pdb").write_text(native)
    fixed = tuple(a for a in records if (a.chain, a.residue_number, a.insertion_code) not in hidden)
    generated = tuple(a for a in records if (a.chain, a.residue_number, a.insertion_code) in hidden)
    request = DiffusionRequest(schema_version=DIFFUSION_SCHEMA_VERSION,
        normalized_pdb=ArtifactReference("source.pdb", hashlib.sha256(source.encode()).hexdigest()),
        target_sequences=tuple(targets), sequence_placements=tuple(placements), gaps=tuple(gaps), fixed_atoms=fixed,
        generated_atoms=generated, retained_explicit_links=(), candidate_count=6, seeds=DEVELOPMENT_SEEDS)
    (output / "evaluation-mask.json").write_text(request.to_json())
    (output / "target.fasta").write_text("".join(f">{t.chain}\n{t.sequence}\n" for t in targets))
    findings = assess_input_geometry(request, source.encode())
    baseline = assessment_adapter(request, source)(native)
    structure = _structure(native, "reference.pdb")
    assert_all_l(structure.pdb.topology, structure.pdb.positions)
    errors = [f for f in findings if f.severity is Severity.ERROR]
    save(output / "source-geometry.json", [asdict(f) | {"severity": f.severity.value} for f in findings])
    save(output / "reference-baseline.json", asdict(baseline))
    # This is mask preparation, NOT an override of public backend admission.
    manifest = {"case_id": spec["id"], "sequence_by_chain": {t.chain: t.sequence for t in targets},
        "source_admissible": not errors, "reference_baseline_accepted": baseline.geometry_passed and not baseline.admission_failures,
        "observed_atoms": len(fixed), "generated_atoms": len(generated), "source_repaired": False,
        "observed_terminal_oxt_preserved": [dict(chain=a.chain, residue_number=a.residue_number, atom_name=a.atom_name)
            for a in fixed if a.atom_name == "OXT"],
        "public_capability_enabled": False, "backend_transport_preflight": "not executed; unsupported cells must remain explicit",
        "reference_coords_supplied_to_sampling": False, "sampling_executed": False,
        "files": {name: digest(output / name) for name in ("source.pdb", "reference-projection-NOT-CONDITIONING.pdb", "evaluation-mask.json")}}
    save(output / "manifest.json", manifest)
    if errors or not manifest["reference_baseline_accepted"]:
        raise ValueError("source/reference baseline failed; retained mask is not sampling-ready")
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    save(output / "protocol.json", {"purpose": "prospective new-mask admission only; zero sampling calls", "cases": CASES,
        "source_sha256": {s["fixture"]: digest(PROJECT / s["fixture"]) for s in CASES},
        "script_sha256": digest(Path(__file__)), "no_case_replacement": True, "qualification": False})
    (output / "preparation-script.py").write_bytes(Path(__file__).read_bytes())
    rows = []
    for spec in CASES:
        try:
            manifest = prepare(spec, output / spec["id"])
            rows.append({"case": spec["id"], "mask_admitted": True, "manifest": manifest})
        except Exception:
            error = traceback.format_exc()
            save(output / spec["id"] / "refusal.json", {"traceback": error, "sampling_executed": False})
            rows.append({"case": spec["id"], "mask_admitted": False, "error": error})
    save(output / "results.json", rows)
    print(json.dumps([{key: value for key, value in row.items() if key != "manifest"} for row in rows], indent=2))


if __name__ == "__main__":
    main()
