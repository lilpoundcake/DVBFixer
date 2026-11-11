#!/usr/bin/env python3
"""New source-bound complete-reference projection, with hidden coordinates withheld.

Explicit projection endpoints and hidden numeric intervals are declared before
sampling. This does not repair source geometry, rewrite occupancies, use a greedy
FASTA placement, or modify an existing frozen request.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import asdict
from pathlib import Path

import numpy as np
from scripts.prepare_rfd3_local_projection import prepare as prepare_rfd3

from dvbfixer.diagnose.report import Severity
from dvbfixer.ffutils.geometry import assert_all_l
from dvbfixer.model.diffusion.anchor_acceptance import _structure
from dvbfixer.model.diffusion.input_geometry import InputGeometryError, assess_input_geometry
from dvbfixer.model.diffusion.protein_atoms import _CANONICAL_ONE_TO_THREE, CANONICAL_HEAVY_ATOMS
from dvbfixer.model.diffusion.staged_runtime import assessment_adapter
from dvbfixer.model.diffusion_cli import build_cli_diffusion_request


def prepare(reference, output, *, chain, first, last, hidden_first, hidden_last, case_id):
    if len(chain) != 1 or not first < hidden_first <= hidden_last < last:
        raise ValueError("declare a single-case-sensitive chain and an internal hidden interval")
    if last-first+1 > 256:
        raise ValueError("this first reference panel is bounded to 256 residues")
    original = reference.read_bytes()
    records, names = {}, {}
    for line in original.decode().splitlines(True):
        if not line.startswith("ATOM  ") or line[21] != chain or not first <= int(line[22:26]) <= last:
            continue
        number, atom_name, residue_name = int(line[22:26]), line[12:16].strip(), line[17:20].strip()
        if line[16].strip() or line[26].strip() or float(line[54:60]) != 1:
            raise ValueError("reference projection contains alternatives, insertion codes or partial occupancy")
        if (number, atom_name) in records or residue_name not in CANONICAL_HEAVY_ATOMS:
            raise ValueError("reference projection has duplicate/noncanonical identities")
        if not np.isfinite([float(line[i:i+8]) for i in (30,38,46)]).all():
            raise ValueError("reference projection has nonfinite coordinates")
        if number in names and names[number] != residue_name:
            raise ValueError("reference residue chemistry is ambiguous")
        names[number] = residue_name
        records[number, atom_name] = line
    if set(names) != set(range(first, last+1)):
        raise ValueError("reference interval is not complete; do not invent native coordinates")
    for number, name in names.items():
        if {atom for residue,atom in records if residue==number} != set(CANONICAL_HEAVY_ATOMS[name]):
            raise ValueError("reference has missing/extra canonical heavy atoms or measured hydrogens")
    symbols = {name:symbol for symbol,name in _CANONICAL_ONE_TO_THREE.items()}
    sequence = "".join(symbols[names[i]] for i in range(first,last+1))
    header = f"REMARK 999 EXPLICIT {case_id} {chain}{first}-{last} REFERENCE PROJECTION\n"
    native = header+"".join(records.values())+"TER\nEND\n"
    structure = _structure(native, "reference-projection.pdb")
    assert_all_l(structure.pdb.topology, structure.pdb.positions)
    source = header+"".join(line for (number,_),line in records.items() if not hidden_first<=number<=hidden_last)+"TER\nEND\n"
    output.mkdir(exist_ok=False, parents=True)
    (output / "reference-entry-NOT-CONDITIONING.pdb").write_bytes(original)
    (output / "reference-projection-NOT-CONDITIONING.pdb").write_text(native)
    source_path = output / "source.pdb"
    source_path.write_text(source)
    request = build_cli_diffusion_request(source_path, {chain:sequence}, seeds=(7,), profile="protpardelle-1c-cc94-mps-v1")
    if len(request.gaps) != 1 or {int(r.residue_number) for r in request.gaps[0].generated_residues} != set(range(hidden_first,hidden_last+1)):
        raise ValueError("shared sequence placement differs from the prospectively hidden interval")
    request_path = output / "evaluation-mask.json"
    request_path.write_text(request.to_json())
    findings = assess_input_geometry(request, source.encode())
    errors = [f for f in findings if f.severity is Severity.ERROR]
    native_assessment = assessment_adapter(request, source)(native)
    native_accepted = native_assessment.geometry_passed and not native_assessment.admission_failures
    (output / "native-reference-assessment.json").write_text(json.dumps(asdict(native_assessment),indent=2)+"\n")
    (output / "source-geometry.json").write_text(json.dumps({"findings":[asdict(f) for f in findings], "errors":len(errors)},
                                                         default=lambda value:value.value, indent=2)+"\n")
    manifest = {"case_id":case_id, "chain":chain, "projection":{"first":first,"last":last},
        "hidden_interval":{"first":hidden_first,"last":hidden_last}, "sequence":sequence,
        "reference_entry_sha256":hashlib.sha256(original).hexdigest(), "generated_atoms":len(request.generated_atoms),
        "observed_atoms":len(request.fixed_atoms), "source_admissible":not errors, "source_repaired":False,
        "reference_baseline_accepted":native_accepted,
        "reference_coords_supplied_to_sampling":False, "new_projection_not_frozen_replay":True,
        "files":{name:hashlib.sha256((output/name).read_bytes()).hexdigest() for name in
                 ("source.pdb","reference-projection-NOT-CONDITIONING.pdb","evaluation-mask.json")}}
    (output / "manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
    if errors:
        raise InputGeometryError(f"reference-case source has {len(errors)} severe findings; bytes and diagnostics retained, no sampling input published")
    if not native_accepted:
        raise InputGeometryError("native reference projection fails baseline geometry/custody; evidence retained, no sampling input published")
    prepare_rfd3(request_path, source_path, (0,), output / "rfd3", flanks=len(sequence))
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--chain", required=True)
    for name in ("first", "last", "hidden-first", "hidden-last"):
        parser.add_argument("--"+name, type=int, required=True)
    parser.add_argument("--case-id", required=True)
    args = parser.parse_args()
    print(json.dumps(prepare(args.reference,args.output,chain=args.chain,first=args.first,last=args.last,
        hidden_first=args.hidden_first,hidden_last=args.hidden_last,case_id=args.case_id),indent=2))
