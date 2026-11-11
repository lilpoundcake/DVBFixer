#!/usr/bin/env python3
"""Bind native RFD3/Boltz output identities to a declared reference-case request.

Only observed coordinates enter fitting. Preserve native drift before source
restoration and store full precision arrays. Reference coordinates are not read.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import re
from pathlib import Path

import numpy as np
from scripts.audit_mac_reconstruction_smoke import inventory

from dvbfixer.model.diffusion.contract import AtomIdentity, DiffusionRequest
from dvbfixer.model.diffusion.geometry import weighted_kabsch
from dvbfixer.model.diffusion.pdb_materialize import materialize_candidate_pdb


def run(case, prediction, output, *, model):
    request_path = case / "evaluation-mask.json"
    request = DiffusionRequest.from_json(request_path.read_text())
    source_bytes = (case / "source.pdb").read_bytes()
    if hashlib.sha256(source_bytes).hexdigest() != request.normalized_pdb.sha256:
        raise ValueError("reference-case source changed")
    manifest = json.loads((case / "rfd3/manifest.json").read_text())
    native = inventory(prediction)
    source = inventory(case / "source.pdb")
    source_by_index = None
    if model == "rfd3":
        metadata_path = prediction.with_suffix("").with_suffix(".json")
        index_map = json.loads(metadata_path.read_text())["diffused_index_map"]
        mapping_custody = {"metadata_sha256":hashlib.sha256(metadata_path.read_bytes()).hexdigest()}
    else:
        helper = Path(__file__).parents[1] / "deploy/boltz-2/atom_mapping.py"
        spec = importlib.util.spec_from_file_location("_boltz_reference_axis", helper)
        if spec is None or spec.loader is None:
            raise ValueError("Boltz canonical atom-mapping helper unavailable")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        source_by_index = module.target_residue_map(request)
        index_map = {f"{r.chain}{r.residue_number}":f"{request.target_sequences[0].chain}{i+1}" for i,r in source_by_index.items()}
        mapping_custody = {"canonical_boltz_mapping_helper_sha256":hashlib.sha256(helper.read_bytes()).hexdigest(),
                           "native_export_convention":"single-chain full target, one-based residue ordinal; complete inventory/sequence checked"}
    axis, names, elements, coordinates = [], [], [], []
    represented = set()
    for item in manifest["axis"]:
        atom = AtomIdentity(**item["atom"])
        label = f"{atom.chain}{atom.residue_number}"
        if label not in index_map:
            raise ValueError("native identity mapping omits a requested source residue")
        match = re.fullmatch(r"(.)(-?\d+)([A-Za-z]?)", index_map[label])
        if match is None:
            raise ValueError("unsupported native residue identity")
        chain, number, code = match.groups()
        mapped = AtomIdentity(chain,number,code,atom.atom_name)
        if mapped in represented or mapped not in native or native[mapped][0] != item["residue_name"]:
            raise ValueError("native source/sequence/atom identity is ambiguous")
        represented.add(mapped)
        axis.append(atom)
        names.append(item["residue_name"])
        elements.append(native[mapped][1])
        coordinates.append(native[mapped][2])
    if represented != native.keys() or set(axis) != set(request.fixed_atoms) | set(request.generated_atoms):
        raise ValueError("native atom inventory differs from complete reference-case mask")
    raw = np.asarray(coordinates)
    index = {atom:i for i,atom in enumerate(axis)}
    fixed = tuple(sorted(source))
    transform = weighted_kabsch(np.array([raw[index[a]] for a in fixed]), np.array([source[a][2] for a in fixed]))
    aligned = transform.apply(raw)
    drift = np.array([np.linalg.norm(aligned[index[a]]-source[a][2]) for a in fixed])
    # Retain native observed drift here; source restoration happens separately in
    # the independent audit, never hidden inside this remapping operation.
    text = materialize_candidate_pdb(source_bytes.decode(), tuple(axis), aligned, tuple(names), request, elements=tuple(elements))
    output.mkdir(exist_ok=False)
    (output / "aligned-native.pdb").write_text(text)
    np.savez_compressed(output / "native-arrays.npz", raw=raw,aligned=aligned,rotation=transform.rotation,
                        translation=transform.translation,atom_name=np.array([a.atom_name for a in axis]),
                        residue_number=np.array([a.residue_number for a in axis]),chain=np.array([a.chain for a in axis]))
    report = {"model":model, "prediction_sha256":hashlib.sha256(prediction.read_bytes()).hexdigest(),
        "case_request_sha256":hashlib.sha256(request_path.read_bytes()).hexdigest(), "source_sha256":request.normalized_pdb.sha256,
        "mapping_custody":mapping_custody, "native_observed_rms_angstrom":float(np.sqrt(np.mean(drift**2))),
        "native_observed_max_angstrom":float(drift.max()), "source_restored":False,"reference_coordinates_read":False,
        "qualification":False}
    (output / "mapping.json").write_text(json.dumps(report,indent=2)+"\n")
    return output / "aligned-native.pdb"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", type=Path)
    parser.add_argument("prediction", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--model", choices=("rfd3","boltz"), required=True)
    args = parser.parse_args()
    print(run(args.case,args.prediction,args.output,model=args.model))
