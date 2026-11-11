#!/usr/bin/env python3
"""Pinned Boltz data/token/feature-axis check, with no weights or remote MSA."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def run(case: Path, cache: Path, output: Path):
    from atom_mapping import model_atom_axis, partition_requested_atoms, validate_feature_axis
    from boltz.data.module.inferencev2 import PredictionDataset, load_input
    from boltz.data.tokenize.boltz2 import Boltz2Tokenizer
    from boltz.data.types import Manifest
    from boltz.main import process_inputs
    from rdkit import Chem

    from dvbfixer.model.diffusion.contract import DiffusionRequest

    request = DiffusionRequest.from_json((case / "evaluation-mask.json").read_text())
    output.mkdir(parents=True, exist_ok=False)
    Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)
    process_inputs(data=[case / "boltz-2/input.yaml"], out_dir=output / "boltz", ccd_path=cache / "ccd.pkl", mol_dir=cache / "mols",
        msa_server_url="https://api.colabfold.com", msa_pairing_strategy="greedy", boltz2=True,
        preprocessing_threads=1, use_msa_server=False)
    processed = output / "boltz/processed"
    manifest = Manifest.load(processed / "manifest.json")
    if len(manifest.records) != 1:
        raise ValueError("preflight expected exactly one target")
    record = manifest.records[0]
    paths = [processed / p for p in ("structures", "msa", "constraints", "templates", "mols")]
    tokenized = Boltz2Tokenizer().tokenize(load_input(record, *paths))
    axis, _, tokens = model_atom_axis(tokenized, request, chain_name_by_asym_id={c.chain_id: c.chain_name for c in record.chains})
    represented, passthrough = partition_requested_atoms(request, axis)
    features = PredictionDataset(manifest, paths[0], paths[1], cache / "mols", *paths[2:])[0]
    validate_feature_axis(axis, tokens, features)
    if any(a.atom_name != "OXT" for a in passthrough) or not set(request.generated_atoms) <= set(axis):
        raise ValueError("unexpected missing native observation or generated atom")
    (output / "preflight.json").write_text(json.dumps({"real_atoms": len(axis), "represented_fixed": len(represented),
        "request_sha256": hashlib.sha256((case / "evaluation-mask.json").read_bytes()).hexdigest(),
        "observed_terminal_passthrough": len(passthrough), "exact_token_atom_feature_axis_verified": True,
        "checkpoint_loaded": False, "remote_msa_used": False, "reserved_seeds_sampled": False}, indent=2)+"\n")
    print(json.dumps({"case": case.name, "real_atoms": len(axis), "checkpoint_loaded": False}), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("case", "cache", "output"):
        parser.add_argument("--"+name, required=True, type=Path)
    args = parser.parse_args()
    run(args.case.resolve(), args.cache.resolve(), args.output.resolve())
