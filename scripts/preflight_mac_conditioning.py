#!/usr/bin/env python3
"""Pinned RFD3 input-builder readback without loading model weights."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def run(inputs: Path, output: Path):
    from rfd3.inference.input_parsing import DesignInputSpecification

    manifest = json.loads((inputs / "manifest.json").read_text())
    spec = json.loads((inputs / "rfd3-input.json").read_text())
    if len(spec) != 1:
        raise ValueError("single declared conditioning input required")
    array = DesignInputSpecification(**next(iter(spec.values()))).build()
    fixed = np.asarray(array.is_motif_atom_with_fixed_coord, dtype=bool)
    sequence = np.asarray(array.is_motif_atom_with_fixed_seq, dtype=bool)
    if len(array) != len(manifest["axis"]) or fixed.sum() != manifest["sampler_fixed_atoms"] or not sequence.all():
        raise ValueError("pinned builder mask/sequence differs from declaration")
    if not np.all(array.coord[~fixed] == 0):
        raise ValueError("pinned builder retained coordinates on unfixed conditioning atoms")
    # Numeric source IDs are transported in src_component, not native output IDs.
    actual = {}
    for component, atom_name, coord, mask in zip(array.src_component, array.atom_name, array.coord, fixed):
        key = str(component).strip(), str(atom_name).strip()
        if key in actual:
            raise ValueError("pinned builder duplicated source atom identity")
        actual[key] = (coord, bool(mask))
    for row in manifest["axis"]:
        atom = row["atom"]
        key = atom["chain"]+atom["residue_number"], atom["atom_name"]
        if key not in actual or actual[key][1] != row["sampler_fixed"]:
            raise ValueError("pinned builder changed an atom-wise conditioning mask")
    output.mkdir(parents=True, exist_ok=False)
    np.savez_compressed(output / "built-input.npz", fixed=fixed, sequence_fixed=sequence, coordinates=array.coord,
        source_component=array.src_component, atom_name=array.atom_name)
    receipt = {"atoms": len(array), "sampler_fixed": int(fixed.sum()), "sampler_unfixed": int((~fixed).sum()),
        "exact_source_identity_masks_checked": True, "all_sequence_fixed": bool(sequence.all()), "checkpoint_loaded": False,
        "unfixed_coordinates_zero": True, "manifest_sha256": hashlib.sha256((inputs / "manifest.json").read_bytes()).hexdigest(),
        "spec_sha256": hashlib.sha256((inputs / "rfd3-input.json").read_bytes()).hexdigest(),
        "sampling": False, "information": manifest["information"], "source_target_spans": manifest["source_target_spans"]}
    (output / "preflight.json").write_text(json.dumps(receipt, indent=2)+"\n")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    run(args.inputs.resolve(), args.output.resolve())
