#!/usr/bin/env python3
"""Private single-chain, whole-gap transport for admitted expansion masks.

The shared smoke builders intentionally admit one gap. This separately labeled
adapter handles a declared atomic pair without changing those public/frozen
builders. It supplies observed coordinates only and fixes the entire sequence.
Checkpoint consumer preflight/execution remains a separate gate.
"""
from __future__ import annotations

import argparse
import json
from dataclasses import replace
from pathlib import Path

import gemmi
from scripts.audit_mac_reconstruction_smoke import inventory
from scripts.mac_reconstruction_comparison import digest, save, target_axis
from scripts.prepare_rfd3_local_projection import prepare as prepare_rfd3

from dvbfixer.model.diffusion.contract import DiffusionRequest
from dvbfixer.model.diffusion.protein_atoms import _CANONICAL_ONE_TO_THREE


def prepare(case: Path, output: Path):
    manifest = json.loads((case / "manifest.json").read_text())
    if not manifest["source_admissible"] or not manifest["reference_baseline_accepted"]:
        raise ValueError("case did not pass source/reference admission")
    for name, expected in manifest["files"].items():
        if digest(case / name) != expected:
            raise ValueError("prepared expansion bytes changed")
    request = DiffusionRequest.from_json((case / "evaluation-mask.json").read_text())
    if len(request.target_sequences) != 1 or request.retained_explicit_links or request.heterogen_contexts:
        raise ValueError("private transport admits single-chain canonical protein only")
    target_axis(request)
    if set(inventory(case / "source.pdb")) != set(request.fixed_atoms):
        raise ValueError("source inventory differs from fixed mask")
    request = replace(request, seeds=(7,), candidate_count=1)
    output.mkdir(parents=True, exist_ok=False)
    for name in ("source.pdb", "reference-projection-NOT-CONDITIONING.pdb"):
        (output / name).write_bytes((case / name).read_bytes())
    (output / "evaluation-mask.json").write_text(request.to_json())
    target, placement = request.target_sequences[0], request.sequence_placements[0]
    manifest = {**manifest, "chain": target.chain, "sequence": target.sequence,
        "preparation_manifest_sha256": digest(case / "manifest.json"),
        "files": {name: digest(output / name) for name in ("source.pdb", "reference-projection-NOT-CONDITIONING.pdb", "evaluation-mask.json")}}
    save(output / "manifest.json", manifest)
    prepare_rfd3(output / "evaluation-mask.json", output / "source.pdb", tuple(range(len(request.gaps))),
        output / "rfd3", flanks=len(target.sequence))
    # No held-out coordinates are read by either template encoder below.
    boltz = output / "boltz-2"
    boltz.mkdir()
    names = [_CANONICAL_ONE_TO_THREE[symbol] for symbol in target.sequence]
    seqres = [f"SEQRES {row:3d} {target.chain} {len(names):4d}  {' '.join(names[start:start+13])}\n"
              for row, start in enumerate(range(0, len(names), 13), 1)]
    source = (output / "source.pdb").read_text()
    (boltz / "template.pdb").write_text("".join(seqres)+"".join(line for line in source.splitlines(keepends=True) if not line.startswith("SEQRES")))
    save(boltz / "input.yaml", {"version": 1, "sequences": [{"protein": {"id": target.chain, "sequence": target.sequence, "msa": "empty"}}],
        "templates": [{"pdb": str(boltz / "template.pdb"), "chain_id": target.chain, "template_id": target.chain+"1"}]})
    protenix = output / "protenix-v1"
    protenix.mkdir()
    structure = gemmi.read_structure(str(output / "source.pdb"))
    structure.setup_entities()
    save(protenix / "template.json", [{"mmcif": structure.make_mmcif_document().as_string(),
        "queryIndices": list(placement.observed_target_indices), "templateIndices": list(range(len(placement.observed_residues)))}])
    save(protenix / "input.json", [{"name": "private_mac_expansion", "sequences": [{"proteinChain": {
        "sequence": target.sequence, "count": 1, "templatesPath": str(protenix / "template.json")}}], "modelSeeds": [7]}])
    receipt = {"case": case.name, "transport": "private-single-chain-whole-gap-expansion-v1",
        "atomic_group": list(range(len(request.gaps))), "seed": 7,
        "fixed_atoms": len(request.fixed_atoms), "generated_atoms": len(request.generated_atoms),
        "reference_coordinates_encoded": False, "remote_msa": False, "checkpoint_loaded": False,
        "consumer_axis_preflight": "required before sampling", "public_capability_enabled": False}
    save(output / "transport.json", receipt)
    return receipt


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    print(json.dumps(prepare(args.case.resolve(), args.output.resolve()), indent=2))
