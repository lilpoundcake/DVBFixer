"""Private independent-cohort inputs and native transport with observed OXT custody."""
from __future__ import annotations

import json
import re
from dataclasses import replace
from pathlib import Path

import gemmi
import numpy as np
from scripts.audit_mac_reconstruction_smoke import inventory
from scripts.mac_conditioning_contracts import prepare as prepare_rfd3
from scripts.mac_reconstruction_comparison import digest, save, target_axis
from scripts.mac_validation_policy import SEEDS

from dvbfixer.model.diffusion.contract import AtomIdentity, DiffusionRequest
from dvbfixer.model.diffusion.geometry import weighted_kabsch
from dvbfixer.model.diffusion.pdb_materialize import materialize_candidate_pdb
from dvbfixer.model.diffusion.protein_atoms import _CANONICAL_ONE_TO_THREE


def prepare(case: Path, output: Path):
    source_manifest = json.loads((case / "manifest.json").read_text())
    if not source_manifest["source_admissible"] or not source_manifest["reference_baseline_accepted"]:
        raise ValueError("independent case source/reference baseline did not pass")
    for name, expected in source_manifest["files"].items():
        if digest(case / name) != expected:
            raise ValueError("independent case bytes changed")
    request = DiffusionRequest.from_json((case / "evaluation-mask.json").read_text())
    if len(request.target_sequences) != 1 or len(request.gaps) != 1:
        raise ValueError("fixed validation slice requires a single target/gap")
    output.mkdir(parents=True, exist_ok=False)
    for name in ("source.pdb", "reference-projection-NOT-CONDITIONING.pdb"):
        (output / name).write_bytes((case / name).read_bytes())
    request = replace(request, seeds=SEEDS, candidate_count=len(SEEDS))
    (output / "evaluation-mask.json").write_text(request.to_json())
    receipt = prepare_rfd3(output, output / "rfd3", gaps=(0,), context="full", information="all-heavy")
    represented = {AtomIdentity(**row["atom"]) for row in receipt["axis"]}
    passthrough = set(request.fixed_atoms)-represented
    ordered = target_axis(request)[request.target_sequences[0].chain]
    last = ordered[-1]
    expected_passthrough = {a for a in request.fixed_atoms if a.atom_name == "OXT" and
                           (a.chain, a.residue_number, a.insertion_code) == (last.chain, last.residue_number, last.insertion_code)}
    if passthrough != expected_passthrough:
        raise ValueError("only an already-observed true-terminal OXT may be outside the canonical RFD3/Boltz axis")
    target, placement = request.target_sequences[0], request.sequence_placements[0]
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
    save(protenix / "input.json", [{"name": "private_fixed_source_validation", "sequences": [{"proteinChain": {
        "sequence": target.sequence, "count": 1, "templatesPath": str(protenix / "template.json")}}], "modelSeeds": list(SEEDS)}])
    save(output / "transport.json", {"source_manifest_sha256": digest(case / "manifest.json"),
        "fixed_observed_terminal_passthrough": [{"chain": a.chain, "residue_number": a.residue_number, "atom_name": a.atom_name} for a in sorted(passthrough)],
        "passthrough_is_not_model_constraint_or_prediction": True, "original_OXT_bytes_kept_in_source_and_final": True,
        "protenix_expected_terminal_OXT_represented": True, "all_other_source_heavy_atoms_conditioned": True,
        "reference_coordinates_encoded": False, "checkpoint_loaded": False, "public_capability_enabled": False})
    return receipt


def transport(case: Path, prediction: Path, output: Path, *, model: str):
    request = DiffusionRequest.from_json((case / "evaluation-mask.json").read_text())
    source = inventory(case / "source.pdb")
    if digest(case / "source.pdb") != request.normalized_pdb.sha256:
        raise ValueError("validation source digest changed")
    native = inventory(prediction)
    if model == "modeller" or model == "protenix":
        if set(native) != set(request.fixed_atoms) | set(request.generated_atoms):
            raise ValueError("complete native requested inventory differs")
        mapped = {a: value[2] for a, value in native.items()}
        mapped_names = {a: value[0] for a, value in native.items()}
        custody = {"identity_transport": "already mapped exact requested IDs"}
    else:
        manifest = json.loads((case / "rfd3/manifest.json").read_text())
        if model == "rfd3":
            metadata = prediction.with_suffix("").with_suffix(".json")
            mapping = json.loads(metadata.read_text())["diffused_index_map"]
            custody = {"metadata_sha256": digest(metadata)}
        elif model == "boltz":
            chains = {a.chain for a in native}
            if len(chains) != 1:
                raise ValueError("Boltz native chain inventory differs")
            native_chain = next(iter(chains))
            target = request.target_sequences[0]
            mapping = {f"{r.chain}{r.residue_number}": f"{native_chain}{i+1}"
                       for i, r in enumerate(target_axis(request)[target.chain])}
            custody = {"identity_transport": "unique native chain, full sequence/atom verified ordinal mapping", "native_chain": native_chain,
                "original_chain_preserved_in_output": target.chain}
        else:
            raise ValueError("unknown private validation backend")
        mapped, mapped_names, seen = {}, {}, set()
        for row in manifest["axis"]:
            atom = AtomIdentity(**row["atom"])
            match = re.fullmatch(r"(.)(-?\d+)([A-Za-z]?)", mapping.get(f"{atom.chain}{atom.residue_number}", ""))
            if match is None:
                raise ValueError("native residue mapping is missing/truncated")
            chain, number, code = match.groups()
            key = AtomIdentity(chain, number, code, atom.atom_name)
            if key in seen or key not in native or native[key][0] != row["residue_name"]:
                raise ValueError("native chemistry/axis differs from known sequence")
            seen.add(key)
            mapped[atom], mapped_names[atom] = native[key][2], native[key][0]
        # A pinned reader may represent the existing terminal oxygen. If present,
        # retain and audit it; never drop it to force a canonical-axis match.
        terminal = target_axis(request)[request.target_sequences[0].chain][-1]
        for atom in source:
            if atom.atom_name != "OXT" or (atom.chain, atom.residue_number, atom.insertion_code) != (terminal.chain, terminal.residue_number, terminal.insertion_code):
                continue
            match = re.fullmatch(r"(.)(-?\d+)([A-Za-z]?)", mapping.get(f"{atom.chain}{atom.residue_number}", ""))
            if match is not None:
                chain, number, code = match.groups()
                key = AtomIdentity(chain, number, code, "OXT")
                if key in native and key not in seen:
                    seen.add(key)
                    mapped[atom], mapped_names[atom] = native[key][2], native[key][0]
        if seen != set(native):
            raise ValueError("unexpected native atoms")
    required = set(request.fixed_atoms) | set(request.generated_atoms)
    missing = required-set(mapped)
    last = target_axis(request)[request.target_sequences[0].chain][-1]
    allowed = {a for a in source if a.atom_name == "OXT" and (a.chain, a.residue_number, a.insertion_code) ==
               (last.chain, last.residue_number, last.insertion_code)}
    if missing and (model not in {"boltz", "rfd3"} or not missing <= allowed):
        raise ValueError("a requested non-passthrough native atom is missing")
    expected_names = {}
    for target in request.target_sequences:
        for residue, symbol in zip(target_axis(request)[target.chain], target.sequence):
            expected_names[residue.chain, residue.residue_number, residue.insertion_code] = _CANONICAL_ONE_TO_THREE[symbol]
    if any(mapped_names[a] != expected_names[a.chain, a.residue_number, a.insertion_code]
           or not np.isfinite(xyz).all() for a, xyz in mapped.items()):
        raise ValueError("native chemistry or finite-coordinate check failed")
    fixed = tuple(sorted(set(source) & set(mapped)))
    transform = weighted_kabsch(np.array([mapped[a] for a in fixed]), np.array([source[a][2] for a in fixed]))
    aligned = {a: transform.apply(np.array([xyz]))[0] for a, xyz in mapped.items()}
    drift = np.array([np.linalg.norm(aligned[a]-source[a][2]) for a in fixed])
    axis = tuple(sorted(required))
    names = tuple(mapped_names[a] if a in mapped_names else source[a][0] for a in axis)
    # Missing OXT is explicitly borrowed from observations, never called native.
    native_with_passthrough = np.array([aligned[a] if a in aligned else source[a][2] for a in axis])
    restored = np.array([source[a][2] if a in source else aligned[a] for a in axis])
    source_text = (case / "source.pdb").read_text()
    output.mkdir(parents=True, exist_ok=False)
    for name, xyz in (("aligned-proposal-with-passthrough.pdb", native_with_passthrough), ("source-restored.pdb", restored)):
        (output / name).write_text(materialize_candidate_pdb(source_text, axis, xyz, names, request, elements=tuple(a.atom_name[0] for a in axis)))
    (output / "source.pdb").write_text(source_text)
    (output / "evaluation-mask.json").write_text(request.to_json())
    np.savez_compressed(output / "aligned-arrays.npz", coordinates=native_with_passthrough,
        native_represented=np.array([a in mapped for a in axis]), rotation=transform.rotation, translation=transform.translation)
    report = {"backend": model, "prediction_sha256": digest(prediction), "source_sha256": digest(case / "source.pdb"),
        "custody": custody, "fixed_atoms_represented_in_native": len(fixed), "fixed_atoms_passed_through_only": len(missing),
        "native_observed_rms_drift": float(np.sqrt(np.mean(drift**2))), "native_observed_max_drift": float(drift.max()),
        "native_inventory_checked_without_heldout_coordinates": True, "source_restoration_separate": True, "qualification": False}
    save(output / "transport.json", report)
    return report
