"""Private source-only RFD3 A/C/D conditioning contracts, never public admission.

Sampling-unfixed atoms and newly generated ownership are distinct: a backbone
information ablation withholds observed sidechains from the model, but they
remain immutable observations in the eventual materialized structure.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
from scripts.audit_mac_reconstruction_smoke import inventory
from scripts.mac_reconstruction_comparison import digest, save, target_axis

from dvbfixer.model.diffusion.contract import AtomIdentity, DiffusionRequest
from dvbfixer.model.diffusion.geometry import weighted_kabsch
from dvbfixer.model.diffusion.pdb_materialize import materialize_candidate_pdb
from dvbfixer.model.diffusion.protein_atoms import _CANONICAL_ONE_TO_THREE, CANONICAL_HEAVY_ATOMS

BACKBONE = frozenset({"N", "CA", "C", "O"})


def prepare(case: Path, output: Path, *, gaps: tuple[int, ...], context: str, information: str):
    case, output = case.resolve(), output.resolve()
    if context not in {"full", "flanks8", "spatial6"} or information not in {"all-heavy", "backbone"}:
        raise ValueError("undeclared context/atom-information policy")
    request = DiffusionRequest.from_json((case / "evaluation-mask.json").read_text())
    if len(request.target_sequences) != 1 or request.retained_explicit_links or request.heterogen_contexts:
        raise ValueError("first conditioning slice is one canonical protein chain, no links")
    if not gaps or len(set(gaps)) != len(gaps) or any(not 0 <= i < len(request.gaps) for i in gaps):
        raise ValueError("declare unique valid whole-gap indices")
    source_path = case / request.normalized_pdb.path
    if digest(source_path) != request.normalized_pdb.sha256:
        raise ValueError("source differs from the original mask")
    source = inventory(source_path)
    source_lines = {AtomIdentity(line[21], line[22:26].strip(), line[26].strip(), line[12:16].strip()): line
                    for line in source_path.read_text().splitlines(keepends=True) if line.startswith("ATOM  ")}
    if set(source) != set(request.fixed_atoms):
        raise ValueError("source inventory differs from declared fixed observations")
    target, placement = request.target_sequences[0], request.sequence_placements[0]
    residues = target_axis(request)[target.chain]
    ordinal = {(r.chain, r.residue_number, r.insertion_code): i for i, r in enumerate(residues)}
    selected_generated = {i for g in gaps for i in range(request.gaps[g].target_interval.start, request.gaps[g].target_interval.stop)}
    observed_indices = set(placement.observed_target_indices)
    flanks = {i for g in gaps for i in observed_indices if
              request.gaps[g].target_interval.start-8 <= i < request.gaps[g].target_interval.start or
              request.gaps[g].target_interval.stop <= i < request.gaps[g].target_interval.stop+8}
    included = set(observed_indices) if context == "full" else set(flanks)
    if context == "spatial6":
        anchors = np.array([value[2] for atom, value in source.items()
                            if ordinal[atom.chain, atom.residue_number, atom.insertion_code] in flanks])
        if len(anchors) < 3:
            raise ValueError("spatial context lacks an observed source frame")
        for atom, value in source.items():
            if np.min(np.linalg.norm(anchors-value[2], axis=1)) <= 6.0:
                included.add(ordinal[atom.chain, atom.residue_number, atom.insertion_code])
    included |= selected_generated
    # Any other missing region stays absent, never a fixed zero-valued obstacle.
    if not included <= observed_indices | selected_generated or len(included) > 256:
        raise ValueError("conditioning includes an unselected gap or exceeds its shape cap")
    ranges = []
    for index in sorted(included):
        if not ranges or ranges[-1][1] != index:
            ranges.append([index, index+1])
        else:
            ranges[-1][1] += 1
    lines, axis, masks = [], [], {}
    for index in sorted(included):
        residue = residues[index]
        if residue.insertion_code:
            raise ValueError("RFD3 insertion-code transport is not admitted here")
        name = _CANONICAL_ONE_TO_THREE[target.sequence[index]]
        generated = index in selected_generated
        masks[f"{residue.chain}{residue.residue_number}"] = "" if generated else "ALL" if information == "all-heavy" else "N,CA,C,O"
        for atom_name in CANONICAL_HEAVY_ATOMS[name]:
            atom = AtomIdentity(residue.chain, residue.residue_number, residue.insertion_code, atom_name)
            fixed = not generated and (information == "all-heavy" or atom_name in BACKBONE)
            if fixed and atom not in source:
                raise ValueError("missing fixed conditioning atom")
            xyz = source[atom][2] if fixed else np.zeros(3)
            field = f" {atom_name:<3}" if len(atom_name) < 4 else atom_name
            if fixed:
                line = source_lines[atom]
                lines.append(line[:6]+f"{len(lines)+1:5d}"+line[11:])
            else:
                lines.append(f"ATOM  {len(lines)+1:5d} {field} {name:3s} {atom.chain}{int(atom.residue_number):4d}    "
                             f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{1:6.2f}{0:6.2f}          {atom_name[0]:>2s}\n")
            axis.append({"atom": {"chain": atom.chain, "residue_number": atom.residue_number, "insertion_code": atom.insertion_code, "atom_name": atom.atom_name},
                "residue_name": name, "target_generated": generated, "sampler_fixed": fixed, "original_observation": atom in source})
    output.mkdir(parents=True, exist_ok=False)
    (output / "conditioning.pdb").write_text("REMARK 999 PRIVATE SOURCE-ONLY CONDITIONING ABLATION\n"+"".join(lines)+"TER\nEND\n")
    contigs = []
    for start, stop in ranges:
        if [int(residues[i].residue_number) for i in range(start, stop)] != list(range(int(residues[start].residue_number), int(residues[stop-1].residue_number)+1)):
            raise ValueError("noncontiguous source numbering cannot be represented by a numeric contig")
        contigs.append(f"{target.chain}{residues[start].residue_number}-{residues[stop-1].residue_number}")
    save(output / "rfd3-input.json", {"conditioning_control": {"input": str(output / "conditioning.pdb"),
        "contig": ",/0,".join(contigs), "select_fixed_atoms": masks, "select_unfixed_sequence": False, "plddt_enhanced": True}})
    receipt = {"contract": "private-source-only-conditioning-v1", "case": str(case), "source_sha256": digest(source_path),
        "request_sha256": digest(case / "evaluation-mask.json"), "gaps": gaps, "context": context, "information": information,
        "source_target_spans": ranges, "axis": axis, "sampler_fixed_atoms": sum(row["sampler_fixed"] for row in axis),
        "sampler_unfixed_atoms": sum(not row["sampler_fixed"] for row in axis), "target_generated_atoms": sum(row["target_generated"] for row in axis),
        "full_source_fixed_atoms": len(source), "observed_sidechains_are_not_new_generated_ownership": True,
        "reference_coordinates_read": False, "omitted_source_not_assumed_irrelevant": True, "public_capability_enabled": False}
    save(output / "manifest.json", receipt)
    return receipt


def native_coordinates(inputs: Path, run: Path):
    manifest = json.loads((inputs / "manifest.json").read_text())
    predictions = list((run / "predictions").glob("*.cif.gz"))
    if len(predictions) != 1:
        raise ValueError("missing/ambiguous native output")
    prediction = predictions[0]
    metadata = prediction.with_suffix("").with_suffix(".json")
    mapping = json.loads(metadata.read_text())["diffused_index_map"]
    native = inventory(prediction)
    source = inventory(Path(manifest["case"]) / "source.pdb")
    if digest(Path(manifest["case"]) / "source.pdb") != manifest["source_sha256"]:
        raise ValueError("source custody changed")
    original, seen = {}, set()
    fixed = []
    target_generated = set()
    for row in manifest["axis"]:
        atom = AtomIdentity(**row["atom"])
        label = f"{atom.chain}{atom.residue_number}"
        match = re.fullmatch(r"(.)(-?\d+)([A-Za-z]?)", mapping.get(label, ""))
        if match is None:
            raise ValueError("invalid/truncated native identity transport")
        chain, number, code = match.groups()
        transported = AtomIdentity(chain, number, code, atom.atom_name)
        if transported in seen or transported not in native or native[transported][0] != row["residue_name"]:
            raise ValueError("native sequence/atom inventory differs")
        seen.add(transported)
        original[atom] = native[transported][2]
        if row["sampler_fixed"]:
            fixed.append(atom)
        if row["target_generated"]:
            target_generated.add(atom)
    if seen != set(native) or len(fixed) < 3:
        raise ValueError("native atom inventory or fixed frame is incomplete")
    transform = weighted_kabsch(np.array([original[a] for a in fixed]), np.array([source[a][2] for a in fixed]))
    aligned = {a: transform.apply(np.array([xyz]))[0] for a, xyz in original.items()}
    represented_observed = sorted(set(aligned) & set(source))
    drift = np.array([np.linalg.norm(aligned[a]-source[a][2]) for a in represented_observed])
    return {a: aligned[a] for a in target_generated}, {"native_sha256": digest(prediction), "metadata_sha256": digest(metadata),
        "manifest_sha256": digest(inputs / "manifest.json"), "source_sha256": manifest["source_sha256"],
        "atomic_gap_indices": manifest["gaps"], "fit_uses_only_sampler_fixed_observations": True,
        "represented_observed_rms_drift": float(np.sqrt(np.mean(drift**2))), "represented_observed_max_drift": float(drift.max()),
        "sampler_fixed_atoms": len(fixed), "new_generated_ownership": len(target_generated),
        "sampling_unfixed_not_equal_to_new_generated_ownership": manifest["sampler_unfixed_atoms"] != len(target_generated)}


def materialize(case: Path, groups):
    request = DiffusionRequest.from_json((case / "evaluation-mask.json").read_text())
    source = inventory(case / "source.pdb")
    generated = {}
    for coordinates in groups:
        if generated.keys() & coordinates.keys():
            raise ValueError("overlapping selected donor ownership")
        generated.update(coordinates)
    if generated.keys() != set(request.generated_atoms):
        raise ValueError("selected groups do not cover all and only target-generated atoms")
    residues = target_axis(request)[request.target_sequences[0].chain]
    names = {r: _CANONICAL_ONE_TO_THREE[s] for r, s in zip(residues, request.target_sequences[0].sequence)}
    atom_names = {(r.chain, r.residue_number, r.insertion_code): name for r, name in names.items()}
    axis = tuple(sorted(set(source) | set(generated)))
    xyz = np.array([source[a][2] if a in source else generated[a] for a in axis])
    return materialize_candidate_pdb((case / "source.pdb").read_text(), axis, xyz,
        tuple(atom_names[a.chain, a.residue_number, a.insertion_code] for a in axis), request,
        elements=tuple(a.atom_name[0] for a in axis))
