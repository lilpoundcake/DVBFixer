#!/usr/bin/env python3
"""Weightless RFD3 full-target transport check for the retained1KTZ construct."""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from dataclasses import asdict
from pathlib import Path

import numpy as np
from scripts.mac_reconstruction_comparison import digest, save, target_axis

from dvbfixer.model.diffusion.contract import AtomIdentity, DiffusionRequest
from dvbfixer.model.diffusion.gap_portfolio import _coordinates
from dvbfixer.model.diffusion.protein_atoms import _CANONICAL_ONE_TO_THREE, CANONICAL_HEAVY_ATOMS
from dvbfixer.renumber import update_conect

PROJECT = Path(__file__).resolve().parents[1]
ROOT = PROJECT / "tmp/mac-real-gap-continuation-2026-11-12"
RUNTIME = PROJECT / ".artifacts/mac-reconstruction-2026-11-11"
ENGINE = RUNTIME / "engines/rfd3"
NATIVE_PYTHON = RUNTIME / "envs/rfd3/bin/python"


def conditioning(request, source_text):
    """All chains/gaps, exact fixed records, zero placeholders, no source rewrite."""
    axis = target_axis(request)
    points, records = _coordinates(source_text)
    fixed, restored = set(request.fixed_atoms), set(request.generated_atoms)
    if set(points) != fixed or fixed & restored:
        raise ValueError("Source differs from complete fixed-observation inventory")
    lines, rows, masks, contigs, passthrough, serial_map = [], [], {}, [], set(), {}
    for target in request.target_sequences:
        residues = axis[target.chain]
        numbers = [int(r.residue_number) for r in residues]
        if any(r.insertion_code for r in residues) or numbers != list(range(numbers[0], numbers[-1] + 1)):
            raise ValueError("This full-target preflight requires contiguous numeric source IDs without insertion codes")
        if any(not -999 <= number <= 9999 for number in numbers):
            raise ValueError("Conditioning residue number exceeds PDB field capacity")
        contigs.append(f"{target.chain}{numbers[0]}-{numbers[-1]}")
        for index, residue in enumerate(residues):
            name = _CANONICAL_ONE_TO_THREE[target.sequence[index]]
            expected = {AtomIdentity(residue.chain, residue.residue_number, residue.insertion_code, atom_name)
                        for atom_name in CANONICAL_HEAVY_ATOMS[name]}
            missing = expected & restored
            observed = expected & fixed
            if expected != missing | observed or (missing and observed):
                raise ValueError("Whole-residue canonical ownership required; no partial-observation reinterpretation")
            masks[f"{residue.chain}{residue.residue_number}"] = "" if missing else "ALL"
            for atom_name in CANONICAL_HEAVY_ATOMS[name]:
                atom = AtomIdentity(residue.chain, residue.residue_number, residue.insertion_code, atom_name)
                if atom in fixed:
                    line = records[atom]
                    if line[17:20].strip() != name or line[16].strip() or float(line[54:60]) != 1:
                        raise ValueError("Observed residue/altloc/occupancy differs from the request")
                    old_serial = int(line[6:11])
                    if old_serial in serial_map:
                        raise ValueError("Ambiguous repeated source atom serial")
                    serial_map[old_serial] = len(rows) + 1
                    line = line[:6] + f"{len(rows) + 1:5d}" + line[11:]
                    xyz = tuple(float(v) for v in points[atom])
                else:
                    field = f" {atom_name:<3}" if len(atom_name) < 4 else atom_name
                    line = (f"ATOM  {len(rows) + 1:5d} {field} {name:3s} {atom.chain}{int(atom.residue_number):4d}    "
                            f"{0:8.3f}{0:8.3f}{0:8.3f}{1:6.2f}{0:6.2f}          {atom_name[0]:>2s}\n")
                    xyz = (0., 0., 0.)
                lines.append(line if line.endswith("\n") else line + "\n")
                rows.append({"atom": asdict(atom), "residue_name": name, "fixed": atom in fixed, "coordinates": xyz})
        lines.append("TER\n")
        terminal = residues[-1]
        passthrough.update(atom for atom in fixed if atom.atom_name == "OXT" and
                           (atom.chain, atom.residue_number, atom.insertion_code) ==
                           (terminal.chain, terminal.residue_number, terminal.insertion_code))
    represented = {AtomIdentity(**row["atom"]) for row in rows}
    if fixed - represented != passthrough or restored - represented or represented != (fixed - passthrough) | restored:
        raise ValueError("Only observed true-terminal OXT may be preserved outside the native canonical axis")
    annotations = []
    for line in source_text.splitlines(keepends=True):
        if line.startswith(("SSBOND", "LINK")):
            annotations.append(line if line.endswith("\n") else line + "\n")
        elif line.startswith("CONECT"):
            mapped = update_conect(line, serial_map)
            if mapped is None:
                raise ValueError("Cannot carry a dangling or pass-through CONECT into conditioning")
            annotations.append(mapped)
    if len(rows) >= 99999:
        raise ValueError("Conditioning serial capacity exceeded")
    text = "REMARK 999 PRIVATE COMPLETE DEPOSITED TARGET; ZERO MISSING-COORDINATE PLACEHOLDERS\n"
    text += "".join(annotations + lines) + "END\n"
    return text, {"axis": rows, "contig": ",/0,".join(contigs), "select_fixed_atoms": masks,
        "preserved_only_atoms": [asdict(atom) for atom in sorted(passthrough)],
        "source_link_count": len(request.retained_explicit_links), "all_chains_and_gaps_included": True,
        "OXT_passthrough_is_not_native_conditioning": True}


def check_array(array, manifest, request):
    """Exact identity/masks/sequence; common translation only; full peptide/link graph."""
    required = ("src_component", "atom_name", "res_name", "coord", "chain_id", "res_id", "ins_code",
                "is_motif_atom_with_fixed_coord", "is_motif_atom_with_fixed_seq", "bonds")
    if any(not hasattr(array, name) for name in required):
        raise ValueError("Native stage lost required transport annotations")
    expected = {(row["atom"]["chain"] + row["atom"]["residue_number"], row["atom"]["atom_name"]): row
                for row in manifest["axis"]}
    actual = {}
    for index in range(len(array)):
        key = str(array.src_component[index]), str(array.atom_name[index])
        if key in actual:
            raise ValueError("Native builder duplicated a source atom")
        actual[key] = index
    if actual.keys() != expected.keys():
        raise ValueError("Native atom inventory differs from complete canonical target")
    fixed = np.asarray(array.is_motif_atom_with_fixed_coord, dtype=bool)
    sequence_fixed = np.asarray(array.is_motif_atom_with_fixed_seq, dtype=bool)
    coords = np.asarray(array.coord)
    if (fixed.shape != (len(expected),) or sequence_fixed.shape != (len(expected),) or
            coords.shape != (len(expected), 3) or not np.isfinite(coords).all() or not sequence_fixed.all()):
        raise ValueError("Native coordinates/known-sequence mask invalid")
    deltas = []
    identity_index = {}
    residue_groups, chain_groups = {}, {}
    for key, row in expected.items():
        i = actual[key]
        atom = AtomIdentity(**row["atom"])
        identity_index[atom] = i
        native_residue = str(array.chain_id[i]), int(array.res_id[i]), str(array.ins_code[i])
        residue_groups.setdefault(key[0], set()).add(native_residue)
        chain_groups.setdefault(atom.chain, set()).add(native_residue[0])
        if str(array.res_name[i]) != row["residue_name"] or bool(fixed[i]) != row["fixed"]:
            raise ValueError("Native known sequence or atom-specific observation mask changed")
        if row["fixed"]:
            deltas.append(coords[i] - np.asarray(row["coordinates"]))
        elif not np.array_equal(coords[i], np.zeros(3)):
            raise ValueError("Missing-coordinate placeholder retained nonzero coordinates")
    if (any(len(group) != 1 for group in residue_groups.values()) or
            len({next(iter(group)) for group in residue_groups.values()}) != len(residue_groups) or
            any(len(group) != 1 for group in chain_groups.values()) or
            len({next(iter(group)) for group in chain_groups.values()}) != len(chain_groups)):
        raise ValueError("Native residue/chain groups collapsed or split the original full-target axes")
    if len(deltas) < 3:
        raise ValueError("Native transport lacks an observed coordinate frame")
    translation = np.mean(deltas, axis=0)
    residual = float(np.max(np.linalg.norm(np.asarray(deltas) - translation, axis=1)))
    if residual > 1e-4:
        raise ValueError("Native builder changed relative observed coordinates beyond float32 representation tolerance")
    if array.bonds is None:
        raise ValueError("Native builder has no covalent graph")
    bonds = {frozenset((int(a), int(b))) for a, b, *_ in array.bonds.as_array()}
    peptide, junctions = set(), set()
    axis = target_axis(request)
    restored = set(request.generated_atoms)
    for residues in axis.values():
        for left, right in zip(residues, residues[1:]):
            a = AtomIdentity(left.chain, left.residue_number, left.insertion_code, "C")
            b = AtomIdentity(right.chain, right.residue_number, right.insertion_code, "N")
            edge = frozenset((identity_index[a], identity_index[b]))
            peptide.add(edge)
            if (a in restored) != (b in restored):
                junctions.add(edge)
    links = {frozenset((identity_index[link.atom1], identity_index[link.atom2])) for link in request.retained_explicit_links}
    missing_peptide, missing_links = peptide - bonds, links - bonds
    reverse = {index: atom for atom, index in identity_index.items()}
    cross_residue = {edge for edge in bonds if len({(reverse[i].chain, reverse[i].residue_number, reverse[i].insertion_code) for i in edge}) > 1}
    extras = cross_residue - peptide - links
    def describe(edges):
        return [[asdict(reverse[index]) for index in sorted(edge)] for edge in sorted(edges, key=lambda e: tuple(sorted(e)))]

    return {"atoms": len(array), "fixed_atoms": int(fixed.sum()), "unfixed_atoms": int((~fixed).sum()),
        "source_components": len({key[0] for key in actual}), "native_chains": len(set(map(str, array.chain_id))),
        "all_sequence_fixed": True, "zero_missing_coordinates": True, "translation_angstrom": translation.tolist(),
        "max_observed_translation_residual_angstrom": residual, "float32_coordinate_tolerance_angstrom": 1e-4,
        "expected_peptide_bonds": len(peptide), "missing_peptide_bonds": len(missing_peptide),
        "expected_observed_missing_junctions": len(junctions), "missing_junction_bonds": len(junctions - bonds),
        "expected_retained_links": len(links), "missing_retained_links": len(missing_links),
        "unexpected_inter_residue_bonds": len(extras),
        "missing_retained_link_pairs": describe(missing_links), "missing_junction_pairs": describe(junctions - bonds),
        "unexpected_inter_residue_pairs": describe(extras),
        "graph_transport_passed": not missing_peptide and not missing_links and not extras}


def declare(output):
    reviews_path = ROOT / "construct-target-reviews-attempt01/construct-target-reviews.json"
    reviews = json.loads(reviews_path.read_text())
    case = next(row for row in reviews["rows"] if row["entry_id"] == "1KTZ")
    for name, expected in reviews["bindings"].items():
        if digest(Path(name)) != expected:
            raise ValueError("Reviewed source evidence changed")
    request_path = ROOT / "source-preparation-attempt01/execution/1KTZ-seed-0/request-source-prepared.json"
    request = DiffusionRequest.from_json(request_path.read_text())
    source_path = request_path.parent / "protein-projection.pdb"
    if digest(source_path) != request.normalized_pdb.sha256 or case["target_residue_instances"] != 468:
        raise ValueError("Bounded468-residue source/request changed")
    projection_path = request_path.parent / "projection.json"
    full_context = request_path.parent / "full-assembly-context.pdb"
    if digest(full_context) != json.loads(projection_path.read_text())["full_context_sha256"]:
        raise ValueError("Original full water context changed")
    text, manifest = conditioning(request, source_path.read_text())
    output.mkdir(parents=True, exist_ok=False)
    (output / "conditioning.pdb").write_text(text)
    save(output / "input.json", {"known_construct_completion": {"input": str(output / "conditioning.pdb"),
        "contig": manifest["contig"], "select_fixed_atoms": manifest["select_fixed_atoms"],
        "select_unfixed_sequence": False, "plddt_enhanced": True}})
    save(output / "manifest.json", {**manifest, "request_path": str(request_path), "request_sha256": digest(request_path),
        "source_path": str(source_path), "source_sha256": digest(source_path),
        "full_context_path": str(full_context), "full_context_sha256": digest(full_context),
        "water_context": case["full_context_omission_counts"], "water_endpoint_treatment_pending": True,
        "recorded_scope_refusal_unchanged": case["recorded_preparation_scope_only_not_new_backend_test"],
        "no_sampling_or_public_admission": True})
    consumed_paths = ("models/rfd3/src/rfd3/inference/input_parsing.py", "models/rfd3/src/rfd3/utils/inference.py")
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ENGINE, text=True).strip()
    (output / "engine-consumer-source.diff").write_bytes(subprocess.check_output(["git", "diff", "--", *consumed_paths], cwd=ENGINE))
    files = [Path(__file__), reviews_path, request_path, source_path, projection_path, full_context, NATIVE_PYTHON,
        output / "conditioning.pdb", output / "input.json", output / "manifest.json", output / "engine-consumer-source.diff",
        ENGINE / "models/rfd3/src/rfd3/inference/input_parsing.py", ENGINE / "models/rfd3/src/rfd3/utils/inference.py"]
    save(output / "plan.json", {"python": sys.executable, "python_sha256": digest(Path(sys.executable)),
        "source_code_and_input_sha256": {str(path): digest(path) for path in files}, "count": 1,
        "native_engine_revision": revision,
        "jobs": [{"case_id": "1KTZ-rfd3-native-data", "case": str(source_path), "seed": 0,
            "command": [str(NATIVE_PYTHON), "-m", "scripts.preflight_retained_rfd3_full_target", "worker", "--plan", str(output),
                "--output", str(output / "execution/1KTZ-rfd3-native-data-seed-0")]}],
        "scope": "Fixed468-residue complete target, data-only builder and pipeline input, not a size probe or checkpoint execution",
        "new_structural_model_calls": 0, "new_optimizer_calls": 0, "public_profile_or_frozen_request_changed": False})
    print(json.dumps({"target_residues": 468, "canonical_atoms": len(manifest["axis"]), "preserved_only_atoms": len(manifest["preserved_only_atoms"])}))


def worker(plan, output):
    if {p.name for p in output.iterdir()} != {"command.json", "output.log"}:
        raise ValueError("Native data worker requires a fresh executor reservation")
    os.environ["PYTORCH_ENABLE_MPS_FALLBACK"] = "0"
    os.environ.update(HF_HUB_OFFLINE="1", HF_DATASETS_OFFLINE="1", WANDB_MODE="disabled")
    import socket

    import torch

    def forbidden(*args, **kwargs):
        raise RuntimeError("Data-only preflight forbids model/checkpoint/network execution")

    torch.load = forbidden
    torch.serialization.load = forbidden
    torch.jit.load = forbidden
    socket.create_connection = forbidden
    socket.socket.connect = forbidden
    socket.socket.connect_ex = forbidden
    handle = torch.nn.modules.module.register_module_forward_pre_hook(forbidden)
    try:
        import safetensors
        import safetensors.torch
    except ImportError:
        safetensors_guard = False
    else:
        safetensors.safe_open = forbidden
        safetensors.torch.load_file = forbidden
        safetensors_guard = True
    try:
        from rfd3.inference import input_parsing

        if not Path(input_parsing.__file__).resolve().is_relative_to(ENGINE.resolve()):
            raise ValueError("Native input parser is not the owned pinned engine")
        manifest = json.loads((plan / "manifest.json").read_text())
        request = DiffusionRequest.from_json(Path(manifest["request_path"]).read_text())
        spec = json.loads((plan / "input.json").read_text())
        builder = input_parsing.DesignInputSpecification(**next(iter(spec.values())))
        built = builder.build()
        stages = {}

        def retain_stage(name, array):
            fields = {"coordinates": "coord", "source_component": "src_component", "atom_name": "atom_name",
                "residue_name": "res_name", "chain_id": "chain_id", "residue_id": "res_id", "insertion_code": "ins_code",
                "fixed": "is_motif_atom_with_fixed_coord", "sequence_fixed": "is_motif_atom_with_fixed_seq"}
            arrays = {key: np.asarray(getattr(array, attr)) for key, attr in fields.items() if hasattr(array, attr)}
            arrays["bonds"] = array.bonds.as_array() if array.bonds is not None else np.empty((0, 3), dtype=int)
            np.savez_compressed(output / f"{name}-input.npz", **arrays)
            try:
                stages[name] = check_array(array, manifest, request)
            except ValueError as exc:
                stages[name] = {"graph_transport_passed": False, "transport_error": str(exc)}
            save(output / f"{name}-readback.json", stages[name])

        retain_stage("builder", built)
        data = input_parsing.prepare_pipeline_input_from_atom_array(built.copy())
        retain_stage("pipeline", data["atom_array"])
        save(output / "result.json", {"case": "1KTZ", "stages": stages,
            "data_transport_passed": all(stage["graph_transport_passed"] for stage in stages.values()),
            "native_parser_path": input_parsing.__file__, "native_python": sys.executable,
            "torch_version": torch.__version__, "numpy_version": np.__version__,
            "checkpoint_loaded": False, "model_forward_executed": False, "network_allowed": False,
            "safetensors_guard_installed": safetensors_guard,
            "recorded_scope_refusal_unchanged": True, "sampling_or_scientific_admission": False,
            "preserved_only_OXT_atoms": manifest["preserved_only_atoms"],
            "OXT_passthrough_is_not_a_native_observation_constraint": True,
            "limits": "Input builder/pipeline data transport only, not final sampler tensors, fixed-coordinate dynamics, MPS execution or whole-context acceptance"})
    finally:
        handle.remove()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("declare", "run", "worker"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--plan", type=Path)
    args = parser.parse_args()
    if args.mode == "declare":
        declare(args.output.resolve())
    elif args.mode == "run":
        # The parent executor needs the full science stack; the isolated native
        # data worker deliberately must not import that dependency tree.
        from scripts.run_mac_modeller_comparison import execute

        execute(args.output.resolve())
    else:
        if args.plan is None:
            parser.error("worker requires --plan")
        worker(args.plan.resolve(), args.output.resolve())
