"""Independent IDPForge transport/native custody readback; optional raw extraction."""
from __future__ import annotations

import argparse
import ast
import json
import re
from dataclasses import asdict
from pathlib import Path

import numpy as np
from scripts.audit_anchor_calibration_grid import diagnose
from scripts.audit_retained_full_context_controls import atoms
from scripts.hybrid_gap_reconstruction import diagnostic_counts
from scripts.mac_reconstruction_comparison import digest, save, target_axis
from scripts.retained_evidence import freeze_code, verify_bindings

from dvbfixer.model.diffusion.contract import DiffusionRequest
from dvbfixer.model.diffusion.pdb_materialize import materialize_candidate_pdb
from dvbfixer.model.diffusion.protein_atoms import _CANONICAL_ONE_TO_THREE, CANONICAL_HEAVY_ATOMS
from dvbfixer.model.diffusion.staged_runtime import assessment_adapter


def atom14_names(path):
    tree = ast.parse(path.read_text())
    assignment = next(n for n in tree.body if isinstance(n, ast.Assign)
                      and any(isinstance(t, ast.Name) and t.id == "restype_name_to_atom14_names" for t in n.targets))
    names = ast.literal_eval(assignment.value)
    for name, canonical in CANONICAL_HEAVY_ATOMS.items():
        if len(names[name]) != 14 or set(filter(None, names[name])) != set(canonical):
            raise ValueError("native atom14 names do not match canonical inventory")
    return names


def register_ca(xyz, moving_ca, source_ca):
    """One proper rigid transform for every atom, not local geometry repair."""
    if moving_ca.shape != source_ca.shape or moving_ca.ndim != 2 or moving_ca.shape[1] != 3:
        raise ValueError("invalid registration correspondence")
    if len(moving_ca) < 3 or not all(np.isfinite(a).all() for a in (xyz, moving_ca, source_ca)):
        raise ValueError("invalid registration coordinates")
    ac, bc = moving_ca.mean(axis=0), source_ca.mean(axis=0)
    u, singular, vt = np.linalg.svd((moving_ca-ac).T @ (source_ca-bc))
    if singular[1] < 1e-8:
        raise ValueError("degenerate registration")
    correction = np.eye(3)
    correction[-1, -1] = np.linalg.det(u @ vt)
    rotation = u @ correction @ vt
    translation = bc - ac @ rotation
    result = xyz @ rotation + translation
    delta = moving_ca @ rotation + translation - source_ca
    return result, {"rotation_row_vector": rotation.tolist(), "translation_A": translation.tolist(),
                    "observed_CA_rms_A": float(np.sqrt(np.mean(np.sum(delta**2, axis=-1)))),
                    "observed_CA_max_A": float(np.linalg.norm(delta, axis=-1).max()),
                    "reflection": False, "per_residue_fit": False}


def audit(root):
    plan_path = root / "plan.json"
    plan = json.loads(plan_path.read_text())
    verify_bindings(plan["source_code_and_input_sha256"], (plan_path,))
    execution = json.loads((root / "execution.json").read_text())
    if execution["returncode"] or execution["guard_stop_reason"] is not None:
        raise ValueError("native worker failed or was stopped")
    request = DiffusionRequest.from_json((root / "request.json").read_text())
    if digest(root / "source.pdb") != request.normalized_pdb.sha256:
        raise ValueError("source hash mismatch")
    source = atoms(root / "source.pdb")
    axis = target_axis(request)
    order = [r for t in request.target_sequences for r in axis[t.chain]]
    generated_residues = {(r.chain, r.residue_number, r.insertion_code) for g in request.gaps for r in g.generated_residues}
    keys = [(r.chain, r.residue_number, r.insertion_code) for r in order]
    fixed = np.array([k not in generated_residues for k in keys])
    seq = "".join(t.sequence for t in request.target_sequences)
    names = atom14_names(Path(plan["runtime_src"]) / "openfold/np/residue_constants.py")
    template = dict(np.load(root / "template.npz", allow_pickle=False))
    initial = dict(np.load(root / "native-initial.npz", allow_pickle=False))
    origin = source[(*keys[int(np.flatnonzero(fixed)[0])], "CA")]["xyz"]
    if str(template["seq"]) != seq or not np.array_equal(template["mask"], fixed):
        raise ValueError("native template sequence or residue mask changed")
    expected = np.zeros((len(order), 14, 3))
    observed_mask = np.zeros((len(order), 14), dtype=bool)
    native_names = []
    for i, (key, aa) in enumerate(zip(keys, seq, strict=True)):
        row = names[_CANONICAL_ONE_TO_THREE[aa]]
        native_names.append(row)
        if fixed[i]:
            for j, name in enumerate(row):
                if name:
                    expected[i, j] = np.asarray(source[(*key, name)]["xyz"]) - origin
                    observed_mask[i, j] = True
    if (not np.array_equal(template["coord"], expected)
            or np.any(template["torsion"][~fixed]) or not np.isfinite(template["torsion"]).all()):
        raise ValueError("template invented source coordinates or unknown torsions")
    used = (expected[None] - template["coord_offset"][0][None, None, None]).astype(np.float32)
    if (not np.array_equal(initial["template_coord"], used)
            or not np.array_equal(initial["template_mask"], fixed[None])
            or not np.array_equal(initial["template_torsion"], template["torsion"][None].astype(np.float32))):
        raise ValueError("native conditioning differs from declared prepared template")
    trace = [json.loads(line) for line in (root / "consumer-trace.jsonl").read_text().splitlines()]
    count = 1 if plan["mode"] == "preflight" else 40
    paths = sorted((root / "consumer-states").glob("*.npz"))
    if len(trace) != count or len(paths) != count:
        raise ValueError("incomplete native model consumer trace")
    aa = np.array([["ARNDCQEGHILKMFPSTWYV".index(a) for a in seq]])
    ss = None
    for i, (row, path) in enumerate(zip(trace, paths, strict=True)):
        state = dict(np.load(path, allow_pickle=False))
        if ss is None:
            ss = state["ss"]
        if (row["call"] != i+1 or row["device"] != "cpu" or row["target_residues"] != len(order)
                or not np.array_equal(state["sequence"], aa)
                or not np.array_equal(state["resi"], np.arange(len(order))[None])
                or not np.array_equal(state["ss"], ss)
                or not np.all(state["mask"] == 1)
                or not np.all(state["t"][0, fixed] == 0)
                or not np.all(state["t"][0, ~fixed] == 199-5*i)
                or not np.array_equal(state["xyz"][0, fixed], used[0, fixed, :5])
                or not np.array_equal(state["tor"][0, fixed], initial["template_torsion"][0, fixed])):
            raise ValueError("saved native consumer sequence/index/mask/conditioning changed")
    if plan.get("template_ss_policy", "stock") == "radian-consistent":
        sec = str(template["sec"])
        if not any(c in sec for c in ("A", "B")):
            raise ValueError("this SS control requires deterministic native SS encoding")
        sec = re.sub(r'[AH]{6,}', lambda m: "H" * len(m.group()), sec)
        expected_ss = np.array([["HEPABCL-".index(c) for c in sec]])
        if not np.array_equal(ss, expected_ss):
            raise ValueError("corrected source secondary structure did not reach native model")
    report = {"plan_sha256": digest(plan_path), "auditor_sha256": digest(Path(__file__)),
              "consumer_passed": True, "native_model_calls": count, "target_residues": len(order),
              "fixed_residues": int(fixed.sum()), "generated_residues": int((~fixed).sum()),
              "source_atoms": len(source), "sampling_completed": plan["mode"] == "sample",
              "full_context_cropped": False, "strict_atomic_preservation_claimed": False,
              "public_qualification": False}
    if plan["mode"] == "sample":
        previous = dict(np.load(Path(plan["preflight"]) / "native-initial.npz", allow_pickle=False))
        if set(previous) != set(initial) or not all(np.array_equal(previous[k], initial[k]) for k in initial):
            raise ValueError("sample/preflight initialization differs")
        endpoint = dict(np.load(root / "native-endpoint.npz", allow_pickle=False))
        xyz = endpoint["positions"]
        expected_exists = np.array([[[bool(a) for a in row] for row in native_names]], dtype=int)
        if (xyz.shape != (1, len(order), 14, 3) or not np.isfinite(xyz).all()
                or not np.array_equal(endpoint["aatype"], aa)
                or not np.array_equal(endpoint["residue_index"], np.arange(len(order))[None])
                or not np.array_equal(endpoint["atom14_atom_exists"], expected_exists)):
            raise ValueError("native endpoint sequence/atom axis changed")
        restored_frame = xyz[0] + np.asarray(origin) + template["coord_offset"][0]
        observed = expected + origin
        delta = np.linalg.norm((restored_frame-observed)[observed_mask], axis=-1)
        report["native_observed_heavy_motion_A"] = {"max": float(delta.max()),
            "rms": float(np.sqrt(np.mean(delta**2))), "mean": float(delta.mean()), "globally_fitted": False}
        report["identical_preflight_initialization"] = True
    return report


def extract(root, frame_policy="native"):
    report = audit(root)
    if not report["sampling_completed"]:
        raise ValueError("no completed raw sample")
    plan = json.loads((root / "plan.json").read_text())
    request = DiffusionRequest.from_json((root / "request.json").read_text())
    source = atoms(root / "source.pdb")
    order = target_axis(request)[request.target_sequences[0].chain]
    seq = request.target_sequences[0].sequence
    names = atom14_names(Path(plan["runtime_src"]) / "openfold/np/residue_constants.py")
    template = dict(np.load(root / "template.npz", allow_pickle=False))
    anchor = request.gaps[0].right_anchor
    origin = source[(anchor.chain, anchor.residue_number, anchor.insertion_code, "CA")]["xyz"]
    xyz = np.load(root / "native-endpoint.npz", allow_pickle=False)["positions"][0] + origin + template["coord_offset"][0]
    fixed = np.array([(r.chain, r.residue_number, r.insertion_code, "CA") in source for r in order])
    registration = None
    if frame_policy == "ca-rigid-fit":
        source_ca = np.array([source[r.chain, r.residue_number, r.insertion_code, "CA"]["xyz"]
                              for r, held in zip(order, fixed, strict=True) if held])
        xyz, registration = register_ca(xyz, xyz[fixed, 1], source_ca)
    elif frame_policy != "native":
        raise ValueError("unknown extraction frame policy")
    folder = root / ("raw-extraction" if frame_policy == "native" else "rigid-registration-extraction")
    folder.mkdir(exist_ok=False)
    save(folder / "plan.json", {"protocol": "private-idpforge-raw-extraction-v1", "frame_policy": frame_policy,
         "source_code_and_input_sha256": {str(p): digest(p) for p in (
             root / "plan.json", root / "native-endpoint.npz", root / "template.npz", root / "source.pdb", root / "request.json",
             Path(__file__).resolve(), Path(plan["runtime_src"]) / "openfold/np/residue_constants.py")},
         "registration": registration, "new_diffusion_calls": 0, "source_restoration_separate": True,
         "local_geometry_repair": False, "selection_or_resampling": False})
    freeze_code(folder / "plan.json")
    coordinates, lines = {}, []
    residue_name_map = {}
    serial = 0
    for i, (r, letter) in enumerate(zip(order, seq, strict=True)):
        resname = _CANONICAL_ONE_TO_THREE[letter]
        residue_name_map[r.chain, r.residue_number, r.insertion_code] = resname
        for j, name in enumerate(names[resname]):
            if not name:
                continue
            coordinates[(r.chain, r.residue_number, r.insertion_code, name)] = xyz[i, j]
            serial += 1
            p = xyz[i, j]
            lines.append(f'ATOM  {serial:5d} {name:^4s} {resname} {r.chain}{int(r.residue_number):4d}{r.insertion_code or " "}   '
                         f'{p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}  1.00  0.00           {name[0]}  \n')
    raw = folder / "native-mapped.pdb"
    raw.write_text("".join(lines) + "TER\nEND\n")
    generated_xyz = [coordinates[(a.chain, a.residue_number, a.insertion_code, a.atom_name)] for a in request.generated_atoms]
    residue_names = tuple(residue_name_map[a.chain, a.residue_number, a.insertion_code] for a in request.generated_atoms)
    text = materialize_candidate_pdb((root / "source.pdb").read_text(), request.generated_atoms,
                                     generated_xyz, residue_names, request)
    candidate = folder / "source-restored-candidate.pdb"
    candidate.write_text(text)
    actual = atoms(candidate)
    generated = {(a.chain, a.residue_number, a.insertion_code, a.atom_name) for a in request.generated_atoms}
    if set(actual) != set(source) | generated or any(actual[a]["line"] != row["line"] for a, row in source.items()):
        raise ValueError("explicit source restoration changed source records")
    if any(actual[a]["xyz"] != tuple(float(f"{v:.3f}") for v in coordinates[a]) for a in generated):
        raise ValueError("materialization changed the generated coordinates")
    assessment = assessment_adapter(request, (root / "source.pdb").read_text())(text)
    save(folder / "assessment.json", asdict(assessment))
    raw_findings = diagnose(raw, folder / "native-diagnose.json")
    restored_findings = diagnose(candidate, folder / "source-restored-diagnose.json")
    summary = {"plan_sha256": digest(root / "plan.json"), "native_endpoint_sha256": digest(root / "native-endpoint.npz"),
               "native_observed_heavy_motion_A": report["native_observed_heavy_motion_A"],
               "native_diagnose": diagnostic_counts(raw_findings),
               "source_restored_diagnose": diagnostic_counts(restored_findings),
               "source_restored_records_exact": True, "source_restored_geometry_passed": assessment.geometry_passed,
               "new_diffusion_calls": 0, "relaxation": False, "global_fit": frame_policy == "ca-rigid-fit",
               "registration": registration,
               "public_output_published": False, "candidate_sha256": digest(candidate),
               "note": "raw mapped neural atom14 versus separate source restoration, not native preservation"}
    save(folder / "summary.json", summary)
    verify_bindings(json.loads((folder / "plan.json").read_text())["source_code_and_input_sha256"])
    return summary


def audit_extraction(folder):
    """Read-only, independently reconstruct both serialized coordinate arms."""
    root = folder.parent
    audit(root)
    plan_path = folder / "plan.json"
    declaration = json.loads(plan_path.read_text())
    verify_bindings(declaration["source_code_and_input_sha256"], (plan_path,))
    request = DiffusionRequest.from_json((root / "request.json").read_text())
    source = atoms(root / "source.pdb")
    target = request.target_sequences[0]
    order = target_axis(request)[target.chain]
    names = atom14_names(Path(json.loads((root / "plan.json").read_text())["runtime_src"]) / "openfold/np/residue_constants.py")
    template = np.load(root / "template.npz", allow_pickle=False)
    anchor = request.gaps[0].right_anchor
    origin = source[anchor.chain, anchor.residue_number, anchor.insertion_code, "CA"]["xyz"]
    xyz = np.load(root / "native-endpoint.npz", allow_pickle=False)["positions"][0] + origin + template["coord_offset"][0]
    if declaration["frame_policy"] == "ca-rigid-fit":
        transform = declaration["registration"]
        rotation, translation = np.array(transform["rotation_row_vector"]), np.array(transform["translation_A"])
        if not np.allclose(rotation.T @ rotation, np.eye(3), atol=1e-10) or not np.isclose(np.linalg.det(rotation), 1., atol=1e-10):
            raise ValueError("extraction is not one proper rigid transform")
        xyz = xyz @ rotation + translation
    elif declaration["frame_policy"] != "native":
        raise ValueError("unrecognized extraction policy")
    coordinates = {(r.chain, r.residue_number, r.insertion_code, name): tuple(float(f"{v:.3f}") for v in xyz[i, j])
                   for i, (r, aa) in enumerate(zip(order, target.sequence, strict=True))
                   for j, name in enumerate(names[_CANONICAL_ONE_TO_THREE[aa]]) if name}
    mapped, restored = atoms(folder / "native-mapped.pdb"), atoms(folder / "source-restored-candidate.pdb")
    generated = {(a.chain, a.residue_number, a.insertion_code, a.atom_name) for a in request.generated_atoms}
    if (set(mapped) != set(coordinates) or set(restored) != set(source) | generated
            or any(mapped[a]["xyz"] != xyz for a, xyz in coordinates.items())
            or any(restored[a]["xyz"] != coordinates[a] for a in generated)
            or any(restored[a]["line"] != row["line"] for a, row in source.items())):
        raise ValueError("extraction changed neural coordinates, ownership or source records")
    assessment = assessment_adapter(request, (root / "source.pdb").read_text())((folder / "source-restored-candidate.pdb").read_text())
    summary = json.loads((folder / "summary.json").read_text())
    if assessment.geometry_passed != summary["source_restored_geometry_passed"]:
        raise ValueError("source-restored geometry decision does not reproduce")
    return {"plan_sha256": digest(plan_path), "native_endpoint_sha256": digest(root / "native-endpoint.npz"),
            "frame_policy": declaration["frame_policy"], "source_records_exact": True,
            "generated_coordinates_match_raw_endpoint_and_declared_transform": True,
            "source_restored_geometry_passed": assessment.geometry_passed,
            "generated_atoms": len(generated), "source_atoms": len(source), "new_diffusion_calls": 0}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("--extract", action="store_true")
    parser.add_argument("--frame-policy", choices=("native", "ca-rigid-fit"), default="native")
    parser.add_argument("--audit-extraction", action="store_true")
    args = parser.parse_args()
    if args.audit_extraction:
        result = audit_extraction(args.root.resolve())
        save(args.root / "independent-extraction-readback.json", result)
    elif args.extract:
        result = extract(args.root.resolve(), args.frame_policy)
    else:
        result = audit(args.root.resolve())
        save(args.root / "independent-readback.json", result)
    print(json.dumps(result))
