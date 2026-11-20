"""Read-only native-first audit; never hide native drift by source copying."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scripts.audit_retained_full_context_controls import atoms
from scripts.mac_reconstruction_comparison import digest, save, target_axis
from scripts.retained_evidence import verify_bindings

from dvbfixer.model.diffusion.contract import AtomIdentity, DiffusionRequest
from dvbfixer.model.diffusion.geometry import weighted_kabsch
from dvbfixer.model.diffusion.protein_atoms import (
    _CANONICAL_ONE_TO_THREE,
    ATOM37_ORDER,
    CANONICAL_HEAVY_ATOMS,
)


def native_coordinates(view, native, native_mask=None, *, expected_tokens=716):
    engine = DiffusionRequest.from_json(json.dumps(view["engine_view"]["request"]))
    axis = target_axis(engine)
    rows, native_names = {}, {}
    cursor = 0
    for target in engine.target_sequences:
        for residue, letter in zip(axis[target.chain], target.sequence, strict=True):
            rows[residue] = cursor
            native_names[residue] = CANONICAL_HEAVY_ATOMS[_CANONICAL_ONE_TO_THREE[letter]]
            cursor += 1
    if native.shape != (cursor, 37, 3) or cursor != expected_tokens or not np.isfinite(native).all():
        raise ValueError("native atom37 axis is not a finite complete tetramer")
    columns = {name: i for i, name in enumerate(ATOM37_ORDER)}
    row_keys = {(r.chain, r.residue_number, r.insertion_code): r for r in rows}
    if native_mask is not None:
        expected = np.zeros((cursor, 37), dtype=bool)
        for residue, row in rows.items():
            for name in native_names[residue]:
                expected[row, columns[name]] = True
        for atom in engine.fixed_atoms:
            residue = row_keys[(atom.chain, atom.residue_number, atom.insertion_code)]
            expected[rows[residue], columns[atom.atom_name]] = True
        if native_mask.dtype != np.bool_ or not np.array_equal(native_mask, expected):
            raise ValueError("native output mask is not canonical predictions plus exact fixed input atoms")
    decoded = {}
    for engine_atom, source_atom in view["engine_view"]["atom_to_source"]:
        e, s = AtomIdentity(**engine_atom), AtomIdentity(**source_atom)
        if e.atom_name not in columns:
            raise ValueError("an observed/generated atom is not representable; no omission allowed")
        residue = next(r for r in rows if (r.chain, r.residue_number, r.insertion_code)
                       == (e.chain, e.residue_number, e.insertion_code))
        # An array column alone is NOT evidence of a predicted atom. The pinned
        # canonical residue output mask omits OXT even though atom37 has a slot.
        # Return absent observations to the caller as missing, never as invented
        # coordinates or silently source-restored atoms.
        present = native_mask[rows[residue], columns[e.atom_name]] if native_mask is not None else e.atom_name in native_names[residue]
        if not present:
            continue
        if s in decoded:
            raise ValueError("duplicate native/source atom mapping")
        decoded[s] = native[rows[residue], columns[e.atom_name]].copy()
    return decoded


def run(root, output):
    output.mkdir(parents=True, exist_ok=False)
    (output / "executed-audit.py").write_bytes(Path(__file__).read_bytes())
    plan = json.loads((root / "plan.json").read_text())
    verify_bindings(plan["source_code_and_input_sha256"], (root / "plan.json",))
    execution_path = root / "guarded/execution.json"
    execution = json.loads(execution_path.read_text()) if execution_path.is_file() else None
    native_dir = root / "native-sampling"
    trace_path = native_dir / "trace.json"
    trace = json.loads(trace_path.read_text()) if trace_path.is_file() else {}
    report = {"plan_sha256": digest(root / "plan.json"), "guarded_execution": execution,
        "audit_code_sha256": digest(Path(__file__)),
        "audit_revision": "native-inpainting-mask-aware-v3" if plan.get("fixed_scaffold_inpainting") else "native-residue-mask-aware-v2",
        "trace": trace, "sampling_success": False, "native_fixed_passed": False,
        "sampling_execution_completed": bool(execution and not execution["returncode"]
            and not execution["guard_stop_reason"] and trace.get("ode_updates") == 500),
        "sampling_success_scope": "accepted sampling output; actual 500-update execution reported separately",
        "protein_mask_geometry_passed": False, "refinement_calls": 0, "source_coordinates_copied_before_gate": False,
        "public_acceptance": False}
    if not (native_dir / "native.npz").is_file():
        report.update(outcome="no_native_endpoint", reason=trace.get("error", execution.get("guard_stop_reason") if execution else "execution unavailable"))
        save(output / "receipt.json", report)
        return
    if digest(native_dir / "source.pdb") != digest(Path(plan["source"])):
        raise ValueError("actual sampled source differs from declared protein-only input")
    if digest(native_dir / "native.npz") != trace.get("native_sha256"):
        raise ValueError("native payload differs from its committed trace")
    if (trace.get("seed") != plan["seed"] or trace.get("tokens") != 716
            or trace.get("generated_gaps") != 4 or trace.get("sample_invocations") != 1):
        raise ValueError("actual diffusion task differs from the finite declaration")
    view = json.loads((native_dir / "engine-view.json").read_text())
    req = DiffusionRequest.from_json(json.dumps(view["request"]))
    if len(view["positions"]) != 716 or req.protein_context.preserved_only or len(req.gaps) != 4:
        raise ValueError("actual sampling used incomplete/local conditioning")
    if {c.atom for c in req.protein_context.fixed_conditioning} != set(req.fixed_atoms):
        raise ValueError("source observations became output-only")
    with np.load(native_dir / "native.npz", allow_pickle=False) as data:
        native = data["atom37"].copy()
        native_mask = data["atom_mask"].copy() if "atom_mask" in data else None
    if plan.get("fixed_scaffold_inpainting"):
        if native_mask is None or trace.get("inpainting_policy") != plan["inpainting_policy"]:
            raise ValueError("missing explicit native inpainting/output-mask evidence")
        from scripts.audit_cc94_inpainting import audit_trajectory
        engine = DiffusionRequest.from_json(json.dumps(view["engine_view"]["request"]))
        axis = target_axis(engine)
        row_index, cursor = {}, 0
        for target in engine.target_sequences:
            for residue in axis[target.chain]:
                row_index[(residue.chain, residue.residue_number, residue.insertion_code)] = cursor
                cursor += 1
        masks = []
        for owned in (engine.fixed_atoms, engine.generated_atoms):
            mask = np.zeros((1, cursor, 37), dtype=bool)
            for atom in owned:
                mask[0, row_index[(atom.chain, atom.residue_number, atom.insertion_code)], ATOM37_ORDER.index(atom.atom_name)] = True
            masks.append(mask)
        proof = audit_trajectory(native_dir, *masks, native,
                                 tolerance=plan["trajectory_float_tolerance_angstrom"])
        save(output / "inpainting-audit.json", proof)
        report["inpainting_audit"] = {k: v for k, v in proof.items() if k != "payload_sha256"}
    elif native_mask is not None:
        raise ValueError("an undeclared output-mask policy cannot replace the canonical v1 mask")
    decoded = native_coordinates(view, native, native_mask)
    source = atoms(Path(plan["source"]))
    missing = [a for a in req.fixed_atoms if a not in decoded]
    fixed = [a for a in req.fixed_atoms if a in decoded]
    target = np.array([source[(a.chain, a.residue_number, a.insertion_code, a.atom_name)]["xyz"] for a in fixed])
    predicted = np.array([decoded[a] for a in fixed])
    transform = weighted_kabsch(predicted, target)
    fitted = transform.apply(predicted)
    errors = np.linalg.norm(fitted-target, axis=1)
    maximum = float(errors.max())
    report.update(native_fixed_rms_angstrom=float(np.sqrt(np.mean(errors**2))), native_fixed_max_angstrom=maximum,
        native_fixed_tolerance_angstrom=plan["native_float_tolerance_angstrom"],
        fixed_atom_count=len(req.fixed_atoms), represented_native_fixed_atoms=len(fixed),
        missing_native_observed_atoms=[{"chain": a.chain, "resid": a.residue_number, "icode": a.insertion_code,
                                        "atom": a.atom_name} for a in missing],
        output_only_observed_atoms=0, native_sha256=digest(native_dir / "native.npz"),
        native_fixed_per_chain={chain: float(max(errors[i] for i, a in enumerate(fixed) if a.chain == chain))
                                for chain in sorted({a.chain for a in fixed})},
        native_backbone_rms_angstrom=float(np.sqrt(np.mean([e*e for e, a in zip(errors, fixed, strict=True) if a.atom_name in {"N", "CA", "C", "O"}]))),
        native_ca_max_angstrom=float(max(e for e, a in zip(errors, fixed, strict=True) if a.atom_name == "CA")),
        audit_revision="native-inpainting-mask-aware-v3" if plan.get("fixed_scaffold_inpainting") else "native-residue-mask-aware-v2",
        audit_code_sha256=digest(Path(__file__)))
    if missing or maximum > plan["native_float_tolerance_angstrom"]:
        reasons = []
        if missing:
            reasons.append("observed_atoms_absent_from_native_residue_mask")
        if maximum > plan["native_float_tolerance_angstrom"]:
            reasons.append("native_observation_motion")
        report.update(outcome="native_observation_contract_refused", reasons=reasons,
                      source_restored_candidate_written=False)
        save(output / "receipt.json", report)
        return
    report["native_fixed_passed"] = True
    if not execution or execution["returncode"] or execution["guard_stop_reason"] or trace.get("ode_updates") != 500:
        report.update(outcome="native_retained_but_execution_incomplete", source_restored_candidate_written=False)
        save(output / "receipt.json", report)
        return
    from scripts.audit_anchor_calibration_grid import diagnose

    from dvbfixer.ffutils.dat import DatRecord
    from dvbfixer.model.diffusion.pdb_materialize import materialize_candidate_pdb
    from dvbfixer.model.diffusion.staged_runtime import assessment_adapter

    generated = tuple(req.generated_atoms)
    residue_names = {}
    for gap in req.gaps:
        sequence = next(t.sequence for t in req.target_sequences if t.chain == gap.chain)
        for i, residue in zip(range(gap.target_interval.start, gap.target_interval.stop), gap.generated_residues, strict=True):
            residue_names[(residue.chain, residue.residue_number, residue.insertion_code)] = _CANONICAL_ONE_TO_THREE[sequence[i]]
    names = tuple(residue_names[(a.chain, a.residue_number, a.insertion_code)] for a in generated)
    coordinates = transform.apply(np.array([decoded[a] for a in generated]))
    source_text = Path(plan["source"]).read_text()
    candidate = materialize_candidate_pdb(source_text, generated, coordinates, names, req,
                                         elements=tuple(a.atom_name[0] for a in generated))
    path = output / "raw-protein-only.pdb"
    path.write_text(candidate)
    final = atoms(path)
    if any(final[a]["line"] != record["line"] for a, record in source.items()):
        raise ValueError("materialization changed observed protein records")
    assessment = assessment_adapter(req, source_text)(candidate)
    from dataclasses import asdict
    save(output / "assessment.json", asdict(assessment))
    cli = diagnose(path, output / "full-cli.json")
    record = DatRecord()
    record.added_atoms = [{"chain": a.chain, "resid": a.residue_number, "icode": a.insertion_code,
                           "resname": name, "atom": a.atom_name, "element": a.atom_name[0]}
                          for a, name in zip(generated, names, strict=True)]
    record.save(output / "generated.dat")
    report.update(sampling_success=True, source_restored_candidate_written=True,
        outcome="raw_geometry_pass" if assessment.geometry_passed and not assessment.admission_failures else "raw_geometry_failure",
        protein_mask_geometry_passed=assessment.geometry_passed and not assessment.admission_failures,
        generated_atoms=len(generated), source_records_exact_after_native_gate=True,
        candidate_sha256=digest(path), full_cli_errors=sum(f["severity"] == "ERROR" for f in cli))
    verify_bindings(plan["source_code_and_input_sha256"], (root / "plan.json",))
    save(output / "receipt.json", report)
