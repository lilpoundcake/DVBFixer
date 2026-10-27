#!/usr/bin/env python3
"""Research-only paired cc89 sampling / local neural denoising / OpenMM matrix.

No production profile, acceptance threshold, or original source is modified.
Run sampling with the pinned Apple environment and evaluation with the scientific
environment. See docs/plans/cc89-local-refinement-ablation.md for frozen policies.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import sys
import time
from dataclasses import asdict
from io import StringIO
from pathlib import Path

import numpy as np

from dvbfixer.model.diffusion.contract import (
    AtomIdentity,
    DiffusionRequest,
    ResidueIdentity,
    SamplerTrace,
)
from dvbfixer.model.diffusion.pdb_materialize import materialize_candidate_pdb
from dvbfixer.model.diffusion.protein_context import (
    engine_context_view,
    independent_gap_requests,
    sampling_entities,
    source_context_position_features,
)
from dvbfixer.model.diffusion.runner_refinement import rewrite_generated_coordinates

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cc89_lrd_protocol import (
    boundary_residues,
    low_noise_schedule,
    observed_frame_center,
    relaxation_residues,
)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path: Path, value) -> None:
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n")


def atom_records(text: str):
    records = {}
    for line in text.splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            atom = AtomIdentity(line[21], line[22:26].strip(), line[26].strip(), line[12:16].strip())
            if atom in records:
                raise ValueError(f"ambiguous atom: {atom}")
            records[atom] = (np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]), line[17:20].strip(), line)
    return records


def load_adapter(path: Path):
    spec = importlib.util.spec_from_file_location("cc89_lrd_frozen_adapter", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def lrd_gap(model, adapter, request, gap_index, text, *, include_observed):
    """Denoise only one explicit boundary mask; never condition on another gap."""
    import torch
    from protpardelle.common import residue_constants

    single = independent_gap_requests(request)[gap_index]
    entities = sampling_entities(single)
    residues = tuple(residue for entity in entities for residue in entity.residues)
    row = {residue: index for index, residue in enumerate(residues)}
    assert len(row) == len(residues) <= 512
    axis = tuple(original for _, original in engine_context_view(single).atom_to_source)
    records = atom_records(text)
    active_residues = boundary_residues(single, single.gaps[0], include_observed=include_observed)
    base = np.zeros((1, len(row), 37, 3), dtype=np.float32)
    present = np.zeros(base.shape[:-1], dtype=bool)
    mobile = np.zeros_like(present)
    observed = np.zeros_like(present)
    original_atoms = set(single.fixed_atoms)
    indices = {}
    for atom in axis:
        residue = ResidueIdentity(atom.chain, atom.residue_number, atom.insertion_code)
        if atom.atom_name not in residue_constants.atom_order:
            continue  # source-only atom names remain untouched in the full PDB
        key = (0, row[residue], residue_constants.atom_order[atom.atom_name])
        base[key] = records[atom][0]
        present[key] = True
        mobile[key] = residue in active_residues
        observed[key] = atom in original_atoms
        indices[atom] = key
    if not mobile.any() or np.count_nonzero(present & ~mobile) < 3:
        raise ValueError("LRD requires mobile atoms and immutable context anchors")
    center = observed_frame_center(base, observed)
    centered = np.where(present[..., None], base - center, 0).astype(np.float32)
    # Same noise array at shared atoms for both mobility policies and base arms.
    noise = np.random.default_rng(7 + gap_index).standard_normal(base.shape).astype(np.float32)
    schedule = low_noise_schedule()
    device = torch.device("mps")
    reference = torch.as_tensor(centered, device=device)
    mobile_t = torch.as_tensor(mobile, device=device)
    present_t = torch.as_tensor(present, device=device)
    state = reference + schedule[0] * torch.as_tensor(noise, device=device) * mobile_t[..., None]
    crop = reference.masked_fill(mobile_t[..., None], 0)
    positions, chains = source_context_position_features(single)
    seq_mask = torch.ones((1, len(row)), device=device)
    residue_index = torch.as_tensor([positions], dtype=torch.long, device=device)
    chain_index = torch.as_tensor([chains], dtype=torch.long, device=device)
    control = adapter._NeutralDenoiserChainPhase()
    probe = adapter._SourcePositionProbe(positions, chains, device)
    rope = model.struct_model.net.transformer.layers[0][0].rope
    steps = []
    self_cond = None
    start = time.perf_counter()
    with torch.inference_mode(), adapter._probe_hook(model.struct_model, control), adapter._rope_probe(rope, probe):
        for index, (sigma, next_sigma) in enumerate(zip(schedule[:-1], schedule[1:], strict=True)):
            static_input_error = float((state[~mobile_t] - reference[~mobile_t]).abs().max().item())
            if static_input_error != 0:
                raise ValueError("LRD moved an inactive atom before denoising")
            denoised, _, next_self_cond, _ = model.forward(
                noisy_coords=state, noise_level=torch.tensor([sigma], device=device, dtype=torch.float32),
                seq_mask=seq_mask, residue_index=residue_index, chain_index=chain_index,
                hotspot_mask=torch.zeros_like(residue_index), struct_crop_cond=crop,
                struct_self_cond=self_cond if model.config.train.self_cond_train_prob > 0.5 else None,
                run_mpnn_model=False)
            proposed = state + 1.2 * (float(next_sigma - sigma) / float(sigma)) * (state - denoised)
            state = torch.where(mobile_t[..., None], proposed, reference)
            state = torch.where(present_t[..., None], state, torch.zeros_like(state))
            self_cond = next_self_cond
            if not bool(torch.isfinite(state).all()):
                raise ValueError("non-finite LRD state")
            steps.append({"index": index, "sigma": float(sigma), "next_sigma": float(next_sigma),
                          "inactive_max_error_angstrom": static_input_error})
    torch.mps.synchronize()
    assert control.calls == 25 and probe.calls == 50 and probe.result is not None
    output = state.cpu().numpy() + center
    updates = {atom: output[key].astype(float) for atom, key in indices.items() if mobile[key]}
    report = {"gap_index": gap_index, "window": 3, "include_observed": include_observed,
        "frame_policy": "same-observed-source-centroid-for-both-mobility-arms", "frame_center": center.tolist(),
        "noise_seed": 7 + gap_index, "steps": steps, "step_scale": 1.2, "extra_churn": 0,
        "neutral_phase_calls": control.calls, "first_layer_rope_calls": probe.calls,
        "mobile_atoms": [asdict(atom) for atom in updates],
        "represented_atoms": len(indices), "wall_seconds": time.perf_counter() - start}
    return updates, report


def sample(args, request):
    import protpardelle.core.models as engine
    import torch

    root = args.workspace
    if (root / "initial").exists():
        raise FileExistsError("do not overwrite an existing sampling experiment")
    adapter = load_adapter(args.adapter)
    source = (root / request.normalized_pdb.path).read_text()
    manifest = {"research_only": True, "seed": 7, "initial_steps": 500,
        "initial_step_scale": 1.2, "initial_s_churn": 200.0,
        "lrd_steps": 25, "lrd_sigma_start": 1.5, "lrd_sigma_min": 0.0064,
        "lrd_step_scale": 1.2, "lrd_churn": 0, "lrd_rho": 7, "window": 3,
        "reinjection_factor": "initial generation only; inactive LRD mask always projected",
        "checkpoint_sha256": digest(args.checkpoint), "config_sha256": digest(args.config),
        "adapter_sha256": digest(args.adapter), "engine_models_sha256": digest(Path(engine.__file__)),
        "script_sha256": digest(Path(__file__)), "request_sha256": digest(root / "request.json"),
        "source_sha256": digest(root / request.normalized_pdb.path), "torch": torch.__version__,
        "device": "mps", "python": sys.executable, "production_publication": False}
    write_json(root / "manifest.json", manifest)
    for reinjection in (False, True):
        arm = "reinjection-on" if reinjection else "reinjection-off"
        directory = root / "initial" / arm
        directory.mkdir(parents=True)
        updates, names, evidence = {}, {}, []
        for index, single in enumerate(independent_gap_requests(request)):
            path = root / f"{arm}-gap-{index}.json"
            path.write_text(single.to_json())
            output = directory / f"gap-{index}"
            output.mkdir()
            result, _ = adapter._run_staged(path, output, args.config, args.checkpoint,
                steps=500, step_scale=1.2, s_churn=200.0, per_step_reinjection=reinjection,
                device_name="mps", mps_profile=False, chain_phase_policy="neutral-v1",
                position_feature_policy="source-physical-v1", corrected_mps_research=True)
            write_json(output / "runner-result.json", result.to_dict())
            trace_path = root / result.candidates[0].sampler_trace_artifact.path
            trace = SamplerTrace.from_json(trace_path.read_text())
            assert trace.callback_update_count == (500 if reinjection else 0)
            assert trace.denoising_update_count == 500
            if reinjection:
                assert trace.max_projection_error_angstrom == 0
            records = atom_records((root / result.candidates[0].coordinate_artifact.path).read_text())
            original_records = atom_records(source)
            assert all(records[atom][2] == original_records[atom][2] for atom in original_records)
            for atom in single.generated_atoms:
                assert atom not in updates
                updates[atom], names[atom] = records[atom][:2]
            evidence.append({"gap_index": index, "callback_updates": trace.callback_update_count,
                             "trace_sha256": digest(trace_path)})
        assert set(updates) == set(request.generated_atoms)
        text = materialize_candidate_pdb(source, request.generated_atoms,
            [updates[a] for a in request.generated_atoms], tuple(names[a] for a in request.generated_atoms), request)
        (directory / "candidate.pdb").write_text(text)
        write_json(directory / "evidence.json", evidence)
        print("INITIAL COMPLETE", arm, flush=True)
    run_lrd(args, request, adapter)


def run_lrd(args, request, adapter=None):
    import protpardelle.core.models as engine
    from protpardelle.core.models import load_model

    root = args.workspace
    if (root / "lrd").exists():
        raise FileExistsError("do not overwrite existing LRD evidence")
    adapter = adapter or load_adapter(args.adapter)
    manifest = json.loads((root / "manifest.json").read_text())
    assert manifest['engine_models_sha256'] == digest(Path(engine.__file__))
    assert manifest['checkpoint_sha256'] == digest(args.checkpoint)
    assert manifest['config_sha256'] == digest(args.config)
    assert manifest['request_sha256'] == digest(root / 'request.json')
    write_json(root / "lrd-manifest.json", {"script_sha256": digest(Path(__file__)),
        "protocol_sha256": digest(Path(__file__).with_name('cc89_lrd_protocol.py')),
        "frame_policy": "same-observed-source-centroid-for-both-mobility-arms",
        "initial_candidates": {arm: digest(root / 'initial' / arm / 'candidate.pdb')
                               for arm in ('reinjection-off', 'reinjection-on')}})
    model = load_model(args.config, args.checkpoint, device="mps")
    model.eval()
    for arm in ("reinjection-off", "reinjection-on"):
        base_text = (root / "initial" / arm / "candidate.pdb").read_text()
        for include_observed in (False, True):
            mobility = "observed-mobile" if include_observed else "generated-only"
            directory = root / "lrd" / arm / mobility
            directory.mkdir(parents=True)
            updates, evidence = {}, []
            for index in range(len(request.gaps)):
                result, trace = lrd_gap(model, adapter, request, index, base_text, include_observed=include_observed)
                if set(result) & set(updates):
                    raise ValueError("overlapping LRD masks require an explicit policy; joint gaps are deferred")
                updates.update(result)
                evidence.append(trace)
            (directory / "candidate.pdb").write_text(rewrite_generated_coordinates(base_text, updates))
            write_json(directory / "lrd-trace.json", evidence)
            print("LRD COMPLETE", arm, mobility, flush=True)


def evaluate(args, request):
    from openmm.app import PDBFile

    from dvbfixer.diagnose.chemistry import (
        check_backbone_bond_angles,
        check_bond_lengths,
        check_peptide_omegas,
    )
    from dvbfixer.diagnose.report import Severity
    from dvbfixer.diagnose.steric import clashes_python
    from dvbfixer.model.diffusion.boundary_refinement import (
        BOUNDARY_REFINEMENT_RESTART_COUNT,
        BOUNDARY_REFINEMENT_REVISION,
        refine_generated_region,
    )
    from dvbfixer.model.diffusion.quality import check_ramachandran, check_sidechain_chi12
    from dvbfixer.model.diffusion.validate import (
        ValidationThresholds,
        _parse_pdb_text,
        _Structure,
        _validate_candidate,
    )

    root = args.workspace
    include_anchors = args.minimize_anchors
    report_path = root / ("comparison-with-anchors.json" if include_anchors else "comparison.json")
    if report_path.exists():
        raise FileExistsError(f"do not overwrite an existing comparison: {report_path}")
    def structure(text, path):
        atoms, residues, links = _parse_pdb_text(text)
        return _Structure(path, PDBFile(StringIO(text)), atoms, residues, frozenset(links), text)

    source_text = (root / request.normalized_pdb.path).read_text()
    source = structure(source_text, root / request.normalized_pdb.path)
    relaxed_residues = relaxation_residues(request, include_anchors=include_anchors)
    anchors = {r for gap in request.gaps for r in (gap.left_anchor, gap.right_anchor) if r is not None}
    local = {(r.chain, r.residue_number + r.insertion_code) for gap in request.gaps
             for r in (*gap.generated_residues, gap.left_anchor, gap.right_anchor) if r is not None}
    report = []

    def assess(path, arm, processing, mobility, refinement=None):
        text = path.read_text()
        candidate = structure(text, path)
        summary, _ = _validate_candidate(request, source, candidate, ValidationThresholds())
        displacements = np.array([np.linalg.norm(candidate.atoms[a].coordinates_angstrom - source.atoms[a].coordinates_angstrom)
                                  for a in request.fixed_atoms])
        anchor_mask = np.array([ResidueIdentity(a.chain, a.residue_number, a.insertion_code) in anchors
                                for a in request.fixed_atoms])
        geometry_failures = [gate for gate in summary.hard_gate_failures
                             if gate not in {"fixed-heavy-atom-rmsd", "fixed-heavy-atom-max-displacement"}]
        # Mobile observed atoms extend beyond the normal generated-only validator
        # neighborhood. Audit the whole released window plus one residue on each
        # side, including its NEW outer junctions, without relaxing normal gates.
        impacted = set()
        for gap in request.gaps:
            selected = boundary_residues(request, gap, include_observed=mobility == "observed-mobile")
            if include_anchors:
                selected = selected | {r for r in (gap.left_anchor, gap.right_anchor) if r is not None}
            placement = next(p for p in request.sequence_placements if p.chain == gap.chain)
            positions = dict(zip(placement.observed_target_indices, placement.observed_residues, strict=True))
            positions.update(zip(range(gap.target_interval.start, gap.target_interval.stop), gap.generated_residues, strict=True))
            selected_indices = {i for i, residue in positions.items() if residue in selected}
            expanded = {j for i in selected_indices for j in (i - 1, i, i + 1)}
            impacted.update((r.chain, r.residue_number + r.insertion_code) for i, r in positions.items() if i in expanded)
        expanded_findings = []
        for check in (check_bond_lengths, check_backbone_bond_angles, check_peptide_omegas,
                      check_ramachandran, check_sidechain_chi12):
            for finding in check(candidate.pdb.topology, candidate.pdb.positions):
                if (finding.chain, finding.resid) not in impacted:
                    continue
                if finding.severity is Severity.ERROR or finding.category == "non_planar_amide":
                    expanded_findings.append(asdict(finding))
        for finding in clashes_python(candidate.pdb.topology, candidate.pdb.positions, clash_warn_a=0.9, clash_error_a=0.9):
            partner = finding.extra.get("clash_partner", {})
            if (finding.chain, finding.resid) in impacted or (partner.get("chain"), partner.get("resid")) in impacted:
                expanded_findings.append(asdict(finding))
        row = {"arm": arm, "processing": processing, "mobility": mobility,
            "minimize_anchors": include_anchors,
            "coordinate_path": str(path.relative_to(root)), "sha256": digest(path),
            "strict_validation": asdict(summary), "geometry_failures_excluding_source_drift": geometry_failures,
            "expanded_mobile_window_geometry_findings": expanded_findings,
            "geometry_only_passed_not_acceptance": not geometry_failures and not expanded_findings,
            "observed_max_displacement_angstrom": float(displacements.max()),
            "observed_rmsd_angstrom": float(np.sqrt(np.mean(displacements ** 2))),
            "observed_atoms_moved_over_0_01a": int(np.count_nonzero(displacements > 0.01)),
            "anchor_max_displacement_angstrom": float(displacements[anchor_mask].max(initial=0)),
            "non_anchor_observed_max_displacement_angstrom": float(displacements[~anchor_mask].max(initial=0)),
            "local_omega_findings": [asdict(f) for f in check_peptide_omegas(candidate.pdb.topology, candidate.pdb.positions)
                                     if (f.chain, f.resid) in local], "refinement": refinement}
        report.append(row)
        write_json(report_path, report)
        print(arm, processing, mobility, "STRICT", summary.passed,
              "GEOMETRY", geometry_failures, "DRIFT", row["observed_max_displacement_angstrom"], flush=True)

    for arm in ("reinjection-off", "reinjection-on"):
        for mobility, path in [("not-applicable", root / "initial" / arm / "candidate.pdb"),
                               *[(mode, root / "lrd" / arm / mode / "candidate.pdb")
                                 for mode in ("generated-only", "observed-mobile")]]:
            if mobility != "not-applicable" and not include_anchors:
                assess(path, arm, "LRD-only", mobility)
            text = path.read_text()
            parsed_atoms, _, _ = _parse_pdb_text(text)
            relaxed_atoms = tuple(dict.fromkeys((*request.generated_atoms, *(
                atom for atom in request.fixed_atoms
                if include_anchors and ResidueIdentity(atom.chain, atom.residue_number, atom.insertion_code) in anchors
                and parsed_atoms[atom].element.upper() not in {"H", "D"}))))
            directory = root / ("processed-with-anchors" if include_anchors else "processed") / arm / mobility
            directory.mkdir(parents=True, exist_ok=False)
            started = time.perf_counter()
            try:
                refined = refine_generated_region(text, generated_residues=relaxed_residues,
                    generated_atoms=relaxed_atoms, max_iterations=500,
                    restart_count=BOUNDARY_REFINEMENT_RESTART_COUNT, perturbation_angstrom=0.25,
                    platform_name="CPU", random_seed=7)
            except Exception as exc:
                failure = {"arm": arm, "processing": "OpenMM-only" if mobility == "not-applicable" else "LRD→OpenMM",
                    "mobility": mobility, "refinement_failed": True, "error": f"{type(exc).__name__}: {exc}",
                    "minimize_anchors": include_anchors,
                    "input_path": str(path.relative_to(root))}
                write_json(directory / "failure.json", failure)
                report.append(failure)
                write_json(report_path, report)
                print("REFINEMENT FAILED", failure, flush=True)
                continue
            result_path = directory / "candidate.pdb"
            result_text = rewrite_generated_coordinates(text, refined.coordinates_angstrom)
            before_records, after_records = atom_records(text), atom_records(result_text)
            allowed = set(relaxed_atoms)
            if any(before_records[a][2] != after_records[a][2] for a in before_records if a not in allowed):
                raise ValueError("OpenMM output changed an atom outside the explicit minimization mask")
            result_path.write_text(result_text)
            params = {"revision": BOUNDARY_REFINEMENT_REVISION, "platform": refined.platform,
                "max_iterations": 500, "restart_count": BOUNDARY_REFINEMENT_RESTART_COUNT,
                "perturbation_angstrom": 0.25, "seed": 7, "wall_seconds": time.perf_counter() - started,
                "initial_energy_kj_mol": refined.initial_energy_kj_mol, "final_energy_kj_mol": refined.final_energy_kj_mol,
                "mobile_residues": [asdict(r) for r in relaxed_residues], "mobile_heavy_atom_count": len(relaxed_atoms),
                "minimize_anchors": include_anchors,
                "observed_atoms_frozen_at": ("except immediate anchors; " if include_anchors else "")
                    + ("LRD output" if mobility != "not-applicable" else "original source")}
            write_json(directory / "refinement.json", params)
            assess(result_path, arm, "OpenMM-only" if mobility == "not-applicable" else "LRD→OpenMM", mobility, params)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("phase", choices=("sample", "lrd", "evaluate"))
    parser.add_argument("--workspace", required=True, type=Path)
    parser.add_argument("--adapter", type=Path)
    parser.add_argument("--config", type=Path)
    parser.add_argument("--checkpoint", type=Path)
    parser.add_argument("--minimize-anchors", action="store_true",
                        help="Research evaluation only: also minimize immediate left/right anchor residues.")
    args = parser.parse_args()
    args.workspace = args.workspace.resolve()
    if args.minimize_anchors and args.phase != "evaluate":
        parser.error("--minimize-anchors applies only to evaluation; sampling remains unchanged")
    request = DiffusionRequest.from_json((args.workspace / "request.json").read_text())
    if request.seeds != (7,) or request.candidate_count != 1:
        parser.error("this frozen experiment requires a single seed-7 candidate")
    if args.phase in {"sample", "lrd"}:
        if any(value is None for value in (args.adapter, args.config, args.checkpoint)):
            parser.error("sampling requires --adapter, --config and --checkpoint")
        (sample if args.phase == "sample" else run_lrd)(args, request)
    else:
        evaluate(args, request)


if __name__ == "__main__":
    main()
