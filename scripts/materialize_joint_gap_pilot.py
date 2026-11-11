#!/usr/bin/env python3
"""Independently rederive a private joint proposal from native arrays, then refine.

The saved full prediction is an output-only graft frame. No observed atom or
other predicted gap can move. Pair-local geometry is not topological/full-model
qualification; all whole-complex findings remain in the retained reports.
"""
from __future__ import annotations

import argparse
import hashlib
from dataclasses import asdict
from pathlib import Path


def main():
    import numpy as np

    from dvbfixer.ffutils.geometry import assert_all_l
    from dvbfixer.model.diffusion.anchor_acceptance import _structure
    from dvbfixer.model.diffusion.boundary_refinement import refine_generated_region
    from dvbfixer.model.diffusion.contract import DiffusionRequest
    from dvbfixer.model.diffusion.gap_portfolio import _coordinates, digest
    from dvbfixer.model.diffusion.geometry import weighted_kabsch
    from dvbfixer.model.diffusion.joint_gap_research import (
        COMPLETE_REVISION,
        MOTIF_COMPONENT_REVISION,
        REVISION,
        STRUCTURE_REVISION,
        build_joint_gap_view,
        build_joint_structure_view,
        complete_context_admission,
        motif_component_admission,
        plan_complete_context,
        plan_motif_component,
    )
    from dvbfixer.model.diffusion.portfolio_evidence import json_bytes, strict_json
    from dvbfixer.model.diffusion.private_artifact import read_private_file_bytes
    from dvbfixer.model.diffusion.protein_atoms import ATOM37_ORDER
    from dvbfixer.model.diffusion.runner_refinement import rewrite_generated_coordinates
    from dvbfixer.model.diffusion.staged_runtime import assessment_adapter

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--spatial-shell", action="store_true",
        help="Explicit full-structure all-generated + bounded observed 5 A shell policy")
    parser.add_argument("--max-iterations", type=int, default=250)
    parser.add_argument("--restart-count", type=int, default=4)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    directory = args.output.resolve(strict=True)
    root = args.sample.resolve(strict=True)

    def write(name, value):
        data = value if isinstance(value, bytes) else json_bytes(value)
        with (directory / name).open("xb") as stream:
            stream.write(data)

    def read(name):
        return read_private_file_bytes(root, name, max_bytes=16_000_000)

    write("run-script.py", Path(__file__).read_bytes())
    write("protocol.json", {"scope": __doc__, "max_iterations": args.max_iterations, "restart_count": args.restart_count,
        "platform": "CPU", "observed_spatial_mobility": args.spatial_shell, "new_sampling": False})
    try:
        trace = strict_json(read("trace.json"))
        if (trace.get("kind") not in {REVISION, STRUCTURE_REVISION, COMPLETE_REVISION, MOTIF_COMPONENT_REVISION} or trace.get("status") != "native-proposal-retained"
                or trace.get("sample_invocations") != 1 or trace.get("ode_updates") != 500
                or trace.get("public_profile_enabled") is not False):
            raise ValueError("joint sampling trace is incomplete or not a single shared trajectory")
        parent_bytes, source, baseline = read("parent-request.json"), read("source.pdb"), read("baseline-output-only.pdb")
        parent = DiffusionRequest.from_dict(strict_json(parent_bytes))
        if trace["kind"] == MOTIF_COMPONENT_REVISION:
            view = plan_motif_component(parent, source, tuple(trace["gap_indices"]), seed=trace["seed"])
            if not motif_component_admission(view)["supported"]:
                raise ValueError("neighbor motif component exceeds its execution admission")
            if args.spatial_shell:
                raise ValueError("component sampling requires assembly before all-generated spatial relaxation")
        elif trace["kind"] == COMPLETE_REVISION:
            view = plan_complete_context(parent, source, tuple(trace["gap_indices"]), seed=trace["seed"])
            if not complete_context_admission(view)["supported"]:
                raise ValueError("complete input exceeds its current execution admission")
        elif trace["kind"] == STRUCTURE_REVISION:
            if trace["gap_indices"] != list(range(len(parent.gaps))):
                raise ValueError("full-structure trace does not cover all gaps")
            view = build_joint_structure_view(parent, source, seed=trace["seed"])
        else:
            if args.spatial_shell:
                raise ValueError("spatial shell policy requires the explicit full-structure research operation")
            view = build_joint_gap_view(parent, source, tuple(trace["gap_indices"]), seed=trace["seed"])
        if (digest(parent.to_json()) != trace["parent_request_sha256"]
                or hashlib.sha256(source).hexdigest() != trace["source_sha256"]
                or hashlib.sha256(baseline).hexdigest() != trace["baseline_output_only_sha256"]
                or read("research-group.json") != view.request.to_json().encode()
                or strict_json(read("generated-axis.json")) != [asdict(a) for a in view.request.generated_atoms]):
            raise ValueError("joint parent/source/output frame/mask custody differs")
        native_bytes = read("native.npz")
        if hashlib.sha256(native_bytes).hexdigest() != trace["native_sha256"]:
            raise ValueError("joint native bytes differ")
        from io import BytesIO

        with np.load(BytesIO(native_bytes), allow_pickle=False) as archive:
            if set(archive.files) != {"atom37"}:
                raise ValueError("joint native archive has unexpected arrays")
            native = archive["atom37"]
        if native.shape != (len(view.positions), 37, 3) or not np.isfinite(native).all():
            raise ValueError("joint native array shape/values differ")
        # Independent physical residue-to-token mapping, no sampler mapping receipt.
        engine = view.engine_view.request
        axis = {}
        cursor = 0
        for target in engine.target_sequences:
            placement = next(p for p in engine.sequence_placements if p.chain == target.chain)
            mapped = dict(zip(placement.observed_target_indices, placement.observed_residues, strict=True))
            for gap in engine.gaps:
                if gap.chain == target.chain:
                    mapped.update(zip(range(gap.target_interval.start, gap.target_interval.stop), gap.generated_residues, strict=True))
            if set(mapped) != set(range(len(target.sequence))):
                raise ValueError("joint engine does not map every token exactly")
            axis.update({(r.chain, r.residue_number, r.insertion_code): cursor + i for i, r in mapped.items()})
            cursor += len(target.sequence)
        source_xyz = _coordinates(source.decode())[0]
        mapping = view.engine_view.atom_to_source

        def point(atom):
            return native[axis[atom.chain, atom.residue_number, atom.insertion_code], ATOM37_ORDER.index(atom.atom_name)]

        # OXT is conditioned/preserved but not part of the model's canonical
        # predicted-heavy mask. Never align its dummy native coordinate.
        fixed = [(e, s) for e, s in mapping if s in set(parent.fixed_atoms) and e.atom_name != "OXT"]
        fit = weighted_kabsch(np.asarray([point(e) for e, _ in fixed]), np.asarray([source_xyz[s] for _, s in fixed]))
        generated = set(view.request.generated_atoms)
        coordinates = {s: fit.apply(np.asarray([point(e)]))[0] for e, s in mapping if s in generated}
        with np.load(BytesIO(read("generated.npz")), allow_pickle=False) as archive:
            if set(archive.files) != {"coordinates"} or not np.allclose(
                    archive["coordinates"], np.asarray([coordinates[a] for a in view.request.generated_atoms]), rtol=0, atol=1e-8):
                raise ValueError("joint generated proposal differs from independent native materialization")
        source_text, baseline_text = source.decode(), baseline.decode()
        assess = assessment_adapter(parent, source_text)
        if assess(baseline_text).admission_failures:
            raise ValueError("joint output graft frame changed source observations")
        before, rows = _coordinates(baseline_text)
        if not generated <= before.keys():
            raise ValueError("joint output frame lacks requested generated atoms")
        raw = rewrite_generated_coordinates(baseline_text, coordinates)
        write("joint_raw_full_complex.pdb", raw.encode())
        write("raw-assessment.json", asdict(assess(raw)))
        write("independent-native-audit.json", {"passed_native_materialization": True,
            "raw_sha256": digest(raw), "sample_invocations": 1, "ode_updates": 500,
            "group_residues": [asdict(r) for g in view.request.gaps for r in g.generated_residues],
            "generated_atoms": len(generated), "other_records_unchanged": True})
        affected = {r for g in view.request.gaps for r in (*g.generated_residues, g.left_anchor, g.right_anchor) if r is not None}
        restraints = ()
        allowed_observed = frozenset()
        mobility = None
        movements = ()
        movement_passed = True
        trials = []
        if args.spatial_shell:
            from dvbfixer.model.diffusion.anchor_policy import (
                AnchorMovementLimits,
                check_anchor_movements,
            )
            from dvbfixer.model.diffusion.spatial_relaxation import (
                audit_spatial_relaxation,
                spatial_mobility,
                spatial_restraints,
            )

            mobility = spatial_mobility(parent, source_text, raw)
            restraints = spatial_restraints(parent, source_text, mobility)
            allowed_observed = frozenset(r.atom for r in restraints)
            write("spatial-mobility.json", {"policy": asdict(mobility), "restraints": [asdict(r) for r in restraints],
                "limits": asdict(AnchorMovementLimits()), "generated_atoms": [asdict(a) for a in parent.generated_atoms]})

        def score(xyz):
            projected = rewrite_generated_coordinates(raw, xyz)
            a = assess(projected, allowed_observed)
            local = sum(i.actionable and bool(set(i.residues) & affected) for i in a.issues)
            write(f"refinement-trial-{len(trials):03d}.pdb", projected.encode())
            trials.append({"assessment": asdict(a), "coordinate_sha256": digest(projected)})
            if args.spatial_shell:
                movement = check_anchor_movements(restraints, _coordinates(projected)[0], AnchorMovementLimits())
                return len(a.admission_failures), sum(not m.within_limits for m in movement), sum(i.actionable for i in a.issues)
            return len(a.admission_failures), local, sum(i.actionable for i in a.issues)

        extra = {"observed_anchor_restraints": restraints, "anchor_movement_limits": AnchorMovementLimits()} if args.spatial_shell else {}
        result = refine_generated_region(raw,
            generated_residues=tuple(r for g in view.request.gaps for r in g.generated_residues),
            generated_atoms=view.request.generated_atoms, max_iterations=args.max_iterations, restart_count=args.restart_count,
            random_seed=trace["seed"], perturbation_angstrom=0.25, platform_name="CPU",
            geometry_preconditioning=True, candidate_score=score, **extra)
        final = rewrite_generated_coordinates(raw, result.coordinates_angstrom)
        final_xyz, final_rows = _coordinates(final)
        mobile = generated | set(allowed_observed)
        if final != rewrite_generated_coordinates(baseline_text, {a: final_xyz[a] for a in mobile}):
            raise ValueError("joint refinement changed records outside the declared group")
        structure = _structure(final, "joint-refined.pdb")
        assert_all_l(structure.pdb.topology, structure.pdb.positions)
        a = assess(final, allowed_observed)
        if args.spatial_shell:
            mobility, movements, a, movement_passed = audit_spatial_relaxation(parent, source_text, raw, final)
            write("spatial-audit.json", {"policy": asdict(mobility), "movements": [asdict(m) for m in movements],
                "assessment": asdict(a), "passed": movement_passed, "global_all_l": True})
        if a.admission_failures:
            raise ValueError("joint refined output has admission failures")
        write("joint_refined_full_complex.pdb", final.encode())
        local = [asdict(i) for i in a.issues if i.actionable and set(i.residues) & affected]
        write("refined-assessment.json", asdict(a))
        write("summary.json", {"scope": ("private full-structure spatial refinement; NOT topology qualification" if args.spatial_shell
            else "private pair-local geometry control, NOT topology or whole-complex acceptance"),
            "pair_geometry_passed": not local, "remaining_pair_findings": local,
            "whole_complex_geometry_passed": a.geometry_passed, "movement_and_geometry_passed": movement_passed and a.geometry_passed,
            "observed_spatial_mobility": args.spatial_shell,
            "observed_shell_residue_count": len(mobility.observed_residues) if mobility else 0,
            "observed_movement_failures": [asdict(m) for m in movements if not m.within_limits], "raw_sha256": digest(raw),
            "refined_sha256": digest(final), "sample_invocations": 1, "ode_updates": 500,
            "generated_atoms": len(generated), "records_outside_mobile_mask_unchanged": True,
            "observed_atoms_unchanged": not args.spatial_shell,
            "global_all_l": True, "topological_qualification": False, "published": False,
            "scientific_qualification": False, "search_diagnostics": result.search_diagnostics, "trials": trials})
        print((directory / "summary.json").read_text(), flush=True)
    except Exception as exc:
        write("failure.json", {"type": type(exc).__name__, "error": str(exc)})
        raise


if __name__ == "__main__":
    main()
