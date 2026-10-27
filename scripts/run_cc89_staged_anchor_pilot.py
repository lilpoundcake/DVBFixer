"""Research-only replay of raw-based geometry-gated refinement on saved 9CQD."""

from __future__ import annotations

import argparse
import hashlib
import shutil
from dataclasses import asdict
from pathlib import Path

from scripts.cc89_anchor_restraints import SOURCE_SHA, assess, residue
from scripts.run_cc89_lrd_ablation import digest, rewrite_generated_coordinates, write_json

from dvbfixer.model.diffusion.anchor_policy import ANCHOR_POLICY_REVISION
from dvbfixer.model.diffusion.boundary_refinement import refine_generated_region
from dvbfixer.model.diffusion.contract import DiffusionRequest, ResidueIdentity
from dvbfixer.model.diffusion.refinement_stages import (
    RepairAssessment,
    RepairIssue,
    run_refinement_stages,
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    request = DiffusionRequest.from_json((args.baseline / "request.json").read_text())
    source_path = args.baseline / request.normalized_pdb.path
    if digest(source_path) != SOURCE_SHA:
        raise ValueError("unexpected source")
    source_text = source_path.read_text()
    args.output.mkdir(parents=True, exist_ok=False)
    shutil.copy2(source_path, args.output / "source.pdb")
    shutil.copy2(args.baseline / "request.json", args.output / "request.json")
    code = args.output / "code"
    code.mkdir()
    from dvbfixer.model.diffusion import anchor_policy, boundary_refinement, refinement_stages
    paths = [Path(__file__), Path(anchor_policy.__file__), Path(boundary_refinement.__file__),
             Path(refinement_stages.__file__), Path("scripts/cc89_anchor_restraints.py")]
    for path in paths:
        shutil.copy2(path, code / path.name)
    write_json(args.output / "manifest.json", {"policy": ANCHOR_POLICY_REVISION,
        "source_sha256": SOURCE_SHA, "request_sha256": digest(args.baseline / "request.json"),
        "code_sha256": {path.name: digest(path) for path in paths}, "platform": "CPU", "threads": 1,
        "solvent": False, "seed": 7, "iterations": 500, "restarts": 8, "perturbation_A": 0.25,
        "publication_enabled": False, "uniform_control": "B weighting disabled; anchors still movable"})
    report = []
    for arm in ("reinjection-off", "reinjection-on"):
        for b_policy in ("disabled", "auto"):
            directory = args.output / arm / b_policy
            directory.mkdir(parents=True)
            raw_path = args.baseline / "initial" / arm / "candidate.pdb"
            raw_text = raw_path.read_text()
            shutil.copy2(raw_path, directory / "raw.pdb")
            state = {"mobile": set(request.generated_atoms)}
            generated_residues = tuple(dict.fromkeys(residue(atom) for atom in request.generated_atoms))
            known = {residue(atom) for atom in request.fixed_atoms + request.generated_atoms}
            lookup = {(r.chain, r.residue_number + r.insertion_code): r for r in known}
            assessment_rows = []

            def assessment(text, directory=directory, assessment_rows=assessment_rows, state=state, lookup=lookup):
                path = directory / f"assessment-{len(assessment_rows):02d}.pdb"
                path.write_text(text)
                metrics = assess(path, request, source_text, state["mobile"])
                assessment_rows.append(metrics)
                issues = []
                for f in metrics["expanded_geometry_findings"]:
                    affected = set()
                    for atom in f.get("extra", {}).get("atoms", []):
                        r = ResidueIdentity(atom["chain"], atom["resid"], atom.get("icode", ""))
                        affected.add(r)
                    for item in (f, f.get("extra", {}).get("clash_partner", {})):
                        r = lookup.get((item.get("chain"), item.get("resid")))
                        if r is not None:
                            affected.add(r)
                    issues.append(RepairIssue(f["category"], tuple(sorted(affected)), f["message"]))
                admission = []
                for gate in metrics["geometry_gates_excluding_source_drift"]:
                    if gate.startswith("generated-or-junction-") or gate in {
                        "junction-peptide-connectivity", "generated-backbone-break", "severe-steric-overlap", "d-ca-chirality"}:
                        # A gate without a detailed finding still blocks success,
                        # but never guesses an observed anchor identity.
                        issues.append(RepairIssue(gate, (), gate))
                    else:
                        admission.append(gate)
                return RepairAssessment(tuple(issues), tuple(admission))

            def minimize(text, restraints, state=state, generated_residues=generated_residues):
                state["mobile"] = set(request.generated_atoms) | {item.atom for item in restraints}
                result = refine_generated_region(text, generated_residues=generated_residues,
                    generated_atoms=request.generated_atoms, observed_anchor_restraints=restraints,
                    max_iterations=500, restart_count=8, perturbation_angstrom=0.25,
                    platform_name="CPU", random_seed=7)
                return rewrite_generated_coordinates(text, result.coordinates_angstrom)

            result = run_refinement_stages(raw_text, source_text, request, assess=assessment,
                minimize=minimize, bfactor_source=b_policy, notify=lambda line: print(line, flush=True))
            attempts = []
            for index, attempt in enumerate(result.attempts):
                row = asdict(attempt)
                text = row.pop("output_text")
                if text is not None:
                    path = directory / f"stage-{index}-{attempt.stage}.pdb"
                    path.write_text(text)
                    row["coordinate_file"] = str(path.relative_to(args.output))
                    row["coordinate_sha256"] = digest(path)
                attempts.append(row)
            outcome = {"arm": arm, "bfactor_policy": b_policy,
                "raw_sha256": hashlib.sha256(raw_text.encode("ascii")).hexdigest(),
                "selected_stage": result.selected_stage,
                "selected_anchors": [asdict(r) for r in result.selected_anchors],
                "requires_anchor_acceptance_contract": result.requires_anchor_acceptance_contract,
                "attempts": attempts, "assessments": assessment_rows, "publication_enabled": False}
            write_json(directory / "report.json", outcome)
            report.append(outcome)
            write_json(args.output / "comparison.json", report)
            print(arm, b_policy, result.selected_stage, "anchors", result.selected_anchors,
                  "geometry", assessment_rows[-1]["geometry_only_pass_not_acceptance"],
                  "RMSD", assessment_rows[-1]["combined_anchor_backbone_rmsd_W_aligned_A"], flush=True)
    write_json(args.output / "completed.json", {"runs": len(report), "publication_enabled": False})


if __name__ == "__main__":
    main()
