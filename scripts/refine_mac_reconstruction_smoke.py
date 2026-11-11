#!/usr/bin/env python3
"""Private shared generated-only repair of an audited proposal, no anchor release."""
from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import asdict
from pathlib import Path

from dvbfixer.model.diffusion.boundary_refinement import refine_generated_region
from dvbfixer.model.diffusion.contract import DiffusionRequest
from dvbfixer.model.diffusion.runner_refinement import rewrite_generated_coordinates
from dvbfixer.model.diffusion.staged_runtime import assessment_adapter


def run(audit_dir, output, *, platform="Reference", geometry_preconditioning=False, force_shared_refinement=False):
    request = DiffusionRequest.from_json((audit_dir / "evaluation-mask.json").read_text())
    text = (audit_dir / "source-restored.pdb").read_bytes().decode()
    source_bytes = (audit_dir / request.normalized_pdb.path).read_bytes()
    source = source_bytes.decode()
    if hashlib.sha256(source_bytes).hexdigest() != request.normalized_pdb.sha256:
        raise ValueError("refinement source differs from the declared request")
    assess = assessment_adapter(request, source)
    before = assess(text)
    output.mkdir(parents=True, exist_ok=False)
    (output / "input.pdb").write_text(text)
    (output / "entry-script.py").write_bytes(Path(__file__).read_bytes())
    implementation = Path(refine_generated_region.__code__.co_filename)
    (output / "boundary-refinement.py").write_bytes(implementation.read_bytes())
    declaration = {"mask": "ALL generated atoms, no observed motion", "seed": 7,
        "input_sha256": hashlib.sha256(text.encode()).hexdigest(), "platform": platform,
        "max_iterations": 1000, "restart_count": 0, "reference_coordinates_supplied": False,
        "geometry_preconditioning": geometry_preconditioning,
        "force_shared_post_assembly_refinement": force_shared_refinement,
        "boundary_refinement_sha256": hashlib.sha256(implementation.read_bytes()).hexdigest(),
        "preconditioning_iteration_caps": [100, 75, 75, 1000] if geometry_preconditioning else None}
    (output / "declaration.json").write_text(json.dumps(declaration, indent=2) + "\n")
    if before.admission_failures:
        (output / "failure.json").write_text(json.dumps({"reason": "raw-admission-failure", "findings": before.admission_failures})+"\n")
        raise ValueError("refinement cannot authorize source/identity admission failures")
    if before.geometry_passed and not before.admission_failures and not force_shared_refinement:
        (output / "status.json").write_text(json.dumps({"status": "raw-clean-no-repair-needed", "resampled": False}) + "\n")
        return

    def score(coordinates):
        result = assess(rewrite_generated_coordinates(text, coordinates))
        return len(result.admission_failures), sum(issue.actionable for issue in result.issues)

    try:
        result = refine_generated_region(text, generated_residues=tuple(r for gap in request.gaps for r in gap.generated_residues),
            generated_atoms=request.generated_atoms, max_iterations=1000, restart_count=0,
            platform_name=platform, random_seed=7, candidate_score=score,
            geometry_preconditioning=geometry_preconditioning)
        if set(result.coordinates_angstrom) != set(request.generated_atoms):
            raise ValueError("generated-only refiner returned a different movable atom mask")
        refined = rewrite_generated_coordinates(text, result.coordinates_angstrom)
        after = assess(refined)
        (output / "refined.pdb").write_text(refined)
        report = {"before": asdict(before), "after": asdict(after), "geometry_passed": after.geometry_passed,
            "source_and_geometry_passed": after.geometry_passed and not after.admission_failures,
            "initial_energy_kj_mol": result.initial_energy_kj_mol, "final_energy_kj_mol": result.final_energy_kj_mol,
            "platform": result.platform, "search_diagnostics": result.search_diagnostics,
            "scientific_qualification": False}
        (output / "status.json").write_text(json.dumps(report, indent=2) + "\n")
        print(json.dumps({"output": str(output), "platform": result.platform,
            "geometry_passed": after.geometry_passed, "admission_failures": list(after.admission_failures),
            "actionable_findings": sum(issue.actionable for issue in after.issues)}), flush=True)
    except Exception as exc:
        (output / "failure.json").write_text(json.dumps({"type": type(exc).__name__, "message": str(exc)}) + "\n")
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--platform", choices=("Reference", "CPU"), default="Reference")
    parser.add_argument("--geometry-preconditioning", action="store_true")
    parser.add_argument("--force-shared-refinement", action="store_true")
    args = parser.parse_args()
    run(args.audit, args.output, platform=args.platform, geometry_preconditioning=args.geometry_preconditioning,
        force_shared_refinement=args.force_shared_refinement)
