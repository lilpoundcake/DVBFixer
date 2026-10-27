"""Frozen research replay of all eight corrected-cc89 pilot raw candidates."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import shutil
import time
import traceback
from collections import Counter
from dataclasses import asdict
from pathlib import Path

from openmm import OpenMMException

from dvbfixer.ffutils.geometry import ChiralityError
from dvbfixer.model.diffusion import (
    anchor_acceptance,
    anchor_policy,
    boundary_refinement,
    refinement_stages,
)
from dvbfixer.model.diffusion.anchor_acceptance import (
    _assessment,
    _observed_lines,
    _structure,
    audit_staged_anchor_repair,
)
from dvbfixer.model.diffusion.anchor_policy import AnchorMovementLimits, NumericalRefinementError
from dvbfixer.model.diffusion.boundary_refinement import (
    BoundaryRefinementError,
    refine_generated_region,
)
from dvbfixer.model.diffusion.contract import DiffusionRequest, RunnerResult
from dvbfixer.model.diffusion.pdb_metadata import reconcile_pdb_metadata
from dvbfixer.model.diffusion.refinement_stages import RepairAssessment, run_refinement_stages
from dvbfixer.model.diffusion.runner_refinement import rewrite_generated_coordinates

CASES = (
    ("prospective-8ek6-H-cdr1-gap5", 7, "8EK6", "89f51d77dda781811e96460563b17e606fcf76133a114a791315cc2140368d3f"),
    ("prospective-8ek6-H-cdr1-gap5", 19, "8EK6", "6ec903751aa847688f17dd70220edad881f360ed2f36b413a676df7da4004f7c"),
    ("prospective-9ww8-gap10", 7, "9WW8", "2e325aff771b01539f41af3d67967dbaafee7417af814c8ade666a61fc5b798a"),
    ("prospective-9ww8-gap10", 19, "9WW8", "bfc629878f276e5b23eb8a6c31d1fb8e1e0696cf98d4656929afd6383c0ae612"),
    ("prospective-9ww8-gap5", 7, "9WW8", "654738da58fd021d55c593d7a246d5f4a4101593eaf0fd71a172e1a634d150e4"),
    ("prospective-9ww8-gap5", 19, "9WW8", "51711f92784987c36a3cd017c686e081764887c628a20d3a681900eb73883fc6"),
    ("prospective-9z0u-H-cdr3-gap10", 7, "9Z0U", "59bea4a2b0d619cd763ad1d2bc353407468b31ecb381a6001d12226986ce51e3"),
    ("prospective-9z0u-H-cdr3-gap10", 19, "9Z0U", "62e44be7c89062a524795eacba420bab1a0c6aaf50aaff97c94b9e527bb80855"),
)


def sha(text: str) -> str:
    return hashlib.sha256(text.encode("ascii")).hexdigest()


def write_json(path, value):
    def finite(item):
        if isinstance(item, float) and not math.isfinite(item):
            return None
        if isinstance(item, dict):
            return {key: finite(value) for key, value in item.items()}
        if isinstance(item, (tuple, list)):
            return [finite(value) for value in item]
        return item
    path.write_text(json.dumps(finite(value), indent=2, ensure_ascii=False, allow_nan=False) + "\n")


def verify_experimental_origin(normalized, original, *, entry_metadata=None, expected_entry=None):
    """Exact observed atom/coordinate/occupancy/B provenance, independent of IDs."""
    headers = [line for line in original.splitlines() if line.startswith("EXPDTA")]
    method = [line[6:].strip() for line in headers]
    method_source = "tracked-PDB-header"
    if not method and entry_metadata is not None:
        if entry_metadata.get("entry", {}).get("id") != expected_entry:
            raise ValueError("experimental method metadata has wrong entry identity")
        method = [item["method"] for item in entry_metadata.get("exptl", [])]
        method_source = "RCSB-entry-metadata"
    if not any("X-RAY DIFFRACTION" in line or "NEUTRON DIFFRACTION" in line for line in method):
        raise ValueError("upstream experimental method unavailable")
    def records(text):
        return Counter((line[12:16].strip(), line[17:20].strip(),
                        *(float(line[i:j]) for i, j in ((30, 38), (38, 46), (46, 54), (54, 60), (60, 66))))
                       for line in text.splitlines() if line.startswith(("ATOM  ", "HETATM")))
    selected = records(normalized)
    if not selected or selected - records(original):
        raise ValueError("normalized observed coordinates/B/occupancy do not match upstream structure")
    return {"experimental_methods": method, "method_source": method_source,
            "upstream_sha256": sha(original), "normalized_sha256": sha(normalized),
            "observed_atom_count": sum(selected.values()), "matching_ignores_only_serial_chain_residue_renumbering": True}


def load_case(root, specification, metadata_directory):
    name, seed, pdb_id, expected_raw = specification
    workspace = root / name / f"seed{seed}"
    request = DiffusionRequest.from_json((workspace / "request.json").read_text())
    if request.seeds != (seed,) or request.candidate_count != 1 or len(request.gaps) != 1:
        raise ValueError("unexpected frozen request")
    source = (workspace / request.normalized_pdb.path).read_bytes().decode("ascii")
    if sha(source) != request.normalized_pdb.sha256:
        raise ValueError("source digest mismatch")
    runner = RunnerResult.from_json((workspace / "raw-result.json").read_text())
    if len(runner.candidates) != 1:
        raise ValueError("unexpected raw candidate count")
    artifact = runner.candidates[0].coordinate_artifact
    raw_path = (workspace / artifact.path).resolve(strict=True)
    if not raw_path.is_relative_to(workspace.resolve()):
        raise ValueError("raw path escapes workspace")
    raw = raw_path.read_bytes().decode("ascii")
    if sha(raw) != artifact.sha256 or sha(raw) != expected_raw:
        raise ValueError("frozen raw digest mismatch")
    original_path = Path("tests/fixtures/cc89_apple_handoff") / f"{pdb_id}.pdb"
    metadata_bytes = (metadata_directory / f"{pdb_id}.json").read_bytes()
    provenance = verify_experimental_origin(source, original_path.read_bytes().decode("ascii"),
        entry_metadata=json.loads(metadata_bytes), expected_entry=pdb_id)
    provenance["entry_metadata_sha256"] = hashlib.sha256(metadata_bytes).hexdigest()
    return request, source, raw, provenance


def run_case(spec, loaded, directory, limits):
    request, source_text, archived_raw, provenance = loaded
    raw = reconcile_pdb_metadata(archived_raw, request).text
    if _observed_lines(raw) != _observed_lines(archived_raw):
        raise ValueError("metadata reconciliation changed atom records")
    directory.mkdir(parents=True)
    (directory / "source.pdb").write_text(source_text)
    (directory / "archived-raw.pdb").write_text(archived_raw)
    (directory / "raw.pdb").write_text(raw)
    (directory / "request.json").write_text(request.to_json())
    write_json(directory / "experimental-origin.json", provenance)
    source = _structure(source_text, "source.pdb")
    source_lines = _observed_lines(source_text)
    observed = source.atoms.keys() - set(request.generated_atoms)
    generated_residues = tuple(dict.fromkeys(anchor_policy.atom_residue(atom) for atom in request.generated_atoms))
    short = all(len(gap.generated_residues) <= 5 for gap in request.gaps)
    native_cache = {}
    native_calls = []
    rows = []
    for b_policy in ("auto", "disabled", "experimental"):
        subdir = directory / b_policy
        subdir.mkdir()
        state = {"restraints": ()}
        assessments = []
        calls = []

        def assess(text, state=state, assessments=assessments):
            candidate = _structure(text, "candidate.pdb")
            assessment, strict = _assessment(request, source, candidate)
            allowed = {item.atom for item in state["restraints"]}
            admission = list(assessment.admission_failures)
            if allowed:
                admission = [gate for gate in admission if gate not in {
                    "fixed-heavy-atom-rmsd", "fixed-heavy-atom-max-displacement"}]
            lines = _observed_lines(text)
            for atom in observed:
                original, current = source_lines[atom], lines.get(atom)
                if current is None or (atom not in allowed and current != original):
                    admission.append("observed-record-outside-mobility-mask-changed")
                elif atom in allowed and (original[:30] != current[:30] or original[54:] != current[54:]):
                    admission.append("observed-anchor-metadata-changed")
            assessments.append({"assessment": asdict(assessment), "historical_strict": asdict(strict)})
            return RepairAssessment(assessment.issues, tuple(dict.fromkeys(admission)))

        def minimize(text, restraints, state=state, calls=calls):
            if text != raw:
                raise ValueError("staged minimizer did not receive raw")
            state["restraints"] = restraints
            key = tuple(sorted((item.atom, item.reference_angstrom, item.k_kj_mol_nm2) for item in restraints))
            reused = key in native_cache
            calls.append({"reused_physical_minimization": reused, "input_sha256": sha(text),
                          "restraints": [asdict(item) for item in restraints]})
            if not reused:
                started = time.perf_counter()
                try:
                    result = refine_generated_region(text, generated_residues=generated_residues,
                        generated_atoms=request.generated_atoms, observed_anchor_restraints=restraints,
                        max_iterations=500 if short else 250, restart_count=8,
                        perturbation_angstrom=0.25 if short else 0.75, platform_name="CPU", random_seed=request.seeds[0])
                    native_cache[key] = rewrite_generated_coordinates(text, result.coordinates_angstrom)
                except (BoundaryRefinementError, ChiralityError, OpenMMException, ValueError) as exc:
                    error_type = NumericalRefinementError if isinstance(exc, NumericalRefinementError) else RuntimeError
                    native_cache[key] = error_type(f"{type(exc).__name__}: {exc}")
                native_calls.append({"stage": "anchors" if restraints else "generated-only",
                    "wall_seconds": time.perf_counter() - started, "restraint_count": len(restraints)})
            if isinstance(native_cache[key], Exception):
                raise native_cache[key]
            return native_cache[key]

        result = run_refinement_stages(raw, source_text, request, assess=assess, minimize=minimize,
            bfactor_source=b_policy, notify=lambda line: print(line, flush=True))
        attempts = []
        for index, attempt in enumerate(result.attempts):
            entry = asdict(attempt)
            text = entry.pop("output_text")
            if text is not None:
                file = subdir / f"stage-{index}-{attempt.stage}.pdb"
                file.write_text(text)
                entry.update(coordinate_file=str(file.relative_to(directory)), coordinate_sha256=sha(text))
            attempts.append(entry)
        outcome = {"case": spec[0], "seed": spec[1], "B_policy": b_policy,
            "archived_raw_sha256": sha(archived_raw), "raw_sha256": sha(raw),
            "selected_stage": result.selected_stage, "selected_anchors": [asdict(r) for r in result.selected_anchors],
            "attempts": attempts, "assessments": assessments, "minimizer_calls": calls,
            "publication_qualified": False, "research_pass": False}
        if result.selected_stage == "generated-and-anchors":
            audit = audit_staged_anchor_repair(request, source_text=source_text, raw_text=raw,
                stage1_text=result.attempts[1].output_text, candidate_text=result.coordinate_text,
                stage1_input_sha256=result.attempts[1].input_sha256,
                stage2_input_sha256=result.attempts[2].input_sha256, limits=limits,
                stage1_failure=result.attempts[1].numerical_failure)
            outcome["audit"] = asdict(audit)
            outcome["research_pass"] = audit.summary.passed
        elif result.selected_stage is not None:
            outcome["research_pass"] = assessments[-1]["historical_strict"]["passed"]
        write_json(subdir / "report.json", outcome)
        rows.append(outcome)
        print(spec[0], spec[1], b_policy, result.selected_stage, "RESEARCH PASS", outcome["research_pass"], flush=True)
    write_json(directory / "native-calls.json", native_calls)
    return {"case": spec[0], "seed": spec[1], "unique_native_minimizations": len(native_calls), "policies": rows}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pilot-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--entry-metadata-dir", required=True, type=Path)
    defaults = AnchorMovementLimits()
    parser.add_argument("--backbone-rmsd-max", default=defaults.backbone_rmsd_angstrom_max, type=float)
    parser.add_argument("--heavy-displacement-max", default=defaults.heavy_displacement_angstrom_max, type=float)
    parser.add_argument("--ca-displacement-max", default=defaults.ca_displacement_angstrom_max, type=float)
    args = parser.parse_args()
    limits = AnchorMovementLimits(args.backbone_rmsd_max, args.heavy_displacement_max, args.ca_displacement_max)
    # Freeze and verify ALL inputs before the first minimization.
    inputs = [load_case(args.pilot_root, spec, args.entry_metadata_dir) for spec in CASES]
    args.output.mkdir(parents=True, exist_ok=False)
    code = args.output / "code"
    code.mkdir()
    paths = [Path(__file__), Path(anchor_policy.__file__), Path(anchor_acceptance.__file__),
             Path(boundary_refinement.__file__), Path(refinement_stages.__file__),
             Path("docs/research/anchor-repair-corrected-cohort-v1.md")]
    for path in paths:
        shutil.copy2(path, code / path.name)
    metadata_output = args.output / "entry-metadata"
    metadata_output.mkdir()
    for pdb_id in sorted({spec[2] for spec in CASES}):
        shutil.copy2(args.entry_metadata_dir / f"{pdb_id}.json", metadata_output / f"{pdb_id}.json")
    write_json(args.output / "manifest.json", {"cases": CASES, "limits": asdict(limits),
        "native_parameters": {"platform": "CPU", "threads": 1, "solvent": False, "restarts": 8,
                              "short_gap_iterations": 500, "short_gap_perturbation_A": 0.25,
                              "long_gap_iterations": 250, "long_gap_perturbation_A": 0.75},
        "code_sha256": {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in paths},
        "publication_qualified": False, "nonfinite_metrics_encoded_as_null": True})
    report = []
    for spec, loaded in zip(CASES, inputs, strict=True):
        directory = args.output / spec[0] / f"seed{spec[1]}"
        try:
            report.append(run_case(spec, loaded, directory, limits))
        except Exception as exc:
            # Preserve failure evidence for every preregistered case, including
            # unexpected errors. Such errors make the cohort command fail.
            directory.mkdir(parents=True, exist_ok=True)
            (directory / "failure.txt").write_text(traceback.format_exc())
            report.append({"case": spec[0], "seed": spec[1], "error": f"{type(exc).__name__}: {exc}"})
        write_json(args.output / "comparison.json", report)
    errors = sum("error" in row for row in report)
    write_json(args.output / "completed.json", {"cases": len(report), "execution_errors": errors, "publication_qualified": False})
    raise SystemExit(bool(errors))


if __name__ == "__main__":
    main()
