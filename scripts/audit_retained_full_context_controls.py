#!/usr/bin/env python3
"""Independent saved-file audit of bounded full-context controls, no modeling."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import subprocess
import sys
from collections import Counter
from dataclasses import asdict
from pathlib import Path

from dvbfixer.ffutils.dat import DatRecord
from dvbfixer.model.diffusion.contract import DiffusionRequest


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def atoms(path):
    """Independent full-context parser; retain fractional HETATM metadata exactly."""
    result = {}
    for line in path.read_text().splitlines(keepends=True):
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        key = line[21], line[22:26].strip(), line[26].strip(), line[12:16].strip()
        occupancy = float(line[54:60])
        if (key in result or line[16].strip() or not 0 < occupancy <= 1
                or (line.startswith("ATOM  ") and occupancy != 1) or line[76:78].strip() in {"H", "D"}):
            raise ValueError("ambiguous or unsupported observation metadata")
        xyz = tuple(float(line[i:i+8]) for i in (30, 38, 46))
        if not all(math.isfinite(value) for value in xyz):
            raise ValueError("nonfinite native coordinates")
        result[key] = {"xyz": xyz, "line": line, "resname": line[17:20].strip(), "element": line[76:78].strip()}
    return result


def ordered_groups(rows):
    result = []
    seen = set()
    for key, row in rows.items():
        if not result or result[-1][0] != key[:3]:
            if key[:3] in seen:
                raise ValueError("noncontiguous saved residue")
            seen.add(key[:3])
            result.append((key[:3], []))
        result[-1][1].append((key, row))
    return result


def native_points(request, source_context, raw):
    target = request.target_sequences[0]
    placement = request.sequence_placements[0]
    mapping = dict(zip(placement.observed_target_indices, placement.observed_residues, strict=True))
    for gap in request.gaps:
        for i, residue in zip(range(gap.target_interval.start, gap.target_interval.stop), gap.generated_residues, strict=True):
            if i in mapping:
                raise ValueError("frozen axis overlap")
            mapping[i] = residue
    if set(mapping) != set(range(len(target.sequence))):
        raise ValueError("incomplete frozen residue axis")
    groups = ordered_groups(raw)
    context = [rows for _, rows in ordered_groups(source_context) if rows[0][1]["line"].startswith("HETATM")]
    if len(groups) != len(mapping) + len(context):
        raise ValueError("native complete context was not retained")
    result = {}
    for i, (_, rows) in enumerate(groups[:len(mapping)]):
        residue = mapping[i]
        for key, row in rows:
            identity = (residue.chain, residue.residue_number, residue.insertion_code, key[3])
            if identity in result:
                raise ValueError("native duplicate protein identity")
            result[identity] = row["xyz"]
    for (_, rows), reference in zip(groups[len(mapping):], context, strict=True):
        by_name = {key[3]: row["xyz"] for key, row in rows}
        if len(reference) == 1 and reference[0][1]["resname"] == "HOH" and set(by_name) == {"OH2"}:
            by_name = {"O": by_name["OH2"]}
        if set(by_name) != {key[3] for key, _ in reference}:
            raise ValueError("native context heavy atom inventory changed")
        for key, _ in reference:
            if key in result:
                raise ValueError("native context atom collides with protein")
            result[key] = by_name[key[3]]
    return result


def check_stage(request, source_context, native, full, dat_path):
    points = native_points(request, source_context, native)
    generated = {(a.chain, a.residue_number, a.insertion_code, a.atom_name) for a in request.generated_atoms}
    fixed = set(source_context)
    extra = set(points) - fixed - generated
    last = request.gaps[-1].generated_residues[-1] if request.gaps[-1].right_anchor is None else request.sequence_placements[0].observed_residues[-1]
    terminal_oxt = (last.chain, last.residue_number, last.insertion_code, "OXT")
    if extra - {terminal_oxt} or not fixed | generated <= points.keys():
        raise ValueError("native atom ownership changed")
    distances = {key: math.dist(row["xyz"], points[key]) for key, row in source_context.items() if row["xyz"] != points[key]}
    if set(full) != fixed | generated or any(full[key]["line"] != row["line"] for key, row in source_context.items()):
        raise ValueError("full export changed source/context records")
    if any(full[key]["xyz"] != points[key] for key in generated):
        raise ValueError("materialized generated atoms differ from native output")
    record = DatRecord.load(dat_path)
    owned = {(a["chain"], a["resid"], a["icode"], a["atom"]) for a in record.added_atoms}
    if owned != generated or len(record.added_atoms) != len(generated) or owned & fixed:
        raise ValueError("canonical sidecar does not identify exactly generated atoms")
    return {"fixed_atoms": len(fixed), "generated_atoms": len(generated), "changed_native_observations": len(distances),
            "max_native_observation_motion_angstrom": max(distances.values(), default=0.),
            "all_context_observation_records_byte_identical": True, "generated_coordinates_match_actual_native_output": True,
            "canonical_dat_ownership_exact": True, "native_engine_only_terminal_OXT_excluded": bool(extra)}


def main(run, output):
    plan_path = run / "plan.json"
    plan = json.loads(plan_path.read_text())
    snapshots_path = run / "executed-code/index.json"
    snapshots = json.loads(snapshots_path.read_text())
    if snapshots["plan_sha256"] != sha(plan_path):
        raise ValueError("executed-code snapshots bind another plan")
    aliases = {row["original_path"]: row for row in snapshots["snapshots"]}
    verified = []
    for name, expected in plan["source_code_and_input_sha256"].items():
        path = Path(name)
        if sha(path) != expected:
            alias = aliases.get(name)
            if path.suffix != ".py" or alias is None or alias["sha256"] != expected or sha(Path(alias["snapshot_path"])) != expected:
                raise ValueError(f"declared input changed: {name}")
            path = Path(alias["snapshot_path"])
        verified.append(path)
    output.mkdir(parents=True, exist_ok=False)
    reports = []
    for job in plan["jobs"]:
        case = Path(job["case"])
        work = run / "execution" / f"{job['case_id']}-seed-{job['seed']}"
        execution_path = work / "execution.json"
        if not execution_path.is_file():
            reports.append({"entry_id": job["case_id"], "status": "not_executed"})
            continue
        execution = json.loads(execution_path.read_text())
        if not execution["completed"] or execution["returncode"] != 0:
            reports.append({"entry_id": job["case_id"], "status": "retained_execution_failure", "execution": execution})
            verified.append(execution_path)
            continue
        request = DiffusionRequest.from_json((case / "request.json").read_text())
        source = atoms(case / "prepared-full-context.pdb")
        declared = json.loads((work / "native/stages.json").read_text())
        calls = json.loads((work / "native/native-interface-calls-final.json").read_text())
        expected = {(a.chain, a.residue_number, a.insertion_code, a.atom_name) for a in request.generated_atoms}
        for call in calls:
            if call["operation"] == "selection":
                selected = {(a["chain"], a["residue_number"], a["insertion_code"], a["atom_name"]) for a in call["atoms"]}
                if selected != expected or len(call["atoms"]) != len(expected):
                    raise ValueError("native optimizer selection included observations or missed generated atoms")
        stages = []
        for row in declared["stages"]:
            label = row["stage"]
            raw_path = work / "native" / ("native.pdb" if label == "raw" else "generated-only-CG-native.pdb")
            full_path = work / "native" / f"{label}-full-context.pdb"
            sidecar = work / "native" / f"{label}.dat"
            audit = check_stage(request, source, atoms(raw_path), atoms(full_path), sidecar)
            if audit["changed_native_observations"] != len(row["native_fixed"]["changed_atoms"]):
                raise ValueError("independent native observation audit disagrees")
            audit.update({"stage": label, "declared_accepted": row["accepted"], "CLI": {}})
            for source_label, path in (("original", case / "original-full-context.pdb"),
                                        ("prepared", case / "prepared-full-context.pdb"), ("candidate", full_path)):
                report = output / f"{job['case_id']}-{label}-{source_label}.json"
                command = [sys.executable, "-m", "dvbfixer.cli", "diagnose", str(path), "--format", "json",
                           "--output", str(report), "--include-water", "--clash-cutoff", "0.6,0.9"]
                env = {**os.environ, "PATH": str(Path(sys.executable).parent) + os.pathsep + os.environ.get("PATH", ""),
                       "OPENMM_CPU_THREADS": "1", "OMP_NUM_THREADS": "1"}
                with (report.with_suffix(".log")).open("xb") as log:
                    result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, timeout=180, env=env)
                if result.returncode not in (0, 1) or not report.is_file():
                    raise RuntimeError("full diagnostics failed, not a pass")
                findings = json.loads(report.read_text())["findings"]
                audit["CLI"][source_label] = {"returncode": result.returncode, "by_severity": dict(Counter(f["severity"] for f in findings)),
                                             "report": str(report), "sha256": sha(report)}
                verified.append(report)
            stages.append(audit)
            verified += [raw_path, full_path, sidecar]
        reports.append({"entry_id": job["case_id"], "seed": job["seed"], "status": "native_stages_independently_read_back",
                        "stages": stages, "native_iupac_calls_skipped": sum(call["operation"] == "to_iupac" for call in calls),
                        "declared_accepted": declared["accepted"], "public_output_published": False})
        verified += [execution_path, work / "native/stages.json", work / "native/native-interface-calls-final.json"]
    verified += [plan_path, snapshots_path, Path(__file__)]
    with (output / "receipt.json").open("x") as stream:
        json.dump({"rows": reports, "new_model_or_optimizer_calls": 0,
                   "bindings": {str(p): sha(p) for p in verified}}, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write("\n")
    print(json.dumps({"rows": [{"entry": row["entry_id"], "status": row["status"], "accepted": row.get("declared_accepted", False)} for row in reports]}))


def audit_saved_analysis(analysis, output):
    from openmm.app import PDBFile

    from dvbfixer.diagnose.chemistry import run_all
    from dvbfixer.diagnose.report import Severity
    from dvbfixer.diagnose.steric import clashes_python
    from dvbfixer.ffutils.geometry import ChiralityError, assert_all_l
    from dvbfixer.model.diffusion.quality import check_ramachandran, check_sidechain_chi12
    from dvbfixer.model.diffusion.staged_runtime import assessment_adapter

    plan_path, completed_path = analysis / "plan.json", analysis / "completed.json"
    plan = json.loads(plan_path.read_text())
    aliases = {}
    if "previous" in plan:
        previous = Path(plan["previous"])
        index = json.loads((previous / "executed-code/index.json").read_text())
        if sha(previous / "plan.json") != index["plan_sha256"]:
            raise ValueError("historical code snapshot binds another declaration")
        aliases = {row["original_path"]: Path(row["snapshot_path"]) for row in index["snapshots"]}
    retained_code = []
    for name, expected in plan["bindings"].items():
        if sha(Path(name)) != expected:
            path = aliases.get(name)
            if Path(name).suffix != ".py" or path is None or sha(path) != expected:
                raise ValueError(f"saved-analysis binding changed: {name}")
            retained_code.append(path)
    output.mkdir(parents=True, exist_ok=False)
    reports, paths, cli_cache = [], [plan_path, completed_path, Path(__file__), *retained_code], {}

    def geometry(path):
        pdb = PDBFile(str(path))
        rows = run_all(pdb.topology, pdb.positions)
        rows += check_ramachandran(pdb.topology, pdb.positions)
        rows += check_sidechain_chi12(pdb.topology, pdb.positions)
        rows += clashes_python(pdb.topology, pdb.positions, include_water=True, clash_warn_a=.6, clash_error_a=.9)
        return [asdict(row) for row in rows if row.severity != Severity.INFO]

    def diagnose(label, path):
        if label in cli_cache:
            return cli_cache[label]
        report = output / f"diagnose-{label}.json"
        command = [sys.executable, "-m", "dvbfixer.cli", "diagnose", str(path), "--format", "json", "--output", str(report),
                   "--include-water", "--clash-cutoff", "0.6,0.9"]
        env = {**os.environ, "PATH": str(Path(sys.executable).parent) + os.pathsep + os.environ.get("PATH", ""),
               "OPENMM_CPU_THREADS": "1", "OMP_NUM_THREADS": "1"}
        with report.with_suffix(".log").open("xb") as log:
            result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, timeout=180, env=env)
        if result.returncode not in (0, 1) or not report.is_file():
            raise RuntimeError("full diagnostic execution failed")
        findings = json.loads(report.read_text())["findings"]
        cli_cache[label] = {"returncode": result.returncode, "report": str(report), "sha256": sha(report),
                            "by_severity": dict(Counter(row["severity"] for row in findings)),
                            "terminal_preparation_findings": [row for row in findings if row["category"] == "missing_terminals"]}
        paths.append(report)
        return cli_cache[label]

    completed = json.loads(completed_path.read_text())
    for row in completed["rows"]:
        case, directory, native = Path(row["case"]), Path(row["directory"]), Path(row["native"])
        request = DiffusionRequest.from_json((case / "request.json").read_text())
        source_path = case / "prepared-full-context.pdb"
        full_path = directory / "full-context.pdb"
        source, full, raw = atoms(source_path), atoms(full_path), atoms(native)
        report = check_stage(request, source, raw, full, directory / "generated.dat")
        stream = iter(full_path.read_text().splitlines(keepends=True))
        if not all(any(candidate == original for candidate in stream) for original in source_path.read_text().splitlines(keepends=True)):
            raise ValueError("original full-context records are not an ordered exact subsequence")
        projection_path = directory / "projection.pdb"
        pdb = PDBFile(str(projection_path))
        all_l = True
        try:
            assert_all_l(pdb.topology, pdb.positions)
        except ChiralityError:
            all_l = False
        assessment = assessment_adapter(request, (case / "source.pdb").read_text())(projection_path.read_text())
        def signature(finding):
            return json.dumps(finding, sort_keys=True, ensure_ascii=False)
        baseline = Counter(signature(finding) for finding in geometry(source_path))
        new_findings = []
        for finding in geometry(full_path):
            key = signature(finding)
            if baseline[key]:
                baseline[key] -= 1
            else:
                new_findings.append(finding)
        raw_job = next(job for job in plan["jobs"] if job["entry_id"] == row["entry_id"] and job["stage"] == "raw")
        original_points = native_points(request, source, atoms(Path(raw_job["native"])))
        current_points = native_points(request, source, raw)
        generated = {(a.chain, a.residue_number, a.insertion_code, a.atom_name) for a in request.generated_atoms}
        max_motion = max(math.dist(original_points[key], current_points[key]) for key in generated)
        accepted = (not report["changed_native_observations"] and all_l and assessment.geometry_passed
                    and not assessment.admission_failures and not new_findings and max_motion <= 3.)
        if accepted != row["accepted"] or len(new_findings) != len(row["full_diagnostics"]["new_or_changed_findings"]):
            raise ValueError("independent full-geometry acceptance differs from saved analysis")
        report.update({"entry_id": row["entry_id"], "stage": row["stage"], "accepted": accepted, "global_all_l": all_l,
             "new_geometry_findings": new_findings, "max_generated_motion_from_original_raw_angstrom": max_motion,
             "original_full_context_record_subsequence_exact": True, "shared_canonical_geometry_detectors": True,
             "CLI": {"original": diagnose(f"{row['entry_id']}-original", case / "original-full-context.pdb"),
                     "prepared": diagnose(f"{row['entry_id']}-prepared", source_path),
                     "candidate": diagnose(f"{row['entry_id']}-{row['stage']}", full_path)},
             "public_output_published": False, "MD_readiness_or_coordinate_accuracy_established": False})
        reports.append(report)
        paths += [native, source_path, case / "request.json", projection_path, full_path, directory / "generated.dat", directory / "result.json"]
    with (output / "receipt.json").open("x") as stream:
        json.dump({"rows": reports, "model_calls": 0, "optimizer_calls": 0, "bindings": {str(p): sha(p) for p in paths}},
                  stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write("\n")
    print(json.dumps({"rows": [{"entry": row["entry_id"], "stage": row["stage"], "accepted": row["accepted"],
                               "native_fixed": row["fixed_atoms"], "new_geometry": len(row["new_geometry_findings"])} for row in reports]}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path)
    parser.add_argument("--saved-analysis", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if (args.run is None) == (args.saved_analysis is None):
        parser.error("supply exactly one of --run or --saved-analysis")
    if args.saved_analysis is not None:
        audit_saved_analysis(args.saved_analysis.resolve(), args.output.resolve())
    else:
        main(args.run.resolve(), args.output.resolve())
