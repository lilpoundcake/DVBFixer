"""Independent native annealing/control readback; never import MODELLER or optimize."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import runpy
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[1]
HELPERS = PROJECT / "tmp/mac-real-gap-continuation-2026-11-12/1KTZ-native-terminal-restoration-attempt01/verify_readback.py"
HANDED_HELPER = PROJECT / "tmp/mac-real-gap-continuation-2026-11-14/1KTZ-inward-wall-refinement-attempt01/verify_readback.py"


def sha(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def main(root):
    helpers = runpy.run_path(str(HELPERS))
    key, inventory, water_contacts = [helpers[n] for n in ("key", "inventory", "water_contacts")]
    handedness = runpy.run_path(str(HANDED_HELPER))["handedness"]
    plan = json.loads((root / "plan.json").read_text())
    for name, expected in plan["source_code_and_input_sha256"].items():
        if sha(Path(name)) != expected:
            raise ValueError("native selected procedure input/code changed")
    if len(plan["jobs"]) != 1:
        raise ValueError("native selected procedure is not a single declared job")
    job = plan["jobs"][0]
    work = root / "execution" / f"{job['case_id']}-seed-{job['seed']}"
    result = json.loads((work / "result.json").read_text())
    execution = json.loads((work / "execution.json").read_text())
    if execution["returncode"] or execution["guard_stop_reason"] or not execution["completed"]:
        raise ValueError("native selected execution failed its unchanged resource guard")
    declaration = json.loads((root / "declaration.json").read_text())
    receipt = {"worker_result": result, "execution": execution, "scientific_qualification": False,
               "native_selected_MD_interface_and_fixation_passed": False, "optimizer_reexecuted": False,
               "scientific_restoration_accepted": False}
    files = [Path(__file__), HELPERS, HANDED_HELPER, root / "declaration.json", root / "plan.json",
             work / "result.json", work / "execution.json"]
    policy_path = work / "native-policy.json"
    rows = []
    if policy_path.exists():
        policy = json.loads(policy_path.read_text())
        active = {key(a) for a in policy["active_atoms"]}
        initial = {key(row["atom"]): tuple(row["xyz_angstrom"]) for row in policy["initial_coordinates"]}
        approved = json.loads((work / "approved-topology-policy.json").read_text())
        targets = {tuple(sorted(key(a) for a in pair)) for pair in approved["target_pairs"]}
        actual = [tuple(sorted(key(a) for a in pair)) for pair in policy["actual_SS_pairs"]]
        source = inventory(work / "source.pdb")
        context = inventory(work / "full-context.pdb")
        waters = {k: row for k, row in context.items() if row[0] == "HOH"}
        parent = inventory(work / "parent-torsion-PROPOSAL.pdb")
        original = inventory(work / "raw-candidate-REJECTED.pdb")
        expected = inventory(work / "native-input-fixed-waters.pdb")
        if (len(active) != 202 or len(waters) != 326 or len(actual) != 21 or set(actual) != targets
                or not expected.keys() <= initial.keys() or len(parent) - len(active) != 3494):
            raise ValueError("native selected graph/mask/full context differs")
        extra = initial.keys() - expected.keys()
        if not extra <= {(chain, "136", "", "OXT") for chain in ("B", "D")}:
            raise ValueError("unexpected temporary native extra atom")
        native_input_residual = max(math.dist(initial[k], row[2]) for k, row in expected.items())
        if native_input_residual > 1e-4:
            raise ValueError("native selected input transport drifted before optimization")
        for number, stage in enumerate(declaration["schedule"], 1):
            coordinate_path = work / f"native-stage-{number:02d}-coordinates.json"
            if not coordinate_path.exists():
                break
            path = work / f"native-stage-{number:02d}-{stage['name']}.pdb"
            native = inventory(path)
            points = {key(row["atom"]): tuple(row["xyz_angstrom"])
                      for row in json.loads(coordinate_path.read_text())["coordinates"]}
            if points.keys() != initial.keys() or native.keys() != initial.keys():
                raise ValueError("native selected stage lost/added identity")
            for k, xyz in points.items():
                if not all(math.isfinite(v) for v in xyz) or native[k][2] != tuple(float(f"{v:8.3f}") for v in xyz):
                    raise ValueError("native selected stage is nonfinite or its PDB differs from raw array")
            fixed = max(math.dist(initial[k], points[k]) for k in points if k not in active)
            raw_max = max(math.dist(original[k][2], points[k]) for k in active)
            active_max = max(math.dist(initial[k], points[k]) for k in active)
            checked, bad = handedness(points)
            graph = []
            for line in path.read_text().splitlines():
                if line.startswith("SSBOND"):
                    graph.append(tuple(sorted(((line[15], line[17:21].strip(), line[21].strip(), "SG"),
                                               (line[29], line[31:35].strip(), line[35].strip(), "SG")))))
            if len(graph) != 21 or set(graph) != targets or fixed > 1e-4 or raw_max > 60.:
                raise ValueError("native selected stage broke graph/fixed/60Ascope")
            rows.append({"ordinal": number, "method": stage["method"], "native_fixed_max_motion_angstrom": fixed,
                         "active_max_motion_from_parent_angstrom": active_max,
                         "active_max_motion_from_original_raw_angstrom": raw_max,
                         "independent_chiral_centers": checked, "D_centers": bad})
            if result["native_procedure_completed"] and bad:
                raise ValueError("native selected procedure passed despite a Dcentre")
        receipt.update(native_stage_readbacks=rows, native_initial_max_residual_angstrom=native_input_residual,
                       actual_native_SS_pairs=21, active_atoms=202, fixed_protein_atoms=3494,
                       fixed_water_oxygens=326, temporary_fixed_extras=len(extra))
        if result["native_procedure_completed"]:
            if len(rows) != len(declaration["schedule"]) or not rows or rows[-1]["active_max_motion_from_parent_angstrom"] <= 1e-6:
                raise ValueError("native selected stages incomplete or no selected motion")
            receipt["native_selected_MD_interface_and_fixation_passed"] = True
        final_path = work / "annealed-terminal-UNQUALIFIED.pdb"
        if declaration["mode"] == "api-control" and final_path.exists():
            raise ValueError("API control unexpectedly materialized a scientific derivative")
        if final_path.exists():
            if not receipt["native_selected_MD_interface_and_fixation_passed"]:
                raise ValueError("native materialization without full native fixation gate")
            final = inventory(final_path)
            if final.keys() != parent.keys() or any(final[k] != row for k, row in source.items()):
                raise ValueError("annealed final inventory or source observations changed")
            for old, new in zip((work / "parent-torsion-PROPOSAL.pdb").read_bytes().splitlines(keepends=True),
                                final_path.read_bytes().splitlines(keepends=True), strict=True):
                k = (old[21:22].decode(), old[22:26].strip().decode(), old[26:27].strip().decode(), old[12:16].strip().decode())
                if old.startswith((b"ATOM  ", b"HETATM")) and k in active:
                    if old[:30] != new[:30] or old[54:] != new[54:]:
                        raise ValueError("annealed terminal metadata changed")
                    if final[k][2] != tuple(float(f"{v:8.3f}") for v in points[k]):
                        raise ValueError("annealed terminal final differs from native coordinates")
                elif old != new:
                    raise ValueError("annealed terminal changed an inactive physical record")
            generated = {key(a) for a in json.loads((work / "request.json").read_text())["generated_atoms"]}
            counts = {"baseline": water_contacts(source, waters),
                      "generated": water_contacts({k: final[k] for k in generated}, waters),
                      "active_terminal": water_contacts({k: final[k] for k in active}, waters)}
            receipt.update(observations_byte_identical=2986, inactive_protein_records_byte_identical=3494,
                           water_contacts=counts, final_sha256=sha(final_path))
        files += [p for p in work.iterdir() if p.is_file() and p.suffix in {".json", ".pdb"} and p not in files]
    receipt["bindings"] = {str(p.resolve()): sha(p) for p in files}
    with (root / "readback-receipt.json").open("x") as stream:
        json.dump(receipt, stream, indent=2)
        stream.write("\n")
    print(json.dumps({k: v for k, v in receipt.items() if k not in {"bindings", "execution", "worker_result"}}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    main(parser.parse_args().root.resolve())
