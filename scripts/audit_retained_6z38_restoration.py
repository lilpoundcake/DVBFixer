#!/usr/bin/env python3
"""Independent readback of retained 6Z38 files; no sampling or coordinate repair."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

from dvbfixer.ffutils.dat import DatRecord
from dvbfixer.model.diffusion.contract import DiffusionRequest


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def atoms(path):
    rows = {}
    for line in path.read_text().splitlines(keepends=True):
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        key = (line[21], line[22:26].strip(), line[26].strip(), line[12:16].strip())
        if key in rows or line[16].strip() or float(line[54:60]) != 1:
            raise ValueError("duplicate or ambiguous retained atom")
        xyz = tuple(float(line[i:i+8]) for i in (30, 38, 46))
        if not all(math.isfinite(x) for x in xyz):
            raise ValueError("nonfinite retained coordinates")
        rows[key] = {"xyz": xyz, "line": line, "resname": line[17:20].strip(), "element": line[76:78].strip()}
    return rows


def motion(source, candidate):
    if not source.keys() <= candidate.keys():
        raise ValueError("missing observations")
    changed = {key: math.dist(row["xyz"], candidate[key]["xyz"]) for key, row in source.items()
               if row["xyz"] != candidate[key]["xyz"]}
    groups = defaultdict(list)
    for key in source:
        groups[key[:3]].append(key)
    permutations = []
    for residue, keys in groups.items():
        if any(key in changed for key in keys) and Counter(source[key]["xyz"] for key in keys) == Counter(candidate[key]["xyz"] for key in keys):
            permutations.append(list(residue))
    return {"changed_named_observations": len(changed), "max_angstrom": max(changed.values(), default=0.),
            "residues_with_exact_coordinate_permutations": permutations,
            "coordinate_permutation_is_not_named_atom_custody": True}


def main(run, output, executed_code_index=None):
    plan = json.loads((run / "plan.json").read_text())
    if plan["seed"] != 102 or plan["count"] != 1:
        raise ValueError("this readback binds only the declared seed102 attempt")
    snapshots = {}
    if executed_code_index is not None:
        for row in json.loads(executed_code_index.read_text())["snapshots"]:
            path = Path(row["snapshot_path"])
            if sha(path) != row["sha256"]:
                raise ValueError("executed-code snapshot changed")
            snapshots[row["original_path"]] = path
    used_snapshots = {}
    for name, expected in plan["source_code_and_input_sha256"].items():
        if sha(Path(name)) != expected:
            snapshot = snapshots.get(name)
            if (Path(name).suffix != ".py" or not Path(name).is_relative_to(Path(__file__).resolve().parent)
                    or snapshot is None or sha(snapshot) != expected):
                raise ValueError(f"declared input changed without exact executed-code snapshot: {name}")
            used_snapshots[name] = str(snapshot)
    case = run / "cases/6Z38"
    work = run / "execution/6Z38-seed-102"
    execution = json.loads((work / "execution.json").read_text())
    if not execution["completed"] or execution["returncode"] != 0 or execution["guard_stop_reason"] is not None:
        raise ValueError("no completed native execution")
    request = DiffusionRequest.from_json((case / "request.json").read_text())
    def key(atom):
        return (atom.chain, atom.residue_number, atom.insertion_code, atom.atom_name)
    fixed, generated = {key(a) for a in request.fixed_atoms}, {key(a) for a in request.generated_atoms}
    source = atoms(case / "source.pdb")
    native = atoms(work / "proposal/mapped-native.pdb")
    raw = atoms(work / "proposal/native.pdb")
    initial = atoms(work / "proposal/target.ini")
    restored = atoms(work / "source-restored.pdb")
    terminal_oxt = ("A", "432", "", "OXT")
    if (set(source) != fixed or set(native) != fixed | generated or set(restored) != fixed | generated
            or set(raw) != fixed | generated | {terminal_oxt} or terminal_oxt in generated
            or any(native[a]["xyz"] != raw[a]["xyz"] for a in native)):
        raise ValueError("native identity/coordinate transport or ownership mismatch")
    if any(source[a]["line"] != restored[a]["line"] for a in source):
        raise ValueError("materialization did not preserve exact source records")
    native_motion = motion(source, native)
    strict = json.loads((work / "strict-readback.json").read_text())
    if (native_motion["changed_named_observations"] != len(strict["native_observation_audit"]["changed_atoms"])
            or not math.isclose(native_motion["max_angstrom"], strict["native_observation_audit"]["max_motion_angstrom"], abs_tol=1e-12)
            or strict["restoration_accepted"]):
        raise ValueError("independent native rejection disagrees with retained report")
    output.mkdir(parents=True, exist_ok=False)
    ownership = DatRecord(description="REJECTED private 6Z38 diagnostic proposal; generated mask only, not an accepted restoration",
                          target_chains={"A": 432})
    for chain, resid, icode, name in sorted(generated):
        row = restored[(chain, resid, icode, name)]
        ownership.added_atoms.append({"chain": chain, "resid": resid, "icode": icode,
                                     "resname": row["resname"], "atom": name, "element": row["element"]})
        summary = ownership.residue_summary.setdefault(f"{chain}/{row['resname']}{resid}{icode}", {"heavy": 0, "hydrogen": 0})
        summary["heavy"] += 1
    ownership.save(output / "rejected-source-restored.dat")
    loaded = DatRecord.load(output / "rejected-source-restored.dat")
    if {(a["chain"], a["resid"], a["icode"], a["atom"]) for a in loaded.added_atoms} != generated:
        raise ValueError("canonical sidecar changed ownership")
    parent = Path(next(name for name in plan["source_code_and_input_sha256"]
                       if name.endswith("6Z38-seed-0/full-assembly-context.pdb"))).parent
    diagnoses = {}
    env = {**os.environ, "PATH": str(Path(sys.executable).parent) + os.pathsep + os.environ.get("PATH", ""),
           "OPENMM_CPU_THREADS": "1", "OMP_NUM_THREADS": "1"}
    for label, path in (("original", parent / "full-assembly-context.pdb"), ("prepared_source", case / "source.pdb"),
                        ("restored_rejected", work / "source-restored.pdb")):
        report = output / f"diagnose-{label}.json"
        command = [sys.executable, "-m", "dvbfixer.cli", "diagnose", str(path), "--format", "json", "--output", str(report),
                   "--clash-cutoff", "0.6,0.9"]
        with (output / f"diagnose-{label}.log").open("xb") as log:
            result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, env=env, timeout=180)
        if result.returncode not in (0, 1) or not report.is_file():
            raise RuntimeError("CLI diagnostics failed; no pass inferred")
        diagnostic = json.loads(report.read_text())
        diagnoses[label] = {"returncode": result.returncode, "path": str(report), "sha256": sha(report),
                            "by_severity": dict(Counter(f["severity"] for f in diagnostic["findings"]))}
    payloads = [case / "source.pdb", case / "request.json", run / "plan.json", work / "execution.json",
                work / "strict-readback.json", work / "proposal/native.pdb", work / "proposal/target.ini",
                work / "proposal/mapped-native.pdb", work / "source-restored.pdb", Path(__file__), output / "rejected-source-restored.dat"]
    if executed_code_index is not None:
        payloads += [executed_code_index, *snapshots.values()]
    report = {"entry_id": "6Z38", "seed": 102, "independent_native_motion": native_motion,
              "initial_native_motion": motion(source, initial), "restored_observation_records_exact": True,
              "generated_atoms": len(generated), "fixed_atoms": len(fixed), "engine_only_terminal_OXT_exported": False,
              "canonical_dat_excludes_observations_and_engine_OXT": True, "CLI_diagnostics": diagnoses,
              "declared_implementation_verified_from_exact_snapshots": used_snapshots,
              "new_sampling_or_refinement_calls": 0, "restoration_accepted": False, "public_output_published": False,
              "bindings": {str(path): sha(path) for path in payloads}}
    with (output / "receipt.json").open("x") as stream:
        json.dump(report, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write("\n")
    print(json.dumps({"native_motion": native_motion, "diagnostics": diagnoses, "accepted": False}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--executed-code-index", type=Path)
    args = parser.parse_args()
    main(args.run.resolve(), args.output.resolve(), args.executed_code_index.resolve() if args.executed_code_index else None)
