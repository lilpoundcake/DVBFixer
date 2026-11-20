"""Read-only independent custody/geometry audit of a bounded terminal comparison."""
from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict
from pathlib import Path

import numpy as np
from scripts.audit_retained_full_context_controls import atoms
from scripts.hybrid_gap_reconstruction import audit_native, key
from scripts.mac_reconstruction_comparison import digest, normalize_modeller, save
from scripts.retained_evidence import verify_bindings
from scripts.terminal_reconstruction_research import prepare_backbone

from dvbfixer.model.diffusion.contract import DiffusionRequest
from dvbfixer.model.diffusion.staged_runtime import assessment_adapter


def audit(root):
    plan = json.loads((root / "plan.json").read_text())
    verify_bindings(plan["source_code_and_input_sha256"], (root / "plan.json",))
    source_path = root / "source.pdb"
    request = DiffusionRequest.from_json((root / "request.json").read_text())
    if digest(source_path) != request.normalized_pdb.sha256:
        raise ValueError("comparison source hash mismatch")
    source = atoms(source_path)
    generated, fixed = set(map(key, request.generated_atoms)), set(map(key, request.fixed_atoms))
    payload_expected = {"fixed": [{"identity": list(a), "xyz": r["xyz"]} for a, r in source.items()],
                        "generated": [list(key(a)) for a in request.generated_atoms]}
    native_modeller = root / "modeller-native"
    mapped, _ = normalize_modeller(request, (native_modeller / "native.pdb").read_text())
    if mapped != (native_modeller / "mapped-native.pdb").read_text():
        raise ValueError("MODELLER identity mapping does not reproduce")
    modeller = atoms(native_modeller / "mapped-native.pdb")
    raw_motion = max(math.dist(modeller[a]["xyz"], row["xyz"]) for a, row in source.items())
    proposals = {"idpforge": atoms(Path(plan["sample"]) / "raw-extraction/native-mapped.pdb"), "modeller": modeller}
    results = {}
    for arm, proposal in proposals.items():
        folder = root / arm
        text, payload, construction = prepare_backbone(request, source, proposal)
        if (text != (folder / "input-backbone.pdb").read_text()
                or json.loads(json.dumps(payload)) != json.loads((folder / "payload.json").read_text())
                or construction != json.loads((folder / "construction.json").read_text())):
            raise ValueError("source-bound reverse construction does not reproduce")
        execution = json.loads((folder / "execution.json").read_text())
        if execution["returncode"] or execution["guard_stop_reason"] is not None:
            results[arm] = {"completed": False}
            continue
        work = folder / "candidate-00"
        packed = json.loads((work / "packed-native.json").read_text())
        final = json.loads((work / "native.json").read_text())
        initial_coords, _ = audit_native(payload_expected, packed)
        final_coords, maximum = audit_native(payload_expected, final)
        mobility = json.loads((work / "cartesian-mobility.json").read_text())
        returned = json.loads((work / "cartesian-returned-state.json").read_text())
        if returned["max_iterations"] != payload["max_iterations"]:
            raise ValueError("minimization budget differs from the declaration")
        before, after = np.array(mobility["before"]), np.array(returned["after"])
        identities = [tuple(r["identity"]) for r in mobility["atoms"]]
        if (len(set(identities)) != len(identities) or not fixed | generated <= set(identities)
                or before.shape != after.shape or before.shape != (len(identities), 3)
                or not np.isfinite(before).all() or not np.isfinite(after).all()):
            raise ValueError("Cartesian native atom axis changed")
        for i, row in enumerate(mobility["atoms"]):
            a = identities[i]
            hydrogen = a[3].lstrip("0123456789").startswith("H")
            if row["mobile"] != (a in generated or hydrogen) or row["auxiliary"] != (a not in fixed | generated):
                raise ValueError("Cartesian mobility differs from generated-heavy/disposable-H policy")
            last = payload["residues"][-1]
            allowed_oxt = (last["chain"], last["residue_number"], last["insertion_code"], "OXT")
            # Rosetta uses an unphysical NV virtual atom for PRO ring closure.
            # It is fixed, not a generated heavy atom or an emitted observation.
            pro_n = source.get((*a[:3], "N"))
            pro_virtual = a[3] == "NV" and pro_n is not None and pro_n["resname"] == "PRO"
            if a not in fixed | generated and not hydrogen and a != allowed_oxt and not pro_virtual:
                raise ValueError("unexpected native auxiliary heavy atom")
            if not row["mobile"] and np.linalg.norm(after[i]-before[i]) > 1e-7:
                raise ValueError("Cartesian optimizer moved an immutable atom")
            if a in fixed | generated and (not np.array_equal(before[i], initial_coords[a]) or not np.array_equal(after[i], final_coords[a])):
                raise ValueError("recorded atom coordinates differ from native optimizer states")
        stages = {}
        for stage, native in (("packed", packed), ("minimized", final)):
            coordinates, _ = audit_native(payload_expected, native)
            path = folder / f"{stage}.pdb"
            actual = atoms(path)
            if (set(actual) != fixed | generated or any(actual[a]["line"] != row["line"] for a, row in source.items())
                    or any(actual[a]["xyz"] != tuple(float(f"{v:.3f}") for v in coordinates[a]) for a in generated)):
                raise ValueError("serialized output changed native generations or source records")
            assessment = assessment_adapter(request, source_path.read_text())(path.read_text())
            from io import StringIO

            from openmm.app import PDBFile

            from dvbfixer.ffutils.geometry import assert_all_l

            structure = PDBFile(StringIO(path.read_text()))
            assert_all_l(structure.topology, structure.positions)
            recorded = json.loads(path.with_suffix(".assessment.json").read_text())
            if asdict(assessment) != recorded:
                # JSON turns tuples into lists; compare its wire representation.
                if json.loads(json.dumps(asdict(assessment))) != recorded:
                    raise ValueError("independent terminal geometry assessment does not reproduce")
            anchor = request.gaps[0].right_anchor
            last = request.gaps[0].generated_residues[-1]
            stages[stage] = {"geometry_passed": assessment.geometry_passed,
                "source_records_exact": True, "global_all_l": True, "candidate_sha256": digest(path),
                "junction_C_N_A": math.dist(actual[last.chain, last.residue_number, last.insertion_code, "C"]["xyz"],
                                            actual[anchor.chain, anchor.residue_number, anchor.insertion_code, "N"]["xyz"])}
        results[arm] = {"completed": True, "native_fixed_max_A": maximum, "stages": stages,
                        "max_iterations": returned["max_iterations"], "convergence_claimed": False}
    return {"plan_sha256": digest(root / "plan.json"), "arms": results, "MODELLER_raw_observed_max_A": raw_motion,
            "source_atoms": len(fixed), "generated_atoms": len(generated), "new_diffusion_calls": 0,
            "public_qualification": False}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    args = parser.parse_args()
    result = audit(args.root.resolve())
    save(args.root / "independent-readback.json", result)
    print(json.dumps(result))
