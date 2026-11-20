"""Isolated existing-PyRosetta worker for private saved-diffusion GenKIC controls.

No dvbfixer/scipy imports, installations, model calls or sequence design. Native
coordinates are recorded before source-coordinate reinsertion in the parent.
"""
from __future__ import annotations

import hashlib
import json
import math
import sys
from pathlib import Path


def save(path, value):
    with path.open("x") as stream:
        json.dump(value, stream, indent=2, allow_nan=False)


def carbonyl_oxygen(alpha, carbon, next_n):
    """Only a generated O is idealized; neither peptide anchor is moved."""
    import numpy as np

    alpha, carbon, next_n = map(lambda v: np.array(v, dtype=float), (alpha, carbon, next_n))
    u = alpha - carbon
    if np.linalg.norm(u) < 1e-8:
        raise ValueError("degenerate carbonyl axis")
    u /= np.linalg.norm(u)
    v = next_n - carbon
    v -= np.dot(v, u) * u
    if np.linalg.norm(v) < 1e-8:
        raise ValueError("degenerate peptide plane")
    v /= np.linalg.norm(v)
    angle = math.radians(120.8)
    return (carbon + 1.23 * (math.cos(angle) * u - math.sin(angle) * v)).tolist()


def closure_xml(gap, mobile, width, attempts):
    lo, hi = gap["start"], gap["stop"]
    cut = gap["cut"]
    selected = ",".join(map(str, mobile))
    perturber = ""
    if width:
        atomsets = "".join(
            f'<AddAtoms res1="{i}" atom1="N" res2="{i}" atom2="CA"/>'
            f'<AddAtoms res1="{i}" atom1="CA" res2="{i}" atom2="C"/>' for i in range(lo, hi + 1))
        perturber = f'<AddPerturber effect="perturb_dihedral">{atomsets}<AddValue value="{width}"/></AddPerturber>'
    return f'''<ROSETTASCRIPTS>
<SCOREFXNS><ScoreFunction name="full" weights="ref2015"/></SCOREFXNS>
<RESIDUE_SELECTORS><Index name="mobile" resnums="{selected}"/><Not name="fixed" selector="mobile"/></RESIDUE_SELECTORS>
<TASKOPERATIONS><RestrictToRepacking name="repack_only"/>
<OperateOnResidueSubset name="freeze" selector="fixed"><PreventRepackingRLT/></OperateOnResidueSubset></TASKOPERATIONS>
<MOVERS><PackRotamersMover name="pack" scorefxn="full" task_operations="repack_only,freeze"/>
<GeneralizedKIC name="close" closure_attempts="{attempts}" stop_when_n_solutions_found="1"
selector="lowest_energy_selector" selector_scorefunction="full" pre_selection_mover="pack"
correct_polymer_dependent_atoms="true">
{''.join(f'<AddResidue res_index="{i}"/>' for i in range(lo, hi + 1))}
<SetPivots res1="{lo}" atom1="CA" res2="{(lo + hi)//2}" atom2="CA" res3="{hi}" atom3="CA"/>
<CloseBond res1="{cut}" atom1="C" res2="{cut + 1}" atom2="N" prioratom_res="{cut}" prioratom="CA" followingatom_res="{cut + 1}" followingatom="CA" bondlength="1.33" angle1="116.2" angle2="121.7" torsion="180" randomize_flanking_torsions="false"/>
{perturber}
</GeneralizedKIC></MOVERS><PROTOCOLS><Add mover="close"/></PROTOCOLS></ROSETTASCRIPTS>'''


def main(output):
    import pyrosetta
    from pyrosetta import rosetta

    payload = json.loads((output / "payload.json").read_text())
    operation = payload.get("operation", "closure-and-pack")
    if operation not in {"closure-and-pack", "pack-only", "pack-and-cartmin"}:
        raise ValueError("unknown private Rosetta operation")
    pyrosetta.init(f'-mute all -constant_seed -jran {payload["seed"]} -multithreading:total_threads 1 -ex1 -ex2aro -ignore_unrecognized_res false')
    save(output / "runtime.json", {"version": pyrosetta.version(), "python": sys.version,
        "device": "CPU", "license": "existing operator-provided PyRosetta; separate license, not redistributed",
        "worker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    pose = pyrosetta.pose_from_pdb(str(output / "input-backbone.pdb"))
    if pose.total_residue() != len(payload["residues"]):
        raise ValueError("Rosetta changed the residue inventory")
    index = {}
    for i, r in enumerate(payload["residues"], 1):
        info = pose.pdb_info()
        if (info.chain(i), str(info.number(i)), info.icode(i).strip()) != (r["chain"], r["residue_number"], r["insertion_code"]):
            raise ValueError("Rosetta changed chain/residue/insertion-code identity")
        index[r["chain"], r["residue_number"], r["insertion_code"]] = i

    fixed = {tuple(r["identity"]): r["xyz"] for r in payload["fixed"]}
    generated = list(map(tuple, payload["generated"]))
    def point(p, a):
        xyz = p.residue(index[a[:3]]).xyz(a[3])
        return [xyz.x, xyz.y, xyz.z]
    def audit(p):
        maximum = max(math.dist(xyz, point(p, a)) for a, xyz in fixed.items())
        if maximum > 1e-7:
            raise ValueError(f"native fixed atoms moved: {maximum}")
        return maximum
    audit(pose)
    save(output / "import-audit.json", {"native_fixed_max_A": audit(pose), "all_observations_represented": True})
    score = (pyrosetta.create_score_function("ref2015_cart") if operation == "pack-and-cartmin"
             else pyrosetta.get_fa_scorefxn())
    mobile = sorted({index[a[:3]] for a in generated})
    trials = []
    for trial in range(payload["candidates"]):
        folder = output / f"candidate-{trial:02d}"
        folder.mkdir()
        current = pose.clone()
        closures = []
        width = (0, 10, 20, 40)[trial % 4]
        for gap_index, gap in enumerate(payload["gaps"] if operation == "closure-and-pack" else []):
            # The preselector packs only this gap; subsequent joint packing also
            # resolves generated/generated contacts between different gaps.
            xml = closure_xml(gap, range(gap["start"], gap["stop"] + 1), width, payload["attempts"])
            (folder / f"gap-{gap_index}.xml").write_text(xml)
            objects = rosetta.protocols.rosetta_scripts.XmlObjects.create_from_string(xml)
            mover = objects.get_mover("close")
            mover.apply(current)
            audit(current)
            closures.append({"gap": gap, "successful": mover.last_run_successful(), "perturbation_degrees": width})
            if not mover.last_run_successful():
                break
            # GenKIC closes the N-CA-C chain, but its native dependent-atom
            # correction did not repair O across our internal cut (retained v2
            # A81:O-C-A82:N was 88.9 degrees). Rebuild generated carbonyls only,
            # then perform final joint packing on this actual materialization.
            oxygen_moves = []
            for i in range(gap["start"], gap["stop"] + 1):
                def xyz(residue, name, p=current):
                    v = p.residue(residue).xyz(name)
                    return [v.x, v.y, v.z]
                oxygen = carbonyl_oxygen(xyz(i, "CA"), xyz(i, "C"), xyz(i + 1, "N"))
                oxygen_moves.append(math.dist(xyz(i, "O"), oxygen))
                atom = rosetta.core.id.AtomID(current.residue(i).atom_index("O"), i)
                current.set_xyz(atom, rosetta.numeric.xyzVector_double_t(*oxygen))
            closures[-1]["generated_carbonyl_count"] = len(oxygen_moves)
            closures[-1]["carbonyl_max_movement_A"] = max(oxygen_moves)
            audit(current)
        if operation == "closure-and-pack" and (not all(c["successful"] for c in closures) or len(closures) != len(payload["gaps"])):
            save(folder / "refusal.json", {"reason": "no-closed-solution", "closure": closures})
            trials.append({"candidate": trial, "status": "closure-refused"})
            continue
        task = rosetta.core.pack.task.TaskFactory.create_packer_task(current)
        task.restrict_to_repacking()
        for i in range(1, current.total_residue() + 1):
            if i not in mobile:
                task.nonconst_residue_task(i).prevent_repacking()
        rosetta.protocols.minimization_packing.PackRotamersMover(score, task).apply(current)
        maximum = audit(current)
        if operation == "pack-and-cartmin":
            save(folder / "packed-native.json", {"score": score(current), "native_fixed_max_A": maximum,
                 "atoms": [{"identity": list(a), "xyz": point(current, a)} for a in [*fixed, *generated]]})
            movemap = rosetta.core.kinematics.MoveMap()
            movemap.set_bb(False)
            movemap.set_chi(False)
            movemap.set_jump(False)
            mobility = []
            before = []
            for i in range(1, current.total_residue() + 1):
                residue = current.residue(i)
                for j in range(1, residue.natoms() + 1):
                    # Only generated heavy atoms and disposable, unobserved H.
                    # Rosetta's auxiliary C-terminal OXT stays fixed/excluded.
                    a = (*next(k for k, v in index.items() if v == i), residue.atom_name(j).strip())
                    variable = a in generated or residue.atom_is_hydrogen(j)
                    movemap.set_atom(rosetta.core.id.AtomID(j, i), variable)
                    v = residue.xyz(j)
                    mobility.append({"identity": list(a), "mobile": variable,
                                     "auxiliary": a not in fixed and a not in generated})
                    before.append([v.x, v.y, v.z])
            save(folder / "cartesian-mobility.json", {"atoms": mobility, "before": before})
            mover = rosetta.protocols.minimization_packing.MinMover(movemap, score, "lbfgs_armijo_nonmonotone", 1e-6, True)
            mover.cartesian(True)
            mover.max_iter(payload["max_iterations"])
            mover.apply(current)
            maximum = audit(current)
            after = []
            for i in range(1, current.total_residue() + 1):
                for j in range(1, current.residue(i).natoms() + 1):
                    v = current.residue(i).xyz(j)
                    after.append([v.x, v.y, v.z])
            if any(math.dist(a, b) > 1e-7 for row, a, b in zip(mobility, before, after, strict=True) if not row["mobile"]):
                raise ValueError("cartesian minimization moved an atom outside its explicit mobile mask")
            save(folder / "cartesian-returned-state.json", {"after": after, "max_iterations": payload["max_iterations"],
                 "native_fixed_max_A": maximum, "convergence_claimed": False,
                 "last_move_status": str(mover.get_last_move_status())})
        save(folder / "native.json", {"score": score(current), "native_fixed_max_A": maximum, "closure": closures, "operation": operation,
            "atoms": [{"identity": list(a), "xyz": point(current, a)} for a in [*fixed, *generated]]})
        trials.append({"candidate": trial, "status": "native-recorded"})
    save(output / "worker-completed.json", {"trials": trials, "new_diffusion_calls": 0})


if __name__ == "__main__":
    main(Path(sys.argv[1]).resolve())
