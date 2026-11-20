"""One declared fixed/local Cartesian Rosetta control, no retries or source clamp."""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path


def save(path, value):
    with path.open("x") as stream:
        json.dump(value, stream, indent=2, allow_nan=False)


def checked_budget(payload):
    resume = bool(payload.get("resume_from"))
    restraint = checked_junction_restraint(payload)
    budget_control = payload.get("budget_control")
    if budget_control is not None and (budget_control != {"passes": 3, "iterations_per_pass": 3000}
                                       or not resume or restraint or payload["policy"] != "local"):
        raise ValueError("undeclared budget control")
    cap = 3000 if budget_control else 1000 if resume else 500
    if (payload["policy"] not in ("fixed", "local") or payload["max_iterations"] != cap
            or payload["seed"] != 101 or (resume and payload["policy"] != "local") or (restraint and not resume)):
        raise ValueError("undeclared refinement policy")
    return resume, cap


def checked_junction_restraint(payload):
    restraint = payload.get("junction_length_restraint")
    expected = {"first": ["B", "90", "", "C"], "second": ["B", "91", "", "N"],
                "target_A": 1.33, "sd_A": .02, "weight": 1.}
    if restraint is not None and restraint != expected:
        raise ValueError("undeclared junction length restraint")
    return restraint


def resume_coordinates(metadata, before, after, *, expected_parent_cap=500):
    if metadata != before["atoms"] or after["max_iterations"] != expected_parent_cap:
        raise ValueError("resume requires the identical native atom/mobility axis and declared parent iteration stage")
    points = after["coordinates"]
    if (len(points) != len(metadata) or any(len(p) != 3 or not all(math.isfinite(v) for v in p) for p in points)):
        raise ValueError("invalid native resume coordinates")
    return points


def main(root):
    import pyrosetta
    from pyrosetta import rosetta as r

    payload = json.loads((root / "payload.json").read_text())
    resume, cap = checked_budget(payload)
    junction_restraint = checked_junction_restraint(payload)
    pyrosetta.init('-mute all -constant_seed -jran 101 -multithreading:total_threads 1 -ex1 -ex2aro -ignore_unrecognized_res false')
    pose = pyrosetta.pose_from_pdb(str(root / "starting.pdb"))
    source = {tuple(row["identity"]): row["xyz"] for row in payload["source"]}
    generated = set(map(tuple, payload["generated"]))
    local = set(map(tuple, payload["observed_mobile"]))
    if not local <= source.keys() or generated & source.keys() or (payload["policy"] == "fixed" and local):
        raise ValueError("invalid generated/observed masks")
    atom_ids = {}
    metadata = []
    for i in range(1, pose.total_residue()+1):
        residue = pose.residue(i)
        rk = pose.pdb_info().chain(i), str(pose.pdb_info().number(i)), pose.pdb_info().icode(i).strip()
        for j in range(1, residue.natoms()+1):
            a = (*rk, residue.atom_name(j).strip())
            if a in atom_ids:
                raise ValueError("duplicate Rosetta atom identity")
            atom_ids[a] = r.core.id.AtomID(j, i)
            metadata.append({"identity": list(a), "hydrogen": residue.atom_is_hydrogen(j),
                "virtual": residue.atom_type(j).is_virtual(), "mobile": not residue.atom_type(j).is_virtual()
                and (a in generated or a in local or residue.atom_is_hydrogen(j))})
    if not source.keys() | generated <= atom_ids.keys():
        raise ValueError("Rosetta import lost requested atoms")
    def xyz(a):
        v = pose.xyz(atom_ids[a])
        return [v.x, v.y, v.z]
    if max(math.dist(p, xyz(a)) for a, p in source.items()) > 1e-7:
        raise ValueError("import changed observations")
    score = pyrosetta.create_score_function("ref2015_cart")
    score.set_weight(r.core.scoring.dihedral_constraint, 1.)
    gen_residues = {a[:3] for a in generated}
    amides = []
    for i in range(1, pose.total_residue()):
        if pose.chain(i) != pose.chain(i+1):
            continue
        first = pose.pdb_info().chain(i), str(pose.pdb_info().number(i)), pose.pdb_info().icode(i).strip()
        second = pose.pdb_info().chain(i+1), str(pose.pdb_info().number(i+1)), pose.pdb_info().icode(i+1).strip()
        if first not in gen_residues and second not in gen_residues:
            continue
        axis = [(*first, "CA"), (*first, "C"), (*second, "N"), (*second, "CA")]
        pose.add_constraint(r.core.scoring.constraints.DihedralConstraint(
            *[atom_ids[a] for a in axis], r.core.scoring.func.CircularHarmonicFunc(math.pi, .15)))
        amides.append([list(a) for a in axis])
    if local:
        score.set_weight(r.core.scoring.coordinate_constraint, 1.)
        reference = next(a for a in source if a not in local and a[3] == "CA")
        for a in sorted(local):
            ref = r.numeric.xyzVector_double_t(*source[a])
            limit = .5 if a[3] in {"N", "CA", "C", "O"} else 1.25
            for func in (r.core.scoring.func.HarmonicFunc(0., .2), r.core.scoring.constraints.BoundFunc(0., limit, .05, "local-bound")):
                pose.add_constraint(r.core.scoring.constraints.CoordinateConstraint(atom_ids[a], atom_ids[reference], ref, func))
    if resume:
        parent = Path(payload["resume_from"])
        original = json.loads((parent / "payload.json").read_text())
        if ({k: v for k, v in payload.items() if k not in {"resume_from", "max_iterations", "junction_length_restraint", "budget_control"}} != {
                k: v for k, v in original.items() if k not in {"max_iterations", "resume_from"}}
                or original.get("junction_length_restraint") is not None):
            raise ValueError("resume changed original references, ownership, mask or seed")
        points = resume_coordinates(metadata, json.loads((parent / "native-before.json").read_text()),
                                    json.loads((parent / "native-after.json").read_text()),
                                    expected_parent_cap=1000 if junction_restraint or payload.get("budget_control") else 500)
        for row, point in zip(metadata, points, strict=True):
            pose.set_xyz(atom_ids[tuple(row["identity"])], r.numeric.xyzVector_double_t(*point))
        recovered = [xyz(tuple(row["identity"])) for row in metadata]
        old_score = json.loads((parent / "native.json").read_text())["score"]
        recovered_score = score(pose)
        if recovered != points or not math.isclose(recovered_score, old_score, rel_tol=1e-10, abs_tol=1e-6):
            raise ValueError("native all-atom state/objective did not resume exactly")
        save(root / "resume-audit.json", {"all_coordinates_identical_before_minimization": True,
             "parent_score": old_score, "resumed_score": recovered_score, "packing_executed": False,
             "original_source_references_retained": True, "LBFGS_history_restored": False})
    else:
        task = r.core.pack.task.TaskFactory.create_packer_task(pose)
        task.restrict_to_repacking()
        for i in range(1, pose.total_residue()+1):
            rk = pose.pdb_info().chain(i), str(pose.pdb_info().number(i)), pose.pdb_info().icode(i).strip()
            if rk not in gen_residues:
                task.nonconst_residue_task(i).prevent_repacking()
        r.protocols.minimization_packing.PackRotamersMover(score, task).apply(pose)
        if max(math.dist(p, xyz(a)) for a, p in source.items()) > 1e-7:
            raise ValueError("packing changed an observed heavy atom")
    base_score = score.clone()
    if junction_restraint:
        first, second = tuple(junction_restraint["first"]), tuple(junction_restraint["second"])
        if first not in generated or second not in local or not resume:
            raise ValueError("targeted junction endpoints are not both authorized mobile atoms")
        i, j = atom_ids[first].rsd(), atom_ids[second].rsd()
        if pose.residue(i).connected_residue_at_upper() != j or pose.residue(j).connected_residue_at_lower() != i:
            raise ValueError("targeted junction has no native polymer connection")
        if score.get_weight(r.core.scoring.atom_pair_constraint) != 0:
            raise ValueError("base score unexpectedly already includes atom-pair constraints")
        score.set_weight(r.core.scoring.atom_pair_constraint, junction_restraint["weight"])
        pose.add_constraint(r.core.scoring.constraints.AtomPairConstraint(atom_ids[first], atom_ids[second],
            r.core.scoring.func.HarmonicFunc(junction_restraint["target_A"], junction_restraint["sd_A"])))
        save(root / "junction-length-restraint.json", {**junction_restraint,
             "native_polymer_connection_verified": True, "both_endpoints_mobile": True,
             "before_A": math.dist(xyz(first), xyz(second)), "base_score_before": base_score(pose),
             "augmented_score_before": score(pose), "type": "harmonic restraint, not an exact distance constraint"})
    movemap = r.core.kinematics.MoveMap()
    movemap.set_bb(False)
    movemap.set_chi(False)
    movemap.set_jump(False)
    for row in metadata:
        movemap.set_atom(atom_ids[tuple(row["identity"])], row["mobile"])
    before = [xyz(tuple(row["identity"])) for row in metadata]
    save(root / "native-before.json", {"atoms": metadata, "coordinates": before, "score": score(pose),
        "trans_amides": amides, "source_packing_motion_max_A": None if resume else 0.,
        "packing_executed": not resume, "version": pyrosetta.version(),
        "base_score": base_score(pose), "junction_length_restraint": junction_restraint})
    mover = r.protocols.minimization_packing.MinMover(movemap, score, 'lbfgs_armijo_nonmonotone', 1e-6, True)
    mover.cartesian(True)
    mover.max_iter(cap)
    if payload.get("budget_control"):
        index = {tuple(row["identity"]): i for i, row in enumerate(metadata)}
        previous = before
        for stage in range(1, 4):
            initial_score = score(pose)
            mover.apply(pose)
            current = [xyz(tuple(row["identity"])) for row in metadata]
            distance = math.dist(current[index['B', '90', '', 'C']], current[index['B', '91', '', 'N']])
            save(root / f"pass-{stage}.json", {"coordinates": current, "score_before": initial_score,
                 "score_after": score(pose), "junction_C_N_A": distance, "iteration_cap": cap,
                 "maximum_atom_step_A": max(math.dist(a, b) for a, b in zip(previous, current, strict=True)),
                 "last_move_status": str(mover.get_last_move_status()),
                 "reported_score_diff": mover.score_diff_after_minimization(),
                 "actual_iterations": None, "native_stop_reason": "not exposed by this MinMover binding",
                 "native_convergence_proven": False, "LBFGS_history_reinitialized": True})
            previous = current
    else:
        mover.apply(pose)
    after = [xyz(tuple(row["identity"])) for row in metadata]
    immutable_error = max(math.dist(a, b) for row, a, b in zip(metadata, before, after, strict=True) if not row["mobile"])
    if immutable_error > 1e-7:
        raise ValueError("native minimizer moved outside the declared Cartesian mask")
    save(root / "native-after.json", {"coordinates": after, "immutable_max_A": immutable_error,
        "max_iterations": cap, "last_move_status": str(mover.get_last_move_status()), "convergence_claimed": False})
    save(root / "native.json", {"atoms": [{"identity": list(a), "xyz": xyz(a)} for a in [*source, *sorted(generated)]],
        "score": score(pose), "native_outside_mask_max_A": max(math.dist(p, xyz(a)) for a, p in source.items() if a not in local),
        "source_restore_after_minimization": False, "base_score": base_score(pose),
        "junction_C_N_A": math.dist(xyz(tuple(junction_restraint["first"])), xyz(tuple(junction_restraint["second"]))) if junction_restraint else None})


if __name__ == "__main__":
    output = Path(sys.argv[1]).resolve()
    try:
        main(output)
    except Exception as exc:
        import traceback
        save(output / "native-failure.json", {"error": str(exc), "traceback": traceback.format_exc()})
        raise
