#!/usr/bin/env python3
"""One standard-stereochemistry CG continuation per saved full-context proposal."""
from __future__ import annotations

import argparse
import json
import math
import os
import shutil
import sys
import traceback
from dataclasses import asdict
from pathlib import Path

from scripts.mac_reconstruction_comparison import digest, save, target_axis
from scripts.retained_modeller_full_context import fixed_audit, native_transport, residue_groups
from scripts.run_mac_modeller_comparison import execute

from dvbfixer.model.diffusion.contract import AtomIdentity, DiffusionRequest

PROJECT = Path(__file__).resolve().parents[1]
OPTIMIZER_SEED = 104


def declare(parent, output, previous=None):
    completed = parent / "completed.json"
    rows = json.loads(completed.read_text())["rows"]
    selected = [next(row for row in rows if row["entry_id"] == entry and row["stage"] == stage)
                for entry, stage in (("1CC6", "generated-only-CG"), ("9M3K", "raw"))]
    if any(row["accepted"] or not row["global_all_l"] or not row["native_fixed"]["exact_at_PDB_precision"] for row in selected):
        raise ValueError("require retained globally-L, fixed-context geometry failures")
    output.mkdir(parents=True, exist_ok=False)
    files = [completed, parent / "plan.json", Path(__file__), PROJECT / "scripts/retained_modeller_full_context.py",
             PROJECT / "scripts/recover_retained_full_context_analysis.py", PROJECT / "scripts/run_mac_modeller_comparison.py",
             PROJECT / "scripts/run_mac_reconstruction_guarded.py", PROJECT / "src/dvbfixer/model/modeller_run.py",
             PROJECT / "src/dvbfixer/ffutils/geometry.py"]
    files += list((PROJECT / "src/dvbfixer/model/diffusion").glob("*.py"))
    if previous is not None:
        work = previous / "execution/1CC6-seed-104"
        receipt = json.loads((work / "execution.json").read_text())
        if (receipt["completed"] or receipt["returncode"] != 1 or (work / "native-policy.json").exists()
                or (work / "stereo-native.pdb").exists() or (previous / "execution/9M3K-seed-104").exists()):
            raise ValueError("only the verified pre-CG initialization failure can continue unchanged optimizer budget")
        files += [previous / "plan.json", work / "execution.json", work / "failure.json", work / "before.pdb"]
        files += list((previous / "executed-code").iterdir())
    jobs = []
    for row in selected:
        case = output / "cases" / row["entry_id"]
        case.mkdir(parents=True)
        raw = next(r for r in rows if r["entry_id"] == row["entry_id"] and r["stage"] == "raw")
        manifest = {"entry_id": row["entry_id"], "original_case": row["case"], "parent_native": row["native"],
                    "original_raw_native": raw["native"], "native_seed": row["seed"], "optimizer_seed": OPTIMIZER_SEED,
                    "parent_analysis": str(Path(row["directory"]) / "result.json")}
        save(case / "manifest.json", manifest)
        files += [case / "manifest.json", Path(row["native"]), Path(raw["native"]), Path(manifest["parent_analysis"])]
        files += [Path(row["case"]) / name for name in ("source.pdb", "prepared-full-context.pdb", "request.json", "alignment.pir")]
        jobs.append({"case_id": row["entry_id"], "case": str(case), "seed": OPTIMIZER_SEED,
                     "command": [sys.executable, "-m", "scripts.run_saved_full_context_stereo", "worker", "--case", str(case),
                                 "--output", str(output / "execution" / f"{row['entry_id']}-seed-{OPTIMIZER_SEED}")]})
    save(output / "plan.json", {"python": sys.executable, "python_sha256": digest(Path(sys.executable)),
         "source_code_and_input_sha256": {str(p): digest(p) for p in files}, "jobs": jobs, "count": 2,
         "native_proposal_seed": 103, "optimizer_seed": OPTIMIZER_SEED, "new_model_make_calls": 0,
         "previous_pre_CG_failure": None if previous is None else str(previous),
         "previous_actual_CG_calls": 0,
         "new_neural_calls": 0, "new_randomization_or_MD_calls": 0, "CG_calls_per_case": 1, "CG_max_iterations": 1000,
         "rationale": "rebuild standard stereo restraints rather than reuse residual loop-scoring restraints; no detector thresholds changed",
         "chirality": "existing native canonical restraints; reject D, never reflect observations",
         "new_peptide_omega": "trans Gaussian 180 degrees / 5-degree stdev for generated peptide C-N bonds only",
         "scope": "all declared generated heavy atoms jointly; every other protein/ligand/cofactor/water atom fixed",
         "movement_reference": "original saved raw proposal, not last derivative; maximum 3 A",
         "solvent_box": False, "new_downloads": 0, "resource_guard": "existing serial1800s/job,RSS8GiB,host-pressure,10GiB evidence",
         "automatic_retries": 0, "public_capability_or_qualification": False})


def worker(case, output):
    from modeller import Alignment, Environ, Model, Selection, features, forms, log, physical
    from modeller.optimizers import ConjugateGradients

    from dvbfixer.ffutils.geometry import assert_all_l
    from dvbfixer.model.diffusion.anchor_acceptance import _structure
    from dvbfixer.model.modeller_run import _add_chirality_restraints

    manifest = json.loads((case / "manifest.json").read_text())
    original_case = Path(manifest["original_case"])
    request = DiffusionRequest.from_json((original_case / "request.json").read_text())
    context = (original_case / "prepared-full-context.pdb").read_text()
    parent = Path(manifest["parent_native"])
    projection, points, _ = native_transport(request, context, parent.read_text())
    if not fixed_audit(context, points)["exact_at_PDB_precision"]:
        raise ValueError("parent native observations are not fixed")
    structure = _structure(projection, "parent.pdb")
    assert_all_l(structure.pdb.topology, structure.pdb.positions)
    generated = set(request.generated_atoms)
    axis = target_axis(request)[request.target_sequences[0].chain]
    previous = Path.cwd()
    try:
        os.chdir(output)
        shutil.copyfile(original_case / "prepared-full-context.pdb", output / "source.pdb")
        shutil.copyfile(original_case / "alignment.pir", output / "alignment.pir")
        log.none()
        env = Environ(rand_seed=-OPTIMIZER_SEED)
        env.io.hetatm = env.io.water = True
        env.io.atom_files_directory = [str(output)]
        env.libs.topology.read(file="$(LIB)/top_heav.lib")
        env.libs.parameters.read(file="$(LIB)/par.lib")
        aln = Alignment(env, file=str(output / "alignment.pir"))
        model = Model(env)
        model.generate_topology(aln["target"])
        model.transfer_xyz(aln)
        parent_groups = residue_groups(parent.read_text())
        if len(model.residues) != len(parent_groups):
            raise ValueError("standard stereo topology changed residue count")
        native_axis, selected = {}, []
        for index, (residue, (_, atoms)) in enumerate(zip(model.residues, parent_groups, strict=True)):
            expected = {atom.identity.atom_name: atom for atom in atoms}
            if {atom.name for atom in residue.atoms} != expected.keys():
                raise ValueError("standard stereo topology changed native atoms")
            for atom in residue.atoms:
                before = expected[atom.name]
                atom.x, atom.y, atom.z = before.xyz
                # transfer_xyz marks absent template coordinates occupancy zero.
                # Load metadata from the retained native proposal too; never
                # replace experimental partial occupancy with a fabricated one.
                atom.occ = float(before.line[54:60])
                atom.biso = float(before.line[60:66])
                if index < len(axis):
                    identity = AtomIdentity(axis[index].chain, axis[index].residue_number, axis[index].insertion_code, atom.name)
                else:
                    identity = before.identity
                native_axis[atom.index] = identity
                if identity in generated:
                    selected.append(atom)
        if {native_axis[atom.index] for atom in selected} != generated or len(selected) != len(generated):
            raise ValueError("stereo optimizer mask is not exactly all generated atoms")
        model.write(file=str(output / "before.pdb"))
        _, loaded, _ = native_transport(request, context, (output / "before.pdb").read_text())
        if loaded.keys() != points.keys() or any(loaded[key].xyz != points[key].xyz for key in loaded):
            raise ValueError("native stereo initialization did not exactly load retained parent")
        in_memory_before = {atom.index: (atom.x, atom.y, atom.z) for atom in model.atoms}
        model.restraints.make(Selection(model), restraint_type="stereo", spline_on_site=False)
        chirality, skipped = _add_chirality_restraints(model)
        trans = 0
        for bond in model.bonds:
            first, last = tuple(bond)
            if first.name == "N" and last.name == "C":
                first, last = last, first
            if (first.name != "C" or last.name != "N" or first.residue.index == last.residue.index
                    or not {native_axis[first.index], native_axis[last.index]} & generated):
                continue
            model.restraints.add(forms.Gaussian(group=physical.omega_dihedral,
                feature=features.Dihedral(first.residue.atoms["CA"], first, last, last.residue.atoms["CA"]),
                mean=math.pi, stdev=math.radians(5.)))
            trans += 1
        model.env.edat.dynamic_sphere = True
        model.env.edat.nonbonded_sel_atoms = 1
        save(output / "native-policy.json", {"active_atoms": [asdict(native_axis[atom.index]) for atom in selected],
             "actual_native_atoms": len(native_axis), "chirality_restraints": chirality, "chirality_skips": skipped,
             "new_peptide_trans_restraints": trans, "parent_exact_at_PDB_precision": True,
             "CG_max_iterations": 1000, "model_make_calls": 0, "randomization_calls": 0, "MD_calls": 0})
        ConjugateGradients(min_atom_shift=.001).optimize(Selection(*selected), max_iterations=1000)
        if any((atom.x, atom.y, atom.z) != in_memory_before[atom.index] for atom in model.atoms if native_axis[atom.index] not in generated):
            raise ValueError("native stereo CG moved an immobile atom")
        model.write(file=str(output / "stereo-native.pdb"))
        _, final, _ = native_transport(request, context, (output / "stereo-native.pdb").read_text())
        fixed = fixed_audit(context, final)
        _, raw, _ = native_transport(request, context, Path(manifest["original_raw_native"]).read_text())
        max_motion = max(math.dist(raw[key].xyz, final[key].xyz) for key in generated)
        save(output / "result.json", {"entry_id": manifest["entry_id"], "native_fixed": fixed,
             "nonactive_in_memory_coordinates_exact": True, "max_generated_motion_from_original_raw_angstrom": max_motion,
             "movement_limit_angstrom": 3., "movement_scope_passed": max_motion <= 3.,
             "actual_CG_calls": 1, "model_make_calls": 0, "public_output_published": False})
    except Exception:
        save(output / "failure.json", {"traceback": traceback.format_exc(), "automatic_retry": False})
        raise
    finally:
        os.chdir(previous)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("declare", "run", "worker"))
    parser.add_argument("--parent", type=Path)
    parser.add_argument("--case", type=Path)
    parser.add_argument("--previous", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.mode == "declare":
        if args.parent is None:
            parser.error("declare requires --parent")
        declare(args.parent.resolve(), args.output.resolve(), args.previous.resolve() if args.previous else None)
    elif args.mode == "run":
        execute(args.output.resolve())
    else:
        if args.case is None:
            parser.error("worker requires --case")
        worker(args.case.resolve(), args.output.resolve())
