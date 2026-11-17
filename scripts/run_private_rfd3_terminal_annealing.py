"""Private selected-atom MODELLER API control/terminal annealing, not public CLI.

All scientific parameters/inputs are frozen by an external guarded plan. This
script never reads external MODELLER source, license material or configuration.
"""
from __future__ import annotations

import argparse
import importlib.abc
import json
import math
import socket
import sys
from dataclasses import asdict
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[1]
PARENT = PROJECT / "tmp/mac-real-gap-continuation-2026-11-14/1KTZ-terminal-torsion-study-attempt01"
PRIOR = PARENT / "execution/1KTZ-terminal-torsions-seed-0"


def save(path, value):
    with path.open("x") as stream:
        json.dump(value, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write("\n")


class NoNeuralImports(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname.split(".")[0] in {"torch", "rfd3", "safetensors"}:
            raise RuntimeError("selected native annealing forbids neural imports")


def forbidden(*args, **kwargs):
    raise RuntimeError("selected native annealing forbids network access")


def main(output, declaration_path):
    sys.meta_path.insert(0, NoNeuralImports())
    socket.create_connection = socket.socket.connect = socket.socket.connect_ex = forbidden
    declaration = json.loads(declaration_path.read_text())
    result = {"mode": declaration["mode"], "neural_calls": 0, "public_publication": False,
              "scientific_qualification": False, "actual_MD_calls": 0, "actual_CG_calls": 0}
    phase = "parent-binding"
    try:
        from scripts.mac_reconstruction_comparison import target_axis

        from dvbfixer.ffutils.geometry import assert_all_l
        from dvbfixer.model.diffusion.anchor_acceptance import _structure
        from dvbfixer.model.diffusion.contract import AtomIdentity, DiffusionRequest
        from dvbfixer.model.diffusion.protein_atoms import _CANONICAL_ONE_TO_THREE
        from dvbfixer.model.diffusion.rfd3_custody_geometry import (
            audit_geometry_with_declared_links,
        )
        from dvbfixer.model.diffusion.rfd3_custody_records import CustodyError, digest, read_records
        from dvbfixer.model.diffusion.rfd3_refinement import build_refinement_input
        from dvbfixer.model.diffusion.rfd3_terminal_restoration import (
            audit_native_terminal_transport,
            audit_terminal_derivative,
        )
        from dvbfixer.model.diffusion.rfd3_water_context import assess_water_context
        from dvbfixer.model.diffusion.runner_refinement import rewrite_generated_coordinates
        from dvbfixer.model.modeller_run import _add_chirality_restraints

        if declaration["mode"] not in {"api-control", "terminal-annealing"}:
            raise CustodyError("unknown private native procedure")
        receipt = json.loads((PARENT / "readback-receipt.json").read_text())
        for name, expected in receipt["bindings"].items():
            if digest(Path(name).read_bytes()) != expected:
                raise CustodyError("selected torsion parent changed")
        if declaration["mode"] == "terminal-annealing":
            control = Path(declaration["required_control_readback"])
            if digest(control.read_bytes()) != declaration["required_control_readback_sha256"]:
                raise CustodyError("native API control binding changed")
            control_report = json.loads(control.read_text())
            if not control_report["native_selected_MD_interface_and_fixation_passed"]:
                raise CustodyError("native selected-MD control did not pass")
            for name, expected in control_report["bindings"].items():
                if digest(Path(name).read_bytes()) != expected:
                    raise CustodyError("native API control evidence changed")
        for name in ("source.pdb", "full-context.pdb", "request.json", "approved-topology-policy.json", "raw-candidate-REJECTED.pdb"):
            with (output / name).open("xb") as stream:
                stream.write((PRIOR / name).read_bytes())
        parent = (PRIOR / "selected-terminal-PROPOSAL.pdb").read_bytes()
        with (output / "parent-torsion-PROPOSAL.pdb").open("xb") as stream:
            stream.write(parent)
        request = DiffusionRequest.from_json((output / "request.json").read_text())
        source, context, raw = [(output / n).read_bytes() for n in
                               ("source.pdb", "full-context.pdb", "raw-candidate-REJECTED.pdb")]
        policy = json.loads((output / "approved-topology-policy.json").read_text())
        active = frozenset(AtomIdentity(**a) for a in policy["active_terminal_atoms"])
        targets = {tuple(sorted(AtomIdentity(**a) for a in pair)) for pair in policy["target_pairs"]}
        if len(active) != 202 or len(targets) != 21:
            raise CustodyError("selected native graph or mask changed")
        ff = build_refinement_input(request, source, context, parent,
                                    context_sha256=digest(context), raw_sha256=digest(parent))
        input_path = output / "native-input-fixed-waters.pdb"
        with input_path.open("xb") as stream:
            stream.write(ff.text.encode())
        ff_atoms, _, _, _ = read_records(ff.text)
        expected = {a.identity: a for a in ff_atoms}
        raw_atoms, _, _, _ = read_records(raw.decode())
        raw_points = {a.identity: a.xyz for a in raw_atoms}
        residues = tuple(r for t in request.target_sequences for r in target_axis(request)[t.chain])
        sequence = "".join(t.sequence for t in request.target_sequences)
        phase = "native-topology"
        from modeller import Environ, Selection, features, forms, log, physical
        from modeller.optimizers import ConjugateGradients, MolecularDynamics
        from modeller.scripts import complete_pdb

        log.none()
        env = Environ(rand_seed=declaration["native_seed"])
        env.io.hetatm = env.io.water = True
        env.libs.topology.read(file="$(LIB)/top_heav.lib")
        env.libs.parameters.read(file="$(LIB)/par.lib")
        env.edat.dynamic_sphere = True
        env.edat.dynamic_lennard = env.edat.dynamic_coulomb = False
        ordinal = {}

        def patches(model):
            native = tuple(model.residues)
            if len(native) != 794:
                raise CustodyError("native full protein/water residue axis changed")
            for n, r, letter in zip(native[:468], residues, sequence, strict=True):
                if _CANONICAL_ONE_TO_THREE.get(n.code, n.code) != _CANONICAL_ONE_TO_THREE[letter]:
                    raise CustodyError("native ordinal sequence changed")
                ordinal[n.index] = r
            by_residue = {(r.chain, r.residue_number, r.insertion_code): n
                          for n, r in ((n, ordinal[n.index]) for n in native[:468])}
            for pair in sorted(targets):
                model.patch(residue_type="DISU", residues=tuple(by_residue[a.chain, a.residue_number, a.insertion_code] for a in pair))

        model = complete_pdb(env, str(input_path), special_patches=patches, transfer_res_num=True)
        axis = {}
        for native in tuple(model.residues)[:468]:
            residue = ordinal[native.index]
            for atom in native.atoms:
                identity = AtomIdentity(residue.chain, residue.residue_number, residue.insertion_code, atom.name)
                if identity in axis:
                    raise CustodyError("ambiguous native protein atom")
                axis[identity] = atom
        waters = [a for a in ff_atoms if a.residue_name == "HOH"]
        for native, water in zip(tuple(model.residues)[468:], waters, strict=True):
            atoms = tuple(native.atoms)
            if len(atoms) != 1 or water.identity in axis:
                raise CustodyError("ambiguous native water transport")
            axis[water.identity] = atoms[0]
        allowed_extra = {AtomIdentity(c, "136", "", "OXT") for c in ("B", "D")}
        if not expected.keys() <= axis.keys() or not axis.keys() - expected.keys() <= allowed_extra:
            raise CustodyError("native original-heavy inventory changed")
        by_index = {a.index: identity for identity, a in axis.items()}
        if len(by_index) != len(axis):
            raise CustodyError("ambiguous native atom index")
        bonds = [tuple(bond) for bond in model.bonds]
        actual = [tuple(sorted(by_index[a.index] for a in pair)) for pair in bonds if all(a.name == "SG" for a in pair)]
        if len(actual) != 21 or set(actual) != targets:
            raise CustodyError("actual native topology differs from approved21pairs")
        selected = Selection(*(axis[a] for a in sorted(active)))
        if {by_index[a.index] for a in selected} != active:
            raise CustodyError("native MD selection is not exactly202approved missingatoms")

        def coordinates():
            return {identity: (a.x, a.y, a.z) for identity, a in axis.items()}

        initial = coordinates()
        audit_native_terminal_transport({a: expected[a].xyz for a in expected}, {a: initial[a] for a in expected},
            active=active, max_displacement_angstrom=60.)
        model.restraints.make(Selection(model), restraint_type="stereo", spline_on_site=False)
        added, skipped = _add_chirality_restraints(model)
        omega_count = 0
        for first, last in bonds:
            if first.name == "N" and last.name == "C":
                first, last = last, first
            if first.name != "C" or last.name != "N" or first.residue.index == last.residue.index:
                continue
            if not {by_index[first.index], by_index[last.index]} & active:
                continue
            model.restraints.add(forms.Gaussian(group=physical.omega_dihedral,
                feature=features.Dihedral(first.residue.atoms["CA"], first, last, last.residue.atoms["CA"]),
                mean=math.pi, stdev=math.radians(5.)))
            omega_count += 1
        if added < 24 or omega_count != 24:
            raise CustodyError("native chirality/terminalamide restraints missing")
        save(output / "native-policy.json", {"actual_SS_pairs": [[asdict(a) for a in p] for p in sorted(actual)],
             "active_atoms": [asdict(a) for a in sorted(active)], "chirality_restraints": added,
             "chirality_skips": skipped, "trans_amide_restraints": omega_count,
             "initial_coordinates": [{"atom": asdict(a), "xyz_angstrom": xyz} for a, xyz in initial.items()],
             "source_sha256": digest(source), "context_sha256": digest(context), "parent_sha256": digest(parent),
             "original_raw_sha256": digest(raw), "water_serial_mapping": ff.water_serial_mapping})
        stage_reports = []
        for ordinal, stage in enumerate(declaration["schedule"], 1):
            phase = stage["name"]
            model.env.edat.dynamic_sphere = stage.get("dynamic_sphere", True)
            if stage["method"] == "MD":
                optimizer = MolecularDynamics(temperature=stage["temperature_kelvin"], md_time_step=2.,
                    cap_atom_shift=.2, init_velocities=True, equilibrate=20, md_return="FINAL")
                result["actual_MD_calls"] += 1
            elif stage["method"] == "CG":
                optimizer = ConjugateGradients(min_atom_shift=.001)
                result["actual_CG_calls"] += 1
            else:
                raise CustodyError("unknown declared native optimizer")
            optimizer.optimize(selected, max_iterations=stage["iterations"])
            path = output / f"native-stage-{ordinal:02d}-{phase}.pdb"
            model.write(file=str(path))
            current = coordinates()
            save(output / f"native-stage-{ordinal:02d}-coordinates.json", {"coordinates": [
                 {"atom": asdict(a), "xyz_angstrom": xyz} for a, xyz in current.items()]})
            transport = audit_native_terminal_transport(initial, current, active=active, max_displacement_angstrom=60.)
            raw_max = max(math.dist(raw_points[a], current[a]) for a in active)
            if raw_max > 60.:
                raise CustodyError("native annealing exceeded original-raw60Aterminalscope")
            structure = _structure(path.read_text(), path.name)
            assert_all_l(structure.pdb.topology, structure.pdb.positions)
            row = {"ordinal": ordinal, "stage": stage, "transport": asdict(transport),
                   "max_active_displacement_from_original_raw_angstrom": raw_max, "global_L": True}
            stage_reports.append(row)
            save(output / f"native-stage-{ordinal:02d}-audit.json", row)
        if not stage_reports or stage_reports[-1]["transport"]["max_active_displacement_angstrom"] <= 1e-6:
            raise CustodyError("native procedure did not actually move selected atoms")
        save(output / "native-before-materialization.json", {"stages": stage_reports,
             "all_nonactive_native_fidelity_passed": True, "global_L_passed": True, "no_posthoc_clamp": True})
        if declaration["mode"] == "terminal-annealing":
            phase = "final-materialization-and-assessment"
            final = rewrite_generated_coordinates(parent.decode(), {a: current[a] for a in active}).encode()
            custody = audit_terminal_derivative(request, source, context, raw, final, active=active,
                raw_sha256=digest(raw), context_sha256=digest(context), max_displacement_angstrom=60.)
            with (output / "annealed-terminal-UNQUALIFIED.pdb").open("xb") as stream:
                stream.write(final)
            save(output / "derivative-custody.json", asdict(custody))
            geometry = audit_geometry_with_declared_links(request, source, final)
            water = assess_water_context(request, source, context, context_sha256=digest(context), candidate=final)
            save(output / "whole-candidate-geometry.json", asdict(geometry))
            save(output / "whole-candidate-water.json", asdict(water))
            result.update(final_sha256=digest(final), geometry_passed=geometry.geometry_passed,
                          water_clear=water.restored_region_clear, full_scientific_candidate_accepted=False)
        result.update(native_procedure_completed=True, native_stage_count=len(stage_reports),
                      actual_native_SS_pairs=21, native_active_atoms=202,
                      max_nonactive_native_motion_angstrom=max(r["transport"]["max_nonactive_residual_angstrom"] for r in stage_reports),
                      global_L=True, modeller_version=__import__("modeller").__version__)
        phase = "completed"
    except Exception as error:
        frames, tb = [], error.__traceback__
        while tb is not None:
            frames.append({"file": tb.tb_frame.f_code.co_filename, "function": tb.tb_frame.f_code.co_name, "line": tb.tb_lineno})
            tb = tb.tb_next
        save(output / "failure.json", {"phase": phase, "type": type(error).__name__, "message": str(error), "frames": frames})
        result.update(native_procedure_completed=False, failure_phase=phase, failure_type=type(error).__name__, failure_message=str(error))
    result["phase"] = phase
    save(output / "result.json", result)
    print(json.dumps(result))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--declaration", type=Path, required=True)
    args = parser.parse_args()
    main(args.output.resolve(), args.declaration.resolve())
