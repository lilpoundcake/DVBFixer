#!/usr/bin/env python3
"""Private, source-bound six-case anchor study inventory and non-study FF preflight.

This entry point cannot execute the calibration grid or sample new proposals.
Temporary PDBFixer/H preparation may perform its own internal minimization;
the research refiner is stopped before its first optimization call.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import traceback
from collections import Counter
from dataclasses import asdict
from pathlib import Path

from scripts.mac_reconstruction_comparison import digest, save, target_axis
from scripts.retained_evidence import freeze_code, verify_bindings
from scripts.retained_modeller_full_context import merge_projection, read_atoms

from dvbfixer.model.diffusion.contract import AtomIdentity, DiffusionRequest, ResidueIdentity

PROJECT = Path(__file__).resolve().parents[1]
RETAINED = PROJECT / "tmp/mac-real-gap-continuation-2026-11-12"
LATEST = PROJECT / "tmp/mac-real-gap-continuation-2026-11-14"
ENTRIES = ("3C7D", "6Z38", "1NON", "2BE1", "1HE8", "1KTZ")
ENVELOPES = ((.10, .10, .50), (.25, .25, .75), (.50, .50, 1.25), (.75, .75, 1.50))
STRENGTHS = (1000., 10000., 100000.)


def variants():
    cells = [{"id": "fixed", "window": 0, "envelope": None, "strength": None}]
    for window in (1, 2):
        for index, envelope in enumerate(ENVELOPES):
            for strength in STRENGTHS:
                cells.append({"id": f"w{window}-e{index+1}-k{int(strength)}", "window": window,
                              "envelope": list(envelope), "strength": strength})
    return cells


def anchor_residues(request, width):
    """Exactly width consecutive observations OUTWARD from internal-gap flanks."""
    if type(width) is not int or width not in (0, 1, 2):
        raise ValueError("undeclared anchor window")
    axes = target_axis(request)
    observed = {r for p in request.sequence_placements for r in p.observed_residues}
    selected = set()
    for gap in request.gaps:
        if gap.left_anchor is None or gap.right_anchor is None:
            continue
        axis = axes[gap.chain]
        for anchor, step in ((gap.left_anchor, -1), (gap.right_anchor, 1)):
            start = axis.index(anchor)
            for offset in range(width):
                index = start + step * offset
                if not 0 <= index < len(axis) or axis[index] not in observed:
                    break
                selected.add(axis[index])
    return tuple(sorted(selected))


def declare_ff_c_termini(fixer, termini):
    """Explicit user-approved free C termini, auxiliary FF copy only.

    PDBFixer otherwise considers a following water/ligand to be the end of the
    chain and misses OXT. Declare exact target-chain ends, never every break or
    every residue before a heterogen. Existing external C bonds fail closed.
    """
    from dvbfixer.model.diffusion.protein_atoms import CANONICAL_HEAVY_ATOMS

    residues = {}
    for residue in fixer.topology.residues():
        identity = ResidueIdentity(residue.chain.id, str(residue.id), residue.insertionCode.strip())
        if identity in residues:
            raise ValueError("ambiguous FF terminal residue identity")
        residues[identity] = residue
    if len(set(termini)) != len(termini):
        raise ValueError("duplicate declared FF C terminus")
    events = []
    for identity in termini:
        if identity not in residues:
            raise ValueError("declared FF C terminus absent")
        residue = residues[identity]
        axis = {a.name: a for a in residue.atoms()}
        if residue.name not in CANONICAL_HEAVY_ATOMS or not {"N", "CA", "C", "O"} <= axis.keys():
            raise ValueError("unsupported/incomplete declared FF C terminus")
        if any((a == axis["C"] and b.residue != residue) or (b == axis["C"] and a.residue != residue)
               for a, b in fixer.topology.bonds()):
            raise ValueError("declared free FF C terminus has an external covalent C bond")
        if "OXT" in axis:
            continue
        previous = tuple(fixer.missingTerminals.get(residue, ()))
        if "OXT" not in previous:
            fixer.missingTerminals[residue] = [*previous, "OXT"]
        events.append({"residue": asdict(identity), "resname": residue.name,
                       "auxiliary_atom": "OXT", "already_planned_by_pdbfixer": "OXT" in previous,
                       "output_ownership_extended": False})
    return events


def selected_inputs(entry):
    review = RETAINED / f"retained-source-review-matrix-attempt01/execution/retained-source-review-seed-0/report/entries/{entry}.json"
    if entry == "3C7D":
        root = LATEST / "3C7D-full-context-restoration-attempt01"
        case = root / "cases/3C7D"
        native = root / "execution/3C7D-seed-104/native"
        return {"request": case / "request.json", "source": case / "source.pdb",
                "context": case / "prepared-full-context.pdb", "proposal": native / "raw-full-context.pdb",
                "lineage": native / "stages.json", "review": review,
                "proposal_kind": "MODELLER seed104 raw full-context; no continuation selected"}
    if entry == "6Z38":
        root = LATEST / "6Z38-modeller-restoration-attempt02"
        row = json.loads(review.read_text())
        return {"request": root / "cases/6Z38/request.json", "source": root / "cases/6Z38/source.pdb",
                "context": Path(row["request_path"]).parent / "full-assembly-context.pdb",
                "proposal": root / "execution/6Z38-seed-102/source-restored.pdb",
                "lineage": root / "outcome.json", "review": review,
                "proposal_kind": "retained source-restored MODELLER derivative; native fixation failure is not erased"}
    if entry in ("1NON", "2BE1"):
        root = RETAINED / "retained-modeller-fixation-seed101-attempt01"
        row = json.loads(review.read_text())
        return {"request": root / f"cases/{entry}/request.json", "source": root / f"cases/{entry}/source.pdb",
                "context": Path(row["request_path"]).parent / "full-assembly-context.pdb",
                "proposal": root / f"execution/{entry}-seed-101/source-restored.pdb",
                "lineage": root / f"execution/{entry}-seed-101/result.json", "review": review,
                "proposal_kind": "MODELLER seed101 source-restored derivative; native fixation failure retained"}
    if entry == "1HE8":
        root = PROJECT / "tmp/mac-reconstruction-2026-11-11/1he8-two-anchor-engineering-feasibility"
        return {"request": root / "request.json", "source": root / "source.pdb",
                "context": RETAINED / "assembly-smoke-1HE8-attempt02/assembly1.pdb",
                "proposal": root / "input.pdb", "lineage": root / "declaration.json",
                "proposal_kind": "retained richer-context patchwork input; full-context mapping must pass independently"}
    if entry == "1KTZ":
        root = LATEST / "1KTZ-restored-heavy-export-attempt02"
        result = root / "execution/1KTZ-restored-heavy-export-seed-0/result"
        return {"request": result / "original-request.json", "source": result / "original-source-protein.pdb",
                "context": result / "original-full-context.pdb",
                "proposal": result / "1KTZ-restored-heavy-with-waters.pdb",
                "lineage": root / "readback-receipt.json", "review": review,
                "proposal_kind": "retained reviewed heavy export; 710 additions; 21 SS target; terminal OXT deferred"}
    raise ValueError("undeclared entry")


def same_observation(first, second):
    # Context insertion may disambiguate serials; no other observed field changes.
    return first[:6] == second[:6] and first[11:] == second[11:]


def prepared_context(source, original):
    """Replay only the explicitly bound historical protein preparation, no repair.

    Keep original context bytes separately. The calibrated reference is the saved
    prepared protein, not a claim that historical preparation kept all observations
    fixed. Original tensors and nonprotein records remain untouched.
    """
    import numpy as np

    from dvbfixer.model.diffusion.runner_refinement import rewrite_generated_coordinates

    protein = read_atoms(source)
    context = read_atoms(original, allow_fractional_heterogens=True)
    original_protein = {a for a in context if context[a].line.startswith("ATOM  ")}
    if not original_protein <= protein.keys():
        raise ValueError("historical preparation removed original protein observations")
    changed = {}
    for atom in original_protein:
        old, new = context[atom], protein[atom]
        if old.line[:6] != new.line[:6] or old.line[11:30] != new.line[11:30] or old.line[54:] != new.line[54:]:
            raise ValueError("historical preparation changed observed noncoordinate metadata")
        if old.xyz != new.xyz:
            changed[atom] = np.asarray(new.xyz)
    replayed = rewrite_generated_coordinates(original, changed)
    merged = merge_projection(source, replayed)
    evidence = {"historical_prepared_coordinate_changes": [
        {"atom": asdict(a), "original_xyz": context[a].xyz, "prepared_xyz": protein[a].xyz,
         "displacement_angstrom": math.dist(context[a].xyz, protein[a].xyz)} for a in sorted(changed)],
        "historical_added_fixed_atoms": [asdict(a) for a in sorted(protein.keys()-original_protein)],
        "new_source_repair_executed": False, "reference": "bound historical prepared source",
        "original_context_retained_separately": True}
    return merged, evidence


def validate_proposal(request, source, context, proposal):
    """No copying may hide changed observations in the selected base proposal."""
    protein = read_atoms(source)
    original = read_atoms(context, allow_fractional_heterogens=True)
    candidate = read_atoms(proposal, allow_fractional_heterogens=True)
    if set(protein) != set(request.fixed_atoms):
        raise ValueError("source observation mask differs from frozen request")
    if not set(protein) <= candidate.keys():
        raise ValueError("proposal lost source observations")
    if any(not same_observation(candidate[a].line, protein[a].line) for a in protein):
        raise ValueError("selected proposal observations are not exact source records")
    generated = set(request.generated_atoms)
    context_only = set(original) - set(protein)
    if generated & (set(protein) | context_only):
        raise ValueError("generated/context ownership overlaps")
    if set(candidate) != set(protein) | context_only | generated:
        raise ValueError("full-context proposal atom inventory differs")
    for atom in context_only:
        if candidate[atom].line != original[atom].line:
            raise ValueError("proposal changed a context observation")
    # Coordinates are the only later mutable fields, including for anchor atoms.
    return {"observed_protein_atoms": len(protein), "context_atoms": len(context_only),
            "generated_atoms": len(generated), "base_observations_exact": True}


def declare(output):
    import openmm.app
    import pdbfixer

    output.mkdir(parents=True, exist_ok=False)
    # Bind immutable protocol bytes, not the living handoff updated after a run.
    protocol = output / "protocol.md"
    protocol.write_bytes((PROJECT / "docs/plans/experimental-structure-completion-handoff.md").read_bytes())
    files = [Path(__file__), protocol, PROJECT / "scripts/audit_anchor_calibration_preflight.py"]
    files += [PROJECT / f"scripts/{name}.py" for name in (
        "mac_reconstruction_comparison", "retained_evidence", "retained_modeller_full_context",
        "run_mac_modeller_comparison", "run_mac_reconstruction_guarded")]
    files += list((PROJECT / "src/dvbfixer/model/diffusion").glob("*.py"))
    files += [PROJECT / "src/dvbfixer/ffutils/geometry.py"]
    files += [PROJECT / "src/dvbfixer/ffutils/__init__.py", PROJECT / "src/dvbfixer/profile_readiness.py"]
    files += sorted((Path(openmm.app.__file__).parent / "data").rglob("*.xml"))
    for directory in ("templates", "soft.xml"):
        resource = Path(pdbfixer.__file__).parent / directory
        files += [resource] if resource.is_file() else sorted(p for p in resource.rglob("*") if p.is_file())
    cases, jobs = [], []
    for entry in ENTRIES:
        row = {"entry": entry, "grid_cells": variants(), "preparation_status": "pending"}
        try:
            inputs = selected_inputs(entry)
            row["proposal_kind"] = inputs.pop("proposal_kind")
            if any(not p.is_file() for p in inputs.values()):
                raise FileNotFoundError([str(p) for p in inputs.values() if not p.is_file()])
            row["input_bindings"] = {str(p): digest(p) for p in inputs.values()}
            files += list(inputs.values())
            case = output / "cases" / entry
            case.mkdir(parents=True)
            request = DiffusionRequest.from_json(inputs["request"].read_text())
            source = inputs["source"].read_text()
            if digest(inputs["source"]) != request.normalized_pdb.sha256:
                raise ValueError("request/source digest mismatch")
            context = inputs["context"].read_text()
            proposed = inputs["proposal"].read_text()
            if "review" in inputs:
                historical_review = json.loads(inputs["review"].read_text())
                if digest(inputs["source"]) != historical_review["source_sha256"]:
                    raise ValueError("historical review/source digest mismatch")
            lineage = json.loads(inputs["lineage"].read_text())
            if "restored_sha256" in lineage and digest(inputs["proposal"]) != lineage["restored_sha256"]:
                raise ValueError("historical derivative digest mismatch")
            if entry not in ("1KTZ", "3C7D"):
                # Only merge after auditing protein observations; never clamp drift.
                protein = read_atoms(source)
                observed = read_atoms(proposed)
                if any(a not in observed or not same_observation(observed[a].line, protein[a].line) for a in protein):
                    raise ValueError("base derivative does not preserve exact protein records")
                context_view, preparation_evidence = prepared_context(source, context)
                proposed = merge_projection(proposed, context_view)
                save(case / "historical-preparation.json", preparation_evidence)
                row["historical_prepared_coordinate_changes"] = len(preparation_evidence["historical_prepared_coordinate_changes"])
            row.update(validate_proposal(request, source, context, proposed))
            for name, text in (("request.json", request.to_json()), ("source.pdb", source),
                               ("context.pdb", context), ("proposal.pdb", proposed)):
                (case / name).write_text(text)
            row["windows"] = {str(w): [asdict(r) for r in anchor_residues(request, w)] for w in (0, 1, 2)}
            row["gaps"] = [{"chain": g.chain, "length": g.target_interval.stop-g.target_interval.start,
                            "internal": g.left_anchor is not None and g.right_anchor is not None}
                           for g in request.gaps]
            row["context_heterogen_atoms"] = dict(Counter(a.name for a in read_atoms(context, allow_fractional_heterogens=True).values()
                                                         if a.line.startswith("HETATM")))
            row["preparation_status"] = "bound; FF preflight pending"
            files += list(case.iterdir())
            jobs.append({"case_id": entry, "case": str(case), "seed": 7,
                         "command": [sys.executable, "-m", "scripts.anchor_calibration_preflight", "worker",
                                     "--case", str(case), "--output", str(output / f"execution/{entry}-seed-7")]})
        except (ValueError, OSError) as exc:
            row.update(preparation_status="input_blocked", error_type=type(exc).__name__, reason=str(exc))
        cases.append(row)
    save(output / "inventory.json", {"cases": cases, "grid_executed": 0,
                                    "no_new_proposals": True, "selection_before_new_optimization": True})
    files.append(output / "inventory.json")
    save(output / "plan.json", {"protocol": "anchor-calibration-preflight-v1", "python": sys.executable,
        "python_sha256": digest(Path(sys.executable)), "source_code_and_input_sha256": {str(p): digest(p) for p in files},
        "jobs": jobs, "count": len(jobs), "study_optimizer_calls": 0, "new_neural_calls": 0,
        "grid_execution_enabled": False, "temporary_chemistry_preparation_may_minimize": True,
        "explicit_FF_only_C_terminal_policy": "user-approved free uncapped target-chain C termini at pH7; auxiliary OXT only; original/output/dat unchanged",
        "resource_policy": "serial1800s/job,RSS8GiB,host-pressure,10GiB evidence"})
    freeze_code(output / "plan.json")


class PreflightCaptured(RuntimeError):
    """Deliberate stop before the first research minimization call."""


def worker(case, output):
    import numpy as np
    from openmm import XmlSerializer
    from openmm.unit import dalton, kilojoule_per_mole, nanometer

    from dvbfixer.ffutils.geometry import assert_all_l
    from dvbfixer.model.diffusion import boundary_refinement as kernel
    from dvbfixer.model.diffusion.anchor_acceptance import _structure
    from dvbfixer.model.diffusion.anchor_policy import AnchorMovementLimits, build_anchor_restraints
    from dvbfixer.profile_readiness import pdbfixer_resources_ready

    if {p.name for p in output.iterdir()} != {"command.json", "output.log"}:
        raise ValueError("preflight requires a fresh guarded reservation")
    row = {"entry": case.name, "study_optimizer_calls": 0, "grid_execution_ready": False,
           "temporary_preparation_may_minimize": True, "system_preflight_passed": False}
    old_simulation = kernel.Simulation
    old_fixer = kernel.PDBFixer
    old_freeze = kernel._freeze_observed_heavy_atoms
    terminal_events = []
    try:
        if not pdbfixer_resources_ready():
            raise ValueError("PDBFixer resources unusable; no piecemeal repair")
        request = DiffusionRequest.from_json((case / "request.json").read_text())
        source, proposal = [(case / n).read_text() for n in ("source.pdb", "proposal.pdb")]
        parsed = _structure(proposal, "retained-full-context.pdb")
        # No unreported chirality repair of the base proposal during preflight.
        assert_all_l(parsed.pdb.topology, parsed.pdb.positions)
        expected = read_atoms(proposal, allow_fractional_heterogens=True)
        termini = tuple(axis[-1] for axis in target_axis(request).values())
        approved_ss_pairs = None
        if case.name == "1KTZ":
            from dvbfixer.model.diffusion.rfd3_custody_records import read_records

            _, _, annotated_pairs, _ = read_records(proposal)
            approved_ss_pairs = {tuple(sorted(pair)) for pair in annotated_pairs
                                 if all(a.atom_name == "SG" for a in pair)}
            if len(approved_ss_pairs) != 21:
                raise ValueError("retained approved candidate annotations do not contain21SSpairs")
        restraints = build_anchor_restraints(source, anchor_residues(request, 2),
            observed_atoms=request.fixed_atoms, base_k_kj_mol_nm2=10000., bfactor_source="disabled")
        mobile = set(request.generated_atoms) | {r.atom for r in restraints}
        captured = []

        class TemporaryFixer(old_fixer):
            def findMissingAtoms(self):
                super().findMissingAtoms()
                terminal_events.extend(declare_ff_c_termini(self, termini))

        def capture_base_system(system, topology, original, movable):
            # Capture before mass freezing and auxiliary/window-dependent forces.
            # The future cell adapter must reuse this chemistry and initial NPZ,
            # not rerun hydrogen placement independently in each cell.
            (output / "base-system.xml").write_text(XmlSerializer.serialize(system))
            old_freeze(system, topology, original, movable)

        class CaptureSimulation(old_simulation):
            def minimizeEnergy(self, *args, **kwargs):
                atoms = list(self.topology.atoms())
                identities = [AtomIdentity(a.residue.chain.id, str(a.residue.id), a.residue.insertionCode.strip(), a.name) for a in atoms]
                if len(set(identities)) != len(identities):
                    raise ValueError("temporary topology has duplicate atom identities")
                state = self.context.getState(getPositions=True, getEnergy=True, getForces=True)
                xyz = state.getPositions(asNumpy=True).value_in_unit(nanometer)
                energy = float(state.getPotentialEnergy().value_in_unit(kilojoule_per_mole))
                if not np.isfinite(xyz).all() or not math.isfinite(energy):
                    raise ValueError("initial full-context energy/coordinates nonfinite")
                positions = dict(zip(identities, xyz * 10, strict=True))
                if not expected.keys() <= positions.keys():
                    raise ValueError("temporary model lost retained heavy atoms")
                for identity, atom in expected.items():
                    if not np.allclose(positions[identity], atom.xyz, rtol=0, atol=1.e-4):
                        raise ValueError(f"temporary preparation changed retained heavy atom {identity}")
                fixed = {identities[a.index] for a in atoms if self.system.getParticleMass(a.index).value_in_unit(dalton) == 0}
                if fixed != expected.keys() - mobile:
                    raise ValueError("actual fixed heavy particles differ from declared window")
                bonds = list(self.topology.bonds())
                ss_pairs = [tuple(sorted((identities[a.index], identities[b.index]))) for a, b in bonds if a.name == b.name == "SG"]
                if approved_ss_pairs is not None and set(ss_pairs) != approved_ss_pairs:
                    raise ValueError("1KTZ actual topology differs from exact confirmed21SSpairs")
                rows = [{"index": a.index, "identity": asdict(identities[a.index]), "resname": a.residue.name,
                         "chain_index": a.residue.chain.index, "residue_index": a.residue.index,
                         "element": a.element.symbol if a.element else None,
                         "mass_dalton": float(self.system.getParticleMass(a.index).value_in_unit(dalton)),
                         "temporary_addition": identities[a.index] not in expected} for a in atoms]
                (output / "temporary-system.xml").write_text(XmlSerializer.serialize(self.system))
                save(output / "temporary-axis.json", {"atoms": rows,
                    "bonds": [[a.index, b.index] for a, b in bonds],
                    "SS_pairs": [[asdict(a), asdict(b)] for a, b in ss_pairs]})
                with (output / "initial-state.npz").open("xb") as stream:
                    np.savez_compressed(stream, positions_nm=xyz,
                        full_force_kj_mol_nm=state.getForces(asNumpy=True).value_in_unit(kilojoule_per_mole/nanometer))
                captured.append({"particles": len(atoms), "fixed_heavy": len(fixed),
                    "temporary_additions": sum(r["temporary_addition"] for r in rows),
                    "SS_pairs": len(ss_pairs), "initial_full_energy_kj_mol": energy,
                    "stopped_before_minimize": True})
                raise PreflightCaptured()

        kernel.Simulation = CaptureSimulation
        kernel.PDBFixer = TemporaryFixer
        kernel._freeze_observed_heavy_atoms = capture_base_system
        try:
            kernel.refine_generated_region(proposal,
                generated_residues=tuple(r for gap in request.gaps for r in gap.generated_residues),
                generated_atoms=request.generated_atoms, max_iterations=1000, restart_count=0,
                random_seed=7, platform_name="CPU", geometry_preconditioning=True,
                observed_anchor_restraints=restraints, anchor_movement_limits=AnchorMovementLimits())
        except PreflightCaptured:
            pass
        if len(captured) != 1:
            raise ValueError("preflight did not reach exactly one pre-minimization state")
        row.update(system_preflight_passed=True, system=captured[0],
                   remaining_gate="shared-state grid adapter and independent readback not implemented/verified yet")
    except Exception as exc:
        row.update(error_type=type(exc).__name__, reason=str(exc))
        (output / "failure.txt").write_text(traceback.format_exc())
    finally:
        kernel.Simulation = old_simulation
        kernel.PDBFixer = old_fixer
        kernel._freeze_observed_heavy_atoms = old_freeze
    row["explicit_FF_only_terminal_additions"] = terminal_events
    save(output / "result.json", row)
    print(json.dumps(row, ensure_ascii=False))


def summarize(output):
    inventory = json.loads((output / "inventory.json").read_text())
    results = []
    for row in inventory["cases"]:
        path = output / f"execution/{row['entry']}-seed-7/result.json"
        result = json.loads(path.read_text()) if path.is_file() else None
        results.append({"entry": row["entry"], "input_status": row["preparation_status"],
                        "input_reason": row.get("reason"), "forcefield_preflight": result})
    save(output / "preflight-summary.json", {"cases": results, "study_variants_executed": 0,
        "scientific_qualification": False, "grid_execution_ready": False})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("declare", "run", "worker", "summarize"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--case", type=Path)
    args = parser.parse_args()
    if args.mode == "declare":
        declare(args.output.resolve())
    elif args.mode == "worker":
        if args.case is None:
            parser.error("worker requires --case")
        worker(args.case.resolve(), args.output.resolve())
    elif args.mode == "run":
        from scripts.run_mac_modeller_comparison import execute
        verify_bindings(json.loads((args.output / "plan.json").read_text())["source_code_and_input_sha256"])
        execute(args.output.resolve())
        summarize(args.output.resolve())
    else:
        summarize(args.output.resolve())
