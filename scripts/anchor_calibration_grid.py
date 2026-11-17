#!/usr/bin/env python3
"""Private paired anchor grid from a single frozen FF initial state per case."""
from __future__ import annotations

import argparse
import json
import math
import sys
import traceback
from dataclasses import asdict
from pathlib import Path

import numpy as np
from scripts.anchor_calibration_preflight import anchor_residues, variants
from scripts.audit_retained_full_context_controls import atoms as independent_atoms
from scripts.mac_reconstruction_comparison import digest, save
from scripts.retained_evidence import freeze_code, verify_bindings
from scripts.retained_modeller_full_context import findings, read_atoms

from dvbfixer.model.diffusion.contract import AtomIdentity, DiffusionRequest, ResidueIdentity


class GridNumericalError(RuntimeError):
    """Finite-state failure, distinct from an implementation/identity error."""


def identity(atom):
    return AtomIdentity(atom.residue.chain.id, str(atom.residue.id), atom.residue.insertionCode.strip(), atom.name)


def read_temporary(snapshot):
    from openmm import XmlSerializer
    from openmm.app import Topology, element

    axis = json.loads((snapshot / "temporary-axis.json").read_text())
    topology = Topology()
    chains, residues, particles = {}, {}, []
    for index, row in enumerate(axis["atoms"]):
        if row["index"] != index:
            raise ValueError("noncanonical frozen particle order")
        atom = AtomIdentity(**row["identity"])
        ci, ri = row["chain_index"], row["residue_index"]
        if ci not in chains:
            chains[ci] = topology.addChain(atom.chain)
        if ri not in residues:
            residues[ri] = topology.addResidue(row["resname"], chains[ci], atom.residue_number,
                                              insertionCode=atom.insertion_code)
        particles.append(topology.addAtom(atom.atom_name, element.get_by_symbol(row["element"]), residues[ri]))
    for first, second in axis["bonds"]:
        topology.addBond(particles[first], particles[second])
    if [identity(a) for a in topology.atoms()] != [AtomIdentity(**r["identity"]) for r in axis["atoms"]]:
        raise ValueError("reconstructed topology particle axis differs")
    system = XmlSerializer.deserialize((snapshot / "base-system.xml").read_text())
    with np.load(snapshot / "initial-state.npz", allow_pickle=False) as data:
        positions = data["positions_nm"].copy()
    if positions.shape != (len(particles), 3) or not np.isfinite(positions).all():
        raise ValueError("invalid common frozen initial state")
    if system.getNumParticles() != len(particles) or system.getNumConstraints():
        raise ValueError("base system axis/constraints differ")
    return topology, system, positions


def cell_by_id(name):
    matches = [cell for cell in variants() if cell["id"] == name]
    if len(matches) != 1:
        raise ValueError("undeclared grid cell")
    return matches[0]


def declared_worker_job(plan, case, snapshot, cell, output):
    matches = [job for job in plan["jobs"] if job["case_id"] == f"{case.name}-{cell['id']}"]
    if len(matches) != 1:
        raise ValueError("worker cell is not in the frozen execution manifest")
    job = matches[0]
    command = job["command"]
    if (Path(job["case"]).resolve() != case.resolve()
            or Path(command[command.index("--snapshot")+1]).resolve() != snapshot.resolve()
            or Path(command[command.index("--output")+1]).resolve() != output.resolve()):
        raise ValueError("worker paths differ from frozen execution manifest")
    return job


def build_cell(request, source, proposal, topology, system, positions, cell):
    from openmm.unit import nanometer

    from dvbfixer.ffutils.geometry import build_ca_chirality_force
    from dvbfixer.model.diffusion import boundary_refinement as kernel
    from dvbfixer.model.diffusion.anchor_policy import AnchorMovementLimits, build_anchor_restraints
    from dvbfixer.model.diffusion.geometry_conditioning import (
        add_anchor_boundary_walls,
        add_geometry_conditioning,
    )

    if cell != cell_by_id(cell["id"]):
        raise ValueError("grid parameters differ from frozen cell")
    restraints = build_anchor_restraints(source, anchor_residues(request, cell["window"]),
        observed_atoms=request.fixed_atoms, base_k_kj_mol_nm2=cell["strength"] or 10000., bfactor_source="disabled")
    limits = None if cell["envelope"] is None else AnchorMovementLimits(
        backbone_rmsd_angstrom_max=cell["envelope"][0], ca_displacement_angstrom_max=cell["envelope"][1],
        heavy_displacement_angstrom_max=cell["envelope"][2])
    mobile = set(request.generated_atoms) | {r.atom for r in restraints}
    original = {a: np.asarray(r.xyz) for a, r in read_atoms(proposal, allow_fractional_heterogens=True).items()}
    kernel._freeze_observed_heavy_atoms(system, topology, original, mobile)
    if restraints:
        system.addForce(kernel._build_anchor_position_force(topology, restraints))
        add_anchor_boundary_walls(system, topology, restraints, limits)
    chirality, count = build_ca_chirality_force(topology, positions*nanometer)
    if count:
        system.addForce(chirality)
    residues = {ResidueIdentity(a.chain, a.residue_number, a.insertion_code) for a in mobile}
    for force, count in (kernel._build_trans_amide_force(topology, residues),
                         kernel._build_peptide_angle_force(topology, residues)):
        if count:
            system.addForce(force)
    add_geometry_conditioning(system)
    return restraints, limits, mobile


def independent_coordinate_audit(before, after, mobile):
    """Physical record audit independent of the coordinate rewriter."""
    old, new = before.splitlines(keepends=True), after.splitlines(keepends=True)
    if len(old) != len(new):
        raise ValueError("grid output changed physical record count")
    seen = set()
    for first, second in zip(old, new, strict=True):
        atom = None
        if first.startswith(("ATOM  ", "HETATM")):
            atom = AtomIdentity(first[21], first[22:26].strip(), first[26].strip(), first[12:16].strip())
        if atom in mobile:
            if first[:30] != second[:30] or first[54:] != second[54:]:
                raise ValueError("mobile atom metadata changed")
            seen.add(atom)
        elif first != second:
            raise ValueError("outside-window/auxiliary record changed")
    if seen != mobile:
        raise ValueError("mobile atom record inventory differs")


def independent_movements(source_path, candidate_path, restraints, limits):
    source, candidate = independent_atoms(source_path), independent_atoms(candidate_path)
    grouped = {}
    for restraint in restraints:
        a = restraint.atom
        atom = (a.chain, a.residue_number, a.insertion_code, a.atom_name)
        displacement = math.dist(source[atom]["xyz"], candidate[atom]["xyz"])
        grouped.setdefault(atom[:3], {})[atom[3]] = displacement
    rows = []
    for residue, values in sorted(grouped.items()):
        backbone = [values[name] for name in ("N", "CA", "C", "O")]
        rmsd = math.sqrt(sum(v*v for v in backbone)/4)
        maximum, ca = max(values.values()), values["CA"]
        rows.append({"residue": residue, "backbone_rmsd": rmsd, "ca_displacement": ca,
                     "heavy_max": maximum, "within_limits": rmsd <= limits.backbone_rmsd_angstrom_max
                     and ca <= limits.ca_displacement_angstrom_max and maximum <= limits.heavy_displacement_angstrom_max})
    return rows


def observed_environment_metrics(topology, before, after, observed):
    """Descriptive proximity and intrasegment C3 changes, never acceptance gates."""
    from MDAnalysis.analysis.dssp import assign
    from scipy.spatial import cKDTree

    def residue(a):
        return a.chain, a.residue_number, a.insertion_code
    bonded_residues = {tuple(sorted((residue(identity(a)), residue(identity(b)))))
                       for a, b in topology.bonds() if a.residue != b.residue}
    axis = sorted(observed)
    def contacts(coordinates):
        points = np.array([coordinates[a] for a in axis])
        pairs = set()
        for first, second in cKDTree(points).query_pairs(4.5, output_type="ndarray"):
            ra, rb = residue(axis[first]), residue(axis[second])
            pair = tuple(sorted((ra, rb)))
            if ra != rb and pair not in bonded_residues:
                pairs.add(pair)
        return pairs
    old, new = contacts(before), contacts(after)
    # Keep chain/break boundaries; this descriptive metric does not claim a full
    # interchain DSSP assignment. Compare the SAME initial completed proposal.
    peptide_edges = {(a.residue, b.residue) if a.name == "C" else (b.residue, a.residue)
                     for a, b in topology.bonds() if {a.name, b.name} == {"C", "N"} and a.residue != b.residue}
    segments, segment = [], []
    for chain in topology.chains():
        for current in chain.residues():
            names = {a.name: identity(a) for a in current.atoms() if a.name in {"N", "CA", "C", "O", "H"}}
            complete = {"N", "CA", "C", "O"} <= names.keys() and all(a in before and a in after for a in names.values())
            if segment and (not complete or (segment[-1][0], current) not in peptide_edges):
                segments.append(segment)
                segment = []
            if complete:
                segment.append((current, names))
        if segment:
            segments.append(segment)
            segment = []
    secondary = []
    for segment in segments:
        if len(segment) < 6:
            continue
        codes = []
        for coordinates in (before, after):
            # Actual auxiliary FF H where present. H=N gives zero donor energy
            # for PRO/termini without a canonical backbone H, rather than
            # fabricating an amide hydrogen on proline.
            xyz = np.array([[*(coordinates[names[name]] for name in ("N", "CA", "C", "O")),
                             coordinates[names.get("H", names["N"])]] for _, names in segment])
            codes.append(np.array(["-", "H", "E"])[np.argmax(assign(xyz), axis=-1)])
        for (_, names), first, second in zip(segment, *codes, strict=True):
            if names["CA"] in observed and first != second:
                secondary.append({"residue": residue(names["CA"]), "before": str(first), "after": str(second)})
    return {"proximity_cutoff_angstrom": 4.5, "proximity_metric": "observed-heavy residue proximity; exclude covalently neighboring residues",
            "lost_contacts": sorted(old-new), "new_contacts": sorted(new-old),
            "secondary_method": "MDAnalysis DSSP C3, auxiliary FF H; missing H donors disabled; continuous peptide segments >=6 residues; no intersegment/interchain sheets",
            "secondary_reference": "same completed initial proposal, not experimental missing-coordinate truth",
            "secondary_changes": secondary}


def worker(case, snapshot, cell, output):
    from openmm import OpenMMException, Platform, VerletIntegrator, XmlSerializer
    from openmm.app import Simulation
    from openmm.unit import dalton, kilojoule_per_mole, nanometer, picosecond
    from scripts.anchor_calibration_preflight import prepared_context
    from scripts.audit_anchor_calibration_grid import evaluation_view
    from scripts.run_retained_6z38_restoration import compare_findings
    from scripts.trace_private_minimization import PassiveReporter

    from dvbfixer.ffutils.dat import DatRecord
    from dvbfixer.ffutils.geometry import ChiralityError, assert_all_l
    from dvbfixer.model.diffusion.anchor_acceptance import _structure
    from dvbfixer.model.diffusion.anchor_policy import check_anchor_movements
    from dvbfixer.model.diffusion.geometry_conditioning import (
        BONDED_GROUP,
        FULL_NONBONDED_GROUP,
        SOFTCORE_GROUP,
    )
    from dvbfixer.model.diffusion.runner_refinement import rewrite_generated_coordinates
    from dvbfixer.model.diffusion.staged_runtime import assessment_adapter

    if {p.name for p in output.iterdir()} != {"command.json", "output.log"}:
        raise ValueError("grid worker requires fresh guarded reservation")
    plan_path = output.parent.parent / "plan.json"
    plan = json.loads(plan_path.read_text())
    verify_bindings(plan["source_code_and_input_sha256"])
    declared_worker_job(plan, case, snapshot, cell, output)
    row = {"entry": case.name, "cell": cell, "completed": False, "public_acceptance": False,
           "initial_state_sha256": digest(snapshot / "initial-state.npz"), "restarts": 0}
    try:
        request = DiffusionRequest.from_json((case / "request.json").read_text())
        source, proposal = [(case / name).read_text() for name in ("source.pdb", "proposal.pdb")]
        topology, system, initial = read_temporary(snapshot)
        restraints, limits, mobile = build_cell(request, source, proposal, topology, system, initial, cell)
        axis = list(topology.atoms())
        identities = [identity(a) for a in axis]
        native_before = dict(zip(identities, initial*10, strict=True))
        observed = set(read_atoms(proposal, allow_fractional_heterogens=True))-set(request.generated_atoms)
        fixed = observed-mobile
        actual_fixed = {identity(a) for a in axis if system.getParticleMass(a.index).value_in_unit(dalton) == 0}
        if actual_fixed != fixed:
            raise ValueError("actual grid fixed-particle mask differs")
        save(output / "mobility.json", {"cell": cell, "restraints": [asdict(r) for r in restraints],
             "mobile": [asdict(a) for a in sorted(mobile)], "fixed_count": len(fixed)})
        (output / "system.xml").write_text(XmlSerializer.serialize(system))
        integrator = VerletIntegrator(.001*picosecond)
        simulation = Simulation(topology, system, integrator, Platform.getPlatformByName("CPU"), {"Threads": "1"})
        simulation.context.setPositions(initial*nanometer)
        bonded, full = 1 << BONDED_GROUP, (1 << BONDED_GROUP) | (1 << FULL_NONBONDED_GROUP)
        calls = []
        for index, (name, groups, scale, alpha, cap) in enumerate((
            ("bonded-geometry", bonded, 0., 1., 100),
            ("soft-steric-1", bonded | (1 << SOFTCORE_GROUP), 1., 1., 75),
            ("soft-steric-2", bonded | (1 << SOFTCORE_GROUP), 1., .2, 75),
            ("full-forcefield", full, 0., .2, 1000))):
            integrator.setIntegrationForceGroups(groups)
            simulation.context.setParameter("soft_scale", scale)
            simulation.context.setParameter("soft_alpha", alpha)
            reporter = PassiveReporter(output / f"phase-{index+1}")
            simulation.minimizeEnergy(tolerance=10*kilojoule_per_mole/nanometer, maxIterations=cap, reporter=reporter)
            state = simulation.context.getState(getPositions=True, getEnergy=True, getForces=True, groups=groups)
            positions = state.getPositions(asNumpy=True).value_in_unit(nanometer)
            forces = state.getForces(asNumpy=True).value_in_unit(kilojoule_per_mole/nanometer)
            energy = float(state.getPotentialEnergy().value_in_unit(kilojoule_per_mole))
            with (output / f"phase-{index+1}/returned-state.npz").open("xb") as stream:
                np.savez_compressed(stream, positions_nm=positions, active_forces=forces)
            if not np.isfinite(positions).all() or not np.isfinite(forces).all() or not math.isfinite(energy):
                raise GridNumericalError("nonfinite optimizer state")
            assert_all_l(topology, state.getPositions())
            calls.append({"phase": name, "max_iterations": cap, "active_groups": groups,
                          "callback_count": reporter.calls, "last_iteration_index": reporter.last_iteration,
                          "energy_kj_mol": energy, "termination_reason": "not exposed by API; convergence not inferred"})
            save(output / f"phase-{index+1}/completed.json", calls[-1])
        native_after = dict(zip(identities, positions*10, strict=True))
        if any(not np.allclose(native_after[a], native_before[a], rtol=0, atol=1.e-4) for a in fixed):
            raise ValueError("native optimizer moved fixed heavy atoms")
        final = rewrite_generated_coordinates(proposal, {a: native_after[a] for a in mobile})
        independent_coordinate_audit(proposal, final, mobile)
        (output / "refined.pdb").write_text(final)
        serialized = _structure(final, "grid-serialized-final.pdb")
        assert_all_l(serialized.pdb.topology, serialized.pdb.positions)
        native_moves = check_anchor_movements(restraints, native_after, limits) if restraints else ()
        serialized_moves = independent_movements(case / "source.pdb", output / "refined.pdb", restraints, limits)
        baseline = (case / "context.pdb").read_text()
        if case.name != "1KTZ":
            baseline, _ = prepared_context(source, baseline)
        baseline_findings = findings(baseline, output / "baseline-diagnostic.pdb")
        final_findings = findings(final, output / "final-diagnostic.pdb")
        full_diagnostics = compare_findings(baseline_findings, final_findings)
        anchor_atoms = {r.atom for r in restraints}
        protein_view = evaluation_view(request, source, output / "refined.pdb", anchor_atoms)
        (output / "protein-evaluation-view.pdb").write_text(protein_view)
        assessment = assessment_adapter(request, source)(protein_view, anchor_atoms)
        record = DatRecord()
        initial_atoms = read_atoms(proposal, allow_fractional_heterogens=True)
        record.added_atoms = [{"chain": a.chain, "resid": a.residue_number, "icode": a.insertion_code,
                               "resname": initial_atoms[a].name, "atom": a.atom_name,
                               "element": initial_atoms[a].element} for a in request.generated_atoms]
        record.save(output / "generated.dat")
        decoded = DatRecord.load(output / "generated.dat")
        owned = {AtomIdentity(a["chain"], a["resid"], a["icode"], a["atom"]) for a in decoded.added_atoms}
        if owned != set(request.generated_atoms) or len(decoded.added_atoms) != len(owned) or owned & observed:
            raise ValueError("independent canonical dat readback differs from generated-only ownership")
        generated_distances = [float(np.linalg.norm(native_after[a]-native_before[a])) for a in request.generated_atoms]
        row.update(completed=True, minimizer_calls=calls, native_fixed_preserved=True,
            outside_window_records_unchanged=True, global_all_l=True,
            native_anchor_movements=[asdict(m) for m in native_moves], serialized_anchor_movements=serialized_moves,
            movement_passed=all(m.within_limits for m in native_moves) and all(m["within_limits"] for m in serialized_moves),
            restoration_mask_assessment=asdict(assessment), geometry_subset_diagnostics=full_diagnostics,
            assessment_scope_revision="protein-request-view-plus-full-context-custody-v2",
            diagnostic_scope="chemistry/Ramachandran/chi12/sterics subset; full CLI postflight requires independent readback",
            environment=observed_environment_metrics(topology, native_before, native_after, observed),
            generated_rms_movement_angstrom=math.sqrt(sum(d*d for d in generated_distances)/len(generated_distances)),
            generated_max_movement_angstrom=max(generated_distances),
            generated_dat_sha256=digest(output / "generated.dat"),
            final_sha256=digest(output / "refined.pdb"))
        row["outcome_flags"] = []
        if not row["movement_passed"]:
            row["outcome_flags"].append("envelope_violation")
        if not assessment.geometry_passed or assessment.admission_failures:
            row["outcome_flags"].append("restoration_geometry_or_admission_failure")
        if full_diagnostics["new_or_changed_findings"]:
            row["outcome_flags"].append("introduced_geometry_subset_findings")
        row["restoration_mask_and_motion_passed"] = row["movement_passed"] and assessment.geometry_passed and not assessment.admission_failures
        row["geometry_subset_no_new_findings"] = not full_diagnostics["new_or_changed_findings"]
    except ChiralityError as exc:
        row.update(error_type=type(exc).__name__, reason=str(exc), outcome_flags=["global_chirality_failure"], unresolved_failure=False)
        (output / "failure.txt").write_text(traceback.format_exc())
    except (OpenMMException, GridNumericalError) as exc:
        row.update(error_type=type(exc).__name__, reason=str(exc), outcome_flags=["numerical_failure"], unresolved_failure=False)
        (output / "failure.txt").write_text(traceback.format_exc())
    except Exception as exc:
        row.update(error_type=type(exc).__name__, reason=str(exc), outcome_flags=["implementation_or_invariant_error"], unresolved_failure=True)
        (output / "failure.txt").write_text(traceback.format_exc())
    save(output / "result.json", row)
    if row.get("unresolved_failure"):
        raise RuntimeError("unresolved grid implementation/invariant error; stop, retain and do not retry")


def declare(preflight, output):
    output.mkdir(parents=True, exist_ok=False)
    parent_plan_path = preflight / "plan.json"
    parent = json.loads(parent_plan_path.read_text())
    resolved = verify_bindings(parent["source_code_and_input_sha256"], (parent_plan_path,))
    historical_bindings = {str(resolved[name]): expected for name, expected in parent["source_code_and_input_sha256"].items()}
    receipt_path = preflight / "independent-readback-v2/receipt.json"
    receipt = json.loads(receipt_path.read_text())
    if receipt["plan_sha256"] != digest(parent_plan_path):
        raise ValueError("preflight readback binds a different plan")
    files = [Path(__file__), Path(__file__).with_name("trace_private_minimization.py"),
             Path(__file__).with_name("audit_retained_full_context_controls.py"),
             Path(__file__).with_name("audit_anchor_calibration_preflight.py"), receipt_path, parent_plan_path,
             Path(__file__).with_name("audit_anchor_calibration_grid.py"),
             Path(__file__).with_name("run_mac_modeller_comparison.py"),
             Path(__file__).resolve().parents[1] / "tests/test_anchor_calibration_grid.py"]
    jobs, admitted, blocked = [], [], []
    for row in receipt["cases"]:
        entry = row["entry"]
        if not row["system_preflight_passed"] or not row.get("readback", {}).get("independent_particle_readback_passed"):
            blocked.append({"entry": entry, "reason": row.get("reason", row.get("readback_error"))})
            continue
        case, snapshot = preflight / "cases" / entry, preflight / f"execution/{entry}-seed-7"
        preparation_seconds = json.loads((snapshot / "execution.json").read_text())["elapsed_seconds"]
        if not 0 < preparation_seconds < 1800:
            raise ValueError("shared preparation exhausted the per-variant time budget")
        for name, expected in row["readback"]["system_and_initial_state_sha256"].items():
            if digest(snapshot / name) != expected:
                raise ValueError("preflight initial state changed")
        files += list(case.iterdir()) + [snapshot / name for name in (
            "base-system.xml", "temporary-axis.json", "initial-state.npz", "result.json", "execution.json")]
        admitted.append(entry)
        for cell in variants():
            case_id = f"{entry}-{cell['id']}"
            jobs.append({"case_id": case_id, "case": str(case), "seed": 7,
                "shared_preparation_seconds": preparation_seconds, "timeout_seconds": 1800-preparation_seconds,
                "command": [sys.executable, "-m", "scripts.anchor_calibration_grid", "worker", "--case", str(case),
                            "--snapshot", str(snapshot), "--cell", cell["id"],
                            "--output", str(output / f"execution/{case_id}-seed-7")]})
    if not 4 <= len(admitted) <= 6:
        raise ValueError("grid requires the agreed4–6admittedcases; no silent reduced panel")
    save(output / "plan.json", {"protocol": "anchor-calibration-grid-v2-scope-corrected", "python": sys.executable,
        "python_sha256": digest(Path(sys.executable)), "source_code_and_input_sha256": {
            **historical_bindings, **{str(p): digest(p) for p in files}},
        "jobs": jobs, "count": len(jobs), "admitted": admitted, "blocked": blocked,
        "common_initial_state_per_case": True, "all25_even_if_zero_passes": True, "restarts": 0,
        "phase_caps": [100, 75, 75, 1000], "new_proposals": 0, "new_source_repairs": 0,
        "contact_metric": "4.5A observed-heavy residue proximity, exclude covalent neighbors; descriptive only",
        "secondary_metric": "MDAnalysis intrasegment DSSP C3, auxiliary FF H (missing H donors disabled), >=6-residue segments, initial proposal reference; no interchain sheets",
        "resource_policy": "serial1800s/job,RSS8GiB,host-pressure,10GiB evidence", "public_qualification": False})
    freeze_code(output / "plan.json")


def summarize(root):
    plan = json.loads((root / "plan.json").read_text())
    rows = []
    for job in plan["jobs"]:
        directory = root / f"execution/{job['case_id']}-seed-{job['seed']}"
        path = directory / "result.json"
        execution_path = directory / "execution.json"
        result = json.loads(path.read_text()) if path.is_file() else None
        execution = json.loads(execution_path.read_text()) if execution_path.is_file() else None
        rows.append({"case_id": job["case_id"], "state": "not_launched" if not directory.exists() else
                     "reported" if result else "resource_stopped" if execution and execution["guard_stop_reason"] else "incomplete",
                     "result_path": str(path) if result else None,
                     "outcome_flags": result.get("outcome_flags", []) if result else [],
                     "restoration_mask_and_motion_passed": result.get("restoration_mask_and_motion_passed") if result else None,
                     "geometry_subset_no_new_findings": result.get("geometry_subset_no_new_findings", result.get("whole_complex_no_new_findings")) if result else None,
                     "guard_stop_reason": execution.get("guard_stop_reason") if execution else None,
                     "elapsed_seconds": execution.get("elapsed_seconds") if execution else None})
    save(root / "summary.json", {"declared": len(rows), "cases": rows, "blocked_inputs": plan["blocked"],
        "no_unresolved_implementation_errors": not any("implementation_or_invariant_error" in r["outcome_flags"] for r in rows),
        "all_cells_reported": all(r["state"] == "reported" for r in rows), "production_policy_selected": False})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("declare", "worker", "run"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--preflight", type=Path)
    parser.add_argument("--case", type=Path)
    parser.add_argument("--snapshot", type=Path)
    parser.add_argument("--cell")
    args = parser.parse_args()
    if args.mode == "declare":
        if args.preflight is None:
            parser.error("declare requires --preflight")
        declare(args.preflight.resolve(), args.output.resolve())
    elif args.mode == "worker":
        if args.case is None or args.snapshot is None or args.cell is None:
            parser.error("worker requires --case/--snapshot/--cell")
        worker(args.case.resolve(), args.snapshot.resolve(), cell_by_id(args.cell), args.output.resolve())
    else:
        from scripts.run_mac_modeller_comparison import execute
        try:
            execute(args.output.resolve())
        finally:
            summarize(args.output.resolve())
