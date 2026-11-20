"""Private research admission/returned-state guards; no repair or minimization.

Not a public completion capability or a production-default change. Covalent X–H
equilibria come only from the captured FF, never the sampled coordinates. OpenMM
hard constraints involving massless observed parents are explicitly unsupported.
"""
from __future__ import annotations

import math

import numpy as np
from openmm import HarmonicBondForce
from openmm.unit import dalton, nanometer
from scipy.spatial import cKDTree

from dvbfixer.diagnose.chemistry import _BOND_LEN_WARN_FRAC, check_stereocentre_geometry
from dvbfixer.diagnose.geometry import atom_identity
from dvbfixer.diagnose.report import findings_to_dict_list
from dvbfixer.diagnose.steric import clashes_python
from dvbfixer.ffutils.geometry import _COINCIDENT_TOL_NM, find_d_residues

REVISION = "private-anchor-refinement-guards-v1"
FORCE_SCALE_KJ_MOL_NM = 10.0


def _axis(topology, system, positions):
    atoms = list(topology.atoms())
    points = np.asarray(positions, dtype=float)
    if (system.getNumParticles() != len(atoms) or points.shape != (len(atoms), 3)
            or any(a.index != i or a.element is None for i, a in enumerate(atoms))):
        raise ValueError("guard requires a complete canonical FF particle axis")
    masses = np.array([system.getParticleMass(i).value_in_unit(dalton) for i in range(len(atoms))])
    if not np.isfinite(masses).all() or (masses < 0).any():
        raise ValueError("invalid FF particle masses")
    return atoms, points, masses


def xh_equilibria(topology, system):
    """Inventory unique actual covalent X–H edges and their FF equilibrium lengths."""
    atoms = list(topology.atoms())
    terms = {}
    for force in system.getForces():
        if isinstance(force, HarmonicBondForce):
            for i in range(force.getNumBonds()):
                first, second, length, _ = force.getBondParameters(i)
                terms.setdefault(frozenset((first, second)), []).append(float(length.value_in_unit(nanometer)))
    parents = {}
    for first, second in topology.bonds():
        if first.element is None or second.element is None:
            raise ValueError("unknown bonded element")
        if first.element.symbol == "H" and second.element.symbol != "H":
            parents.setdefault(first.index, []).append(second.index)
        elif second.element.symbol == "H" and first.element.symbol != "H":
            parents.setdefault(second.index, []).append(first.index)
    rows = []
    for atom in atoms:
        if atom.element is None or atom.element.symbol != "H":
            continue
        neighbors = parents.get(atom.index, [])
        if len(neighbors) != 1:
            raise ValueError("hydrogen requires exactly one supplied heavy parent")
        parent = neighbors[0]
        lengths = terms.get(frozenset((atom.index, parent)), [])
        if len(lengths) != 1 or not math.isfinite(lengths[0]) or lengths[0] <= 0:
            raise ValueError("covalent X–H requires one positive FF harmonic equilibrium")
        massless = [i for i in (atom.index, parent) if system.getParticleMass(i).value_in_unit(dalton) == 0]
        rows.append({"hydrogen_index": atom.index, "parent_index": parent,
            "hydrogen": atom_identity(atom).to_dict(), "parent": atom_identity(atoms[parent]).to_dict(),
            "equilibrium_nm": lengths[0], "massless_indices": massless,
            "openmm_hard_constraint_supported": not massless})
    return rows


def immutable_context_admission(topology, system, positions_nm):
    """Severe immutable heavy contacts cannot be repaired by moving generated atoms."""
    atoms, points, masses = _axis(topology, system, positions_nm)
    if not np.isfinite(points).all():
        return {"admitted": False, "reasons": ["nonfinite_input"], "immutable_severe_clashes": []}
    provenance = {}
    findings = clashes_python(topology, points*nanometer, include_water=True,
                              clash_warn_a=.6, clash_error_a=.9, provenance=provenance)
    immutable = {(a.residue.chain.id, str(a.residue.id), a.residue.insertionCode.strip(), a.name) for a in atoms
                 if masses[a.index] == 0 and a.element.symbol != "H"}
    severe = []
    for finding in findings_to_dict_list(findings):
        if finding["category"] != "clash" or finding["severity"] != "ERROR":
            continue
        identities = finding["details"]["atoms"]
        if all((a["chain"], a["resid"], a["icode"], a["atom"]) in immutable for a in identities):
            severe.append(finding)
    reasons = []
    if provenance.get("status") != "complete":
        reasons.append("unresolved_steric_connectivity")
    if severe:
        reasons.append("incompatible_immutable_context")
    return {"admitted": not reasons, "reasons": reasons, "steric": provenance,
            "immutable_severe_clashes": severe, "threshold_angstrom": .9,
            "policy": "report/refuse; no water deletion, source motion or chemical reinterpretation"}


def _rms(values):
    """Scaled Cartesian-component RMS, finite even for enormous saved forces."""
    if not values.size:
        return None
    maximum = float(np.max(np.abs(values)))
    return 0. if maximum == 0 else maximum*float(np.sqrt(np.mean((values/maximum)**2)))


def returned_state_guard(topology, system, initial_nm, positions_nm, forces, generated_indices):
    """Independent unconstrained Cartesian stationarity + H/identity checks.

    The active-component RMS uses ONLY positive-mass DOFs. Generated-heavy RMS
    is separately reported, not used to retroactively invent a per-mask threshold.
    Constraint-projected stationarity is deliberately unsupported in this version.
    """
    atoms, points, masses = _axis(topology, system, positions_nm)
    initial = np.asarray(initial_nm, dtype=float)
    gradient = np.asarray(forces, dtype=float)
    if initial.shape != points.shape or gradient.shape != points.shape:
        raise ValueError("returned-state axes differ")
    if system.getNumConstraints():
        raise ValueError("guard v1 does not interpret constraint-projected forces")
    generated = sorted(generated_indices)
    if any(i < 0 or i >= len(atoms) or atoms[i].element.symbol == "H" or masses[i] <= 0 for i in generated):
        raise ValueError("generated-heavy force axis differs from mobile FF ownership")
    if not np.isfinite(initial).all():
        raise ValueError("nonfinite initial observations")
    if not np.isfinite(points).all() or not np.isfinite(gradient).all():
        return {"passed": False, "reasons": ["nonfinite_returned_state"], "revision": REVISION}
    mobile = masses > 0
    if not mobile.any():
        raise ValueError("returned-state guard requires active particles")
    reasons = []
    fixed_exact = np.array_equal(initial[~mobile], points[~mobile])
    if not fixed_exact:
        reasons.append("native_fixed_observation_motion")
    active_rms = _rms(gradient[mobile])
    if active_rms > FORCE_SCALE_KJ_MOL_NM:
        reasons.append("nonstationary_active_cartesian_state")
    bonds = xh_equilibria(topology, system)
    bad_bonds = []
    for row in bonds:
        distance = float(np.linalg.norm(points[row["hydrogen_index"]]-points[row["parent_index"]]))
        relative = abs(distance-row["equilibrium_nm"])/row["equilibrium_nm"]
        if relative >= _BOND_LEN_WARN_FRAC-1.e-12:
            bad_bonds.append({**row, "distance_nm": distance, "relative_deviation": relative})
    if bad_bonds:
        reasons.append("covalent_xh_geometry_failure")
    collapses = []
    for first, second in sorted(cKDTree(points).query_pairs(_COINCIDENT_TOL_NM)):
        if "H" not in (atoms[first].element.symbol, atoms[second].element.symbol):
            continue
        distance = float(np.linalg.norm(points[first]-points[second]))
        if distance < _COINCIDENT_TOL_NM:
            collapses.append({"atoms": [atom_identity(atoms[i]).to_dict() for i in (first, second)],
                              "distance_nm": distance})
    if collapses:
        reasons.append("hydrogen_atom_collapse")
    d_ca = find_d_residues(topology, points*nanometer)
    if d_ca:
        reasons.append("global_d_ca")
    stereo = findings_to_dict_list(check_stereocentre_geometry(topology, points*nanometer))
    generated_residue_keys = {(atoms[i].residue.chain.id, str(atoms[i].residue.id),
                               atoms[i].residue.insertionCode.strip()) for i in generated}
    generated_stereo = [f for f in stereo if any((a["chain"], a["resid"], a["icode"])
                       in generated_residue_keys for a in f["details"]["atoms"])]
    if generated_stereo:
        reasons.append("generated_stereocentre_geometry")
    return {"passed": not reasons, "reasons": reasons, "revision": REVISION,
        "native_fixed_exact": fixed_exact, "active_cartesian_component_rms": active_rms,
        "active_particles": int(mobile.sum()), "total_particles": len(atoms),
        "massless_zeroed_all_component_rms": active_rms*math.sqrt(float(mobile.sum())/len(atoms)),
        "generated_heavy_component_rms": _rms(gradient[generated]),
        "force_scale_kj_mol_nm": FORCE_SCALE_KJ_MOL_NM,
        "stationarity_scope": "independent active-Cartesian test; not an OpenMM termination code or geometry pass",
        "xh_warning_relative_deviation": _BOND_LEN_WARN_FRAC,
        "xh_equilibrium_source": "captured HarmonicBondForce, not raw coordinates",
        "xh_bond_count": len(bonds), "xh_geometry_failures": bad_bonds,
        "openmm_xh_constraints_blocked_by_massless": sum(not b["openmm_hard_constraint_supported"] for b in bonds),
        "hydrogen_collapse_threshold_nm": _COINCIDENT_TOL_NM, "hydrogen_collapses": collapses,
        "d_ca": d_ca, "generated_stereocentre_findings": generated_stereo}
