"""Offline, geometry-only reconstruction for known isolated nonprotein components.

This application service intentionally does not produce force-field parameters.
Its public integration is typed Python: callers provide one exact component
instance and one exclusive authority object, then may atomically publish a
successful geometry bundle.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import re
import shutil
import tempfile
from dataclasses import asdict, is_dataclass, replace
from enum import Enum
from pathlib import Path
from typing import Any

import numpy as np

from dvbfixer.domain.nonprotein_reconstruction import (
    AtomCoordinate,
    AtomMovement,
    AuthorityProvenance,
    ChemicalAtom,
    ChemicalBond,
    ChemicalGraph,
    ComponentClass,
    CoordinateSource,
    ExternalLink,
    LocalPinnedCcdAuthority,
    MicrostateProvenance,
    ObservedAtom,
    OnlineEnrichmentAuthority,
    ParameterizationDecision,
    ParameterizationRequest,
    ReconstructionRequest,
    ReconstructionResult,
    ReconstructionStatus,
    UserMappedAuthority,
    ValidationFinding,
)
from dvbfixer.domain.parameterization import COMPLEX_COFACTORS, classify_parameterization
from dvbfixer.domain.structure_identity import ExactAtomRef

CCD_FORMAT = "dvbfixer-ccd-snapshot-v1"
MICROSTATE_POLICY = "declared-ph-range-v1"
MAX_AUTHORITY_BYTES = 5 * 1024 * 1024
MAX_ATOMS = 512
MAX_BONDS = 2048
MAX_STATES = 64
MAX_JSON_DEPTH = 24
MAX_TEXT_LENGTH = 256
_SAFE_ARTIFACT_NAME = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]{0,127}\Z")
_METALS = frozenset(
    {"LI", "NA", "K", "RB", "CS", "MG", "CA", "SR", "BA", "AL", "MN", "FE", "CO", "NI", "CU", "ZN", "CD", "HG"}
)
_CLASS_A_ELEMENTS = frozenset({"H", "B", "C", "N", "O", "F", "SI", "P", "S", "CL", "SE", "BR", "I"})


class _OutcomeError(Exception):
    def __init__(
        self,
        status: ReconstructionStatus,
        code: str,
        message: str,
        component_class: ComponentClass = ComponentClass.D,
        microstate: MicrostateProvenance | None = None,
    ) -> None:
        super().__init__(message)
        self.status = status
        self.finding = ValidationFinding(code, message)
        self.component_class = component_class
        self.microstate = microstate


def _json_value(value: Any) -> Any:
    if is_dataclass(value) and not isinstance(value, type):
        return _json_value(asdict(value))
    if isinstance(value, Enum):
        return value.value
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, tuple):
        return [_json_value(item) for item in value]
    if isinstance(value, list):
        return [_json_value(item) for item in value]
    if isinstance(value, dict):
        return {str(key): _json_value(item) for key, item in value.items()}
    return value


def _canonical_bytes(value: Any) -> bytes:
    return json.dumps(
        _json_value(value), sort_keys=True, separators=(",", ":"), ensure_ascii=False
    ).encode("utf-8")


def _digest(value: Any) -> str:
    return hashlib.sha256(_canonical_bytes(value)).hexdigest()


def chemical_graph_digest(graph: ChemicalGraph) -> str:
    """Return the canonical digest required by ``UserMappedAuthority``."""
    return _digest(graph)


def _json_depth(value: Any, depth: int = 0) -> int:
    if depth > MAX_JSON_DEPTH:
        return depth
    if isinstance(value, dict):
        return max((_json_depth(item, depth + 1) for item in value.values()), default=depth)
    if isinstance(value, list):
        return max((_json_depth(item, depth + 1) for item in value), default=depth)
    return depth


def _number(value: Any, field: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise _OutcomeError(ReconstructionStatus.INVALID_INPUT, "invalid-authority", f"{field} must be numeric")
    result = float(value)
    if not math.isfinite(result):
        raise _OutcomeError(ReconstructionStatus.INVALID_INPUT, "invalid-authority", f"{field} must be finite")
    return result


def _graph_from_record(component_id: str, state: Any) -> ChemicalGraph:
    if not isinstance(state, dict):
        raise _OutcomeError(ReconstructionStatus.INVALID_INPUT, "invalid-authority", "state must be an object")
    atoms_raw = state.get("atoms")
    bonds_raw = state.get("bonds")
    if not isinstance(atoms_raw, list) or not isinstance(bonds_raw, list):
        raise _OutcomeError(ReconstructionStatus.INVALID_INPUT, "invalid-authority", "state atoms and bonds must be arrays")
    if not 1 <= len(atoms_raw) <= MAX_ATOMS or len(bonds_raw) > MAX_BONDS:
        raise _OutcomeError(ReconstructionStatus.INVALID_INPUT, "authority-bounds", "authority atom or bond count exceeds limits")
    atoms: list[ChemicalAtom] = []
    for raw in atoms_raw:
        if not isinstance(raw, dict) or not isinstance(raw.get("ideal_position"), list) or len(raw["ideal_position"]) != 3:
            raise _OutcomeError(ReconstructionStatus.INVALID_INPUT, "invalid-authority", "every atom needs a three-value ideal_position")
        position = raw["ideal_position"]
        atoms.append(
            ChemicalAtom(
                name=str(raw.get("name", "")),
                element=str(raw.get("element", "")),
                formal_charge=int(raw.get("formal_charge", 0)),
                ideal_position=(
                    _number(position[0], "ideal_position"),
                    _number(position[1], "ideal_position"),
                    _number(position[2], "ideal_position"),
                ),
                aromatic=bool(raw.get("aromatic", False)),
                stereo=str(raw["stereo"]) if raw.get("stereo") is not None else None,
            )
        )
    bonds: list[ChemicalBond] = []
    for raw in bonds_raw:
        if not isinstance(raw, dict):
            raise _OutcomeError(ReconstructionStatus.INVALID_INPUT, "invalid-authority", "every bond must be an object")
        bonds.append(
            ChemicalBond(
                atom1=str(raw.get("atom1", "")),
                atom2=str(raw.get("atom2", "")),
                order=_number(raw.get("order"), "bond order"),
                ideal_length=_number(raw.get("ideal_length"), "ideal bond length"),
                aromatic=bool(raw.get("aromatic", False)),
            )
        )
    return ChemicalGraph(
        component_id=component_id,
        state_id=str(state.get("id", "")),
        atoms=tuple(atoms),
        bonds=tuple(bonds),
        net_charge=int(state.get("net_charge", 0)),
        ph_min=_number(state["ph_min"], "ph_min") if state.get("ph_min") is not None else None,
        ph_max=_number(state["ph_max"], "ph_max") if state.get("ph_max") is not None else None,
        priority=int(state.get("priority", 0)),
        redox_state=str(state["redox_state"]) if state.get("redox_state") is not None else None,
        coordination_model=(
            str(state["coordination_model"])
            if state.get("coordination_model") is not None
            else None
        ),
        component_count=int(state.get("component_count", 1)),
    )


def _load_local_authority(authority: LocalPinnedCcdAuthority) -> tuple[list[ChemicalGraph], AuthorityProvenance]:
    path = authority.path
    try:
        if path.is_symlink() or not path.is_file():
            raise _OutcomeError(ReconstructionStatus.INVALID_INPUT, "invalid-authority-path", "local CCD snapshot must be a regular non-symlink file")
        size = path.stat().st_size
        if size > MAX_AUTHORITY_BYTES:
            raise _OutcomeError(ReconstructionStatus.INVALID_INPUT, "authority-bounds", "local CCD snapshot exceeds the 5 MiB limit")
        content = path.read_bytes()
        if len(content) > MAX_AUTHORITY_BYTES:
            raise _OutcomeError(ReconstructionStatus.INVALID_INPUT, "authority-bounds", "local CCD snapshot exceeds the 5 MiB limit")
    except _OutcomeError:
        raise
    except OSError as exc:
        raise _OutcomeError(ReconstructionStatus.FAILED, "authority-read-failed", str(exc)) from exc
    actual_digest = hashlib.sha256(content).hexdigest()
    if not re.fullmatch(r"[0-9a-f]{64}", authority.sha256) or actual_digest != authority.sha256:
        raise _OutcomeError(ReconstructionStatus.INVALID_INPUT, "authority-digest-mismatch", "local CCD snapshot SHA-256 does not match its required pin")
    try:
        data = json.loads(content)
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise _OutcomeError(ReconstructionStatus.INVALID_INPUT, "invalid-authority", f"invalid CCD snapshot JSON: {exc}") from exc
    if _json_depth(data) > MAX_JSON_DEPTH or not isinstance(data, dict):
        raise _OutcomeError(ReconstructionStatus.INVALID_INPUT, "authority-bounds", "CCD snapshot nesting exceeds limits")
    if data.get("format") != CCD_FORMAT or not isinstance(data.get("version"), str):
        raise _OutcomeError(ReconstructionStatus.INVALID_INPUT, "invalid-authority", f"CCD snapshot must declare format {CCD_FORMAT!r} and a version")
    components = data.get("components")
    if not isinstance(components, dict) or authority.component_id not in components:
        raise _OutcomeError(ReconstructionStatus.UNSUPPORTED, "component-not-in-snapshot", f"component {authority.component_id!r} is absent from the pinned snapshot")
    component = components[authority.component_id]
    if not isinstance(component, dict) or not isinstance(component.get("states"), list):
        raise _OutcomeError(ReconstructionStatus.INVALID_INPUT, "invalid-authority", "component states must be an array")
    states = component["states"]
    if not 1 <= len(states) <= MAX_STATES:
        raise _OutcomeError(ReconstructionStatus.INVALID_INPUT, "authority-bounds", "authority state count exceeds limits")
    return (
        [_graph_from_record(authority.component_id, state) for state in states],
        AuthorityProvenance(authority.mode, str(path), actual_digest, data["version"]),
    )


def _validate_graph(graph: ChemicalGraph) -> None:
    if not graph.component_id or not graph.state_id:
        raise _OutcomeError(ReconstructionStatus.INVALID_INPUT, "incomplete-graph", "component and state identifiers are required")
    if not 1 <= len(graph.atoms) <= MAX_ATOMS or len(graph.bonds) > MAX_BONDS:
        raise _OutcomeError(ReconstructionStatus.INVALID_INPUT, "authority-bounds", "graph atom or bond count exceeds limits")
    names = [atom.name for atom in graph.atoms]
    if any(not name or len(name) > MAX_TEXT_LENGTH for name in names) or len(names) != len(set(names)):
        raise _OutcomeError(ReconstructionStatus.AMBIGUOUS, "ambiguous-atom-names", "authority atom names must be nonempty and unique")
    if graph.component_count < 1 or graph.component_count > MAX_ATOMS:
        raise _OutcomeError(ReconstructionStatus.INVALID_INPUT, "authority-bounds", "component count is outside supported bounds")
    if (
        (graph.ph_min is None) != (graph.ph_max is None)
        or graph.ph_min is not None
        and graph.ph_max is not None
        and graph.ph_min > graph.ph_max
    ):
        raise _OutcomeError(ReconstructionStatus.INVALID_INPUT, "invalid-ph-range", "microstate pH bounds must be both absent or an ordered range")
    if sum(atom.formal_charge for atom in graph.atoms) != graph.net_charge:
        raise _OutcomeError(ReconstructionStatus.INVALID_INPUT, "charge-mismatch", "atom formal charges do not sum to the declared net charge")
    known = set(names)
    seen_bonds: set[frozenset[str]] = set()
    for atom in graph.atoms:
        if not atom.element or any(not math.isfinite(value) for value in atom.ideal_position):
            raise _OutcomeError(ReconstructionStatus.INVALID_INPUT, "incomplete-graph", f"atom {atom.name!r} lacks element or finite ideal coordinates")
    for bond in graph.bonds:
        pair = frozenset((bond.atom1, bond.atom2))
        if bond.atom1 == bond.atom2 or not {bond.atom1, bond.atom2} <= known or pair in seen_bonds:
            raise _OutcomeError(ReconstructionStatus.INVALID_INPUT, "invalid-bond-graph", "bond endpoints must be distinct, known, and unique")
        if bond.order <= 0 or bond.order > 3 or bond.ideal_length <= 0:
            raise _OutcomeError(ReconstructionStatus.INVALID_INPUT, "invalid-bond-graph", "bond order and ideal length must be positive")
        seen_bonds.add(pair)
    if len(graph.atoms) > 1:
        adjacency: dict[str, set[str]] = {name: set() for name in names}
        for bond in graph.bonds:
            adjacency[bond.atom1].add(bond.atom2)
            adjacency[bond.atom2].add(bond.atom1)
        visited: set[str] = set()
        pending = [names[0]]
        while pending:
            current = pending.pop()
            if current in visited:
                continue
            visited.add(current)
            pending.extend(adjacency[current] - visited)
        if visited != known:
            raise _OutcomeError(ReconstructionStatus.UNSUPPORTED, "disconnected-graph", "Class A requires one connected molecule")


def _resolve_authority(request: ReconstructionRequest) -> tuple[list[ChemicalGraph], AuthorityProvenance]:
    authority = request.authority
    if isinstance(authority, OnlineEnrichmentAuthority):
        raise _OutcomeError(ReconstructionStatus.UNSUPPORTED, "online-enrichment-disabled", "online enrichment has no approved adapter and is explicitly refused")
    if isinstance(authority, LocalPinnedCcdAuthority):
        return _load_local_authority(authority)
    if isinstance(authority, UserMappedAuthority):
        if not authority.source_label or len(authority.source_label) > MAX_TEXT_LENGTH or authority.source_digest != chemical_graph_digest(authority.graph):
            raise _OutcomeError(ReconstructionStatus.INVALID_INPUT, "authority-digest-mismatch", "user-mapped graph requires its exact canonical digest and source label")
        return [authority.graph], AuthorityProvenance(authority.mode, authority.source_label, authority.source_digest, None)
    raise _OutcomeError(ReconstructionStatus.INVALID_INPUT, "authority-mode", "exactly one recognized authority mode is required")


def _select_microstate(request: ReconstructionRequest, graphs: list[ChemicalGraph]) -> tuple[ChemicalGraph, MicrostateProvenance]:
    if not math.isfinite(request.ph) or not 0.0 <= request.ph <= 14.0:
        raise _OutcomeError(ReconstructionStatus.INVALID_INPUT, "invalid-ph", "explicit pH must be finite and between 0 and 14")
    by_id = {graph.state_id: graph for graph in graphs}
    if len(by_id) != len(graphs):
        raise _OutcomeError(ReconstructionStatus.AMBIGUOUS, "duplicate-state", "authority state identifiers are not unique")
    if request.microstate_override is not None:
        selected = by_id.get(request.microstate_override)
        if selected is None:
            provenance = MicrostateProvenance(
                MICROSTATE_POLICY,
                request.ph,
                None,
                tuple(sorted(by_id)),
                tuple((state_id, by_id[state_id].priority) for state_id in sorted(by_id)),
                True,
                "requested override is absent",
            )
            raise _OutcomeError(ReconstructionStatus.INVALID_INPUT, "unknown-state-override", "microstate override is absent from the selected authority", microstate=provenance)
        provenance = MicrostateProvenance(MICROSTATE_POLICY, request.ph, selected.state_id, tuple(sorted(by_id)), tuple((state_id, by_id[state_id].priority) for state_id in sorted(by_id)), True, "explicit user override")
        return selected, provenance
    if isinstance(request.authority, UserMappedAuthority):
        selected = graphs[0]
        provenance = MicrostateProvenance(MICROSTATE_POLICY, request.ph, selected.state_id, (selected.state_id,), ((selected.state_id, selected.priority),), False, "user graph is a locked microspecies")
        return selected, provenance
    eligible = [
        graph for graph in graphs
        if graph.ph_min is not None and graph.ph_max is not None and graph.ph_min <= request.ph <= graph.ph_max
    ]
    ranked = tuple(graph.state_id for graph in sorted(eligible, key=lambda item: (item.priority, item.state_id)))
    candidate_scores = tuple(
        (graph.state_id, graph.priority if graph in eligible else None)
        for graph in sorted(graphs, key=lambda item: item.state_id)
    )
    if not eligible:
        provenance = MicrostateProvenance(MICROSTATE_POLICY, request.ph, None, (), candidate_scores, False, "no declared pH range covers the request")
        raise _OutcomeError(ReconstructionStatus.UNSUPPORTED, "microstate-unresolved", "no declared authority microstate covers the explicit pH", microstate=provenance)
    best_priority = min(graph.priority for graph in eligible)
    best = [graph for graph in eligible if graph.priority == best_priority]
    if len(best) != 1:
        provenance = MicrostateProvenance(MICROSTATE_POLICY, request.ph, None, ranked, candidate_scores, False, "equal-priority states overlap at pH")
        raise _OutcomeError(ReconstructionStatus.AMBIGUOUS, "microstate-ambiguous", "equally authoritative microstates overlap at the explicit pH", microstate=provenance)
    selected = best[0]
    return selected, MicrostateProvenance(MICROSTATE_POLICY, request.ph, selected.state_id, ranked, candidate_scores, False, "lowest declared priority covering pH")


def _classify(request: ReconstructionRequest, graph: ChemicalGraph) -> ComponentClass:
    elements = {atom.element.upper() for atom in graph.atoms}
    if (
        any(atom.element.upper() in _METALS for atom in graph.atoms)
        or graph.redox_state is not None
        or graph.coordination_model is not None
        or graph.component_count > 1
        or graph.component_id.upper() in COMPLEX_COFACTORS
    ):
        actual = ComponentClass.C
    elif request.external_links:
        actual = ComponentClass.B
    elif "C" not in elements or not elements <= _CLASS_A_ELEMENTS:
        actual = ComponentClass.D
    else:
        actual = ComponentClass.A
    if request.requested_class not in (ComponentClass.AUTO, actual):
        raise _OutcomeError(ReconstructionStatus.UNSUPPORTED, "class-mismatch", f"evidence classifies the component as Class {actual.value}, not requested Class {request.requested_class.value}", actual)
    if actual is ComponentClass.B:
        raise _OutcomeError(ReconstructionStatus.UNSUPPORTED, "class-b-unsupported", "linked carbohydrate, PTM, and covalent-ligand reconstruction has no approved linkage-aware backend", actual)
    if actual is ComponentClass.C:
        raise _OutcomeError(ReconstructionStatus.UNSUPPORTED, "class-c-unsupported", "metal, redox, cofactor, and multi-component reconstruction has no approved state/coordination backend", actual)
    if actual is ComponentClass.D:
        raise _OutcomeError(ReconstructionStatus.UNSUPPORTED, "class-d-unsupported", "component is not a complete supported isolated organic Class A graph", actual)
    return actual


def _fit_coordinates(request: ReconstructionRequest, graph: ChemicalGraph) -> tuple[tuple[AtomCoordinate, ...], tuple[AtomMovement, ...], float]:
    if len(request.observed_atoms) > MAX_ATOMS:
        raise _OutcomeError(ReconstructionStatus.INVALID_INPUT, "request-bounds", "observed atom count exceeds limits")
    observed_by_name: dict[str, ObservedAtom] = {}
    graph_by_name = {atom.name: atom for atom in graph.atoms}
    for observed in request.observed_atoms:
        if observed.identity.component != request.component:
            raise _OutcomeError(ReconstructionStatus.INVALID_INPUT, "identity-mismatch", "every observed atom must belong to the exact requested component instance")
        if observed.identity.name in observed_by_name:
            raise _OutcomeError(ReconstructionStatus.AMBIGUOUS, "duplicate-observed-atom", f"observed atom {observed.identity.name!r} is not unique")
        authority_atom = graph_by_name.get(observed.identity.name)
        if authority_atom is None or authority_atom.element.upper() != observed.element.upper():
            raise _OutcomeError(ReconstructionStatus.INVALID_INPUT, "atom-map-mismatch", f"observed atom {observed.identity.name!r} does not map exactly by name and element")
        if any(not math.isfinite(value) for value in observed.position):
            raise _OutcomeError(ReconstructionStatus.INVALID_INPUT, "invalid-coordinate", "observed coordinates must be finite")
        observed_by_name[observed.identity.name] = observed
    missing = [atom for atom in graph.atoms if atom.name not in observed_by_name]
    rotation = np.eye(3)
    translation = np.zeros(3)
    anchor_rmsd = 0.0
    if missing:
        if len(observed_by_name) < 3:
            raise _OutcomeError(ReconstructionStatus.UNSUPPORTED, "insufficient-anchors", "reconstruction requires at least three observed anchors")
        ordered = sorted(observed_by_name)
        ideal = np.asarray([graph_by_name[name].ideal_position for name in ordered], dtype=float)
        target = np.asarray([observed_by_name[name].position for name in ordered], dtype=float)
        ideal_centered = ideal - ideal.mean(axis=0)
        if np.linalg.matrix_rank(ideal_centered) < 2:
            raise _OutcomeError(ReconstructionStatus.AMBIGUOUS, "collinear-anchors", "observed atom mapping does not define a unique rigid frame")
        target_centered = target - target.mean(axis=0)
        u_matrix, _singular, vt_matrix = np.linalg.svd(ideal_centered.T @ target_centered)
        rotation = vt_matrix.T @ u_matrix.T
        if np.linalg.det(rotation) < 0:
            vt_matrix[-1, :] *= -1
            rotation = vt_matrix.T @ u_matrix.T
        translation = target.mean(axis=0) - ideal.mean(axis=0) @ rotation.T
        fitted = ideal @ rotation.T + translation
        anchor_rmsd = float(np.sqrt(np.mean(np.sum((fitted - target) ** 2, axis=1))))
        if anchor_rmsd > request.geometry_policy.max_anchor_rmsd:
            raise _OutcomeError(ReconstructionStatus.FAILED, "anchor-rmsd", f"authority geometry anchor RMSD {anchor_rmsd:.3f} A exceeds {request.geometry_policy.max_anchor_rmsd:.3f} A")
    coordinates: list[AtomCoordinate] = []
    movements: list[AtomMovement] = []
    for atom in graph.atoms:
        existing = observed_by_name.get(atom.name)
        identity = ExactAtomRef(request.component, atom.name)
        if existing is not None:
            position = existing.position
            source = CoordinateSource.OBSERVED
            movement = 0.0
            movements.append(AtomMovement(identity, movement))
        else:
            generated = np.asarray(atom.ideal_position) @ rotation.T + translation
            position = (float(generated[0]), float(generated[1]), float(generated[2]))
            source = CoordinateSource.SOURCE_TEMPLATE
        coordinates.append(AtomCoordinate(identity, atom.element, position, source, atom.name))
    return tuple(coordinates), tuple(movements), anchor_rmsd


def _validate_geometry(request: ReconstructionRequest, graph: ChemicalGraph, coordinates: tuple[AtomCoordinate, ...]) -> tuple[ValidationFinding, ...]:
    positions = {atom.authority_atom: np.asarray(atom.position) for atom in coordinates}
    bonded = {frozenset((bond.atom1, bond.atom2)) for bond in graph.bonds}
    findings: list[ValidationFinding] = []
    for bond in graph.bonds:
        distance = float(np.linalg.norm(positions[bond.atom1] - positions[bond.atom2]))
        if abs(distance - bond.ideal_length) > request.geometry_policy.max_bond_deviation:
            findings.append(ValidationFinding("bond-geometry", f"bond {bond.atom1}-{bond.atom2} length {distance:.3f} A differs from authority ideal {bond.ideal_length:.3f} A", (bond.atom1, bond.atom2)))
    adjacency: dict[str, list[str]] = {atom.name: [] for atom in graph.atoms}
    for bond in graph.bonds:
        adjacency[bond.atom1].append(bond.atom2)
        adjacency[bond.atom2].append(bond.atom1)
    ideal_positions = {atom.name: np.asarray(atom.ideal_position) for atom in graph.atoms}
    for atom in graph.atoms:
        neighbors = sorted(adjacency[atom.name])
        if atom.stereo is None or len(neighbors) < 3:
            continue
        first, second, third = neighbors[:3]
        ideal_volume = float(
            np.dot(
                np.cross(ideal_positions[first] - ideal_positions[atom.name], ideal_positions[second] - ideal_positions[atom.name]),
                ideal_positions[third] - ideal_positions[atom.name],
            )
        )
        output_volume = float(
            np.dot(
                np.cross(positions[first] - positions[atom.name], positions[second] - positions[atom.name]),
                positions[third] - positions[atom.name],
            )
        )
        if abs(ideal_volume) < 1e-8 or abs(output_volume) < 1e-8 or ideal_volume * output_volume < 0:
            findings.append(ValidationFinding("stereochemistry", f"stereocentre {atom.name} is degenerate or inverted relative to authority state {atom.stereo}", (atom.name, first, second, third)))
    names = sorted(positions)
    for index, atom1 in enumerate(names):
        for atom2 in names[index + 1 :]:
            if frozenset((atom1, atom2)) in bonded:
                continue
            distance = float(np.linalg.norm(positions[atom1] - positions[atom2]))
            if distance < request.geometry_policy.min_nonbonded_distance:
                findings.append(ValidationFinding("internal-clash", f"nonbonded atoms {atom1}-{atom2} are {distance:.3f} A apart", (atom1, atom2)))
    return tuple(findings)


def _preflight_request(request: ReconstructionRequest) -> None:
    identity_fields = (
        request.component.chain_id,
        request.component.sequence_number,
        request.component.insertion_code,
        request.component.alternate_location,
    )
    if any(len(value) > MAX_TEXT_LENGTH for value in identity_fields):
        raise _OutcomeError(ReconstructionStatus.INVALID_INPUT, "request-bounds", "component identity field exceeds limits")
    if len(request.observed_atoms) > MAX_ATOMS or len(request.external_links) > MAX_BONDS:
        raise _OutcomeError(ReconstructionStatus.INVALID_INPUT, "request-bounds", "request atom or link count exceeds limits")
    policy = request.geometry_policy
    if (
        policy.max_anchor_rmsd < 0
        or policy.max_bond_deviation < 0
        or policy.min_nonbonded_distance < 0
        or policy.max_observed_movement != 0
    ):
        raise _OutcomeError(ReconstructionStatus.INVALID_INPUT, "geometry-policy", "V1 requires nonnegative tolerances and fixed observed atoms")
    if isinstance(request.authority, UserMappedAuthority):
        _validate_graph(request.authority.graph)


def _geometry_digest(
    component: Any,
    graph: ChemicalGraph,
    links: tuple[ExternalLink, ...],
    coordinates: tuple[AtomCoordinate, ...],
) -> str:
    return _digest(
        {
            "component": component,
            "graph": graph,
            "external_links": links,
            "coordinates": coordinates,
        }
    )


def reconstruct_nonprotein(request: ReconstructionRequest) -> ReconstructionResult:
    """Reconstruct one Class A component without publishing or parameterizing it."""
    input_digest: str | None = None
    provenance: AuthorityProvenance | None = None
    microstate: MicrostateProvenance | None = None
    graph: ChemicalGraph | None = None
    try:
        _preflight_request(request)
        input_digest = _digest(request)
        graphs, provenance = _resolve_authority(request)
        for graph in graphs:
            _validate_graph(graph)
        graph, microstate = _select_microstate(request, graphs)
        component_class = _classify(request, graph)
        coordinates, movements, anchor_rmsd = _fit_coordinates(request, graph)
        findings = _validate_geometry(request, graph, coordinates)
        if findings:
            return ReconstructionResult(
                status=ReconstructionStatus.FAILED,
                component_class=component_class,
                component=request.component,
                resolved_graph=graph,
                external_links=request.external_links,
                findings=findings, authority=provenance, microstate=microstate,
                input_digest=input_digest,
            )
        output_digest = _geometry_digest(request.component, graph, request.external_links, coordinates)
        success_finding = ValidationFinding("geometry-validated", f"fixed-anchor geometry passed; fitted anchor RMSD {anchor_rmsd:.3f} A")
        return ReconstructionResult(
            status=ReconstructionStatus.SUCCEEDED,
            component_class=component_class,
            component=request.component,
            resolved_graph=graph,
            external_links=request.external_links,
            coordinates=coordinates,
            movements=movements,
            findings=(success_finding,),
            authority=provenance,
            microstate=microstate,
            input_digest=input_digest,
            output_digest=output_digest,
            geometry_approved=False,
            md_ready=False,
        )
    except _OutcomeError as exc:
        return ReconstructionResult(
            status=exc.status,
            component_class=exc.component_class,
            component=request.component,
            resolved_graph=graph,
            external_links=request.external_links,
            findings=(exc.finding,),
            authority=provenance,
            microstate=exc.microstate or microstate,
            input_digest=input_digest,
        )
    except Exception as exc:
        return ReconstructionResult(
            status=ReconstructionStatus.FAILED,
            component_class=ComponentClass.D,
            component=request.component,
            findings=(ValidationFinding("internal-error", f"{type(exc).__name__}: {exc}"),),
            authority=provenance,
            microstate=microstate,
            input_digest=input_digest,
        )


def decide_parameterization(request: ParameterizationRequest) -> ParameterizationDecision:
    """Choose an MD route without running it or claiming MD readiness."""
    geometry = request.geometry
    if geometry.status is not ReconstructionStatus.SUCCEEDED or not geometry.output_digest:
        return ParameterizationDecision(ReconstructionStatus.INVALID_INPUT, None, "parameterization requires a successful geometry result", geometry.output_digest, request.force_field_family, request.policy)
    if geometry.resolved_graph is None or geometry.output_digest != _geometry_digest(
        geometry.component,
        geometry.resolved_graph,
        geometry.external_links,
        geometry.coordinates,
    ):
        return ParameterizationDecision(ReconstructionStatus.INVALID_INPUT, None, "geometry digest does not match its graph, links, identity, and coordinates", geometry.output_digest, request.force_field_family, request.policy)
    if not request.geometry_approved:
        return ParameterizationDecision(ReconstructionStatus.UNSUPPORTED, None, "geometry requires an independent approval before MD routing", geometry.output_digest, request.force_field_family, request.policy)
    if geometry.component_class is not ComponentClass.A and not request.explicit_strip:
        return ParameterizationDecision(ReconstructionStatus.UNSUPPORTED, None, "only approved Class A geometry can enter the isolated-ligand route", geometry.output_digest, request.force_field_family, request.policy)
    route = classify_parameterization(
        "",
        exact_native_match=request.exact_native_match,
        user_template=request.validated_user_template,
        strip=request.explicit_strip,
    )
    if route.value == "gaff-candidate" and not request.allow_isolated_ligand_candidate:
        return ParameterizationDecision(ReconstructionStatus.UNSUPPORTED, None, "no explicit compatible parameterization evidence or isolated-ligand candidate policy", geometry.output_digest, request.force_field_family, request.policy)
    return ParameterizationDecision(ReconstructionStatus.SUCCEEDED, route, "route candidate selected; MD parameters and system validation have not run", geometry.output_digest, request.force_field_family, request.policy, md_ready=False)


def _write_json(path: Path, value: Any) -> None:
    with path.open("wb") as handle:
        handle.write(_canonical_bytes(value) + b"\n")
        handle.flush()
        os.fsync(handle.fileno())


def publish_reconstruction(result: ReconstructionResult, root: Path, artifact_name: str) -> Path:
    """Atomically publish a complete successful bundle as one directory rename."""
    if result.status is not ReconstructionStatus.SUCCEEDED or not result.output_digest:
        raise ValueError("only a successful reconstruction result can be published")
    if result.resolved_graph is None or result.output_digest != _geometry_digest(
        result.component,
        result.resolved_graph,
        result.external_links,
        result.coordinates,
    ):
        raise ValueError("reconstruction output digest does not match its payload")
    if not _SAFE_ARTIFACT_NAME.fullmatch(artifact_name):
        raise ValueError("artifact_name must be a bounded portable name")
    if root.is_symlink():
        raise ValueError("publication root must be an existing non-symlink directory")
    root = root.resolve(strict=True)
    if not root.is_dir():
        raise ValueError("publication root must be an existing non-symlink directory")
    destination = root / artifact_name
    if destination.exists() or destination.is_symlink():
        raise FileExistsError(f"publication destination already exists: {destination}")
    staging = Path(tempfile.mkdtemp(prefix=".dvbfixer-reconstruction-", dir=root))
    committed = False
    try:
        _write_json(staging / "geometry.json", {"output_digest": result.output_digest, "coordinates": result.coordinates})
        _write_json(staging / "provenance.json", result)
        directory_fd = os.open(staging, os.O_RDONLY)
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
        os.replace(staging, destination)
        committed = True
        root_fd = os.open(root, os.O_RDONLY)
        try:
            os.fsync(root_fd)
        finally:
            os.close(root_fd)
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        if committed:
            shutil.rmtree(destination, ignore_errors=True)
        raise
    return destination


def reconstruct_and_publish(request: ReconstructionRequest, root: Path, artifact_name: str) -> tuple[ReconstructionResult, Path | None]:
    """Run geometry reconstruction and publish only a complete success."""
    result = reconstruct_nonprotein(request)
    if result.status is not ReconstructionStatus.SUCCEEDED:
        return result, None
    published = publish_reconstruction(result, root, artifact_name)
    return replace(result, geometry_approved=False), published
