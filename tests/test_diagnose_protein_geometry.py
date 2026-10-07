"""Focused contracts for policy-free protein geometry evidence."""

from __future__ import annotations

import json
from typing import Any

import pytest

from dvbfixer.diagnose.geometry import (
    AtomIdentity,
    BoundaryIdentity,
    GeometryMeasurement,
    ResidueIdentity,
)
from dvbfixer.diagnose.protein_geometry import (
    GeometryLookup,
    MeasurementResult,
    ProteinAtom,
    ProteinResidue,
    ProteinStructure,
    backbone_torsions,
    canonical_heavy_bonds,
    connectivity_measurements,
    peptide_geometry,
    ramachandran_class,
    rotamer_reference_result,
    sidechain_torsions,
    steric_measurements,
)


def _atom(
    residue: ResidueIdentity,
    name: str,
    position: tuple[float, float, float],
    element: str = "C",
) -> ProteinAtom:
    return ProteinAtom(
        AtomIdentity(residue.chain, residue.resid, residue.icode, name),
        position,
        element,
    )


def _structure(
    residues: tuple[tuple[ResidueIdentity, str], ...],
    atoms: tuple[ProteinAtom, ...],
    bonds: tuple[tuple[AtomIdentity, AtomIdentity], ...] = (),
) -> ProteinStructure:
    return ProteinStructure(
        tuple(ProteinResidue(identity, name) for identity, name in residues),
        atoms,
        frozenset(frozenset(pair) for pair in bonds),
    )


def test_measurement_result_requires_consistent_metadata() -> None:
    atom = AtomIdentity("A", "1", "", "CA")
    measurement = GeometryMeasurement(
        "example.v1",
        (atom,),
        1.0,
        "angstrom",
        metadata={"source": "one"},
    )

    with pytest.raises(ValueError, match="metadata must agree"):
        MeasurementResult(
            "example.v1",
            (atom,),
            measurement,
            metadata={"source": "two"},
        )


def test_lookup_and_reference_failures_are_explicit() -> None:
    residue = ResidueIdentity("A", "1")
    first = AtomIdentity("A", "1", "", "CA")
    second = AtomIdentity("A", "1", "", "CB")
    third = AtomIdentity("A", "1", "", "CG")

    missing = GeometryLookup(_structure(((residue, "ALA"),), ())).measure(
        "example.v1", (first, second), "distance"
    )
    assert missing.reason == "missing_atom"

    duplicate_residue = GeometryLookup(
        _structure(
            ((residue, "ALA"), (residue, "ALA")),
            (_atom(residue, "CA", (0.0, 0.0, 0.0)),),
        )
    ).measure("example.v1", (first, second), "distance")
    assert duplicate_residue.reason == "duplicate_residue"

    duplicate_atom = _atom(residue, "CA", (0.0, 0.0, 0.0))
    duplicate = GeometryLookup(
        _structure(((residue, "ALA"),), (duplicate_atom, duplicate_atom))
    ).measure("example.v1", (first, second), "distance")
    assert duplicate.reason == "duplicate_atom"

    non_finite = GeometryLookup(
        _structure(
            ((residue, "ALA"),),
            (
                _atom(residue, "CA", (float("nan"), 0.0, 0.0)),
                _atom(residue, "CB", (1.0, 0.0, 0.0)),
            ),
        )
    ).measure("example.v1", (first, second), "distance")
    assert non_finite.reason == "non_finite_coordinates"

    degenerate = GeometryLookup(
        _structure(
            ((residue, "ALA"),),
            (
                _atom(residue, "CA", (0.0, 0.0, 0.0)),
                _atom(residue, "CB", (0.0, 0.0, 0.0)),
                _atom(residue, "CG", (1.0, 0.0, 0.0)),
            ),
        )
    ).measure("example.v1", (first, second, third), "angle")
    assert degenerate.reason == "degenerate_geometry"

    unsupported = sidechain_torsions(_structure(((residue, "MSE"),), ()), residue)[0]
    assert unsupported.reason == "unsupported_residue"
    assert rotamer_reference_result(residue, "ALA").reason == "unsupported_reference_data"


def test_ramachandran_classes_and_pre_pro_torsions() -> None:
    assert ramachandran_class("ALA") == "general"
    assert ramachandran_class("ALA", "PRO") == "pre-pro"
    assert ramachandran_class("GLY", "PRO") == "gly"
    assert ramachandran_class("PRO", "PRO") == "pro"
    assert ramachandran_class("MSE") is None

    previous = ResidueIdentity("D", "81")
    current = ResidueIdentity("D", "82", "A")
    following = ResidueIdentity("D", "83")
    previous_c = _atom(previous, "C", (-1.0, 0.0, 1.0))
    current_n = _atom(current, "N", (0.0, 0.0, 0.0), "N")
    current_ca = _atom(current, "CA", (1.0, 0.0, 0.0))
    current_c = _atom(current, "C", (1.0, 1.0, 0.0))
    following_n = _atom(following, "N", (2.0, 1.0, 1.0), "N")
    structure = _structure(
        ((previous, "ALA"), (current, "ALA"), (following, "PRO")),
        (previous_c, current_n, current_ca, current_c, following_n),
        (
            (previous_c.identity, current_n.identity),
            (current_c.identity, following_n.identity),
        ),
    )

    phi, psi = backbone_torsions(
        structure,
        current,
        (BoundaryIdentity(previous, current), BoundaryIdentity(current, following)),
    )

    assert phi.measurement is not None
    assert psi.measurement is not None
    assert phi.metadata["ramachandran_class"] == "pre-pro"
    assert phi.atoms[1].residue == current
    serialized_atoms = phi.to_dict()["atoms"]
    assert isinstance(serialized_atoms, list)
    assert serialized_atoms[1] == {
        "chain": "D",
        "resid": "82",
        "icode": "A",
        "atom": "N",
    }
    assert '"chain": "D"' in json.dumps(phi.to_dict(), sort_keys=True)


def test_pre_pro_class_is_undefined_when_declared_neighbor_is_missing() -> None:
    previous = ResidueIdentity("A", "1")
    current = ResidueIdentity("A", "2")
    missing_following = ResidueIdentity("A", "3")
    previous_c = _atom(previous, "C", (-1.0, 0.0, 1.0))
    current_n = _atom(current, "N", (0.0, 0.0, 0.0), "N")
    current_ca = _atom(current, "CA", (1.0, 0.0, 0.0))
    current_c = _atom(current, "C", (1.0, 1.0, 0.0))
    structure = _structure(
        ((previous, "ALA"), (current, "ALA")),
        (previous_c, current_n, current_ca, current_c),
        ((previous_c.identity, current_n.identity),),
    )

    phi, psi = backbone_torsions(
        structure,
        current,
        (
            BoundaryIdentity(previous, current),
            BoundaryIdentity(current, missing_following),
        ),
    )

    assert phi.reason == "missing_residue"
    assert psi.reason == "missing_residue"
    assert phi.metadata["ramachandran_class"] is None


def test_peptide_geometry_requires_explicit_same_chain_bond() -> None:
    left = ResidueIdentity("A", "1")
    right = ResidueIdentity("A", "2")
    left_ca = _atom(left, "CA", (0.0, 0.0, 0.0))
    left_c = _atom(left, "C", (1.0, 0.0, 0.0))
    left_o = _atom(left, "O", (1.0, 1.0, 0.0), "O")
    right_n = _atom(right, "N", (2.0, 0.0, 1.0), "N")
    right_ca = _atom(right, "CA", (3.0, 1.0, 1.0))
    residues = ((left, "ALA"), (right, "GLY"))
    atoms = (left_ca, left_c, left_o, right_n, right_ca)
    boundary = BoundaryIdentity(left, right)

    unbonded = peptide_geometry(_structure(residues, atoms), boundary)
    assert {result.reason for result in unbonded} == {"unsupported_connectivity"}

    bonded = peptide_geometry(
        _structure(residues, atoms, ((left_c.identity, right_n.identity),)),
        boundary,
        boundary_kind="junction",
    )
    assert all(result.measurement is not None for result in bonded)
    assert all(result.metadata["boundary_kind"] == "junction" for result in bonded)
    assert all(
        result.metadata == result.measurement.metadata
        for result in bonded
        if result.measurement is not None
    )
    omega_planarity = next(
        result for result in bonded if result.metric == "peptide_omega_planarity.v1"
    )
    improper_planarity = next(
        result for result in bonded if result.metric == "peptide_carbonyl_improper_planarity.v1"
    )
    assert "distance_to_cis_degrees" in omega_planarity.metadata
    assert "distance_to_trans_degrees" in omega_planarity.metadata
    assert "distance_to_cis_degrees" not in improper_planarity.metadata
    assert "distance_to_nearest_plane_degrees" in improper_planarity.metadata

    cross_chain_right = ResidueIdentity("a", "2")
    cross_atoms = atoms[:-2] + (
        _atom(cross_chain_right, "N", (2.0, 0.0, 1.0), "N"),
        _atom(cross_chain_right, "CA", (3.0, 1.0, 1.0)),
    )
    cross = peptide_geometry(
        _structure(((left, "ALA"), (cross_chain_right, "GLY")), cross_atoms),
        BoundaryIdentity(left, cross_chain_right),
    )
    assert {result.reason for result in cross} == {"unsupported_connectivity"}


def test_sidechain_chi_coverage_and_atomless_identity() -> None:
    arg = ResidueIdentity("H", "82", "A")
    coordinates = {
        "N": (0.0, 0.0, 0.0),
        "CA": (1.0, 0.0, 0.0),
        "CB": (1.0, 1.0, 0.0),
        "CG": (2.0, 1.0, 1.0),
        "CD": (3.0, 2.0, 1.0),
        "NE": (4.0, 2.0, 2.0),
        "CZ": (5.0, 3.0, 2.0),
        "NH1": (6.0, 3.0, 3.0),
    }
    structure = _structure(
        ((arg, "ARG"),),
        tuple(
            _atom(arg, name, point, "N" if name in {"N", "NE", "NH1"} else "C")
            for name, point in coordinates.items()
        ),
    )

    chis = sidechain_torsions(structure, arg)
    assert [result.metric for result in chis] == [
        "sidechain_chi1.v1",
        "sidechain_chi2.v1",
        "sidechain_chi3.v1",
        "sidechain_chi4.v1",
        "sidechain_chi5.v1",
    ]
    assert all(result.measurement is not None for result in chis)

    missing = ResidueIdentity("h", "82", "A")
    undefined = sidechain_torsions(structure, missing)[0].to_dict()
    assert undefined["reason"] == "missing_residue"
    assert undefined["metadata"] == {"residue": missing.to_dict()}


def test_connectivity_supports_cyx_junctions_and_declared_terminals() -> None:
    first = ResidueIdentity("A", "1")
    second = ResidueIdentity("A", "2")
    first_sg = _atom(first, "SG", (0.0, 0.0, 0.0), "S")
    second_sg = _atom(second, "SG", (2.0, 0.0, 0.0), "S")
    disulfide = BoundaryIdentity(first, second)
    cyx = _structure(
        ((first, "CYX"), (second, "CYX")),
        (first_sg, second_sg),
        ((first_sg.identity, second_sg.identity),),
    )

    results = connectivity_measurements(cyx, disulfides=(disulfide,))
    disulfide_presence = [
        result
        for result in results
        if result.metric == "canonical_connectivity.v1"
        and result.metadata.get("bond_kind") == "disulfide"
    ]
    assert len(disulfide_presence) == 1
    disulfide_measurement = disulfide_presence[0].measurement
    assert disulfide_measurement is not None
    assert disulfide_measurement.value == 1.0
    assert disulfide_presence[0].metadata["definition_source"] == "declared_disulfide"

    missing_terminal = ResidueIdentity("D", "82", "A")
    terminal = connectivity_measurements(
        ProteinStructure((), ()),
        c_termini=frozenset((missing_terminal,)),
    )[0]
    assert terminal.metric == "terminal_atom_presence.v1"
    assert terminal.reason == "missing_residue"
    assert terminal.atoms == (AtomIdentity("D", "82", "A", "OXT"),)

    left = ResidueIdentity("J", "10")
    right = ResidueIdentity("J", "11")
    left_c = _atom(left, "C", (0.0, 0.0, 0.0))
    right_n = _atom(right, "N", (1.3, 0.0, 0.0), "N")
    junction_results = connectivity_measurements(
        _structure(
            ((left, "ALA"), (right, "GLY")),
            (left_c, right_n),
            ((left_c.identity, right_n.identity),),
        ),
        junctions=(BoundaryIdentity(left, right),),
    )
    junction = [
        result
        for result in junction_results
        if result.metric == "canonical_connectivity.v1"
        and result.metadata.get("bond_kind") == "junction"
    ]
    assert len(junction) == 1
    junction_measurement = junction[0].measurement
    assert junction_measurement is not None
    assert junction_measurement.value == 1.0
    assert junction[0].metadata["definition_source"] == "declared_junction"

    terminal_residue = ResidueIdentity("T", "1")
    terminal_atoms = (
        _atom(terminal_residue, "N", (0.0, 0.0, 0.0), "N"),
        _atom(terminal_residue, "CA", (1.0, 0.0, 0.0)),
        _atom(terminal_residue, "C", (2.0, 0.0, 0.0)),
        _atom(terminal_residue, "O", (2.5, 1.0, 0.0), "O"),
        _atom(terminal_residue, "OXT", (2.5, -1.0, 0.0), "O"),
        _atom(terminal_residue, "CB", (1.0, 1.0, 0.0)),
    )
    terminal_results = connectivity_measurements(
        _structure(((terminal_residue, "ALA"),), terminal_atoms),
        n_termini=frozenset((terminal_residue,)),
        c_termini=frozenset((terminal_residue,)),
    )
    terminal_bond = next(
        result
        for result in terminal_results
        if result.metric == "canonical_connectivity.v1"
        and {atom.atom for atom in result.atoms} == {"C", "OXT"}
    )
    template_bond = next(
        result
        for result in terminal_results
        if result.metric == "canonical_connectivity.v1"
        and {atom.atom for atom in result.atoms} == {"CA", "CB"}
    )
    assert terminal_bond.metadata["definition_source"] == "declared_terminal"
    assert template_bond.metadata["definition_source"] == "openmm_residue_template"


def test_canonical_connectivity_reference_is_immutable() -> None:
    bonds: Any = canonical_heavy_bonds()
    original = bonds["ALA"]

    with pytest.raises(TypeError):
        bonds["ALA"] = ()

    assert canonical_heavy_bonds()["ALA"] == original


def test_steric_search_covers_large_radii_and_topology_exclusions() -> None:
    first = ResidueIdentity("A", "1")
    second = ResidueIdentity("A", "2")
    first_k = _atom(first, "K", (0.0, 0.0, 0.0), "K")
    second_k = _atom(second, "K", (5.2, 0.0, 0.0), "K")
    structure = _structure(((first, "UNK"), (second, "UNK")), (first_k, second_k))

    overlaps = steric_measurements(
        structure,
        generated=frozenset((first_k.identity,)),
        fixed=frozenset((second_k.identity,)),
    )

    assert len(overlaps) == 1
    assert overlaps[0].measurement is not None
    assert overlaps[0].measurement.value > 0.0
    assert overlaps[0].metadata["pair_class"] == "generated/fixed"
    assert overlaps[0].metadata["cutoff_angstrom"] > 5.2

    bonded = _structure(
        ((first, "UNK"), (second, "UNK")),
        (first_k, second_k),
        ((first_k.identity, second_k.identity),),
    )
    assert (
        steric_measurements(
            bonded,
            generated=frozenset((first_k.identity, second_k.identity)),
            fixed=frozenset(),
        )
        == ()
    )
    assert (
        len(
            steric_measurements(
                bonded,
                generated=frozenset((first_k.identity, second_k.identity)),
                fixed=frozenset(),
                exclude_bond_hops=0,
            )
        )
        == 1
    )


def test_steric_roles_and_hydrogen_selection_are_explicit() -> None:
    generated_residue = ResidueIdentity("A", "1")
    fixed_residue = ResidueIdentity("A", "2")
    hydrogen = _atom(generated_residue, "H", (0.0, 0.0, 0.0), "H")
    oxygen = _atom(fixed_residue, "O", (1.0, 0.0, 0.0), "O")
    structure = _structure(
        ((generated_residue, "UNK"), (fixed_residue, "UNK")),
        (hydrogen, oxygen),
    )

    assert (
        steric_measurements(
            structure,
            generated=frozenset((hydrogen.identity,)),
            fixed=frozenset((oxygen.identity,)),
        )
        == ()
    )
    hydrogen_aware = steric_measurements(
        structure,
        generated=frozenset((hydrogen.identity,)),
        fixed=frozenset((oxygen.identity,)),
        include_hydrogens=True,
    )
    assert len(hydrogen_aware) == 1
    assert hydrogen_aware[0].metadata["hydrogen_pair"] is True
    assert hydrogen_aware[0].metadata["pair_class"] == "generated/fixed"

    junction = steric_measurements(
        structure,
        generated=frozenset(),
        fixed=frozenset((oxygen.identity,)),
        junction=frozenset((hydrogen.identity,)),
        include_hydrogens=True,
    )
    assert len(junction) == 1
    assert junction[0].metadata["pair_class"] == "junction"
