"""Conditional chemistry/transport regressions, not actual-polymer qualification."""
from __future__ import annotations

import hashlib
import json
from dataclasses import replace
from pathlib import Path

import pytest
from scripts.chemical_reconstruction_peptide_controls import POLICY, inventory, load_ccd, run_probes

from dvbfixer.chemical_reconstruction.ccd_authority import CcdAuthorityError, verify_cip
from dvbfixer.chemical_reconstruction.construction import propose
from dvbfixer.chemical_reconstruction.contracts import Operation
from dvbfixer.chemical_reconstruction.observed_peptide import (
    compile_peptide_plan,
    prepare_observed_probe,
    select_interior_target,
)
from dvbfixer.chemical_reconstruction.peptide_state import (
    PeptideStateError,
    assemble_witness,
    derive_local_peptide,
)
from dvbfixer.chemical_reconstruction.proposals import PlacementKind
from dvbfixer.domain.structure_identity import ComponentInstanceRef

TARGET = ComponentInstanceRef(1, "A", "2", "B")


def local(name="ALA"):
    return derive_local_peptide(load_ccd(name), load_ccd("ALA"), TARGET, inventory(name))


@pytest.mark.parametrize("name", ["ALA", "THR", "LYS", "PHE"])
def test_explicit_substitutions_prove_chemistry_without_witness_geometry(name: str) -> None:
    a = local(name)
    assert a.removed_atoms == {"OXT", "HXT", "H2"}
    assert len(a.witness.extra_bonds) == 2
    assert a.state.graph.net_charge == (1 if name == "LYS" else 0)
    assert a.parent.state.graph.net_charge == a.state.graph.net_charge
    assert {x.name for x in a.state.graph.atoms} == a.inventory.atoms
    assert {x.name for x in a.parent.state.graph.atoms} >= {"OXT", "HXT", "H2"}
    positions = {x.name: x.ideal_position for x in a.witness.flat_fragments.atoms}
    with pytest.raises(CcdAuthorityError, match="radical or incomplete"):
        verify_cip(a.witness.flat_fragments, positions)
    assert verify_cip(a.witness.flat_fragments, positions, extra_bonds=a.witness.extra_bonds) == a.witness.verified_cip


def test_unexpected_ile_naming_refusal_and_explicit_lyn_gate_are_retained() -> None:
    with pytest.raises(PeptideStateError, match="inventory-mismatch"):
        local("ILE")
    with pytest.raises(PeptideStateError, match="variant-projection-unqualified"):
        derive_local_peptide(load_ccd("LYS"), load_ccd("ALA"), TARGET, inventory("LYN"))


@pytest.mark.parametrize("mutation", ["nonleaving", "port-reuse", "partial-branch", "disconnected"])
def test_no_blanket_leaving_atom_deletion_or_incomplete_assembly(mutation: str) -> None:
    witness = local().witness
    links = witness.substitutions
    if mutation == "nonleaving":
        source = links[0].carbon.component
        from dvbfixer.domain.structure_identity import ExactAtomRef
        links = (replace(links[0], removed=links[0].removed | {ExactAtomRef(source, "CB")}), links[1])
    elif mutation == "port-reuse":
        links = (links[0], links[0])
    elif mutation == "partial-branch":
        links = (replace(links[0], removed=frozenset(a for a in links[0].removed if a.name != "HXT")), links[1])
    else:
        links = ()
    with pytest.raises(ValueError):
        assemble_witness(witness.source_units, links)


def _line(serial, name, residue, number, xyz, element, *, chain="A", code=" ", record="ATOM"):
    x, y, z = xyz
    atom_field = f" {name:<3}" if len(name) < 4 else name
    return (f"{record:<6}{serial:5d} {atom_field} {residue:>3} {chain}{number:4d}{code}   "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00          {element:>2}\n")


def synthetic_observed(name="ALA", *, rotate_carbonyl=False):
    """Incomplete opaque neighbours deliberately do NOT qualify whole geometry."""
    a = local(name)
    points = {x.name: x.ideal_position for x in a.state.graph.atoms}
    n, c = points["N"], points["C"]
    lines = ["HEADER    SYNTHETIC SOURCE TRANSPORT ONLY\n",
             _line(1, "N", "GLY", 1, (-20., 0., 0.), "N"),
             _line(2, "C", "GLY", 1, (n[0], n[1], n[2]+1.3), "C")]
    serial = 3
    for atom in a.state.graph.atoms:
        if atom.element == "H":
            continue
        point = atom.ideal_position
        if rotate_carbonyl and atom.name == "O":
            # Rotate O by pi about CA-C, without moving CA/N/C or any withheld CB.
            import math

            from dvbfixer.chemical_reconstruction.measurement import _dot, _sub, _unit
            axis = _unit(_sub(points["CA"], c))
            delta = _sub(point, c)
            point = tuple(c[i]+2*_dot(delta, axis)*axis[i]-delta[i] for i in range(3))
            assert math.dist(point, c) == pytest.approx(math.dist(atom.ideal_position, c))
        lines.append(_line(serial, atom.name, name, 2, point, atom.element, code="B"))
        serial += 1
    lines.extend((_line(serial, "N", "GLY", 3, (c[0], c[1], c[2]+1.3), "N"),
                  _line(serial+1, "C", "GLY", 3, (20., 0., 0.), "C"),
                  _line(serial+2, "C1", "LIG", 99, (100., 100., 100.), "C", chain="a", record="HETATM"),
                  _line(serial+3, "O1", "LIG", 99, (101.4, 100., 100.), "O", chain="a", record="HETATM"),
                  f"CONECT{serial+2:5d}{serial+3:5d}\n", "END\n"))
    return "".join(lines).encode("ascii"), a


@pytest.mark.parametrize("name,mask", [("ALA", frozenset({"CB"})), ("THR", frozenset({"CB", "OG1", "CG2"}))])
def test_whole_source_context_survives_and_witnesses_never_become_observations(name, mask) -> None:
    source, a = synthetic_observed(name)
    probe = prepare_observed_probe(source, a, withheld=mask)
    assert probe.source_atom_count-len(probe.withheld) == len(probe.request.observations)
    assert probe.request.source_sha256 == hashlib.sha256(source).hexdigest()
    assert len(probe.unassessed_internal_connections) == 1
    assert {o.identity.component.chain_id for o in probe.request.observations} == {"A", "a"}
    assert all("witness" not in c.identity.chain_id for c in probe.request.components)
    plan = compile_peptide_plan(probe, POLICY)
    assert plan.placements[0].kind == PlacementKind.CONE
    result = propose(probe.request, plan)
    assert result.status == "declared-geometry-passed", result
    assert not result.geometry_approved and not result.md_ready and not result.source_coverage_verified
    assert {x.identity for x in result.raw.additions} == probe.request.ownership.generated


def test_cb_cone_is_not_driven_by_carbonyl_phase() -> None:
    candidates = []
    for rotated in (False, True):
        source, a = synthetic_observed(rotate_carbonyl=rotated)
        probe = prepare_observed_probe(source, a, withheld=frozenset({"CB"}))
        result = propose(probe.request, compile_peptide_plan(probe, POLICY))
        assert result.status == "declared-geometry-passed", result
        candidates.append(result.raw.additions[0].position)
    assert candidates[0] == candidates[1]


@pytest.mark.parametrize("mutation", ["ter", "dangling", "target-h", "mask-backbone", "h-operation"])
def test_source_or_geometry_gates_refuse_without_source_repair(mutation: str) -> None:
    source, a = synthetic_observed()
    mask = frozenset({"CB"})
    operation = Operation.HEAVY
    if mutation == "ter":
        lines = source.decode().splitlines(keepends=True)
        source = "".join((*lines[:3], "TER\n", *lines[3:])).encode()
    elif mutation == "dangling":
        source = source.replace(b"END\n", b"CONECT99999    1\nEND\n")
    elif mutation == "target-h":
        lines = source.decode().splitlines(keepends=True)
        index = next(i for i, line in enumerate(lines) if line.startswith("ATOM") and line[22:26].strip() == "3")
        lines.insert(index, _line(999, "H", "ALA", 2, (0., 0., 0.), "H", code="B"))
        source = "".join(lines).encode()
    elif mutation == "mask-backbone":
        mask = frozenset({"CA"})
    else:
        operation = Operation.HYDROGENS
    before = source
    with pytest.raises(ValueError):
        prepare_observed_probe(source, a, withheld=mask, operation=operation)
    assert source == before


def test_file_order_selection_preserves_case_and_insertion_code() -> None:
    source, _ = synthetic_observed("THR")
    assert select_interior_target(source, "THR", chain="A") == TARGET
    with pytest.raises(PeptideStateError):
        select_interior_target(source, "THR", chain="a")


def test_frozen_probe_regression_does_not_count_admission_as_coordinate_success(tmp_path: Path) -> None:
    output = tmp_path / "probes"
    s = run_probes(output)
    assert s["declared"] == s["executed_entries"] == 12
    chemistry = [c for c in s["controls"] if c["kind"] == "conditional-chemistry"]
    assert sum(c["status"] == "conditional-chemistry-passed" for c in chemistry) == 4
    observed = [c for c in s["controls"] if c["kind"] == "observed-probe"]
    assert all(c["status"] == "refused" for c in observed)
    for c in observed:
        evidence = json.loads((output / f"{c['id']}.json").read_bytes())
        assert "result" not in evidence  # No coordinate constructor executed.
    with pytest.raises(FileExistsError):
        run_probes(output)
