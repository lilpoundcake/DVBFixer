"""Close contacts must never manufacture their own steric exclusions."""

from __future__ import annotations

from io import StringIO
from types import SimpleNamespace

import pytest

pytest.importorskip("openmm")
pytest.importorskip("scipy")

from openmm.app import Element, PDBFile, Topology  # noqa: E402
from openmm.unit import Quantity, angstrom  # noqa: E402

from dvbfixer.diagnose import steric  # noqa: E402
from dvbfixer.diagnose.report import Severity  # noqa: E402


def _structure(specs, coordinates, bonds=()):
    top = Topology()
    chain = top.addChain("A")
    atoms = []
    for index, (resname, names) in enumerate(specs):
        residue = top.addResidue(resname, chain, id=str(index + 1))
        for name, element in names:
            atoms.append(top.addAtom(name, Element.getBySymbol(element), residue))
    for i, j in bonds:
        top.addBond(atoms[i], atoms[j])
    return top, Quantity(coordinates, angstrom)


@pytest.mark.parametrize("distance", [2.0, 1.9, 1.8, 0.0])
def test_shorter_nonbonded_contact_never_disappears(distance):
    top, pos = _structure([("L1", [("C1", "C")]), ("L2", [("C2", "C")])],
                          [(0, 0, 0), (distance, 0, 0)])
    findings = steric.clashes_python(top, pos)
    measurements = steric.steric_overlap_measurements_python(top, pos)
    assert len(findings) == len(measurements) == 1
    assert findings[0].severity is Severity.ERROR
    assert measurements[0].value == pytest.approx(3.4 - distance)
    assert measurements[0].metadata["bond_exclusion_policy"] == steric.BOND_EXCLUSION_POLICY


@pytest.mark.parametrize("include_water, count", [(False, 0), (True, 1)])
def test_carbon_water_contact(include_water, count):
    top, pos = _structure([("ALA", [("CB", "C")]), ("HOH", [("O", "O")])],
                          [(0, 0, 0), (1.8, 0, 0)])
    findings = steric.clashes_python(top, pos, include_water=include_water)
    assert len(findings) == count
    if count:
        assert findings[0].severity is Severity.ERROR


def test_ignored_water_cannot_bridge_nonwater_atoms():
    top, pos = _structure(
        [("L1", [("C1", "C")]), ("L2", [("C2", "C")]), ("HOH", [("O", "O")])],
        [(-1.2, 0, 0), (1.2, 0, 0), (0, 0, 0)],
    )
    findings = steric.clashes_python(top, pos)
    assert len(findings) == 1
    assert findings[0].extra["overlap_angstrom"] == pytest.approx(1.0)


def test_false_edge_cannot_suppress_neighbors():
    top, pos = _structure(
        [("LIG", [("C1", "C"), ("C2", "C")]), ("HOH", [("O", "O")])],
        [(0, 0, 0), (0, 1.5, 0), (1.8, 0, 0)], [(0, 1)],
    )
    findings = steric.clashes_python(top, pos, include_water=True)
    assert len(findings) == 2
    assert {f.atom for f in findings} == {"C1", "C2"}


def test_very_short_h_acceptor_contact_is_not_a_bond():
    top, pos = _structure(
        [("ALA", [("N", "N"), ("H", "H")]), ("HOH", [("O", "O")])],
        [(0, 0, 0), (1, 0, 0), (2, 0, 0)], [(0, 1)],
    )
    findings = steric.clashes_python(top, pos, include_water=True)
    assert any(f.atom == "H" and f.category == "clash" for f in findings)


def test_only_template_edges_can_assign_a_hydrogen_parent():
    top, pos = _structure(
        [("ALA", [("CB", "C"), ("HB1", "H")]), ("SER", [("OG", "O")])],
        [(0, 0, 0), (1.09, 0, 0), (2.94, 0, 0)],
    )
    graph, unresolved = steric._build_neighbor_graph(top, False)
    assert graph == {0: {1}, 1: {0}}
    assert not unresolved
    assert any(f.atom == "HB1" for f in steric.clashes_python(top, pos))


def test_duplicate_atom_names_are_not_selected_for_template_edges():
    top, pos = _structure([("ALA", [("CB", "C"), ("CB", "C"), ("HB1", "H")])],
                          [(0, 0, 0), (1, 0, 0), (2, 0, 0)])
    graph, unresolved = steric._build_neighbor_graph(top, False)
    assert graph == {}
    assert len(unresolved) == 1
    assert [f.category for f in steric.clashes_python(top, pos)] == ["steric_connectivity"]


def test_templates_preserve_case_sensitive_chains_and_insertion_codes():
    top = Topology()
    for chain_id in ("D", "d"):
        chain = top.addChain(chain_id)
        for icode in ("", "A"):
            residue = top.addResidue("HIE", chain, id="82", insertionCode=icode)
            top.addAtom("CB", Element.getBySymbol("C"), residue)
            top.addAtom("CG", Element.getBySymbol("C"), residue)
    graph, unresolved = steric._build_neighbor_graph(top, False)
    assert not unresolved
    assert graph == {i: {i + 1 if i % 2 == 0 else i - 1} for i in range(8)}


def test_empty_measurement_set_has_explicit_provenance():
    provenance = {}
    assert steric.steric_overlap_measurements_python(Topology(), [], provenance=provenance) == []
    assert provenance["status"] == "complete"
    assert provenance["unresolved_residues"] == []


@pytest.mark.parametrize("resname", ["HIS", "HIE", "HSD", "GLH", "GLUP"])
def test_template_graph_is_coordinate_independent_and_does_not_mutate(resname):
    names = [("CB", "C"), ("CG", "C")]
    top, pos = _structure([(resname, names), ("LIG", [("C1", "C")])],
                          [(0, 0, 0), (0, 1.5, 0), (1.8, 0, 0)])
    graph, unresolved = steric._build_neighbor_graph(top, include_water=False)
    assert graph == {0: {1}, 1: {0}}
    assert unresolved == []
    assert top.getNumBonds() == 0
    assert next(top.residues()).name == resname
    assert len(steric.clashes_python(top, pos)) == 2


@pytest.mark.parametrize("resname", ["LIG", "NAG", "0GA"])
def test_missing_connectivity_warns_even_without_neighbor_pairs(resname):
    top, pos = _structure([(resname, [("C1", "C"), ("C2", "C")])],
                          [(0, 0, 0), (10, 0, 0)])
    provenance = {}
    findings = steric.clashes_python(top, pos, provenance=provenance)
    assert [f.category for f in findings] == ["steric_connectivity"]
    assert findings[0].severity is Severity.WARNING
    assert provenance["status"] == "partial"
    assert provenance["unresolved_residues"][0]["resname"] == resname
    measurement_provenance = {}
    assert steric.steric_overlap_measurements_python(top, pos, provenance=measurement_provenance) == []
    assert measurement_provenance == provenance


def test_missing_intra_residue_graph_cannot_create_external_exclusions():
    top, pos = _structure(
        [("LIG", [("C1", "C"), ("C2", "C")]), ("ALA", [("CB", "C")])],
        [(0, 0, 0), (0, 1.5, 0), (1.8, 0, 0)],
    )
    findings = steric.clashes_python(top, pos)
    assert sum(f.category == "steric_connectivity" for f in findings) == 1
    assert sum(f.category == "clash" for f in findings) == 2


def test_partial_connectivity_status_survives_dispatch(monkeypatch):
    top, pos = _structure([("LIG", [("C1", "C"), ("C2", "C")])],
                          [(0, 0, 0), (1.5, 0, 0)])
    monkeypatch.setattr(steric, "_find_binary", lambda _: None)
    provenance = {}
    findings = steric.run_all(top, pos, "unused.pdb", provenance=provenance)
    assert [f.category for f in findings] == ["steric_connectivity"]
    assert provenance["fallback"] is True
    assert provenance["status"] == "partial"


def test_close_intra_residue_atoms_cannot_shortcut_a_template_path():
    top, pos = _structure(
        [("LYS", [("CB", "C"), ("CG", "C"), ("CD", "C"), ("CE", "C"), ("NZ", "N")]),
         ("LIG", [("C1", "C")])],
        [(0, 0, 0), (-1.5, 0, 0), (-1.5, 1.5, 0), (0, 1.5, 0), (0, 0, 1.5), (1.8, 0, 0)],
        [(4, 5)],
    )
    graph, unresolved = steric._build_neighbor_graph(top, False)
    assert 4 not in graph[0]  # folded CB/NZ geometry is not a new bond
    assert not unresolved
    findings = steric.clashes_python(top, pos)
    assert any(f.atom == "CB" and f.extra["clash_partner"]["atom"] == "C1" for f in findings)


@pytest.mark.parametrize("left,right", [
    (("ALA", "C", "C"), ("GLY", "N", "N")),
    (("CYS", "SG", "S"), ("CYS", "SG", "S")),
    (("NAG", "C1", "C"), ("NAG", "O4", "O")),
])
def test_supplied_inter_residue_bonds_and_one_three_one_four_survive(left, right):
    left_res, first, first_element = left
    right_res, second, second_element = right
    top, pos = _structure(
        [(left_res, [(first, first_element), ("X", "C")]),
         (right_res, [(second, second_element), ("Y", "C")])],
        [(0, 0, 0), (0, 1.5, 0), (1.5, 0, 0), (1.5, 1.5, 0)],
        [(0, 1), (0, 2), (2, 3)],
    )
    assert steric.clashes_python(top, pos) == []
    assert top.getNumBonds() == 3


@pytest.mark.parametrize("mode", ["python", "missing_probe", "failed_probe", "probe"])
def test_loaded_pdb_dispatch_and_fallback(monkeypatch, mode):
    pdb = PDBFile(StringIO(
        "HETATM    1  C1  LG1 A   1       0.000   0.000   0.000  1.00 20.00           C\n"
        "HETATM    2  C2  LG2 A   2       1.800   0.000   0.000  1.00 20.00           C\nEND\n"
    ))
    monkeypatch.setattr(steric, "_find_binary", lambda _: None if mode == "missing_probe" else "probe")
    from dvbfixer.diagnose import probe
    monkeypatch.setattr(probe, "analyze", lambda *a, **kw: SimpleNamespace(
        findings=[] if mode == "probe" else None, status="complete" if mode == "probe" else "failed",
        version="2.18", reason="test failure",
    ))
    provenance = {}
    findings = steric.run_all(pdb.topology, pdb.positions, "unused.pdb",
                              prefer_probe=mode != "python", provenance=provenance)
    assert len(findings) == (0 if mode == "probe" else 1)
    assert provenance["actual_engine"] == ("molprobity-probe" if mode == "probe" else steric.DETERMINISTIC_STERIC_ENGINE)
    if mode != "probe":
        assert findings[0].severity is Severity.ERROR
        assert provenance["bond_exclusion_policy"] == steric.BOND_EXCLUSION_POLICY
