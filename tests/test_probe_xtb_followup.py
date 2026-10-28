"""Identity and numerical contracts for the research-only follow-up tools."""
import importlib
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("openmm")
pytest.importorskip("rdkit")


@pytest.fixture
def probe(monkeypatch):
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[1] / "scripts"))
    return importlib.import_module("audit_probe_adapter")


def record(kind="bo", source=None, target=None, gap="-0.848"):
    def field(chain, resid, icode, name, atom):
        return f"{chain:>2}{resid:>4}{icode:1}{name:>3} {atom:<4} "
    source = source or field("H", "82", "a", "0FA", "H1")
    target = target or field("h", "-1", "B", "LIG", "C1")
    return ":".join(["", "1->1", kind, source, target, "11", gap, "-0.476",
                     "0", "0", "0", "0", "0", "H", "C", "0", "0", "0", "0", "0"]) + "\n"


def inventory():
    return {("H", "82", "a", "h1", ""): {"resname": "0fA"},
            ("h", "-1", "B", "C1", ""): {"resname": "LIG"}}


def test_probe_resolves_lossy_names_without_normalizing_chain_or_icode(probe):
    contacts, findings = probe.parse_contacts(record() + record("hb"), inventory())
    assert len(contacts) == 2 and len(findings) == 1
    assert findings[0]["atoms"] == tuple(sorted(inventory()))
    assert findings[0]["contact_types"] == ["bo", "hb"]
    assert findings[0]["overlap_angstrom"] == pytest.approx(.848)


@pytest.mark.parametrize("raw", [
    "unsupported\n", record() + "unsupported\n", record(gap="nan"),
    record().replace("1->1", "2->1"), record().replace(" H  82a", " H  82A"),
    record().replace(":11:", ":0:"), record().replace(":-0.476:", ":inf:"),
])
def test_probe_rejects_malformed_or_unresolved_records(probe, raw):
    with pytest.raises(ValueError):
        probe.parse_contacts(raw, inventory())


def test_probe_rejects_noninjective_atom_aliases_even_if_stdout_empty(probe):
    atoms = inventory() | {("H", "82", "a", "H1", ""): {"resname": "0fA"}}
    with pytest.raises(ValueError, match="ambiguous"):
        probe.parse_contacts("", atoms)


def test_probe_clean_result_and_threshold_are_not_failures(probe):
    assert probe.parse_contacts("", inventory()) == ([], [])
    assert probe.parse_contacts(record(gap="-0.399"), inventory())[1] == []
    assert len(probe.parse_contacts(record(gap="-0.400"), inventory())[1]) == 1


def test_harmonic_restraint_gradient_units_and_mask(monkeypatch):
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[1] / "scripts"))
    context = importlib.import_module("audit_xtb_context")
    x = np.array([[.2, .3, -.1], [.1, -.2, .4], [-.2, .1, .2]])
    reference = np.zeros_like(x)
    heavy = np.array([0, 2])
    original_gradient = x.copy()

    def evaluate(coords):
        return context.restrained_energy(.5 * np.sum(coords**2), coords, coords, reference, heavy, 5.)

    _, gradient = evaluate(x)
    assert np.array_equal(x, original_gradient)
    assert np.array_equal(gradient[1], x[1])
    for atom in range(3):
        for axis in range(3):
            plus, minus = x.copy(), x.copy()
            plus[atom, axis] += 1e-5
            minus[atom, axis] -= 1e-5
            finite_difference = (evaluate(plus)[0] - evaluate(minus)[0]) / 2e-5
            assert gradient[atom, axis] == pytest.approx(finite_difference, abs=1e-10)


def test_clash_severity_change_is_not_a_new_physical_pair(monkeypatch):
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[1] / "scripts"))
    summary = importlib.import_module("summarize_probe_xtb_followup")
    pair = (("H", "82", "A", "CA"), ("h", "82", "", "C1"))
    other = (("H", "83", "", "CB"), ("h", "82", "", "C1"))
    lost = [("clash", "WARNING", pair, "H", "82", "CA")]
    new = [("clash", "ERROR", pair, "H", "82", "CA"), ("clash", "WARNING", other, "H", "83", "CB")]
    assert summary.physical_clash_changes(new, lost) == (1, 0)


def test_context_optimizer_never_exposes_frozen_coordinates_as_variables(tmp_path, monkeypatch):
    import json
    from types import SimpleNamespace

    from openmm.app import Element, Topology

    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[1] / "scripts"))
    context = importlib.import_module("audit_xtb_context")
    top, sub = Topology(), Topology()
    top.addAtom("CA", Element.getBySymbol("C"), top.addResidue("ALA", top.addChain("P"), "1"))
    top.addAtom("C1", Element.getBySymbol("C"), top.addResidue("LIG", top.addChain("a"), "82", insertionCode="B"))
    sub.addAtom("C1", Element.getBySymbol("C"), sub.addResidue("LIG", sub.addChain("a"), "82", insertionCode="B"))
    reference = np.array([[0., 0., 0.], [3., 0., 0.]])
    monkeypatch.setattr(context, "load_case", lambda *a: (SimpleNamespace(topology=top), sub, reference,
                                                       reference[1:], {1: 0}, [1]))
    evaluated = []

    class Quadratic:
        def __init__(self, *args):
            self.numeric_libraries = []

        def evaluate(self, coordinates):
            evaluated.append(coordinates.copy())
            # The backend actively wants to move the protein too. Only
            # projection, not a zero protein gradient, keeps it stationary.
            delta = coordinates - np.array([[10., 0., 0.], [4., 0., 0.]])
            return .5 * float(np.sum(delta**2)), delta

        def close(self):
            pass

    monkeypatch.setattr(context, "GFNFF", Quadratic)
    monkeypatch.setattr(context, "finding_sets", lambda *a: (set(), set()))
    monkeypatch.setattr(context, "geometry_snapshot", lambda *a: {})
    monkeypatch.setattr(context, "pocket_snapshot", lambda *a: [])
    monkeypatch.setattr(context, "geometry_deltas", lambda *a: {
        "bonds_changed_over_0_1_angstrom": [], "angles_changed_over_10_degrees": [], "potential_stereovolume_flips": [],
    })
    baseline = tmp_path / "input"
    baseline.write_text("test reference identity")
    context.worker(SimpleNamespace(baseline=baseline, glycam=False, scope="context", k=0., maxiter=50,
                                   output=tmp_path, library=tmp_path / "unused.so"))
    assert evaluated and all(np.array_equal(coords[0], reference[0]) for coords in evaluated)
    final = np.load(tmp_path / "raw-final.npy", allow_pickle=False)
    assert final[1] == pytest.approx([4, 0, 0])
    assert json.loads((tmp_path / "report.json").read_text())["converged"]
