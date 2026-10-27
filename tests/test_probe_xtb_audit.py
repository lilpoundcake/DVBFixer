"""Dependency audit contracts; known safety gaps remain strict xfails.

These tests do not change the scientific adapters or require either executable.
"""

from subprocess import CompletedProcess, TimeoutExpired

import pytest

pytest.importorskip("openmm")

from openmm import Vec3
from openmm.app import Element, Topology
from openmm.unit import Quantity, nanometer

from dvbfixer.diagnose import steric
from dvbfixer.minimize import refine


@pytest.fixture
def molecule():
    top = Topology()
    residue = top.addResidue("LIG", top.addChain("a"), "1")
    top.addAtom("C1", Element.getBySymbol("C"), residue)
    top.addAtom("O1", Element.getBySymbol("O"), residue)
    return top, Quantity([Vec3(0, 0, 0), Vec3(0.14, 0, 0)], nanometer)


@pytest.mark.parametrize("failure", ["missing", "timeout", "nonzero"])
def test_probe_failure_falls_back(monkeypatch, molecule, failure):
    def run(*args, **kwargs):
        if failure == "missing":
            raise FileNotFoundError("probe")
        if failure == "timeout":
            raise TimeoutExpired("probe", 60)
        return CompletedProcess(args, 1, "", "failure")

    monkeypatch.setattr(steric, "_find_binary", lambda name: "probe")
    monkeypatch.setattr(steric.subprocess, "run", run)
    sentinel = [object()]
    monkeypatch.setattr(steric, "clashes_python", lambda *a, **kw: sentinel)
    assert steric.run_all(*molecule, "input.pdb") is sentinel


@pytest.mark.xfail(strict=True, reason="incompatible Probe stdout is accepted as a clean result")
def test_probe_incompatible_output_not_clean(monkeypatch):
    monkeypatch.setattr(steric, "_run_probe", lambda *a: "unsupported output format\n")
    assert steric.clashes_probe("input.pdb", "probe") is None


@pytest.mark.xfail(strict=True, reason="default self mode accepts one selection, not two")
def test_probe_command_uses_valid_selection_count(monkeypatch):
    def run(command, **kwargs):
        selections = command[command.index("-het") + 1:-1]
        assert selections == ["ALL"] or selections in (
            ["-self", "ALL"], ["-both", "ALL", "ALL"], ["-once", "ALL", "ALL"]
        )
        return CompletedProcess(command, 0, "", "")

    monkeypatch.setattr(steric.subprocess, "run", run)
    steric._run_probe("input.pdb", "probe")


@pytest.mark.xfail(strict=True, reason="condensed field 5 is dot count; mingap is field 6")
def test_probe_condensed_overlap(monkeypatch):
    raw = ":1->1:bo: A   1 LIG  C   : a   1 LIG  H3  :11:-0.848:-0.476:0:0:0:0:0:C:H:0:0:0:0:0\n"
    monkeypatch.setattr(steric, "_run_probe", lambda *a: raw)
    findings = steric.clashes_probe("input.pdb", "probe", clash_warn_a=0.4)
    assert findings is not None and len(findings) == 1
    assert findings[0].extra["overlap_angstrom"] == pytest.approx(0.848)


@pytest.mark.parametrize("failure", ["timeout", "nonzero", "no-output"])
def test_xtb_failures_return_none(monkeypatch, molecule, failure):
    import subprocess

    def run(*args, **kwargs):
        if failure == "timeout":
            raise TimeoutExpired("xtb", 3600)
        return CompletedProcess(args, int(failure == "nonzero"), "", "failure")

    monkeypatch.setattr(subprocess, "run", run)
    assert refine._run_xtb(*molecule, "xtb", 2, False) is None


def test_missing_xtb_preserves_positions(monkeypatch, molecule):
    monkeypatch.setattr(refine, "_find_binary", lambda name: None)
    assert refine.refine_with_xtb(*molecule) is molecule[1]


def test_large_protein_without_heterogens_switches_without_running_xtb(monkeypatch, molecule, capsys):
    from types import SimpleNamespace

    empty_top = Topology()
    top = SimpleNamespace(atoms=lambda: iter(range(5001)))
    monkeypatch.setattr(refine, "_find_binary", lambda name: "xtb")
    monkeypatch.setattr(refine, "_extract_heterogen_subsystem",
                        lambda *a: (empty_top, molecule[1], {}, set()))
    monkeypatch.setattr(refine, "_run_xtb", lambda *a, **kw: pytest.fail("empty subsystem must not run"))
    assert refine.refine_with_xtb(top, molecule[1]) is molecule[1]
    assert "auto-switching" in capsys.readouterr().out


def test_xtb_merge_restores_only_the_frozen_mask(monkeypatch, molecule):
    top, positions = molecule
    shifted = Quantity([Vec3(.1, .1, .1), Vec3(.24, .1, .1)], nanometer)
    monkeypatch.setattr(refine, "_find_binary", lambda name: "xtb")
    monkeypatch.setattr(refine, "_extract_heterogen_subsystem",
                        lambda *a: (top, positions, {0: 0, 1: 1}, {0}))
    monkeypatch.setattr(refine, "_run_xtb", lambda *a, **kw: shifted)
    result = refine.refine_with_xtb(top, positions, heterogens_only=True)
    assert result[0] == positions[0]
    assert result[1] == shifted[1]


@pytest.mark.parametrize("xyz", [
    "1\nwrong count\nC 0 0 0\n",
    "2\nswapped elements\nO 0 0 0\nC 1.4 0 0\n",
    "2\nnonfinite\nC nan 0 0\nO 1.4 0 0\n",
])
@pytest.mark.xfail(strict=True, reason="xTB output count/elements/finiteness are not validated")
def test_xtb_rejects_invalid_output(monkeypatch, molecule, xyz):
    import subprocess
    from pathlib import Path

    def run(*args, **kwargs):
        Path("xtbopt.xyz").write_text(xyz)
        return CompletedProcess(args, 0, "", "")

    monkeypatch.setattr(subprocess, "run", run)
    assert refine._run_xtb(*molecule, "xtb", 2, False) is None


@pytest.mark.xfail(strict=True, reason="exit zero and xtbopt.xyz suffice without convergence evidence")
def test_xtb_rejects_unconverged_output(monkeypatch, molecule):
    import subprocess
    from pathlib import Path

    def run(*args, **kwargs):
        Path("xtbopt.xyz").write_text("2\nunconverged\nC 0 0 0\nO 1.4 0 0\n")
        return CompletedProcess(args, 0, "FAILED TO CONVERGE", "")

    monkeypatch.setattr(subprocess, "run", run)
    assert refine._run_xtb(*molecule, "xtb", 2, False) is None


@pytest.mark.parametrize("count,heterogen", [(5000, False), (5001, False), (5001, True)])
def test_switch_boundary_with_real_topology(monkeypatch, capsys, count, heterogen):
    top = Topology()
    chain = top.addChain("a")
    protein = top.addResidue("ALA", chain, "1")
    for i in range(count - int(heterogen)):
        top.addAtom(f"C{i}", Element.getBySymbol("C"), protein)
    if heterogen:
        ligand = top.addResidue("LIG", chain, "2", insertionCode="A")
        top.addAtom("C1", Element.getBySymbol("C"), ligand)
    positions = Quantity([Vec3(i * .2, 0, 0) for i in range(count)], nanometer)
    monkeypatch.setattr(refine, "_find_binary", lambda name: "xtb")
    calls = []

    def run(subtop, subpos, *args, **kwargs):
        calls.append(subtop.getNumAtoms())
        return subpos

    monkeypatch.setattr(refine, "_run_xtb", run)
    result = refine.refine_with_xtb(top, positions)
    assert len(result) == count
    assert calls == ([count] if count == 5000 else [1] if heterogen else [])
    assert ("auto-switching" in capsys.readouterr().out) == (count > 5000)
