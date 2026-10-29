"""Unavailable engines, failed calculations and legitimate negative findings differ."""

from __future__ import annotations

import subprocess
from pathlib import Path
from types import SimpleNamespace

import pytest

from dvbfixer import antibody

PDB = ("ATOM      1  N   ALA a   1       0.000   0.000   0.000  1.00  0.00           N\n"
       "ATOM      2  CA  ALA a   1       1.450   0.000   0.000  1.00  0.00           C\n"
       "ATOM      3  C   ALA a   1       2.000   1.400   0.000  1.00  0.00           C\nEND\n")


@pytest.fixture
def anarci_ready(tmp_path, monkeypatch):
    database = tmp_path / "HMMs"
    database.mkdir()
    for suffix in ("", ".h3f", ".h3i", ".h3m", ".h3p"):
        (database / ("ALL.hmm" + suffix)).write_text("fixture database")
    module = SimpleNamespace(HMM_path=str(database), anarci=lambda *a, **kw: ([None], [None], []))
    monkeypatch.setattr(antibody, "importlib", SimpleNamespace(import_module=lambda name: module))
    monkeypatch.setattr(antibody, "shutil", SimpleNamespace(which=lambda name: "/env/bin/hmmscan"))
    return module


def test_anarci_absence_is_not_no_antibody_domain(monkeypatch):
    def unavailable(name):
        raise ImportError("missing")
    monkeypatch.setattr(antibody, "importlib", SimpleNamespace(import_module=unavailable))
    with pytest.raises(antibody.AntibodyToolError, match="ANARCI is unavailable"):
        antibody._number_v_domain("ACDEFG", "imgt")


def test_hmmer_and_pressed_database_are_required(anarci_ready, monkeypatch):
    monkeypatch.setattr(antibody.shutil, "which", lambda name: None)
    with pytest.raises(antibody.AntibodyToolError, match="hmmscan"):
        antibody.require_anarci()
    monkeypatch.setattr(antibody.shutil, "which", lambda name: "/env/bin/hmmscan")
    (Path(anarci_ready.HMM_path) / "ALL.hmm.h3i").unlink()
    with pytest.raises(antibody.AntibodyToolError, match="database is missing or incomplete"):
        antibody.require_anarci()


def test_anarci_execution_failure_does_not_become_negative_domain(anarci_ready):
    def failed(*args, **kwargs):
        raise RuntimeError("hmmscan exited 2")
    anarci_ready.anarci = failed
    with pytest.raises(antibody.AntibodyToolError, match="execution failed.*no antibody-domain decision"):
        antibody._number_v_domain("ACDEFG", "imgt")


def test_successful_no_domain_is_still_a_legitimate_negative(anarci_ready):
    assert antibody._number_v_domain("ACDEFG", "imgt") is None


def test_explicit_antibody_scheme_failure_writes_no_seqres_fallback(tmp_path, monkeypatch):
    from dvbfixer import renumber

    source, output = tmp_path / "input.pdb", tmp_path / "output.pdb"
    source.write_text(PDB)

    def failed():
        raise antibody.AntibodyToolError("hmmscan unavailable")

    monkeypatch.setattr(antibody, "require_anarci", failed)
    with pytest.raises(SystemExit, match="hmmscan unavailable"):
        renumber.main([str(source), "--scheme", "imgt", "-o", str(output)])
    assert not output.exists()
    assert source.read_text() == PDB


@pytest.mark.parametrize("propka_state,reduce_state", [
    ("disabled", "disabled"), ("completed", "completed"),
    ("unavailable", "unavailable"), ("failed", "failed"),
])
def test_prepare_records_disabled_unavailable_failed_and_completed_separately(
    tmp_path, monkeypatch, capsys, propka_state, reduce_state,
):
    pytest.importorskip("openmm")
    from dvbfixer import prep_backend, protonate
    from dvbfixer.prepare.pipeline import _run_propka_reduce_variants

    source = tmp_path / "input.pdb"
    source.write_text(PDB)

    def propka(path):
        if propka_state == "unavailable":
            raise protonate.PropkaUnavailableError("not installed")
        if propka_state == "failed":
            raise ValueError("pKa calculation error")
        if propka_state == "disabled":
            raise AssertionError("disabled PROPKA was called")
        return object()

    def reduce(input_path, output_path, **kwargs):
        if reduce_state == "unavailable":
            raise prep_backend.MissingPrepToolError("reduce absent")
        if reduce_state == "failed":
            raise prep_backend.ReduceError("reduce calculation error")
        if reduce_state == "disabled":
            raise AssertionError("disabled Reduce was called")
        Path(output_path).write_text(Path(input_path).read_text())

    monkeypatch.setattr(protonate, "run_propka", propka)
    monkeypatch.setattr(protonate, "get_pka_results", lambda container: [])
    monkeypatch.setattr(prep_backend, "run_reduce", reduce)
    statuses = {}
    assert _run_propka_reduce_variants(source, 7.0, set(), use_propka=propka_state != "disabled",
                                      use_reduce=reduce_state != "disabled", stage_status=statuses) == {}
    assert statuses == {"propka": propka_state, "reduce": reduce_state}
    output = capsys.readouterr().out
    assert ("WARNING:" in output) == (propka_state in {"unavailable", "failed"})


def test_tleap_requires_reduce_before_starting_tleap(tmp_path, monkeypatch):
    from dvbfixer import prep_backend

    source, output = tmp_path / "input.pdb", tmp_path / "output.pdb"
    source.write_text(PDB)
    monkeypatch.setattr(prep_backend.shutil, "which", lambda name: "/env/bin/tleap" if name == "tleap" else None)
    monkeypatch.setattr(prep_backend, "run_tleap", lambda *a, **kw: pytest.fail("tleap started without Reduce"))
    with pytest.raises(prep_backend.MissingPrepToolError, match="MolProbity Reduce separately"):
        prep_backend.run_prep(source, output)
    assert not output.exists()


def test_noncanonical_tleap_input_fails_before_external_work(tmp_path, monkeypatch):
    from dvbfixer import prep_backend

    source, output = tmp_path / "input.pdb", tmp_path / "output.pdb"
    source.write_text(PDB.replace("ALA", "NAG"))
    monkeypatch.setattr(prep_backend, "_require_binary", lambda *a, **kw: pytest.fail("tool lookup before chemistry check"))
    with pytest.raises(prep_backend.PrepBackendError, match="a/NAG1.*--backend legacy"):
        prep_backend.run_prep(source, output)
    assert not output.exists()


def test_standalone_protonate_requires_reduce_before_propka_calculation(tmp_path, monkeypatch):
    from dvbfixer import protonate

    source, output = tmp_path / "input.pdb", tmp_path / "output.pdb"
    source.write_text(PDB)
    monkeypatch.setattr(protonate, "_find_reduce_binary", lambda *a: None)
    monkeypatch.setattr(protonate, "run_propka", lambda *a: pytest.fail("PROPKA started without required Reduce"))
    with pytest.raises(SystemExit, match="Reduce is unavailable"):
        protonate.main([str(source), "-o", str(output)])
    assert not output.exists()


@pytest.mark.parametrize("code,contents,success", [(1, PDB, True), (2, PDB, False), (0, "usage only\n", False)])
def test_reduce_exit_and_pdb_output_both_matter(tmp_path, monkeypatch, code, contents, success):
    import os
    import tempfile

    from dvbfixer import protonate

    fd, output = tempfile.mkstemp(dir=tmp_path)
    monkeypatch.setattr(tempfile, "mkstemp", lambda **kw: (fd, output))

    def run(argv, **kwargs):
        kwargs["stdout"].write(contents)
        return subprocess.CompletedProcess(argv, code, stderr="native diagnostic")

    monkeypatch.setattr(subprocess, "run", run)
    if success:
        assert protonate._run_reduce("input.pdb", "/env/reduce") == output
        assert Path(output).read_text() == PDB
    else:
        with pytest.raises(RuntimeError):
            protonate._run_reduce("input.pdb", "/env/reduce")
        assert not Path(output).exists()
    with pytest.raises(OSError):
        os.fstat(fd)


def test_trajectory_mismatch_is_actionable_before_analysis(tmp_path):
    mda = pytest.importorskip("MDAnalysis")
    from dvbfixer.trajectory_input import TrajectoryInputError, load_universe

    topology = tmp_path / "three.pdb"
    topology.write_text(PDB)
    other = tmp_path / "two.pdb"
    other.write_text("\n".join(PDB.splitlines()[:2]) + "\nEND\n")
    universe = mda.Universe(str(other))
    universe.dimensions = [10, 10, 10, 90, 90, 90]
    trajectory = tmp_path / "two atoms.dcd"
    with mda.Writer(str(trajectory), n_atoms=2) as writer:
        writer.write(universe.atoms)
    with pytest.raises(TrajectoryInputError, match="matching atom counts/order"):
        load_universe(topology, trajectory)
    valid = load_universe(other, trajectory)
    assert valid.atoms.n_atoms == 2
    valid.trajectory.close()


def test_unsupported_tpr_has_reader_version_and_alternative_topology_advice(tmp_path):
    pytest.importorskip("MDAnalysis")
    from dvbfixer.trajectory_input import TrajectoryInputError, load_universe

    topology = tmp_path / "unsupported.tpr"
    topology.write_bytes(b"not a supported TPR header")
    with pytest.raises(TrajectoryInputError, match="TPR writer version.*GRO/PDB"):
        load_universe(topology)


def test_psi4_environment_command_preserves_named_route_and_supports_prefix(tmp_path):
    from dvbfixer.parametrize import _psi4_environment_command

    assert _psi4_environment_command('conda', 'psi4', ['python', 'worker.py']) == [
        'conda', 'run', '-n', 'psi4', 'python', 'worker.py']
    prefix = tmp_path / 'QM environment'
    assert _psi4_environment_command('mamba', str(prefix), ['python']) == [
        'mamba', 'run', '-p', str(prefix), 'python']


def test_missing_psi4_environment_fails_before_parameterization(monkeypatch):
    from dvbfixer import parametrize

    monkeypatch.setattr(parametrize, '_find_env_runner', lambda: 'micromamba')
    monkeypatch.setattr(subprocess, 'run', lambda *a, **kw: subprocess.CompletedProcess(a[0], 1, '', 'private env failure'))
    with pytest.raises(RuntimeError, match='PSI4/psiresp cannot be imported'):
        parametrize.preflight_qm_engine('psi4', psi4_env='missing')
