"""Production Probe contracts: validated input/output and visible fallback."""
from subprocess import CompletedProcess, TimeoutExpired

import pytest

pytest.importorskip("openmm")

from openmm import Vec3
from openmm.app import Element, PDBFile, Topology
from openmm.unit import Quantity, nanometer

from dvbfixer.diagnose import probe, steric
from dvbfixer.diagnose.geometry import AtomIdentity


@pytest.fixture
def structure(tmp_path):
    top = Topology()
    for chain, resid, icode, name, atom in [("H", "82", "a", "0fA", "c1"), ("h", "-1", "B", "LIG", "C2")]:
        residue = top.addResidue(name, top.addChain(chain), resid, insertionCode=icode)
        top.addAtom(atom, Element.getBySymbol("C"), residue)
    positions = Quantity([Vec3(0, 0, 0), Vec3(.2, 0, 0)], nanometer)
    path = tmp_path / "original.pdb"
    with path.open("w") as handle:
        PDBFile.writeFile(top, positions, handle, keepIds=True)
    return top, positions, path


def record(kind="bo", gap="-0.848", source=None):
    source = source or " H  82a0FA  C1  "
    return ":".join(["", "1->1", kind, source, " h  -1BLIG  C2  ", "11", gap, "-0.476",
                     "0", "0", "0", "0", "0", "C", "C", "0", "0", "0", "0", "0"]) + "\n"


def native(monkeypatch, raw="", failure=None, version=probe.SUPPORTED_VERSION):
    calls = []

    def run(command, **kwargs):
        if command[-1] == "-version":
            return CompletedProcess(command, 0, version + "\n", "")
        from pathlib import Path
        calls.append((command, Path(command[-1]).read_text()))
        if failure == "missing":
            raise FileNotFoundError("probe")
        if failure == "timeout":
            raise TimeoutExpired(command, 60)
        if failure == "encoding":
            raise UnicodeError("invalid UTF-8")
        return CompletedProcess(command, int(failure == "nonzero"), raw, "")

    monkeypatch.setattr(probe.subprocess, "run", run)
    monkeypatch.setattr(steric, "_find_binary", lambda name: "probe")
    return calls


def test_invocation_normalization_identity_and_pair_aggregation(monkeypatch, structure):
    top, positions, path = structure
    path.write_text("".join(line[:12] + f"{line[12:16].strip():<4}" + line[16:]
                            if line.startswith(("ATOM  ", "HETATM")) else line
                            for line in path.read_text().splitlines(keepends=True)))
    original = path.read_bytes()
    calls = native(monkeypatch, record() + record("so", "-0.7") + record("hb", "-0.9"))
    provenance = {}
    findings = steric.run_all(top, positions, path, clash_warn_a=.4, clash_error_a=.5,
                              probe_radii="nuclear", provenance=provenance)
    assert len(findings) == 1
    finding = findings[0]
    assert finding.extra["overlap_angstrom"] == pytest.approx(.848)
    assert finding.extra["contact_types"] == ["bo", "hb", "so"]
    assert finding.extra["atoms"] == [AtomIdentity("H", "82", "a", "c1").to_dict(),
                                      AtomIdentity("h", "-1", "B", "C2").to_dict()]
    assert finding.resname == "0fA" and finding.atom == "c1" and finding.resid == "82a"
    command, normalized = calls[0]
    assert command[command.index("-self") + 1:-1] == ["ALL"]
    assert "-nuclear" in command and "-nowaters" in command
    assert normalized != path.read_text()  # lowercase atom spelling survives; field alignment changes
    assert path.read_bytes() == original
    assert provenance["actual_engine"] == "molprobity-probe" and not provenance["fallback"]
    assert provenance["probe_version"] == probe.SUPPORTED_VERSION


@pytest.mark.parametrize("raw", ["unsupported\n", record() + "bad\n", record(gap="nan"),
                                  record(gap="inf"), record(source=" H  82A0FA  C1  "),
                                  record().replace(":11:", ":0:")])
def test_invalid_output_falls_back_instead_of_becoming_clean(monkeypatch, structure, raw, capsys):
    native(monkeypatch, raw)
    metadata = {}
    findings = steric.run_all(*structure, provenance=metadata)
    assert len(findings) == 1 and findings[0].extra["engine"]["name"] == "dvbfixer-python-steric"
    assert metadata["actual_engine"] == "dvbfixer-python-steric" and metadata["fallback"]
    assert metadata["fallback_reason"]
    assert "WARNING" in capsys.readouterr().err


@pytest.mark.parametrize("failure", ["missing", "timeout", "nonzero", "encoding"])
def test_native_failures_fall_back(monkeypatch, structure, failure):
    native(monkeypatch, failure=failure)
    metadata = {}
    assert len(steric.run_all(*structure, provenance=metadata)) == 1
    assert metadata["probe_status"] == "failed" and metadata["fallback"]


def test_unknown_version_falls_back(monkeypatch, structure):
    calls = native(monkeypatch, version="probe.9.99")
    metadata = {}
    steric.run_all(*structure, provenance=metadata)
    assert not calls and metadata["fallback"]
    assert metadata["probe_version"] == "probe.9.99"


def test_clean_result_keeps_engine_metadata_without_warning(monkeypatch, structure, capsys):
    native(monkeypatch)
    metadata = {}
    assert steric.run_all(*structure, provenance=metadata) == []
    assert metadata["actual_engine"] == "molprobity-probe" and metadata["status"] == "complete"
    assert capsys.readouterr().err == ""


def test_blank_chain_matches_loaded_pdb_identity(monkeypatch, structure):
    top, positions, path = structure
    path.write_text("".join(line[:21] + " " + line[22:]
                            if line.startswith(("ATOM  ", "HETATM")) and line[21] == "H" else line
                            for line in path.read_text().splitlines(keepends=True)))
    native(monkeypatch, record().replace(" H  82a", "    82a"))
    loaded = PDBFile(str(path))
    findings = steric.run_all(loaded.topology, loaded.positions, path, clash_warn_a=.4)
    assert findings[0].chain == next(loaded.topology.chains()).id == " "
    assert findings[0].extra["atoms"][0]["chain"] == " "


def test_optional_missing_binary_is_visible_but_not_a_warning(monkeypatch, structure, capsys):
    monkeypatch.setattr(steric, "_find_binary", lambda name: None)
    metadata = {}
    steric.run_all(*structure, provenance=metadata)
    assert metadata["probe_status"] == "unavailable" and metadata["actual_engine"] == "dvbfixer-python-steric"
    assert capsys.readouterr().err == ""


def test_water_selection_and_conect_integrity(monkeypatch, structure):
    top, positions, path = structure
    water = top.addResidue("SOL", top.addChain("W"), "3")
    top.addAtom("O", Element.getBySymbol("O"), water)
    positions = Quantity([Vec3(0, 0, 0), Vec3(.2, 0, 0), Vec3(.4, 0, 0)], nanometer)
    with path.open("w") as handle:
        PDBFile.writeFile(top, positions, handle, keepIds=True)
    atom_lines = [line for line in path.read_text().splitlines() if line.startswith(("ATOM  ", "HETATM"))]
    serial = atom_lines[-1][6:11]
    with path.open("a") as handle:
        handle.write(f"CONECT{atom_lines[0][6:11]}{serial}\n")
    calls = native(monkeypatch)
    steric.run_all(top, positions, path)
    assert "SOL" not in calls[0][1] and "CONECT" not in calls[0][1]
    steric.run_all(top, positions, path, include_water=True)
    assert "SOL" in calls[1][1] and "-wat2wat" in calls[1][0] and "-nowaters" not in calls[1][0]


def test_ambiguous_case_mapping_rejected_before_native_call(monkeypatch, structure):
    top, positions, path = structure
    lines = path.read_text().splitlines(keepends=True)
    first = next(line for line in lines if line.startswith(("ATOM  ", "HETATM")))
    path.write_text("".join(lines) + first[:12] + " C1 " + first[16:])
    calls = native(monkeypatch)
    result = probe.analyze(path, "probe", include_water=False, radii="electron", warn_a=.4, error_a=.5)
    assert result.status == "unsupported_input" and result.findings is None and not calls


@pytest.mark.parametrize("fmt", ["json", "text"])
def test_clean_cli_report_exposes_actual_engine(monkeypatch, structure, capsys, fmt):
    import json

    from dvbfixer.diagnose import pipeline
    native(monkeypatch)
    top, positions, path = structure
    monkeypatch.setattr(pipeline, "_load_topology_and_fixer", lambda *a: (top, positions, None))
    with pytest.raises(SystemExit) as exc:
        pipeline.main([str(path), "--only", "steric", "--format", fmt])
    assert exc.value.code == 0
    output = capsys.readouterr().out
    if fmt == "json":
        report = json.loads(output)
        assert report["steric"]["actual_engine"] == "molprobity-probe" and report["findings"] == []
    else:
        assert "Steric engine: molprobity-probe" in output and "No findings" in output
