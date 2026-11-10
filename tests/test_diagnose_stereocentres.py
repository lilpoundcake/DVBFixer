"""Warning-only stereo/shape contracts, independent of reconstruction research."""

from __future__ import annotations

import json
import math
import os
import subprocess
import sys
from pathlib import Path

import pytest

from dvbfixer.diagnose.geometry import stereocentre_geometry


@pytest.mark.parametrize("angle", [-55.0, -10.01, -9.99, 0.0, 9.99, 10.01, 55.0])
def test_ordered_angle_and_volume(angle):
    radians = math.radians(angle)
    result = stereocentre_geometry(
        (0, 0, 0),
        (2, 0, 0),
        (0, 3, 0),
        (4 * math.cos(radians), 0, 4 * math.sin(radians)),
    )
    assert result is not None
    assert result.out_of_plane_degrees == pytest.approx(angle)
    assert result.triple_product_angstrom3 == pytest.approx(24 * math.sin(radians))


def test_tetrahedral_control_rigid_transform_scaling_and_reflection():
    points = [(0, 0, 0), (1, 1, 1), (1, -1, -1), (-1, 1, -1)]
    ideal = stereocentre_geometry(*points)
    assert ideal is not None
    assert ideal.out_of_plane_degrees == pytest.approx(math.degrees(math.asin(math.sqrt(2 / 3))))
    transformed = [(3 - 2 * y, 4 + 2 * x, 5 + 2 * z) for x, y, z in points]
    measured = stereocentre_geometry(*transformed)
    assert measured is not None
    assert measured.out_of_plane_degrees == pytest.approx(ideal.out_of_plane_degrees)
    assert measured.triple_product_angstrom3 == pytest.approx(8 * ideal.triple_product_angstrom3)
    reflected = stereocentre_geometry(*[(-x, y, z) for x, y, z in points])
    assert reflected is not None
    assert reflected.out_of_plane_degrees == pytest.approx(-ideal.out_of_plane_degrees)


@pytest.mark.parametrize(
    "points",
    [
        [(0, 0, 0), (0, 0, 0), (0, 1, 0), (0, 0, 1)],
        [(0, 0, 0), (1, 0, 0), (1, 0.01, 0), (0, 0, 1)],
        [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 0)],
        [(math.nan, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)],
        [(0, 0, 0), (math.inf, 0, 0), (0, 1, 0), (0, 0, 1)],
        [(0, 0, 0), (1e200, 0, 0), (0, 1e200, 0), (0, 0, 1e200)],
    ],
)
def test_undefined_stereocentre_is_not_assigned(points):
    assert stereocentre_geometry(*points) is None


def _topology(points, name="THR", chain="d", icode="A"):
    pytest.importorskip("openmm")
    from openmm.app import Element, Topology
    from openmm.unit import Quantity, angstrom

    topology = Topology()
    residue = topology.addResidue(name, topology.addChain(chain), id="82", insertionCode=icode)
    for atom in points:
        topology.addAtom(atom, Element.getBySymbol(atom[0]), residue)
    return topology, Quantity(list(points.values()), angstrom)


def _thr(sign=1):
    # CA stays clearly L while the sign of CB alone changes.
    return dict(
        CB=(0, 0, 0),
        CA=(1.5, 0, 0),
        OG1=(-0.45, 1.3, 0),
        CG2=(-0.5, -0.7, sign * 1.1),
        N=(1.5, 0, 1.45),
        C=(1.5, 1.5, 0),
    )


def _alpha(angle):
    radius = 1.53
    return dict(
        CA=(0, 0, 0),
        N=(1.45, 0, 0),
        C=(-0.5, 1.4, 0),
        CB=(-radius * math.cos(math.radians(angle)), 0, radius * math.sin(math.radians(angle))),
    )


def test_thr_beta_warning_has_exact_identity_units_and_no_ca_error():
    from dvbfixer.diagnose.chemistry import check_ca_chirality, check_stereocentre_geometry
    from dvbfixer.diagnose.report import Severity

    assert check_stereocentre_geometry(*_topology(_thr())) == []
    top, positions = _topology(_thr(-1))
    assert check_ca_chirality(top, positions) == []
    before = positions[:]
    (finding,) = check_stereocentre_geometry(top, positions)
    assert finding.category == "thr_beta_chirality" and finding.severity == Severity.WARNING
    assert (finding.chain, finding.resid, finding.atom) == ("d", "82A", "CB")
    assert finding.extra["observed_angstrom3"] < 0
    assert finding.extra["observed_degrees"] < -10
    assert finding.extra["atoms"][0] == dict(chain="d", resid="82", icode="A", atom="CB")
    assert positions == before


@pytest.mark.parametrize(
    "angle,warn",
    [(-10.01, False), (-9.99, True), (0, True), (9.99, True), (10.01, False), (55, False)],
)
def test_ca_planarity_threshold_is_sign_symmetric_and_warning_only(angle, warn):
    from dvbfixer.diagnose.chemistry import check_ca_chirality, check_stereocentre_geometry
    from dvbfixer.diagnose.report import Severity

    top, positions = _topology(_alpha(angle), name="ALA")
    findings = check_stereocentre_geometry(top, positions)
    assert bool(findings) is warn
    assert all(f.category == "ca_planarity" and f.severity == Severity.WARNING for f in findings)
    assert bool(check_ca_chirality(top, positions)) is (angle < 0)


@pytest.mark.parametrize("name", ["VAL", "ILE", "TPO", "UNK", "LIG"])
def test_noncanonical_thr_and_nonstereogenic_val_are_not_assigned_beta(name):
    from dvbfixer.diagnose.chemistry import check_stereocentre_geometry

    assert check_stereocentre_geometry(*_topology(_thr(-1), name=name)) == []


@pytest.mark.parametrize(
    "name,expected",
    [
        ("ALA", True),
        ("HIE", True),
        ("GLUP", True),
        ("MSE", True),
        ("TPO", True),
        ("GLY", False),
        ("UNK", False),
        ("LIG", False),
    ],
)
def test_planarity_recognizes_only_protein_and_shared_variants(name, expected):
    from dvbfixer.diagnose.chemistry import check_stereocentre_geometry

    assert bool(check_stereocentre_geometry(*_topology(_alpha(2), name=name))) is expected


@pytest.mark.parametrize(
    "fault",
    ["missing", "duplicate", "element", "nan", "coincident", "collinear", "duplicate_residue"],
)
def test_ambiguous_or_undefined_required_atoms_are_not_guessed(fault):
    pytest.importorskip("openmm")
    from openmm.app import Element
    from openmm.unit import angstrom

    from dvbfixer.diagnose.chemistry import check_stereocentre_geometry

    points = _thr(-1)
    if fault == "missing":
        del points["OG1"]
    top, positions = _topology(points)
    residue = next(top.residues())
    if fault == "duplicate":
        top.addAtom("OG1", Element.getBySymbol("O"), residue)
        positions.append((1, 0, 0) * angstrom)
    elif fault == "duplicate_residue":
        top.addResidue("THR", residue.chain, id="82", insertionCode="A")
    elif fault == "element":
        next(a for a in top.atoms() if a.name == "OG1").element = Element.getBySymbol("C")
    elif fault in {"nan", "coincident", "collinear"}:
        atom = next(a for a in top.atoms() if a.name == "OG1")
        point = {"nan": (math.nan, 0, 0), "coincident": (0, 0, 0), "collinear": (1, 0, 0)}[fault]
        positions[atom.index] = point * angstrom
    assert check_stereocentre_geometry(top, positions) == []


@pytest.mark.parametrize("identifier", ["1UBQ", "1A3N", "1TIM"])
def test_native_fixture_controls_have_no_new_stereo_warnings(identifier):
    pytest.importorskip("openmm")
    from openmm.app import PDBFile

    from dvbfixer.diagnose.chemistry import check_stereocentre_geometry

    path = Path(__file__).parent / "fixtures" / "protein_context" / f"{identifier}.pdb"
    original = path.read_bytes()
    pdb = PDBFile(str(path))
    assert check_stereocentre_geometry(pdb.topology, pdb.positions) == []
    assert path.read_bytes() == original


@pytest.mark.parametrize("angle,scale,warn", [(55, 0.01, False), (3, 100, True)])
def test_planarity_is_not_an_absolute_volume_gate(angle, scale, warn):
    from dvbfixer.diagnose.chemistry import check_stereocentre_geometry

    points = {name: tuple(scale * x for x in xyz) for name, xyz in _alpha(angle).items()}
    assert bool(check_stereocentre_geometry(*_topology(points, name="ALA"))) is warn


@pytest.mark.parametrize("angle,warn", [(-10.01, True), (-9.99, False), (0, False), (55, False)])
def test_nearly_flat_thr_beta_is_not_assigned_an_inversion(angle, warn):
    from dvbfixer.diagnose.chemistry import check_stereocentre_geometry

    points = _thr()
    points["CG2"] = (-math.cos(math.radians(angle)), 0, math.sin(math.radians(angle)))
    findings = check_stereocentre_geometry(*_topology(points))
    assert bool(findings) is warn
    assert all(f.category == "thr_beta_chirality" for f in findings)


@pytest.mark.parametrize("fault", ["missing", "duplicate", "nan", "collinear"])
def test_ca_geometry_does_not_guess_missing_or_ambiguous_atoms(fault):
    pytest.importorskip("openmm")
    from openmm.app import Element
    from openmm.unit import angstrom

    from dvbfixer.diagnose.chemistry import check_stereocentre_geometry

    points = _alpha(2)
    if fault == "missing":
        del points["N"]
    elif fault == "nan":
        points["CB"] = (math.nan, 0, 0)
    elif fault == "collinear":
        points["C"] = (-1, 0, 0)
    top, positions = _topology(points, name="ALA")
    if fault == "duplicate":
        top.addAtom("CA", Element.getBySymbol("C"), next(top.residues()))
        positions.append((0, 0, 0) * angstrom)
    assert check_stereocentre_geometry(top, positions) == []


def test_native_thr_single_branch_reflection_detected_without_moving_ca():
    pytest.importorskip("openmm")
    from openmm.app import PDBFile
    from openmm.unit import Quantity, angstrom

    from dvbfixer.diagnose.chemistry import check_ca_chirality, check_stereocentre_geometry

    pdb = PDBFile(str(Path(__file__).parent / "fixtures/protein_context/1UBQ.pdb"))
    residue = next(r for r in pdb.topology.residues() if r.name == "THR")
    atoms = {a.name: a for a in residue.atoms()}
    points = [tuple(p) for p in pdb.positions.value_in_unit(angstrom)]
    centre, first, second, branch = [points[atoms[n].index] for n in ("CB", "CA", "OG1", "CG2")]
    a, b, c = [tuple(x - y for x, y in zip(p, centre)) for p in (first, second, branch)]
    normal = (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])
    coefficient = 2 * sum(x * y for x, y in zip(normal, c)) / sum(x * x for x in normal)
    points[atoms["CG2"].index] = tuple(x - coefficient * y for x, y in zip(branch, normal))
    positions = Quantity(points, angstrom)
    (finding,) = check_stereocentre_geometry(pdb.topology, positions)
    assert finding.category == "thr_beta_chirality"
    assert finding.resid == residue.id
    assert check_ca_chirality(pdb.topology, positions) == check_ca_chirality(
        pdb.topology, pdb.positions
    )
    assert tuple(positions[atoms["CA"].index].value_in_unit(angstrom)) == tuple(
        pdb.positions[atoms["CA"].index].value_in_unit(angstrom)
    )


def _write_pdb(path, top, positions):
    from openmm.app import PDBFile

    with path.open("w") as stream:
        PDBFile.writeFile(top, positions, stream, keepIds=True)


@pytest.mark.parametrize("kind", ["beta", "planarity"])
def test_text_json_severity_scope_and_input_preservation(tmp_path, capsys, kind):
    pytest.importorskip("pdbfixer")
    from dvbfixer.diagnose.pipeline import main

    source = tmp_path / "стерео.pdb"
    top, positions = _topology(_thr(-1)) if kind == "beta" else _topology(_alpha(3), name="ALA")
    _write_pdb(source, top, positions)
    original = source.read_bytes()
    category = "thr_beta_chirality" if kind == "beta" else "ca_planarity"

    def run(options):
        try:
            main([str(source), "--only", "chemistry", *options])
        except SystemExit as error:
            assert error.code == 0
        return capsys.readouterr().out

    text = run([])
    assert "D-isomer error: NO" in text and "WARNING" in text
    assert "Chemistry / bond geometry" in text and "no automatic repair" in text
    data = json.loads(
        run(["--format", "json", "--scope-residue", "d:82:A", "--scope-residue", "D:82:A"])
    )
    assert data["chirality"]["d_isomer_error"] is False
    warning = next(f for f in data["findings"] if f["category"] == category)
    assert warning["severity"] == "WARNING"
    details = warning["details"]
    assert details["geometry_schema"] == "dvbfixer.geometry.v1"
    assert details["planarity_threshold_degrees"] == 10
    assert math.isfinite(details["observed_degrees"])
    assert math.isfinite(details["observed_angstrom3"])
    assert all(a["chain"] == "d" and a["icode"] == "A" for a in details["atoms"])
    assert data["scope"]["finding_aggregates"][category]["count"] == 1
    maxima = data["scope"]["finding_aggregates"][category]["maxima"]
    assert maxima["absolute_triple_product_angstrom3"] == abs(details["observed_angstrom3"])
    assert maxima["observed_absolute_degrees"] == abs(details["observed_degrees"])
    other = json.loads(
        run(["--format", "json", "--scope-residue", "D:82:A", "--scope-residue", "d:82"])
    )
    assert other["scope"]["finding_aggregates"][category]["count"] == 0
    filtered = json.loads(
        run(["--format", "json", "--severity", "ERROR", "--scope-residue", "d:82:A"])
    )
    assert not filtered["findings"]
    assert filtered["scope"]["finding_aggregates"][category]["count"] == 0
    assert source.read_bytes() == original


@pytest.mark.parametrize("output_format", ["text", "json"])
def test_negative_nearplanarity_does_not_suppress_existing_d_error(tmp_path, capsys, output_format):
    pytest.importorskip("pdbfixer")
    from dvbfixer.diagnose.pipeline import main

    source = tmp_path / "negative.pdb"
    _write_pdb(source, *_topology(_alpha(-3), name="ALA"))
    with pytest.raises(SystemExit) as error:
        main([str(source), "--only", "chemistry", "--format", output_format])
    assert error.value.code == 1
    output = capsys.readouterr().out
    if output_format == "text":
        assert "D-isomer error: YES" in output
        assert "near-planar Cα" in output
    else:
        payload = json.loads(output)
        assert payload["chirality"]["d_isomer_error"] is True
        assert {f["category"] for f in payload["findings"]} >= {"chirality", "ca_planarity"}


@pytest.mark.parametrize("family", ["structural", "steric"])
def test_stereo_warnings_are_not_run_in_other_families(tmp_path, capsys, family):
    pytest.importorskip("pdbfixer")
    from dvbfixer.diagnose.pipeline import main

    source = tmp_path / "wrong-thr.pdb"
    _write_pdb(source, *_topology(_thr(-1)))
    with pytest.raises(SystemExit) as error:
        main([str(source), "--only", family, "--format", "json"])
    assert error.value.code in (0, 1)
    payload = json.loads(capsys.readouterr().out)
    assert payload["chirality"]["checked"] is False
    assert not any(
        f["category"] in {"ca_planarity", "thr_beta_chirality"} for f in payload["findings"]
    )


@pytest.mark.parametrize("output_format", ["text", "json"])
@pytest.mark.parametrize(
    "kind,severity,expected_category,expected_exit",
    [
        ("beta", "WARNING", "thr_beta_chirality", 0),
        ("planarity", "WARNING", "ca_planarity", 0),
        ("negative", "WARNING", "ca_planarity", 1),
        ("beta", "ERROR", None, 0),
    ],
)
def test_public_cli_stereo_logging_and_exit_contract(
    tmp_path, output_format, kind, severity, expected_category, expected_exit
):
    pytest.importorskip("pdbfixer")
    points = _thr(-1) if kind == "beta" else _alpha(-3 if kind == "negative" else 3)
    source = tmp_path / "стерео.pdb"
    _write_pdb(source, *_topology(points, name="THR" if kind == "beta" else "ALA"))
    original = source.read_bytes()
    log = tmp_path / "diagnose.log"
    # Exercise the actual unified entry point and fd-level capture in a child,
    # always importing this checkout rather than a different editable install.
    src = str(Path(__file__).resolve().parents[1] / "src")
    env = dict(os.environ, PYTHONPATH=src + os.pathsep + os.environ.get("PYTHONPATH", ""))
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "dvbfixer.cli",
            "diagnose",
            str(source),
            "--only",
            "chemistry",
            "--format",
            output_format,
            "--severity",
            severity,
            "--log-file",
            str(log),
        ],
        env=env,
        capture_output=True,
        text=True,
        timeout=60,
        check=False,
    )
    assert result.returncode == expected_exit, result.stderr + result.stdout
    assert source.read_bytes() == original
    logged = log.read_text()
    assert "\x1b[" not in logged
    assert "DIAGNOSTIC SUMMARY" not in result.stdout
    if expected_exit == 0:
        assert "!!! ERROR !!!" not in result.stdout + result.stderr + logged
    if expected_category is None:
        assert "DIAGNOSTIC SUMMARY" not in result.stderr + logged
    elif output_format == "text" and expected_exit == 0:
        assert "DIAGNOSTIC SUMMARY: 0 error(s)," in result.stderr
    if output_format == "json":
        payload = json.loads(result.stdout)
        assert payload["chirality"]["d_isomer_error"] is (kind == "negative")
        if expected_category is None:
            assert payload["findings"] == []
        else:
            assert any(f["category"] == expected_category for f in payload["findings"])
        # Logging retains the report, not an ASCII-escaped approximation.
        assert "стерео.pdb" in logged
        if expected_category is not None:
            assert "observed_angstrom3" in logged
            assert ("Cβ" if kind == "beta" else "Cα") in logged
    else:
        assert f"D-isomer error: {'YES' if kind == 'negative' else 'NO'}" in result.stdout
        assert ("WARNING" in result.stdout) is (expected_category is not None)
        if expected_category is not None:
            assert "no automatic repair" in logged


def test_public_batch_cli_keeps_json_reports_and_error_diagnostics(tmp_path):
    pytest.importorskip("pdbfixer")
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    _write_pdb(inputs / "beta.pdb", *_topology(_thr(-1)))
    _write_pdb(inputs / "negative.pdb", *_topology(_alpha(-3), name="ALA"))
    originals = {p.name: p.read_bytes() for p in inputs.iterdir()}
    output = tmp_path / "reports"
    log = tmp_path / "batch.log"
    src = str(Path(__file__).resolve().parents[1] / "src")
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "dvbfixer.cli",
            "diagnose",
            "--input-dir",
            str(inputs),
            "--output-dir",
            str(output),
            "--only",
            "chemistry",
            "--format",
            "json",
            "--log-file",
            str(log),
        ],
        env=dict(os.environ, PYTHONPATH=src + os.pathsep + os.environ.get("PYTHONPATH", "")),
        capture_output=True,
        text=True,
        timeout=60,
        check=False,
    )
    assert result.returncode == 1, result.stderr + result.stdout
    assert "0 execution failures" in result.stdout
    assert "DIAGNOSTIC SUMMARY" in result.stderr
    beta = json.loads((output / "beta_diagnose.json").read_text())
    assert beta["chirality"]["d_isomer_error"] is False
    assert any(f["category"] == "thr_beta_chirality" for f in beta["findings"])
    negative = json.loads((output / "negative_diagnose.json").read_text())
    assert negative["chirality"]["d_isomer_error"] is True
    assert any(f["category"] == "ca_planarity" for f in negative["findings"])
    assert {p.name: p.read_bytes() for p in inputs.iterdir()} == originals
    assert "\x1b[" not in log.read_text()
