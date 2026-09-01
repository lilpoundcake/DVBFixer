from __future__ import annotations

import tomllib
from pathlib import Path

import pytest


def test_sugar_consumers_share_the_canonical_registry() -> None:
    from dvbfixer.acpype_export import _is_glycam_sugar
    from dvbfixer.cluster import _is_sugar
    from dvbfixer.prepare.cli import SUGAR_RESNAMES
    from dvbfixer.residue_registry import (
        CHARMM_SUGAR_RESNAMES,
        PDB_SUGAR_RESNAMES,
        is_pdb_or_glycam_sugar,
    )

    assert SUGAR_RESNAMES == set(PDB_SUGAR_RESNAMES | CHARMM_SUGAR_RESNAMES)
    for name in PDB_SUGAR_RESNAMES:
        assert _is_sugar(name)
        assert _is_glycam_sugar(name) == is_pdb_or_glycam_sugar(name)
    for name in CHARMM_SUGAR_RESNAMES:
        assert _is_sugar(name)


def test_wheel_configuration_installs_bundled_gromacs_force_fields() -> None:
    root = Path(__file__).parents[1]
    config = tomllib.loads((root / "pyproject.toml").read_text())
    data_files = config["tool"]["setuptools"]["data-files"]
    assert data_files["share/dvbfixer/FF/amber99sb-ildn-lipid21.ff"]
    assert data_files["share/dvbfixer/FF/charmm36_ljpme-jul2022.ff"]


def test_numpy_constraint_matches_python_311_typecheck_baseline() -> None:
    root = Path(__file__).parents[1]
    config = tomllib.loads((root / "pyproject.toml").read_text())
    dependencies = config["project"]["dependencies"]
    environment = (root / "environment.yml").read_text()

    assert "numpy<2.5" in dependencies
    assert "- numpy <2.5" in environment


def test_strip_solvent_can_retain_only_input_waters() -> None:
    openmm = pytest.importorskip("openmm")
    from openmm.app import PDBFile

    from dvbfixer.minimize.pipeline import _water_residue_key, strip_solvent

    pdb_text = """\
HETATM    1  O   HOH A 101       0.000   0.000   0.000  1.00  0.00           O
HETATM    2  O   HOH B 202       5.000   0.000   0.000  1.00  0.00           O
END
"""
    from io import StringIO

    pdb = PDBFile(StringIO(pdb_text))
    waters = list(pdb.topology.residues())
    topology, _ = strip_solvent(
        pdb.topology, pdb.positions, {_water_residue_key(waters[0])},
    )
    kept = list(topology.residues())
    assert len(kept) == 1
    assert _water_residue_key(kept[0]) == _water_residue_key(waters[0])
    assert openmm is not None


def test_zbs_forwards_keep_water_to_prepare_and_minimize(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    import dvbfixer.minimize
    import dvbfixer.prepare
    from dvbfixer.zbs import _run_pipeline, parse_args

    input_path = tmp_path / "input.pdb"
    input_path.write_text("END\n")
    calls: dict[str, list[str]] = {}

    def fake_step(name: str):
        def run(argv: list[str]) -> None:
            calls[name] = argv
            output = Path(argv[argv.index("-o") + 1])
            output.write_text("END\n")
        return run

    monkeypatch.setattr(dvbfixer.prepare, "main", fake_step("prepare"))
    monkeypatch.setattr(dvbfixer.minimize, "main", fake_step("minimize"))
    args = parse_args([
        str(input_path), "-o", str(tmp_path / "output.pdb"),
        "--skip-renumber", "--skip-model", "--keep-water",
        "--no-postflight", "--no-align-to-input",
    ])
    _run_pipeline(args, input_path)

    assert "--keep-water" in calls["prepare"]
    assert "--keep-water" in calls["minimize"]
