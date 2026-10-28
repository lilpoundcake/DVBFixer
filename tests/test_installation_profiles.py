"""Installation composition contracts; no scientific imports needed."""

import subprocess
import sys
from pathlib import Path

import pytest

from dvbfixer.installation_profiles import composition, environment_spec, load_catalog


def test_desktop_inherits_science_and_only_adds_gui():
    standard = composition("standard")
    desktop = composition("desktop")
    assert [p for p in desktop if p["group"] == "science"] == standard
    assert [p["spec"] for p in desktop if p not in standard] == ["nodejs >=22"]
    assert {"pyscf", "mdanalysis", "plotly", "acpype", "modeller"} <= {
        p["spec"] for p in standard
    }


def test_single_openbabel_provider_and_no_ml_metapackage():
    packages = composition("standard")
    assert [p["spec"] for p in packages if p["spec"].startswith("openbabel")] == [
        "openbabel-wheel >=3.1.1.23"
    ]
    assert next(p for p in packages if p["spec"] == "acpype")["source"] == "pip"
    assert "openff-toolkit-base" in {p["spec"] for p in packages}
    assert not {"openff-toolkit", "torch", "pytorch", "xtb", "nodejs >=22"} & {
        p["spec"] for p in packages
    }


def test_addons_are_explicit_and_external_runtimes_are_separate():
    assert "muscle >=5" not in environment_spec("standard")
    assert "muscle >=5" in environment_spec("standard", ("msa-alternatives",))
    assert load_catalog()["external"]["psi4"]["kind"] == "separate-environment"
    with pytest.raises(ValueError, match="Unknown"):
        composition("standard", ("protenix-cuda",))
    with pytest.raises(ValueError, match="Unknown"):
        composition("full")


def test_generated_composition_has_no_drift():
    root = Path(__file__).resolve().parents[1]
    subprocess.run([sys.executable, str(root / "scripts/gen_installation_profiles.py"), "--check"],
                   check=True, cwd=root)
