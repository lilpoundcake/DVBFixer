"""Installation composition contracts; no scientific imports needed."""

import subprocess
import sys
from pathlib import Path

import pytest

from dvbfixer.installation_profiles import (
    composition,
    dependency_matrix,
    environment_spec,
    load_catalog,
)


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
    alternatives = environment_spec("standard", ("msa-alternatives",))
    assert "bioconda::muscle >=5,<6" in alternatives
    assert "bioconda::clustalo" in alternatives
    assert load_catalog()["external"]["psi4"]["kind"] == "separate-environment"
    assert {"pydantic<2", "libint=2.9.0"} <= set(load_catalog()["external"]["psi4"]["dependencies"])
    with pytest.raises(ValueError, match="Unknown"):
        composition("standard", ("protenix-cuda",))
    with pytest.raises(ValueError, match="Unknown"):
        composition("full")


def test_generated_composition_has_no_drift():
    root = Path(__file__).resolve().parents[1]
    subprocess.run([sys.executable, str(root / "scripts/gen_installation_profiles.py"), "--check"],
                   check=True, cwd=root)


@pytest.mark.parametrize("profile", ["standard", "desktop"])
def test_ambertools_compat_is_explicit_and_keeps_all_science(profile):
    base = composition(profile)
    combined = composition(profile, ("ambertools-compat",))
    assert all(package in combined for package in base)
    specs = {p["spec"] for p in combined if p["group"] == "ambertools-compat"}
    assert specs == {"python <3.13", "numpy >=1.26,<2", "biopython >=1.83,<1.86", "pdb2pqr"}
    assert not any(p["group"] == "ambertools-compat" for p in base)
    assert "python >=3.11,<3.14" in environment_spec(profile, ("ambertools-compat",))
    assert any(p["spec"] == "pyscf" for p in combined)


def test_managed_probe_platforms_do_not_hide_operator_tool_readiness():
    addon = load_catalog()["addons"]["probe"]
    assert addon["managed_platforms"] == ["linux-64"]
    assert addon["manual_route"] == "docs/installation.md#native-apple-probe-alternative"
    assert any(p.get("check") == "probe-version" for p in composition("standard", ("probe",)))
    matrix = dependency_matrix()
    assert "| `probe` | linux-64 |" in matrix
    assert "(installation.md#native-apple-probe-alternative)" in matrix
