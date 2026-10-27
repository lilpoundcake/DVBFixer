"""Dependency-policy and native-library checks for the scientific environment."""

from __future__ import annotations

import subprocess
import sys
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_conda_blas_and_lapack_select_one_backend() -> None:
    dependencies = {
        line.split("#", 1)[0].strip().removeprefix("- ").strip()
        for line in (ROOT / "environment.yml").read_text().splitlines()
        if line.startswith("  - ")
    }
    assert {"libblas=*=*openblas", "libcblas=*=*openblas", "liblapack=*=*openblas"} <= dependencies
    assert "python >=3.11,<3.14" in dependencies
    assert "numpy <2.5" in dependencies


def test_numpy_scipy_native_linear_algebra_smoke() -> None:
    subprocess.run([sys.executable, str(ROOT / "scripts/check_scientific_linalg.py")], check=True, timeout=30)


def test_ter_naming_regressions_run_in_fork_safe_ci() -> None:
    workflow = (ROOT / ".github/workflows/ci.yml").read_text()
    fast_job = workflow.split("\n  lint-and-fast-tests:\n", 1)[1].split("\n  package-distributions:\n", 1)[0]
    pytest_step = fast_job.split("- name: pytest (fast lane, pure naming and registry contracts)", 1)[1]
    pytest_step = pytest_step.split("\n      - name:", 1)[0]
    assert "tests/test_atom_names_ter_integration.py" in pytest_step
    assert "tests/test_environment_policy.py" in pytest_step


def test_current_release_examples_match_authoritative_version() -> None:
    version = tomllib.loads((ROOT / "pyproject.toml").read_text())["project"]["version"]
    assert f"Current release: **{version}**." in (ROOT / "README.md").read_text()
    assert f"# dvbfixer {version}\n" in (ROOT / "docs/installation.md").read_text()
