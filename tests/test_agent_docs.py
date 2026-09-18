from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def _load_checker():
    path = ROOT / "scripts/check_agent_docs.py"
    spec = importlib.util.spec_from_file_location("check_agent_docs", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_agent_documentation_maps_are_valid() -> None:
    checker = _load_checker()
    assert checker.validate_agent_docs(ROOT) == []


def test_reference_path_parser_ignores_external_consumers() -> None:
    checker = _load_checker()
    assert checker._path_part("src/dvbfixer/cli.py::main") == "src/dvbfixer/cli.py"
    assert checker._path_part("external GROMACS tooling") is None


def test_invalid_status_is_reported() -> None:
    checker = _load_checker()
    errors: list[str] = []
    checker._check_status({"status": "finished"}, "tasks.example", errors)
    assert errors == ["tasks.example: invalid status 'finished'"]
