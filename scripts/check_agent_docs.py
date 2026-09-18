#!/usr/bin/env python
"""Validate the lightweight, machine-readable agent knowledge map."""

from __future__ import annotations

import re
import tomllib
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
AGENT_DOCS = Path("docs/agent")
MAP_FILES = ("tasks.toml", "contracts.toml", "invariants.toml")
ALLOWED_STATUSES = {
    "implemented",
    "partial",
    "missing",
    "known-gap",
    "proposed",
    "research",
    "deprecated",
}
REQUIRED_CONTEXT_HEADINGS = {
    "Purpose",
    "Scope",
    "Capabilities",
    "Entry Points",
    "Contracts",
    "Invariants",
    "Callers",
    "Adapters",
    "Side Effects",
    "Known Divergences",
    "Proposed Work",
    "Focused Verification",
}
MARKDOWN_LINK_RE = re.compile(r"\[[^]]*]\(([^)]+)\)")


def _load_toml(root: Path, relative: Path, errors: list[str]) -> dict[str, Any]:
    path = root / relative
    try:
        return tomllib.loads(path.read_text())
    except FileNotFoundError:
        errors.append(f"missing map: {relative}")
    except tomllib.TOMLDecodeError as exc:
        errors.append(f"invalid TOML in {relative}: {exc}")
    return {}


def _path_part(reference: str) -> str | None:
    value = reference.split("::", 1)[0].strip()
    if not value or value.startswith("external "):
        return None
    if "/" not in value and not value.endswith((".py", ".ts", ".md", ".toml")):
        return None
    return value


def _check_references(
    root: Path, references: list[str], location: str, errors: list[str]
) -> None:
    for reference in references:
        path = _path_part(reference)
        if path is not None and not (root / path).exists():
            errors.append(f"{location}: missing referenced path {path!r}")


def _as_strings(value: Any, location: str, errors: list[str]) -> list[str]:
    if value is None:
        return []
    if not isinstance(value, list) or not all(isinstance(item, str) for item in value):
        errors.append(f"{location}: expected a list of strings")
        return []
    return value


def _check_status(record: dict[str, Any], location: str, errors: list[str]) -> None:
    status = record.get("status")
    if status not in ALLOWED_STATUSES:
        errors.append(f"{location}: invalid status {status!r}")


def _check_markdown(
    root: Path, relative: Path, errors: list[str], *, context_page: bool = True
) -> None:
    path = root / relative
    if not path.exists():
        errors.append(f"missing context page: {relative}")
        return
    text = path.read_text()
    if context_page:
        status_pattern = r"^Status: (?:" + "|".join(sorted(ALLOWED_STATUSES)) + r")\s*$"
        if not re.search(status_pattern, text, re.MULTILINE):
            errors.append(f"{relative}: missing or invalid Status line")
        if not re.search(r"^Verified on: \d{4}-\d{2}-\d{2}\s*$", text, re.MULTILINE):
            errors.append(f"{relative}: missing or invalid Verified on line")
        if not re.search(
            r"^Verified at commit: `[0-9a-f]{40}`\s*$", text, re.MULTILINE
        ):
            errors.append(f"{relative}: missing or invalid Verified at commit line")
        headings = set(re.findall(r"^## (.+)$", text, re.MULTILINE))
        missing = sorted(REQUIRED_CONTEXT_HEADINGS - headings)
        if missing:
            errors.append(f"{relative}: missing headings: {', '.join(missing)}")
    for target in MARKDOWN_LINK_RE.findall(text):
        if target.startswith(("http://", "https://", "#", "mailto:")):
            continue
        target_path = target.split("#", 1)[0]
        if target_path and not (path.parent / target_path).resolve().exists():
            errors.append(f"{relative}: broken link {target!r}")


def validate_agent_docs(root: Path = ROOT) -> list[str]:
    """Return validation errors for the agent documentation below ``root``."""
    errors: list[str] = []
    base = AGENT_DOCS
    loaded = {
        name: _load_toml(root, base / name, errors)
        for name in MAP_FILES
    }
    if errors:
        return errors

    for name, data in loaded.items():
        if data.get("schema_version") != 1:
            errors.append(f"{base / name}: schema_version must be 1")
        if not isinstance(data.get("verified_on"), str):
            errors.append(f"{base / name}: verified_on must be a date string")
        commit = data.get("verified_at_commit")
        if not isinstance(commit, str) or not re.fullmatch(r"[0-9a-f]{40}", commit):
            errors.append(f"{base / name}: verified_at_commit must be a full git SHA")

    tasks_data = loaded["tasks.toml"]
    contracts_data = loaded["contracts.toml"]
    invariants_data = loaded["invariants.toml"]
    tasks = tasks_data.get("tasks", {})
    groups = tasks_data.get("change_groups", {})
    contracts = contracts_data.get("contracts", {})
    invariants = invariants_data.get("invariants", {})

    for table_name, table in (
        ("tasks", tasks),
        ("change_groups", groups),
        ("contracts", contracts),
        ("invariants", invariants),
    ):
        if not isinstance(table, dict) or not table:
            errors.append(f"{table_name}: expected a non-empty table")

    context_paths: set[str] = set()
    for task_id, record in tasks.items():
        location = f"tasks.{task_id}"
        _check_status(record, location, errors)
        context = record.get("context")
        if not isinstance(context, str):
            errors.append(f"{location}: context must be a path")
        else:
            context_paths.add(context)
            _check_references(root, [context], location, errors)

        starts = _as_strings(record.get("start"), f"{location}.start", errors)
        callers = _as_strings(record.get("callers"), f"{location}.callers", errors)
        tests = _as_strings(record.get("test_files"), f"{location}.test_files", errors)
        verification = _as_strings(
            record.get("verification"), f"{location}.verification", errors
        )
        _check_references(root, starts + callers + tests, location, errors)

        for contract_id in _as_strings(
            record.get("contracts"), f"{location}.contracts", errors
        ):
            if contract_id not in contracts:
                errors.append(f"{location}: unknown contract {contract_id!r}")
        for invariant_id in _as_strings(
            record.get("invariants"), f"{location}.invariants", errors
        ):
            if invariant_id not in invariants:
                errors.append(f"{location}: unknown invariant {invariant_id!r}")
        group = record.get("change_group")
        if group not in groups:
            errors.append(f"{location}: unknown change group {group!r}")
        if record.get("status") == "implemented" and (not starts or not tests):
            errors.append(f"{location}: implemented tasks require start and test_files")
        if record.get("status") in {"implemented", "partial"} and not verification:
            errors.append(f"{location}: active tasks require verification commands")

    for group_id, record in groups.items():
        location = f"change_groups.{group_id}"
        if record.get("parallel_policy") not in {"exclusive-owner", "parallel-safe"}:
            errors.append(f"{location}: invalid parallel_policy")
        files = _as_strings(record.get("files"), f"{location}.files", errors)
        _check_references(root, files, location, errors)

    for contract_id, record in contracts.items():
        location = f"contracts.{contract_id}"
        _check_status(record, location, errors)
        owner = record.get("owner")
        tests = _as_strings(record.get("test_files"), f"{location}.test_files", errors)
        refs = list(tests)
        if isinstance(owner, str):
            refs.append(owner)
        else:
            errors.append(f"{location}: owner must be a reference")
        refs += _as_strings(record.get("producers"), f"{location}.producers", errors)
        refs += _as_strings(record.get("consumers"), f"{location}.consumers", errors)
        _check_references(root, refs, location, errors)
        for field in (
            "description",
            "identity_keys",
            "source_mutation",
            "commit_behavior",
            "failure_behavior",
            "known_divergences",
        ):
            if not record.get(field):
                errors.append(f"{location}: missing required field {field}")
        if record.get("status") == "implemented" and not tests:
            errors.append(f"{location}: implemented contracts require test_files")

    for invariant_id, record in invariants.items():
        location = f"invariants.{invariant_id}"
        _check_status(record, location, errors)
        owner = record.get("owner")
        tests = _as_strings(record.get("test_files"), f"{location}.test_files", errors)
        refs = list(tests)
        if isinstance(owner, str):
            refs.append(owner)
        else:
            errors.append(f"{location}: owner must be a reference")
        refs += _as_strings(record.get("enforcement"), f"{location}.enforcement", errors)
        _check_references(root, refs, location, errors)
        for field in ("statement", "failure_mode", "known_exceptions"):
            if field not in record:
                errors.append(f"{location}: missing required field {field}")
        if record.get("status") == "implemented" and not tests:
            errors.append(f"{location}: implemented invariants require test_files")

    _check_markdown(root, base / "README.md", errors, context_page=False)
    for context in sorted(context_paths):
        if context != str(base / "README.md"):
            _check_markdown(root, Path(context), errors)
    indexed_contexts = {
        Path(context)
        for context in context_paths
        if context != str(base / "README.md")
    }
    actual_contexts = {
        path.relative_to(root)
        for path in (root / base / "contexts").glob("*.md")
    }
    for context in sorted(actual_contexts - indexed_contexts):
        errors.append(f"unindexed context page: {context}")
    return errors


def main() -> int:
    errors = validate_agent_docs()
    if errors:
        for error in errors:
            print(f"ERROR: {error}")
        return 1
    print("agent documentation maps are valid")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
