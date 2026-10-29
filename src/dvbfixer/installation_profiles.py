"""Dependency composition shared by installers, generated specs and readiness checks.

Standard-library only: usable before installing the scientific stack.
"""

from __future__ import annotations

import json
import tomllib
from importlib.resources import files
from typing import Any


def load_catalog() -> dict[str, Any]:
    catalog = tomllib.loads(
        files("dvbfixer").joinpath("data/installation-profiles.toml").read_text(encoding="utf-8")
    )
    for package in catalog["packages"]:
        package.setdefault("channel", "PyPI" if package.get("source") == "pip"
                           else catalog["default_conda_channel"])
    return catalog


def composition(profile: str, addons: tuple[str, ...] = ()) -> list[dict[str, Any]]:
    catalog = load_catalog()
    groups: set[str] = set()
    seen: set[str] = set()
    current: str | None = profile
    while current is not None:
        if current in seen:
            raise ValueError("Cyclic installation profile inheritance")
        seen.add(current)
        if current not in catalog["profiles"]:
            raise ValueError(f"Unknown installation profile: {current}")
        entry = catalog["profiles"][current]
        groups.update(entry.get("groups", []))
        current = entry.get("inherits")
    for addon in addons:
        if addon not in catalog["addons"]:
            raise ValueError(f"Unknown in-environment add-on: {addon}")
        groups.update(catalog["addons"][addon]["groups"])
    return [package for package in catalog["packages"] if package["group"] in groups]


def environment_spec(profile: str, addons: tuple[str, ...] = ()) -> str:
    """Render portable YAML; JSON quoting avoids YAML/path ambiguity."""
    packages = composition(profile, addons)
    lines = ["# Generated from src/dvbfixer/data/installation-profiles.toml.",
             f"name: dvbfixer-{profile}", "channels:"]
    lines.extend(f"  - {channel}" for channel in load_catalog()["channels"])
    lines.append("  - nodefaults")
    lines.append("dependencies:")
    lines.extend(f"  - {json.dumps(p['channel'] + '::' + p['spec'])}"
                 for p in packages if p.get("source") != "pip")
    pip = [p["spec"] for p in packages if p.get("source") == "pip"]
    if pip:
        lines.append("  - pip:")
        lines.extend(f"      - {json.dumps(spec)}" for spec in pip)
    return "\n".join(lines) + "\n"


def dependency_matrix() -> str:
    lines = ["# Installation dependency matrix", "",
             "Generated from `src/dvbfixer/data/installation-profiles.toml`.", "",
             "Desktop inherits every Standard scientific dependency. `science` means both profiles;",
             "`gui` means Desktop only. Other groups are explicit optional add-ons.", "",
             "| Package / constraint | Source | Group | Consumers | Runtime / manual check |",
             "|---|---|---|---|---|"]
    for p in load_catalog()["packages"]:
        checks = [p["module"]] if "module" in p else []
        checks.extend(p.get("executables", []))
        checks.extend(p[k] for k in ("check", "manual") if k in p)
        lines.append(f"| `{p['spec']}` | {p['channel']} | "
                     f"{p['group']} | {p['consumer']} | {'; '.join(checks)} |")
    lines.extend(["", "Channels, in priority order: " + ", ".join(load_catalog()["channels"]) + ".",
                  "", "## Managed add-on platforms", "",
                  "These restrictions apply to provisioning, not to doctor checks for an operator-provided tool.",
                  "", "| Add-on | Managed platforms | Manual alternative |", "|---|---|---|"])
    catalog = load_catalog()
    for name, addon in catalog["addons"].items():
        platforms = ", ".join(addon.get("managed_platforms", catalog["platforms"]))
        route = addon.get("manual_route")
        alternative = f"[Installation guide]({route.removeprefix('docs/')})" if route else "—"
        lines.append(f"| `{name}` | {platforms} | {alternative} |")
    lines.extend(["", "Version/build and execution evidence: [profile acceptance](installation-profile-acceptance.md).",
                  "Licenses/notices remain governed by upstream distributions and `THIRD_PARTY_NOTICES/`.", ""])
    return "\n".join(lines)
