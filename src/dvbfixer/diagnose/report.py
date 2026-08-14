"""Findings model + plain-text formatter for ``dvbfixer diagnose``."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import StrEnum
from typing import Any


class Severity(StrEnum):
    """Ordered severity levels. ERROR > WARNING > INFO."""

    ERROR = "ERROR"      # would break downstream tools if left as-is
    WARNING = "WARNING"  # tolerated by pipeline but suspect
    INFO = "INFO"        # noted for user review


# Precedence order for the report (higher-severity findings first).
_SEV_ORDER = {Severity.ERROR: 0, Severity.WARNING: 1, Severity.INFO: 2}


@dataclass
class Finding:
    """A single per-residue (or per-atom) diagnostic finding."""

    severity: Severity
    category: str          # 'coincident_atoms', 'clash', 'valence', ...
    chain: str
    resid: str             # includes iCode if any, e.g. '100A'
    resname: str
    message: str
    atom: str | None = None       # None for residue-level findings
    fix_hint: str | None = None    # optional suggested `dvbfixer X` command
    extra: dict[str, Any] = field(default_factory=dict)  # engine-specific detail

    def sort_key(self) -> tuple:
        """Sort by severity → chain → numeric resid → iCode → atom."""
        try:
            rs_num = int("".join(c for c in self.resid if c.isdigit() or c == "-") or "0")
        except ValueError:
            rs_num = 0
        rs_icode = "".join(c for c in self.resid if c.isalpha())
        return (
            _SEV_ORDER[self.severity],
            self.chain,
            rs_num,
            rs_icode,
            self.atom or "",
        )


# ---------------------------------------------------------------------------
# Category → human-readable heading
# ---------------------------------------------------------------------------
_CATEGORY_HEADINGS: dict[str, str] = {
    # Meta
    "multi_model": "Input",
    # Structural
    "missing_residues": "Structural integrity",
    "missing_atoms": "Structural integrity",
    "missing_terminals": "Structural integrity",
    "coincident_atoms": "Structural integrity",
    "misplaced_hydrogen": "Structural integrity",
    "altloc_conflict": "Structural integrity",
    "chain_break": "Structural integrity",
    "duplicate_chain_coordinates": "Structural integrity",
    "insertion_codes": "Structural integrity",
    "seqres_gap": "Structural integrity",
    # Chemistry
    "valence": "Chemistry / bond geometry",
    "bond_length": "Chemistry / bond geometry",
    "bond_angle": "Chemistry / bond geometry",
    "cis_peptide": "Chemistry / bond geometry",
    "non_planar_amide": "Chemistry / bond geometry",
    "chirality": "Chemistry / bond geometry",
    "disulfide_geometry": "Chemistry / bond geometry",
    # Steric
    "clash": "Steric analysis",
}


def _heading_for(category: str) -> str:
    return _CATEGORY_HEADINGS.get(category, "Other findings")


# Which dvbfixer subcommand fixes each category, used to synthesize the
# "Suggested next step" hint at the bottom of the report.
_FIX_TOOL: dict[str, str] = {
    "multi_model": "dvbfixer split",
    "missing_residues": "dvbfixer model",
    "missing_atoms": "dvbfixer prepare",
    "missing_terminals": "dvbfixer prepare",
    "coincident_atoms": "dvbfixer prepare",
    "misplaced_hydrogen": "dvbfixer prepare",
    "altloc_conflict": "manual: pick one altLoc",
    "chain_break": "dvbfixer model  (if it's a real gap)",
    "duplicate_chain_coordinates": "manual: restore MODEL/ENDMDL records or remove the duplicate chain",
    "insertion_codes": "dvbfixer renumber  (if renumbering is desired)",
    "seqres_gap": "dvbfixer model",
    "valence": "manual: inspect input; may need CONECT cleanup",
    "bond_length": "dvbfixer minimize",
    "bond_angle": "dvbfixer minimize",
    "cis_peptide": "manual: verify against reference structure",
    "non_planar_amide": "dvbfixer minimize",
    "chirality": "manual: rebuild the residue",
    "disulfide_geometry": "dvbfixer minimize",
    "clash": "dvbfixer minimize",
}


def sev_at_least(sev: Severity, threshold: Severity) -> bool:
    """True if ``sev`` should be included under a ``threshold`` filter."""
    return _SEV_ORDER[sev] <= _SEV_ORDER[threshold]


def findings_to_dict_list(findings: list[Finding]) -> list[dict[str, Any]]:
    """Serialise findings for the ``--format json`` output.

    Sorted by severity → chain → resid → atom so the JSON output
    mirrors the text report's ordering.
    """
    ordered = sorted(findings, key=Finding.sort_key)
    out: list[dict[str, Any]] = []
    for f in ordered:
        item: dict[str, Any] = {
            "severity": f.severity.value,
            "category": f.category,
            "chain": f.chain,
            "resid": f.resid,
            "resname": f.resname,
            "atom": f.atom,
            "message": f.message,
            "fix_hint": f.fix_hint,
        }
        if f.extra:
            item["details"] = f.extra
        out.append(item)
    return out


def format_report(
    input_path: str,
    n_atoms: int,
    n_residues: int,
    n_chains: int,
    findings: list[Finding],
    *,
    chirality_checked: bool = True,
    chirality_repairs: list[dict[str, str]] | None = None,
) -> str:
    """Render the plain-text report."""
    lines: list[str] = []
    bar = "=" * 80
    dash = "-" * 40
    lines.append(bar)
    lines.append(f"dvbfixer diagnose — {input_path}")
    lines.append(bar)
    lines.append(f"Loaded: {n_atoms} atoms, {n_residues} residues, {n_chains} chains")
    d_findings = [f for f in findings if f.category == "chirality"]
    if chirality_checked:
        status = "YES" if d_findings else "NO"
    else:
        status = "NOT CHECKED"
    lines.append(f"D-isomer error: {status}")
    repairs = chirality_repairs or []
    if repairs:
        labels = ", ".join(
            f"{item['chain']}/{item['resname']}{item['resid']} ({item['stage']})"
            for item in repairs
        )
        lines.append(f"Forced D→L repair history: YES — {labels}")
        lines.append("WARNING: inspect hydrogen angles and local geometry around repaired residues.")
    else:
        lines.append("Forced D→L repair history: NO")
    lines.append("")

    if not findings:
        lines.append("No findings. Structure looks clean.")
        lines.append("")
        lines.append(bar)
        return "\n".join(lines)

    # Group by heading, then sort within each group.
    by_heading: dict[str, list[Finding]] = {}
    for f in findings:
        by_heading.setdefault(_heading_for(f.category), []).append(f)

    # Preserve a deterministic heading order.
    for heading in ["Input", "Structural integrity", "Chemistry / bond geometry",
                    "Steric analysis", "Other findings"]:
        group = by_heading.get(heading)
        if not group:
            continue
        lines.append(heading)
        lines.append(dash)
        for f in sorted(group, key=Finding.sort_key):
            loc = f"{f.chain}/{f.resname}{f.resid}"
            if f.atom:
                loc = f"{loc}:{f.atom}"
            lines.append(f"  {f.severity.value:<7s} {loc:<22s} {f.message}")
            if f.fix_hint:
                lines.append(f"                                 fix: {f.fix_hint}")
        lines.append("")

    # Summary
    counts = {sev: 0 for sev in Severity}
    for f in findings:
        counts[f.severity] += 1
    residues_with_error = {(f.chain, f.resid) for f in findings if f.severity == Severity.ERROR}
    residues_with_warning = {(f.chain, f.resid) for f in findings if f.severity == Severity.WARNING}

    lines.append("Summary")
    lines.append(dash)
    lines.append(f"  ERROR:   {counts[Severity.ERROR]} findings across "
                 f"{len(residues_with_error)} residues")
    lines.append(f"  WARNING: {counts[Severity.WARNING]} findings across "
                 f"{len(residues_with_warning)} residues")
    lines.append(f"  INFO:    {counts[Severity.INFO]} findings")

    # Suggested next step — pick the tool that fixes the most-common
    # ERROR category (falls back to prepare).
    if counts[Severity.ERROR] > 0:
        cat_counts: dict[str, int] = {}
        for f in findings:
            if f.severity == Severity.ERROR:
                cat_counts[f.category] = cat_counts.get(f.category, 0) + 1
        top_cat = max(cat_counts, key=cat_counts.get)
        tool = _FIX_TOOL.get(top_cat, "dvbfixer prepare")
        lines.append("")
        lines.append(f"Suggested next step: `{tool} {input_path}`")

    lines.append("")
    lines.append(bar)
    return "\n".join(lines)
