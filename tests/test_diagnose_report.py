"""Unit tests for the diagnose report layer (Severity, Finding, formatter)."""

from __future__ import annotations

from dvbfixer.diagnose.report import (
    Finding,
    Severity,
    format_report,
    sev_at_least,
)


def _f(sev: Severity, chain: str = "A", resid: str = "10", **kw) -> Finding:
    return Finding(
        severity=sev,
        category=kw.get("category", "coincident_atoms"),
        chain=chain,
        resid=resid,
        resname=kw.get("resname", "SER"),
        message=kw.get("message", "test finding"),
        atom=kw.get("atom"),
        fix_hint=kw.get("fix_hint"),
    )


class TestSeverityOrder:
    def test_at_least_error_threshold(self) -> None:
        assert sev_at_least(Severity.ERROR, Severity.ERROR)
        assert not sev_at_least(Severity.WARNING, Severity.ERROR)
        assert not sev_at_least(Severity.INFO, Severity.ERROR)

    def test_at_least_warning_threshold(self) -> None:
        assert sev_at_least(Severity.ERROR, Severity.WARNING)
        assert sev_at_least(Severity.WARNING, Severity.WARNING)
        assert not sev_at_least(Severity.INFO, Severity.WARNING)

    def test_at_least_info_threshold(self) -> None:
        assert sev_at_least(Severity.ERROR, Severity.INFO)
        assert sev_at_least(Severity.WARNING, Severity.INFO)
        assert sev_at_least(Severity.INFO, Severity.INFO)


class TestFindingSortKey:
    def test_severity_first(self) -> None:
        errors = _f(Severity.ERROR).sort_key()
        warnings = _f(Severity.WARNING).sort_key()
        infos = _f(Severity.INFO).sort_key()
        assert errors < warnings < infos

    def test_chain_then_numeric_resid(self) -> None:
        a = _f(Severity.ERROR, chain="A", resid="5").sort_key()
        b = _f(Severity.ERROR, chain="A", resid="10").sort_key()
        c = _f(Severity.ERROR, chain="B", resid="5").sort_key()
        # Numeric ordering: A:5 < A:10 < B:5
        assert a < b < c

    def test_icode_after_numeric(self) -> None:
        base = _f(Severity.ERROR, chain="A", resid="100").sort_key()
        icode = _f(Severity.ERROR, chain="A", resid="100A").sort_key()
        assert base < icode


class TestFormatReport:
    def test_clean_structure(self) -> None:
        text = format_report("clean.pdb", 100, 20, 2, [])
        assert "No findings" in text
        assert "clean.pdb" in text

    def test_findings_grouped_by_heading(self) -> None:
        findings = [
            _f(Severity.ERROR, category="coincident_atoms"),
            _f(Severity.ERROR, category="valence"),
            _f(Severity.ERROR, category="clash"),
        ]
        text = format_report("test.pdb", 100, 20, 2, findings)
        # Deterministic heading order.
        pos_structural = text.index("Structural integrity")
        pos_chemistry = text.index("Chemistry / bond geometry")
        pos_steric = text.index("Steric analysis")
        assert pos_structural < pos_chemistry < pos_steric

    def test_summary_counts(self) -> None:
        findings = [
            _f(Severity.ERROR),
            _f(Severity.ERROR),
            _f(Severity.WARNING),
        ]
        text = format_report("t.pdb", 100, 20, 2, findings)
        assert "ERROR:   2 findings" in text
        assert "WARNING: 1 findings" in text

    def test_suggested_next_step_appears_on_error(self) -> None:
        findings = [_f(Severity.ERROR, category="coincident_atoms")]
        text = format_report("t.pdb", 10, 1, 1, findings)
        assert "Suggested next step" in text
        assert "dvbfixer prepare" in text

    def test_no_suggestion_when_only_warnings(self) -> None:
        findings = [_f(Severity.WARNING, category="clash")]
        text = format_report("t.pdb", 10, 1, 1, findings)
        assert "Suggested next step" not in text
