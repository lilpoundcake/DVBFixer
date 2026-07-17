"""Unit tests for the ``--acpype`` spec parsers in ``dvbfixer.top.acpype``.

These parsers are pure argument-transformation helpers with no OpenMM /
Modeller / ACPYPE dependency, so they run cleanly in the fast lane.
"""

from __future__ import annotations

import pytest

from dvbfixer.top.acpype import (
    _parse_protonation_overrides,
    _parse_ss_pairs,
)


class TestParseSsPairs:
    def test_single_pair(self) -> None:
        assert _parse_ss_pairs(["A:22:A:96"]) == {("A", 22), ("A", 96)}

    def test_multiple_pairs(self) -> None:
        result = _parse_ss_pairs(["A:22:A:96", "B:10:B:80"])
        assert result == {("A", 22), ("A", 96), ("B", 10), ("B", 80)}

    def test_inter_chain(self) -> None:
        assert _parse_ss_pairs(["A:22:B:96"]) == {("A", 22), ("B", 96)}

    def test_empty(self) -> None:
        assert _parse_ss_pairs([]) == set()

    def test_malformed_dropped(self) -> None:
        """Historical behaviour: 3-part or 5-part specs silently skipped."""
        assert _parse_ss_pairs(["A:22", "A:22:A:96", "A:B:C:D:E"]) == {("A", 22), ("A", 96)}


class TestParseProtonationOverrides:
    def test_none_and_empty(self) -> None:
        assert _parse_protonation_overrides(None, []) == {}
        assert _parse_protonation_overrides("all", []) == {}

    def test_single_protonate_spec(self) -> None:
        assert _parse_protonation_overrides("A:66:GLH", []) == {("A", 66): "GLH"}

    def test_charmm_names_normalised(self) -> None:
        """CHARMM ASPP/GLUP/HSP names get mapped to AMBER equivalents."""
        assert _parse_protonation_overrides("A:66:GLUP", []) == {("A", 66): "GLH"}
        assert _parse_protonation_overrides("A:66:HSD", []) == {("A", 66): "HID"}
        assert _parse_protonation_overrides("A:66:ASPP", []) == {("A", 66): "ASH"}

    def test_multiple_specs(self) -> None:
        result = _parse_protonation_overrides("A:66:GLH,B:20:HID", [])
        assert result == {("A", 66): "GLH", ("B", 20): "HID"}

    def test_unknown_variant_silently_dropped(self) -> None:
        """Anything not in AMBER_VARIANTS (after normalisation) is dropped."""
        assert _parse_protonation_overrides("A:66:BOGUS", []) == {}

    def test_his_specs_merge(self) -> None:
        result = _parse_protonation_overrides("A:66:GLH", ["B:20:HID", "C:30:HIP"])
        assert result == {("A", 66): "GLH", ("B", 20): "HID", ("C", 30): "HIP"}

    def test_his_specs_alone(self) -> None:
        result = _parse_protonation_overrides(None, ["A:66:HIE"])
        assert result == {("A", 66): "HIE"}

    def test_malformed_protonate_without_colon_exits(self) -> None:
        """Historical behaviour: invalid --protonate values print an error and sys.exit(1)."""
        with pytest.raises(SystemExit) as excinfo:
            _parse_protonation_overrides("bogus", [])
        assert excinfo.value.code == 1
