"""Unit tests for the shared variant helpers in ``dvbfixer.ffutils.variants``."""

from __future__ import annotations

from pathlib import Path

import pytest

from dvbfixer.ffutils.variants import (
    ALL_VARIANTS,
    AMBER_VARIANT_TO_PARENT,
    CHARMM_VARIANT_TO_PARENT,
    VARIANT_TO_PARENT,
    restore_variants_post_addhydrogens,
    scan_variant_names,
    text_rename_variants_to_parent,
)


class _FakeAtom:
    def __init__(self, name: str) -> None:
        self.name = name


class _FakeChain:
    def __init__(self, chain_id: str) -> None:
        self.id = chain_id


class _FakeResidue:
    def __init__(self, chain_id: str, res_id: str, name: str, icode: str = "") -> None:
        self.chain = _FakeChain(chain_id)
        self.id = res_id
        self.name = name
        self.insertionCode = icode

    def atoms(self):
        return []


class _FakeTopology:
    def __init__(self, residues: list[_FakeResidue]) -> None:
        self._residues = residues

    def residues(self):
        return list(self._residues)


def test_variant_tables_cover_expected_names() -> None:
    for name in ("HID", "HIE", "HIP", "ASH", "GLH", "CYX", "CYM", "LYN"):
        assert AMBER_VARIANT_TO_PARENT[name]
    for name in ("HSD", "HSE", "HSP", "ASPP", "GLUP", "LSN"):
        assert CHARMM_VARIANT_TO_PARENT[name]
    # Combined lookup is the union.
    for name in ALL_VARIANTS:
        assert VARIANT_TO_PARENT[name] in {"HIS", "ASP", "GLU", "CYS", "LYS"}


def test_variant_tables_do_not_overlap() -> None:
    """No name appears in both AMBER and CHARMM tables (would create ambiguity)."""
    assert not set(AMBER_VARIANT_TO_PARENT) & set(CHARMM_VARIANT_TO_PARENT)


def test_scan_variant_names_finds_cyx(default_pdb: Path) -> None:
    """default.pdb has CYX residues; scan should locate them."""
    saved = scan_variant_names(default_pdb)
    assert saved, "expected at least one variant"
    assert any(name == "CYX" for name in saved.values()), (
        f"expected CYX in scan, got {set(saved.values())}"
    )


def test_scan_variant_names_ignores_standard_residues(small_pdb: Path) -> None:
    """ASN.pdb has one ASN — no variants."""
    saved = scan_variant_names(small_pdb)
    assert saved == {}


def test_text_rename_variants_to_parent_rewrites_cyx(
    tmp_workdir: Path, default_pdb: Path
) -> None:
    """CYX residues in default.pdb → CYS in the rewritten file."""
    out_path, saved = text_rename_variants_to_parent(default_pdb)
    assert saved, "expected renames"
    assert out_path != str(default_pdb), "should have written a temp file"

    # Rewritten file has no CYX left.
    rewritten_names = {
        ln[17:20].strip()
        for ln in Path(out_path).read_text().splitlines()
        if ln.startswith(("ATOM  ", "HETATM"))
    }
    assert "CYX" not in rewritten_names
    assert "CYS" in rewritten_names

    # Clean up temp file.
    Path(out_path).unlink(missing_ok=True)


def test_text_rename_variants_to_parent_no_op_on_standard(
    small_pdb: Path,
) -> None:
    """ASN.pdb has no variants — rename returns the input path unchanged."""
    out_path, saved = text_rename_variants_to_parent(small_pdb)
    assert out_path == str(small_pdb)
    assert saved == {}


def test_text_rename_variants_preserves_columns(
    tmp_workdir: Path, default_pdb: Path
) -> None:
    """Rewritten line must keep atom coord columns byte-identical."""
    orig_lines = {
        ln[6:11].strip(): ln[30:54]
        for ln in default_pdb.read_text().splitlines()
        if ln.startswith(("ATOM  ", "HETATM"))
    }
    out_path, _saved = text_rename_variants_to_parent(default_pdb)
    if out_path == str(default_pdb):
        pytest.skip("no rewrite happened; nothing to compare")

    for ln in Path(out_path).read_text().splitlines():
        if ln.startswith(("ATOM  ", "HETATM")):
            serial = ln[6:11].strip()
            assert ln[30:54] == orig_lines[serial], f"coord drift at serial {serial}"
    Path(out_path).unlink(missing_ok=True)


def test_restore_variants_topology_shape_only() -> None:
    """2-tuple `(chain_id, res_id)` keys — the shape used when the topology
    was pre-renamed via `rename_variants_to_parent_in_topology`.
    """
    top = _FakeTopology([
        _FakeResidue("A", "5", "HIS"),
        _FakeResidue("A", "7", "LYS"),
    ])
    saved = {("A", "5"): "HIE", ("A", "7"): "LYN"}
    restore_variants_post_addhydrogens(top, saved)
    names = [r.name for r in top.residues()]
    assert names == ["HIE", "LYN"]


def test_restore_variants_text_shape_only() -> None:
    """3-tuple `(chain, resseq, icode)` keys — the shape used when the
    file was pre-renamed via `text_rename_variants_to_parent`.
    """
    top = _FakeTopology([_FakeResidue("A", "5", "HIS")])
    saved = {("A", "5", ""): "HIE"}
    restore_variants_post_addhydrogens(top, saved)
    assert top.residues()[0].name == "HIE"


def test_restore_variants_mixed_shapes_no_crash() -> None:
    """Regression for the tuple-arity bug that crashed ``zbs`` step 5.

    protonate.py builds ``_saved`` by merging text-shape (3-tuple) keys
    into the topology-shape (2-tuple) map produced by
    ``rename_variants_to_parent_in_topology``. The old fallback
    destructured every key as a 3-tuple unconditionally and crashed with
    ``ValueError: not enough values to unpack (expected 3, got 2)`` the
    moment ``saved`` contained a 2-tuple entry alongside a 3-tuple one.

    Ensure the fallback tolerates a dict containing BOTH shapes without
    raising — and, per the icode-collapse audit fix, that a 3-tuple entry
    recorded under a DIFFERENT insertion code than the residue actually
    has is correctly NOT applied (residue A:7 here has no icode; the
    saved entry is for icode "B" — a different, hypothetical sibling
    residue — so A:7 must stay unrestored).
    """
    top = _FakeTopology([
        _FakeResidue("A", "5", "HIS"),
        _FakeResidue("A", "7", "LYS"),
    ])
    # Mixed: A:5 keyed by 2-tuple (topology shape, direct match).
    # A:7's *icode* ("") doesn't appear in `saved` at all — only a
    # same-resSeq, different-icode ("B") entry does, which must NOT match.
    saved = {
        ("A", "5"): "HIE",
        ("A", "7", "B"): "LYN",
    }
    # Should not raise, and must not cross-contaminate icode "" with icode "B".
    restore_variants_post_addhydrogens(top, saved)
    names = [r.name for r in top.residues()]
    assert names == ["HIE", "LYS"]


def test_restore_variants_icode_exact_match() -> None:
    """Two residues sharing a resSeq via insertion code (Kabat CDR-loop
    style, e.g. ``H:82``/``H:82A``) must each get their OWN variant, not
    whichever one a broad chain+resSeq scan happens to find first."""
    top = _FakeTopology([
        _FakeResidue("H", "82", "HIS", icode=""),
        _FakeResidue("H", "82", "LYS", icode="A"),
    ])
    saved = {
        ("H", "82", ""): "HIE",
        ("H", "82", "A"): "LYN",
    }
    restore_variants_post_addhydrogens(top, saved)
    names = [r.name for r in top.residues()]
    assert names == ["HIE", "LYN"]


def test_restore_variants_no_op_when_no_match() -> None:
    """Residues in the topology that aren't in `saved` keep their name."""
    top = _FakeTopology([
        _FakeResidue("A", "5", "HIS"),
        _FakeResidue("A", "999", "ARG"),
    ])
    saved = {("A", "5"): "HIE"}
    restore_variants_post_addhydrogens(top, saved)
    names = [r.name for r in top.residues()]
    assert names == ["HIE", "ARG"]
