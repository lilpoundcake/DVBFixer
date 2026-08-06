"""Regression tests for `dvbfixer conect`.

CONECT inference depends on OpenBabel; if OpenBabel isn't importable the
tests are skipped rather than failing (some CI envs won't have it).
"""

from __future__ import annotations

from pathlib import Path

import pytest

openbabel = pytest.importorskip("openbabel", reason="conect inference needs OpenBabel")

from dvbfixer.conect import main  # noqa: E402


def _count_conect(pdb: Path) -> int:
    return sum(1 for ln in pdb.read_text().splitlines() if ln.startswith("CONECT"))


def test_conect_writes_records(tmp_workdir: Path, hinge_ch3_glycosylated_pdb: Path) -> None:
    """The tracked hinge fixture has cysteines — an SS bond is emitted."""
    out = tmp_workdir / "conect.pdb"
    main([str(hinge_ch3_glycosylated_pdb), "-o", str(out)])
    n = _count_conect(out)
    assert n > 0, "expected at least one CONECT record on a Cys-containing PDB"


def test_conect_is_idempotent(tmp_workdir: Path, hinge_ch3_glycosylated_pdb: Path) -> None:
    """Running conect twice produces the same output."""
    out1 = tmp_workdir / "pass1.pdb"
    out2 = tmp_workdir / "pass2.pdb"
    main([str(hinge_ch3_glycosylated_pdb), "-o", str(out1)])
    main([str(out1), "-o", str(out2)])
    # Compare CONECT sets to be robust to line order.
    c1 = sorted(ln.rstrip() for ln in out1.read_text().splitlines() if ln.startswith("CONECT"))
    c2 = sorted(ln.rstrip() for ln in out2.read_text().splitlines() if ln.startswith("CONECT"))
    assert c1 == c2


def test_conect_refuses_in_place_without_force(tmp_workdir: Path, hinge_ch3_glycosylated_pdb: Path) -> None:
    """--output == input should error unless --force is passed."""
    inp = tmp_workdir / "in.pdb"
    inp.write_bytes(hinge_ch3_glycosylated_pdb.read_bytes())
    with pytest.raises(SystemExit):
        main([str(inp), "-o", str(inp)])
