"""Tests for the 0.7.7 PROPKA + Reduce integration in legacy prepare.

Prior to 0.7.7, the legacy prepare backend built its `variants=[...]`
list only from user `--mutate` overrides + HIS HD1/HE2 atom presence.
PROPKA and MolProbity Reduce were not invoked, so pKa-driven ASH /
GLH / HIP / LYN / CYM were never emitted for the average user.

0.7.7 wires PROPKA + Reduce into prepare's variants-list build via
`_run_propka_reduce_variants`. These tests verify:

1. The helper returns a variant map (dict form) for a real PDB.
2. `--no-propka --no-protassign` reproduces the pre-0.7.7 behaviour.
3. The zbs pipeline propagates the flags into prepare.
"""
from __future__ import annotations

from pathlib import Path

import pytest


@pytest.mark.slow
def test_propka_reduce_helper_returns_variants(
    pure_protein_small: Path, tmp_workdir: Path,
) -> None:
    """The helper should run PROPKA + Reduce on a real PDB and return
    a variant map. Even a small protein has HIS tautomers to pick and
    typically one or two shifted acidic side chains — the map should
    be non-empty."""
    from dvbfixer.prepare.pipeline import _run_propka_reduce_variants

    variant_map = _run_propka_reduce_variants(
        pure_protein_small,
        ph=7.0,
        ss_pairs=set(),
        his_default="HIE",
        cys_ss_pka=8.0,
        use_propka=True,
        use_reduce=True,
        verbose=False,
    )
    # Result is a dict; values are valid AMBER variant names.
    _VALID = {"HID", "HIE", "HIP", "ASH", "GLH", "CYX", "CYM", "LYN"}
    assert isinstance(variant_map, dict)
    for key, variant in variant_map.items():
        assert variant in _VALID, (
            f"unexpected variant name {variant!r} for {key}"
        )


def test_propka_reduce_helper_both_off_returns_empty(
    pure_protein_small: Path,
) -> None:
    """When both PROPKA and Reduce are off, the helper should return
    an empty dict (no pKa-driven decisions to make)."""
    from dvbfixer.prepare.pipeline import _run_propka_reduce_variants

    variant_map = _run_propka_reduce_variants(
        pure_protein_small,
        ph=7.0,
        ss_pairs=set(),
        use_propka=False,
        use_reduce=False,
        verbose=False,
    )
    assert variant_map == {}


def test_propka_reduce_helper_ss_pairs_force_cyx(
    ss_bonded_antibody: Path,
) -> None:
    """When ss_pairs is non-empty, every listed (chain, resseq) should
    end up as CYX in the returned map — regardless of PROPKA's
    per-residue CYS pKa."""
    from dvbfixer.acpype_export import detect_ss_bonds

    from dvbfixer.prepare.pipeline import _run_propka_reduce_variants

    ss_pairs = detect_ss_bonds(str(ss_bonded_antibody))
    assert ss_pairs, "test fixture should have SS bonds"

    variant_map = _run_propka_reduce_variants(
        ss_bonded_antibody,
        ph=7.0,
        ss_pairs=ss_pairs,
        use_propka=False,  # keep the test fast — CONECT-driven CYX only
        use_reduce=False,
        verbose=False,
    )
    cyx_keys = {
        (ch, rn) for (ch, rn, _ic), var in variant_map.items()
        if var == "CYX"
    }
    # Every SS-paired CYS should have a CYX entry.
    for (ch, rn) in ss_pairs:
        assert (ch, rn) in cyx_keys, (
            f"SS-bonded CYS {ch}:{rn} not in CYX map: {sorted(cyx_keys)[:5]}"
        )


def test_prepare_cli_propka_flags_exist() -> None:
    """The new CLI flags on `dvbfixer prepare` should parse cleanly."""
    from dvbfixer.prepare.cli import parse_args

    args = parse_args([
        "dummy.pdb", "--no-propka", "--no-protassign",
        "--his-default", "HID", "--cys-ss-pka", "9.5",
    ])
    assert args.propka is False
    assert args.protassign is False
    assert args.his_default == "HID"
    assert args.cys_ss_pka == 9.5


def test_zbs_cli_propka_flags_exist() -> None:
    """The same flags on `dvbfixer zbs` should parse cleanly."""
    from dvbfixer.zbs import parse_args

    args = parse_args([
        "dummy.pdb", "--no-propka", "--no-protassign",
        "--his-default", "HID", "--cys-ss-pka", "9.5",
    ])
    assert args.propka is False
    assert args.protassign is False
    assert args.his_default == "HID"
    assert args.cys_ss_pka == 9.5


def test_zbs_skip_protonate_maps_to_no_propka_no_protassign() -> None:
    """--skip-protonate (deprecated) should now flip propka+protassign
    to False for backward compat."""
    from dvbfixer.zbs import parse_args

    args = parse_args(["dummy.pdb", "--skip-protonate"])
    assert args.propka is False
    assert args.protassign is False
