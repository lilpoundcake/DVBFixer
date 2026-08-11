from __future__ import annotations

import argparse

import pytest

from dvbfixer.cli_types import (
    nonnegative_float,
    nonnegative_int,
    positive_float,
    positive_int,
)


@pytest.mark.parametrize(
    ("validator", "bad"),
    [
        (positive_int, "0"),
        (positive_int, "-1"),
        (nonnegative_int, "-1"),
        (positive_float, "0"),
        (positive_float, "-0.1"),
        (positive_float, "nan"),
        (positive_float, "inf"),
        (nonnegative_float, "-0.1"),
    ],
)
def test_shared_numeric_validators_reject_invalid_values(validator, bad: str) -> None:
    with pytest.raises(argparse.ArgumentTypeError):
        validator(bad)


def test_model_counts_fail_during_argument_parsing() -> None:
    from dvbfixer.model.cli import parse_args

    for option in ("--num-models", "--num-loops", "--num-output"):
        with pytest.raises(SystemExit) as caught:
            parse_args(["input.pdb", option, "0"])
        assert caught.value.code == 2


@pytest.mark.parametrize(
    "argv",
    [
        ["input.pdb", "--ss", "A:not-a-number:B:2"],
        ["input.pdb", "--ss", "A:1:B"],
        ["input.pdb", "--his", "A:10:UNKNOWN"],
        ["input.pdb", "--his", "A:x:HIE"],
        ["input.pdb", "--protonate", "A:not-a-number"],
        ["input.pdb", "--protonate", "A:10:UNKNOWN"],
    ],
)
def test_top_structured_specs_fail_during_argument_parsing(argv: list[str]) -> None:
    from dvbfixer.top.cli import parse_args

    with pytest.raises(SystemExit) as caught:
        parse_args(argv)
    assert caught.value.code == 2


def test_top_structured_specs_preserve_valid_strings() -> None:
    from dvbfixer.top.cli import parse_args

    args = parse_args([
        "input.pdb", "--ss", "A:1:B:2", "--his", "A:3:HIE",
        "--protonate", "A:4:ASH,B:5",
    ])
    assert args.ss == ["A:1:B:2"]
    assert args.his == ["A:3:HIE"]
    assert args.protonate == "A:4:ASH,B:5"


def test_pull_atom_specs_fail_during_argument_parsing() -> None:
    pytest.importorskip("openmm")
    from dvbfixer.pull import parse_args

    with pytest.raises(SystemExit) as caught:
        parse_args(["input.pdb", "--bond", "A:not-a-number:SG", "B:2:SG"])
    assert caught.value.code == 2


def test_removed_zbs_skip_protonate_option_is_rejected() -> None:
    from dvbfixer.zbs import parse_args

    with pytest.raises(SystemExit) as caught:
        parse_args(["input.pdb", "--skip-protonate"])
    assert caught.value.code == 2


@pytest.mark.parametrize(
    "option",
    ["--psi4-method", "--psi4-nthreads", "--psi4-memory"],
)
def test_removed_parametrize_psi4_aliases_are_rejected(option: str) -> None:
    from dvbfixer.parametrize import parse_args

    with pytest.raises(SystemExit) as caught:
        parse_args(["input.pdb", option, "value"])
    assert caught.value.code == 2
