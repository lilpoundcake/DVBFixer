"""Removed xTB options must fail rather than silently choosing another engine."""
import pytest

pytest.importorskip("openmm")


@pytest.mark.parametrize("arguments", [["--xtb-refine"], ["--xtb-cycles", "200"]])
def test_minimize_rejects_removed_options(arguments):
    from dvbfixer.minimize.cli import parse_args
    with pytest.raises(SystemExit) as exc:
        parse_args(["input.pdb", "--no-solvent", *arguments])
    assert exc.value.code == 2


def test_zbs_rejects_removed_refinement_choice():
    from dvbfixer.zbs import parse_args
    with pytest.raises(SystemExit) as exc:
        parse_args(["input.pdb", "--no-solvent", "--refine", "xtb"])
    assert exc.value.code == 2


def test_minimize_has_no_xtb_api():
    from dvbfixer import minimize
    from dvbfixer.minimize import refine
    assert not hasattr(minimize, "refine_with_xtb")
    assert not hasattr(refine, "refine_with_xtb")
    assert not hasattr(refine, "_run_xtb")
