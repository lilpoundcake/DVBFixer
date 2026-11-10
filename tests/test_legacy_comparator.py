"""Declared production comparator preserves .dat, no-solvent and raw naming."""
from pathlib import Path

from scripts.research_legacy_comparator import commands


def test_comparator_keeps_production_physics_and_explicit_dev_solvent_convention():
    prepare,minimize=commands('/python',Path('/owned/case'))
    assert prepare[prepare.index('--backend')+1]=='legacy'
    assert '--no-solvent' in minimize
    assert prepare[prepare.index('--dat')+1]==minimize[minimize.index('--dat')+1]
    for command in (prepare,minimize):
        assert command[command.index('--atom-naming')+1]=='standard'
        assert not {'--platform','--no-propka','--no-protassign','--restraint-k','--weak-k','--max-iter'} & set(command)
