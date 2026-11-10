"""Stage observer uses an isolated child and cannot pick a scientific platform."""
from pathlib import Path

from scripts.research_legacy_stereo_trace import worker_argv


def test_child_runner_keeps_exact_case_and_plan_without_optimizer_overrides():
    argv=worker_argv('/python',Path('/plan.json'),'case-A',Path('/owned/case'))
    assert argv[0]=='/python' and argv[2]=='worker'
    assert argv[argv.index('--case')+1]=='case-A'
    assert not {'--platform','--weak-k','--restraint-k','--max-iter'} & set(argv)
