"""Native capture failures must not silently become absence of diagnostics."""
import pytest
from scripts.research_chi_termination_flush import flush_checked


def test_failed_native_flush_is_a_blocker():
    with pytest.raises(OSError,match='fflush failed'):
        flush_checked(lambda _:1)


def test_successful_flush_does_not_change_an_optimizer_setting():
    calls=[]
    def flush(value):
        calls.append(value)
        return 0
    flush_checked(flush)
    assert calls==[None]
