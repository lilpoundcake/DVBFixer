"""Borrowed stream lifetime must not break late batch-capture cleanup."""

import io

import pytest

from dvbfixer.batch import _Mirror


def test_mirror_can_close_after_borrowed_capture_target():
    target = io.StringIO()
    mirror = _Mirror(target)
    mirror.write("WARNING: retained diagnostic\n")
    target.close()
    mirror.flush()
    mirror.close()
    assert mirror.closed
    assert mirror.captured() == "WARNING: retained diagnostic\n"


def test_mirror_closing_does_not_close_its_borrowed_target():
    target = io.StringIO()
    mirror = _Mirror(target)
    mirror.write("forwarded\n")
    mirror.close()
    assert not target.closed
    target.write("still usable\n")
    assert target.getvalue() == "forwarded\nstill usable\n"


@pytest.mark.parametrize("exception", [ValueError, OSError])
def test_mirror_does_not_hide_flush_failures_from_open_targets(exception):
    class FailingTarget(io.StringIO):
        fail = True

        def flush(self):
            if self.fail:
                raise exception("real I/O failure")
            return super().flush()

    target = FailingTarget()
    mirror = _Mirror(target)
    try:
        with pytest.raises(exception, match="real I/O failure"):
            mirror.flush()
    finally:
        target.fail = False
        mirror.close()


def test_mirror_handles_target_closing_during_flush():
    class ClosingTarget(io.StringIO):
        def flush(self):
            self.close()
            raise ValueError("I/O operation on closed file")

    target = ClosingTarget()
    mirror = _Mirror(target)
    mirror.flush()
    mirror.close()
    assert target.closed
