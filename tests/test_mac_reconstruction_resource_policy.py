from scripts.run_mac_reconstruction_guarded import resource_stop_reason


def test_normal_pressure_does_not_imply_a_swap_growth_failure():
    assert resource_stop_reason(elapsed=20, timeout=1800, peak_rss=2*1024**3,
        pressure=1, critical_seconds=0, low_headroom_seconds=0) is None


def test_transient_pressure_and_low_headroom_are_not_sustained():
    assert resource_stop_reason(elapsed=20, timeout=1800, peak_rss=2*1024**3,
        pressure=4, critical_seconds=4.9, low_headroom_seconds=4.9) is None


def test_critical_pressure_and_resource_caps_still_fail_closed():
    base = dict(elapsed=20, timeout=1800, peak_rss=2*1024**3,
        pressure=4, critical_seconds=5, low_headroom_seconds=0)
    assert resource_stop_reason(**base) == "sustained-critical-host-memory-pressure"
    assert resource_stop_reason(**{**base, "pressure": 1, "critical_seconds": 0,
        "low_headroom_seconds": 5}) == "sustained-host-memory-headroom"
    assert resource_stop_reason(**{**base, "peak_rss": 9*1024**3}) == "aggregate-process-rss-limit"
    assert resource_stop_reason(**{**base, "elapsed": 1801}) == "wall-time-limit"
