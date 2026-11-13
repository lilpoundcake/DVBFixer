"""Packaging controls only; compact wrapper is NOT a new native scientific run."""
import subprocess
import sys

from scripts.openmm_cpu_sign_minimal import make_system, measure
from scripts.research_openmm_cpu_sign import analytic


def test_compact_system_is_byte_identical_to_frozen_synthetic_model():
    import openmm
    from scripts.research_openmm_cpu_sign import make_system as original

    assert openmm.XmlSerializer.serialize(make_system()) == openmm.XmlSerializer.serialize(original(5e7, True))


def test_compact_scalar_analysis_matches_independent_math_without_native_call():
    import numpy as np

    p = [[-0.0500025, 0., 0.], [0.0500025, 0., 0.], [10., 10., 0.]]
    physical, gp = analytic(p, 5e7, True, 1e6)
    wrong = (physical-gp).astype(np.float32).astype(np.float64).flatten().tolist()
    report = measure([v for point in p for v in point], wrong, {'restraint strength': 1e6}, True)
    assert not report['physical']['original_constraint_allowance_passed']
    assert report['source_sign_float32']['maximum_error'] == 0


def test_compact_entrypoint_is_nonexecuting_by_default():
    result = subprocess.run([sys.executable, 'scripts/openmm_cpu_sign_minimal.py'],
                            check=True, capture_output=True, text=True)
    assert 'no Context or minimizer call' in result.stdout
