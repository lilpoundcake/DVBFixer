"""Optional real separate-process RESP acceptance, never inferred from an import."""

import os

import pytest


@pytest.mark.slow
def test_separate_psi4_runtime_computes_two_stage_resp(tmp_path):
    prefix = os.environ.get("DVBFIXER_TEST_PSI4_PREFIX")
    if not prefix:
        pytest.skip("not-tested: provide DVBFIXER_TEST_PSI4_PREFIX for the separate optional QM runtime")
    import numpy as np

    from dvbfixer.parametrize import _compute_resp_charges_psi4, preflight_qm_engine

    water = tmp_path / "water.pdb"
    water.write_text(
        "HETATM    1  O   HOH A   1       0.000   0.000   0.000  1.00  0.00           O\n"
        "HETATM    2  H1  HOH A   1       0.758   0.000   0.504  1.00  0.00           H\n"
        "HETATM    3  H2  HOH A   1      -0.758   0.000   0.504  1.00  0.00           H\nEND\n")
    preflight_qm_engine("psi4", psi4_env=prefix)
    charges = _compute_resp_charges_psi4(water, psi4_env=prefix, nthreads=1, memory="1GB")
    assert len(charges) == 3 and np.isfinite(charges).all()
    assert abs(sum(charges)) < 1e-5 and charges[0] < 0 and min(charges[1:]) > 0
