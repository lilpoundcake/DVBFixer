"""Tests for backend-neutral constrained-sampler conformance."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pytest

from dvbfixer.model.diffusion.contract import AtomIdentity
from dvbfixer.model.diffusion.sampler import (
    SamplerCapabilities,
    SamplerConformance,
    SamplingAblationMode,
    assess_sampler_conformance,
)

REPO_ROOT = Path(__file__).resolve().parent.parent


def test_template_conditioning_does_not_claim_reinjection_control() -> None:
    capabilities = SamplerCapabilities(False, False, False, False)

    report = assess_sampler_conformance(
        capabilities,
        SamplingAblationMode.TEMPLATE_ONLY,
    )

    assert report.supported
    assert not capabilities.supports_per_step_reinjection()
    assert capabilities.supports_mode(SamplingAblationMode.TEMPLATE_ONLY)


def test_reinjection_requires_every_per_step_control() -> None:
    capabilities = SamplerCapabilities(
        mutable_state_each_step=True,
        identity_mapping_each_step=False,
        fixed_coordinate_overwrite_each_step=True,
        localized_boundary_refinement=False,
    )

    report = assess_sampler_conformance(
        capabilities,
        SamplingAblationMode.REINJECTION,
    )

    assert not report.supported
    assert report.missing_capabilities == ("identity-mapping-each-step",)
    assert not capabilities.supports_mode(SamplingAblationMode.REINJECTION)


def test_boundary_refinement_is_a_separate_required_capability() -> None:
    reinjection_only = SamplerCapabilities(True, True, True, False)
    full = SamplerCapabilities(True, True, True, True)

    missing = assess_sampler_conformance(
        reinjection_only,
        SamplingAblationMode.REINJECTION_BOUNDARY_REFINEMENT,
    )
    supported = assess_sampler_conformance(
        full,
        SamplingAblationMode.REINJECTION_BOUNDARY_REFINEMENT,
    )

    assert missing.missing_capabilities == ("localized-boundary-refinement",)
    assert supported.supported
    assert full.supports_mode(
        SamplingAblationMode.REINJECTION_BOUNDARY_REFINEMENT
    )


def test_conformance_record_rejects_inconsistent_decisions() -> None:
    with pytest.raises(ValueError, match="supported"):
        SamplerConformance(
            SamplingAblationMode.REINJECTION,
            supported=True,
            missing_capabilities=("mutable-state-each-step",),
        )
    with pytest.raises(ValueError, match="unsupported"):
        SamplerConformance(
            SamplingAblationMode.REINJECTION,
            supported=False,
        )


def test_protenix_callback_synchronizes_and_reinjects_on_device() -> None:
    torch = pytest.importorskip("torch")
    module_path = REPO_ROOT / "deploy/protenix-v1/reinjection.py"
    spec = importlib.util.spec_from_file_location("protenix_reinjection", module_path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    identities = tuple(AtomIdentity("A", str(index + 1), "", "CA") for index in range(5))
    target = np.asarray(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.5, 0.5, 0.5],
        ]
    )
    rotation = np.asarray([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    sampled = target @ rotation.T + np.asarray([4.0, -3.0, 2.0])
    state = torch.as_tensor(np.stack([sampled, sampled]), dtype=torch.float64)
    fixed = {identities[index]: target[index] for index in range(4)}

    reinjector = module.FixedAtomReinjector(identities, fixed)
    synchronized = reinjector(state, step_index=0, step_count=2)

    expected = torch.as_tensor(target, dtype=torch.float64)
    assert torch.equal(synchronized[:, :4], expected[:4].expand(2, -1, -1))
    assert torch.allclose(synchronized[:, 4], expected[4].expand(2, -1), atol=1e-12)
