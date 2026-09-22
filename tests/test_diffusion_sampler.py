"""Tests for backend-neutral constrained-sampler conformance."""

from __future__ import annotations

import pytest

from dvbfixer.model.diffusion.sampler import (
    SamplerCapabilities,
    SamplerConformance,
    SamplingAblationMode,
    assess_sampler_conformance,
)


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
