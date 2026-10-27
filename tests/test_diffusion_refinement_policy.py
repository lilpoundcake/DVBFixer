"""Pure request-policy compatibility and fail-closed decoding tests."""

from dataclasses import replace

import pytest

from dvbfixer.model.diffusion.contract import (
    BackendOption,
    DiffusionContractError,
    DiffusionRequest,
)
from dvbfixer.model.diffusion.refinement_policy import RefinementPolicy, read_refinement_policy
from tests.test_diffusion_contract import _request


def test_historical_request_roundtrip_keeps_bytes_and_does_not_enable_motion():
    request = _request()
    original = request.to_json()
    assert read_refinement_policy(request.backend_options) is None
    assert DiffusionRequest.from_json(original).to_json() == original
    assert "geometry_refinement" not in original


@pytest.mark.parametrize("anchors", [True, False])
@pytest.mark.parametrize("bfactors", [True, False])
def test_explicit_switches_roundtrip_and_bind_limits(anchors, bfactors):
    policy = RefinementPolicy(anchors, bfactors)
    request = _request()
    updated = replace(request, backend_options=(*request.backend_options, *policy.backend_options()))
    decoded = read_refinement_policy(DiffusionRequest.from_json(updated.to_json()).backend_options)
    assert decoded == policy
    assert updated.to_json() != request.to_json()
    assert policy.movement_limits.backbone_rmsd_angstrom_max == 0.5
    assert policy.movement_limits.ca_displacement_angstrom_max == 0.5
    assert policy.movement_limits.heavy_displacement_angstrom_max == 1.25


def test_new_policy_defaults_enable_stage_two_and_b_weighting():
    assert RefinementPolicy().relax_anchors and RefinementPolicy().bfactor_weighting


@pytest.mark.parametrize("mutation", [
    ("policy", "future-v2"), ("relax_anchors", "yes"),
    ("bfactor_weighting", "TRUE"), ("unknown", "true"),
])
def test_request_rejects_unknown_versions_keys_and_noncanonical_booleans(mutation):
    name, value = mutation
    options = tuple(option for option in RefinementPolicy().backend_options()
                    if option.name != "geometry_refinement." + name)
    options += (BackendOption("geometry_refinement." + name, value),)
    with pytest.raises(DiffusionContractError):
        replace(_request(), backend_options=options)


def test_partial_and_duplicate_policy_cannot_authorize_motion():
    options = RefinementPolicy().backend_options()
    for malformed in (options[:1], options[1:], (*options, options[0])):
        with pytest.raises(DiffusionContractError):
            read_refinement_policy(malformed)


@pytest.mark.parametrize("value", [1, "false", None])
def test_policy_constructor_rejects_truthy_non_booleans(value):
    with pytest.raises(DiffusionContractError):
        RefinementPolicy(relax_anchors=value)
