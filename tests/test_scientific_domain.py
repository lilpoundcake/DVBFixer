from __future__ import annotations

import pytest

from dvbfixer.domain.parameterization import (
    ParameterizationRoute,
    classify_parameterization,
)
from dvbfixer.domain.structure_identity import allocate_chain_ids


def test_chain_allocator_reserves_coordinate_and_metadata_ids() -> None:
    assert allocate_chain_ids(4, alphabet=list("ABCDEFGHIL"), reserved={"A", "H", "L"}) == (
        "B", "C", "D", "E"
    )


def test_chain_allocator_fails_instead_of_reusing_an_identity() -> None:
    with pytest.raises(ValueError, match="unique molecule chain IDs"):
        allocate_chain_ids(2, alphabet=["A", "B"], reserved={"A"})


def test_complex_cofactors_do_not_fall_into_generic_gaff() -> None:
    assert classify_parameterization("HEM") is ParameterizationRoute.UNSUPPORTED_COMPLEX_COFACTOR
    assert classify_parameterization("FAD") is ParameterizationRoute.UNSUPPORTED_COMPLEX_COFACTOR
    assert classify_parameterization("PLM") is ParameterizationRoute.GAFF_CANDIDATE
    assert classify_parameterization("HEM", user_template=True) is ParameterizationRoute.USER_TEMPLATE
