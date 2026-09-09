"""Pure policy vocabulary for choosing a retained component's parameter route."""

from __future__ import annotations

from enum import StrEnum


class ParameterizationRoute(StrEnum):
    EXACT_NATIVE_TEMPLATE = "exact-native-template"
    GAFF_CANDIDATE = "gaff-candidate"
    USER_TEMPLATE = "user-template"
    UNSUPPORTED_COMPLEX_COFACTOR = "unsupported-complex-cofactor"
    EXPLICIT_STRIP = "explicit-strip"


COMPLEX_COFACTORS = frozenset({"HEM", "HEME", "FAD", "FMN", "NAD", "NAP"})


def classify_parameterization(
    residue_name: str,
    *,
    exact_native_match: bool = False,
    user_template: bool = False,
    strip: bool = False,
) -> ParameterizationRoute:
    """Classify a component without guessing compatibility from its name alone."""
    if strip:
        return ParameterizationRoute.EXPLICIT_STRIP
    if user_template:
        return ParameterizationRoute.USER_TEMPLATE
    if exact_native_match:
        return ParameterizationRoute.EXACT_NATIVE_TEMPLATE
    if residue_name.upper() in COMPLEX_COFACTORS:
        return ParameterizationRoute.UNSUPPORTED_COMPLEX_COFACTOR
    return ParameterizationRoute.GAFF_CANDIDATE
