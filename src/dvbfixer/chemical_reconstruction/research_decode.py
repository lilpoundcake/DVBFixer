"""Bounded private known-DTO JSON transport, not a public/pickle/.dat loader.

Only trusted Python annotations reachable from the explicit root DTOs below are
used. JSON never names a module/class. File/hash/source provenance belongs to the
caller; successful decoding establishes neither scientific acceptance nor trust.
"""
from __future__ import annotations

import json
import math
import types
from dataclasses import fields, is_dataclass
from enum import Enum
from typing import Literal, TypeVar, Union, cast, get_args, get_origin, get_type_hints

from .conformer_search import ConformerSearchPlan, ConformerSearchResult
from .contracts import Observation, PreparedContextRequest, ReconstructionRequest
from .declared_peptide import DeclaredPeptideAuthority
from .heavy_conditioning import HeavyConditionedRequest
from .joint_hydrogen_geometry import JointMetricGeometryPlan, JointSourceMetricGeometryPlan
from .metric_geometry import GeneratedMetricGeometryPlan, SourceMetricGeometryPlan
from .observed_peptide import ObservedPeptideProbe
from .peptide_state import LocalPeptideAuthority
from .planar_port import PlanarPortAuthority, PlanarPortDeclaration
from .projection_contracts import ProjectionPlan, ProjectionResult
from .proposals import GeometryPlan, ProposalResult
from .staged_contracts import StagedResult

MAX_BYTES = 64 * 1024 * 1024
MAX_NODES = 2_000_000
MAX_DEPTH = 64
T = TypeVar("T")
ROOTS = (ReconstructionRequest, DeclaredPeptideAuthority, ConformerSearchResult, ObservedPeptideProbe, StagedResult,
         ProjectionResult, PlanarPortAuthority, PlanarPortDeclaration, Observation, GeometryPlan, ProposalResult)
VARIANTS: dict[type, tuple[type, ...]] = {
    ReconstructionRequest: (ReconstructionRequest, PreparedContextRequest, HeavyConditionedRequest),
    LocalPeptideAuthority: (LocalPeptideAuthority, DeclaredPeptideAuthority),
    GeometryPlan: (GeometryPlan, GeneratedMetricGeometryPlan, SourceMetricGeometryPlan,
                   JointMetricGeometryPlan, JointSourceMetricGeometryPlan, ConformerSearchPlan, ProjectionPlan),
    GeneratedMetricGeometryPlan: (GeneratedMetricGeometryPlan, SourceMetricGeometryPlan,
                                  JointMetricGeometryPlan, JointSourceMetricGeometryPlan),
}


class ResearchDecodeError(ValueError):
    """Malformed, unknown, ambiguous or over-budget private evidence."""


def _pairs(items: list[tuple[str, object]]) -> dict[str, object]:
    result: dict[str, object] = {}
    for key, value in items:
        if key in result:
            raise ResearchDecodeError("duplicate-JSON-key")
        result[key] = value
    return result


def _constant(value: str) -> object:
    raise ResearchDecodeError("nonfinite-JSON-number:" + value)


def parse_research_json(data: bytes) -> object:
    """Bound bytes before JSON parsing, then bound every node before DTO decoding."""
    if type(data) is not bytes or len(data) > MAX_BYTES:
        raise ResearchDecodeError("research-JSON-byte-budget-or-type")
    try:
        value: object = json.loads(data.decode("utf-8"), object_pairs_hook=_pairs, parse_constant=_constant)
    except (UnicodeError, json.JSONDecodeError, RecursionError) as exc:
        raise ResearchDecodeError("malformed-research-JSON") from exc
    pending = [(value, 0)]
    count = 0
    while pending:
        node, depth = pending.pop()
        count += 1
        if count > MAX_NODES or depth > MAX_DEPTH:
            raise ResearchDecodeError("research-JSON-node-or-depth-budget")
        if type(node) is dict:
            pending.extend((x, depth + 1) for x in node.values())
        elif type(node) is list:
            pending.extend((x, depth + 1) for x in node)
        elif type(node) in {int, float}:
            try:
                finite = math.isfinite(cast(float, node))
            except OverflowError as exc:
                raise ResearchDecodeError("out-of-range-JSON-number") from exc
            if not finite:
                raise ResearchDecodeError("nonfinite-JSON-number")
    return value


def _known_types() -> frozenset[type]:
    known: set[type] = set()
    pending: list[object] = [*ROOTS, *(t for row in VARIANTS.values() for t in row)]
    while pending:
        annotation = pending.pop()
        if get_origin(annotation) is not None:
            pending.extend(get_args(annotation))
        elif isinstance(annotation, type) and annotation not in known:
            known.add(annotation)
            if is_dataclass(annotation):
                pending.extend(get_type_hints(annotation).values())
    return frozenset(known)


KNOWN = _known_types()


def _restore(annotation: object, value: object) -> object:
    origin, args = get_origin(annotation), get_args(annotation)
    if origin in {Union, types.UnionType}:
        for option in args:
            try:
                return _restore(option, value)
            except (ValueError, TypeError):
                continue
        raise ResearchDecodeError("no-exact-DTO-union-variant")
    if origin is Literal:
        if not any(type(value) is type(option) and value == option for option in args):
            raise ResearchDecodeError("wrong-DTO-literal")
        return value
    if origin in {tuple, frozenset}:
        if type(value) is not list:
            raise ResearchDecodeError("DTO-sequence-needs-JSON-array")
        if origin is frozenset:
            result = tuple(_restore(args[0], item) for item in value)
            if len(set(result)) != len(result):
                raise ResearchDecodeError("duplicate-DTO-set-item")
            return frozenset(result)
        if len(args) == 2 and args[1] is Ellipsis:
            return tuple(_restore(args[0], item) for item in value)
        if len(value) != len(args):
            raise ResearchDecodeError("wrong-DTO-tuple-length")
        return tuple(_restore(t, item) for t, item in zip(args, value, strict=True))
    if not isinstance(annotation, type) or annotation not in KNOWN:
        raise ResearchDecodeError("unknown-DTO-type")
    if issubclass(annotation, Enum):
        if type(value) is not str:
            raise ResearchDecodeError("wrong-DTO-enum-scalar")
        return annotation(value)
    if is_dataclass(annotation):
        if type(value) is not dict:
            raise ResearchDecodeError("DTO-record-needs-JSON-object")
        variants = VARIANTS.get(annotation, (annotation,))
        matches = [t for t in variants if set(value) == {f.name for f in fields(t)}]
        if len(matches) != 1:
            raise ResearchDecodeError("unknown-or-ambiguous-DTO-fields")
        selected = matches[0]
        hints = get_type_hints(selected)
        return selected(**{f.name: _restore(hints[f.name], value[f.name]) for f in fields(selected)})
    if annotation is type(None) and value is None:
        return None
    if annotation in {str, int, bool} and type(value) is annotation:
        return value
    # Preserve original int-vs-float representation in coordinate/digest records;
    # bool is explicitly NOT a numeric value, and nonfinite floats are forbidden.
    if annotation is float and type(value) in {int, float} and math.isfinite(cast(float, value)):
        return value
    raise ResearchDecodeError("wrong-DTO-scalar-type")


def decode_research_record(data: bytes, expected: type[T]) -> T:
    """Strict exact-field transport for a known root DTO; no source/digest bypass."""
    if expected not in ROOTS:
        raise ResearchDecodeError("unknown-research-root-DTO")
    value = _restore(expected, parse_research_json(data))
    if not isinstance(value, expected):
        raise ResearchDecodeError("research-root-type-mismatch")
    return value
