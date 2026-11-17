"""Private native-disulfide target policy; never rewrite frozen source links.

Chemical target selection must be explicitly approved and source-bound by the
caller. This module validates atom ownership and graph consistency, not the
biochemical truth of an annotation or qualification of restored coordinates.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass

from dvbfixer.model.diffusion.contract import AtomIdentity
from dvbfixer.model.diffusion.rfd3_custody_records import CustodyError, Pair


@dataclass(frozen=True)
class NativeGraphAudit:
    explicit_source_pairs: int
    source_inferred_observed_pairs: int
    target_pairs: int
    additional_active_pairs: tuple[Pair, ...]
    immutable_observed_pairs: tuple[Pair, ...]
    cysteine_endpoints: int
    source_annotations_rewritten: bool = False
    scientific_acceptance: bool = False


def validate_native_graph(
    *, cysteines: frozenset[AtomIdentity], observed: frozenset[AtomIdentity],
    active: frozenset[AtomIdentity], explicit_source_pairs: tuple[Pair, ...],
    source_inferred_pairs: tuple[Pair, ...], target_pairs: tuple[Pair, ...],
) -> NativeGraphAudit:
    """Every CYS has one partner; preserve source links and exact mobile ownership.

    Source-inferred pairs must have two actual observed endpoints. An additional
    target pair may involve an active missing atom, never an unapproved inactive
    proposal atom. Input pairs are identity-based and canonicalized, not guessed
    from distances or cysteine order. No graph edges are added here.
    """
    if not cysteines or observed & active or not active:
        raise CustodyError("native graph requires disjoint observed/active identities")
    if not cysteines <= observed | active or any(a.atom_name != "SG" for a in cysteines):
        raise CustodyError("native graph cysteine inventory is incomplete or not SG")

    def checked(pairs):
        if any(len(p) != 2 or p[0] == p[1] or not set(p) <= cysteines for p in pairs):
            raise CustodyError("native graph has invalid or non-CYS endpoints")
        normalized = tuple(tuple(sorted(p)) for p in pairs)
        if len(set(normalized)) != len(normalized):
            raise CustodyError("native graph has duplicate endpoint pairs")
        return set(normalized)

    explicit, inferred, target = map(checked, (explicit_source_pairs, source_inferred_pairs, target_pairs))
    if not explicit <= inferred <= target:
        raise CustodyError("native target lost an explicit or source-inferred bond")
    if any(not set(p) <= observed for p in inferred):
        raise CustodyError("source-only inferred graph includes a missing-coordinate atom")
    extra = target - inferred
    if any(not set(p) & active for p in extra):
        raise CustodyError("new target bond between observations lacks source-only support")
    counts = Counter(atom for p in target for atom in p)
    if set(counts) != cysteines or any(count != 1 for count in counts.values()):
        raise CustodyError("native target must assign exactly one partner to every CYS")
    return NativeGraphAudit(len(explicit), len(inferred), len(target), tuple(sorted(extra)),
                            tuple(sorted(inferred)), len(cysteines))
