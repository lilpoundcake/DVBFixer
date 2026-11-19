"""Private exact inventory transport; naming alone is not chemical authority."""
from __future__ import annotations

import re
from dataclasses import dataclass

from dvbfixer.domain.force_field_naming import ForceFieldTarget, NamingProfile, atom_name_policy

from .contracts import digest_record
from .peptide_state import PeptideStateError, TemplateInventory


def _name(value: str) -> None:
    if not isinstance(value, str) or re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9'*]{0,15}", value) is None:
        raise PeptideStateError("invalid-inventory-atom-name")


def _values(raw: TemplateInventory, target: ForceFieldTarget, profile: NamingProfile
            ) -> tuple[TemplateInventory, tuple[tuple[str, str, str], ...]]:
    if target is not ForceFieldTarget.AMBER19SB or profile is not NamingProfile.STANDARD:
        raise PeptideStateError("inventory-naming-policy-unqualified")
    if not isinstance(raw, TemplateInventory):
        raise PeptideStateError("explicit-raw-inventory-required")
    if not isinstance(raw.source, str) or not raw.source.strip() or not isinstance(raw.state_name, str):
        raise PeptideStateError("invalid-inventory-source-metadata")
    _name(raw.state_name)
    if raw.aliases:
        raise PeptideStateError("already-aliased-inventory-not-raw")
    for name in raw.atoms:
        _name(name)
    policy = atom_name_policy(source_residue_name=raw.state_name, target_residue_name=raw.state_name,
                              target=target, profile=profile, atom_names=set(raw.atoms))
    rows = []
    for name in sorted(raw.atoms):
        if name in policy:
            destination, rule = policy[name]
            rows.append((name, destination, rule.value))
        else:
            rows.append((name, name, "identity"))
    names = {source: destination for source, destination, _ in rows}
    for destination in names.values():
        _name(destination)
    if len(set(names.values())) != len(names):
        raise PeptideStateError("inventory-naming-collision")
    if not isinstance(raw.external_tokens, tuple):
        raise PeptideStateError("immutable-external-inventory-rows-required")
    external = []
    for row in raw.external_tokens:
        if not isinstance(row, tuple) or len(row) not in {2, 4, 5}:
            raise PeptideStateError("malformed-external-inventory-row")
        mapped = []
        has_external = False
        for token in row:
            if not isinstance(token, str) or not token:
                raise PeptideStateError("invalid-inventory-reference")
            if token[0] in "+-":
                # Neighbour identity/state is unknown. Only invariant backbone
                # ports are qualified here; never apply target sidechain aliases.
                if token[1:] not in {"C", "N"}:
                    raise PeptideStateError("external-inventory-name-unqualified")
                has_external = True
                mapped.append(token)
            elif token in names:
                mapped.append(names[token])
            else:
                raise PeptideStateError("unknown-local-inventory-reference")
        if not has_external:
            raise PeptideStateError("external-inventory-row-without-external-port")
        external.append(tuple(mapped))
    canonical = TemplateInventory(raw.source, raw.source_sha256, raw.state_name,
                                  frozenset(names.values()),
                                  frozenset(frozenset(names[n] for n in bond) for bond in raw.bonds),
                                  tuple((a, b) for a, b, _ in rows if a != b), tuple(external))
    return canonical, tuple(rows)


@dataclass(frozen=True)
class InventoryNameTransport:
    raw: TemplateInventory
    canonical: TemplateInventory
    target: ForceFieldTarget
    profile: NamingProfile
    mapping: tuple[tuple[str, str, str], ...]

    def __post_init__(self) -> None:
        canonical, mapping = _values(self.raw, self.target, self.profile)
        if self.canonical != canonical or self.mapping != mapping:
            raise PeptideStateError("inventory-transport-binding-mismatch")

    @property
    def digest(self) -> str:
        return digest_record(("inventory-name-transport-v1", self))


def transport_inventory(raw: TemplateInventory, *, target: ForceFieldTarget,
                        profile: NamingProfile) -> InventoryNameTransport:
    """Apply one explicit canonical owner policy simultaneously; never choose a state."""
    canonical, mapping = _values(raw, target, profile)
    return InventoryNameTransport(raw, canonical, target, profile, mapping)
