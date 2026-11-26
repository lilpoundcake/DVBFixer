"""Shared metadata-only adapter for assert_all_l; no bonds or coordinates built."""
from __future__ import annotations

from collections.abc import Mapping

from dvbfixer.domain.structure_identity import ComponentInstanceRef

from .contracts import Observation, _text, _tuple


class MetadataTopologyError(ValueError):
    """Unrepresentable metadata, not a request to infer or repair chemistry."""


def assert_metadata_all_l(observations: tuple[Observation, ...],
                          residue_names: Mapping[ComponentInstanceRef, str]) -> None:
    """Group exact residues before adding atoms; generated atoms may arrive last.

    Grouping changes no coordinates. Models/chain case/insertion/occurrence remain
    distinct even though OpenMM's detector prints only a subset of that identity.
    The existing detector and its skip/deadband rules are the ONLY chirality rule.
    """
    from openmm import Vec3
    from openmm.app import Topology, element
    from openmm.unit import angstrom

    from dvbfixer.ffutils.geometry import assert_all_l

    _tuple(observations, Observation, "metadata-only coordinates", maximum=100000)
    if (len({o.identity for o in observations}) != len(observations)
            or set(residue_names) != {o.identity.component for o in observations}):
        raise MetadataTopologyError("metadata-topology-exact-identity-inventory")
    groups: dict[tuple[int, str], dict[ComponentInstanceRef, list[Observation]]] = {}
    for observation in observations:
        component = observation.identity.component
        groups.setdefault((component.model, component.chain_id), {}).setdefault(component, []).append(observation)
    topology = Topology()
    positions = []
    for (_, chain_id), residues in groups.items():
        chain = topology.addChain(chain_id)
        for component, atoms in residues.items():
            name = residue_names[component]
            _text(name, "actual source residue name")
            residue = topology.addResidue(name, chain, id=component.sequence_number, insertionCode=component.insertion_code)
            for observation in atoms:
                try:
                    atom_element = element.get_by_symbol(observation.element)
                except KeyError as exc:
                    raise MetadataTopologyError("metadata-topology-unknown-element:"+observation.element) from exc
                topology.addAtom(observation.identity.name, atom_element, residue)
                positions.append(Vec3(*observation.angstrom)*angstrom)
    assert_all_l(topology, positions)
