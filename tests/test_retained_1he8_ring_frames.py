from copy import deepcopy

import numpy as np
import pytest
from scripts.audit_retained_1he8_ring_frames import focus_probe, summarize_frames

from dvbfixer.model.diffusion.contract import AtomIdentity


def fixture():
    ring = tuple(AtomIdentity("A", "972", "B", f"C{i}") for i in range(6))
    witness = (AtomIdentity("a", "1044", "", "CA"), AtomIdentity("a", "1044", "", "CB"))
    points = {atom: (float(np.cos(i * np.pi / 3)), float(np.sin(i * np.pi / 3)), 0.) for i, atom in enumerate(ring)}
    points.update({witness[0]: (.2, .1, -1.), witness[1]: (.2, .1, 1.)})
    bonds = frozenset([frozenset((a, b)) for a, b in zip(ring, ring[1:] + ring[:1])] + [frozenset(witness)])
    return points, bonds, ring, witness


def test_known_pair_is_diagnosed_without_coordinate_mutation_or_acceptance_gate():
    points, bonds, ring, witness = fixture()
    original = deepcopy(points)
    result = focus_probe(points, bonds, ring, witness, frozenset(ring))
    assert points == original
    assert result["status"] == "defined-surface"
    assert result["witness_piercings"] == 1
    assert result["probe"]["piercings"][0]["bond"][0]["chain"] == "a"
    assert result["probe"]["ring"][0]["insertion_code"] == "B"
    assert not result["publication_gate"]


def test_missing_ring_is_visible_only_for_explicit_partial_source():
    points, bonds, ring, witness = fixture()
    del points[ring[0]]
    with pytest.raises(ValueError, match="lacks topology"):
        focus_probe(points, bonds, ring, witness, frozenset())
    result = focus_probe(points, bonds, ring, witness, frozenset(), partial_source=True)
    assert result["status"] == "missing-observed-ring"
    assert result["probe"] is None
    assert len(result["missing_ring_atoms"]) == 1


def test_nonplanar_ring_is_undefined_not_a_clean_surface():
    points, bonds, ring, witness = fixture()
    points[ring[0]] = (1., 0., 2.)
    result = focus_probe(points, bonds, ring, witness, frozenset())
    assert result["status"] == "undefined-surface"
    assert result["probe"]["undefined_reason"]
    assert result["witness_piercings"] == 0


def test_source_with_missing_witness_is_not_claimed_as_tested_pair():
    points, bonds, ring, witness = fixture()
    del points[witness[1]]
    result = focus_probe(points, bonds, ring, witness, frozenset(), partial_source=True)
    assert not result["witness_atoms_present"]
    assert result["missing_topology_coordinate_atoms"] == 1
    assert result["available_coordinate_edges"] == len(bonds) - 1


@pytest.mark.parametrize("removed", ["ring", "witness"])
def test_missing_explicit_topology_is_not_inferred_from_nearby_atoms(removed):
    points, bonds, ring, witness = fixture()
    edge = frozenset(ring[:2]) if removed == "ring" else frozenset(witness)
    with pytest.raises(ValueError, match="edges absent"):
        focus_probe(points, bonds - {edge}, ring, witness, frozenset())


def test_all_missing_source_does_not_fabricate_coordinates():
    _, bonds, ring, witness = fixture()
    result = focus_probe({}, bonds, ring, witness, frozenset(), partial_source=True)
    assert result["status"] == "missing-observed-ring"
    assert not result["witness_atoms_present"]
    assert result["available_coordinate_edges"] == 0


def test_summary_keeps_undefined_missing_and_defined_separate():
    points, bonds, ring, witness = fixture()
    positive = focus_probe(points, bonds, ring, witness, frozenset())
    points[ring[0]] = (1., 0., 2.)
    undefined = focus_probe(points, bonds, ring, witness, frozenset())
    missing = focus_probe({}, bonds, ring, witness, frozenset(), partial_source=True)
    result = summarize_frames([{"screen": row} for row in (positive, undefined, missing)])
    assert (result["frames"], result["defined_surfaces"], result["undefined_surfaces"], result["missing_observed_rings"]) == (3, 1, 1, 1)
    assert result["frames_with_piercings"] == result["frames_with_witness_piercings"] == 1


def test_invalid_coordinate_and_topology_refused():
    points, bonds, ring, witness = fixture()
    points[witness[0]] = (np.nan, 0., 0.)
    with pytest.raises(ValueError, match="finite"):
        focus_probe(points, bonds, ring, witness, frozenset())
    with pytest.raises(ValueError, match="Invalid topology"):
        focus_probe(points, bonds | {frozenset((ring[0],))}, ring, witness, frozenset())
