from types import SimpleNamespace

from dvbfixer.model.diffusion.contract import AtomIdentity, ResidueIdentity
from dvbfixer.model.diffusion.gap_neighbors import source_gap_components


def test_neighbor_graph_is_transitive_but_keeps_distant_case_sensitive_chain(monkeypatch):
    import dvbfixer.model.diffusion.gap_neighbors as module

    chains = ("A", "B", "C", "a")
    gaps = tuple(SimpleNamespace(chain=c, target_interval=SimpleNamespace(start=1, stop=2)) for c in chains)
    placements = tuple(SimpleNamespace(chain=c, observed_target_indices=(0,),
        observed_residues=(ResidueIdentity(c, "1", ""),)) for c in chains)
    points = {AtomIdentity(c, "1", "", name): (x, y, 0.0)
        for c, x in zip(chains, (0.0, 5.0, 10.0, 30.0), strict=True)
        for name, y in (("N", 0.0), ("CA", 0.2), ("C", 0.4))}
    monkeypatch.setattr(module, "source_coordinates", lambda *_: points)
    groups, edges = source_gap_components(SimpleNamespace(gaps=gaps, sequence_placements=placements), "")
    assert groups == ((0, 1, 2), (3,))
    assert {e["gaps"] for e in edges} == {(0, 1), (1, 2)}
