"""Source-bound spatial gap grouping, without refinement/runtime imports."""
from __future__ import annotations

import numpy as np

from dvbfixer.model.diffusion.contract import DiffusionContractError, ResidueIdentity
from dvbfixer.model.diffusion.protein_context import source_coordinates


def source_gap_components(request, source_text):
    """Transitive components of 8-residue flanks with heavy contacts <=6 Å.

    Never split a component to fit a token budget. Components are evidence of
    proximity, not a guarantee that initially separated generated loops cannot
    meet. Use original observations, never predictions, to define this graph.
    """
    source = source_coordinates(request, source_text.encode())
    points = []
    for gap in request.gaps:
        placement = next(p for p in request.sequence_placements if p.chain == gap.chain)
        residues = {r for i, r in zip(placement.observed_target_indices, placement.observed_residues, strict=True)
            if gap.target_interval.start - 8 <= i < gap.target_interval.start
            or gap.target_interval.stop <= i < gap.target_interval.stop + 8}
        axis = np.asarray([xyz for a, xyz in source.items()
            if ResidueIdentity(a.chain, a.residue_number, a.insertion_code) in residues])
        if len(axis) < 3:
            raise DiffusionContractError("spatial grouping lacks observed flank evidence")
        points.append(axis)
    groups = [{i} for i in range(len(points))]
    edges = []
    for i, first in enumerate(points):
        for j in range(i):
            distance = float(np.linalg.norm(first[:, None, :] - points[j][None, :, :], axis=2).min())
            if distance <= 6.0:
                edges.append({"gaps": (j, i), "minimum_flank_distance_angstrom": distance})
                left = next(g for g in groups if i in g)
                right = next(g for g in groups if j in g)
                if left is not right:
                    left.update(right)
                    groups.remove(right)
    return tuple(sorted(tuple(sorted(g)) for g in groups)), tuple(edges)
