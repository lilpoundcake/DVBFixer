"""Backbone torsions from explicit observed peptide links, never residue numbering."""
from __future__ import annotations

import math

from dvbfixer.diagnose.geometry import signed_dihedral_degrees
from dvbfixer.domain.structure_identity import ExactAtomRef
from dvbfixer.protein_completion.contracts import CompletionRequest
from dvbfixer.protein_completion.rotamers import RotamerLibraryError


def backbone_angles(request: CompletionRequest) -> tuple[float, float]:
    observed = {o.atom: o for o in request.observations}
    n, ca, c = (ExactAtomRef(request.residue, name) for name in ("N", "CA", "C"))
    if any(atom not in observed for atom in (n, ca, c)):
        raise RotamerLibraryError("backbone atoms unavailable")

    def partner(atom: ExactAtomRef, name: str) -> ExactAtomRef:
        partners = [b if a == atom else a for a, b in request.links if atom in (a, b)]
        if (len(partners) != 1 or partners[0].name != name
                or partners[0].component == request.residue or partners[0] not in observed):
            raise RotamerLibraryError("phi/psi require unambiguous explicit observed peptide neighbours")
        other = partners[0]
        if (observed[other].element != name or other.component.model != request.residue.model
                or other.component.alternate_location or observed[other].occupancy != 1):
            raise RotamerLibraryError("invalid observed peptide neighbour")
        if not 1.0 <= math.dist(observed[atom].xyz, observed[other].xyz) <= 1.8:
            raise RotamerLibraryError("invalid observed peptide bond length")
        return other

    previous, following = partner(n, "C"), partner(c, "N")
    if previous.component == following.component:
        raise RotamerLibraryError("cyclic or self-joined backbone requires separate support")
    phi = signed_dihedral_degrees(*(observed[a].xyz for a in (previous, n, ca, c)))
    psi = signed_dihedral_degrees(*(observed[a].xyz for a in (n, ca, c, following)))
    if phi is None or psi is None:
        raise RotamerLibraryError("degenerate backbone torsions")
    return phi, psi
