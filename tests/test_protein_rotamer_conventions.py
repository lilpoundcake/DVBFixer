"""Independent Biopython sign/identity checks; no real rotamer table is used."""
from __future__ import annotations

import math
from pathlib import Path

import pytest
from Bio.PDB.vectors import Vector, calc_dihedral

from dvbfixer.diagnose.geometry import signed_dihedral_degrees
from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef
from dvbfixer.protein_completion.backbone import backbone_angles
from dvbfixer.protein_completion.contracts import CompletionRequest, Observation
from dvbfixer.protein_completion.geometry import extend, l_cb

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize("resid", (20, 57, 65))
def test_fixture_phi_psi_and_chi_signs_agree_with_independent_biopython(resid):
    source = (ROOT/"tests/fixtures/protein_context/1UBQ.pdb").read_text()
    observations = []
    for line in source.splitlines():
        if line.startswith("ATOM  ") and line[21] == "A" and int(line[22:26]) in (resid-1, resid, resid+1):
            component = ComponentInstanceRef(1, "A", line[22:26].strip())
            observations.append(Observation(ExactAtomRef(component, line[12:16].strip()),
                                            line[76:78].strip(), tuple(float(line[i:i+8]) for i in (30, 38, 46))))
    residue = ComponentInstanceRef(1, "A", str(resid))
    n, ca, c, cb, og = (ExactAtomRef(residue, name) for name in ("N", "CA", "C", "CB", "OG"))
    previous = ExactAtomRef(ComponentInstanceRef(1, "A", str(resid-1)), "C")
    following = ExactAtomRef(ComponentInstanceRef(1, "A", str(resid+1)), "N")
    req = CompletionRequest("a"*64, residue, "SER", tuple(observations), (), ((previous, n), (c, following)))
    positions = {o.atom: o.xyz for o in observations}
    phi, psi = backbone_angles(req)
    for angle, atoms in ((phi, (previous, n, ca, c)), (psi, (n, ca, c, following)),
                         (signed_dihedral_degrees(*(positions[a] for a in (n, ca, cb, og))), (n, ca, cb, og))):
        bio = math.degrees(calc_dihedral(*(Vector(*positions[a]) for a in atoms)))
        assert angle == pytest.approx(bio, abs=1e-9)


@pytest.mark.parametrize("chi", (-179., -66., 62., 178.))
def test_own_coordinate_builder_obeys_independent_signed_chi(chi):
    n, ca, c = (-1.46, 0., 0.), (0., 0., 0.), (0.55, 1.43, 0.)
    cb = l_cb(n, ca, c)
    og = extend(n, ca, cb, 1.417, 110.8, chi)
    independent = math.degrees(calc_dihedral(*(Vector(*p) for p in (n, ca, cb, og))))
    assert abs((independent-chi+180) % 360-180) < 1e-9
