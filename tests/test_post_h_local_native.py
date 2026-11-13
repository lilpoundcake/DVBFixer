"""Native enrollment transport controls, no scientific native execution."""
import pytest
from scripts.research_post_h_local_native import chemistry

from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef


def test_shared_chemistry_preserves_case_insertion_and_resolves_only_known_aliases():
    from scripts.hydrogen_initialization_guard import Atom

    a = Atom(ExactAtomRef(ComponentInstanceRef(1, 'A', '2', 'B'), 'CA'), 'GLH', 'C')
    parent = Atom(a.identity, 'GLU', 'C')
    other = Atom(ExactAtomRef(ComponentInstanceRef(1, 'a', '2', 'B'), 'CA'), 'GLU', 'C')
    assert chemistry([a]) == chemistry([parent])
    assert chemistry([a]) != chemistry([other])


def test_phase_transition_rejects_anchor_retargeting():
    import openmm
    from scripts.research_thr2_hswap import replace_restraint_parameters

    def system(k, x):
        value = openmm.System()
        value.addParticle(12.)
        f = openmm.CustomExternalForce('0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)')
        for name in ('k', 'x0', 'y0', 'z0'):
            f.addPerParticleParameter(name)
        f.addParticle(0, [k, x, 0., 0.])
        value.addForce(f)
        return value

    with pytest.raises(ValueError, match='same-anchor'):
        replace_restraint_parameters(system(100., 0.), system(10., 1.))
