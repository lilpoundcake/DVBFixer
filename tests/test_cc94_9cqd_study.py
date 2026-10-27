"""Research metric controls; do not require neural weights or ignored inputs."""

from types import SimpleNamespace

import pytest

pytest.importorskip("numpy")
pytest.importorskip("openmm")

from scripts.cc94_9cqd_metrics import score, segment_distance

from dvbfixer.model.diffusion.contract import AtomIdentity


@pytest.mark.parametrize(
    "p0,p1,q0,q1,expected",
    [
        ((-1, 0, 0), (1, 0, 0), (0, -1, 0), (0, 1, 0), 0),
        ((-1, 0, 0), (1, 0, 0), (0, -1, 2), (0, 1, 2), 2),
        ((0, 0, 0), (1, 0, 0), (0, 3, 0), (1, 3, 0), 3),
        ((0, 0, 0), (1, 0, 0), (2, 0, 0), (3, 0, 0), 1),
        ((0, 0, 0), (2, 0, 0), (1, 0, 0), (3, 0, 0), 0),
        ((0, 0, 0), (0, 0, 0), (-1, 1, 0), (1, 1, 0), 1),
        ((0, 0, 0), (0, 0, 0), (0, 0, 2), (0, 0, 2), 2),
    ],
)
def test_analytic_segments(p0, p1, q0, q1, expected):
    assert segment_distance(p0, p1, q0, q1) == pytest.approx(expected)
    assert segment_distance(q1, q0, p1, p0) == pytest.approx(expected)


def atom_line(serial, name, chain, resid, x):
    return f"ATOM  {serial:5d} {name:>4} ALA {chain}{resid:4d}    {x:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00           C\n"


def test_proximity_to_other_chain_does_not_create_bond_exclusion():
    generated = AtomIdentity("M", "155", "", "CA")
    fixed = AtomIdentity("X", "1", "", "CA")
    source = atom_line(1, "CA", "X", 1, 0.2) + "END\n"
    raw = atom_line(1, "CA", "M", 155, 0.0) + source
    request = SimpleNamespace(
        generated_atoms=(generated,),
        fixed_atoms=(fixed,),
        retained_explicit_links=(),
        sequence_placements=(),
        gaps=(),
    )
    result = score(raw, source, request)
    assert result["generated_interchain_covalent_edges"] == 0
    assert result["authoritative_severe_clash_count"] == 1
    assert result["X"]["pairs_below_1_9_angstrom"] == 1
    assert result["X"]["minimum_distance_angstrom"] == pytest.approx(0.2)
