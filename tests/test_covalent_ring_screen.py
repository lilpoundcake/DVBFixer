from pathlib import Path

from scripts.screen_covalent_rings import canonical_cycle, screen

from dvbfixer.model.diffusion.contract import AtomIdentity


def test_cycle_normalization_does_not_depend_on_rotation_or_direction():
    cycle = tuple(AtomIdentity('a', '7', 'B', n) for n in ('CG','CD1','CE1','CZ','CE2','CD2'))
    expected = canonical_cycle(cycle)
    assert canonical_cycle(cycle[2:]+cycle[:2]) == expected
    assert canonical_cycle(cycle[::-1]) == expected
    assert all(a.chain == 'a' and a.insertion_code == 'B' for a in expected)


def test_native_ubiquitin_cycle_basis_is_a_diagnostic_control_not_qualification():
    fixture = Path(__file__).parent / 'fixtures/protein_context/1UBQ.pdb'
    text = ''.join(line for line in fixture.read_text().splitlines(True)
                   if line.startswith('ATOM  ') and int(line[22:26]) <= 72)+'TER\nEND\n'
    result = screen(text)
    assert result['cycles_considered'] == 7
    assert result['piercings'] == 0
    assert not result['qualified'] and not result['publication_gate']
