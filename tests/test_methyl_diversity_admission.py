"""Shared source admission retains heavy-only/occupancy/metadata refusals."""
import pytest
from scripts.freeze_methyl_diversity import admission


def source():
    atoms = [('N', 'N', (1., 0., 0.)), ('CA', 'C', (0., 0., 0.)), ('C', 'C', (0., 1., 0.)),
             ('O', 'O', (0., 1., 1.)), ('CB', 'C', (0., 0., 1.))]
    lines = []
    for serial, (name, element, (x, y, z)) in enumerate(atoms, 1):
        lines.append(f'ATOM  {serial:5d} {name:^4s} ALA A   1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2s}\n')
    return (''.join(lines)+'TER\nEND\n').encode()


def test_heavy_source_admission_is_a_read_only_scope_check():
    raw = source()
    original = raw[:]
    result = admission(raw)
    assert result['atoms'] == 5 and result['residues'] == 1 and result['heavy_only']
    assert raw == original


@pytest.mark.parametrize('fault', ['hydrogen', 'occupancy', 'altloc', 'dangling_conect', 'duplicate_serial'])
def test_source_refusals_are_not_cleaned_up_to_enter_the_cohort(fault):
    lines = source().decode().splitlines(keepends=True)
    if fault == 'hydrogen':
        line = lines[0]
        lines[0] = line[:12]+' H  '+line[16:76]+' H'+line[78:]
    elif fault == 'occupancy':
        lines[0] = lines[0][:54]+'  0.50'+lines[0][60:]
    elif fault == 'altloc':
        lines[0] = lines[0][:16]+'A'+lines[0][17:]
    elif fault == 'dangling_conect':
        lines.insert(-1, 'CONECT    1  999\n')
    else:
        lines[1] = lines[1][:6]+'    1'+lines[1][11:]
    raw = ''.join(lines).encode()
    with pytest.raises(ValueError):
        admission(raw)


def test_chain_case_is_not_normalized_by_admission():
    lines = source().decode().splitlines(keepends=True)
    for i in range(5):
        lines[i] = lines[i][:21]+'a'+lines[i][22:]
    assert admission(''.join(lines).encode())['atoms'] == 5
