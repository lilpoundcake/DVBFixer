import pytest
from scripts.export_rfd3_context_ensemble import atom_records, model_lines

from dvbfixer.model.diffusion.contract import AtomIdentity


def line(chain, number, serial):
    return f'ATOM  {serial:5d}  CA  ALA {chain}{number:4d}    {1:8.3f}{2:8.3f}{3:8.3f}{1:6.2f}{20:6.2f}           C'


def test_visualization_preserves_identity_coordinates_and_noncontiguous_breaks():
    text = '\n'.join([line('A', 436, 90), line('A', 437, 95), line('A', 489, 101), line('a', 489, 105)])
    records = atom_records(text)
    axis = tuple(records)
    result = model_lines(axis, records, model_index=2)
    assert result.count('TER\n') == 3
    readback = atom_records(result)
    assert set(readback) == set(records)
    assert {a.chain for a in readback} == {'A', 'a'}
    assert all(readback[a][30:] == records[a][30:] for a in axis)
    assert [int(l[6:11]) for l in result.splitlines() if l.startswith('ATOM  ')] == [1,2,3,4]


def test_visualization_refuses_missing_or_duplicate_atoms():
    row = line('A', 1, 1)
    with pytest.raises(ValueError, match='duplicate'):
        atom_records(row+'\n'+row)
    with pytest.raises(ValueError, match='omits'):
        model_lines((AtomIdentity('A', '1', '', 'CA'),), {}, model_index=1)
