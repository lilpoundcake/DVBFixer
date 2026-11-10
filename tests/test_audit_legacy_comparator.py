"""Real observation motion remains visible; diagnostic remasking never restores it."""
from scripts.audit_legacy_comparator import atoms, diagnostic_projection, motion

from dvbfixer.domain.structure_identity import ComponentInstanceRef


def pdb(name,x,serial=1,element='C'):
    return f'ATOM  {serial:5d} {name:^4s} VAL A   1    {x:8.3f}{0.:8.3f}{0.:8.3f}  1.00  0.00          {element:>2s}\n'.encode()


def test_moved_original_is_measured_in_original_frame_not_aligned_or_restored():
    result=motion(pdb('CA',0.),pdb('CA',2.))
    assert result['maximum_displacement_angstrom']==2.
    assert result['rms_displacement_angstrom']==2.
    assert result['moved_at_pdb_precision']==1
    assert not result['aligned'] and not result['observations_restored']


def test_missing_observation_is_retained_not_silently_excluded():
    result=motion(pdb('CA',0.)+pdb('N',1.,2,'N'),pdb('CA',0.))
    assert result['original_heavy_count']==2 and result['matched_heavy_count']==1
    assert len(result['missing_identities'])==1


def test_remasking_uses_real_heavy_positions_and_retains_raw_input():
    target=ComponentInstanceRef(1,'A','1','')
    raw=pdb('CA',2.)+pdb('CB',3.,2)+pdb('H',4.,3,'H')
    projected=diagnostic_projection(raw,target,('CB',))
    assert len(atoms(projected))==1
    assert next(iter(atoms(projected).values())).position==(2.,0.,0.)
    assert b' H  ' in raw or len(raw.splitlines())==3
    assert len(raw.splitlines())==3
