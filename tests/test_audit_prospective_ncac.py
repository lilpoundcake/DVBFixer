"""Independent snapshot reader retains real import blockers without scoring them."""
import gzip
import json

import pytest
from scripts.audit_prospective_ncac import canonical, snapshot_check
from scripts.research_rosetta_joint import sha
from scripts.research_rosetta_repack import observed_atoms
from scripts.rosetta_fixed_contract import import_audit


def fixture(tmp_path,blocked=False):
    plan=tmp_path/'plan'
    run=tmp_path/'run'
    plan.mkdir()
    folder=run/'case'/'own-0'
    folder.mkdir(parents=True)
    pdb=b'ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C  \nEND\n'
    (plan/'input.pdb').write_bytes(pdb)
    expected=observed_atoms(pdb)
    physical=list(expected)
    if blocked:
        physical.append(dict(expected[0],atom='CB',xyz=[1.,0.,0.]))
    s=dict(heavy_atoms=physical,virtual_atoms=[],residues=[dict(chain='A',resid=1,icode=' ',name='ALA')],
           hydrogens=[],external=[],target_internal=[])
    raw=gzip.compress(json.dumps(s).encode(),mtime=0)
    (folder/'pose-snapshot.json.gz').write_bytes(raw)
    row=dict(id='own-0',pdb='input.pdb',snapshot_sha256=sha(raw),
        import_audit=canonical(import_audit(expected,physical,[],[],[],[],[])),
        score_status='import_contract_blocker' if blocked else 'scored',total_score=None if blocked else 2.,
        score_terms={} if blocked else {'term':dict(raw=1.,weight=2.)},virtual_coordinates_unchanged_by_scoring=True)
    (folder/'score.json').write_text(json.dumps(row))
    case=dict(id='case',expected_external=[],target_internal=[],target_key=['A',1,' ','ALA'])
    return plan,run,case,row,folder


@pytest.mark.parametrize('blocked',[False,True])
def test_fixed_heavy_snapshot_accepts_clean_or_preserves_real_blocker(tmp_path,blocked):
    plan,run,case,row,_=fixture(tmp_path,blocked)
    result=snapshot_check(plan,run,case,row)
    assert result==dict(extras=[],virtual=[],context_h=[])


def test_blocked_import_cannot_publish_a_usable_score(tmp_path):
    plan,run,case,row,folder=fixture(tmp_path,True)
    row['total_score']=5.
    (folder/'score.json').write_text(json.dumps(row))
    with pytest.raises(ValueError,match='blocked score'):
        snapshot_check(plan,run,case,row)


def test_score_term_mismatch_is_not_a_scientific_pass(tmp_path):
    plan,run,case,row,folder=fixture(tmp_path)
    row['total_score']=4.
    (folder/'score.json').write_text(json.dumps(row))
    with pytest.raises(ValueError,match='score-term'):
        snapshot_check(plan,run,case,row)
