"""New-system enrollment and graph transport do not select by outcome."""
import json
from pathlib import Path

import pytest
from scripts import prospective_ncac_cohort as module
from scripts.compare_branched_completion import measure

from dvbfixer.domain.structure_identity import ComponentInstanceRef
from dvbfixer.protein_completion.branched_source import plan_branched_source
from dvbfixer.protein_completion.service import complete_sidechain


def test_whole_mask_is_fixed_and_unsupported_chemistry_refuses():
    assert module.whole_mask('VAL')==['CB','CG1','CG2']
    assert module.whole_mask('THR')==['CB','OG1','CG2']
    with pytest.raises(ValueError):
        module.whole_mask('ILE')


def test_new_systems_are_predeclared_not_success_selected():
    assert list(module.SYSTEMS)==['1PGB','1BPI','1CSP']


def test_score_transport_retains_all_proposals_but_only_scores_feasible(tmp_path):
    reference,_=module.project_first(Path('tests/fixtures/protein_context/independent_followup/2CI2.pdb').read_bytes())
    source=b''.join(line for line in reference.splitlines(keepends=True) if not (
        line.startswith(b'ATOM  ') and line[21:22]==b'I' and line[22:26].strip()==b'22'
        and line[12:16].strip() in {b'CB',b'OG1',b'CG2'}))
    request=plan_branched_source(source,ComponentInstanceRef(1,'I','22')).request
    own=complete_sidechain(request)
    rows=[]
    for i,attempt in enumerate(own.attempts):
        audit,pdb=measure(source,request,attempt.candidate.additions,anchored=False)
        rows.append((i,audit,pdb))
    artifacts={}
    case=module.score_case('synthetic','baseline',source,request,rows,tmp_path,artifacts)
    assert len(list((tmp_path/'synthetic-baseline').glob('own-*.pdb')))==len(rows)
    assert [c['index'] for c in case['cells']]==[i for i,audit,_ in rows if audit.accepted]
    assert (tmp_path/'synthetic-baseline/source.pdb').read_bytes()==source
    assert not any('reference' in k or 'native' in k for k in case)
    assert len(case['generated'])==3
    assert json.loads(json.dumps(case))['target_key']==['I',22,' ','THR']
