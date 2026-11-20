"""Rank the same 32 A73–75 backbones after packing, against model 12."""
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree
from scripts.audit_anchor_calibration_grid import diagnose
from scripts.audit_retained_full_context_controls import atoms
from scripts.local_1non_closure import run
from scripts.mac_reconstruction_comparison import digest, save

from dvbfixer.model.diffusion.contract import AtomIdentity, DiffusionRequest
from dvbfixer.model.diffusion.geometry import weighted_kabsch
from dvbfixer.model.diffusion.runner_refinement import rewrite_generated_coordinates
from dvbfixer.model.diffusion.staged_runtime import assessment_adapter


def select(raw,src,out,opts,chain,lo,hi):
    assert (chain,lo,hi)==('A',73,75)
    model,source=atoms(raw),atoms(src)
    bb=('N','CA','C','O'); window={k for k in model if k[0]=='A' and lo<=int(k[1])<=hi}
    assert window.isdisjoint(source) and all(model[k]['line']==v['line'] for k,v in source.items())
    req=DiffusionRequest.from_json((src.parent/'context-request.json').read_text()); assess=assessment_adapter(req,src.read_text())
    library={}
    for i in range(lo,hi+1):
        rg=('A',str(i),'');side=sorted(k for k in window if k[:3]==rg and k[3] not in bb)
        donors=[]
        for d in source:
            if d[3]!='CA' or source[d]['resname']!=model[rg+('CA',)]['resname']:continue
            if all(d[:3]+(k[3],) in source for k in side):
                donors.append((np.array([source[d[:3]+(n,)]['xyz'] for n in bb[:3]]),np.array([source[d[:3]+(k[3],)]['xyz'] for k in side])))
        library[i]=side,donors
    def score(path,directory):
        d=diagnose(path,directory/'full-cli.json');a=assess(path.read_text());issues=Counter((x.category,tuple((t.chain,t.residue_number,t.insertion_code) for t in x.residues)) for x in a.issues if x.residues and x.category!='clash');issues.update((f['category'],f.get('chain'),f.get('resid'),f.get('atom')) for f in d if f['category'] in ('ca_planarity','chirality','bond_length','bond_angle','non_planar_amide','coincident_atoms'))
        return sum(f['category']=='clash' and f['severity']=='ERROR' for f in d),issues,a.admission_failures
    n0,i0,ad=score(raw,out); assert n0==7 and not ad
    rows=[];base=raw.read_text()
    for j,(cost,angles,x) in enumerate(opts):
        folder=out/f'{j:02d}';folder.mkdir();xyz={k:np.array(v['xyz']) for k,v in model.items()}
        for i in library:
            for t,n in enumerate(bb):xyz['A',str(i),'',n]=x[5*(i-lo)+t]
        for _ in range(3):
            for i,(side,donors) in library.items():
                rg=('A',str(i),'');target=np.array([xyz[rg+(n,)] for n in bb[:3]]);tree=cKDTree([v for k,v in xyz.items() if k[:3]!=rg]);positions=np.stack([weighted_kabsch(b,target).apply(s) for b,s in donors]);ds=tree.query(positions.reshape(-1,3),k=5)[0].reshape(len(donors),len(side),5);best=np.argmin((np.minimum(ds-3.35,0)**2).sum(axis=(1,2)));xyz.update(zip(side,positions[best],strict=True))
        text=rewrite_generated_coordinates(base,{AtomIdentity(*k):xyz[k] for k in window});path=folder/'candidate.pdb';path.write_text(text);p=atoms(path);assert set(p)==set(model) and all(v['line']==p[k]['line'] for k,v in model.items() if k not in window)
        n,issues,ad=score(path,folder);new=issues-i0;rows.append(dict(index=j,backbone_cost=cost,clashes=n,new_findings={str(k):v for k,v in new.items()},eligible=n<n0 and not new and not ad))
    valid=[r for r in rows if r['eligible']];winner=min(valid,key=lambda r:(r['clashes'],r['backbone_cost'])) if valid else None
    save(out/'comparison.json',dict(parent_sha256=digest(raw),baseline=n0,window='A73–75',rows=rows,winner=winner));print('clash distribution:',Counter(r['clashes'] for r in rows),'eligible:',len(valid),'winner:',winner)
    if winner:
        import shutil
        shutil.copyfile(out/f"{winner['index']:02d}"/'candidate.pdb',Path('tmp/1NON-inspection/15_A73_packed_selection.pdb'))


if __name__=='__main__':run(packed_selector=select)
