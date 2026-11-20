"""Private RFD3-103 backbone/observed-rotamer test on protein-only 1NON."""
import itertools
import math
from collections import defaultdict
from dataclasses import asdict
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.spatial import cKDTree
from scripts.audit_anchor_calibration_grid import diagnose
from scripts.audit_retained_full_context_controls import atoms
from scripts.mac_reconstruction_comparison import digest, save

from dvbfixer.ffutils.geometry import assert_all_l
from dvbfixer.model.diffusion.contract import AtomIdentity, DiffusionRequest
from dvbfixer.model.diffusion.geometry import weighted_kabsch
from dvbfixer.model.diffusion.runner_refinement import rewrite_generated_coordinates
from dvbfixer.model.diffusion.staged_runtime import assessment_adapter


def run():
    from openmm.app import PDBFile
    b = Path('tmp/1non-backends/attempt01')
    src = b/'protein-only.pdb'
    raw = Path('tmp/1NON-inspection/08_RFD3_seed103_source_records.pdb')
    out = b/'hybrid-backbone-103'
    out.mkdir(exist_ok=False)
    source, model = atoms(src), atoms(raw)
    bb = {'N','CA','C','O'}
    generated = set(model)-set(source)
    mobile = sorted(k for k in generated if k[3] in bb)
    keys = mobile+list(source)
    index = {k:i for i,k in enumerate(keys)}
    ref = np.array([model[k]['xyz'] for k in keys]); m=len(mobile)
    assert (m,len(generated),len(source))==(288,556,5040)
    def aid(c,i,n):
        return c,str(i),'',n
    # Measured median trans peptide geometry from intact peptide units of this source.
    units = [[(0,'N'),(0,'CA'),(0,'C')],[(0,'CA'),(0,'C'),(0,'O'),(1,'N'),(1,'CA')]]
    constraints={}
    for unit in units:
        ds=defaultdict(list)
        for c in 'ABCD':
            for i in range(1,179):
                ks=[aid(c,i+j,n) for j,n in unit]
                if all(k in source for k in ks):
                    for l,h in itertools.combinations(range(len(ks)),2):
                        ds[l,h].append(math.dist(source[ks[l]]['xyz'],source[ks[h]]['xyz']))
        for c in 'ABCD':
            for i in range(72,91):
                ids=[index[aid(c,i+j,n)] for j,n in unit]
                for (l,h),v in ds.items():
                    if min(ids[l],ids[h])<m:
                        constraints[tuple(sorted([ids[l],ids[h]]))]=np.median(v)
    pairs=cKDTree(ref).query_pairs(9,output_type='ndarray'); near=[]
    for i,j in pairs:
        if min(i,j)>=m:continue
        a,z=keys[i],keys[j]
        if a[0]==z[0] and abs(int(a[1])-int(z[1]))<=1:continue
        near.append((i,j))
    l,h=np.array(list(constraints)).T; target=np.array(list(constraints.values())); u,v=np.array(near).T
    def objective(q):
        x=ref.copy(); x[:m]=q.reshape(m,3); grad=np.zeros_like(x);delta=x[:m]-ref[:m]; E=.5*np.sum(delta**2); grad[:m]=delta
        for i,j,t,w,wall in [(l,h,target,1600,False),(u,v,3.45,35,True)]:
            dx=x[i]-x[j]; dist=np.linalg.norm(dx,axis=1); residual=dist-t
            if wall:residual=np.minimum(residual,0)
            E+=.5*w*(residual@residual);f=(w*residual/np.maximum(dist,1.e-10))[:,None]*dx
            np.add.at(grad,i,f);np.add.at(grad,j,-f)
        return E,grad[:m].ravel()
    opt=minimize(objective,ref[:m].ravel(),jac=True,method='L-BFGS-B',bounds=[(a-2.5,a+2.5) for a in ref[:m].ravel()],options={'maxiter':1500,'ftol':1e-12})
    placed=dict(zip(mobile,opt.x.reshape(m,3),strict=True));background={k:np.array(a['xyz']) for k,a in source.items()};background.update(placed)
    groups=defaultdict(list)
    for k in source:groups[k[:3]].append(k)
    library=defaultdict(list)
    for g,ks in groups.items():library[source[ks[0]]['resname']].append(g)
    catalog={};selected={};donors={}
    for g in sorted({k[:3] for k in generated}):
        side=sorted(k for k in generated if k[:3]==g and k[3] not in bb)
        if not side:continue
        names=['N','CA','C']; tgt=np.array([placed[g+(n,)] for n in names]); opts=[]
        for d in library[model[side[0]]['resname']]:
            if not all(d+(k[3],) in source for k in side):continue
            base=np.array([background[d+(n,)] for n in names]); tf=weighted_kabsch(base-base[1],tgt-tgt[1]); xyz=tf.apply(np.array([background[d+(k[3],)] for k in side])-base[1])+tgt[1];opts.append((d,xyz))
        catalog[g]=(side,opts)
    for sweep in range(4):
        for g,(side,opts) in catalog.items():
            cloud=[p for k,p in background.items() if k[:3]!=g];tree=cKDTree(cloud)
            xyz=np.stack([x for _,x in opts]);ds=tree.query(xyz.reshape(-1,3),k=4)[0].reshape(len(opts),len(side),4);scores=np.minimum(ds-3.45,0)**2;best=int(np.argmin(scores.sum(axis=(1,2))));d,x=opts[best];background.update(zip(side,x));selected.update(zip(side,x));donors[str(g)]=d
    placed.update(selected);assert set(placed)==generated
    text=rewrite_generated_coordinates(raw.read_text(),{AtomIdentity(*k):v for k,v in placed.items()});p=out/'packed.pdb';p.write_text(text);parsed=atoms(p);assert set(parsed)==set(model) and all(parsed[k]['line']==a['line'] for k,a in source.items());pdb=PDBFile(str(p));assert_all_l(pdb.topology,pdb.positions)
    findings=diagnose(p,out/'full-cli.json');q=DiffusionRequest.from_json((b/'context-request.json').read_text());ass=assessment_adapter(q,src.read_text())(text);save(out/'assessment.json',asdict(ass));save(out/'donors.json',donors);save(out/'receipt.json',{'source':digest(src),'raw':digest(raw),'packed':digest(p),'backbone_count':m,'rebuilt_side_atoms':len(selected),'backbone_rms_A':float(np.sqrt(np.mean(np.sum((opt.x.reshape(m,3)-ref[:m])**2,axis=1)))),'solver':str(opt.message),'iterations':int(opt.nit),'geometry_passed':ass.geometry_passed,'clashes':sum(f['category']=='clash' and f['severity']=='ERROR' for f in findings),'junctions':{c:[math.dist(parsed[aid(c,l,'C')]['xyz'],parsed[aid(c,h,'N')]['xyz']) for l,h in ((72,73),(90,91))] for c in 'ABCD'}})


if __name__=='__main__':run()
