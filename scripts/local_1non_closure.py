"""Private bounded torsional window closure on saved RFD3-derived model 12."""
import math
from dataclasses import asdict
from pathlib import Path

import numpy as np
from numba import njit
from scipy.optimize import least_squares
from scipy.spatial import cKDTree
from openmm.app import PDBFile
from scripts.audit_anchor_calibration_grid import diagnose
from scripts.audit_retained_full_context_controls import atoms
from scripts.mac_reconstruction_comparison import save

from dvbfixer.ffutils.geometry import assert_all_l
from dvbfixer.model.diffusion.contract import AtomIdentity, DiffusionRequest
from dvbfixer.model.diffusion.geometry import weighted_kabsch
from dvbfixer.model.diffusion.runner_refinement import rewrite_generated_coordinates
from dvbfixer.model.diffusion.staged_runtime import assessment_adapter


@njit
def atom(a,b,c,length,angle,t):
    u=c-b;u=u/np.linalg.norm(u);v=np.cross(b-a,u);v=v/np.linalg.norm(v);w=np.cross(v,u);ang=angle*np.pi/180
    return c+length*(-np.cos(ang)*u+np.sin(ang)*(w*np.cos(t)+v*np.sin(t)))


@njit
def build(left,q):
    n=(len(q)-1)//2;arr=np.zeros((n*5+2,3));ca,c,o=left
    N=atom(o,ca,c,1.33,116.2,np.pi);CA=atom(ca,c,N,1.46,121.7,np.pi+q[-1]);prev=c
    for i in range(n):
        C=atom(prev,N,CA,1.525,111.2,q[2*i]);NN=atom(N,CA,C,1.33,116.2,q[2*i+1]);O=atom(N,CA,C,1.23,120.8,q[2*i+1]+np.pi)
        CB=atom(prev,N,CA,1.53,110.5,q[2*i]-2.12)
        arr[5*i]=N;arr[5*i+1]=CA;arr[5*i+2]=C;arr[5*i+3]=O;arr[5*i+4]=CB
        CA=atom(CA,C,NN,1.46,121.7,np.pi);prev=C;N=NN
    arr[-2]=N;arr[-1]=CA;return arr


def run(expanded=False, packed_selector=None):
    base=Path('tmp/1non-backends/attempt01');insp=Path('tmp/1NON-inspection');raw=insp/'12_RFD3_junction_control.pdb';src=base/'protein-only.pdb';out=base/('A73-packed-selection' if packed_selector else 'torsion-windows-expanded' if expanded else 'torsion-windows');out.mkdir(exist_ok=False)
    model=atoms(raw);source=atoms(src);xyz={k:np.array(a['xyz']) for k,a in model.items()};bb=['N','CA','C','O','CB'];rng=np.random.default_rng(103);updates={};history=[]
    for chain,lo,hi in ([('A',73,77),('B',73,77),('B',86,90)] if expanded else [('A',73,75),('B',73,75),('B',88,90)]):
        n=hi-lo+1;key=lambda i,name:(chain,str(i),'',name)
        left=np.array([xyz[key(lo-1,name)] for name in ['CA','C','O']]);right=np.array([xyz[key(hi+1,name)] for name in ['N','CA']]);ref=np.array([xyz[key(i,'CA')] for i in range(lo,hi+1)])
        envkeys=[k for k,p in xyz.items() if not(k[0]==chain and lo<=int(k[1])<=hi) and np.min(np.linalg.norm(ref-p,axis=1))<11];env=np.array([xyz[k] for k in envkeys]);inc=[]
        for i in range(lo,hi+1):
            for name in bb:
                inc.append([not(k[0]==chain and abs(int(k[1])-i)==1 and k[3] in bb[:4] and name in bb[:4]) and key(i,name) in model for k in envkeys])
        inc=np.array(inc);tree=cKDTree(env);J,K=np.where(inc & (np.linalg.norm(np.repeat(ref,5,axis=0)[:,None]-env,axis=2)<8.5));
        def resid(q):
            x=build(left,q);dist=np.linalg.norm(x[J]-env[K],axis=1);return np.r_[(140 if expanded else 30)*(x[-2:]-right).ravel(),5*np.minimum(dist-3.35,0),.35*(x[1:-2:5]-ref).ravel(),(80 if expanded else 4)*q[-1]]
        opts=[]
        for rep in range(32):
            q=np.r_[np.deg2rad(np.array([[-65,-35],[-95,90],[-125,125],[-80,145]])[rng.integers(4,size=n)].ravel())+rng.normal(0,.35,2*n),0.];res=least_squares(resid,q,max_nfev=100,ftol=1e-6);x=build(left,res.x);opts.append((float(res.cost),res.x.tolist(),x))
        if packed_selector:
            return packed_selector(raw,src,out,opts,chain,lo,hi)
        cost,angles,x=min(opts,key=lambda opt:opt[0]);history.append({'chain':chain,'window':[lo,hi],'cost':cost,'terminal_max_A':float(np.linalg.norm(x[-2:]-right,axis=1).max()),'q':angles});
        for i in range(lo,hi+1):
            for j,name in enumerate(bb[:4]):updates[key(i,name)]=x[5*(i-lo)+j];xyz[key(i,name)]=updates[key(i,name)]
        # Pack side chains from observed same-residue conformers in source, using corrected backbone.
        for _ in range(3):
            for i in range(lo,hi+1):
                side=[k for k in xyz if k[:3]==key(i,'')[:3] and k[3] not in bb[:4]]
                if not side:continue
                names=['N','CA','C'];b0=np.array([xyz[key(i,name)] for name in names]);tree=cKDTree([p for k,p in xyz.items() if k[:3]!=key(i,'')[:3]]);rot=[]
                for donor in source:
                    if donor[3]!='CA' or source[donor]['resname']!=model[key(i,'CA')]['resname']:continue
                    ks=[donor[:3]+(k[3],) for k in side];
                    if not all(k in source for k in ks):continue
                    srcbb=np.array([source[donor[:3]+(name,)]['xyz'] for name in names]);t=weighted_kabsch(srcbb,b0);sc=t.apply(np.array([source[k]['xyz'] for k in ks]));ds=tree.query(sc,k=5)[0];rot.append((np.square(np.minimum(ds-3.35,0)).sum(),sc))
                chosen=min(rot,key=lambda z:z[0])[1];xyz.update(zip(side,chosen));updates.update(zip(side,chosen))
    text=rewrite_generated_coordinates(raw.read_text(),{AtomIdentity(*k):v for k,v in updates.items()});path=out/'candidate.pdb';path.write_text(text);parsed=atoms(path);assert all(v['line']==parsed[k]['line'] for k,v in model.items() if k not in updates);assert not (updates.keys()&source.keys());p=PDBFile(str(path));assert_all_l(p.topology,p.positions);q=DiffusionRequest.from_json((base/'context-request.json').read_text());a=assessment_adapter(q,src.read_text())(text);save(out/'assessment.json',asdict(a));diagnose(path,out/'full-cli.json');save(out/'windows.json',history)


if __name__=='__main__':
    import sys
    run('--expanded' in sys.argv)
