"""Read-only per-call diagnostic on four fixed-sequence ASP73 copies."""
import json

import numpy as np


def install(module, features, example, root):
    names = ['N', 'CA', 'C', 'O', 'CB', 'CG', 'OD1', 'OD2']
    keys = list(zip(example['src_component'].astype(str), example['gt_atom_name'].astype(str), strict=True))
    ids = np.array([[keys.index((c+'73', n)) for n in names] for c in 'ABCD'])
    idx = ids.ravel().tolist()
    token = features['atom_to_token_map'][ids[:, 1].tolist()].detach().cpu().tolist()
    details = {'atom_order': names, 'chain_order': 'ABCD', 'features': {}}
    for k,v in features.items():
        if getattr(v,'ndim',0) and v.shape[0] in (5592,716) and (v.ndim==1 or v.shape[1] not in (716,5592)):
            arr=v[idx if v.shape[0]==5592 else token].detach().cpu().numpy().reshape((4,-1)+tuple(v.shape[1:]))
            details['features'][k] = arr.tolist()
    (root/'asp-input.json').write_text(json.dumps(details))
    count=0
    def observe(mod,args,kw,out):
        nonlocal count
        count+=1
        row={'call':count,'sigma':float(kw['t'].item())}
        for label, tensor in [('noisy',kw['X_noisy_L']),('denoised',out['X_L'])]:
            x=tensor[0,idx].detach().cpu().numpy().reshape(4,8,3)
            row[label]={label:np.linalg.norm(x[:,i]-x[:,j],axis=1).tolist() for label,i,j in
                        [('CA_CB',1,4),('CB_CG',4,5),('CG_OD1',5,6),('CG_OD2',5,7)]}
        logits=out['sequence_logits_I'][0,token].detach().cpu().numpy()
        row['sequence_head_argmax']=logits.argmax(-1).tolist()
        with (root/'asp-trace.jsonl').open('a') as f:f.write(json.dumps(row,allow_nan=False)+'\n')
    return module.register_forward_hook(observe,with_kwargs=True)
