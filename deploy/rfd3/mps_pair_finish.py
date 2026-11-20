"""Inference-only full pair initializer: stream additions and projection by tiles."""
import ast

import torch
from torch.nn import functional as F


def finish(module, pair, single, tokens, idx, rows=64):
    if module.training or torch.is_grad_enabled() or pair.ndim != 3 or rows < 1:
        raise ValueError('requires owned unbatched inference pair tensor')
    l, r = module.process_single_l(single), module.process_single_m(single)
    z = module.process_z(tokens)
    for i in range(0,len(idx),rows):
        sl = slice(i,i+rows)
        tile = pair[sl]
        tile.add_(l[sl,None,:]+r[None,:,:])
        tile.add_(z[idx[sl,None],idx[None,:]])
        tile.add_(module.pair_mlp(tile))
    # Tile the uncontracted key axis, retaining the native two reductions.
    count = int(idx.max().item())+1
    one = F.one_hot(idx.long(),count).to(pair.dtype)
    temp = pair.new_empty(count,len(idx),pair.shape[-1])
    for i in range(0,len(idx),rows):
        sl = slice(i,i+rows)
        projected = module.process_pll(pair[:,sl].contiguous())
        temp[:,sl] = torch.einsum('ia,acd->icd',one.T,projected)
    total = torch.einsum('cj,icd->ijd',one,temp)
    n = one.sum(0)
    return total/torch.outer(n,n).clamp(min=1).unsqueeze(-1)


def install(code):
    tree = ast.parse(code)
    matches = []
    for node in ast.walk(tree):
        body = getattr(node,'orelse',None)
        if not isinstance(body,list):
            continue
        for i,stmt in enumerate(body):
            if isinstance(stmt,ast.Assign) and ast.unparse(stmt.value).startswith('P_LL + (self.process_single_l'):
                originals = [ast.unparse(x) for x in body[i:i+5]]
                if ('self.process_z(' not in originals[1] or 'self.pair_mlp(' not in originals[2]
                        or 'P_LL.contiguous()' not in originals[3] or 'pairwise_mean_pool(' not in originals[4]):
                    raise ValueError('unexpected pair initializer')
                body[i:i+5] = ast.parse('pooled_atom_level_features = _finish(self,P_LL,C_L,Z_init_II,tok_idx)').body
                matches.append(i)
                break
    if len(matches)!=1:
        raise ValueError('could not identify unique pair initializer block')
    return ast.unparse(tree)+'\n'


def gather_pair(pair, indices):
    query = torch.arange(indices.shape[1], device=indices.device)[None,:,None]
    if pair.ndim == 3:
        return pair[query,indices]
    batch = torch.arange(indices.shape[0],device=indices.device)[:,None,None]
    return pair[batch,query,indices]


def sparse_rows(Q,K,V,B,indices,H,gather_bias=False,G=None,full=False,rows=64):
    if full or gather_bias or torch.is_grad_enabled():
        raise ValueError('requires native pre-gathered sparse inference')
    batch,length,c=Q.shape
    bidx=torch.arange(batch,device=Q.device)[:,None,None]
    result=torch.empty_like(Q)
    for i in range(0,length,rows):
        sl=slice(i,i+rows)
        idx=indices[:,sl]
        q=Q[:,sl].unflatten(-1,(H,c//H)).transpose(1,2)
        k=K[bidx,idx].unflatten(-1,(H,c//H)).permute(0,3,1,2,4)
        logits=torch.einsum('...ld,...lkd->...lk',q,k)/(c//H)**.5+B[:,sl].permute(0,3,1,2)
        v=V[bidx,idx].unflatten(-1,(H,c//H)).permute(0,3,1,2,4)
        out=torch.einsum('...ij,...ijc->...ic',logits.softmax(-1),v).transpose(1,2).reshape(batch,idx.shape[1],c)
        result[:,sl]=out if G is None else out*G[:,sl]
    return result


def install_local_bias():
    import inspect
    import textwrap

    from rfd3.model.layers.attention import LocalAttentionPairBias
    native = LocalAttentionPairBias.forward
    code = textwrap.dedent(inspect.getsource(native))
    marker = 'P_LL_sparse = None\n        use_sparse_pll = False'
    if code.count(marker)!=1:
        raise ValueError('unknown local attention branch')
    code=code.replace(marker,'use_sparse_pll = not full and not self.training and not torch.is_grad_enabled() and P_LL.device.type == "mps"\n        P_LL_sparse = _gather(P_LL, indices) if use_sparse_pll else None')
    ns=dict(native.__globals__,_gather=gather_pair)
    import os
    if os.environ.get('DVBFIXER_PRIVATE_SPARSE_ROWS')=='1':
        code=code.replace('attn_out = sparse_pairbias_attention(\n                    Q=q,','attn_out = _rows(\n                    Q=q,',1)
        if 'attn_out = _rows(' not in code:
            raise ValueError('unknown sparse attention branch')
        ns['_rows']=sparse_rows
    exec(compile(code,'<private-local-bias>','exec'),ns)
    LocalAttentionPairBias.forward=ns['forward']
