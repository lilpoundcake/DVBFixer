"""One paired last-state denoiser call with identical initializer outputs."""
import json
import os
import runpy
from pathlib import Path

import numpy as np
import torch
from rfd3.model.layers.attention import LocalAttentionPairBias
from rfd3.model.layers.layer_utils import Transition

original = {cls: cls.forward for cls in (LocalAttentionPairBias, Transition)}
root = Path(os.environ['DVBFIXER_RFD3_TRANSPORT_OUTPUT']).parent


def probe(owner, module, kw):
    parent = Path(os.environ['PAIR_PARENT'])/'transport'
    with np.load(parent/'native-transport.npz') as saved:
        x = torch.as_tensor(saved['last_noisy'].copy(),device='mps')
        endpoint = saved['final_X_L'].copy()
        free = ~saved['is_motif_atom_with_fixed_coord']
    sampler = owner.inference_sampler.sampler
    sched = sampler._construct_inference_noise_schedule(device=x.device,partial_t=None)
    t = sched[-2]*(1+(sampler.gamma_0 if sched[-1]>sampler.gamma_min else 0))
    args = dict(kw,X_noisy_L=x,t=t.tile(1))
    assert not owner.use_classifier_free_guidance and not sampler.allow_realignment
    cpu, mps = torch.get_rng_state(),torch.mps.get_rng_state()
    def call():
        torch.set_rng_state(cpu)
        torch.mps.set_rng_state(mps)
        return module.forward(**{**args,'f':args['f'].copy()})['X_L'].detach()
    base=call()
    reconstructed = x+sampler.step_scale*(sched[-1]-t)*((x-base)/t)
    report={'t':float(t),'blocked_reconstructed_endpoint_max_A':float(np.abs(reconstructed.cpu().numpy()-endpoint).max())}
    np.save(root/'blocked.npy',base.cpu().numpy())
    (root/'probe.json').write_text(json.dumps(report,indent=2))
    if report['blocked_reconstructed_endpoint_max_A']>1e-4:
        raise RuntimeError('failed archived endpoint replay')
    torch.mps.synchronize();torch.mps.empty_cache()
    for cls,fn in original.items(): cls.forward=fn
    for cls in (LocalAttentionPairBias,Transition): assert cls.forward is original[cls]
    if os.environ.get('PAIR_NATIVE_TRIM'):
        def trim(*unused):
            torch.mps.synchronize()
            torch.mps.empty_cache()
        for sub in module.modules():
            sub.register_forward_pre_hook(trim)
            sub.register_forward_hook(trim)
    if os.environ.get('PAIR_CPU_REFERENCE'):
        def host(v):
            if isinstance(v,torch.Tensor): return v.cpu()
            if isinstance(v,dict): return {k:host(t) for k,t in v.items()}
            if isinstance(v,list): return [host(t) for t in v]
            if isinstance(v,tuple): return tuple(host(t) for t in v)
            return v
        hostargs=host(args)
        torch.set_num_threads(2)
        torch.set_rng_state(cpu)
        native=module.cpu().forward(**hostargs)['X_L'].detach().numpy()
        report['reference_device']='cpu'
    else:
        native=call().cpu().numpy()
        report['reference_device']='mps'
    np.save(root/'native.npy',native)
    diff=np.linalg.norm(native-base.cpu().numpy(),axis=-1)[0]
    ex=np.load(parent/'native-example.npz'); ids=[(str(c),str(n)) for c,n in zip(ex['src_component'],ex['gt_atom_name'])]
    ix=[ids.index(('A73',n)) for n in ('CA','CB')]
    report.update(generated_rms_A=float(np.sqrt(np.mean(diff[free]**2))),generated_max_A=float(diff[free].max()),
        A73_CA_CB_blocked=float(torch.linalg.norm(base[0,ix[0]]-base[0,ix[1]])),A73_CA_CB_native=float(np.linalg.norm(native[0,ix[0]]-native[0,ix[1]])))
    (root/'probe.json').write_text(json.dumps(report,indent=2)+'\n')
    raise SystemExit(0)


if __name__=='__main__':
    import sys
    import types
    exported = types.ModuleType('last_step_probe')
    exported.probe = probe
    sys.modules['last_step_probe']=exported
    runpy.run_module('mps_accumulate_entry',run_name='__main__')
