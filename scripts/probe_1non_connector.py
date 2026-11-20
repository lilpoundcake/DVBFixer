"""Inspect first real cc94 denoiser call, then stop before any ODE update."""
import argparse
import json
import runpy
import sys
from pathlib import Path
from unittest.mock import patch

import numpy as np
import torch
from protpardelle.common import residue_constants as rc
from protpardelle.core import models
from scripts.run_1non_diffusion_full import native_arguments


class ProbeComplete(Exception):
    pass


def run(root, out):
    plan = json.loads((root/'plan.json').read_text())
    source = np.zeros((716,37,3))
    known = np.zeros((716,37), bool)
    for line in Path(plan['source']).read_text().splitlines():
        if line.startswith('ATOM  '):
            i = 'ABCD'.index(line[21])*179+int(line[22:26])-1
            j = rc.atom_order[line[12:16].strip()]
            known[i,j] = True
            source[i,j] = [float(line[k:k+8]) for k in (30,38,46)]
    assert known.sum() == 5040
    saved = root/'native-sampling/inpainting'
    with np.load(saved/'reference.npz') as data:
        np.testing.assert_array_equal(known,data['known_mask'][0])
        gen = data['generated_mask'][0]
    out.mkdir(parents=True,exist_ok=False)
    load = models.load_model

    def instrument(*args, **kwargs):
        model = load(*args, **kwargs)
        original = model.struct_model.forward

        def forward(*args, **kw):
            crop = kw['struct_crop_cond']
            assert args[2].shape == (1,716) and bool((args[2]==1).all())
            capture = []
            handle = model.struct_model.net.register_forward_pre_hook(lambda mod, inputs: capture.append(inputs[0].detach().cpu().numpy()) or None)
            base = original(*args,**kw)
            handle.remove()
            array = crop[0].cpu().numpy()
            np.testing.assert_array_equal(np.any(array!=0,axis=-1),known)
            np.testing.assert_array_equal(capture[0][0,...,-3:], (crop/model.sigma_data)[0].cpu().numpy())
            x,y = array[known].astype(float),source[known]
            u,_,v = np.linalg.svd((x-x.mean(0)).T@(y-y.mean(0)))
            rot = u@np.diag([1.,1.,np.linalg.det(u@v)])@v
            fitted = (x-x.mean(0))@rot+y.mean(0)
            error = np.linalg.norm(fitted-y,axis=1)
            assert error.max()<1.e-4
            compare = {}
            with np.load(saved/'projection-clean-prediction-0000.npz') as data:
                archived_error = float(np.max(np.abs(base[0].cpu().numpy()[gen]-data['generated_before'])))
            for label in ('repeat','hide_B','hide_all'):
                control = crop.clone()
                if label=='hide_B':
                    control[:,179:358]=0
                if label=='hide_all':
                    control.zero_()
                delta=(original(*args,**{**kw,'struct_crop_cond':control})-base)[0].cpu().numpy()
                compare[label]={c: float(np.sqrt(np.mean(np.sum(delta[i*179:(i+1)*179][gen[i*179:(i+1)*179]]**2,axis=-1)))) for i,c in enumerate('ABCD')}
            np.savez_compressed(out/'network-input.npz',embedding=capture[0],known=known,generated=gen,source=source,crop=array)
            report={'probe':'first denoising call, same noisy coordinates in all four forwards',
                    'observed_per_chain':known.reshape(4,179,37).sum(axis=(1,2)).tolist(),
                    'seq_mask_count':716,'sigma':float(args[1].item()),'sigma_data':float(model.sigma_data),
                    'embedding_shape':capture[0].shape,'crop_channels_match_normalized_coordinates':True,
                    'proper_rigid_fit_max_angstrom':float(error.max()), 'rms_angstrom':float(np.sqrt(np.mean(error**2))),
                    'archived_step0_generated_max_abs_diff':archived_error,'generated_rms_change_angstrom':compare}
            (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
            raise ProbeComplete('intentional stop after first-step connector diagnostic')
        model.struct_model.forward=forward
        return model
    argv=native_arguments(plan,out/'probe-input')
    try:
        with patch.object(models,'load_model',instrument), patch.object(sys,'argv',argv):
            runpy.run_path(argv[0],run_name='__main__')
    except ProbeComplete:
        return
    raise RuntimeError('probe failed to intercept')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path)
    p.add_argument('output',type=Path)
    args=p.parse_args()
    run(args.root.resolve(),args.output.resolve())
