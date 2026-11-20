"""Private four-component cc94 comparison using the frozen 1NON source and ownership."""
import argparse
import json
import runpy
import subprocess
import sys
from pathlib import Path

from scripts.mac_reconstruction_comparison import digest, save
from scripts.run_1non_diffusion_full import PROJECT, environment


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode', choices=['declare','run','worker'])
    p.add_argument('--root', required=True, type=Path)
    p.add_argument('--index', type=int, choices=range(4))
    p.add_argument('--full-own-chain',action='store_true')
    args = p.parse_args()
    root = args.root.resolve()
    if args.mode == 'declare':
        from scripts.retained_evidence import freeze_code

        from dvbfixer.model.diffusion.contract import DiffusionRequest
        from dvbfixer.model.diffusion.joint_gap_research import motif_component_admission, plan_motif_component
        base = PROJECT/'tmp/1non-fixed-scaffold-inpainting-2026-11-17/attempt01/study/plan.json'
        b = json.loads(base.read_text())
        root.mkdir(parents=True,exist_ok=False)
        source = Path(b['source']).read_bytes()
        parent = DiffusionRequest.from_json(Path(b['parent_request']).read_text())
        jobs = []
        for i in range(4):
            view = plan_motif_component(parent,source,(i,),seed=23,full_own_chain=args.full_own_chain)
            admission = motif_component_admission(view)
            assert admission['supported'] and len(view.request.generated_atoms)==139
            jobs.append(admission)
        files = [Path(path) for path in b['source_code_and_input_sha256']]+[base,Path(__file__)]
        for name in ('source','checkpoint','parent_request'):
            assert digest(Path(b[name]))==b['source_code_and_input_sha256'][b[name]]
        plan={**b,'protocol':'private-1non-local-components-v1','input':'one frozen spatial motif per gap',
              'comparison':'four separate local calls vs one full-tetramer inpainting call, seed23, no placement',
              'model_calls':4,'candidate_count':1,'jobs':jobs,'baseline_manifest':str(base),
              'source_code_and_input_sha256':{str(f.resolve()):digest(f) for f in files}}
        if args.full_own_chain:
            plan.update(full_own_chain=True,protocol='private-1non-full-own-chain-v1', input='full own chain + unchanged partner motif',comparison='same four-call local control')
        save(root/'plan.json',plan)
        freeze_code(root/'plan.json')
        return
    from scripts.run_1non_diffusion_full import load_plan
    plan = load_plan(root)
    if args.mode=='worker':
        i=args.index
        script=str(PROJECT/'deploy/protpardelle-1c/joint_gap_research.py')
        sys.argv=[script,'--parent-request',plan['parent_request'],'--source',plan['source'],'--baseline',plan['source'],
                  '--gap-indices',str(i),'--local-inpainting','--fixed-scaffold-inpainting','--seed','23',
                  '--checkpoint',plan['checkpoint'],'--refinement-python',plan['refinement_python'],
                  '--output',str(root/f'native-{i}'),'--acknowledge-private-joint']
        if plan.get('full_own_chain'):
            sys.argv.append('--full-own-chain')
        runpy.run_path(script,run_name='__main__')
    else:
        (root/'reservation').touch(exist_ok=False)
        for i in range(4):
            subprocess.run([sys.executable,'-m','scripts.run_mac_reconstruction_guarded','--python',plan['python'],
                '--module','scripts.run_1non_motif_control','--output',str(root/f'guard-{i}'),'--timeout','1800',
                '--resource-policy','pressure-v2','--','worker','--root',str(root),'--index',str(i)],
                env=environment(plan),check=True)


if __name__=='__main__':
    main()
