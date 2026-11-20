"""Single matched-budget 1NON replay: add a signed sine wall on generated L centres."""
from __future__ import annotations

import argparse
import json
import math
import shutil
import sys
from pathlib import Path

from openmm import CustomCompoundBondForce, XmlSerializer
from openmm.unit import dalton
from scripts.anchor_calibration_grid import identity, read_temporary
from scripts.audit_1non_generated_refinement import audit_prepared
from scripts.mac_reconstruction_comparison import digest, save
from scripts.refine_1non_generated_only import worker
from scripts.retained_evidence import freeze_code, verify_bindings


def shape_force():
    # Sine of signed elevation of CB above N–CA–C. Regularize normalization
    # at collapsed bonds/plane, avoid angular singularities (no asin).
    f = CustomCompoundBondForce(4,
        '5000*max(0,shape_floor-s)^2;s=(nx*cx+ny*cy+nz*cz)/sqrt((nx^2+ny^2+nz^2)*(cx^2+cy^2+cz^2)+1e-24);'
        'nx=ay*bz-az*by;ny=az*bx-ax*bz;nz=ax*by-ay*bx;'
        'ax=x1-x2;ay=y1-y2;az=z1-z2;bx=x3-x2;by=y3-y2;bz=z3-z2;cx=x4-x2;cy=y4-y2;cz=z4-z2')
    f.addGlobalParameter('shape_floor', math.sin(math.radians(20)))
    f.setName('Generated CA signed elevation wall')
    return f


def declare(parent, root):
    path = parent / 'plan.json'
    old = json.loads(path.read_text())
    proof = parent / 'independent-audit/refinement-receipt.json'
    audit = json.loads(proof.read_text())
    if (old['protocol'] != 'private-1non-diffusion-generated-only-reduced-xh-v1'
            or audit['plan_sha256'] != digest(path) or not audit['preparation']['raw_to_prepared_heavy_bits_exact']):
        raise ValueError('requires audited baseline')
    verify_bindings(old['source_code_and_input_sha256'], (path,))
    job = old['jobs'][0].copy()
    snapshot = Path(job['snapshot'])
    system = snapshot.parent / 'system.xml'
    history = parent / 'failure-analysis-2026-11-18/analysis.json'
    verify_bindings({str(p): json.loads(history.read_text())['input_sha256'][str(p)] for p in
                     (system, *(snapshot / n for n in ('base-system.xml','temporary-axis.json','initial-state.npz')))})
    root.mkdir(parents=True, exist_ok=False)
    work = root / 'execution/1NON-generated-only-seed-23'
    job.update(command=[sys.executable, '-m', 'scripts.shape_1non_control', 'worker', '--output', str(work)],
               snapshot=str(work / 'prepared'), timeout_seconds=1800-10)
    files = [path, proof, history, system, *snapshot.iterdir(), *Path(job['case']).iterdir(),
             *Path('scripts').glob('*.py'), *Path('src/dvbfixer').rglob('*.py'), Path(sys.executable),
             Path('tests/test_1non_shape_force.py')]
    plan = {**old, 'jobs':[job], 'protocol':'private-1non-generated-ca-shape-v1', 'parent_snapshot':str(snapshot),
            'parent_sha256':digest(path),'preparation':'reuse exact audited H/FF snapshot; charge 10s for shared preparation',
            'intervention':'group 0, 0.5*10000*max(0,sin(20deg)-signed_elevation_sine)^2 kJ/mol; all generated non-Gly CA only',
            'source_code_and_input_sha256':{str(p.resolve()):digest(p) for p in files}}
    save(root/'plan.json',plan)
    freeze_code(root/'plan.json')


def replay(output):
    plan = json.loads((output.parent.parent/'plan.json').read_text())
    def prepare(request, source, proposal, directory):
        snapshot = Path(plan['parent_snapshot'])
        shutil.copytree(snapshot, directory)
        audit_prepared(Path(plan['jobs'][0]['case']),directory)
        topology, _, points = read_temporary(directory)
        system = XmlSerializer.deserialize((snapshot.parent/'system.xml').read_text())
        gen_ids = set(request.generated_atoms)
        generated = {a.index for a in topology.atoms() if identity(a) in gen_ids}
        if len(generated) != 556 or {identity(a) for a in topology.atoms() if system.getParticleMass(a.index).value_in_unit(dalton)==0} != set(request.fixed_atoms):
            raise ValueError('snapshot mobile/fixed inventory differs')
        f = shape_force()
        centres = []
        for r in topology.residues():
            names = {a.name:a.index for a in r.atoms()}
            if names['CA'] in generated and r.name != 'GLY':
                ids = [names[n] for n in ('N','CA','C','CB')]
                if not set(ids) <= generated:
                    raise ValueError('partial centre not allowed')
                f.addBond(ids, [])
                centres.append([r.chain.id,r.id,r.insertionCode.strip()])
        if not centres:
            raise ValueError('no generated centres')
        system.addForce(f)
        save(output/'shape-centres.json',centres)
        return topology,system,points,generated,0
    worker(output,prepare=prepare)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=('declare','worker','run'))
    parser.add_argument('--parent',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    root=args.output.resolve()
    if args.mode=='declare':
        declare(args.parent.resolve(),root)
    elif args.mode=='worker':
        replay(root)
    else:
        from scripts.run_mac_modeller_comparison import execute
        execute(root)
