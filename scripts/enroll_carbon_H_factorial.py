"""Enroll existing native H starts for a fixed carbon-CH2/CH3 factorial study.

Metadata only: no coordinates, groups, Contexts, minimization or new inputs.
Keep all reconstruction-stress conditions and inherited source refusals.
"""
from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

from scripts.enroll_methyl_diversity import SOURCES
from scripts.research_dunbrack_prior import sha, write

ARMS = ('unchanged', 'all_methyl_zero', 'all_methylene', 'combined_zero')
STRESS = ('1SHG-VAL9', '1SHG-VAL23', '1SHG-THR24', '1SHG-THR32')


def enroll(legacy_plan, capture_plan, capture_root, output):
    legacy_raw = legacy_plan.read_bytes()
    legacy = json.loads(legacy_raw)
    expected = [(label, seed, arm) for label in STRESS for seed in (1, 2) for arm in ('raw', 'beta_only', 'full_generated')]
    if [(c['id'], c['seed'], c['source_arm']) for c in legacy['cases']] != expected:
        raise ValueError('complete original ALL24 reconstruction-stress cohort required')
    capture_raw = capture_plan.read_bytes()
    capture = json.loads(capture_raw)
    run_path = capture_root/'evidence.json'
    run_raw = run_path.read_bytes()
    run = json.loads(run_raw)
    if (capture['schema'] != 'dvbfixer.methyl-diversity-capture-plan.v2'
            or run['plan_sha256'] != sha(capture_raw)
            or [c['id'] for c in capture['cases']] != [label for label, _ in SOURCES]
            or [(c['id'], c['seed']) for c in run['cases']] != [(label, seed) for label, _ in SOURCES for seed in (1, 2)]
            or Counter(c['status'] for c in run['cases']) != {'captured':10, 'refused_source':6}):
        raise ValueError('complete fixed eight-input/two-seed capture cohort required')
    cases = []
    for c in legacy['cases']:
        source = Path(c['source'])
        cases.append(dict(id=c['id'], seed=c['seed'], conditioning=c['source_arm'],
            context_kind='geometry-reconstruction-stress', pdb='1SHG', status='native_start_retained',
            system=str(source/'h-internal-system.xml'), atoms=str(source/'h-internal-before.json.gz'),
            pre_H=str(source/'h-input.json.gz'), topology=str(source/'topology.json'),
            baseline_path=c['path'], transport_contract='legacy-H-capture; verify native bond graph/index/complete identities before replay'))
    for c in run['cases']:
        source = capture_root/c['id']/f"seed-{c['seed']}"
        cases.append(dict(id=c['id'], seed=c['seed'], conditioning='prepared-view',
            context_kind='prepared-protein-view', pdb=c['id'], status=c['status'],
            system=str(source/'H-system.xml'), atoms=str(source/'H-topology.json'),
            pre_H=str(source/'H-input.json.gz'), topology=str(source/'H-topology.json'),
            baseline_path=str(source/'H-path.json.gz'), ownership=str(source/'H-ownership.json'),
            transport_contract='actual native ForceField topology/object lineage from V2 capture; inherited refusal never executed'))
    inputs = {str(legacy_plan):sha(legacy_raw), str(capture_plan):sha(capture_raw), str(run_path):sha(run_raw)}
    for name in ('scripts/enroll_carbon_H_factorial.py', 'scripts/enroll_methyl_diversity.py',
                 'scripts/graph_methyl_geometry.py', 'scripts/research_methylene_geometry.py'):
        inputs[name] = sha(Path(name).read_bytes())
    record = dict(schema='dvbfixer.carbon-H-factorial-enrollment.v1', cases=cases, inputs=inputs,
        summary=dict(contexts=40, native_starts=34, inherited_source_seed_refusals=6,
            arm_cells=160, potential_replays=136, inherited_refused_arm_cells=24),
        policy=dict(arms=ARMS, H_seeds=[1,2], platform='Reference', tolerance_kj_mol_nm=1., maximum_iterations=50,
            per_cell_seconds=900, total_seconds=14400,
            timing='existing-input qualification fixed before factorial adapter/initializers are evaluated; previously studied cases, not blinded proteins or temporal holdout',
            enrollment='ALL24 original 1SHG reconstruction-stress starts plus ALL16 fixed eight-input source/seed capture entries; no replacement, source cleanup or new seed',
            scope='24 stress conditions on ONE 1SHG structure plus five admitted prepared protein views/two seeds; inherited three source refusals/two seeds remain',
            stress='intentionally poor generated-heavy source conditions remain explicit stress strata, not newly admitted healthy structures or a claim H initialization repairs heavy geometry',
            ownership='actual native bond graph, full case-sensitive identities and pre-H/newly mobile H ownership; no measured/kept H rewrite or guessed label aliases',
            construction='reuse frozen graph CH3 phase-zero second-shell frame and heavy-compatible CH2 signed heavy-plane formula; each original C-H radius and all nonowned coordinates retained',
            combination='disjoint CH2/CH3 H ownership; apply both fixed constructors without energy/phase search, revisits, H regeneration or force/mass/parameter edits',
            admission='preserve source/runtime and undefined/infeasible construction refusals; report full and matched denominators; never repair/crop to gain entry',
            controls='unchanged nm endpoints must reproduce each original native H endpoint; known positive/intersection replays are computational controls, NOT independent new samples',
            audit='all callbacks/native mobile gradients, focal and global H geometry, input heavy-shape findings, radii/slots/nonownership and all adverse contexts; no result-based arm selection',
            limits='H-stage only; no whole-system acceptance, native accuracy, public candidate/PDB, new chemistry/platform qualification or production capability'),
        coordinates_read=False, groups_discovered=False, arms_evaluated=False, production_changed=False)
    write(output, record)
    print(json.dumps(record['summary'], indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--legacy-plan', type=Path, required=True)
    parser.add_argument('--capture-plan', type=Path, required=True)
    parser.add_argument('--capture-root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    enroll(args.legacy_plan, args.capture_plan, args.capture_root, args.output)
