"""Read-only source/preparation audit; no new H controls or optimizer calls."""
from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

from scripts.branched_template_baseline import preflight
from scripts.enroll_methyl_diversity import SOURCES
from scripts.freeze_methyl_diversity import admission
from scripts.prospective_independent_cohort import project_first
from scripts.research_dunbrack_prior import sha, write
from scripts.research_ha_minimizer_path import verify

from dvbfixer.ffutils.dat import DatRecord


def check_controls(controls, label, digest):
    if (controls['id'] != label or controls['plan_sha256'] != digest
            or controls['hooks_restored'] is not True or controls['native_coordinates_read'] is not False
            or not controls['fixer_instances'] or not controls['H_calls']):
        raise ValueError('preparation observer contract differs')
    if any(c != dict(platform='Reference') for c in controls['fixer_instances']):
        raise ValueError('preparation fixer platform differs')
    if any(c != dict(platform='Reference', seed=1, state_restored=True) for c in controls['H_calls']):
        raise ValueError('preparation H seed/platform/state differs')


def audit(plan_path, root, output):
    raw = plan_path.read_bytes()
    plan = json.loads(raw)
    verify(plan)
    if [(c['id'], c['source']) for c in plan['cases']] != list(SOURCES):
        raise ValueError('fixed eight-source enrollment differs')
    if preflight() != plan['runtime']:
        raise ValueError('preparation runtime/resource binding differs')
    run_raw = (root/'evidence.json').read_bytes()
    run = json.loads(run_raw)
    contract_raw = (root/'execution-contract.json').read_bytes()
    contract = json.loads(contract_raw)
    if (run['source_plan_sha256'] != sha(raw) or (root/'executed-source-plan.json').read_bytes() != raw
            or run['execution_contract_sha256'] != sha(contract_raw)
            or contract['source_plan_sha256'] != sha(raw)
            or contract['orchestrator_sha256'] != sha(Path('scripts/prepare_methyl_diversity.py').read_bytes())
            or [c['id'] for c in run['cases']] != [c['id'] for c in plan['cases']]):
        raise ValueError('preparation execution/implementation binding differs')
    reports = []
    for source, result in zip(plan['cases'], run['cases'], strict=True):
        folder = root/source['id']
        for name, digest in result['artifacts'].items():
            if sha((folder/name).read_bytes()) != digest:
                raise ValueError('preparation artifact changed')
        if json.loads((folder/'result.json').read_bytes()) != result:
            raise ValueError('per-source result differs from complete run')
        original = Path(source['source']).read_bytes()
        projected, chain = project_first(original)
        if sha(original) != source['source_sha256'] or sha(projected) != source['projection_sha256'] or chain != source['chain']:
            raise ValueError('independently projected source differs')
        projection = Path(source['projection'])
        if projection.read_bytes() != projected or (projection.parent/'original.pdb').read_bytes() != original:
            raise ValueError('retained source/projection bytes differ')
        try:
            inventory = admission(projected)
            refusal = None
        except ValueError as error:
            inventory, refusal = None, str(error)
        if (inventory is None) != (source['status'] == 'refused_source') or result['source_status'] != source['status']:
            raise ValueError('independent source admission differs')
        report = dict(id=source['id'], chain=chain, source_sha256=sha(original),
            projection_sha256=sha(projected), status=result['status'], artifacts=result['artifacts'],
            returncode=result['returncode'], source_inventory=inventory, refusal=refusal)
        if refusal is not None:
            if refusal != source['refusal'] or result['status'] != 'refused_source' or result['artifacts']:
                raise ValueError('source refusal was changed or executed')
        elif result['status'] == 'prepared':
            if result['returncode'] != 0 or (folder/'input.pdb').read_bytes() != projected:
                raise ValueError('prepared input ownership/return code differs')
            controls = json.loads((folder/'prepare-controls.json').read_bytes())
            check_controls(controls, source['id'], sha(raw))
            dat = DatRecord.load(folder/'prepared.dat')
            # Canonical schema decode; do not interpret raw sidecar JSON here.
            if len(dat.added_keys()) != len(dat.added_atoms):
                raise ValueError('duplicate canonical added-atom identities')
            report.update(controls=controls, dat_added_atoms=len(dat.added_atoms),
                dat_added_heavy=sum(a['element'] not in {'H', 'D'} for a in dat.added_atoms),
                dat_variant_overrides=dat.variant_overrides,
                prepared_coordinate_records=sum(l.startswith((b'ATOM  ', b'HETATM')) for l in (folder/'prepared.pdb').read_bytes().splitlines()),
                borrowed_input_exact=True)
        reports.append(report)
    verify(plan)
    summary = dict(enrolled=len(reports), statuses=dict(Counter(c['status'] for c in reports)),
        scope='source admission and legacy preparation only; no graph-H trial or scientific qualification')
    write(output, dict(schema='dvbfixer.methyl-diversity-admission-audit.v1', cases=reports,
        summary=summary, source_plan_sha256=sha(raw), enrollment_sha256=plan['enrollment_sha256'],
        run_sha256=sha(run_raw), execution_contract=contract, runtime=plan['runtime'],
        auditor_sha256=sha(Path(__file__).read_bytes()), production_changed=False,
        new_H_control_executed=False, scientific_qualification=False))
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    audit(args.plan, args.run, args.output)
