"""Offline index/view only; preserve a failed exclusive-write orchestration."""
from __future__ import annotations

import argparse
from pathlib import Path

from scripts.audit_thr2_native_capture import read
from scripts.research_dunbrack_prior import sha, write
from scripts.research_hydrogen_guard_native import check


def assemble_case(index, case, folder, recorded):
    artifacts = {str(p.relative_to(folder)):sha(p.read_bytes()) for p in sorted(folder.rglob('*')) if p.is_file()}
    if recorded is not None:
        if recorded['artifacts'] != artifacts:
            raise ValueError('previously indexed artifacts changed')
        return recorded
    returned = True
    for number in (1, 2, 3):
        phase = folder/f'phase-{number}'
        if not all((phase/name).is_file() for name in ('result.json', 'path.json.gz', 'returned-state.xml')):
            returned = False
            continue
        result, path = read(phase/'result.json'), read(phase/'path.json.gz')
        if not result['returned'] or not path['returned'] or result['callbacks'] != len(path['callbacks']):
            returned = False
    return dict(index=index, **case, status='phases-returned-exit-unrecorded' if returned else 'partial-exit-unrecorded',
                returncode=None, artifacts=artifacts)


def recover(plan_path, root, view, log, output):
    import openmm
    from openmm.unit import nanometer

    plan = read(plan_path)
    check(plan)
    partial = read(root/'evidence.json')
    if partial['plan_sha256'] != sha(plan_path.read_bytes()) or (root/'executed-plan.json').read_bytes() != plan_path.read_bytes():
        raise ValueError('original plan binding mismatch')
    known = {c['index']:c for c in partial['cases']}
    cases = []
    for index, case in enumerate(plan['cases']):
        folder = root/f'path-{index:02d}'
        row = assemble_case(index, case, folder, known.get(index))
        for number in (1, 2, 3):
            phase = folder/f'phase-{number}'
            if (phase/'returned-state.xml').exists() and (phase/'path.json.gz').exists():
                state = openmm.XmlSerializer.deserialize((phase/'returned-state.xml').read_text())
                positions = state.getPositions(asNumpy=True).value_in_unit(nanometer).tolist()
                if positions != read(phase/'path.json.gz')['final_nm']:
                    raise ValueError('serialized native return differs from path')
        cases.append(row)
    view.mkdir(parents=True, exist_ok=False)
    (view/'executed-plan.json').write_bytes(plan_path.read_bytes())
    for index in range(len(plan['cases'])):
        (view/f'path-{index:02d}').symlink_to(root/f'path-{index:02d}', target_is_directory=True)
    write(view/'evidence.json', dict(plan_sha256=sha(plan_path.read_bytes()), cases=cases))
    write(output, dict(schema='dvbfixer.private-guard-offline-summary-recovery.v1',
        plan_sha256=sha(plan_path.read_bytes()),
        original_partial_index_sha256=sha((root/'evidence.json').read_bytes()),
        original_orchestrator_log_sha256=sha(log.read_bytes()),
        recovered_index_sha256=sha((view/'evidence.json').read_bytes()),
        source_sha256=sha(Path(__file__).read_bytes()),
        cases=[{k:v for k, v in c.items() if k != 'artifacts'} for c in cases],
        original_orchestrator_failed=True, original_files_modified=False,
        new_contexts=0, new_minimizations=0,
        limits='Offline accounting only. Unrecorded process exit is unknown, not zero; native return flags/States are separately checked.'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('plan', 'root', 'view', 'log', 'output'):
        parser.add_argument('--'+name, type=Path, required=True)
    args = parser.parse_args()
    recover(args.plan.resolve(), args.root.resolve(), args.view.resolve(), args.log.resolve(), args.output.resolve())
