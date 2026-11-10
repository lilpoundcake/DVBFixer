"""Private mechanism probes, not an observation-preserving completion backend.

Plan first; each probe runs in a fresh PyRosetta process. Never restore moved
observations to conceal the residue-level packer's different ownership contract.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
from pathlib import Path

CASES = ('1TIM-VAL39-CG1', '1TIM-VAL39-CB-CG1-CG2',
         '1EMV-THR11-CG2', '1EMV-THR11-CB-OG1-CG2')


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def write_json(path: Path, value: object) -> None:
    with path.open('x') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')


def observed_atoms(source: bytes) -> list[dict]:
    rows = []
    seen = set()
    for line in source.decode('ascii').splitlines():
        if not line.startswith('ATOM  '):
            continue
        if line[16].strip() or float(line[54:60]) != 1:
            raise ValueError('probe requires explicit unit-occupancy no-altloc projection')
        row = dict(chain=line[21], resid=int(line[22:26]), icode=line[26],
                   resname=line[17:20], atom=line[12:16].strip(),
                   xyz=[float(line[i:i+8]) for i in (30, 38, 46)])
        key = identity(row)
        if key in seen:
            raise ValueError('duplicate observed identity')
        seen.add(key)
        rows.append(row)
    return rows


def identity(row: dict) -> tuple:
    return tuple(row[k] for k in ('chain', 'resid', 'icode', 'resname', 'atom'))


def movement(observations: list[dict], snapshot: list[dict]) -> dict:
    points = {identity(a): a['xyz'] for a in snapshot}
    missing, moved = [], []
    maximum = 0.
    for atom in observations:
        key = identity(atom)
        if key not in points:
            missing.append(atom)
            continue
        distance = math.dist(atom['xyz'], points[key])
        maximum = max(maximum, distance)
        if distance > 1e-7:
            moved.append(dict(atom, distance_angstrom=distance, after_xyz=points[key]))
    return dict(maximum_angstrom=maximum, missing_observations=missing,
                moved_above_1e_7_angstrom=moved, observations_preserved=not missing and not moved)


def plan(frozen: Path, output: Path) -> None:
    record = json.loads((frozen/'plan.json').read_bytes())
    cases = {c['id']: c for c in record['cases']}
    payloads = {}
    for name in CASES:
        case = cases[name]
        for artifact in (case['source'], case['reference']):
            data = (frozen/artifact).read_bytes()
            if sha(data) != record['artifacts'][artifact]:
                raise ValueError('frozen comparison input changed')
            observed_atoms(data)
            payloads[artifact] = data
    output.mkdir(parents=True, exist_ok=False)
    for name, data in payloads.items():
        (output/name).write_bytes(data)
    write_json(output/'plan.json', dict(
        schema='dvbfixer.rosetta-repack-mechanism-plan.v1',
        prior_plan_sha256=sha((frozen/'plan.json').read_bytes()),
        harness_sha256=sha(Path(__file__).read_bytes()),
        cases=[cases[n] for n in CASES],
        artifacts={n: sha(data) for n, data in payloads.items()},
        seeds=[7], auto_pack_on_import=[False, True], preserve_input_cb=[False, True],
        score_function='ref2015', explicit_task='repack only selected target; IncludeCurrent; no extra chi',
        import_h_optimization=False, scientific_accuracy_benchmark=False,
        reference_use='post-proposal reporting only; never supplied to PyRosetta',
        output_policy='raw diagnostic poses and snapshots only; no production materialization'))


def pose_snapshot(pose) -> list[dict]:
    rows = []
    info = pose.pdb_info()
    for i in range(1, pose.total_residue()+1):
        residue = pose.residue(i)
        for a in range(1, residue.nheavyatoms()+1):
            p = residue.xyz(a)
            rows.append(dict(chain=info.chain(i), resid=info.number(i), icode=info.icode(i),
                             resname=residue.name3(), atom=residue.atom_name(a).strip(), xyz=[p.x, p.y, p.z]))
    return rows


def probe(plan_root: Path, case_id: str, output: Path, *, auto: bool, preserve_cb: bool, seed: int) -> None:
    raw_plan = (plan_root/'plan.json').read_bytes()
    record = json.loads(raw_plan)
    if sha(Path(__file__).read_bytes()) != record['harness_sha256']:
        raise ValueError('probe implementation changed after planning')
    if seed not in record['seeds']:
        raise ValueError('seed not predeclared')
    for name, digest in record['artifacts'].items():
        if Path(name).name != name or sha((plan_root/name).read_bytes()) != digest:
            raise ValueError('planned input changed')
    case = next(c for c in record['cases'] if c['id'] == case_id)
    source = (plan_root/case['source']).read_bytes()
    observations = observed_atoms(source)
    output.mkdir(parents=True, exist_ok=False)
    (output/'input.pdb').write_bytes(source)
    import pyrosetta

    flags = (f'-constant_seed -jran {seed} -mute all -multithreading:total_threads 1 '
             f'-packing:pack_missing_sidechains {str(auto).lower()} -packing:no_optH true '
             f'-packing:preserve_c_beta {str(preserve_cb).lower()} -packing:repack_only')
    pyrosetta.init(flags)
    r = pyrosetta.rosetta
    pose = pyrosetta.pose_from_pdb(str(output/'input.pdb'))
    initial = pose_snapshot(pose)
    write_json(output/'loaded-snapshot.json', initial)
    pose.dump_pdb(str(output/'loaded-diagnostic.pdb'))
    target = case['target']
    index = pose.pdb_info().pdb2pose(target['chain_id'], int(target['sequence_number']), target['insertion_code'] or ' ')
    if not index or pose.residue(index).name3() != case['residue_name']:
        raise ValueError('target identity changed at import')
    sf = pyrosetta.create_score_function(record['score_function'])
    before_score = sf(pose)
    task = r.core.pack.task.TaskFactory.create_packer_task(pose)
    task.initialize_from_command_line()
    task.restrict_to_repacking()
    task.or_include_current(True)
    for i in range(1, pose.total_residue()+1):
        if i != index:
            task.nonconst_residue_task(i).prevent_repacking()
    # Build a diagnostic candidate set independently; the mover rebuilds its own set.
    graph = r.core.pack.create_packer_graph(pose, sf, task)
    rotamers = r.core.pack.rotamer_set.RotamerSet_()
    rotamers.set_resid(index)
    rotamers.build_rotamers(pose, sf, task, graph)
    candidates = []
    info = pose.pdb_info()
    for j in range(1, rotamers.num_rotamers()+1):
        residue = rotamers.rotamer(j)
        rows = []
        for a in range(1, residue.nheavyatoms()+1):
            p = residue.xyz(a)
            rows.append(dict(chain=info.chain(index), resid=info.number(index), icode=info.icode(index),
                             resname=residue.name3(), atom=residue.atom_name(a).strip(), xyz=[p.x, p.y, p.z]))
        target_observed = [o for o in observations if o['chain'] == target['chain_id']
                           and o['resid'] == int(target['sequence_number'])
                           and o['icode'] == (target['insertion_code'] or ' ')]
        candidates.append(dict(index=j, name=residue.name(), chis=[residue.chi(k) for k in range(1, residue.nchi()+1)],
                               target_observed_movement=movement(target_observed, rows), atoms=rows))
    write_json(output/'diagnostic-rotamers.json', candidates)
    mover = r.protocols.minimization_packing.PackRotamersMover(sf, task)
    mover.apply(pose)
    final = pose_snapshot(pose)
    write_json(output/'repacked-snapshot.json', final)
    pose.dump_pdb(str(output/'repacked-diagnostic.pdb'))
    weights = {str(t): sf.get_weight(t) for t in sf.get_nonzero_weighted_scoretypes()}
    write_json(output/'result.json', dict(
        schema='dvbfixer.rosetta-repack-mechanism-result.v1', plan_sha256=sha(raw_plan), case=case,
        python=platform.python_version(), pyrosetta_version=pyrosetta.version(), init_flags=flags,
        auto_pack=auto, preserve_input_cb=preserve_cb, seed=seed, score_weights=weights,
        pose_residues=pose.total_residue(), task=str(task), target_pose_index=index,
        task_preserve_c_beta=task.residue_task(index).preserve_c_beta(),
        diagnostic_rotamer_count=len(candidates),
        load_vs_source=movement(observations, initial), repack_vs_source=movement(observations, final),
        repack_vs_loaded=movement(initial, final), score_before=before_score, score_after=sf(pose),
        scientific_qualification=False, observation_restoration=False,
        files={p.name: sha(p.read_bytes()) for p in sorted(output.iterdir())}))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('plan')
    p.add_argument('--frozen', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p = sub.add_parser('probe')
    p.add_argument('--plan', type=Path, required=True)
    p.add_argument('--case', required=True, choices=CASES)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--auto', action='store_true')
    p.add_argument('--preserve-cb', action='store_true')
    p.add_argument('--seed', type=int, default=7)
    args = parser.parse_args()
    if args.command == 'plan':
        plan(args.frozen, args.output)
    else:
        probe(args.plan, args.case, args.output, auto=args.auto, preserve_cb=args.preserve_cb, seed=args.seed)


if __name__ == '__main__':
    main()
