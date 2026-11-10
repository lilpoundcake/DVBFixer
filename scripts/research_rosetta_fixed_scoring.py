"""Plan-before-score research: import completed candidates without repacking.

Scoring is not a new geometry gate. All input heavy coordinates, target graph
and external links must survive import; only documented terminal OXT is allowed
as extra heavy context. Rebuild missing H separately for every completed input.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path

from scripts.research_rosetta_repack import identity, observed_atoms, pose_snapshot, write_json
from scripts.rosetta_fixed_contract import choose_scored_own, import_audit


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def plan(comparison: Path, output: Path) -> None:
    # Import the normal scientific stack only for planning, never in the PyRosetta runner.
    from dataclasses import asdict

    from scripts.compare_branched_completion import measure

    from dvbfixer.diagnose.geometry import VDW_RADII_ANGSTROM
    from dvbfixer.domain.structure_identity import ComponentInstanceRef, ExactAtomRef
    from dvbfixer.pdb_component_reconstruction import _component_ref, _parse_atom_lines
    from dvbfixer.protein_completion.anchored_source import plan_anchored_source
    from dvbfixer.protein_completion.branched_source import plan_branched_source
    from dvbfixer.protein_completion.contracts import Addition
    from dvbfixer.protein_completion.validation import sidechain_bonds

    raw = comparison.read_bytes()
    prior = json.loads(raw)
    prior_root, frozen_root = Path(prior['retained_directory']), Path(prior['plan_root'])
    frozen_raw = (frozen_root/'plan.json').read_bytes()
    if sha(frozen_raw) != prior['plan_sha256']:
        raise ValueError('prior frozen plan changed')
    frozen = json.loads(frozen_raw)
    for name, digest in frozen['source_hashes'].items():
        if sha(Path(name).read_bytes()) != digest:
            raise ValueError('prior frozen implementation changed')
    for name, digest in frozen['artifacts'].items():
        if sha((frozen_root/name).read_bytes()) != digest:
            raise ValueError('prior frozen artifact changed')
    cases, payloads = [], {}
    for declaration in frozen['cases']:
        name = declaration['id']
        if Path(name).name != name:
            raise ValueError('flat case identity required')
        row = prior['results'][name]
        for filename, digest in row['retained_files'].items():
            if sha((prior_root/name/filename).read_bytes()) != digest:
                raise ValueError('retained comparison artifact changed')
        source = (prior_root/name/'input.pdb').read_bytes()
        if sha(source) != frozen['artifacts'][declaration['source']]:
            raise ValueError('comparison source differs from frozen plan')
        payloads[name+'/source.pdb'] = source
        if row['status'] == 'shared_admission_refused':
            cases.append(dict(declaration, status=row['status'], reason=row['reason']))
            continue
        target = ComponentInstanceRef(**declaration['target'])
        anchored = declaration['mode'] == 'anchored'
        source_plan = plan_anchored_source(source, target) if anchored else plan_branched_source(source, target)
        request = source_plan.request
        attempts = json.loads((prior_root/name/'own-attempts.json').read_bytes())

        def rebuild(values):
            return tuple(Addition(ExactAtomRef(ComponentInstanceRef(**a['atom']['component']), a['atom']['name']),
                                  a['element'], tuple(a['xyz'])) for a in values)

        reference = (frozen_root/declaration['reference']).read_bytes()
        truth = {a.name: a.position for a in _parse_atom_lines(reference.decode().splitlines(keepends=True))
                 if _component_ref(a) == target}
        proposals = [('own', i, rebuild(a['candidate']['additions']), a['audit']['accepted'])
                     for i, a in enumerate(attempts['attempts'])]
        template = json.loads((prior_root/name/'template-proposal.json').read_bytes())
        proposals += [('template', 0, rebuild(template['additions']), False),
                      ('native_control', 0, tuple(Addition(ExactAtomRef(target, n), n[0], truth[n])
                                                  for n in request.missing_atoms), False)]
        cells = []
        for arm, index, additions, old_accepted in proposals:
            audit, pdb = measure(source, request, additions, anchored=anchored)
            filename = f'{name}/{arm}-{index}.pdb'
            payloads[filename] = pdb
            # Native reference error is deliberately absent from runner-facing cells.
            cells.append(dict(id=f'{arm}-{index}', arm=arm, index=index, pdb=filename,
                              old_source_accepted=old_accepted, common_accepted=audit.accepted,
                              common_audit=asdict(audit)))
        lookup = {(a.chain_id, a.sequence_number, a.insertion_code, a.name):
                  identity(dict(chain=a.chain_id, resid=int(a.sequence_number), icode=a.insertion_code or ' ',
                                resname=a.residue_name, atom=a.name))
                  for a in _parse_atom_lines(payloads[cells[0]['pdb']].decode().splitlines(keepends=True))}

        def key(atom, lookup=lookup):
            c = atom.component
            if c.model != 1 or c.occurrence != 1 or c.alternate_location:
                raise ValueError('scorer supports unambiguous single-model projections only')
            return lookup[(c.chain_id, c.sequence_number, c.insertion_code, atom.name)]

        external = [(key(a), key(b)) for a, b in request.links]
        internal = [(key(ExactAtomRef(target, a)), key(ExactAtomRef(target, b)))
                    for a, b in sidechain_bonds(request.residue_name)
                    if (target.chain_id, target.sequence_number, target.insertion_code, a) in lookup
                    and (target.chain_id, target.sequence_number, target.insertion_code, b) in lookup]
        reference_name = name+'/reference-reporting-only.pdb'
        payloads[reference_name] = reference
        cases.append(dict(declaration, status='admitted', cells=cells, expected_external=external,
                          target_internal=internal, old_selected_index=attempts['selected_index'],
                          reference_reporting_only=reference_name,
                          target_key=(target.chain_id, int(target.sequence_number), target.insertion_code or ' ', request.residue_name),
                          generated=[key(ExactAtomRef(target, n)) for n in request.missing_atoms]))
    output.mkdir(parents=True, exist_ok=False)
    for filename, data in payloads.items():
        destination = output/filename
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(data)
    write_json(output/'plan.json', dict(
        schema='dvbfixer.rosetta-fixed-scoring-plan.v2', cases=cases,
        prior_comparison_sha256=sha(raw), prior_plan_sha256=sha(frozen_raw),
        artifacts={p: sha(data) for p, data in payloads.items()},
        source_hashes={p: sha(Path(p).read_bytes()) for p in (
            'scripts/research_rosetta_fixed_scoring.py', 'scripts/rosetta_fixed_contract.py',
            'scripts/research_rosetta_repack.py', 'scripts/compare_branched_completion.py')},
        expected_rosetta_revision='1a56185c2592611dec4c9c75ddc9468cd2227c1f',
        score_function='ref2015', seed=7, radii=dict(VDW_RADII_ANGSTROM),
        init_flags='-constant_seed -jran 7 -mute all -multithreading:total_threads 1 -packing:pack_missing_sidechains false -packing:no_optH true',
        hydrogen_protocol='fresh import for each rounded heavy candidate; no H optimization; no repack or minimization',
        extra_heavy_policy='only unobserved OXT on Rosetta upper termini without declared outgoing C peptide link',
        virtual_atom_policy='typed virtual atoms are separate scorer-owned pseudoatoms, retained in snapshots; never exempt an input atom by name',
        selection_policy='minimum total score over complete old-source-and-common-gate-feasible own pool; stable index tie; any missing feasible score blocks selection',
        native_policy='native controls scored after own selection committed; never eligible for selection',
        topology_scope='source-derived external links and canonical target internal bonds; not a complete independent context internal bond-order audit',
        retrospective=True, public_capability=False))


def pose_edges(pose, target_key: tuple) -> tuple[list, list]:
    external, internal = set(), set()
    info = pose.pdb_info()

    def key(i, a):
        residue = pose.residue(i)
        return (info.chain(i), info.number(i), info.icode(i), residue.name3(), residue.atom_name(a).strip())

    for i in range(1, pose.total_residue()+1):
        residue = pose.residue(i)
        if key(i, 1)[:4] == target_key:
            for a in range(1, residue.nheavyatoms()+1):
                for b in residue.bonded_neighbor(a):
                    if b <= residue.nheavyatoms():
                        internal.add(tuple(sorted((key(i, a), key(i, b)))))
        for c in range(1, residue.n_possible_residue_connections()+1):
            j = residue.connected_residue_at_resconn(c)
            if not j:
                continue
            a = residue.residue_connect_atom_index(c)
            other = pose.residue(j)
            b = other.residue_connect_atom_index(residue.residue_connection_conn_id(c))
            if a <= residue.nheavyatoms() and b <= other.nheavyatoms():
                external.add(tuple(sorted((key(i, a), key(j, b)))))
    return sorted(external), sorted(internal)


def physical_snapshot(pose) -> tuple[list, list]:
    """Separate type-declared pseudoatoms, without ever whitelisting an atom name."""
    virtual = []
    info = pose.pdb_info()
    for i in range(1, pose.total_residue()+1):
        residue = pose.residue(i)
        for a in range(1, residue.natoms()+1):
            if residue.atom_type(a).is_virtual():
                xyz = residue.xyz(a)
                virtual.append(dict(chain=info.chain(i), resid=info.number(i), icode=info.icode(i),
                                    resname=residue.name3(), atom=residue.atom_name(a).strip(),
                                    xyz=[xyz.x, xyz.y, xyz.z], virtual_by_rosetta_type=True,
                                    rosetta_atom_type=residue.atom_type(a).name()))
    keys = {identity(a) for a in virtual}
    return [a for a in pose_snapshot(pose) if identity(a) not in keys], virtual


def score_cell(pyrosetta, sf, plan_root: Path, output: Path, case: dict, cell: dict) -> dict:
    source = (plan_root/cell['pdb']).read_bytes()
    expected = observed_atoms(source)
    pose = pyrosetta.pose_from_pdb(str(plan_root/cell['pdb']))
    before, virtual_before = physical_snapshot(pose)
    external, internal = pose_edges(pose, tuple(case['target_key']))
    expected_keys = {identity(a) for a in expected}
    has_outgoing = {tuple(a if a[-1] == 'C' else b)[:4] for a, b in case['expected_external']
                    if {a[-1], b[-1]} == {'C', 'N'}}
    allowed = []
    info = pose.pdb_info()
    for i in range(1, pose.total_residue()+1):
        residue = pose.residue(i)
        key = (info.chain(i), info.number(i), info.icode(i), residue.name3(), 'OXT')
        if residue.is_upper_terminus() and key[:4] not in has_outgoing and key not in expected_keys:
            allowed.append(key)
    audit = import_audit(expected, before, case['expected_external'], external,
                         case['target_internal'], internal, allowed)
    score = sf(pose) if audit['accepted'] else None
    after, virtual_after = physical_snapshot(pose)
    if before != after:
        raise ValueError('score function moved heavy atoms')
    if score is not None and not math.isfinite(score):
        raise ValueError('nonfinite Rosetta score')
    snapshot = dict(heavy_atoms=after, virtual_atoms=virtual_after, external=external, target_internal=internal, allowed_extra=allowed,
                    residues=[], hydrogens=[])
    for i in range(1, pose.total_residue()+1):
        residue = pose.residue(i)
        snapshot['residues'].append(dict(chain=info.chain(i), resid=info.number(i), icode=info.icode(i),
                                          name=residue.name(), atoms=[residue.atom_name(a).strip() for a in range(1, residue.natoms()+1)]))
        for a in range(residue.nheavyatoms()+1, residue.natoms()+1):
            if residue.atom_type(a).is_virtual():
                continue
            xyz = residue.xyz(a)
            snapshot['hydrogens'].append(dict(chain=info.chain(i), resid=info.number(i), icode=info.icode(i),
                                               atom=residue.atom_name(a).strip(), xyz=[xyz.x, xyz.y, xyz.z]))
    snapshot_bytes = json.dumps(snapshot, sort_keys=True, allow_nan=False).encode()
    (output/'pose-snapshot.json.gz').write_bytes(gzip.compress(snapshot_bytes, mtime=0))
    terms = {str(t): dict(raw=pose.energies().total_energies()[t], weight=sf.get_weight(t))
             for t in sf.get_nonzero_weighted_scoretypes()} if score is not None else {}
    fingerprint = sha(json.dumps(dict(residues=snapshot['residues'], external=external, target_internal=internal),
                                 sort_keys=True).encode())
    return dict(score_status='scored' if audit['accepted'] else 'import_contract_blocker',
                total_score=score, import_audit=audit, score_terms=terms,
                type_topology_fingerprint=fingerprint, hydrogen_count=len(snapshot['hydrogens']),
                virtual_atom_count=len(virtual_after), virtual_coordinates_unchanged_by_scoring=virtual_before == virtual_after,
                snapshot_sha256=sha((output/'pose-snapshot.json.gz').read_bytes()))


def run(plan_root: Path, output: Path) -> None:
    raw = (plan_root/'plan.json').read_bytes()
    record = json.loads(raw)
    if record['schema'] != 'dvbfixer.rosetta-fixed-scoring-plan.v2':
        raise ValueError('separate version-2 physical/virtual inventory plan required')
    for paths, root in ((record['source_hashes'], Path('.')), (record['artifacts'], plan_root)):
        for filename, digest in paths.items():
            if sha((root/filename).read_bytes()) != digest:
                raise ValueError('frozen scoring input or implementation changed')
    import pyrosetta

    if record['expected_rosetta_revision'] not in pyrosetta.version():
        raise ValueError('unexpected installed Rosetta build')
    pyrosetta.init(record['init_flags'])
    sf = pyrosetta.create_score_function(record['score_function'])
    output.mkdir(parents=True, exist_ok=False)
    write_json(output/'environment.json', dict(pyrosetta_version=pyrosetta.version(),
                                              weights={str(t): sf.get_weight(t) for t in sf.get_nonzero_weighted_scoretypes()}))
    (output/'executed-plan.json').write_bytes(raw)
    outcomes = {}
    for case in record['cases']:
        name = case['id']
        folder = output/name
        folder.mkdir()
        if case['status'] != 'admitted':
            outcomes[name] = dict(status=case['status'], reason=case['reason'])
            write_json(folder/'result.json', outcomes[name])
            continue
        scored, own, fingerprint = [], [], None

        def execute(cell, folder=folder, case=case, scored=scored, own=own):
            nonlocal fingerprint
            destination = folder/cell['id']
            destination.mkdir()
            try:
                result = score_cell(pyrosetta, sf, plan_root, destination, case, cell)
                if result['score_status'] == 'scored':
                    if fingerprint is None:
                        fingerprint = result['type_topology_fingerprint']
                    elif fingerprint != result['type_topology_fingerprint']:
                        result.update(score_status='inconsistent_type_topology_blocker', total_score=None)
            except Exception as error:
                result = dict(score_status='execution_blocker', total_score=None,
                              error_type=type(error).__name__, error=str(error))
            row = dict(cell, **result)
            write_json(destination/'score.json', row)
            scored.append(row)
            if cell['arm'] == 'own':
                own.append(row)

        for cell in case['cells']:
            if cell['arm'] != 'native_control':
                execute(cell)
        # Commit selection BEFORE importing any native-reference control.
        selection = choose_scored_own(own)
        write_json(folder/'selection-before-native.json', selection)
        for cell in case['cells']:
            if cell['arm'] == 'native_control':
                execute(cell)
        outcomes[name] = dict(status='evaluated', selection=selection, cells=scored)
        write_json(folder/'result.json', outcomes[name])
    write_json(output/'evidence.json', dict(schema='dvbfixer.rosetta-fixed-scoring-raw.v2',
                                          plan_sha256=sha(raw), plan_root=str(plan_root),
                                          results=outcomes, production_enabled=False))
    print(json.dumps(dict(cases=len(outcomes), admitted=sum(r['status'] == 'evaluated' for r in outcomes.values()),
                         cells=sum(len(r.get('cells', [])) for r in outcomes.values()),
                         scored=sum(c['score_status'] == 'scored' for r in outcomes.values() for c in r.get('cells', []))), indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('plan')
    p.add_argument('--comparison', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p = sub.add_parser('run')
    p.add_argument('--plan', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.command == 'plan':
        plan(args.comparison, args.output)
    else:
        run(args.plan, args.output)


if __name__ == '__main__':
    main()
