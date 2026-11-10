"""Independent saved-snapshot audit; no PyRosetta import or new proposals."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

from scripts.branched_stereo_oracle import rdkit_labels
from scripts.research_rosetta_repack import movement, observed_atoms, write_json


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def summarize(root: Path, plan_root: Path, contact_root: Path, source_root: Path, output: Path) -> None:
    plan = json.loads((plan_root/'plan.json').read_bytes())
    index = json.loads((root/'execution-index.json').read_bytes())
    assert len(index['outcomes']) == 16
    cases = {c['id']: c for c in plan['cases']}
    results = {}
    for execution in index['outcomes']:
        if execution['returncode'] != 0 or not execution['result_present']:
            results[execution['id']] = dict(status='execution_blocker', execution=execution)
            continue
        folder = root/execution['id']
        row = json.loads((folder/'result.json').read_bytes())
        assert row['plan_sha256'] == sha(plan_root/'plan.json')
        for name, digest in row['files'].items():
            assert sha(folder/name) == digest, (folder, name)
        case = cases[row['case']['id']]
        assert sha(folder/'input.pdb') == plan['artifacts'][case['source']]
        source = observed_atoms((folder/'input.pdb').read_bytes())
        loaded = json.loads((folder/'loaded-snapshot.json').read_bytes())
        final = json.loads((folder/'repacked-snapshot.json').read_bytes())
        assert movement(source, loaded) == row['load_vs_source']
        assert movement(source, final) == row['repack_vs_source']
        assert movement(loaded, final) == row['repack_vs_loaded']
        rotamers = json.loads((folder/'diagnostic-rotamers.json').read_bytes())
        target = case['target']

        def is_target(atom, target=target):
            return (atom['chain'] == target['chain_id'] and atom['resid'] == int(target['sequence_number'])
                    and atom['icode'] == (target['insertion_code'] or ' '))

        target_observed = [a for a in source if is_target(a)]
        for candidate in rotamers:
            assert movement(target_observed, candidate['atoms']) == candidate['target_observed_movement']
        points = {a['atom']: tuple(a['xyz']) for a in final if is_target(a)}
        reference = {a['atom']: tuple(a['xyz']) for a in observed_atoms((plan_root/case['reference']).read_bytes()) if is_target(a)}
        cip = rdkit_labels(points, case['residue_name'])
        backbone = [a for a in source if a['atom'] in ('N', 'CA', 'C', 'O')]
        identity_preserved = {(a['chain'], a['resid'], a['icode'], a['resname']) for a in source} == {
            (a['chain'], a['resid'], a['icode'], a['resname']) for a in final}
        assert identity_preserved
        assert movement(backbone, final)['observations_preserved']
        requested = case['mask']
        generated_rmsd = math.sqrt(sum(math.dist(points[n], reference[n])**2 for n in requested)/len(requested))
        results[execution['id']] = dict(
            status='completed', case_id=case['id'], auto_pack=row['auto_pack'], preserve_cb=row['preserve_input_cb'],
            task_preserve_c_beta=row['task_preserve_c_beta'],
            target_rdkit_cip=cip, target_expected_cip=({'CA': 'S', 'CB': 'R'} if case['residue_name'] == 'THR' else {'CA': 'S'}),
            residue_identity_preserved=identity_preserved, backbone_preserved=True,
            load_vs_source=row['load_vs_source'], repack_vs_source=row['repack_vs_source'],
            diagnostic_rotamer_count=len(rotamers),
            diagnostic_observation_preserving_rotamers=sum(c['target_observed_movement']['observations_preserved'] for c in rotamers),
            diagnostic_unique_chi1_count=len({round(c['chis'][0], 5) for c in rotamers}),
            generated_exact_name_rmsd_reporting_only_angstrom=generated_rmsd,
            ownership_comparable_to_generated_only_method=row['repack_vs_source']['observations_preserved'],
            raw_total_score_before=row['score_before'], raw_total_score_after=row['score_after'],
            result_sha256=sha(folder/'result.json'), retained_files=row['files'])
    contact = json.loads((contact_root/'evidence.json').read_bytes())
    manifest = json.loads((source_root/'source-manifest.json').read_bytes())
    for name, row in manifest['files'].items():
        assert sha(source_root/name) == row['sha256']
    summary = dict(
        schema='dvbfixer.rosetta-repack-mechanism-evidence.v1',
        installed_rosetta_revision=manifest['revision'],
        pyrosetta_version=json.loads((root/index['outcomes'][0]['id']/'result.json').read_bytes())['pyrosetta_version'],
        source_manifest=manifest, retained_source_directory=str(source_root),
        plan=str(plan_root/'plan.json'), plan_sha256=sha(plan_root/'plan.json'),
        retained_run_directory=str(root), execution_index_sha256=sha(root/'execution-index.json'),
        harness_sha256=sha(Path('scripts/research_rosetta_repack.py')),
        summarizer_sha256=sha(Path(__file__)), contact_evidence_sha256=sha(contact_root/'evidence.json'),
        contact=contact, results=results,
        aggregate=dict(planned=16, completed=sum(r['status'] == 'completed' for r in results.values()),
                       load_preserves_observations=sum(r.get('load_vs_source', {}).get('observations_preserved', False) for r in results.values()),
                       repack_preserves_observations=sum(r.get('repack_vs_source', {}).get('observations_preserved', False) for r in results.values())),
        protocol_failures_retained=[
            'tmp/pdbfixer-research/rosetta-mechanism-run-01/1TIM-VAL39-CG1-auto0-cb0/execution-error.json',
            'tmp/pdbfixer-research/rosetta-contact-01/execution-error.json'],
        scientific_qualification=False, matched_accuracy_comparison=False,
        limitations=['Four retrospective masks, one seed, one installed build; no broad or exact-repeat qualification.',
                     'Diagnostic rotamer set is built separately from mover; it is not a captured annealer trace.',
                     'PyRosetta uses its bundled empirical data; no table export, redistribution or own-reader authority established.',
                     'Target-only CIP checks are not an independent whole-protein stereo qualification.',
                     'Contact contraction controls are unphysical pair diagnostics, not minimized structures or total energies.'])
    write_json(output, summary)
    print(json.dumps(summary['aggregate'], indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--contact', type=Path, required=True)
    parser.add_argument('--sources', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    summarize(args.run, args.plan, args.contact, args.sources, args.output)


if __name__ == '__main__':
    main()
