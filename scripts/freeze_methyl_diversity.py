"""Source admission for the immutable graph-methyl convenience enrollment.

No H-control evaluation or source cleanup. Preserve original/projection bytes,
all failures and the first-model/first-ATOM-chain policy already enrolled.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from scripts.enroll_methyl_diversity import SOURCES
from scripts.prospective_independent_cohort import project_first
from scripts.research_dunbrack_prior import sha, write
from scripts.research_ha_minimizer_path import verify

from dvbfixer.protein_completion.branched_source import _thr_beta_checks
from dvbfixer.protein_completion.links import parse_declared_links
from dvbfixer.protein_completion.source import _global_ca, _layout


def admission(projected):
    lines, atoms, _ = _layout(projected)
    _global_ca(atoms, allow_missing_cb=None)
    _thr_beta_checks(atoms, pending_target=None)
    serials = {a.serial for a in atoms}
    if len(serials) != len(atoms):
        raise ValueError('duplicate source atom serial')
    by_ref = {(a.chain_id, a.sequence_number, a.insertion_code, a.name): a for a in atoms}
    links = parse_declared_links(lines, atoms)
    for link in links:
        def point(ref):
            component = ref.component
            return by_ref[(component.chain_id, component.sequence_number, component.insertion_code, ref.name)].position
        if link.reported_distance is not None and abs(math.dist(point(link.left), point(link.right))-link.reported_distance) > .1:
            raise ValueError('declared link distance disagrees with retained source atoms')
    conect = 0
    for line in lines:
        if not line.startswith('CONECT'):
            continue
        fields = line[6:].rstrip('\n')
        values = [int(fields[i:i+5]) for i in range(0, len(fields), 5) if fields[i:i+5].strip()]
        if len(values) < 2 or any(value not in serials for value in values):
            raise ValueError('dangling/empty source CONECT after the enrolled projection')
        if values[0] in values[1:] or len(set(values[1:])) != len(values[1:]):
            raise ValueError('self or unsupported multiplicity source CONECT')
        conect += 1
    return dict(atoms=len(atoms), residues=len({(a.chain_id, a.sequence_number, a.insertion_code) for a in atoms}),
        declared_links=len(links), conect_records=conect, heavy_only=True,
        scope='shared canonical source layout, global CA/THR-beta and exact declared-link endpoint checks; not complete-complex/chemistry qualification')


def freeze(enrollment_path, runtime_path, root, output):
    raw = enrollment_path.read_bytes()
    enrollment = json.loads(raw)
    if [(c['id'], c['source']) for c in enrollment['sources']] != list(SOURCES):
        raise ValueError('source cohort changed from the fixed enrollment definition')
    for name, digest in enrollment['inputs'].items():
        if sha(Path(name).read_bytes()) != digest:
            raise ValueError(f'enrolled source/policy changed: {name}')
    runtime_raw = runtime_path.read_bytes()
    runtime = json.loads(runtime_raw)
    root.mkdir(parents=True, exist_ok=False)
    inputs = {**enrollment['inputs'], str(enrollment_path): sha(raw), str(runtime_path): sha(runtime_raw)}
    cases = []
    for entry in enrollment['sources']:
        folder = root/entry['id']
        folder.mkdir(exist_ok=False)
        original = Path(entry['source']).read_bytes()
        if sha(original) != entry['source_sha256'] or enrollment['inputs'].get(entry['source']) != entry['source_sha256']:
            raise ValueError('source bytes differ from enrolled source identity')
        (folder/'original.pdb').write_bytes(original)
        inputs[str(folder/'original.pdb')] = sha(original)
        row = dict(entry, status='refused_source', chain=None, projection=None)
        try:
            projected, chain = project_first(original)
            path = folder/'input.pdb'
            path.write_bytes(projected)
            inputs[str(path)] = sha(projected)
            row.update(chain=chain, projection=str(path), projection_sha256=sha(projected))
            result = admission(projected)
            row.update(status='admitted', source_inventory=result)
        except (ValueError, UnicodeDecodeError) as error:
            row.update(error_type=type(error).__name__, refusal=str(error))
        cases.append(row)
        write(folder/'admission.json', row)
        inputs[str(folder/'admission.json')] = sha((folder/'admission.json').read_bytes())
    sources = {str(p): sha(p.read_bytes()) for p in Path('src/dvbfixer').rglob('*.py')}
    for name in ('scripts/freeze_methyl_diversity.py', 'scripts/enroll_methyl_diversity.py',
                 'scripts/prospective_independent_cohort.py', 'scripts/prospective_geometry_cohort.py',
                 'scripts/branched_template_baseline.py', 'scripts/research_dunbrack_prior.py'):
        sources[name] = sha(Path(name).read_bytes())
    write(output, dict(schema='dvbfixer.methyl-diversity-source-plan.v1', cases=cases, inputs=inputs,
        source_hashes=sources, enrollment_sha256=sha(raw), runtime=runtime,
        policy=enrollment['policy'], production_changed=False, native_coordinates_read=False))
    verify(json.loads(output.read_bytes()))
    print(json.dumps(dict(enrolled=len(cases), admitted=sum(c['status']=='admitted' for c in cases),
        refused=sum(c['status']!='admitted' for c in cases)), indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--enrollment', type=Path, required=True)
    parser.add_argument('--runtime', type=Path, required=True)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    freeze(args.enrollment, args.runtime, args.root, args.output)
