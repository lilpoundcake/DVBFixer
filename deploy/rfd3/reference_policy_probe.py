"""Report effective native feature policy on the supplied example before sampling."""
import hashlib
import json
import os
import runpy
from pathlib import Path

import numpy as np
from rfd3.transforms.design_transforms import CreateDesignReferenceFeatures

original = CreateDesignReferenceFeatures.forward


def inspect_features(self, data):
    result = original(self, data)
    a, f = result['atom_array'], result['feats']
    groups = {}
    for role, mask in [('observed', a.is_motif_atom_with_fixed_coord),
                       ('generated', ~a.is_motif_atom_with_fixed_coord)]:
        groups[role] = {'atoms': int(mask.sum()), 'proteins': int(a.is_protein[mask].sum()),
                        'sequence_fixed': int(a.is_motif_atom_with_fixed_seq[mask].sum()),
                        'ref_mask_true': int(f['ref_mask'][mask].sum()),
                        'nonzero_ref_pos': int(np.any(f['ref_pos'][mask] != 0, axis=1).sum())}
    assert [groups[k]['atoms'] for k in ('observed','generated')] == [5036,556]
    import inspect
    path = Path(inspect.getsourcefile(original))
    report = {'native_transform_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
              'effective': {k:getattr(self,k) for k in ('generate_conformers',
                'generate_conformers_for_non_protein_only','provide_reference_conformer_when_unmasked')},
              'groups':groups}
    with Path(os.environ['REFERENCE_PROBE_REPORT']).open('x') as stream:
        json.dump(report,stream,indent=2)
    raise SystemExit(0)


if __name__ == '__main__':
    CreateDesignReferenceFeatures.forward = inspect_features
    runpy.run_module('rfd3.cli', run_name='__main__')
