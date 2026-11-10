"""Retain exact upstream minimizer sources at the revision reported by runtime.

Read-only forensic provenance, not a binary attestation, patched runtime or
license permission for unrelated empirical tables. Upstream notices retained.
"""
from __future__ import annotations

import argparse
import json
import platform
import re
import urllib.request
from pathlib import Path

from scripts.research_dunbrack_prior import sha, write

FILES = (
    'openmmapi/src/LocalEnergyMinimizer.cpp',
    'openmmapi/src/ContextImpl.cpp',
    'platforms/reference/src/ReferenceKernels.cpp',
    'platforms/reference/src/SimTKReference/ReferenceMinimize.cpp',
    'libraries/lbfgs/src/lbfgs.cpp',
)


def capture(output):
    import openmm
    import openmm.version

    revision = openmm.version.git_revision
    if not re.fullmatch('[0-9a-f]{40}', revision):
        raise ValueError('runtime does not report an exact upstream revision')
    output.mkdir(parents=True, exist_ok=False)
    records = []
    for name in FILES:
        url = f'https://raw.githubusercontent.com/openmm/openmm/{revision}/{name}'
        with urllib.request.urlopen(url, timeout=60) as response:
            raw = response.read()
        if len(raw) > 2_000_000 or b'Copyright' not in raw and b'copyright' not in raw:
            raise ValueError('unexpected upstream source payload')
        path = output/name
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open('xb') as stream:
            stream.write(raw)
        records.append(dict(path=name, url=url, sha256=sha(raw), bytes=len(raw)))
    record = dict(schema='dvbfixer.openmm-minimizer-source.v1', files=records,
        runtime_version=openmm.__version__, runtime_full_version=openmm.version.version,
        runtime_declared_git_revision=revision, platform=platform.system(), machine=platform.machine(),
        upstream_notices_retained=True, runtime_binary_attestation=False,
        installed_runtime_modified=False, production_changed=False,
        capture_script_sha256=sha(Path(__file__).read_bytes()))
    write(output/'manifest.json', record)
    print(json.dumps(record, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    capture(args.output)
