# Contributing to dvbfixer

Bug reports, documentation improvements and code contributions are welcome from
anyone. Keep changes focused, describe their purpose and include relevant tests.
For a substantial feature or dependency change, discuss the scope in an issue
before implementing it.

## Licensing your contribution

Unless explicitly agreed otherwise, contributions intentionally submitted for
inclusion in dvbfixer's original code or documentation are provided under the
project's [BSD-3-Clause license](LICENSE). You retain your copyright; submitting
a contribution does not assign it to the maintainer. Only submit work you have
the right to license, including any necessary employer or co-author permission.

Preserve existing copyright and license notices. Identify copied/adapted code,
datasets, parameters, images and other third-party material, and provide its
source, exact version and license. Do not assume BSD-3-Clause covers dependencies,
force fields, model weights or fixtures. Contributions to separately licensed
material must respect its applicable terms. See
[third-party notices](THIRD_PARTY_NOTICES/README.md) and the
[licensing audit](docs/research/third-party-license-audit.md).

Do not submit credentials, MODELLER keys, restricted software or model weights.
The project's BSD license does not grant permission to redistribute those items.
This policy governs new submissions; it is not evidence of retroactive consent
for contributions predating the license declaration.

## Development checks

Follow the setup in [README](README.md) and the domain invariants in
[AGENTS.md](AGENTS.md). Python changes should run relevant tests first, then
`pytest -m 'not slow' -q` when the scientific environment is available. GUI changes
should run `npm run typecheck` and `npm test -- --run` from `gui/`.

For license or dependency declaration changes, also run:

```bash
python scripts/gen_license_inventory.py
python scripts/gen_license_inventory.py --check
python -m pytest tests/test_license_inventory.py tests/test_license_evidence.py -q
```

Generating an inventory is not legal clearance of a combined distribution.
