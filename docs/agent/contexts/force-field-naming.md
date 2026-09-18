# Force-field Naming

Status: partial

Verified on: 2026-09-18

Verified at commit: `5478403ee139cb8e4245e84eef1ed3927417e863`

## Purpose

This context owns representation-level residue and atom names required by
AMBER, CHARMM, and GROMACS consumers. It restores explicit protonation variants
lost by toolkit adapters and applies supported target-profile atom shifts.

It does not choose protonation states, parameterize molecules, validate complete
force-field compatibility, or generate a topology.

## Scope

Domain vocabulary and naming policy live in
`src/dvbfixer/domain/force_field_naming.py`. The pure application boundary is
`src/dvbfixer/force_field_naming.py::convert_force_field_naming`: PDB text and
typed options in, transformed text plus a structured report out. Existing
pipelines continue through
`src/dvbfixer/ffutils/ff_names.py::apply_variants_to_pdb_text`, an atomic
in-place compatibility adapter over that service.

`dvbfixer rename` is outside this operation. It collapses variant residue names
to canonical PDB parents and does not produce GROMACS atom naming.

## Capabilities

| Capability | Status | Owner | Evidence |
|---|---|---|---|
| Restore explicit AMBER variants after OpenMM normalization | implemented | `force_field_naming.py::convert_force_field_naming` | `tests/test_force_field_naming.py`, `tests/test_ff_names.py` |
| AMBER/GROMACS protein, cap, and terminal atom shifts | implemented | `domain/force_field_naming.py::atom_name_policy` | `tests/test_force_field_naming.py`, `tests/test_ff_names.py`, `tests/test_terminal_caps.py` |
| AMBER/GROMACS nucleic-acid sugar hydrogen shifts | partial | `domain/force_field_naming.py::GROMACS_AMBER_NA_ATOM_RENAMES` | mappings exist; terminal hydroxyl coverage remains unverified |
| CHARMM residue variants, cap names, and backbone `H` to `HN` | implemented | `domain/force_field_naming.py::atom_name_policy` | `tests/test_force_field_naming.py`, `tests/test_ff_names.py` |
| Pure text-in/result-out naming service | implemented | `force_field_naming.py::convert_force_field_naming` | `tests/test_force_field_naming.py` |
| Structured per-change report and stable failures | implemented | `NamingConversionResult`, `NamingConversionError` | `tests/test_force_field_naming.py` |
| Collision rejection | implemented | `force_field_naming.py::convert_force_field_naming` | collision regressions in pure-service and adapter tests |
| Multi-model V1 policy | implemented | `force_field_naming.py::convert_force_field_naming` | more than one `MODEL` is rejected before rendering |

## Entry Points

| Task | Start symbol | Required contracts |
|---|---|---|
| Change force-field naming policy/service | `force_field_naming.py::convert_force_field_naming` | `pdb-force-field-naming`, `dat-record-lifecycle` |
| Change legacy pipeline file adaptation | `ff_names.py::apply_variants_to_pdb_text` | `pdb-force-field-naming` |
| Change protonation variant identity | `variants.py::scan_variant_names` | `dat-record-lifecycle` |
| Add a public naming command | `command_registry.py::COMMAND_REGISTRY` | public command surface in `tasks.toml` |

The complete task record and focused commands are in
[`../tasks.toml`](../tasks.toml).

## Contracts

The authoritative contract records are in
[`../contracts.toml`](../contracts.toml).

`pdb-force-field-naming` is pure and suitable for a future CLI/HTTP adapter.
The legacy wrapper still mutates a path for existing pipeline callers, but it
validates first, performs no write on a no-op or failure, and commits changed
text with a same-directory atomic replacement. An HTTP adapter must call the
pure service through the planned command/application boundary and create a new
workspace artifact rather than use the in-place wrapper.

`dat-record-lifecycle` carries explicit variants across pipeline stages using
`chain:resid:icode`. New code must not collapse that identity to two fields.

## Invariants

The relevant records in [`../invariants.toml`](../invariants.toml) are:

- `INV-NAMING-PRESERVE-STRUCTURE`;
- `INV-NAMING-EXPLICIT-VARIANTS`;
- `INV-IDENTITY-ICODE`;
- `INV-NAMING-IDEMPOTENT`;
- `INV-NAMING-NO-COLLISION`.
- `INV-NAMING-MODEL-POLICY`.

The naming-specific invariants are implemented by the pure service. The wider
insertion-code invariant remains partial because unrelated pipeline-local maps
have not all migrated to typed three-field identities.

## Callers

- `prepare/pipeline.py` applies the helper near final output on the legacy
  backend.
- `minimize/pipeline.py` restores variants and target naming after writing the
  minimized structure.
- `protonate.py` applies it on its legacy output paths.
- `zbs.py` forwards atom-naming options to prepare and minimize rather than
  owning naming rules.

The `tleap-reduce` early-return paths do not consistently apply the final naming
helper. Treat this as a known backend discrepancy.

## Adapters

- PDB fixed-column parsing/rendering is contained in the pure application
  service; domain policy does not depend on file I/O or external toolkits.
- Atomic path replacement and legacy key coercion are adapter concerns in
  `ffutils/ff_names.py`.
- OpenMM is an upstream adapter that can canonicalize explicit variant names.
- GROMACS force-field files are evidence for target spellings, not domain code.
- `.dat` is the published pipeline language for preserving variant choices.
- The future Node HTTP layer is a transport/workspace adapter and must not
  duplicate naming tables or scientific decisions.

## Side Effects

- The pure service has no side effects.
- Its result includes transformed text, per-coordinate changes, summary counts,
  diagnostics, model identity, and stable rule/error codes.
- The compatibility wrapper reads and may atomically replace its supplied path.
- The wrapper retains the historical changed-`ATOM`/`HETATM` integer return.

## Known Divergences

- V1 rejects multiple `MODEL` blocks rather than converting them independently.
- Legacy two-tuple keys remain accepted only by the compatibility adapter and
  mean an exact blank insertion code, never an insertion-code wildcard.
- Nucleic-acid terminal hydroxyl naming is not fully validated.
- `--atom-naming standard` suppresses shifts; it is not a reverse converter.
- The `tleap-reduce` early-return paths still bypass the final naming adapter.
- No dedicated public CLI or HTTP adapter exists yet.

## Proposed Work

The API design is in [`../../plans/api-roadmap.md`](../../plans/api-roadmap.md).
Phase 1's pure typed boundary and safety validation are implemented. The next
slice is a dedicated non-destructive CLI adapter with JSON reporting, followed
by the workspace-scoped HTTP route. Do not expose the in-place compatibility
wrapper as an HTTP handler.

## Focused Verification

```bash
pytest -q tests/test_force_field_naming.py tests/test_ff_names.py \
  tests/test_variants_gromacs_lyn.py tests/test_prepare_icode_variants.py \
  tests/test_terminal_caps.py tests/test_scientific_domain.py
python scripts/check_agent_docs.py
```
