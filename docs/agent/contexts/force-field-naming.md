# Force-field Naming

Status: partial

Verified on: 2026-09-18

Verified at commit: `425f290eb85760246766f1d1500e51672c640b2d`

## Purpose

This context owns representation-level residue and atom names required by
AMBER, CHARMM, and GROMACS consumers. It restores explicit protonation variants
lost by toolkit adapters and applies supported target-profile atom shifts.

It does not choose protonation states, parameterize molecules, validate complete
force-field compatibility, or generate a topology.

## Scope

Current production behavior is centered on
`src/dvbfixer/ffutils/ff_names.py::apply_variants_to_pdb_text`. It rewrites PDB
text in place after OpenMM or another pipeline stage has produced user-visible
output.

`dvbfixer rename` is outside this operation. It collapses variant residue names
to canonical PDB parents and does not produce GROMACS atom naming.

## Capabilities

| Capability | Status | Owner | Evidence |
|---|---|---|---|
| Restore explicit AMBER variants after OpenMM normalization | implemented | `ff_names.py::apply_variants_to_pdb_text` | `tests/test_ff_names.py` |
| AMBER/GROMACS protein, cap, and terminal atom shifts | implemented | `ff_names.py::_effective_amber_rename_map` | `tests/test_ff_names.py`, `tests/test_terminal_caps.py` |
| AMBER/GROMACS nucleic-acid sugar hydrogen shifts | partial | `ff_names.py::GROMACS_AMBER_NA_ATOM_RENAMES` | table assertions exist; complete end-to-end coverage does not |
| CHARMM residue variants, cap names, and backbone `H` to `HN` | implemented | `ff_names.py::apply_variants_to_pdb_text` | `tests/test_ff_names.py` |
| Pure text-in/result-out naming service | proposed | API roadmap | not implemented |
| Structured per-change diagnostics | proposed | API roadmap | not implemented |
| Collision rejection | missing | none | current helper does not validate collisions |
| Model-aware terminal classification | missing | none | current classifier groups only by chain |

## Entry Points

| Task | Start symbol | Required contracts |
|---|---|---|
| Change force-field naming | `ff_names.py::apply_variants_to_pdb_text` | `pdb-force-field-naming`, `dat-record-lifecycle` |
| Change protonation variant identity | `variants.py::scan_variant_names` | `dat-record-lifecycle` |
| Add a public naming command | `command_registry.py::COMMAND_REGISTRY` | public command surface in `tasks.toml` |

The complete task record and focused commands are in
[`../tasks.toml`](../tasks.toml).

## Contracts

The authoritative contract records are in
[`../contracts.toml`](../contracts.toml).

`pdb-force-field-naming` currently mutates a path in place. This is unsuitable
as a direct HTTP boundary. A future API must call a pure conversion service and
commit a new artifact only after validation succeeds.

`dat-record-lifecycle` carries explicit variants across pipeline stages using
`chain:resid:icode`. New code must not collapse that identity to two fields.

## Invariants

The relevant records in [`../invariants.toml`](../invariants.toml) are:

- `INV-NAMING-PRESERVE-STRUCTURE`;
- `INV-NAMING-EXPLICIT-VARIANTS`;
- `INV-IDENTITY-ICODE`;
- `INV-NAMING-IDEMPOTENT`;
- `INV-NAMING-NO-COLLISION`.

Only the first two are fully implemented. Do not infer that the current helper
already guarantees all proposed API invariants.

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

- PDB fixed-column parsing/rendering is an adapter concern in `ff_names.py`.
- OpenMM is an upstream adapter that can canonicalize explicit variant names.
- GROMACS force-field files are evidence for target spellings, not domain code.
- `.dat` is the published pipeline language for preserving variant choices.
- The future Node HTTP layer is a transport/workspace adapter and must not
  duplicate naming tables or scientific decisions.

## Side Effects

- The current helper reads and may rewrite the supplied file path in place.
- It writes only when at least one `ATOM` or `HETATM` line changes.
- The write is not an atomic temporary-file replacement.
- It returns a changed-line count, not changed residues or a structured report.
- Validation does not run as a complete pre-commit phase.

## Known Divergences

- Terminal classification is not `MODEL`-aware.
- Legacy two-tuple variant keys remain accepted and can affect insertion-code
  siblings through fallback lookup.
- Target atom-name collisions are not rejected.
- Corresponding `ANISOU` and `TER` labels are not updated.
- Four-character CHARMM residue names are not cleanly idempotent under the
  current parser slices.
- CHARMM `LYN`/`LSN` hydrogen-pair conversion is not safely repeatable.
- Nucleic-acid terminal hydroxyl naming is not fully validated.
- `--atom-naming standard` suppresses shifts; it is not a reverse converter.

## Proposed Work

The proposed API design is in
[`../../plans/api-roadmap.md`](../../plans/api-roadmap.md). Its first code phase
must extract a pure conversion service with typed request/result, collision
validation, explicit model policy, and structured changes. Do not begin by
wrapping the in-place helper in an HTTP handler.

## Focused Verification

```bash
pytest -q tests/test_ff_names.py tests/test_variants_gromacs_lyn.py \
  tests/test_prepare_icode_variants.py tests/test_terminal_caps.py
python scripts/check_agent_docs.py
```
