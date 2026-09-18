# Structure Diagnostics

Status: implemented

Verified on: 2026-09-18

Verified at commit: `425f290eb85760246766f1d1500e51672c640b2d`

## Purpose

This context owns read-only structure-quality analysis and deterministic text or
JSON reports. It identifies structural, chemistry, and steric findings; it does
not repair coordinates, rename residues, rewrite connectivity, or decide whether
a pipeline may publish a structure.

Runtime warning/error highlighting is a neighboring CLI concern. It observes
process output but does not define scientific findings or report serialization.

## Scope

Production behavior starts at `src/dvbfixer/diagnose/pipeline.py::main` and uses
the `Finding` model in `src/dvbfixer/diagnose/report.py`. Structural checks combine
raw fixed-column PDB inspection, OpenMM topology/positions, PDBFixer discovery,
and shared geometry and duplicate-chain helpers.

The unified CLI normalizes CIF once to temporary PDB before calling diagnose.
Direct calls to `dvbfixer.diagnose.main` receive the PDB-oriented implementation.

## Capabilities

| Capability | Status | Owner | Evidence |
|---|---|---|---|
| Missing pieces, atom placement, chain breaks, duplicate frames | implemented | `diagnose/structural.py::run_all` | `tests/test_diagnose_pipeline.py`, `tests/test_duplicate_chain_coordinates.py` |
| Valence, bond geometry, peptide geometry, chirality, disulfides | implemented | `diagnose/chemistry.py::run_all` | `tests/test_diagnose_chemistry.py` |
| Probe-first clash checks with pure-Python fallback | implemented | `diagnose/steric.py::run_all` | `tests/test_diagnose_steric.py` |
| Deterministic text and machine-readable JSON reports | implemented | `diagnose/report.py`, `diagnose/pipeline.py` | `tests/test_diagnose_report.py`, `tests/test_diagnose_pipeline.py` |
| Repair or mutation of the diagnosed structure | missing by design | none | report-only CLI contract |

## Entry Points

| Task | Start symbol | Primary tests |
|---|---|---|
| Change orchestration, filtering, output, or status | `diagnose/pipeline.py::main` | `tests/test_diagnose_pipeline.py` |
| Add or change a finding category | `diagnose/structural.py::run_all`, `chemistry.py::run_all`, `steric.py::run_all` | focused diagnose module tests |
| Change report shape or ordering | `diagnose/report.py::Finding`, `format_report`, `findings_to_dict_list` | `tests/test_diagnose_report.py` |
| Change duplicate-chain policy | `pdbutils/duplicates.py::duplicate_protein_chain_coordinates` | `tests/test_duplicate_chain_coordinates.py` |
| Change CLI output capture or summary | `runtime.py::tee_output` | `tests/test_cli_runtime.py` |

## Contracts

`Finding` carries severity, category, case-sensitive chain, residue identifier,
residue name, message, optional atom and fix hint, and optional details. The
residue identifier includes its insertion code, such as `100A`.

Severity filtering occurs before serialization and exit-status selection. Text
and JSON findings sort by severity, chain, numeric residue, insertion code, and
atom. JSON exposes input/count metadata, findings, severity counts, and chirality
audit data; user-facing Unicode is emitted with `ensure_ascii=False`.

Exit statuses are part of the public contract:

- `0`: no ERROR finding remains after the selected checks and severity filter;
- `1`: at least one reported ERROR finding remains; this is a completed analysis,
  not an execution failure;
- `2`: argument, missing-input, or load failure.

## Invariants

- Diagnose is report-only: input structure bytes and coordinates are unchanged.
- A clean run emits a report but no empty fd-level diagnostic-summary banner.
- Chain IDs remain exact and case-sensitive; chains `D` and `d` are distinct.
- Insertion codes remain part of finding residue identity.
- Explicit PDB chain transitions are boundaries, never internal chain breaks.
- Multi-model input is analyzed as MODEL 1 only and receives a WARNING finding.
- Coordinate-identical complete protein chains are suspicious but non-fatal.
- JSON retains readable Unicode in paths, units, arrows, and messages.

## Callers

- `src/dvbfixer/cli.py::main` applies CIF normalization, runtime headers, batch
  dispatch, fd-level teeing, and optional log capture around diagnose.
- `src/dvbfixer/batch.py::run_directory` treats status `1` as findings, writes a
  report per structure, continues by default, and distinguishes execution
  failures in its aggregate summary.
- `src/dvbfixer/zbs.py::_run_pipeline` writes postflight JSON. ERROR findings warn
  by default and fail only under `--strict-postflight`; statuses above `1` fail.
- Humans and CI consume standalone text or JSON reports directly.

## Adapters

- OpenMM supplies topology, bonds, positions, and chain/residue objects.
- PDBFixer discovers missing residues, atoms, and terminal atoms without adding
  them.
- Raw PDB parsing preserves records not represented reliably by OpenMM, including
  MODEL, SEQRES, altLoc, insertion-code, and audit REMARK information.
- MolProbity `probe` is optional; failure falls back to SciPy-based clash checks.
- `structure_input.py` is the sole CIF-to-PDB boundary and preserves compatible
  single-character chain IDs while mapping incompatible ones.

## Side Effects

- `-o/--output` creates or replaces the requested report; otherwise the report is
  printed to stdout.
- Multi-model analysis creates and removes a temporary MODEL-1 PDB directory.
- Probe may run as a subprocess. No check writes to the diagnosed structure.
- Verbose family timings go to stderr and are not part of the report schema.
- The unified CLI may append captured stdout/stderr to `--log-file`.

## Known Divergences

- The text summary counts affected residues by `(chain, resid)` rather than a
  richer atom/residue key; JSON retains each complete finding separately.
- Duplicate-chain detection intentionally ignores residue numbers. It requires
  equal protein-atom counts and matching residue names, atom names, file order,
  and coordinates within PDB precision.
- Multi-model files do not receive per-frame analysis; only MODEL 1 is checked.
- Direct module invocation bypasses unified CIF normalization and fd-level runtime
  capture.
- Some structural helper failures are suppressed to keep reporting resilient, so
  absence of a category is not proof that every optional heuristic ran.

## Proposed Work

No repair behavior belongs in this context. Future service/API work may extract a
typed analysis request/result that accepts a normalized structure artifact and
returns findings without `sys.exit`, while preserving the current CLI report and
status contract. Such work must not duplicate scientific checks in Node or a GUI.

## Focused Verification

```bash
pytest -q tests/test_diagnose_report.py tests/test_diagnose_pipeline.py \
  tests/test_diagnose_chemistry.py tests/test_diagnose_steric.py \
  tests/test_duplicate_chain_coordinates.py tests/test_cli_runtime.py \
  tests/test_batch.py tests/test_zbs_postflight.py
python scripts/check_agent_docs.py
```

The fd-level integration boundary is deliberate: `runtime.tee_output` captures
Python and inherited child stdout/stderr, emphasizes recognized WARNING/ERROR
lines, deduplicates them, and emits a summary during normal or exceptional
cleanup. It must not parse report files, alter JSON, invent `Finding` objects, or
change diagnose exit status. Reports written with `-o` therefore remain the
authoritative diagnostic artifact even when no report lines cross the captured
file descriptors.
