# Private bounded multi-target completion controls

Status: executable research slice, **not a public replacement backend**. Follows
the [source adapter](protein-completion-source-adapter.md) and
[declared-link admission](protein-completion-link-adapter.md). Production
PDBFixer consumers, geometry thresholds and defaults remain unchanged.

Continuation: [private replay-verified evidence records](protein-completion-publication.md)
now retain accepted and rejected controls using manifest-last commit. They remain
disjoint from public bundles; the earlier source hashes below stay historical.

## Why independent completion is insufficient

Two proposals can each pass against all source observations, yet their **new**
atoms can clash with each other. Writing successful single-target proposals one
at a time also risks treating earlier additions as observed source or exposing a
partial successful output when a later target fails.

`protein_completion/cohort.py` now provides a separate immutable `CohortPlan`,
retained per-target and joint attempts, bounded selection, and an all-or-nothing
in-memory materializer/readback audit. It does not sequentially modify the input.

## Contract and limits

- `plan_cohort(source, exact_targets)` admits **one to four distinct targets**.
  Every target must have missing supported sidechain atoms; no-op targets refuse.
  Targets are canonicalized into source encounter order, not case-normalized.
  Full-source SHA, masks, observations, links and policy bind the plan digest.
- Only explicitly selected missing CB centres are exempted before generation.
  Every target is independently admitted using the existing source rules; other
  missing/undefined or non-L protein Cα centres still refuse. The public
  single-target planner never inherits cohort exemptions.
- Source/link chemistry remains ALA/SER/free-CYS reconstruction with the existing
  canonical-heavy protein/water scope. Linked targets and unsupported chemistry
  still refuse. No broader force-field or multi-chi authority was added.
- This first cohort slice is **staggered-control-only**. It has no empirical
  library parameter or silent authority fallback. Authentic Dunbrack and coupled
  empirical packing remain separate gates.
- Per-target candidates are generated/audited against unchanged source bytes.
  Only individually passing candidates enter the joint Cartesian product.
  At most **81** combinations are enumerated; context/combination budgets fail
  closed, without truncating or claiming unsampled combinations were rejected.
- For each combination, every target is independently re-audited with all other
  targets' additions included as immutable context. A temporary contextual
  digest binds those atoms; original source-bound candidates remain in evidence.
  Both new–source and new–new clashes are therefore checked.
- Select among jointly passing combinations by capped worst nonbonded clearance,
  with stable enumeration as tie break. This is a small feasibility search,
  **not calibrated energy, rotamer probability, native recovery or general packing**.

## All-or-nothing materialization

Requests and candidates are rebound/re-audited before a single insertion pass.
Serials are allocated in output encounter order while respecting source atom/TER
serials. The single-target and cohort adapters share insertion/physical-readback
helpers; there are no intermediate accepted structures.

Readback requires exact original records, exact combined generated identities,
canonical generated records and bounded PDB coordinate rounding. Every target is
then audited again against **rounded additions from all targets**, followed by
the global protein Cα gate. Any failure raises without returning accepted PDB
bytes. Other raw-feasible combinations are not retried after a selected
combination fails rounding; such refusal is retained, not hidden.

The returned `MaterializedCohort` produces a fresh canonical `DatRecord` owning
exactly the generated atoms across all targets, preserving chain case and
insertion codes. No original atom becomes generated. The service writes no files
and implements no public bundle, CLI, registry or filesystem publisher.

## Executed evidence

[Portable record](protein-completion-cohort-evidence.json), raw inputs/plans/
attempts/outputs under `tmp/pdbfixer-research/source-cohort-01/`:

| Declared case | Outcome |
|---|---|
| 1UBQ A1–72 projection: ALA28 + ALA46 CB | Completed and independently read back |
| Same projection: SER20 + SER57 + SER65 CB/OG | Completed and independently read back |
| Same projection: SER20 + ALA28 + ALA46 + SER57 sidechains | Completed and independently read back |
| Synthetic two-SER interacting proposals | Six combinations; five joint passes, one retained clash rejection |
| Synthetic two-ALA proposals | Each passes alone; the only joint combination clashes and is refused |
| Synthetic two-ALA plus boundary water | Raw joint pass; PDB rounding fails one target and the whole output refuses |

Real controls preserve an explicitly documented **ATOM A1–72 projection**, not
full-1UBQ acceptance, occupancy rewriting or native-quality superiority. Original
fixture bytes and pre-execution declarations are retained. The three real
requests generate 2/6/6 atoms; the synthetic success generates four.

Automated tests additionally exercise all **81/81** feasible combinations of four
widely separated synthetic SER targets, source/order/context tampering, missing
and duplicate target selections, resource limits, changed/missing output atoms,
chain `A` versus `a`, insertion codes, single-target equivalence, and retained
individual-target failures. They are software/geometry controls, not empirical
rotamer benchmarks.

Combined focused suite: **288 passed in 12.49 s**; package mypy passed for ten
files and targeted Ruff passed. Eight prior single-target controls were rerun
under `tmp/pdbfixer-research/source-adapter-03/`, retaining separate evidence.
Earlier source hashes remain historical; they are not overwritten after shared
helper extraction.

```bash
PYTHONPATH=src python -m pytest tests/test_protein_completion_cohort.py -q
PYTHONPATH=src python scripts/probe_protein_completion_cohort.py --output tmp/pdbfixer-research/new-cohort-controls --summary tmp/pdbfixer-research/new-cohort-summary.json
```

## Remaining work

Evidence/publication contracts, broader chemistry and multi-target empirical
search need separate review and qualification. Actual licensed rotamer data,
independent scientific comparisons and all nine production consumer adapters
remain open. Do not infer public readiness from these bounded controls.
