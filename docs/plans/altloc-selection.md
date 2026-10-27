# Alternate-location selection

## Agreed behavior

- Exactly one labeled state per residue; shared blank-altloc atoms survive.
- Manual, occupancy, first-encountered, preferred-label, and automatic plus
  explicit per-residue overrides all use one Python implementation.
- GUI checkboxes are mutually exclusive per residue. Preview is not publication.
- Write a new structure with cleared altloc labels, never change the source.
- Preserve coordinates and occupancy; reconcile ANISOU, CONECT and LINK without
  inventing bonds. Refuse ambiguous identities or conflicting shared atoms.
- Python owns scientific selection; GUI supplies a versioned, source-digest-bound
  selection plan and publishes results only in the owning workspace.

## Implementation sequence

- [x] Python analysis/selection service, deterministic policies, explicit plan.
- [x] CLI analysis, dry-run, report and separate output; registry/reference/schema.
- [x] Tests: identities, policies, ties, missing preferred labels, records, input safety.
- [x] GUI analysis/review with mutually exclusive checkboxes and Mol* previews.
- [x] Workspace-scoped execution/publication with provenance and stale-input protection.
- [x] Early diffusion diagnostics directing users to resolve altlocs first.
- [x] User documentation and focused/full available verification.

## Policy details

Occupancy is the mean occupancy of labeled heavy atoms, not an atomwise mix.
Ties resolve deterministically with a visible warning; missing occupancies and
different atom inventories are reported. First means input encounter order, not
alphabetical A. Preferred-label absence defaults to an error unless an explicit
fallback is selected. Manual requires one selection for every alternate residue.
An override always wins. States with different residue names are selected whole.
Keys retain model ordinal, chain case, residue number, insertion code and record
occurrence. No implicit cross-residue correlation is inferred from matching labels.

## Validation boundary

Selecting a state is not chemistry repair, parameterization or optimization.
Preview must display the exact Python-selected atom identities. No automatic
altloc cleanup inside diffusion; unresolved alternatives must fail before the
external sampler is launched, with actionable diagnostics.

## Progress and evidence

Implemented in `src/dvbfixer/altloc_selection.py` and `altlocs.py`, with the
workspace adapter in `gui/server/altloc-api.ts`, the checkbox editor in
`gui/src/components/AltlocPanel.tsx` and isolated Mol* viewer in
`gui/src/components/AltlocPreview.tsx`. No main-viewer state is replaced by preview.
The report uses model-specific serial masks from Python, not JS atom selection
heuristics. The endpoint is internal GUI/workspace API, not a new V1 OpenAPI claim.

Verification:

- Python policy/CLI tests cover all modes, manual overrides, source digest,
  insertion codes, case-distinct chains, MODEL/occurrence identity, microheterogeneity,
  missing occupancy, ties, partial inventories, ANISOU, LINK, CONECT, dry run and
  no-replace publication. Tracked 1A3N/1DQJ/2VLN/2VLQ selections preserve coordinates.
- Mol* parser/query test verifies both A/B are retained and the exact Python serial
  mask selects B, without requiring WebGL. UI state tests verify mutual exclusion.
- API tests cover analysis without registration, reviewed apply, stale inputs,
  revoked access, malformed outputs, symlinks/hardlinks, manifest rebasing and
  failure cleanup. GUI suite: 349 passing tests.
- Real Node adapter → installed CLI → report validation → reviewed apply passed
  on tracked 2VLQ: six alternate residues, 27 removed atoms; source digest unchanged.
- GUI typecheck, ESLint, build and standalone server smoke pass. Existing bundle-size
  and mixed static/dynamic import warnings remain non-failing.
- Ruff, strict API-surface mypy, documentation, agent maps, release metadata,
  generated CLI reference and GUI-spec checks pass.
- Fast Python suite initially passed 1,373 tests with two optional skips and
  31 slow tests deselected. Final rerun after additional regressions is pending.

Visual browser/WebGL acceptance remains a manual follow-up: the data/query path
is tested, but no screenshot or interactive browser result is claimed. A real ML
inference run is not part of these checks. Review drafts are transient; closing
the editor requires re-analysis, whereas applied plans/reports are retained.

Earlier diffusion warning/fixed-link changes in this working tree are separate
work and have not been discarded. Nothing has been committed or pushed in this
implementation task.
