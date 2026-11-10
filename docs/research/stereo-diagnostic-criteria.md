# Warning-only stereocentre diagnostics

The `diagnose` chemistry family adds two **detection-only WARNING categories**:
`ca_planarity` and `thr_beta_chirality`. Neither changes coordinates, reconstructs
atoms, assigns general CIP stereochemistry, or changes `prepare`/`minimize`.
The existing D-Cα detector and its ERROR semantics are unchanged. See the
[command guide](../commands/diagnose.md) for reporting and filtering.

## Measurement and conservative criteria

For a centre X and three ordered neighbors A, B, C, measure
`v = ((A-X) × (B-X)) · (C-X)` in Å³ (a scalar triple product, **not** volume
divided by six). The signed out-of-plane angle is
`asin(v / (|(A-X) × (B-X)| |C-X|))` in degrees. Unlike an absolute volume gate,
this angle does not depend on bond lengths. Ideal tetrahedral three-neighbor
geometry has an absolute angle of approximately **54.74°**.

| Category | Ordered atoms X, A, B, C | WARNING criterion |
|---|---|---|
| `ca_planarity` | CA, N, C, CB | Absolute angle **≤10°**, of either sign |
| `thr_beta_chirality` | CB, CA, OG1, CG2 | Signed angle **<−10°**; canonical THR R requires positive |

Ten degrees is a deliberately narrow severe-flattening screen, far from ideal
tetrahedral geometry, **not** an empirically optimized acceptance boundary or a
residue-specific force-field reference. Nearly flat THR-Cβ is deliberately not
assigned an inversion from a small negative sign. This check does not cover every
possible malformed THR geometry. There is no VAL-Cβ chirality finding: its methyl
branches are chemically equivalent. ILE, modified THR such as TPO, ligands and
arbitrary residue names do not receive canonical THR-Cβ assignments.

Cα screening recognizes canonical amino acids, the shared AMBER/CHARMM
variant-to-parent mapping, and diagnose's existing nonstandard-amino-acid whitelist;
GLY is excluded. Each measurement requires one unambiguous residue identity and
one correctly named, correctly element-typed atom per required name. Nonfinite
coordinates, coincident vectors (length ≤1e−12 Å), and ill-defined first/second
planes (sine <0.1) are skipped rather than assigned handedness. These are checks
of the loaded topology's selected coordinates, not alternative-conformer or
occupancy validation; structural checks remain separate. Missing topology bonds
do not disable the name-based stereo screen.

This policy does **not** promote the research 0.5 Å³ gate into production. A
positive near-planar Cα can warn without a D assignment. A negative near-planar
Cα can have both a shape WARNING and the unchanged D-Cα ERROR, where the existing
detector's own deadband and plane checks allow it. No existing error is suppressed.

## Validation evidence (2026-11-10)

Portable tests in `tests/test_diagnose_stereocentres.py` cover analytic signed
controls, ideal tetrahedra, rigid transforms, scaling, reflection, both sides of
the 10° boundary, zero/nonfinite/ill-defined geometry, missing/duplicate atoms and
residues, element mismatches, recognized variants, unsupported chemistry, exact
case-sensitive chains/insertion codes, text/JSON/severity/scope behavior, and
unchanged input bytes. A reflected branch from native 1UBQ THR is detected while
its Cα coordinates and existing chirality findings stay unchanged.

Existing [tracked fixtures](../../tests/fixtures/README.md), without cleaning or
rebuilding, give the following bounded native controls:

| Fixture | Measurable Cα | THR-Cβ | Minimum absolute Cα angle | Minimum THR-Cβ angle | New warnings |
|---|---:|---:|---:|---:|---:|
| 1UBQ | 70 | 7 | 40.242° | 51.000° | 0 |
| 1A3N | 532 | 32 | 44.758° | 49.006° | 0 |
| 1TIM | 440 | 20 | 34.389° | 49.620° | 0 |

These **1,042 Cα and 59 THR-Cβ** measurements validate sign convention and provide
a small no-new-warning control, not a population false-positive rate. Fixture
provenance/checksums remain in their existing registry; no new structural assets
or generated research outputs are required by the tests.

A separate read-only sweep loaded all **32 PDB fixtures** in this checkout and
reported zero new stereo warnings, with unchanged input bytes and no load failures.
That broader development smoke is not an independent sensitivity/specificity
benchmark, and does not imply that the fixtures are clean in other categories.

Separate manual readback of four retained research endpoints, using only the
production detector (no research-module imports), found:

| Retained endpoint | Finding | Signed angle in final PDB |
|---|---|---:|
| Legacy 1CRN THR2 | THR-Cβ opposite canonical R | −51.973° |
| Legacy 1ENH THR6 | THR-Cβ opposite canonical R | −48.392° |
| 1SHG THR24 seed-1 raw | Target Cα near-planarity | +0.196° |
| 1SHG THR24 seed-1 full-generated | Observed-context GLN16 Cα near-planarity | +8.910° |

Inputs were byte-identical after readback. Rounded final-PDB GLN16 triple product
is +0.479232 Å³; its positive sign is not a D-Cα finding. These are retrospective
smokes, not independent threshold qualification or proof of inversion causality.
Private endpoint files are deliberately **not** portable test dependencies.
No coordinate repair, empirical rotamer tables, PyRosetta dependency, scientific
backend/platform changes or public completion capability is introduced.

## Reproduce the portable controls

```bash
pytest -q tests/test_diagnose_stereocentres.py
dvbfixer diagnose tests/fixtures/protein_context/1UBQ.pdb \
  --only chemistry --format json
```

Other existing chemistry categories can still report findings. A clean new
stereo screen is not a guarantee that a structure is physically valid or safe
for a downstream workflow; warnings request source/identity review, not an
automatic repair command.

## Unified CLI follow-up

Subprocess testing after the first implementation found a presentation bug that
direct `diagnose.main()` tests could not expose: fd-level runtime annotation
inserted `!!! WARNING !!!` inside a JSON message containing the word "warning".
The unified CLI now opts single-input diagnose JSON into verbatim stdout capture.
JSON fields are data, not free-text diagnostics; stderr remains fully captured,
emphasized and summarized, including inherited child-process stderr and exception
cleanup. Other commands and human-readable batch progress retain normal handling.

The same tests exposed false runtime ERROR markers on `D-isomer error: NO` and
the count-table row `ERROR: 0 findings across 0 residues`. Exact negative status
and count-table forms are no longer classified as extra diagnostics. Actual
per-residue warnings/errors and positive D-isomer status remain visible.

The additional controls exercise public text/JSON commands, ERROR filtering,
Unicode log content, exit codes, mixed warning/error batch report publication,
format-option spelling/last-value semantics, child stdout beyond the pipe read
chunk size without a final newline, real child stderr and exceptional exit.
None of these changes affects scientific geometry thresholds or coordinates.

Initial branch verification at `039e1ba`: **2,231 fast tests passed**, three optional Torch/
Biotite skips and 33 slow tests deselected. This includes 61 stereocentre controls.
Targeted Ruff lint/format, generated CLI/GUI checks and agent documentation maps
pass. The pure geometry module passes mypy; a broader four-module diagnose mypy
check reports the same six pre-existing errors reproduced from `main`, with no
new errors in this extension. This is not full slow-suite/platform acceptance.

Follow-up verification: **2,245 fast tests passed**, the same three optional skips
and 33 slow deselections; **164 diagnose/runtime tests passed** in the focused
run, including its glycoprotein slow smoke. This includes 70 stereocentre tests
and nine CLI/runtime tests. Full-source Ruff, CLI/runtime/geometry mypy,
generated CLI/GUI checks, agent maps, version synchronization and license-inventory
freshness pass. Inventory freshness does not establish legal clearance. Retained
before/after readback against `039e1ba` reproduces both presentation defects and
confirms their removal without scientific execution.
