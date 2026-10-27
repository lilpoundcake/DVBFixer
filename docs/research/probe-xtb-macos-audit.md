# Probe / xTB: macOS audit, 2026-10-28

Subsequent execution: [Linux continuation](probe-xtb-linux-audit.md) reproduces
the failed acceptance gates with fresh Linux baselines and records additional
1VCU results. Linux deferral below describes the original macOS session.

**Outcome: bounded macOS execution complete; production scientific acceptance
failed. Linux execution is deferred by the user, not passed.** Neither production
adapters nor installed scientific dependencies nor packaging were changed.
Recommendations below require a separate decision before implementation.

## Evidence and scope

Branch: `test/probe-xtb`. Primary evidence is
`.artifacts/probe-xtb/complete-macos/`: 24 tasks returned zero, fresh preparation
and CPU OpenMM baselines, complete native stdout/stderr, binary fingerprints,
input digests, extraction maps and script snapshots. Its report explicitly says
`execution_complete=true`, `no_xtb_discovery_verified=true`,
`adapter_contracts_met=false`. Successful audit-script execution is not successful
optimization. The guard covers preparation and ordinary minimization without an
xTB request; it is not a global proof about every CLI path.

Supplemental evidence: `geometry-final/`, `separate-{fixed,installed}-{0,1}/`,
`results-final.json`, `package-provenance/`, upstream unit logs, native contract
fixtures and `contracts-without-xfail.xml`. These supersede the partial conclusions
in `results.json` and `results-stage3.json`; earlier pilot logs remain historical.
The transferable archive includes these items and the source repository snapshot.
See [Linux handoff](../plans/probe-xtb-linux-handoff.md).

The final portable extended runner was also replayed on those frozen macOS
baselines: **41/41 tasks returned zero**, including all component/geometry arms.
This is a replay-tool check, not another independent preparation. Its report
correctly leaves the fresh-preparation discovery gate false.

This implements the Probe/xTB investigation in B.12–B.13 of
`research/9cqd-protpardelle-multichain:docs/plans/environment-dependency-audit.md`.
It does **not** implement the unrelated chemical-stack or packaging migration.
Inputs are tracked, checksummed ASN, glycosylated hinge/CH3, 7X35/PLM and 1UBQ,
plus generated seeded controls. All development minimizations used `--no-solvent`.
ALPB-water is implicit solvent in isolated xTB calculations, not a solvent box.

## Builds and mechanism controls

Host: macOS 26.6.2 arm64, scientific Python 3.12.13. One numerical thread.

| Component | Exact identity / observation |
|---|---|
| Installed xTB | 6.7.1 (`edcfbbe`), conda-forge `gfortran_hc17bbfb_5` |
| Package SHA-256 | `df7623a53b06a503a377f4c021471991129d107a3047b688946fe7ae8bcf626b` |
| Recipe revision | `811e76f7850e403de1b564d2dd449fe0a0d8d9cb`; bundled patches do not backport #1429 |
| Fixed upstream xTB | `a45f05d397268fd683b5f9f824ed59912c5490b5` (#1429 merge) |
| Successful build | Meson 1.7.1, Ninja 1.13.0, GNU Fortran 16.2.0, Netlib BLAS/LAPACK, OpenMP off |
| Upstream unit suite | 30/30 suites passed; optional tblite/CPCM-X disabled, not qualified |
| Probe | 2.18, `091d117f5bceef87389002dccfe26df78516fdda`, C/Make build |

Version strings alone do not establish backports. Current release checks and
recipe inspection found the installed v6.7.1 line, not a released #1429 build.
The comparison is package-versus-pinned-upstream, not a one-change causal proof:
other upstream changes/build differences also exist.

The 540-atom separated-ethanol and 542-atom connected-alkane controls actually
activate fragmented diagonalization under the default optimizer (logged).
Fixed upstream converges with max anchor drift below `1.4e-11 Å`; installed
default fails before useful optimization. Installed `rf` reaches partial
trajectories with drift up to 6.17 Å, without convergence. Smaller installed
controls converge yet drift 0.041–0.165 Å. The predeclared drift tolerance is
`1e-5 Å`; near-zero movement in an initial frame is not a pass.

CMake/Accelerate numerical failures are **not** evidence against #1429. The same
CMake object build relinked to Netlib converged on the neutral alkane and PLM;
Meson/Netlib passed upstream suites. Linear algebra configuration matters.
Explicit small-system L-ANC trials also failed/drifted; no universal optimizer
workaround is established.

## Matched fresh complex results

Clashes below use the same Python detector and coordinates before/after each
arm; they are not Probe-vs-Python quality rankings.

| Arm | Convergence / geometry after frozen-coordinate restoration |
|---|---|
| Installed GFN-FF, glycan, default and `rf` | Converges at cycle 94; raw frozen max/RMS drift 1.244/0.810 Å; clashes **34 → 49** |
| Fixed GFN-FF, glycan, 200-cycle request | No convergence; partial frame only, not accepted |
| Fixed GFN-FF, glycan, 1000-cycle request | Default/`rf` converge at cycle 582; raw max drift `8.1e-10 Å`; clashes **34 → 66** |
| Installed GFN-FF, PLM, production arguments | Converges; heavy-atom RMS/max displacement 1.041/2.422 Å; **17/29** original 4 Å contacts lost; clashes **56 → 62** |
| Installed GFN2, charge −1, UHF 0, vacuum | Converges; RMS/max 0.717/1.978 Å; **12/29** contacts lost; clashes **56 → 63** |
| Installed GFN2, charge −1, UHF 0, ALPB water | Converges; RMS/max 0.647/1.395 Å; **11/29** contacts lost; clashes **56 → 62** |

GFN-FF charge/spin/water sensitivity and all fixed-build PLM arms were also run.
Fixed PLM arms lose 17–18/29 contacts and end with 60–62 clashes. These geometric
contacts are not binding affinity. Energies from unlike methods are not ranked.
GFN2 costs roughly 3 seconds on this 49-atom ligand, versus roughly 1–2 seconds
for GFN-FF; concurrent timings are observations, not calibrated benchmarks or
evidence of scalability/support for arbitrary complexes.

The installed glycan's two frozen/mobile boundary bonds change from
1.516/1.517 Å to **1.657/1.642 Å** after anchor restoration, despite raw xTB
values near 1.50 Å. Restoring anchors conceals drift, not its geometric effect.
The fixed build keeps these bonds near 1.494 Å.

Geometry-delta analysis measures all topology bonds/angles and ordered carbon
tetrahedral volumes. Installed glycan: two bonds change >0.1 Å and four angles
>10°. Fixed extended glycan and PLM arms: none exceed those relative thresholds;
no candidate tetrahedral sign flips were observed. These are change screens, not
absolute ideal-geometry/CIP validation. Whole-system RMS is diluted by unchanged
protein atoms. Chemistry diagnostics retain pre-existing hinge amide/angle/
disulfide findings; this is not an initially perfect structure. No broad
stereochemical or glycan planarity qualification is claimed.

The two hinge components (247 atoms each) were optimized separately in both
binaries. Fixed components converge at cycles 422 and 380 with negligible anchor
drift. Recombining both separate results gives **67 clashes**, versus 66 jointly
and 34 initially. Joint-versus-separate mobile coordinates differ by RMS 1.57 Å,
maximum 7.65 Å. This demonstrates sensitivity to grouping/optimization basin,
not that coupling alone caused every change. Separation is not a pocket model.
Both approaches omit the surrounding protein; a context-aware alternative needs
its own boundary/capping/electrostatics design and acceptance study.

1VCU strict ligand parameterization fails on EPE A/382 with N explicit valence 4.
That blocker is retained in `1vcu-openmm.log`; no ligand stripping or protein-only
substitution was used to manufacture a successful comparator.

## Probe versus Python

Identical explicit hydrogens/coordinates were used in both radius arms. The
production command supplies `ALL ALL` without a two-selection mode and Probe
returns `ERROR: could not open file 1: ALL`. Its condensed parser also reads
field 5 (dot count) instead of field 6 (minimum gap). Direct research invocations
use `-self ALL`, preserving this distinction from the unchanged production path.

Ground-truth controls establish:

- Both engines: clean separated methanes and ASN yield no clashes; inter-methane
  clashes are found; the bonded ethane C–C contact is excluded.
- Probe adds useful coverage: octane C1–C5 at **2.62675 Å**, separated by four
  graph bonds, is reported; Python misses it because it excludes same-residue
  pairs. Control selection uses geometry/connectivity, not detector outcomes.
- Methylamine/acetone H–O is 1.8 Å at 180° or 90°. Both controls retain Probe
  `hb` classification and Python excludes the target pair. This does **not**
  demonstrate an angular-selectivity advantage. The bent control has other
  clashes and is not an otherwise-clean negative.

Real-input discordant pairs include resolved identities, distances, Probe contact
classes and Python exclusions. Raw count differences do not measure precision:
Probe radii/classification differ, and one pair may have both `hb` and `so` dot
classes. Summary counts retain `bo`/`so`; do not treat every such pair as an
independently adjudicated false positive/true clash. Atom-field alignment alone
changed an earlier glycoprotein count from 3478 to 47 without moving coordinates.
Preserve format/radius controls rather than interpreting that as improvement.
Probe runtime on fresh hinge/PLM was about 0.22/0.58 seconds versus Python
0.08/0.23 seconds; building Probe is small, but integration validation is real work.

## Coverage and contract disposition

| Required area | Disposition |
|---|---|
| Actual native failure/timeout/missing executable | Exercised; native fixtures bounded and validated |
| Probe malformed stdout / valid overlap | **Fail**: malformed becomes clean; valid overlap is missed |
| xTB count/elements/finite coordinates/convergence | **Fail**: invalid values accepted; truncation raises instead of safe skip |
| Atom ordering | Element sequence checked by research reader; same-element permutations are not detectable in XYZ alone |
| >5000 / no heterogen / boundary 5000 | Real-topology control tests verify switch and no native call for empty subsystem; not a 5001-atom QM benchmark |
| Engine/status/fallback visibility | Findings carry engine metadata; clean Probe results and fallback lack complete requested/actual/status provenance. xTB skip warnings exist, but empty subsystem can be silent |
| Joint/separate heterogens, pocket loss, interface | Matched comparisons completed; scientific gate fails |
| Python regression suite | **1493 passed, 2 skipped, 31 deselected, 7 strict xfailed**; `--runxfail` gives seven failures |
| Linux x86_64 | Deliberately deferred, replay required |

The two skips require unavailable `biotite.structure` and Torch. An initial
repeat had six subprocess failures because the scientific environment's `bin`
was absent from PATH; correcting only PATH produced the result above. Both logs
are retained. Focused audit/helper tests: 27 passed, seven strict xfailed. Ruff,
strict-target mypy, agent-doc validation, CLI reference and GUI schema sync passed.

## Recommendation, not implementation

**Probe:** retain as optional; a demonstrated intramolecular benefit rules out
removal merely because Python exists. Repair command/parser/format validation and
explicit engine/fallback reporting before trusting automatic selection. Do not
make it an obligatory installation dependency on this study.

**xTB:** current integration is not accepted as a safe whole-complex polishing
step. A package upgrade/#1429 alone is insufficient: fixation improves but pocket
clashes worsen. Keep available pending review, propose a separately qualified
opt-in profile with strict output/convergence/identity validation, raw drift and
post-merge geometry gates, explicit skipped/failed status and retained evidence.
GFN2 must be a separate charge/spin/solvent-aware mode, not an automatic substitute.

No package removal, default change, release bump or production repair is authorized
by these results. General benefit, broader clinical/structural cohorts and a
context-aware optimizer remain future research, not unfinished macOS executions.
