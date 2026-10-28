# Probe / xTB: Linux continuation, 2026-10-28

Completed continuation: [Probe adapter, full-context and restraint-ablation
study with final decisions](probe-xtb-followup.md).

**The fresh Linux audit completed its 46-task execution matrix. Scientific
acceptance still fails.** This continues the [macOS audit](probe-xtb-macos-audit.md)
and [Linux handoff](../plans/probe-xtb-linux-handoff.md) on branch
`test/probe-xtb`, starting at `5566065fb6f3bd5517110b8e57ae866fa8a4e720`.

## Evidence and environment

Evidence root: `.artifacts/probe-xtb/` (ignored by git).

| Evidence | Contents |
|---|---|
| `linux-fresh/` | Fresh preparation, CPU OpenMM baselines, 46 tasks, native logs, extraction maps, geometry reports and checksums |
| `linux-summary.json` | Original bounded matrix index, including the joint-glycan timeout |
| `linux-summary-longtimeout.json` | Original matrix plus explicitly selected supplemental joint/component comparison |
| `linux-longtimeout/` | Supplemental joint-glycan calculation with a 300-second native limit |
| `linux-1vcu/` | Additional strict comparator, empty-cache GAFF2 run, atom-identity/chirality checks, Probe and xTB comparisons |
| `linux-provenance/` | Package/build inventories, native digests, `ldd`, source commits, downloaded xTB package and its embedded recipe |
| `linux-build/` | Pinned source checkouts, native build, compiler/build configuration and upstream unit logs |
| `linux-fast-suite.{log,xml}` | Full fast-lane regression results |
| `linux-contracts-without-xfail.{log,xml}` | Seven explicitly failing adapter contracts |

Host: Linux x86_64, kernel `5.15.0-102-generic`, glibc 2.39. Scientific
environment: existing local `dvbfixer` environment, Python 3.13.15, NumPy 2.4.6,
OpenMM 8.6.1, PDBFixer 1.12.0, RDKit 2026.3.1, PROPKA 3.5.1. This is independent
Linux **coordinate preparation**, not a fresh solve of `environment.yml` or an
identical scientific dependency stack to macOS.

The environment's editable installation initially resolved a different worktree.
All retained scientific runs and test reports explicitly use
`PYTHONPATH="$PWD/src"`; module paths and verbose CLI import evidence confirm the
requested checkout. Installed distribution metadata still says 0.8.6, whereas the
executed source and CLI header are 0.10.0. Scientific subprocesses use the
environment's `bin` on `PATH`.

All development minimizations use `--no-solvent`, CPU OpenMM and 200 iterations
per phase. The runner sets OMP/OpenBLAS/MKL numerical threads to one. Inputs match
the tracked fixture manifest. PLM preparation uses authoritative
`PLM=CCCCCCCCCCCCCCCC(=O)[O-]` and strict ligand parameterization. Baseline
coordinates and clash/contact counts can differ from fresh macOS preparation.
The macOS evidence archive was unavailable here, so no frozen-macOS-coordinate
cross-platform arm is claimed.

## Independently built tools

| Tool | Identity |
|---|---|
| Probe | 2.18, `091d117f5bceef87389002dccfe26df78516fdda`, system C/Make build |
| Installed xTB | 6.7.1 (`edcfbbe`), conda-forge Linux `gfortran_ha02e54d_5` |
| Fixed xTB | `a45f05d397268fd683b5f9f824ed59912c5490b5`, #1429 merge |
| Fixed-build toolchain | GNU C/Fortran 13.3.0, Meson 1.7.1, Ninja 1.13.0 |
| Fixed-build numerical libraries | Netlib BLAS/LAPACK 3.11.0, OpenMP off; optional tblite/CPCM-X disabled |
| Fixed-build upstream tests | 30/30 unit suites passed; disabled optional backends are not qualified |

The isolated build environment lives at `/tmp/opencode/probe-xtb-build`.
`ldd` resolves the fixed binary's BLAS/LAPACK to that environment's Netlib
libraries. The package arm uses its independently installed numerical libraries;
this remains a package-versus-upstream comparison, not a one-change causal proof.
Required source dependencies resolve to mctc-lib v0.5.1, jonquil v0.3.0,
toml-f v0.4.3 and test-drive v0.5.0; their exact commits are recorded.

SHA-256:

```text
Probe binary:   46317c8c6bb743b6f9ab5081ac2bfec578f423b870050e0e7e0d16d5bfa5b3c6
Installed xTB:  271dd6cdd5da127a8e39919e7e32d924790fb346f0057b4183769ba8b3411887
Fixed xTB:      60bf6c5928e570148dbc0142bfcc3b8144dd4f946fa542feb772619649757ce4
xTB .conda:     9fe81c43fbe3909e239841d8fe27a7d719f4e6441cc5a1eade784f12b7459edf
```

The downloaded `.conda` digest matches installed package metadata. Its embedded
recipe and patches are preserved and match the inspected feedstock revision
`811e76f7850e403de1b564d2dd449fe0a0d8d9cb`; patches do not backport #1429's
`frag_hess.f90` / `relaxation_engine.f90` changes.

## Execution and adapter contracts

The matrix report says `execution_complete=true`,
`no_xtb_discovery_verified=true`, `adapter_contracts_met=false`. The first flag
means the audit scripts finished, including evidence collection for native
failures/timeouts; it is not a claim that all native optimizations converged.
The discovery guard exercised fresh preparation and ordinary minimization.

| Check | Result |
|---|---|
| Focused audit/helper tests, final | 28 passed, 7 strict xfailed |
| Contracts with `--runxfail` | 12 passed, **7 failed**, as expected on unchanged adapters |
| `pytest -m 'not slow' -q` | **1488 passed, 7 skipped, 31 deselected, 7 strict xfailed**, 233.85 s |

One skip needs `biotite.structure`; six need the unavailable GROMACS `gmx`.
The full suite also reports two SWIG deprecations and four unraisable
`batch._Mirror.flush` warnings on closed streams; the log preserves them.
The final focused run additionally verifies supplemental-joint selection and
baseline-digest rejection in the summary helper; this test was added after the
full-suite run. Ruff on changed Python files, `git diff --check` and agent-doc
validation pass.
Seven contract failures reproduce the invalid Probe selection count, malformed
output treated as clean, incorrect condensed-field parsing, and absent xTB
count/element/finiteness/convergence validation. Missing executable, native
failure, timeout, empty subsystem and the 5000/5001 switch controls were exercised.
These boundary controls are not a large-system QM performance benchmark.

## xTB mechanism controls

The predeclared raw frozen-atom maximum-drift tolerance is `1e-5 Å`.

- Installed small controls converge yet drift **0.04095 Å** (ethanol) and
  **0.16482 Å** (two ethanols).
- Installed default 540-atom separated ethanol and 542-atom connected alkane
  controls activate fragmented diagonalization and terminate with **SIGSEGV
  (`returncode=-11`)** before useful optimization. A zero-drift initial frame is
  not a successful result.
- Installed `rf` returns zero but does not converge; final partial coordinate
  drift is **2.660 Å** and **3.039 Å** respectively.
- Fixed default and `rf` converge on all four controls. Largest anchor drift
  across the eight runs is **1.64e-11 Å**. Default large controls explicitly
  activate fragmented diagonalization; `rf` uses its separate optimizer path.

## Fresh matched complex results

Clashes use the same Python detector and coordinates before/after each arm.
They are not rankings of Probe against Python. All PLM arms converge.

| Arm | Observation |
|---|---|
| Installed glycan, default / `rf`, 200 or 1000 requested cycles | Converges at cycle 122, raw frozen max/RMS drift **1.623/0.971 Å**, clashes **33 → 49** after restoring anchors |
| Fixed glycan, 200 requested cycles | Returns zero without convergence; cycle banner 219, diagnostic partial geometry **33 → 60** clashes |
| Fixed glycan, 1000 requested cycles, original 90-second limit | Both native runs time out; no accepted final geometry |
| Fixed glycan, 1000 requested cycles, supplemental 300-second limit | Default / `rf` converge at cycle 582 in 105.9 / 104.3 s; raw max drift **6.23e-10 Å**, clashes **33 → 65** |
| Installed PLM, production GFN-FF | Heavy-atom RMS/max displacement **0.998/2.363 Å**, **16/28** original 4 Å contacts lost, clashes **57 → 61** |
| Installed PLM, GFN2, charge −1, UHF 0, vacuum | **0.707/1.962 Å**, **11/28** contacts lost, clashes **57 → 63** |
| Installed PLM, GFN2, charge −1, UHF 0, ALPB water | **0.640/1.446 Å**, **9/28** contacts lost, clashes **57 → 61** |
| Fixed PLM, production GFN-FF | **0.985/2.326 Å**, **17/28** contacts lost, clashes **57 → 61** |
| Fixed PLM, GFN2 vacuum / ALPB water | RMS **1.070/0.983 Å**, **17/28 / 16/28** contacts lost, clashes **57 → 63 / 61** |

The remaining GFN-FF charge/spin/water arms also converge but end at 62–64
clashes and lose 16–17/28 original contacts. Unlike total energies are not
ranked. Contact retention is geometric, not binding affinity.

Installed joint glycan linkage bonds change from **1.516/1.516 Å** to
**1.666/1.661 Å** after anchor restoration, although raw optimized bonds are
near 1.50 Å. Geometry deltas detect two changes over 0.1 Å and four angle
changes over 10°. Separate installed components still drift 0.813/0.828 Å.
Fixed separate 247-atom components converge at cycles 420/380 with raw drift
below 7e-10 Å, but each individually increases full-complex clashes to 49.
The original index correctly refuses a joint/separate comparison because the
joint calculation timed out.

The completed longer-timeout joint calculation keeps both linkage bonds near
1.494 Å and triggers no relative bond/angle or tetrahedral-sign screen, yet
nearly doubles pocket clashes. Combining the two separate fixed-component
results also gives **65 clashes**. Joint-versus-separate mobile coordinates
differ by RMS **1.579 Å**, maximum **7.639 Å**, over 476 mobile atoms.
This demonstrates grouping/basin sensitivity; neither arm contains a protein
pocket model. The initial supplemental launch was interrupted by the execution
harness and retained separately; `fixed-glycan-extended-v2` is the completed run.

Relative bond/angle and ordered tetrahedral-volume screens are not absolute
chemical or CIP validation. No candidate sign flips were found in the measured
arms. Existing hinge amide/angle/disulfide findings remain in the diagnostics.

## Additional 1VCU comparator

**The macOS strict-parameterization blocker did not reproduce on this Linux
stack.** Fresh preparation and strict whole-system minimization succeed with
one DAN and two EPE residues present. A second strict minimization uses a new
empty `DVBFIXER_LIG_CACHE` below `linux-1vcu/fresh-ligand-cache`; verbose output
confirms GAFF2 template construction for DAN and EPE. Both outputs preserve all
50 ligand heavy-atom identities, have finite coordinates and zero detected D-Cα.
The empty-cache output moves ligand heavy atoms by RMS/max 0.371/0.680 Å.
Preparation's exhausted seeded-rebuild warning is retained. This establishes
local execution, not the cause of the macOS/Linux parameterization difference.

An additional matched GFN-FF comparison uses that empty-cache OpenMM output.
Both xTB binaries and both default/`rf` optimizers converge, but full-complex
clashes increase **89 → 113**. DAN loses **42/90** original 4 Å contacts;
EPE 382 loses **10/48**, EPE 383 loses **21/34**. Heavy-atom RMS displacement
is approximately 1.067 Å for DAN and 0.448/1.60 Å for the two EPE residues.
This is joint isolated-heterogen optimization, without a protein pocket or
covalent anchors; successful native convergence does not preserve the interface.

## Probe

Both electron/nuclear radius arms use identical coordinates and explicit H,
with atom-field normalization. The production invocation still supplies
`ALL ALL` in self mode and fails; direct research calls use `-self ALL`.

| Input | Python findings | Probe overlap pairs, electron | Probe overlap pairs, nuclear |
|---|---:|---:|---:|
| ASN | 0 | 0 | 0 |
| Hinge OpenMM baseline | 33 | 55 | 28 |
| 7X35 OpenMM baseline | 57 | 88 | 45 |
| Prepared ubiquitin | 22 | 35 | 10 |
| Additional 1VCU baseline | 89 | 130 | 47 |

Raw counts reflect different radii, exclusions and contact classes, not
adjudicated precision. Controls reproduce clean separated methanes, inter-methane
clashes and exclusion of the bonded ethane contact. Probe detects octane C1–C5
at 2.62675 Å / four graph bonds, which Python misses because the atoms share a
residue. Both 180° and 90° H-bond controls retain `hb`; an angular-selectivity
advantage is not established.

## Reproduction

Run from the requested checkout in the scientific environment. Prepend its
`bin` to `PATH`; select the checkout explicitly even if an editable installation
points elsewhere:

```bash
export PYTHONPATH="$PWD/src"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python scripts/run_probe_xtb_audit.py \
  --output .artifacts/probe-xtb/linux-fresh-replay \
  --probe .artifacts/probe-xtb/linux-build/probe-source/probe \
  --xtb installed="$CONDA_PREFIX/bin/xtb" \
  --xtb fixed=.artifacts/probe-xtb/linux-build/xtb-build/xtb \
  --extended-analysis

# Review the completed supplemental run while retaining the bounded matrix.
python scripts/summarize_probe_xtb_audit.py \
  --audit .artifacts/probe-xtb/linux-fresh \
  --fixed-component-run .artifacts/probe-xtb/linux-fresh/fixed-glycan-separate-0 \
  --fixed-component-run .artifacts/probe-xtb/linux-fresh/fixed-glycan-separate-1 \
  --fixed-joint-run .artifacts/probe-xtb/linux-longtimeout/fixed-glycan-extended-v2 \
  --output .artifacts/probe-xtb/linux-summary-longtimeout-replay.json
```

The handoff documents pinned builds and test commands. Output directories must
be new. Preserve bounded failures and use separate supplemental directories for
longer-timeout replays. Build artifacts depend on the isolated Netlib environment
remaining available. Linux reproduces the adapter failures, fixation mechanism
improvement and pocket-geometry degradation; an xTB upgrade alone does not
qualify whole-complex polishing.
