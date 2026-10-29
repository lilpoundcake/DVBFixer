# Installation profiles: native Apple acceptance handoff

Date: 2026-10-29. Branch: `feature/installation-profiles`.

## Latest native checkpoint and Linux continuation

The installed-package finalization passed Standard/Desktop licensed science
(28/14 checks), 407 installed GUI tests and real Chrome acceptance. The expanded
cc94 MPS cohort published **25/25** candidates across repaired 1TIM model/Homology,
two internal 1UBQ gaps and both terminal cases, with independently verified
publications and retained private per-gap traces. Installed Protenix input
construction passed full/motif controls at 152/512/513 target residues; this
does not imply native Protenix inference. See
[the finalization record](../installation-profile-acceptance.md#native-installed-package-finalization-2026-10-29)
for exact outcomes, evidence and limitations. The instructions below retain the
original bootstrap protocol; statements of what was initially untested are
historical and superseded by the recorded results.

On Linux:

1. Fetch `feature/installation-profiles` into a separate clean checkout; record
   the actual commit and platform/driver/device versions. Preserve this Mac's
   environments, input originals, failed attempts and local evidence. Runtime
   evidence under `.artifacts/` is not transferred by Git; do not assume it is
   available on the new host or overwrite an existing Linux prefix/checkpoint.
2. Use the approved [installation-profile handoff](installation-profiles-handoff.md)
   and catalog-managed new prefixes/runners. Verify installed imports with
   `PYTHONPATH`/`PYTHONHOME` unset, preflight real CUDA tensors and the prescribed
   refinement process layout. Keep protected MODELLER license input out of argv
   and logs. Never use clobber or disabled conda safety settings.
3. Exercise **actual Protenix checkpoint inference** through public model and
   Homology paths, including `auto` full input at/below the 512-residue adapter
   budget, motif selection above it, explicit `full/local`, incomplete-input
   refusal, and source repair before context selection. The budget is not a
   free-VRAM estimate; invalid chemistry must never trigger silent cropping.
4. Re-run the current cc94 CUDA model/Homology, multigap, terminal, batch/resume
   and repair/publication paths after the latest changes. Use
   `scripts/smoke_diffusion_workflows.py` with a suitable frozen cohort and
   `scripts/check_protenix_hardware_acceptance.py` for its supported hardware
   checks. Declare seeds before sampling and preserve every outcome, original,
   prepared source, request, trace and private workspace. Independently read
   published bundles; CLI exit 0 alone does not mean every seed passed.
5. For the user-supplied 9CQD/5J13 cases use their tracked PDB **and FASTA** under
   `tests/fixtures/protein_context/`. Preserve all admission/geometry refusals;
   do not bypass the archived partial 9CQD coverage gate or call an unexecuted
   seed a scientific rejection. New production Protpardelle requests use cc94,
   never silently cc89. Keep `--no-solvent` on development minimize/ZBS runs.
   Fresh native 9CQD with all FASTA tails passed 4/5 seeds (101 rejected by
   junction/anchor geometry); accepted variants used bounded observed-anchor
   refinement, not a zero-movement baseline. Original 5J13 was refused before
   sampling because of SO4 B301–B303, and ordinary diagnose also found incomplete
   heavy atoms/hydrogen defects. Preserve that original refusal. Any future
   prepared protein-only derivative needs a separate, explicitly approved scope
   and provenance; simply removing sulfate is not established as sufficient.
6. Run focused and fast tests, GUI checks if affected, generation/version/docs
   checks and installed science smokes; record actual outcomes before claiming
   the cross-platform plan complete. Native conda package collisions and managed
   osx-arm64 Probe packaging are upstream Mac gates, not solved by Linux success.

## Executed native checkpoint

Micromamba 2.9.0 and mamba 2.10.0 Standard/Desktop licensed smokes passed;
installed Desktop passed 405 tests and a real Chrome workspace/import/diagnose
job. Standard extended checks, native GROMACS, MSA alternatives and separate
PSI4/RESP passed. Corrected cc94 MPS provisioning/preflight and one public 1UBQ
bundle also passed, with 500 MPS updates and fallback disabled. This is an
installation case, not general cc94 MPS scientific qualification.
After checkpoint `71f99e9`, the pinned native Probe source build also passed
exact-version doctor and actual electron/nuclear diagnostics with a fallback
control; see the manual alternative in the installation guide. The managed
conda-based Probe add-on remains unavailable on arm64.
Public native MPS follow-up also passed all five predeclared 1UBQ four-residue
model seeds and a Homology template-plan case, including real CPU refinement
and verified publication. These controlled installation cases do not establish
the wider multi-chain/context research gate.
The controlled 1TIM interface deletion also passed all five declared MPS seeds
after automatic bounded source repair (29 severe source findings removed; 109
warnings retained). Original/repaired input and independent movement/publication
audits are retained. The earlier unprepared cohort lost seed-101 private evidence
through partial-success cleanup; this is explicitly recorded, not replaced by
the new cohort. The unsupported `none` attempt never sampled. Wider qualification
remains open; see the detailed acceptance record.

Open gates: conda 26.9.4 macOS upstream package-file collisions; no declared
bioconda Probe 2.18 osx-arm64 build; broader cc94 MPS model/Homology/context
qualification. The archived selected-gap 9CQD request fails current public
coverage admission before sampling; preserve that result and do not bypass the
gate. Exact versions, warnings, failures and local evidence locations are in
[the acceptance record](../installation-profile-acceptance.md#native-apple-acceptance-bootstrap-2026-10-29).
The whole plan remains incomplete. No commit/push is implied by this checkpoint.

## Why the next host must be Apple Silicon

WP1–WP4 implementation and the feasible Linux WP5 acceptance are complete.
The Standard/Desktop, diffusion and optional PSI4 inputs and Linux CI/manual-path
migration are ready for native platform acceptance. Linux cannot establish
macOS binary compatibility, Desktop behavior or cc94 MPS execution. A conda
cross-solve, Docker Linux container or mocked MPS availability is not that proof.

This is a native **macOS arm64** task. Keep the Linux checkout/evidence and the
existing Mac environments. Create new native prefixes; do not copy Linux envs
or overwrite an existing engine checkout. Use ordinary prefix paths without
whitespace. Keep weights and MODELLER credentials outside Git.

## 1. Checkout and bootstrap

Fetch the pushed branch into a separate checkout and record `git rev-parse HEAD`
and `git status --short`. A Python 3.11+ bootstrap and one supported manager must
be available. The scripts prepend the current checkout's `src`; independently
check that imports do not resolve to an old editable worktree:

```bash
PYTHONPATH="$PWD/src" python3 -c 'import platform,dvbfixer; print(platform.machine(),dvbfixer.__version__,dvbfixer.__file__)'
```

Start with micromamba, then repeat the declared manager routes with mamba and
conda. Linux tested versions were 2.9.0, 2.10.0 and 26.9.4 respectively. Record
actual Mac versions. Run manager acceptance sequentially or use distinct package
caches: a mixed concurrent Linux matrix exposed a shared cache `.partial` race;
owned-prefix `--resume` completed those attempts without removing environments.

## 2. Standard and Desktop

Example new prefixes (their parent must already exist):

```bash
python3 scripts/install_dvbfixer.py --profile standard --manager micromamba \
  --prefix /path/to/new-standard --non-interactive --skip-modeller
python3 scripts/install_dvbfixer.py --profile desktop --manager micromamba \
  --prefix /path/to/new-desktop --non-interactive --skip-modeller
```

Skipping the license returns `action-required` (exit 2), not full product
readiness. Finish from a user-owned 0600 credential file, without placing its
value in argv or a report:

```bash
python3 scripts/install_dvbfixer.py --profile standard --manager micromamba \
  --prefix /path/to/new-standard --non-interactive --resume --configure-only \
  --modeller-key-fd 3 3< "$KEY_FILE"
```

Repeat for Desktop. The helper checks the real MODELLER import and real
model/Homology execution. It runs scientific smokes and records versions,
native conda builds, commands and logs below each prefix. Keep the protected
source file local; never retrieve the key from another installed environment.

Run extended acceptance in the new Standard environment:

```bash
micromamba run -p /path/to/new-standard python scripts/smoke_installation_profile.py \
  --with-modeller --extended --output /path/to/new-evidence-directory
micromamba run -p /path/to/new-standard python -m pytest -m 'not slow' -q
```

Install dev extras only in this new test environment if pytest is needed. Run
the slow suite as resources permit; every dev minimize/ZBS invocation must use
`--no-solvent` (the tracked slow tests and smoke script already do).

Desktop must pass fresh lockfile install, typecheck, tests, client/server build
and standalone HTTP smoke. Check both the installed host and a real browser
workspace/scientific job. Source GUI tests alone do not establish the installed
Desktop route. Record versions, warnings and actual result states.

## 3. Add-ons and formats

- Use a new Standard test prefix with `--addon msa-alternatives --addon probe`.
  Verify MAFFT auto plus explicit MUSCLE 5 and Clustal Omega, not just PATH.
  Check the exact supported Probe version and actual diagnostic engine/fallback.
- Run the extended script with `--msa-alternatives` in that prefix. Exercise
  PDBx/mmCIF, small-molecule CIF, SDF, native Open Babel tools, GAFF/OpenFF,
  authoritative SMILES and automatic chemistry, GLYCAM/CHARMM and topology paths.
- Verify DCD/GRO/NetCDF reading/writing, the application structural fit and
  actual cluster/Plotly output. For a supported GROMACS installation, supply
  `--gmx /path/to/gmx` and check generated TPR plus a matching coordinate input.
  Unsupported TPR writer versions must produce the explicit reader/alternative
  topology message, never an implicit topology substitution.
- PySCF/RESP is **part of Standard/Desktop** and must compute real charges.
- Optional PSI4 stays separate. `--runtime psi4` creates `<root>/environment`;
  select it using `--psi4-env /path/to/root/environment`. The validated Linux
  recipe pins PSI4 1.9.1, psiresp 0.4.2, Pydantic <2 and Libint 2.9.0. Native
  Apple availability and the actual kernel/RESP path are still **not-tested**.
  If available, run `tests/test_optional_psi4.py` with
  `DVBFIXER_TEST_PSI4_PREFIX` set. Do not install DVBFixer in that Python 3.9 env
  or change the main Python constraint to accommodate it.

Gaussian input generation was verified on Linux; running licensed Gaussian is
an external operator step, not an Apple or Standard dependency gate.

## 4. cc94 MPS provisioning and execution

Provide the official cc94 weights/config under their terms, retaining the
`model_params/{weights,configs}` layout and existing SHA256 identities. Do not
substitute cc89 or download weights implicitly. Choose a new root and contained
runner directory:

```bash
python3 scripts/install_dvbfixer.py --runtime protpardelle-mps --manager micromamba \
  --prefix /path/to/new-cc94-mps --non-interactive \
  --checkpoint /path/to/model_params/weights/cc94_epoch3100.pth \
  --runner-directory /path/to/test-runners
```

Check native Python 3.12/PyTorch 2.6.0 sampling, only the maintained Apple patch,
MPS fallback **disabled**, and the separate Python 3.12/OpenMM 8.5.2 CPU refiner.
The existing runner must pass its native tensor and refinement kernel probes.
No CPU sampling, cc89 intervention or backend substitution is an accepted fallback.

Put the contained runner directory on PATH and run a real public model request
from Standard. The helper's licensed smoke outputs include a one-residue 1UBQ
deletion and target FASTA for an installation/execution check:

```bash
PATH="/path/to/test-runners:$PATH" micromamba run -p /path/to/new-standard \
  dvbfixer model /path/to/gap-input.pdb --fasta /path/to/target.fasta \
  --backend diffusion --diffusion-model protpardelle \
  --diffusion-work-parent /path/to/native-work -o /path/to/new-native-bundle
```

A finished CLI or preflight is not candidate acceptance. Verify the bundle with
`dvbfixer.model.diffusion.publication.read_published_bundle`; `.selected()` must
only succeed for `status == 'success'`. Preserve rejected candidates, raw arrays,
traces and failure diagnostics. Record native chain/position features, 500 ODE
updates, checkpoint/source identity, actual MPS placement, CPU-refinement identity,
fixed-coordinate preservation and all independent geometry/chirality gates.

Repeat representative model/Homology and multi-seed/context cases from the
[cc94 research record](../research/9cqd-cc94-multichain-comparison.md) before making
MPS scientific claims. Linux's cc94 pilot, the 1UBQ installation case and historical
cc89 Apple cohorts do not qualify cc94 MPS. Record accepted/rejected seeds without
selecting a favorable output or relaxing admission/geometry gates.

## 5. Completion and publication

Retain a matrix for each manager/profile/add-on with **passed**, **failed**,
**action-required** or **not-tested**, plus actual commands, versions/builds and
digests. Successful solve/import, empty skips and Linux proxies are not Mac passes.
Never put credentials, checkpoints, raw private structures or environment trees
in the commit. Update the acceptance record and canonical checklist only from
executed evidence. Finish remaining WP5/WP6 platform claims, rerun required
checks, and commit/push the Mac results as a separate increment.
