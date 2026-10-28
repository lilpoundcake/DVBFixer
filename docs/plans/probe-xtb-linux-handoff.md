# Probe / xTB Linux handoff

After this historical replay, [ADR 0012](../adr/0012-remove-xtb-refinement-and-repair-probe.md)
removed production xTB refinement and repaired Probe. The pre-removal contract
tests/commands below belong to commit `3b4ac0b` or the transferred source snapshot.
Current production tests are `tests/test_probe_integration.py` and
`tests/test_xtb_removed.py`; independent native research tools remain available.

The [macOS audit](../research/probe-xtb-macos-audit.md) completed its bounded
execution matrix and failed acceptance gates. The subsequent
[Linux continuation](../research/probe-xtb-linux-audit.md) completed fresh
preparation and its 46-task matrix; scientific acceptance also failed. This
document retains the original replay protocol.
Reproduce the failures as well as successes; do not interpret zero runner exit
as chemical acceptance. Production fixes and packaging decisions are separate.

## Transfer and verify

Transfer `.artifacts/probe-xtb/probe-xtb-linux-handoff.tar.gz` and its adjacent
`.sha256` file. The archive contains a repository snapshot plus uncommitted audit
scripts/tests/docs, selected checksummed inputs, macOS evidence, source archives
and build provenance. It excludes macOS executables/build objects and unrelated
local work. Keep the original git checkout if you want branch history: the
portable source tree is not a git repository. The archive was built before the
audit commit; the `test/probe-xtb` branch carries the audit scripts/tests/docs.
Ignored `.artifacts/` evidence is not included in a git push: transfer the archive
and its digest separately.

```bash
sha256sum -c probe-xtb-linux-handoff.tar.gz.sha256
mkdir probe-xtb-handoff
tar -xzf probe-xtb-linux-handoff.tar.gz -C probe-xtb-handoff
cd probe-xtb-handoff
sha256sum -c MANIFEST.sha256
cd repository
```

The manifest verifies archived regular files, including baseline coordinates and
raw logs. macOS paths in historical commands are provenance, not Linux commands.
Native-fixture scripts in evidence have macOS shebangs: **regenerate them with the
runner**, never execute the archived fixtures on Linux.

## Scientific environment

Use the repository's [installation instructions](../installation.md) and a new
Linux scientific environment from `environment.yml`; preserve Python
`>=3.11,<3.14`, NumPy `<2.5` and MODELLER license handling. Do not overwrite an
existing environment or transfer the macOS environment as Linux executables.
Record the Linux solver output/build strings independently. Installing build
tools below requires an isolated build environment, not edits to production deps.

On Docker/VirtioFS use a native container-overlay micromamba root, not a macOS
bind mount (see [known issues](../known-issues.md)). Keep the scientific env's
`bin` directory on PATH so child `dvbfixer`/AmberTools calls resolve correctly.

## Build the independently pinned tools

`sources/` contains source tar files with full upstream license files and git
commit provenance in `provenance.json`. Build in fresh writable directories.
No sudo/global installation is required. Probe source is commit
`091d117f5bceef87389002dccfe26df78516fdda` (2.18):

```bash
mkdir -p ../build/probe
tar -xf ../sources/probe.tar -C ../build/probe
make -C ../build/probe
```

Fixed xTB is `a45f05d397268fd683b5f9f824ed59912c5490b5`, the #1429 merge.
Source snapshots include mctc-lib v0.5.1, jonquil v0.3.0, toml-f v0.4.3 and
test-drive v0.5.0 as their checked-out exact commits; provenance records them.
With GNU C/Fortran and Netlib development BLAS/LAPACK available in the isolated
build environment, use Meson 1.7.1 / Ninja 1.13.0:

```bash
mkdir -p ../build/xtb-source
tar -xf ../sources/xtb.tar -C ../build/xtb-source
FC=gfortran meson setup ../build/xtb-build ../build/xtb-source \
  --buildtype=release --wrap-mode=nodownload \
  -Dbuild_name=a45f05d-audit \
  -Dlapack=netlib -Dopenmp=false -Dtblite=disabled -Dcpcmx=disabled
meson compile -C ../build/xtb-build -j 4
meson test -C ../build/xtb-build --no-rebuild --suite unit --print-errorlogs
```

These are Linux replay instructions, not a claim the Linux build already passed.
Record compiler, Meson, Ninja, BLAS/LAPACK implementation, `ldd` output, executable
digests and upstream test logs. Reject an unintended Accelerate/MKL/OpenBLAS
substitution when reproducing the Netlib arm; other backends require labelled
sensitivity arms. Never reuse macOS SHA-256 as an expected Linux binary digest.

The installed/package arm should use independently resolved Linux xTB 6.7.1,
with its package digest, recipe revision and patches recorded. Its build string
will differ from macOS. If #1429 has been backported, label the arm accordingly;
do not silently call it an unpatched baseline. A version string is insufficient.

## Replay

Run in the scientific environment. Set `OMP_NUM_THREADS=1`,
`OPENBLAS_NUM_THREADS=1`, `MKL_NUM_THREADS=1`; the runner also enforces these.
If an editable installation resolves another worktree, set
`PYTHONPATH="$PWD/src"` and verify `dvbfixer.__file__` before running.

```bash
python scripts/run_probe_xtb_audit.py \
  --output .artifacts/probe-xtb/linux-fresh \
  --probe ../build/probe/probe \
  --xtb installed="$CONDA_PREFIX/bin/xtb" \
  --xtb fixed=../build/xtb-build/xtb --extended-analysis

python scripts/summarize_probe_xtb_audit.py \
  --audit .artifacts/probe-xtb/linux-fresh \
  --fixed-component-run .artifacts/probe-xtb/linux-fresh/fixed-glycan-separate-0 \
  --fixed-component-run .artifacts/probe-xtb/linux-fresh/fixed-glycan-separate-1 \
  --output .artifacts/probe-xtb/linux-summary.json

pytest tests/test_probe_xtb_audit.py tests/test_probe_xtb_research_metrics.py -q
pytest tests/test_probe_xtb_audit.py --runxfail -q
pytest -m 'not slow' -q
```

The `--runxfail` command is **expected to fail** on unchanged adapters. Preserve
its XML/logs. Seven strict xfails must not be described as passing safety gates.
Extended mode adds component and geometry tasks to the base 24-task protocol.
Directories must not exist. Fresh baseline preparation uses authoritative PLM
SMILES and strict parameterization; preserve any failures, never strip ligands
to create a substitute successful run.

For a separate numerical-platform comparison, repeat in another output directory
with `--baselines ../evidence/complete-macos/baselines`. This reuses macOS frozen
coordinates; it is **not** independent preparation or a new platform qualification.
The report marks `no_xtb_discovery_verified=false` in frozen-baseline mode because
preparation/minimization were not exercised.

## Required review gates

1. Inspect nested reports, native return codes, convergence banners and final XYZ,
   not merely runner return codes. Partial trajectories are diagnostics only.
2. Validate count, element sequence, finite coordinates and stable extraction maps.
   Same-element permutations cannot be proven absent using XYZ alone; preserve
   engine order guarantees and reject unsupported identity changes.
3. Require raw frozen max drift ≤`1e-5 Å` **and** acceptable geometry after merge.
   Restoring anchors cannot turn a drifting optimization into a valid one.
4. Compare matched OpenMM-only bonds/angles, tetrahedral signs, planarity findings,
   clash pairs, ligand displacement and 4 Å pocket contacts. Relative screens
   (>0.1 Å bond change, >10° angle change) are not absolute chemistry validation.
5. Keep GFN-FF and explicit GFN2 (charge −1, UHF 0, vacuum/ALPB-water) separate.
   Do not rank unlike total energies or infer affinity from contact retention.
6. Recheck Probe command/parser/format/fallback gaps. Retain radius arms, format
   normalization, known bonded and intramolecular controls and `hb` ambiguity.
7. Record the 5000/5001 boundary and empty-heterogen behavior as control-flow tests,
   not a large-system QM performance benchmark.

Expected qualitative macOS findings: Probe has useful intramolecular coverage but
broken production invocation/parser; installed xTB drifts; #1429/Netlib fixes the
tested fixation mechanism but still worsens pocket clashes. Exact fresh-prep
coordinates/counts are not bitwise acceptance criteria across platforms.
Only after reviewing Linux evidence should dependency/profile or production
algorithm changes be proposed.
