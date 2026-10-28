# ADR 0012: Remove xTB refinement and repair optional Probe diagnostics

Date: 2026-10-28

Status: Accepted. Implementation authorized by the project owner.

## Context

The [macOS audit](../research/probe-xtb-macos-audit.md),
[Linux reproduction](../research/probe-xtb-linux-audit.md), and
[completed follow-up](../research/probe-xtb-followup.md) distinguish successful
execution from useful structure preparation. The follow-up includes 60 valid
Probe executions, eight failure controls and 31 xTB optimization jobs.

The installed xTB 6.7.1 violates tested frozen-coordinate constraints and crashes
on fragmented-Hessian controls. Upstream #1429 repairs that mechanism, but does
not make isolated-subsystem polishing preserve the surrounding complex:
glycan clashes increase 33 → 65 in the matched fixed-build native run, and PLM
and DAN/EPE lose pocket contacts and gain clashes. Explicit GFN2 charge/spin/
solvent arms do not establish a passing replacement.

The full-context, restrained hinge experiment converges with zero frozen drift
and all original pocket contacts retained, but increases clashes 33 → 35.
The same positional restraint without the protein context gives identical
clash pairs and contact retention in about 14 seconds rather than 28.5 minutes.
Restraint-only PLM and DAN/EPE also increase clashes; DAN/EPE additionally fails
the relative angle-change screen. No tested optimization arm passes all declared
comparative screens. These findings do not prove xTB unsuitable for unrelated
isolated-molecule applications; they do not justify maintaining this integration
as a supported complex-preparation feature.

The production adapter additionally lacks sufficient count, element, finite-
coordinate and convergence validation. Its advertised higher-quality polishing
claim is unsupported. Keeping a dependency and a user-facing route without a
qualified use case adds installation and maintenance cost and can degrade
already prepared complexes.

Probe has a different disposition: its additional intramolecular clash coverage
is demonstrated, and the research adapter validates a bounded integration repair.
Its production command/parser are broken. Atom-field formatting can change
thousands of reported contacts without moving atoms; Probe uppercases residue/
atom names but preserves case-sensitive chain IDs. These are integration
requirements, not evidence for removing this diagnostic backend.

## Decision

1. Remove xTB from supported production refinement, CLI/ZBS/GUI options,
   dependency/installation profiles and environment diagnostics. Do not silently
   substitute another engine when an obsolete xTB option is requested.
2. Retain OpenMM minimization and the separately explicit OpenBabel refinement
   route. Removing xTB does not qualify OpenBabel or justify comparative
   higher/lower-quality claims.
3. Preserve historical research, checksummed evidence and reproducible research
   tools. Decouple research helpers from deleted production xTB functions as
   necessary; the evidence must remain inspectable and must not be rewritten as
   successful scientific acceptance.
4. Keep Probe optional and repair its invocation, formatting, strict output
   validation, atom identity mapping, pair aggregation, water policy and explicit
   engine/fallback reporting. Preserve case-sensitive chains, insertion codes and
   original residue/atom spellings. Reject ambiguous lossy name mappings.
5. Treat clean native output as clean only after successful supported invocation
   and validation. Native failures and incompatible output fall back to Python
   with visible provenance. Finding-count differences and `hb` labels alone do
   not establish comparative precision or angular H-bond validation.

## Consequences

Existing `--xtb-refine`, `--xtb-cycles` and `zbs --refine xtb` invocations cease
to be accepted. Users should rerun ordinary OpenMM minimization and inspect its
diagnostics; no replacement polishing engine is selected automatically.
Managed environments no longer request xTB. Existing operator-installed binaries
are not uninstalled from a user's machine.

Historical audit failures are replaced by relevant production regression tests
for removed-option rejection and repaired Probe behavior. Probe's supported
version/format, radii convention, requested/actual engine and fallback reason
must remain explicit. A future xTB or other quantum-assisted backend would need
a new design and independent scientific acceptance, not restoration of the
removed isolated-subsystem adapter.

## Implementation verification

- Production `src/`, `environment.yml` and the generated GUI command schema
  contain no xTB route. Regression tests reject the retired CLI options and
  verify absence of the former Python API and doctor capability.
- Production Probe matches all **30** normalized native reference arms from
  the frozen research corpus, including exact pairs, overlaps and contact
  classes; every original input remains byte-identical.
- **Eight** real-subprocess protocol cases pass: six failure/fallback cases,
  a successful clean result and a valid overlap. xTB is explicitly recorded
  as removed in the current contract tool; historical failures remain archived.
- Native CLI checks with Probe's nuclear radii produce a clean ASN report,
  nine clashing-methane findings and a 67-finding 1VCU report with waters enabled.
  All three report the actual Probe engine and preserve the input. The clean
  report has no diagnostic warning banner.
- Fast Python suite: **1513 passed, 7 skipped, 31 deselected**. GUI typecheck
  passes; **400 GUI tests** pass. Ruff, strict-target mypy, agent/user document
  checks and release-metadata checks pass. The seven old strict xfails are
  replaced by passing Probe/removal regressions rather than reclassified as
  successful historical xTB runs.
- CLI reference pages are regenerated and checked with **Python 3.11**, the
  canonical CI argparse baseline; GUI schema synchronization also passes.
