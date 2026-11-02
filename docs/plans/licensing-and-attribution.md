# Licensing and third-party attribution plan

Updated: **2026-11-02**, repository baseline **`4a328a6` / 0.10.1**.
Status: **preliminary audit refreshed for implemented installation profiles;
declaration inventory, native installed evidence and preliminary artifact checks
implemented; BSD-3-Clause and Stanislav Skripkin attribution implemented on
2026-11-02; full third-party clearance pending**.

[Audit findings and sources](../research/third-party-license-audit.md) ·
[Implementation handoff](licensing-audit-handoff.md) ·
[Owner decision gate](../research/licensing-decisions-2026-10-31.md) ·
[Plan index](README.md) · [Dependency audit](environment-dependency-audit.md) ·
[Existing notices](../../THIRD_PARTY_NOTICES/README.md)

## Goal and boundaries

Give DVBFixer's own code an explicit, owner-approved license and prepare accurate
attribution for dependencies, copied/modified material, force fields, model
weights and datasets. Verify what is actually distributed in each artifact.

This plan does not select MIT, Apache or GPL by default, approve commercial use,
authorize redistribution of MODELLER/weights, or change scientific backends.
Do not download large artifacts merely to write notices. Do not remove existing
notices or relabel upstream material as DVBFixer-owned.

## Completed discovery

- [x] Inspect Python declarations/imports, subprocess tools, GUI lockfile and
  deployment/model identities.
- [x] Document code-versus-weight licensing, including MIT / CC-BY-4.0 for
  Protpardelle-1c and the Protenix v1/v2 distinction.
- [x] Identify direct GPL imports, missing project-license declaration and
  incomplete bundled-force-field notice coverage.
- [x] Record source links, evidence limits and unresolved release gates in the
  audit report. This is preliminary discovery, not a complete SBOM or legal review.
- [x] Reconcile the plan with the implemented Standard/Desktop catalog, add-ons,
  separate PSI4/ML runtimes and current packaging declarations (2026-10-31).
- [x] Recount GUI lockfile metadata and add previously omitted direct Python
  dependencies to the preliminary audit. This does not close WP1 or WP2.

## Current release decision and audit coverage

At baseline `4a328a6` there was no root project license or Python/GUI license
declaration. On 2026-11-02 the owner selected BSD-3-Clause and confirmed the name
Stanislav Skripkin; LICENSE, synchronized metadata, contribution terms and project
citation were added. Third-party notices still cover only Amber19SB/GROMACS.
Do not describe the complete 0.10.1 distribution as license-cleared or infer
third-party rights from passing installation/scientific tests.

The audit must cover **all used dependencies**, not just direct Python imports:

| Boundary | Current implementation | Required audit unit |
|---|---|---|
| Source/sdist and wheel | 20 direct Python requirements; three packaged FF trees and notices | Original/copied code, package data, exact built archives and build tooling |
| Standard | Scientific stack includes PySCF/RESP, MODELLER, AmberTools, ANARCI/HMMER and MAFFT | Exact conda builds, pip distributions, native libraries and installed data on Linux and Apple |
| Desktop | Inherits Standard; adds Node/npm and GUI | Standard plus browser assets, standalone server, external npm runtime packages and Node notices |
| Add-ons | `msa-alternatives`, `probe`, `ambertools-compat` | Actual selected solve; include `pdb2pqr` and version changes, not merely added package names |
| Separate integrations | PSI4/psiresp; operator-owned Gaussian | Separate environment or operator-supplied software and its intended use |
| Managed ML runtimes | Helper clones/pins/patches Protenix v1 and cc94 CUDA/MPS sources; separate refinement runtimes | Original and patched sources, checkpoint/auxiliary assets, transitive native/GPU stack and any distributed runtime archive/image |
| Historical/research paths | cc89 replay, Boltz-2, RFdiffusion and legacy environment specs | Distinct scope; do not silently treat these as Standard or as cleared for redistribution |

Use `src/dvbfixer/data/installation-profiles.toml` and the implemented
[installation handoff](installation-profiles-handoff.md), not the superseded
core/full/extras proposal. PySCF stays in Standard; xTB is removed from current
profiles. An installation helper downloading upstream software onto a user's
machine is not the same artifact as a redistributed prebuilt environment; record
both modes without assuming either grants blanket permission.

**At the initial refresh:** preliminary direct-component and model observations existed;
complete resolved transitive/native inventories, archived exact-version legal
evidence, compatibility decisions and built-artifact checks were absent. Unknown
coverage remained an open gate, not a clean result.

**Continuation:** exact installed metadata/legal evidence for eight Apple prefixes,
local npm evidence, 15 upstream source captures and actual wheel/sdist snapshot
inspections are now retained in the [evidence index](../research/licensing-evidence-2026-10-31.json).
This supersedes the initial no-artifact/no-installed-evidence state, not the open
gates: complete cross-platform/dynamic-link closure, all manual legal reviews,
lock-clean GUI/image artifacts and final release checks remain incomplete.

## Work packages

### WP1 — Define and inventory distribution boundaries

- [x] Implement a stdlib declaration inventory with drift checking and tests:
  [generated register](../research/license-inventory.json), regenerated by
  `python3 scripts/gen_license_inventory.py` (`--check` checks drift only).
  Its 804 records include repeated uses/lock entries/assets, not 804 unique
  dependencies. Runtime manifests remain pending review inputs; resolved
  installed/transitive coverage and legal clearance are still open below.
- [x] Capture eight existing Apple prefix inventories and available legal files,
  local npm installation differences, and hash-bound upstream evidence; preserve
  raw captures under project-local `tmp/licensing-audit-2026-10-31/`.
- [x] Inspect a staged source snapshot's actual wheel/sdist and compare all 120
  declared packaged files to source digests. This is a partial WP5 result, not
  approval of all intended artifacts or the later working tree.
- [ ] Inventory source checkout/sdist, Python wheel, browser bundle, standalone
  server and each separately distributed runner/container.
- [ ] Record exact package/build versions and distinguish required, optional,
  research, build/test and operator-provided components.
- [ ] Inventory copied code, local patches, embedded assets, parameter files,
  fixtures and auxiliary model data; do not rely only on imports.
- [ ] Capture actual Standard/Desktop/add-on solves on `linux-64` and
  `osx-arm64`, including pip-over-conda packages, channels/build strings and
  distribution hashes. A composition spec is not a resolved lock or SBOM.
- [ ] Compare the current supported catalog with historical/research specs;
  mark removed tools as historical rather than auditing them as current defaults.

Each component record must include name, version/commit, official source URL,
local use/files, source and data/weight licenses separately, copyright holders,
modifications, distribution boundary, required notices, citation and verification
status. Also record integration mode (import/link/subprocess/data), platform,
package/build digest, whether it ships or is fetched separately, evidence URL/hash,
review date and decision owner. Preserve license evidence with an immutable
upstream reference where possible; metadata alone is not final proof.

Use explicit states: `declared` (metadata only), `text-verified` (exact upstream
terms inspected), `artifact-verified` (actual shipped files/notices inspected),
and `blocked`/`unknown`. Track compatibility decisions separately: a verified
license text is not itself approval of the combined distribution.

**Acceptance:** every intended release artifact has an explicit inventory and
unknown entries remain visible rather than receiving guessed SPDX identifiers.

### WP2 — Resolve compatibility and redistribution blockers

- [ ] Review exact Open Babel GPL v2 and ACPYPE GPL v3 terms and their direct
  Python integration, including any GPL-version compatibility issue.
- [ ] Determine `-only` versus `-or-later` from actual Open Babel/wheel source
  and bundled-file terms; do not invent an SPDX expression from a generic GPL
  label. Include other copyleft and proprietary components in the same analysis.
- [ ] Determine supported project/distribution license combinations; obtain
  specialist review where needed. An optional import or subprocess refactor is
  not automatically a legal solution and must be separately justified.
- [ ] Obtain and inspect actual MODELLER terms for redistribution and intended
  academic/commercial/server use; do not distribute registration keys.
- [x] Inspect/capture the public MODELLER EULA and non-academic policy; record
  its internal/non-commercial and written-permission restrictions. Individual
  agreements and intended deployment rights still require owner evidence.
- [x] Inspect exact Open Babel wheel wrapper/fork sources and ACPYPE release
  source, distinguishing the optional amb2gmx route from other Open Babel calls.
  Record ACPYPE's bundled Amber payload and correct selected Reduce 4.15 terms;
  this source/terms evidence does not close the compatibility tasks above.
- [ ] Establish provenance and terms for the AMBER99SB-ILDN/Lipid21 and CHARMM
  bundled trees, CHARMM caps and relevant parameter datasets.
- [x] Compare the tracked CHARMM tree to a hash-bound official content reference:
  59 identical, three BUF-modified files, one omitted reference file. Preserve
  the owner's GROMACS acquisition report separately; exact package and rights
  remain open in [the provenance note](../provenance/charmm36-jul2022.md).
- [ ] Make an explicit ship/exclude/permission-required decision for each
  unresolved FF tree: all three trees are currently wheel data, not merely
  optional dependencies downloaded by the user. Do not remove them without a
  separately approved packaging change.
- [ ] Verify model auxiliary data, ANARCI HMM profiles and non-pristine fixtures.
- [x] Pin RFdiffusion's later weight-license clarification to immutable commit
  `820bfdfaded8c260b962dc40a3171eae316b6ce0`, separately from execution source.

**Acceptance:** each blocker has a documented permission/compatibility finding
or an explicit restriction on the affected artifact. Unresolved combinations
are not advertised as release-cleared.

### WP3 — Give DVBFixer itself a license

- [ ] Confirm copyright holder(s), contributor rights and any employer or
  third-party rights in the original code.
- [x] Obtain the owner's policy and license-family choice: BSD-3-Clause for
  genuinely original code. Combined-distribution compatibility remains in WP2.
- [x] Record the approved decision and its scope. A license for original code
  does not relicense dependencies, force fields, datasets or weights.
- [x] Record the owner's explicit 2026-11-02 choice of BSD-3-Clause for genuinely
  original code, followed by confirmation of
  `Copyright (c) 2026, Stanislav Skripkin`. Prior contribution-rights review
  remains open and is not closed by confirming the maintainer's attribution.
- [x] Add root `LICENSE` with approved text and copyright information.
- [x] Add synchronized Python and GUI license metadata as appropriate to their
  actual licensing scopes; ensure wheel/sdist license-file inclusion.
- [x] Add a concise README explanation of own-code and third-party terms.
- [x] Add CONTRIBUTING terms for new submissions under BSD-3-Clause without
  copyright assignment; do not treat them as retroactive consent.

**Acceptance:** the owner has approved the actual license; source and built
artifacts declare it consistently and explicitly preserve third-party terms.

### WP4 — Publish notices, model terms and scientific citations

- [ ] Expand `THIRD_PARTY_NOTICES/README.md` and include required full upstream
  license/copyright/NOTICE texts, not merely links.
- [ ] Preserve the existing GROMACS notice and Amber19SB provenance.
- [x] Add a dedicated [model-license document](../model-licenses.md) with linked commits, checkpoint names,
  SHA-256, source URLs, separate code/weight/data terms and redistribution status.
- [ ] Add provenance for older FF bundles and adapted CHARMM cap parameters.
- [ ] Add a validated scientific bibliography (for example,
  `docs/citations.bib`) for actual tools, models and parameter sets.
- [x] Add `CITATION.cff` for how to cite DVBFixer itself after owner metadata is
  confirmed; it must not substitute for third-party citations or legal notices.
- [ ] Link final documents from README/user documentation and deployment guides.

**Acceptance:** legal requirements and recommended scientific citations are
clearly distinguished; no claim assigns a wrapper's license to all its data.

### WP5 — Verify artifacts and runtime provenance

- [ ] Generate release-specific SBOMs for Python, GUI and distributed images;
  include system/BLAS/GPU dependencies where present.
- [ ] Inspect notices in wheel/sdist, actual browser bundle/server package and
  runner images. Do not count every lockfile platform variant as shipped code.
- [ ] Inspect the standalone server's external npm dependency closure:
  `--packages=external` means its JS bundle alone is not the full runtime.
- [ ] Check source/patch notices and data terms in helper-provisioned runtime
  trees and any exported archives, independently of checkpoint integrity checks.
- [ ] Satisfy applicable source-availability, modification and other
  GPL/LGPL/MPL obligations for each distribution boundary.
- [ ] Verify model weights and restricted software are not accidentally bundled.
- [x] Populate Protpardelle checkpoint provenance as CC-BY-4.0 separately from
  MIT source provenance; add narrow tests without altering scientific results.
- [ ] Verify archived/patched engine sources retain upstream notices.

**Acceptance:** artifact contents match the inventory and provenance identifies
the licenses of the actual model and checkpoint, not a generic model-family label.

### WP6 — Maintain the evidence

- [x] Add declaration-drift checking and pure inventory/evidence tests to the
  existing lightweight CI lane. These checks do not assert license clearance.
- [ ] Add proportionate checks for missing license metadata/notices and their
  inclusion in release artifacts; scanners do not replace manual review.
- [ ] Add release-checklist items for dependency, checkpoint, FF and asset changes.
- [ ] Revalidate changed versions and preserve unresolved/new terms explicitly.
- [ ] Validate local documentation links and agent-document maps after edits.

**Acceptance:** future updates cannot silently inherit approval for different
weights, versions or parameter archives.

## Recommended execution order and owner decisions

Execute WP1 and WP2 first; project-license approval in WP3 depends on them.
Draft WP4 from verified entries while blockers are resolved. WP5 validates real
builds, and WP6 makes that repeatable. This is not authority to perform an
unrequested release or packaging/backend migration.

Owner decisions: copyright attribution, desired licensing policy, supported
commercial/service distribution, which ready-made containers will be published,
and whether any weights or restricted tools will ever be redistributed.

### Immediate next deliverables, in priority order

1. **P0 / WP1:** a versioned component register plus per-artifact resolved
   inventories. Account for every direct requirement, catalog entry, GUI package,
   runner dependency and bundled/copied asset; enumerate transitive closure.
2. **P0 / WP2:** a written decision matrix for Open Babel/ACPYPE, MODELLER and
   the two older shipped FF trees. Each row needs exact evidence, intended use,
   reviewer and disposition; unresolved distribution modes remain blocked.
3. **Owner gate / WP3:** confirm rights and intended permissive/copyleft and
   commercial/service policy; approve actual license text only after the
   relevant compatibility findings. Do not add a guessed root LICENSE.
4. **P1 / WP4:** complete legal notices and separate model/data/parameter
   provenance; add scientific citations without treating them as legal clearance.
5. **P1 / WP5–WP6:** inspect real release artifacts, verify applicable source
   obligations, then automate evidence/notice drift checks for future releases.

The initial refresh completed documentation reconciliation only. It did not claim a
new full upstream-license review, grant redistribution rights, approve public
partial-atom diffusion, or change any scientific/install capability.

Implementation has since progressed beyond that initial refresh. See the
[decision report](../research/licensing-decisions-2026-10-31.md) before continuing:
the owner has clarified that license selection is the purpose of this work.
On 2026-11-02 the owner explicitly withdrew the briefly considered non-commercial
restriction. The current goal is genuine open source with public contributions
and commercial use permitted. MIT and BSD-3-Clause are suitable permissive
candidates; neither requires author notification of commercial use. After the
questionnaire and comparison, the owner explicitly chose **BSD-3-Clause** for
original code and confirmed Stanislav Skripkin as the copyright attribution.
LICENSE, metadata and contribution documentation are implemented. Do not pursue NC/PolyForm terms or
reopen the license-family decision. The combined stack is not release-approved.
Continue independent evidence work;
seek approval before changing distribution scope or removing dependencies.
No dependency was removed and no disputed FF material was relabeled as BSD.

## Definition of done

- [x] DVBFixer original code has an owner-approved explicit license declaration;
  existing third-party/contributor rights still require the reviews above.
- [ ] Every shipped third-party component has verified source/terms/notices.
- [ ] Source, checkpoints and auxiliary data are licensed separately where needed.
- [ ] Compatibility and redistribution blockers are resolved or enforced as
  artifact-specific exclusions/restrictions.
- [ ] Release artifacts include required legal material and source access.
- [ ] Scientific citations and project citation metadata are complete and linked.
- [ ] SBOM/evidence and maintenance checks match the actual supported releases.
