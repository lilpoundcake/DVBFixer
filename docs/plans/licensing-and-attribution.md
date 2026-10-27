# Licensing and third-party attribution plan

Date: **2026-10-27**. Status: **findings documented; implementation and project
license selection pending**.

[Audit findings and sources](../research/third-party-license-audit.md) ·
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

## Work packages

### WP1 — Define and inventory distribution boundaries

- [ ] Inventory source checkout/sdist, Python wheel, browser bundle, standalone
  server and each separately distributed runner/container.
- [ ] Record exact package/build versions and distinguish required, optional,
  research, build/test and operator-provided components.
- [ ] Inventory copied code, local patches, embedded assets, parameter files,
  fixtures and auxiliary model data; do not rely only on imports.
- [ ] Align with the existing environment/dependency plan without assuming its
  proposed packaging changes have already happened.

Each component record must include name, version/commit, official source URL,
local use/files, source and data/weight licenses separately, copyright holders,
modifications, distribution boundary, required notices, citation and verification
status. Preserve license evidence with an immutable upstream reference where
possible; metadata alone is not final proof.

**Acceptance:** every intended release artifact has an explicit inventory and
unknown entries remain visible rather than receiving guessed SPDX identifiers.

### WP2 — Resolve compatibility and redistribution blockers

- [ ] Review exact Open Babel GPL v2 and ACPYPE GPL v3 terms and their direct
  Python integration, including any GPL-version compatibility issue.
- [ ] Determine supported project/distribution license combinations; obtain
  specialist review where needed. An optional import or subprocess refactor is
  not automatically a legal solution and must be separately justified.
- [ ] Obtain and inspect actual MODELLER terms for redistribution and intended
  academic/commercial/server use; do not distribute registration keys.
- [ ] Establish provenance and terms for the AMBER99SB-ILDN/Lipid21 and CHARMM
  bundled trees, CHARMM caps and relevant parameter datasets.
- [ ] Verify model auxiliary data, ANARCI HMM profiles and non-pristine fixtures.
- [ ] Pin evidence for RFdiffusion's later weight-license clarification.

**Acceptance:** each blocker has a documented permission/compatibility finding
or an explicit restriction on the affected artifact. Unresolved combinations
are not advertised as release-cleared.

### WP3 — Give DVBFixer itself a license

- [ ] Confirm copyright holder(s), contributor rights and any employer or
  third-party rights in the original code.
- [ ] Ask the owner to choose permissive versus copyleft policy, then select an
  actual license compatible with WP2. Assess MIT, Apache-2.0 or appropriate GPL
  options rather than recommending one without the dependency review.
- [ ] Record the approved decision and its scope. A license for original code
  does not relicense dependencies, force fields, datasets or weights.
- [ ] Add root `LICENSE` with approved text and copyright information.
- [ ] Add synchronized Python and GUI license metadata as appropriate to their
  actual licensing scopes; ensure wheel/sdist license-file inclusion.
- [ ] Add a concise README explanation of own-code and third-party terms.

**Acceptance:** the owner has approved the actual license; source and built
artifacts declare it consistently and explicitly preserve third-party terms.

### WP4 — Publish notices, model terms and scientific citations

- [ ] Expand `THIRD_PARTY_NOTICES/README.md` and include required full upstream
  license/copyright/NOTICE texts, not merely links.
- [ ] Preserve the existing GROMACS notice and Amber19SB provenance.
- [ ] Add a dedicated model-license document with commits, checkpoint names,
  SHA-256, source URLs, separate code/weight/data terms and redistribution status.
- [ ] Add provenance for older FF bundles and adapted CHARMM cap parameters.
- [ ] Add a validated scientific bibliography (for example,
  `docs/citations.bib`) for actual tools, models and parameter sets.
- [ ] Add `CITATION.cff` for how to cite DVBFixer itself after owner metadata is
  confirmed; it must not substitute for third-party citations or legal notices.
- [ ] Link final documents from README/user documentation and deployment guides.

**Acceptance:** legal requirements and recommended scientific citations are
clearly distinguished; no claim assigns a wrapper's license to all its data.

### WP5 — Verify artifacts and runtime provenance

- [ ] Generate release-specific SBOMs for Python, GUI and distributed images;
  include system/BLAS/GPU dependencies where present.
- [ ] Inspect notices in wheel/sdist, actual browser bundle/server package and
  runner images. Do not count every lockfile platform variant as shipped code.
- [ ] Satisfy applicable source-availability, modification and other
  GPL/LGPL/MPL obligations for each distribution boundary.
- [ ] Verify model weights and restricted software are not accidentally bundled.
- [ ] Populate Protpardelle checkpoint provenance as CC-BY-4.0 separately from
  MIT source provenance; add narrow tests without altering scientific results.
- [ ] Verify archived/patched engine sources retain upstream notices.

**Acceptance:** artifact contents match the inventory and provenance identifies
the licenses of the actual model and checkpoint, not a generic model-family label.

### WP6 — Maintain the evidence

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

## Definition of done

- [ ] DVBFixer original code has an owner-approved explicit license.
- [ ] Every shipped third-party component has verified source/terms/notices.
- [ ] Source, checkpoints and auxiliary data are licensed separately where needed.
- [ ] Compatibility and redistribution blockers are resolved or enforced as
  artifact-specific exclusions/restrictions.
- [ ] Release artifacts include required legal material and source access.
- [ ] Scientific citations and project citation metadata are complete and linked.
- [ ] SBOM/evidence and maintenance checks match the actual supported releases.
