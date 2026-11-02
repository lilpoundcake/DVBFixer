# Licensing audit implementation handoff

**Closed for the agreed source-plus-installer scope (2026-11-02).** Start with
[the final distribution policy](../licensing-status.md) and
[completion record](../research/licensing-completion-2026-11-02.json), not the
historical full-environment audit sequence below. No further owner decisions
are required for this task. Keep force fields, dependencies and installation
behavior unchanged; no prebuilt environment/container is a publication target.
Remaining third-party legal uncertainties are retained, not marked resolved.
Final checks: 63 narrow tests passed; 42 user-documentation pages and agent maps
validated, version metadata remains 0.10.1. Eight notices match fresh bounded
snapshot sdist/wheel and isolated no-deps install; archive screens contain no
native candidates or screened sensitive names. Artifact/input hashes and logs
are retained in `notice-final-readback/` below the existing audit evidence root,
with report identity in the completion record. This does not publish a release
or assert scientific/whole-environment/legal qualification.

Date: **2026-10-31**. Starting point: **main `4a328a6`, 0.10.1**.
Working branch: **`feature/licensing-audit`**.

[Policy and work packages](licensing-and-attribution.md) ·
[Preliminary findings](../research/third-party-license-audit.md) ·
[Continuation / owner gate](../research/licensing-decisions-2026-10-31.md) ·
[Plan index](README.md)

## Mandate and non-goals

Implement the repository-local licensing plan, beginning with a reproducible
inventory. This document is the continuation contract, not a harness task list.
Do not pick a project license, infer owner/contributor rights, remove scientific
dependencies, change packaging or approve commercial redistribution implicitly.
No licenses are granted by an inventory or a successful test.

The working tree initially contained three modified licensing documents from
the preceding review and untracked `FF/amber14sb.ff/`,
`FF/charmm36-mar2019-cphmd.ff/`, `FF/old_amber19sb.ff/` trees. Preserve them;
the untracked trees are not approved release inputs. Do not add them to git or
the generated inventory merely because they exist locally.

## Deliverables and sequence

### A — Reproducible declaration inventory (WP1, first implementation slice)

- [x] Add a stdlib-only generator and checked-in JSON inventory.
- [x] Enumerate Python runtime/dev/build requirements, profile/add-on entries,
  external integrations, GUI direct declarations and every npm lock entry.
- [x] Hash the tracked packaged FF trees, notices and Python package data;
  record environment/runner manifests and patches as pending review inputs.
- [x] Keep declared license metadata separate from verified text, distribution
  inclusion and compatibility. Never convert a scanner result to legal approval.
- [x] Add deterministic/drift/failure tests without importing scientific modules.

This slice is **not** an installed SBOM. YAML runtime declarations are initially
hashed review inputs, not parsed/resolved transitive dependency inventories.
The generator must say so explicitly. Source fixture provenance, copied-code
search and runtime binary/library closure remain work, not hidden exclusions.

### B — Resolved inventories and evidence register (remaining WP1)

1. Select exact Linux and Apple Standard/Desktop prefixes and add-on solves.
   Preserve conda explicit package URLs/builds, pip distribution versions and
   file hashes; do not install/rebuild environments just for notices without
   agreeing the execution scope.
2. Inventory separate PSI4, Protenix sampling/refinement and cc94 CUDA/MPS
   runtimes. Include actual native/BLAS/compiler/GPU libraries, auxiliary data,
   npm server closure and browser bundled assets. Declare historical/research
   scope separately.
3. Create reviewed records keyed to package/build/artifact digests: upstream
   URL/commit, full license/copyright/NOTICE text hashes, local use, integration
   mode, modifications, distribution mode, reviewer/date and unresolved issues.
4. Audit copied/adapted code, fixtures, downloaded CCD/HMM/model data and force
   fields independently of wrapper-code licensing. Reconcile every inventory
   item to an evidence record or an explicit unknown.

**Gate:** every intended artifact has a resolved closure and all omissions are
visible. Declared dependency enumeration alone cannot close WP1.

### C — Compatibility decisions (WP2, release blockers)

| Priority | Issue | Required finding |
|---|---|---|
| P0 | Open Babel/wheel + ACPYPE direct imports | Exact source/build terms, GPL editions and any `or-later` grant; analysis of the actual combined distribution, not a guessed optional-import exemption |
| P0 | MODELLER | Actual terms for academic/commercial/service use and redistribution; operator-provided and bundled modes distinguished; no keys in evidence |
| P0 | Packaged AMBER99SB-ILDN/Lipid21 and CHARMM trees/caps | Original archives/port history, modifications, rights and notices for every part; explicit permission/restriction or separately approved exclusion |
| P1 | Other copyleft/native/runtime/data material | Exact LGPL/MPL and model/data obligations for each distribution boundary; source-availability and modification obligations where applicable |

Request specialist review for unresolved compatibility questions. Neither a GPL
dependency nor a permissive original-code grant alone settles the licensing of
all combined artifacts. Keep ship/exclude/permission-required decisions explicit.

### D — Owner decision and publication (WP3–WP4)

After relevant C findings: confirm copyright/contributor/employer rights and
intended permissive/copyleft, commercial/service and redistribution policy.
Obtain approval of actual license text and scope. Only then add root LICENSE,
Python/GUI metadata and user-facing explanation. Complete third-party legal
texts, source/patch notices, model/data terms and separately validated citations.
Preserve Protpardelle MIT **source** versus CC-BY-4.0 **weights** and Protenix
v1 versus v2 distinctions.

### E — Real artifact verification and maintenance (WP5–WP6)

Build and inspect wheel/sdist, browser/server packages and any intended
environment/container exports. Compare contents to inventories, check notices
and applicable source delivery, verify restricted artifacts/credentials are not
bundled, then add release checks. Do not add CI claims of legal clearance based
on metadata presence. Packaging changes require their own approval.

## First-slice execution and validation

Expected new files: this handoff, `scripts/gen_license_inventory.py`,
`docs/research/license-inventory.json`, `tests/test_license_inventory.py`.
Existing licensing plan, audit and index will link the handoff and describe the
limited implemented coverage. No future links are considered complete until
the referenced files exist.

Run:

```bash
python3 scripts/gen_license_inventory.py
python3 scripts/gen_license_inventory.py --check
python3 -m pytest tests/test_license_inventory.py -q
PYTHONPATH=src python3 -m pytest tests/test_license_inventory.py tests/test_installation_profiles.py -q
python3 scripts/check_agent_docs.py
git diff --check
```

Validate local links in changed Markdown. Record counts/test results below,
and update only completed A checkboxes. Do not mark B–E complete from A results.
No scientific outputs or environment mutations are needed for A.

## Continuation status

Slice A implemented on the named branch: 804 records, including 558 npm lock
entries, 40 catalog entries, 120 tracked packaged files and 24 pending review
inputs. Counts are observations of this baseline, not unique dependency counts.
Further work now includes eight retained Apple installed inventories, local npm
evidence, upstream legal sources, actual wheel/sdist snapshot inspections, a
portable evidence index and a concrete compatibility/owner decision report.
Tools: `collect_license_evidence.py`, `fetch_license_sources.py` and
`summarize_license_evidence.py` under `scripts/`. Full retained reports/texts are
in `tmp/licensing-audit-2026-10-31/`; final installed/npm/archive captures use
`installed-v4`, `npm-v3`, `build/inspection-v3-*`, and upstream uses `upstream-v2`.
Earlier failed/partial captures remain preserved. Run the summarizer only with
those retained inputs; it verifies blob digests and does not fetch missing data.

Full B/C/E acceptance and D's owner/legal gates remain open. The owner clarified
that they created DVBFixer and need assistance choosing its license. The linked
report initially recommended MIT conditionally. On 2026-11-02 the owner briefly
considered an NC restriction, then explicitly withdrew it: the current goal is
genuine open source, public contributions and commercial use permitted.
The owner questionnaire favored permissive terms. After comparing alternatives,
the owner explicitly selected **BSD-3-Clause** on 2026-11-02 for original code.
The owner confirmed **Stanislav Skripkin**; the standard LICENSE, synchronized
metadata/README, inbound contribution terms and CITATION.cff are now implemented.
Wheel/sdist license inclusion, clean installed wheel smoke and actual GUI client
LICENSE copying passed; evidence is in the decision report's 2026-11-02 section.
Review existing contribution rights and preserve third-party scopes. Do not
reopen the license-family choice or ask for the already-confirmed author name.
BSD does not require commercial-use notification; its additional clause concerns
endorsement, not ordinary use. Do not continue the NC/PolyForm route. GPL
combined-distribution and other third-party constraints remain separate.
Do not demand a policy choice
before doing independent audit work; approval is needed before making the grant
or changing distribution/architecture, not before explaining the options.
Missing Linux/GUI/native evidence and manual text review remain explicit
follow-up work, not completed tasks. No commit, merge or release is
implied. No environments were installed or modified. Ruff is absent from the current `python3`; report that limitation
rather than installing into the user's base interpreter.

Validation on 2026-10-31: the new inventory and installation-profile tests passed
**15/15** using `PYTHONPATH=src` (the base interpreter has no installed dvbfixer).
Declaration drift check, agent-document maps and `git diff --check` passed;
all 52 local file links in the four changed Markdown documents resolved.
Those were initial slice-A results. Continuation added pure evidence tests and
narrow scientific/provenance tests, CI drift checks, model terms and an initial
three-entry bibliography. Final test results and artifact limitations are in the
continuation report. The only runtime-code change is two provenance license
fields; sampling settings and coordinates are unchanged. Legal clearance remains
unperformed. The full scientific suite and checkpoint sampling were not rerun.

### Exact-source dependency continuation (2026-11-02)

Foundation commit **`6dae875`** was pushed to `origin/feature/licensing-audit` at
the owner's request. The next source-review slice is documented in
[dependency review](../research/dependency-license-review-2026-11-02.md) and
[CHARMM provenance](../provenance/charmm36-jul2022.md). Full captures are retained
under the original ignored evidence root, with a new
`dependency-review-2026-11-02/` subtree. Do not overwrite historical snapshots.

Key corrections: inspect the actual Open Babel wheel recipe's fork, not only
the official upstream tag; ACPYPE's platform wheel contains additional Amber
native/data material absent from its small sdist; the default `amb2gmx=True`
constructor bypasses its chirality/Open Babel call, so co-installation alone is
not a proof of combined-work incompatibility. Selected Reduce 4.15 is not simply
Apache-2.0. The collector now reads RECORD-owned package-internal legal files.

The owner withdrew the GROMACS acquisition hypothesis: original source is
unknown; do not keep requesting a remembered package/archive. Independent
comparison finds 59 current matches, three BUF changes and omitted
`watermodels.dat`. All 62 initial-import files at `c7714fc` match the official
reference; BUF edits are in `8bb94b4` and `1bbd3bf`, without established original
value derivation. The quoted 2021 converter header differs from the tracked
2022 header. Both resolved converter revisions have BSD-3-Clause LICENSE texts;
parameter/port rights remain separate. Continue authoritative upstream research.

No scientific pipeline, dependency composition or force-field bytes were changed
in this slice. Review outcomes still do not establish redistribution clearance.
Validation for this slice: 34 narrow tests passed; declaration drift, user docs,
agent maps, changed-document local links and diff checks passed. Supplemental
source blobs were independently hash-verified. No scientific suite/sampling
acceptance is inferred from these checks.

### Native wheel follow-up

Exact ACPYPE 2026.10.15/macOS arm64 and Open Babel wheel 3.1.1.23/CPython 3.13
macOS arm64 archives now have PyPI-digest and internal RECORD verification.
All 77/168 native candidates match the selected installed Apple prefix. The
collector retains named archive legal texts and detects extensionless binaries
by magic; 47 narrow tests passed. Full reports/recipe captures are retained in
the existing dependency-review evidence root and bound by the review JSON.

ACPYPE's pinned vendoring recipe omits per-package conda legal/source metadata,
uses a non-build-locked AmberTools prefix and separately adds old charmmgen.
Its slim sdist deliberately omits the native payload. This establishes concrete
follow-up targets, not upstream noncompliance or legal clearance. The notices
README now contains a scoped external-dependency index, not a complete export
notice/source package. Continue per-binary source/build/notice mapping, the
GPL combined-work question and Linux evidence; do not infer resolution from
the 245 installed-byte matches or four sampled otool load-command reports.

Linux follow-up now includes archive-only inspection of the exact x86_64 wheels:
ACPYPE 601 files/53 native candidates, Open Babel CPython 3.13 346/165; all hashed
RECORD members and PyPI identities verified. No Linux installation/execution or
installed-prefix correspondence is implied. Original CHARMM-format July-2022
archive search is also recorded separately in the provenance note, with its
bounded-search limitations and no inferred grant.
Final narrow suite: 49 passed after predeclared archive-digest and unlisted
RECORD regression controls. Declaration inventory was regenerated for the
changed notice README; no dependency or scientific behavior changed.

### Component/source continuation after `311e596`

Pinned ACPYPE charmmgen containers match both wheel binaries. A build script and
September 2026 rebuild commits supersede the suggestion of necessarily unchanged
2016-era binaries; the component is legacy, but these builds were refreshed.
The script shallow-clones unpinned AmberClassic and applies an environment-name
patch, so exact source revision/closure and grant applicability remain open.
Do not claim the historical GPLv2-or-later root-license candidate proves a
particular wheel's rights or inherit ACPYPE's wrapper grant automatically.

Both Open Babel wheel InChI headers identify 1.04/custom Licence No. 1.0 or later.
The correct source path is `include/inchi_api.h`, not the different
`include/inchi/inchi_api.h`; explicit matching path evidence is retained. The
primary five-page license PDF was read: preserve sections 1/4/6 notice/source/
linking requirements, section 3's unexercised GPL-conversion option and section
15's modified-output naming condition. No binary/source proof or relicensing is
inferred. Cairo's Linux recipe-named 1.15.14 archive supplies implementation
LGPL2.1-or-MPL1.1 terms; auxiliary GPL files are separate and Apple provenance
remains open. Ten pinned source/history captures, the InChI PDF and the Cairo
archive/legal texts are retained; five report hashes are in the existing JSON.
The first charmmgen extraction assumption and the nonmatching InChI header
comparison remain preserved as partial evidence; use the v2 comparison reports.
Validation for this continuation: 51 narrow tests and user documentation checks
passed; notice inventory regenerated. No complete corresponding-source delivery,
closed-combination permission, scientific execution or legal clearance is claimed.

### Reduce continuation after `c6cb14f`

Selected osx-arm64 build `hbdacb55_1` archive and its embedded recipe are retained.
The source archive matches the recipe digest; tag v4.15 resolves to `e7485302`.
Historical bioconda `098d2ea5` template/build script match the package bytes.
Four directory/root notices plus additional inspected source blocks now ship as
notice-only files; preserve their text, source hashes and appended-final-LF
convention. UCSF acknowledgement is included, but the modified-versions-free
condition is not resolved by a generic BSD label. Dictionary equality across
source/package/installed copies does not clear underlying data rights.

Installed Reduce differs from the package executable. Binary prefix relocation
is declared, but not proven to account for all changes; do not claim equality
or silently patch/reinstall to make it match. Eight pinned sources and three
retained report hashes are recorded. No independent rebuild, scientific execution,
Linux selected-build mapping or full distribution clearance is inferred.
Validation: 53 narrow tests, user docs, declaration drift and diff checks pass;
three retained report hashes, captured blobs and nine supplemental source notice
blocks checked. New notice files are included in the regenerated declaration
inventory. Fresh wheel/readback is pending because base Python lacks `build`;
no package was installed into the user's environment to bypass this limitation.

### Dictionary and notice packaging after `78192ad`

Reduce's pinned updater and wwPDB primary usage/CCD pages now establish the
declared legacy-HET input and explicit transformation steps. Current PDB archive
policy is CC0; this is positive upstream-data evidence, not a blanket grant for
Reduce edits/software. Observed dictionary: 39,488 residues, 9,431 xxCONECT and
2,139,516 CRLF records. Historical source snapshot and exact regeneration remain
unestablished; the updater was not executed.

Fresh notice readback is no longer blocked by the missing `build` frontend:
existing setuptools.build_meta built a bounded tracked-source snapshot without
tool installation. Seven notice files match sdist, wheel and isolated no-deps
pip target. Source snapshot/input hashes, artifacts and logs are retained under
`notice-artifact-readback/`, helper in `build-notice-readback.py`, below the audit
root. This is not a full-repository release sdist or scientific smoke. No native
payload/screened sensitive filename found. Do not infer environment export or
legal clearance from this packaging check; source-data/combined-work gates remain.

### Owner steering: practical scope, not exhaustive forensics

The owner explicitly said the audit was going too deep. The primary continuation
is now [practical publication status](../licensing-status.md): own BSD grant,
source/wheel force-field rights, one focused Open Babel/ACPYPE combined-work
review, and operator-licensed MODELLER. Do not keep downloading or mapping every
native library merely to improve a research inventory. Preserve existing work;
resume per-library export compliance only when a concrete prebuilt artifact is
actually planned. No architecture, dependency, data-exclusion or relicensing
change is authorized by this scope reduction. Upstream questions are proposed,
not sent or answered. Unknowns remain explicit rather than called violations.
Practical-scope validation: 53 narrow tests, user/agent documentation checks,
64 changed/new local links, JSON and diff checks pass. Dictionary/build evidence
already gathered is retained and indexed, not expanded into a new full-platform
qualification. Three retained reports and seven current installed notice texts
were independently read back; publication checklist is not legal approval.
