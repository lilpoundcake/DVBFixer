# Licensing audit implementation handoff

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
