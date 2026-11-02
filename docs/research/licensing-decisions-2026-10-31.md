# Licensing continuation: evidence and owner decision gate

Date: **2026-10-31**. Branch: `feature/licensing-audit`, baseline `4a328a6`.
Status: **substantial inventory/evidence work completed; not release clearance**.

[Handoff](../plans/licensing-audit-handoff.md) ·
[Evidence index](licensing-evidence-2026-10-31.json) ·
[Source register](licensing-upstream-sources.json) ·
[Original audit](third-party-license-audit.md) ·
[Model terms](../model-licenses.md)

## Evidence retained

All full captures are below `tmp/licensing-audit-2026-10-31/` **inside this
project**, not solely in system temp. The tracked index binds 12 reports and
1,244 component identities to hashes. These are package/build/metadata identities,
not 1,244 unique projects or completed legal reviews. Conda and Python records
overlap. Legal-file blobs are deduplicated by SHA-256 within each capture;
the index generator independently checks every referenced blob. Preserve the
ignored evidence directory when transferring this investigation: Git alone
does not transfer full texts or locally installed environments.

| Captured existing environment | Conda records | Python distribution records |
|---|---:|---:|
| Apple Standard (`460cce9` prefix) | 245 | 110 |
| Apple Desktop (`8f0cd33` return prefix) | 248 | 110 |
| Apple ambertools-compat Standard | 246 | 117 |
| Apple ambertools-compat Desktop | 249 | 114 |
| Apple MSA add-on | 247 | 110 |
| cc94 MPS sampling | 37 | 96 |
| cc94 CPU refinement | 157 | 51 |
| Separate Apple PSI4 | 200 | 103 |

These are read-only metadata snapshots, **not fresh solves or environment
usability checks**. Known PDBFixer asset damage in older prefixes is not repaired
or waived by capturing their metadata. Gaussian, native Probe, system frameworks,
complete dynamic linkage and Linux/CUDA/Protenix installed closures remain open.
The Linux acceptance record is historical evidence, not a substitute for its
unavailable package inventory. No SSH credentials or remote access were assumed.
Some Python 3.12 prefixes expose `lib/python3.1` as a filesystem alias. The
collector now resolves/deduplicates site-package roots, records those aliases,
and tests that they do not double-count distributions. Superseded captures are
retained; the final installed evidence is `installed-v4`, not earlier versions.

The local GUI contains 497 of 558 lockfile entries and 486 installed entries
with top-level/legal-directory legal texts. Missing entries include platform
variants; this is not automatically 61 broken dependencies. Four installed
versions **do differ** from the lock: body-parser 2.3.0 vs 2.4.0, iconv-lite
0.7.2 vs 0.7.3, qs 6.16.0 vs 6.17.0, zod 4.3.6 vs 4.6.7. The tree was preserved,
not repaired with `npm ci`. Its contents are not a lock-clean release closure.
Access to the retained installed Desktop share directory was denied; that
directory was not subsequently accessed through another tool.

Fifteen of sixteen registered upstream sources were captured in `upstream-v2`.
The CHARMM policy page was readable through the web tool, but raw archival fetch
timed out at TLS handshake; this unavailable capture remains explicit. Earlier
`upstream` captures also retain a failed guessed ParmEd root `LICENSE` path:
the pinned ParmEd tree has no such root file. Its conda package-cache license
was captured separately. Do not fabricate a missing upstream license URL.

## High-priority findings

### 1. Open Babel and ACPYPE: GPL-version compatibility is a real blocker

For Open Babel 3.1.1, tag `openbabel-3-1-1` resolves to commit
`cbd4db43f8908b874864280fdc03bf92569eebc1`. Its README and the grants in
`src/mol.cpp` / `include/openbabel/mol.h` expressly name **GPL version 2**,
without an `or later` grant. Installed `openbabel-wheel` 3.1.1.23 includes GPLv2
text (SHA-256 `2b0c8f433d672cf5c5e3e42c0fcef6920ed0388641cc0cbb443859da03edd96d`).
This is evidence of v2-only terms in essential source files, not an exhaustive
file-level bill for the wheel or proof of an unseen relicensing permission.

Installed ACPYPE 2026.10.15 declares **GPL-3.0-or-later**, independently confirmed
by that exact PyPI release. Its captured GPLv3 text hashes to
`589ed823e9a84c56feb95ac58e7cf384626b9cbf4fda2a907bc36e103de1bad2`.
DVBFixer directly imports Open Babel in chemistry paths and
`acpype.topol.MolTopol` in `acpype_export.py`. This is not a subprocess-only
architecture. The GNU compatibility FAQ explains the GPLv2-only/GPLv3 conflict
for a combined work. Whether each intended DVBFixer artifact constitutes such
a combined work requires review of the actual integration/distribution.

**Disposition:** do not approve the combined distribution by merely choosing
GPLv3 for DVBFixer. An MIT/Apache grant for original code does not erase
third-party obligations either. Do not treat optional imports, separate pip
downloads or a proposed subprocess rewrite as automatic legal cures. Any
architecture/packaging change requires a separately approved plan and, where
necessary, specialist or upstream permission review.

### 2. MODELLER: actual public EULA is now inspected

The public academic EULA grants non-transferable use. Section 5 limits use to
internal non-commercial purposes and forbids transfer/distribution of the
program or improvements without prior written permission. Section 7 further
limits copying; section 9 specifies a five-year term. The registration page
limits free availability to academic non-profit institutions; the non-academic
page directs commercial entities **and government labs** to BIOVIA.

**Disposition:** exclude MODELLER binaries/keys from any distribution decision
unless appropriate written rights are provided. User-side conda retrieval is
not redistribution by a prebuilt-image publisher, but it does not establish
eligibility of every user or authorize a commercial/hosted service. The inspected
public EULA is not evidence of the owner's individual agreement or its dates.
Do not contact licensors, register users, accept terms or request keys on the
owner's behalf.

### 3. Bundled FF rights are a present packaging issue

Both actual wheel and sdist contain the three declared FF trees. Existing
Amber19SB/GROMACS provenance is preserved. AMBER99SB-ILDN/Lipid21 names its
authors/port but does not establish the mixed tree's exact origin and grants.
The CHARMM tree names July 2022 parameters, CGenFF 4.6 and the charmm2gmx
conversion revision; the official page links the corresponding December 2022
corrected archive. An archive link or citations do **not** establish full
redistribution rights. CHARMM cap XML is adapted parameter material, not made
DVBFixer-owned by its XML wrapper.

**Disposition:** retain the current files unchanged for review. Release of
affected artifacts remains blocked pending provenance/permissions or an
owner-approved exclusion/refactoring. Do not assign GROMACS's LGPL to all
CHARMM data merely because it uses GROMACS format.

### 4. Other verified distinctions

- ParmEd 4.3.1's checked conda build declares LGPL-2.1-or-later and its cache
  includes LGPL 2.1 text. This narrows the earlier generic LGPL observation;
  the generic license's example `or later` paragraph is not itself a grant
  for every upstream file or different build.
- Installed MDAnalysis is **2.11.0**, not the older 2.10.0 preliminary example.
  Its conda declaration is LGPL-3.0-or-later; its full legal inventory includes
  separately licensed components/assets. Do not label all contents LGPL alone.
- OpenMM 8.6.1's actual `Licenses.txt` distinguishes MIT API/Reference/CPU,
  LGPL GPU platforms and additional SFMT/Rice/AsmJit notices. A conda package's
  single `LGPL-3.0-or-later` field is not a full file-level inventory.
- Protpardelle source MIT and checkpoint CC-BY-4.0 were reconfirmed separately.
  The shared adapter now populates both existing provenance fields after its
  existing supported-checkpoint digest checks. No sampling settings, request
  schema, frozen requests, expected output hashes or coordinates were changed.
- RFdiffusion's later weight clarification is now pinned to immutable commit
  `820bfdfaded8c260b962dc40a3171eae316b6ce0`, distinct from the frozen runner's
  source revision. Other model-data/HMM and fixture rights remain open.
- A narrow attribution-marker search of Python sources found no clear external
  copied-code markers; this is **not** a similarity/provenance audit and does
  not prove original ownership. Fixture derivatives still need individual review.

## Actual wheel/sdist inspection

Built a deliberate source snapshot with the existing base setuptools 75.1.0
and wheel 0.44.0, with no dependency installation and no changes to live prefixes.
No local untracked FF trees were adopted. Outputs and manifest/hash reports are
retained under `build/`; these are unpublished audit artifacts, not a release.
The snapshot predates the later two-field Protpardelle provenance edit.

- Wheel: 284 files; sdist: 462 files.
- All **120 declared packaged files** match source SHA-256 in both archives.
- Existing two notices files are included; no root project LICENSE or project
  license metadata was supplied. The old FF permission gaps are unchanged.
- Filename screening found no registered model-weight extensions or protected
  key/config filenames. This is not a comprehensive secret-content scan.
- Neither archive contains `tests/fixtures/` or deployment files. The sdist
  contains Python tests without their structural inputs; source-provenance
  documentation is also absent. Packaging remediation is a separate approved
  change, not silently performed during this audit.
- GUI browser/server bundles and runtime images are **not artifact-verified**.
  npm drift and unavailable installed GUI evidence must not be relabelled a pass.

## License-selection recommendation

**Owner clarification:** the user created DVBFixer and is asking for help choosing
its license. Asking the user to choose a licensing policy before providing an
analysis was not a valid reason to stop independent audit work. Contributor/title
verification remains a publication check, not an assumption that the user lacks
authority to discuss licensing. Git signatures alone do not establish ownership
or employer rights.

**Current owner requirement (2026-11-02):** a genuine open-source project with
public participation and commercial use permitted. The owner explicitly withdrew
the previously considered non-commercial restriction after discussing its scope.
Do not continue the NC/PolyForm route or ask the owner to define commercial use
as a licensing exclusion. No non-commercial license was adopted.

MIT and BSD-3-Clause are both suitable permissive candidates for genuinely
original code. Neither requires the author to be notified of commercial use;
both require preservation of applicable notices on redistribution.
BSD-3-Clause additionally requires prior written permission to use the
copyright-holder/contributor names to endorse or promote derived products.
That is not a general permission or notification requirement for commercial use.
The owner subsequently answered the selection questionnaire: closed commercial
derivatives are acceptable; modified online services need not disclose source;
embedding code in closed products is very important; maximum adoption takes
priority over copyleft; no additional conditions are requested.

**Owner selection (2026-11-02): BSD-3-Clause for genuinely original DVBFixer
code.** The questionnaire favored a permissive grant; after comparing MIT and
BSD-3-Clause, the owner explicitly chose BSD-3-Clause. This supersedes the
earlier MIT recommendation. Do not reopen the license-family decision or add
NC/notification/royalty conditions. The owner subsequently confirmed the full
name **Stanislav Skripkin**; the attribution is
`Copyright (c) 2026, Stanislav Skripkin`. This choice does not clear third-party
integrations or relicense their material.

Open participation also needs clear inbound contribution terms and contributor
authority; public pull requests are not an unlimited grant of relicensing rights.
The current comparison is:

| Candidate for original code | Assessment against current findings |
|---|---|
| MIT | Simple permissive grant, compatible with inclusion under either GPLv2 or GPLv3 terms where their respective conditions apply; does not require downstream changes to be open |
| BSD-3-Clause | Likewise permits commercial/closed derivatives and is GPLv2/GPLv3-compatible; adds an explicit non-endorsement condition, not author notification |
| Apache-2.0 | Explicit patent provisions, but incompatible with GPLv2-only for a combined work; therefore not the preferred default while that integration remains unresolved |
| GPL-2.0-or-later | Consider if mandatory source reciprocity is an explicit goal; licensing DVBFixer this way still does not make third-party GPLv2-only and GPLv3 components mutually compatible |
| GPL-3.0-or-later | Consider only for an intended GPLv3 distribution whose dependency boundary is cleared; not a blanket remedy for the current Open Babel issue |

A permissive grant on original code and the terms applicable to a combined
distribution are different questions. BSD-3-Clause does not relicense Open Babel,
ACPYPE, MODELLER, parameters or weights. Nor can original code that is actually
a derivative/copy of covered third-party code simply be declared BSD-3-Clause. Preserve
separate scopes and complete the provenance/combined-work review.

## Owner decisions now required

The owner has withdrawn the non-commercial requirement; preserve that decision.
BSD-3-Clause is now selected for original code; further license-family approval
is not needed. The maintainer's copyright attribution is now confirmed and the
root LICENSE has been added. Review rights to pre-existing contributions
separately from combined-distribution review before representing all material
as covered by the owner's grant.

Implementation order: add the standard, unmodified BSD-3-Clause text with
confirmed attribution; synchronize Python/GUI metadata and README; document
inbound contributions under the same license without automatic copyright
assignment; preserve third-party scopes/notices; regenerate declaration evidence
and test source/wheel inclusion. Do not condition all original-code licensing
work on finishing every transitive inventory, but do not declare combined
artifacts release-cleared while their independent gates remain open.

Distribution choices are a **separate** gate: whether to publish prebuilt
environments/containers, restricted tools or weights; how applicable third-party
terms constrain commercial/hosted deployments; and whether to authorize packaging changes or seek third-party
permissions. Original Lipid21/CHARMM archives or written grants, if available,
would help resolve their evidence gaps. No dependency has been removed or
architecture changed as a legal workaround.

Missing Linux/native/browser closures, manual evidence reconciliation and the
bibliography can continue without prior policy approval. Full WP1/WP2/WP4/WP5
acceptance remains open. The original-code BSD declaration is now implemented;
altering other distribution boundaries still requires the relevant approval.

## BSD implementation and verification (2026-11-02)

Added the standard BSD-3-Clause text with `Copyright (c) 2026, Stanislav Skripkin`.
Python metadata uses the supported setuptools>=64-compatible license table and
explicit license-file inclusion; GUI package and root lock entry declare the
same SPDX identifier. README distinguishes original code from third-party
terms. CONTRIBUTING defines same-license new submissions without assignment
and explicitly does not establish retroactive consent. CITATION.cff records the
confirmed given/family names without inventing a release date or email address.

The GUI retains matching LICENSE copies at its package root and in public
assets, so the existing installer copy and Vite client build carry the text.
Tests require byte equality with the canonical root LICENSE. These are original
project notices, not complete dependency notices.

New retained evidence:
`tmp/licensing-audit-2026-10-31/bsd-implementation-2026-11-02/`.
The wheel and sdist were built from an explicit new source snapshot, excluding
untracked local FF trees. Their LICENSE bytes match the repository; the sdist
also includes matching CONTRIBUTING.md and CITATION.cff. The wheel passed an
installed metadata/LICENSE/CLI/FF smoke in a new dependency-free venv, without
installing scientific dependencies. `wheel-smoke.json` retains hashes/results;
`build.log` retains build diagnostics. This is not third-party clearance.

Validation: **29** pure inventory/evidence/profile tests passed, GUI typecheck
passed and **407** GUI tests passed. An actual Vite client build was retained
under `gui-client/` and its LICENSE matches the canonical text. The local npm
tree's previously recorded lock drift is not repaired or waived by that build.
The earlier 2026-10-31 evidence/index remains historical and is not overwritten
with these new artifacts.

## Validation and limits

- Pure inventory/evidence/profile tests: **27 passed** with the base Python.
- Narrow provenance/cc94/research-script/inventory/evidence tests: **62 passed,
  one existing Torch-dependent test skipped** because Torch is absent from the
  selected Standard test prefix. No new skip was introduced. These tests do not
  establish usability of all resources in that older prefix.
- The intact `460cce9` Standard interpreter has no pytest installed; nothing was
  installed into it. Base Python lacks OpenMM, so the scientific checks used the
  existing `8f0cd33` return Standard interpreter with `PYTHONPATH=src` instead.
- Agent maps, user-doc checks, declaration drift and `git diff --check` passed.
- Full retained-report/blob readback passed when generating the 12-report index.
  Both archived artifacts matched all 120 source-data digests.
- Ruff is unavailable in the base interpreter and was not installed there.
  CI now includes the pure licensing tests and drift check, not a fake legal
  approval gate. Full scientific suite, fresh checkpoint sampling, lock-clean GUI
  build and final release-artifact acceptance were not performed.
