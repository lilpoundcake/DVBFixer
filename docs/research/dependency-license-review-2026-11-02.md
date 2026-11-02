# Dependency licensing: exact sources and distribution boundaries

Date: **2026-11-02**. Baseline: `6dae875`, branch `feature/licensing-audit`.
Status: **engineering evidence review, not a legal opinion or release clearance**.

[Evidence register](dependency-license-review-2026-11-02.json) ·
[Original audit](third-party-license-audit.md) ·
[CHARMM provenance](../provenance/charmm36-jul2022.md) ·
[Handoff](../plans/licensing-audit-handoff.md)

## What this slice adds

Seven small exact/pinned sources were captured with full-byte hashes, including
two source archives checked against PyPI's published SHA-256 values. The
GROMACS reference listing and a supplemental installed-prefix capture are also
retained under
`tmp/licensing-audit-2026-10-31/dependency-review-2026-11-02/`.
The official CHARMM archive is retained at the parent evidence directory.
No source archive was installed or executed and no scientific environment or
force-field parameter was modified.

The collector now includes Python RECORD-owned legal files **inside package
directories**, not just `.dist-info`. This found additional bundled AmberTools
notices inside ACPYPE. It explicitly excludes license-key filenames and keeps
native filenames as candidates, not proof of full dynamic-link/source closure.
The earlier 2026-10-31 inventory remains an immutable historical capture; it is
not silently rewritten to pretend it already covered these notices.

## 1. Open Babel: inspect the fork actually named by the wheel recipe

The exact `openbabel-wheel` **3.1.1.23** source archive is only 13,256 bytes and
is a build wrapper, **not the complete Open Babel corresponding source**.
Its pyproject pins engine commit
`77993b9a3b96fb9bd86249098beb97ab0fcbafc6`; its CMake recipe fetches
`https://github.com/njzjz/openbabel`, not the official upstream repository.
The wrapper tag resolves to `c6b2731dbd0a559ee56b8084b6d9997df1beb16f`.

The fork's pinned README, `include/openbabel/mol.h` and `src/mol.cpp` were
inspected and captured. They retain grants naming **GPL version 2**, without
`or later` wording. This confirms the earlier observation for these essential
files in the actual recipe source, rather than merely assuming the official
3.1.1 tag is identical to the wheel's fork. It is not an exhaustive file-level
license inventory or proof that the installed binary was built from that recipe.

Before redistributing the binary wheel/environment, audit the actual wheel's
vendored libraries, build modifications and corresponding-source route. The
recipe enables InChI and its platform build settings reference other native
libraries; the wrapper license alone cannot clear all of them.

## 2. ACPYPE: direct API route, optional Open Babel calls, additional binaries

Exact release **2026.10.15** declares `GPL-3.0-or-later` and requires
`openbabel-wheel>=3.1.1.23`. Its source archive hash is
`58fab1ab715cfd81e942c1c61fd386526af9d179735ec9e17101de33dc1a7e6b`.
The archived and selected installed `src/acpype/topol.py` match byte-for-byte
(SHA-256 `81470baacde3e6c8868f381c40beb8f2d9c5a9d268db336c28d8475537579a10`).
The release requires Python >=3.12; do not transfer this exact-release finding
to older ACPYPE releases selected by a Python 3.11 resolver.

DVBFixer directly constructs `acpype.topol.MolTopol(amb2gmx=True)` and calls
`writeGromacsTopolFiles` in `acpype_export.export_gromacs`. Public entry points
include `top --acpype` and `transplant --gromacs`. ACPYPE is imported **inside that
function**. Other pipelines importing SG-bond/glycan helpers from the same module
do not thereby execute that ACPYPE import.

ACPYPE itself has direct Open Babel/pybel calls for SMILES and chirality.
However, the inspected MolTopol constructor's `amb2gmx=True` branch bypasses
`getChirals`; the default PDB/prmtop export route is not proven to execute an
Open Babel call just because both packages are installed. This is static source
inspection, not a runtime trace or a legal finding of separation. DVBFixer has
other Open Babel API paths; invocation/CLI input normalization can matter.

The installed ACPYPE distribution also includes `acpype/amber_macos/`, with
**69 native-file candidates** plus executables and parameter data. Its bundled
LICENSE describes a mixed AmberTools distribution (GPLv3 default, LGPL and
other exceptions) and a public-domain statement for the specified `dat/leap`
force-field files. This is not a blanket public-domain grant for every local
AMBER/GROMACS port, Lipid21 modification, or all ACPYPE data. The small ACPYPE
sdist does not contain these platform Amber files. Therefore:

- Audit the wheel's bundled software/data independently from the Python wrapper.
- Do not count a separately installed conda AmberTools package as the only
  Amber distribution present in the prefix.
- Static `set_for_pip` setup prefers an existing antechamber and only wires the
  bundled tools if it is absent. That affects execution, not what is redistributed.
- Do not claim that co-installation alone proves a GPLv2/GPLv3 conflict, or that
  skipping one call automatically resolves combined-work compatibility.

## 3. Reduce: selected 4.15 terms differ from the previous generic observation

Selected build: bioconda `reduce-4.15-hbdacb55_1` on osx-arm64; package digest and
the captured license digest are recorded in the evidence register. Its metadata
says `BSD-4-Clause-UC`, but the actual license file contains **both**:

- UCSF copyright/acknowledgement and non-endorsement language, including an
  acknowledgement requirement for documentation/advertising materials.
- The J. Michael Word / ReliaSolve notice requiring that its notice remain
  unchanged and that modified program versions also be "available for free".

It is not sufficient to label this selected distribution Apache-2.0 or ordinary
BSD-3-Clause. Preserve the full text, establish the affected source portions and
resolve the meaning/requirements for any intended modified or binary export.
The inspected text is not, by itself, a blanket prohibition of commercial use
of unmodified Reduce; do not impose a new non-commercial restriction on dvbfixer.
Reduce is invoked via subprocess, which is a different boundary from the
Open Babel/ACPYPE Python APIs. Aggregation and each component's own notices still
need review for prebuilt environments.

## 4. CHARMM: acquisition unknown, content and edit history identified

The owner withdrew the GROMACS hypothesis: acquisition is **unknown**. Independent
comparison identifies 59 current matches, three BUF-modified files and omitted
`watermodels.dat`. All 62 initial-import files at `c7714fc` match the official
reference; BUF changes are in `8bb94b4` and `1bbd3bf`. This identifies an edit
timeline, not original derivation/permission for the values.

The quoted 2021/`077b6a1` converter header differs from the tracked nonempty
2022/`7b82040` file. Both resolved converter revisions have BSD-3-Clause LICENSE
texts, copyright 2021 András Wacha. That license does not automatically cover
input parameters or every output. See the provenance note and evidence register.

## Distribution decisions remain separate

### Exact macOS wheel/native follow-up (baseline `2aa74fa`)

Two actual PyPI wheels were downloaded without installation, checked against
PyPI SHA-256 identities and inspected without extraction or package execution:

| Wheel | Archive files | SHA-256 RECORD checks | Native magic candidates | Selected installed native matches |
|---|---:|---:|---:|---:|
| ACPYPE 2026.10.15, macOS 11 arm64 | 625 | 624 | 77 | 77 |
| Open Babel wheel 3.1.1.23, CPython 3.13/macOS 14 arm64 | 348 | 347 | 168 | 168 |

Only each RECORD itself is unhashed; neither wheel has unlisted files. Native
counts include versioned/duplicate copies, not distinct projects. Previous
suffix-only candidates missed extensionless executables: the collector now
recognizes binary magic as well and retains full named legal texts from archives.
It rejects missing/duplicate RECORDs, hash/size mismatches and malformed rows.
The generic zip-inspection fixture is now explicitly `.zip`, with separate real
wheel-contract tests; missing RECORDs are not silently tolerated.

The selected installed prefix is `dvbfixer-apple-standard-460cce9` from the
earlier Apple evidence. All 245 native candidates match their wheel bytes.
`otool -L` samples confirm `_openbabel.so` directly links bundled
`libopenbabel.7.dylib`, and ACPYPE's teLeap links bundled netCDF. This reads load
commands, not application execution; it is not a full dependency-resolution or
source-correspondence proof.

ACPYPE tag `2026.10.15` resolves to
`8c6626a0b1a1c71205f0784f800d4611e32b92bc`. Four pinned recipes were captured:
vendor script, macOS update script, distribution builder and release workflow.
They explain concrete remaining source/notice gaps:

- `vendor_amber.py` copies selected programs, transitive libraries and declared
  data from a conda prefix. Its legal-file list includes root `LICENSE` and
  `GNU_LGPL_v3`, but not the per-package conda license/build/source metadata.
- The macOS update defaults to AmberTools 26, without an exact conda package/build
  lock. It can use an existing prefix. This does not identify every wheel binary's
  exact upstream revision, build scripts or patches.
- An old `charmmgen_macos.tgz` supplies legacy charmmgen outside the normal
  vendoring closure. `otool` shows both x86_64 and arm64 slices in this binary.
  Source/build provenance needs its own record, not attribution to modern AT26.
- The distribution builder deliberately strips both Amber trees from the sdist
  and keeps one in each platform wheel. The small sdist therefore cannot by
  itself supply corresponding source for the vendored native payload.

ACPYPE has three named legal files; Open Babel has one. Open Babel also contains
InChI and a `.dylibs` tree with Cairo, font/X11 and gettext-related library names;
ACPYPE carries numerous native libraries. These names identify additional review
targets, not a verified license/version assignment for each. Neither a root GPL
text nor a generic upstream repository URL substitutes for all required notices
and a suitable source-delivery route. This is an evidence gap, not a finding that
upstream has violated a license; other source offers may exist.

Full manifests, legal-text blobs, recipe captures and the installed comparison
remain in the existing evidence root. Their hashes and artifact identities are
in the linked JSON register. The notice README now distinguishes these external
requirements from files bundled in dvbfixer's own wheel. No dependency has been
removed, relicensed or approved for a closed combined distribution.

Two Linux x86_64 PyPI wheels were subsequently inspected on this Mac without
installing or executing them. ACPYPE has **601 files / 53 native candidates**;
Open Babel CPython 3.13 has **346 / 165**. PyPI digests and all 600/345 hashed
RECORD entries pass; only RECORD itself is unhashed, and no files are unlisted.
The same three/one named legal texts are present. Open Babel's Linux `.libs`
payload differs from Apple's `.dylibs`; the individual component/source review
must remain platform-specific. These archive inspections do **not** establish
the contents of a selected installed Linux prefix or Linux runtime acceptance.

For the intended commercial/service policy, distinguish obligations from
prohibitions: GPLv2 and GPLv3 permit commercial use and are not AGPL. Ordinary
server use without distributing copies is not, by itself, an AGPL-style duty to
publish server source. Distributing a combined application or environment is a
different boundary and can require copyleft terms, notices and source delivery.
This does not clear MODELLER, parameter/data rights, or any particular combined
work. No source offer on behalf of dvbfixer has been made by this audit.

| Boundary | What is established | Still required |
|---|---|---|
| Original dvbfixer code | Owner-approved BSD-3-Clause declaration | Existing contribution/copy provenance; not a relicense of dependencies |
| Inspected dvbfixer wheel | Own namespace and declared FF data; third-party Python libraries are external requirements | Combined-work/API review and the bundled FF/data rights; absence of bundled GPL binaries is not an automatic exemption |
| User-created scientific prefix | Exact selected versions and notices can be identified | Operator eligibility for restricted tools and applicable use terms; co-installation is not proof of redistribution or incompatibility |
| Exported environment/container | Would include native libraries, ACPYPE's Amber payload and possibly MODELLER | Component notices/source delivery and explicit MODELLER redistribution permission; no image is cleared here |
| Generated structures/topologies | Program output is not automatically GPL/BSD | Any actually copied parameter/template material and independent input/data rights |

BSD-3-Clause is compatible with GPLv2 and GPLv3 individually. If a particular
combined work is legally subject to both GPLv2-only and GPLv3 obligations,
changing dvbfixer's root license does not provide a common license. Equally,
neither optional imports nor a hypothetical subprocess refactor is a legal
clearance decision. Record the concrete boundary and seek specialist/upstream
review where necessary, without modifying scientific architecture as a shortcut.

MODELLER's previously captured academic EULA remains a separate use/distribution
gate. Its restrictions do not make all original dvbfixer code non-commercial,
and BSD on dvbfixer does not grant rights to MODELLER. No commercial agreement,
redistribution consent or license key has been supplied or assumed by this review.

## Next work

1. Review CHARMM parameter/port terms and BUF-value derivation independently of
   the unknown acquisition; recover Lipid21 port archive/history as well.
2. Review the binary wheel/native/data closure and actual corresponding-source
   route for Open Babel and ACPYPE's bundled Amber, not just sdists/metadata.
3. Resolve Reduce 4.15's source-specific notices; do not transfer a different
   upstream revision's license to this package.
4. Continue the broader manual dependency register, then prepare artifact-specific
   notices/source delivery. Do not publish a generic all-dependencies-BSD claim.

## Verification

Native follow-up final validation: **49 narrow tests passed**. Actual wheel
RECORD checks passed for all four archived wheels, and all 245 selected Apple
installed native files matched. The collector now also accepts
`--archive ... --expected-sha256 ...` to reject a wrong artifact before reading
members or capturing legal texts; unlisted/unhashed RECORD entries stay explicit.
These checks do not imply complete native source or legal clearance.

Follow-up provenance correction: both historical converter LICENSE captures
were retained under `converter-licenses/` and independently hash-verified; both
contain the same 1,513-byte BSD text. The 62 initial-import files were re-read
from Git and their hashes retained in `initial-import-verification.json` below
the evidence root. The same 34-test narrow suite and user documentation checks
passed again. This correction does not resolve parameter redistribution rights
or the remaining native/source/compatibility audit.

The narrow inventory/evidence/profile suite passed **34 tests**. Declaration
drift checking, user documentation checks, agent maps and `git diff --check`
passed. All 42 local links in the five changed/new review/plan documents
resolved, and all seven supplemental captured-source blobs were independently
read back and SHA-256 checked. The tracked CHARMM comparison is checked against
the current local file bytes in the test suite; it remains explicitly
`redistribution_clearance: unknown`.

This slice changes audit tools, tests and documentation only. It does not rerun
checkpoint sampling or the full scientific suite and does not establish native
binary/source correspondence, complete third-party notices or legal clearance.
