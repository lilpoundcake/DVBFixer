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

## 4. CHARMM: owner-reported GROMACS source, independent content comparison

The owner confirmed acquisition from **GROMACS**. The comparison with an official
MacKerell archive establishes a content reference, not a different acquisition
history. It identifies 59 identical files, three BUF-modified files and an
untracked/omitted reference `watermodels.dat`. The exact original GROMACS
package/release and permission for modifications remain open. The checked
GROMACS v2026.3 source lists charmm27, not this Jul-2022 tree; that single check
does not disprove acquisition through another GROMACS package or collection.
See the separate provenance note and machine-readable manifests.

## Distribution decisions remain separate

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

1. Identify the owner's exact GROMACS source/package for CHARMM and BUF changes;
   recover the corresponding Lipid21 port archive/history independently.
2. Review the binary wheel/native/data closure and actual corresponding-source
   route for Open Babel and ACPYPE's bundled Amber, not just sdists/metadata.
3. Resolve Reduce 4.15's source-specific notices; do not transfer a different
   upstream revision's license to this package.
4. Continue the broader manual dependency register, then prepare artifact-specific
   notices/source delivery. Do not publish a generic all-dependencies-BSD claim.

## Verification

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
