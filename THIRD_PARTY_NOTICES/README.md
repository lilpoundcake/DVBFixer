# Third-Party Notices

This notice set is not yet a complete component inventory or release license
clearance. See the [preliminary audit](../docs/research/third-party-license-audit.md)
and [licensing plan](../docs/plans/licensing-and-attribution.md) for additional
dependencies, model-weight terms and unresolved attribution work. Those documents
do not change the licenses of any bundled material.

`FF/amber19sb.ff` is an unmodified copy of the GROMACS `v2026.3`
Amber19SB force-field directory. GROMACS is distributed under the GNU Lesser
General Public License version 2.1. The full upstream license is included as
`GROMACS-LGPL-2.1.txt`; force-field attribution and citations remain in
`FF/amber19sb.ff/forcefield.doc`.

Full source, commit, checksums, comparison evidence, and compatibility audit:
[`docs/provenance/amber19sb-gromacs.md`](../docs/provenance/amber19sb-gromacs.md).

## External scientific dependencies — partial notice index

These components are external Python/environment requirements, not binaries
copied into dvbfixer's own wheel. An exported environment can include them and
their additional native libraries. This index is not the complete license-text
or corresponding-source set needed for such an export.

| Component reviewed | Identified terms and scope | Open distribution work |
|---|---|---|
| Open Babel wheel 3.1.1.23 | Inspected core source grants name GPLv2 without an or-later grant; the selected wheel includes a GPLv2 text | Combined-work review, vendored-component notices and exact corresponding-source delivery |
| ACPYPE 2026.10.15 | Python wrapper declares GPL-3.0-or-later; platform wheel includes a separately licensed mixed AmberTools/native/data payload | Individual bundled-component terms/source, legacy charmmgen provenance and API compatibility |
| Reduce 4.15, selected bioconda build | UCSF acknowledgement/non-endorsement and Word/ReliaSolve notices; not generic Apache-2.0 | Preserve full selected notices and assess modified/binary distribution requirements |
| MODELLER | Separately licensed; project BSD does not grant MODELLER use or redistribution rights | Operator eligibility and any necessary commercial/service/redistribution permission |

Exact artifact identities, retained full-text hashes and unresolved questions:
[dependency review](../docs/research/dependency-license-review-2026-11-02.md).
Versions/builds selected on another platform or Python version require their own
review. Listing a component here does not relicense it or constitute a source
offer. Original dvbfixer code remains under its own BSD-3-Clause grant; do not
describe the entire scientific environment as BSD-only.
