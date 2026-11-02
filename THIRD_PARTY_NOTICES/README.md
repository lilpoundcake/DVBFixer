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

### Reduce 4.15 acknowledgement and source notices

Reduce includes software developed by the **University of California,
San Francisco**. This acknowledgement is not an endorsement by the University.
The identified source release also credits J. Michael Word, ReliaSolve,
Richardson Lab at Duke University, and Ralf W. Grosse-Kunstleve in the respective
source portions; do not collapse their notices into a generic BSD label.

The four complete upstream directory/root license texts are reproduced in
`REDUCE-4.15-ROOT.txt`, `REDUCE-4.15-LIBPDB.txt`, `REDUCE-4.15-SOURCE.txt` and
`REDUCE-4.15-TOOLCLASSES.txt`. Each preserves the original text with one terminal
LF appended for repository formatting; no notice wording is changed.
`REDUCE-4.15-ADDITIONAL.txt` preserves inspected source-specific notice blocks.
These notices document an external dependency; no Reduce executable, dictionary
or source archive is added to dvbfixer's own wheel by these files.

Identified source: Reduce tag `v4.15`, commit
`e748530212b104d2e3328bcb98d6f6972ff2997b`, archive SHA-256
`f2f993e3f86ded38135d6433e0a7c2ed10fbe5da37f232c04d7316702582ed06`.
These texts are not permission to disregard the "modified versions ... available
for free" condition or a source offer. Exact condition interpretation, any
intended modified build, complete copyright closure and separately derived
dictionary-data rights remain under review.

Identified subcomponents of those external wheels (not additional software
bundled into dvbfixer's own wheel):

- **InChI 1.04:** the reviewed Open Babel wheel headers name IUPAC/InChI-Trust
  InChI Licence No. 1.0 or later, not an assumed current MIT grant. The historical
  text has notice/source/linking and modified-output naming requirements.
- **Cairo 1.15.14 source reference:** the Linux recipe names this release; its
  implementation offers LGPL 2.1 or MPL 1.1. Auxiliary tools have separate terms.
  Exact wheel build/source correspondence and Apple's Cairo version remain open.
- **ACPYPE charmmgen:** separate legacy component with September 2026 rebuild
  history and a patched AmberClassic build recipe. Matching binary containers
  have been identified; exact source revision and applicable source closure are
  not established. Do not infer the component's terms from ACPYPE's wrapper alone.

Exact artifact identities, retained full-text hashes and unresolved questions:
[dependency review](../docs/research/dependency-license-review-2026-11-02.md).
Versions/builds selected on another platform or Python version require their own
review. Listing a component here does not relicense it or constitute a source
offer. Original dvbfixer code remains under its own BSD-3-Clause grant; do not
describe the entire scientific environment as BSD-only.
