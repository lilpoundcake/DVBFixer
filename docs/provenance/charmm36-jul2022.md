# CHARMM36 LJ-PME July 2022: source identification and open rights review

Updated **2026-11-02**, review baseline `6dae875`.

[Content comparison](charmm36-jul2022-comparison.json) ·
[Dependency review](../research/dependency-license-review-2026-11-02.md) ·
[Licensing plan](../plans/licensing-and-attribution.md)

## Acquisition versus content reference

The original acquisition source is **unknown**. The owner withdrew the earlier
GROMACS hypothesis as uncertain on 2026-11-02. Neither GROMACS acquisition nor
direct MacKerell download is confirmed. Continue independent research without
requiring the owner to remember an archive or release.

An independent **content reference** was retrieved from the official MacKerell
[force-field page](https://mackerell.umaryland.edu/charmm_ff.shtml):

- [July 2022 LJ-PME archive](https://mackerell.umaryland.edu/download.php?filename=CHARMM_ff_params_files/charmm36_ljpme-jul2022.ff.tgz)
- Retrieved: 2026-11-02, using curl after Python's TLS handshake timed out.
- Bytes: `1103301`.
- SHA-256: `87f95a03b83d57e6b1965eb4ab798136bc7b78270f8fb85a77f82596d496c15c`.
- Retained locally: `tmp/licensing-audit-2026-10-31/charmm36_ljpme-jul2022.ff.tgz`.

The download URL is mutable; the hash identifies the compared bytes, not every
past or future archive served there. Neither matching content nor public
download access establishes the owner's acquisition chain or a redistribution
grant.

## Byte comparison

The reference archive contains 63 regular files; the repository tracks 62 files
in this tree. After removing their different parent-directory prefixes:

- **59 files are byte-identical**.
- **Three files differ:**
  - `atomtypes.atp`: adds the `BUF` buffer-ion entry.
  - `ffnonbonded.itp`: adds the `BUF` atom type and BUF/HT, BUF/OT pair entries.
  - `ions.itp`: adds the `BUF` molecule definition and optional position restraint.
- `watermodels.dat` exists in the reference archive but is not tracked locally.
- No tracked local filename is absent from the reference archive.

The machine-readable comparison retains both manifests and each differing
file's hashes. Changes are observations, not attribution to a particular author
or a finding about scientific correctness. No parameters were changed, replaced
or added during this review. In particular, the original author/source and
permission for the BUF additions remain unverified.

## Recovered history and converter headers

At initial import `c7714fc61c1dfeb0fa16514764bf909cc0fde6a8` (2026-03-10),
all 62 imported files match the corresponding official reference bytes.
`watermodels.dat` was already absent. Commit
`8bb94b4f6e7514a31de74b458e7086da66393982` (2026-03-12) introduces dummy BUF
definitions; `1bbd3bffd75d85eb9a7070b57bcfdfc85f32253a` later that day replaces
them and adds the current BUF atom/pair values. Both have author signature
`lilpoundcake`. The latter message calls the values the user's parameters, but
does not establish their original derivation or rights.

The quoted header names `0.6.dev1+g077b6a1`, dated 2021-10-05, and an empty
missing-dihedral file. The tracked `ffmissingdihedrals.itp` instead names
`0.7.dev45+g7b82040.d20221208`, dated 2022-12-08, and is nonempty. It matches the
reference archive. These are different versions, not interchangeable evidence.

Converter revisions resolve to `077b6a13ae3233d5ad08c9910a3f33d1ea1bdf05` and
`7b820405b929f1bd7d1e749ffa742061f87213a9`. Both historical LICENSE texts are
BSD-3-Clause, copyright 2021 András Wacha. The
[converter documentation](https://awacha.gitlab.io/charmm2gmx/) identifies the
tool as producing the canonical GROMACS CHARMM port since July 2021; the
Wacha/Lemkul paper (DOI `10.1021/acs.jcim.3c00860`) describes the conversion.
**Converter licensing does not automatically license the input parameters or
all generated output.** Parameter/port permissions remain a separate review.

## Reproducing the current comparison

Reproduce without extraction or modification (use a fresh output path):

```bash
python scripts/compare_forcefield_archive.py \
  --archive tmp/licensing-audit-2026-10-31/charmm36_ljpme-jul2022.ff.tgz \
  --expected-sha256 87f95a03b83d57e6b1965eb4ab798136bc7b78270f8fb85a77f82596d496c15c \
  --local-tree FF/charmm36_ljpme-jul2022.ff \
  --archive-prefix charmm36_ljpme-jul2022.ff \
  --output tmp/licensing-audit-2026-10-31/charmm-comparison-repeat.json
```

## The checked GROMACS reference is not proof of the owner's package

The official GROMACS `v2026.3` source at commit
`121090014570a53a17ea391bcddae45e5ea05eb4` was inspected separately. Its
`share/top` directory lists `charmm27.ff`, but not this Jul-2022 tree. The listing
was captured from an immutable commit reference in the upstream source register.

This is **one checked release**, not a statement that no GROMACS package,
third-party force-field collection, custom build or installation can contain
the tree. An installed `share/top` location establishes where the files were
found, but not necessarily who assembled or licensed every file. The owner's
GROMACS hypothesis is withdrawn; do not transfer that release's LGPL to a
different bundle.

## Remaining permission evidence

The **original CHARMM-format** `toppar_c36_jul22.tgz` was also retrieved from
the same official download page on 2026-11-02, separately from the GROMACS port.
It is 11,998,126 bytes, SHA-256
`dcb5c266bc35b1b17829a48b1fb8e058822262916b6d28bde603c2067b1b38c4`, retained
under `dependency-review-2026-11-02/toppar_c36_jul22-original.tgz` below the audit
evidence root. A bounded text search of regular members up to 4 MB did not locate
a general redistribution grant; copyright matches in example outputs/binaries
are not such a grant. Nested archives and complete file-level permissions still
need review; do not present this search as exhaustive clearance or a prohibition.
The official page links the upstream
[original parameter repository](https://github.com/Lemkul-Lab/toppar_c36_releases),
whose currently inspected root has no LICENSE and whose README identifies
maintainers but supplies no general grant. That current repository observation
is not a retroactive license assignment to July 2022.

No file named LICENSE/COPYING/README/copyright was found in the reference tar,
and a bounded text search found no explicit license/redistribution grant there.
Its `forcefield.doc` supplies converter, parameter-version and citation details,
not a complete grant for every portion. This search is not a legal determination
that no permission exists elsewhere.

Continue reviewing upstream parameter/port permissions and BUF-value derivation
using the identified content and edit history. Unknown acquisition does not
prevent this independent review.
Until then redistribution clearance is **unknown**, not BSD-3-Clause by virtue
of the project license and not LGPL merely because the data are in GROMACS format.
