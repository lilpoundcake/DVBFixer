# CHARMM36 LJ-PME July 2022: source identification and open rights review

Updated **2026-11-02**, review baseline `6dae875`.

[Content comparison](charmm36-jul2022-comparison.json) ·
[Dependency review](../research/dependency-license-review-2026-11-02.md) ·
[Licensing plan](../plans/licensing-and-attribution.md)

## Acquisition versus content reference

The owner reports that `FF/charmm36_ljpme-jul2022.ff/` was obtained **from
GROMACS**. Its exact GROMACS release, package, archive or installed source path
has not yet been identified. Preserve that report: this investigation does not
claim the owner downloaded the files directly from MacKerell.

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
found, but not necessarily who assembled or licensed every file. Do not reject
the owner's report or silently transfer the checked release's LGPL to a
different bundle.

## Remaining permission evidence

No file named LICENSE/COPYING/README/copyright was found in the reference tar,
and a bounded text search found no explicit license/redistribution grant there.
Its `forcefield.doc` supplies converter, parameter-version and citation details,
not a complete grant for every portion. This search is not a legal determination
that no permission exists elsewhere.

Next evidence needed: the owner's actual GROMACS package/archive/version or
source path, its associated notices, and the provenance of BUF modifications.
Then evaluate the upstream parameters, conversion and local changes separately.
Until then redistribution clearance is **unknown**, not BSD-3-Clause by virtue
of the project license and not LGPL merely because the data are in GROMACS format.
