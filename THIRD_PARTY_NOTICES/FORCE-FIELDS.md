# Bundled force-field attribution

These scientific data retain their upstream rights and attribution. The
DVBFixer BSD-3-Clause grant does not relicense the data. The project retains the
existing sets for its source-plus-installer distribution; this decision is not
a new permission grant from their authors. No force-field parameters were
changed as part of licensing work.

| Material | Attribution and known changes |
|---|---|
| `FF/amber19sb.ff/` | Unmodified GROMACS v2026.3 copy. Preserve `forcefield.doc` and `GROMACS-LGPL-2.1.txt`; exact source/checksums are in the repository's `docs/provenance/amber19sb-gromacs.md`. |
| `FF/amber99sb-ildn-lipid21.ff/` | Existing mixed AMBER99SB-ILDN protein / AMBER94 nucleic / Lipid21 data. Its `forcefield.doc` credits Kresten Lindorff-Larsen and coworkers at D. E. Shaw Research, and Dickson et al. for Lipid21. The exact mixed-port acquisition/permission chain is not established; do not attribute the entire mixture to the project's BSD license or assume it is an untouched upstream release. |
| `FF/charmm36_ljpme-jul2022.ff/` | CHARMM36/CGenFF 4.6 GROMACS conversion by András Wacha, based on E. Prabhu Raman, Justin A. Lemkul, Robert Best and Alexander D. MacKerell, Jr. Preserve the per-file headers and `forcefield.doc`, including its extensive source/citation list. |
| `src/dvbfixer/data/charmm-terminal-caps.xml` | Explicit-residue adaptations of CHARMM36 ACE/CT3 terminal patches to DVBFixer's separate ACE/NME representation, as documented in the XML comment. CHARMM-derived parameter attribution remains separate from the surrounding application code. |

CHARMM's initial 62-file import matches the corresponding official archive
files. Later local BUF changes affect `atomtypes.atp`, `ffnonbonded.itp` and
`ions.itp`: commits `8bb94b4` and `1bbd3bf`, dated 2026-03-12, author signature
`lilpoundcake`. The reference also contains `watermodels.dat`, absent at initial
import. Original acquisition and original derivation of BUF values remain
unknown. See `docs/provenance/charmm36-jul2022.md` in the source repository for
full commit identities, archive comparison and evidence limitations.

The identified charmm2gmx converter revisions use BSD-3-Clause; licensing the
conversion program does not automatically assign that license to its input
parameters or every output. No new license is inferred from GROMACS/OpenMM
format, public download access, a publication citation or common usage alone.

Recommended scientific citations (not substitutes for legal terms):

- Lindorff-Larsen et al., *Improved side-chain torsion potentials for the Amber
  ff99SB protein force field*, Proteins 78:1950–1958 (2010), as recorded upstream.
- Dickson et al., JCTC (2022), Lipid21, as identified by this bundle's
  `forcefield.doc`; consult the upstream material for the full reference.
- The component-specific CHARMM references listed in its `forcefield.doc`.
- Wacha and Lemkul, charmm2gmx conversion method, DOI `10.1021/acs.jcim.3c00860`.

Redistribution of a new modified set or an independently assembled scientific
environment requires assessment of that actual artifact. This attribution page
does not declare any unidentified material public domain or establish a ban on
normal calculations.
