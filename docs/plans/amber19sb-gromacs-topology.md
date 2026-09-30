# Proposed `amber19sb` GROMACS topology backend

**Status:** implemented on 2026-09-30 after maintainer approval of all gate
items. Implementation evidence, the selected `v2026.3` provenance, checksums,
metadata audit, exact water/ion matrix, redistribution assessment, and external
acceptance status are recorded in
[`../provenance/amber19sb-gromacs.md`](../provenance/amber19sb-gromacs.md).

## Objective and scope

Add a distinct, **protein-only** `dvbfixer top --ff amber19sb` RTP backend.
It must coexist with, and never change the meaning, default, files, water/ion
selection, or results of the current `top --ff amber` backend
(`amber99sb-ildn-lipid21.ff`).  `amber19sb` is not an OpenMM alias and is not
a replacement for the separate OpenMM `amber19` / `amber` aliases.

The backend must reject unsupported non-protein chemistry before writing a
partial topology.  In particular, it does not promise Lipid21, GLYCAM,
carbohydrate, nucleic-acid, arbitrary-ligand, or mixed-complex support.  Those
systems retain their current backend choices (including `amber`, `charmm`, or
`--acpype`) until separately designed and validated.

## Source provenance and update rule

The candidate local source tree is `FF/amber19sb.ff`.  The prior byte-for-byte
comparison found it identical to the official GROMACS source at
[`github.com/gromacs/gromacs` tag `v2026.1`](https://github.com/gromacs/gromacs/tree/v2026.1),
commit `f41e7fae…`.  At plan writing, that is the newest stable `v2026` tag.
This is provenance for the candidate only; it is not a claim that DVBfixer
currently ships the tree.

At implementation time, the owner must query GitHub releases/tags and identify
the **latest stable** GROMACS tag (exclude release candidates, beta releases,
and other prereleases).  Do not blindly import `v2026.1` merely because this
plan names it.  Record the selected tag, fully resolved commit, retrieval date,
upstream file manifest with SHA-256 values, and the exact comparison command in
the implementation evidence.  If the current latest stable tag differs from
`v2026.1`, repeat the complete compatibility, water/ion, license, and GROMACS
validation below against that tag.  Any mismatch from an official checkout is
a blocker until reviewed; do not hand-edit upstream parameter files to make the
comparison pass.

## Why this is a separate backend

| Aspect | Existing `top --ff amber` | Proposed `top --ff amber19sb` |
|---|---|---|
| Bundled directory | `amber99sb-ildn-lipid21.ff` | official `amber19sb.ff` from the latest stable GROMACS tag |
| Intended coverage | Current AMBER/RTP route, including its existing Lipid21 and water/ion policy | Canonical protein-only upstream distribution |
| Compatibility promise | Preserved exactly | New, explicitly bounded capability |
| Water and ions | Current DVBfixer water selector plus maintained ion-set substitution | Only upstream-supported water/ion combinations, enumerated and validated pair-by-pair |
| Default | `amber` | Never default; explicit opt-in |

## Capability descriptor and admission policy

Introduce one declarative force-field capability descriptor, consumed by CLI
validation, pipeline admission, writers, help generation, and tests rather
than scattered `args.ff == "amber"` / `"charmm"` tests.  It should describe at
least:

- CLI name, bundled directory, family, provenance identifier, and whether the
  backend is selectable as a default;
- accepted molecule classes (`protein` and explicit supported protein caps only
  for `amber19sb`), rejected classes, supported protonation/terminal handling,
  and whether CMAP, ARN, R2B, and TDB inputs are available;
- exact upstream water models and exact ion molecule/parameter sources;
- the allowed **water × ion** pairs, including whether each is supplied by
  upstream files or a separately reviewed DVBfixer transformation;
- whether custom `--ff-dir`, `--ion-set`, `--merge`, `--acpype`, and
  topology-matched PDB output are meaningful for that descriptor.

The descriptor must drive a preflight classification before output files are
opened.  For `amber19sb`, reject mixed systems and unsupported residues with a
clear error naming the component and a supported alternative.  It must not
silently omit chains or inherit `amber`'s Lipid21/glycan behavior.  Preserve
case-sensitive chain IDs and `(chain, resid, icode)` identity in the resulting
diagnostics and topology path.

## Implementation sequence

1. **Revalidate and stage the upstream distribution.** Fetch the latest stable
   tag as described above, compare every staged file byte-for-byte, retain
   upstream notices, and add the complete directory only after approval.  Add
   it to source and wheel package data alongside—not in place of—the existing
   AMBER and CHARMM directories.  Ensure source checkouts and installed wheels
   resolve the same descriptor root.
2. **Parse the upstream metadata without AMBER99 assumptions.** Audit the
   parser/builder against the candidate's `forcefield.itp`, `ffnonbonded.itp`,
   `ffbonded.itp`, `aminoacids.rtp`, `.r2b`, `.arn`, terminal databases and
   hydrogen databases.  Add a descriptor-controlled file inventory check.
   Absence of an optional file must be represented explicitly, not guessed.
3. **Handle ARN, termini, and CMAP correctly.** Confirm all ARN directions
   needed for input/output atom spelling; test N- and C-terminal blocks and
   supported ACE/NME caps.  Determine from the upstream files whether CMAP is
   present and required.  If it is, generalize the writer so AMBER19SB emits
   its upstream CMAP parameters and the builder emits matching `[ cmap ]`
   terms; if not, test and document that no CMAP path is selected.  Do not use
   the current `ff_type == "charmm"` condition as an implicit CMAP policy.
4. **Validate every upstream water/ion pair.** Derive the matrix from the
   selected upstream directory, including all water models it actually ships
   and every compatible ion definition/parameter pairing.  For every allowed
   pair, test that `top` accepts it, emits self-contained water and ion ITPs,
   has no missing atom type or molecule type, and passes GROMACS preprocessing.
   For every unsupported cross-pair, reject before writing.  Do not apply the
   existing hand-maintained AMBER99 `ION_PARAMS` substitutions to Amber19SB
   unless an approved provenance-preserving design explicitly establishes that
   they are compatible; no fallback to a merely similarly named water model.
5. **Keep current routing intact.** Add `amber19sb` as a third explicit CLI
   value and route directory selection, diagnostics, builder settings, and
   writers exclusively through the descriptor.  Keep `amber` mapped to
   `amber99sb-ildn-lipid21.ff`, default `amber`, and its present ion behavior.
   `--acpype` remains its documented AMBER14+GLYCAM path and must not acquire
   Amber19SB semantics.
6. **Update user-facing and package evidence.** Regenerate CLI reference and
   GUI command schema after argparse changes; update the force-field matrix and
   `top` command documentation to distinguish the two GROMACS AMBER choices,
   protein-only limits, provenance, and the validated water/ion matrix.  Add
   required license/copyright notices and attribution exactly as required by
   the selected upstream distribution and this repository's distribution
   policy; legal review must confirm redistribution terms before publishing a
   wheel.

## Verification plan

- Add unit tests for descriptor validation, CLI choices/default preservation,
  bundle resolution in source and installed-wheel contexts, file inventory,
  ARN mappings, terminal/cap handling, and any CMAP branch.
- Add protein fixtures spanning standard residues, disulfide handling,
  histidine/protonation variants, termini, ACE/NME caps, and an unsupported
  non-protein input that fails before artifacts are created.
- For each allowed upstream water/ion pair, generate a minimal protein-plus-
  water/ion topology and run `gmx grompp` using an appropriate minimal MDP and
  coordinates.  Assert zero missing includes/types/moleculetypes and retain
  the command output as CI diagnostics on failure.  Run negative tests for
  every disallowed pair.
- Run regression tests proving existing `top --ff amber` output selection and
  its current supported matrix remain unchanged; test `charmm` separately to
  protect its current restrictions.
- Run the topology task's focused tests, the new parameterized tests, generated
  CLI/GUI checks, `python scripts/check_agent_docs.py`, and `git diff --check`.

## Explicit approval gate

**No implementation may begin** until a maintainer approves all of the
following evidence in one reviewable change proposal:

1. latest-stable-tag lookup, full commit, byte-identical manifest comparison,
   and provenance record;
2. the completed capability descriptor and unsupported-chemistry policy;
3. the complete upstream-derived water × ion matrix and a decision for every
   pair (supported or preflight-rejected);
4. ARN, terminal, CMAP, and ion-file compatibility audit results;
5. packaging and license/redistribution review; and
6. the proposed fixture and `gmx grompp` CI strategy, including availability of
   a pinned GROMACS executable.

Failure of any gate item means retain the present `amber` route unchanged and
do not add `--ff amber19sb` or ship `FF/amber19sb.ff`.
