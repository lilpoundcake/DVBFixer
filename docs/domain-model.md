# Scientific domain model

DVBfixer uses a focused domain model where file-format details must not
silently change molecular meaning. It is intentionally not a rewrite of the
CLI, GUI, or external-tool adapters.

## Structure identity

`dvbfixer.domain.structure_identity` defines immutable atom, residue, and
molecular-component references, component kinds, and chain-ID allocation. A
generated PDB chain ID must be unique and must not reuse a polymer or
metadata-reserved ID. PDB has only 62 portable one-character IDs, so exhaustion
is reported rather than producing ambiguous output.

`split --unique-molecule-chains` applies that policy. After empirical chain
assignment, it preserves the resulting polymer chain IDs, assigns each non-solvent heterogen component a new ID, rewrites
coordinate identity records, leaves CONECT records unchanged, and records the
operation in `REMARK 999`. Waters and monatomic ions are not separately named
molecules.

`--keep-heterogens` retains all heterogen records, including solvent, ions,
and buffers. With this option atom serials remain unchanged and solvent stays
inside its source MODEL. Solvent/ion records do not consume unique molecule IDs.

For legacy PDB, a residue is currently the component boundary. Empirical Split preserves
atom serials in all modes so unchanged CONECT records retain their endpoints.
PDBx/mmCIF `label_asym_id` remains the authoritative identity whenever an input
adapter can carry it losslessly.

## Parameterization policy

`dvbfixer.domain.parameterization` expresses five routes: exact native template,
isolated organic GAFF candidate, user template, unsupported complex cofactor,
and explicit stripping.

A residue name alone is never proof of force-field compatibility. Exact atom
and bond compatibility is checked by OpenMM. Complex cofactors such as HEM and
FAD are not sent through generic isolated-ligand GAFF because coordination,
covalent attachment, and redox state require a validated model. Supply one with
repeatable `--extra-ff FILE`, or explicitly use `--strip-heterogens`.

Small organic ligands remain GAFF candidates. Charge corrections are based on
connectivity, including tetracoordinate quaternary ammonium and ionizable
phosphate groups; they do not depend on LBN/POPC atom names. Chemically similar
deposited lipids can have incompatible naming and must not be blindly renamed
to a force-field template.

Gemmi, OpenMM, Open Babel, OpenFF, and AmberTools are adapters around these
decisions, not part of the domain model.

## Integration limits

DDD is partially integrated: chain allocation and the complex-cofactor routing
guard are active callers of the domain policies. The five route names are a
policy vocabulary, not five completed backend implementations. Unknown-residue
screening currently uses the loaded template names; OpenMM subsequently checks
atom/bond compatibility. The presence of an XML template name therefore does
not establish that an incomplete or differently bonded residue can be used.

There is no general component-completeness validator, metal/redox-state resolver,
or domain service for whole-complex geometry regularization. PDB multi-residue
heterogens are not yet grouped by connected-component analysis for naming.
Engine setup and fallback behavior still live in `minimize/pipeline.py` and
`lig_params.py`. See [relaxation research](research/whole-complex-relaxation.md)
for a proposed architecture and validation criteria; none of those proposed
backends are exposed by this release.
