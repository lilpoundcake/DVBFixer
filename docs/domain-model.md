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

`split --unique-molecule-chains` applies that policy. It preserves polymer
chain IDs, assigns each non-solvent heterogen component a new ID, rewrites
coordinate identity records, leaves CONECT records unchanged, and records the
operation in `REMARK 999`. Waters and monatomic ions are not separately named
molecules.

For legacy PDB, a residue is the conservative component boundary because the
existing split writer reserializes coordinates while preserving CONECT text.
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
