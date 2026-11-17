# Project goals and scope

## Purpose

**dvbfixer makes preparation of biomolecular structures for downstream
computational analysis convenient, reproducible and independently checkable.**
It identifies structural defects, completes missing structural information and
reconciles chemical representations while preserving data provenance and making
uncertainty explicit.

This page defines project-level goals, not a claim that every input class is
supported or that every research method is qualified for production. Current
capabilities and limitations are documented in the [command guides](commands/),
[domain model](domain-model.md) and [known issues](known-issues.md).

## Goals

1. **Make structural preparation accessible.** Provide consistent CLI workflows
   and a GUI for inspecting inputs, selecting operations and reviewing results.
2. **Detect and explain defects.** Report missing atoms/residues, connectivity,
   stereochemical and steric problems, and distinguish inherited findings from
   defects introduced during preparation.
3. **Complete structures from explicit evidence.** Rebuild missing atoms and
   known-sequence residue gaps, and support comparative modeling with explicit
   target sequences and template provenance. Never silently invent sequence to
   resolve an ambiguous construct or mapping.
4. **Prepare chemical representations for downstream tools.** Support protonation,
   hydrogen placement, naming conversion, restrained minimization and topology
   export within the supported chemistry and force-field contracts.
5. **Preserve molecular meaning and traceability.** Retain chain/residue/atom
   identities, relevant context and connectivity; distinguish observations from
   modeled additions and record any explicitly permitted changes.
6. **Make results reproducible and reviewable.** Record inputs, settings, tool
   versions, computational ownership and validation findings. Preserve failures
   and explain unsupported cases rather than silently changing the task.

## Experimental structure completion research

The current reconstruction study uses experimentally determined structures from
public databases, principally the PDB. It estimates absent coordinates for
residues whose identities and order are established by the deposited construct's
sequence information. The observed scaffold and relevant context are preserved
under explicit task contracts.

This is **coordinate completion at a known sequence, not new protein sequence
design**. Modeling tools can generate alternative coordinate proposals, but those
coordinates remain computational estimates, not experimental measurements.
Geometry checks alone cannot establish the true conformation of an unresolved
region.

Test-set curation and model execution are separate operations. Structural coverage,
data quality, construct identity and reproducible grouping/splitting guide the
test set; it is not selected to optimize biological activity. Public availability
establishes source accessibility, not automatic suitability for every operation.
The current study's scope and execution boundaries are in the
[experimental structure completion handoff](plans/experimental-structure-completion-handoff.md).
They do not restrict every other supported project workflow to that study's inputs.

## Non-goals and limitations

- De novo sequence design, optimization of biological activity and experimental
  synthesis/expression/testing are not project goals. An explicitly supplied
  target sequence is an input to modeling, not an instruction to optimize function.
- Preparation does not certify experimental truth, biological function, correct
  topology or MD readiness for arbitrary structures, ligands and cofactors.
- A clean geometry report, successful software test or isolated research example
  does not establish broad scientific qualification.
- Research prototypes do not automatically become public capabilities or change
  production defaults. New scope and execution require their own explicit decisions.
- Existing tools are reused rather than assumed interchangeable. Open-source code,
  model-weight terms and redistribution rights are distinct; MODELLER and other
  dependencies retain their own licensing requirements.

## Quality criteria

A successful workflow has an explicit input and task, preserves required identity
and context, accounts for every modeled addition and permitted change, and passes
the applicable independent checks. Remaining warnings, unsupported chemistry,
uncertainty and incomplete preparation remain visible to the user.

See [architecture](../ARCHITECTURE.md) for implementation boundaries,
[pipelines](pipelines.md) for usage recipes and the
[plan status index](plans/README.md) for implementation and research progress.
