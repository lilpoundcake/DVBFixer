# Diffusion Confirmatory Cohort Materialization

`scripts/materialize_diffusion_confirmatory_cohort.py` converts the metadata-locked
cases in `diffusion-confirmatory-cohort.json` into screened, runnable benchmark
workspaces. The default output is below the ignored `.artifacts/` root.

```bash
python scripts/materialize_diffusion_confirmatory_cohort.py
python scripts/materialize_diffusion_confirmatory_cohort.py --case-id 36hb
```

The script downloads only the case's exact
`https://files.rcsb.org/download/<PDB>.cif` URL. It records and verifies the
download SHA256 but does not treat transport metadata as coordinate evidence.
Gemmi performs case-sensitive label/auth/entity selection. The selected chain is
then converted through `dvbfixer.structure_input.normalize_structure`, and the
resulting request must pass the existing diffusion scope gate.

## Stable Layout

Layout version `dvbfixer-diffusion-confirmatory-workspace-v4` is:

```text
<root>/
  downloads/<PDB>.cif
  downloads/<PDB>.source.json
  cases/<pdb-id>/
    .complete.json
    accepted.json | rejected.json
    workspace/                    # accepted cases only
      input/normalized.pdb
      reference.pdb
      target.fasta
      request.json
      source.json
  report.json
```

Case IDs are lowercase PDB IDs. Chain IDs remain case-sensitive everywhere.
`report.json` gives accepted workspaces as root-relative paths such as
`cases/36hb/workspace`; a batch runner must enumerate accepted records rather
than infer workspaces from directory names.

Workspace PDBs are heavy-atom-only. Explicit hydrogens present in an experimental
mmCIF are removed deterministically after the shared CIF boundary and before
reference/source/request construction. `input/normalized.pdb` therefore has
exactly the identities in `DiffusionRequest.fixed_atoms`, while `reference.pdb`
has exactly the union of `fixed_atoms` and `generated_atoms`. Backends place
hydrogens after consuming this contract, matching the established pilot policy.

## Publication And Resume Rules

A case is complete only when `.complete.json` exists and verifies its case
record. Accepted records also carry SHA256 and byte length for every workspace
artifact. Resume re-hashes all artifacts before reusing a case. Scientific
rejections receive an atomic marker; download, filesystem, and unexpected
software failures are `operational-error` report entries without a marker and
are retried on the next invocation.

Filtered `--case-id` runs preserve completed records for other cases and report
unselected unfinished cases as `pending`. `report.json` schema version 1 has:

- `layout_version`, `cohort_sha256`, mask seed, and context radius
- `summary` counts for `accepted`, `rejected`, `operational-error`, and `pending`
- ordered `cases` records matching cohort order
- accepted-case source identity/digests, artifact digests, and exact mask identities
- rejected-case stable `reason_codes`, detail, and structured evidence
- minimum-final-accepted, screening-complete, and confirmatory-sample-valid decisions

The mask rank is a SHA256 ordering over immutable case metadata and a fixed
selection seed. Candidates must be internal, have the planned length (5 or 10),
have two complete observed anchors, and have no unsupported local context.

The RCSB entity sequence is retained only for clustering and provenance. The
request target and FASTA use `observed-coordinate-sequence`: canonical residues
with unique, strictly increasing `label_seq_id` placement. Natural missing
entity positions are omitted and are never generated. A mask and its two
anchors must form a contiguous entity-sequence run, so the artificial withheld
region cannot span a natural coordinate gap.

The frozen metadata pool contains 500 ordered sequence-cluster/PDB
representatives, balanced as 250 planned 5-residue masks and 250 planned
10-residue masks. All 500 receive accepted or rejected screening records before
inference. In locked `screening_index` order, the first at most 300 accepted
cases are selected for inference; later accepted cases are retained as
`reserve-not-selected` and are not runnable. At least 180 selected independent
groups are required. Screening does not stop when that minimum is reached.

## Live Regeneration Handoff

The repository currently retains the superseded 200-case JSON until an agent
with network access regenerates it. Run from the repository root:

```bash
python scripts/build_diffusion_confirmatory_cohort.py \
  docs/research/diffusion-confirmatory-cohort.json \
  --groups 500 --candidates 800
sha256sum docs/research/diffusion-confirmatory-cohort.json
```

Replace `REGENERATE_AFTER_LIVE_RCSB_BUILD` in
`diffusion-gap-reconstruction-inventory.toml` with that digest, change the
inventory status to `candidate-metadata-locked-coordinate-screening-pending`,
and synchronize `matching_sequence_cluster_count` and
`returned_representative_count` from the regenerated JSON. Then run the full
blind coordinate screen without `--case-id`:

```bash
python scripts/materialize_diffusion_confirmatory_cohort.py \
  --output-root .artifacts/diffusion-confirmatory-cohort-500-v4
```

Do not begin inference unless `report.json` says `screening_complete: true` and
`confirmatory_sample_valid: true`. The report must contain all 500 ordered case
records even when the accepted count reaches 180 earlier. The runner verifies
`selected_for_inference: true`; reserve workspaces and completion markers may
exist, but runner discovery excludes them.
