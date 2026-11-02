# Licensing and distribution policy

Updated 2026-11-02. **The agreed licensing/documentation task is complete.**
This records the project's distribution choice, not a legal opinion or blanket
third-party license clearance.

## What DVBFixer distributes

- **Original code:** BSD-3-Clause, copyright 2026 Stanislav Skripkin.
  Commercial use and proprietary derivatives are permitted under that grant.
- **Primary delivery:** the source repository and existing Standard/Desktop
  installer. The installer obtains dependencies from their respective sources
  into the user's environment; it does not grant licenses to those dependencies.
- **Force fields:** the existing tracked sets remain included. Upstream headers,
  citations and notices are preserved; known local changes are identified in
  [bundled-data attribution](../THIRD_PARTY_NOTICES/FORCE-FIELDS.md).
  They are **not relicensed under the project's BSD grant**.
- **Wheel:** supported as a Python packaging mechanism, not a separate end-user
  installation solution or an additional publication target for this task.
  Existing wheel support and data inclusion are unchanged.
- **No prebuilt scientific environment or container is published as part of
  this task.** The Desktop profile means local installation/build, not delivery
  of a pre-populated environment. No new Docker/container pipeline is introduced.
- **MODELLER and model weights:** not copied into project release artifacts.
  Installation/provisioning and use remain subject to the selected component's
  own terms. Users obtain any required licenses themselves.

## What users need to know

| Component | Practical rule |
|---|---|
| DVBFixer original code | Retain the BSD copyright, license and disclaimer. No author notification is required for ordinary commercial use. |
| MODELLER | Obtain a license appropriate to academic/commercial/service use. A configured key is not permission to redistribute MODELLER. |
| Open Babel / ACPYPE | GPL obligations remain relevant to distributing combined applications. Installing dependencies separately is not an automatic legal exemption, and mere co-installation is not proof of incompatibility. Do not promise a BSD-only or freely closable complete stack. |
| Force fields and other data | Preserve their attribution and applicable terms. Widespread use in GROMACS/OpenMM supports normal scientific practice, but does not itself establish a license for every modified file. Unknown permission is not a finding of a prohibition. |
| Reduce | Preserve selected notices and the UCSF acknowledgement. Its custom condition for modified versions is not replaced by a generic BSD/Apache label. |
| Diffusion models | Use the [exact code/weight terms](model-licenses.md), not the project's BSD license or another model version's terms. |

GPLv2/v3 permit commercial use and are not AGPL. Ordinary server use alone does
not create an AGPL-style server-source publication obligation. Local container
construction and redistribution of a ready-made image are different actions;
containers are not categorically forbidden by GPL. Restricted software and
data terms still apply in the selected use/distribution context.

## Completion and limits

The owner has selected the license, attribution and source-plus-installer route;
no further owner decision or recovery of old download history is needed for this
task. License metadata, contributor terms, notices, model terms, provenance,
inventory maintenance and packaging checks are in place. Scientific code,
dependencies, force-field bytes and installer behavior were not changed to
work around a license question.

Unresolved historical FF permissions and GPL combined-work interpretation remain
recorded in the [evidence review](research/dependency-license-review-2026-11-02.md).
Completion does **not** mean those third-party rights have been adjudicated.
No upstream permission, legal opinion, completed exhaustive SBOM or source offer
is invented. No emails were sent. Further forensic research is not an active
deliverable or a request for action from the owner.

If the project later distributes a prebuilt environment, restricted binaries,
checkpoints or a proprietary combined application, assess that concrete artifact
and its notice/source/permission requirements before publishing it. Reassess
changed dependency/model/data versions; do not inherit clearance from old hashes.

## Maintenance checklist

For a dependency, FF, model or asset update: preserve upstream notices, record
the actual version/source and known modifications, update model/data terms where
needed, regenerate `scripts/gen_license_inventory.py`, and run its `--check`
plus the license tests. Never label a scanner result as legal approval. A future
release artifact must carry its declared notices; an optional wheel is not an
exception to attribution requirements. No PyPI/GitHub release is made by closing
this task, and no merge into `main` is implied.

See [LICENSE](../LICENSE), [notices](../THIRD_PARTY_NOTICES/README.md),
[installation](installation.md) and the
[completion record](research/licensing-completion-2026-11-02.json).
