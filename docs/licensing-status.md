# Licensing: practical publication status

Updated 2026-11-02. This is an engineering release checklist, not a legal opinion.
**The audit is paused at the practical release boundary, not complete.**

## Short answer

Original DVBFixer code has an owner-approved **BSD-3-Clause** license. That grant
permits commercial use, modification and redistribution, including proprietary
derivatives, subject to its terms. It does not grant rights to third-party code,
force fields, model weights or restricted scientific software.

There is no evidence here that the entire current distribution is cleared for
all commercial/proprietary uses. Equally, unresolved questions are **not findings
of infringement** and do not make original DVBFixer code non-commercial.

## What can be published, and what needs attention?

| Material / action | Practical status | What to do |
|---|---|---|
| Original code contributed by the owner | **BSD grant in place** | Retain LICENSE, copyright and disclaimer. Contributor/copied-code rights remain their own question. |
| Entire Git repository or source archive | **Some bundled data needs clarification** | The repository includes force fields. Do not describe the whole tree as BSD-only or fully cleared merely because external program binaries are absent. |
| Current Python wheel | **Packaging checked; rights review incomplete** | It includes three FF trees and notices, but not external Python/native dependencies. Resolve bundled FF permissions and the GPL integration/distribution question before making blanket proprietary-distribution claims. |
| Installation instructions / dependency manifests | **Different from shipping a prebuilt environment** | Users obtain dependencies under their respective terms. This reduces what DVBFixer itself delivers; it does not automatically resolve direct GPL API integration or operator eligibility. |
| Complete environment, container or bundled desktop installer | **Not cleared by this audit** | Do not add a new prebuilt distribution on the strength of the root BSD license. Restricted software permission and relevant notices/source delivery must be resolved for that artifact. |
| Commercial calculations / hosted service | **No blanket prohibition from BSD or GPL alone** | GPLv2/v3 are not AGPL; ordinary server use alone does not impose AGPL-style server-source publication. MODELLER and other restricted software still need applicable use rights; distributing clients/binaries is a separate action. |

## Only the substantial open questions

| Component | What is known | Needed decision / clarification |
|---|---|---|
| **Open Babel / ACPYPE** | Direct Python API integrations; inspected core Open Babel grants name GPLv2 without or-later wording, ACPYPE declares GPLv3+. Mere co-installation does not establish an incompatibility. | One focused legal/upstream review of the actual combined distribution and the desired proprietary-derivative policy. Do not change backends or relicense the owner's code as a speculative workaround. |
| **MODELLER** | Separate academic/commercial terms; no redistribution agreement has been assumed. | Users need appropriate licenses. **Do not bundle it without applicable permission.** A DVBFixer license key setup is not redistribution permission. |
| **Bundled force fields** | Amber19SB has recorded GROMACS provenance/notices. CHARMM's official reference and local BUF edit history are identified; original acquisition remains unknown. | Clarify redistribution/modification rights for CHARMM and the AMBER99SB-ILDN/Lipid21 mixture, caps and BUF additions. Do not require the owner to remember the old download. |
| **Reduce** | Exact selected source/build identified; full notices and UCSF acknowledgement included. | Clarify the custom modified-versions-free condition if distributing a modified build; keep separate provenance for its transformed dictionary. This is not an established blanket commercial-use ban. |
| **Model weights** | Source-code and weight licenses differ; runtime pins matter. | Consult [model terms](model-licenses.md) for the exact backend. Do not ship checkpoints under the project BSD license or infer a later model version's rights. |

## Recommended next actions — not yet performed

1. **Finish the existing BSD source/wheel release review first.** Do not expand
   the task into certifying a complete container or scientific environment.
2. Send one concise request to the CHARMM/port maintainers: may the identified
   July-2022 files be redistributed in an open-source package and its commercial
   derivatives, and modified with attribution? Record the applicable terms.
3. Obtain a focused opinion on the Open Babel/ACPYPE combined-distribution
   boundary. The project owner's BSD choice alone cannot settle that question.
4. Keep MODELLER operator-licensed unless explicit bundling permission is obtained.
5. If a prebuilt environment is actually planned, resume per-library notices and
   corresponding-source work **for that selected artifact**, not every possible
   platform/build in advance.

No email has been sent, no upstream permission obtained, and no scientific
dependency, force-field content or distribution architecture changed by this
checklist. Excluding or downloading data separately would be a future approved
packaging decision, not an action already taken or an automatic legal solution.

## Evidence, if needed

- [Project license](../LICENSE) and [partial notices](../THIRD_PARTY_NOTICES/README.md).
- [Detailed dependency findings](research/dependency-license-review-2026-11-02.md)
  and [CHARMM provenance](provenance/charmm36-jul2022.md).
- Four selected dependency wheels inspected; selected Apple native-byte matches
  verified. Source correspondence/rights clearance is not inferred from hashes.
- All seven notice files checked in a freshly built bounded-snapshot sdist,
  wheel and isolated no-dependency install. This is packaging evidence, not
  certification of a full-repository source release or combined environment.
