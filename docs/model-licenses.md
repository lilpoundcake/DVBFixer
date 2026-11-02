# Model source, weights and data terms

Updated **2026-11-02**. This document records third-party terms; it does not
replace the original-code [BSD-3-Clause license](../LICENSE) or approve
redistribution of a complete runtime.

| Integration | Source | Weights | Operational boundary |
|---|---|---|---|
| Protpardelle-1c cc94 | MIT at pinned engine revision | CC-BY-4.0, official Zenodo record 16817230 | New requests; CUDA/MPS runtimes are separate from Standard |
| Protpardelle-1c cc89 | Same pinned MIT source | CC-BY-4.0 for the corresponding official archive | Explicit historical/frozen-request replay, not the new default profile |
| Protenix v1 | Apache-2.0 at pinned engine revision | Pinned README includes v1 parameters under Apache-2.0 | Separate optional CUDA runtime; do not extend this finding to v2 weights |
| Boltz-2 | MIT in inspected pinned repository | Pinned README/model card includes weights under MIT | Research-only; `mols.tar` and auxiliary data require separate review |
| RFdiffusion v1 | BSD-3-Clause at pinned runner source | Later upstream clarification includes linked weights | Clarification commit differs from pinned execution source; research scope |

Exact checkpoint names, SHA-256 identities and pinned source links live in the
[audit's model inventory](research/third-party-license-audit.md#model-source-weights-and-auxiliary-artifacts).
The [upstream source register](research/licensing-upstream-sources.json) and
[captured evidence index](research/licensing-evidence-2026-10-31.json) bind the
current legal observations. The [decision report](research/licensing-decisions-2026-10-31.md)
states remaining release gates. Do not infer terms for different checkpoints,
future versions or all contents of an upstream cache.

## Preserve separate provenance

Protpardelle's shared adapter now reports `source_license="MIT"` and
`checkpoint_license="CC-BY-4.0"` separately, after existing digest verification.
This metadata correction does not change sampling or confer scientific acceptance.
Existing frozen requests/results are not rewritten. Modified engine source must
retain its notices and describe local patches if distributed.

If CC-BY-4.0 weights or adaptations are shared, preserve attribution identifying
the creators/source, license link and applicable change indications. A paper
citation alone is not necessarily sufficient attribution, and the weight license
does not automatically apply to generated structures. Check other applicable
terms instead of asserting blanket ownership or restrictions on all outputs.

CCD/RDKit caches, cluster/release tables, HMM profiles, molecular archives,
native libraries and CUDA components are separate materials. Download capability,
checkpoint digest correctness and a model family's source license do not clear
their terms. No weights are added to DVBFixer distributions by this work.

## Citations

Keep legal notices separate from scientific citations. The
[initial bibliography](citations.bib) contains a small verified starting set,
not every tool/model/parameter reference. Additional model citation requirements
and links remain in the audit. The owner confirmed Stanislav Skripkin and the
initial project [CITATION.cff](../CITATION.cff) is now present; the complete third-party
bibliography and rights review remain open.
