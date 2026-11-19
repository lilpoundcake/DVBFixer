# WP2E: explicit observed-source preparation without atom edits

Date: 2026-11-19. Prospective protocol, baseline `c0e4d8b`; continue on persistent
`feature/chemical-reconstruction`. Freeze before executing the new source controls.

## Scope and documentation plan

Reuse the whole-PDB parser, exact identities, original observed-peptide adapter,
canonical inventory transport, peptide witnesses, constructor and independent
auditor. This slice is private heavy-atom research, not public source repair or
general H placement. Documentation: this protocol, source validation JSON, active
handoff and agent task map. Pending links: none. Preserve prior protocols/results.

Create a separate immutable, source-bound preparation record containing original
and prepared ASCII/LF content, both digests and exact removed metadata records
(one-based line index, original line, reason). Its validation must recompute the
projection. Atom records, coordinates, occupancy, names, order, chain case,
insertion codes, other headers and all valid linkage records remain byte-identical.
The scientific request receives only retained observations, not withheld
experimental coordinates or the complete original/prepared source payloads.

## Declared metadata and occupancy policy

1. CONECT: use the existing `renumber.update_conect` owner with an identity serial
   map to decide dangling-record rejection, never a partial rewrite. Drop a whole
   record if any endpoint serial is absent, or if it has fewer than two serials.
   Record every removal. Preserve valid CONECT bytes exactly. Strictly refuse
   malformed/non-numeric/truncated nonblank fields and self-links rather than
   deleting them as nuisance metadata. Do not infer missing bonds or reinterpret
   valid unknown LINK chemistry. No input file is overwritten.
2. Altloc/multiple models: no selection is authorized here. Keep the original
   source and refuse explicitly before atom projection; do not strip 1A3N's
   hemes, choose occupancy winners or silently delete alternate atoms.
3. Fractional occupancy: a separate versioned request/policy may explicitly retain
   the exact fractional-occupancy atoms in unselected **opaque** context. Preserve
   their occupancy and coordinates unchanged and include them fully/conservatively
   in clash measurements. They remain chemically unassessed, not repaired or
   probabilistically discounted. Exact scope, prepared-source SHA and preparation
   digest must be bound into the new request. Original request schema/default
   digests and original occupancy refusal remain unchanged.
4. The allowed fractional set must exactly equal non-unit observations, excluding
   every target, resolved component, owned/mobile atom, constraint atom and
   connection endpoint. Validate target occupancy BEFORE withholding so a mask
   cannot hide uncertain selected observations. Zero/invalid occupancy remains
   unsupported. Fractional coordinates are never geometry anchors. No altloc
   exemption and no weakening of the existing observed H guard.

Preparation is not scientific admission. Bind its evidence into a separate
prepared-observed request, then use the existing admission/compiler/proposal path.
Any remaining source/chemistry/geometry refusal survives. Do not bypass admission,
crop context, widen bounds, rename observed atoms or try an alternative engine.
RTP inventory canonicalization is the explicit WP2C policy, not a new fallback.
The original observed-peptide adapter and WP2B runner retain their old behavior.

## Three prospective real-source controls

Use the original tracked sources/pins and the SAME masks selected in WP2B:

| Source | Target | Withhold |
|---|---|---|
| `protein_context/1UBQ.pdb` | A/PHE45 | CB CG CD1 CD2 CE1 CE2 CZ |
| `regressions/1FR2.pdb` | first interior chain-A THR (previously A/THR11) | CB OG1 CG2 |
| `protein_context/1A3N.pdb` | first interior chain-A ALA (previously A/ALA5) | CB |

Three declared entries; preparation, request admission, coordinate construction
and geometry are separate denominators. Hypotheses: 1UBQ can retain fractional
opaque context without falsifying occupancy; 1FR2 can remove only dangling/empty
CONECT metadata; 1A3N still refuses unresolved altloc. Do NOT predict geometry
success or change any bound after seeing outcomes. Retain all parse, planning,
constructor and raw/rounded failures. Use the unchanged WP2B observed-heavy
policy: 0.10 A bond, 12-degree angle, 0.7 A minimum separation, 0.02 A anchor bound;
stereo minimum 0.05 and planar tolerance 0.05 A. No held-out RMSD ranking/scoring.

Check protocol/source pins before exclusive output reservation. Retain original
source bytes for every entry and prepared bytes/report when available, all inside
the project evidence directory. Re-read original files after execution and verify
unchanged hashes. Evidence must prove that only declared metadata lines changed,
all other records and retained observations survived, and withholding remains an
explicit separate operation. No PDB candidate/.dat/public bundle is materialized.

## Regression and acceptance gates

Test preparation recomputation/binding, mixed valid/dangling records with complete
record dropping, malformed/self links, exact occupancy retention, conservative
fractional-context clash rejection, attempts to hide fractional target atoms in a
mask, fractional constraint/link endpoints, wrong source/evidence digest, altloc/
multi-model refusal, no original-file write and old request/digest/refusal behavior.
Source preparation must not become a generic escape hatch from admission.

Run focused source/heavy/H/identity tests, Ruff and scoped strict mypy, then full
fast suite. Report inherited transitive typing errors separately and verify any
new dependency reachability against unchanged archived source. Replay old pinned
controls byte-for-byte as regressions, not new scientific evidence. Preserve
original fixture/FF bytes and final source/report hashes; update handoff/task map.

Actual neighbour chemistry, full-polymer state/CIP, general real-complex heavy/H
qualification, H operation on observed sources and public/default/dependency or
production/data-license permission remain OPEN. No prospective pass count closes
those gates.
