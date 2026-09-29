# Diffusion Confirmatory Cohort Results

- Status: complete.
- Evaluation date: 2026-09-29.
- Cohort: 231 selected sequence-cluster representatives from 500 screened groups.
- Eligible candidate: Protenix v1.
- Comparator: Salilab MODELLER 10.8.
- Proxy-only backend: Boltz-2.
- Aggregate schema: 1.
- Aggregate artifact: `.artifacts/diffusion-confirmatory-cohort-500-v5/aggregate.json`.
- Aggregate SHA-256: `55e4cbc6c16ee6bf34a14aabf61ec71afcddfb436f15edb5195c3d8b3f717ad7`.
- Decision: `candidate-preferred` for Protenix v1 within the tested scope.

## Scope And Provenance

- [x] Screen all 500 metadata-locked independence groups before inference.
- [x] Select 231 accepted groups, above the preregistered minimum of 180 and below the maximum of 300.
- [x] Use one withheld internal gap per sequence-cluster independence group.
- [x] Evaluate 128 five-residue gaps and 103 ten-residue gaps.
- [x] Preserve case-sensitive chain IDs and insertion codes.
- [x] Require complete canonical generated heavy atoms.
- [x] Preserve fixed heavy coordinates exactly after final serialization.
- [x] Count hard-gate failures against the paired pass-rate endpoint.
- [x] Keep Boltz-2 outside confirmatory selection because its training leakage remains unresolved.

The v5 materialization corrected the generated-atom contract for 13 cases whose
deposited withheld residues had incomplete side chains. The normalized input,
reference, target FASTA, and fixed-atom sets are byte- or identity-equivalent
between v4 and v5 for all 231 selected cases. The generated-atom sets are
equivalent for 218 cases; their completed v4 inference results are therefore
reused. The remaining 13 cases add 103 canonical heavy atoms and were rerun
under the v5 contract. The final aggregate contains 218 reused results and 13
v5 reruns.

## Primary Results

| Backend | Role | Passed | Pass rate | Wilson 95% CI | Median valid gap-backbone RMSD |
|---|---|---:|---:|---:|---:|
| Protenix v1 | Eligible candidate | 220/231 | 95.24% | 91.68-97.32% | 0.408 A |
| MODELLER 10.8 | Comparator | 169/231 | 73.16% | 67.10-78.46% | 5.330 A |
| Boltz-2 | Proxy only | 220/231 | 95.24% | 91.68-97.32% | 0.423 A |

Correcting the 13 incomplete generated-atom masks increased both all-atom
backends from 208/231 to 220/231 passing cases. Twelve corrected cases pass all
hard gates; `9ejk` remains a hard-gate failure for both all-atom backends.

## Protenix Versus MODELLER

- [x] Paired case count: 231.
- [x] Both passed: 160.
- [x] Protenix only passed: 60.
- [x] MODELLER only passed: 9.
- [x] Neither passed: 2.
- [x] Pass-rate difference: +22.08 percentage points for Protenix.
- [x] Pass-rate-difference 95% CI: +15.58 to +28.57 percentage points.
- [x] Paired noninferiority margin: 5 percentage points.
- [x] Noninferiority demonstrated: yes.
- [x] Median paired RMSD difference: -4.667 A for Protenix.
- [x] Cluster-bootstrap RMSD-difference 95% CI: -5.036 to -4.181 A.
- [x] Both-valid independence groups: 160, above the preregistered minimum of 30.

Protenix satisfies the preregistered selection rule and is preferred to the
MODELLER comparator for the tested one-chain, internal-gap scope.

## Protenix Versus Boltz-2

- [x] Both passed: 216.
- [x] Protenix only passed: 4.
- [x] Boltz-2 only passed: 4.
- [x] Neither passed: 7.
- [x] Median paired RMSD difference among common passing cases: -0.005 A for Boltz-2.
- [x] Boltz-2 confirmatory eligibility: false.
- [x] Exclusion reason: unresolved training leakage.

The two all-atom backends are descriptively similar on this cohort. Boltz-2 is
not eligible for selection regardless of observed quality until an immutable
training manifest or an author-provided held-out split resolves membership and
homology leakage.

## Gap-Length Strata

| Gap length | Backend | Passed | Pass rate | Median valid RMSD |
|---:|---|---:|---:|---:|
| 5, n=128 | Protenix v1 | 125 | 97.66% | 0.326 A |
| 5, n=128 | Boltz-2 | 125 | 97.66% | 0.328 A |
| 10, n=103 | Protenix v1 | 95 | 92.23% | 0.623 A |
| 10, n=103 | Boltz-2 | 95 | 92.23% | 0.585 A |

The expected difficulty increase is visible for ten-residue gaps, but both
all-atom backends remain above 92% hard-gate pass rate in that stratum.

## Remaining Hard-Gate Failures

- Shared Protenix/Boltz-2 failures:
  - `22uy`
  - `24rc`
  - `8xyp`
  - `9ddr`
  - `9e16`
  - `9ejk`
  - `9jjb`
- Protenix-only failures:
  - `13ib`
  - `9g9t`
  - `9jco`
  - `9jqa`
- Boltz-2-only failures:
  - `9ed2`
  - `9hab`
  - `9i9l`
  - `9j4f`
- Observed failure classes:
  - junction peptide connectivity;
  - generated-region backbone breaks;
  - bond-length and bond-angle outliers;
  - amide-planarity outliers;
  - Ramachandran outliers;
  - severe steric overlap.

## Quality Distribution

| Backend | Q1 | Median | Q3 | P90 | Maximum among passing cases |
|---|---:|---:|---:|---:|---:|
| Protenix v1 | 0.286 A | 0.408 A | 0.838 A | 1.757 A | 10.768 A |
| Boltz-2 | 0.268 A | 0.423 A | 0.751 A | 1.444 A | 6.194 A |
| MODELLER 10.8 | 3.796 A | 5.330 A | 7.204 A | 9.059 A | 29.209 A |

Hard-gate success is not an absolute RMSD threshold. The high-RMSD passing tail
represents geometrically valid alternative conformations and must remain visible
rather than being silently removed from the endpoint.

## Runtime And Memory

| Backend | Median complete wall time | Aggregate sequential wall time | Peak RAM | Peak VRAM |
|---|---:|---:|---:|---:|
| Protenix v1 | 209.0 s | 22.30 h | 3.68 GiB | 6.61 GiB |
| Boltz-2 | 160.0 s | 19.63 h | 5.38 GiB | 3.93 GiB |
| MODELLER 10.8 | 16.2 s | 1.50 h | 0.15 GiB | Not applicable |

Runtime and memory are operational evidence only and do not override scientific
eligibility or the preregistered selection rule.

## Decision And Limits

- [x] Select Protenix v1 for implementation of an opt-in experimental
  `model --backend diffusion` path.
- [x] Keep MODELLER as the default and supported production baseline.
- [x] Prohibit automatic fallback between diffusion and MODELLER.
- [x] Keep Boltz-2 as proxy-only evidence.
- [x] Preserve every frozen hard gate in the production path.
- [ ] Complete the production runner environment and checkpoint distribution review.
- [ ] Add public CLI dispatch and production request construction.
- [ ] Add production PDB/`.dat` publication without partial outputs.
- [ ] Add representative multichain/interface evidence before claiming that scope.
- [ ] Add retained ligand, glycan, PTM, cofactor, metal, or covalent-link evidence before claiming those scopes.
- [ ] Make production support a separate release decision after the experimental integration passes deployment acceptance.

The result supports implementation of the narrow experimental backend. It does
not authorize changing the default backend, shipping model weights, claiming
unsupported chemistry, or promoting diffusion-backed homology modeling.

## Reproduction

```bash
python scripts/materialize_diffusion_confirmatory_cohort.py \
  --output-root .artifacts/diffusion-confirmatory-cohort-500-v5

python scripts/run_diffusion_confirmatory_cohort.py aggregate \
  .artifacts/diffusion-confirmatory-cohort-500-v5 \
  --output .artifacts/diffusion-confirmatory-cohort-500-v5/aggregate.json \
  --bootstrap-samples 10000
```

The ignored aggregate is guarded by
`.artifacts/diffusion-confirmatory-cohort-500-v5/aggregate.json.complete.json`,
which records every aggregation input and the output digest.
