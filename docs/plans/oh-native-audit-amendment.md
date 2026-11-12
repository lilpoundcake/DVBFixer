# Offline audit transport amendment (no new native execution)

The [frozen cohort](oh-native-qualification.md), committed at `5d99286`, executed
all twelve declared attempts once. Eleven returned success; case 11 stopped
after CHARMM-template fallback with `NamingConversionError: Duplicate source atom
identity`. Preserve that failed preparation and every saved H boundary; do not
repeat it or change its scientific configuration.

The first frozen auditor stopped at case 2 with `ValueError: All residues within
a chain must be contiguous`. Its reused small-graph reconstruction helper merged
all identical chain IDs. Real PDB topologies can have ordered A/B/C protein
blocks followed by separate A/B/C water blocks. This is an auditor transport
failure, not a production repair failure. Original source and failure status are
retained at `tmp/hydrogen-postflight-contract-01/native-audit-v1-failure/`.

## Prospective offline-only correction

Before re-auditing the retained cohort:

1. Keep the original frozen observer, auditor, native artifacts and plan bytes
   unchanged. A separate `scripts/audit_oh_native_v2.py` adapter may replace only
   the auditor's topology restoration function during an offline call.
2. Create a new chain whenever the ordered chain-ID run changes; preserve exact
   case-sensitive IDs, residue/atom encounter order and insertion codes. Repeated
   ID runs are distinct chain objects, not renamed chains.
3. Require each residue identity group to be contiguous and globally unique in
   the saved state. Reject ambiguous repeated residue identities instead of
   merging them. All twelve saved pre-repair states pass this admission check.
4. Restore every bond endpoint and exact saved nm float unchanged. Verify
   atom/residue identities, graph and count before matched replay. Do not fill,
   optimize, reflect, rename, re-seed or sample anything.
5. Retain the original frozen auditor's full scalar result, then add adapter and
   amendment digests plus transport inventory. Keep all original refusal and
   coverage gates. Even if matched repair succeeds, 11/12 preparation success
   cannot become a complete-cohort qualification pass.
6. Test separated repeated-ID runs, insertion codes and chain case, ambiguous
   repeated-residue refusal, and exact matched replay before offline execution.

## Limitation

The original snapshot did not serialize native chain-object occurrence indices.
Ordered label runs are enough to preserve unique residue ownership and the
explicit bond graph required by this local repair; they do not establish exact
reconstruction of adjacent native chain objects that share the same label.
Report this transport limitation. Exact candidate-coordinate replay is still
required; any mismatch remains a blocker. Future observer schemas may serialize
chain/residue occurrence indices explicitly, but this cohort is not resampled.
