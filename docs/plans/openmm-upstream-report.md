# Upstream report preparation and required human opening

## Scope and mandatory publication gate

After `32be60e`, the user instructed continuation in response to the explicit
upstream-issue proposal. Prepare a report on `research/openmm-upstream-report`,
but do NOT publish automatically: upstream's AI_POLICY.md requires every issue
to BEGIN with a human-written summary in the submitter's own words. AI-generated
technical material must be clearly labeled and verified. A generic continued-
work instruction is not a human-authored technical summary and must never be
passed off as one. Ask the user for that opening and review before publication.
Upstream also asks that replies to maintainer questions come from the human,
not an AI; do not auto-reply later. No patch, new environment, production/main
change, external comment, issue or PR is authorized past that unmet policy gate.

Documentation plan: this plan and frozen single packaging-check JSON; a clearly
AI-generated technical appendix; complete readiness/evidence JSON; then canonical
handoff. Pending links none; create referenced files before linking and validate
all links/hashes/public payloads. Preserve prior scientific records unchanged.

## Duplicate and contribution checks

Read OpenMM bug-report template, CONTRIBUTING.md and AI_POLICY.md. Search open
and closed issues/PRs for evaluateCpu, constraint+gradient, penalty+gradient,
constraint+sign and overflow+minimizer. Related #5394 concerns unconstrained CUDA
fallback performance; do not hijack it. #5434/#5435 cover known Apple M3 out-of-
range float-to-long overflow, with a merged saturation fix. Keep the intended
constraint-gradient sign report distinct: the successful trigger is INSIDE
signed fixed-point range, above the detector threshold. Report only tested
revision/platform; do not imply current master or all releases are affected.

Bug template requests testInstallation diagnostics; record if not executed,
not a fabricated pass. Existing available-platform/plugin metadata and completed
CPU/Reference controls are not a substitute for that exact command's output.
Native library directory discovery remains permission-denied/not bypassed;
runtime strings are not installed-binary provenance.

## Single prospective portable-entry packaging check

The compact script's native entry was not executed during the prior factorial
study (its System and raw-callback offline replay passed). Upstream specifically
asks that an AI-generated reproducer be checked to run and reproduce the bug.
Before presenting it as runnable evidence, execute its EXACT unchanged portable
entry ONCE, as a separately frozen packaging verification, NOT another scientific
case or qualification seed. This is not retry-until-green: the factorial primary
already succeeded; no constants, source, tolerances, model or point selection
may change. Retain failure and stop if this one invocation fails.

Bind script SHA256, original controlled evidence/validation, exact interpreter
metadata and command before execution. One child process, one native minimizer
call,30-second outer limit, kill process group on timeout. Project-owned exclusive
output/temp/cache directories; PYTHONDONTWRITEBYTECODE=1; no package installation,
structural input, extra force query or allowed external library read. Capture
stdout/stderr and actual exit status. Require returned=true, at least one callback
and a witness rejecting physical gradient even with original allowance while
source-sign-float32 is exact. Missing/negative results block issue readiness;
no retry, parameter edits or hidden extra calls. Retain source hash before/after.

Prior scientific accounting stays24 calls/174 callbacks; this adds ONE explicit
portable-entry packaging call only. It cannot change first-series failures,
structural qualification or installed-binary/internal-event attestation.

## Draft and final handoff

Prepare a technical appendix with accurate measured scalar witness, tested
environment, source hypothesis, control matrix, pinned code/data links and
limitations. No credentials, local usernames/paths, biological inputs or measured
coordinate/RNG arrays. Fixed synthetic starting constants in source are not
private observations. State whether compact entry verification passed; do not
claim a human review or opening exists until provided. Publish NOTHING until
human summary/review is supplied; never fill it with an AI impersonation.
