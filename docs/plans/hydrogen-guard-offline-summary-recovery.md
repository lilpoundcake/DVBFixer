# Offline-only summary recovery; no native retry

The two native paths frozen at `a430c86` produced all six saved returned phases.
The parent runner then exited with `FileExistsError`: the exclusive `write`
helper was incorrectly called on the same `evidence.json` inside the loop.
The first path's summary is intact; the second path's worker exit code was not
durably recorded. This is an orchestration failure, not a passed whole run.

Before offline audit, create a **new** private audit-view directory. Bind the
original frozen plan, partial summary, parent error log and every retained path
artifact. Preserve the original directory and all its bytes. Reuse its path
directories through read-only-by-convention symlinks; copy the exact frozen plan
and create one new complete accounting index in the view. Never modify the
original partial index or frozen runner, and execute no new Context/minimization
in recovery. The subsequent independent physical audit may create its normal
readback Contexts; these are not new trajectories.

Use the original recorded status/exit code for any path already indexed. For an
unindexed path, verify all three phase result/path/returned-state files: true
native-return flags, callback counts and exact final-coordinate equality with
the serialized returned State. Report `phases-returned-exit-unrecorded`, with
null worker exit code; never invent `completed` or exit zero. Missing/failed
phases remain partial. The audit-view index binds hashes of all original files.
The frozen independent physics/identity auditor then reads that new index and
the untouched original payloads. Preserve parent failure, missing exit-code
status and every gradient/convergence limitation alongside any favorable CIP.

Documentation scope: this explicit amendment, recovery/audit scalar records,
then update the existing result/handoff. No native enrollment, force field,
proposal, transport tolerance, stopping rule or scientific assertion changes.
Pending links: none. Public coordinate/RNG/state payloads remain forbidden.
