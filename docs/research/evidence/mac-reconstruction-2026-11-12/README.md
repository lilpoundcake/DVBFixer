# Retained Mac reconstruction research evidence

Snapshot date:2026-11-12. These are research artifacts, **not accepted public
structure bundles**, model weights or runnable environments. Rejected proposals,
technical failures, original source inputs, masks, declarations, source snapshots,
native/restored/final structures and audit results remain distinguishable.

`manifest.json` records each original repository-relative path, size and SHA256,
plus ordered archive-part sizes/hashes and explicit exclusions. Parts are at most
32MiB so no individual artifact exceeds normal Git-hosting file limits. Each study
is one gzip-compressed tar stream, split into ordered `.partNNN` files.

To inspect a study, first verify the part hashes against the manifest. Concatenate
its parts in lexical order and extract **into a separate empty directory**, not
over a working checkout:

```bash
cat mac-reconstruction-comparison-2026-11-11.tar.gz.part* > /path/to/review/study.tar.gz
tar -tzf /path/to/review/study.tar.gz
tar -xzf /path/to/review/study.tar.gz -C /path/to/empty-review-directory
```

Then check extracted files against the manifest. Original paths begin with `tmp/`.
No symlinks are archived. Runtime/cache metadata and browser profiles are excluded;
owned engine environments/checkpoints under `.artifacts/`, local installations,
unrelated force-field directories and licensing-audit downloads are not included.
Host paths inside historical receipts are retained provenance, not portable runtime
locations. Reprovisioning still requires licensed dependencies and separately
obtained checkpoints under their own terms. Never execute archived source snapshots
merely because they are evidence; use reviewed project code.

The full archives intentionally include duplicate native/derivative records so
historical failures and source custody are not rewritten into a success-only view.
For readable conclusions, start with the
[comparison report](../../mac-reconstruction-comparison-2026-11-11.md) and
[real-gap follow-up](../../../plans/mac-reconstruction-real-gaps.md).
