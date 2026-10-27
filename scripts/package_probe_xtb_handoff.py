"""Package the audited source snapshot and selected evidence, never local binaries.

Run only after all audit/test processes finish. Existing output is refused.
No commits, dependency changes, or recursive collection of unrelated work.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import subprocess
import tarfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / ".artifacts/probe-xtb"
OVERLAY = [
    "scripts/audit_probe_contacts.py", "scripts/audit_probe_xtb_contracts.py",
    "scripts/audit_xtb_controls.py", "scripts/audit_xtb_geometry_deltas.py",
    "scripts/run_probe_xtb_audit.py", "scripts/summarize_probe_xtb_audit.py",
    "scripts/package_probe_xtb_handoff.py", "tests/test_probe_xtb_audit.py",
    "tests/test_probe_xtb_research_metrics.py", "docs/research/probe-xtb-macos-audit.md",
    "docs/plans/probe-xtb-linux-handoff.md",
]
EVIDENCE = [
    "complete-macos", "extended-macos", "geometry-final", "package-provenance",
    "separate-fixed-0", "separate-fixed-1", "separate-installed-0", "separate-installed-1",
    "native-contracts-final", "cmake-netlib-neutral-control", "cmake-netlib-calibration",
    "upstream-fixed-controls", "upstream-fixed-glycan", "upstream-fixed-plm",
    "glycan-lanc-v1", "meson-netlib-glycan", "meson-netlib-controls",
    "probe-aligned-v1", "probe-aligned-v2", "results-final.json", "results-extended.json",
    "contracts-without-xfail.xml", "fast-final.log", "fast-final.xml",
    "fast-final-path.log", "fast-final-path.xml", "focused-final.log", "quality-checks.log",
    "1vcu-prepare.log", "1vcu-openmm.log", "docs-plan.md", "docs-links-todo.md", "docs-quality-report.json",
    "environment-provenance.json",
    "audit-plan-source.md",
    "upstream-metadata",
]


def git(repo, *arguments):
    return subprocess.run(["git", "-C", str(repo), *arguments], check=True, capture_output=True).stdout


def source_tar(repo, dependencies=()):
    output = io.BytesIO()
    with tarfile.open(fileobj=output, mode="w") as archive:
        for location, prefix in [(repo, ""), *dependencies]:
            data = git(location, "archive", "--format=tar", "HEAD")
            with tarfile.open(fileobj=io.BytesIO(data)) as source:
                for member in source:
                    if not member.isfile():
                        if member.issym() or member.islnk():
                            raise ValueError(f"Review source link before packaging: {member.name}")
                        continue
                    member.name = prefix + member.name
                    archive.addfile(member, source.extractfile(member))
    return output.getvalue()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or Path(str(args.output) + ".sha256").exists():
        parser.error("Output and digest must not exist")
    missing = [str(AUDIT / name) for name in EVIDENCE if not (AUDIT / name).exists()]
    if missing:
        parser.error(f"Evidence incomplete: {missing}")
    for name in ("complete-macos", "extended-macos"):
        if not json.loads((AUDIT / name / "report.json").read_text())["execution_complete"]:
            parser.error(f"Replay incomplete: {name}")
    upstream = AUDIT / "xtb-upstream-fixed"
    dependency_names = ["mctc-lib", "jonquil", "toml-f", "test-drive"]
    dependencies = [(upstream / "subprojects" / name, f"subprojects/{name}/") for name in dependency_names]
    repositories = {"dvbfixer": ROOT, "xtb": upstream, "probe": AUDIT / "probe-source",
                    **{name: path for name, (path, _) in zip(dependency_names, dependencies)}}
    provenance = {name: {"commit": git(path, "rev-parse", "HEAD").decode().strip(),
                         "tracked_diff": git(path, "diff", "HEAD").decode()}
                  for name, path in repositories.items()}
    provenance["overlay"] = OVERLAY
    provenance["note"] = "Source snapshot, not a git checkout; native fixtures are historical evidence only"
    digests = {}
    with args.output.open("xb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as compressed, \
            tarfile.open(fileobj=compressed, mode="w|") as archive:

        def add(name, data, mode=0o644):
            if name in digests:
                raise ValueError(f"Duplicate archive entry: {name}")
            if data[:4] in (b"\xcf\xfa\xed\xfe", b"\xfe\xed\xfa\xcf", b"\xca\xfe\xba\xbe", b"\x7fELF"):
                raise ValueError(f"Native executable/object excluded: {name}")
            info = tarfile.TarInfo(name)
            info.size, info.mode = len(data), mode
            archive.addfile(info, io.BytesIO(data))
            digests[name] = hashlib.sha256(data).hexdigest()

        with tarfile.open(fileobj=io.BytesIO(git(ROOT, "archive", "--format=tar", "HEAD"))) as source:
            for member in source:
                if member.isfile() and member.name not in OVERLAY:
                    add("repository/" + member.name, source.extractfile(member).read(), member.mode)
                elif member.issym() or member.islnk():
                    raise ValueError(f"Review repository link: {member.name}")
        for name in OVERLAY:
            add("repository/" + name, (ROOT / name).read_bytes())
        for name in EVIDENCE:
            path = AUDIT / name
            files = sorted(p for p in path.rglob("*") if p.is_file()) if path.is_dir() else [path]
            for file in files:
                if "__pycache__" not in file.parts and file.suffix != ".pyc":
                    add("evidence/" + str(file.relative_to(AUDIT)), file.read_bytes())
        for name in ("testlog.txt", "meson-log.txt"):
            path = AUDIT / "xtb-upstream-meson/meson-logs" / name
            add("evidence/upstream-meson-logs/" + name, path.read_bytes())
        for name, path in [("buildoptions.json", "intro-buildoptions.json"), ("dependencies.json", "intro-dependencies.json")]:
            add("evidence/upstream-meson-config/" + name,
                (AUDIT / "xtb-upstream-meson/meson-info" / path).read_bytes())
        add("sources/xtb.tar", source_tar(upstream, dependencies))
        add("sources/probe.tar", source_tar(AUDIT / "probe-source"))
        add("provenance.json", (json.dumps(provenance, indent=2) + "\n").encode())
        add("README.md", b"Start with repository/docs/plans/probe-xtb-linux-handoff.md.\nVerify MANIFEST.sha256 before use. Linux acceptance is not yet tested.\n")
        manifest = "".join(f"{digest}  {name}\n" for name, digest in sorted(digests.items()))
        add("MANIFEST.sha256", manifest.encode())
    with tarfile.open(args.output, "r:gz") as archive:
        count = 0
        for member in archive:
            if member.isfile():
                assert hashlib.sha256(archive.extractfile(member).read()).hexdigest() == digests[member.name]
                count += 1
    digest = hashlib.sha256(args.output.read_bytes()).hexdigest()
    Path(str(args.output) + ".sha256").write_text(f"{digest}  {args.output.name}\n")
    print(json.dumps({"archive": str(args.output), "sha256": digest, "verified_files": count,
                      "bytes": args.output.stat().st_size}, indent=2))


if __name__ == "__main__":
    main()
