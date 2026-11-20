"""Copy pinned research dependencies with explicit CPU/Python import fixes only."""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
from pathlib import Path

IDPFORGE = "e8da98924e956fd307a6df903432e5a271637076"
OPENFOLD = "be2ec1841f16c966c65ae0e7599ebbadc725757d"


def replacements(name):
    if name == "esm/esmfold/trunk.py":
        return [("from dataclasses import dataclass", "from dataclasses import dataclass, field"),
                ("structure_module: StructureModuleConfig = StructureModuleConfig()",
                 "structure_module: StructureModuleConfig = field(default_factory=StructureModuleConfig)")]
    if name in {"openfold/model/structure_module.py", "openfold/utils/kernel/attention_core.py"}:
        return [('attn_core_inplace_cuda = importlib.import_module("attn_core_inplace_cuda")',
                 'attn_core_inplace_cuda = (importlib.import_module("attn_core_inplace_cuda")\n'
                 '                          if torch.cuda.is_available() else None)')]
    return []


def patch_text(name, text):
    for before, after in replacements(name):
        if text.count(before) != 1:
            raise ValueError(f"unexpected pinned source for {name}")
        text = text.replace(before, after)
    return text


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--upstream", type=Path, required=True)
    parser.add_argument("--openfold", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    for source, revision in ((args.upstream, IDPFORGE), (args.openfold, OPENFOLD)):
        actual = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=source, text=True).strip()
        dirty = subprocess.check_output(["git", "status", "--porcelain", "--untracked-files=no"], cwd=source, text=True).strip()
        if actual != revision or dirty:
            raise ValueError("require clean pinned research sources")
    args.output.mkdir(parents=True, exist_ok=False)
    for name, parent in (("idpforge", args.upstream), ("esm", args.upstream), ("openfold", args.openfold)):
        shutil.copytree(parent / name, args.output / name, ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
    shutil.copyfile(args.upstream / "LICENSE", args.output / "IDPFORGE-LICENSE")
    shutil.copyfile(args.openfold / "LICENSE", args.output / "OPENFOLD-LICENSE")
    rows = []
    for name in ("esm/esmfold/trunk.py", "openfold/model/structure_module.py", "openfold/utils/kernel/attention_core.py"):
        path = args.output / name
        before = path.read_bytes()
        after = patch_text(name, before.decode()).encode()
        path.write_bytes(after)
        rows.append({"path": name, "before_sha256": hashlib.sha256(before).hexdigest(),
                     "after_sha256": hashlib.sha256(after).hexdigest(),
                     "replacements": replacements(name)})
    (args.output / "compatibility.json").write_text(json.dumps({
        "idpforge_revision": IDPFORGE, "openfold_revision": OPENFOLD, "patches": rows,
        "inference_math_changed": False, "cuda_kernels_available": False,
        "note": "CPU default attention only; requesting the unavailable CUDA kernel fails, never substitutes"}, indent=2) + "\n")


if __name__ == "__main__":
    main()
