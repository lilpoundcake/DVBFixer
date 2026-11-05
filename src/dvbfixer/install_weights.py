"""Opt-in official model downloads into a verified, repository-local cache."""

from __future__ import annotations

import hashlib
import os
import shutil
import tarfile
import tempfile
import time
import urllib.parse
import urllib.request
from pathlib import Path

from dvbfixer.install_progress import stage
from dvbfixer.installation_profiles import load_catalog
from dvbfixer.runtime_install import asset_requirements, checked_file

CC94_URL = "https://zenodo.org/api/records/16817230/files/protpardelle-1c.tar.gz/content"
CC94_SIZE = 4294232741
CC94_MD5 = "59fc043a9d981d6c74baa2d051ea164a"
PROTENIX_URL = "https://protenix.tos-cn-beijing.volces.com/checkpoint/protenix_base_default_v1.0.0.pt"
PROTENIX_SIZE = 1475950125


def description(name: str) -> tuple[str, int, str, str]:
    if name in {"protpardelle-cuda", "protpardelle-mps"}:
        return CC94_URL, CC94_SIZE, "CC BY 4.0", "https://zenodo.org/records/16817230"
    if name == "protenix-cuda":
        return PROTENIX_URL, PROTENIX_SIZE, "Apache 2.0 (v1 weights only)", "https://github.com/bytedance/Protenix/blob/85767b811c40ed46e73a9b39519cf6bfca8701ba/LICENSE"
    raise ValueError("No managed weights download is declared for this backend")


def cache_directory(repository: Path, name: str) -> Path:
    identity = "protpardelle-cc94" if name.startswith("protpardelle-") else "protenix-v1"
    description(name)  # Reject undeclared profiles, including historical cc89.
    path = repository / "models" / identity
    for item in (repository / "models", path):
        if item.is_symlink():
            raise ValueError(f"Model cache must not be a symlink: {item}")
    return path


def cached_checkpoint(repository: Path, name: str) -> Path | None:
    directory = cache_directory(repository, name)
    if not directory.exists():
        return None
    requirements = asset_requirements(load_catalog()["external"][name], repository)
    checkpoint = directory / "model_params/weights" / requirements["checkpoint"]["name"]
    for parent in (checkpoint.parent.parent, checkpoint.parent, directory / "model_params/configs"):
        if parent.is_symlink():
            raise ValueError(f"Model cache must not contain linked directories: {parent}")
    checked_file(checkpoint, requirements["checkpoint"]["sha256"], "Cached checkpoint")
    if "config" in requirements:
        checked_file(directory / "model_params/configs/cc94.yaml", requirements["config"], "Cached cc94 config")
    return checkpoint


class SecureRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        if urllib.parse.urlsplit(newurl).scheme != "https":
            raise ValueError("Refusing a non-HTTPS model download redirect")
        return super().redirect_request(req, fp, code, msg, headers, newurl)


def download(url: str, target: Path, size: int) -> None:
    """Bounded transfer into our private staging area, with visible byte progress."""
    opener = urllib.request.build_opener(SecureRedirect())
    count = 0
    notice = time.monotonic()
    with opener.open(url, timeout=60) as response, target.open("xb") as stream:
        if urllib.parse.urlsplit(response.geturl()).scheme != "https":
            raise ValueError("Model downloads require HTTPS")
        while chunk := response.read(min(1024 * 1024, size - count + 1)):
            count += len(chunk)
            if count > size:
                raise ValueError("Model download exceeds the declared size")
            stream.write(chunk)
            if time.monotonic() - notice >= 5:
                print(f"  Downloaded {count / 1e9:.2f}/{size / 1e9:.2f} GB ({100 * count / size:.1f}%)", flush=True)
                notice = time.monotonic()
        if count != size:
            raise ValueError(f"Incomplete model download: expected {size} bytes, received {count}; retry the download")
        stream.flush()
        os.fsync(stream.fileno())
    print(f"  Downloaded {count:,} bytes; verifying frozen assets.", flush=True)


def extract_cc94(archive: Path, directory: Path) -> None:
    """Copy only the two exact regular members; never extract arbitrary paths/links."""
    members = {
        "protpardelle-1c/model_params/weights/cc94_epoch3100.pth": 301344516,
        "protpardelle-1c/model_params/configs/cc94.yaml": 1024 * 1024,
    }
    seen = set()
    with tarfile.open(archive, "r|gz") as stream:
        for member in stream:
            if member.name not in members:
                continue
            if member.name in seen or not member.isfile() or member.size > members[member.name]:
                raise ValueError("Invalid or duplicate cc94 archive member")
            seen.add(member.name)
            target = directory / member.name.removeprefix("protpardelle-1c/")
            target.parent.mkdir(parents=True, exist_ok=True)
            source = stream.extractfile(member)
            if source is None:
                raise ValueError("Missing cc94 archive member bytes")
            with source, target.open("xb") as output:
                shutil.copyfileobj(source, output)
    if seen != set(members):
        raise ValueError("Official archive does not contain both required cc94 files")


def obtain(repository: Path, name: str, *, dry_run: bool = False) -> Path | None:
    """Caller must obtain consent first. Existing verified weights require no network."""
    existing = cached_checkpoint(repository, name)
    if existing is not None:
        print(f"Reusing SHA256-verified model weights: {existing}", flush=True)
        return existing
    url, size, _, _ = description(name)
    destination = cache_directory(repository, name)
    if dry_run:
        print(f"Preview: would download {size / 1e9:.2f} GB from {url} into {destination}; no download applied.", flush=True)
        return None
    destination.parent.mkdir(parents=True, exist_ok=True)
    if shutil.disk_usage(destination.parent).free < size + 400_000_000:
        raise ValueError("Insufficient free disk space for model download and extraction")
    requirements = asset_requirements(load_catalog()["external"][name], repository)
    with tempfile.TemporaryDirectory(prefix=".weights-", dir=destination.parent) as temporary:
        root = Path(temporary)
        payload = root / "payload"
        payload.mkdir()
        with stage("Download and verify model weights"):
            if "config" in requirements:
                archive = root / "archive.tar.gz"
                download(url, archive, size)
                with archive.open("rb") as stream:
                    if hashlib.file_digest(stream, "md5").hexdigest() != CC94_MD5:
                        raise ValueError("Official cc94 archive checksum mismatch")
                extract_cc94(archive, payload)
            else:
                weights = payload / "model_params/weights"
                weights.mkdir(parents=True)
                download(url, weights / requirements["checkpoint"]["name"], size)
            checkpoint = payload / "model_params/weights" / requirements["checkpoint"]["name"]
            checked_file(checkpoint, requirements["checkpoint"]["sha256"], "Downloaded checkpoint")
            if "config" in requirements:
                checked_file(payload / "model_params/configs/cc94.yaml", requirements["config"], "Downloaded cc94 config")
            # Exclusive reservation plus no-replace files: never overwrite a cache
            # created by another process, nor silently repair a damaged cache.
            try:
                destination.mkdir()
            except FileExistsError:
                return cached_checkpoint(repository, name)
            for source in payload.rglob("*"):
                target = destination / source.relative_to(payload)
                if source.is_dir():
                    target.mkdir(exist_ok=True)
                else:
                    os.link(source, target)
    return cached_checkpoint(repository, name)


def select(args, repository: Path, *, interactive: bool) -> None:
    """Select shared cache, consented download, local file or deferred setup."""
    from dvbfixer.install_ui import ask, choose

    if args.checkpoint is not None:
        if args.download_weights:
            raise ValueError("Choose --checkpoint or --download-weights, not both")
        return
    if args.runtime == "psi4":
        if args.download_weights:
            raise ValueError("PSI4 does not use model weights")
        return
    if not interactive and not args.download_weights:
        return
    existing = cached_checkpoint(repository, args.runtime)
    if existing is not None:
        print(f"Reusing SHA256-verified model weights: {existing}", flush=True)
        args.checkpoint = existing
        return
    url, size, license_name, terms = description(args.runtime)
    print(f"\nManaged weights: {size / 1e9:.2f} GB transfer from {url}\n"
          f"Terms: {license_name} — {terms}\n"
          f"Shared storage (excluded from Git): {cache_directory(repository, args.runtime)}\n"
          "Only the required pinned weights/config are retained. Protenix data/Kalign remain separate inputs.", flush=True)
    if args.download_weights:
        if not args.accept_model_terms:
            if not interactive:
                raise ValueError("Review the model terms above, then pass --accept-model-terms with --download-weights")
            from dvbfixer.install_ui import confirm

            if not confirm("Accept the displayed model terms and download weights?"):
                raise EOFError()
        selection = "download"
    else:
        selection = choose("Model weights", [
            ("download", "Accept displayed model terms and download into shared models/"),
            ("local", "Use an existing local weights file"),
            ("later", "Configure later (backend will not be ready)"),
        ], default="later")
    if selection == "download":
        args.checkpoint = obtain(repository, args.runtime, dry_run=args.dry_run)
    elif selection == "local":
        path = ask("Checkpoint weights path (empty to configure later): ").strip()
        args.checkpoint = Path(path).expanduser() if path else None
