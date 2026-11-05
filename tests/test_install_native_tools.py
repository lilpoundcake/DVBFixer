"""Actual filesystem controls for the pip/native obminimize collision."""

import csv
import json
import os

import pytest

from dvbfixer.install_native_tools import preserve_native_openbabel_tool


def native_tool(tmp_path):
    prefix = tmp_path / "env"
    (prefix / "conda-meta").mkdir(parents=True)
    (prefix / "conda-meta/openbabel.json").write_text(json.dumps({
        "name": "openbabel", "files": ["bin/obminimize"],
    }))
    executable = prefix / "bin/obminimize"
    executable.parent.mkdir()
    executable.write_bytes(b"native binary")
    executable.chmod(0o755)
    record = prefix / "lib/python3.12/site-packages/dvbfixer-1.dist-info/RECORD"
    record.parent.mkdir(parents=True)
    record.write_text("../../../bin/obminimize,hash,12\ndvbfixer/__init__.py,hash,20\n")
    return prefix, executable, record


@pytest.mark.parametrize("fail", [False, True])
def test_hardlinked_cache_survives_pip_write_and_failure(tmp_path, fail):
    prefix, executable, record = native_tool(tmp_path)
    cache = tmp_path / "cached-binary"
    os.link(executable, cache)

    def apply():
        with preserve_native_openbabel_tool(prefix):
            assert not executable.exists()
            executable.write_text("from dvbfixer.openbabel_cli import obminimize\n")
            if fail:
                raise RuntimeError("pip failed")

    if fail:
        with pytest.raises(RuntimeError, match="pip failed"):
            apply()
    else:
        apply()
    assert executable.read_bytes() == cache.read_bytes() == b"native binary"
    assert executable.stat().st_mode & 0o777 == 0o755
    with record.open(newline="") as stream:
        assert list(csv.reader(stream)) == [["dvbfixer/__init__.py", "hash", "20"]]


def test_relative_native_symlink_is_restored_without_touching_target(tmp_path):
    prefix, executable, _ = native_tool(tmp_path)
    target = executable.with_name("obminimize-native")
    executable.rename(target)
    executable.symlink_to(target.name)
    with preserve_native_openbabel_tool(prefix):
        executable.write_text("from dvbfixer.openbabel_cli import obminimize\n")
    assert executable.is_symlink()
    assert os.readlink(executable) == target.name
    assert target.read_bytes() == b"native binary"


def test_unexpected_concurrent_file_is_not_overwritten(tmp_path, capsys):
    prefix, executable, record = native_tool(tmp_path)
    with pytest.raises(ValueError, match="changed unexpectedly"):
        with preserve_native_openbabel_tool(prefix):
            executable.write_bytes(b"user replacement")
    assert executable.read_bytes() == b"user replacement"
    assert "bin/obminimize" in record.read_text()
    backup = capsys.readouterr().out.split("backup retained for recovery: ")[1].strip()
    from pathlib import Path

    assert Path(backup).read_bytes() == b"native binary"
    Path(backup).unlink()
    Path(backup).parent.rmdir()
