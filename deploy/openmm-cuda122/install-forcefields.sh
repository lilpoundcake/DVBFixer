#!/usr/bin/env bash
set -euo pipefail

prefix=${1:-${CONDA_PREFIX:-}}
if [[ -z "$prefix" || ! -x "$prefix/bin/python" ]]; then
    echo "usage: $0 /path/to/conda-environment" >&2
    exit 2
fi

data_dir=$(
    "$prefix/bin/python" -c \
        'from pathlib import Path; import openmm.app; print(Path(openmm.app.__file__).parent / "data")'
)
target_dir="$data_dir/amber19"
mkdir -p "$target_dir"

tmp_dir=$(mktemp -d)
trap 'rm -rf "$tmp_dir"' EXIT
base_url=https://raw.githubusercontent.com/openmm/openmm/8.3.0/wrappers/python/openmm/app/data/amber19

curl --fail --location --silent --show-error \
    "$base_url/protein.ff19SB.xml" --output "$tmp_dir/protein.ff19SB.xml"
curl --fail --location --silent --show-error \
    "$base_url/tip3p.xml" --output "$tmp_dir/tip3p.xml"

(
    cd "$tmp_dir"
    printf '%s  %s\n' \
        086bfdb1c05e7d5cb1d330e6c39fab5056a95e89492cdbd3044c746db13e9d09 \
        protein.ff19SB.xml \
        3f4b188dbcb6c02863230eaca231e927fb6bf3307ce947d8a50d0f46f6dd83d9 \
        tip3p.xml | sha256sum --check --strict
)

install -m 0644 "$tmp_dir/protein.ff19SB.xml" "$target_dir/protein.ff19SB.xml"
install -m 0644 "$tmp_dir/tip3p.xml" "$target_dir/tip3p.xml"
