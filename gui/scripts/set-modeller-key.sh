#!/usr/bin/env bash
# Sets the Modeller academic license key inside the active dvbfixer
# micromamba env. Modeller reads this from <env>/lib/modeller-*/modlib/modeller/config.py.
#
# Usage:
#   micromamba activate dvbfixer
#   MODELLER_LICENSE_KEY=<your-license-key> bash scripts/set-modeller-key.sh
#
# Or with an explicit env prefix:
#   MODELLER_LICENSE_KEY=<your-license-key> bash scripts/set-modeller-key.sh /opt/conda/envs/dvbfixer
#
# Get a free academic license at https://salilab.org/modeller/registration.html

set -euo pipefail

if [ -z "${MODELLER_LICENSE_KEY:-}" ]; then
  echo "error: MODELLER_LICENSE_KEY env var is required."
  echo "       Register at https://salilab.org/modeller/registration.html"
  echo "       Then: MODELLER_LICENSE_KEY=<your-license-key> bash scripts/set-modeller-key.sh"
  exit 1
fi

ENV_PREFIX="${1:-${CONDA_PREFIX:-}}"
if [ -z "${ENV_PREFIX}" ]; then
  echo "error: no env prefix found. Activate the dvbfixer env first or pass it as an argument."
  echo "       e.g. bash scripts/set-modeller-key.sh /opt/conda/envs/dvbfixer"
  exit 1
fi

# Find Modeller config.py inside the env
ENV_PYTHON="${ENV_PREFIX}/bin/python"
if [ ! -x "${ENV_PYTHON}" ]; then
  echo "error: environment Python is not executable: ${ENV_PYTHON}"
  exit 1
fi

CONFIG=$(find "${ENV_PREFIX}/lib" -maxdepth 4 -path "*/modeller*/modlib/modeller/config.py" 2>/dev/null | head -1)

if [ -z "${CONFIG}" ]; then
  echo "error: could not find modeller/config.py under ${ENV_PREFIX}/lib"
  echo "       is Modeller installed in this environment?"
  exit 1
fi

CONFIG_DIR=$(dirname "${CONFIG}")
NEW_CONFIG=$(mktemp "${CONFIG_DIR}/.config.py.new.XXXXXX")
OLD_CONFIG=$(mktemp "${CONFIG_DIR}/.config.py.old.XXXXXX")
trap 'rm -f "${NEW_CONFIG}" "${OLD_CONFIG}"' EXIT
cp "${CONFIG}" "${OLD_CONFIG}"

CONFIG_SOURCE="${CONFIG}" CONFIG_TARGET="${NEW_CONFIG}" "${ENV_PYTHON}" - <<'PY'
import os
import re
from pathlib import Path

source = Path(os.environ["CONFIG_SOURCE"]).read_text()
key = os.environ["MODELLER_LICENSE_KEY"]
line = "license = " + repr(key) + "\n"
updated, count = re.subn(r"^license\s*=.*$", line.rstrip(), source, flags=re.MULTILINE)
if count == 0:
    updated = source.rstrip() + "\n\n" + line
Path(os.environ["CONFIG_TARGET"]).write_text(updated)
PY
chmod --reference="${CONFIG}" "${NEW_CONFIG}" 2>/dev/null || true
mv "${NEW_CONFIG}" "${CONFIG}"

if ! "${ENV_PYTHON}" -c 'import modeller; modeller.environ()' >/dev/null 2>&1; then
  mv "${OLD_CONFIG}" "${CONFIG}"
  echo "error: Modeller rejected the supplied license; previous configuration restored."
  exit 1
fi
rm -f "${OLD_CONFIG}"
echo "Modeller license configured and validated."
