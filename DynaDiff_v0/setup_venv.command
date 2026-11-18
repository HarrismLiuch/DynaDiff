#!/bin/bash
set -euo pipefail

project_dir="$(cd -- "$(dirname -- "$0")" && pwd -P)"
venv_dir="$project_dir/.venv"

if [[ ! -x "$venv_dir/bin/python" ]]; then
    candidates=(python3 python3.14 python3.13 python3.12 python3.11 python3.10
                /opt/homebrew/bin/python3 /usr/local/bin/python3 /opt/anaconda3/bin/python3)
    if [[ -n "${DYNADIFF_PYTHON:-}" ]]; then
        candidates=("$DYNADIFF_PYTHON")
    fi
    interpreter=""
    for candidate in "${candidates[@]}"; do
        if command -v "$candidate" >/dev/null 2>&1 &&
            "$candidate" -c 'import sys, venv; sys.exit(sys.version_info < (3, 10))' 2>/dev/null; then
            interpreter="$candidate"
            break
        fi
    done
    if [[ -z "$interpreter" ]]; then
        echo "DynaDiff requires Python 3.10 or later with the venv module." >&2
        echo "Install Python from python.org, then run setup_venv.command again." >&2
        exit 1
    fi
    echo "Creating the DynaDiff environment with $interpreter..."
    "$interpreter" -m venv "$venv_dir"
fi

"$venv_dir/bin/python" -m pip --disable-pip-version-check install -r "$project_dir/requirements.txt"
"$venv_dir/bin/python" -c 'import ase, gemmi, numpy, scipy.linalg, matplotlib, cmasher, cmcrameri.cm; from PySide6 import QtWidgets'
touch "$venv_dir/.dynadiff-ready"
echo "DynaDiff is ready. Open DynaDiff.command to launch."
