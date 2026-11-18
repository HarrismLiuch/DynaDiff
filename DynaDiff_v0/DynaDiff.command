#!/bin/bash
set -euo pipefail

project_dir="$(cd -- "$(dirname -- "$0")" && pwd -P)"
venv_dir="$project_dir/.venv"
if [[ ! -x "$venv_dir/bin/python" || ! -f "$venv_dir/.dynadiff-ready" ||
      "$project_dir/requirements.txt" -nt "$venv_dir/.dynadiff-ready" ]]; then
    "$project_dir/setup_venv.command"
fi
exec "$venv_dir/bin/python" "$project_dir/DynaDiff.py" "$@"
