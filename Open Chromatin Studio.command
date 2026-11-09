#!/bin/sh
set -eu
studio_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
cd "$studio_dir"
if [ -x "$studio_dir/.venv/bin/python" ]; then
    studio_python="$studio_dir/.venv/bin/python"
elif [ -x "$studio_dir/../../work/sim-env/bin/python" ]; then
    studio_python="$studio_dir/../../work/sim-env/bin/python"
else
    studio_python=python3
fi
if ! "$studio_python" -c 'import numpy, scipy' >/dev/null 2>&1; then
    echo "Install the local environment first; see README.md for setup instructions."
    exit 1
fi
exec "$studio_python" run.py --open "$@"
