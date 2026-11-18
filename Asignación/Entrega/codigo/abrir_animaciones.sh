#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
if [[ -n "${PYTHON:-}" ]]; then
    interprete="$PYTHON"
elif [[ -x .venv/bin/python ]]; then
    interprete=.venv/bin/python
elif [[ "$PWD" == */Entrega/codigo && -x ../../.venv/bin/python ]]; then
    interprete="$(cd ../.. && pwd -P)/.venv/bin/python"
else
    echo 'Cree el entorno: python -m venv .venv' >&2
    echo 'Instale: .venv/bin/python -m pip install -r requirements-vispy.txt' >&2
    exit 1
fi
ulimit -c 0
exec "$interprete" animaciones.py "$@"
