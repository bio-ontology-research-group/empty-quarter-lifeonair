#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
bash lifeonair/run_lifeonair.sh --verify
(cd lifeonair/manuscript && latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex > /dev/null)
if grep -Eq 'Citation .* undefined|Reference .* undefined|There were undefined references|There were undefined citations' lifeonair/manuscript/main.log; then
  echo 'Unresolved manuscript references' >&2
  exit 1
fi
echo 'REPRODUCIBILITY: frozen-evidence checks and manuscript build PASS'
