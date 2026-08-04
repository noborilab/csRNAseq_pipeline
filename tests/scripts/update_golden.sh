#!/usr/bin/env bash
# Refresh the committed golden QC tables from a finished e2e run.
#
#   ./tests/scripts/update_golden.sh <e2e output dir>
#
# Only do this when the new values are known to be right: a golden mismatch means either a
# regression or a deliberate change (a metric redefinition, or an upgrade of HOMER, bwa or
# bfqutils). Inspect the diff first, and say which it was in the commit message.
set -euo pipefail
if [[ $# -ne 1 || ! -d "$1" ]]; then
    echo "usage: $0 <e2e output dir>" >&2
    exit 1
fi
REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
for f in qc_final_cs.txt qc_final_in.txt; do
    if [[ ! -f "$1/$f" ]]; then
        echo "ERROR: $1/$f not found; is that an e2e output directory?" >&2
        exit 1
    fi
    diff -u "$REPO_ROOT/tests/data/golden/$f" "$1/$f" || true
    cp "$1/$f" "$REPO_ROOT/tests/data/golden/$f"
    echo "updated tests/data/golden/$f"
done
