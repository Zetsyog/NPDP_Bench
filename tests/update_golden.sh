#!/usr/bin/env bash
# Regenerate the reference checksums in tests/golden/.
#
# Only run this when a change of results is intended (e.g. a kernel bug fix),
# and review the diff: the golden files are what `make test` compares against.

set -eu

ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
source "$ROOT/tests/common.sh"

BUILD=$(mktemp -d)
trap 'rm -rf "$BUILD"' EXIT
mkdir -p "$ROOT/tests/golden"

for b in "${NPDP_BENCHES[@]}"; do
    golden_hash "$b" "$BUILD" -O0 >"$ROOT/tests/golden/$b.sha256"
    echo "$b: $(cat "$ROOT/tests/golden/$b.sha256")"
done
