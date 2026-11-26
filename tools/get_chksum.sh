#!/usr/bin/env bash
# This script is used to get the checksum of a benchmark result.
# It must be run from inside a benchmark directory (e.g. ./counting).
# Usage: ./get_chksum.sh [extra_flags] <target>
#
# If no -DNPDP_RSEED=... is given in extra_flags, a fixed seed is used so that
# two invocations (e.g. original vs pluto) compute on the same input.

usage() {
    echo "Usage: $0 [extra_flags] <target>" >&2
    echo "Example: $0 \"-DN=200 -DNPDP_RSEED=21\" counting_pluto" >&2
    exit 1
}

if [[ $# -eq 1 ]]; then
    EXTRA_FLAGS=""
    target=$1
elif [[ $# -eq 2 ]]; then
    EXTRA_FLAGS=$1
    target=$2
else
    usage
fi

if [[ -z "$target" ]]; then
    usage
fi

if [[ "$EXTRA_FLAGS" != *"-DNPDP_RSEED"* ]]; then
    EXTRA_FLAGS="$EXTRA_FLAGS -DNPDP_RSEED=42"
fi
EXTRA_FLAGS="$EXTRA_FLAGS -DNPDP_DUMP"

dump_file=$(mktemp)
trap 'rm -f "$dump_file"' EXIT

# delete the old binary so that it is rebuilt with the new flags
rm -f "$target"

# compile the benchmark with dump option.
# PLUTO_QUIET is used instead of PLUTO_EXTRA_FLAGS so that per-benchmark
# pluto flags defined in the Makefiles are kept.
make "$target" PLUTO_QUIET="-q" EXTRA_FLAGS="$EXTRA_FLAGS" 1>&2 || exit 1

# run the benchmark (result is dumped on stderr)
"./$target" 2>"$dump_file" >/dev/null || exit 1

# get the checksum
sha256sum "$dump_file" | awk '{print $1}'
