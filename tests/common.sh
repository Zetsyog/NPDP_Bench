#!/usr/bin/env bash
# Shared definitions for the test scripts (sourced, not executed).

NPDP_BENCHES=(counting counting_im knuth mcc mcm mea nussinov nw sw sw3d triang zuker)
# Benchmarks using NPDP_PAIRED (base pairing, optional lookup table)
NPDP_PAIR_BENCHES=(counting counting_im mcc mea nussinov)

# Golden checksums are computed at this size and seed.
NPDP_GOLDEN_N=64
NPDP_GOLDEN_SEED=1

# golden_hash <bench> <build_dir> <optimization flags...>
# Print the sha256 of the result dump of <bench> built with the given flags.
golden_hash() {
    local bench=$1 build_dir=$2
    shift 2
    local bin="$build_dir/$bench.golden"
    ${CC:-gcc} "$@" -fopenmp -ffp-contract=off -I "$ROOT/include" -I "$ROOT/$bench" \
        -DN="$NPDP_GOLDEN_N" -DNPDP_RSEED="$NPDP_GOLDEN_SEED" -DNPDP_DUMP \
        -o "$bin" "$ROOT/$bench/$bench.c" -lm || return 1
    (cd "$build_dir" && "$bin" 2>&1 >/dev/null) | sha256sum | awk '{print $1}'
}
