#!/usr/bin/env bash
# Compare two result dumps (NPDP_DUMP output) with a numerical tolerance.
#
# Usage: compare_dumps.sh <reference> <candidate> [rtol] [atol]
#
# Both files are split into whitespace separated tokens. Tokens that are both
# numbers must satisfy |a - b| <= atol + rtol * max(|a|, |b|); any other token
# must be identical. Exit status 0 if the dumps match, 1 otherwise (the first
# mismatch is printed on stderr), 2 on usage error.
#
# Defaults: rtol = 1e-6, atol = 1e-5. Dumps print floating point values with
# 3 to 5 decimals, so a transformation that only reassociates floating point
# operations can flip the last printed digit; use rtol=0 atol=0 to require
# an exact match.

if [[ $# -lt 2 || $# -gt 4 ]]; then
    echo "Usage: $0 <reference> <candidate> [rtol=1e-6] [atol=1e-5]" >&2
    exit 2
fi

ref=$1
cand=$2
rtol=${3:-1e-6}
atol=${4:-1e-5}

for f in "$ref" "$cand"; do
    if [[ ! -r "$f" ]]; then
        echo "compare_dumps: cannot read $f" >&2
        exit 2
    fi
done

awk -v rtol="$rtol" -v atol="$atol" -v ref="$ref" -v cand="$cand" '
function isnum(s) { return s ~ /^[-+]?([0-9]+\.?[0-9]*|\.[0-9]+)([eE][-+]?[0-9]+)?$/ }
function abs(x) { return x < 0 ? -x : x }
BEGIN {
    n = 0
    while ((getline line < ref) > 0) {
        k = split(line, t)
        for (i = 1; i <= k; i++) a[++n] = t[i]
    }
    m = 0
    while ((getline line < cand) > 0) {
        k = split(line, t)
        for (i = 1; i <= k; i++) b[++m] = t[i]
    }
    if (n != m) {
        printf("compare_dumps: %d tokens in %s, %d in %s\n", n, ref, m, cand) > "/dev/stderr"
        exit 1
    }
    for (i = 1; i <= n; i++) {
        if (isnum(a[i]) && isnum(b[i])) {
            x = a[i] + 0; y = b[i] + 0
            tol = atol + rtol * (abs(x) > abs(y) ? abs(x) : abs(y))
            if (abs(x - y) > tol) {
                printf("compare_dumps: token %d differs: %s vs %s (tolerance %g)\n", i, a[i], b[i], tol) > "/dev/stderr"
                exit 1
            }
        } else if (a[i] != b[i]) {
            printf("compare_dumps: token %d differs: \"%s\" vs \"%s\"\n", i, a[i], b[i]) > "/dev/stderr"
            exit 1
        }
    }
    exit 0
}'
