#!/usr/bin/env bash
# NPDP Bench test suite. Run from anywhere: `make test` or `bash tests/run_tests.sh`.
#
#  1. warnings   every benchmark builds with -Wall -Wextra -Werror
#  2. sanitizers ASan + UBSan (+ LSan) report nothing at N=40
#  3. golden     result dumps at -O0 and -O3 -march=native match the reference
#                checksums in tests/golden/ (regenerate with `make golden`)
#  4. answers    known-answer tests on hand-checkable inputs (known_answers.c)
#  5. pair table the RNA kernels built with -DNPDP_PAIR_TABLE give the same
#                golden checksums and known answers
#  6. valgrind   no invalid or uninitialized accesses at N=16 (skipped if
#                valgrind is not installed)
#
# Environment: CC (default gcc), ASAN_OPTIONS (default detect_leaks=1; LSan
# does not work in some sandboxes, set detect_leaks=0 there).

set -u

ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
source "$ROOT/tests/common.sh"

CC=${CC:-gcc}
export ASAN_OPTIONS=${ASAN_OPTIONS:-detect_leaks=1}
BUILD=$(mktemp -d)
trap 'rm -rf "$BUILD"' EXIT

failures=0
pass() { printf '  %-12s ok\n' "$1"; }
fail() {
    printf '  %-12s FAIL %s\n' "$1" "$2"
    failures=$((failures + 1))
}

# build <bench> <output> <flags...>
build() {
    local bench=$1 out=$2
    shift 2
    $CC "$@" -fopenmp -ffp-contract=off -I "$ROOT/include" -I "$ROOT/$bench" -o "$out" "$ROOT/$bench/$bench.c" -lm
}

echo "== warnings (-Wall -Wextra -Werror)"
for b in "${NPDP_BENCHES[@]}"; do
    if log=$(build "$b" "$BUILD/$b.w" -O2 -Wall -Wextra -Werror -Wno-unknown-pragmas -DN=40 -DNPDP_DUMP -DNPDP_TIME 2>&1); then
        pass "$b"
    else
        fail "$b" "$(echo "$log" | grep -m1 -E 'error|warning')"
    fi
done

echo "== sanitizers (ASan/UBSan, N=40)"
for b in "${NPDP_BENCHES[@]}"; do
    if ! build "$b" "$BUILD/$b.san" -O1 -g -fsanitize=address,undefined -fno-sanitize-recover=all -DN=40 -DNPDP_RSEED=1 -DNPDP_DUMP 2>"$BUILD/$b.cerr"; then
        fail "$b" "build failed: $(head -1 "$BUILD/$b.cerr")"
        continue
    fi
    if (cd "$BUILD" && "./$b.san" >/dev/null 2>"$b.san.log"); then
        if grep -qE "runtime error|ERROR: (Address|Leak)Sanitizer" "$BUILD/$b.san.log"; then
            fail "$b" "$(grep -m1 -E 'runtime error|ERROR:' "$BUILD/$b.san.log")"
        else
            pass "$b"
        fi
    else
        fail "$b" "exit $?: $(grep -m1 -E 'runtime error|ERROR:' "$BUILD/$b.san.log")"
    fi
done

echo "== golden checksums (N=$NPDP_GOLDEN_N, seed $NPDP_GOLDEN_SEED, -O0 and -O3)"
for b in "${NPDP_BENCHES[@]}"; do
    golden_file="$ROOT/tests/golden/$b.sha256"
    if [[ ! -f "$golden_file" ]]; then
        fail "$b" "missing $golden_file (run: make golden)"
        continue
    fi
    expected=$(cat "$golden_file")
    ok=1
    for opt in "-O0" "-O3 -march=native"; do
        # shellcheck disable=SC2086
        got=$(golden_hash "$b" "$BUILD" $opt)
        if [[ "$got" != "$expected" ]]; then
            fail "$b" "$opt: $got != $expected"
            ok=0
            break
        fi
    done
    [[ $ok -eq 1 ]] && pass "$b"
done

echo "== known answers"
for t in MCM NUSSINOV SW NW TRIANG; do
    name=$(echo "$t" | tr 'A-Z' 'a-z')
    if ! $CC -O1 -g -fsanitize=address,undefined -fno-sanitize-recover=all -I "$ROOT/include" -DN=16 -DNPDP_NO_MAIN -DTEST_$t \
        -o "$BUILD/ka_$name" "$ROOT/tests/known_answers.c" -lm 2>"$BUILD/ka_$name.cerr"; then
        fail "$name" "build failed: $(head -1 "$BUILD/ka_$name.cerr")"
        continue
    fi
    if out=$("$BUILD/ka_$name" 2>&1); then
        pass "$name"
    else
        fail "$name" "$(echo "$out" | head -1)"
    fi
done

echo "== pair lookup table (-DNPDP_PAIR_TABLE gives the same results)"
for b in "${NPDP_PAIR_BENCHES[@]}"; do
    expected=$(cat "$ROOT/tests/golden/$b.sha256" 2>/dev/null)
    ok=1
    for opt in "-O0" "-O3 -march=native"; do
        # shellcheck disable=SC2086
        got=$(golden_hash "$b" "$BUILD" $opt -DNPDP_PAIR_TABLE)
        if [[ "$got" != "$expected" ]]; then
            fail "$b" "$opt -DNPDP_PAIR_TABLE: $got != $expected"
            ok=0
            break
        fi
    done
    [[ $ok -eq 1 ]] && pass "$b"
done
if $CC -O1 -g -fsanitize=address,undefined -fno-sanitize-recover=all -I "$ROOT/include" -DN=16 -DNPDP_NO_MAIN -DTEST_NUSSINOV -DNPDP_PAIR_TABLE \
    -o "$BUILD/ka_nussinov_table" "$ROOT/tests/known_answers.c" -lm 2>"$BUILD/ka_table.cerr" && out=$("$BUILD/ka_nussinov_table" 2>&1); then
    pass "answers"
else
    fail "answers" "nussinov known answers with -DNPDP_PAIR_TABLE: $(head -1 "$BUILD/ka_table.cerr") ${out:-}"
fi

if command -v valgrind >/dev/null 2>&1; then
    echo "== valgrind (N=16)"
    for b in "${NPDP_BENCHES[@]}"; do
        build "$b" "$BUILD/$b.vg" -O0 -g -DN=16 -DNPDP_RSEED=1 -DNPDP_DUMP 2>/dev/null
        if (cd "$BUILD" && valgrind -q --error-exitcode=9 "./$b.vg" >/dev/null 2>"$b.vg.log"); then
            pass "$b"
        else
            fail "$b" "$(grep -m1 -E 'Invalid|uninitialised' "$BUILD/$b.vg.log")"
        fi
    done
else
    echo "== valgrind: skipped (not installed)"
fi

echo
if [[ $failures -eq 0 ]]; then
    echo "All tests passed."
else
    echo "$failures test(s) failed."
    exit 1
fi
