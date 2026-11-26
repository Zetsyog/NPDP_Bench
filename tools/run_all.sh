#!/usr/bin/env bash
# Run every benchmark x variant x thread count and write one tidy CSV.
#
# Usage: tools/run_all.sh [-o results.csv]
#
# For each benchmark and size:
#   1. the original program is built with NPDP_DUMP to get a reference result;
#   2. each variant is built with NPDP_DUMP and its result compared to the
#      reference (tools/compare_dumps.sh) -> column `correct`;
#   3. each variant is built with NPDP_TIME and, for each thread count, run
#      once as warm-up then NB_RUNS times; the median of the NB_MEDIAN_RUNS
#      central runs and their standard deviation are reported.
#
# Everything is built in a temporary copy of the sources, the source tree is
# left untouched.
#
# Configuration (environment variables):
#   BENCHES         benchmarks to run (default: all)
#   VARIANTS        "orig" and/or "pluto" (default: orig, plus pluto if PLUTO is set)
#   THREADS         thread counts for parallel variants (default: "1 2 4 ... nproc")
#   SIZES           per-benchmark sizes, e.g. "counting=2000 sw3d=150"
#                   (default: see DEFAULT_SIZES below)
#   NB_RUNS         timed runs per measurement (default 5)
#   NB_MEDIAN_RUNS  central runs used for the median/stddev (default 3)
#   STDDEV_TRESHOLD abort a measurement if stddev exceeds this (default 1e9 = never)
#   TIMEOUT         timeout of one run, in seconds (default 600)
#   SEED            NPDP_RSEED used for every build (default 42)
#   FLUSH_CACHE     1 to build with -DNPDP_FLUSH_CACHE (default 0)
#   EXTRA_DEFS      extra -D flags for every build, recorded in the CSV
#                   (e.g. EXTRA_DEFS=-DNPDP_PAIR_TABLE)
#   CFLAGS          compiler flags (default: same as common.mk)
#   PLUTO           path to polycc (enables the pluto variant)
#   PLUTO_FLAGS     pluto flags (default: same as common.mk)
#   OMP_PROC_BIND, OMP_PLACES  thread pinning (default: close, cores)

set -u

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
ROOT=$(cd -- "$SCRIPT_DIR/.." && pwd)
source "$SCRIPT_DIR/utils.sh"

OUTPUT="results.csv"
while getopts "o:h" opt; do
    case $opt in
    o) OUTPUT=$OPTARG ;;
    *)
        sed -n '2,34p' "$0" >&2
        exit 1
        ;;
    esac
done

# Sizes giving sequential runs of roughly 1-5 s (measured on an Intel Core
# i5-1235U laptop, gcc 16 -O3). Override with SIZES.
declare -A DEFAULT_SIZES=(
    [counting]=2000 [counting_im]=2000 [knuth]=2000 [mcc]=1500
    [mcm]=2000 [mea]=350 [nussinov]=1500 [nw]=2000 [sw]=2000
    [sw3d]=220 [triang]=1800 [zuker]=600
)

BENCHES=${BENCHES:-"counting counting_im knuth mcc mcm mea nussinov nw sw sw3d triang zuker"}
if [[ -z "${VARIANTS:-}" ]]; then
    VARIANTS="orig"
    [[ -n "${PLUTO:-}" ]] && VARIANTS="orig pluto"
fi
if [[ -z "${THREADS:-}" ]]; then
    THREADS=""
    nproc=$(nproc)
    for ((t = 1; t < nproc; t *= 2)); do THREADS="$THREADS $t"; done
    THREADS="$THREADS $nproc"
fi
NB_RUNS=${NB_RUNS:-5}
NB_MEDIAN_RUNS=${NB_MEDIAN_RUNS:-3}
STDDEV_TRESHOLD=${STDDEV_TRESHOLD:-1e9}
TIMEOUT=${TIMEOUT:-600}
SEED=${SEED:-42}
FLUSH_CACHE=${FLUSH_CACHE:-0}
CFLAGS=${CFLAGS:-"-march=native -O3 -fopenmp -ffp-contract=off"}
PLUTO_FLAGS=${PLUTO_FLAGS:-"--tile --parallel"}
export OMP_PROC_BIND=${OMP_PROC_BIND:-close}
export OMP_PLACES=${OMP_PLACES:-cores}

declare -A SIZE
for b in $BENCHES; do SIZE[$b]=${DEFAULT_SIZES[$b]:-1000}; done
for kv in ${SIZES:-}; do SIZE[${kv%%=*}]=${kv##*=}; done

if [[ " $VARIANTS " == *" pluto "* && -z "${PLUTO:-}" ]]; then
    log "VARIANTS contains pluto but PLUTO is not set (export PLUTO=/path/to/polycc)"
    exit 1
fi

# ---------------------------------------------------------------------------
# Metadata recorded in every row
# ---------------------------------------------------------------------------
csv() { # quote a CSV field
    local s=${1//\"/\'}
    printf '"%s"' "$s"
}
META_TIMESTAMP=$(date -Iseconds)
META_VERSION=$(cat "$ROOT/VERSION" 2>/dev/null || echo unknown)
META_COMMIT=$(git -C "$ROOT" rev-parse --short HEAD 2>/dev/null || echo unknown)
META_DIRTY=$([[ -n "$(git -C "$ROOT" status --porcelain 2>/dev/null)" ]] && echo yes || echo no)
META_HOST=$(uname -n)
META_CPU=$(grep -m1 "model name" /proc/cpuinfo 2>/dev/null | cut -d: -f2- | sed 's/^ *//')
META_CPU=${META_CPU:-$(uname -m)}
META_GOVERNOR=$(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor 2>/dev/null || echo n/a)
META_CC=$(gcc --version | head -1)
META_PLUTO=n/a
[[ -n "${PLUTO:-}" ]] && META_PLUTO=$("$PLUTO" --version 2>&1 | head -1)

if [[ "$META_DIRTY" == yes ]]; then
    log "warning: the source tree has uncommitted changes (recorded as git_dirty=yes)"
fi
if [[ "$META_GOVERNOR" != performance && "$META_GOVERNOR" != n/a ]]; then
    log "warning: CPU frequency governor is '$META_GOVERNOR' (use 'performance' for stable timings)"
fi

echo "timestamp,npdp_version,git_commit,git_dirty,host,cpu_model,governor,cc_version,pluto_version,cflags,pluto_flags,flush_cache,seed,extra_defs,bench,variant,n,threads,omp_proc_bind,omp_places,median_time,stddev,nb_runs,nb_median_runs,correct" >"$OUTPUT"
log "Writing results to $OUTPUT"

# ---------------------------------------------------------------------------
# Build helpers (in a temporary copy of the sources)
# ---------------------------------------------------------------------------
BUILD=$(mktemp -d)
trap 'rm -rf "$BUILD"' EXIT
cp -r "$ROOT/include" "$ROOT/common.mk" "$BUILD/"

# make_variant <bench> <variant> <extra flags> -> prints the binary path
make_variant() {
    local bench=$1 variant=$2 flags=$3 target=$1
    [[ $variant == pluto ]] && target="${bench}_pluto"
    rm -f "$BUILD/$bench/$target"
    if ! make -C "$BUILD/$bench" "$target" CFLAGS="$CFLAGS" PLUTO="${PLUTO:-}" PLUTO_FLAGS="$PLUTO_FLAGS" PLUTO_QUIET=-q \
        EXTRA_FLAGS="$flags" >"$BUILD/make.log" 2>&1; then
        log "build failed: $bench $variant"
        cat "$BUILD/make.log" >&2
        return 1
    fi
    echo "$BUILD/$bench/$target"
}

# run_dump <binary> <output file>
run_dump() {
    (cd "$(dirname "$1")" && timeout "$TIMEOUT" "$1" 2>"$2" >/dev/null)
}

# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------
for bench in $BENCHES; do
    n=${SIZE[$bench]}
    rm -rf "${BUILD:?}/$bench"
    cp -r "$ROOT/$bench" "$BUILD/"
    base_flags="-DN=$n -DNPDP_RSEED=$SEED ${EXTRA_DEFS:-}"
    log "=== $bench (n=$n)"

    ref_bin=$(make_variant "$bench" orig "$base_flags -DNPDP_DUMP") || continue
    if ! run_dump "$ref_bin" "$BUILD/ref.dump"; then
        log "reference run failed or timed out: $bench"
        continue
    fi

    for variant in $VARIANTS; do
        correct=no
        if bin=$(make_variant "$bench" "$variant" "$base_flags -DNPDP_DUMP") && run_dump "$bin" "$BUILD/var.dump" &&
            bash "$SCRIPT_DIR/compare_dumps.sh" "$BUILD/ref.dump" "$BUILD/var.dump"; then
            correct=yes
        fi
        [[ $correct == no ]] && log "warning: $bench $variant result differs from the reference"

        time_flags="$base_flags -DNPDP_TIME"
        [[ "$FLUSH_CACHE" == 1 ]] && time_flags="$time_flags -DNPDP_FLUSH_CACHE"
        bin=$(make_variant "$bench" "$variant" "$time_flags") || continue

        # the original program is sequential: thread count is irrelevant
        threads_list=$THREADS
        [[ $variant == orig ]] && threads_list=1

        for t in $threads_list; do
            export OMP_NUM_THREADS=$t
            # warm-up run (discarded)
            (cd "$BUILD/$bench" && timeout "$TIMEOUT" "$bin" >/dev/null 2>&1)
            result=$(run_bench "$bin" "$NB_RUNS" "$NB_MEDIAN_RUNS" "$STDDEV_TRESHOLD" "$TIMEOUT" "$BUILD/$bench")
            status=$?
            if [[ $status -ne 0 ]]; then
                log "measurement failed: $bench $variant threads=$t"
                continue
            fi
            read -r median stddev <<<"$result"
            log "$(printf '%-12s %-6s n=%-6s threads=%-3s %10s s  (stddev %s)' "$bench" "$variant" "$n" "$t" "$median" "$stddev")"

            {
                csv "$META_TIMESTAMP"; printf ','; csv "$META_VERSION"; printf ','; csv "$META_COMMIT"; printf ','
                csv "$META_DIRTY"; printf ','; csv "$META_HOST"; printf ','; csv "$META_CPU"; printf ','
                csv "$META_GOVERNOR"; printf ','; csv "$META_CC"; printf ','; csv "$META_PLUTO"; printf ','
                csv "$CFLAGS"; printf ','; csv "$PLUTO_FLAGS"; printf ','; csv "$FLUSH_CACHE"; printf ','
                csv "$SEED"; printf ','; csv "${EXTRA_DEFS:-}"; printf ','; csv "$bench"; printf ','; csv "$variant"; printf ','
                printf '%s,%s,' "$n" "$t"
                csv "$OMP_PROC_BIND"; printf ','; csv "$OMP_PLACES"; printf ','
                printf '%s,%s,%s,%s,' "$median" "$stddev" "$NB_RUNS" "$NB_MEDIAN_RUNS"
                csv "$correct"; printf '\n'
            } >>"$OUTPUT"
        done
    done
done

log "Done: $OUTPUT"
