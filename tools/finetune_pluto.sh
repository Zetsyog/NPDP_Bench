#!/usr/bin/env bash
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
source "$SCRIPT_DIR/utils.sh"

# This script is used to fine-tune a given program optimized by Pluto.
# It will do a exhaustive search of the best tile sizes
# for the given program

usage() {
    echo "Usage: $0 <source.c> <output_csv> <config>"
    exit 1
}

if [[ $# -ne 3 ]]; then
    usage
fi

SOURCE=$1
SRC_DIR=$(dirname "$SOURCE")
CSV_OUTPUT=$(realpath $2)
CONFIG=$3

BENCH_FILENAME=$(basename -- "$SOURCE")
BENCH_EXT="${BENCH_FILENAME##*.}"
BENCH_FILENAME="${BENCH_FILENAME%.*}"
BENCH_NAME="$BENCH_FILENAME"

# pluto environment
if [ -z ${POLYCC+x} ]; then
    log "POLYCC is not set"
    exit 1
fi

if [ -z ${PLUTO_EXTRA_FLAGS+x} ]; then
    log "PLUTO_EXTRA_FLAGS is not set"
    exit 1
fi

# Tile sizes
if [ -z ${TILESIZE_LVL_1+x} ]; then
    log "TILESIZE_LVL_1 is not set"
    exit 1
fi
if [ -z ${TILESIZE_LVL_2+x} ]; then
    log "TILESIZE_LVL_2 is not set"
    exit 1
fi
if [ -z ${TILESIZE_LVL_3+x} ]; then
    log "TILESIZE_LVL_3 is not set"
    exit 1
fi

if [ -z ${CFLAGS+x} ]; then
    log "CFLAGS is not set"
    exit 1
fi

if [ -z ${EXTRA_FLAGS+x} ]; then
    log "EXTRA_FLAGS is not set"
    exit 1
fi
EXTRA_FLAGS="$EXTRA_FLAGS -lm -DMETRICS_TIME"

if [ -z ${INCLUDE_DIR+x} ]; then
    log "INCLUDE_DIR is not set"
    exit 1
fi

if [ -z ${NB_MEDIAN_RUNS+x} ]; then
    log "NB_MEDIAN_RUNS is not set"
    exit 1
fi

if [ -z ${STDDEV_TRESHOLD+x} ]; then
    log "STDDEV_TRESHOLD is not set"
    exit 1
fi

if [ -z ${TIMEOUT+x} ]; then
    log "TIMEOUT is not set"
    exit 1
fi

# Check if output file exists
log "Writing result data to $CSV_OUTPUT"
if [[ -f $CSV_OUTPUT ]]; then
    # if it exists, check if we must override it
    if [[ $OVERRIDE_RESULTS -eq 0 ]]; then
        log "Output data file already exists. Aborting."
        exit 1
    else
        log "Output directory already contains result file. Overriding."
        echo "" >"$CSV_OUTPUT"
    fi
fi

echo "DIV1,DIV2,DIV3,EXEC_TYPE,THREAD_NB,TIME,DATASET_SIZE,CONFIG,CFLAGS,EXTRA_FLAGS,STDDEV" >$CSV_OUTPUT
log "==============================="
log "Finetuning $SOURCE"
log "==============================="
log "Polybench dir is $POLYBENCH_DIR"
log "Configuration name $CONFIG"
log "CONFIG DUMP"
log "==============================="
print_config
log "==============================="

BENCH_SRC=$(realpath $SOURCE)
EXEC="./bench.out"

printf "%-40s%-8s%-8s%-8s%-12s%-12s%-12s\n" "NAME" "DIV0" "DIV1" "DIV2" "TIME" "STD_DEV" "BEST_TIME" 1>&2

best_exec_time=999
best_tiles=""

for d1 in ${TILESIZE_LVL_1}; do
    for d2 in ${TILESIZE_LVL_2}; do
        for d3 in ${TILESIZE_LVL_3}; do
            pluto_src="${SRC_DIR}/bench.pluto.c"
            echo "$d1" >tile.sizes
            echo "$d2" >>tile.sizes
            echo "$d3" >>tile.sizes
            # Tile program
            ${POLYCC} -q ${PLUTO_EXTRA_FLAGS} "${BENCH_SRC}" -o "$pluto_src" >/dev/null 2>&1
            if [[ $? -ne 0 ]]; then
                log "Error while generating tiled program"
                log "${POLYCC} ${PLUTO_EXTRA_FLAGS} ${BENCH_SRC} -o \"$pluto_src\""
                exit 1
            fi
            rm -f tile.sizes ./*.cloog ./srcfilename ./regtile

            # Compile the program
            gcc ${CFLAGS} -I "$SRC_DIR" -I "${INCLUDE_DIR}" "$pluto_src" -o "${EXEC}" ${EXTRA_FLAGS}
            if [[ $? -ne 0 ]]; then
                log "Error while compiling program"
                log "gcc ${CFLAGS} -I $SRC_DIR -I ${INCLUDE_DIR} $pluto_src -o ${EXEC} ${EXTRA_FLAGS}"
                exit 1
            fi
            rm -f "$pluto_src"

            bench_result=$(run_bench "${EXEC}" "$NB_RUNS" "$NB_MEDIAN_RUNS" "$STDDEV_TRESHOLD" "$TIMEOUT" "$SRC_DIR")
            status=$?

            rm -f "$EXEC"

            read -r exec_time std_dev <<<"$bench_result"

            printf "%-40s%-8u%-8u%-8u" "${BENCH_NAME}" "$d1" "$d2" "$d3" 1>&2

            # if execution failed, we print error and continue
            if [[ $status -ne 0 ]]; then
                tmp=$(echo "$bench_result" | head -n 1)
                echo -n -e "\terror($status) $tmp" 1>&2
                log ""
                continue
            fi
            printf "%-12f%-12f" "${exec_time}" "${std_dev}" 1>&2
            better=$(echo "${exec_time} < ${best_exec_time}" | bc -l)
            # if execution time is better than previous best, update best
            if [[ "$better" -eq 1 ]]; then
                best_exec_time=$exec_time
                best_tiles="$d1 x $d2 x $d3"
                printf "%-12f - new best!" "${best_exec_time}" 1>&2
            else
                printf "%-12f" "${best_exec_time}" 1>&2
            fi
            log ""

            echo "$d1,$d2,$d3,$src_name,$thread_nb,$exec_time,$DATASET_SIZE,$CONFIG,$CFLAGS,$EXTRA_FLAGS,$std_dev" >>$CSV_OUTPUT

        done
    done
done

log "** RESULTS for $BENCH_NAME with $OMP_NUM_THREADS threads"
printf "**** time = %-10s - TILE_SIZES = %-15s ****\n" "${best_exec_time}" "${best_tiles}" 1>&2
log ""
