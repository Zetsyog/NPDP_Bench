#!/usr/bin/env bash

log() { echo "$@" 1>&2; }

# usage: compute_mean_exec_time [result file] [total number of runs] [number of median runs to use] [stddev treshold]
compute_mean_exec_time() {
    local total_runs="$2"
    local median_runs="$3"
    local treshold="$4"

    # sort the result file
    sort -n $1 >$1.sorted

    local first_line="$(expr $total_runs / 2 - $median_runs / 2)"
    local last_line="$(expr $total_runs / 2 + $median_runs / 2)"

    # extract the selected lines
    sed -n "${first_line},${last_line}p" "$1.sorted" >"$1.median"

    # compute mean
    local mean=$(awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' "$1.median")

    # compute standard deviation
    local stddev=$(awk '{delta = $1 - avg; avg += delta / NR; mean2 += delta * ($1 - avg); } END { print sqrt(mean2 / NR); }' "$1.median")

    # check if the standard deviation is too high
    if ! awk "BEGIN{ exit ($stddev > $treshold) }"; then
        log "Standard deviation is too high: ${stddev} > ${treshold}"
        log "Aborting experiment"
        return 1
    fi

    # cleanup
    rm "$1.sorted" "$1.median"

    echo "$mean" "$stddev"
}

# usage: run_bench [program] [number of runs] [number of median runs to use] [stddev treshold] [timeout]
# returns 0 if ok
# returns 1 if an error occured
# echo the execution time
run_bench() {
    local prog=$1
    local nb_runs=$2
    local nb_median_runs=$3
    local stddev_treshold=$4
    local timeout=$5
    local debug_dir=${6:-"./"}

    local result_file="${debug_dir}/____exec_times.tmp"
    local debug_file="${debug_dir}/____debug.tmp"
    local error_file="${debug_dir}/error.log"
    local exec_time

    rm -f $result_file $debug_file

    if [[ "$nb_runs" -eq 1 ]]; then
        timeout $timeout $prog 2>$debug_file >>$result_file
        prog_status=$?
        if [[ "$prog_status" -eq 124 ]]; then
            log "Program timed out"
            echo "$timeout" "0.0"
            return 0
        fi
        if [[ "$prog_status" -ne 0 ]]; then
            cat $debug_file
            echo "======================================" >$error_file
            cat $debug_file >$error_file
            return 1
        fi
        exec_time=$(cat $result_file)
        rm -f $result_file $debug_file

        echo "$exec_time" "0.0"
    else
        # run program nb_runs times
        # store execution times in a temp file
        for ((i = 0; i < nb_runs; i++)); do
            timeout $timeout $prog 2>$debug_file >>$result_file
            prog_status=$?
            if [[ "$prog_status" -eq 124 ]]; then
                log "Program timed out"
                echo "$timeout" "0.0"
                return 0
            fi
            if [[ "$prog_status" -ne 0 ]]; then
                cat $debug_file
                echo "======================================" >$error_file
                cat $debug_file >$error_file
                return 1
            fi
        done
        compute_mean_exec_time "$result_file" "$nb_runs" "$nb_median_runs" "$stddev_treshold"
        if [[ "$?" -ne 0 ]]; then
            return 1
        fi
    fi
    rm -f $result_file $debug_file
    return 0
}

# unset all env variables used by the scripts
reset_config() {
    unset GEN_POLYCC
    unset POLYCC
    unset PLUTO_EXTRA_FLAGS
    unset GEN_OMP_SCHEDULES
    unset DATASET_SIZE
    unset TILESIZE_LVL_1
    unset TILESIZE_LVL_2
    unset TILESIZE_LVL_3
    unset THREAD_NB_LIST
    unset CFLAGS
    unset EXTRA_FLAGS
    unset CC
    unset DRYRUN
    unset OVERRIDE_RESULTS
    unset POLYBENCH_DIR

    unset NB_MEDIAN_RUNS
    unset NB_RUNS
    unset STDDEV_TRESHOLD
    unset TIMEOUT

    # OpenMP environment variables
    unset OMP_PLACES
    unset OMP_WAIT_POLICY
    unset OMP_DYNAMIC
    unset OMP_STATIC
    unset OMP_PROC_BIND

}

print_config() {
    echo -e "# Pluto configuration"
    echo -e "GEN_POLYCC=$GEN_POLYCC"
    echo -e "\tPOLYCC=$POLYCC"
    echo -e "\tPLUTO_EXTRA_FLAGS=$PLUTO_EXTRA_FLAGS"

    echo -e "# Polybench configration"
    echo -e "DATASET_SIZE=$DATASET_SIZE"
    echo -e "THREAD_NB_LIST=$THREAD_NB_LIST"
    echo -e "TILESIZE_LVL_1=$TILESIZE_LVL_1"
    echo -e "TILESIZE_LVL_2=$TILESIZE_LVL_2"
    echo -e "TILESIZE_LVL_3=$TILESIZE_LVL_3"

    echo -e "# Compiler configuration"
    echo -e "CC=$CC"
    echo -e "\tCFLAGS=$CFLAGS"
    echo -e "\tEXTRA_FLAGS=$EXTRA_FLAGS"

    echo -e "# Misc"
    echo -e "GEN_OMP_SCHEDULES=$GEN_OMP_SCHEDULES"
    echo -e "DRYRUN=$DRYRUN"
    echo -e "OVERRIDE_RESULTS=$OVERRIDE_RESULTS"
    echo -e "POLYBENCH_DIR=$POLYBENCH_DIR"

    echo -e "# Experiment parameters"
    echo -e "NB_MEDIAN_RUNS=$NB_MEDIAN_RUNS"
    echo -e "NB_RUNS=$NB_RUNS"
    echo -e "STDDEV_TRESHOLD=$STDDEV_TRESHOLD"
    echo -e "TIMEOUT=$TIMEOUT"

    # OpenMP environment variables
    echo -e "# OpenMP environment:"
    echo -e "OMP_PLACES=$OMP_PLACES"
    echo -e "OMP_WAIT_POLICY=$OMP_WAIT_POLICY"
    echo -e "OMP_DYNAMIC=$OMP_DYNAMIC"
    echo -e "OMP_PROC_BIND=$OMP_PROC_BIND"
    echo -e "OMP_SCHEDULE=$OMP_SCHEDULE"

}
