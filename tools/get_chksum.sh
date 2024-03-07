#!/usr/bin/bash
# This script is used to get the checksum of a benchmark
# Usage: ./get_chksum.sh <target> [extra_flags]

# check arguments
if [ $# -lt 1 ]; then
    echo "Usage: ./get_chksum.sh <target> [extra_flags]"
    exit 1
fi
target=$1
EXTRA_FLAGS="$2 -DRESULT_DUMP"

# cleanup
make clean

# compile the benchmark with dump option
make "$target" PLUTO_EXTRA_FLAGS="-q" EXTRA_FLAGS="$EXTRA_FLAGS"

# run the benchmark
"./$target" 2>__temp_dump

# get the checksum
chksum=$(sha256sum __temp_dump | awk '{print $1}')

# remove the temporary files
rm -f __temp_dump

make clean

echo "$chksum"
