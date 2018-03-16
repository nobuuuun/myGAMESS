#!/bin/bash
# Abort on Error
ncpus=$1
folder=$2
file=$3

set -e

BUILDING="TEST SET LIBCCHEM VALIDATION : $folder"
export PING_SLEEP=30s
export WORKDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export BUILD_OUTPUT=$WORKDIR/build.out

touch $BUILD_OUTPUT

dump_output() {
   echo $(date) - building $BUILDING ... completed
}

dump_output_error() {
   echo Lines of $BUILDING build output:
   cat $BUILD_OUTPUT
}

error_handler() {
  echo ERROR: An error was encountered with the build.
  dump_output_error
  exit 1
}
# If an error occurs, run our error handler to output a tail of the build
trap 'error_handler' ERR

# Set up a repeating loop to send some output to Travis.
bash -c "while true; do echo \$(date) - building $BUILDING ...; sleep $PING_SLEEP; done" &
PING_LOOP_PID=$!

tests/runtest.py --ncpus=$ncpus --folder=$folder --file=$file

# The build finished without returning an error so dump a tail of the output
dump_output

# nicely terminate the ping output loop
kill $PING_LOOP_PID
