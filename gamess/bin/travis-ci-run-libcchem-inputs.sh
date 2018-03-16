#!/bin/bash
# Abort on Error
set -e

BUILDING="LIBCCHEM TEST VALIDATION"
export PING_SLEEP=60s
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

tests/runtest.py --folder=libcchem/ --ncpus=1 > $BUILD_OUTPUT 2>&1

# The build finished without returning an error so dump a tail of the output
dump_output

# nicely terminate the ping output loop
kill $PING_LOOP_PID
