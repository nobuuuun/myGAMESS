#!/bin/bash
# Abort on Error
set -e

BUILDING="EFMO w64_r0.3 RUN"
export PING_SLEEP=60s
export WORKDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export BUILD_OUTPUT=$WORKDIR/build.out

touch $BUILD_OUTPUT

dump_output() {
   echo $(date) - building $BUILDING ... completed
   cat $BUILD_OUTPUT
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

tests/runtest.py --folder=efmo --file=w64_r0.3 --ncpus=1 > $BUILD_OUTPUT 2>&1
cd tests
./checkgms.py --verbose_parsing --verbose_validation --exit_on_fail >> $BUILD_OUTPUT 2>&1
cd ../
rm tests/efmo/w64_r0.3.*log
tests/runtest.py --folder=efmo --file=w64_r0.3 --ncpus=2 >> $BUILD_OUTPUT 2>&1
cd tests
./checkgms.py --verbose_parsing --verbose_validation --exit_on_fail >> $BUILD_OUTPUT 2>&1
cd ../
rm tests/efmo/w64_r0.3.*log

# The build finished without returning an error so dump a tail of the output
dump_output

# nicely terminate the ping output loop
kill $PING_LOOP_PID

