#!/bin/bash
# Abort on Error
set -e

BUILDING="EIGEN 3.3.3"
export PING_SLEEP=60s
export WORKDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export BUILD_OUTPUT=$WORKDIR/build.out

touch $BUILD_OUTPUT

dump_output() {
   echo $(date) - extracting $BUILDING ... completed
}

dump_output_error() {
   echo $BUILDING extraction output:
   cat $BUILD_OUTPUT
}

error_handler() {
  echo ERROR: An error was encountered with the extraction.
  dump_output_error
  exit 1
}
# If an error occurs, run our error handler to output a tail of the extraction
trap 'error_handler' ERR

# Set up a repeating loop to send some output to Travis.
bash -c "while true; do echo \$(date) - extracting $BUILDING ...; sleep $PING_SLEEP; done" &
PING_LOOP_PID=$!

cd /home/travis
wget --no-check-certificate http://bitbucket.org/eigen/eigen/get/3.3.3.tar.gz > $BUILD_OUTPUT 2>&1
tar -xf 3.3.3.tar.gz >> $BUILD_OUTPUT 2>&1
rm -rf 3.3.3.tar.gz

# The extraction finished without returning an error so dump the output
dump_output

# nicely terminate the ping output loop
kill $PING_LOOP_PID
