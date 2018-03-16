#!/bin/bash
# Abort on Error
set -e

BUILDING="MPICH 3.2"
export PING_SLEEP=120s
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

cd /home/travis
wget --no-check-certificate http://www.mpich.org/static/downloads/3.2/mpich-3.2.tar.gz > $BUILD_OUTPUT 2>&1
tar -xf mpich-3.2.tar.gz >> $BUILD_OUTPUT 2>&1
rm -rf mpich-3.2.tar.gz
cd mpich-3.2
mkdir build-mpich
cd build-mpich
../configure \
  --enable-fast=all,O2 \
  --enable-fast \
  --enable-shared \
  --prefix=/usr/local >> $BUILD_OUTPUT 2>&1
make -s -j $NUM_CPU_CORES >> $BUILD_OUTPUT 2>&1
sudo make -j $NUM_CPU_CORES install >> $BUILD_OUTPUT 2>&1

# The build finished without returning an error so dump a tail of the output
dump_output

# nicely terminate the ping output loop
kill $PING_LOOP_PID
