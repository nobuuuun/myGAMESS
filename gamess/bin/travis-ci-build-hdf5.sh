#!/bin/bash
# Abort on Error
set -e

BUILDING="HDF5 1.10.1"
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
tar -xf hdf5-1.10.1.tar.gz >> $BUILD_OUTPUT 2>&1
rm -rf hdf5-1.10.1.tar.gz
cd hdf5-1.10.1
mkdir build-hdf5
cd build-hdf5
../configure \
  --enable-shared \
  --enable-static \
  --enable-fortran \
  --enable-cxx \
  --with-szlib=/opt/szip \
  --prefix=/usr/local >> $BUILD_OUTPUT 2>&1
make -s -j $NUM_CPU_CORES >> $BUILD_OUTPUT 2>&1
sudo make -j $NUM_CPU_CORES install >> $BUILD_OUTPUT 2>&1

# The build finished without returning an error so dump a tail of the output
dump_output

# nicely terminate the ping output loop
kill $PING_LOOP_PID
