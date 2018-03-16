#!/bin/bash
# Abort on Error
set -e

BUILDING="LIBINT"
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
git clone https://github.com/evaleev/libint.git libint
cd libint
git checkout 1f82c27a6b29a863f35f19f28180d3f05e91bb0a
./autogen.sh
cd /home/travis
mkdir build-libbint
cd build-libbint
../libint/configure \
  --prefix=/opt/libint \
  --enable-1body=1 \
  --enable-eri=1 \
  --with-max-am=4 \
  --with-cartgauss-ordering=gamess \
  --with-cxx-optflags=-O2 \
  CXXFLAGS=-O3 >> $BUILD_OUTPUT 2>&1
make -j $NUM_CPU_CORES >> $BUILD_OUTPUT 2>&1
sudo make -j $NUM_CPU_CORES install >> $BUILD_OUTPUT 2>&1

# The build finished without returning an error so dump a tail of the output
dump_output

# nicely terminate the ping output loop
kill $PING_LOOP_PID
