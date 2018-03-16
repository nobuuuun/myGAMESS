## LIBCCHEM BUILD INSTRUCTIONS

### Prequisites:
0.  GNU compilers v. 4.8.5 or higher, CUDA (optional), Intel Math Kernel Library (MKL)
1.  HDF5 - file system
    * Home Page : https://support.hdfgroup.org/
    * Download : https://support.hdfgroup.org/HDF5/release/obtain518.html
    * Instructions :
      * Download : https://support.hdfgroup.org/ftp/HDF5/current18/src/hdf5-1.8.18.tar.gz
      * Extract hdf5-1.8.18.tar.gz
      * ``cd hdf5-1.8.18``
      * ```
         ./configure \
         --prefix=$HOME/opt/hdf5 \
         --enable-shared \
         --enable-static \
         --enable-fortran \
         --enable-cxx
        ```
      * ``make -j6``
      * ``make -j6 install``
      * Note the full path to ``$HOME/opt/hdf5``. You will need this for the GAMESS config process.
2.  EIGEN - C++ template library for linear algebra:matrices, vectors, numerical solvers, and related algorithms
    * Home Page : http://eigen.tuxfamily.org/index.php?title=Main_Page
    * Download : http://bitbucket.org/eigen/eigen/get/3.3.3.tar.gz
    * Instructions :
      * Download http://bitbucket.org/eigen/eigen/get/3.3.3.tar.gz
      * Extract 3.3.3.tar.gz
      * ``mv eigen-eigen-* $HOME/opt/eigen``
      * Note the full path to ``$HOME/opt/eigen``. You will need this for the GAMESS config process.
3.  BOOST - C++ libraries
    * Home Page : http://www.boost.org/
    * Download : https://sourceforge.net/projects/boost/files/boost/1.63.0/boost_1_63_0.tar.gz
      * Instructions :
        * Download https://sourceforge.net/projects/boost/files/boost/1.63.0/boost_1_63_0.tar.gz
        * Extract boost_1_63_0.tar.gz
        * ``cd boost_1_63_0  ``
        * ``./bootstrap.sh --prefix=$HOME/opt/boost --with-toolset=gcc stage``
        * ``./b2 install --prefix=$HOME/opt/boost stage``
4.  GMP - GNU multiple precision arithmetic library
    * Home Page : https://gmplib.org/
    * Download : https://ftp.gnu.org/gnu/gmp/gmp-6.1.2.tar.xz
    * Instructions :
      * Download https://ftp.gnu.org/gnu/gmp/gmp-6.1.2.tar.xz
      * Extract gmp-6.1.2.tar.xz
      * ``cd gmp-6.1.2``
      * ``./configure --prefix=$HOME/opt/gmp --enable-cxx``
      * ``make -j6``
      * ``make -j6 install``
5.  LIBINT - C/C++ integral library (optional - needed for RI-MP2 or ZAPT gradients)
    * Home Page : https://github.com/evaleev/libint
    * Download : via GIT
    * Instructions :
      * ``export INCLUDEPATH=$HOME/opt/boost/include:$HOME/opt/gmp/include:$INCLUDEPATH``
      * ``export CPATH=$HOME/opt/boost/include:$HOME/opt/gmp/include:$CPATH``
      * ``export LD_RUN_PATH=$HOME/opt/boost/include:$HOME/opt/gmp/include:$LD_RUN_PATH``
      * ``export LD_LIBRARY_PATH=$HOME/opt/boost/include:$HOME/opt/gmp/include:$LD_LIBRARY_PATH``
      * ``git clone https://github.com/evaleev/libint.git libint``
      * ``cd libint``
      * ``git checkout 1f82c27a6b29a863f35f19f28180d3f05e91bb0a``
      * ``./autogen.sh``
      * ``mkdir build``
      * ``cd build``
      * ```
        ../configure \
        --prefix=$HOME/opt/libint \
        --enable-1body=1 \
        --enable-eri=1 \
        --with-max-am=4 \
        --with-cartgauss-ordering=gamess \
        --with-incdirs="-I$HOME/opt/boost/include/ -I$HOME/opt/gmp/include/" \
        --with-libdirs="-L$HOME/opt/boost/lib/ -L$HOME/opt/gmp/lib/" \
        --with-cxx-optflags=-O2 \
        CXXFLAGS=-O3 \
        CPPFLAGS=-I$HOME/opt/gmp
        ```
      * ``make -j6``
      * ``make -j6 install``
      * Comment out the first instance of ``renorm();`` in $HOME/opt/libint/include/libint2/shell.h
        *  Meaning change ``renorm();`` to ``//renorm();``
      * Note the full path to ``$HOME/opt/libint``. You will need this for the GAMESS config process.

### Run GAMESS ./config
