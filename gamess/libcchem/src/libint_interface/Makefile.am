#AM_CPPFLAGS = -I$(top_srcdir)/src @AM_CPPFLAGS@ -std=c++11
#AM_CPPFLAGS = -I$(top_srcdir)/src -I/opt/local/include/eigen3/ -I/Users/luke/devel/install/libint_gamess_1-6-16/include/ -I/Users/luke/devel/install/libint_gamess_1-6-16/include/libint2  @AM_CPPFLAGS@ -std=c++11

#AM_CPPFLAGS = -I$(top_srcdir)/src -I/opt/local/include/eigen3/ -I/Users/luke/devel/install/libint_gamess_12-10-15/include/ -I/Users/luke/devel/install/libint_gamess_12-10-15/include/libint2  @AM_CPPFLAGS@ -std=c++11

#AM_CPPFLAGS = -DHAVE_CONFIG_H -I. -I../../config -I../../src -I/Users/luke/devel/install/libint_gamess_12-1-16/include/ -I/Users/luke/devel/install/libint_gamess_12-1-16/include/libint2  -I /Users/luke/devel/build/ri-mp2-git/libcchem/src -I /Users/luke/devel/build/ri-mp2-git/libcchem/rysq/src -I /Users/luke/devel/build/ri-mp2-git/libcchem/src/externals -std=c++11 -DNDEBUG -I/opt/local/include/boost -O3 -msse3 -fopenmp -O3 -MT
include $(top_srcdir)/config/flymake.mk

AM_CPPFLAGS = -I$(top_srcdir)/src @LIBINT_FLAGS@ @AM_CPPFLAGS@ -std=c++11

noinst_LTLIBRARIES = liblibintinterface.la
liblibintinterface_la_SOURCES = libcchem_libint.cpp

