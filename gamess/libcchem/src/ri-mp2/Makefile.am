AM_CPPFLAGS = -I$(top_srcdir)/src @AM_CPPFLAGS@  @LIBINT_FLAGS@

include $(top_srcdir)/config/flymake.mk

noinst_LTLIBRARIES = librimp2.la
librimp2_la_SOURCES = metric.cpp
librimp2_la_SOURCES += ri-energy.cpp
librimp2_la_SOURCES += three-center-work.cpp
librimp2_la_SOURCES += ri-zvector.cpp
librimp2_la_SOURCES += ri-lagrangian.cpp
librimp2_la_SOURCES += ri-gradient.cpp
librimp2_la_SOURCES += JK.cpp
librimp2_la_SOURCES += gradient-helpers.cpp
librimp2_la_SOURCES += ri-async.cpp
librimp2_la_SOURCES += DFIntGradient.cpp
librimp2_la_SOURCES += deviceGradientEngine.cpp 

if HAVE_CUBLAS
SUBDIRS = kernel
librimp2_la_LIBADD = kernel/librimp2kernel.la
endif
