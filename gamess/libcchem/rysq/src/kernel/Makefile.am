include $(top_srcdir)/config/flymake.am
include $(top_srcdir)/config/cheetah.mk

BUILT_SOURCES = quadrature1-impl.hpp quadrature2-impl.hpp
BUILT_SOURCES += derivative-overlap-quadrature-impl.hpp
BUILT_SOURCES += derivative-twoc-quadrature-impl.hpp
BUILT_SOURCES += derivative-threec-quadrature-impl.hpp
BUILT_SOURCES += twoc-quadrature-impl.hpp

AM_CPPFLAGS = -I$(top_srcdir)/src @AM_CPPFLAGS@

quadrature2.py: $(top_srcdir)/lib/python/rysq-generate.py
	$(top_srcdir)/lib/python/rysq-generate.py > $@

quadrature2-impl.hpp: quadrature2.py
