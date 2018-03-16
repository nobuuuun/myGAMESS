include $(top_srcdir)/config/flymake.mk
AM_CPPFLAGS =  -I$(top_srcdir)/src -I$(CUDA_ROOT)/include @AM_CPPFLAGS@

AM_CUDAFLAGS = -rdc=true

noinst_LTLIBRARIES = librimp2kernel.la

librimp2kernel_o_OBJECTS = ri-energy-kernel.o ri-gradient-kernel.o rimp2kernel.o
rimp2kernel_o_OBJECTS = ri-gradient-kernel.o ri-energy-kernel.o
librimp2kernel_la_DEPENDENCIES = rimp2kernel.o

librimp2kernel_la_SOURCES =
librimp2kernel_la_SOURCES += ri-energy-kernel.cu
librimp2kernel_la_SOURCES += ri-gradient-kernel.cu

librimp2kernel_la_LIBADD = rimp2kernel.lo


rimp2kernel.o: $(rimp2kernel_o_OBJECTS)
	$(LIBTOOL) --tag=CXX --mode=compile $(CUDA) $(CUDAFLAGS) -dlink -lib -lcublas_device -lcudadevrt -o rimp2kernel.o $(rimp2kernel_o_OBJECTS)




librimp2kernel.la: $(librimp2kernel_la_OBJECTS) $(librimp2kernel_la_DEPENDENCIES)
	$(CUDALINK)  $(librimp2kernel_la_OBJECTS) $(librimp2kernel_la_LIBADD) $(LIBS)


