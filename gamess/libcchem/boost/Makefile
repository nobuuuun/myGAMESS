ifndef PREFIX
    PREFIX := $(PWD)
endif

VERSION=1.46.1
BOOST=$(shell echo "boost_$(VERSION)" | sed "s/\./_/g")
LIBRARIES := date_time,thread #,math

TOOLSET=

ifeq ($(TOOLSET),)
TOOLSET := $(shell $(CXX) -v 2>&1 | grep -q "gcc version" && echo "gcc")
endif

ifeq ($(TOOLSET),)
TOOLSET := $(shell $(CXX) --version | grep -q "icc" && echo "intel")
endif

ifeq ($(TOOLSET),)
TOOLSET := $(shell basename $(CXX))
endif

boost: $(BOOST).tar.gz
	tar -xzf $(BOOST).tar.gz
	cd $(BOOST) && patch -p0 < ../$(BOOST).nvcc.patch
	#Uncomment the two lines below if you run into issues with TIME_UTC
	#sed -i 's/TIME_UTC/TIME_UTC_/g' boost_1_46_1/boost/thread/xtime.hpp
	#sed -i 's/TIME_UTC/TIME_UTC_/g' boost_1_46_1/libs/thread/src/pthread/timeconv.inl
	cd $(BOOST) && patch -p0 < ../boost-gcc-4.7.patch
	echo "using $(TOOLSET) : : $(CXX) ;" > $(BOOST)/tools/build/v2/user-config.jam 
	cd $(BOOST) && \
          sh bootstrap.sh --with-libraries=$(LIBRARIES) && \
          ./bjam --toolset=$(TOOLSET) --prefix=$(PREFIX) install
	rm -fr $(BOOST)

clean:
	rm -fr $(BOOST) lib include

$(BOOST).tar.gz:
	wget http://downloads.sourceforge.net/project/boost/boost/$(VERSION)/$(BOOST).tar.gz
