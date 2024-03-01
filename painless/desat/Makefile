UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	CXX = 'xcrun g++-mp-10'
	CC = 'xcrun gcc-mp-10'
	INCLUDES = '-I/opt/local/include/'
	LIBS = '-L/opt/local/lib/'
	P_OPTS = '--with-gmp=/opt/local/include/'
endif

##################################################
###                 DeSAT                      ###
##################################################
all:  build-cadical
	+ $(MAKE) -C Desat release  		&& \
	mv Desat/desat .



##################################################
###              	  CaDiCaL                     ###
##################################################
build-cadical:
	cd cadical	                        &&\
	./configure                             &&\
	make

clean:
	+ $(MAKE) clean -C Desat
	+ $(MAKE) clean -C cadical
	rm -rf desat
