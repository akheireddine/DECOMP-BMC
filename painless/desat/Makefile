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
all: build-minisat
	+ $(MAKE) -C Desat release  		&& \
	mv Desat/desat .

##################################################
###                MiniSat                     ###
##################################################
build-minisat:
	+ $(MAKE) -C MiniSat r

clean:
	+ $(MAKE) -C MiniSat clean
	+ $(MAKE) clean -C Desat
	rm -rf desat
