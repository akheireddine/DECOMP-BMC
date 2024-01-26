UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	CXX = 'xcrun g++-mp-10'
	CC = 'xcrun gcc-mp-10'
	INCLUDES = '-I/opt/local/include/'
	LIBS = '-L/opt/local/lib/'
	P_OPTS = '--with-gmp=/opt/local/include/'
endif

##################################################
###                PaInleSS                    ###
##################################################
all: build-minisat build-maplecomsps build-desat 
	+ cd painless-src 						&& \
	make release							&& \
	mv painless ../painless

##################################################
###                Minisat                     ###
##################################################
build-minisat:
	+ $(MAKE) -C minisat

##################################################
###                MapleCOMSPS                 ###
##################################################
build-maplecomsps:
	+ $(MAKE) -C mapleCOMSPS r

##################################################
###                   DeSAT                    ###
##################################################
build-desat:
	+ $(MAKE) -C desat						&& \
        cd desat &&  sh libdesat.sh					&& \
	cd ..

docs:
	rm -rf documents
	doxygen doxygen.config
	mkdir -p documents
	mv html latex documents

clean:
	+ $(MAKE) -C minisat clean
	+ $(MAKE) -C mapleCOMSPS clean
	+ $(MAKE) -C desat clean
	+ $(MAKE) clean -C painless-src
	rm -rf painless documents
