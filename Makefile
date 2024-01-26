# Environnment variables.
CC = g++
INCLUDES=
LIBS = -L/opt/local/lib/

# This part is for supporting MacOS.
# Currently this is not very felxible since paths
# follows the standard macport installation.
# One should note that
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	CC = xcrun g++-mp-10
	INCLUDES = -I/opt/local/include/
	LIBS = -L/opt/local/lib/
endif

SRCS := $(shell find ./src -name "*.cc")

OBJS := $(addsuffix .o, $(basename $(SRCS)))

EXEC = decomp-bmc

CXXFLAGS = -Wall -Wextra  -ansi -std=c++17 -O3

LDFLAGS = -L./painless/lib \
          -L./painless/lib -lpainless \
          -L./NuSMV-2.6.0/NuSMV/lib -lnusmv \
          -lxml2 -lpthread -lz -lreadline \
          -lm -lgmpxx -lgmp -fopenmp $(LIBS)

CXXFLAGS += -I./painless/mapleCOMSPS \
            -I./painless/kissat \
            -I./painless \
            -I./src \
            -I./NuSMV-2.6.0/NuSMV/code \
            -I./NuSMV-2.6.0/NuSMV/build/build-cudd/include \
            -I./NuSMV-2.6.0/NuSMV/build \
             -DPRODUCE_PROOF $(INCLUDES)
            #  #-DFULL_LABELING

all: build-nusmv build-painless build
	echo Compilation done.

build: $(EXEC)

check:
	./tests/test.sh

$(EXEC): decomp-bmc.o $(OBJS)
	$(CC) -o $@ $^ $(CXXFLAGS) $(LDFLAGS)

%.o: %.cc
	$(CC) -o $@ -c $< $(CXXFLAGS) $(LDFLAGS)

build-nusmv:
	cd NuSMV-2.6.0/NuSMV && mkdir -p build    		 	&& \
	cd build                                                 	&& \
	cmake .. -DPYTHON_EXECUTABLE:FILEPATH=/usr/bin/python2.7 	&& \
	make                                                     	&& \
	cd .. 							 	&& \
	mkdir -p lib						 	&& \
	$(CC) -shared -o lib/libnusmv.so -lxml2 -lpthread    	           \
	-lz -lreadline  -lm -lgmpxx $(LIBS) -lgmp               	   \
	`find ./build/code/nusmv/shell/cmd/ ./build/code/nusmv/addons_core -type f -wholename "*.o"`   `find ./build/code/nusmv/core   -type f -wholename "*.o"` ./build/build-cudd/lib/*.a `find ./build/build-MiniSat/ -type f -wholename "*.a"`


build-painless:
	cd painless                                               && \
	make                                                      && \
	sh libpainless.sh

clean:
	rm -f $(OBJS) decomp-bmc.o $(EXEC)

cleanall: clean
	cd  painless && make clean && cd .. && rm -rf NuSMV-2.6.0/NuSMV/build && rm -Rf NuSMV-2.6.0/NuSMV/lib  

mrproper: clean
	rm -rf $(EXEC) decomp-bmc.o

.PHONY: clean mrproper
