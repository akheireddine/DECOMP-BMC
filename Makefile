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

EXEC = bmc-decomp

CXXFLAGS = -Wall -Wextra  -ansi -std=c++17 -O3

LDFLAGS = -L./painless/lib -l_standard \
          -L./painless/lib -lpainless \
          -lxml2 -lpthread -lz -lreadline \
          -lm -lgmpxx -lgmp -fopenmp $(LIBS)

CXXFLAGS += -I./painless/mapleCOMSPS \
            -I./painless/kissat \
            -I./painless \
            -I./src \
             -DPRODUCE_PROOF $(INCLUDES)
            #  #-DFULL_LABELING

all: build-painless build
	echo Compilation done.

build: $(EXEC)

check:
	./tests/test.sh

$(EXEC): bmcDecompMain.o $(OBJS)
	$(CC) -o $@ $^ $(CXXFLAGS) $(LDFLAGS)

%.o: %.cc
	$(CC) -o $@ -c $< $(CXXFLAGS) $(LDFLAGS)

build-painless:
	cd painless                                               && \
	make                                                      && \
	sh libpainless.sh

clean:
	rm -f $(OBJS) bmcDecompMain.o $(EXEC)

cleanall: clean
	cd  painless && make clean && cd .. 

mrproper: clean
	rm -rf $(EXEC) bmcDecompMain.o

.PHONY: clean mrproper
