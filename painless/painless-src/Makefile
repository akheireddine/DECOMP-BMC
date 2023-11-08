SRCS = $(shell find . -name "*.cpp")
SRCS += $(shell find . -name "*.cc")

OBJS = $(addsuffix .o, $(basename $(SRCS)))

EXEC = painless

## Detect OS
UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Darwin)
	INCLUDES=-I/opt/local/include/
endif

CXX = g++
ifeq ($(UNAME_S),Darwin)
	CXX = xcrun g++-mp-10
endif

LIBS = -lminisat -L../minisat/build/release/lib/                 \
       -lglucose -L../glucose/parallel                           \
       -llgl -L../lingeling/                                     \
       -lmapleCOMSPS -L../mapleCOMSPS/build/release/lib/         \
       -ldesat -L../desat/lib/					 \
       -lpthread -lz -lm -lgmpxx -L/opt/local/lib/ -lgmp -fopenmp

CXXFLAGS = -I../minisat -I../glucose -I../lingeling              \
           -I../mapleCOMSPS          			           \
           -I../desat/Desat -I../desat/Desat/libExpression	 \
           $(INCLUDES) -DPRODUCE_PROOF                           \
           -D __STDC_LIMIT_MACROS -D __STDC_FORMAT_MACROS        \
           -std=c++11 -O3


release: CXXFLAGS += -D NDEBUG -Wno-deprecated
release: $(EXEC)

debug: CXXFLAGS += -DEBUG -g -Wno-deprecated -ggdb
debug: $(EXEC)

$(EXEC): $(OBJS)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LIBS)

%.o: %.cpp
	$(CXX) -c $< -o $@ $(CXXFLAGS) $(LIBS)

clean:
	rm -f $(OBJS) $(EXEC)
