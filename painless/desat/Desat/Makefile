CC = gcc
SRCS  = $(shell find libExpression/ -name "*.cpp")
SRCS += $(shell find libRBC/ -name "*.cpp")
SRCS += $(shell find libDeSAT/ -name "*.cpp")
SRCS += $(shell find DeSAT/ -name "*.cpp")

OBJS = $(addsuffix .o, $(basename $(SRCS)))
OBJ += ../cadical/build/*.o

EXEC = desat

## Detect OS
UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Darwin)
	INCLUDES=-I/opt/local/include/
endif

CXX = g++
ifeq ($(UNAME_S),Darwin)
	CXX = xcrun g++-mp-10
endif

LIBS += -lpthread -lz -lm -lgmpxx -L/opt/local/lib/ -lgmp -fopenmp \
	-L../cadical/build 					   \
	-L../cadical/build -lcadical

CXXFLAGS =  -I../cadical/ -I../cadical/src 			\
	   -IlibDeSAT -IlibExpression -IlibRBC  -Ilibcadical	\
           $(INCLUDES) 						\
           -D __STDC_LIMIT_MACROS -D __STDC_FORMAT_MACROS       \
           -std=c++17 -O3

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
