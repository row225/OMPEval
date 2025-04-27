CXXFLAGS += -O3 -std=c++11 -Wall -Wpedantic

# TBB configuration (using local version)
TBB_ROOT = oneapi-tbb-2022.0.0
TBB_INCLUDE = -I$(TBB_ROOT)/include
TBB_LIB = -L$(TBB_ROOT)/lib/intel64/gcc4.8
TBB_LIBS = $(TBB_LIB) -ltbb -ltbbmalloc

ifdef SYSTEMROOT
    CXXFLAGS += -lpthread
else
    CXXFLAGS += -pthread
endif

ifeq ($(SSE4),1)
    CXXFLAGS += -msse4.2
endif

SRCS := $(wildcard omp/*.cpp)
OBJS := ${SRCS:.cpp=.o}

all: lib/ompeval.a test main

lib:
	mkdir -p lib

lib/ompeval.a: $(OBJS) | lib
	ar rcs $@ $^

test: test.cpp benchmark.cpp lib/ompeval.a
	$(CXX) $(CXXFLAGS) $(TBB_INCLUDE) -o $@ $^ $(TBB_LIBS) -Wl,-rpath,$(TBB_ROOT)/lib/intel64/gcc4.8

main: main.cpp lib/ompeval.a
	$(CXX) $(CXXFLAGS) $(TBB_INCLUDE) -o $@ $^ $(TBB_LIBS) -Wl,-rpath,$(TBB_ROOT)/lib/intel64/gcc4.8

clean:
	$(RM) test test.exe main main.exe lib/ompeval.a $(OBJS)
