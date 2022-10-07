# This makefile will build the C++ components of Porechop.
# EDIT: qbonenfant
# Modified the makefile so it also build approx_counter tool

# Example commands:
#   make (build in release mode)
#   make debug (build in debug mode)
#   make clean (deletes *.o files, which aren't required to run the aligner)
#   make distclean (deletes *.o files and the *.so file, which is required to run the aligner)
#   make CXX=g++-5 (build with a particular compiler)
#   make CXXFLAGS="-Werror -g3" (build with particular compiler flags)


# CXX and CXXFLAGS can be overridden by the user.
CXX         ?= g++
CXXFLAGS    ?= -Wall -Wextra -pedantic -mtune=native

# These flags are required for the build to work.
FLAGS        = -std=c++17 -Iporechop_abi/include -fPIC
SHARED       = -shared  # trying to fix LDFLAG issues by splitting shared / not shared
LDFLAGS     ?= -Lporechop_abi/include

# Different debug/optimisation levels for debug/release builds.
DEBUGFLAGS   = -DSEQAN_ENABLE_DEBUG=1 -g
RELEASEFLAGS = -O3 -D NDEBUG

TARGET       = porechop_abi/cpp_functions.so
SHELL        = /bin/sh
# modified SOURCES / HEADER to only search cpp/h file in src folder
SOURCES      = $(shell find porechop_abi/src -name "*.cpp")
HEADERS      = $(shell find porechop_abi/src -name "*.h")
OBJECTS      = $(SOURCES:.cpp=.o)


# Adapt finder (k-mer approximate counter)
APPROXCOUNTER_TGT  = porechop_abi/approx_counter
APPROXCOUNTER_SRC  = porechop_abi/ab_initio_src/approx_counter.cpp
APPROXCOUNTER_FLAG = -DSEQAN_HAS_ZLIB=1

# Compatibility flag library
COMPAT_TGT = porechop_abi/compatibility.so
COMPAT_SRC = porechop_abi/ab_initio_src/compatibility.cpp
COMPAT_HDR = porechop_abi/ab_initio_src/compatibility.h
COMPAT_OBJ = $(COMPAT_SRC:.cpp=.o)

# Adding Compatibility library header to headers list
HEADERS+=$(COMPAT_HDR)


# MSA consensus
MSA_TGT = porechop_abi/msa_consensus
MSA_SRC = porechop_abi/ab_initio_src/msa_consensus.cpp

# Linux needs '-soname' while Mac needs '-install_name'
PLATFORM     = $(shell uname)
ifeq ($(PLATFORM), Darwin)
SONAME	= -install_name
OMP 	= 
LRT 	= 
LZ      = -lz
FLAGS   += -D_LIBCPP_DISABLE_AVAILABILITY

else
SONAME	= -soname
OMP 	= -fopenmp
LRT 	= -lrt
LZ      = -lz
endif


all: $(TARGET) $(COMPAT_TGT) $(MSA_TGT) $(APPROXCOUNTER_TGT)

.PHONY: release
release: FLAGS+=$(RELEASEFLAGS)
release: all

.PHONY: debug
debug: FLAGS+=$(DEBUGFLAGS)
debug: all


$(TARGET): $(OBJECTS)
	$(CXX) $(FLAGS) $(CXXFLAGS) $(SHARED) $(LDFLAGS) -Wl,$(SONAME),$(TARGET) -o $(TARGET) $(OBJECTS)

$(APPROXCOUNTER_TGT):
	$(CXX) $(FLAGS) $(APPROXCOUNTER_FLAG) $(OMP) $(CXXFLAGS) $(LDFLAGS) $(APPROXCOUNTER_SRC) $(LRT) $(LZ) -o $(APPROXCOUNTER_TGT)

$(COMPAT_TGT): $(COMPAT_OBJ)
	$(CXX) $(FLAGS) $(CXXFLAGS) $(SHARED) $(LDFLAGS) -Wl,$(SONAME),$(COMPAT_TGT) -o $(COMPAT_TGT) $(COMPAT_OBJ)

$(MSA_TGT):
	$(CXX) $(FLAGS) $(OMP) $(CXXFLAGS) $(LDFLAGS) $(MSA_SRC) $(LRT) -o $(MSA_TGT)
	
clean:
	$(RM) $(OBJECTS) $(COMPAT_OBJ)

distclean: clean
	$(RM) $(TARGET) $(APPROXCOUNTER_TGT) $(COMPAT_TGT) $(MSA_TGT)

%.o: %.cpp $(HEADERS)
	$(CXX) $(FLAGS) $(CXXFLAGS) -c -o $@ $<
