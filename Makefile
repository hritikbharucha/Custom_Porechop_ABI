# This makefile will build the C++ components of Porechop.
# EDIT: qbonenfant
# Modified the makefile so it also build adaptFinder tool

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
FLAGS        = -std=c++14 -Iporechop/include -fPIC
LDFLAGS      = -shared

# Different debug/optimisation levels for debug/release builds.
DEBUGFLAGS   = -DSEQAN_ENABLE_DEBUG=1 -g
RELEASEFLAGS = -O3 -D NDEBUG

TARGET       = porechop/cpp_functions.so
SHELL        = /bin/sh
# modified SOURCES / HEADER to only search cpp/h file in src folder
SOURCES      = $(shell find porechop/src -name "*.cpp")
HEADERS      = $(shell find porechop/src -name "*.h")
OBJECTS      = $(SOURCES:.cpp=.o)


# Adapt finder (k-mer approximate counter)
ADAPTFINDER_SRC = adaptFinder/adaptFinder.cpp
ADAPTFINDER_TGT = porechop/adaptFinder

# Compatibility flag library
COMPAT_TGT = porechop/compatibility.so
COMPAT_SRC = porechop/ab_initio_src/compatibility.cpp
COMPAT_HDR = porechop/ab_initio_src/compatibility.h
COMPAT_OBJ = $(COMPAT_SRC:.cpp=.o)

# Adding Compatibility library header to headers list
HEADERS+=$(COMPAT_HDR)


# MSA consensus
MSA_TGT = porechop/msa_consensus
MSA_SRC = porechop/ab_initio_src/msa_consensus.cpp

# Linux needs '-soname' while Mac needs '-install_name'
PLATFORM     = $(shell uname)
ifeq ($(PLATFORM), Darwin)
SONAME	= -install_name
OMP 	= 
LRT 	= 

else
SONAME	= -soname
OMP 	= -fopenmp
LRT 	= -lrt
endif


all: $(TARGET) $(COMPAT_TGT) $(MSA_TGT) $(ADAPTFINDER_TGT)

.PHONY: release
release: FLAGS+=$(RELEASEFLAGS)
release: all

.PHONY: debug
debug: FLAGS+=$(DEBUGFLAGS)
debug: all


$(TARGET): $(OBJECTS)
	$(CXX) $(FLAGS) $(CXXFLAGS) $(LDFLAGS) -Wl,$(SONAME),$(TARGET) -o $(TARGET) $(OBJECTS)

$(ADAPTFINDER_TGT):
	$(CXX) $(FLAGS) $(OMP) $(CXXFLAGS)  $(ADAPTFINDER_SRC) $(LRT) -o $(ADAPTFINDER_TGT)

$(COMPAT_TGT): $(COMPAT_OBJ)
	$(CXX) $(FLAGS) $(CXXFLAGS) $(LDFLAGS) -Wl,$(SONAME),$(COMPAT_TGT) -o $(COMPAT_TGT) $(COMPAT_OBJ)

$(MSA_TGT):
	$(CXX) $(FLAGS) $(OMP) $(CXXFLAGS) $(MSA_SRC) $(LRT) -o $(MSA_TGT)
	
clean:
	$(RM) $(OBJECTS) $(COMPAT_OBJ)

distclean: clean
	$(RM) $(TARGET) $(ADAPTFINDER_TGT) $(COMPAT_TGT) $(MSA_TGT)

%.o: %.cpp $(HEADERS)
	$(CXX) $(FLAGS) $(CXXFLAGS) -c -o $@ $<
