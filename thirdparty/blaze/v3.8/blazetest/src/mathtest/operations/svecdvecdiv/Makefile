#==================================================================================================
#
#  Makefile for the svecdvecdiv module of the Blaze test suite
#
#  Copyright (C) 2012-2020 Klaus Iglberger - All Rights Reserved - All Rights Reserved
#
#  This file is part of the Blaze library. You can redistribute it and/or modify it under
#  the terms of the New (Revised) BSD License. Redistribution and use in source and binary
#  forms, with or without modification, are permitted provided that the following conditions
#  are met:
#
#  1. Redistributions of source code must retain the above copyright notice, this list of
#     conditions and the following disclaimer.
#  2. Redistributions in binary form must reproduce the above copyright notice, this list
#     of conditions and the following disclaimer in the documentation and/or other materials
#     provided with the distribution.
#  3. Neither the names of the Blaze development group nor the names of its contributors
#     may be used to endorse or promote products derived from this software without specific
#     prior written permission.
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
#  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
#  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
#  SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
#  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
#  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
#  BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
#  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
#  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
#  DAMAGE.
#
#==================================================================================================


# Including the compiler and library settings
ifneq ($(MAKECMDGOALS),reset)
ifneq ($(MAKECMDGOALS),clean)
-include ../../../Makeconfig
endif
endif


# Setting the source, object and dependency files
SRC = $(wildcard ./*.cpp)
DEP = $(SRC:.cpp=.d)
OBJ = $(SRC:.cpp=.o)
BIN = $(SRC:.cpp=)


# General rules
default: VCaV3a VCaV3b VCbV3a VCbV3b VCaVHa VCaVHb VCbVHa VCbVHb \
         VCaVDa VCaVDb VCbVDa VCbVDb VCaVUa VCaVUb VCbVUa VCbVUb AliasingTest
all: $(BIN)
essential: VCaV3a VCaVHa VCaVDa VCaVUa AliasingTest
single: VCaVDa


# Build rules
VCaV2a: VCaV2a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCaV2b: VCaV2b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCbV2a: VCbV2a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCbV2b: VCbV2b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

VCaV3a: VCaV3a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCaV3b: VCaV3b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCbV3a: VCbV3a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCbV3b: VCbV3b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

VCaV4a: VCaV4a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCaV4b: VCaV4b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCbV4a: VCbV4a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCbV4b: VCbV4b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

VCaV5a: VCaV5a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCaV5b: VCaV5b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCbV5a: VCbV5a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCbV5b: VCbV5b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

VCaV6a: VCaV6a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCaV6b: VCaV6b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCbV6a: VCbV6a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCbV6b: VCbV6b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

VCaVDa: VCaVDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCaVDb: VCaVDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCbVDa: VCbVDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCbVDb: VCbVDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

VCaVHa: VCaVHa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCaVHb: VCaVHb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCbVHa: VCbVHa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCbVHb: VCbVHb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

VCaVUa: VCaVUa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCaVUb: VCaVUb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCbVUa: VCbVUa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCbVUb: VCbVUb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

AliasingTest: AliasingTest.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)


# Cleanup
reset:
	@$(RM) $(OBJ) $(BIN)
clean:
	@$(RM) $(OBJ) $(BIN) $(DEP)


# Makefile includes
ifneq ($(MAKECMDGOALS),reset)
ifneq ($(MAKECMDGOALS),clean)
-include $(DEP)
endif
endif


# Makefile generation
%.d: %.cpp
	@$(CXX) -MM -MP -MT "$*.o $*.d" -MF $@ $(CXXFLAGS) $<


# Setting the independent commands
.PHONY: default all essential single reset clean
