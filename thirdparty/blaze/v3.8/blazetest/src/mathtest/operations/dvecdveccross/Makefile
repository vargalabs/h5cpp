#==================================================================================================
#
#  Makefile for the dvecdveccross module of the Blaze test suite
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
default: V3aV3a V3aV3b V3bV3a V3bV3b VHaVHa VHaVHb VHbVHa VHbVHb \
         VDaVDa VDaVDb VDbVDa VDbVDb VUaVUa VUaVUb VUbVUa VUbVUb AliasingTest
all: $(BIN)
essential: V3aV3a VHaVHa VDaVDa VUaVUa AliasingTest
single: VDaVDa


# Build rules
V3aV3a: V3aV3a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V3aV3b: V3aV3b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V3aVDa: V3aVDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V3aVDb: V3aVDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V3aVHa: V3aVHa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V3aVHb: V3aVHb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V3aVUa: V3aVUa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V3aVUb: V3aVUb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V3bV3a: V3bV3a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V3bV3b: V3bV3b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V3bVDa: V3bVDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V3bVDb: V3bVDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V3bVHa: V3bVHa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V3bVHb: V3bVHb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V3bVUa: V3bVUa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V3bVUb: V3bVUb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

VDaV3a: VDaV3a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaV3b: VDaV3b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaVDa: VDaVDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaVDb: VDaVDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaVHa: VDaVHa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaVHb: VDaVHb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaVUa: VDaVUa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaVUb: VDaVUb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbV3a: VDbV3a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbV3b: VDbV3b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbVDa: VDbVDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbVDb: VDbVDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbVHa: VDbVHa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbVHb: VDbVHb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbVUa: VDbVUa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbVUb: VDbVUb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

VHaV3a: VHaV3a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHaV3b: VHaV3b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHaVDa: VHaVDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHaVDb: VHaVDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHaVHa: VHaVHa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHaVHb: VHaVHb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHaVUa: VHaVUa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHaVUb: VHaVUb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHbV3a: VHbV3a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHbV3b: VHbV3b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHbVDa: VHbVDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHbVDb: VHbVDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHbVHa: VHbVHa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHbVHb: VHbVHb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHbVUa: VHbVUa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHbVUb: VHbVUb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

VUaV3a: VUaV3a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUaV3b: VUaV3b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUaVDa: VUaVDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUaVDb: VUaVDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUaVHa: VUaVHa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUaVHb: VUaVHb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUaVUa: VUaVUa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUaVUb: VUaVUb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUbV3a: VUbV3a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUbV3b: VUbV3b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUbVDa: VUbVDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUbVDb: VUbVDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUbVHa: VUbVHa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUbVHb: VUbVHb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUbVUa: VUbVUa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUbVUb: VUbVUb.o
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
