#==================================================================================================
#
#  Makefile for the dvecdvecouter module of the Blaze test suite
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
V2aV2a: V2aV2a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V2aV2b: V2aV2b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V2aV3a: V2aV3a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V2aV3b: V2aV3b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V2aV4a: V2aV4a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V2aV4b: V2aV4b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V2aV5a: V2aV5a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V2aV5b: V2aV5b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V2aV6a: V2aV6a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V2aV6b: V2aV6b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V2aVDa: V2aVDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V2aVDb: V2aVDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V2bV2a: V2bV2a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V2bV2b: V2bV2b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V2bV3a: V2bV3a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V2bV3b: V2bV3b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V2bV4a: V2bV4a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V2bV4b: V2bV4b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V2bV5a: V2bV5a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V2bV5b: V2bV5b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V2bV6a: V2bV6a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V2bV6b: V2bV6b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V2bVDa: V2bVDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V2bVDb: V2bVDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

V3aV2a: V3aV2a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V3aV2b: V3aV2b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V3aV3a: V3aV3a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V3aV3b: V3aV3b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V3aV4a: V3aV4a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V3aV4b: V3aV4b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V3aV5a: V3aV5a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V3aV5b: V3aV5b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V3aV6a: V3aV6a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V3aV6b: V3aV6b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V3aVDa: V3aVDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V3aVDb: V3aVDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V3bV2a: V3bV2a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V3bV2b: V3bV2b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V3bV3a: V3bV3a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V3bV3b: V3bV3b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V3bV4a: V3bV4a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V3bV4b: V3bV4b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V3bV5a: V3bV5a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V3bV5b: V3bV5b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V3bV6a: V3bV6a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V3bV6b: V3bV6b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V3bVDa: V3bVDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V3bVDb: V3bVDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

V4aV2a: V4aV2a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V4aV2b: V4aV2b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V4aV3a: V4aV3a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V4aV3b: V4aV3b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V4aV4a: V4aV4a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V4aV4b: V4aV4b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V4aV5a: V4aV5a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V4aV5b: V4aV5b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V4aV6a: V4aV6a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V4aV6b: V4aV6b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V4aVDa: V4aVDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V4aVDb: V4aVDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V4bV2a: V4bV2a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V4bV2b: V4bV2b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V4bV3a: V4bV3a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V4bV3b: V4bV3b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V4bV4a: V4bV4a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V4bV4b: V4bV4b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V4bV5a: V4bV5a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V4bV5b: V4bV5b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V4bV6a: V4bV6a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V4bV6b: V4bV6b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V4bVDa: V4bVDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V4bVDb: V4bVDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

V5aV2a: V5aV2a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V5aV2b: V5aV2b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V5aV3a: V5aV3a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V5aV3b: V5aV3b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V5aV4a: V5aV4a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V5aV4b: V5aV4b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V5aV5a: V5aV5a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V5aV5b: V5aV5b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V5aV6a: V5aV6a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V5aV6b: V5aV6b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V5aVDa: V5aVDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V5aVDb: V5aVDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V5bV2a: V5bV2a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V5bV2b: V5bV2b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V5bV3a: V5bV3a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V5bV3b: V5bV3b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V5bV4a: V5bV4a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V5bV4b: V5bV4b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V5bV5a: V5bV5a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V5bV5b: V5bV5b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V5bV6a: V5bV6a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V5bV6b: V5bV6b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V5bVDa: V5bVDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V5bVDb: V5bVDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

V6aV2a: V6aV2a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V6aV2b: V6aV2b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V6aV3a: V6aV3a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V6aV3b: V6aV3b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V6aV4a: V6aV4a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V6aV4b: V6aV4b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V6aV5a: V6aV5a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V6aV5b: V6aV5b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V6aV6a: V6aV6a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V6aV6b: V6aV6b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V6aVDa: V6aVDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V6aVDb: V6aVDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V6aVHa: V6aVHa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V6aVHb: V6aVHb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V6aVUa: V6aVUa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V6aVUb: V6aVUb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V6bV2a: V6bV2a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V6bV2b: V6bV2b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V6bV3a: V6bV3a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V6bV3b: V6bV3b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V6bV4a: V6bV4a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V6bV4b: V6bV4b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V6bV5a: V6bV5a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V6bV5b: V6bV5b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V6bV6a: V6bV6a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V6bV6b: V6bV6b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V6bVDa: V6bVDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V6bVDb: V6bVDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V6bVHa: V6bVHa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V6bVHb: V6bVHb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V6bVUa: V6bVUa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V6bVUb: V6bVUb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

VDaV2a: VDaV2a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaV2b: VDaV2b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaV3a: VDaV3a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaV3b: VDaV3b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaV4a: VDaV4a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaV4b: VDaV4b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaV5a: VDaV5a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaV5b: VDaV5b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaV6a: VDaV6a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaV6b: VDaV6b.o
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
VDbV2a: VDbV2a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbV2b: VDbV2b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbV3a: VDbV3a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbV3b: VDbV3b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbV4a: VDbV4a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbV4b: VDbV4b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbV5a: VDbV5a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbV5b: VDbV5b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbV6a: VDbV6a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbV6b: VDbV6b.o
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

VHaV6a: VHaV6a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHaV6b: VHaV6b.o
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
VHbV6a: VHbV6a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHbV6b: VHbV6b.o
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

VUaV6a: VUaV6a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUaV6b: VUaV6b.o
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
VUbV6a: VUbV6a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUbV6b: VUbV6b.o
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
