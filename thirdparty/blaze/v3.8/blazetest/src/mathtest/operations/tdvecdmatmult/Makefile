#==================================================================================================
#
#  Makefile for the tdvecdmatmult module of the Blaze test suite
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
default: V3aM3x3a V3aM3x3b V3bM3x3a V3bM3x3b \
         VHaMHa VHaMHb VHbMHa VHbMHb \
         VDaMDa VDaMDb VDbMDa VDbMDb \
         VUaMUa VUaMUb VUbMUa VUbMUb \
         VDaSDa VDaSDb VDbSDa VDbSDb \
         VDaHDa VDaHDb VDbHDa VDbHDb \
         VDaLDa VDaLDb VDbLDa VDbLDb \
         VDaUDa VDaUDb VDbUDa VDbUDb \
         VDaDDa VDaDDb VDbDDa VDbDDb \
         AliasingTest
all: $(BIN)
essential: V3aM3x3a VHaMHa VDaMDa VUaMUa VDaSDa VDaHDa VDaLDa VDaUDa VDaDDa AliasingTest
single: VDaMDa


# Build rules
V2aM2x2a: V2aM2x2a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V2aM2x2b: V2aM2x2b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V2aMDa: V2aMDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V2aMDb: V2aMDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V2bM2x2a: V2bM2x2a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V2bM2x2b: V2bM2x2b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V2bMDa: V2bMDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V2bMDb: V2bMDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

V3aM3x3a: V3aM3x3a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V3aM3x3b: V3aM3x3b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V3aMDa: V3aMDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V3aMDb: V3aMDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V3bM3x3a: V3bM3x3a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V3bM3x3b: V3bM3x3b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V3bMDa: V3bMDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V3bMDb: V3bMDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

V4aM4x4a: V4aM4x4a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V4aM4x4b: V4aM4x4b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V4aMDa: V4aMDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V4aMDb: V4aMDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V4bM4x4a: V4bM4x4a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V4bM4x4b: V4bM4x4b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V4bMDa: V4bMDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V4bMDb: V4bMDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

V5aM5x5a: V5aM5x5a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V5aM5x5b: V5aM5x5b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V5aMDa: V5aMDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V5aMDb: V5aMDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V5bM5x5a: V5bM5x5a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V5bM5x5b: V5bM5x5b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V5bMDa: V5bMDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V5bMDb: V5bMDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

V6aM6x6a: V6aM6x6a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V6aM6x6b: V6aM6x6b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V6aMDa: V6aMDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V6aMDb: V6aMDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V6aMHa: V6aMHa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V6aMHb: V6aMHb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V6aMUa: V6aMUa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V6aMUb: V6aMUb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V6bM6x6a: V6bM6x6a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V6bM6x6b: V6bM6x6b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V6bMDa: V6bMDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V6bMDb: V6bMDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V6bMHa: V6bMHa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V6bMHb: V6bMHb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V6bMUa: V6bMUa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V6bMUb: V6bMUb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

V7aM7x13a: V7aM7x13a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V7aM7x13b: V7aM7x13b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V7bM7x13a: V7bM7x13a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V7bM7x13b: V7bM7x13b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

V16aM16x8a: V16aM16x8a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V16aM16x8b: V16aM16x8b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V16bM16x8a: V16bM16x8a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V16bM16x8b: V16bM16x8b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

VDaDDa: VDaDDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaDDb: VDaDDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaHDa: VDaHDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaHDb: VDaHDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaLDa: VDaLDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaLDb: VDaLDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaM2x2a: VDaM2x2a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaM2x2b: VDaM2x2b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaM3x3a: VDaM3x3a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaM3x3b: VDaM3x3b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaM4x4a: VDaM4x4a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaM4x4b: VDaM4x4b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaM5x5a: VDaM5x5a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaM5x5b: VDaM5x5b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaM6x6a: VDaM6x6a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaM6x6b: VDaM6x6b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaM7x13a: VDaM7x13a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaM7x13b: VDaM7x13b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaM16x8a: VDaM16x8a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaM16x8b: VDaM16x8b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaMDa: VDaMDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaMDb: VDaMDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaMHa: VDaMHa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaMHb: VDaMHb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaMUa: VDaMUa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaMUb: VDaMUb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaSDa: VDaSDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaSDb: VDaSDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaSLDa: VDaSLDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaSLDb: VDaSLDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaSUDa: VDaSUDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaSUDb: VDaSUDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaUDa: VDaUDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaUDb: VDaUDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbDDa: VDbDDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbDDb: VDbDDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbHDa: VDbHDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbHDb: VDbHDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbLDa: VDbLDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbLDb: VDbLDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbM2x2a: VDbM2x2a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbM2x2b: VDbM2x2b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbM3x3a: VDbM3x3a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbM3x3b: VDbM3x3b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbM4x4a: VDbM4x4a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbM4x4b: VDbM4x4b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbM5x5a: VDbM5x5a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbM5x5b: VDbM5x5b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbM6x6a: VDbM6x6a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbM6x6b: VDbM6x6b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbM7x13a: VDbM7x13a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbM7x13b: VDbM7x13b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbM16x8a: VDbM16x8a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbM16x8b: VDbM16x8b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbMDa: VDbMDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbMDb: VDbMDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbMHa: VDbMHa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbMHb: VDbMHb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbMUa: VDbMUa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbMUb: VDbMUb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbSDa: VDbSDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbSDb: VDbSDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbSLDa: VDbSLDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbSLDb: VDbSLDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbSUDa: VDbSUDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbSUDb: VDbSUDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbUDa: VDbUDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbUDb: VDbUDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

VHaM2x2a: VHaM2x2a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHaM2x2b: VHaM2x2b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHaM3x3a: VHaM3x3a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHaM3x3b: VHaM3x3b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHaM4x4a: VHaM4x4a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHaM4x4b: VHaM4x4b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHaM5x5a: VHaM5x5a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHaM5x5b: VHaM5x5b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHaM6x6a: VHaM6x6a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHaM6x6b: VHaM6x6b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHaM7x13a: VHaM7x13a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHaM7x13b: VHaM7x13b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHaM16x8a: VHaM16x8a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHaM16x8b: VHaM16x8b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHaMDa: VHaMDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHaMDb: VHaMDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHaMHa: VHaMHa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHaMHb: VHaMHb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHaMUa: VHaMUa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHaMUb: VHaMUb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHbM2x2a: VHbM2x2a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHbM2x2b: VHbM2x2b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHbM3x3a: VHbM3x3a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHbM3x3b: VHbM3x3b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHbM4x4a: VHbM4x4a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHbM4x4b: VHbM4x4b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHbM5x5a: VHbM5x5a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHbM5x5b: VHbM5x5b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHbM6x6a: VHbM6x6a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHbM6x6b: VHbM6x6b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHbM7x13a: VHbM7x13a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHbM7x13b: VHbM7x13b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHbM16x8a: VHbM16x8a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHbM16x8b: VHbM16x8b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHbMDa: VHbMDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHbMDb: VHbMDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHbMHa: VHbMHa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHbMHb: VHbMHb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHbMUa: VHbMUa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHbMUb: VHbMUb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

VUaM2x2a: VUaM2x2a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUaM2x2b: VUaM2x2b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUaM3x3a: VUaM3x3a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUaM3x3b: VUaM3x3b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUaM4x4a: VUaM4x4a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUaM4x4b: VUaM4x4b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUaM5x5a: VUaM5x5a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUaM5x5b: VUaM5x5b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUaM6x6a: VUaM6x6a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUaM6x6b: VUaM6x6b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUaM7x13a: VUaM7x13a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUaM7x13b: VUaM7x13b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUaM16x8a: VUaM16x8a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUaM16x8b: VUaM16x8b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUaMDa: VUaMDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUaMDb: VUaMDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUaMHa: VUaMHa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUaMHb: VUaMHb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUaMUa: VUaMUa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUaMUb: VUaMUb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUbM2x2a: VUbM2x2a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUbM2x2b: VUbM2x2b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUbM3x3a: VUbM3x3a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUbM3x3b: VUbM3x3b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUbM4x4a: VUbM4x4a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUbM4x4b: VUbM4x4b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUbM5x5a: VUbM5x5a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUbM5x5b: VUbM5x5b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUbM6x6a: VUbM6x6a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUbM6x6b: VUbM6x6b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUbM7x13a: VUbM7x13a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUbM7x13b: VUbM7x13b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUbM16x8a: VUbM16x8a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUbM16x8b: VUbM16x8b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUbMDa: VUbMDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUbMDb: VUbMDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUbMHa: VUbMHa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUbMHb: VUbMHb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUbMUa: VUbMUa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUbMUb: VUbMUb.o
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
