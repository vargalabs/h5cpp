#==================================================================================================
#
#  Makefile for the tsvecdmatmult module of the Blaze test suite
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
default: VCaM3x3a VCaM3x3b VCbM3x3a VCbM3x3b \
         VCaMHa VCaMHb VCbMHa VCbMHb \
         VCaMDa VCaMDb VCbMDa VCbMDb \
         VCaMUa VCaMUb VCbMUa VCbMUb \
         VCaSDa VCaSDb VCbSDa VCbSDb \
         VCaHDa VCaHDb VCbHDa VCbHDb \
         VCaLDa VCaLDb VCbLDa VCbLDb \
         VCaUDa VCaUDb VCbUDa VCbUDb \
         VCaDDa VCaDDb VCbDDa VCbDDb \
         AliasingTest
all: $(BIN)
essential: VCaM3x3a VCaMHa VCaMDa VCaMUa VCaSDa VCaHDa VCaLDa VCaUDa VCaDDa AliasingTest
single: VCaMDa


# Build rules
VCaDDa: VCaDDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCaDDb: VCaDDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCaHDa: VCaHDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCaHDb: VCaHDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCaLDa: VCaLDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCaLDb: VCaLDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCaM2x2a: VCaM2x2a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCaM2x2b: VCaM2x2b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCaM3x3a: VCaM3x3a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCaM3x3b: VCaM3x3b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCaM4x4a: VCaM4x4a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCaM4x4b: VCaM4x4b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCaM5x5a: VCaM5x5a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCaM5x5b: VCaM5x5b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCaM6x6a: VCaM6x6a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCaM6x6b: VCaM6x6b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCaM7x13a: VCaM7x13a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCaM7x13b: VCaM7x13b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCaM16x8a: VCaM16x8a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCaM16x8b: VCaM16x8b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCaMDa: VCaMDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCaMDb: VCaMDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCaMHa: VCaMHa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCaMHb: VCaMHb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCaMUa: VCaMUa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCaMUb: VCaMUb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCaSDa: VCaSDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCaSDb: VCaSDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCaSLDa: VCaSLDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCaSLDb: VCaSLDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCaSUDa: VCaSUDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCaSUDb: VCaSUDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCaUDa: VCaUDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCaUDb: VCaUDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

VCbDDa: VCbDDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCbDDb: VCbDDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCbHDa: VCbHDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCbHDb: VCbHDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCbLDa: VCbLDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCbLDb: VCbLDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCbM2x2a: VCbM2x2a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCbM2x2b: VCbM2x2b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCbM3x3a: VCbM3x3a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCbM3x3b: VCbM3x3b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCbM4x4a: VCbM4x4a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCbM4x4b: VCbM4x4b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCbM5x5a: VCbM5x5a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCbM5x5b: VCbM5x5b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCbM6x6a: VCbM6x6a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCbM6x6b: VCbM6x6b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCbM7x13a: VCbM7x13a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCbM7x13b: VCbM7x13b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCbM16x8a: VCbM16x8a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCbM16x8b: VCbM16x8b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCbMDa: VCbMDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCbMDb: VCbMDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCbMHa: VCbMHa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCbMHb: VCbMHb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCbMUa: VCbMUa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCbMUb: VCbMUb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCbSDa: VCbSDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCbSDb: VCbSDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCbSLDa: VCbSLDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCbSLDb: VCbSLDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCbSUDa: VCbSUDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCbSUDb: VCbSUDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCbUDa: VCbUDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCbUDb: VCbUDb.o
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
