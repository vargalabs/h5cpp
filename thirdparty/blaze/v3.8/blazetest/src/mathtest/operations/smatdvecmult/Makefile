#==================================================================================================
#
#  Makefile for the smatdvecmult module of the Blaze test suite
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
default: MCaV3a MCaV3b MCbV3a MCbV3b \
         MCaVHa MCaVHb MCbVHa MCbVHb \
         MCaVDa MCaVDb MCbVDa MCbVDb \
         MCaVUa MCaVUb MCbVUa MCbVUb \
         SCaVDa SCaVDb SCbVDa SCbVDb \
         HCaVDa HCaVDb HCbVDa HCbVDb \
         LCaVDa LCaVDb LCbVDa LCbVDb \
         UCaVDa UCaVDb UCbVDa UCbVDb \
         DCaVDa DCaVDb DCbVDa DCbVDb \
         AliasingTest
all: $(BIN)
essential: MCaV3a MCaVHa MCaVDa MCaVUa SCaVDa HCaVDa LCaVDa UCaVDa DCaVDa AliasingTest
single: MCaVDa


# Build rules
DCaVDa: DCaVDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
DCaVDb: DCaVDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
DCbVDa: DCbVDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
DCbVDb: DCbVDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

HCaVDa: HCaVDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
HCaVDb: HCaVDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
HCbVDa: HCbVDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
HCbVDb: HCbVDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

LCaVDa: LCaVDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
LCaVDb: LCaVDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
LCbVDa: LCbVDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
LCbVDb: LCbVDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

MCaV2a: MCaV2a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MCaV2b: MCaV2b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MCaV3a: MCaV3a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MCaV3b: MCaV3b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MCaV4a: MCaV4a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MCaV4b: MCaV4b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MCaV5a: MCaV5a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MCaV5b: MCaV5b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MCaV6a: MCaV6a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MCaV6b: MCaV6b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MCaVDa: MCaVDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MCaVDb: MCaVDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MCaVHa: MCaVHa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MCaVHb: MCaVHb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MCaVUa: MCaVUa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MCaVUb: MCaVUb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

MCbV2a: MCbV2a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MCbV2b: MCbV2b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MCbV3a: MCbV3a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MCbV3b: MCbV3b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MCbV4a: MCbV4a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MCbV4b: MCbV4b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MCbV5a: MCbV5a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MCbV5b: MCbV5b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MCbV6a: MCbV6a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MCbV6b: MCbV6b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MCbVDa: MCbVDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MCbVDb: MCbVDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MCbVHa: MCbVHa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MCbVHb: MCbVHb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MCbVUa: MCbVUa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MCbVUb: MCbVUb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

MIaV3a: MIaV3a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MIaV3b: MIaV3b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MIaVDa: MIaVDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MIaVDb: MIaVDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MIaVHa: MIaVHa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MIaVHb: MIaVHb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MIaVUa: MIaVUa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MIaVUb: MIaVUb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MIbV3a: MIbV3a.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MIbV3b: MIbV3b.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MIbVDa: MIbVDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MIbVDb: MIbVDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MIbVHa: MIbVHa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MIbVHb: MIbVHb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MIbVUa: MIbVUa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MIbVUb: MIbVUb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

SCaVDa: SCaVDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SCaVDb: SCaVDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SCbVDa: SCbVDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SCbVDb: SCbVDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

UCaVDa: UCaVDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UCaVDb: UCaVDb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UCbVDa: UCbVDa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UCbVDb: UCbVDb.o
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
