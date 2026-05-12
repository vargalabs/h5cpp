#==================================================================================================
#
#  Makefile for the dmatsvecmult module of the Blaze test suite
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
default: M3x3aVCa M3x3aVCb M3x3bVCa M3x3bVCb \
         MHaVCa MHaVCb MHbVCa MHbVCb \
         MDaVCa MDaVCb MDbVCa MDbVCb \
         MUaVCa MUaVCb MUbVCa MUbVCb \
         SDaVCa SDaVCb SDbVCa SDbVCb \
         HDaVCa HDaVCb HDbVCa HDbVCb \
         LDaVCa LDaVCb LDbVCa LDbVCb \
         UDaVCa UDaVCb UDbVCa UDbVCb \
         DDaVCa DDaVCb DDbVCa DDbVCb \
         AliasingTest
all: $(BIN)
essential: M3x3aVCa MHaVCa MDaVCa MUaVCa SDaVCa HDaVCa LDaVCa UDaVCa DDaVCa AliasingTest
single: MDaVCa


# Build rules
DDaVCa: DDaVCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
DDaVCb: DDaVCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
DDbVCa: DDbVCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
DDbVCb: DDbVCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

LDaVCa: LDaVCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
LDaVCb: LDaVCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
LDbVCa: LDbVCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
LDbVCb: LDbVCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

HDaVCa: HDaVCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
HDaVCb: HDaVCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
HDbVCa: HDbVCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
HDbVCb: HDbVCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

M2x2aVCa: M2x2aVCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M2x2aVCb: M2x2aVCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M2x2bVCa: M2x2bVCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M2x2bVCb: M2x2bVCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

M3x3aVCa: M3x3aVCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M3x3aVCb: M3x3aVCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M3x3bVCa: M3x3bVCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M3x3bVCb: M3x3bVCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

M4x4aVCa: M4x4aVCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M4x4aVCb: M4x4aVCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M4x4bVCa: M4x4bVCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M4x4bVCb: M4x4bVCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

M5x5aVCa: M5x5aVCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M5x5aVCb: M5x5aVCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M5x5bVCa: M5x5bVCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M5x5bVCb: M5x5bVCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

M6x6aVCa: M6x6aVCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M6x6aVCb: M6x6aVCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M6x6bVCa: M6x6bVCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M6x6bVCb: M6x6bVCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

M7x13aVCa: M7x13aVCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M7x13aVCb: M7x13aVCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M7x13bVCa: M7x13bVCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M7x13bVCb: M7x13bVCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

M16x8aVCa: M16x8aVCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M16x8aVCb: M16x8aVCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M16x8bVCa: M16x8bVCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M16x8bVCb: M16x8bVCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

MDaVCa: MDaVCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDaVCb: MDaVCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDbVCa: MDbVCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDbVCb: MDbVCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

MHaVCa: MHaVCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MHaVCb: MHaVCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MHbVCa: MHbVCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MHbVCb: MHbVCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

MUaVCa: MUaVCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MUaVCb: MUaVCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MUbVCa: MUbVCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MUbVCb: MUbVCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

SDaVCa: SDaVCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SDaVCb: SDaVCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SDbVCa: SDbVCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SDbVCb: SDbVCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

SLDaVCa: SLDaVCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SLDaVCb: SLDaVCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SLDbVCa: SLDbVCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SLDbVCb: SLDbVCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

SUDaVCa: SUDaVCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SUDaVCb: SUDaVCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SUDbVCa: SUDbVCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SUDbVCb: SUDbVCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

UDaVCa: UDaVCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UDaVCb: UDaVCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UDbVCa: UDbVCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UDbVCb: UDbVCb.o
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
