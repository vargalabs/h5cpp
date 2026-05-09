#==================================================================================================
#
#  Makefile for the dmatsmatschur module of the Blaze test suite
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
default: M3x3aMCa M3x3aMCb M3x3bMCa M3x3bMCb \
         MHaMCa MHaMCb MHbMCa MHbMCb \
         MDaMCa MDaMCb MDbMCa MDbMCb \
         MUaMCa MUaMCb MUbMCa MUbMCb \
         SDaSCa SDaSCb SDbSCa SDbSCb \
         HDaHCa HDaHCb HDbHCa HDbHCb \
         LDaLCa LDaLCb LDbLCa LDbLCb \
         UDaUCa UDaUCb UDbUCa UDbUCb \
         DDaDCa DDaDCb DDbDCa DDbDCb \
         AliasingTest
all: $(BIN)
essential: M3x3aMCa MHaMCa MDaMCa MUaMCa SDaSCa HDaHCa LDaLCa UDaUCa DDaDCa AliasingTest
single: MDaMCa


# Build rules
D3x3aDCa: D3x3aDCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
D3x3aDCb: D3x3aDCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
D3x3bDCa: D3x3bDCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
D3x3bDCb: D3x3bDCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

DDaDCa: DDaDCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
DDaDCb: DDaDCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
DDaHCa: DDaHCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
DDaHCb: DDaHCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
DDaLCa: DDaLCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
DDaLCb: DDaLCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
DDaMCa: DDaMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
DDaMCb: DDaMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
DDaMIa: DDaMIa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
DDaMIb: DDaMIb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
DDaSCa: DDaSCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
DDaSCb: DDaSCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
DDaUCa: DDaUCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
DDaUCb: DDaUCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
DDbDCa: DDbDCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
DDbDCb: DDbDCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
DDbHCa: DDbHCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
DDbHCb: DDbHCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
DDbLCa: DDbLCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
DDbLCb: DDbLCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
DDbMCa: DDbMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
DDbMCb: DDbMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
DDbMIa: DDbMIa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
DDbMIb: DDbMIb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
DDbSCa: DDbSCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
DDbSCb: DDbSCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
DDbUCa: DDbUCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
DDbUCb: DDbUCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

DHaDCa: DHaDCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
DHaDCb: DHaDCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
DHbDCa: DHbDCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
DHbDCb: DHbDCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

H3x3aHCa: H3x3aHCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
H3x3aHCb: H3x3aHCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
H3x3bHCa: H3x3bHCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
H3x3bHCb: H3x3bHCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

HDaDCa: HDaDCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
HDaDCb: HDaDCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
HDaHCa: HDaHCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
HDaHCb: HDaHCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
HDaLCa: HDaLCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
HDaLCb: HDaLCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
HDaMCa: HDaMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
HDaMCb: HDaMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
HDaMIa: HDaMIa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
HDaMIb: HDaMIb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
HDaSCa: HDaSCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
HDaSCb: HDaSCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
HDaUCa: HDaUCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
HDaUCb: HDaUCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
HDbDCa: HDbDCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
HDbDCb: HDbDCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
HDbHCa: HDbHCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
HDbHCb: HDbHCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
HDbLCa: HDbLCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
HDbLCb: HDbLCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
HDbMCa: HDbMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
HDbMCb: HDbMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
HDbMIa: HDbMIa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
HDbMIb: HDbMIb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
HDbSCa: HDbSCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
HDbSCb: HDbSCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
HDbUCa: HDbUCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
HDbUCb: HDbUCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

HHaHCa: HHaHCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
HHaHCb: HHaHCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
HHbHCa: HHbHCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
HHbHCb: HHbHCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

L3x3aLCa: L3x3aLCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
L3x3aLCb: L3x3aLCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
L3x3bLCa: L3x3bLCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
L3x3bLCb: L3x3bLCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

LDaDCa: LDaDCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
LDaDCb: LDaDCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
LDaHCa: LDaHCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
LDaHCb: LDaHCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
LDaLCa: LDaLCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
LDaLCb: LDaLCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
LDaMCa: LDaMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
LDaMCb: LDaMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
LDaMIa: LDaMIa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
LDaMIb: LDaMIb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
LDaSCa: LDaSCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
LDaSCb: LDaSCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
LDaUCa: LDaUCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
LDaUCb: LDaUCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
LDbDCa: LDbDCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
LDbDCb: LDbDCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
LDbHCa: LDbHCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
LDbHCb: LDbHCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
LDbLCa: LDbLCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
LDbLCb: LDbLCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
LDbMCa: LDbMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
LDbMCb: LDbMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
LDbMIa: LDbMIa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
LDbMIb: LDbMIb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
LDbSCa: LDbSCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
LDbSCb: LDbSCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
LDbUCa: LDbUCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
LDbUCb: LDbUCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

LHaLCa: LHaLCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
LHaLCb: LHaLCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
LHbLCa: LHbLCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
LHbLCb: LHbLCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

M2x2aMCa: M2x2aMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M2x2aMCb: M2x2aMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M2x2bMCa: M2x2bMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M2x2bMCb: M2x2bMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

M3x3aMCa: M3x3aMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M3x3aMCb: M3x3aMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M3x3aMIa: M3x3aMIa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M3x3aMIb: M3x3aMIb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M3x3bMCa: M3x3bMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M3x3bMCb: M3x3bMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M3x3bMIa: M3x3bMIa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M3x3bMIb: M3x3bMIb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

M4x4aMCa: M4x4aMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M4x4aMCb: M4x4aMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M4x4bMCa: M4x4bMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M4x4bMCb: M4x4bMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

M5x5aMCa: M5x5aMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M5x5aMCb: M5x5aMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M5x5bMCa: M5x5bMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M5x5bMCb: M5x5bMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

M6x6aMCa: M6x6aMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M6x6aMCb: M6x6aMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M6x6bMCa: M6x6bMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M6x6bMCb: M6x6bMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

M7x13aMCa: M7x13aMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M7x13aMCb: M7x13aMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M7x13bMCa: M7x13bMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M7x13bMCb: M7x13bMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

M16x8aMCa: M16x8aMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M16x8aMCb: M16x8aMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M16x8bMCa: M16x8bMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
M16x8bMCb: M16x8bMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

MDaDCa: MDaDCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDaDCb: MDaDCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDaHCa: MDaHCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDaHCb: MDaHCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDaLCa: MDaLCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDaLCb: MDaLCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDaMCa: MDaMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDaMCb: MDaMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDaMIa: MDaMIa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDaMIb: MDaMIb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDaSCa: MDaSCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDaSCb: MDaSCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDaUCa: MDaUCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDaUCb: MDaUCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDbDCa: MDbDCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDbDCb: MDbDCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDbHCa: MDbHCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDbHCb: MDbHCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDbLCa: MDbLCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDbLCb: MDbLCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDbMCa: MDbMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDbMCb: MDbMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDbMIa: MDbMIa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDbMIb: MDbMIb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDbSCa: MDbSCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDbSCb: MDbSCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDbUCa: MDbUCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MDbUCb: MDbUCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

MHaMCa: MHaMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MHaMCb: MHaMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MHaMIa: MHaMIa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MHaMIb: MHaMIb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MHbMCa: MHbMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MHbMCb: MHbMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MHbMIa: MHbMIa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MHbMIb: MHbMIb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

MUaMCa: MUaMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MUaMCb: MUaMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MUaMIa: MUaMIa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MUaMIb: MUaMIb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MUbMCa: MUbMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MUbMCb: MUbMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MUbMIa: MUbMIa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MUbMIb: MUbMIb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

S3x3aSCa: S3x3aSCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
S3x3aSCb: S3x3aSCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
S3x3bSCa: S3x3bSCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
S3x3bSCb: S3x3bSCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

SDaDCa: SDaDCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SDaDCb: SDaDCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SDaHCa: SDaHCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SDaHCb: SDaHCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SDaLCa: SDaLCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SDaLCb: SDaLCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SDaMCa: SDaMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SDaMCb: SDaMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SDaMIa: SDaMIa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SDaMIb: SDaMIb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SDaSCa: SDaSCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SDaSCb: SDaSCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SDaUCa: SDaUCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SDaUCb: SDaUCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SDbDCa: SDbDCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SDbDCb: SDbDCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SDbHCa: SDbHCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SDbHCb: SDbHCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SDbLCa: SDbLCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SDbLCb: SDbLCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SDbMCa: SDbMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SDbMCb: SDbMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SDbMIa: SDbMIa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SDbMIb: SDbMIb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SDbSCa: SDbSCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SDbSCb: SDbSCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SDbUCa: SDbUCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SDbUCb: SDbUCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

SHaSCa: SHaSCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SHaSCb: SHaSCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SHbSCa: SHbSCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SHbSCb: SHbSCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

U3x3aUCa: U3x3aUCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
U3x3aUCb: U3x3aUCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
U3x3bUCa: U3x3bUCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
U3x3bUCb: U3x3bUCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

UDaDCa: UDaDCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UDaDCb: UDaDCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UDaHCa: UDaHCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UDaHCb: UDaHCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UDaLCa: UDaLCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UDaLCb: UDaLCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UDaMCa: UDaMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UDaMCb: UDaMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UDaMIa: UDaMIa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UDaMIb: UDaMIb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UDaSCa: UDaSCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UDaSCb: UDaSCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UDaUCa: UDaUCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UDaUCb: UDaUCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UDbDCa: UDbDCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UDbDCb: UDbDCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UDbHCa: UDbHCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UDbHCb: UDbHCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UDbLCa: UDbLCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UDbLCb: UDbLCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UDbMCa: UDbMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UDbMCb: UDbMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UDbMIa: UDbMIa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UDbMIb: UDbMIb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UDbSCa: UDbSCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UDbSCb: UDbSCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UDbUCa: UDbUCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UDbUCb: UDbUCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

UHaUCa: UHaUCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UHaUCb: UHaUCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UHbUCa: UHbUCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UHbUCb: UHbUCb.o
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
