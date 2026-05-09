#==================================================================================================
#
#  Makefile for the smatsmatmult module of the Blaze test suite
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
default: MCaMCa MCaMCb MCbMCa MCbMCb \
         SCaSCa SCaSCb SCbSCa SCbSCb \
         HCaHCa HCaHCb HCbHCa HCbHCb \
         LCaLCa LCaLCb LCbLCa LCbLCb \
         UCaUCa UCaUCb UCbUCa UCbUCb \
         DCaDCa DCaDCb DCbDCa DCbDCb \
         AliasingTest
all: $(BIN)
essential: MCaMCa SCaSCa HCaHCa LCaLCa UCaUCa DCaDCa AliasingTest
single: MCaMCa


# Build rules
DCaDCa: DCaDCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
DCaDCb: DCaDCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
DCaHCa: DCaHCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
DCaHCb: DCaHCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
DCaLCa: DCaLCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
DCaLCb: DCaLCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
DCaMCa: DCaMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
DCaMCb: DCaMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
DCaMIa: DCaMIa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
DCaMIb: DCaMIb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
DCaSCa: DCaSCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
DCaSCb: DCaSCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
DCaUCa: DCaUCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
DCaUCb: DCaUCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
DCbDCa: DCbDCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
DCbDCb: DCbDCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
DCbHCa: DCbHCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
DCbHCb: DCbHCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
DCbLCa: DCbLCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
DCbLCb: DCbLCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
DCbMCa: DCbMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
DCbMCb: DCbMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
DCbMIa: DCbMIa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
DCbMIb: DCbMIb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
DCbSCa: DCbSCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
DCbSCb: DCbSCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
DCbUCa: DCbUCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
DCbUCb: DCbUCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

HCaDCa: HCaDCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
HCaDCb: HCaDCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
HCaHCa: HCaHCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
HCaHCb: HCaHCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
HCaLCa: HCaLCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
HCaLCb: HCaLCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
HCaMCa: HCaMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
HCaMCb: HCaMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
HCaMIa: HCaMIa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
HCaMIb: HCaMIb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
HCaSCa: HCaSCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
HCaSCb: HCaSCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
HCaUCa: HCaUCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
HCaUCb: HCaUCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
HCbDCa: HCbDCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
HCbDCb: HCbDCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
HCbHCa: HCbHCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
HCbHCb: HCbHCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
HCbLCa: HCbLCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
HCbLCb: HCbLCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
HCbMCa: HCbMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
HCbMCb: HCbMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
HCbMIa: HCbMIa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
HCbMIb: HCbMIb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
HCbSCa: HCbSCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
HCbSCb: HCbSCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
HCbUCa: HCbUCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
HCbUCb: HCbUCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

LCaDCa: LCaDCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
LCaDCb: LCaDCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
LCaHCa: LCaHCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
LCaHCb: LCaHCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
LCaLCa: LCaLCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
LCaLCb: LCaLCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
LCaMCa: LCaMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
LCaMCb: LCaMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
LCaMIa: LCaMIa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
LCaMIb: LCaMIb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
LCaSCa: LCaSCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
LCaSCb: LCaSCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
LCaUCa: LCaUCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
LCaUCb: LCaUCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
LCbDCa: LCbDCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
LCbDCb: LCbDCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
LCbHCa: LCbHCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
LCbHCb: LCbHCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
LCbLCa: LCbLCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
LCbLCb: LCbLCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
LCbMCa: LCbMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
LCbMCb: LCbMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
LCbMIa: LCbMIa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
LCbMIb: LCbMIb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
LCbSCa: LCbSCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
LCbSCb: LCbSCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
LCbUCa: LCbUCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
LCbUCb: LCbUCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

MCaDCa: MCaDCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MCaDCb: MCaDCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MCaHCa: MCaHCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MCaHCb: MCaHCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MCaLCa: MCaLCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MCaLCb: MCaLCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MCaMCa: MCaMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MCaMCb: MCaMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MCaMIa: MCaMIa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MCaMIb: MCaMIb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MCaSCa: MCaSCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MCaSCb: MCaSCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MCaUCa: MCaUCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MCaUCb: MCaUCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MCbDCa: MCbDCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MCbDCb: MCbDCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MCbHCa: MCbHCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MCbHCb: MCbHCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MCbLCa: MCbLCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MCbLCb: MCbLCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MCbMCa: MCbMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MCbMCb: MCbMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MCbMIa: MCbMIa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MCbMIb: MCbMIb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MCbSCa: MCbSCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MCbSCb: MCbSCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MCbUCa: MCbUCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MCbUCb: MCbUCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

MIaDCa: MIaDCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MIaDCb: MIaDCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MIaHCa: MIaHCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MIaHCb: MIaHCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MIaLCa: MIaLCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MIaLCb: MIaLCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MIaMCa: MIaMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MIaMCb: MIaMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MIaSCa: MIaSCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MIaSCb: MIaSCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MIaUCa: MIaUCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MIaUCb: MIaUCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MIbDCa: MIbDCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MIbDCb: MIbDCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MIbHCa: MIbHCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MIbHCb: MIbHCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MIbLCa: MIbLCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MIbLCb: MIbLCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MIbMCa: MIbMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MIbMCb: MIbMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MIbSCa: MIbSCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MIbSCb: MIbSCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MIbUCa: MIbUCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
MIbUCb: MIbUCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

SCaDCa: SCaDCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SCaDCb: SCaDCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SCaHCa: SCaHCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SCaHCb: SCaHCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SCaLCa: SCaLCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SCaLCb: SCaLCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SCaMCa: SCaMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SCaMCb: SCaMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SCaMIa: SCaMIa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SCaMIb: SCaMIb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SCaSCa: SCaSCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SCaSCb: SCaSCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SCaUCa: SCaUCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SCaUCb: SCaUCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SCbDCa: SCbDCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SCbDCb: SCbDCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SCbHCa: SCbHCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SCbHCb: SCbHCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SCbLCa: SCbLCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SCbLCb: SCbLCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SCbMCa: SCbMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SCbMCb: SCbMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SCbMIa: SCbMIa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SCbMIb: SCbMIb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SCbSCa: SCbSCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SCbSCb: SCbSCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SCbUCa: SCbUCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SCbUCb: SCbUCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

UCaDCa: UCaDCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UCaDCb: UCaDCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UCaHCa: UCaHCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UCaHCb: UCaHCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UCaLCa: UCaLCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UCaLCb: UCaLCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UCaMCa: UCaMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UCaMCb: UCaMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UCaMIa: UCaMIa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UCaMIb: UCaMIb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UCaSCa: UCaSCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UCaSCb: UCaSCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UCaUCa: UCaUCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UCaUCb: UCaUCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UCbDCa: UCbDCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UCbDCb: UCbDCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UCbHCa: UCbHCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UCbHCb: UCbHCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UCbLCa: UCbLCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UCbLCb: UCbLCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UCbMCa: UCbMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UCbMCb: UCbMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UCbMIa: UCbMIa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UCbMIb: UCbMIb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UCbSCa: UCbSCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UCbSCb: UCbSCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UCbUCa: UCbUCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UCbUCb: UCbUCb.o
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
