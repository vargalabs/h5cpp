#==================================================================================================
#
#  Makefile for the tdvecsmatmult module of the Blaze test suite
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
default: V3aMCa V3aMCb V3bMCa V3bMCb \
         VHaMCa VHaMCb VHbMCa VHbMCb \
         VDaMCa VDaMCb VDbMCa VDbMCb \
         VUaMCa VUaMCb VUbMCa VUbMCb \
         VDaSCa VDaSCb VDbSCa VDbSCb \
         VDaHCa VDaHCb VDbHCa VDbHCb \
         VDaLCa VDaLCb VDbLCa VDbLCb \
         VDaUCa VDaUCb VDbUCa VDbUCb \
         VDaDCa VDaDCb VDbDCa VDbDCb \
         AliasingTest
all: $(BIN)
essential: V3aMCa VHaMCa VDaMCa VUaMCa VDaSCa VDaHCa VDaLCa VDaUCa VDaDCa AliasingTest
single: VDaMCa


# Build rules
V2aMCa: V2aMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V2aMCb: V2aMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V2bMCa: V2bMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V2bMCb: V2bMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

V3aMCa: V3aMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V3aMCb: V3aMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V3aMIa: V3aMIa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V3aMIb: V3aMIb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V3bMCa: V3bMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V3bMCb: V3bMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V3bMIa: V3bMIa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V3bMIb: V3bMIb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

V4aMCa: V4aMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V4aMCb: V4aMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V4bMCa: V4bMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V4bMCb: V4bMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

V5aMCa: V5aMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V5aMCb: V5aMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V5bMCa: V5bMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V5bMCb: V5bMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

V6aMCa: V6aMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V6aMCb: V6aMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V6bMCa: V6bMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V6bMCb: V6bMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

VDaDCa: VDaDCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaDCb: VDaDCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaHCa: VDaHCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaHCb: VDaHCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaLCa: VDaLCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaLCb: VDaLCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaMCa: VDaMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaMCb: VDaMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaMIa: VDaMIa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaMIb: VDaMIb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaSCa: VDaSCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaSCb: VDaSCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaUCa: VDaUCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaUCb: VDaUCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbDCa: VDbDCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbDCb: VDbDCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbHCa: VDbHCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbHCb: VDbHCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbLCa: VDbLCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbLCb: VDbLCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbMCa: VDbMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbMCb: VDbMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbMIa: VDbMIa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbMIb: VDbMIb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbSCa: VDbSCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbSCb: VDbSCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbUCa: VDbUCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbUCb: VDbUCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

VHaMCa: VHaMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHaMCb: VHaMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHaMIa: VHaMIa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHaMIb: VHaMIb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHbMCa: VHbMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHbMCb: VHbMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHbMIa: VHbMIa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHbMIb: VHbMIb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

VUaMCa: VUaMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUaMCb: VUaMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUaMIa: VUaMIa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUaMIb: VUaMIb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUbMCa: VUbMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUbMCb: VUbMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUbMIa: VUbMIa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUbMIb: VUbMIb.o
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
