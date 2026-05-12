#==================================================================================================
#
#  Makefile for the tsvecsmatmult module of the Blaze test suite
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
default: VCaMCa VCaMCb VCbMCa VCbMCb \
         VCaSCa VCaSCb VCbSCa VCbSCb \
         VCaHCa VCaHCb VCbHCa VCbHCb \
         VCaLCa VCaLCb VCbLCa VCbLCb \
         VCaUCa VCaUCb VCbUCa VCbUCb \
         VCaDCa VCaDCb VCbDCa VCbDCb \
         AliasingTest
all: $(BIN)
essential: VCaMCa VCaSCa VCaHCa VCaLCa VCaUCa VCaDCa AliasingTest
single: VCaMCa


# Build rules
VCaDCa: VCaDCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCaDCb: VCaDCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCaHCa: VCaHCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCaHCb: VCaHCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCaLCa: VCaLCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCaLCb: VCaLCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCaMCa: VCaMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCaMCb: VCaMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCaMIa: VCaMIa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCaMIb: VCaMIb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCaSCa: VCaSCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCaSCb: VCaSCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCaUCa: VCaUCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCaUCb: VCaUCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

VCbDCa: VCbDCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCbDCb: VCbDCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCbHCa: VCbHCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCbHCb: VCbHCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCbLCa: VCbLCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCbLCb: VCbLCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCbMCa: VCbMCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCbMCb: VCbMCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCbMIa: VCbMIa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCbMIb: VCbMIb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCbSCa: VCbSCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCbSCb: VCbSCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCbUCa: VCbUCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VCbUCb: VCbUCb.o
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
