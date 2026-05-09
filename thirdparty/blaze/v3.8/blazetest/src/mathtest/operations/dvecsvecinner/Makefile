#==================================================================================================
#
#  Makefile for the dvecsvecinner module of the Blaze test suite
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
default: V3aVCa V3aVCb V3bVCa V3bVCb VHaVCa VHaVCb VHbVCa VHbVCb \
         VDaVCa VDaVCb VDbVCa VDbVCb VUaVCa VUaVCb VUbVCa VUbVCb
all: $(BIN)
essential: V3aVCa VHaVCa VDaVCa VUaVCa
single: VDaVCa


# Build rules
V2aVCa: V2aVCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V2aVCb: V2aVCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V2bVCa: V2bVCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V2bVCb: V2bVCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

V3aVCa: V3aVCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V3aVCb: V3aVCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V3bVCa: V3bVCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V3bVCb: V3bVCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

V4aVCa: V4aVCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V4aVCb: V4aVCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V4bVCa: V4bVCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V4bVCb: V4bVCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

V5aVCa: V5aVCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V5aVCb: V5aVCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V5bVCa: V5bVCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V5bVCb: V5bVCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

V6aVCa: V6aVCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V6aVCb: V6aVCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V6bVCa: V6bVCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
V6bVCb: V6bVCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

VDaVCa: VDaVCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDaVCb: VDaVCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbVCa: VDbVCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VDbVCb: VDbVCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

VHaVCa: VHaVCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHaVCb: VHaVCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHbVCa: VHbVCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VHbVCb: VHbVCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)

VUaVCa: VUaVCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUaVCb: VUaVCb.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUbVCa: VUbVCa.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
VUbVCb: VUbVCb.o
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
