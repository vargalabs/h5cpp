#==================================================================================================
#
#  Makefile for the SIMD module of the Blaze test suite
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
-include ../../Makeconfig
endif
endif


# Setting the source, object and dependency files
SRC = $(wildcard ./*.cpp)
DEP = $(SRC:.cpp=.d)
OBJ = $(SRC:.cpp=.o)
BIN = $(SRC:.cpp=)


# General rules
default: all
all: $(BIN)
essential: $(BIN)
single: $(BIN)


# Build rules
Char: Char.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
ComplexChar: ComplexChar.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
ComplexFloat: ComplexFloat.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
ComplexDouble: ComplexDouble.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
ComplexInt: ComplexInt.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
ComplexLong: ComplexLong.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
ComplexShort: ComplexShort.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
ComplexSignedChar: ComplexSignedChar.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
ComplexUnsignedChar: ComplexUnsignedChar.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
ComplexUnsignedInt: ComplexUnsignedInt.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
ComplexUnsignedLong: ComplexUnsignedLong.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
ComplexUnsignedShort: ComplexUnsignedShort.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
Double: Double.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
Float: Float.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
IncludeTest: IncludeTest.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
Int: Int.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
Long: Long.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
Short: Short.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
SignedChar: SignedChar.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UnsignedChar: UnsignedChar.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UnsignedInt: UnsignedInt.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UnsignedLong: UnsignedLong.o
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBRARIES)
UnsignedShort: UnsignedShort.o
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
