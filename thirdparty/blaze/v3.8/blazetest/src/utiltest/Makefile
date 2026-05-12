#==================================================================================================
#
#  Makefile for the utiltest module of the Blaze test suite
#
#  Copyright (C) 2012-2020 Klaus Iglberger - All Rights Reserved
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


# Build rules
default: all

all: constraints alignedallocator memory numericcast smallarray typetraits valuetraits

essential: all

single: all


# Internal rules
alignedallocator:
	@echo
	@echo "Building the aligned allocator tests..."
	@$(MAKE) --no-print-directory -C ./alignedallocator $(MAKECMDGOALS)

memory:
	@echo
	@echo "Building the memory tests..."
	@$(MAKE) --no-print-directory -C ./memory $(MAKECMDGOALS)

numericcast:
	@echo
	@echo "Building the numeric cast tests..."
	@$(MAKE) --no-print-directory -C ./numericcast $(MAKECMDGOALS)

smallarray:
	@echo
	@echo "Building the small array tests..."
	@$(MAKE) --no-print-directory -C ./smallarray $(MAKECMDGOALS)

typetraits:
	@echo
	@echo "Building the type traits tests..."
	@$(MAKE) --no-print-directory -C ./typetraits $(MAKECMDGOALS)

valuetraits:
	@echo
	@echo "Building the value traits tests..."
	@$(MAKE) --no-print-directory -C ./valuetraits $(MAKECMDGOALS)


# Cleanup
reset:
	@$(MAKE) --no-print-directory -C ./alignedallocator reset
	@$(MAKE) --no-print-directory -C ./memory reset
	@$(MAKE) --no-print-directory -C ./numericcast reset
	@$(MAKE) --no-print-directory -C ./smallarray reset
	@$(MAKE) --no-print-directory -C ./typetraits reset
	@$(MAKE) --no-print-directory -C ./valuetraits reset

clean:
	@$(MAKE) --no-print-directory -C ./alignedallocator clean
	@$(MAKE) --no-print-directory -C ./memory clean
	@$(MAKE) --no-print-directory -C ./numericcast clean
	@$(MAKE) --no-print-directory -C ./smallarray clean
	@$(MAKE) --no-print-directory -C ./typetraits clean
	@$(MAKE) --no-print-directory -C ./valuetraits clean


# Setting the independent commands
.PHONY: default all essential single reset clean \
        alignedallocator memory numericcast smallarray typetraits valuetraits
