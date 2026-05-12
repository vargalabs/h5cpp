#==================================================================================================
#
#  Makefile for the mathtest module of the Blaze test suite
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

all: shims simd blas lapack typetraits traits constraints functors \
     vectors matrices views adaptors operations

essential: all

single: all


# Internal rules
shims:
	@echo
	@echo "Building the mathematical shims operation tests..."
	@$(MAKE) --no-print-directory -C ./shims $(MAKECMDGOALS)

simd:
	@echo
	@echo "Building the SIMD operation tests..."
	@$(MAKE) --no-print-directory -C ./simd $(MAKECMDGOALS)

blas:
	@echo
	@echo "Building the BLAS operation tests..."
	@$(MAKE) --no-print-directory -C ./blas $(MAKECMDGOALS)

lapack:
	@echo
	@echo "Building the LAPACK operation tests..."
	@$(MAKE) --no-print-directory -C ./lapack $(MAKECMDGOALS)

typetraits:
	@echo
	@echo "Building the typetraits operation tests..."
	@$(MAKE) --no-print-directory -C ./typetraits $(MAKECMDGOALS)

traits:
	@$(MAKE) --no-print-directory -C ./traits $(MAKECMDGOALS)

constraints:
	@echo
	@echo "Building the constraints operation tests..."
	@$(MAKE) --no-print-directory -C ./constraints $(MAKECMDGOALS)

functors:
	@echo
	@echo "Building the functors operation tests..."
	@$(MAKE) --no-print-directory -C ./functors $(MAKECMDGOALS)

vectors:
	@$(MAKE) --no-print-directory -C ./vectors $(MAKECMDGOALS)

matrices:
	@$(MAKE) --no-print-directory -C ./matrices $(MAKECMDGOALS)

views:
	@$(MAKE) --no-print-directory -C ./views $(MAKECMDGOALS)

adaptors:
	@$(MAKE) --no-print-directory -C ./adaptors $(MAKECMDGOALS)

operations:
	@$(MAKE) --no-print-directory -C ./operations $(MAKECMDGOALS)


# Cleanup
reset:
	@$(MAKE) --no-print-directory -C ./shims reset
	@$(MAKE) --no-print-directory -C ./simd reset
	@$(MAKE) --no-print-directory -C ./blas reset
	@$(MAKE) --no-print-directory -C ./lapack reset
	@$(MAKE) --no-print-directory -C ./typetraits reset
	@$(MAKE) --no-print-directory -C ./traits reset
	@$(MAKE) --no-print-directory -C ./constraints reset
	@$(MAKE) --no-print-directory -C ./functors reset
	@$(MAKE) --no-print-directory -C ./vectors reset
	@$(MAKE) --no-print-directory -C ./matrices reset
	@$(MAKE) --no-print-directory -C ./views reset
	@$(MAKE) --no-print-directory -C ./adaptors reset
	@$(MAKE) --no-print-directory -C ./operations reset

clean:
	@$(MAKE) --no-print-directory -C ./shims clean
	@$(MAKE) --no-print-directory -C ./simd clean
	@$(MAKE) --no-print-directory -C ./blas clean
	@$(MAKE) --no-print-directory -C ./lapack clean
	@$(MAKE) --no-print-directory -C ./typetraits clean
	@$(MAKE) --no-print-directory -C ./traits clean
	@$(MAKE) --no-print-directory -C ./constraints clean
	@$(MAKE) --no-print-directory -C ./functors clean
	@$(MAKE) --no-print-directory -C ./vectors clean
	@$(MAKE) --no-print-directory -C ./matrices clean
	@$(MAKE) --no-print-directory -C ./views clean
	@$(MAKE) --no-print-directory -C ./adaptors clean
	@$(MAKE) --no-print-directory -C ./operations clean


# Setting the independent commands
.PHONY: default all essential single reset clean \
        shims simd blas lapack typetraits traits constraints functors \
        vectors matrices views adaptors operations
