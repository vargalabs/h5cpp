#==================================================================================================
#
#  Makefile for the vectors module of the Blaze test suite
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


# Build rules
default: all

all: densevector staticvector hybridvector dynamicvector customvector uniformvector initializervector \
     sparsevector compressedvector zerovector \
     vectorserializer

essential: all

single: all


# Internal rules
densevector:
	@echo
	@echo "Building the DenseVector operation tests..."
	@$(MAKE) --no-print-directory -C ./densevector $(MAKECMDGOALS)

staticvector:
	@echo
	@echo "Building the StaticVector tests..."
	@$(MAKE) --no-print-directory -C ./staticvector $(MAKECMDGOALS)

hybridvector:
	@echo
	@echo "Building the HybridVector tests..."
	@$(MAKE) --no-print-directory -C ./hybridvector $(MAKECMDGOALS)

dynamicvector:
	@echo
	@echo "Building the DynamicVector tests..."
	@$(MAKE) --no-print-directory -C ./dynamicvector $(MAKECMDGOALS)

customvector:
	@echo
	@echo "Building the CustomVector tests..."
	@$(MAKE) --no-print-directory -C ./customvector $(MAKECMDGOALS)

uniformvector:
	@echo
	@echo "Building the UniformVector tests..."
	@$(MAKE) --no-print-directory -C ./uniformvector $(MAKECMDGOALS)

initializervector:
	@echo
	@echo "Building the InitializerVector tests..."
	@$(MAKE) --no-print-directory -C ./initializervector $(MAKECMDGOALS)

sparsevector:
	@echo
	@echo "Building the SparseVector operation tests..."
	@$(MAKE) --no-print-directory -C ./sparsevector $(MAKECMDGOALS)

compressedvector:
	@echo
	@echo "Building the CompressedVector tests..."
	@$(MAKE) --no-print-directory -C ./compressedvector $(MAKECMDGOALS)

zerovector:
	@echo
	@echo "Building the ZeroVector tests..."
	@$(MAKE) --no-print-directory -C ./zerovector $(MAKECMDGOALS)

vectorserializer:
	@echo
	@echo "Building the VectorSerializer class tests..."
	@$(MAKE) --no-print-directory -C ./vectorserializer $(MAKECMDGOALS)


# Cleanup
reset:
	@$(MAKE) --no-print-directory -C ./densevector reset
	@$(MAKE) --no-print-directory -C ./staticvector reset
	@$(MAKE) --no-print-directory -C ./hybridvector reset
	@$(MAKE) --no-print-directory -C ./dynamicvector reset
	@$(MAKE) --no-print-directory -C ./customvector reset
	@$(MAKE) --no-print-directory -C ./uniformvector reset
	@$(MAKE) --no-print-directory -C ./initializervector reset
	@$(MAKE) --no-print-directory -C ./sparsevector reset
	@$(MAKE) --no-print-directory -C ./compressedvector reset
	@$(MAKE) --no-print-directory -C ./zerovector reset
	@$(MAKE) --no-print-directory -C ./vectorserializer reset

clean:
	@$(MAKE) --no-print-directory -C ./densevector clean
	@$(MAKE) --no-print-directory -C ./staticvector clean
	@$(MAKE) --no-print-directory -C ./hybridvector clean
	@$(MAKE) --no-print-directory -C ./dynamicvector clean
	@$(MAKE) --no-print-directory -C ./customvector clean
	@$(MAKE) --no-print-directory -C ./uniformvector clean
	@$(MAKE) --no-print-directory -C ./initializervector clean
	@$(MAKE) --no-print-directory -C ./sparsevector clean
	@$(MAKE) --no-print-directory -C ./compressedvector clean
	@$(MAKE) --no-print-directory -C ./zerovector clean
	@$(MAKE) --no-print-directory -C ./vectorserializer clean


# Setting the independent commands
.PHONY: default all essential single reset clean \
        densevector staticvector hybridvector dynamicvector customvector uniformvector initializervector \
        sparsevector compressedvector zerovector \
        vectorserializer
