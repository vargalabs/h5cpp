#==================================================================================================
#
#  Makefile for the matrices module of the Blaze test suite
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

all: densematrix staticmatrix hybridmatrix dynamicmatrix custommatrix uniformmatrix initializermatrix \
     sparsematrix compressedmatrix identitymatrix zeromatrix \
     matrixserializer

essential: all

single: all


# Internal rules
densematrix:
	@echo
	@echo "Building the DenseMatrix operation tests..."
	@$(MAKE) --no-print-directory -C ./densematrix $(MAKECMDGOALS)

staticmatrix:
	@echo
	@echo "Building the StaticMatrix tests..."
	@$(MAKE) --no-print-directory -C ./staticmatrix $(MAKECMDGOALS)

hybridmatrix:
	@echo
	@echo "Building the HybridMatrix tests..."
	@$(MAKE) --no-print-directory -C ./hybridmatrix $(MAKECMDGOALS)

dynamicmatrix:
	@echo
	@echo "Building the DynamicMatrix tests..."
	@$(MAKE) --no-print-directory -C ./dynamicmatrix $(MAKECMDGOALS)

custommatrix:
	@echo
	@echo "Building the CustomMatrix tests..."
	@$(MAKE) --no-print-directory -C ./custommatrix $(MAKECMDGOALS)

uniformmatrix:
	@echo
	@echo "Building the UniformMatrix tests..."
	@$(MAKE) --no-print-directory -C ./uniformmatrix $(MAKECMDGOALS)

initializermatrix:
	@echo
	@echo "Building the InitializerMatrix tests..."
	@$(MAKE) --no-print-directory -C ./initializermatrix $(MAKECMDGOALS)

sparsematrix:
	@echo
	@echo "Building the SparseMatrix operation tests..."
	@$(MAKE) --no-print-directory -C ./sparsematrix $(MAKECMDGOALS)

compressedmatrix:
	@echo
	@echo "Building the CompressedMatrix tests..."
	@$(MAKE) --no-print-directory -C ./compressedmatrix $(MAKECMDGOALS)

identitymatrix:
	@echo
	@echo "Building the IdentityMatrix tests..."
	@$(MAKE) --no-print-directory -C ./identitymatrix $(MAKECMDGOALS)

zeromatrix:
	@echo
	@echo "Building the ZeroMatrix tests..."
	@$(MAKE) --no-print-directory -C ./zeromatrix $(MAKECMDGOALS)

matrixserializer:
	@echo
	@echo "Building the MatrixSerializer class tests..."
	@$(MAKE) --no-print-directory -C ./matrixserializer $(MAKECMDGOALS)


# Cleanup
reset:
	@$(MAKE) --no-print-directory -C ./densematrix reset
	@$(MAKE) --no-print-directory -C ./staticmatrix reset
	@$(MAKE) --no-print-directory -C ./hybridmatrix reset
	@$(MAKE) --no-print-directory -C ./dynamicmatrix reset
	@$(MAKE) --no-print-directory -C ./custommatrix reset
	@$(MAKE) --no-print-directory -C ./uniformmatrix reset
	@$(MAKE) --no-print-directory -C ./initializermatrix reset
	@$(MAKE) --no-print-directory -C ./sparsematrix reset
	@$(MAKE) --no-print-directory -C ./compressedmatrix reset
	@$(MAKE) --no-print-directory -C ./identitymatrix reset
	@$(MAKE) --no-print-directory -C ./zeromatrix reset
	@$(MAKE) --no-print-directory -C ./matrixserializer reset

clean:
	@$(MAKE) --no-print-directory -C ./densematrix clean
	@$(MAKE) --no-print-directory -C ./staticmatrix clean
	@$(MAKE) --no-print-directory -C ./hybridmatrix clean
	@$(MAKE) --no-print-directory -C ./dynamicmatrix clean
	@$(MAKE) --no-print-directory -C ./custommatrix clean
	@$(MAKE) --no-print-directory -C ./uniformmatrix clean
	@$(MAKE) --no-print-directory -C ./initializermatrix clean
	@$(MAKE) --no-print-directory -C ./sparsematrix clean
	@$(MAKE) --no-print-directory -C ./compressedmatrix clean
	@$(MAKE) --no-print-directory -C ./identitymatrix clean
	@$(MAKE) --no-print-directory -C ./zeromatrix clean
	@$(MAKE) --no-print-directory -C ./matrixserializer clean


# Setting the independent commands
.PHONY: default all essential single reset clean \
        densematrix staticmatrix hybridmatrix dynamicmatrix custommatrix uniformmatrix initializermatrix \
        sparsematrix compressedmatrix identitymatrix zeromatrix \
        matrixserializer
