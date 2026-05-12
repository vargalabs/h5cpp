#==================================================================================================
#
#  Makefile for the views module of the Blaze test suite
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

all: subvector elements submatrix row rows column columns band

essential: all

single: all


# Internal rules
subvector:
	@echo
	@echo "Building the Subvector tests..."
	@$(MAKE) --no-print-directory -C ./subvector $(MAKECMDGOALS)

elements:
	@echo
	@echo "Building the Elements tests..."
	@$(MAKE) --no-print-directory -C ./elements $(MAKECMDGOALS)

submatrix:
	@echo
	@echo "Building the Submatrix tests..."
	@$(MAKE) --no-print-directory -C ./submatrix $(MAKECMDGOALS)

row:
	@echo
	@echo "Building the Row tests..."
	@$(MAKE) --no-print-directory -C ./row $(MAKECMDGOALS)

rows:
	@echo
	@echo "Building the Rows tests..."
	@$(MAKE) --no-print-directory -C ./rows $(MAKECMDGOALS)

column:
	@echo
	@echo "Building the Column tests..."
	@$(MAKE) --no-print-directory -C ./column $(MAKECMDGOALS)

columns:
	@echo
	@echo "Building the Columns tests..."
	@$(MAKE) --no-print-directory -C ./columns $(MAKECMDGOALS)

band:
	@echo
	@echo "Building the Band tests..."
	@$(MAKE) --no-print-directory -C ./band $(MAKECMDGOALS)


# Cleanup
reset:
	@$(MAKE) --no-print-directory -C ./subvector reset
	@$(MAKE) --no-print-directory -C ./elements reset
	@$(MAKE) --no-print-directory -C ./submatrix reset
	@$(MAKE) --no-print-directory -C ./row reset
	@$(MAKE) --no-print-directory -C ./rows reset
	@$(MAKE) --no-print-directory -C ./column reset
	@$(MAKE) --no-print-directory -C ./columns reset
	@$(MAKE) --no-print-directory -C ./band reset

clean:
	@$(MAKE) --no-print-directory -C ./subvector clean
	@$(MAKE) --no-print-directory -C ./elements clean
	@$(MAKE) --no-print-directory -C ./submatrix clean
	@$(MAKE) --no-print-directory -C ./row clean
	@$(MAKE) --no-print-directory -C ./rows clean
	@$(MAKE) --no-print-directory -C ./column clean
	@$(MAKE) --no-print-directory -C ./columns clean
	@$(MAKE) --no-print-directory -C ./band clean


# Setting the independent commands
.PHONY: default all essential single reset clean \
        subvector elements submatrix row rows column columns band
