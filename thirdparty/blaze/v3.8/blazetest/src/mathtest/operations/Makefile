#==================================================================================================
#
#  Makefile for the operations module of the Blaze test suite
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

all: dvecdvecadd dvecsvecadd svecdvecadd svecsvecadd \
     dvecdvecsub dvecsvecsub svecdvecsub svecsvecsub \
     dvecdvecmult dvecsvecmult svecdvecmult svecsvecmult \
     dvecdveccross dvecsveccross svecdveccross svecsveccross \
     dvecdvecinner dvecsvecinner svecdvecinner svecsvecinner \
     dvecdvecouter dvecsvecouter svecdvecouter svecsvecouter \
     dvecdveckron dvecsveckron svecdveckron svecsveckron \
     dvecdvecdiv svecdvecdiv \
     dvecdvecmin dvecdvecmax \
     dvecreduce svecreduce \
     dvecexpand svecexpand \
     dvecrepeat svecrepeat \
     dmatdvecmult dmatsvecmult smatdvecmult smatsvecmult \
     tdvecdmatmult tsvecdmatmult tdvecsmatmult tsvecsmatmult \
     dmatdmatadd dmatsmatadd smatdmatadd smatsmatadd \
     dmatdmatsub dmatsmatsub smatdmatsub smatsmatsub \
     dmatdmatschur dmatsmatschur smatdmatschur smatsmatschur \
     dmatdmatmult dmatsmatmult smatdmatmult smatsmatmult \
     dmatdmatkron dmatsmatkron smatdmatkron smatsmatkron \
     dmatdmatmin dmatdmatmax \
     dmatreduce smatreduce \
     dmatrepeat smatrepeat \
     determinant lu llh qr rq ql lq inversion lse exponential eigen svd

essential: all

single: all


# Internal rules
dvecdvecadd:
	@echo
	@echo "Building the dense vector/dense vector addition tests..."
	@$(MAKE) --no-print-directory -C ./dvecdvecadd $(MAKECMDGOALS)

dvecsvecadd:
	@echo
	@echo "Building the dense vector/sparse vector addition tests..."
	@$(MAKE) --no-print-directory -C ./dvecsvecadd $(MAKECMDGOALS)

svecdvecadd:
	@echo
	@echo "Building the sparse vector/dense vector addition tests..."
	@$(MAKE) --no-print-directory -C ./svecdvecadd $(MAKECMDGOALS)

svecsvecadd:
	@echo
	@echo "Building the sparse vector/sparse vector addition tests..."
	@$(MAKE) --no-print-directory -C ./svecsvecadd $(MAKECMDGOALS)

dvecdvecsub:
	@echo
	@echo "Building the dense vector/dense vector subtraction tests..."
	@$(MAKE) --no-print-directory -C ./dvecdvecsub $(MAKECMDGOALS)

dvecsvecsub:
	@echo
	@echo "Building the dense vector/sparse vector subtraction tests..."
	@$(MAKE) --no-print-directory -C ./dvecsvecsub $(MAKECMDGOALS)

svecdvecsub:
	@echo
	@echo "Building the sparse vector/dense vector subtraction tests..."
	@$(MAKE) --no-print-directory -C ./svecdvecsub $(MAKECMDGOALS)

svecsvecsub:
	@echo
	@echo "Building the sparse vector/sparse vector subtraction tests..."
	@$(MAKE) --no-print-directory -C ./svecsvecsub $(MAKECMDGOALS)

dvecdvecmult:
	@echo
	@echo "Building the dense vector/dense vector multiplication tests..."
	@$(MAKE) --no-print-directory -C ./dvecdvecmult $(MAKECMDGOALS)

dvecsvecmult:
	@echo
	@echo "Building the dense vector/sparse vector multiplication tests..."
	@$(MAKE) --no-print-directory -C ./dvecsvecmult $(MAKECMDGOALS)

svecdvecmult:
	@echo
	@echo "Building the sparse vector/dense vector multiplication tests..."
	@$(MAKE) --no-print-directory -C ./svecdvecmult $(MAKECMDGOALS)

svecsvecmult:
	@echo
	@echo "Building the sparse vector/sparse vector multiplication tests..."
	@$(MAKE) --no-print-directory -C ./svecsvecmult $(MAKECMDGOALS)

dvecdveccross:
	@echo
	@echo "Building the dense vector/dense vector cross product tests..."
	@$(MAKE) --no-print-directory -C ./dvecdveccross $(MAKECMDGOALS)

dvecsveccross:
	@echo
	@echo "Building the dense vector/sparse vector cross product tests..."
	@$(MAKE) --no-print-directory -C ./dvecsveccross $(MAKECMDGOALS)

svecdveccross:
	@echo
	@echo "Building the sparse vector/dense vector cross product tests..."
	@$(MAKE) --no-print-directory -C ./svecdveccross $(MAKECMDGOALS)

svecsveccross:
	@echo
	@echo "Building the sparse vector/sparse vector cross product tests..."
	@$(MAKE) --no-print-directory -C ./svecsveccross $(MAKECMDGOALS)

dvecdvecinner:
	@echo
	@echo "Building the dense vector/dense vector inner product tests..."
	@$(MAKE) --no-print-directory -C ./dvecdvecinner $(MAKECMDGOALS)

dvecsvecinner:
	@echo
	@echo "Building the dense vector/sparse vector inner product tests..."
	@$(MAKE) --no-print-directory -C ./dvecsvecinner $(MAKECMDGOALS)

svecdvecinner:
	@echo
	@echo "Building the sparse vector/dense vector inner product tests..."
	@$(MAKE) --no-print-directory -C ./svecdvecinner $(MAKECMDGOALS)

svecsvecinner:
	@echo
	@echo "Building the sparse vector/sparse vector inner product tests..."
	@$(MAKE) --no-print-directory -C ./svecsvecinner $(MAKECMDGOALS)

dvecdvecouter:
	@echo
	@echo "Building the dense vector/dense vector outer product tests..."
	@$(MAKE) --no-print-directory -C ./dvecdvecouter $(MAKECMDGOALS)

dvecsvecouter:
	@echo
	@echo "Building the dense vector/sparse vector outer product tests..."
	@$(MAKE) --no-print-directory -C ./dvecsvecouter $(MAKECMDGOALS)

svecdvecouter:
	@echo
	@echo "Building the sparse vector/dense vector outer product tests..."
	@$(MAKE) --no-print-directory -C ./svecdvecouter $(MAKECMDGOALS)

svecsvecouter:
	@echo
	@echo "Building the sparse vector/sparse vector outer product tests..."
	@$(MAKE) --no-print-directory -C ./svecsvecouter $(MAKECMDGOALS)

dvecdveckron:
	@echo
	@echo "Building the dense vector/dense vector Kronecker product tests..."
	@$(MAKE) --no-print-directory -C ./dvecdveckron $(MAKECMDGOALS)

dvecsveckron:
	@echo
	@echo "Building the dense vector/sparse vector Kronecker product tests..."
	@$(MAKE) --no-print-directory -C ./dvecsveckron $(MAKECMDGOALS)

svecdveckron:
	@echo
	@echo "Building the sparse vector/dense vector Kronecker product tests..."
	@$(MAKE) --no-print-directory -C ./svecdveckron $(MAKECMDGOALS)

svecsveckron:
	@echo
	@echo "Building the sparse vector/sparse vector Kronecker product tests..."
	@$(MAKE) --no-print-directory -C ./svecsveckron $(MAKECMDGOALS)

dvecdvecdiv:
	@echo
	@echo "Building the dense vector/dense vector division tests..."
	@$(MAKE) --no-print-directory -C ./dvecdvecdiv $(MAKECMDGOALS)

svecdvecdiv:
	@echo
	@echo "Building the sparse vector/dense vector division tests..."
	@$(MAKE) --no-print-directory -C ./svecdvecdiv $(MAKECMDGOALS)

dvecdvecmin:
	@echo
	@echo "Building the dense vector/dense vector minimum tests..."
	@$(MAKE) --no-print-directory -C ./dvecdvecmin $(MAKECMDGOALS)

dvecdvecmax:
	@echo
	@echo "Building the dense vector/dense vector maximum tests..."
	@$(MAKE) --no-print-directory -C ./dvecdvecmax $(MAKECMDGOALS)

dvecreduce:
	@echo
	@echo "Building the dense vector reduction tests..."
	@$(MAKE) --no-print-directory -C ./dvecreduce $(MAKECMDGOALS)

svecreduce:
	@echo
	@echo "Building the sparse vector reduction tests..."
	@$(MAKE) --no-print-directory -C ./svecreduce $(MAKECMDGOALS)

dvecexpand:
	@echo
	@echo "Building the dense vector expansion tests..."
	@$(MAKE) --no-print-directory -C ./dvecexpand $(MAKECMDGOALS)

svecexpand:
	@echo
	@echo "Building the sparse vector expansion tests..."
	@$(MAKE) --no-print-directory -C ./svecexpand $(MAKECMDGOALS)

dvecrepeat:
	@echo
	@echo "Building the dense vector repeat tests..."
	@$(MAKE) --no-print-directory -C ./dvecrepeat $(MAKECMDGOALS)

svecrepeat:
	@echo
	@echo "Building the sparse vector repeat tests..."
	@$(MAKE) --no-print-directory -C ./svecrepeat $(MAKECMDGOALS)

dmatdvecmult:
	@echo
	@echo "Building the dense matrix/dense vector multiplication tests..."
	@$(MAKE) --no-print-directory -C ./dmatdvecmult $(MAKECMDGOALS)

dmatsvecmult:
	@echo
	@echo "Building the dense matrix/sparse vector multiplication tests..."
	@$(MAKE) --no-print-directory -C ./dmatsvecmult $(MAKECMDGOALS)

smatdvecmult:
	@echo
	@echo "Building the sparse matrix/dense vector multiplication tests..."
	@$(MAKE) --no-print-directory -C ./smatdvecmult $(MAKECMDGOALS)

smatsvecmult:
	@echo
	@echo "Building the sparse matrix/sparse vector multiplication tests..."
	@$(MAKE) --no-print-directory -C ./smatsvecmult $(MAKECMDGOALS)

tdvecdmatmult:
	@echo
	@echo "Building the dense vector/dense matrix multiplication tests..."
	@$(MAKE) --no-print-directory -C ./tdvecdmatmult $(MAKECMDGOALS)

tsvecdmatmult:
	@echo
	@echo "Building the sparse vector/dense matrix multiplication tests..."
	@$(MAKE) --no-print-directory -C ./tsvecdmatmult $(MAKECMDGOALS)

tdvecsmatmult:
	@echo
	@echo "Building the dense vector/sparse matrix multiplication tests..."
	@$(MAKE) --no-print-directory -C ./tdvecsmatmult $(MAKECMDGOALS)

tsvecsmatmult:
	@echo
	@echo "Building the sparse vector/sparse matrix multiplication tests..."
	@$(MAKE) --no-print-directory -C ./tsvecsmatmult $(MAKECMDGOALS)

dmatdmatadd:
	@echo
	@echo "Building the dense matrix/dense matrix addition tests..."
	@$(MAKE) --no-print-directory -C ./dmatdmatadd $(MAKECMDGOALS)

dmatsmatadd:
	@echo
	@echo "Building the dense matrix/sparse matrix addition tests..."
	@$(MAKE) --no-print-directory -C ./dmatsmatadd $(MAKECMDGOALS)

smatdmatadd:
	@echo
	@echo "Building the sparse matrix/dense matrix addition tests..."
	@$(MAKE) --no-print-directory -C ./smatdmatadd $(MAKECMDGOALS)

smatsmatadd:
	@echo
	@echo "Building the sparse matrix/sparse matrix addition tests..."
	@$(MAKE) --no-print-directory -C ./smatsmatadd $(MAKECMDGOALS)

dmatdmatsub:
	@echo
	@echo "Building the dense matrix/dense matrix subtraction tests..."
	@$(MAKE) --no-print-directory -C ./dmatdmatsub $(MAKECMDGOALS)

dmatsmatsub:
	@echo
	@echo "Building the dense matrix/sparse matrix subtraction tests..."
	@$(MAKE) --no-print-directory -C ./dmatsmatsub $(MAKECMDGOALS)

smatdmatsub:
	@echo
	@echo "Building the sparse matrix/dense matrix subtraction tests..."
	@$(MAKE) --no-print-directory -C ./smatdmatsub $(MAKECMDGOALS)

smatsmatsub:
	@echo
	@echo "Building the sparse matrix/sparse matrix subtraction tests..."
	@$(MAKE) --no-print-directory -C ./smatsmatsub $(MAKECMDGOALS)

dmatdmatschur:
	@echo
	@echo "Building the dense matrix/dense matrix Schur product tests..."
	@$(MAKE) --no-print-directory -C ./dmatdmatschur $(MAKECMDGOALS)

dmatsmatschur:
	@echo
	@echo "Building the dense matrix/sparse matrix Schur product tests..."
	@$(MAKE) --no-print-directory -C ./dmatsmatschur $(MAKECMDGOALS)

smatdmatschur:
	@echo
	@echo "Building the sparse matrix/dense matrix Schur product tests..."
	@$(MAKE) --no-print-directory -C ./smatdmatschur $(MAKECMDGOALS)

smatsmatschur:
	@echo
	@echo "Building the sparse matrix/sparse matrix Schur product tests..."
	@$(MAKE) --no-print-directory -C ./smatsmatschur $(MAKECMDGOALS)

dmatdmatmult:
	@echo
	@echo "Building the dense matrix/dense matrix multiplication tests..."
	@$(MAKE) --no-print-directory -C ./dmatdmatmult $(MAKECMDGOALS)

dmatsmatmult:
	@echo
	@echo "Building the dense matrix/sparse matrix multiplication tests..."
	@$(MAKE) --no-print-directory -C ./dmatsmatmult $(MAKECMDGOALS)

smatdmatmult:
	@echo
	@echo "Building the sparse matrix/dense matrix multiplication tests..."
	@$(MAKE) --no-print-directory -C ./smatdmatmult $(MAKECMDGOALS)

smatsmatmult:
	@echo
	@echo "Building the sparse matrix/sparse matrix multiplication tests..."
	@$(MAKE) --no-print-directory -C ./smatsmatmult $(MAKECMDGOALS)

dmatdmatkron:
	@echo
	@echo "Building the dense matrix/dense matrix Kronecker product tests..."
	@$(MAKE) --no-print-directory -C ./dmatdmatkron $(MAKECMDGOALS)

dmatsmatkron:
	@echo
	@echo "Building the dense matrix/sparse matrix Kronecker product tests..."
	@$(MAKE) --no-print-directory -C ./dmatsmatkron $(MAKECMDGOALS)

smatdmatkron:
	@echo
	@echo "Building the sparse matrix/dense matrix Kronecker product tests..."
	@$(MAKE) --no-print-directory -C ./smatdmatkron $(MAKECMDGOALS)

smatsmatkron:
	@echo
	@echo "Building the sparse matrix/sparse matrix Kronecker product tests..."
	@$(MAKE) --no-print-directory -C ./smatsmatkron $(MAKECMDGOALS)

dmatdmatmin:
	@echo
	@echo "Building the dense matrix/dense matrix minimum tests..."
	@$(MAKE) --no-print-directory -C ./dmatdmatmin $(MAKECMDGOALS)

dmatdmatmax:
	@echo
	@echo "Building the dense matrix/dense matrix maximum tests..."
	@$(MAKE) --no-print-directory -C ./dmatdmatmax $(MAKECMDGOALS)

dmatreduce:
	@$(MAKE) --no-print-directory -C ./dmatreduce $(MAKECMDGOALS)

smatreduce:
	@$(MAKE) --no-print-directory -C ./smatreduce $(MAKECMDGOALS)

dmatrepeat:
	@echo
	@echo "Building the dense matrix repeat tests..."
	@$(MAKE) --no-print-directory -C ./dmatrepeat $(MAKECMDGOALS)

smatrepeat:
	@echo
	@echo "Building the sparse matrix repeat tests..."
	@$(MAKE) --no-print-directory -C ./smatrepeat $(MAKECMDGOALS)

determinant:
	@echo
	@echo "Building the matrix determinant tests..."
	@$(MAKE) --no-print-directory -C ./determinant $(MAKECMDGOALS)

lu:
	@echo
	@echo "Building the matrix LU decomposition tests..."
	@$(MAKE) --no-print-directory -C ./lu $(MAKECMDGOALS)

llh:
	@echo
	@echo "Building the matrix LLH decomposition tests..."
	@$(MAKE) --no-print-directory -C ./llh $(MAKECMDGOALS)

qr:
	@echo
	@echo "Building the matrix QR decomposition tests..."
	@$(MAKE) --no-print-directory -C ./qr $(MAKECMDGOALS)

rq:
	@echo
	@echo "Building the matrix RQ decomposition tests..."
	@$(MAKE) --no-print-directory -C ./rq $(MAKECMDGOALS)

ql:
	@echo
	@echo "Building the matrix QL decomposition tests..."
	@$(MAKE) --no-print-directory -C ./ql $(MAKECMDGOALS)

lq:
	@echo
	@echo "Building the matrix LQ decomposition tests..."
	@$(MAKE) --no-print-directory -C ./lq $(MAKECMDGOALS)

inversion:
	@echo
	@echo "Building the matrix inversion tests..."
	@$(MAKE) --no-print-directory -C ./inversion $(MAKECMDGOALS)

lse:
	@echo
	@echo "Building the matrix LSE tests..."
	@$(MAKE) --no-print-directory -C ./lse $(MAKECMDGOALS)

exponential:
	@echo
	@echo "Building the matrix exponential tests..."
	@$(MAKE) --no-print-directory -C ./exponential $(MAKECMDGOALS)

eigen:
	@echo
	@echo "Building the matrix eigenvalue/eigenvector tests..."
	@$(MAKE) --no-print-directory -C ./eigen $(MAKECMDGOALS)

svd:
	@echo
	@echo "Building the matrix singular value/vector tests..."
	@$(MAKE) --no-print-directory -C ./svd $(MAKECMDGOALS)


# Cleanup
reset:
	@$(MAKE) --no-print-directory -C ./dvecdvecadd reset
	@$(MAKE) --no-print-directory -C ./dvecsvecadd reset
	@$(MAKE) --no-print-directory -C ./svecdvecadd reset
	@$(MAKE) --no-print-directory -C ./svecsvecadd reset
	@$(MAKE) --no-print-directory -C ./dvecdvecsub reset
	@$(MAKE) --no-print-directory -C ./dvecsvecsub reset
	@$(MAKE) --no-print-directory -C ./svecdvecsub reset
	@$(MAKE) --no-print-directory -C ./svecsvecsub reset
	@$(MAKE) --no-print-directory -C ./dvecdvecmult reset
	@$(MAKE) --no-print-directory -C ./dvecsvecmult reset
	@$(MAKE) --no-print-directory -C ./svecdvecmult reset
	@$(MAKE) --no-print-directory -C ./svecsvecmult reset
	@$(MAKE) --no-print-directory -C ./dvecdveccross reset
	@$(MAKE) --no-print-directory -C ./dvecsveccross reset
	@$(MAKE) --no-print-directory -C ./svecdveccross reset
	@$(MAKE) --no-print-directory -C ./svecsveccross reset
	@$(MAKE) --no-print-directory -C ./dvecdvecinner reset
	@$(MAKE) --no-print-directory -C ./dvecsvecinner reset
	@$(MAKE) --no-print-directory -C ./svecdvecinner reset
	@$(MAKE) --no-print-directory -C ./svecsvecinner reset
	@$(MAKE) --no-print-directory -C ./dvecdvecouter reset
	@$(MAKE) --no-print-directory -C ./dvecsvecouter reset
	@$(MAKE) --no-print-directory -C ./svecdvecouter reset
	@$(MAKE) --no-print-directory -C ./svecsvecouter reset
	@$(MAKE) --no-print-directory -C ./dvecdveckron reset
	@$(MAKE) --no-print-directory -C ./dvecsveckron reset
	@$(MAKE) --no-print-directory -C ./svecdveckron reset
	@$(MAKE) --no-print-directory -C ./svecsveckron reset
	@$(MAKE) --no-print-directory -C ./dvecdvecdiv reset
	@$(MAKE) --no-print-directory -C ./svecdvecdiv reset
	@$(MAKE) --no-print-directory -C ./dvecdvecmin reset
	@$(MAKE) --no-print-directory -C ./dvecdvecmax reset
	@$(MAKE) --no-print-directory -C ./dvecreduce reset
	@$(MAKE) --no-print-directory -C ./svecreduce reset
	@$(MAKE) --no-print-directory -C ./dvecexpand reset
	@$(MAKE) --no-print-directory -C ./svecexpand reset
	@$(MAKE) --no-print-directory -C ./dvecrepeat reset
	@$(MAKE) --no-print-directory -C ./svecrepeat reset
	@$(MAKE) --no-print-directory -C ./dmatdvecmult reset
	@$(MAKE) --no-print-directory -C ./dmatsvecmult reset
	@$(MAKE) --no-print-directory -C ./smatdvecmult reset
	@$(MAKE) --no-print-directory -C ./smatsvecmult reset
	@$(MAKE) --no-print-directory -C ./tdvecdmatmult reset
	@$(MAKE) --no-print-directory -C ./tsvecdmatmult reset
	@$(MAKE) --no-print-directory -C ./tdvecsmatmult reset
	@$(MAKE) --no-print-directory -C ./tsvecsmatmult reset
	@$(MAKE) --no-print-directory -C ./dmatdmatadd reset
	@$(MAKE) --no-print-directory -C ./dmatsmatadd reset
	@$(MAKE) --no-print-directory -C ./smatdmatadd reset
	@$(MAKE) --no-print-directory -C ./smatsmatadd reset
	@$(MAKE) --no-print-directory -C ./dmatdmatsub reset
	@$(MAKE) --no-print-directory -C ./dmatsmatsub reset
	@$(MAKE) --no-print-directory -C ./smatdmatsub reset
	@$(MAKE) --no-print-directory -C ./smatsmatsub reset
	@$(MAKE) --no-print-directory -C ./dmatdmatschur reset
	@$(MAKE) --no-print-directory -C ./dmatsmatschur reset
	@$(MAKE) --no-print-directory -C ./smatdmatschur reset
	@$(MAKE) --no-print-directory -C ./smatsmatschur reset
	@$(MAKE) --no-print-directory -C ./dmatdmatmult reset
	@$(MAKE) --no-print-directory -C ./dmatsmatmult reset
	@$(MAKE) --no-print-directory -C ./smatdmatmult reset
	@$(MAKE) --no-print-directory -C ./smatsmatmult reset
	@$(MAKE) --no-print-directory -C ./dmatdmatkron reset
	@$(MAKE) --no-print-directory -C ./dmatsmatkron reset
	@$(MAKE) --no-print-directory -C ./smatdmatkron reset
	@$(MAKE) --no-print-directory -C ./smatsmatkron reset
	@$(MAKE) --no-print-directory -C ./dmatdmatmin reset
	@$(MAKE) --no-print-directory -C ./dmatdmatmax reset
	@$(MAKE) --no-print-directory -C ./dmatreduce reset
	@$(MAKE) --no-print-directory -C ./smatreduce reset
	@$(MAKE) --no-print-directory -C ./dmatrepeat reset
	@$(MAKE) --no-print-directory -C ./smatrepeat reset
	@$(MAKE) --no-print-directory -C ./determinant reset
	@$(MAKE) --no-print-directory -C ./lu reset
	@$(MAKE) --no-print-directory -C ./llh reset
	@$(MAKE) --no-print-directory -C ./qr reset
	@$(MAKE) --no-print-directory -C ./rq reset
	@$(MAKE) --no-print-directory -C ./ql reset
	@$(MAKE) --no-print-directory -C ./lq reset
	@$(MAKE) --no-print-directory -C ./inversion reset
	@$(MAKE) --no-print-directory -C ./lse reset
	@$(MAKE) --no-print-directory -C ./exponential reset
	@$(MAKE) --no-print-directory -C ./eigen reset
	@$(MAKE) --no-print-directory -C ./svd reset

clean:
	@$(MAKE) --no-print-directory -C ./dvecdvecadd clean
	@$(MAKE) --no-print-directory -C ./dvecsvecadd clean
	@$(MAKE) --no-print-directory -C ./svecdvecadd clean
	@$(MAKE) --no-print-directory -C ./svecsvecadd clean
	@$(MAKE) --no-print-directory -C ./dvecdvecsub clean
	@$(MAKE) --no-print-directory -C ./dvecsvecsub clean
	@$(MAKE) --no-print-directory -C ./svecdvecsub clean
	@$(MAKE) --no-print-directory -C ./svecsvecsub clean
	@$(MAKE) --no-print-directory -C ./dvecdvecmult clean
	@$(MAKE) --no-print-directory -C ./dvecsvecmult clean
	@$(MAKE) --no-print-directory -C ./svecdvecmult clean
	@$(MAKE) --no-print-directory -C ./svecsvecmult clean
	@$(MAKE) --no-print-directory -C ./dvecdveccross clean
	@$(MAKE) --no-print-directory -C ./dvecsveccross clean
	@$(MAKE) --no-print-directory -C ./svecdveccross clean
	@$(MAKE) --no-print-directory -C ./svecsveccross clean
	@$(MAKE) --no-print-directory -C ./dvecdvecinner clean
	@$(MAKE) --no-print-directory -C ./dvecsvecinner clean
	@$(MAKE) --no-print-directory -C ./svecdvecinner clean
	@$(MAKE) --no-print-directory -C ./svecsvecinner clean
	@$(MAKE) --no-print-directory -C ./dvecdvecouter clean
	@$(MAKE) --no-print-directory -C ./dvecsvecouter clean
	@$(MAKE) --no-print-directory -C ./svecdvecouter clean
	@$(MAKE) --no-print-directory -C ./svecsvecouter clean
	@$(MAKE) --no-print-directory -C ./dvecdveckron clean
	@$(MAKE) --no-print-directory -C ./dvecsveckron clean
	@$(MAKE) --no-print-directory -C ./svecdveckron clean
	@$(MAKE) --no-print-directory -C ./svecsveckron clean
	@$(MAKE) --no-print-directory -C ./dvecdvecdiv clean
	@$(MAKE) --no-print-directory -C ./svecdvecdiv clean
	@$(MAKE) --no-print-directory -C ./dvecdvecmin clean
	@$(MAKE) --no-print-directory -C ./dvecdvecmax clean
	@$(MAKE) --no-print-directory -C ./dvecreduce clean
	@$(MAKE) --no-print-directory -C ./svecreduce clean
	@$(MAKE) --no-print-directory -C ./dvecexpand clean
	@$(MAKE) --no-print-directory -C ./svecexpand clean
	@$(MAKE) --no-print-directory -C ./dvecrepeat clean
	@$(MAKE) --no-print-directory -C ./svecrepeat clean
	@$(MAKE) --no-print-directory -C ./dmatdvecmult clean
	@$(MAKE) --no-print-directory -C ./dmatsvecmult clean
	@$(MAKE) --no-print-directory -C ./smatdvecmult clean
	@$(MAKE) --no-print-directory -C ./smatsvecmult clean
	@$(MAKE) --no-print-directory -C ./tdvecdmatmult clean
	@$(MAKE) --no-print-directory -C ./tsvecdmatmult clean
	@$(MAKE) --no-print-directory -C ./tdvecsmatmult clean
	@$(MAKE) --no-print-directory -C ./tsvecsmatmult clean
	@$(MAKE) --no-print-directory -C ./dmatdmatadd clean
	@$(MAKE) --no-print-directory -C ./dmatsmatadd clean
	@$(MAKE) --no-print-directory -C ./smatdmatadd clean
	@$(MAKE) --no-print-directory -C ./smatsmatadd clean
	@$(MAKE) --no-print-directory -C ./dmatdmatsub clean
	@$(MAKE) --no-print-directory -C ./dmatsmatsub clean
	@$(MAKE) --no-print-directory -C ./smatdmatsub clean
	@$(MAKE) --no-print-directory -C ./smatsmatsub clean
	@$(MAKE) --no-print-directory -C ./dmatdmatschur clean
	@$(MAKE) --no-print-directory -C ./dmatsmatschur clean
	@$(MAKE) --no-print-directory -C ./smatdmatschur clean
	@$(MAKE) --no-print-directory -C ./smatsmatschur clean
	@$(MAKE) --no-print-directory -C ./dmatdmatmult clean
	@$(MAKE) --no-print-directory -C ./dmatsmatmult clean
	@$(MAKE) --no-print-directory -C ./smatdmatmult clean
	@$(MAKE) --no-print-directory -C ./smatsmatmult clean
	@$(MAKE) --no-print-directory -C ./dmatdmatkron clean
	@$(MAKE) --no-print-directory -C ./dmatsmatkron clean
	@$(MAKE) --no-print-directory -C ./smatdmatkron clean
	@$(MAKE) --no-print-directory -C ./smatsmatkron clean
	@$(MAKE) --no-print-directory -C ./dmatdmatmin clean
	@$(MAKE) --no-print-directory -C ./dmatdmatmax clean
	@$(MAKE) --no-print-directory -C ./dmatreduce clean
	@$(MAKE) --no-print-directory -C ./smatreduce clean
	@$(MAKE) --no-print-directory -C ./dmatrepeat clean
	@$(MAKE) --no-print-directory -C ./smatrepeat clean
	@$(MAKE) --no-print-directory -C ./determinant clean
	@$(MAKE) --no-print-directory -C ./lu clean
	@$(MAKE) --no-print-directory -C ./llh clean
	@$(MAKE) --no-print-directory -C ./qr clean
	@$(MAKE) --no-print-directory -C ./rq clean
	@$(MAKE) --no-print-directory -C ./ql clean
	@$(MAKE) --no-print-directory -C ./lq clean
	@$(MAKE) --no-print-directory -C ./inversion clean
	@$(MAKE) --no-print-directory -C ./lse clean
	@$(MAKE) --no-print-directory -C ./exponential clean
	@$(MAKE) --no-print-directory -C ./eigen clean
	@$(MAKE) --no-print-directory -C ./svd clean


# Setting the independent commands
.PHONY: default all essential single reset clean \
        dvecdvecadd dvecsvecadd svecdvecadd svecsvecadd \
        dvecdvecsub dvecsvecsub svecdvecsub svecsvecsub \
        dvecdvecmult dvecsvecmult svecdvecmult svecsvecmult \
        dvecdveccross dvecsveccross svecdveccross svecsveccross \
        dvecdvecinner dvecsvecinner svecdvecinner svecsvecinner \
        dvecdvecouter dvecsvecouter svecdvecouter svecsvecouter \
        dvecdveckron dvecsveckron svecdveckron svecsveckron \
        dvecdvecdiv svecdvecdiv \
        dvecdvecmin dvecdvecmax \
        dvecreduce svecreduce \
        dvecexpand svecexpand \
        dvecrepeat svecrepeat \
        dmatdvecmult dmatsvecmult smatdvecmult smatsvecmult \
        tdvecdmatmult tsvecdmatmult tdvecsmatmult tsvecsmatmult \
        dmatdmatadd dmatsmatadd smatdmatadd smatsmatadd \
        dmatdmatsub dmatsmatsub smatdmatsub smatsmatsub \
        dmatdmatschur dmatsmatschur smatdmatschur smatsmatschur \
        dmatdmatmult dmatsmatmult smatdmatmult smatsmatmult \
        dmatdmatkron dmatsmatkron smatdmatkron smatsmatkron \
        dmatdmatmin dmatdmatmax \
        dmatreduce smatreduce \
        dmatrepeat smatrepeat \
        determinant lu llh qr rq ql lq inversion lse exponential eigen svd
