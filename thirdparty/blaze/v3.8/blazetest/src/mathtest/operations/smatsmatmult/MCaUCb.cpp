//=================================================================================================
/*!
//  \file src/mathtest/operations/smatsmatmult/MCaUCb.cpp
//  \brief Source file for the MCaUCb sparse matrix/sparse matrix multiplication math test
//
//  Copyright (C) 2012-2020 Klaus Iglberger - All Rights Reserved
//
//  This file is part of the Blaze library. You can redistribute it and/or modify it under
//  the terms of the New (Revised) BSD License. Redistribution and use in source and binary
//  forms, with or without modification, are permitted provided that the following conditions
//  are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this list of
//     conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice, this list
//     of conditions and the following disclaimer in the documentation and/or other materials
//     provided with the distribution.
//  3. Neither the names of the Blaze development group nor the names of its contributors
//     may be used to endorse or promote products derived from this software without specific
//     prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
//  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
//  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
//  SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
//  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
//  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
//  BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
//  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
//  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
//  DAMAGE.
*/
//=================================================================================================


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <cstdlib>
#include <iostream>
#include <blaze/math/CompressedMatrix.h>
#include <blaze/math/UpperMatrix.h>
#include <blazetest/mathtest/Creator.h>
#include <blazetest/mathtest/operations/smatsmatmult/OperationTest.h>
#include <blazetest/system/MathTest.h>

#ifdef BLAZE_USE_HPX_THREADS
#  include <hpx/hpx_main.hpp>
#endif


//=================================================================================================
//
//  MAIN FUNCTION
//
//=================================================================================================

//*************************************************************************************************
int main()
{
   std::cout << "   Running 'MCaUCb'..." << std::endl;

   using blazetest::mathtest::TypeA;
   using blazetest::mathtest::TypeB;

   try
   {
      // Matrix type definitions
      using MCa = blaze::CompressedMatrix<TypeA>;
      using UCb = blaze::UpperMatrix< blaze::CompressedMatrix<TypeB> >;

      // Creator type definitions
      using CMCa = blazetest::Creator<MCa>;
      using CUCb = blazetest::Creator<UCb>;

      // Running tests with small matrices
      for( size_t i=0UL; i<=6UL; ++i ) {
         for( size_t j=0UL; j<=6UL; ++j ) {
            RUN_SMATSMATMULT_OPERATION_TEST( CMCa( i, j,     0UL ), CUCb( j,     0UL ) );
            RUN_SMATSMATMULT_OPERATION_TEST( CMCa( i, j,     0UL ), CUCb( j, 0.2*j*j ) );
            RUN_SMATSMATMULT_OPERATION_TEST( CMCa( i, j,     0UL ), CUCb( j, 0.5*j*j ) );
            RUN_SMATSMATMULT_OPERATION_TEST( CMCa( i, j, 0.3*i*j ), CUCb( j,     0UL ) );
            RUN_SMATSMATMULT_OPERATION_TEST( CMCa( i, j, 0.3*i*j ), CUCb( j, 0.2*j*j ) );
            RUN_SMATSMATMULT_OPERATION_TEST( CMCa( i, j, 0.3*i*j ), CUCb( j, 0.5*j*j ) );
            RUN_SMATSMATMULT_OPERATION_TEST( CMCa( i, j,     i*j ), CUCb( j,     0UL ) );
            RUN_SMATSMATMULT_OPERATION_TEST( CMCa( i, j,     i*j ), CUCb( j, 0.2*j*j ) );
            RUN_SMATSMATMULT_OPERATION_TEST( CMCa( i, j,     i*j ), CUCb( j, 0.5*j*j ) );
         }
      }

      // Running tests with large matrices
      RUN_SMATSMATMULT_OPERATION_TEST( CMCa( 37UL, 15UL,  7UL ), CUCb( 15UL,  7UL ) );
      RUN_SMATSMATMULT_OPERATION_TEST( CMCa( 37UL, 37UL,  7UL ), CUCb( 37UL,  7UL ) );
      RUN_SMATSMATMULT_OPERATION_TEST( CMCa( 37UL, 63UL, 13UL ), CUCb( 63UL, 13UL ) );
      RUN_SMATSMATMULT_OPERATION_TEST( CMCa( 32UL, 16UL,  8UL ), CUCb( 16UL,  8UL ) );
      RUN_SMATSMATMULT_OPERATION_TEST( CMCa( 32UL, 32UL,  8UL ), CUCb( 32UL,  8UL ) );
      RUN_SMATSMATMULT_OPERATION_TEST( CMCa( 32UL, 64UL, 16UL ), CUCb( 64UL, 16UL ) );
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during sparse matrix/sparse matrix multiplication:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
