//=================================================================================================
/*!
//  \file src/mathtest/operations/smatdmatmult/HCaUDb.cpp
//  \brief Source file for the HCaUDb sparse matrix/dense matrix multiplication math test
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
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/HermitianMatrix.h>
#include <blaze/math/UpperMatrix.h>
#include <blazetest/mathtest/Creator.h>
#include <blazetest/mathtest/operations/smatdmatmult/OperationTest.h>
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
   std::cout << "   Running 'HCaUDb'..." << std::endl;

   using blazetest::mathtest::ScalarA;
   using blazetest::mathtest::ScalarB;

   try
   {
      // Matrix type definitions
      using HCa = blaze::HermitianMatrix< blaze::CompressedMatrix<ScalarA> >;
      using UDb = blaze::UpperMatrix< blaze::DynamicMatrix<ScalarB> >;

      // Creator type definitions
      using CHCa = blazetest::Creator<HCa>;
      using CUDb = blazetest::Creator<UDb>;

      // Running tests with small matrices
      for( size_t i=0UL; i<=6UL; ++i ) {
         for( size_t j=0UL; j<=i*i; ++j ) {
            RUN_SMATDMATMULT_OPERATION_TEST( CHCa( i, j ), CUDb( i ) );
         }
      }

      // Running tests with large matrices
      RUN_SMATDMATMULT_OPERATION_TEST( CHCa(  31UL,  7UL ), CUDb(  31UL ) );
      RUN_SMATDMATMULT_OPERATION_TEST( CHCa(  67UL,  7UL ), CUDb(  67UL ) );
      RUN_SMATDMATMULT_OPERATION_TEST( CHCa( 127UL, 13UL ), CUDb( 127UL ) );
      RUN_SMATDMATMULT_OPERATION_TEST( CHCa(  32UL,  8UL ), CUDb(  32UL ) );
      RUN_SMATDMATMULT_OPERATION_TEST( CHCa(  64UL,  8UL ), CUDb(  64UL ) );
      RUN_SMATDMATMULT_OPERATION_TEST( CHCa( 128UL, 16UL ), CUDb( 128UL ) );
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during sparse matrix/dense matrix multiplication:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
