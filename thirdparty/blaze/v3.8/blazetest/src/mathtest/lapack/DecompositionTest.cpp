//=================================================================================================
/*!
//  \file src/mathtest/lapack/DecompositionTest.cpp
//  \brief Source file for the LAPACK decomposition test
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
#include <blaze/util/Complex.h>
#include <blazetest/mathtest/lapack/DecompositionTest.h>

#ifdef BLAZE_USE_HPX_THREADS
#  include <hpx/hpx_main.hpp>
#endif


namespace blazetest {

namespace mathtest {

namespace lapack {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the DecompositionTest class test.
//
// \exception std::runtime_error Decomposition error detected.
*/
DecompositionTest::DecompositionTest()
{
   using blaze::complex;


   //=====================================================================================
   // Single precision tests
   //=====================================================================================

   //testGetrf< float >();
   //testSytrf< float >();
   //testPotrf< float >();
   //testGeqrf< float >();
   //testOrgqr< float >();
   //testOrg2r< float >();
   //testOrmqr< float >();
   //testGerqf< float >();
   //testOrgrq< float >();
   //testOrgr2< float >();
   //testOrmrq< float >();
   //testGeqlf< float >();
   //testOrgql< float >();
   //testOrg2l< float >();
   //testOrmql< float >();
   //testGelqf< float >();
   //testOrglq< float >();
   //testOrgl2< float >();
   //testOrmlq< float >();


   //=====================================================================================
   // Double precision tests
   //=====================================================================================

   testGetrf< double >();
   testSytrf< double >();
   testPotrf< double >();
   testGeqrf< double >();
   testOrgqr< double >();
   testOrg2r< double >();
   testOrmqr< double >();
   testGerqf< double >();
   testOrgrq< double >();
   testOrgr2< double >();
   testOrmrq< double >();
   testGeqlf< double >();
   testOrgql< double >();
   testOrg2l< double >();
   testOrmql< double >();
   testGelqf< double >();
   testOrglq< double >();
   testOrgl2< double >();
   testOrmlq< double >();


   //=====================================================================================
   // Single precision complex tests
   //=====================================================================================

   //testGetrf< complex<float> >();
   //testSytrf< complex<float> >();
   //testHetrf< complex<float> >();
   //testPotrf< complex<float> >();
   //testGeqrf< complex<float> >();
   //testUngqr< complex<float> >();
   //testUng2r< complex<float> >();
   //testUnmqr< complex<float> >();
   //testGerqf< complex<float> >();
   //testUngrq< complex<float> >();
   //testUngr2< complex<float> >();
   //testUnmrq< complex<float> >();
   //testGeqlf< complex<float> >();
   //testUngql< complex<float> >();
   //testUng2l< complex<float> >();
   //testUnmql< complex<float> >();
   //testGelqf< complex<float> >();
   //testUnglq< complex<float> >();
   //testUngl2< complex<float> >();
   //testUnmlq< complex<float> >();


   //=====================================================================================
   // Double precision complex tests
   //=====================================================================================

   testGetrf< complex<double> >();
   testSytrf< complex<double> >();
   testHetrf< complex<double> >();
   testPotrf< complex<double> >();
   testGeqrf< complex<double> >();
   testUngqr< complex<double> >();
   testUng2r< complex<double> >();
   testUnmqr< complex<double> >();
   testGerqf< complex<double> >();
   testUngrq< complex<double> >();
   testUngr2< complex<double> >();
   testUnmrq< complex<double> >();
   testGeqlf< complex<double> >();
   testUngql< complex<double> >();
   testUng2l< complex<double> >();
   testUnmql< complex<double> >();
   testGelqf< complex<double> >();
   testUnglq< complex<double> >();
   testUngl2< complex<double> >();
   testUnmlq< complex<double> >();
}
//*************************************************************************************************

} // namespace lapack

} // namespace mathtest

} // namespace blazetest




//=================================================================================================
//
//  MAIN FUNCTION
//
//=================================================================================================

//*************************************************************************************************
int main()
{
   std::cout << "   Running LAPACK decomposition test..." << std::endl;

   try
   {
      RUN_LAPACK_DECOMPOSITION_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during LAPACK decomposition test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
