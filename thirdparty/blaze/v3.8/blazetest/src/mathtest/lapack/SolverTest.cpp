//=================================================================================================
/*!
//  \file src/mathtest/lapack/SolverTest.cpp
//  \brief Source file for the LAPACK solver test
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
#include <blazetest/mathtest/lapack/SolverTest.h>

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
/*!\brief Constructor for the SolverTest class test.
//
// \exception std::runtime_error Solver error detected.
*/
SolverTest::SolverTest()
{
   using blaze::complex;


   //=====================================================================================
   // Single precision tests
   //=====================================================================================

   //testGesv< float >();
   //testSysv< float >();
   //testPosv< float >();
   //testTrsv< float >();


   //=====================================================================================
   // Double precision tests
   //=====================================================================================

   testGesv< double >();
   testSysv< double >();
   testPosv< double >();
   testTrsv< double >();


   //=====================================================================================
   // Single precision complex tests
   //=====================================================================================

   //testGesv< complex<float> >();
   //testSysv< complex<float> >();
   //testHesv< complex<float> >();
   //testPosv< complex<float> >();
   //testTrsv< complex<float> >();


   //=====================================================================================
   // Double precision complex tests
   //=====================================================================================

   testGesv< complex<double> >();
   testSysv< complex<double> >();
   testHesv< complex<double> >();
   testPosv< complex<double> >();
   testTrsv< complex<double> >();
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
   std::cout << "   Running LAPACK solver test..." << std::endl;

   try
   {
      RUN_LAPACK_SOLVER_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during LAPACK solver test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
