//=================================================================================================
/*!
//  \file src/mathtest/operations/lse/DenseTest.cpp
//  \brief Source file for the dense matrix LSE test
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
#include <blazetest/mathtest/operations/lse/DenseTest.h>

#ifdef BLAZE_USE_HPX_THREADS
#  include <hpx/hpx_main.hpp>
#endif


namespace blazetest {

namespace mathtest {

namespace operations {

namespace lse {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the DenseTest LSE test.
//
// \exception std::runtime_error LSE error detected.
*/
DenseTest::DenseTest()
{
   for( size_t i=0UL; i<=12UL; ++i )
   {
      //testGeneral  < float >( i );
      //testSymmetric< float >( i );
      //testHermitian< float >( i );
      //testLower    < float >( i );
      //testUniLower < float >( i );
      //testUpper    < float >( i );
      //testUniUpper < float >( i );
      //testDiagonal < float >( i );

      testGeneral  < double >( i );
      testSymmetric< double >( i );
      testHermitian< double >( i );
      testLower    < double >( i );
      testUniLower < double >( i );
      testUpper    < double >( i );
      testUniUpper < double >( i );
      testDiagonal < double >( i );

      //testGeneral  < complex<float> >( i );
      //testSymmetric< complex<float> >( i );
      //testHermitian< complex<float> >( i );
      //testLower    < complex<float> >( i );
      //testUniLower < complex<float> >( i );
      //testUpper    < complex<float> >( i );
      //testUniUpper < complex<float> >( i );
      //testDiagonal < complex<float> >( i );

      testGeneral  < complex<double> >( i );
      testSymmetric< complex<double> >( i );
      testHermitian< complex<double> >( i );
      testLower    < complex<double> >( i );
      testUniLower < complex<double> >( i );
      testUpper    < complex<double> >( i );
      testUniUpper < complex<double> >( i );
      testDiagonal < complex<double> >( i );
   }
}
//*************************************************************************************************

} // namespace lse

} // namespace operations

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
   std::cout << "   Running dense matrix LSE test..." << std::endl;

   try
   {
      RUN_LSE_DENSE_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during dense matrix LSE test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
