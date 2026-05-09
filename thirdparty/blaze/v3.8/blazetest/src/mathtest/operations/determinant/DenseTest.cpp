//=================================================================================================
/*!
//  \file src/mathtest/operations/determinant/DenseTest.cpp
//  \brief Source file for the dense matrix determinant test
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
#include <blaze/math/DiagonalMatrix.h>
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/HermitianMatrix.h>
#include <blaze/math/LowerMatrix.h>
#include <blaze/math/SymmetricMatrix.h>
#include <blaze/math/UniLowerMatrix.h>
#include <blaze/math/UniUpperMatrix.h>
#include <blaze/math/UpperMatrix.h>
#include <blaze/util/Complex.h>
#include <blazetest/mathtest/operations/determinant/DenseTest.h>

#ifdef BLAZE_USE_HPX_THREADS
#  include <hpx/hpx_main.hpp>
#endif


namespace blazetest {

namespace mathtest {

namespace operations {

namespace determinant {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the DenseTest determinant test.
//
// \exception std::runtime_error Error during determinant computation error detected.
*/
DenseTest::DenseTest()
{
   using blaze::DynamicMatrix;
   using blaze::SymmetricMatrix;
   using blaze::HermitianMatrix;
   using blaze::LowerMatrix;
   using blaze::UniLowerMatrix;
   using blaze::UpperMatrix;
   using blaze::UniUpperMatrix;
   using blaze::DiagonalMatrix;
   using blaze::rowMajor;
   using blaze::columnMajor;

   using cfloat  = blaze::complex<float>;
   using cdouble = blaze::complex<double>;


   //=====================================================================================
   // Specific matrix tests
   //=====================================================================================

   testSpecific();


   //=====================================================================================
   // Random 2x2 matrix tests
   //=====================================================================================

   testRandom2x2< DynamicMatrix<float  ,rowMajor> >();
   testRandom2x2< DynamicMatrix<double ,rowMajor> >();
   testRandom2x2< DynamicMatrix<cfloat ,rowMajor> >();
   testRandom2x2< DynamicMatrix<cdouble,rowMajor> >();

   testRandom2x2< DynamicMatrix<float  ,columnMajor> >();
   testRandom2x2< DynamicMatrix<double ,columnMajor> >();
   testRandom2x2< DynamicMatrix<cfloat ,columnMajor> >();
   testRandom2x2< DynamicMatrix<cdouble,columnMajor> >();

   testRandom2x2< SymmetricMatrix< DynamicMatrix<double,rowMajor> > >();
   testRandom2x2< HermitianMatrix< DynamicMatrix<double,rowMajor> > >();
   testRandom2x2< LowerMatrix< DynamicMatrix<double,rowMajor>     > >();
   testRandom2x2< UniLowerMatrix< DynamicMatrix<double,rowMajor>  > >();
   testRandom2x2< UpperMatrix< DynamicMatrix<double,rowMajor>     > >();
   testRandom2x2< UniUpperMatrix< DynamicMatrix<double,rowMajor>  > >();
   testRandom2x2< DiagonalMatrix< DynamicMatrix<double,rowMajor>  > >();

   testRandom2x2< SymmetricMatrix< DynamicMatrix<double,columnMajor> > >();
   testRandom2x2< HermitianMatrix< DynamicMatrix<double,columnMajor> > >();
   testRandom2x2< LowerMatrix< DynamicMatrix<double,columnMajor>     > >();
   testRandom2x2< UniLowerMatrix< DynamicMatrix<double,columnMajor>  > >();
   testRandom2x2< UpperMatrix< DynamicMatrix<double,columnMajor>     > >();
   testRandom2x2< UniUpperMatrix< DynamicMatrix<double,columnMajor>  > >();
   testRandom2x2< DiagonalMatrix< DynamicMatrix<double,columnMajor>  > >();


   //=====================================================================================
   // Random 3x3 matrix tests
   //=====================================================================================

   testRandom3x3< DynamicMatrix<float  ,rowMajor> >();
   testRandom3x3< DynamicMatrix<double ,rowMajor> >();
   testRandom3x3< DynamicMatrix<cfloat ,rowMajor> >();
   testRandom3x3< DynamicMatrix<cdouble,rowMajor> >();

   testRandom3x3< DynamicMatrix<float  ,columnMajor> >();
   testRandom3x3< DynamicMatrix<double ,columnMajor> >();
   testRandom3x3< DynamicMatrix<cfloat ,columnMajor> >();
   testRandom3x3< DynamicMatrix<cdouble,columnMajor> >();

   testRandom3x3< SymmetricMatrix< DynamicMatrix<double,rowMajor> > >();
   testRandom3x3< HermitianMatrix< DynamicMatrix<double,rowMajor> > >();
   testRandom3x3< LowerMatrix< DynamicMatrix<double,rowMajor>     > >();
   testRandom3x3< UniLowerMatrix< DynamicMatrix<double,rowMajor>  > >();
   testRandom3x3< UpperMatrix< DynamicMatrix<double,rowMajor>     > >();
   testRandom3x3< UniUpperMatrix< DynamicMatrix<double,rowMajor>  > >();
   testRandom3x3< DiagonalMatrix< DynamicMatrix<double,rowMajor>  > >();

   testRandom3x3< SymmetricMatrix< DynamicMatrix<double,columnMajor> > >();
   testRandom3x3< HermitianMatrix< DynamicMatrix<double,columnMajor> > >();
   testRandom3x3< LowerMatrix< DynamicMatrix<double,columnMajor>     > >();
   testRandom3x3< UniLowerMatrix< DynamicMatrix<double,columnMajor>  > >();
   testRandom3x3< UpperMatrix< DynamicMatrix<double,columnMajor>     > >();
   testRandom3x3< UniUpperMatrix< DynamicMatrix<double,columnMajor>  > >();
   testRandom3x3< DiagonalMatrix< DynamicMatrix<double,columnMajor>  > >();


   //=====================================================================================
   // Random 4x4 matrix tests
   //=====================================================================================

   testRandom4x4< DynamicMatrix<float  ,rowMajor> >();
   testRandom4x4< DynamicMatrix<double ,rowMajor> >();
   testRandom4x4< DynamicMatrix<cfloat ,rowMajor> >();
   testRandom4x4< DynamicMatrix<cdouble,rowMajor> >();

   testRandom4x4< DynamicMatrix<float  ,columnMajor> >();
   testRandom4x4< DynamicMatrix<double ,columnMajor> >();
   testRandom4x4< DynamicMatrix<cfloat ,columnMajor> >();
   testRandom4x4< DynamicMatrix<cdouble,columnMajor> >();

   testRandom4x4< SymmetricMatrix< DynamicMatrix<double,rowMajor> > >();
   testRandom4x4< HermitianMatrix< DynamicMatrix<double,rowMajor> > >();
   testRandom4x4< LowerMatrix< DynamicMatrix<double,rowMajor>     > >();
   testRandom4x4< UniLowerMatrix< DynamicMatrix<double,rowMajor>  > >();
   testRandom4x4< UpperMatrix< DynamicMatrix<double,rowMajor>     > >();
   testRandom4x4< UniUpperMatrix< DynamicMatrix<double,rowMajor>  > >();
   testRandom4x4< DiagonalMatrix< DynamicMatrix<double,rowMajor>  > >();

   testRandom4x4< SymmetricMatrix< DynamicMatrix<double,columnMajor> > >();
   testRandom4x4< HermitianMatrix< DynamicMatrix<double,columnMajor> > >();
   testRandom4x4< LowerMatrix< DynamicMatrix<double,columnMajor>     > >();
   testRandom4x4< UniLowerMatrix< DynamicMatrix<double,columnMajor>  > >();
   testRandom4x4< UpperMatrix< DynamicMatrix<double,columnMajor>     > >();
   testRandom4x4< UniUpperMatrix< DynamicMatrix<double,columnMajor>  > >();
   testRandom4x4< DiagonalMatrix< DynamicMatrix<double,columnMajor>  > >();


   //=====================================================================================
   // Random 5x5 matrix tests
   //=====================================================================================

   testRandom5x5< DynamicMatrix<float  ,rowMajor> >();
   testRandom5x5< DynamicMatrix<double ,rowMajor> >();
   testRandom5x5< DynamicMatrix<cfloat ,rowMajor> >();
   testRandom5x5< DynamicMatrix<cdouble,rowMajor> >();

   testRandom5x5< DynamicMatrix<float  ,columnMajor> >();
   testRandom5x5< DynamicMatrix<double ,columnMajor> >();
   testRandom5x5< DynamicMatrix<cfloat ,columnMajor> >();
   testRandom5x5< DynamicMatrix<cdouble,columnMajor> >();

   testRandom5x5< SymmetricMatrix< DynamicMatrix<double,rowMajor> > >();
   testRandom5x5< HermitianMatrix< DynamicMatrix<double,rowMajor> > >();
   testRandom5x5< LowerMatrix< DynamicMatrix<double,rowMajor>     > >();
   testRandom5x5< UniLowerMatrix< DynamicMatrix<double,rowMajor>  > >();
   testRandom5x5< UpperMatrix< DynamicMatrix<double,rowMajor>     > >();
   testRandom5x5< UniUpperMatrix< DynamicMatrix<double,rowMajor>  > >();
   testRandom5x5< DiagonalMatrix< DynamicMatrix<double,rowMajor>  > >();

   testRandom5x5< SymmetricMatrix< DynamicMatrix<double,columnMajor> > >();
   testRandom5x5< HermitianMatrix< DynamicMatrix<double,columnMajor> > >();
   testRandom5x5< LowerMatrix< DynamicMatrix<double,columnMajor>     > >();
   testRandom5x5< UniLowerMatrix< DynamicMatrix<double,columnMajor>  > >();
   testRandom5x5< UpperMatrix< DynamicMatrix<double,columnMajor>     > >();
   testRandom5x5< UniUpperMatrix< DynamicMatrix<double,columnMajor>  > >();
   testRandom5x5< DiagonalMatrix< DynamicMatrix<double,columnMajor>  > >();


   //=====================================================================================
   // Random 6x6 matrix tests
   //=====================================================================================

   testRandom6x6< DynamicMatrix<float  ,rowMajor> >();
   testRandom6x6< DynamicMatrix<double ,rowMajor> >();
   testRandom6x6< DynamicMatrix<cfloat ,rowMajor> >();
   testRandom6x6< DynamicMatrix<cdouble,rowMajor> >();

   testRandom6x6< DynamicMatrix<float  ,columnMajor> >();
   testRandom6x6< DynamicMatrix<double ,columnMajor> >();
   testRandom6x6< DynamicMatrix<cfloat ,columnMajor> >();
   testRandom6x6< DynamicMatrix<cdouble,columnMajor> >();

   testRandom6x6< SymmetricMatrix< DynamicMatrix<double,rowMajor> > >();
   testRandom6x6< HermitianMatrix< DynamicMatrix<double,rowMajor> > >();
   testRandom6x6< LowerMatrix< DynamicMatrix<double,rowMajor>     > >();
   testRandom6x6< UniLowerMatrix< DynamicMatrix<double,rowMajor>  > >();
   testRandom6x6< UpperMatrix< DynamicMatrix<double,rowMajor>     > >();
   testRandom6x6< UniUpperMatrix< DynamicMatrix<double,rowMajor>  > >();
   testRandom6x6< DiagonalMatrix< DynamicMatrix<double,rowMajor>  > >();

   testRandom6x6< SymmetricMatrix< DynamicMatrix<double,columnMajor> > >();
   testRandom6x6< HermitianMatrix< DynamicMatrix<double,columnMajor> > >();
   testRandom6x6< LowerMatrix< DynamicMatrix<double,columnMajor>     > >();
   testRandom6x6< UniLowerMatrix< DynamicMatrix<double,columnMajor>  > >();
   testRandom6x6< UpperMatrix< DynamicMatrix<double,columnMajor>     > >();
   testRandom6x6< UniUpperMatrix< DynamicMatrix<double,columnMajor>  > >();
   testRandom6x6< DiagonalMatrix< DynamicMatrix<double,columnMajor>  > >();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the determinant functionality with specific, predetermined matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function computes determinants for specific, predetermined matrices. In case an error is
// detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testSpecific()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      {
         test_ = "Row-major det() function (0x0)";

         blaze::DynamicMatrix<double,blaze::rowMajor> A;

         const double determinant( det( A ) );

         if( determinant != 0.0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid determinant evaluation\n"
                << " Details:\n"
                << "   Result: " << determinant << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         test_ = "Row-major det() function (1x1)";

         blaze::DynamicMatrix<double,blaze::rowMajor> A( 1UL, 1UL );
         randomize( A );

         const double determinant( det( A ) );

         if( determinant != A(0,0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid determinant evaluation\n"
                << " Details:\n"
                << "   Result: " << determinant << "\n"
                << "   Expected result: " << A(0,0) << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         test_ = "Row-major det() function (unilower)";

         blaze::UniLowerMatrix< blaze::DynamicMatrix<double,blaze::rowMajor> > A( 9UL );
         randomize( A );

         const double determinant( det( A ) );

         if( determinant != 1.0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid determinant evaluation\n"
                << " Details:\n"
                << "   Result: " << determinant << "\n"
                << "   Expected result: 1\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         test_ = "Row-major det() function (uniupper)";

         blaze::UniUpperMatrix< blaze::DynamicMatrix<double,blaze::rowMajor> > A( 9UL );
         randomize( A );

         const double determinant( det( A ) );

         if( determinant != 1.0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid determinant evaluation\n"
                << " Details:\n"
                << "   Result: " << determinant << "\n"
                << "   Expected result: 1\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         test_ = "Row-major det() function (non-square)";

         blaze::DynamicMatrix<double,blaze::rowMajor> A( 2UL, 3UL );

         try {
            const double determinant( det( A ) );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Computing the determinant for a non-square matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << determinant << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      {
         test_ = "Column-major det() function (0x0)";

         blaze::DynamicMatrix<double,blaze::columnMajor> A;

         const double determinant( det( A ) );

         if( determinant != 0.0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid determinant evaluation\n"
                << " Details:\n"
                << "   Result: " << determinant << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         test_ = "Column-major det() function (1x1)";

         blaze::DynamicMatrix<double,blaze::columnMajor> A( 1UL, 1UL );
         randomize( A );

         const double determinant( det( A ) );

         if( determinant != A(0,0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid determinant evaluation\n"
                << " Details:\n"
                << "   Result: " << determinant << "\n"
                << "   Expected result: " << A(0,0) << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         test_ = "Column-major det() function (unilower)";

         blaze::UniLowerMatrix< blaze::DynamicMatrix<double,blaze::columnMajor> > A( 9UL );
         randomize( A );

         const double determinant( det( A ) );

         if( determinant != 1.0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid determinant evaluation\n"
                << " Details:\n"
                << "   Result: " << determinant << "\n"
                << "   Expected result: 1\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         test_ = "Column-major det() function (uniupper)";

         blaze::UniUpperMatrix< blaze::DynamicMatrix<double,blaze::columnMajor> > A( 9UL );
         randomize( A );

         const double determinant( det( A ) );

         if( determinant != 1.0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid determinant evaluation\n"
                << " Details:\n"
                << "   Result: " << determinant << "\n"
                << "   Expected result: 1\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         test_ = "Column-major det() function (non-square)";

         blaze::DynamicMatrix<double,blaze::columnMajor> A( 2UL, 3UL );

         try {
            const double determinant( det( A ) );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Computing the determinant for a non-square matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << determinant << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

#endif
}
//*************************************************************************************************

} // namespace determinant

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
   std::cout << "   Running dense matrix determinant test..." << std::endl;

   try
   {
      RUN_DETERMINANT_DENSE_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during dense matrix determinant test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
