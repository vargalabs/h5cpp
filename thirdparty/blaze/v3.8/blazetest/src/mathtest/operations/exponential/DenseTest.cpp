//=================================================================================================
/*!
//  \file src/mathtest/operations/exponential/DenseTest.cpp
//  \brief Source file for the dense matrix exponential test
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
#include <blaze/math/IdentityMatrix.h>
#include <blaze/math/LowerMatrix.h>
#include <blaze/math/StrictlyLowerMatrix.h>
#include <blaze/math/StrictlyUpperMatrix.h>
#include <blaze/math/SymmetricMatrix.h>
#include <blaze/math/UniLowerMatrix.h>
#include <blaze/math/UniUpperMatrix.h>
#include <blaze/math/UpperMatrix.h>
#include <blaze/math/ZeroMatrix.h>
#include <blaze/util/Complex.h>
#include <blazetest/mathtest/operations/exponential/DenseTest.h>

#ifdef BLAZE_USE_HPX_THREADS
#  include <hpx/hpx_main.hpp>
#endif


namespace blazetest {

namespace mathtest {

namespace operations {

namespace exponential {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the DenseTest exponential test.
//
// \exception std::runtime_error Matrix exponential error detected.
*/
DenseTest::DenseTest()
{
   using blaze::DynamicMatrix;
   using blaze::SymmetricMatrix;
   using blaze::HermitianMatrix;
   using blaze::LowerMatrix;
   using blaze::UniLowerMatrix;
   using blaze::StrictlyLowerMatrix;
   using blaze::UpperMatrix;
   using blaze::UniUpperMatrix;
   using blaze::StrictlyUpperMatrix;
   using blaze::DiagonalMatrix;

   using cplx = blaze::complex<double>;


   //=====================================================================================
   // Specific matrix tests
   //=====================================================================================

   testSpecific();


   //=====================================================================================
   // Random matrix tests
   //=====================================================================================

   for( size_t i=0UL; i<8UL; ++i )
   {
      testRandom< DynamicMatrix<double> >( i );
      testRandom< DynamicMatrix<cplx  > >( i );

      testRandom< SymmetricMatrix    < DynamicMatrix<double> > >( i );
      testRandom< SymmetricMatrix    < DynamicMatrix<cplx  > > >( i );
      testRandom< HermitianMatrix    < DynamicMatrix<double> > >( i );
      testRandom< HermitianMatrix    < DynamicMatrix<cplx  > > >( i );
      testRandom< LowerMatrix        < DynamicMatrix<double> > >( i );
      testRandom< LowerMatrix        < DynamicMatrix<cplx  > > >( i );
      testRandom< UniLowerMatrix     < DynamicMatrix<double> > >( i );
      testRandom< UniLowerMatrix     < DynamicMatrix<cplx  > > >( i );
      testRandom< StrictlyLowerMatrix< DynamicMatrix<double> > >( i );
      testRandom< StrictlyLowerMatrix< DynamicMatrix<cplx  > > >( i );
      testRandom< UpperMatrix        < DynamicMatrix<double> > >( i );
      testRandom< UpperMatrix        < DynamicMatrix<cplx  > > >( i );
      testRandom< UniUpperMatrix     < DynamicMatrix<double> > >( i );
      testRandom< UniUpperMatrix     < DynamicMatrix<cplx  > > >( i );
      testRandom< StrictlyUpperMatrix< DynamicMatrix<double> > >( i );
      testRandom< StrictlyUpperMatrix< DynamicMatrix<cplx  > > >( i );
      testRandom< DiagonalMatrix     < DynamicMatrix<double> > >( i );
      testRandom< DiagonalMatrix     < DynamicMatrix<cplx  > > >( i );
   }
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
         test_ = "Row-major dense matrix exponential (0x0)";

         blaze::DynamicMatrix<double,blaze::rowMajor> A;
         blaze::DynamicMatrix<double,blaze::rowMajor> B( matexp( A ) );

         if( B.rows() != 0UL || B.columns() != 0UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix exponential failed\n"
                << " Details:\n"
                << "   matexp(A):\n" << B << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         test_ = "Row-major dense matrix exponential ( matexp(Z) )";

         blaze::DynamicMatrix<double,blaze::rowMajor> A( blaze::ZeroMatrix<double>( 9UL, 9UL ) );
         blaze::DynamicMatrix<double,blaze::rowMajor> B( matexp( A ) );

         if( B.rows() != 9UL || B.columns() != 9UL || !isIdentity( B ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix exponential failed\n"
                << " Details:\n"
                << "   matexp(A):\n" << B << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         test_ = "Row-major dense matrix exponential ( matexp(I) )";

         blaze::DynamicMatrix<double,blaze::rowMajor> A( blaze::IdentityMatrix<double>( 9UL ) );
         blaze::DynamicMatrix<double,blaze::rowMajor> B( matexp( A ) );

         if( B.rows() != 9UL || B.columns() != 9UL ||
             !isDiagonal( B ) ||
             diagonal( B ) != blaze::uniform( 9UL, std::exp(1.0) ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix exponential failed\n"
                << " Details:\n"
                << "   matexp(A):\n" << B << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         test_ = "Row-major dense matrix exponential ( matexp(trans(A)) == trans(matexp(A)) )";

         blaze::DynamicMatrix<double,blaze::rowMajor> A( 9UL, 9UL );
         randomize( A, -1.0, 1.0 );

         blaze::DynamicMatrix<double,blaze::rowMajor> B( matexp( trans( A ) ) );
         blaze::DynamicMatrix<double,blaze::rowMajor> C( trans( matexp( A ) ) );

         if( B != C ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix exponential failed\n"
                << " Details:\n"
                << "   matexp( trans(A) ):\n" << B << "\n"
                << "   trans( matexp(A) ):\n" << C << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         test_ = "Row-major dense matrix exponential ( matexp(ctrans(A)) == ctrans(matexp(A)) )";

         blaze::DynamicMatrix<complex<double>,blaze::rowMajor> A( 9UL, 9UL );
         randomize( A, -1.0, 1.0 );

         blaze::DynamicMatrix<complex<double>,blaze::rowMajor> B( matexp( ctrans( A ) ) );
         blaze::DynamicMatrix<complex<double>,blaze::rowMajor> C( ctrans( matexp( A ) ) );

         if( B != C ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix exponential failed\n"
                << " Details:\n"
                << "   matexp( ctrans(A) ):\n" << B << "\n"
                << "   ctrans( matexp(A) ):\n" << C << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         test_ = "Row-major dense matrix exponential ( det(matexp(A)) == exp(trace(A)) )";

         blaze::DynamicMatrix<double,blaze::rowMajor> A( 9UL, 9UL );
         randomize( A, -1.0, 1.0 );

         blaze::DynamicMatrix<double,blaze::rowMajor> B( matexp( A ) );

         if( !isEqual( det( B ), exp( trace( A ) ) ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix exponential failed\n"
                << " Details:\n"
                << "   det( matexp(A) ):\n" << det( B ) << "\n"
                << "   exp( trace(A) ):\n" << exp( trace(A) ) << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         test_ = "Row-major dense matrix exponential ( matexp(-A) * matexp(A) == I )";

         blaze::DynamicMatrix<double,blaze::rowMajor> A( 9UL, 9UL );
         randomize( A, -1.0, 1.0 );

         blaze::DynamicMatrix<double,blaze::rowMajor> B( matexp( -A ) * matexp( A ) );
         blaze::DynamicMatrix<double,blaze::rowMajor> C( matexp( A ) * matexp( -A ) );

         if( B.rows() != 9UL || B.columns() != 9UL || !isIdentity( B ) ||
             C.rows() != 9UL || C.columns() != 9UL || !isIdentity( C ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix exponential failed\n"
                << " Details:\n"
                << "   matexp(-A)*matexp(A):\n" << B << "\n"
                << "   matexp(A)*matexp(-A):\n" << C << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         test_ = "Row-major dense matrix exponential ( matexp(aA)*matexp(bA) == matexp((a+b)A) )";

         blaze::DynamicMatrix<double,blaze::rowMajor> A( 9UL, 9UL );
         randomize( A, -1.0, 1.0 );

         const double a( blaze::rand<double>( -1.0, 1.0 ) );
         const double b( blaze::rand<double>( -1.0, 1.0 ) );

         blaze::DynamicMatrix<double,blaze::rowMajor> B( matexp( a*A ) * matexp( b*A ) );
         blaze::DynamicMatrix<double,blaze::rowMajor> C( matexp( (a+b)*A ) );

         if( B.rows() != 9UL || B.columns() != 9UL || C.rows() != 9UL || C.columns() != 9UL ||
             B != C ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix exponential failed\n"
                << " Details:\n"
                << "   matexp(a*A) * matexp(b*A):\n" << B << "\n"
                << "   matexp((a+b)*A):\n" << C << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         test_ = "Row-major dense matrix exponential ( matexp(B*A*inv(B)) == B*matexp(A)*inv(B) )";

         blaze::DynamicMatrix<double,blaze::rowMajor> A( 9UL, 9UL );
         randomize( A, -1.0, 1.0 );

         blaze::DynamicMatrix<double,blaze::rowMajor> B( 9UL, 9UL );
         do {
            randomize( B, -1.0, 1.0 );
         }
         while( blaze::isZero( det( B ) ) );

         blaze::DynamicMatrix<double,blaze::rowMajor> C( matexp( B*A*inv(B) ) );
         blaze::DynamicMatrix<double,blaze::rowMajor> D( B*matexp(A)*inv(B) );

         if( C != D ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix exponential failed\n"
                << " Details:\n"
                << "   matexp( B*A*inv(B) ):\n" << C << "\n"
                << "   B*matexp(A)*inv(B):\n" << D << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         test_ = "Row-major dense matrix exponential (non-square)";

         blaze::DynamicMatrix<double,blaze::rowMajor> A( 2UL, 3UL );

         try {
            matexp( A );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Exponential of a non-square matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << A << "\n";
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
         test_ = "Column-major dense matrix exponential (0x0)";

         blaze::DynamicMatrix<double,blaze::columnMajor> A;
         blaze::DynamicMatrix<double,blaze::columnMajor> B( matexp( A ) );

         if( B.rows() != 0UL || B.columns() != 0UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix exponential failed\n"
                << " Details:\n"
                << "   matexp(A):\n" << B << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         test_ = "Column-major dense matrix exponential ( matexp(Z) )";

         blaze::DynamicMatrix<double,blaze::columnMajor> A( blaze::ZeroMatrix<double>( 9UL, 9UL ) );
         blaze::DynamicMatrix<double,blaze::columnMajor> B( matexp( A ) );

         if( B.rows() != 9UL || B.columns() != 9UL || !isIdentity( B ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix exponential failed\n"
                << " Details:\n"
                << "   matexp(A):\n" << B << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         test_ = "Column-major dense matrix exponential ( matexp(I) )";

         blaze::DynamicMatrix<double,blaze::columnMajor> A( blaze::IdentityMatrix<double>( 9UL ) );
         blaze::DynamicMatrix<double,blaze::columnMajor> B( matexp( A ) );

         if( B.rows() != 9UL || B.columns() != 9UL ||
             !isDiagonal( B ) ||
             diagonal( B ) != blaze::uniform( 9UL, std::exp(1.0) ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix exponential failed\n"
                << " Details:\n"
                << "   matexp(A):\n" << B << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         test_ = "Column-major dense matrix exponential ( matexp(trans(A)) == trans(matexp(A)) )";

         blaze::DynamicMatrix<double,blaze::columnMajor> A( 9UL, 9UL );
         randomize( A, -1.0, 1.0 );

         blaze::DynamicMatrix<double,blaze::columnMajor> B( matexp( trans( A ) ) );
         blaze::DynamicMatrix<double,blaze::columnMajor> C( trans( matexp( A ) ) );

         if( B != C ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix exponential failed\n"
                << " Details:\n"
                << "   matexp( trans(A) ):\n" << B << "\n"
                << "   trans( matexp(A) ):\n" << C << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         test_ = "Column-major dense matrix exponential ( matexp(ctrans(A)) == ctrans(matexp(A)) )";

         blaze::DynamicMatrix<complex<double>,blaze::columnMajor> A( 9UL, 9UL );
         randomize( A, -1.0, 1.0 );

         blaze::DynamicMatrix<complex<double>,blaze::columnMajor> B( matexp( ctrans( A ) ) );
         blaze::DynamicMatrix<complex<double>,blaze::columnMajor> C( ctrans( matexp( A ) ) );

         if( B != C ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix exponential failed\n"
                << " Details:\n"
                << "   matexp( ctrans(A) ):\n" << B << "\n"
                << "   ctrans( matexp(A) ):\n" << C << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         test_ = "Column-major dense matrix exponential ( det(matexp(A)) == exp(trace(A)) )";

         blaze::DynamicMatrix<double,blaze::columnMajor> A( 9UL, 9UL );
         randomize( A, -1.0, 1.0 );

         blaze::DynamicMatrix<double,blaze::columnMajor> B( matexp( A ) );

         if( !isEqual( det( B ), exp( trace( A ) ) ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix exponential failed\n"
                << " Details:\n"
                << "   det( matexp(A) ):\n" << det( B ) << "\n"
                << "   exp( trace(A) ):\n" << exp( trace(A) ) << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         test_ = "Column-major dense matrix exponential ( matexp(-A) * matexp(A) == I )";

         blaze::DynamicMatrix<double,blaze::columnMajor> A( 9UL, 9UL );
         randomize( A, -1.0, 1.0 );

         blaze::DynamicMatrix<double,blaze::columnMajor> B( matexp( -A ) * matexp( A ) );
         blaze::DynamicMatrix<double,blaze::columnMajor> C( matexp( A ) * matexp( -A ) );

         if( B.rows() != 9UL || B.columns() != 9UL || !isIdentity( B ) ||
             C.rows() != 9UL || C.columns() != 9UL || !isIdentity( C ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix exponential failed\n"
                << " Details:\n"
                << "   matexp(-A)*matexp(A):\n" << B << "\n"
                << "   matexp(A)*matexp(-A):\n" << C << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         test_ = "Column-major dense matrix exponential ( matexp(aA)*matexp(bA) == matexp((a+b)A) )";

         blaze::DynamicMatrix<double,blaze::columnMajor> A( 9UL, 9UL );
         randomize( A, -1.0, 1.0 );

         const double a( blaze::rand<double>( -1.0, 1.0 ) );
         const double b( blaze::rand<double>( -1.0, 1.0 ) );

         blaze::DynamicMatrix<double,blaze::columnMajor> B( matexp( a*A ) * matexp( b*A ) );
         blaze::DynamicMatrix<double,blaze::columnMajor> C( matexp( (a+b)*A ) );

         if( B.rows() != 9UL || B.columns() != 9UL || C.rows() != 9UL || C.columns() != 9UL ||
             B != C ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix exponential failed\n"
                << " Details:\n"
                << "   matexp(a*A) * matexp(b*A):\n" << B << "\n"
                << "   matexp((a+b)*A):\n" << C << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         test_ = "Column-major dense matrix exponential ( matexp(B*A*inv(B)) == B*matexp(A)*inv(B) )";

         blaze::DynamicMatrix<double,blaze::columnMajor> A( 9UL, 9UL );
         randomize( A, -1.0, 1.0 );

         blaze::DynamicMatrix<double,blaze::columnMajor> B( 9UL, 9UL );
         do {
            randomize( B, -1.0, 1.0 );
         }
         while( blaze::isZero( det( B ) ) );

         blaze::DynamicMatrix<double,blaze::columnMajor> C( matexp( B*A*inv(B) ) );
         blaze::DynamicMatrix<double,blaze::columnMajor> D( B*matexp(A)*inv(B) );

         if( C != D ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix exponential failed\n"
                << " Details:\n"
                << "   matexp( B*A*inv(B) ):\n" << C << "\n"
                << "   B*matexp(A)*inv(B):\n" << D << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         test_ = "Column-major dense matrix exponential (non-square)";

         blaze::DynamicMatrix<double,blaze::columnMajor> A( 2UL, 3UL );

         try {
            matexp( A );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Exponential of a non-square matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << A << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

#endif
}
//*************************************************************************************************

} // namespace exponential

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
   std::cout << "   Running dense matrix exponential test..." << std::endl;

   try
   {
      RUN_EXPONENTIAL_DENSE_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during dense matrix exponential test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
