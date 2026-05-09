//=================================================================================================
/*!
//  \file blazetest/mathtest/operations/lse/DenseTest.h
//  \brief Header file for the dense matrix LSE test
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

#ifndef _BLAZETEST_MATHTEST_OPERATIONS_LSE_DENSETEST_H_
#define _BLAZETEST_MATHTEST_OPERATIONS_LSE_DENSETEST_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <sstream>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <blaze/math/DiagonalMatrix.h>
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/HermitianMatrix.h>
#include <blaze/math/LowerMatrix.h>
#include <blaze/math/Submatrix.h>
#include <blaze/math/SymmetricMatrix.h>
#include <blaze/math/UniLowerMatrix.h>
#include <blaze/math/UniUpperMatrix.h>
#include <blaze/math/UpperMatrix.h>
#include <blaze/util/Random.h>
#include <blazetest/system/LAPACK.h>
#include <blazetest/system/Types.h>


namespace blazetest {

namespace mathtest {

namespace operations {

namespace lse {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class for all dense matrix LSE tests.
//
// This class represents a test suite for the dense matrix LSE kernels. It solves a series of
// LSEs with various sizes and on all dense matrix types of the Blaze library.
*/
class DenseTest
{
 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit DenseTest();
   // No explicitly declared copy constructor.
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   // No explicitly declared destructor.
   //**********************************************************************************************

 private:
   //**Test functions******************************************************************************
   /*!\name Test functions */
   //@{
   template< typename Type > void testGeneral  ( size_t N );
   template< typename Type > void testSymmetric( size_t N );
   template< typename Type > void testHermitian( size_t N );
   template< typename Type > void testLower    ( size_t N );
   template< typename Type > void testUniLower ( size_t N );
   template< typename Type > void testUpper    ( size_t N );
   template< typename Type > void testUniUpper ( size_t N );
   template< typename Type > void testDiagonal ( size_t N );
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   std::string test_;  //!< Label of the currently performed test.
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the LSE kernels with random \f$ N \times N \f$ general matrices.
//
// \param N The number of rows and columns of the matrix.
// \return void
// \exception std::runtime_error Error detected.
//
// This function tests the dense matrix LSE kernels for random \f$ N \times N \f$ general
// matrices. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void DenseTest::testGeneral( size_t N )
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   using blaze::DynamicMatrix;
   using blaze::DynamicVector;
   using blaze::rowMajor;
   using blaze::columnMajor;
   using blaze::rowVector;
   using blaze::solve;
   using blaze::isDefault;


   //=====================================================================================
   // Single right-hand side
   //=====================================================================================

   {
      test_ = "General LSE (single rhs, automatic)";

      DynamicMatrix<Type> A( N, N );
      DynamicVector<Type> b( N );

      if( N != 0UL ) {
         do {
            randomize( A );
         }
         while( isDefault( det( A ) ) );

         randomize( b );
      }

      const DynamicMatrix<Type,rowMajor>    A1( A );
      const DynamicMatrix<Type,columnMajor> A2( A );

      DynamicVector<Type> x1( N );
      DynamicVector<Type> x2( N );

      x1 = solve( A1, b );
      x2 = solve( A2, b );

      if( A*x1 != b || A*x2 != b || x1 != x2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   Row-major solution (x1):\n" << x1 << "\n"
             << "   Column-major solution (x2):\n" << x2 << "\n"
             << "   A * x1 =\n" << ( A * x1 ) << "\n"
             << "   A * x2 =\n" << ( A * x2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "General LSE (single rhs, automatic, transpose)";

      DynamicMatrix<Type> A( N, N );
      DynamicVector<Type> b( N );

      if( N != 0UL ) {
         do {
            randomize( A );
         }
         while( isDefault( det( A ) ) );

         randomize( b );
      }

      const DynamicMatrix<Type,rowMajor>    A1( A );
      const DynamicMatrix<Type,columnMajor> A2( A );

      DynamicVector<Type> x1( N );
      DynamicVector<Type> x2( N );

      x1 = solve( trans(A1), b );
      x2 = solve( trans(A2), b );

      if( trans(A)*x1 != b || trans(A)*x2 != b || x1 != x2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << trans(A) << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   Row-major solution (x1):\n" << x1 << "\n"
             << "   Column-major solution (x2):\n" << x2 << "\n"
             << "   A * x1 =\n" << ( trans(A) * x1 ) << "\n"
             << "   A * x2 =\n" << ( trans(A) * x2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "General LSE (single rhs, inv(A)*b)";

      DynamicMatrix<Type> A( N, N );
      DynamicVector<Type> b( N );

      if( N != 0UL ) {
         do {
            randomize( A );
         }
         while( isDefault( det( A ) ) );

         randomize( b );
      }

      const DynamicMatrix<Type,rowMajor>    A1( A );
      const DynamicMatrix<Type,columnMajor> A2( A );

      DynamicVector<Type> x1( N );
      DynamicVector<Type> x2( N );

      x1 = inv(A1) * b;
      x2 = inv(A2) * b;

      if( A*x1 != b || A*x2 != b || x1 != x2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   Row-major solution (x1):\n" << x1 << "\n"
             << "   Column-major solution (x2):\n" << x2 << "\n"
             << "   A * x1 =\n" << ( A * x1 ) << "\n"
             << "   A * x2 =\n" << ( A * x2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "General LSE (single rhs, b*inv(A))";

      DynamicMatrix<Type> A( N, N );
      DynamicVector<Type,rowVector> b( N );

      if( N != 0UL ) {
         do {
            randomize( A );
         }
         while( isDefault( det( A ) ) );

         randomize( b );
      }

      const DynamicMatrix<Type,rowMajor>    A1( A );
      const DynamicMatrix<Type,columnMajor> A2( A );

      DynamicVector<Type,rowVector> x1( N );
      DynamicVector<Type,rowVector> x2( N );

      x1 = b * inv(A1);
      x2 = b * inv(A2);

      if( x1*A != b || x2*A != b || x1 != x2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   Row-major solution (x1):\n" << x1 << "\n"
             << "   Column-major solution (x2):\n" << x2 << "\n"
             << "   x1 * A =\n" << ( x1 * A ) << "\n"
             << "   x2 * A =\n" << ( x2 * A ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "General LSE (single rhs, non-square)";

      const DynamicMatrix<Type,rowMajor> A( 2UL, 3UL );
      const DynamicVector<Type> b( 2UL );
      DynamicVector<Type> x;

      try {
         x = solve( A, b );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE with non-square system matrix succeeded\n"
             << " Details:\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   Solution (x):\n" << x << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   {
      test_ = "General LSE (single rhs, non-square)";

      const DynamicMatrix<Type,rowMajor> A( 3UL, 2UL );
      const DynamicVector<Type> b( 2UL );
      DynamicVector<Type> x;

      try {
         x = solve( A, b );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE with non-square system matrix succeeded\n"
             << " Details:\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   Solution (x):\n" << x << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   {
      test_ = "General LSE (single rhs, non-matching right-hand side)";

      const DynamicMatrix<Type,rowMajor> A( 2UL, 2UL );
      const DynamicVector<Type> b( 3UL );
      DynamicVector<Type> x;

      try {
         x = solve( A, b );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE with invalid right-hand side succeeded\n"
             << " Details:\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   Solution (x):\n" << x << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Multiple right-hand sides
   //=====================================================================================

   {
      test_ = "General LSE (multiple rhs, automatic)";

      DynamicMatrix<Type> A( N, N );
      DynamicMatrix<Type> B( N, 3UL );

      if( N != 0UL ) {
         do {
            randomize( A );
         }
         while( isDefault( det( A ) ) );

         randomize( B );
      }

      const DynamicMatrix<Type,rowMajor>    A1( A );
      const DynamicMatrix<Type,columnMajor> A2( A );

      const DynamicMatrix<Type,rowMajor>    B1( B );
      const DynamicMatrix<Type,columnMajor> B2( B );

      DynamicMatrix<Type,rowMajor>    X1;
      DynamicMatrix<Type,columnMajor> X2;

      X1 = solve( A1, B1 );
      X2 = solve( A2, B2 );

      if( A*X1 != B || A*X2 != B || X1 != X2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Right-hand sides (B):\n" << B << "\n"
             << "   Row-major solutions (X1):\n" << X1 << "\n"
             << "   Column-major solutions (X2):\n" << X2 << "\n"
             << "   A * X1 =\n" << ( A * X1 ) << "\n"
             << "   A * X2 =\n" << ( A * X2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "General LSE (multiple rhs, automatic, transpose)";

      DynamicMatrix<Type> A( N, N );
      DynamicMatrix<Type> B( 3UL, N );

      if( N != 0UL ) {
         do {
            randomize( A );
         }
         while( isDefault( det( A ) ) );

         randomize( B );
      }

      const DynamicMatrix<Type,rowMajor>    A1( A );
      const DynamicMatrix<Type,columnMajor> A2( A );

      const DynamicMatrix<Type,rowMajor>    B1( B );
      const DynamicMatrix<Type,columnMajor> B2( B );

      DynamicMatrix<Type,rowMajor>    X1;
      DynamicMatrix<Type,columnMajor> X2;

      X1 = solve( trans(A1), trans(B1) );
      X2 = solve( trans(A2), trans(B2) );

      if( trans(A)*X1 != trans(B) || trans(A)*X2 != trans(B) || X1 != X2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << trans(A) << "\n"
             << "   Right-hand sides (B):\n" << trans(B) << "\n"
             << "   Row-major solutions (X1):\n" << X1 << "\n"
             << "   Column-major solutions (X2):\n" << X2 << "\n"
             << "   A * X1 =\n" << ( trans(A) * X1 ) << "\n"
             << "   A * X2 =\n" << ( trans(A) * X2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "General LSE (multiple rhs, inv(A)*B)";

      DynamicMatrix<Type> A( N, N );
      DynamicMatrix<Type> B( N, 3UL );

      if( N != 0UL ) {
         do {
            randomize( A );
         }
         while( isDefault( det( A ) ) );

         randomize( B );
      }

      const DynamicMatrix<Type,rowMajor>    A1( A );
      const DynamicMatrix<Type,columnMajor> A2( A );

      const DynamicMatrix<Type,rowMajor>    B1( B );
      const DynamicMatrix<Type,columnMajor> B2( B );

      DynamicMatrix<Type,rowMajor>    X1;
      DynamicMatrix<Type,columnMajor> X2;

      X1 = inv(A1) * B1;
      X2 = inv(A2) * B2;

      if( A*X1 != B || A*X2 != B || X1 != X2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Right-hand sides (B):\n" << B << "\n"
             << "   Row-major solutions (X1):\n" << X1 << "\n"
             << "   Column-major solutions (X2):\n" << X2 << "\n"
             << "   A * X1 =\n" << ( A * X1 ) << "\n"
             << "   A * X2 =\n" << ( A * X2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "General LSE (multiple rhs, B*inv(A))";

      DynamicMatrix<Type> A( N, N );
      DynamicMatrix<Type> B( 3UL, N );

      if( N != 0UL ) {
         do {
            randomize( A );
         }
         while( isDefault( det( A ) ) );

         randomize( B );
      }

      const DynamicMatrix<Type,rowMajor>    A1( A );
      const DynamicMatrix<Type,columnMajor> A2( A );

      const DynamicMatrix<Type,rowMajor>    B1( B );
      const DynamicMatrix<Type,columnMajor> B2( B );

      DynamicMatrix<Type,rowMajor>    X1;
      DynamicMatrix<Type,columnMajor> X2;

      X1 = B1 * inv(A1);
      X2 = B2 * inv(A2);

      if( X1*A != B || X2*A != B || X1 != X2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Right-hand sides (B):\n" << B << "\n"
             << "   Row-major solutions (X1):\n" << X1 << "\n"
             << "   Column-major solutions (X2):\n" << X2 << "\n"
             << "   X1 * A =\n" << ( X1 * A ) << "\n"
             << "   X2 * A =\n" << ( X2 * A ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "General LSE (multiple rhs, non-square)";

      DynamicMatrix<Type,rowMajor> A( 2UL, 3UL );
      DynamicMatrix<Type> B( 2UL, 5UL );
      DynamicMatrix<Type> X;

      try {
         X = solve( A, B );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE with non-square system matrix succeeded\n"
             << " Details:\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Right-hand side (B):\n" << B << "\n"
             << "   Solution (X):\n" << X << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   {
      test_ = "General LSE (multiple rhs, non-square)";

      DynamicMatrix<Type,rowMajor> A( 3UL, 2UL );
      DynamicMatrix<Type> B( 2UL, 5UL );
      DynamicMatrix<Type> X;

      try {
         X = solve( A, B );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE with non-square system matrix succeeded\n"
             << " Details:\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Right-hand side (B):\n" << B << "\n"
             << "   Solution (X):\n" << X << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   {
      test_ = "General LSE (multiple rhs, non-matching right-hand side)";

      DynamicMatrix<Type,rowMajor> A( 2UL, 2UL );
      DynamicMatrix<Type> B( 3UL, 5UL );
      DynamicMatrix<Type> X;

      try {
         X = solve( A, B );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE with invalid right-hand side succeeded\n"
             << " Details:\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Right-hand side (B):\n" << B << "\n"
             << "   Solution (X):\n" << X << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the LSE kernels with random \f$ N \times N \f$ symmetric matrices.
//
// \param N The number of rows and columns of the matrix.
// \return void
// \exception std::runtime_error Error detected.
//
// This function tests the dense matrix LSE kernels for random \f$ N \times N \f$ symmetric
// matrices. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void DenseTest::testSymmetric( size_t N )
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   using blaze::SymmetricMatrix;
   using blaze::DynamicMatrix;
   using blaze::DynamicVector;
   using blaze::rowMajor;
   using blaze::columnMajor;
   using blaze::solve;
   using blaze::isDefault;


   //=====================================================================================
   // Single right-hand side
   //=====================================================================================

   {
      test_ = "Symmetric LSE (single rhs, automatic)";

      SymmetricMatrix< DynamicMatrix<Type> > A( N );
      DynamicVector<Type> b( N );

      if( N != 0UL ) {
         do {
            randomize( A );
         }
         while( isDefault( det( A ) ) );

         randomize( b );
      }

      const SymmetricMatrix< DynamicMatrix<Type,rowMajor>    > A1( A );
      const SymmetricMatrix< DynamicMatrix<Type,columnMajor> > A2( A );

      DynamicVector<Type> x1( N );
      DynamicVector<Type> x2( N );

      x1 = solve( A1, b );
      x2 = solve( A2, b );

      if( A*x1 != b || A*x2 != b || x1 != x2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   Row-major solution (x1):\n" << x1 << "\n"
             << "   Column-major solution (x2):\n" << x2 << "\n"
             << "   A * x1 =\n" << ( A * x1 ) << "\n"
             << "   A * x2 =\n" << ( A * x2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Symmetric LSE (single rhs, automatic, transpose)";

      SymmetricMatrix< DynamicMatrix<Type> > A( N );
      DynamicVector<Type> b( N );

      if( N != 0UL ) {
         do {
            randomize( A );
         }
         while( isDefault( det( A ) ) );

         randomize( b );
      }

      const SymmetricMatrix< DynamicMatrix<Type,rowMajor>    > A1( A );
      const SymmetricMatrix< DynamicMatrix<Type,columnMajor> > A2( A );

      DynamicVector<Type> x1( N );
      DynamicVector<Type> x2( N );

      x1 = solve( trans(A1), b );
      x2 = solve( trans(A2), b );

      if( trans(A)*x1 != b || trans(A)*x2 != b || x1 != x2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << trans(A) << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   Row-major solution (x1):\n" << x1 << "\n"
             << "   Column-major solution (x2):\n" << x2 << "\n"
             << "   A * x1 =\n" << ( trans(A) * x1 ) << "\n"
             << "   A * x2 =\n" << ( trans(A) * x2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Symmetric LSE (single rhs, declsym)";

      SymmetricMatrix< DynamicMatrix<Type> > A( N );
      DynamicVector<Type> b( N );

      if( N != 0UL ) {
         do {
            randomize( A );
         }
         while( isDefault( det( A ) ) );

         randomize( b );
      }

      const DynamicMatrix<Type,rowMajor>    A1( A );
      const DynamicMatrix<Type,columnMajor> A2( A );

      DynamicVector<Type> x1( N );
      DynamicVector<Type> x2( N );

      solve( declsym( A1 ), x1, b );
      solve( declsym( A2 ), x2, b );

      if( A*x1 != b || A*x2 != b || x1 != x2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   Row-major solution (x1):\n" << x1 << "\n"
             << "   Column-major solution (x2):\n" << x2 << "\n"
             << "   A * x1 =\n" << ( A * x1 ) << "\n"
             << "   A * x2 =\n" << ( A * x2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Symmetric LSE (single rhs, declsym, transpose)";

      SymmetricMatrix< DynamicMatrix<Type> > A( N );
      DynamicVector<Type> b( N );

      if( N != 0UL ) {
         do {
            randomize( A );
         }
         while( isDefault( det( A ) ) );

         randomize( b );
      }

      const DynamicMatrix<Type,rowMajor>    A1( A );
      const DynamicMatrix<Type,columnMajor> A2( A );

      DynamicVector<Type> x1( N );
      DynamicVector<Type> x2( N );

      solve( declsym( trans(A1) ), x1, b );
      solve( declsym( trans(A2) ), x2, b );

      if( A*x1 != b || A*x2 != b || x1 != x2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << trans(A) << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   Row-major solution (x1):\n" << x1 << "\n"
             << "   Column-major solution (x2):\n" << x2 << "\n"
             << "   A * x1 =\n" << ( trans(A) * x1 ) << "\n"
             << "   A * x2 =\n" << ( trans(A) * x2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "General LSE (single rhs, non-matching right-hand side)";

      const SymmetricMatrix< DynamicMatrix<Type> > A( 2UL );
      const DynamicVector<Type> b( 3UL );
      DynamicVector<Type> x;

      try {
         x = solve( A, b );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE with invalid right-hand side succeeded\n"
             << " Details:\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   Solution (x):\n" << x << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Multiple right-hand sides
   //=====================================================================================

   {
      test_ = "Symmetric LSE (multiple rhs, automatic)";

      SymmetricMatrix< DynamicMatrix<Type> > A( N );
      DynamicMatrix<Type> B( N, 3UL );

      if( N != 0UL ) {
         do {
            randomize( A );
         }
         while( isDefault( det( A ) ) );

         randomize( B );
      }

      const SymmetricMatrix< DynamicMatrix<Type,rowMajor> >    A1( A );
      const SymmetricMatrix< DynamicMatrix<Type,columnMajor> > A2( A );

      const DynamicMatrix<Type,rowMajor>    B1( B );
      const DynamicMatrix<Type,columnMajor> B2( B );

      DynamicMatrix<Type,rowMajor>    X1;
      DynamicMatrix<Type,columnMajor> X2;

      X1 = solve( A1, B1 );
      X2 = solve( A2, B2 );

      if( A*X1 != B || A*X2 != B || X1 != X2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Right-hand sides (B):\n" << B << "\n"
             << "   Row-major solutions (X1):\n" << X1 << "\n"
             << "   Column-major solutions (X2):\n" << X2 << "\n"
             << "   A * X1 =\n" << ( A * X1 ) << "\n"
             << "   A * X2 =\n" << ( A * X2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Symmetric LSE (multiple rhs, automatic, transpose)";

      SymmetricMatrix< DynamicMatrix<Type> > A( N );
      DynamicMatrix<Type> B( 3UL, N );

      if( N != 0UL ) {
         do {
            randomize( A );
         }
         while( isDefault( det( A ) ) );

         randomize( B );
      }

      const SymmetricMatrix< DynamicMatrix<Type,rowMajor> >    A1( A );
      const SymmetricMatrix< DynamicMatrix<Type,columnMajor> > A2( A );

      const DynamicMatrix<Type,rowMajor>    B1( B );
      const DynamicMatrix<Type,columnMajor> B2( B );

      DynamicMatrix<Type,rowMajor>    X1;
      DynamicMatrix<Type,columnMajor> X2;

      X1 = solve( trans(A1), trans(B1) );
      X2 = solve( trans(A2), trans(B2) );

      if( trans(A)*X1 != trans(B) || trans(A)*X2 != trans(B) || X1 != X2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << trans(A) << "\n"
             << "   Right-hand sides (B):\n" << trans(B) << "\n"
             << "   Row-major solutions (X1):\n" << X1 << "\n"
             << "   Column-major solutions (X2):\n" << X2 << "\n"
             << "   A * X1 =\n" << ( trans(A) * X1 ) << "\n"
             << "   A * X2 =\n" << ( trans(A) * X2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Symmetric LSE (multiple rhs, declsym)";

      SymmetricMatrix< DynamicMatrix<Type> > A( N );
      DynamicMatrix<Type> B( N, 3UL );

      if( N != 0UL ) {
         do {
            randomize( A );
         }
         while( isDefault( det( A ) ) );

         randomize( B );
      }

      const DynamicMatrix<Type,rowMajor>    A1( A );
      const DynamicMatrix<Type,columnMajor> A2( A );

      const DynamicMatrix<Type,rowMajor>    B1( B );
      const DynamicMatrix<Type,columnMajor> B2( B );

      DynamicMatrix<Type,rowMajor>    X1;
      DynamicMatrix<Type,columnMajor> X2;

      solve( declsym( A1 ), X1, B1 );
      solve( declsym( A2 ), X2, B2 );

      if( A*X1 != B || A*X2 != B || X1 != X2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Right-hand sides (B):\n" << B << "\n"
             << "   Row-major solutions (X1):\n" << X1 << "\n"
             << "   Column-major solutions (X2):\n" << X2 << "\n"
             << "   A * X1 =\n" << ( A * X1 ) << "\n"
             << "   A * X2 =\n" << ( A * X2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Symmetric LSE (multiple rhs, declsym, transpose)";

      SymmetricMatrix< DynamicMatrix<Type> > A( N );
      DynamicMatrix<Type> B( 3UL, N );

      if( N != 0UL ) {
         do {
            randomize( A );
         }
         while( isDefault( det( A ) ) );

         randomize( B );
      }

      const DynamicMatrix<Type,rowMajor>    A1( A );
      const DynamicMatrix<Type,columnMajor> A2( A );

      const DynamicMatrix<Type,rowMajor>    B1( B );
      const DynamicMatrix<Type,columnMajor> B2( B );

      DynamicMatrix<Type,rowMajor>    X1;
      DynamicMatrix<Type,columnMajor> X2;

      solve( declsym( trans(A1) ), X1, trans(B1) );
      solve( declsym( trans(A2) ), X2, trans(B2) );

      if( trans(A)*X1 != trans(B) || trans(A)*X2 != trans(B) || X1 != X2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << trans(A) << "\n"
             << "   Right-hand sides (B):\n" << trans(B) << "\n"
             << "   Row-major solutions (X1):\n" << X1 << "\n"
             << "   Column-major solutions (X2):\n" << X2 << "\n"
             << "   A * X1 =\n" << ( trans(A) * X1 ) << "\n"
             << "   A * X2 =\n" << ( trans(A) * X2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Symmetric LSE (multiple rhs, non-matching right-hand side)";

      SymmetricMatrix< DynamicMatrix<Type> > A( 2UL );
      DynamicMatrix<Type> B( 3UL, 5UL );
      DynamicMatrix<Type> X;

      try {
         X = solve( A, B );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE with invalid right-hand side succeeded\n"
             << " Details:\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Right-hand side (B):\n" << B << "\n"
             << "   Solution (X):\n" << X << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the LSE kernels with random \f$ N \times N \f$ Hermitian matrices.
//
// \param N The number of rows and columns of the matrix.
// \return void
// \exception std::runtime_error Error detected.
//
// This function tests the dense matrix LSE kernels for random \f$ N \times N \f$ Hermitian
// matrices. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void DenseTest::testHermitian( size_t N )
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   using blaze::HermitianMatrix;
   using blaze::DynamicMatrix;
   using blaze::DynamicVector;
   using blaze::rowMajor;
   using blaze::columnMajor;
   using blaze::solve;
   using blaze::isDefault;


   //=====================================================================================
   // Single right-hand side
   //=====================================================================================

   {
      test_ = "Hermitian LSE (single rhs, automatic)";

      HermitianMatrix< DynamicMatrix<Type> > A( N );
      DynamicVector<Type> b( N );

      if( N != 0UL ) {
         do {
            randomize( A );
         }
         while( isDefault( det( A ) ) );

         randomize( b );
      }

      const HermitianMatrix< DynamicMatrix<Type,rowMajor>    > A1( A );
      const HermitianMatrix< DynamicMatrix<Type,columnMajor> > A2( A );

      DynamicVector<Type> x1( N );
      DynamicVector<Type> x2( N );

      x1 = solve( A1, b );
      x2 = solve( A2, b );

      if( A*x1 != b || A*x2 != b || x1 != x2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   Row-major solution (x1):\n" << x1 << "\n"
             << "   Column-major solution (x2):\n" << x2 << "\n"
             << "   A * x1 =\n" << ( A * x1 ) << "\n"
             << "   A * x2 =\n" << ( A * x2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Hermitian LSE (single rhs, automatic, transpose)";

      HermitianMatrix< DynamicMatrix<Type> > A( N );
      DynamicVector<Type> b( N );

      if( N != 0UL ) {
         do {
            randomize( A );
         }
         while( isDefault( det( A ) ) );

         randomize( b );
      }

      const HermitianMatrix< DynamicMatrix<Type,rowMajor>    > A1( A );
      const HermitianMatrix< DynamicMatrix<Type,columnMajor> > A2( A );

      DynamicVector<Type> x1( N );
      DynamicVector<Type> x2( N );

      x1 = solve( trans(A1), b );
      x2 = solve( trans(A2), b );

      if( trans(A)*x1 != b || trans(A)*x2 != b || x1 != x2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << trans(A) << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   Row-major solution (x1):\n" << x1 << "\n"
             << "   Column-major solution (x2):\n" << x2 << "\n"
             << "   A * x1 =\n" << ( trans(A) * x1 ) << "\n"
             << "   A * x2 =\n" << ( trans(A) * x2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Hermitian LSE (single rhs, declherm)";

      HermitianMatrix< DynamicMatrix<Type> > A( N );
      DynamicVector<Type> b( N );

      if( N != 0UL ) {
         do {
            randomize( A );
         }
         while( isDefault( det( A ) ) );

         randomize( b );
      }

      const DynamicMatrix<Type,rowMajor>    A1( A );
      const DynamicMatrix<Type,columnMajor> A2( A );

      DynamicVector<Type> x1( N );
      DynamicVector<Type> x2( N );

      solve( declherm( A1 ), x1, b );
      solve( declherm( A2 ), x2, b );

      if( A*x1 != b || A*x2 != b || x1 != x2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   Row-major solution (x1):\n" << x1 << "\n"
             << "   Column-major solution (x2):\n" << x2 << "\n"
             << "   A * x1 =\n" << ( A * x1 ) << "\n"
             << "   A * x2 =\n" << ( A * x2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Hermitian LSE (single rhs, declherm, transpose)";

      HermitianMatrix< DynamicMatrix<Type> > A( N );
      DynamicVector<Type> b( N );

      if( N != 0UL ) {
         do {
            randomize( A );
         }
         while( isDefault( det( A ) ) );

         randomize( b );
      }

      const DynamicMatrix<Type,rowMajor>    A1( A );
      const DynamicMatrix<Type,columnMajor> A2( A );

      DynamicVector<Type> x1( N );
      DynamicVector<Type> x2( N );

      solve( declherm( trans(A1) ), x1, b );
      solve( declherm( trans(A2) ), x2, b );

      if( trans(A)*x1 != b || trans(A)*x2 != b || x1 != x2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << trans(A) << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   Row-major solution (x1):\n" << x1 << "\n"
             << "   Column-major solution (x2):\n" << x2 << "\n"
             << "   A * x1 =\n" << ( trans(A) * x1 ) << "\n"
             << "   A * x2 =\n" << ( trans(A) * x2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "General LSE (single rhs, non-matching right-hand side)";

      const HermitianMatrix< DynamicMatrix<Type> > A( 2UL );
      const DynamicVector<Type> b( 3UL );
      DynamicVector<Type> x;

      try {
         x = solve( A, b );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE with invalid right-hand side succeeded\n"
             << " Details:\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   Solution (x):\n" << x << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Multiple right-hand sides
   //=====================================================================================

   {
      test_ = "Hermitian LSE (multiple rhs, automatic)";

      HermitianMatrix< DynamicMatrix<Type> > A( N );
      DynamicMatrix<Type> B( N, 3UL );

      if( N != 0UL ) {
         do {
            randomize( A );
         }
         while( isDefault( det( A ) ) );

         randomize( B );
      }

      const HermitianMatrix< DynamicMatrix<Type,rowMajor> >    A1( A );
      const HermitianMatrix< DynamicMatrix<Type,columnMajor> > A2( A );

      const DynamicMatrix<Type,rowMajor>    B1( B );
      const DynamicMatrix<Type,columnMajor> B2( B );

      DynamicMatrix<Type,rowMajor>    X1;
      DynamicMatrix<Type,columnMajor> X2;

      X1 = solve( A1, B1 );
      X2 = solve( A2, B2 );

      if( A*X1 != B || A*X2 != B || X1 != X2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Right-hand sides (B):\n" << B << "\n"
             << "   Row-major solutions (X1):\n" << X1 << "\n"
             << "   Column-major solutions (X2):\n" << X2 << "\n"
             << "   A * X1 =\n" << ( A * X1 ) << "\n"
             << "   A * X2 =\n" << ( A * X2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Hermitian LSE (multiple rhs, automatic, transpose)";

      HermitianMatrix< DynamicMatrix<Type> > A( N );
      DynamicMatrix<Type> B( 3UL, N );

      if( N != 0UL ) {
         do {
            randomize( A );
         }
         while( isDefault( det( A ) ) );

         randomize( B );
      }

      const HermitianMatrix< DynamicMatrix<Type,rowMajor> >    A1( A );
      const HermitianMatrix< DynamicMatrix<Type,columnMajor> > A2( A );

      const DynamicMatrix<Type,rowMajor>    B1( B );
      const DynamicMatrix<Type,columnMajor> B2( B );

      DynamicMatrix<Type,rowMajor>    X1;
      DynamicMatrix<Type,columnMajor> X2;

      X1 = solve( trans(A1), trans(B1) );
      X2 = solve( trans(A2), trans(B2) );

      if( trans(A)*X1 != trans(B) || trans(A)*X2 != trans(B) || X1 != X2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << trans(A) << "\n"
             << "   Right-hand sides (B):\n" << trans(B) << "\n"
             << "   Row-major solutions (X1):\n" << X1 << "\n"
             << "   Column-major solutions (X2):\n" << X2 << "\n"
             << "   A * X1 =\n" << ( trans(A) * X1 ) << "\n"
             << "   A * X2 =\n" << ( trans(A) * X2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Hermitian LSE (multiple rhs, declherm)";

      HermitianMatrix< DynamicMatrix<Type> > A( N );
      DynamicMatrix<Type> B( N, 3UL );

      if( N != 0UL ) {
         do {
            randomize( A );
         }
         while( isDefault( det( A ) ) );

         randomize( B );
      }

      const DynamicMatrix<Type,rowMajor>    A1( A );
      const DynamicMatrix<Type,columnMajor> A2( A );

      const DynamicMatrix<Type,rowMajor>    B1( B );
      const DynamicMatrix<Type,columnMajor> B2( B );

      DynamicMatrix<Type,rowMajor>    X1;
      DynamicMatrix<Type,columnMajor> X2;

      solve( declherm( A1 ), X1, B1 );
      solve( declherm( A2 ), X2, B2 );

      if( A*X1 != B || A*X2 != B || X1 != X2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Right-hand sides (B):\n" << B << "\n"
             << "   Row-major solutions (X1):\n" << X1 << "\n"
             << "   Column-major solutions (X2):\n" << X2 << "\n"
             << "   A * X1 =\n" << ( A * X1 ) << "\n"
             << "   A * X2 =\n" << ( A * X2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Hermitian LSE (multiple rhs, declherm, transpose)";

      HermitianMatrix< DynamicMatrix<Type> > A( N );
      DynamicMatrix<Type> B( 3UL, N );

      if( N != 0UL ) {
         do {
            randomize( A );
         }
         while( isDefault( det( A ) ) );

         randomize( B );
      }

      const DynamicMatrix<Type,rowMajor>    A1( A );
      const DynamicMatrix<Type,columnMajor> A2( A );

      const DynamicMatrix<Type,rowMajor>    B1( B );
      const DynamicMatrix<Type,columnMajor> B2( B );

      DynamicMatrix<Type,rowMajor>    X1;
      DynamicMatrix<Type,columnMajor> X2;

      solve( declherm( trans(A1) ), X1, trans(B1) );
      solve( declherm( trans(A2) ), X2, trans(B2) );

      if( trans(A)*X1 != trans(B) || trans(A)*X2 != trans(B) || X1 != X2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << trans(A) << "\n"
             << "   Right-hand sides (B):\n" << trans(B) << "\n"
             << "   Row-major solutions (X1):\n" << X1 << "\n"
             << "   Column-major solutions (X2):\n" << X2 << "\n"
             << "   A * X1 =\n" << ( trans(A) * X1 ) << "\n"
             << "   A * X2 =\n" << ( trans(A) * X2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Hermitian LSE (multiple rhs, non-matching right-hand side)";

      HermitianMatrix< DynamicMatrix<Type> > A( 2UL );
      DynamicMatrix<Type> B( 3UL, 5UL );
      DynamicMatrix<Type> X;

      try {
         X = solve( A, B );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE with invalid right-hand side succeeded\n"
             << " Details:\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Right-hand side (B):\n" << B << "\n"
             << "   Solution (X):\n" << X << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the LSE kernels with random \f$ N \times N \f$ lower triangular matrices.
//
// \param N The number of rows and columns of the matrix.
// \return void
// \exception std::runtime_error Error detected.
//
// This function tests the dense matrix LSE kernels for random \f$ N \times N \f$ lower triangular
// matrices. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void DenseTest::testLower( size_t N )
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   using blaze::LowerMatrix;
   using blaze::UpperMatrix;
   using blaze::DynamicMatrix;
   using blaze::DynamicVector;
   using blaze::rowMajor;
   using blaze::columnMajor;
   using blaze::solve;
   using blaze::isDefault;


   //=====================================================================================
   // Single right-hand side
   //=====================================================================================

   {
      test_ = "Lower LSE (single rhs, automatic)";

      LowerMatrix< DynamicMatrix<Type> > A( N );
      DynamicVector<Type> b( N );

      if( N != 0 ) {
         do {
            randomize( A );
         }
         while( isDefault( det( A ) ) );

         randomize( b );
      }

      const LowerMatrix< DynamicMatrix<Type,rowMajor>    > A1( A );
      const LowerMatrix< DynamicMatrix<Type,columnMajor> > A2( A );

      DynamicVector<Type> x1( N );
      DynamicVector<Type> x2( N );

      x1 = solve( A1, b );
      x2 = solve( A2, b );

      if( A*x1 != b || A*x2 != b || x1 != x2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   Row-major solution (x1):\n" << x1 << "\n"
             << "   Column-major solution (x2):\n" << x2 << "\n"
             << "   A * x1 =\n" << ( A * x1 ) << "\n"
             << "   A * x2 =\n" << ( A * x2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Lower LSE (single rhs, automatic, transpose)";

      UpperMatrix< DynamicMatrix<Type> > A( N );
      DynamicVector<Type> b( N );

      if( N != 0 ) {
         do {
            randomize( A );
         }
         while( isDefault( det( A ) ) );

         randomize( b );
      }

      const UpperMatrix< DynamicMatrix<Type,rowMajor>    > A1( A );
      const UpperMatrix< DynamicMatrix<Type,columnMajor> > A2( A );

      DynamicVector<Type> x1( N );
      DynamicVector<Type> x2( N );

      x1 = solve( trans(A1), b );
      x2 = solve( trans(A2), b );

      if( trans(A)*x1 != b || trans(A)*x2 != b || x1 != x2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << trans(A) << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   Row-major solution (x1):\n" << x1 << "\n"
             << "   Column-major solution (x2):\n" << x2 << "\n"
             << "   A * x1 =\n" << ( trans(A) * x1 ) << "\n"
             << "   A * x2 =\n" << ( trans(A) * x2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Lower LSE (single rhs, decllow)";

      LowerMatrix< DynamicMatrix<Type> > A( N );
      DynamicVector<Type> b( N );

      if( N != 0UL ) {
         do {
            randomize( A );
         }
         while( isDefault( det( A ) ) );

         randomize( b );
      }

      const DynamicMatrix<Type,rowMajor>    A1( A );
      const DynamicMatrix<Type,columnMajor> A2( A );

      DynamicVector<Type> x1( N );
      DynamicVector<Type> x2( N );

      solve( decllow( A1 ), x1, b );
      solve( decllow( A2 ), x2, b );

      if( A*x1 != b || A*x2 != b || x1 != x2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   Row-major solution (x1):\n" << x1 << "\n"
             << "   Column-major solution (x2):\n" << x2 << "\n"
             << "   A * x1 =\n" << ( A * x1 ) << "\n"
             << "   A * x2 =\n" << ( A * x2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Lower LSE (single rhs, decllow, transpose)";

      UpperMatrix< DynamicMatrix<Type> > A( N );
      DynamicVector<Type> b( N );

      if( N != 0UL ) {
         do {
            randomize( A );
         }
         while( isDefault( det( A ) ) );

         randomize( b );
      }

      const DynamicMatrix<Type,rowMajor>    A1( A );
      const DynamicMatrix<Type,columnMajor> A2( A );

      DynamicVector<Type> x1( N );
      DynamicVector<Type> x2( N );

      solve( decllow( trans(A1) ), x1, b );
      solve( decllow( trans(A2) ), x2, b );

      if( trans(A)*x1 != b || trans(A)*x2 != b || x1 != x2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << trans(A) << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   Row-major solution (x1):\n" << x1 << "\n"
             << "   Column-major solution (x2):\n" << x2 << "\n"
             << "   A * x1 =\n" << ( trans(A) * x1 ) << "\n"
             << "   A * x2 =\n" << ( trans(A) * x2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "General LSE (single rhs, non-matching right-hand side)";

      const LowerMatrix< DynamicMatrix<Type> > A( 2UL );
      const DynamicVector<Type> b( 3UL );
      DynamicVector<Type> x;

      try {
         x = solve( A, b );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE with invalid right-hand side succeeded\n"
             << " Details:\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   Solution (x):\n" << x << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Multiple right-hand sides
   //=====================================================================================

   {
      test_ = "Lower LSE (multiple rhs, automatic)";

      LowerMatrix< DynamicMatrix<Type> > A( N );
      DynamicMatrix<Type> B( N, 3UL );

      if( N != 0UL ) {
         do {
            randomize( A );
         }
         while( isDefault( det( A ) ) );

         randomize( B );
      }

      const LowerMatrix< DynamicMatrix<Type,rowMajor> >    A1( A );
      const LowerMatrix< DynamicMatrix<Type,columnMajor> > A2( A );

      const DynamicMatrix<Type,rowMajor>    B1( B );
      const DynamicMatrix<Type,columnMajor> B2( B );

      DynamicMatrix<Type,rowMajor>    X1;
      DynamicMatrix<Type,columnMajor> X2;

      X1 = solve( A1, B1 );
      X2 = solve( A2, B2 );

      if( A*X1 != B || A*X2 != B || X1 != X2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Right-hand sides (B):\n" << B << "\n"
             << "   Row-major solutions (X1):\n" << X1 << "\n"
             << "   Column-major solutions (X2):\n" << X2 << "\n"
             << "   A * X1 =\n" << ( A * X1 ) << "\n"
             << "   A * X2 =\n" << ( A * X2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Lower LSE (multiple rhs, automatic, transpose)";

      UpperMatrix< DynamicMatrix<Type> > A( N );
      DynamicMatrix<Type> B( 3UL, N );

      if( N != 0UL ) {
         do {
            randomize( A );
         }
         while( isDefault( det( A ) ) );

         randomize( B );
      }

      const UpperMatrix< DynamicMatrix<Type,rowMajor> >    A1( A );
      const UpperMatrix< DynamicMatrix<Type,columnMajor> > A2( A );

      const DynamicMatrix<Type,rowMajor>    B1( B );
      const DynamicMatrix<Type,columnMajor> B2( B );

      DynamicMatrix<Type,rowMajor>    X1;
      DynamicMatrix<Type,columnMajor> X2;

      X1 = solve( trans(A1), trans(B1) );
      X2 = solve( trans(A2), trans(B2) );

      if( trans(A)*X1 != trans(B) || trans(A)*X2 != trans(B) || X1 != X2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << trans(A) << "\n"
             << "   Right-hand sides (B):\n" << trans(B) << "\n"
             << "   Row-major solutions (X1):\n" << X1 << "\n"
             << "   Column-major solutions (X2):\n" << X2 << "\n"
             << "   A * X1 =\n" << ( trans(A) * X1 ) << "\n"
             << "   A * X2 =\n" << ( trans(A) * X2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Lower LSE (multiple rhs, decllow)";

      LowerMatrix< DynamicMatrix<Type> > A( N );
      DynamicMatrix<Type> B( N, 3UL );

      if( N != 0UL ) {
         do {
            randomize( A );
         }
         while( isDefault( det( A ) ) );

         randomize( B );
      }

      const DynamicMatrix<Type,rowMajor>    A1( A );
      const DynamicMatrix<Type,columnMajor> A2( A );

      const DynamicMatrix<Type,rowMajor>    B1( B );
      const DynamicMatrix<Type,columnMajor> B2( B );

      DynamicMatrix<Type,rowMajor>    X1;
      DynamicMatrix<Type,columnMajor> X2;

      solve( decllow( A1 ), X1, B1 );
      solve( decllow( A2 ), X2, B2 );

      if( A*X1 != B || A*X2 != B || X1 != X2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Right-hand sides (B):\n" << B << "\n"
             << "   Row-major solutions (X1):\n" << X1 << "\n"
             << "   Column-major solutions (X2):\n" << X2 << "\n"
             << "   A * X1 =\n" << ( A * X1 ) << "\n"
             << "   A * X2 =\n" << ( A * X2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Lower LSE (multiple rhs, decllow, transpose)";

      UpperMatrix< DynamicMatrix<Type> > A( N );
      DynamicMatrix<Type> B( 3UL, N );

      if( N != 0UL ) {
         do {
            randomize( A );
         }
         while( isDefault( det( A ) ) );

         randomize( B );
      }

      const DynamicMatrix<Type,rowMajor>    A1( A );
      const DynamicMatrix<Type,columnMajor> A2( A );

      const DynamicMatrix<Type,rowMajor>    B1( B );
      const DynamicMatrix<Type,columnMajor> B2( B );

      DynamicMatrix<Type,rowMajor>    X1;
      DynamicMatrix<Type,columnMajor> X2;

      solve( decllow( trans(A1) ), X1, trans(B1) );
      solve( decllow( trans(A2) ), X2, trans(B2) );

      if( trans(A)*X1 != trans(B) || trans(A)*X2 != trans(B) || X1 != X2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << trans(A) << "\n"
             << "   Right-hand sides (B):\n" << trans(B) << "\n"
             << "   Row-major solutions (X1):\n" << X1 << "\n"
             << "   Column-major solutions (X2):\n" << X2 << "\n"
             << "   A * X1 =\n" << ( trans(A) * X1 ) << "\n"
             << "   A * X2 =\n" << ( trans(A) * X2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Lower LSE (multiple rhs, non-matching right-hand side)";

      LowerMatrix< DynamicMatrix<Type> > A( 2UL );
      DynamicMatrix<Type> B( 3UL, 5UL );
      DynamicMatrix<Type> X;

      try {
         X = solve( A, B );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE with invalid right-hand side succeeded\n"
             << " Details:\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Right-hand side (B):\n" << B << "\n"
             << "   Solution (X):\n" << X << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the LSE kernels with random \f$ N \times N \f$ lower unitriangular matrices.
//
// \param N The number of rows and columns of the matrix.
// \return void
// \exception std::runtime_error Error detected.
//
// This function tests the dense matrix LSE kernels for random \f$ N \times N \f$ lower
// unitriangular matrices. In case an error is detected, a \a std::runtime_error exception is
// thrown.
*/
template< typename Type >
void DenseTest::testUniLower( size_t N )
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   using blaze::UniLowerMatrix;
   using blaze::UniUpperMatrix;
   using blaze::DynamicMatrix;
   using blaze::DynamicVector;
   using blaze::rowMajor;
   using blaze::columnMajor;
   using blaze::solve;
   using blaze::isDefault;


   //=====================================================================================
   // Single right-hand side
   //=====================================================================================

   {
      test_ = "UniLower LSE (single rhs, automatic)";

      UniLowerMatrix< DynamicMatrix<Type> > A( N );
      DynamicVector<Type> b( N );

      if( N != 0UL ) {
         do {
            randomize( A );
         }
         while( isDefault( det( A ) ) );

         randomize( b );
      }

      const UniLowerMatrix< DynamicMatrix<Type,rowMajor>    > A1( A );
      const UniLowerMatrix< DynamicMatrix<Type,columnMajor> > A2( A );

      DynamicVector<Type> x1( N );
      DynamicVector<Type> x2( N );

      x1 = solve( A1, b );
      x2 = solve( A2, b );

      if( A*x1 != b || A*x2 != b || x1 != x2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   Row-major solution (x1):\n" << x1 << "\n"
             << "   Column-major solution (x2):\n" << x2 << "\n"
             << "   A * x1 =\n" << ( A * x1 ) << "\n"
             << "   A * x2 =\n" << ( A * x2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "UniLower LSE (single rhs, automatic, transpose)";

      UniUpperMatrix< DynamicMatrix<Type> > A( N );
      DynamicVector<Type> b( N );

      if( N != 0UL ) {
         do {
            randomize( A );
         }
         while( isDefault( det( A ) ) );

         randomize( b );
      }

      const UniUpperMatrix< DynamicMatrix<Type,rowMajor>    > A1( A );
      const UniUpperMatrix< DynamicMatrix<Type,columnMajor> > A2( A );

      DynamicVector<Type> x1( N );
      DynamicVector<Type> x2( N );

      x1 = solve( trans(A1), b );
      x2 = solve( trans(A2), b );

      if( trans(A)*x1 != b || trans(A)*x2 != b || x1 != x2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << trans(A) << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   Row-major solution (x1):\n" << x1 << "\n"
             << "   Column-major solution (x2):\n" << x2 << "\n"
             << "   A * x1 =\n" << ( trans(A) * x1 ) << "\n"
             << "   A * x2 =\n" << ( trans(A) * x2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "UniLower LSE (single rhs, declunilow)";

      UniLowerMatrix< DynamicMatrix<Type> > A( N );
      DynamicVector<Type> b( N );

      if( N != 0UL ) {
         do {
            randomize( A );
         }
         while( isDefault( det( A ) ) );

         randomize( b );
      }

      const DynamicMatrix<Type,rowMajor>    A1( A );
      const DynamicMatrix<Type,columnMajor> A2( A );

      DynamicVector<Type> x1( N );
      DynamicVector<Type> x2( N );

      solve( declunilow( A1 ), x1, b );
      solve( declunilow( A2 ), x2, b );

      if( A*x1 != b || A*x2 != b || x1 != x2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   Row-major solution (x1):\n" << x1 << "\n"
             << "   Column-major solution (x2):\n" << x2 << "\n"
             << "   A * x1 =\n" << ( A * x1 ) << "\n"
             << "   A * x2 =\n" << ( A * x2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "UniLower LSE (single rhs, declunilow, transpose)";

      UniUpperMatrix< DynamicMatrix<Type> > A( N );
      DynamicVector<Type> b( N );

      if( N != 0UL ) {
         do {
            randomize( A );
         }
         while( isDefault( det( A ) ) );

         randomize( b );
      }

      const DynamicMatrix<Type,rowMajor>    A1( A );
      const DynamicMatrix<Type,columnMajor> A2( A );

      DynamicVector<Type> x1( N );
      DynamicVector<Type> x2( N );

      solve( declunilow( trans(A1) ), x1, b );
      solve( declunilow( trans(A2) ), x2, b );

      if( trans(A)*x1 != b || trans(A)*x2 != b || x1 != x2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << trans(A) << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   Row-major solution (x1):\n" << x1 << "\n"
             << "   Column-major solution (x2):\n" << x2 << "\n"
             << "   A * x1 =\n" << ( trans(A) * x1 ) << "\n"
             << "   A * x2 =\n" << ( trans(A) * x2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "General LSE (single rhs, non-matching right-hand side)";

      const UniLowerMatrix< DynamicMatrix<Type> > A( 2UL );
      const DynamicVector<Type> b( 3UL );
      DynamicVector<Type> x;

      try {
         x = solve( A, b );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE with invalid right-hand side succeeded\n"
             << " Details:\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   Solution (x):\n" << x << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Multiple right-hand sides
   //=====================================================================================

   {
      test_ = "UniLower LSE (multiple rhs, automatic)";

      UniLowerMatrix< DynamicMatrix<Type> > A( N );
      DynamicMatrix<Type> B( N, 3UL );

      if( N != 0UL ) {
         do {
            randomize( A );
         }
         while( isDefault( det( A ) ) );

         randomize( B );
      }

      const UniLowerMatrix< DynamicMatrix<Type,rowMajor> >    A1( A );
      const UniLowerMatrix< DynamicMatrix<Type,columnMajor> > A2( A );

      const DynamicMatrix<Type,rowMajor>    B1( B );
      const DynamicMatrix<Type,columnMajor> B2( B );

      DynamicMatrix<Type,rowMajor>    X1;
      DynamicMatrix<Type,columnMajor> X2;

      X1 = solve( A1, B1 );
      X2 = solve( A2, B2 );

      if( A*X1 != B || A*X2 != B || X1 != X2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Right-hand sides (B):\n" << B << "\n"
             << "   Row-major solutions (X1):\n" << X1 << "\n"
             << "   Column-major solutions (X2):\n" << X2 << "\n"
             << "   A * X1 =\n" << ( A * X1 ) << "\n"
             << "   A * X2 =\n" << ( A * X2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "UniLower LSE (multiple rhs, automatic, transpose)";

      UniUpperMatrix< DynamicMatrix<Type> > A( N );
      DynamicMatrix<Type> B( 3UL, N );

      if( N != 0UL ) {
         do {
            randomize( A );
         }
         while( isDefault( det( A ) ) );

         randomize( B );
      }

      const UniUpperMatrix< DynamicMatrix<Type,rowMajor> >    A1( A );
      const UniUpperMatrix< DynamicMatrix<Type,columnMajor> > A2( A );

      const DynamicMatrix<Type,rowMajor>    B1( B );
      const DynamicMatrix<Type,columnMajor> B2( B );

      DynamicMatrix<Type,rowMajor>    X1;
      DynamicMatrix<Type,columnMajor> X2;

      X1 = solve( trans(A1), trans(B1) );
      X2 = solve( trans(A2), trans(B2) );

      if( trans(A)*X1 != trans(B) || trans(A)*X2 != trans(B) || X1 != X2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << trans(A) << "\n"
             << "   Right-hand sides (B):\n" << trans(B) << "\n"
             << "   Row-major solutions (X1):\n" << X1 << "\n"
             << "   Column-major solutions (X2):\n" << X2 << "\n"
             << "   A * X1 =\n" << ( trans(A) * X1 ) << "\n"
             << "   A * X2 =\n" << ( trans(A) * X2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "UniLower LSE (multiple rhs, declunilow)";

      UniLowerMatrix< DynamicMatrix<Type> > A( N );
      DynamicMatrix<Type> B( N, 3UL );

      if( N != 0UL ) {
         do {
            randomize( A );
         }
         while( isDefault( det( A ) ) );

         randomize( B );
      }

      const DynamicMatrix<Type,rowMajor>    A1( A );
      const DynamicMatrix<Type,columnMajor> A2( A );

      const DynamicMatrix<Type,rowMajor>    B1( B );
      const DynamicMatrix<Type,columnMajor> B2( B );

      DynamicMatrix<Type,rowMajor>    X1;
      DynamicMatrix<Type,columnMajor> X2;

      solve( declunilow( A1 ), X1, B1 );
      solve( declunilow( A2 ), X2, B2 );

      if( A*X1 != B || A*X2 != B || X1 != X2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Right-hand sides (B):\n" << B << "\n"
             << "   Row-major solutions (X1):\n" << X1 << "\n"
             << "   Column-major solutions (X2):\n" << X2 << "\n"
             << "   A * X1 =\n" << ( A * X1 ) << "\n"
             << "   A * X2 =\n" << ( A * X2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "UniLower LSE (multiple rhs, declunilow, transpose)";

      UniUpperMatrix< DynamicMatrix<Type> > A( N );
      DynamicMatrix<Type> B( 3UL, N );

      if( N != 0UL ) {
         do {
            randomize( A );
         }
         while( isDefault( det( A ) ) );

         randomize( B );
      }

      const DynamicMatrix<Type,rowMajor>    A1( A );
      const DynamicMatrix<Type,columnMajor> A2( A );

      const DynamicMatrix<Type,rowMajor>    B1( B );
      const DynamicMatrix<Type,columnMajor> B2( B );

      DynamicMatrix<Type,rowMajor>    X1;
      DynamicMatrix<Type,columnMajor> X2;

      solve( declunilow( trans(A1) ), X1, trans(B1) );
      solve( declunilow( trans(A2) ), X2, trans(B2) );

      if( trans(A)*X1 != trans(B) || trans(A)*X2 != trans(B) || X1 != X2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << trans(A) << "\n"
             << "   Right-hand sides (B):\n" << trans(B) << "\n"
             << "   Row-major solutions (X1):\n" << X1 << "\n"
             << "   Column-major solutions (X2):\n" << X2 << "\n"
             << "   A * X1 =\n" << ( trans(A) * X1 ) << "\n"
             << "   A * X2 =\n" << ( trans(A) * X2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "UniLower LSE (multiple rhs, non-matching right-hand side)";

      UniLowerMatrix< DynamicMatrix<Type> > A( 2UL );
      DynamicMatrix<Type> B( 3UL, 5UL );
      DynamicMatrix<Type> X;

      try {
         X = solve( A, B );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE with invalid right-hand side succeeded\n"
             << " Details:\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Right-hand side (B):\n" << B << "\n"
             << "   Solution (X):\n" << X << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the LSE kernels with random \f$ N \times N \f$ upper triangular matrices.
//
// \param N The number of rows and columns of the matrix.
// \return void
// \exception std::runtime_error Error detected.
//
// This function tests the dense matrix LSE kernels for random \f$ N \times N \f$ upper triangular
// matrices. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void DenseTest::testUpper( size_t N )
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   using blaze::UpperMatrix;
   using blaze::LowerMatrix;
   using blaze::DynamicMatrix;
   using blaze::DynamicVector;
   using blaze::rowMajor;
   using blaze::columnMajor;
   using blaze::solve;
   using blaze::isDefault;


   //=====================================================================================
   // Single right-hand side
   //=====================================================================================

   {
      test_ = "Upper LSE (single rhs, automatic)";

      UpperMatrix< DynamicMatrix<Type> > A( N );
      DynamicVector<Type> b( N );

      if( N != 0 ) {
         do {
            randomize( A );
         }
         while( isDefault( det( A ) ) );

         randomize( b );
      }

      const UpperMatrix< DynamicMatrix<Type,rowMajor>    > A1( A );
      const UpperMatrix< DynamicMatrix<Type,columnMajor> > A2( A );

      DynamicVector<Type> x1( N );
      DynamicVector<Type> x2( N );

      x1 = solve( A1, b );
      x2 = solve( A2, b );

      if( A*x1 != b || A*x2 != b || x1 != x2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   Row-major solution (x1):\n" << x1 << "\n"
             << "   Column-major solution (x2):\n" << x2 << "\n"
             << "   A * x1 =\n" << ( A * x1 ) << "\n"
             << "   A * x2 =\n" << ( A * x2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Upper LSE (single rhs, automatic, transpose)";

      LowerMatrix< DynamicMatrix<Type> > A( N );
      DynamicVector<Type> b( N );

      if( N != 0 ) {
         do {
            randomize( A );
         }
         while( isDefault( det( A ) ) );

         randomize( b );
      }

      const LowerMatrix< DynamicMatrix<Type,rowMajor>    > A1( A );
      const LowerMatrix< DynamicMatrix<Type,columnMajor> > A2( A );

      DynamicVector<Type> x1( N );
      DynamicVector<Type> x2( N );

      x1 = solve( trans(A1), b );
      x2 = solve( trans(A2), b );

      if( trans(A)*x1 != b || trans(A)*x2 != b || x1 != x2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << trans(A) << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   Row-major solution (x1):\n" << x1 << "\n"
             << "   Column-major solution (x2):\n" << x2 << "\n"
             << "   A * x1 =\n" << ( trans(A) * x1 ) << "\n"
             << "   A * x2 =\n" << ( trans(A) * x2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Upper LSE (single rhs, declupp)";

      UpperMatrix< DynamicMatrix<Type> > A( N );
      DynamicVector<Type> b( N );

      if( N != 0UL ) {
         do {
            randomize( A );
         }
         while( isDefault( det( A ) ) );

         randomize( b );
      }

      const DynamicMatrix<Type,rowMajor>    A1( A );
      const DynamicMatrix<Type,columnMajor> A2( A );

      DynamicVector<Type> x1( N );
      DynamicVector<Type> x2( N );

      solve( declupp( A1 ), x1, b );
      solve( declupp( A2 ), x2, b );

      if( A*x1 != b || A*x2 != b || x1 != x2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   Row-major solution (x1):\n" << x1 << "\n"
             << "   Column-major solution (x2):\n" << x2 << "\n"
             << "   A * x1 =\n" << ( A * x1 ) << "\n"
             << "   A * x2 =\n" << ( A * x2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Upper LSE (single rhs, declupp, transpose)";

      LowerMatrix< DynamicMatrix<Type> > A( N );
      DynamicVector<Type> b( N );

      if( N != 0UL ) {
         do {
            randomize( A );
         }
         while( isDefault( det( A ) ) );

         randomize( b );
      }

      const DynamicMatrix<Type,rowMajor>    A1( A );
      const DynamicMatrix<Type,columnMajor> A2( A );

      DynamicVector<Type> x1( N );
      DynamicVector<Type> x2( N );

      solve( declupp( trans(A1) ), x1, b );
      solve( declupp( trans(A2) ), x2, b );

      if( trans(A)*x1 != b || trans(A)*x2 != b || x1 != x2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << trans(A) << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   Row-major solution (x1):\n" << x1 << "\n"
             << "   Column-major solution (x2):\n" << x2 << "\n"
             << "   A * x1 =\n" << ( trans(A) * x1 ) << "\n"
             << "   A * x2 =\n" << ( trans(A) * x2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "General LSE (single rhs, non-matching right-hand side)";

      const UpperMatrix< DynamicMatrix<Type> > A( 2UL );
      const DynamicVector<Type> b( 3UL );
      DynamicVector<Type> x;

      try {
         x = solve( A, b );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE with invalid right-hand side succeeded\n"
             << " Details:\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   Solution (x):\n" << x << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Multiple right-hand sides
   //=====================================================================================

   {
      test_ = "Upper LSE (multiple rhs, automatic)";

      UpperMatrix< DynamicMatrix<Type> > A( N );
      DynamicMatrix<Type> B( N, 3UL );

      if( N != 0UL ) {
         do {
            randomize( A );
         }
         while( isDefault( det( A ) ) );

         randomize( B );
      }

      const UpperMatrix< DynamicMatrix<Type,rowMajor> >    A1( A );
      const UpperMatrix< DynamicMatrix<Type,columnMajor> > A2( A );

      const DynamicMatrix<Type,rowMajor>    B1( B );
      const DynamicMatrix<Type,columnMajor> B2( B );

      DynamicMatrix<Type,rowMajor>    X1;
      DynamicMatrix<Type,columnMajor> X2;

      X1 = solve( A1, B1 );
      X2 = solve( A2, B2 );

      if( A*X1 != B || A*X2 != B || X1 != X2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Right-hand sides (B):\n" << B << "\n"
             << "   Row-major solutions (X1):\n" << X1 << "\n"
             << "   Column-major solutions (X2):\n" << X2 << "\n"
             << "   A * X1 =\n" << ( A * X1 ) << "\n"
             << "   A * X2 =\n" << ( A * X2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Upper LSE (multiple rhs, automatic, transpose)";

      LowerMatrix< DynamicMatrix<Type> > A( N );
      DynamicMatrix<Type> B( 3UL, N );

      if( N != 0UL ) {
         do {
            randomize( A );
         }
         while( isDefault( det( A ) ) );

         randomize( B );
      }

      const LowerMatrix< DynamicMatrix<Type,rowMajor> >    A1( A );
      const LowerMatrix< DynamicMatrix<Type,columnMajor> > A2( A );

      const DynamicMatrix<Type,rowMajor>    B1( B );
      const DynamicMatrix<Type,columnMajor> B2( B );

      DynamicMatrix<Type,rowMajor>    X1;
      DynamicMatrix<Type,columnMajor> X2;

      X1 = solve( trans(A1), trans(B1) );
      X2 = solve( trans(A2), trans(B2) );

      if( trans(A)*X1 != trans(B) || trans(A)*X2 != trans(B) || X1 != X2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << trans(A) << "\n"
             << "   Right-hand sides (B):\n" << trans(B) << "\n"
             << "   Row-major solutions (X1):\n" << X1 << "\n"
             << "   Column-major solutions (X2):\n" << X2 << "\n"
             << "   A * X1 =\n" << ( trans(A) * X1 ) << "\n"
             << "   A * X2 =\n" << ( trans(A) * X2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Upper LSE (multiple rhs, declupp)";

      UpperMatrix< DynamicMatrix<Type> > A( N );
      DynamicMatrix<Type> B( N, 3UL );

      if( N != 0UL ) {
         do {
            randomize( A );
         }
         while( isDefault( det( A ) ) );

         randomize( B );
      }

      const DynamicMatrix<Type,rowMajor>    A1( A );
      const DynamicMatrix<Type,columnMajor> A2( A );

      const DynamicMatrix<Type,rowMajor>    B1( B );
      const DynamicMatrix<Type,columnMajor> B2( B );

      DynamicMatrix<Type,rowMajor>    X1;
      DynamicMatrix<Type,columnMajor> X2;

      solve( declupp( A1 ), X1, B1 );
      solve( declupp( A2 ), X2, B2 );

      if( A*X1 != B || A*X2 != B || X1 != X2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Right-hand sides (B):\n" << B << "\n"
             << "   Row-major solutions (X1):\n" << X1 << "\n"
             << "   Column-major solutions (X2):\n" << X2 << "\n"
             << "   A * X1 =\n" << ( A * X1 ) << "\n"
             << "   A * X2 =\n" << ( A * X2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Upper LSE (multiple rhs, declupp, transpose)";

      LowerMatrix< DynamicMatrix<Type> > A( N );
      DynamicMatrix<Type> B( 3UL, N );

      if( N != 0UL ) {
         do {
            randomize( A );
         }
         while( isDefault( det( A ) ) );

         randomize( B );
      }

      const DynamicMatrix<Type,rowMajor>    A1( A );
      const DynamicMatrix<Type,columnMajor> A2( A );

      const DynamicMatrix<Type,rowMajor>    B1( B );
      const DynamicMatrix<Type,columnMajor> B2( B );

      DynamicMatrix<Type,rowMajor>    X1;
      DynamicMatrix<Type,columnMajor> X2;

      solve( declupp( trans(A1) ), X1, trans(B1) );
      solve( declupp( trans(A2) ), X2, trans(B2) );

      if( trans(A)*X1 != trans(B) || trans(A)*X2 != trans(B) || X1 != X2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << trans(A) << "\n"
             << "   Right-hand sides (B):\n" << trans(B) << "\n"
             << "   Row-major solutions (X1):\n" << X1 << "\n"
             << "   Column-major solutions (X2):\n" << X2 << "\n"
             << "   A * X1 =\n" << ( trans(A) * X1 ) << "\n"
             << "   A * X2 =\n" << ( trans(A) * X2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Upper LSE (multiple rhs, non-matching right-hand side)";

      UpperMatrix< DynamicMatrix<Type> > A( 2UL );
      DynamicMatrix<Type> B( 3UL, 5UL );
      DynamicMatrix<Type> X;

      try {
         X = solve( A, B );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE with invalid right-hand side succeeded\n"
             << " Details:\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Right-hand side (B):\n" << B << "\n"
             << "   Solution (X):\n" << X << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the LSE kernels with random \f$ N \times N \f$ upper unitriangular matrices.
//
// \param N The number of rows and columns of the matrix.
// \return void
// \exception std::runtime_error Error detected.
//
// This function tests the dense matrix LSE kernels for random \f$ N \times N \f$ upper
// unitriangular matrices. In case an error is detected, a \a std::runtime_error exception is
// thrown.
*/
template< typename Type >
void DenseTest::testUniUpper( size_t N )
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   using blaze::UniUpperMatrix;
   using blaze::UniLowerMatrix;
   using blaze::DynamicMatrix;
   using blaze::DynamicVector;
   using blaze::rowMajor;
   using blaze::columnMajor;
   using blaze::solve;
   using blaze::isDefault;


   //=====================================================================================
   // Single right-hand side
   //=====================================================================================

   {
      test_ = "UniUpper LSE (single rhs, automatic)";

      UniUpperMatrix< DynamicMatrix<Type> > A( N );
      DynamicVector<Type> b( N );

      if( N != 0UL ) {
         do {
            randomize( A );
         }
         while( isDefault( det( A ) ) );

         randomize( b );
      }

      const UniUpperMatrix< DynamicMatrix<Type,rowMajor>    > A1( A );
      const UniUpperMatrix< DynamicMatrix<Type,columnMajor> > A2( A );

      DynamicVector<Type> x1( N );
      DynamicVector<Type> x2( N );

      x1 = solve( A1, b );
      x2 = solve( A2, b );

      if( A*x1 != b || A*x2 != b || x1 != x2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   Row-major solution (x1):\n" << x1 << "\n"
             << "   Column-major solution (x2):\n" << x2 << "\n"
             << "   A * x1 =\n" << ( A * x1 ) << "\n"
             << "   A * x2 =\n" << ( A * x2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "UniUpper LSE (single rhs, automatic, transpose)";

      UniLowerMatrix< DynamicMatrix<Type> > A( N );
      DynamicVector<Type> b( N );

      if( N != 0UL ) {
         do {
            randomize( A );
         }
         while( isDefault( det( A ) ) );

         randomize( b );
      }

      const UniLowerMatrix< DynamicMatrix<Type,rowMajor>    > A1( A );
      const UniLowerMatrix< DynamicMatrix<Type,columnMajor> > A2( A );

      DynamicVector<Type> x1( N );
      DynamicVector<Type> x2( N );

      x1 = solve( trans(A1), b );
      x2 = solve( trans(A2), b );

      if( trans(A)*x1 != b || trans(A)*x2 != b || x1 != x2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << trans(A) << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   Row-major solution (x1):\n" << x1 << "\n"
             << "   Column-major solution (x2):\n" << x2 << "\n"
             << "   A * x1 =\n" << ( trans(A) * x1 ) << "\n"
             << "   A * x2 =\n" << ( trans(A) * x2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "UniUpper LSE (single rhs, decluniupp)";

      UniUpperMatrix< DynamicMatrix<Type> > A( N );
      DynamicVector<Type> b( N );

      if( N != 0UL ) {
         do {
            randomize( A );
         }
         while( isDefault( det( A ) ) );

         randomize( b );
      }

      const DynamicMatrix<Type,rowMajor>    A1( A );
      const DynamicMatrix<Type,columnMajor> A2( A );

      DynamicVector<Type> x1( N );
      DynamicVector<Type> x2( N );

      solve( decluniupp( A1 ), x1, b );
      solve( decluniupp( A2 ), x2, b );

      if( A*x1 != b || A*x2 != b || x1 != x2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   Row-major solution (x1):\n" << x1 << "\n"
             << "   Column-major solution (x2):\n" << x2 << "\n"
             << "   A * x1 =\n" << ( A * x1 ) << "\n"
             << "   A * x2 =\n" << ( A * x2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "UniUpper LSE (single rhs, decluniupp, transpose)";

      UniLowerMatrix< DynamicMatrix<Type> > A( N );
      DynamicVector<Type> b( N );

      if( N != 0UL ) {
         do {
            randomize( A );
         }
         while( isDefault( det( A ) ) );

         randomize( b );
      }

      const DynamicMatrix<Type,rowMajor>    A1( A );
      const DynamicMatrix<Type,columnMajor> A2( A );

      DynamicVector<Type> x1( N );
      DynamicVector<Type> x2( N );

      solve( decluniupp( trans(A1) ), x1, b );
      solve( decluniupp( trans(A2) ), x2, b );

      if( trans(A)*x1 != b || trans(A)*x2 != b || x1 != x2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << trans(A) << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   Row-major solution (x1):\n" << x1 << "\n"
             << "   Column-major solution (x2):\n" << x2 << "\n"
             << "   A * x1 =\n" << ( trans(A) * x1 ) << "\n"
             << "   A * x2 =\n" << ( trans(A) * x2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "General LSE (single rhs, non-matching right-hand side)";

      const UniUpperMatrix< DynamicMatrix<Type> > A( 2UL );
      const DynamicVector<Type> b( 3UL );
      DynamicVector<Type> x;

      try {
         x = solve( A, b );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE with invalid right-hand side succeeded\n"
             << " Details:\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   Solution (x):\n" << x << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Multiple right-hand sides
   //=====================================================================================

   {
      test_ = "UniUpper LSE (multiple rhs, automatic)";

      UniUpperMatrix< DynamicMatrix<Type> > A( N );
      DynamicMatrix<Type> B( N, 3UL );

      if( N != 0UL ) {
         do {
            randomize( A );
         }
         while( isDefault( det( A ) ) );

         randomize( B );
      }

      const UniUpperMatrix< DynamicMatrix<Type,rowMajor> >    A1( A );
      const UniUpperMatrix< DynamicMatrix<Type,columnMajor> > A2( A );

      const DynamicMatrix<Type,rowMajor>    B1( B );
      const DynamicMatrix<Type,columnMajor> B2( B );

      DynamicMatrix<Type,rowMajor>    X1;
      DynamicMatrix<Type,columnMajor> X2;

      X1 = solve( A1, B1 );
      X2 = solve( A2, B2 );

      if( A*X1 != B || A*X2 != B || X1 != X2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Right-hand sides (B):\n" << B << "\n"
             << "   Row-major solutions (X1):\n" << X1 << "\n"
             << "   Column-major solutions (X2):\n" << X2 << "\n"
             << "   A * X1 =\n" << ( A * X1 ) << "\n"
             << "   A * X2 =\n" << ( A * X2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "UniUpper LSE (multiple rhs, automatic, transpose)";

      UniLowerMatrix< DynamicMatrix<Type> > A( N );
      DynamicMatrix<Type> B( 3UL, N );

      if( N != 0UL ) {
         do {
            randomize( A );
         }
         while( isDefault( det( A ) ) );

         randomize( B );
      }

      const UniLowerMatrix< DynamicMatrix<Type,rowMajor> >    A1( A );
      const UniLowerMatrix< DynamicMatrix<Type,columnMajor> > A2( A );

      const DynamicMatrix<Type,rowMajor>    B1( B );
      const DynamicMatrix<Type,columnMajor> B2( B );

      DynamicMatrix<Type,rowMajor>    X1;
      DynamicMatrix<Type,columnMajor> X2;

      X1 = solve( trans(A1), trans(B1) );
      X2 = solve( trans(A2), trans(B2) );

      if( trans(A)*X1 != trans(B) || trans(A)*X2 != trans(B) || X1 != X2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << trans(A) << "\n"
             << "   Right-hand sides (B):\n" << trans(B) << "\n"
             << "   Row-major solutions (X1):\n" << X1 << "\n"
             << "   Column-major solutions (X2):\n" << X2 << "\n"
             << "   A * X1 =\n" << ( trans(A) * X1 ) << "\n"
             << "   A * X2 =\n" << ( trans(A) * X2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "UniUpper LSE (multiple rhs, decluniupp)";

      UniUpperMatrix< DynamicMatrix<Type> > A( N );
      DynamicMatrix<Type> B( N, 3UL );

      if( N != 0UL ) {
         do {
            randomize( A );
         }
         while( isDefault( det( A ) ) );

         randomize( B );
      }

      const DynamicMatrix<Type,rowMajor>    A1( A );
      const DynamicMatrix<Type,columnMajor> A2( A );

      const DynamicMatrix<Type,rowMajor>    B1( B );
      const DynamicMatrix<Type,columnMajor> B2( B );

      DynamicMatrix<Type,rowMajor>    X1;
      DynamicMatrix<Type,columnMajor> X2;

      solve( decluniupp( A1 ), X1, B1 );
      solve( decluniupp( A2 ), X2, B2 );

      if( A*X1 != B || A*X2 != B || X1 != X2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Right-hand sides (B):\n" << B << "\n"
             << "   Row-major solutions (X1):\n" << X1 << "\n"
             << "   Column-major solutions (X2):\n" << X2 << "\n"
             << "   A * X1 =\n" << ( A * X1 ) << "\n"
             << "   A * X2 =\n" << ( A * X2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "UniUpper LSE (multiple rhs, decluniupp, transpose)";

      UniLowerMatrix< DynamicMatrix<Type> > A( N );
      DynamicMatrix<Type> B( 3UL, N );

      if( N != 0UL ) {
         do {
            randomize( A );
         }
         while( isDefault( det( A ) ) );

         randomize( B );
      }

      const DynamicMatrix<Type,rowMajor>    A1( A );
      const DynamicMatrix<Type,columnMajor> A2( A );

      const DynamicMatrix<Type,rowMajor>    B1( B );
      const DynamicMatrix<Type,columnMajor> B2( B );

      DynamicMatrix<Type,rowMajor>    X1;
      DynamicMatrix<Type,columnMajor> X2;

      solve( decluniupp( trans(A1) ), X1, trans(B1) );
      solve( decluniupp( trans(A2) ), X2, trans(B2) );

      if( trans(A)*X1 != trans(B) || trans(A)*X2 != trans(B) || X1 != X2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << trans(A) << "\n"
             << "   Right-hand sides (B):\n" << trans(B) << "\n"
             << "   Row-major solutions (X1):\n" << X1 << "\n"
             << "   Column-major solutions (X2):\n" << X2 << "\n"
             << "   A * X1 =\n" << ( trans(A) * X1 ) << "\n"
             << "   A * X2 =\n" << ( trans(A) * X2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "UniUpper LSE (multiple rhs, non-matching right-hand side)";

      UniUpperMatrix< DynamicMatrix<Type> > A( 2UL );
      DynamicMatrix<Type> B( 3UL, 5UL );
      DynamicMatrix<Type> X;

      try {
         X = solve( A, B );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE with invalid right-hand side succeeded\n"
             << " Details:\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Right-hand side (B):\n" << B << "\n"
             << "   Solution (X):\n" << X << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the LSE kernels with random \f$ N \times N \f$ diagonal matrices.
//
// \param N The number of rows and columns of the matrix.
// \return void
// \exception std::runtime_error Error detected.
//
// This function tests the dense matrix LSE kernels for random \f$ N \times N \f$ diagonal
// matrices. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void DenseTest::testDiagonal( size_t N )
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   using blaze::DiagonalMatrix;
   using blaze::DynamicMatrix;
   using blaze::DynamicVector;
   using blaze::rowMajor;
   using blaze::columnMajor;
   using blaze::solve;
   using blaze::isDefault;


   //=====================================================================================
   // Single right-hand side
   //=====================================================================================

   {
      test_ = "Diagonal LSE (single rhs, automatic)";

      DiagonalMatrix< DynamicMatrix<Type> > A( N );
      DynamicVector<Type> b( N );

      if( N != 0 ) {
         do {
            randomize( A );
         }
         while( isDefault( det( A ) ) );

         randomize( b );
      }

      const DiagonalMatrix< DynamicMatrix<Type,rowMajor>    > A1( A );
      const DiagonalMatrix< DynamicMatrix<Type,columnMajor> > A2( A );

      DynamicVector<Type> x1( N );
      DynamicVector<Type> x2( N );

      x1 = solve( A1, b );
      x2 = solve( A2, b );

      if( A*x1 != b || A*x2 != b || x1 != x2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   Row-major solution (x1):\n" << x1 << "\n"
             << "   Column-major solution (x2):\n" << x2 << "\n"
             << "   A * x1 =\n" << ( A * x1 ) << "\n"
             << "   A * x2 =\n" << ( A * x2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Diagonal LSE (single rhs, automatic, transpose)";

      DiagonalMatrix< DynamicMatrix<Type> > A( N );
      DynamicVector<Type> b( N );

      if( N != 0 ) {
         do {
            randomize( A );
         }
         while( isDefault( det( A ) ) );

         randomize( b );
      }

      const DiagonalMatrix< DynamicMatrix<Type,rowMajor>    > A1( A );
      const DiagonalMatrix< DynamicMatrix<Type,columnMajor> > A2( A );

      DynamicVector<Type> x1( N );
      DynamicVector<Type> x2( N );

      x1 = solve( trans(A1), b );
      x2 = solve( trans(A2), b );

      if( trans(A)*x1 != b || trans(A)*x2 != b || x1 != x2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << trans(A) << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   Row-major solution (x1):\n" << x1 << "\n"
             << "   Column-major solution (x2):\n" << x2 << "\n"
             << "   A * x1 =\n" << ( trans(A) * x1 ) << "\n"
             << "   A * x2 =\n" << ( trans(A) * x2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Diagonal LSE (single rhs, decldiag)";

      DiagonalMatrix< DynamicMatrix<Type> > A( N );
      DynamicVector<Type> b( N );

      if( N != 0UL ) {
         do {
            randomize( A );
         }
         while( isDefault( det( A ) ) );

         randomize( b );
      }

      const DynamicMatrix<Type,rowMajor>    A1( A );
      const DynamicMatrix<Type,columnMajor> A2( A );

      DynamicVector<Type> x1( N );
      DynamicVector<Type> x2( N );

      solve( decldiag( A1 ), x1, b );
      solve( decldiag( A2 ), x2, b );

      if( A*x1 != b || A*x2 != b || x1 != x2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   Row-major solution (x1):\n" << x1 << "\n"
             << "   Column-major solution (x2):\n" << x2 << "\n"
             << "   A * x1 =\n" << ( A * x1 ) << "\n"
             << "   A * x2 =\n" << ( A * x2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Diagonal LSE (single rhs, decldiag, transpose)";

      DiagonalMatrix< DynamicMatrix<Type> > A( N );
      DynamicVector<Type> b( N );

      if( N != 0UL ) {
         do {
            randomize( A );
         }
         while( isDefault( det( A ) ) );

         randomize( b );
      }

      const DynamicMatrix<Type,rowMajor>    A1( A );
      const DynamicMatrix<Type,columnMajor> A2( A );

      DynamicVector<Type> x1( N );
      DynamicVector<Type> x2( N );

      solve( decldiag( trans(A1) ), x1, b );
      solve( decldiag( trans(A2) ), x2, b );

      if( trans(A)*x1 != b || trans(A)*x2 != b || x1 != x2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << trans(A) << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   Row-major solution (x1):\n" << x1 << "\n"
             << "   Column-major solution (x2):\n" << x2 << "\n"
             << "   A * x1 =\n" << ( trans(A) * x1 ) << "\n"
             << "   A * x2 =\n" << ( trans(A) * x2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "General LSE (single rhs, non-matching right-hand side)";

      const DiagonalMatrix< DynamicMatrix<Type> > A( 2UL );
      const DynamicVector<Type> b( 3UL );
      DynamicVector<Type> x;

      try {
         x = solve( A, b );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE with invalid right-hand side succeeded\n"
             << " Details:\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   Solution (x):\n" << x << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Multiple right-hand sides
   //=====================================================================================

   {
      test_ = "Diagonal LSE (multiple rhs, automatic)";

      DiagonalMatrix< DynamicMatrix<Type> > A( N );
      DynamicMatrix<Type> B( N, 3UL );

      if( N != 0UL ) {
         do {
            randomize( A );
         }
         while( isDefault( det( A ) ) );

         randomize( B );
      }

      const DiagonalMatrix< DynamicMatrix<Type,rowMajor> >    A1( A );
      const DiagonalMatrix< DynamicMatrix<Type,columnMajor> > A2( A );

      const DynamicMatrix<Type,rowMajor>    B1( B );
      const DynamicMatrix<Type,columnMajor> B2( B );

      DynamicMatrix<Type,rowMajor>    X1;
      DynamicMatrix<Type,columnMajor> X2;

      X1 = solve( A1, B1 );
      X2 = solve( A2, B2 );

      if( A*X1 != B || A*X2 != B || X1 != X2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Right-hand sides (B):\n" << B << "\n"
             << "   Row-major solutions (X1):\n" << X1 << "\n"
             << "   Column-major solutions (X2):\n" << X2 << "\n"
             << "   A * X1 =\n" << ( A * X1 ) << "\n"
             << "   A * X2 =\n" << ( A * X2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Diagonal LSE (multiple rhs, automatic, transpose)";

      DiagonalMatrix< DynamicMatrix<Type> > A( N );
      DynamicMatrix<Type> B( 3UL, N );

      if( N != 0UL ) {
         do {
            randomize( A );
         }
         while( isDefault( det( A ) ) );

         randomize( B );
      }

      const DiagonalMatrix< DynamicMatrix<Type,rowMajor> >    A1( A );
      const DiagonalMatrix< DynamicMatrix<Type,columnMajor> > A2( A );

      const DynamicMatrix<Type,rowMajor>    B1( B );
      const DynamicMatrix<Type,columnMajor> B2( B );

      DynamicMatrix<Type,rowMajor>    X1;
      DynamicMatrix<Type,columnMajor> X2;

      X1 = solve( trans(A1), trans(B1) );
      X2 = solve( trans(A2), trans(B2) );

      if( trans(A)*X1 != trans(B) || trans(A)*X2 != trans(B) || X1 != X2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << trans(A) << "\n"
             << "   Right-hand sides (B):\n" << trans(B) << "\n"
             << "   Row-major solutions (X1):\n" << X1 << "\n"
             << "   Column-major solutions (X2):\n" << X2 << "\n"
             << "   A * X1 =\n" << ( trans(A) * X1 ) << "\n"
             << "   A * X2 =\n" << ( trans(A) * X2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Diagonal LSE (multiple rhs, decldiag)";

      DiagonalMatrix< DynamicMatrix<Type> > A( N );
      DynamicMatrix<Type> B( N, 3UL );

      if( N != 0UL ) {
         do {
            randomize( A );
         }
         while( isDefault( det( A ) ) );

         randomize( B );
      }

      const DynamicMatrix<Type,rowMajor>    A1( A );
      const DynamicMatrix<Type,columnMajor> A2( A );

      const DynamicMatrix<Type,rowMajor>    B1( B );
      const DynamicMatrix<Type,columnMajor> B2( B );

      DynamicMatrix<Type,rowMajor>    X1;
      DynamicMatrix<Type,columnMajor> X2;

      solve( decldiag( A1 ), X1, B1 );
      solve( decldiag( A2 ), X2, B2 );

      if( A*X1 != B || A*X2 != B || X1 != X2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Right-hand sides (B):\n" << B << "\n"
             << "   Row-major solutions (X1):\n" << X1 << "\n"
             << "   Column-major solutions (X2):\n" << X2 << "\n"
             << "   A * X1 =\n" << ( A * X1 ) << "\n"
             << "   A * X2 =\n" << ( A * X2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Diagonal LSE (multiple rhs, decldiag, transpose)";

      DiagonalMatrix< DynamicMatrix<Type> > A( N );
      DynamicMatrix<Type> B( 3UL, N );

      if( N != 0UL ) {
         do {
            randomize( A );
         }
         while( isDefault( det( A ) ) );

         randomize( B );
      }

      const DynamicMatrix<Type,rowMajor>    A1( A );
      const DynamicMatrix<Type,columnMajor> A2( A );

      const DynamicMatrix<Type,rowMajor>    B1( B );
      const DynamicMatrix<Type,columnMajor> B2( B );

      DynamicMatrix<Type,rowMajor>    X1;
      DynamicMatrix<Type,columnMajor> X2;

      solve( decldiag( trans(A1) ), X1, trans(B1) );
      solve( decldiag( trans(A2) ), X2, trans(B2) );

      if( trans(A)*X1 != trans(B) || trans(A)*X2 != trans(B) || X1 != X2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << trans(A) << "\n"
             << "   Right-hand sides (B):\n" << trans(B) << "\n"
             << "   Row-major solutions (X1):\n" << X1 << "\n"
             << "   Column-major solutions (X2):\n" << X2 << "\n"
             << "   A * X1 =\n" << ( trans(A) * X1 ) << "\n"
             << "   A * X2 =\n" << ( trans(A) * X2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Diagonal LSE (multiple rhs, non-matching right-hand side)";

      DiagonalMatrix< DynamicMatrix<Type> > A( 2UL );
      DynamicMatrix<Type> B( 3UL, 5UL );
      DynamicMatrix<Type> X;

      try {
         X = solve( A, B );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving LSE with invalid right-hand side succeeded\n"
             << " Details:\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Right-hand side (B):\n" << B << "\n"
             << "   Solution (X):\n" << X << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

#endif
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Testing the dense matrix LSE kernels.
//
// \return void
*/
void runTest()
{
   DenseTest();
}
//*************************************************************************************************




//=================================================================================================
//
//  MACRO DEFINITIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Macro for the execution of the dense matrix LSE test.
*/
#define RUN_LSE_DENSE_TEST \
   blazetest::mathtest::operations::lse::runTest()
/*! \endcond */
//*************************************************************************************************

} // namespace lse

} // namespace operations

} // namespace mathtest

} // namespace blazetest

#endif
