//=================================================================================================
/*!
//  \file blazetest/mathtest/lapack/SolverTest.h
//  \brief Header file for the LAPACK solver test
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

#ifndef _BLAZETEST_MATHTEST_LAPACK_SOLVERTEST_H_
#define _BLAZETEST_MATHTEST_LAPACK_SOLVERTEST_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <sstream>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <blaze/math/HermitianMatrix.h>
#include <blaze/math/LAPACK.h>
#include <blaze/math/LowerMatrix.h>
#include <blaze/math/shims/Conjugate.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/StaticMatrix.h>
#include <blaze/math/StaticVector.h>
#include <blaze/math/SymmetricMatrix.h>
#include <blaze/math/UniLowerMatrix.h>
#include <blaze/math/UniUpperMatrix.h>
#include <blaze/math/UpperMatrix.h>
#include <blazetest/system/LAPACK.h>


namespace blazetest {

namespace mathtest {

namespace lapack {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class for all tests of the LAPACK functionality.
//
// This class represents a test suite for LAPACK functionality wrapped by the Blaze library.
*/
class SolverTest
{
 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit SolverTest();
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
   template< typename Type > void testGesv();
   template< typename Type > void testSysv();
   template< typename Type > void testHesv();
   template< typename Type > void testPosv();
   template< typename Type > void testTrsv();
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
/*!\brief Test of the general linear system solver functions (gesv).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the general linear system solver functions for various
// data types. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void SolverTest::testGesv()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major general LSE (single right-hand side)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A, LU;

      do {
         randomize( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      blaze::StaticVector<blaze::blas_int_t,3UL,blaze::columnVector> ipiv;

      LU = A;
      x = b;

      blaze::gesv( LU, x, ipiv.data() );

      if( ( trans( A ) * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   A * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major general LSE (multiple right-hand sides)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A, LU;

      do {
         randomize( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> B, X;
      randomize( B );

      blaze::StaticVector<blaze::blas_int_t,3UL,blaze::columnVector> ipiv;

      LU = A;
      X = B;

      blaze::gesv( LU, X, ipiv.data() );

      if( ( trans( A ) * trans( X ) ) != trans( B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (X):\n" << X << "\n"
             << "   Right-hand side (B):\n" << B << "\n"
             << "   A * X:\n" << ( trans( A ) * trans( X ) ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major general LSE (single right-hand side)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> A, LU;

      do {
         randomize( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      blaze::StaticVector<blaze::blas_int_t,3UL,blaze::columnVector> ipiv;

      LU = A;
      x = b;

      blaze::gesv( LU, x, ipiv.data() );

      if( ( A * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   A * x:\n" << ( A * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major general LSE (multiple right-hand sides)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> A, LU;

      do {
         randomize( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B, X;
      randomize( B );

      blaze::StaticVector<blaze::blas_int_t,3UL,blaze::columnVector> ipiv;

      LU = A;
      X = B;

      blaze::gesv( LU, X, ipiv.data() );

      if( ( A * X ) != B ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (X):\n" << X << "\n"
             << "   Right-hand side (B):\n" << B << "\n"
             << "   A * X:\n" << ( A * X ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the symmetric indefinite linear system solver functions (sysv).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the symmetric indefinite linear system solver functions for
// various data types. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void SolverTest::testSysv()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major symmetric indefinite LSE (single right-hand side, lower part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A, LU;

      do {
         randomize( A );
         A *= trans( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      blaze::StaticVector<blaze::blas_int_t,3UL,blaze::columnVector> ipiv;

      LU = A;
      x = b;

      blaze::sysv( LU, x, 'L', ipiv.data() );

      if( ( trans( A ) * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   A * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major symmetric indefinite LSE (single right-hand side, upper part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A, LU;

      do {
         randomize( A );
         A *= trans( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      blaze::StaticVector<blaze::blas_int_t,3UL,blaze::columnVector> ipiv;

      LU = A;
      x = b;

      blaze::sysv( LU, x, 'U', ipiv.data() );

      if( ( trans( A ) * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   A * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major symmetric indefinite LSE (multiple right-hand sides, lower part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A, LU;

      do {
         randomize( A );
         A *= trans( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> B, X;
      randomize( B );

      blaze::StaticVector<blaze::blas_int_t,3UL,blaze::columnVector> ipiv;

      LU = A;
      X = B;

      blaze::sysv( LU, X, 'L', ipiv.data() );

      if( ( trans( A ) * trans( X ) ) != trans( B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (X):\n" << X << "\n"
             << "   Right-hand side (B):\n" << B << "\n"
             << "   A * X:\n" << ( trans( A ) * trans( X ) ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major symmetric indefinite LSE (multiple right-hand sides, upper part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A, LU;

      do {
         randomize( A );
         A *= trans( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> B, X;
      randomize( B );

      blaze::StaticVector<blaze::blas_int_t,3UL,blaze::columnVector> ipiv;

      LU = A;
      X = B;

      blaze::sysv( LU, X, 'U', ipiv.data() );

      if( ( trans( A ) * trans( X ) ) != trans( B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (X):\n" << X << "\n"
             << "   Right-hand side (B):\n" << B << "\n"
             << "   A * X:\n" << ( trans( A ) * trans( X ) ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major symmetric indefinite LSE (single right-hand side, lower part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> A, LU;

      do {
         randomize( A );
         A *= trans( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      blaze::StaticVector<blaze::blas_int_t,3UL,blaze::columnVector> ipiv;

      LU = A;
      x = b;

      blaze::sysv( LU, x, 'L', ipiv.data() );

      if( ( A * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   A * x:\n" << ( A * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major symmetric indefinite LSE (single right-hand side, upper part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> A, LU;

      do {
         randomize( A );
         A *= trans( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      blaze::StaticVector<blaze::blas_int_t,3UL,blaze::columnVector> ipiv;

      LU = A;
      x = b;

      blaze::sysv( LU, x, 'U', ipiv.data() );

      if( ( A * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   A * x:\n" << ( A * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major symmetric indefinite LSE (multiple right-hand sides, lower part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> A, LU;

      do {
         randomize( A );
         A *= trans( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B, X;
      randomize( B );

      blaze::StaticVector<blaze::blas_int_t,3UL,blaze::columnVector> ipiv;

      LU = A;
      X = B;

      blaze::sysv( LU, X, 'L', ipiv.data() );

      if( ( A * X ) != B ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (X):\n" << X << "\n"
             << "   Right-hand side (B):\n" << B << "\n"
             << "   A * X:\n" << ( A * X ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major symmetric indefinite LSE (multiple right-hand sides, upper part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> A, LU;

      do {
         randomize( A );
         A *= trans( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B, X;
      randomize( B );

      blaze::StaticVector<blaze::blas_int_t,3UL,blaze::columnVector> ipiv;

      LU = A;
      X = B;

      blaze::sysv( LU, X, 'U', ipiv.data() );

      if( ( A * X ) != B ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (X):\n" << X << "\n"
             << "   Right-hand side (B):\n" << B << "\n"
             << "   A * X:\n" << ( A * X ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Hermitian indefinite linear system solver functions (hesv).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the Hermitian indefinite linear system solver functions for
// various data types. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void SolverTest::testHesv()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major Hermitian indefinite LSE (single right-hand side, lower part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A, LU;

      do {
         randomize( A );
         A *= ctrans( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      blaze::StaticVector<blaze::blas_int_t,3UL,blaze::columnVector> ipiv;

      LU = A;
      x = b;

      blaze::hesv( LU, x, 'L', ipiv.data() );

      if( ( trans( A ) * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   A * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major Hermitian indefinite LSE (single right-hand side, upper part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A, LU;

      do {
         randomize( A );
         A *= ctrans( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      blaze::StaticVector<blaze::blas_int_t,3UL,blaze::columnVector> ipiv;

      LU = A;
      x = b;

      blaze::hesv( LU, x, 'U', ipiv.data() );

      if( ( trans( A ) * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   A * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major Hermitian indefinite LSE (multiple right-hand sides, lower part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A, LU;

      do {
         randomize( A );
         A *= ctrans( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> B, X;
      randomize( B );

      blaze::StaticVector<blaze::blas_int_t,3UL,blaze::columnVector> ipiv;

      LU = A;
      X = B;

      blaze::hesv( LU, X, 'L', ipiv.data() );

      if( ( trans( A ) * trans( X ) ) != trans( B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (X):\n" << X << "\n"
             << "   Right-hand side (B):\n" << B << "\n"
             << "   A * X:\n" << ( trans( A ) * trans( X ) ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major Hermitian indefinite LSE (multiple right-hand sides, upper part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A, LU;

      do {
         randomize( A );
         A *= ctrans( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> B, X;
      randomize( B );

      blaze::StaticVector<blaze::blas_int_t,3UL,blaze::columnVector> ipiv;

      LU = A;
      X = B;

      blaze::hesv( LU, X, 'U', ipiv.data() );

      if( ( trans( A ) * trans( X ) ) != trans( B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (X):\n" << X << "\n"
             << "   Right-hand side (B):\n" << B << "\n"
             << "   A * X:\n" << ( trans( A ) * trans( X ) ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major Hermitian indefinite LSE (single right-hand side, lower part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> A, LU;

      do {
         randomize( A );
         A *= ctrans( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      blaze::StaticVector<blaze::blas_int_t,3UL,blaze::columnVector> ipiv;

      LU = A;
      x = b;

      blaze::hesv( LU, x, 'L', ipiv.data() );

      if( ( A * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   A * x:\n" << ( A * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major Hermitian indefinite LSE (single right-hand side, upper part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> A, LU;

      do {
         randomize( A );
         A *= ctrans( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      blaze::StaticVector<blaze::blas_int_t,3UL,blaze::columnVector> ipiv;

      LU = A;
      x = b;

      blaze::hesv( LU, x, 'U', ipiv.data() );

      if( ( A * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   A * x:\n" << ( A * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major Hermitian indefinite LSE (multiple right-hand sides, lower part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> A, LU;

      do {
         randomize( A );
         A *= ctrans( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B, X;
      randomize( B );

      blaze::StaticVector<blaze::blas_int_t,3UL,blaze::columnVector> ipiv;

      LU = A;
      X = B;

      blaze::hesv( LU, X, 'L', ipiv.data() );

      if( ( A * X ) != B ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (X):\n" << X << "\n"
             << "   Right-hand side (B):\n" << B << "\n"
             << "   A * X:\n" << ( A * X ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major Hermitian indefinite LSE (multiple right-hand sides, upper part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> A, LU;

      do {
         randomize( A );
         A *= ctrans( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B, X;
      randomize( B );

      blaze::StaticVector<blaze::blas_int_t,3UL,blaze::columnVector> ipiv;

      LU = A;
      X = B;

      blaze::hesv( LU, X, 'U', ipiv.data() );

      if( ( A * X ) != B ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (X):\n" << X << "\n"
             << "   Right-hand side (B):\n" << B << "\n"
             << "   A * X:\n" << ( A * X ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the positive definite linear system solver functions (posv).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the positive definite linear system solver functions for
// various data types. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void SolverTest::testPosv()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major positive definite LSE (single right-hand side, lower part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A, LU;

      do {
         randomize( A );
         A *= ctrans( A );
         A(0,0) += Type(3);
         A(1,1) += Type(3);
         A(2,2) += Type(3);
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      LU = A;
      x = b;

      blaze::posv( LU, x, 'L' );

      if( ( trans( A ) * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   A * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major positive definite LSE (single right-hand side, upper part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A, LU;

      do {
         randomize( A );
         A *= ctrans( A );
         A(0,0) += Type(3);
         A(1,1) += Type(3);
         A(2,2) += Type(3);
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      LU = A;
      x = b;

      blaze::posv( LU, x, 'U' );

      if( ( trans( A ) * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   A * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major positive definite LSE (multiple right-hand sides, lower part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A, LU;

      do {
         randomize( A );
         A *= ctrans( A );
         A(0,0) += Type(3);
         A(1,1) += Type(3);
         A(2,2) += Type(3);
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> B, X;
      randomize( B );

      LU = A;
      X = B;

      blaze::posv( LU, X, 'L' );

      if( ( trans( A ) * trans( X ) ) != trans( B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (X):\n" << X << "\n"
             << "   Right-hand side (B):\n" << B << "\n"
             << "   A * X:\n" << ( trans( A ) * trans( X ) ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major positive definite LSE (multiple right-hand sides, upper part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A, LU;

      do {
         randomize( A );
         A *= ctrans( A );
         A(0,0) += Type(3);
         A(1,1) += Type(3);
         A(2,2) += Type(3);
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> B, X;
      randomize( B );

      LU = A;
      X = B;

      blaze::posv( LU, X, 'U' );

      if( ( trans( A ) * trans( X ) ) != trans( B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (X):\n" << X << "\n"
             << "   Right-hand side (B):\n" << B << "\n"
             << "   A * X:\n" << ( trans( A ) * trans( X ) ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major positive definite LSE (single right-hand side, lower part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> A, LU;

      do {
         randomize( A );
         A *= ctrans( A );
         A(0,0) += Type(3);
         A(1,1) += Type(3);
         A(2,2) += Type(3);
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      LU = A;
      x = b;

      blaze::posv( LU, x, 'L' );

      if( ( A * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   A * x:\n" << ( A * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major positive definite LSE (single right-hand side, upper part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> A, LU;

      do {
         randomize( A );
         A *= ctrans( A );
         A(0,0) += Type(3);
         A(1,1) += Type(3);
         A(2,2) += Type(3);
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      LU = A;
      x = b;

      blaze::posv( LU, x, 'U' );

      if( ( A * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   A * x:\n" << ( A * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major positive definite LSE (multiple right-hand sides, lower part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> A, LU;

      do {
         randomize( A );
         A *= ctrans( A );
         A(0,0) += Type(3);
         A(1,1) += Type(3);
         A(2,2) += Type(3);
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B, X;
      randomize( B );

      LU = A;
      X = B;

      blaze::posv( LU, X, 'L' );

      if( ( A * X ) != B ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (X):\n" << X << "\n"
             << "   Right-hand side (B):\n" << B << "\n"
             << "   A * X:\n" << ( A * X ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major positive definite LSE (multiple right-hand sides, upper part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> A, LU;

      do {
         randomize( A );
         A *= ctrans( A );
         A(0,0) += Type(3);
         A(1,1) += Type(3);
         A(2,2) += Type(3);
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B, X;
      randomize( B );

      LU = A;
      X = B;

      blaze::posv( LU, X, 'U' );

      if( ( A * X ) != B ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (X):\n" << X << "\n"
             << "   Right-hand side (B):\n" << B << "\n"
             << "   A * X:\n" << ( A * X ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the triangular linear system solver functions (trsv).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the triangular linear system solver functions for various
// data types. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void SolverTest::testTrsv()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major triangular LSE (single right-hand side, lower part, not transposed)";

      blaze::LowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> LU( A );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trsv( LU, x, 'L', 'N', 'N' );

      if( ( trans( A ) * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   A * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major triangular LSE (single right-hand side, lower part, transposed)";

      blaze::LowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> LU( trans( A ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trsv( LU, x, 'U', 'T', 'N' );

      if( ( trans( A ) * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   A * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major triangular LSE (single right-hand side, lower part, conjugate transposed)";

      blaze::LowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> LU( ctrans( A ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trsv( LU, x, 'U', 'C', 'N' );

      if( ( trans( A ) * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   A * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major unitriangular LSE (single right-hand side, lower part, not transposed)";

      blaze::UniLowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> LU( A );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trsv( LU, x, 'L', 'N', 'U' );

      if( ( trans( A ) * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   A * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major unitriangular LSE (single right-hand side, lower part, transposed)";

      blaze::UniLowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> LU( trans( A ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trsv( LU, x, 'U', 'T', 'U' );

      if( ( trans( A ) * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   A * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major unitriangular LSE (single right-hand side, lower part, conjugate transposed)";

      blaze::UniLowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> LU( ctrans( A ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trsv( LU, x, 'U', 'C', 'U' );

      if( ( trans( A ) * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   A * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major triangular LSE (single right-hand side, upper part, not transposed)";

      blaze::UpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> LU( A );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trsv( LU, x, 'U', 'N', 'N' );

      if( ( trans( A ) * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   A * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major triangular LSE (single right-hand side, upper part, transposed)";

      blaze::UpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> LU( trans( A ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trsv( LU, x, 'L', 'T', 'N' );

      if( ( trans( A ) * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   A * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major triangular LSE (single right-hand side, upper part, conjugate transposed)";

      blaze::UpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> LU( ctrans( A ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trsv( LU, x, 'L', 'C', 'N' );

      if( ( trans( A ) * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   A * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major unitriangular LSE (single right-hand side, upper part, not transposed)";

      blaze::UniUpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> LU( A );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trsv( LU, x, 'U', 'N', 'U' );

      if( ( trans( A ) * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   A * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major unitriangular LSE (single right-hand side, upper part, transposed)";

      blaze::UniUpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> LU( trans( A ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trsv( LU, x, 'L', 'T', 'U' );

      if( ( trans( A ) * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   A * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major unitriangular LSE (single right-hand side, upper part, conjugate transposed)";

      blaze::UniUpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> LU( ctrans( A ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trsv( LU, x, 'L', 'C', 'U' );

      if( ( trans( A ) * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   A * x:\n" << ( trans( A ) * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major triangular LSE (single right-hand side, lower part, not transposed)";

      blaze::LowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> LU( A );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trsv( LU, x, 'L', 'N', 'N' );

      if( ( A * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   A * x:\n" << ( A * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major triangular LSE (single right-hand side, lower part, transposed)";

      blaze::LowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> LU( trans( A ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trsv( LU, x, 'U', 'T', 'N' );

      if( ( A * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   A * x:\n" << ( A * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major triangular LSE (single right-hand side, lower part, conjugate transposed)";

      blaze::LowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> LU( ctrans( A ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trsv( LU, x, 'U', 'C', 'N' );

      if( ( A * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   A * x:\n" << ( A * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major unitriangular LSE (single right-hand side, lower part, not transposed)";

      blaze::UniLowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> LU( A );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trsv( LU, x, 'L', 'N', 'U' );

      if( ( A * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   A * x:\n" << ( A * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major unitriangular LSE (single right-hand side, lower part, transposed)";

      blaze::UniLowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> LU( trans( A ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trsv( LU, x, 'U', 'T', 'U' );

      if( ( A * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   A * x:\n" << ( A * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major unitriangular LSE (single right-hand side, lower part, conjugate transposed)";

      blaze::UniLowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> LU( ctrans( A ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trsv( LU, x, 'U', 'C', 'U' );

      if( ( A * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   A * x:\n" << ( A * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major triangular LSE (single right-hand side, upper part, not transposed)";

      blaze::UpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> LU( A );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trsv( LU, x, 'U', 'N', 'N' );

      if( ( A * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   A * x:\n" << ( A * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major triangular LSE (single right-hand side, upper part, transposed)";

      blaze::UpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> LU( trans( A ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trsv( LU, x, 'L', 'T', 'N' );

      if( ( A * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   A * x:\n" << ( A * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major triangular LSE (single right-hand side, upper part, conjugate transposed)";

      blaze::UpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> LU( ctrans( A ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trsv( LU, x, 'L', 'C', 'N' );

      if( ( A * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   A * x:\n" << ( A * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major unitriangular LSE (single right-hand side, upper part, not transposed)";

      blaze::UniUpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> LU( A );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trsv( LU, x, 'U', 'N', 'U' );

      if( ( A * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   A * x:\n" << ( A * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major unitriangular LSE (single right-hand side, upper part, transposed)";

      blaze::UniUpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> LU( trans( A ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trsv( LU, x, 'L', 'T', 'U' );

      if( ( A * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   A * x:\n" << ( A * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major unitriangular LSE (single right-hand side, upper part, conjugate transposed)";

      blaze::UniUpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> LU( ctrans( A ) );

      blaze::StaticVector<Type,3UL,blaze::columnVector> b, x;
      randomize( b );

      x = b;

      blaze::trsv( LU, x, 'L', 'C', 'U' );

      if( ( A * x ) != b ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Solving the LSE failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   System matrix (A):\n" << A << "\n"
             << "   Result (x):\n" << x << "\n"
             << "   Right-hand side (b):\n" << b << "\n"
             << "   A * x:\n" << ( A * x ) << "\n";
         throw std::runtime_error( oss.str() );
      }
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
/*!\brief Testing the LAPACK functionality.
//
// \return void
*/
void runTest()
{
   SolverTest();
}
//*************************************************************************************************




//=================================================================================================
//
//  MACRO DEFINITIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Macro for the execution of the LAPACK solver test.
*/
#define RUN_LAPACK_SOLVER_TEST \
   blazetest::mathtest::lapack::runTest()
/*! \endcond */
//*************************************************************************************************

} // namespace lapack

} // namespace mathtest

} // namespace blazetest

#endif
