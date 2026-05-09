//=================================================================================================
/*!
//  \file blazetest/mathtest/lapack/InversionTest.h
//  \brief Header file for the LAPACK inversion test
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

#ifndef _BLAZETEST_MATHTEST_LAPACK_INVERSIONTEST_H_
#define _BLAZETEST_MATHTEST_LAPACK_INVERSIONTEST_H_


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
class InversionTest
{
 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit InversionTest();
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
   template< typename Type > void testGetri();
   template< typename Type > void testSytri();
   template< typename Type > void testHetri();
   template< typename Type > void testPotri();
   template< typename Type > void testTrtri();
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
/*!\brief Test of the LU-based matrix inversion functions (getri).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the LU-based matrix inversion functions for various data
// types. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void InversionTest::testGetri()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major LU-based matrix inversion";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A;

      do {
         randomize( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> Ainv( A );
      blaze::StaticVector<blaze::blas_int_t,3UL,blaze::rowVector> ipiv;

      blaze::getrf( Ainv, ipiv.data() );
      blaze::getri( Ainv, ipiv.data() );

      if( !blaze::isIdentity( Ainv * A ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: LU-based matrix inversion failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Result:\n" << Ainv << "\n"
             << "   Ainv * A = " << ( Ainv * A ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major LU-based matrix inversion";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> A;

      do {
         randomize( A );
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> Ainv( A );
      blaze::StaticVector<blaze::blas_int_t,3UL,blaze::rowVector> ipiv;

      blaze::getrf( Ainv, ipiv.data() );
      blaze::getri( Ainv, ipiv.data() );

      if( !blaze::isIdentity( Ainv * A ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: LU-based matrix inversion failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Result:\n" << Ainv << "\n"
             << "   Ainv * A = " << ( Ainv * A ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Bunch-Kaufman-based matrix inversion functions for symmetric matrices (sytri).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the Bunch-Kaufman-based matrix inversion functions for
// symmetric indefinite matrices for various data types. In case an error is detected, a
// \a std::runtime_error exception is thrown.
*/
template< typename Type >
void InversionTest::testSytri()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major symmetric matrix inversion (lower part)";

      blaze::SymmetricMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> B( A );
      blaze::StaticVector<blaze::blas_int_t,3UL,blaze::rowVector> ipiv;

      blaze::sytrf( B, 'L', ipiv.data() );
      blaze::sytri( B, 'L', ipiv.data() );

      B(0,1) = B(1,0);
      B(0,2) = B(2,0);
      B(1,2) = B(2,1);

      if( !blaze::isIdentity( A * B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Symmetric matrix inversion failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Result:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major symmetric matrix inversion (upper part)";

      blaze::SymmetricMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> B( A );
      blaze::StaticVector<blaze::blas_int_t,3UL,blaze::rowVector> ipiv;

      blaze::sytrf( B, 'U', ipiv.data() );
      blaze::sytri( B, 'U', ipiv.data() );

      B(1,0) = B(0,1);
      B(2,0) = B(0,2);
      B(2,1) = B(1,2);

      if( !blaze::isIdentity( A * B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Symmetric matrix inversion failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Result:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major symmetric matrix inversion (lower part)";

      blaze::SymmetricMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B( A );
      blaze::StaticVector<blaze::blas_int_t,3UL,blaze::rowVector> ipiv;

      blaze::sytrf( B, 'L', ipiv.data() );
      blaze::sytri( B, 'L', ipiv.data() );

      B(0,1) = B(1,0);
      B(0,2) = B(2,0);
      B(1,2) = B(2,1);

      if( !blaze::isIdentity( A * B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Symmetric matrix inversion failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Result:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major symmetric matrix inversion (upper part)";

      blaze::SymmetricMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B( A );
      blaze::StaticVector<blaze::blas_int_t,3UL,blaze::rowVector> ipiv;

      blaze::sytrf( B, 'U', ipiv.data() );
      blaze::sytri( B, 'U', ipiv.data() );

      B(1,0) = B(0,1);
      B(2,0) = B(0,2);
      B(2,1) = B(1,2);

      if( !blaze::isIdentity( A * B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Symmetric matrix inversion failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Result:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Bunch-Kaufman-based matrix inversion functions for Hermitian matrices (hetri).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the Bunch-Kaufman-based matrix inversion functions for
// Hermitian indefinite matrices for various data types. In case an error is detected, a
// \a std::runtime_error exception is thrown.
*/
template< typename Type >
void InversionTest::testHetri()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   using blaze::conj;


   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major Hermitian matrix inversion (lower part)";

      blaze::HermitianMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> B( A );
      blaze::StaticVector<blaze::blas_int_t,3UL,blaze::rowVector> ipiv;

      blaze::hetrf( B, 'L', ipiv.data() );
      blaze::hetri( B, 'L', ipiv.data() );

      B(0,1) = conj( B(1,0) );
      B(0,2) = conj( B(2,0) );
      B(1,2) = conj( B(2,1) );

      if( !blaze::isIdentity( A * B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Hermitian matrix inversion failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Result:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major Hermitian matrix inversion (upper part)";

      blaze::HermitianMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> B( A );
      blaze::StaticVector<blaze::blas_int_t,3UL,blaze::rowVector> ipiv;

      blaze::hetrf( B, 'U', ipiv.data() );
      blaze::hetri( B, 'U', ipiv.data() );

      B(1,0) = conj( B(0,1) );
      B(2,0) = conj( B(0,2) );
      B(2,1) = conj( B(1,2) );

      if( !blaze::isIdentity( A * B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Hermitian matrix inversion failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Result:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major Hermitian matrix inversion (lower part)";

      blaze::HermitianMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B( A );
      blaze::StaticVector<blaze::blas_int_t,3UL,blaze::rowVector> ipiv;

      blaze::hetrf( B, 'L', ipiv.data() );
      blaze::hetri( B, 'L', ipiv.data() );

      B(0,1) = conj( B(1,0) );
      B(0,2) = conj( B(2,0) );
      B(1,2) = conj( B(2,1) );

      if( !blaze::isIdentity( A * B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Hermitian matrix inversion failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Result:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major Hermitian matrix inversion (upper part)";

      blaze::HermitianMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B( A );
      blaze::StaticVector<blaze::blas_int_t,3UL,blaze::rowVector> ipiv;

      blaze::hetrf( B, 'U', ipiv.data() );
      blaze::hetri( B, 'U', ipiv.data() );

      B(1,0) = conj( B(0,1) );
      B(2,0) = conj( B(0,2) );
      B(2,1) = conj( B(1,2) );

      if( !blaze::isIdentity( A * B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Hermitian matrix inversion failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Result:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Cholesky-based matrix inversion functions (potri).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the Cholesky-based matrix inversion functions for various
// data types. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void InversionTest::testPotri()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   using blaze::conj;


   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major Cholesky-based matrix inversion (lower part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A;

      do {
         randomize( A );
         A *= ctrans( A );
         A(0,0) += Type(3);
         A(1,1) += Type(3);
         A(2,2) += Type(3);
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> B( A );

      blaze::potrf( A, 'L' );
      blaze::potri( A, 'L' );

      A(0,1) = conj( A(1,0) );
      A(0,2) = conj( A(2,0) );
      A(1,2) = conj( A(2,1) );

      if( !blaze::isIdentity( A * B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cholesky-based matrix inversion failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Result:\n" << A << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major Cholesky-based matrix inversion (upper part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A;

      do {
         randomize( A );
         A *= ctrans( A );
         A(0,0) += Type(3);
         A(1,1) += Type(3);
         A(2,2) += Type(3);
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> B( A );

      blaze::potrf( A, 'U' );
      blaze::potri( A, 'U' );

      A(1,0) = conj( A(0,1) );
      A(2,0) = conj( A(0,2) );
      A(2,1) = conj( A(1,2) );

      if( !blaze::isIdentity( A * B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cholesky-based matrix inversion failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Result:\n" << A << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major Cholesky-based matrix inversion (lower part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> A;

      do {
         randomize( A );
         A *= ctrans( A );
         A(0,0) += Type(3);
         A(1,1) += Type(3);
         A(2,2) += Type(3);
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B( A );

      blaze::potrf( A, 'L' );
      blaze::potri( A, 'L' );

      A(0,1) = conj( A(1,0) );
      A(0,2) = conj( A(2,0) );
      A(1,2) = conj( A(2,1) );

      if( !blaze::isIdentity( A * B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cholesky-based matrix inversion failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Result:\n" << A << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major Cholesky-based matrix inversion (upper part)";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> A;

      do {
         randomize( A );
         A *= ctrans( A );
         A(0,0) += Type(3);
         A(1,1) += Type(3);
         A(2,2) += Type(3);
      }
      while( blaze::isDefault( det( A ) ) );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B( A );

      blaze::potrf( A, 'U' );
      blaze::potri( A, 'U' );

      A(1,0) = conj( A(0,1) );
      A(2,0) = conj( A(0,2) );
      A(2,1) = conj( A(1,2) );

      if( !blaze::isIdentity( A * B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cholesky-based matrix inversion failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Result:\n" << A << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the triangular matrix inversion functions (trtri).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the triangular matrix inversion functions for various data
// types. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void InversionTest::testTrtri()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major lower triangular matrix inversion";

      blaze::LowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> B( A );

      blaze::trtri( B, 'L', 'N' );

      if( !blaze::isIdentity( A * B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Lower triangular matrix inversion failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Result:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major lower unitriangular matrix inversion";

      blaze::UniLowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> B( A );

      blaze::trtri( B, 'L', 'U' );

      if( !blaze::isIdentity( A * B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Lower unitriangular matrix inversion failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Result:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major upper triangular matrix inversion";

      blaze::UpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> B( A );

      blaze::trtri( B, 'U', 'N' );

      if( !blaze::isIdentity( A * B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Upper triangular matrix inversion failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Result:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major upper unitriangular matrix inversion";

      blaze::UniUpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> B( A );

      blaze::trtri( B, 'U', 'U' );

      if( !blaze::isIdentity( A * B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Upper unitriangular matrix inversion failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Result:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major lower triangular matrix inversion";

      blaze::LowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B( A );

      blaze::trtri( B, 'L', 'N' );

      if( !blaze::isIdentity( A * B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Lower triangular matrix inversion failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Result:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major lower unitriangular matrix inversion";

      blaze::UniLowerMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B( A );

      blaze::trtri( B, 'L', 'U' );

      if( !blaze::isIdentity( A * B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Lower unitriangular matrix inversion failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Result:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major upper triangular matrix inversion";

      blaze::UpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B( A );

      blaze::trtri( B, 'U', 'N' );

      if( !blaze::isIdentity( A * B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Upper triangular matrix inversion failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Result:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major upper unitriangular matrix inversion";

      blaze::UniUpperMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> > A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B( A );

      blaze::trtri( B, 'U', 'U' );

      if( !blaze::isIdentity( A * B ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Upper unitriangular matrix inversion failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Result:\n" << B << "\n";
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
   InversionTest();
}
//*************************************************************************************************




//=================================================================================================
//
//  MACRO DEFINITIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Macro for the execution of the LAPACK inversion test.
*/
#define RUN_LAPACK_INVERSION_TEST \
   blazetest::mathtest::lapack::runTest()
/*! \endcond */
//*************************************************************************************************

} // namespace lapack

} // namespace mathtest

} // namespace blazetest

#endif
