//=================================================================================================
/*!
//  \file blazetest/mathtest/lapack/DecompositionTest.h
//  \brief Header file for the LAPACK decomposition test
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

#ifndef _BLAZETEST_MATHTEST_LAPACK_DECOMPOSITIONTEST_H_
#define _BLAZETEST_MATHTEST_LAPACK_DECOMPOSITIONTEST_H_


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
class DecompositionTest
{
 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit DecompositionTest();
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
   template< typename Type > void testGetrf();
   template< typename Type > void testSytrf();
   template< typename Type > void testHetrf();
   template< typename Type > void testPotrf();

   template< typename Type > void testGeqrf();
   template< typename Type > void testOrgqr();
   template< typename Type > void testOrg2r();
   template< typename Type > void testUngqr();
   template< typename Type > void testUng2r();
   template< typename Type > void testOrmqr();
   template< typename Type > void testUnmqr();

   template< typename Type > void testGerqf();
   template< typename Type > void testOrgrq();
   template< typename Type > void testOrgr2();
   template< typename Type > void testUngrq();
   template< typename Type > void testUngr2();
   template< typename Type > void testOrmrq();
   template< typename Type > void testUnmrq();

   template< typename Type > void testGeqlf();
   template< typename Type > void testOrgql();
   template< typename Type > void testOrg2l();
   template< typename Type > void testUngql();
   template< typename Type > void testUng2l();
   template< typename Type > void testOrmql();
   template< typename Type > void testUnmql();

   template< typename Type > void testGelqf();
   template< typename Type > void testOrglq();
   template< typename Type > void testOrgl2();
   template< typename Type > void testUnglq();
   template< typename Type > void testUngl2();
   template< typename Type > void testOrmlq();
   template< typename Type > void testUnmlq();
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
/*!\brief Test of the LU decomposition functions (getrf).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the LU decomposition functions for various data types. In
// case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void DecompositionTest::testGetrf()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   test_ = "LU decomposition";

   {
      blaze::StaticMatrix<Type,2UL,5UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,5UL,2UL,blaze::columnMajor> B( trans( A ) );

      blaze::StaticVector<blaze::blas_int_t,2UL,blaze::columnVector> ipivA;
      blaze::StaticVector<blaze::blas_int_t,2UL,blaze::columnVector> ipivB;

      blaze::getrf( A, ipivA.data() );
      blaze::getrf( B, ipivB.data() );

      if( A != trans( B ) || ipivA != ipivB ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: LU decomposition failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A << "\n"
             << "   Row-major pivot elements:\n" << ipivA << "\n"
             << "   Column-major decomposition:\n" << B << "\n"
             << "   Column-major pivot elements:\n" << ipivB << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      blaze::StaticMatrix<Type,5UL,2UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,2UL,5UL,blaze::columnMajor> B( trans( A ) );

      blaze::StaticVector<blaze::blas_int_t,2UL,blaze::columnVector> ipivA;
      blaze::StaticVector<blaze::blas_int_t,2UL,blaze::columnVector> ipivB;

      blaze::getrf( A, ipivA.data() );
      blaze::getrf( B, ipivB.data() );

      if( A != trans( B ) || ipivA != ipivB ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: LU decomposition failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A << "\n"
             << "   Row-major pivot elements:\n" << ipivA << "\n"
             << "   Column-major decomposition:\n" << B << "\n"
             << "   Column-major pivot elements:\n" << ipivB << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Bunch-Kaufman decomposition functions for symmetric matrices (sytrf).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the Bunch-Kaufman decomposition functions for symmetric
// indefinite matrices for various data types. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename Type >
void DecompositionTest::testSytrf()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   test_ = "Symmetric matrix decomposition";

   {
      blaze::SymmetricMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > S;
      randomize( S );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor>    A( S );
      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B( S );

      blaze::StaticVector<blaze::blas_int_t,3UL,blaze::rowVector> ipivA;
      blaze::StaticVector<blaze::blas_int_t,3UL,blaze::rowVector> ipivB;

      blaze::sytrf( A, 'L', ipivA.data() );
      blaze::sytrf( B, 'U', ipivB.data() );

      if( A != trans( B ) || ipivA != ipivB ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Symmetric matrix decomposition failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A << "\n"
             << "   Row-major pivot elements:\n" << ipivA << "\n"
             << "   Column-major decomposition:\n" << B << "\n"
             << "   Column-major pivot elements:\n" << ipivB << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Bunch-Kaufman decomposition functions for Hermitian matrices (hetrf).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the Bunch-Kaufman decomposition functions for Hermitian
// indefinite matrices for various data types. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename Type >
void DecompositionTest::testHetrf()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   test_ = "Hermitian matrix decomposition";

   {
      blaze::HermitianMatrix< blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> > H;
      randomize( H );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor>    A( H );
      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B( H );

      blaze::StaticVector<blaze::blas_int_t,3UL,blaze::rowVector> ipivA;
      blaze::StaticVector<blaze::blas_int_t,3UL,blaze::rowVector> ipivB;

      blaze::hetrf( A, 'L', ipivA.data() );
      blaze::hetrf( B, 'U', ipivB.data() );

      if( A != ctrans( B ) || ipivA != ipivB ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Hermitian matrix decomposition failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A << "\n"
             << "   Row-major pivot elements:\n" << ipivA << "\n"
             << "   Column-major decomposition:\n" << B << "\n"
             << "   Column-major pivot elements:\n" << ipivB << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Cholesky decomposition functions (potrf).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the Cholesky decomposition functions for various data types.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void DecompositionTest::testPotrf()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   test_ = "Cholesky decomposition";

   {
      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A;

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
      blaze::potrf( B, 'L' );

      if( A != B ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cholesky decomposition failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A << "\n"
             << "   Column-major decomposition:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A;

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
      blaze::potrf( B, 'U' );

      if( A != B ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cholesky decomposition failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A << "\n"
             << "   Column-major decomposition:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the QR decomposition functions (geqrf).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the QR decomposition functions for various data types. In
// case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void DecompositionTest::testGeqrf()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   test_ = "QR decomposition";

   {
      blaze::StaticMatrix<Type,2UL,5UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,2UL,5UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,2UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,2UL,blaze::rowVector> tauB;

      blaze::geqrf( A, tauA.data() );
      blaze::geqrf( B, tauB.data() );

      if( A != B || tauA != conj( tauB ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: QR decomposition failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A << "\n"
             << "   Row-major projectors:\n" << tauA << "\n"
             << "   Column-major decomposition:\n" << B << "\n"
             << "   Column-major projectors:\n" << tauB << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      blaze::StaticMatrix<Type,5UL,2UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,5UL,2UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,2UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,2UL,blaze::rowVector> tauB;

      blaze::geqrf( A, tauA.data() );
      blaze::geqrf( B, tauB.data() );

      if( A != B || tauA != conj( tauB ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: QR decomposition failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A << "\n"
             << "   Row-major projectors:\n" << tauA << "\n"
             << "   Column-major decomposition:\n" << B << "\n"
             << "   Column-major projectors:\n" << tauB << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Q reconstruction from a QR decomposition (orgqr).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the Q reconstruction from a QR decomposition for various data
// types. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void DecompositionTest::testOrgqr()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   test_ = "Reconstruction of Q from a QR decomposition (orgqr)";

   {
      blaze::StaticMatrix<Type,2UL,5UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,2UL,5UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,2UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,2UL,blaze::rowVector> tauB;

      blaze::geqrf( A, tauA.data() );
      blaze::geqrf( B, tauB.data() );

      blaze::orgqr( A, tauA.data() );
      blaze::orgqr( B, tauB.data() );

      if( A != B ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q reconstruction failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major reconstruction:\n" << A << "\n"
             << "   Column-major reconstruction:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      blaze::StaticMatrix<Type,5UL,2UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,5UL,2UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,2UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,2UL,blaze::rowVector> tauB;

      blaze::geqrf( A, tauA.data() );
      blaze::geqrf( B, tauB.data() );

      blaze::orgqr( A, tauA.data() );
      blaze::orgqr( B, tauB.data() );

      if( A != B ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q reconstruction failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major reconstruction:\n" << A << "\n"
             << "   Column-major reconstruction:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Q reconstruction from a QR decomposition (org2r).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the Q reconstruction from a QR decomposition for various data
// types. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void DecompositionTest::testOrg2r()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   test_ = "Reconstruction of Q from a QR decomposition (org2r)";

   {
      blaze::StaticMatrix<Type,2UL,5UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,2UL,5UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,2UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,2UL,blaze::rowVector> tauB;

      blaze::geqrf( A, tauA.data() );
      blaze::geqrf( B, tauB.data() );

      blaze::org2r( A, tauA.data() );
      blaze::org2r( B, tauB.data() );

      if( A != B ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q reconstruction failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major reconstruction:\n" << A << "\n"
             << "   Column-major reconstruction:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      blaze::StaticMatrix<Type,5UL,2UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,5UL,2UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,2UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,2UL,blaze::rowVector> tauB;

      blaze::geqrf( A, tauA.data() );
      blaze::geqrf( B, tauB.data() );

      blaze::org2r( A, tauA.data() );
      blaze::org2r( B, tauB.data() );

      if( A != B ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q reconstruction failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major reconstruction:\n" << A << "\n"
             << "   Column-major reconstruction:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Q reconstruction from a QR decomposition (ungqr).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the Q reconstruction from a QR decomposition for various data
// types. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void DecompositionTest::testUngqr()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   test_ = "Reconstruction of Q from a QR decomposition (ungqr)";

   {
      blaze::StaticMatrix<Type,2UL,5UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,2UL,5UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,2UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,2UL,blaze::rowVector> tauB;

      blaze::geqrf( A, tauA.data() );
      blaze::geqrf( B, tauB.data() );

      blaze::ungqr( A, tauA.data() );
      blaze::ungqr( B, tauB.data() );

      if( A != B ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q reconstruction failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major reconstruction:\n" << A << "\n"
             << "   Column-major reconstruction:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      blaze::StaticMatrix<Type,5UL,2UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,5UL,2UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,2UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,2UL,blaze::rowVector> tauB;

      blaze::geqrf( A, tauA.data() );
      blaze::geqrf( B, tauB.data() );

      blaze::ungqr( A, tauA.data() );
      blaze::ungqr( B, tauB.data() );

      if( A != B ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q reconstruction failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major reconstruction:\n" << A << "\n"
             << "   Column-major reconstruction:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Q reconstruction from a QR decomposition (ung2r).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the Q reconstruction from a QR decomposition for various data
// types. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void DecompositionTest::testUng2r()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   test_ = "Reconstruction of Q from a QR decomposition (ung2r)";

   {
      blaze::StaticMatrix<Type,2UL,5UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,2UL,5UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,2UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,2UL,blaze::rowVector> tauB;

      blaze::geqrf( A, tauA.data() );
      blaze::geqrf( B, tauB.data() );

      blaze::ung2r( A, tauA.data() );
      blaze::ung2r( B, tauB.data() );

      if( A != B ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q reconstruction failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major reconstruction:\n" << A << "\n"
             << "   Column-major reconstruction:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      blaze::StaticMatrix<Type,5UL,2UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,5UL,2UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,2UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,2UL,blaze::rowVector> tauB;

      blaze::geqrf( A, tauA.data() );
      blaze::geqrf( B, tauB.data() );

      blaze::ung2r( A, tauA.data() );
      blaze::ung2r( B, tauB.data() );

      if( A != B ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q reconstruction failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major reconstruction:\n" << A << "\n"
             << "   Column-major reconstruction:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the multiplication of Q from a QR decomposition with a matrix (ormqr).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the multiplication of Q from a QR decomposition with a matrix
// for various data types. In case an error is detected, a \a std::runtime_error exception is
// thrown.
*/
template< typename Type >
void DecompositionTest::testOrmqr()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   {
      test_ = "Multiplication of Q from a QR decomposition with a matrix ('L', 'N')";

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,3UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,3UL,blaze::rowVector> tauB;

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> C1;
      randomize( C1 );

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor>    C2( C1 );
      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> C3( C1 );
      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> C4( C1 );

      blaze::geqrf( A, tauA.data() );
      blaze::geqrf( B, tauB.data() );

      blaze::ormqr( C1, A, 'L', 'N', tauA.data() );
      blaze::ormqr( C2, B, 'L', 'N', tauB.data() );
      blaze::ormqr( C3, A, 'L', 'N', tauA.data() );
      blaze::ormqr( C4, B, 'L', 'N', tauB.data() );

      if( C1 != C2 || C1 != C3 || C1 != C4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q multiplication failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major/row-major multiplication:\n" << C1 << "\n"
             << "   Row-major/column-major multiplication:\n" << C2 << "\n"
             << "   Column-major/row-major multiplication:\n" << C3 << "\n"
             << "   Column-major/column-major multiplication:\n" << C4 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Multiplication of Q from a QR decomposition with a matrix ('L', 'T')";

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,3UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,3UL,blaze::rowVector> tauB;

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> C1;
      randomize( C1 );

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor>    C2( C1 );
      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> C3( C1 );
      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> C4( C1 );

      blaze::geqrf( A, tauA.data() );
      blaze::geqrf( B, tauB.data() );

      blaze::ormqr( C1, A, 'L', 'T', tauA.data() );
      blaze::ormqr( C2, B, 'L', 'T', tauB.data() );
      blaze::ormqr( C3, A, 'L', 'T', tauA.data() );
      blaze::ormqr( C4, B, 'L', 'T', tauB.data() );

      if( C1 != C2 || C1 != C3 || C1 != C4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q multiplication failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major/row-major multiplication:\n" << C1 << "\n"
             << "   Row-major/column-major multiplication:\n" << C2 << "\n"
             << "   Column-major/row-major multiplication:\n" << C3 << "\n"
             << "   Column-major/column-major multiplication:\n" << C4 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Multiplication of Q from a QR decomposition with a matrix ('R', 'N')";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,3UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,3UL,blaze::rowVector> tauB;

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> C1;
      randomize( C1 );

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor>    C2( C1 );
      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> C3( C1 );
      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> C4( C1 );

      blaze::geqrf( A, tauA.data() );
      blaze::geqrf( B, tauB.data() );

      blaze::ormqr( C1, A, 'R', 'N', tauA.data() );
      blaze::ormqr( C2, B, 'R', 'N', tauB.data() );
      blaze::ormqr( C3, A, 'R', 'N', tauA.data() );
      blaze::ormqr( C4, B, 'R', 'N', tauB.data() );

      if( C1 != C2 || C1 != C3 || C1 != C4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q multiplication failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major/row-major multiplication:\n" << C1 << "\n"
             << "   Row-major/column-major multiplication:\n" << C2 << "\n"
             << "   Column-major/row-major multiplication:\n" << C3 << "\n"
             << "   Column-major/column-major multiplication:\n" << C4 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Multiplication of Q from a QR decomposition with a matrix ('R', 'T')";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,3UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,3UL,blaze::rowVector> tauB;

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> C1;
      randomize( C1 );

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor>    C2( C1 );
      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> C3( C1 );
      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> C4( C1 );

      blaze::geqrf( A, tauA.data() );
      blaze::geqrf( B, tauB.data() );

      blaze::ormqr( C1, A, 'R', 'T', tauA.data() );
      blaze::ormqr( C2, B, 'R', 'T', tauB.data() );
      blaze::ormqr( C3, A, 'R', 'T', tauA.data() );
      blaze::ormqr( C4, B, 'R', 'T', tauB.data() );

      if( C1 != C2 || C1 != C3 || C1 != C4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q multiplication failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major/row-major multiplication:\n" << C1 << "\n"
             << "   Row-major/column-major multiplication:\n" << C2 << "\n"
             << "   Column-major/row-major multiplication:\n" << C3 << "\n"
             << "   Column-major/column-major multiplication:\n" << C4 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the multiplication of Q from a QR decomposition with a matrix (unmqr).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the multiplication of Q from a QR decomposition with a matrix
// for various data types. In case an error is detected, a \a std::runtime_error exception is
// thrown.
*/
template< typename Type >
void DecompositionTest::testUnmqr()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   {
      test_ = "Multiplication of Q from a QR decomposition with a matrix ('L', 'N')";

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,3UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,3UL,blaze::rowVector> tauB;

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> C1;
      randomize( C1 );

      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> C2( C1 );

      blaze::geqrf( A, tauA.data() );
      blaze::geqrf( B, tauB.data() );

      blaze::unmqr( C1, A, 'L', 'N', tauA.data() );
      blaze::unmqr( C2, B, 'L', 'N', tauB.data() );

      if( C1 != C2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q multiplication failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major multiplication:\n" << C1 << "\n"
             << "   Column-major multiplication:\n" << C2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Multiplication of Q from a QR decomposition with a matrix ('L', 'C')";

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,3UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,3UL,blaze::rowVector> tauB;

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> C1;
      randomize( C1 );

      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> C2( C1 );

      blaze::geqrf( A, tauA.data() );
      blaze::geqrf( B, tauB.data() );

      blaze::unmqr( C1, A, 'L', 'C', tauA.data() );
      blaze::unmqr( C2, B, 'L', 'C', tauB.data() );

      if( C1 != C2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q multiplication failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major multiplication:\n" << C1 << "\n"
             << "   Column-major multiplication:\n" << C2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Multiplication of Q from a QR decomposition with a matrix ('R', 'N')";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,3UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,3UL,blaze::rowVector> tauB;

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> C1;
      randomize( C1 );

      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> C2( C1 );

      blaze::geqrf( A, tauA.data() );
      blaze::geqrf( B, tauB.data() );

      blaze::unmqr( C1, A, 'R', 'N', tauA.data() );
      blaze::unmqr( C2, B, 'R', 'N', tauB.data() );

      if( C1 != C2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q multiplication failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major multiplication:\n" << C1 << "\n"
             << "   Column-major multiplication:\n" << C2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Multiplication of Q from a QR decomposition with a matrix ('R', 'C')";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,3UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,3UL,blaze::rowVector> tauB;

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> C1;
      randomize( C1 );

      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> C2( C1 );

      blaze::geqrf( A, tauA.data() );
      blaze::geqrf( B, tauB.data() );

      blaze::unmqr( C1, A, 'R', 'C', tauA.data() );
      blaze::unmqr( C2, B, 'R', 'C', tauB.data() );

      if( C1 != C2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q multiplication failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major multiplication:\n" << C1 << "\n"
             << "   Column-major multiplication:\n" << C2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the RQ decomposition functions (gerqf).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the RQ decomposition functions for various data types. In
// case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void DecompositionTest::testGerqf()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   test_ = "RQ decomposition";

   {
      blaze::StaticMatrix<Type,2UL,5UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,2UL,5UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,2UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,2UL,blaze::rowVector> tauB;

      blaze::gerqf( A, tauA.data() );
      blaze::gerqf( B, tauB.data() );

      if( A != B || tauA != conj( tauB ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: RQ decomposition failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A << "\n"
             << "   Row-major projectors:\n" << tauA << "\n"
             << "   Column-major decomposition:\n" << B << "\n"
             << "   Column-major projectors:\n" << tauB << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      blaze::StaticMatrix<Type,5UL,2UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,5UL,2UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,2UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,2UL,blaze::rowVector> tauB;

      blaze::gerqf( A, tauA.data() );
      blaze::gerqf( B, tauB.data() );

      if( A != B || tauA != conj( tauB ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: RQ decomposition failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A << "\n"
             << "   Row-major projectors:\n" << tauA << "\n"
             << "   Column-major decomposition:\n" << B << "\n"
             << "   Column-major projectors:\n" << tauB << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Q reconstruction from a RQ decomposition (orgrq).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the Q reconstruction from a RQ decomposition for various data
// types. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void DecompositionTest::testOrgrq()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   test_ = "Reconstruction of Q from a RQ decomposition (orgrq)";

   {
      blaze::StaticMatrix<Type,2UL,5UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,2UL,5UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,2UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,2UL,blaze::rowVector> tauB;

      blaze::gerqf( A, tauA.data() );
      blaze::gerqf( B, tauB.data() );

      blaze::orgrq( A, tauA.data() );
      blaze::orgrq( B, tauB.data() );

      if( A != B ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q reconstruction failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major reconstruction:\n" << A << "\n"
             << "   Column-major reconstruction:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      blaze::StaticMatrix<Type,5UL,2UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,5UL,2UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,2UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,2UL,blaze::rowVector> tauB;

      blaze::gerqf( A, tauA.data() );
      blaze::gerqf( B, tauB.data() );

      blaze::orgrq( A, tauA.data() );
      blaze::orgrq( B, tauB.data() );

      if( A != B ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q reconstruction failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major reconstruction:\n" << A << "\n"
             << "   Column-major reconstruction:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Q reconstruction from a RQ decomposition (orgr2).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the Q reconstruction from a RQ decomposition for various data
// types. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void DecompositionTest::testOrgr2()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   test_ = "Reconstruction of Q from a RQ decomposition (orgr2)";

   {
      blaze::StaticMatrix<Type,2UL,5UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,2UL,5UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,2UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,2UL,blaze::rowVector> tauB;

      blaze::gerqf( A, tauA.data() );
      blaze::gerqf( B, tauB.data() );

      blaze::orgr2( A, tauA.data() );
      blaze::orgr2( B, tauB.data() );

      if( A != B ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q reconstruction failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major reconstruction:\n" << A << "\n"
             << "   Column-major reconstruction:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      blaze::StaticMatrix<Type,5UL,2UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,5UL,2UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,2UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,2UL,blaze::rowVector> tauB;

      blaze::gerqf( A, tauA.data() );
      blaze::gerqf( B, tauB.data() );

      blaze::orgr2( A, tauA.data() );
      blaze::orgr2( B, tauB.data() );

      if( A != B ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q reconstruction failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major reconstruction:\n" << A << "\n"
             << "   Column-major reconstruction:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Q reconstruction from a RQ decomposition (ungrq).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the Q reconstruction from a RQ decomposition for various data
// types. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void DecompositionTest::testUngrq()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   test_ = "Reconstruction of Q from a RQ decomposition (ungrq)";

   {
      blaze::StaticMatrix<Type,2UL,5UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,2UL,5UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,2UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,2UL,blaze::rowVector> tauB;

      blaze::gerqf( A, tauA.data() );
      blaze::gerqf( B, tauB.data() );

      blaze::ungrq( A, tauA.data() );
      blaze::ungrq( B, tauB.data() );

      if( A != B ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q reconstruction failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major reconstruction:\n" << A << "\n"
             << "   Column-major reconstruction:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      blaze::StaticMatrix<Type,5UL,2UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,5UL,2UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,2UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,2UL,blaze::rowVector> tauB;

      blaze::gerqf( A, tauA.data() );
      blaze::gerqf( B, tauB.data() );

      blaze::ungrq( A, tauA.data() );
      blaze::ungrq( B, tauB.data() );

      if( A != B ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q reconstruction failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major reconstruction:\n" << A << "\n"
             << "   Column-major reconstruction:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Q reconstruction from a RQ decomposition (ungr2).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the Q reconstruction from a RQ decomposition for various data
// types. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void DecompositionTest::testUngr2()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   test_ = "Reconstruction of Q from a RQ decomposition (ungr2)";

   {
      blaze::StaticMatrix<Type,2UL,5UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,2UL,5UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,2UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,2UL,blaze::rowVector> tauB;

      blaze::gerqf( A, tauA.data() );
      blaze::gerqf( B, tauB.data() );

      blaze::ungr2( A, tauA.data() );
      blaze::ungr2( B, tauB.data() );

      if( A != B ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q reconstruction failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major reconstruction:\n" << A << "\n"
             << "   Column-major reconstruction:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      blaze::StaticMatrix<Type,5UL,2UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,5UL,2UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,2UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,2UL,blaze::rowVector> tauB;

      blaze::gerqf( A, tauA.data() );
      blaze::gerqf( B, tauB.data() );

      blaze::ungr2( A, tauA.data() );
      blaze::ungr2( B, tauB.data() );

      if( A != B ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q reconstruction failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major reconstruction:\n" << A << "\n"
             << "   Column-major reconstruction:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the multiplication of Q from a RQ decomposition with a matrix (ormrq).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the multiplication of Q from a RQ decomposition with a matrix
// for various data types. In case an error is detected, a \a std::runtime_error exception is
// thrown.
*/
template< typename Type >
void DecompositionTest::testOrmrq()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   {
      test_ = "Multiplication of Q from a RQ decomposition with a matrix ('L', 'N')";

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,3UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,3UL,blaze::rowVector> tauB;

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> C1;
      randomize( C1 );

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor>    C2( C1 );
      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> C3( C1 );
      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> C4( C1 );

      blaze::gerqf( A, tauA.data() );
      blaze::gerqf( B, tauB.data() );

      blaze::ormrq( C1, A, 'L', 'N', tauA.data() );
      blaze::ormrq( C2, B, 'L', 'N', tauB.data() );
      blaze::ormrq( C3, A, 'L', 'N', tauA.data() );
      blaze::ormrq( C4, B, 'L', 'N', tauB.data() );

      if( C1 != C2 || C1 != C3 || C1 != C4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q multiplication failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major/row-major multiplication:\n" << C1 << "\n"
             << "   Row-major/column-major multiplication:\n" << C2 << "\n"
             << "   Column-major/row-major multiplication:\n" << C3 << "\n"
             << "   Column-major/column-major multiplication:\n" << C4 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Multiplication of Q from a RQ decomposition with a matrix ('L', 'T')";

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,3UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,3UL,blaze::rowVector> tauB;

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> C1;
      randomize( C1 );

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor>    C2( C1 );
      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> C3( C1 );
      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> C4( C1 );

      blaze::gerqf( A, tauA.data() );
      blaze::gerqf( B, tauB.data() );

      blaze::ormrq( C1, A, 'L', 'T', tauA.data() );
      blaze::ormrq( C2, B, 'L', 'T', tauB.data() );
      blaze::ormrq( C3, A, 'L', 'T', tauA.data() );
      blaze::ormrq( C4, B, 'L', 'T', tauB.data() );

      if( C1 != C2 || C1 != C3 || C1 != C4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q multiplication failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major/row-major multiplication:\n" << C1 << "\n"
             << "   Row-major/column-major multiplication:\n" << C2 << "\n"
             << "   Column-major/row-major multiplication:\n" << C3 << "\n"
             << "   Column-major/column-major multiplication:\n" << C4 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Multiplication of Q from a RQ decomposition with a matrix ('R', 'N')";

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,3UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,3UL,blaze::rowVector> tauB;

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> C1;
      randomize( C1 );

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor>    C2( C1 );
      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> C3( C1 );
      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> C4( C1 );

      blaze::gerqf( A, tauA.data() );
      blaze::gerqf( B, tauB.data() );

      blaze::ormrq( C1, A, 'R', 'N', tauA.data() );
      blaze::ormrq( C2, B, 'R', 'N', tauB.data() );
      blaze::ormrq( C3, A, 'R', 'N', tauA.data() );
      blaze::ormrq( C4, B, 'R', 'N', tauB.data() );

      if( C1 != C2 || C1 != C3 || C1 != C4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q multiplication failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major/row-major multiplication:\n" << C1 << "\n"
             << "   Row-major/column-major multiplication:\n" << C2 << "\n"
             << "   Column-major/row-major multiplication:\n" << C3 << "\n"
             << "   Column-major/column-major multiplication:\n" << C4 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Multiplication of Q from a RQ decomposition with a matrix ('R', 'T')";

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,3UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,3UL,blaze::rowVector> tauB;

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> C1;
      randomize( C1 );

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor>    C2( C1 );
      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> C3( C1 );
      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> C4( C1 );

      blaze::gerqf( A, tauA.data() );
      blaze::gerqf( B, tauB.data() );

      blaze::ormrq( C1, A, 'R', 'T', tauA.data() );
      blaze::ormrq( C2, B, 'R', 'T', tauB.data() );
      blaze::ormrq( C3, A, 'R', 'T', tauA.data() );
      blaze::ormrq( C4, B, 'R', 'T', tauB.data() );

      if( C1 != C2 || C1 != C3 || C1 != C4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q multiplication failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major/row-major multiplication:\n" << C1 << "\n"
             << "   Row-major/column-major multiplication:\n" << C2 << "\n"
             << "   Column-major/row-major multiplication:\n" << C3 << "\n"
             << "   Column-major/column-major multiplication:\n" << C4 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the multiplication of Q from a RQ decomposition with a matrix (unmrq).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the multiplication of Q from a RQ decomposition with a matrix
// for various data types. In case an error is detected, a \a std::runtime_error exception is
// thrown.
*/
template< typename Type >
void DecompositionTest::testUnmrq()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   {
      test_ = "Multiplication of Q from a RQ decomposition with a matrix ('L', 'N')";

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,3UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,3UL,blaze::rowVector> tauB;

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> C1;
      randomize( C1 );

      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> C2( C1 );

      blaze::gerqf( A, tauA.data() );
      blaze::gerqf( B, tauB.data() );

      blaze::unmrq( C1, A, 'L', 'N', tauA.data() );
      blaze::unmrq( C2, B, 'L', 'N', tauB.data() );

      if( C1 != C2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q multiplication failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major multiplication:\n" << C1 << "\n"
             << "   Column-major multiplication:\n" << C2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Multiplication of Q from a RQ decomposition with a matrix ('L', 'C')";

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,3UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,3UL,blaze::rowVector> tauB;

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> C1;
      randomize( C1 );

      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> C2( C1 );

      blaze::gerqf( A, tauA.data() );
      blaze::gerqf( B, tauB.data() );

      blaze::unmrq( C1, A, 'L', 'C', tauA.data() );
      blaze::unmrq( C2, B, 'L', 'C', tauB.data() );

      if( C1 != C2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q multiplication failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major multiplication:\n" << C1 << "\n"
             << "   Column-major multiplication:\n" << C2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Multiplication of Q from a RQ decomposition with a matrix ('R', 'N')";

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,3UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,3UL,blaze::rowVector> tauB;

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> C1;
      randomize( C1 );

      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> C2( C1 );

      blaze::gerqf( A, tauA.data() );
      blaze::gerqf( B, tauB.data() );

      blaze::unmrq( C1, A, 'R', 'N', tauA.data() );
      blaze::unmrq( C2, B, 'R', 'N', tauB.data() );

      if( C1 != C2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q multiplication failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major multiplication:\n" << C1 << "\n"
             << "   Column-major multiplication:\n" << C2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Multiplication of Q from a RQ decomposition with a matrix ('R', 'C')";

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,3UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,3UL,blaze::rowVector> tauB;

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> C1;
      randomize( C1 );

      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> C2( C1 );

      blaze::gerqf( A, tauA.data() );
      blaze::gerqf( B, tauB.data() );

      blaze::unmrq( C1, A, 'R', 'C', tauA.data() );
      blaze::unmrq( C2, B, 'R', 'C', tauB.data() );

      if( C1 != C2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q multiplication failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major multiplication:\n" << C1 << "\n"
             << "   Column-major multiplication:\n" << C2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the QL decomposition functions (geqlf).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the QL decomposition functions for various data types. In
// case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void DecompositionTest::testGeqlf()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   test_ = "QL decomposition";

   {
      blaze::StaticMatrix<Type,2UL,5UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,2UL,5UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,2UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,2UL,blaze::rowVector> tauB;

      blaze::geqlf( A, tauA.data() );
      blaze::geqlf( B, tauB.data() );

      if( A != B || tauA != conj( tauB ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: QL decomposition failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A << "\n"
             << "   Row-major projectors:\n" << tauA << "\n"
             << "   Column-major decomposition:\n" << B << "\n"
             << "   Column-major projectors:\n" << tauB << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      blaze::StaticMatrix<Type,5UL,2UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,5UL,2UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,2UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,2UL,blaze::rowVector> tauB;

      blaze::geqlf( A, tauA.data() );
      blaze::geqlf( B, tauB.data() );

      if( A != B || tauA != conj( tauB ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: QL decomposition failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A << "\n"
             << "   Row-major projectors:\n" << tauA << "\n"
             << "   Column-major decomposition:\n" << B << "\n"
             << "   Column-major projectors:\n" << tauB << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Q reconstruction from a QL decomposition (orgql).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the Q reconstruction from a QL decomposition for various data
// types. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void DecompositionTest::testOrgql()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   test_ = "Reconstruction of Q from a QL decomposition (orgql)";

   {
      blaze::StaticMatrix<Type,2UL,5UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,2UL,5UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,2UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,2UL,blaze::rowVector> tauB;

      blaze::geqlf( A, tauA.data() );
      blaze::geqlf( B, tauB.data() );

      blaze::orgql( A, tauA.data() );
      blaze::orgql( B, tauB.data() );

      if( A != B ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q reconstruction failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major reconstruction:\n" << A << "\n"
             << "   Column-major reconstruction:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      blaze::StaticMatrix<Type,5UL,2UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,5UL,2UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,2UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,2UL,blaze::rowVector> tauB;

      blaze::geqlf( A, tauA.data() );
      blaze::geqlf( B, tauB.data() );

      blaze::orgql( A, tauA.data() );
      blaze::orgql( B, tauB.data() );

      if( A != B ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q reconstruction failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major reconstruction:\n" << A << "\n"
             << "   Column-major reconstruction:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Q reconstruction from a QL decomposition (org2l).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the Q reconstruction from a QL decomposition for various data
// types. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void DecompositionTest::testOrg2l()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   test_ = "Reconstruction of Q from a QL decomposition (org2l)";

   {
      blaze::StaticMatrix<Type,2UL,5UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,2UL,5UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,2UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,2UL,blaze::rowVector> tauB;

      blaze::geqlf( A, tauA.data() );
      blaze::geqlf( B, tauB.data() );

      blaze::org2l( A, tauA.data() );
      blaze::org2l( B, tauB.data() );

      if( A != B ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q reconstruction failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major reconstruction:\n" << A << "\n"
             << "   Column-major reconstruction:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      blaze::StaticMatrix<Type,5UL,2UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,5UL,2UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,2UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,2UL,blaze::rowVector> tauB;

      blaze::geqlf( A, tauA.data() );
      blaze::geqlf( B, tauB.data() );

      blaze::org2l( A, tauA.data() );
      blaze::org2l( B, tauB.data() );

      if( A != B ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q reconstruction failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major reconstruction:\n" << A << "\n"
             << "   Column-major reconstruction:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Q reconstruction from a QL decomposition (ungql).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the Q reconstruction from a QL decomposition for various data
// types. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void DecompositionTest::testUngql()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   test_ = "Reconstruction of Q from a QL decomposition (ungql)";

   {
      blaze::StaticMatrix<Type,2UL,5UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,2UL,5UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,2UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,2UL,blaze::rowVector> tauB;

      blaze::geqlf( A, tauA.data() );
      blaze::geqlf( B, tauB.data() );

      blaze::ungql( A, tauA.data() );
      blaze::ungql( B, tauB.data() );

      if( A != B ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q reconstruction failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major reconstruction:\n" << A << "\n"
             << "   Column-major reconstruction:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      blaze::StaticMatrix<Type,5UL,2UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,5UL,2UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,2UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,2UL,blaze::rowVector> tauB;

      blaze::geqlf( A, tauA.data() );
      blaze::geqlf( B, tauB.data() );

      blaze::ungql( A, tauA.data() );
      blaze::ungql( B, tauB.data() );

      if( A != B ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q reconstruction failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major reconstruction:\n" << A << "\n"
             << "   Column-major reconstruction:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Q reconstruction from a QL decomposition (ung2l).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the Q reconstruction from a QL decomposition for various data
// types. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void DecompositionTest::testUng2l()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   test_ = "Reconstruction of Q from a QL decomposition (ung2l)";

   {
      blaze::StaticMatrix<Type,2UL,5UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,2UL,5UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,2UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,2UL,blaze::rowVector> tauB;

      blaze::geqlf( A, tauA.data() );
      blaze::geqlf( B, tauB.data() );

      blaze::ung2l( A, tauA.data() );
      blaze::ung2l( B, tauB.data() );

      if( A != B ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q reconstruction failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major reconstruction:\n" << A << "\n"
             << "   Column-major reconstruction:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      blaze::StaticMatrix<Type,5UL,2UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,5UL,2UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,2UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,2UL,blaze::rowVector> tauB;

      blaze::geqlf( A, tauA.data() );
      blaze::geqlf( B, tauB.data() );

      blaze::ung2l( A, tauA.data() );
      blaze::ung2l( B, tauB.data() );

      if( A != B ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q reconstruction failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major reconstruction:\n" << A << "\n"
             << "   Column-major reconstruction:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the multiplication of Q from a QL decomposition with a matrix (ormql).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the multiplication of Q from a QL decomposition with a matrix
// for various data types. In case an error is detected, a \a std::runtime_error exception is
// thrown.
*/
template< typename Type >
void DecompositionTest::testOrmql()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   {
      test_ = "Multiplication of Q from a QL decomposition with a matrix ('L', 'N')";

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,3UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,3UL,blaze::rowVector> tauB;

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> C1;
      randomize( C1 );

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor>    C2( C1 );
      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> C3( C1 );
      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> C4( C1 );

      blaze::geqlf( A, tauA.data() );
      blaze::geqlf( B, tauB.data() );

      blaze::ormql( C1, A, 'L', 'N', tauA.data() );
      blaze::ormql( C2, B, 'L', 'N', tauB.data() );
      blaze::ormql( C3, A, 'L', 'N', tauA.data() );
      blaze::ormql( C4, B, 'L', 'N', tauB.data() );

      if( C1 != C2 || C1 != C3 || C1 != C4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q multiplication failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major/row-major multiplication:\n" << C1 << "\n"
             << "   Row-major/column-major multiplication:\n" << C2 << "\n"
             << "   Column-major/row-major multiplication:\n" << C3 << "\n"
             << "   Column-major/column-major multiplication:\n" << C4 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Multiplication of Q from a QL decomposition with a matrix ('L', 'T')";

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,3UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,3UL,blaze::rowVector> tauB;

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> C1;
      randomize( C1 );

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor>    C2( C1 );
      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> C3( C1 );
      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> C4( C1 );

      blaze::geqlf( A, tauA.data() );
      blaze::geqlf( B, tauB.data() );

      blaze::ormql( C1, A, 'L', 'T', tauA.data() );
      blaze::ormql( C2, B, 'L', 'T', tauB.data() );
      blaze::ormql( C3, A, 'L', 'T', tauA.data() );
      blaze::ormql( C4, B, 'L', 'T', tauB.data() );

      if( C1 != C2 || C1 != C3 || C1 != C4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q multiplication failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major/row-major multiplication:\n" << C1 << "\n"
             << "   Row-major/column-major multiplication:\n" << C2 << "\n"
             << "   Column-major/row-major multiplication:\n" << C3 << "\n"
             << "   Column-major/column-major multiplication:\n" << C4 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Multiplication of Q from a QL decomposition with a matrix ('R', 'N')";

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,3UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,3UL,blaze::rowVector> tauB;

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> C1;
      randomize( C1 );

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor>    C2( C1 );
      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> C3( C1 );
      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> C4( C1 );

      blaze::geqlf( A, tauA.data() );
      blaze::geqlf( B, tauB.data() );

      blaze::ormql( C1, A, 'R', 'N', tauA.data() );
      blaze::ormql( C2, B, 'R', 'N', tauB.data() );
      blaze::ormql( C3, A, 'R', 'N', tauA.data() );
      blaze::ormql( C4, B, 'R', 'N', tauB.data() );

      if( C1 != C2 || C1 != C3 || C1 != C4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q multiplication failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major/row-major multiplication:\n" << C1 << "\n"
             << "   Row-major/column-major multiplication:\n" << C2 << "\n"
             << "   Column-major/row-major multiplication:\n" << C3 << "\n"
             << "   Column-major/column-major multiplication:\n" << C4 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Multiplication of Q from a QL decomposition with a matrix ('R', 'T')";

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,3UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,3UL,blaze::rowVector> tauB;

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> C1;
      randomize( C1 );

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor>    C2( C1 );
      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> C3( C1 );
      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> C4( C1 );

      blaze::geqlf( A, tauA.data() );
      blaze::geqlf( B, tauB.data() );

      blaze::ormql( C1, A, 'R', 'T', tauA.data() );
      blaze::ormql( C2, B, 'R', 'T', tauB.data() );
      blaze::ormql( C3, A, 'R', 'T', tauA.data() );
      blaze::ormql( C4, B, 'R', 'T', tauB.data() );

      if( C1 != C2 || C1 != C3 || C1 != C4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q multiplication failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major/row-major multiplication:\n" << C1 << "\n"
             << "   Row-major/column-major multiplication:\n" << C2 << "\n"
             << "   Column-major/row-major multiplication:\n" << C3 << "\n"
             << "   Column-major/column-major multiplication:\n" << C4 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the multiplication of Q from a QL decomposition with a matrix (unmql).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the multiplication of Q from a QL decomposition with a matrix
// for various data types. In case an error is detected, a \a std::runtime_error exception is
// thrown.
*/
template< typename Type >
void DecompositionTest::testUnmql()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   {
      test_ = "Multiplication of Q from a QL decomposition with a matrix ('L', 'N')";

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,3UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,3UL,blaze::rowVector> tauB;

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> C1;
      randomize( C1 );

      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> C2( C1 );

      blaze::geqlf( A, tauA.data() );
      blaze::geqlf( B, tauB.data() );

      blaze::unmql( C1, A, 'L', 'N', tauA.data() );
      blaze::unmql( C2, B, 'L', 'N', tauB.data() );

      if( C1 != C2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q multiplication failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major multiplication:\n" << C1 << "\n"
             << "   Column-major multiplication:\n" << C2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Multiplication of Q from a QL decomposition with a matrix ('L', 'C')";

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,3UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,3UL,blaze::rowVector> tauB;

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> C1;
      randomize( C1 );

      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> C2( C1 );

      blaze::geqlf( A, tauA.data() );
      blaze::geqlf( B, tauB.data() );

      blaze::unmql( C1, A, 'L', 'C', tauA.data() );
      blaze::unmql( C2, B, 'L', 'C', tauB.data() );

      if( C1 != C2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q multiplication failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major/row-major multiplication:\n" << C1 << "\n"
             << "   Column-major/row-major multiplication:\n" << C2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Multiplication of Q from a QL decomposition with a matrix ('R', 'N')";

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,3UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,3UL,blaze::rowVector> tauB;

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> C1;
      randomize( C1 );

      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> C2( C1 );

      blaze::geqlf( A, tauA.data() );
      blaze::geqlf( B, tauB.data() );

      blaze::unmql( C1, A, 'R', 'N', tauA.data() );
      blaze::unmql( C2, B, 'R', 'N', tauB.data() );

      if( C1 != C2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q multiplication failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major multiplication:\n" << C1 << "\n"
             << "   Column-major multiplication:\n" << C2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Multiplication of Q from a QL decomposition with a matrix ('R', 'C')";

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,3UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,3UL,blaze::rowVector> tauB;

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> C1;
      randomize( C1 );

      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> C2( C1 );

      blaze::geqlf( A, tauA.data() );
      blaze::geqlf( B, tauB.data() );

      blaze::unmql( C1, A, 'R', 'C', tauA.data() );
      blaze::unmql( C2, B, 'R', 'C', tauB.data() );

      if( C1 != C2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q multiplication failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major multiplication:\n" << C1 << "\n"
             << "   Column-major multiplication:\n" << C2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the LQ decomposition functions (gelqf).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the LQ decomposition functions for various data types. In
// case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void DecompositionTest::testGelqf()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   test_ = "LQ decomposition";

   {
      blaze::StaticMatrix<Type,2UL,5UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,2UL,5UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,2UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,2UL,blaze::rowVector> tauB;

      blaze::gelqf( A, tauA.data() );
      blaze::gelqf( B, tauB.data() );

      if( A != B || tauA != conj( tauB ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: LQ decomposition failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A << "\n"
             << "   Row-major projectors:\n" << tauA << "\n"
             << "   Column-major decomposition:\n" << B << "\n"
             << "   Column-major projectors:\n" << tauB << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      blaze::StaticMatrix<Type,5UL,2UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,5UL,2UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,2UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,2UL,blaze::rowVector> tauB;

      blaze::gelqf( A, tauA.data() );
      blaze::gelqf( B, tauB.data() );

      if( A != B || tauA != conj( tauB ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: LQ decomposition failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A << "\n"
             << "   Row-major projectors:\n" << tauA << "\n"
             << "   Column-major decomposition:\n" << B << "\n"
             << "   Column-major projectors:\n" << tauB << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Q reconstruction from a LQ decomposition (orglq).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the Q reconstruction from a LQ decomposition for various data
// types. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void DecompositionTest::testOrglq()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   test_ = "Reconstruction of Q from a LQ decomposition (orglq)";

   {
      blaze::StaticMatrix<Type,2UL,5UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,2UL,5UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,2UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,2UL,blaze::rowVector> tauB;

      blaze::gelqf( A, tauA.data() );
      blaze::gelqf( B, tauB.data() );

      blaze::orglq( A, tauA.data() );
      blaze::orglq( B, tauB.data() );

      if( A != B ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q reconstruction failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major reconstruction:\n" << A << "\n"
             << "   Column-major reconstruction:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      blaze::StaticMatrix<Type,5UL,2UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,5UL,2UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,2UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,2UL,blaze::rowVector> tauB;

      blaze::gelqf( A, tauA.data() );
      blaze::gelqf( B, tauB.data() );

      blaze::orglq( A, tauA.data() );
      blaze::orglq( B, tauB.data() );

      if( A != B ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q reconstruction failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major reconstruction:\n" << A << "\n"
             << "   Column-major reconstruction:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Q reconstruction from a LQ decomposition (orgl2).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the Q reconstruction from a LQ decomposition for various data
// types. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void DecompositionTest::testOrgl2()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   test_ = "Reconstruction of Q from a LQ decomposition (orgl2)";

   {
      blaze::StaticMatrix<Type,2UL,5UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,2UL,5UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,2UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,2UL,blaze::rowVector> tauB;

      blaze::gelqf( A, tauA.data() );
      blaze::gelqf( B, tauB.data() );

      blaze::orgl2( A, tauA.data() );
      blaze::orgl2( B, tauB.data() );

      if( A != B ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q reconstruction failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major reconstruction:\n" << A << "\n"
             << "   Column-major reconstruction:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      blaze::StaticMatrix<Type,5UL,2UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,5UL,2UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,2UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,2UL,blaze::rowVector> tauB;

      blaze::gelqf( A, tauA.data() );
      blaze::gelqf( B, tauB.data() );

      blaze::orgl2( A, tauA.data() );
      blaze::orgl2( B, tauB.data() );

      if( A != B ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q reconstruction failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major reconstruction:\n" << A << "\n"
             << "   Column-major reconstruction:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Q reconstruction from a LQ decomposition (unglq).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the Q reconstruction from a LQ decomposition for various data
// types. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void DecompositionTest::testUnglq()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   test_ = "Reconstruction of Q from a LQ decomposition (unglq)";

   {
      blaze::StaticMatrix<Type,2UL,5UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,2UL,5UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,2UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,2UL,blaze::rowVector> tauB;

      blaze::gelqf( A, tauA.data() );
      blaze::gelqf( B, tauB.data() );

      blaze::unglq( A, tauA.data() );
      blaze::unglq( B, tauB.data() );

      if( A != B ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q reconstruction failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major reconstruction:\n" << A << "\n"
             << "   Column-major reconstruction:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      blaze::StaticMatrix<Type,5UL,2UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,5UL,2UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,2UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,2UL,blaze::rowVector> tauB;

      blaze::gelqf( A, tauA.data() );
      blaze::gelqf( B, tauB.data() );

      blaze::unglq( A, tauA.data() );
      blaze::unglq( B, tauB.data() );

      if( A != B ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q reconstruction failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major reconstruction:\n" << A << "\n"
             << "   Column-major reconstruction:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Q reconstruction from a LQ decomposition (ungl2).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the Q reconstruction from a LQ decomposition for various data
// types. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void DecompositionTest::testUngl2()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   test_ = "Reconstruction of Q from a LQ decomposition (ungl2)";

   {
      blaze::StaticMatrix<Type,2UL,5UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,2UL,5UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,2UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,2UL,blaze::rowVector> tauB;

      blaze::gelqf( A, tauA.data() );
      blaze::gelqf( B, tauB.data() );

      blaze::ungl2( A, tauA.data() );
      blaze::ungl2( B, tauB.data() );

      if( A != B ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q reconstruction failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major reconstruction:\n" << A << "\n"
             << "   Column-major reconstruction:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      blaze::StaticMatrix<Type,5UL,2UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,5UL,2UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,2UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,2UL,blaze::rowVector> tauB;

      blaze::gelqf( A, tauA.data() );
      blaze::gelqf( B, tauB.data() );

      blaze::ungl2( A, tauA.data() );
      blaze::ungl2( B, tauB.data() );

      if( A != B ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q reconstruction failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major reconstruction:\n" << A << "\n"
             << "   Column-major reconstruction:\n" << B << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the multiplication of Q from a LQ decomposition with a matrix (ormlq).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the multiplication of Q from a LQ decomposition with a matrix
// for various data types. In case an error is detected, a \a std::runtime_error exception is
// thrown.
*/
template< typename Type >
void DecompositionTest::testOrmlq()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   {
      test_ = "Multiplication of Q from a LQ decomposition with a matrix ('L', 'N')";

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,3UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,3UL,blaze::rowVector> tauB;

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> C1;
      randomize( C1 );

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor>    C2( C1 );
      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> C3( C1 );
      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> C4( C1 );

      blaze::gelqf( A, tauA.data() );
      blaze::gelqf( B, tauB.data() );

      blaze::ormlq( C1, A, 'L', 'N', tauA.data() );
      blaze::ormlq( C2, B, 'L', 'N', tauB.data() );
      blaze::ormlq( C3, A, 'L', 'N', tauA.data() );
      blaze::ormlq( C4, B, 'L', 'N', tauB.data() );

      if( C1 != C2 || C1 != C3 || C1 != C4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q multiplication failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major/row-major multiplication:\n" << C1 << "\n"
             << "   Row-major/column-major multiplication:\n" << C2 << "\n"
             << "   Column-major/row-major multiplication:\n" << C3 << "\n"
             << "   Column-major/column-major multiplication:\n" << C4 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Multiplication of Q from a LQ decomposition with a matrix ('L', 'T')";

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,3UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,3UL,blaze::rowVector> tauB;

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> C1;
      randomize( C1 );

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor>    C2( C1 );
      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> C3( C1 );
      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> C4( C1 );

      blaze::gelqf( A, tauA.data() );
      blaze::gelqf( B, tauB.data() );

      blaze::ormlq( C1, A, 'L', 'T', tauA.data() );
      blaze::ormlq( C2, B, 'L', 'T', tauB.data() );
      blaze::ormlq( C3, A, 'L', 'T', tauA.data() );
      blaze::ormlq( C4, B, 'L', 'T', tauB.data() );

      if( C1 != C2 || C1 != C3 || C1 != C4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q multiplication failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major/row-major multiplication:\n" << C1 << "\n"
             << "   Row-major/column-major multiplication:\n" << C2 << "\n"
             << "   Column-major/row-major multiplication:\n" << C3 << "\n"
             << "   Column-major/column-major multiplication:\n" << C4 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Multiplication of Q from a LQ decomposition with a matrix ('R', 'N')";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,3UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,3UL,blaze::rowVector> tauB;

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> C1;
      randomize( C1 );

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor>    C2( C1 );
      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> C3( C1 );
      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> C4( C1 );

      blaze::gelqf( A, tauA.data() );
      blaze::gelqf( B, tauB.data() );

      blaze::ormlq( C1, A, 'R', 'N', tauA.data() );
      blaze::ormlq( C2, B, 'R', 'N', tauB.data() );
      blaze::ormlq( C3, A, 'R', 'N', tauA.data() );
      blaze::ormlq( C4, B, 'R', 'N', tauB.data() );

      if( C1 != C2 || C1 != C3 || C1 != C4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q multiplication failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major/row-major multiplication:\n" << C1 << "\n"
             << "   Row-major/column-major multiplication:\n" << C2 << "\n"
             << "   Column-major/row-major multiplication:\n" << C3 << "\n"
             << "   Column-major/column-major multiplication:\n" << C4 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Multiplication of Q from a LQ decomposition with a matrix ('R', 'T')";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,3UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,3UL,blaze::rowVector> tauB;

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> C1;
      randomize( C1 );

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor>    C2( C1 );
      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> C3( C1 );
      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> C4( C1 );

      blaze::gelqf( A, tauA.data() );
      blaze::gelqf( B, tauB.data() );

      blaze::ormlq( C1, A, 'R', 'T', tauA.data() );
      blaze::ormlq( C2, B, 'R', 'T', tauB.data() );
      blaze::ormlq( C3, A, 'R', 'T', tauA.data() );
      blaze::ormlq( C4, B, 'R', 'T', tauB.data() );

      if( C1 != C2 || C1 != C3 || C1 != C4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q multiplication failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major/row-major multiplication:\n" << C1 << "\n"
             << "   Row-major/column-major multiplication:\n" << C2 << "\n"
             << "   Column-major/row-major multiplication:\n" << C3 << "\n"
             << "   Column-major/column-major multiplication:\n" << C4 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the multiplication of Q from a LQ decomposition with a matrix (unmlq).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the multiplication of Q from a LQ decomposition with a matrix
// for various data types. In case an error is detected, a \a std::runtime_error exception is
// thrown.
*/
template< typename Type >
void DecompositionTest::testUnmlq()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   {
      test_ = "Multiplication of Q from a LQ decomposition with a matrix ('L', 'N')";

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,3UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,3UL,blaze::rowVector> tauB;

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> C1;
      randomize( C1 );

      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> C2( C1 );

      blaze::gelqf( A, tauA.data() );
      blaze::gelqf( B, tauB.data() );

      blaze::unmlq( C1, A, 'L', 'N', tauA.data() );
      blaze::unmlq( C2, B, 'L', 'N', tauB.data() );

      if( C1 != C2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q multiplication failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major multiplication:\n" << C1 << "\n"
             << "   Column-major multiplication:\n" << C2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Multiplication of Q from a LQ decomposition with a matrix ('L', 'C')";

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,3UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,3UL,blaze::rowVector> tauB;

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> C1;
      randomize( C1 );

      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> C2( C1 );

      blaze::gelqf( A, tauA.data() );
      blaze::gelqf( B, tauB.data() );

      blaze::unmlq( C1, A, 'L', 'C', tauA.data() );
      blaze::unmlq( C2, B, 'L', 'C', tauB.data() );

      if( C1 != C2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q multiplication failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major multiplication:\n" << C1 << "\n"
             << "   Column-major multiplication:\n" << C2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Multiplication of Q from a LQ decomposition with a matrix ('R', 'N')";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,3UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,3UL,blaze::rowVector> tauB;

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> C1;
      randomize( C1 );

      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> C2( C1 );

      blaze::gelqf( A, tauA.data() );
      blaze::gelqf( B, tauB.data() );

      blaze::unmlq( C1, A, 'R', 'N', tauA.data() );
      blaze::unmlq( C2, B, 'R', 'N', tauB.data() );

      if( C1 != C2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q multiplication failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major multiplication:\n" << C1 << "\n"
             << "   Column-major multiplication:\n" << C2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Multiplication of Q from a LQ decomposition with a matrix ('R', 'C')";

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> A;
      randomize( A );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> B( A );

      blaze::StaticVector<Type,3UL,blaze::rowVector> tauA;
      blaze::StaticVector<Type,3UL,blaze::rowVector> tauB;

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> C1;
      randomize( C1 );

      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> C2( C1 );

      blaze::gelqf( A, tauA.data() );
      blaze::gelqf( B, tauB.data() );

      blaze::unmlq( C1, A, 'R', 'C', tauA.data() );
      blaze::unmlq( C2, B, 'R', 'C', tauB.data() );

      if( C1 != C2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Q multiplication failed\n"
             << " Details:\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major multiplication:\n" << C1 << "\n"
             << "   Column-major multiplication:\n" << C2 << "\n";
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
   DecompositionTest();
}
//*************************************************************************************************




//=================================================================================================
//
//  MACRO DEFINITIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Macro for the execution of the LAPACK decomposition test.
*/
#define RUN_LAPACK_DECOMPOSITION_TEST \
   blazetest::mathtest::lapack::runTest()
/*! \endcond */
//*************************************************************************************************

} // namespace lapack

} // namespace mathtest

} // namespace blazetest

#endif
