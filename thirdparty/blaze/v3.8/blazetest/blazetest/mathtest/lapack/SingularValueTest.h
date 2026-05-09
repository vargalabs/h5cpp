//=================================================================================================
/*!
//  \file blazetest/mathtest/lapack/SingularValueTest.h
//  \brief Header file for the LAPACK singular value test
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

#ifndef _BLAZETEST_MATHTEST_LAPACK_SINGULARVALUETEST_H_
#define _BLAZETEST_MATHTEST_LAPACK_SINGULARVALUETEST_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <sstream>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <blaze/math/LAPACK.h>
#include <blaze/math/StaticMatrix.h>
#include <blaze/math/StaticVector.h>
#include <blaze/math/typetraits/UnderlyingElement.h>
#include <blaze/util/Random.h>
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
class SingularValueTest
{
 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit SingularValueTest();
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
   template< typename Type > void testGesvd();
   template< typename Type > void testGesdd();
   template< typename Type > void testGesvdx();
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
/*!\brief Test of the singular value decomposition functions (gesvd).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the singular value decomposition functions for various data
// types. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void SingularValueTest::testGesvd()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   //=====================================================================================
   // gesvd( DenseMatrix, DenseVector, char, char )
   //=====================================================================================

   {
      test_ = "gesvd( DenseMatrix, DenseVector, char, char ) (3x5, 'N', 'N')";

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> A1;
      randomize( A1 );

      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> A2( A1 );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s1;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s2;

      blaze::gesvd( A1, s1, 'N', 'N' );
      blaze::gesvd( A2, s2, 'N', 'N' );

      if( s1 != s2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Singular value decomposition failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A1 << "\n"
             << "   Row-major singular values:\n" << s1 << "\n"
             << "   Column-major decomposition:\n" << A2 << "\n"
             << "   Column-major singular values:\n" << s2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "gesvd( DenseMatrix, DenseVector, char, char ) (5x3, 'N', 'N')";

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> A1;
      randomize( A1 );

      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> A2( A1 );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s1;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s2;

      blaze::gesvd( A1, s1, 'N', 'N' );
      blaze::gesvd( A2, s2, 'N', 'N' );

      if( s1 != s2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Singular value decomposition failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A1 << "\n"
             << "   Row-major singular values:\n" << s1 << "\n"
             << "   Column-major decomposition:\n" << A2 << "\n"
             << "   Column-major singular values:\n" << s2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "gesvd( DenseMatrix, DenseVector, char, char ) (3x5, 'N', 'O')";

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> A1;
      randomize( A1 );

      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> A2( A1 );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s1;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s2;

      blaze::gesvd( A1, s1, 'N', 'O' );
      blaze::gesvd( A2, s2, 'N', 'O' );

      if( s1 != s2 || abs( A1 ) != abs( A2 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Singular value decomposition failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A1 << "\n"
             << "   Row-major singular values:\n" << s1 << "\n"
             << "   Column-major decomposition:\n" << A2 << "\n"
             << "   Column-major singular values:\n" << s2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "gesvd( DenseMatrix, DenseVector, char, char ) (5x3, 'N', 'O')";

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> A1;
      randomize( A1 );

      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> A2( A1 );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s1;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s2;

      blaze::gesvd( A1, s1, 'N', 'O' );
      blaze::gesvd( A2, s2, 'N', 'O' );

      if( s1 != s2 || abs( A1 ) != abs( A2 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Singular value decomposition failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A1 << "\n"
             << "   Row-major singular values:\n" << s1 << "\n"
             << "   Column-major decomposition:\n" << A2 << "\n"
             << "   Column-major singular values:\n" << s2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "gesvd( DenseMatrix, DenseVector, char, char ) (3x5, 'O', 'N')";

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> A1;
      randomize( A1 );

      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> A2( A1 );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s1;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s2;

      blaze::gesvd( A1, s1, 'O', 'N' );
      blaze::gesvd( A2, s2, 'O', 'N' );

      if( s1 != s2 || abs( A1 ) != abs( A2 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Singular value decomposition failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A1 << "\n"
             << "   Row-major singular values:\n" << s1 << "\n"
             << "   Column-major decomposition:\n" << A2 << "\n"
             << "   Column-major singular values:\n" << s2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "gesvd( DenseMatrix, DenseVector, char, char ) (5x3, 'O', 'N')";

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> A1;
      randomize( A1 );

      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> A2( A1 );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s1;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s2;

      blaze::gesvd( A1, s1, 'O', 'N' );
      blaze::gesvd( A2, s2, 'O', 'N' );

      if( s1 != s2 || abs( A1 ) != abs( A2 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Singular value decomposition failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A1 << "\n"
             << "   Row-major singular values:\n" << s1 << "\n"
             << "   Column-major decomposition:\n" << A2 << "\n"
             << "   Column-major singular values:\n" << s2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // gesvd( DenseMatrix, DenseMatrix, DenseVector, char, char )
   //=====================================================================================

   {
      test_ = "gesvd( DenseMatrix, DenseMatrix, DenseVector, char, char ) (3x5, 'N', 'N')";

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> A1;
      randomize( A1 );
      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> A2( A1 );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s1;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s2;

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> U1;

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> U2;

      blaze::gesvd( A1, U1, s1, 'N', 'N' );
      blaze::gesvd( A2, U2, s2, 'N', 'N' );

      if( s1 != s2 || !isDefault( U1 ) || !isDefault( U2 ) ) {
         std::ostringstream oss;
         oss.precision( 30 );
         oss << " Test: " << test_ << "\n"
             << " Error: Singular value decomposition failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A1 << "\n"
             << "   Row-major singular values:\n" << s1 << "\n"
             << "   Row-major left singular values:\n" << U1 << "\n"
             << "   Column-major decomposition:\n" << A2 << "\n"
             << "   Column-major singular values:\n" << s2 << "\n"
             << "   Column-major left singular values:\n" << U2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "gesvd( DenseMatrix, DenseMatrix, DenseVector, char, char ) (5x3, 'N', 'N')";

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> A1;
      randomize( A1 );
      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> A2( A1 );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s1;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s2;

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> U1;

      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> U2;

      blaze::gesvd( A1, U1, s1, 'N', 'N' );
      blaze::gesvd( A2, U2, s2, 'N', 'N' );

      if( s1 != s2 || !isDefault( U1 ) || !isDefault( U2 ) ) {
         std::ostringstream oss;
         oss.precision( 30 );
         oss << " Test: " << test_ << "\n"
             << " Error: Singular value decomposition failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A1 << "\n"
             << "   Row-major singular values:\n" << s1 << "\n"
             << "   Row-major left singular values:\n" << U1 << "\n"
             << "   Column-major decomposition:\n" << A2 << "\n"
             << "   Column-major singular values:\n" << s2 << "\n"
             << "   Column-major left singular values:\n" << U2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "gesvd( DenseMatrix, DenseMatrix, DenseVector, char, char ) (3x5, 'S', 'O')";

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> A1;
      randomize( A1 );
      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> A2( A1 );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s1;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s2;

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> U1;

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> U2;

      blaze::gesvd( A1, U1, s1, 'S', 'O' );
      blaze::gesvd( A2, U2, s2, 'S', 'O' );

      if( s1 != s2 || abs( U1 ) != abs( U2 ) || abs( A1 ) != abs( A2 ) ) {
         std::ostringstream oss;
         oss.precision( 30 );
         oss << " Test: " << test_ << "\n"
             << " Error: Singular value decomposition failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A1 << "\n"
             << "   Row-major singular values:\n" << s1 << "\n"
             << "   Row-major left singular values:\n" << U1 << "\n"
             << "   Column-major decomposition:\n" << A2 << "\n"
             << "   Column-major singular values:\n" << s2 << "\n"
             << "   Column-major left singular values:\n" << U2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "gesvd( DenseMatrix, DenseMatrix, DenseVector, char, char ) (5x3, 'S', 'O')";

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> A1;
      randomize( A1 );
      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> A2( A1 );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s1;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s2;

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> U1;

      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> U2;

      blaze::gesvd( A1, U1, s1, 'S', 'O' );
      blaze::gesvd( A2, U2, s2, 'S', 'O' );

      if( s1 != s2 || abs( U1 ) != abs( U2 ) || abs( A1 ) != abs( A2 ) ) {
         std::ostringstream oss;
         oss.precision( 30 );
         oss << " Test: " << test_ << "\n"
             << " Error: Singular value decomposition failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A1 << "\n"
             << "   Row-major singular values:\n" << s1 << "\n"
             << "   Row-major left singular values:\n" << U1 << "\n"
             << "   Column-major decomposition:\n" << A2 << "\n"
             << "   Column-major singular values:\n" << s2 << "\n"
             << "   Column-major left singular values:\n" << U2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "gesvd( DenseMatrix, DenseMatrix, DenseVector, char, char ) (3x5, 'A', 'O')";

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> A1;
      randomize( A1 );
      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> A2( A1 );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s1;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s2;

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> U1;

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> U2;

      blaze::gesvd( A1, U1, s1, 'A', 'O' );
      blaze::gesvd( A2, U2, s2, 'A', 'O' );

      if( s1 != s2 || abs( U1 ) != abs( U2 ) || abs( A1 ) != abs( A2 ) ) {
         std::ostringstream oss;
         oss.precision( 30 );
         oss << " Test: " << test_ << "\n"
             << " Error: Singular value decomposition failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A1 << "\n"
             << "   Row-major singular values:\n" << s1 << "\n"
             << "   Row-major left singular values:\n" << U1 << "\n"
             << "   Column-major decomposition:\n" << A2 << "\n"
             << "   Column-major singular values:\n" << s2 << "\n"
             << "   Column-major left singular values:\n" << U2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "gesvd( DenseMatrix, DenseMatrix, DenseVector, char, char ) (5x3, 'A', 'O')";

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> A1;
      randomize( A1 );
      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> A2( A1 );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s1;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s2;

      blaze::StaticMatrix<Type,5UL,5UL,blaze::rowMajor> U1;

      blaze::StaticMatrix<Type,5UL,5UL,blaze::columnMajor> U2;

      blaze::gesvd( A1, U1, s1, 'A', 'O' );
      blaze::gesvd( A2, U2, s2, 'A', 'O' );

      if( s1 != s2 || abs( U1 ) != abs( U2 ) || abs( A1 ) != abs( A2 ) ) {
         std::ostringstream oss;
         oss.precision( 30 );
         oss << " Test: " << test_ << "\n"
             << " Error: Singular value decomposition failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A1 << "\n"
             << "   Row-major singular values:\n" << s1 << "\n"
             << "   Row-major left singular values:\n" << U1 << "\n"
             << "   Column-major decomposition:\n" << A2 << "\n"
             << "   Column-major singular values:\n" << s2 << "\n"
             << "   Column-major left singular values:\n" << U2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // gesvd( DenseMatrix, DenseVector, DenseMatrix, char, char )
   //=====================================================================================

   {
      test_ = "gesvd( DenseMatrix, DenseVector, DenseMatrix, char, char ) (3x5, 'N', 'N')";

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> A1;
      randomize( A1 );
      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> A2( A1 );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s1;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s2;

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> V1;

      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> V2;

      blaze::gesvd( A1, s1, V1, 'N', 'N' );
      blaze::gesvd( A2, s2, V2, 'N', 'N' );

      if( s1 != s2 || !isDefault( V1 ) || !isDefault( V2 ) ) {
         std::ostringstream oss;
         oss.precision( 30 );
         oss << " Test: " << test_ << "\n"
             << " Error: Singular value decomposition failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A1 << "\n"
             << "   Row-major singular values:\n" << s1 << "\n"
             << "   Row-major right singular values:\n" << V1 << "\n"
             << "   Column-major decomposition:\n" << A2 << "\n"
             << "   Column-major singular values:\n" << s2 << "\n"
             << "   Column-major right singular values:\n" << V2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "gesvd( DenseMatrix, DenseVector, DenseMatrix, char, char ) (5x3, 'N', 'N')";

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> A1;
      randomize( A1 );
      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> A2( A1 );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s1;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s2;

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> V1;

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> V2;

      blaze::gesvd( A1, s1, V1, 'N', 'N' );
      blaze::gesvd( A2, s2, V2, 'N', 'N' );

      if( s1 != s2 || !isDefault( V1 ) || !isDefault( V2 ) ) {
         std::ostringstream oss;
         oss.precision( 30 );
         oss << " Test: " << test_ << "\n"
             << " Error: Singular value decomposition failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A1 << "\n"
             << "   Row-major singular values:\n" << s1 << "\n"
             << "   Row-major right singular values:\n" << V1 << "\n"
             << "   Column-major decomposition:\n" << A2 << "\n"
             << "   Column-major singular values:\n" << s2 << "\n"
             << "   Column-major right singular values:\n" << V2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "gesvd( DenseMatrix, DenseVector, DenseMatrix, char, char ) (3x5, 'O', 'S')";

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> A1;
      randomize( A1 );
      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> A2( A1 );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s1;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s2;

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> V1;

      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> V2;

      blaze::gesvd( A1, s1, V1, 'O', 'S' );
      blaze::gesvd( A2, s2, V2, 'O', 'S' );

      if( s1 != s2 || abs( A1 ) != abs( A2 ) || abs( V1 ) != abs( V2 ) ) {
         std::ostringstream oss;
         oss.precision( 30 );
         oss << " Test: " << test_ << "\n"
             << " Error: Singular value decomposition failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A1 << "\n"
             << "   Row-major singular values:\n" << s1 << "\n"
             << "   Row-major right singular values:\n" << V1 << "\n"
             << "   Column-major decomposition:\n" << A2 << "\n"
             << "   Column-major singular values:\n" << s2 << "\n"
             << "   Column-major right singular values:\n" << V2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "gesvd( DenseMatrix, DenseVector, DenseMatrix, char, char ) (5x3, 'O', 'S')";

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> A1;
      randomize( A1 );
      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> A2( A1 );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s1;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s2;

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> V1;

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> V2;

      blaze::gesvd( A1, s1, V1, 'O', 'S' );
      blaze::gesvd( A2, s2, V2, 'O', 'S' );

      if( s1 != s2 || abs( A1 ) != abs( A2 ) || abs( V1 ) != abs( V2 ) ) {
         std::ostringstream oss;
         oss.precision( 30 );
         oss << " Test: " << test_ << "\n"
             << " Error: Singular value decomposition failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A1 << "\n"
             << "   Row-major singular values:\n" << s1 << "\n"
             << "   Row-major right singular values:\n" << V1 << "\n"
             << "   Column-major decomposition:\n" << A2 << "\n"
             << "   Column-major singular values:\n" << s2 << "\n"
             << "   Column-major right singular values:\n" << V2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "gesvd( DenseMatrix, DenseVector, DenseMatrix, char, char ) (3x5, 'O', 'A')";

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> A1;
      randomize( A1 );
      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> A2( A1 );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s1;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s2;

      blaze::StaticMatrix<Type,5UL,5UL,blaze::rowMajor> V1;

      blaze::StaticMatrix<Type,5UL,5UL,blaze::columnMajor> V2;

      blaze::gesvd( A1, s1, V1, 'O', 'A' );
      blaze::gesvd( A2, s2, V2, 'O', 'A' );

      if( s1 != s2 || abs( A1 ) != abs( A2 ) || abs( V1 ) != abs( V2 ) ) {
         std::ostringstream oss;
         oss.precision( 30 );
         oss << " Test: " << test_ << "\n"
             << " Error: Singular value decomposition failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A1 << "\n"
             << "   Row-major singular values:\n" << s1 << "\n"
             << "   Row-major right singular values:\n" << V1 << "\n"
             << "   Column-major decomposition:\n" << A2 << "\n"
             << "   Column-major singular values:\n" << s2 << "\n"
             << "   Column-major right singular values:\n" << V2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "gesvd( DenseMatrix, DenseVector, DenseMatrix, char, char ) (5x3, 'O', 'A')";

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> A1;
      randomize( A1 );
      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> A2( A1 );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s1;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s2;

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> V1;

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> V2;

      blaze::gesvd( A1, s1, V1, 'O', 'A' );
      blaze::gesvd( A2, s2, V2, 'O', 'A' );

      if( s1 != s2 || abs( A1 ) != abs( A2 ) || abs( V1 ) != abs( V2 ) ) {
         std::ostringstream oss;
         oss.precision( 30 );
         oss << " Test: " << test_ << "\n"
             << " Error: Singular value decomposition failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A1 << "\n"
             << "   Row-major singular values:\n" << s1 << "\n"
             << "   Row-major right singular values:\n" << V1 << "\n"
             << "   Column-major decomposition:\n" << A2 << "\n"
             << "   Column-major singular values:\n" << s2 << "\n"
             << "   Column-major right singular values:\n" << V2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // gesvd( DenseMatrix, DenseMatrix, DenseVector, DenseMatrix, char, char )
   //=====================================================================================

   {
      test_ = "gesvd( DenseMatrix, DenseMatrix, DenseVector, DenseMatrix, char, char ) (3x5, 'N', 'N')";

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> A1;
      randomize( A1 );
      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> A2( A1 );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s1;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s2;

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> U1;
      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> V1;

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> U2;
      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> V2;

      blaze::gesvd( A1, U1, s1, V1, 'N', 'N' );
      blaze::gesvd( A2, U2, s2, V2, 'N', 'N' );

      if( s1 != s2 || !isDefault( U1 ) || !isDefault( U2 ) || !isDefault( V1 ) || !isDefault( V2 ) ) {
         std::ostringstream oss;
         oss.precision( 30 );
         oss << " Test: " << test_ << "\n"
             << " Error: Singular value decomposition failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A1 << "\n"
             << "   Row-major singular values:\n" << s1 << "\n"
             << "   Row-major left singular values:\n" << U1 << "\n"
             << "   Row-major right singular values:\n" << V1 << "\n"
             << "   Column-major decomposition:\n" << A2 << "\n"
             << "   Column-major singular values:\n" << s2 << "\n"
             << "   Column-major left singular values:\n" << U2 << "\n"
             << "   Column-major right singular values:\n" << V2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "gesvd( DenseMatrix, DenseMatrix, DenseVector, DenseMatrix, char, char ) (5x3, 'N', 'N')";

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> A1;
      randomize( A1 );
      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> A2( A1 );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s1;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s2;

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> U1;
      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> V1;

      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> U2;
      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> V2;

      blaze::gesvd( A1, U1, s1, V1, 'N', 'N' );
      blaze::gesvd( A2, U2, s2, V2, 'N', 'N' );

      if( s1 != s2 || !isDefault( U1 ) || !isDefault( U2 ) || !isDefault( V1 ) || !isDefault( V2 ) ) {
         std::ostringstream oss;
         oss.precision( 30 );
         oss << " Test: " << test_ << "\n"
             << " Error: Singular value decomposition failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A1 << "\n"
             << "   Row-major singular values:\n" << s1 << "\n"
             << "   Row-major left singular values:\n" << U1 << "\n"
             << "   Row-major right singular values:\n" << V1 << "\n"
             << "   Column-major decomposition:\n" << A2 << "\n"
             << "   Column-major singular values:\n" << s2 << "\n"
             << "   Column-major left singular values:\n" << U2 << "\n"
             << "   Column-major right singular values:\n" << V2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "gesvd( DenseMatrix, DenseMatrix, DenseVector, DenseMatrix, char, char ) (3x5, 'S', 'S')";

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> A1;
      randomize( A1 );
      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> A2( A1 );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s1;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s2;

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> U1;
      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> V1;

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> U2;
      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> V2;

      blaze::gesvd( A1, U1, s1, V1, 'S', 'S' );
      blaze::gesvd( A2, U2, s2, V2, 'S', 'S' );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> S1;
      S1(0,0) = s1[0];
      S1(1,1) = s1[1];
      S1(2,2) = s1[2];

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> S2;
      S2(0,0) = s2[0];
      S2(1,1) = s2[1];
      S2(2,2) = s2[2];

      if( s1 != s2 || ( U1 * S1 * V1 ) != ( U2 * S2 * V2 ) ) {
         std::ostringstream oss;
         oss.precision( 30 );
         oss << " Test: " << test_ << "\n"
             << " Error: Singular value decomposition failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A1 << "\n"
             << "   Row-major singular values:\n" << s1 << "\n"
             << "   Row-major left singular values:\n" << U1 << "\n"
             << "   Row-major right singular values:\n" << V1 << "\n"
             << "   Column-major decomposition:\n" << A2 << "\n"
             << "   Column-major singular values:\n" << s2 << "\n"
             << "   Column-major left singular values:\n" << U2 << "\n"
             << "   Column-major right singular values:\n" << V2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "gesvd( DenseMatrix, DenseMatrix, DenseVector, DenseMatrix, char, char ) (5x3, 'S', 'S')";

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> A1;
      randomize( A1 );
      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> A2( A1 );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s1;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s2;

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> U1;
      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> V1;

      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> U2;
      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> V2;

      blaze::gesvd( A1, U1, s1, V1, 'S', 'S' );
      blaze::gesvd( A2, U2, s2, V2, 'S', 'S' );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> S1;
      S1(0,0) = s1[0];
      S1(1,1) = s1[1];
      S1(2,2) = s1[2];

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> S2;
      S2(0,0) = s2[0];
      S2(1,1) = s2[1];
      S2(2,2) = s2[2];

      if( s1 != s2 || ( U1 * S1 * V1 ) != ( U2 * S2 * V2 ) ) {
         std::ostringstream oss;
         oss.precision( 30 );
         oss << " Test: " << test_ << "\n"
             << " Error: Singular value decomposition failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A1 << "\n"
             << "   Row-major singular values:\n" << s1 << "\n"
             << "   Row-major left singular values:\n" << U1 << "\n"
             << "   Row-major right singular values:\n" << V1 << "\n"
             << "   Column-major decomposition:\n" << A2 << "\n"
             << "   Column-major singular values:\n" << s2 << "\n"
             << "   Column-major left singular values:\n" << U2 << "\n"
             << "   Column-major right singular values:\n" << V2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "gesvd( DenseMatrix, DenseMatrix, DenseVector, DenseMatrix, char, char ) (3x5, 'A', 'A')";

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> A1;
      randomize( A1 );
      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> A2( A1 );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s1;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s2;

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> U1;
      blaze::StaticMatrix<Type,5UL,5UL,blaze::rowMajor> V1;

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> U2;
      blaze::StaticMatrix<Type,5UL,5UL,blaze::columnMajor> V2;

      blaze::gesvd( A1, U1, s1, V1, 'A', 'A' );
      blaze::gesvd( A2, U2, s2, V2, 'A', 'A' );

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> S1;
      S1(0,0) = s1[0];
      S1(1,1) = s1[1];
      S1(2,2) = s1[2];

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> S2;
      S2(0,0) = s2[0];
      S2(1,1) = s2[1];
      S2(2,2) = s2[2];

      if( s1 != s2 || ( U1 * S1 * V1 ) != ( U2 * S2 * V2 ) ) {
         std::ostringstream oss;
         oss.precision( 30 );
         oss << " Test: " << test_ << "\n"
             << " Error: Singular value decomposition failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A1 << "\n"
             << "   Row-major singular values:\n" << s1 << "\n"
             << "   Row-major left singular values:\n" << U1 << "\n"
             << "   Row-major right singular values:\n" << V1 << "\n"
             << "   Column-major decomposition:\n" << A2 << "\n"
             << "   Column-major singular values:\n" << s2 << "\n"
             << "   Column-major left singular values:\n" << U2 << "\n"
             << "   Column-major right singular values:\n" << V2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "gesvd( DenseMatrix, DenseMatrix, DenseVector, DenseMatrix, char, char ) (5x3, 'A', 'A')";

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> A1;
      randomize( A1 );
      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> A2( A1 );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s1;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s2;

      blaze::StaticMatrix<Type,5UL,5UL,blaze::rowMajor> U1;
      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> V1;

      blaze::StaticMatrix<Type,5UL,5UL,blaze::columnMajor> U2;
      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> V2;

      blaze::gesvd( A1, U1, s1, V1, 'A', 'A' );
      blaze::gesvd( A2, U2, s2, V2, 'A', 'A' );

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> S1;
      S1(0,0) = s1[0];
      S1(1,1) = s1[1];
      S1(2,2) = s1[2];

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> S2;
      S2(0,0) = s2[0];
      S2(1,1) = s2[1];
      S2(2,2) = s2[2];

      if( s1 != s2 || ( U1 * S1 * V1 ) != ( U2 * S2 * V2 ) ) {
         std::ostringstream oss;
         oss.precision( 30 );
         oss << " Test: " << test_ << "\n"
             << " Error: Singular value decomposition failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A1 << "\n"
             << "   Row-major singular values:\n" << s1 << "\n"
             << "   Row-major left singular values:\n" << U1 << "\n"
             << "   Row-major right singular values:\n" << V1 << "\n"
             << "   Column-major decomposition:\n" << A2 << "\n"
             << "   Column-major singular values:\n" << s2 << "\n"
             << "   Column-major left singular values:\n" << U2 << "\n"
             << "   Column-major right singular values:\n" << V2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the singular value decomposition functions (gesdd).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the singular value decomposition functions for various data
// types. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void SingularValueTest::testGesdd()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   //=====================================================================================
   // gesdd( DenseMatrix, DenseVector )
   //=====================================================================================

   {
      test_ = "gesdd( DenseMatrix, DenseVector, char ) (3x5)";

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> A1;
      randomize( A1 );

      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> A2( A1 );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s1;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s2;

      blaze::gesdd( A1, s1 );
      blaze::gesdd( A2, s2 );

      if( s1 != s2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Singular value decomposition failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A1 << "\n"
             << "   Row-major singular values:\n" << s1 << "\n"
             << "   Column-major decomposition:\n" << A2 << "\n"
             << "   Column-major singular values:\n" << s2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "gesdd( DenseMatrix, DenseVector, char ) (5x3)";

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> A1;
      randomize( A1 );

      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> A2( A1 );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s1;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s2;

      blaze::gesdd( A1, s1 );
      blaze::gesdd( A2, s2 );

      if( s1 != s2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Singular value decomposition failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A1 << "\n"
             << "   Row-major singular values:\n" << s1 << "\n"
             << "   Column-major decomposition:\n" << A2 << "\n"
             << "   Column-major singular values:\n" << s2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // gesdd( DenseMatrix, DenseMatrix, DenseVector, char )
   //=====================================================================================

   {
      test_ = "gesdd( DenseMatrix, DenseMatrix, DenseVector, char ) (3x5, 'N')";

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> A1;
      randomize( A1 );
      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> A2( A1 );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s1;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s2;

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> U1;

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> U2;

      blaze::gesdd( A1, U1, s1, 'N' );
      blaze::gesdd( A2, U2, s2, 'N' );

      if( s1 != s2 || !isDefault( U1 ) || !isDefault( U2 ) ) {
         std::ostringstream oss;
         oss.precision( 30 );
         oss << " Test: " << test_ << "\n"
             << " Error: Singular value decomposition failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A1 << "\n"
             << "   Row-major singular values:\n" << s1 << "\n"
             << "   Row-major left singular values:\n" << U1 << "\n"
             << "   Column-major decomposition:\n" << A2 << "\n"
             << "   Column-major singular values:\n" << s2 << "\n"
             << "   Column-major left singular values:\n" << U2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "gesdd( DenseMatrix, DenseMatrix, DenseVector, char ) (5x3, 'N')";

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> A1;
      randomize( A1 );
      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> A2( A1 );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s1;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s2;

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> U1;

      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> U2;

      blaze::gesdd( A1, U1, s1, 'N' );
      blaze::gesdd( A2, U2, s2, 'N' );

      if( s1 != s2 || !isDefault( U1 ) || !isDefault( U2 ) ) {
         std::ostringstream oss;
         oss.precision( 30 );
         oss << " Test: " << test_ << "\n"
             << " Error: Singular value decomposition failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A1 << "\n"
             << "   Row-major singular values:\n" << s1 << "\n"
             << "   Row-major left singular values:\n" << U1 << "\n"
             << "   Column-major decomposition:\n" << A2 << "\n"
             << "   Column-major singular values:\n" << s2 << "\n"
             << "   Column-major left singular values:\n" << U2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "gesdd( DenseMatrix, DenseMatrix, DenseVector, char ) (3x5, 'O')";

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> A1;
      randomize( A1 );
      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> A2( A1 );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s1;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s2;

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> U1;

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> U2;

      blaze::gesdd( A1, U1, s1, 'O' );
      blaze::gesdd( A2, U2, s2, 'O' );

      if( s1 != s2 || abs( U1 ) != abs( U2 ) || abs( A1 ) != abs( A2 ) ) {
         std::ostringstream oss;
         oss.precision( 30 );
         oss << " Test: " << test_ << "\n"
             << " Error: Singular value decomposition failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A1 << "\n"
             << "   Row-major singular values:\n" << s1 << "\n"
             << "   Row-major left singular values:\n" << U1 << "\n"
             << "   Column-major decomposition:\n" << A2 << "\n"
             << "   Column-major singular values:\n" << s2 << "\n"
             << "   Column-major left singular values:\n" << U2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // gesdd( DenseMatrix, DenseVector, DenseMatrix, char )
   //=====================================================================================

   {
      test_ = "gesdd( DenseMatrix, DenseVector, DenseMatrix, char ) (3x5, 'N')";

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> A1;
      randomize( A1 );
      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> A2( A1 );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s1;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s2;

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> V1;

      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> V2;

      blaze::gesdd( A1, s1, V1, 'N' );
      blaze::gesdd( A2, s2, V2, 'N' );

      if( s1 != s2 || !isDefault( V1 ) || !isDefault( V2 ) ) {
         std::ostringstream oss;
         oss.precision( 30 );
         oss << " Test: " << test_ << "\n"
             << " Error: Singular value decomposition failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A1 << "\n"
             << "   Row-major singular values:\n" << s1 << "\n"
             << "   Row-major right singular values:\n" << V1 << "\n"
             << "   Column-major decomposition:\n" << A2 << "\n"
             << "   Column-major singular values:\n" << s2 << "\n"
             << "   Column-major right singular values:\n" << V2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "gesdd( DenseMatrix, DenseVector, DenseMatrix, char ) (5x3, 'N')";

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> A1;
      randomize( A1 );
      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> A2( A1 );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s1;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s2;

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> V1;

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> V2;

      blaze::gesdd( A1, s1, V1, 'N' );
      blaze::gesdd( A2, s2, V2, 'N' );

      if( s1 != s2 || !isDefault( V1 ) || !isDefault( V2 ) ) {
         std::ostringstream oss;
         oss.precision( 30 );
         oss << " Test: " << test_ << "\n"
             << " Error: Singular value decomposition failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A1 << "\n"
             << "   Row-major singular values:\n" << s1 << "\n"
             << "   Row-major right singular values:\n" << V1 << "\n"
             << "   Column-major decomposition:\n" << A2 << "\n"
             << "   Column-major singular values:\n" << s2 << "\n"
             << "   Column-major right singular values:\n" << V2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "gesdd( DenseMatrix, DenseVector, DenseMatrix, char ) (5x3, 'O')";

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> A1;
      randomize( A1 );
      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> A2( A1 );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s1;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s2;

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> V1;

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> V2;

      blaze::gesdd( A1, s1, V1, 'O' );
      blaze::gesdd( A2, s2, V2, 'O' );

      if( s1 != s2 || abs( A1 ) != abs( A2 ) || abs( V1 ) != abs( V2 ) ) {
         std::ostringstream oss;
         oss.precision( 30 );
         oss << " Test: " << test_ << "\n"
             << " Error: Singular value decomposition failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A1 << "\n"
             << "   Row-major singular values:\n" << s1 << "\n"
             << "   Row-major right singular values:\n" << V1 << "\n"
             << "   Column-major decomposition:\n" << A2 << "\n"
             << "   Column-major singular values:\n" << s2 << "\n"
             << "   Column-major right singular values:\n" << V2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // gesdd( DenseMatrix, DenseMatrix, DenseVector, DenseMatrix, char )
   //=====================================================================================

   {
      test_ = "gesdd( DenseMatrix, DenseMatrix, DenseVector, DenseMatrix, char ) (3x5, 'N')";

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> A1;
      randomize( A1 );
      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> A2( A1 );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s1;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s2;

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> U1;
      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> V1;

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> U2;
      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> V2;

      blaze::gesdd( A1, U1, s1, V1, 'N' );
      blaze::gesdd( A2, U2, s2, V2, 'N' );

      if( s1 != s2 || !isDefault( U1 ) || !isDefault( U2 ) || !isDefault( V1 ) || !isDefault( V2 ) ) {
         std::ostringstream oss;
         oss.precision( 30 );
         oss << " Test: " << test_ << "\n"
             << " Error: Singular value decomposition failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A1 << "\n"
             << "   Row-major singular values:\n" << s1 << "\n"
             << "   Row-major left singular values:\n" << U1 << "\n"
             << "   Row-major right singular values:\n" << V1 << "\n"
             << "   Column-major decomposition:\n" << A2 << "\n"
             << "   Column-major singular values:\n" << s2 << "\n"
             << "   Column-major left singular values:\n" << U2 << "\n"
             << "   Column-major right singular values:\n" << V2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "gesdd( DenseMatrix, DenseMatrix, DenseVector, DenseMatrix, char ) (5x3, 'N')";

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> A1;
      randomize( A1 );
      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> A2( A1 );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s1;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s2;

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> U1;
      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> V1;

      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> U2;
      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> V2;

      blaze::gesdd( A1, U1, s1, V1, 'N' );
      blaze::gesdd( A2, U2, s2, V2, 'N' );

      if( s1 != s2 || !isDefault( U1 ) || !isDefault( U2 ) || !isDefault( V1 ) || !isDefault( V2 ) ) {
         std::ostringstream oss;
         oss.precision( 30 );
         oss << " Test: " << test_ << "\n"
             << " Error: Singular value decomposition failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A1 << "\n"
             << "   Row-major singular values:\n" << s1 << "\n"
             << "   Row-major left singular values:\n" << U1 << "\n"
             << "   Row-major right singular values:\n" << V1 << "\n"
             << "   Column-major decomposition:\n" << A2 << "\n"
             << "   Column-major singular values:\n" << s2 << "\n"
             << "   Column-major left singular values:\n" << U2 << "\n"
             << "   Column-major right singular values:\n" << V2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "gesdd( DenseMatrix, DenseMatrix, DenseVector, DenseMatrix, char ) (3x5, 'S')";

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> A1;
      randomize( A1 );
      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> A2( A1 );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s1;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s2;

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> U1;
      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> V1;

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> U2;
      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> V2;

      blaze::gesdd( A1, U1, s1, V1, 'S' );
      blaze::gesdd( A2, U2, s2, V2, 'S' );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> S1;
      S1(0,0) = s1[0];
      S1(1,1) = s1[1];
      S1(2,2) = s1[2];

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> S2;
      S2(0,0) = s2[0];
      S2(1,1) = s2[1];
      S2(2,2) = s2[2];

      if( s1 != s2 || ( U1 * S1 * V1 ) != ( U2 * S2 * V2 ) ) {
         std::ostringstream oss;
         oss.precision( 30 );
         oss << " Test: " << test_ << "\n"
             << " Error: Singular value decomposition failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A1 << "\n"
             << "   Row-major singular values:\n" << s1 << "\n"
             << "   Row-major left singular values:\n" << U1 << "\n"
             << "   Row-major right singular values:\n" << V1 << "\n"
             << "   Column-major decomposition:\n" << A2 << "\n"
             << "   Column-major singular values:\n" << s2 << "\n"
             << "   Column-major left singular values:\n" << U2 << "\n"
             << "   Column-major right singular values:\n" << V2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "gesdd( DenseMatrix, DenseMatrix, DenseVector, DenseMatrix, char ) (5x3, 'S')";

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> A1;
      randomize( A1 );
      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> A2( A1 );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s1;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s2;

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> U1;
      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> V1;

      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> U2;
      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> V2;

      blaze::gesdd( A1, U1, s1, V1, 'S' );
      blaze::gesdd( A2, U2, s2, V2, 'S' );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> S1;
      S1(0,0) = s1[0];
      S1(1,1) = s1[1];
      S1(2,2) = s1[2];

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> S2;
      S2(0,0) = s2[0];
      S2(1,1) = s2[1];
      S2(2,2) = s2[2];

      if( s1 != s2 || ( U1 * S1 * V1 ) != ( U2 * S2 * V2 ) ) {
         std::ostringstream oss;
         oss.precision( 30 );
         oss << " Test: " << test_ << "\n"
             << " Error: Singular value decomposition failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A1 << "\n"
             << "   Row-major singular values:\n" << s1 << "\n"
             << "   Row-major left singular values:\n" << U1 << "\n"
             << "   Row-major right singular values:\n" << V1 << "\n"
             << "   Column-major decomposition:\n" << A2 << "\n"
             << "   Column-major singular values:\n" << s2 << "\n"
             << "   Column-major left singular values:\n" << U2 << "\n"
             << "   Column-major right singular values:\n" << V2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "gesdd( DenseMatrix, DenseMatrix, DenseVector, DenseMatrix, char ) (3x5, 'A')";

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> A1;
      randomize( A1 );
      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> A2( A1 );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s1;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s2;

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> U1;
      blaze::StaticMatrix<Type,5UL,5UL,blaze::rowMajor> V1;

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> U2;
      blaze::StaticMatrix<Type,5UL,5UL,blaze::columnMajor> V2;

      blaze::gesdd( A1, U1, s1, V1, 'A' );
      blaze::gesdd( A2, U2, s2, V2, 'A' );

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> S1;
      S1(0,0) = s1[0];
      S1(1,1) = s1[1];
      S1(2,2) = s1[2];

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> S2;
      S2(0,0) = s2[0];
      S2(1,1) = s2[1];
      S2(2,2) = s2[2];

      if( s1 != s2 || ( U1 * S1 * V1 ) != ( U2 * S2 * V2 ) ) {
         std::ostringstream oss;
         oss.precision( 30 );
         oss << " Test: " << test_ << "\n"
             << " Error: Singular value decomposition failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A1 << "\n"
             << "   Row-major singular values:\n" << s1 << "\n"
             << "   Row-major left singular values:\n" << U1 << "\n"
             << "   Row-major right singular values:\n" << V1 << "\n"
             << "   Column-major decomposition:\n" << A2 << "\n"
             << "   Column-major singular values:\n" << s2 << "\n"
             << "   Column-major left singular values:\n" << U2 << "\n"
             << "   Column-major right singular values:\n" << V2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "gesdd( DenseMatrix, DenseMatrix, DenseVector, DenseMatrix, char ) (5x3, 'A')";

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> A1;
      randomize( A1 );
      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> A2( A1 );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s1;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s2;

      blaze::StaticMatrix<Type,5UL,5UL,blaze::rowMajor> U1;
      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> V1;

      blaze::StaticMatrix<Type,5UL,5UL,blaze::columnMajor> U2;
      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> V2;

      blaze::gesdd( A1, U1, s1, V1, 'A' );
      blaze::gesdd( A2, U2, s2, V2, 'A' );

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> S1;
      S1(0,0) = s1[0];
      S1(1,1) = s1[1];
      S1(2,2) = s1[2];

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> S2;
      S2(0,0) = s2[0];
      S2(1,1) = s2[1];
      S2(2,2) = s2[2];

      if( s1 != s2 || ( U1 * S1 * V1 ) != ( U2 * S2 * V2 ) ) {
         std::ostringstream oss;
         oss.precision( 30 );
         oss << " Test: " << test_ << "\n"
             << " Error: Singular value decomposition failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A1 << "\n"
             << "   Row-major singular values:\n" << s1 << "\n"
             << "   Row-major left singular values:\n" << U1 << "\n"
             << "   Row-major right singular values:\n" << V1 << "\n"
             << "   Column-major decomposition:\n" << A2 << "\n"
             << "   Column-major singular values:\n" << s2 << "\n"
             << "   Column-major left singular values:\n" << U2 << "\n"
             << "   Column-major right singular values:\n" << V2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the singular value decomposition functions (gesvdx).
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the singular value decomposition functions for various data
// types. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void SingularValueTest::testGesvdx()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE && BLAZETEST_MATHTEST_LAPACK_SUPPORTS_GESVDX

   //=====================================================================================
   // gesvdx( DenseMatrix, DenseVector )
   //=====================================================================================

   {
      test_ = "gesvdx( DenseMatrix, DenseVector ) (3x5)";

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> A1;
      randomize( A1 );

      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> A2( A1 );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s1;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s2;

      blaze::gesvdx( A1, s1 );
      blaze::gesvdx( A2, s2 );

      if( s1 != s2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Singular value decomposition failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A1 << "\n"
             << "   Row-major singular values:\n" << s1 << "\n"
             << "   Column-major decomposition:\n" << A2 << "\n"
             << "   Column-major singular values:\n" << s2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "gesvdx( DenseMatrix, DenseVector ) (5x3)";

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> A1;
      randomize( A1 );

      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> A2( A1 );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s1;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s2;

      blaze::gesvdx( A1, s1 );
      blaze::gesvdx( A2, s2 );

      if( s1 != s2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Singular value decomposition failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A1 << "\n"
             << "   Row-major singular values:\n" << s1 << "\n"
             << "   Column-major decomposition:\n" << A2 << "\n"
             << "   Column-major singular values:\n" << s2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // gesvdx( DenseMatrix, DenseVector, double, double )
   //=====================================================================================

   {
      test_ = "gesvdx( DenseMatrix, DenseVector, double, double ) (3x5)";

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> A1;
      randomize( A1 );

      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> A2( A1 );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s1;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s2;

      blaze::gesvdx( A1, s1, 0.0, 5.0 );
      blaze::gesvdx( A2, s2, 0.0, 5.0 );

      if( s1 != s2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Singular value decomposition failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A1 << "\n"
             << "   Row-major singular values:\n" << s1 << "\n"
             << "   Column-major decomposition:\n" << A2 << "\n"
             << "   Column-major singular values:\n" << s2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "gesvdx( DenseMatrix, DenseVector, double, double ) (3x5)";

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> A1;
      randomize( A1 );

      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> A2( A1 );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s1;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s2;

      blaze::gesvdx( A1, s1, 0.0, 5.0 );
      blaze::gesvdx( A2, s2, 0.0, 5.0 );

      if( s1 != s2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Singular value decomposition failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A1 << "\n"
             << "   Row-major singular values:\n" << s1 << "\n"
             << "   Column-major decomposition:\n" << A2 << "\n"
             << "   Column-major singular values:\n" << s2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // gesvdx( DenseMatrix, DenseVector, int, int )
   //=====================================================================================

   {
      test_ = "gesvdx( DenseMatrix, DenseVector, int, int ) (3x5)";

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> A1;
      randomize( A1 );

      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> A2( A1 );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,2UL,blaze::columnVector> s1;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,2UL,blaze::columnVector> s2;

      blaze::gesvdx( A1, s1, 0, 1 );
      blaze::gesvdx( A2, s2, 0, 1 );

      if( s1 != s2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Singular value decomposition failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A1 << "\n"
             << "   Row-major singular values:\n" << s1 << "\n"
             << "   Column-major decomposition:\n" << A2 << "\n"
             << "   Column-major singular values:\n" << s2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "gesvdx( DenseMatrix, DenseVector, int, int ) (5x3)";

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> A1;
      randomize( A1 );

      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> A2( A1 );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,2UL,blaze::columnVector> s1;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,2UL,blaze::columnVector> s2;

      blaze::gesvdx( A1, s1, 0, 1 );
      blaze::gesvdx( A2, s2, 0, 1 );

      if( s1 != s2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Singular value decomposition failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A1 << "\n"
             << "   Row-major singular values:\n" << s1 << "\n"
             << "   Column-major decomposition:\n" << A2 << "\n"
             << "   Column-major singular values:\n" << s2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // gesvdx( DenseMatrix, DenseMatrix, DenseVector )
   //=====================================================================================

   {
      test_ = "gesvdx( DenseMatrix, DenseMatrix, DenseVector ) (3x5)";

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> A1;
      randomize( A1 );

      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> A2( A1 );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s1;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s2;

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> U1;

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> U2;

      blaze::gesvdx( A1, U1, s1 );
      blaze::gesvdx( A2, U2, s2 );

      if( s1 != s2 || abs( U1 ) != abs( U2 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Singular value decomposition failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A1 << "\n"
             << "   Row-major singular values:\n" << s1 << "\n"
             << "   Row-major left singular values:\n" << U1 << "\n"
             << "   Column-major decomposition:\n" << A2 << "\n"
             << "   Column-major singular values:\n" << s2 << "\n"
             << "   Column-major left singular values:\n" << U2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "gesvdx( DenseMatrix, DenseMatrix, DenseVector ) (5x3)";

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> A1;
      randomize( A1 );

      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> A2( A1 );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s1;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s2;

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> U1;

      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> U2;

      blaze::gesvdx( A1, U1, s1 );
      blaze::gesvdx( A2, U2, s2 );

      if( s1 != s2 || abs( U1 ) != abs( U2 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Singular value decomposition failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A1 << "\n"
             << "   Row-major singular values:\n" << s1 << "\n"
             << "   Row-major left singular values:\n" << U1 << "\n"
             << "   Column-major decomposition:\n" << A2 << "\n"
             << "   Column-major singular values:\n" << s2 << "\n"
             << "   Column-major left singular values:\n" << U2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // gesvdx( DenseMatrix, DenseMatrix, DenseVector, double, double )
   //=====================================================================================

   {
      test_ = "gesvdx( DenseMatrix, DenseMatrix, DenseVector, double, double ) (3x5)";

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> A1;
      randomize( A1 );

      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> A2( A1 );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s1;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s2;

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> U1;

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> U2;

      blaze::gesvdx( A1, U1, s1, 0.0, 5.0 );
      blaze::gesvdx( A2, U2, s2, 0.0, 5.0 );

      if( s1 != s2 || abs( U1 ) != abs( U2 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Singular value decomposition failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A1 << "\n"
             << "   Row-major singular values:\n" << s1 << "\n"
             << "   Row-major left singular values:\n" << U1 << "\n"
             << "   Column-major decomposition:\n" << A2 << "\n"
             << "   Column-major singular values:\n" << s2 << "\n"
             << "   Column-major left singular values:\n" << U2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "gesvdx( DenseMatrix, DenseMatrix, DenseVector, double, double ) (5x3)";

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> A1;
      randomize( A1 );

      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> A2( A1 );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s1;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s2;

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> U1;

      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> U2;

      blaze::gesvdx( A1, U1, s1, 0.0, 5.0 );
      blaze::gesvdx( A2, U2, s2, 0.0, 5.0 );

      if( s1 != s2 || abs( U1 ) != abs( U2 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Singular value decomposition failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A1 << "\n"
             << "   Row-major singular values:\n" << s1 << "\n"
             << "   Row-major left singular values:\n" << U1 << "\n"
             << "   Column-major decomposition:\n" << A2 << "\n"
             << "   Column-major singular values:\n" << s2 << "\n"
             << "   Column-major left singular values:\n" << U2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // gesvdx( DenseMatrix, DenseMatrix, DenseVector, int, int )
   //=====================================================================================

   {
      test_ = "gesvdx( DenseMatrix, DenseMatrix, DenseVector, int, int ) (3x5)";

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> A1;
      randomize( A1 );

      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> A2( A1 );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,2UL,blaze::columnVector> s1;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,2UL,blaze::columnVector> s2;

      blaze::StaticMatrix<Type,3UL,2UL,blaze::rowMajor> U1;

      blaze::StaticMatrix<Type,3UL,2UL,blaze::columnMajor> U2;

      blaze::gesvdx( A1, U1, s1, 0, 1 );
      blaze::gesvdx( A2, U2, s2, 0, 1 );

      if( s1 != s2 || abs( U1 ) != abs( U2 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Singular value decomposition failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A1 << "\n"
             << "   Row-major singular values:\n" << s1 << "\n"
             << "   Row-major left singular values:\n" << U1 << "\n"
             << "   Column-major decomposition:\n" << A2 << "\n"
             << "   Column-major singular values:\n" << s2 << "\n"
             << "   Column-major left singular values:\n" << U2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "gesvdx( DenseMatrix, DenseMatrix, DenseVector, int, int ) (5x3)";

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> A1;
      randomize( A1 );

      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> A2( A1 );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,2UL,blaze::columnVector> s1;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,2UL,blaze::columnVector> s2;

      blaze::StaticMatrix<Type,5UL,2UL,blaze::rowMajor> U1;

      blaze::StaticMatrix<Type,5UL,2UL,blaze::columnMajor> U2;

      blaze::gesvdx( A1, U1, s1, 0, 1 );
      blaze::gesvdx( A2, U2, s2, 0, 1 );

      if( s1 != s2 || abs( U1 ) != abs( U2 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Singular value decomposition failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A1 << "\n"
             << "   Row-major singular values:\n" << s1 << "\n"
             << "   Row-major left singular values:\n" << U1 << "\n"
             << "   Column-major decomposition:\n" << A2 << "\n"
             << "   Column-major singular values:\n" << s2 << "\n"
             << "   Column-major left singular values:\n" << U2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // gesvdx( DenseMatrix, DenseVector, DenseMatrix )
   //=====================================================================================

   {
      test_ = "gesvdx( DenseMatrix, DenseVector, DenseMatrix ) (3x5)";

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> A1;
      randomize( A1 );

      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> A2( A1 );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s1;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s2;

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> V1;

      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> V2;

      blaze::gesvdx( A1, s1, V1 );
      blaze::gesvdx( A2, s2, V2 );

      if( s1 != s2 || abs( V1 ) != abs( V2 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Singular value decomposition failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A1 << "\n"
             << "   Row-major singular values:\n" << s1 << "\n"
             << "   Row-major right singular values:\n" << V1 << "\n"
             << "   Column-major decomposition:\n" << A2 << "\n"
             << "   Column-major singular values:\n" << s2 << "\n"
             << "   Column-major right singular values:\n" << V2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "gesvdx( DenseMatrix, DenseMatrix, DenseVector, DenseMatrix ) (5x3)";

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> A1;
      randomize( A1 );

      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> A2( A1 );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s1;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s2;

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> V1;

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> V2;

      blaze::gesvdx( A1, s1, V1 );
      blaze::gesvdx( A2, s2, V2 );

      if( s1 != s2 || abs( V1 ) != abs( V2 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Singular value decomposition failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A1 << "\n"
             << "   Row-major singular values:\n" << s1 << "\n"
             << "   Row-major right singular values:\n" << V1 << "\n"
             << "   Column-major decomposition:\n" << A2 << "\n"
             << "   Column-major singular values:\n" << s2 << "\n"
             << "   Column-major right singular values:\n" << V2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // gesvdx( DenseMatrix, DenseVector, DenseMatrix, double, double )
   //=====================================================================================

   {
      test_ = "gesvdx( DenseMatrix, DenseVector, DenseMatrix, double, double ) (3x5)";

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> A1;
      randomize( A1 );

      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> A2( A1 );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s1;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s2;

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> V1;

      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> V2;

      blaze::gesvdx( A1, s1, V1, 0.0, 5.0 );
      blaze::gesvdx( A2, s2, V2, 0.0, 5.0 );

      if( s1 != s2 || abs( V1 ) != abs( V2 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Singular value decomposition failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A1 << "\n"
             << "   Row-major singular values:\n" << s1 << "\n"
             << "   Row-major right singular values:\n" << V1 << "\n"
             << "   Column-major decomposition:\n" << A2 << "\n"
             << "   Column-major singular values:\n" << s2 << "\n"
             << "   Column-major right singular values:\n" << V2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "gesvdx( DenseMatrix, DenseVector, DenseMatrix, double, double ) (5x3)";

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> A1;
      randomize( A1 );

      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> A2( A1 );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s1;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s2;

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> V1;

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> V2;

      blaze::gesvdx( A1, s1, V1, 0.0, 5.0 );
      blaze::gesvdx( A2, s2, V2, 0.0, 5.0 );

      if( s1 != s2 || abs( V1 ) != abs( V2 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Singular value decomposition failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A1 << "\n"
             << "   Row-major singular values:\n" << s1 << "\n"
             << "   Row-major right singular values:\n" << V1 << "\n"
             << "   Column-major decomposition:\n" << A2 << "\n"
             << "   Column-major singular values:\n" << s2 << "\n"
             << "   Column-major right singular values:\n" << V2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // gesvdx( DenseMatrix, DenseVector, DenseMatrix, int, int )
   //=====================================================================================

   {
      test_ = "gesvdx( DenseMatrix, DenseVector, DenseMatrix, int, int ) (3x5)";

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> A1;
      randomize( A1 );

      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> A2( A1 );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,2UL,blaze::columnVector> s1;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,2UL,blaze::columnVector> s2;

      blaze::StaticMatrix<Type,2UL,5UL,blaze::rowMajor> V1;

      blaze::StaticMatrix<Type,2UL,5UL,blaze::columnMajor> V2;

      blaze::gesvdx( A1, s1, V1, 0, 1 );
      blaze::gesvdx( A2, s2, V2, 0, 1 );

      if( s1 != s2 || abs( V1 ) != abs( V2 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Singular value decomposition failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A1 << "\n"
             << "   Row-major singular values:\n" << s1 << "\n"
             << "   Row-major right singular values:\n" << V1 << "\n"
             << "   Column-major decomposition:\n" << A2 << "\n"
             << "   Column-major singular values:\n" << s2 << "\n"
             << "   Column-major right singular values:\n" << V2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "gesvdx( DenseMatrix, DenseVector, DenseMatrix, int, int ) (5x3)";

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> A1;
      randomize( A1 );

      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> A2( A1 );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,2UL,blaze::columnVector> s1;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,2UL,blaze::columnVector> s2;

      blaze::StaticMatrix<Type,2UL,3UL,blaze::rowMajor> V1;

      blaze::StaticMatrix<Type,2UL,3UL,blaze::columnMajor> V2;

      blaze::gesvdx( A1, s1, V1, 0, 1 );
      blaze::gesvdx( A2, s2, V2, 0, 1 );

      if( s1 != s2 || abs( V1 ) != abs( V2 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Singular value decomposition failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A1 << "\n"
             << "   Row-major singular values:\n" << s1 << "\n"
             << "   Row-major right singular values:\n" << V1 << "\n"
             << "   Column-major decomposition:\n" << A2 << "\n"
             << "   Column-major singular values:\n" << s2 << "\n"
             << "   Column-major right singular values:\n" << V2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // gesvdx( DenseMatrix, DenseMatrix, DenseVector, DenseMatrix )
   //=====================================================================================

   {
      test_ = "gesvdx( DenseMatrix, DenseMatrix, DenseVector, DenseMatrix ) (3x5)";

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> A1;
      randomize( A1 );

      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> A2( A1 );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s1;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s2;

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> U1;
      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> V1;

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> U2;
      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> V2;

      blaze::gesvdx( A1, U1, s1, V1 );
      blaze::gesvdx( A2, U2, s2, V2 );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> S1;
      S1(0,0) = s1[0];
      S1(1,1) = s1[1];
      S1(2,2) = s1[2];

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> S2;
      S2(0,0) = s2[0];
      S2(1,1) = s2[1];
      S2(2,2) = s2[2];

      if( s1 != s2 || ( U1 * S1 * V1 ) != ( U2 * S2 * V2 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Singular value decomposition failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A1 << "\n"
             << "   Row-major singular values:\n" << s1 << "\n"
             << "   Row-major left singular values:\n" << U1 << "\n"
             << "   Row-major right singular values:\n" << V1 << "\n"
             << "   Column-major decomposition:\n" << A2 << "\n"
             << "   Column-major singular values:\n" << s2 << "\n"
             << "   Column-major left singular values:\n" << U2 << "\n"
             << "   Column-major right singular values:\n" << V2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "gesvdx( DenseMatrix, DenseMatrix, DenseVector, DenseMatrix ) (5x3)";

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> A1;
      randomize( A1 );

      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> A2( A1 );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s1;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s2;

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> U1;
      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> V1;

      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> U2;
      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> V2;

      blaze::gesvdx( A1, U1, s1, V1 );
      blaze::gesvdx( A2, U2, s2, V2 );

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> S1;
      S1(0,0) = s1[0];
      S1(1,1) = s1[1];
      S1(2,2) = s1[2];

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> S2;
      S2(0,0) = s2[0];
      S2(1,1) = s2[1];
      S2(2,2) = s2[2];

      if( s1 != s2 || ( U1 * S1 * V1 ) != ( U2 * S2 * V2 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Singular value decomposition failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A1 << "\n"
             << "   Row-major singular values:\n" << s1 << "\n"
             << "   Row-major left singular values:\n" << U1 << "\n"
             << "   Row-major right singular values:\n" << V1 << "\n"
             << "   Column-major decomposition:\n" << A2 << "\n"
             << "   Column-major singular values:\n" << s2 << "\n"
             << "   Column-major left singular values:\n" << U2 << "\n"
             << "   Column-major right singular values:\n" << V2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // gesvdx( DenseMatrix, DenseMatrix, DenseVector, DenseMatrix, double, double )
   //=====================================================================================

   {
      test_ = "gesvdx( DenseMatrix, DenseMatrix, DenseVector, DenseMatrix, double, double ) (3x5)";

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> A1;
      randomize( A1 );

      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> A2( A1 );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s1;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s2;

      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> U1;
      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> V1;

      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> U2;
      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> V2;

      blaze::gesvdx( A1, U1, s1, V1, 0.0, 5.0 );
      blaze::gesvdx( A2, U2, s2, V2, 0.0, 5.0 );

      if( s1 != s2 || abs( U1 ) != abs( U2 ) || abs( V1 ) != abs( V2 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Singular value decomposition failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A1 << "\n"
             << "   Row-major singular values:\n" << s1 << "\n"
             << "   Row-major left singular values:\n" << U1 << "\n"
             << "   Row-major right singular values:\n" << V1 << "\n"
             << "   Column-major decomposition:\n" << A2 << "\n"
             << "   Column-major singular values:\n" << s2 << "\n"
             << "   Column-major left singular values:\n" << U2 << "\n"
             << "   Column-major right singular values:\n" << V2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "gesvdx( DenseMatrix, DenseMatrix, DenseVector, DenseMatrix, double, double ) (5x3)";

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> A1;
      randomize( A1 );

      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> A2( A1 );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s1;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,3UL,blaze::columnVector> s2;

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> U1;
      blaze::StaticMatrix<Type,3UL,3UL,blaze::rowMajor> V1;

      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> U2;
      blaze::StaticMatrix<Type,3UL,3UL,blaze::columnMajor> V2;

      blaze::gesvdx( A1, U1, s1, V1, 0.0, 5.0 );
      blaze::gesvdx( A2, U2, s2, V2, 0.0, 5.0 );

      if( s1 != s2 || abs( U1 ) != abs( U2 ) || abs( V1 ) != abs( V2 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Singular value decomposition failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A1 << "\n"
             << "   Row-major singular values:\n" << s1 << "\n"
             << "   Row-major left singular values:\n" << U1 << "\n"
             << "   Row-major right singular values:\n" << V1 << "\n"
             << "   Column-major decomposition:\n" << A2 << "\n"
             << "   Column-major singular values:\n" << s2 << "\n"
             << "   Column-major left singular values:\n" << U2 << "\n"
             << "   Column-major right singular values:\n" << V2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // gesvdx( DenseMatrix, DenseMatrix, DenseVector, DenseMatrix, int, int )
   //=====================================================================================

   {
      test_ = "gesvdx( DenseMatrix, DenseMatrix, DenseVector, DenseMatrix, int, int ) (3x5)";

      blaze::StaticMatrix<Type,3UL,5UL,blaze::rowMajor> A1;
      randomize( A1 );

      blaze::StaticMatrix<Type,3UL,5UL,blaze::columnMajor> A2( A1 );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,2UL,blaze::columnVector> s1;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,2UL,blaze::columnVector> s2;

      blaze::StaticMatrix<Type,3UL,2UL,blaze::rowMajor> U1;
      blaze::StaticMatrix<Type,2UL,5UL,blaze::rowMajor> V1;

      blaze::StaticMatrix<Type,3UL,2UL,blaze::columnMajor> U2;
      blaze::StaticMatrix<Type,2UL,5UL,blaze::columnMajor> V2;

      blaze::gesvdx( A1, U1, s1, V1, 0, 1 );
      blaze::gesvdx( A2, U2, s2, V2, 0, 1 );

      if( s1 != s2 || abs( U1 ) != abs( U2 ) || abs( V1 ) != abs( V2 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Singular value decomposition failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A1 << "\n"
             << "   Row-major singular values:\n" << s1 << "\n"
             << "   Row-major left singular values:\n" << U1 << "\n"
             << "   Row-major right singular values:\n" << V1 << "\n"
             << "   Column-major decomposition:\n" << A2 << "\n"
             << "   Column-major singular values:\n" << s2 << "\n"
             << "   Column-major left singular values:\n" << U2 << "\n"
             << "   Column-major right singular values:\n" << V2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "gesvdx( DenseMatrix, DenseMatrix, DenseVector, DenseMatrix, int, int ) (5x3)";

      blaze::StaticMatrix<Type,5UL,3UL,blaze::rowMajor> A1;
      randomize( A1 );

      blaze::StaticMatrix<Type,5UL,3UL,blaze::columnMajor> A2( A1 );

      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,2UL,blaze::columnVector> s1;
      blaze::StaticVector<blaze::UnderlyingElement_t<Type>,2UL,blaze::columnVector> s2;

      blaze::StaticMatrix<Type,5UL,2UL,blaze::rowMajor> U1;
      blaze::StaticMatrix<Type,2UL,3UL,blaze::rowMajor> V1;

      blaze::StaticMatrix<Type,5UL,2UL,blaze::columnMajor> U2;
      blaze::StaticMatrix<Type,2UL,3UL,blaze::columnMajor> V2;

      blaze::gesvdx( A1, U1, s1, V1, 0, 1 );
      blaze::gesvdx( A2, U2, s2, V2, 0, 1 );

      if( s1 != s2 || abs( U1 ) != abs( U2 ) || abs( V1 ) != abs( V2 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Singular value decomposition failed\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Element type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Row-major decomposition:\n" << A1 << "\n"
             << "   Row-major singular values:\n" << s1 << "\n"
             << "   Row-major left singular values:\n" << U1 << "\n"
             << "   Row-major right singular values:\n" << V1 << "\n"
             << "   Column-major decomposition:\n" << A2 << "\n"
             << "   Column-major singular values:\n" << s2 << "\n"
             << "   Column-major left singular values:\n" << U2 << "\n"
             << "   Column-major right singular values:\n" << V2 << "\n";
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
   SingularValueTest();
}
//*************************************************************************************************




//=================================================================================================
//
//  MACRO DEFINITIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Macro for the execution of the LAPACK singular value test.
*/
#define RUN_LAPACK_SINGULAR_VALUE_TEST \
   blazetest::mathtest::lapack::runTest()
/*! \endcond */
//*************************************************************************************************

} // namespace lapack

} // namespace mathtest

} // namespace blazetest

#endif
