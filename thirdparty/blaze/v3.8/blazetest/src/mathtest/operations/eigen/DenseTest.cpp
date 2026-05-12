//=================================================================================================
/*!
//  \file src/mathtest/operations/eigen/DenseTest.cpp
//  \brief Source file for the dense matrix eigenvalue test
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

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <blaze/math/Column.h>
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/DynamicVector.h>
#include <blaze/math/DiagonalMatrix.h>
#include <blaze/math/HermitianMatrix.h>
#include <blaze/math/LowerMatrix.h>
#include <blaze/math/UpperMatrix.h>
#include <blaze/math/Row.h>
#include <blaze/math/shims/Equal.h>
#include <blaze/math/SymmetricMatrix.h>
#include <blaze/util/Complex.h>
#include <blazetest/mathtest/operations/eigen/DenseTest.h>
#include <blazetest/system/LAPACK.h>

#ifdef BLAZE_USE_HPX_THREADS
#  include <hpx/hpx_main.hpp>
#endif


namespace blazetest {

namespace mathtest {

namespace operations {

namespace eigen {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the DenseTest test.
//
// \exception std::runtime_error Error during eigenvalue/eigenvector computation detected.
*/
DenseTest::DenseTest()
{
   testGeneral();
   testSymmetric();
   testHermitian();
   testLower();
   testUpper();
   testDiagonal();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the eigenvalue/eigenvector evaluation for general matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function tests the dense matrix eigenvalue/eigenvector evaluation for general matrices.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testGeneral()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   using blaze::DynamicMatrix;
   using blaze::DynamicVector;
   using blaze::complex;
   using blaze::eigen;
   using blaze::rowMajor;
   using blaze::columnMajor;
   using blaze::rowVector;

   const auto comparator = []( const auto& v1, const auto& v2 ) {
      return blaze::equal( v1, v2 );
   };


   //=====================================================================================
   // eigen( DenseMatrix, DenseVector )
   //=====================================================================================

   {
      test_ = "eigen( DenseMatrix, DenseVector ) (general, double)";

      DynamicMatrix<double,rowMajor> A1( 5UL, 5UL );
      randomize( A1 );
      DynamicMatrix<double,columnMajor> A2( A1 );

      {
         DynamicVector<complex<double>,rowVector> w1;
         DynamicVector<complex<double>,rowVector> w2;

         eigen( A1, w1 );
         eigen( A2, w2 );

         if( !std::is_permutation( w1.begin(), w1.end(), w2.begin(), comparator ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Eigenvalue computation failed\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Row-major eigenvalues:\n" << w1 << "\n"
                << "   Column-major eigenvalues:\n" << w2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         DynamicVector<complex<double>,rowVector> w1;
         DynamicVector<complex<double>,rowVector> w2;

         eigen( trans( A1 ), w1 );
         eigen( trans( A2 ), w2 );

         if( !std::is_permutation( w1.begin(), w1.end(), w2.begin(), comparator ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Transpose eigenvalue computation failed\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Row-major eigenvalues:\n" << w1 << "\n"
                << "   Column-major eigenvalues:\n" << w2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   {
      test_ = "eigen( DenseMatrix, DenseVector ) (general, complex<double>)";

      DynamicMatrix<complex<double>,rowMajor> A1( 5UL, 5UL );
      randomize( A1 );
      DynamicMatrix<complex<double>,columnMajor> A2( A1 );

      {
         DynamicVector<complex<double>,rowVector> w1;
         DynamicVector<complex<double>,rowVector> w2;

         eigen( A1, w1 );
         eigen( A2, w2 );

         if( !std::is_permutation( w1.begin(), w1.end(), w2.begin(), comparator ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Eigenvalue computation failed\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Row-major eigenvalues:\n" << w1 << "\n"
                << "   Column-major eigenvalues:\n" << w2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         DynamicVector<complex<double>,rowVector> w1;
         DynamicVector<complex<double>,rowVector> w2;

         eigen( trans( A1 ), w1 );
         eigen( trans( A2 ), w2 );

         if( !std::is_permutation( w1.begin(), w1.end(), w2.begin(), comparator ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Transpose eigenvalue computation failed\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Row-major eigenvalues:\n" << w1 << "\n"
                << "   Column-major eigenvalues:\n" << w2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // eigen( DenseMatrix, DenseVector, DenseMatrix )
   //=====================================================================================

   {
      test_ = "eigen( DenseMatrix, DenseVector, DenseMatrix ) (general, double)";

      DynamicMatrix<double,rowMajor> A( 5UL, 5UL );
      randomize( A );

      DynamicMatrix<double,rowMajor>    A1( A );
      DynamicMatrix<double,columnMajor> A2( A );

      {
         DynamicVector<complex<double>,rowVector> w1;
         DynamicVector<complex<double>,rowVector> w2;

         DynamicMatrix<complex<double>,rowMajor>    V1;
         DynamicMatrix<complex<double>,columnMajor> V2;

         eigen( A1, w1, V1 );
         eigen( A2, w2, V2 );

         if( !std::is_permutation( w1.begin(), w1.end(), w2.begin(), comparator ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Eigenvalue computation failed\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Row-major eigenvalues:\n" << w1 << "\n"
                << "   Column-major eigenvalues:\n" << w2 << "\n";
            throw std::runtime_error( oss.str() );
         }

         for( size_t i=0UL; i<V1.rows(); ++i ) {
            checkEigenvector( row( V1, i ), A, w1[i] );
         }

         for( size_t i=0UL; i<V2.columns(); ++i ) {
            checkEigenvector( column( V2, i ), A, w2[i] );
         }
      }

      {
         DynamicVector<complex<double>,rowVector> w1;
         DynamicVector<complex<double>,rowVector> w2;

         DynamicMatrix<complex<double>,rowMajor>    V1;
         DynamicMatrix<complex<double>,columnMajor> V2;

         eigen( trans( A1 ), w1, V1 );
         eigen( trans( A2 ), w2, V2 );

         if( !std::is_permutation( w1.begin(), w1.end(), w2.begin(), comparator ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Transpose eigenvalue computation failed\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Row-major eigenvalues:\n" << w1 << "\n"
                << "   Column-major eigenvalues:\n" << w2 << "\n";
            throw std::runtime_error( oss.str() );
         }

         for( size_t i=0UL; i<V1.rows(); ++i ) {
            checkEigenvector( ctrans( row( V1, i ) ), trans( A ), w1[i] );
         }

         for( size_t i=0UL; i<V2.columns(); ++i ) {
            checkEigenvector( ctrans( column( V2, i ) ), trans( A ), w2[i] );
         }
      }
   }

   {
      test_ = "eigen( DenseMatrix, DenseVector, DenseMatrix ) (general, complex<double>)";

      DynamicMatrix<complex<double>,rowMajor> A( 5UL, 5UL );
      randomize( A );

      DynamicMatrix<complex<double>,rowMajor>    A1( A );
      DynamicMatrix<complex<double>,columnMajor> A2( A );

      {
         DynamicVector<complex<double>,rowVector> w1;
         DynamicVector<complex<double>,rowVector> w2;

         DynamicMatrix<complex<double>,rowMajor>    V1;
         DynamicMatrix<complex<double>,columnMajor> V2;

         eigen( A1, w1, V1 );
         eigen( A2, w2, V2 );

         if( !std::is_permutation( w1.begin(), w1.end(), w2.begin(), comparator ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Eigenvalue computation failed\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Row-major eigenvalues:\n" << w1 << "\n"
                << "   Column-major eigenvalues:\n" << w2 << "\n";
            throw std::runtime_error( oss.str() );
         }

         for( size_t i=0UL; i<V1.rows(); ++i ) {
            checkEigenvector( row( V1, i ), A, w1[i] );
         }

         for( size_t i=0UL; i<V2.columns(); ++i ) {
            checkEigenvector( column( V2, i ), A, w2[i] );
         }
      }

      {
         DynamicVector<complex<double>,rowVector> w1;
         DynamicVector<complex<double>,rowVector> w2;

         DynamicMatrix<complex<double>,rowMajor>    V1;
         DynamicMatrix<complex<double>,columnMajor> V2;

         eigen( trans( A1 ), w1, V1 );
         eigen( trans( A2 ), w2, V2 );

         if( !std::is_permutation( w1.begin(), w1.end(), w2.begin(), comparator ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Transpose eigenvalue computation failed\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Row-major eigenvalues:\n" << w1 << "\n"
                << "   Column-major eigenvalues:\n" << w2 << "\n";
            throw std::runtime_error( oss.str() );
         }

         for( size_t i=0UL; i<V1.rows(); ++i ) {
            checkEigenvector( trans( row( V1, i ) ), trans( A ), w1[i] );
         }

         for( size_t i=0UL; i<V2.columns(); ++i ) {
            checkEigenvector( trans( column( V2, i ) ), trans( A ), w2[i] );
         }
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the eigenvalue/eigenvector evaluation for symmetric matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function tests the dense matrix eigenvalue/eigenvector evaluation for symmetric matrices.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testSymmetric()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   using blaze::SymmetricMatrix;
   using blaze::DynamicMatrix;
   using blaze::DynamicVector;
   using blaze::complex;
   using blaze::eigen;
   using blaze::rowMajor;
   using blaze::columnMajor;
   using blaze::rowVector;

   const auto comparator = []( const auto& v1, const auto& v2 ) {
      return blaze::equal( v1, v2 );
   };


   //=====================================================================================
   // eigen( DenseMatrix, DenseVector )
   //=====================================================================================

   {
      test_ = "eigen( DenseMatrix, DenseVector ) (symmetric, double)";

      SymmetricMatrix< DynamicMatrix<double,rowMajor> > A1( 5UL );
      randomize( A1 );
      SymmetricMatrix< DynamicMatrix<double,columnMajor> > A2( A1 );

      {
         DynamicVector<double,rowVector> w1;
         DynamicVector<double,rowVector> w2;

         eigen( A1, w1 );
         eigen( A2, w2 );

         if( w1 != w2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Eigenvalue computation failed\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Row-major eigenvalues:\n" << w1 << "\n"
                << "   Column-major eigenvalues:\n" << w2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         DynamicVector<double,rowVector> w1;
         DynamicVector<double,rowVector> w2;

         eigen( trans( A1 ), w1 );
         eigen( trans( A2 ), w2 );

         if( w1 != w2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Transpose eigenvalue computation failed\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Row-major eigenvalues:\n" << w1 << "\n"
                << "   Column-major eigenvalues:\n" << w2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   {
      test_ = "eigen( DenseMatrix, DenseVector ) (symmetric, complex<double>)";

      SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> > A1( 5UL );
      randomize( A1 );
      SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> > A2( A1 );

      {
         DynamicVector<complex<double>,rowVector> w1;
         DynamicVector<complex<double>,rowVector> w2;

         eigen( A1, w1 );
         eigen( A2, w2 );

         if( w1 != w2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Eigenvalue computation failed\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Row-major eigenvalues:\n" << w1 << "\n"
                << "   Column-major eigenvalues:\n" << w2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         DynamicVector<complex<double>,rowVector> w1;
         DynamicVector<complex<double>,rowVector> w2;

         eigen( trans( A1 ), w1 );
         eigen( trans( A2 ), w2 );

         if( w1 != w2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Transpose eigenvalue computation failed\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Row-major eigenvalues:\n" << w1 << "\n"
                << "   Column-major eigenvalues:\n" << w2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // eigen( DenseMatrix, DenseVector, DenseMatrix )
   //=====================================================================================

   {
      test_ = "eigen( DenseMatrix, DenseVector, DenseMatrix ) (symmetric, double)";

      SymmetricMatrix< DynamicMatrix<double,rowMajor> > A( 5UL );
      randomize( A );

      SymmetricMatrix< DynamicMatrix<double,rowMajor> >    A1( A );
      SymmetricMatrix< DynamicMatrix<double,columnMajor> > A2( A );

      {
         DynamicVector<double,rowVector> w1;
         DynamicVector<double,rowVector> w2;

         DynamicMatrix<double,rowMajor>    V1;
         DynamicMatrix<double,columnMajor> V2;

         eigen( A1, w1, V1 );
         eigen( A2, w2, V2 );

         if( !std::is_permutation( w1.begin(), w1.end(), w2.begin(), comparator ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Eigenvalue computation failed\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Row-major eigenvalues:\n" << w1 << "\n"
                << "   Column-major eigenvalues:\n" << w2 << "\n";
            throw std::runtime_error( oss.str() );
         }

         for( size_t i=0UL; i<V1.rows(); ++i ) {
            checkEigenvector( row( V1, i ), A, w1[i] );
         }

         for( size_t i=0UL; i<V2.columns(); ++i ) {
            checkEigenvector( column( V2, i ), A, w2[i] );
         }
      }

      {
         DynamicVector<double,rowVector> w1;
         DynamicVector<double,rowVector> w2;

         DynamicMatrix<double,rowMajor>    V1;
         DynamicMatrix<double,columnMajor> V2;

         eigen( trans( A1 ), w1, V1 );
         eigen( trans( A2 ), w2, V2 );

         if( !std::is_permutation( w1.begin(), w1.end(), w2.begin(), comparator ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Transpose eigenvalue computation failed\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Row-major eigenvalues:\n" << w1 << "\n"
                << "   Column-major eigenvalues:\n" << w2 << "\n";
            throw std::runtime_error( oss.str() );
         }

         for( size_t i=0UL; i<V1.rows(); ++i ) {
            checkEigenvector( trans( row( V1, i ) ), trans( A ), w1[i] );
         }

         for( size_t i=0UL; i<V2.columns(); ++i ) {
            checkEigenvector( trans( column( V2, i ) ), trans( A ), w2[i] );
         }
      }
   }

   {
      test_ = "eigen( DenseMatrix, DenseVector, DenseMatrix ) (symmetric, complex<double>)";

      SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> > A( 5UL );
      randomize( A );

      SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >    A1( A );
      SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> > A2( A );

      {
         DynamicVector<complex<double>,rowVector> w1;
         DynamicVector<complex<double>,rowVector> w2;

         DynamicMatrix<complex<double>,rowMajor>    V1;
         DynamicMatrix<complex<double>,columnMajor> V2;

         eigen( A1, w1, V1 );
         eigen( A2, w2, V2 );

         if( !std::is_permutation( w1.begin(), w1.end(), w2.begin(), comparator ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Eigenvalue computation failed\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Row-major eigenvalues:\n" << w1 << "\n"
                << "   Column-major eigenvalues:\n" << w2 << "\n";
            throw std::runtime_error( oss.str() );
         }

         for( size_t i=0UL; i<V1.rows(); ++i ) {
            checkEigenvector( row( V1, i ), A, w1[i] );
         }

         for( size_t i=0UL; i<V2.columns(); ++i ) {
            checkEigenvector( column( V2, i ), A, w2[i] );
         }
      }

      {
         DynamicVector<complex<double>,rowVector> w1;
         DynamicVector<complex<double>,rowVector> w2;

         DynamicMatrix<complex<double>,rowMajor>    V1;
         DynamicMatrix<complex<double>,columnMajor> V2;

         eigen( trans( A1 ), w1, V1 );
         eigen( trans( A2 ), w2, V2 );

         if( !std::is_permutation( w1.begin(), w1.end(), w2.begin(), comparator ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Transpose eigenvalue computation failed\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Row-major eigenvalues:\n" << w1 << "\n"
                << "   Column-major eigenvalues:\n" << w2 << "\n";
            throw std::runtime_error( oss.str() );
         }

         for( size_t i=0UL; i<V1.rows(); ++i ) {
            checkEigenvector( trans( row( V1, i ) ), trans( A ), w1[i] );
         }

         for( size_t i=0UL; i<V2.columns(); ++i ) {
            checkEigenvector( trans( column( V2, i ) ), trans( A ), w2[i] );
         }
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the eigenvalue/eigenvector evaluation for Hermitian matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function tests the dense matrix eigenvalue/eigenvector evaluation for Hermitian matrices.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testHermitian()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   using blaze::HermitianMatrix;
   using blaze::DynamicMatrix;
   using blaze::DynamicVector;
   using blaze::complex;
   using blaze::eigen;
   using blaze::rowMajor;
   using blaze::columnMajor;
   using blaze::rowVector;

   const auto comparator = []( const auto& v1, const auto& v2 ) {
      return blaze::equal( v1, v2 );
   };


   //=====================================================================================
   // eigen( DenseMatrix, DenseVector )
   //=====================================================================================

   {
      test_ = "eigen( DenseMatrix, DenseVector ) (Hermitian, double)";

      HermitianMatrix< DynamicMatrix<double,rowMajor> > A1( 5UL );
      randomize( A1 );
      HermitianMatrix< DynamicMatrix<double,columnMajor> > A2( A1 );

      {
         DynamicVector<double,rowVector> w1;
         DynamicVector<double,rowVector> w2;

         eigen( A1, w1 );
         eigen( A2, w2 );

         if( w1 != w2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Eigenvalue computation failed\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Row-major eigenvalues:\n" << w1 << "\n"
                << "   Column-major eigenvalues:\n" << w2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         DynamicVector<double,rowVector> w1;
         DynamicVector<double,rowVector> w2;

         eigen( trans( A1 ), w1 );
         eigen( trans( A2 ), w2 );

         if( w1 != w2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Transpose eigenvalue computation failed\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Row-major eigenvalues:\n" << w1 << "\n"
                << "   Column-major eigenvalues:\n" << w2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   {
      test_ = "eigen( DenseMatrix, DenseVector ) (Hermitian, complex<double>)";

      HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> > A1( 5UL );
      randomize( A1 );
      HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> > A2( A1 );

      {
         DynamicVector<double,rowVector> w1;
         DynamicVector<double,rowVector> w2;

         eigen( A1, w1 );
         eigen( A2, w2 );

         if( w1 != w2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Eigenvalue computation failed\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Row-major eigenvalues:\n" << w1 << "\n"
                << "   Column-major eigenvalues:\n" << w2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         DynamicVector<double,rowVector> w1;
         DynamicVector<double,rowVector> w2;

         eigen( trans( A1 ), w1 );
         eigen( trans( A2 ), w2 );

         if( w1 != w2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Transpose eigenvalue computation failed\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Row-major eigenvalues:\n" << w1 << "\n"
                << "   Column-major eigenvalues:\n" << w2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // eigen( DenseMatrix, DenseVector, DenseMatrix )
   //=====================================================================================

   {
      test_ = "eigen( DenseMatrix, DenseVector, DenseMatrix ) (Hermitian, double)";

      HermitianMatrix< DynamicMatrix<double,rowMajor> > A( 5UL );
      randomize( A );

      HermitianMatrix< DynamicMatrix<double,rowMajor> >    A1( A );
      HermitianMatrix< DynamicMatrix<double,columnMajor> > A2( A );

      {
         DynamicVector<double,rowVector> w1;
         DynamicVector<double,rowVector> w2;

         DynamicMatrix<double,rowMajor>    V1;
         DynamicMatrix<double,columnMajor> V2;

         eigen( A1, w1, V1 );
         eigen( A2, w2, V2 );

         if( !std::is_permutation( w1.begin(), w1.end(), w2.begin(), comparator ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Eigenvalue computation failed\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Row-major eigenvalues:\n" << w1 << "\n"
                << "   Column-major eigenvalues:\n" << w2 << "\n";
            throw std::runtime_error( oss.str() );
         }

         for( size_t i=0UL; i<V1.rows(); ++i ) {
            checkEigenvector( row( V1, i ), A, w1[i] );
         }

         for( size_t i=0UL; i<V2.columns(); ++i ) {
            checkEigenvector( column( V2, i ), A, w2[i] );
         }
      }

      {
         DynamicVector<double,rowVector> w1;
         DynamicVector<double,rowVector> w2;

         DynamicMatrix<double,rowMajor>    V1;
         DynamicMatrix<double,columnMajor> V2;

         eigen( trans( A1 ), w1, V1 );
         eigen( trans( A2 ), w2, V2 );

         if( !std::is_permutation( w1.begin(), w1.end(), w2.begin(), comparator ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Transpose eigenvalue computation failed\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Row-major eigenvalues:\n" << w1 << "\n"
                << "   Column-major eigenvalues:\n" << w2 << "\n";
            throw std::runtime_error( oss.str() );
         }

         for( size_t i=0UL; i<V1.rows(); ++i ) {
            checkEigenvector( trans( row( V1, i ) ), trans( A ), w1[i] );
         }

         for( size_t i=0UL; i<V2.columns(); ++i ) {
            checkEigenvector( trans( column( V2, i ) ), trans( A ), w2[i] );
         }
      }
   }

   {
      test_ = "eigen( DenseMatrix, DenseVector, DenseMatrix ) (Hermitian, complex<double>)";

      HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> > A( 5UL );
      randomize( A );

      HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >    A1( A );
      HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> > A2( A );

      {
         DynamicVector<double,rowVector> w1;
         DynamicVector<double,rowVector> w2;

         DynamicMatrix<complex<double>,rowMajor>    V1;
         DynamicMatrix<complex<double>,columnMajor> V2;

         eigen( A1, w1, V1 );
         eigen( A2, w2, V2 );

         if( !std::is_permutation( w1.begin(), w1.end(), w2.begin(), comparator ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Eigenvalue computation failed\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Row-major eigenvalues:\n" << w1 << "\n"
                << "   Column-major eigenvalues:\n" << w2 << "\n";
            throw std::runtime_error( oss.str() );
         }

         for( size_t i=0UL; i<V1.rows(); ++i ) {
            checkEigenvector( row( V1, i ), A, w1[i] );
         }

         for( size_t i=0UL; i<V2.columns(); ++i ) {
            checkEigenvector( column( V2, i ), A, w2[i] );
         }
      }

      {
         DynamicVector<double,rowVector> w1;
         DynamicVector<double,rowVector> w2;

         DynamicMatrix<complex<double>,rowMajor>    V1;
         DynamicMatrix<complex<double>,columnMajor> V2;

         eigen( trans( A1 ), w1, V1 );
         eigen( trans( A2 ), w2, V2 );

         if( !std::is_permutation( w1.begin(), w1.end(), w2.begin(), comparator ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Transpose eigenvalue computation failed\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Row-major eigenvalues:\n" << w1 << "\n"
                << "   Column-major eigenvalues:\n" << w2 << "\n";
            throw std::runtime_error( oss.str() );
         }

         for( size_t i=0UL; i<V1.rows(); ++i ) {
            checkEigenvector( trans( row( V1, i ) ), trans( A ), w1[i] );
         }

         for( size_t i=0UL; i<V2.columns(); ++i ) {
            checkEigenvector( trans( column( V2, i ) ), trans( A ), w2[i] );
         }
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the eigenvalue/eigenvector evaluation for lower matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function tests the dense matrix eigenvalue/eigenvector evaluation for lower matrices.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testLower()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   using blaze::LowerMatrix;
   using blaze::DynamicMatrix;
   using blaze::DynamicVector;
   using blaze::complex;
   using blaze::eigen;
   using blaze::rowMajor;
   using blaze::columnMajor;
   using blaze::rowVector;

   const auto comparator = []( const auto& v1, const auto& v2 ) {
      return blaze::equal( v1, v2 );
   };


   //=====================================================================================
   // eigen( DenseMatrix, DenseVector )
   //=====================================================================================

   {
      test_ = "eigen( DenseMatrix, DenseVector ) (lower, double)";

      LowerMatrix< DynamicMatrix<double,rowMajor> > A1( 5UL );
      randomize( A1 );
      LowerMatrix< DynamicMatrix<double,columnMajor> > A2( A1 );

      {
         DynamicVector<double,rowVector> w1;
         DynamicVector<double,rowVector> w2;

         eigen( A1, w1 );
         eigen( A2, w2 );

         if( w1 != w2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Eigenvalue computation failed\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Row-major eigenvalues:\n" << w1 << "\n"
                << "   Column-major eigenvalues:\n" << w2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         DynamicVector<double,rowVector> w1;
         DynamicVector<double,rowVector> w2;

         eigen( trans( A1 ), w1 );
         eigen( trans( A2 ), w2 );

         if( w1 != w2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Transpose eigenvalue computation failed\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Row-major eigenvalues:\n" << w1 << "\n"
                << "   Column-major eigenvalues:\n" << w2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   {
      test_ = "eigen( DenseMatrix, DenseVector ) (lower, complex<double>)";

      LowerMatrix< DynamicMatrix<complex<double>,rowMajor> > A1( 5UL );
      randomize( A1 );
      LowerMatrix< DynamicMatrix<complex<double>,columnMajor> > A2( A1 );

      {
         DynamicVector<complex<double>,rowVector> w1;
         DynamicVector<complex<double>,rowVector> w2;

         eigen( A1, w1 );
         eigen( A2, w2 );

         if( w1 != w2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Eigenvalue computation failed\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Row-major eigenvalues:\n" << w1 << "\n"
                << "   Column-major eigenvalues:\n" << w2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         DynamicVector<complex<double>,rowVector> w1;
         DynamicVector<complex<double>,rowVector> w2;

         eigen( trans( A1 ), w1 );
         eigen( trans( A2 ), w2 );

         if( w1 != w2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Transpose eigenvalue computation failed\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Row-major eigenvalues:\n" << w1 << "\n"
                << "   Column-major eigenvalues:\n" << w2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // eigen( DenseMatrix, DenseVector, DenseMatrix )
   //=====================================================================================

   {
      test_ = "eigen( DenseMatrix, DenseVector, DenseMatrix ) (lower, double)";

      LowerMatrix< DynamicMatrix<double,rowMajor> > A( 5UL );
      randomize( A );

      LowerMatrix< DynamicMatrix<double,rowMajor> >    A1( A );
      LowerMatrix< DynamicMatrix<double,columnMajor> > A2( A );

      {
         DynamicVector<complex<double>,rowVector> w1;
         DynamicVector<complex<double>,rowVector> w2;

         DynamicMatrix<complex<double>,rowMajor>    V1;
         DynamicMatrix<complex<double>,columnMajor> V2;

         eigen( A1, w1, V1 );
         eigen( A2, w2, V2 );

         if( !std::is_permutation( w1.begin(), w1.end(), w2.begin(), comparator ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Eigenvalue computation failed\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Row-major eigenvalues:\n" << w1 << "\n"
                << "   Column-major eigenvalues:\n" << w2 << "\n";
            throw std::runtime_error( oss.str() );
         }

         for( size_t i=0UL; i<V1.rows(); ++i ) {
            checkEigenvector( row( V1, i ), A, w1[i] );
         }

         for( size_t i=0UL; i<V2.columns(); ++i ) {
            checkEigenvector( column( V2, i ), A, w2[i] );
         }
      }

      {
         DynamicVector<complex<double>,rowVector> w1;
         DynamicVector<complex<double>,rowVector> w2;

         DynamicMatrix<complex<double>,rowMajor>    V1;
         DynamicMatrix<complex<double>,columnMajor> V2;

         eigen( trans( A1 ), w1, V1 );
         eigen( trans( A2 ), w2, V2 );

         if( !std::is_permutation( w1.begin(), w1.end(), w2.begin(), comparator ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Transpose eigenvalue computation failed\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Row-major eigenvalues:\n" << w1 << "\n"
                << "   Column-major eigenvalues:\n" << w2 << "\n";
            throw std::runtime_error( oss.str() );
         }

         for( size_t i=0UL; i<V1.rows(); ++i ) {
            checkEigenvector( ctrans( row( V1, i ) ), trans( A ), w1[i] );
         }

         for( size_t i=0UL; i<V2.columns(); ++i ) {
            checkEigenvector( ctrans( column( V2, i ) ), trans( A ), w2[i] );
         }
      }
   }

   {
      test_ = "eigen( DenseMatrix, DenseVector, DenseMatrix ) (lower, complex<double>)";

      LowerMatrix< DynamicMatrix<complex<double>,rowMajor> > A( 5UL );
      randomize( A );

      LowerMatrix< DynamicMatrix<complex<double>,rowMajor> >    A1( A );
      LowerMatrix< DynamicMatrix<complex<double>,columnMajor> > A2( A );

      {
         DynamicVector<complex<double>,rowVector> w1;
         DynamicVector<complex<double>,rowVector> w2;

         DynamicMatrix<complex<double>,rowMajor>    V1;
         DynamicMatrix<complex<double>,columnMajor> V2;

         eigen( A1, w1, V1 );
         eigen( A2, w2, V2 );

         if( !std::is_permutation( w1.begin(), w1.end(), w2.begin(), comparator ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Eigenvalue computation failed\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Row-major eigenvalues:\n" << w1 << "\n"
                << "   Column-major eigenvalues:\n" << w2 << "\n";
            throw std::runtime_error( oss.str() );
         }

         for( size_t i=0UL; i<V1.rows(); ++i ) {
            checkEigenvector( row( V1, i ), A, w1[i] );
         }

         for( size_t i=0UL; i<V2.columns(); ++i ) {
            checkEigenvector( column( V2, i ), A, w2[i] );
         }
      }

      {
         DynamicVector<complex<double>,rowVector> w1;
         DynamicVector<complex<double>,rowVector> w2;

         DynamicMatrix<complex<double>,rowMajor>    V1;
         DynamicMatrix<complex<double>,columnMajor> V2;

         eigen( trans( A1 ), w1, V1 );
         eigen( trans( A2 ), w2, V2 );

         if( !std::is_permutation( w1.begin(), w1.end(), w2.begin(), comparator ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Transpose eigenvalue computation failed\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Row-major eigenvalues:\n" << w1 << "\n"
                << "   Column-major eigenvalues:\n" << w2 << "\n";
            throw std::runtime_error( oss.str() );
         }

         for( size_t i=0UL; i<V1.rows(); ++i ) {
            checkEigenvector( trans( row( V1, i ) ), trans( A ), w1[i] );
         }

         for( size_t i=0UL; i<V2.columns(); ++i ) {
            checkEigenvector( trans( column( V2, i ) ), trans( A ), w2[i] );
         }
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the eigenvalue/eigenvector evaluation for upper matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function tests the dense matrix eigenvalue/eigenvector evaluation for upper matrices.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testUpper()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   using blaze::UpperMatrix;
   using blaze::DynamicMatrix;
   using blaze::DynamicVector;
   using blaze::complex;
   using blaze::eigen;
   using blaze::rowMajor;
   using blaze::columnMajor;
   using blaze::rowVector;

   const auto comparator = []( const auto& v1, const auto& v2 ) {
      return blaze::equal( v1, v2 );
   };


   //=====================================================================================
   // eigen( DenseMatrix, DenseVector )
   //=====================================================================================

   {
      test_ = "eigen( DenseMatrix, DenseVector ) (upper, double)";

      UpperMatrix< DynamicMatrix<double,rowMajor> > A1( 5UL );
      randomize( A1 );
      UpperMatrix< DynamicMatrix<double,columnMajor> > A2( A1 );

      {
         DynamicVector<double,rowVector> w1;
         DynamicVector<double,rowVector> w2;

         eigen( A1, w1 );
         eigen( A2, w2 );

         if( w1 != w2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Eigenvalue computation failed\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Row-major eigenvalues:\n" << w1 << "\n"
                << "   Column-major eigenvalues:\n" << w2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         DynamicVector<double,rowVector> w1;
         DynamicVector<double,rowVector> w2;

         eigen( trans( A1 ), w1 );
         eigen( trans( A2 ), w2 );

         if( w1 != w2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Transpose eigenvalue computation failed\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Row-major eigenvalues:\n" << w1 << "\n"
                << "   Column-major eigenvalues:\n" << w2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   {
      test_ = "eigen( DenseMatrix, DenseVector ) (upper, complex<double>)";

      UpperMatrix< DynamicMatrix<complex<double>,rowMajor> > A1( 5UL );
      randomize( A1 );
      UpperMatrix< DynamicMatrix<complex<double>,columnMajor> > A2( A1 );

      {
         DynamicVector<complex<double>,rowVector> w1;
         DynamicVector<complex<double>,rowVector> w2;

         eigen( A1, w1 );
         eigen( A2, w2 );

         if( w1 != w2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Eigenvalue computation failed\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Row-major eigenvalues:\n" << w1 << "\n"
                << "   Column-major eigenvalues:\n" << w2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         DynamicVector<complex<double>,rowVector> w1;
         DynamicVector<complex<double>,rowVector> w2;

         eigen( trans( A1 ), w1 );
         eigen( trans( A2 ), w2 );

         if( w1 != w2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Transpose eigenvalue computation failed\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Row-major eigenvalues:\n" << w1 << "\n"
                << "   Column-major eigenvalues:\n" << w2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // eigen( DenseMatrix, DenseVector, DenseMatrix )
   //=====================================================================================

   {
      test_ = "eigen( DenseMatrix, DenseVector, DenseMatrix ) (upper, double)";

      UpperMatrix< DynamicMatrix<double,rowMajor> > A( 5UL );
      randomize( A );

      UpperMatrix< DynamicMatrix<double,rowMajor> >    A1( A );
      UpperMatrix< DynamicMatrix<double,columnMajor> > A2( A );

      {
         DynamicVector<complex<double>,rowVector> w1;
         DynamicVector<complex<double>,rowVector> w2;

         DynamicMatrix<complex<double>,rowMajor>    V1;
         DynamicMatrix<complex<double>,columnMajor> V2;

         eigen( A1, w1, V1 );
         eigen( A2, w2, V2 );

         if( !std::is_permutation( w1.begin(), w1.end(), w2.begin(), comparator ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Eigenvalue computation failed\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Row-major eigenvalues:\n" << w1 << "\n"
                << "   Column-major eigenvalues:\n" << w2 << "\n";
            throw std::runtime_error( oss.str() );
         }

         for( size_t i=0UL; i<V1.rows(); ++i ) {
            checkEigenvector( row( V1, i ), A, w1[i] );
         }

         for( size_t i=0UL; i<V2.columns(); ++i ) {
            checkEigenvector( column( V2, i ), A, w2[i] );
         }
      }

      {
         DynamicVector<complex<double>,rowVector> w1;
         DynamicVector<complex<double>,rowVector> w2;

         DynamicMatrix<complex<double>,rowMajor>    V1;
         DynamicMatrix<complex<double>,columnMajor> V2;

         eigen( trans( A1 ), w1, V1 );
         eigen( trans( A2 ), w2, V2 );

         if( !std::is_permutation( w1.begin(), w1.end(), w2.begin(), comparator ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Transpose eigenvalue computation failed\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Row-major eigenvalues:\n" << w1 << "\n"
                << "   Column-major eigenvalues:\n" << w2 << "\n";
            throw std::runtime_error( oss.str() );
         }

         for( size_t i=0UL; i<V1.rows(); ++i ) {
            checkEigenvector( ctrans( row( V1, i ) ), trans( A ), w1[i] );
         }

         for( size_t i=0UL; i<V2.columns(); ++i ) {
            checkEigenvector( ctrans( column( V2, i ) ), trans( A ), w2[i] );
         }
      }
   }

   {
      test_ = "eigen( DenseMatrix, DenseVector, DenseMatrix ) (upper, complex<double>)";

      UpperMatrix< DynamicMatrix<complex<double>,rowMajor> > A( 5UL );
      randomize( A );

      UpperMatrix< DynamicMatrix<complex<double>,rowMajor> >    A1( A );
      UpperMatrix< DynamicMatrix<complex<double>,columnMajor> > A2( A );

      {
         DynamicVector<complex<double>,rowVector> w1;
         DynamicVector<complex<double>,rowVector> w2;

         DynamicMatrix<complex<double>,rowMajor>    V1;
         DynamicMatrix<complex<double>,columnMajor> V2;

         eigen( A1, w1, V1 );
         eigen( A2, w2, V2 );

         if( !std::is_permutation( w1.begin(), w1.end(), w2.begin(), comparator ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Eigenvalue computation failed\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Row-major eigenvalues:\n" << w1 << "\n"
                << "   Column-major eigenvalues:\n" << w2 << "\n";
            throw std::runtime_error( oss.str() );
         }

         for( size_t i=0UL; i<V1.rows(); ++i ) {
            checkEigenvector( row( V1, i ), A, w1[i] );
         }

         for( size_t i=0UL; i<V2.columns(); ++i ) {
            checkEigenvector( column( V2, i ), A, w2[i] );
         }
      }

      {
         DynamicVector<complex<double>,rowVector> w1;
         DynamicVector<complex<double>,rowVector> w2;

         DynamicMatrix<complex<double>,rowMajor>    V1;
         DynamicMatrix<complex<double>,columnMajor> V2;

         eigen( trans( A1 ), w1, V1 );
         eigen( trans( A2 ), w2, V2 );

         if( !std::is_permutation( w1.begin(), w1.end(), w2.begin(), comparator ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Transpose eigenvalue computation failed\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Row-major eigenvalues:\n" << w1 << "\n"
                << "   Column-major eigenvalues:\n" << w2 << "\n";
            throw std::runtime_error( oss.str() );
         }

         for( size_t i=0UL; i<V1.rows(); ++i ) {
            checkEigenvector( trans( row( V1, i ) ), trans( A ), w1[i] );
         }

         for( size_t i=0UL; i<V2.columns(); ++i ) {
            checkEigenvector( trans( column( V2, i ) ), trans( A ), w2[i] );
         }
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the eigenvalue/eigenvector evaluation for diagonal matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function tests the dense matrix eigenvalue/eigenvector evaluation for diagonal matrices.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testDiagonal()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   using blaze::DiagonalMatrix;
   using blaze::DynamicMatrix;
   using blaze::DynamicVector;
   using blaze::complex;
   using blaze::eigen;
   using blaze::rowMajor;
   using blaze::columnMajor;
   using blaze::rowVector;

   const auto comparator = []( const auto& v1, const auto& v2 ) {
      return blaze::equal( v1, v2 );
   };


   //=====================================================================================
   // eigen( DenseMatrix, DenseVector )
   //=====================================================================================

   {
      test_ = "eigen( DenseMatrix, DenseVector ) (diagonal, double)";

      DiagonalMatrix< DynamicMatrix<double,rowMajor> > A1( 5UL );
      randomize( A1 );
      DiagonalMatrix< DynamicMatrix<double,columnMajor> > A2( A1 );

      {
         DynamicVector<double,rowVector> w1;
         DynamicVector<double,rowVector> w2;

         eigen( A1, w1 );
         eigen( A2, w2 );

         if( w1 != w2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Eigenvalue computation failed\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Row-major eigenvalues:\n" << w1 << "\n"
                << "   Column-major eigenvalues:\n" << w2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         DynamicVector<double,rowVector> w1;
         DynamicVector<double,rowVector> w2;

         eigen( trans( A1 ), w1 );
         eigen( trans( A2 ), w2 );

         if( w1 != w2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Transpose eigenvalue computation failed\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Row-major eigenvalues:\n" << w1 << "\n"
                << "   Column-major eigenvalues:\n" << w2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   {
      test_ = "eigen( DenseMatrix, DenseVector ) (diagonal, complex<double>)";

      DiagonalMatrix< DynamicMatrix<complex<double>,rowMajor> > A1( 5UL );
      randomize( A1 );
      DiagonalMatrix< DynamicMatrix<complex<double>,columnMajor> > A2( A1 );

      {
         DynamicVector<complex<double>,rowVector> w1;
         DynamicVector<complex<double>,rowVector> w2;

         eigen( A1, w1 );
         eigen( A2, w2 );

         if( w1 != w2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Eigenvalue computation failed\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Row-major eigenvalues:\n" << w1 << "\n"
                << "   Column-major eigenvalues:\n" << w2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         DynamicVector<complex<double>,rowVector> w1;
         DynamicVector<complex<double>,rowVector> w2;

         eigen( trans( A1 ), w1 );
         eigen( trans( A2 ), w2 );

         if( w1 != w2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Transpose eigenvalue computation failed\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Row-major eigenvalues:\n" << w1 << "\n"
                << "   Column-major eigenvalues:\n" << w2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // eigen( DenseMatrix, DenseVector, DenseMatrix )
   //=====================================================================================

   {
      test_ = "eigen( DenseMatrix, DenseVector, DenseMatrix ) (diagonal, double)";

      DiagonalMatrix< DynamicMatrix<double,rowMajor> > A( 5UL );
      randomize( A );

      DiagonalMatrix< DynamicMatrix<double,rowMajor> >    A1( A );
      DiagonalMatrix< DynamicMatrix<double,columnMajor> > A2( A );

      {
         DynamicVector<double,rowVector> w1;
         DynamicVector<double,rowVector> w2;

         DynamicMatrix<complex<double>,rowMajor>    V1;
         DynamicMatrix<complex<double>,columnMajor> V2;

         eigen( A1, w1, V1 );
         eigen( A2, w2, V2 );

         if( !std::is_permutation( w1.begin(), w1.end(), w2.begin(), comparator ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Eigenvalue computation failed\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Row-major eigenvalues:\n" << w1 << "\n"
                << "   Column-major eigenvalues:\n" << w2 << "\n";
            throw std::runtime_error( oss.str() );
         }

         for( size_t i=0UL; i<V1.rows(); ++i ) {
            checkEigenvector( row( V1, i ), A, w1[i] );
         }

         for( size_t i=0UL; i<V2.columns(); ++i ) {
            checkEigenvector( column( V2, i ), A, w2[i] );
         }
      }

      {
         DynamicVector<double,rowVector> w1;
         DynamicVector<double,rowVector> w2;

         DynamicMatrix<complex<double>,rowMajor>    V1;
         DynamicMatrix<complex<double>,columnMajor> V2;

         eigen( trans( A1 ), w1, V1 );
         eigen( trans( A2 ), w2, V2 );

         if( !std::is_permutation( w1.begin(), w1.end(), w2.begin(), comparator ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Transpose eigenvalue computation failed\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Row-major eigenvalues:\n" << w1 << "\n"
                << "   Column-major eigenvalues:\n" << w2 << "\n";
            throw std::runtime_error( oss.str() );
         }

         for( size_t i=0UL; i<V1.rows(); ++i ) {
            checkEigenvector( ctrans( row( V1, i ) ), trans( A ), w1[i] );
         }

         for( size_t i=0UL; i<V2.columns(); ++i ) {
            checkEigenvector( ctrans( column( V2, i ) ), trans( A ), w2[i] );
         }
      }
   }

   {
      test_ = "eigen( DenseMatrix, DenseVector, DenseMatrix ) (diagonal, complex<double>)";

      DiagonalMatrix< DynamicMatrix<complex<double>,rowMajor> > A( 5UL );
      randomize( A );

      DiagonalMatrix< DynamicMatrix<complex<double>,rowMajor> >    A1( A );
      DiagonalMatrix< DynamicMatrix<complex<double>,columnMajor> > A2( A );

      {
         DynamicVector<complex<double>,rowVector> w1;
         DynamicVector<complex<double>,rowVector> w2;

         DynamicMatrix<complex<double>,rowMajor>    V1;
         DynamicMatrix<complex<double>,columnMajor> V2;

         eigen( A1, w1, V1 );
         eigen( A2, w2, V2 );

         if( !std::is_permutation( w1.begin(), w1.end(), w2.begin(), comparator ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Eigenvalue computation failed\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Row-major eigenvalues:\n" << w1 << "\n"
                << "   Column-major eigenvalues:\n" << w2 << "\n";
            throw std::runtime_error( oss.str() );
         }

         for( size_t i=0UL; i<V1.rows(); ++i ) {
            checkEigenvector( row( V1, i ), A, w1[i] );
         }

         for( size_t i=0UL; i<V2.columns(); ++i ) {
            checkEigenvector( column( V2, i ), A, w2[i] );
         }
      }

      {
         DynamicVector<complex<double>,rowVector> w1;
         DynamicVector<complex<double>,rowVector> w2;

         DynamicMatrix<complex<double>,rowMajor>    V1;
         DynamicMatrix<complex<double>,columnMajor> V2;

         eigen( trans( A1 ), w1, V1 );
         eigen( trans( A2 ), w2, V2 );

         if( !std::is_permutation( w1.begin(), w1.end(), w2.begin(), comparator ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Transpose eigenvalue computation failed\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Row-major eigenvalues:\n" << w1 << "\n"
                << "   Column-major eigenvalues:\n" << w2 << "\n";
            throw std::runtime_error( oss.str() );
         }

         for( size_t i=0UL; i<V1.rows(); ++i ) {
            checkEigenvector( trans( row( V1, i ) ), trans( A ), w1[i] );
         }

         for( size_t i=0UL; i<V2.columns(); ++i ) {
            checkEigenvector( trans( column( V2, i ) ), trans( A ), w2[i] );
         }
      }
   }

#endif
}
//*************************************************************************************************

} // namespace eigen

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
   std::cout << "   Running dense matrix eigenvalue test..." << std::endl;

   try
   {
      RUN_DENSE_EIGEN_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during dense matrix eigenvalue test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
