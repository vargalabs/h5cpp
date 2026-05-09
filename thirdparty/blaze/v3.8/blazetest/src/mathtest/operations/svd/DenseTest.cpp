//=================================================================================================
/*!
//  \file src/mathtest/operations/svd/DenseTest.cpp
//  \brief Source file for the dense matrix singular value test
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
#include <sstream>
#include <stdexcept>
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/DynamicVector.h>
#include <blaze/math/shims/Equal.h>
#include <blaze/util/Complex.h>
#include <blaze/util/Random.h>
#include <blazetest/mathtest/operations/svd/DenseTest.h>
#include <blazetest/system/LAPACK.h>

#ifdef BLAZE_USE_HPX_THREADS
#  include <hpx/hpx_main.hpp>
#endif


namespace blazetest {

namespace mathtest {

namespace operations {

namespace svd {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the DenseTest test.
//
// \exception std::runtime_error Error during singular value computation detected.
*/
DenseTest::DenseTest()
{
   testGeneral();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the LU decomposition functionality for general matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function tests the dense matrix LU decomposition for general matrices. In case an error
// is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testGeneral()
{
   using blaze::DynamicMatrix;
   using blaze::DynamicVector;
   using blaze::complex;
   using blaze::svd;
   using blaze::rowMajor;
   using blaze::columnMajor;
   using blaze::rowVector;


#if BLAZETEST_MATHTEST_LAPACK_MODE
   {
      //=====================================================================================
      // svd( DenseMatrix, DenseVector )
      //=====================================================================================

      {
         test_ = "svd( DenseMatrix, DenseVector ) (double)";

         DynamicMatrix<double,columnMajor> A1( 8UL, 5UL );
         randomize( A1 );
         DynamicMatrix<double,rowMajor> A2( A1 );

         DynamicVector<double,rowVector> s1;
         DynamicVector<double,rowVector> s2;

         svd( A1, s1 );
         svd( A2, s2 );

         if( s1 != s2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Singular value computation failed\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Row-major singular values:\n" << s1 << "\n"
                << "   Column-major singular values:\n" << s2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         test_ = "svd( DenseMatrix, DenseVector ) (complex<double>)";

         DynamicMatrix<complex<double>,columnMajor> A1( 8UL, 5UL );
         randomize( A1 );
         DynamicMatrix<complex<double>,rowMajor> A2( A1 );

         DynamicVector<double,rowVector> s1;
         DynamicVector<double,rowVector> s2;

         svd( A1, s1 );
         svd( A2, s2 );

         if( s1 != s2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Singular value computation failed\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Row-major singular values:\n" << s1 << "\n"
                << "   Column-major singular values:\n" << s2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }


      //=====================================================================================
      // svd( DenseMatrix, DenseMatrix, DenseVector, DenseMatrix )
      //=====================================================================================

      {
         test_ = "svd( DenseMatrix, DenseMatrix, DenseVector, DenseMatrix ) (double)";

         DynamicMatrix<double,columnMajor> A1( 8UL, 5UL );
         randomize( A1 );
         DynamicMatrix<double,rowMajor> A2( A1 );

         DynamicVector<double,rowVector> s1;
         DynamicVector<double,rowVector> s2;

         DynamicMatrix<double,columnMajor> U1;
         DynamicMatrix<double,columnMajor> V1;

         DynamicMatrix<double,rowMajor> U2;
         DynamicMatrix<double,rowMajor> V2;

         svd( A1, U1, s1, V1 );
         svd( A2, U2, s2, V2 );

         const blaze::DynamicMatrix<double,blaze::rowMajor> S1{ { s1[0],   0.0,   0.0,   0.0,   0.0 },
                                                                {   0.0, s1[1],   0.0,   0.0,   0.0 },
                                                                {   0.0,   0.0, s1[2],   0.0,   0.0 },
                                                                {   0.0,   0.0,   0.0, s1[3],   0.0 },
                                                                {   0.0,   0.0,   0.0,   0.0, s1[4] } };
         const blaze::DynamicMatrix<double,blaze::rowMajor> S2{ { s2[0],   0.0,   0.0,   0.0,   0.0 },
                                                                {   0.0, s2[1],   0.0,   0.0,   0.0 },
                                                                {   0.0,   0.0, s2[2],   0.0,   0.0 },
                                                                {   0.0,   0.0,   0.0, s2[3],   0.0 },
                                                                {   0.0,   0.0,   0.0,   0.0, s2[4] } };

         if( s1 != s2 || ( U1 * S1 * V1 ) != ( U2 * S2 * V2 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Singular value computation failed\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Row-major singular values:\n" << s1 << "\n"
                << "   Row-major left singular vectors:\n" << U1 << "\n"
                << "   Row-major right singular vectors:\n" << V1 << "\n"
                << "   Column-major singular values:\n" << s2 << "\n"
                << "   Column-major left singular vectors:\n" << U2 << "\n"
                << "   Column-major right singular vectors:\n" << V2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         test_ = "svd( DenseMatrix, DenseMatrix, DenseVector, DenseMatrix ) (complex<double>)";

         DynamicMatrix<complex<double>,columnMajor> A1( 8UL, 5UL );
         randomize( A1 );
         DynamicMatrix<complex<double>,rowMajor> A2( A1 );

         DynamicVector<double,rowVector> s1;
         DynamicVector<double,rowVector> s2;

         DynamicMatrix<complex<double>,columnMajor> U1;
         DynamicMatrix<complex<double>,columnMajor> V1;

         DynamicMatrix<complex<double>,rowMajor> U2;
         DynamicMatrix<complex<double>,rowMajor> V2;

         svd( A1, U1, s1, V1 );
         svd( A2, U2, s2, V2 );

         const blaze::DynamicMatrix<double,blaze::rowMajor> S1{ { s1[0],   0.0,   0.0,   0.0,   0.0 },
                                                                {   0.0, s1[1],   0.0,   0.0,   0.0 },
                                                                {   0.0,   0.0, s1[2],   0.0,   0.0 },
                                                                {   0.0,   0.0,   0.0, s1[3],   0.0 },
                                                                {   0.0,   0.0,   0.0,   0.0, s1[4] } };
         const blaze::DynamicMatrix<double,blaze::rowMajor> S2{ { s2[0],   0.0,   0.0,   0.0,   0.0 },
                                                                {   0.0, s2[1],   0.0,   0.0,   0.0 },
                                                                {   0.0,   0.0, s2[2],   0.0,   0.0 },
                                                                {   0.0,   0.0,   0.0, s2[3],   0.0 },
                                                                {   0.0,   0.0,   0.0,   0.0, s2[4] } };

         if( s1 != s2 || ( U1 * S1 * V1 ) != ( U2 * S2 * V2 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Singular value computation failed\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Row-major singular values:\n" << s1 << "\n"
                << "   Row-major left singular vectors:\n" << U1 << "\n"
                << "   Row-major right singular vectors:\n" << V1 << "\n"
                << "   Column-major singular values:\n" << s2 << "\n"
                << "   Column-major left singular vectors:\n" << U2 << "\n"
                << "   Column-major right singular vectors:\n" << V2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
#endif


#if BLAZETEST_MATHTEST_LAPACK_MODE && BLAZETEST_MATHTEST_LAPACK_SUPPORTS_GESVDX
   {
      //=====================================================================================
      // svd( DenseMatrix, DenseVector, double, double )
      //=====================================================================================

      {
         test_ = "svd( DenseMatrix, DenseVector, double, double ) (double)";

         DynamicMatrix<double,columnMajor> A1( 8UL, 5UL );
         randomize( A1 );
         DynamicMatrix<double,rowMajor> A2( A1 );

         DynamicVector<double,rowVector> s1;
         DynamicVector<double,rowVector> s2;

         svd( A1, s1, 0.0, 0.5 );
         svd( A2, s2, 0.0, 0.5 );

         if( s1 != s2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Singular value computation failed\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Row-major singular values:\n" << s1 << "\n"
                << "   Column-major singular values:\n" << s2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         test_ = "svd( DenseMatrix, DenseVector, double, double ) (complex<double>)";

         DynamicMatrix<complex<double>,columnMajor> A1( 8UL, 5UL );
         randomize( A1 );
         DynamicMatrix<complex<double>,rowMajor> A2( A1 );

         DynamicVector<double,rowVector> s1;
         DynamicVector<double,rowVector> s2;

         svd( A1, s1, 0.0, 0.5 );
         svd( A2, s2, 0.0, 0.5 );

         if( s1 != s2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Singular value computation failed\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Row-major singular values:\n" << s1 << "\n"
                << "   Column-major singular values:\n" << s2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }


      //=====================================================================================
      // svd( DenseMatrix, DenseVector, int, int )
      //=====================================================================================

      {
         test_ = "svd( DenseMatrix, DenseVector, int, int ) (double)";

         DynamicMatrix<double,columnMajor> A1( 8UL, 5UL );
         randomize( A1 );
         DynamicMatrix<double,rowMajor> A2( A1 );

         DynamicVector<double,rowVector> s1;
         DynamicVector<double,rowVector> s2;

         svd( A1, s1, 0, 1 );
         svd( A2, s2, 0, 1 );

         if( s1 != s2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Singular value computation failed\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Row-major singular values:\n" << s1 << "\n"
                << "   Column-major singular values:\n" << s2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         test_ = "svd( DenseMatrix, DenseVector, int, int ) (complex<double>)";

         DynamicMatrix<complex<double>,columnMajor> A1( 8UL, 5UL );
         randomize( A1 );
         DynamicMatrix<complex<double>,rowMajor> A2( A1 );

         DynamicVector<double,rowVector> s1;
         DynamicVector<double,rowVector> s2;

         svd( A1, s1, 0, 1 );
         svd( A2, s2, 0, 1 );

         if( s1 != s2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Singular value computation failed\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Row-major singular values:\n" << s1 << "\n"
                << "   Column-major singular values:\n" << s2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }


      //=====================================================================================
      // svd( DenseMatrix, DenseMatrix, DenseVector, DenseMatrix, double, double )
      //=====================================================================================

      {
         test_ = "svd( DenseMatrix, DenseMatrix, DenseVector, DenseMatrix, double, double ) (double)";

         DynamicMatrix<double,columnMajor> A1( 8UL, 5UL );
         randomize( A1 );
         DynamicMatrix<double,rowMajor> A2( A1 );

         DynamicVector<double,rowVector> s1;
         DynamicVector<double,rowVector> s2;

         DynamicMatrix<double,columnMajor> U1;
         DynamicMatrix<double,columnMajor> V1;

         DynamicMatrix<double,rowMajor> U2;
         DynamicMatrix<double,rowMajor> V2;

         svd( A1, U1, s1, V1, 0.0, 0.5 );
         svd( A2, U2, s2, V2, 0.0, 0.5 );

         if( s1 != s2 || abs( U1 ) != abs( U2 ) || abs( V1 ) != abs( V2 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Singular value computation failed\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Row-major singular values:\n" << s1 << "\n"
                << "   Row-major left singular vectors:\n" << U1 << "\n"
                << "   Row-major right singular vectors:\n" << V1 << "\n"
                << "   Column-major singular values:\n" << s2 << "\n"
                << "   Column-major left singular vectors:\n" << U2 << "\n"
                << "   Column-major right singular vectors:\n" << V2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         test_ = "svd( DenseMatrix, DenseMatrix, DenseVector, DenseMatrix, double, double ) (complex<double>)";

         DynamicMatrix<complex<double>,columnMajor> A1( 8UL, 5UL );
         randomize( A1 );
         DynamicMatrix<complex<double>,rowMajor> A2( A1 );

         DynamicVector<double,rowVector> s1;
         DynamicVector<double,rowVector> s2;

         DynamicMatrix<complex<double>,columnMajor> U1;
         DynamicMatrix<complex<double>,columnMajor> V1;

         DynamicMatrix<complex<double>,rowMajor> U2;
         DynamicMatrix<complex<double>,rowMajor> V2;

         svd( A1, U1, s1, V1, 0.0, 0.5 );
         svd( A2, U2, s2, V2, 0.0, 0.5 );

         if( s1 != s2 || abs( U1 ) != abs( U2 ) || abs( V1 ) != abs( V2 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Singular value computation failed\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Row-major singular values:\n" << s1 << "\n"
                << "   Row-major left singular vectors:\n" << U1 << "\n"
                << "   Row-major right singular vectors:\n" << V1 << "\n"
                << "   Column-major singular values:\n" << s2 << "\n"
                << "   Column-major left singular vectors:\n" << U2 << "\n"
                << "   Column-major right singular vectors:\n" << V2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }


      //=====================================================================================
      // svd( DenseMatrix, DenseMatrix, DenseVector, DenseMatrix, int, int )
      //=====================================================================================

      {
         test_ = "svd( DenseMatrix, DenseMatrix, DenseVector, DenseMatrix, int, int ) (double)";

         DynamicMatrix<double,columnMajor> A1( 8UL, 5UL );
         randomize( A1 );
         DynamicMatrix<double,rowMajor> A2( A1 );

         DynamicVector<double,rowVector> s1;
         DynamicVector<double,rowVector> s2;

         DynamicMatrix<double,columnMajor> U1;
         DynamicMatrix<double,columnMajor> V1;

         DynamicMatrix<double,rowMajor> U2;
         DynamicMatrix<double,rowMajor> V2;

         svd( A1, U1, s1, V1, 0.0, 0.5 );
         svd( A2, U2, s2, V2, 0.0, 0.5 );

         if( s1 != s2 || abs( U1 ) != abs( U2 ) || abs( V1 ) != abs( V2 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Singular value computation failed\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Row-major singular values:\n" << s1 << "\n"
                << "   Row-major left singular vectors:\n" << U1 << "\n"
                << "   Row-major right singular vectors:\n" << V1 << "\n"
                << "   Column-major singular values:\n" << s2 << "\n"
                << "   Column-major left singular vectors:\n" << U2 << "\n"
                << "   Column-major right singular vectors:\n" << V2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         test_ = "svd( DenseMatrix, DenseMatrix, DenseVector, DenseMatrix, int, int ) (complex<double>)";

         DynamicMatrix<complex<double>,columnMajor> A1( 8UL, 5UL );
         randomize( A1 );
         DynamicMatrix<complex<double>,rowMajor> A2( A1 );

         DynamicVector<double,rowVector> s1;
         DynamicVector<double,rowVector> s2;

         DynamicMatrix<complex<double>,columnMajor> U1;
         DynamicMatrix<complex<double>,columnMajor> V1;

         DynamicMatrix<complex<double>,rowMajor> U2;
         DynamicMatrix<complex<double>,rowMajor> V2;

         svd( A1, U1, s1, V1, 0, 1 );
         svd( A2, U2, s2, V2, 0, 1 );

         if( s1 != s2 || abs( U1 ) != abs( U2 ) || abs( V1 ) != abs( V2 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Singular value computation failed\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Row-major singular values:\n" << s1 << "\n"
                << "   Row-major left singular vectors:\n" << U1 << "\n"
                << "   Row-major right singular vectors:\n" << V1 << "\n"
                << "   Column-major singular values:\n" << s2 << "\n"
                << "   Column-major left singular vectors:\n" << U2 << "\n"
                << "   Column-major right singular vectors:\n" << V2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
#endif
}
//*************************************************************************************************

} // namespace svd

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
   std::cout << "   Running dense matrix singular value test..." << std::endl;

   try
   {
      RUN_DENSE_SVD_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during dense matrix singular value test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
