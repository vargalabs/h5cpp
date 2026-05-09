//=================================================================================================
/*!
//  \file src/mathtest/matrices/densematrix/GeneralTest.cpp
//  \brief Source file for the general DenseMatrix operation test
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

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <blaze/math/dense/DenseMatrix.h>
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/DynamicVector.h>
#include <blazetest/mathtest/matrices/densematrix/GeneralTest.h>
#include <blazetest/mathtest/IsEqual.h>

#ifdef BLAZE_USE_HPX_THREADS
#  include <hpx/hpx_main.hpp>
#endif


namespace blazetest {

namespace mathtest {

namespace matrices {

namespace densematrix {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the GeneralTest class test.
//
// \exception std::runtime_error Operation error detected.
*/
GeneralTest::GeneralTest()
{
   testIsNan();
   testIsInf();
   testIsFinite();
   testIsSquare();
   testIsSymmetric();
   testIsHermitian();
   testIsUniform();
   testIsZero();
   testIsLower();
   testIsUniLower();
   testIsStrictlyLower();
   testIsUpper();
   testIsUniUpper();
   testIsStrictlyUpper();
   testIsDiagonal();
   testIsIdentity();
   testIsPositiveDefinite();
   testMinimum();
   testMaximum();
   testTrace();
   testRank();
   testL1Norm();
   testL2Norm();
   testL3Norm();
   testL4Norm();
   testLpNorm();
   testLinfNorm();
   testMean();
   testVar();
   testStdDev();
   testSoftmax();
   testLeftShift();
   testRightShift();
   testBitand();
   testBitor();
   testBitxor();
   testNot();
   testAnd();
   testOr();
   testGenerate();
   testUniform();
   testZero();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the \c isnan() function for dense matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isnan() function for dense matrices. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testIsNan()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major isnan()";

      // isnan with 0x0 matrix (non-NAN)
      {
         blaze::DynamicMatrix<float,blaze::rowMajor> mat;

         checkRows    ( mat, 0UL );
         checkColumns ( mat, 0UL );
         checkNonZeros( mat, 0UL );

         if( blaze::isnan( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isnan evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isnan with 3x5 matrix (non-NAN)
      {
         blaze::DynamicMatrix<float,blaze::rowMajor> mat
            { {  1.2F, -2.1F,  3.2F, -4.1F,  5.2F },
              { -2.1F,  3.2F, -4.1F,  5.2F, -6.1F },
              {  3.2F, -4.1F,  5.2F, -6.1F,  7.2F } };

         checkRows    ( mat,  3UL );
         checkColumns ( mat,  5UL );
         checkNonZeros( mat, 15UL );

         if( blaze::isnan( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isnan evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isnan with 3x5 matrix (NAN)
      {
         blaze::DynamicMatrix<float,blaze::rowMajor> mat
            { {  1.2F, -2.1F,  3.2F, -4.1F,  5.2F },
              { -2.1F,  3.2F, -4.1F,  5.2F, -6.1F },
              {  3.2F, -4.1F,   NAN, -6.1F,  7.2F } };

         checkRows    ( mat,  3UL );
         checkColumns ( mat,  5UL );
         checkNonZeros( mat, 15UL );

         if( blaze::isnan( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isnan evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major isnan()";

      // isnan with 0x0 matrix (non-NAN)
      {
         blaze::DynamicMatrix<float,blaze::columnMajor> mat;

         checkRows    ( mat, 0UL );
         checkColumns ( mat, 0UL );
         checkNonZeros( mat, 0UL );

         if( blaze::isnan( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isnan evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isnan with 3x5 matrix (non-NAN)
      {
         blaze::DynamicMatrix<float,blaze::columnMajor> mat
            { {  1.2F, -2.1F,  3.2F, -4.1F,  5.2F },
              { -2.1F,  3.2F, -4.1F,  5.2F, -6.1F },
              {  3.2F, -4.1F,  5.2F, -6.1F,  7.2F } };

         checkRows    ( mat,  3UL );
         checkColumns ( mat,  5UL );
         checkNonZeros( mat, 15UL );

         if( blaze::isnan( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isnan evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isnan with 3x5 matrix (NAN)
      {
         blaze::DynamicMatrix<float,blaze::columnMajor> mat
            { {  1.2F, -2.1F,  3.2F, -4.1F,  5.2F },
              { -2.1F,  3.2F, -4.1F,  5.2F, -6.1F },
              {  3.2F, -4.1F,   NAN, -6.1F,  7.2F } };

         checkRows    ( mat,  3UL );
         checkColumns ( mat,  5UL );
         checkNonZeros( mat, 15UL );

         if( blaze::isnan( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isnan evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c isinf() function for dense matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isinf() function for dense matrices. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testIsInf()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major isinf()";

      // isinf with 0x0 matrix (non-inf)
      {
         blaze::DynamicMatrix<float,blaze::rowMajor> mat;

         checkRows    ( mat, 0UL );
         checkColumns ( mat, 0UL );
         checkNonZeros( mat, 0UL );

         if( blaze::isinf( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isinf evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isinf with 3x5 matrix (non-inf)
      {
         blaze::DynamicMatrix<float,blaze::rowMajor> mat
            { {  1.2F, -2.1F,  3.2F, -4.1F,  5.2F },
              { -2.1F,  3.2F, -4.1F,  5.2F, -6.1F },
              {  3.2F, -4.1F,  5.2F, -6.1F,  7.2F } };

         checkRows    ( mat,  3UL );
         checkColumns ( mat,  5UL );
         checkNonZeros( mat, 15UL );

         if( blaze::isinf( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isinf evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isinf with 3x5 matrix (HUGE_VALF)
      {
         blaze::DynamicMatrix<float,blaze::rowMajor> mat
            { {  1.2F, -2.1F,      3.2F, -4.1F,  5.2F },
              { -2.1F,  3.2F,     -4.1F,  5.2F, -6.1F },
              {  3.2F, -4.1F, HUGE_VALF, -6.1F,  7.2F } };

         checkRows    ( mat,  3UL );
         checkColumns ( mat,  5UL );
         checkNonZeros( mat, 15UL );

         if( blaze::isinf( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isinf evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major isinf()";

      // isinf with 0x0 matrix (non-inf)
      {
         blaze::DynamicMatrix<float,blaze::columnMajor> mat;

         checkRows    ( mat, 0UL );
         checkColumns ( mat, 0UL );
         checkNonZeros( mat, 0UL );

         if( blaze::isinf( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isinf evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isinf with 3x5 matrix (non-inf)
      {
         blaze::DynamicMatrix<float,blaze::columnMajor> mat
            { {  1.2F, -2.1F,  3.2F, -4.1F,  5.2F },
              { -2.1F,  3.2F, -4.1F,  5.2F, -6.1F },
              {  3.2F, -4.1F,  5.2F, -6.1F,  7.2F } };

         checkRows    ( mat,  3UL );
         checkColumns ( mat,  5UL );
         checkNonZeros( mat, 15UL );

         if( blaze::isinf( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isinf evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isinf with 3x5 matrix (HUGE_VALF)
      {
         blaze::DynamicMatrix<float,blaze::columnMajor> mat
            { {  1.2F, -2.1F,      3.2F, -4.1F,  5.2F },
              { -2.1F,  3.2F,     -4.1F,  5.2F, -6.1F },
              {  3.2F, -4.1F, HUGE_VALF, -6.1F,  7.2F } };

         checkRows    ( mat,  3UL );
         checkColumns ( mat,  5UL );
         checkNonZeros( mat, 15UL );

         if( blaze::isinf( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isnan evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c isfinite() function for dense matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isfinite() function for dense matrices. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testIsFinite()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major isfinite()";

      // isinf with 0x0 matrix (finite)
      {
         blaze::DynamicMatrix<float,blaze::rowMajor> mat;

         checkRows    ( mat, 0UL );
         checkColumns ( mat, 0UL );
         checkNonZeros( mat, 0UL );

         if( blaze::isfinite( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isfinite evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isfinite with 3x5 matrix (finite)
      {
         blaze::DynamicMatrix<float,blaze::rowMajor> mat
            { {  1.2F, -2.1F,  3.2F, -4.1F,  5.2F },
              { -2.1F,  3.2F, -4.1F,  5.2F, -6.1F },
              {  3.2F, -4.1F,  5.2F, -6.1F,  7.2F } };

         checkRows    ( mat,  3UL );
         checkColumns ( mat,  5UL );
         checkNonZeros( mat, 15UL );

         if( blaze::isfinite( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isfinite evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isfinite with 3x5 matrix (NAN)
      {
         blaze::DynamicMatrix<float,blaze::rowMajor> mat
            { {  1.2F, -2.1F,  3.2F, -4.1F,  5.2F },
              { -2.1F,  3.2F, -4.1F,  5.2F, -6.1F },
              {  3.2F, -4.1F,   NAN, -6.1F,  7.2F } };

         checkRows    ( mat,  3UL );
         checkColumns ( mat,  5UL );
         checkNonZeros( mat, 15UL );

         if( blaze::isfinite( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isfinite evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isfinite with 3x5 matrix (HUGE_VALF)
      {
         blaze::DynamicMatrix<float,blaze::rowMajor> mat
            { {  1.2F, -2.1F,      3.2F, -4.1F,  5.2F },
              { -2.1F,  3.2F,     -4.1F,  5.2F, -6.1F },
              {  3.2F, -4.1F, HUGE_VALF, -6.1F,  7.2F } };

         checkRows    ( mat,  3UL );
         checkColumns ( mat,  5UL );
         checkNonZeros( mat, 15UL );

         if( blaze::isfinite( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isfinite evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major isfinite()";

      // isinf with 0x0 matrix (finite)
      {
         blaze::DynamicMatrix<float,blaze::columnMajor> mat;

         checkRows    ( mat, 0UL );
         checkColumns ( mat, 0UL );
         checkNonZeros( mat, 0UL );

         if( blaze::isfinite( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isfinite evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isfinite with 3x5 matrix (finite)
      {
         blaze::DynamicMatrix<float,blaze::columnMajor> mat
            { {  1.2F, -2.1F,  3.2F, -4.1F,  5.2F },
              { -2.1F,  3.2F, -4.1F,  5.2F, -6.1F },
              {  3.2F, -4.1F,  5.2F, -6.1F,  7.2F } };

         checkRows    ( mat,  3UL );
         checkColumns ( mat,  5UL );
         checkNonZeros( mat, 15UL );

         if( blaze::isfinite( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isfinite evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isfinite with 3x5 matrix (NAN)
      {
         blaze::DynamicMatrix<float,blaze::columnMajor> mat
            { {  1.2F, -2.1F,  3.2F, -4.1F,  5.2F },
              { -2.1F,  3.2F, -4.1F,  5.2F, -6.1F },
              {  3.2F, -4.1F,   NAN, -6.1F,  7.2F } };

         checkRows    ( mat,  3UL );
         checkColumns ( mat,  5UL );
         checkNonZeros( mat, 15UL );

         if( blaze::isfinite( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isfinite evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isfinite with 3x5 matrix (HUGE_VALF)
      {
         blaze::DynamicMatrix<float,blaze::columnMajor> mat
            { {  1.2F, -2.1F,      3.2F, -4.1F,  5.2F },
              { -2.1F,  3.2F,     -4.1F,  5.2F, -6.1F },
              {  3.2F, -4.1F, HUGE_VALF, -6.1F,  7.2F } };

         checkRows    ( mat,  3UL );
         checkColumns ( mat,  5UL );
         checkNonZeros( mat, 15UL );

         if( blaze::isfinite( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isfinite evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c isSquare() function for dense matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isSquare() function for dense matrices. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testIsSquare()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major isSquare()";

      // Square matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );

         checkRows   ( mat, 3UL );
         checkColumns( mat, 3UL );

         if( isSquare( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSquare evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-square matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 3UL, 0 );

         checkRows   ( mat, 2UL );
         checkColumns( mat, 3UL );

         if( isSquare( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSquare evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major isSquare()";

      // Square matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );

         checkRows   ( mat, 3UL );
         checkColumns( mat, 3UL );

         if( isSquare( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSquare evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-square matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 2UL, 0 );

         checkRows   ( mat, 3UL );
         checkColumns( mat, 2UL );

         if( isSquare( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSquare evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c isSymmetric() function for dense matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isSymmetric() function for dense matrices. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testIsSymmetric()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major isSymmetric()";

      // Non-square matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 3UL, 0 );

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 6UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );

         if( isSymmetric( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSymmetric evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Default initialized matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 0UL );

         if( isSymmetric( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSymmetric evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Diagonal matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(1,1) = 2;
         mat(2,2) = 3;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 3UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( isSymmetric( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSymmetric evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-symmetric matrix (additional element in the lower part)
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(1,1) = 2;
         mat(2,0) = 4;
         mat(2,2) = 3;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 4UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 2UL );

         if( isSymmetric( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSymmetric evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-symmetric matrix (additional element in the upper part)
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(0,2) = 4;
         mat(1,1) = 2;
         mat(2,2) = 3;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 4UL );
         checkNonZeros( mat, 0UL, 2UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( isSymmetric( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSymmetric evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Symmetric matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(0,2) = 4;
         mat(1,1) = 2;
         mat(2,0) = 4;
         mat(2,2) = 3;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 5UL );
         checkNonZeros( mat, 0UL, 2UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 2UL );

         if( isSymmetric( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSymmetric evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major isSymmetric()";

      // Non-square matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 2UL, 3UL, 0 );

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 6UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 0UL );

         if( isSymmetric( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSymmetric evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Default initialized matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 0UL );

         if( isSymmetric( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSymmetric evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Diagonal matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(1,1) = 2;
         mat(2,2) = 3;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 3UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( isSymmetric( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSymmetric evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-symmetric matrix (additional element in the lower part)
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(1,1) = 2;
         mat(2,0) = 4;
         mat(2,2) = 3;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 4UL );
         checkNonZeros( mat, 0UL, 2UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( isSymmetric( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSymmetric evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-symmetric matrix (additional element in the upper part)
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(0,2) = 4;
         mat(1,1) = 2;
         mat(2,2) = 3;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 4UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 2UL );

         if( isSymmetric( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSymmetric evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Symmetric matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(0,2) = 4;
         mat(1,1) = 2;
         mat(2,0) = 4;
         mat(2,2) = 3;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 5UL );
         checkNonZeros( mat, 0UL, 2UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 2UL );

         if( isSymmetric( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSymmetric evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c isHermitian() function for dense matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isHermitian() function for dense matrices. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testIsHermitian()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major isHermitian()";

      // Non-square matrix
      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 2UL, 3UL, 0 );

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 6UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );

         if( isHermitian( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isHermitian evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Default initialized matrix
      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL, 0.0 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 0UL );

         if( isHermitian( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isHermitian evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-real diagonal element
      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL, 0.0 );
         mat(1,1).imag( 1 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 1UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 0UL );

         if( isHermitian( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isHermitian evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-Hermitian matrix (additional element in the lower part)
      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL, 0.0 );
         mat(0,0).real( 1 );
         mat(1,1).real( 2 );
         mat(2,0).real( 4 );
         mat(2,2).real( 3 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 4UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 2UL );

         if( isHermitian( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isHermitian evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-Hermitian matrix (additional element in the upper part)
      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL, 0.0 );
         mat(0,0).real( 1 );
         mat(0,2).real( 4 );
         mat(1,1).real( 2 );
         mat(2,2).real( 3 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 4UL );
         checkNonZeros( mat, 0UL, 2UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( isHermitian( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isHermitian evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-Hermitian matrix (invalid pair of elements)
      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL, 0.0 );
         mat(0,0).real( 1 );
         mat(0,2).imag( 4 );
         mat(1,1).real( 2 );
         mat(2,0).imag( 4 );
         mat(2,2).real( 3 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 5UL );
         checkNonZeros( mat, 0UL, 2UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 2UL );

         if( isHermitian( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isHermitian evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Hermitian matrix
      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL, 0.0 );
         mat(0,0).real(  1 );
         mat(0,2).imag(  4 );
         mat(1,1).real(  2 );
         mat(2,0).imag( -4 );
         mat(2,2).real(  3 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 5UL );
         checkNonZeros( mat, 0UL, 2UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 2UL );

         if( isHermitian( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isHermitian evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major isHermitian()";

      // Non-square matrix
      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 2UL, 3UL, 0 );

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 6UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );

         if( isHermitian( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isHermitian evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Default initialized matrix
      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL, 0.0 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 0UL );

         if( isHermitian( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isHermitian evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-real diagonal element
      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL, 0.0 );
         mat(1,1).imag( 1 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 1UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 0UL );

         if( isHermitian( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isHermitian evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-Hermitian matrix (additional element in the lower part)
      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL, 0.0 );
         mat(0,0).real( 1 );
         mat(1,1).real( 2 );
         mat(2,0).real( 4 );
         mat(2,2).real( 3 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 4UL );
         checkNonZeros( mat, 0UL, 2UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( isHermitian( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isHermitian evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-Hermitian matrix (additional element in the upper part)
      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL, 0.0 );
         mat(0,0).real( 1 );
         mat(0,2).real( 4 );
         mat(1,1).real( 2 );
         mat(2,2).real( 3 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 4UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 2UL );

         if( isHermitian( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isHermitian evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-Hermitian matrix (invalid pair of elements)
      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL, 0.0 );
         mat(0,0).real( 1 );
         mat(0,2).imag( 4 );
         mat(1,1).real( 2 );
         mat(2,0).imag( 4 );
         mat(2,2).real( 3 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 5UL );
         checkNonZeros( mat, 0UL, 2UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 2UL );

         if( isHermitian( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isHermitian evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Hermitian matrix
      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL, 0.0 );
         mat(0,0).real(  1 );
         mat(0,2).imag(  4 );
         mat(1,1).real(  2 );
         mat(2,0).imag( -4 );
         mat(2,2).real(  3 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 5UL );
         checkNonZeros( mat, 0UL, 2UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 2UL );

         if( isHermitian( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isHermitian evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c isUniform() function for dense matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isUniform() function for dense matrices. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testIsUniform()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major isUniform()";

      // Uniform matrix (0x3)
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 0UL, 3UL, 5 );

         checkRows    ( mat, 0UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 0UL );
         checkNonZeros( mat, 0UL );

         if( isUniform( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUniform evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Uniform matrix (3x0)
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 0UL, 5 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 0UL );
         checkCapacity( mat, 0UL );
         checkNonZeros( mat, 0UL );

         if( isUniform( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUniform evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Uniform matrix (1x3)
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 1UL, 3UL, 5 );

         checkRows    ( mat, 1UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 3UL );
         checkNonZeros( mat, 3UL );
         checkNonZeros( mat, 0UL, 3UL );

         if( isUniform( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUniform evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Uniform matrix (3x1)
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 1UL, 5 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 1UL );
         checkCapacity( mat, 3UL );
         checkNonZeros( mat, 3UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( isUniform( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUniform evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Uniform matrix (3x5)
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 5UL, 5 );

         checkRows    ( mat,  3UL );
         checkColumns ( mat,  5UL );
         checkCapacity( mat, 15UL );
         checkNonZeros( mat, 15UL );
         checkNonZeros( mat,  0UL, 5UL );
         checkNonZeros( mat,  1UL, 5UL );
         checkNonZeros( mat,  2UL, 5UL );

         if( isUniform( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUniform evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Uniform matrix (5x3)
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 5UL, 3UL, 5 );

         checkRows    ( mat,  5UL );
         checkColumns ( mat,  3UL );
         checkCapacity( mat, 15UL );
         checkNonZeros( mat, 15UL );
         checkNonZeros( mat,  0UL, 3UL );
         checkNonZeros( mat,  1UL, 3UL );
         checkNonZeros( mat,  2UL, 3UL );
         checkNonZeros( mat,  3UL, 3UL );
         checkNonZeros( mat,  4UL, 3UL );

         if( isUniform( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUniform evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-uniform matrix (3x3)
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 5 );
         mat(2,2) = 3;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 9UL );
         checkNonZeros( mat, 0UL, 3UL );
         checkNonZeros( mat, 1UL, 3UL );
         checkNonZeros( mat, 2UL, 3UL );

         if( isUniform( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUniform evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major isUniform()";

      // Uniform matrix (0x3)
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 0UL, 3UL, 5 );

         checkRows    ( mat, 0UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 0UL );
         checkNonZeros( mat, 0UL );

         if( isUniform( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUniform evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Uniform matrix (3x0)
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 0UL, 5 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 0UL );
         checkCapacity( mat, 0UL );
         checkNonZeros( mat, 0UL );

         if( isUniform( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUniform evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Uniform matrix (1x3)
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 1UL, 3UL, 5 );

         checkRows    ( mat, 1UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 3UL );
         checkNonZeros( mat, 3UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( isUniform( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUniform evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Uniform matrix (3x1)
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 1UL, 5 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 1UL );
         checkCapacity( mat, 3UL );
         checkNonZeros( mat, 3UL );
         checkNonZeros( mat, 0UL, 3UL );

         if( isUniform( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUniform evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Uniform matrix (3x5)
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 5UL, 5 );

         checkRows    ( mat,  3UL );
         checkColumns ( mat,  5UL );
         checkCapacity( mat, 15UL );
         checkNonZeros( mat, 15UL );
         checkNonZeros( mat,  0UL, 3UL );
         checkNonZeros( mat,  1UL, 3UL );
         checkNonZeros( mat,  2UL, 3UL );
         checkNonZeros( mat,  3UL, 3UL );
         checkNonZeros( mat,  4UL, 3UL );

         if( isUniform( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUniform evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Uniform matrix (5x3)
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 5UL, 3UL, 5 );

         checkRows    ( mat,  5UL );
         checkColumns ( mat,  3UL );
         checkCapacity( mat, 15UL );
         checkNonZeros( mat, 15UL );
         checkNonZeros( mat,  0UL, 5UL );
         checkNonZeros( mat,  1UL, 5UL );
         checkNonZeros( mat,  2UL, 5UL );

         if( isUniform( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUniform evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-uniform matrix (3x3)
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 5 );
         mat(2,2) = 3;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 9UL );
         checkNonZeros( mat, 0UL, 3UL );
         checkNonZeros( mat, 1UL, 3UL );
         checkNonZeros( mat, 2UL, 3UL );

         if( isUniform( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUniform evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c isZero() function for dense matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isZero() function for dense matrices. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testIsZero()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major isZero()";

      // Zero matrix (0x3)
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 0UL, 3UL, 5 );

         checkRows    ( mat, 0UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 0UL );
         checkNonZeros( mat, 0UL );

         if( isZero( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isZero evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Zero matrix (3x0)
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 0UL, 5 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 0UL );
         checkCapacity( mat, 0UL );
         checkNonZeros( mat, 0UL );

         if( isZero( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isZero evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Zero matrix (1x3)
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 1UL, 3UL, 0 );

         checkRows    ( mat, 1UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 3UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );

         if( isZero( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isZero evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Zero matrix (3x1)
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 1UL, 0 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 1UL );
         checkCapacity( mat, 3UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 0UL );

         if( isZero( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isZero evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Zero matrix (3x5)
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 5UL, 0 );

         checkRows    ( mat,  3UL );
         checkColumns ( mat,  5UL );
         checkCapacity( mat, 15UL );
         checkNonZeros( mat,  0UL );
         checkNonZeros( mat,  0UL, 0UL );
         checkNonZeros( mat,  1UL, 0UL );
         checkNonZeros( mat,  2UL, 0UL );

         if( isZero( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isZero evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Zero matrix (5x3)
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 5UL, 3UL, 0 );

         checkRows    ( mat,  5UL );
         checkColumns ( mat,  3UL );
         checkCapacity( mat, 15UL );
         checkNonZeros( mat,  0UL );
         checkNonZeros( mat,  0UL, 0UL );
         checkNonZeros( mat,  1UL, 0UL );
         checkNonZeros( mat,  2UL, 0UL );
         checkNonZeros( mat,  3UL, 0UL );
         checkNonZeros( mat,  4UL, 0UL );

         if( isZero( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isZero evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-Zero matrix (3x3)
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat{ { 0, 0, 0 },
                                                        { 0, 0, 0 },
                                                        { 0, 0, 3 } };

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 1UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( isZero( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isZero evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major isZero()";

      // Zero matrix (0x3)
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 0UL, 3UL, 5 );

         checkRows    ( mat, 0UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 0UL );
         checkNonZeros( mat, 0UL );

         if( isZero( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isZero evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Zero matrix (3x0)
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 0UL, 5 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 0UL );
         checkCapacity( mat, 0UL );
         checkNonZeros( mat, 0UL );

         if( isZero( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isZero evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Zero matrix (1x3)
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 1UL, 3UL, 0 );

         checkRows    ( mat, 1UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 3UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 0UL );

         if( isZero( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isZero evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Zero matrix (3x1)
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 1UL, 0 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 1UL );
         checkCapacity( mat, 3UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );

         if( isZero( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isZero evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Zero matrix (3x5)
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 5UL, 0 );

         checkRows    ( mat,  3UL );
         checkColumns ( mat,  5UL );
         checkCapacity( mat, 15UL );
         checkNonZeros( mat,  0UL );
         checkNonZeros( mat,  0UL, 0UL );
         checkNonZeros( mat,  1UL, 0UL );
         checkNonZeros( mat,  2UL, 0UL );
         checkNonZeros( mat,  3UL, 0UL );
         checkNonZeros( mat,  4UL, 0UL );

         if( isZero( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isZero evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Zero matrix (5x3)
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 5UL, 3UL, 0 );

         checkRows    ( mat,  5UL );
         checkColumns ( mat,  3UL );
         checkCapacity( mat, 15UL );
         checkNonZeros( mat,  0UL );
         checkNonZeros( mat,  0UL, 0UL );
         checkNonZeros( mat,  1UL, 0UL );
         checkNonZeros( mat,  2UL, 0UL );

         if( isZero( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isZero evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-Zero matrix (3x3)
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat{ { 0, 0, 0 },
                                                           { 0, 0, 0 },
                                                           { 0, 0, 3 } };

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 1UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( isZero( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isZero evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c isLower() function for dense matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isLower() function for dense matrices. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testIsLower()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major isLower()";

      // Non-square matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 3UL, 0 );

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 6UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );

         if( isLower( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isLower evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Default initialized matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 0UL );

         if( isLower( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isLower evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Diagonal matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(1,1) = 2;
         mat(2,2) = 3;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 3UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( isLower( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isLower evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-lower triangular matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(0,2) = 2;
         mat(1,0) = 3;
         mat(1,1) = 4;
         mat(2,2) = 5;
         mat(2,0) = 6;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 6UL );
         checkNonZeros( mat, 0UL, 2UL );
         checkNonZeros( mat, 1UL, 2UL );
         checkNonZeros( mat, 2UL, 2UL );

         if( isLower( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isLower evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Lower triangular matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(1,0) = 2;
         mat(1,1) = 3;
         mat(2,2) = 4;
         mat(2,0) = 5;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 5UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 2UL );
         checkNonZeros( mat, 2UL, 2UL );

         if( isLower( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isLower evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major isLower()";

      // Non-square matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 2UL, 0 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 2UL );
         checkCapacity( mat, 6UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );

         if( isLower( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isLower evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Default initialized matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 0UL );

         if( isLower( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isLower evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Diagonal matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(1,1) = 2;
         mat(2,2) = 3;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 3UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( isLower( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isLower evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-lower triangular matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(0,2) = 2;
         mat(1,0) = 3;
         mat(1,1) = 4;
         mat(2,2) = 5;
         mat(2,0) = 6;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 6UL );
         checkNonZeros( mat, 0UL, 3UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 2UL );

         if( isLower( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isLower evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Lower triangular matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(1,0) = 2;
         mat(1,1) = 3;
         mat(2,2) = 4;
         mat(2,0) = 5;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 5UL );
         checkNonZeros( mat, 0UL, 3UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( isLower( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isLower evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c isUniLower() function for dense matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isUniLower() function for dense matrices. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testIsUniLower()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major isUniLower()";

      // Non-square matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 3UL, 0 );

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 6UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );

         if( isUniLower( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUniLower evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Default initialized matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 0UL );

         if( isUniLower( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUniLower evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Identity matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(1,1) = 1;
         mat(2,2) = 1;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 3UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( isUniLower( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUniLower evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Diagonal matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(1,1) = 2;
         mat(2,2) = 3;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 3UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( isUniLower( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUniLower evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Lower unitriangular matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(1,0) = 2;
         mat(1,1) = 1;
         mat(2,2) = 1;
         mat(2,0) = 3;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 5UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 2UL );
         checkNonZeros( mat, 2UL, 2UL );

         if( isUniLower( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUniLower evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Lower triangular matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(1,0) = 2;
         mat(1,1) = 3;
         mat(2,2) = 4;
         mat(2,0) = 5;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 5UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 2UL );
         checkNonZeros( mat, 2UL, 2UL );

         if( isUniLower( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUniLower evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-lower unitriangular matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(0,2) = 2;
         mat(1,0) = 3;
         mat(1,1) = 1;
         mat(2,2) = 1;
         mat(2,0) = 4;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 6UL );
         checkNonZeros( mat, 0UL, 2UL );
         checkNonZeros( mat, 1UL, 2UL );
         checkNonZeros( mat, 2UL, 2UL );

         if( isUniLower( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUniLower evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major isUniLower()";

      // Non-square matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 2UL, 3UL, 0 );

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 6UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 0UL );

         if( isUniLower( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUniLower evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Default initialized matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 0UL );

         if( isUniLower( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUniLower evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Identity matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(1,1) = 1;
         mat(2,2) = 1;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 3UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( isUniLower( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUniLower evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Diagonal matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(1,1) = 2;
         mat(2,2) = 3;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 3UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( isUniLower( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUniLower evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Lower unitriangular matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(1,0) = 2;
         mat(1,1) = 1;
         mat(2,2) = 1;
         mat(2,0) = 3;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 5UL );
         checkNonZeros( mat, 0UL, 3UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( isUniLower( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUniLower evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Lower triangular matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(1,0) = 2;
         mat(1,1) = 3;
         mat(2,2) = 4;
         mat(2,0) = 5;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 5UL );
         checkNonZeros( mat, 0UL, 3UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( isUniLower( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUniLower evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-lower unitriangular matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(0,2) = 2;
         mat(1,0) = 3;
         mat(1,1) = 1;
         mat(2,2) = 1;
         mat(2,0) = 4;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 6UL );
         checkNonZeros( mat, 0UL, 3UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 2UL );

         if( isUniLower( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUniLower evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c isStrictlyLower() function for dense matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isStrictlyLower() function for dense matrices. In case
// an error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testIsStrictlyLower()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major isStrictlyLower()";

      // Non-square matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 3UL, 0 );

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 6UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );

         if( isStrictlyLower( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isStrictlyLower evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Default initialized matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 0UL );

         if( isStrictlyLower( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isStrictlyLower evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Diagonal matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(1,1) = 2;
         mat(2,2) = 3;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 3UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( isStrictlyLower( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isStrictlyLower evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Strictly lower triangular matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
         mat(1,0) = 2;
         mat(2,0) = 5;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 2UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( isStrictlyLower( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isStrictlyLower evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Lower triangular matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(1,0) = 2;
         mat(1,1) = 3;
         mat(2,2) = 4;
         mat(2,0) = 5;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 5UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 2UL );
         checkNonZeros( mat, 2UL, 2UL );

         if( isStrictlyLower( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isStrictlyLower evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-strictly lower triangular matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
         mat(0,2) = 2;
         mat(1,0) = 3;
         mat(2,0) = 4;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 3UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( isStrictlyLower( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isStrictlyLower evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major isStrictlyLower()";

      // Non-square matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 2UL, 3UL, 0 );

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 6UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );

         if( isStrictlyLower( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isStrictlyLower evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Default initialized matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 0UL );

         if( isStrictlyLower( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isStrictlyLower evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Diagonal matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(1,1) = 2;
         mat(2,2) = 3;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 3UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( isStrictlyLower( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isStrictlyLower evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Strictly lower triangular matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
         mat(1,0) = 2;
         mat(2,0) = 5;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 2UL );
         checkNonZeros( mat, 0UL, 2UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 0UL );

         if( isStrictlyLower( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isStrictlyLower evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Lower triangular matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(1,0) = 2;
         mat(1,1) = 3;
         mat(2,2) = 4;
         mat(2,0) = 5;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 5UL );
         checkNonZeros( mat, 0UL, 3UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( isStrictlyLower( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isStrictlyLower evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-strictly lower triangular matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
         mat(0,2) = 2;
         mat(1,0) = 3;
         mat(2,0) = 4;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 3UL );
         checkNonZeros( mat, 0UL, 2UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( isStrictlyLower( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isStrictlyLower evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c isUpper() function for dense matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isUpper() function for dense matrices. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testIsUpper()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major isUpper()";

      // Non-square matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 3UL, 0 );

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 6UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );

         if( isUpper( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUpper evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Default initialized matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 0UL );

         if( isUpper( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUpper evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Diagonal matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(1,1) = 2;
         mat(2,2) = 3;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 3UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( isUpper( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUpper evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-upper triangular matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(0,2) = 2;
         mat(1,1) = 3;
         mat(1,2) = 4;
         mat(2,0) = 5;
         mat(2,2) = 6;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 6UL );
         checkNonZeros( mat, 0UL, 2UL );
         checkNonZeros( mat, 1UL, 2UL );
         checkNonZeros( mat, 2UL, 2UL );

         if( isUpper( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUpper evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Upper triangular matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(0,2) = 2;
         mat(1,1) = 3;
         mat(1,2) = 4;
         mat(2,2) = 5;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 5UL );
         checkNonZeros( mat, 0UL, 2UL );
         checkNonZeros( mat, 1UL, 2UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( isUpper( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUpper evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major isUpper()";

      // Non-square matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 2UL, 0 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 2UL );
         checkCapacity( mat, 6UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );

         if( isUpper( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUpper evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Default initialized matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 0UL );

         if( isUpper( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUpper evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Diagonal matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(1,1) = 2;
         mat(2,2) = 3;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 3UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( isUpper( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUpper evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-upper triangular matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(0,2) = 2;
         mat(1,1) = 3;
         mat(1,2) = 4;
         mat(2,0) = 5;
         mat(2,2) = 6;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 6UL );
         checkNonZeros( mat, 0UL, 2UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 3UL );

         if( isUpper( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUpper evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Upper triangular matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(0,2) = 2;
         mat(1,1) = 3;
         mat(1,2) = 4;
         mat(2,2) = 5;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 5UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 3UL );

         if( isUpper( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUpper evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c isUniUpper() function for dense matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isUniUpper() function for dense matrices. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testIsUniUpper()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major isUniUpper()";

      // Non-square matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 3UL, 0 );

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 6UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );

         if( isUniUpper( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUniUpper evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Default initialized matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 0UL );

         if( isUniUpper( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUniUpper evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Identity matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(1,1) = 1;
         mat(2,2) = 1;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 3UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( isUniUpper( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUniUpper evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Diagonal matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(1,1) = 2;
         mat(2,2) = 3;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 3UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( isUniUpper( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUniUpper evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Upper unitriangular matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(0,2) = 2;
         mat(1,1) = 1;
         mat(1,2) = 3;
         mat(2,2) = 1;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 5UL );
         checkNonZeros( mat, 0UL, 2UL );
         checkNonZeros( mat, 1UL, 2UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( isUniUpper( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUniUpper evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Upper triangular matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(0,2) = 2;
         mat(1,1) = 3;
         mat(1,2) = 4;
         mat(2,2) = 5;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 5UL );
         checkNonZeros( mat, 0UL, 2UL );
         checkNonZeros( mat, 1UL, 2UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( isUniUpper( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUniUpper evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-upper triangular matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(0,2) = 2;
         mat(1,1) = 1;
         mat(1,2) = 3;
         mat(2,0) = 4;
         mat(2,2) = 1;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 6UL );
         checkNonZeros( mat, 0UL, 2UL );
         checkNonZeros( mat, 1UL, 2UL );
         checkNonZeros( mat, 2UL, 2UL );

         if( isUniUpper( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUniUpper evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major isUniUpper()";

      // Non-square matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 2UL, 3UL, 0 );

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 6UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 0UL );

         if( isUniUpper( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUniUpper evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Default initialized matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 0UL );

         if( isUniUpper( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUniUpper evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Identity matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(1,1) = 1;
         mat(2,2) = 1;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 3UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( isUniUpper( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUniUpper evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Diagonal matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(1,1) = 2;
         mat(2,2) = 3;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 3UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( isUniUpper( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUniUpper evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Upper unitriangular matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(0,2) = 2;
         mat(1,1) = 1;
         mat(1,2) = 3;
         mat(2,2) = 1;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 5UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 3UL );

         if( isUniUpper( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUniUpper evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Upper triangular matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(0,2) = 2;
         mat(1,1) = 3;
         mat(1,2) = 4;
         mat(2,2) = 5;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 5UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 3UL );

         if( isUniUpper( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUniUpper evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-upper triangular matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(0,2) = 2;
         mat(1,1) = 1;
         mat(1,2) = 3;
         mat(2,0) = 4;
         mat(2,2) = 1;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 6UL );
         checkNonZeros( mat, 0UL, 2UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 3UL );

         if( isUniUpper( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUniUpper evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c isStrictlyUpper() function for dense matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isStrictlyUpper() function for dense matrices. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testIsStrictlyUpper()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major isStrictlyUpper()";

      // Non-square matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 3UL, 0 );

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 6UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );

         if( isStrictlyUpper( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isStrictlyUpper evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Default initialized matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 0UL );

         if( isStrictlyUpper( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isStrictlyUpper evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Diagonal matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(1,1) = 2;
         mat(2,2) = 3;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 3UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( isStrictlyUpper( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isStrictlyUpper evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Strictly upper triangular matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
         mat(0,2) = 2;
         mat(1,2) = 4;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 2UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 0UL );

         if( isStrictlyUpper( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isStrictlyUpper evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Upper triangular matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(0,2) = 2;
         mat(1,1) = 3;
         mat(1,2) = 4;
         mat(2,2) = 5;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 5UL );
         checkNonZeros( mat, 0UL, 2UL );
         checkNonZeros( mat, 1UL, 2UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( isStrictlyUpper( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isStrictlyUpper evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-strictly upper triangular matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
         mat(0,2) = 2;
         mat(1,2) = 3;
         mat(2,0) = 4;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 3UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( isStrictlyUpper( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isStrictlyUpper evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major isStrictlyUpper()";

      // Non-square matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 2UL, 3UL, 0 );

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 6UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );

         if( isStrictlyUpper( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isStrictlyUpper evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Default initialized matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 0UL );

         if( isStrictlyUpper( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isStrictlyUpper evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Diagonal matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(1,1) = 2;
         mat(2,2) = 3;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 3UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( isStrictlyUpper( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isStrictlyUpper evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Strictly upper triangular matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
         mat(0,2) = 2;
         mat(1,2) = 4;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 2UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 2UL );

         if( isStrictlyUpper( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isStrictlyUpper evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Upper triangular matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(0,2) = 2;
         mat(1,1) = 3;
         mat(1,2) = 4;
         mat(2,2) = 5;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 5UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 3UL );

         if( isStrictlyUpper( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isStrictlyUpper evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-strictly upper triangular matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
         mat(0,2) = 2;
         mat(1,2) = 3;
         mat(2,0) = 4;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 3UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 2UL );

         if( isStrictlyUpper( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isStrictlyUpper evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c isDiagonal() function for dense matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isDiagonal() function for dense matrices. In case
// an error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testIsDiagonal()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major isDiagonal()";

      // Non-square matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 3UL, 0 );

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 6UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );

         if( isDiagonal( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDiagonal evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Default initialized matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 0UL );

         if( isDiagonal( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDiagonal evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Diagonal matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(1,1) = 2;
         mat(2,2) = 3;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 3UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( isDiagonal( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDiagonal evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Lower matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(1,1) = 2;
         mat(2,0) = 4;
         mat(2,2) = 3;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 4UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 2UL );

         if( isDiagonal( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDiagonal evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Upper matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(0,2) = 4;
         mat(1,1) = 2;
         mat(2,2) = 3;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 4UL );
         checkNonZeros( mat, 0UL, 2UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( isDiagonal( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDiagonal evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major isDiagonal()";

      // Non-square matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 2UL, 3UL, 0 );

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 6UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 0UL );

         if( isDiagonal( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDiagonal evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Default initialized matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 0UL );

         if( isDiagonal( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDiagonal evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Diagonal matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(1,1) = 2;
         mat(2,2) = 3;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 3UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( isDiagonal( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDiagonal evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Lower matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(1,1) = 2;
         mat(2,0) = 4;
         mat(2,2) = 3;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 4UL );
         checkNonZeros( mat, 0UL, 2UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( isDiagonal( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDiagonal evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Upper matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(0,2) = 4;
         mat(1,1) = 2;
         mat(2,2) = 3;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 4UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 2UL );

         if( isDiagonal( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDiagonal evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c isIdentity() function for dense matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isIdentity() function for dense matrices. In case
// an error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testIsIdentity()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major isIdentity()";

      // Non-square matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 3UL, 0 );

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 6UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );

         if( isIdentity( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isIdentity evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Default initialized matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 0UL );

         if( isIdentity( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isIdentity evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Identity matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(1,1) = 1;
         mat(2,2) = 1;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 3UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( isIdentity( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isIdentity evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Incomplete identity matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(1,1) = 0;
         mat(2,2) = 1;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 2UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( isIdentity( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isIdentity evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Diagonal matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(1,1) = 2;
         mat(2,2) = 3;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 3UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( isIdentity( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isIdentity evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Lower matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(1,1) = 1;
         mat(2,0) = 2;
         mat(2,2) = 1;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 4UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 2UL );

         if( isIdentity( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isIdentity evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Upper matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(0,2) = 2;
         mat(1,1) = 1;
         mat(2,2) = 1;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 4UL );
         checkNonZeros( mat, 0UL, 2UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( isIdentity( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isIdentity evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major isIdentity()";

      // Non-square matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 2UL, 3UL, 0 );

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 6UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 0UL );

         if( isIdentity( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isIdentity evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Default initialized matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 0UL );

         if( isIdentity( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isIdentity evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Identity matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(1,1) = 1;
         mat(2,2) = 1;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 3UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( isIdentity( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isIdentity evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Incomplete identity matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(1,1) = 0;
         mat(2,2) = 1;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 2UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( isIdentity( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isIdentity evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Diagonal matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(1,1) = 2;
         mat(2,2) = 3;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 3UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( isIdentity( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isIdentity evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Lower matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(1,1) = 1;
         mat(2,0) = 2;
         mat(2,2) = 1;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 4UL );
         checkNonZeros( mat, 0UL, 2UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( isIdentity( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isIdentity evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Upper matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(0,2) = 2;
         mat(1,1) = 1;
         mat(2,2) = 1;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 4UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 2UL );

         if( isIdentity( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isIdentity evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c isPositiveDefinite() function for dense matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isPositiveDefinite() function for dense matrices. In
// case an error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testIsPositiveDefinite()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major isPositiveDefinite()";

      // 0x0 matrix
      {
         blaze::DynamicMatrix<double,blaze::rowMajor> mat( 0UL, 0UL );

         checkRows    ( mat, 0UL );
         checkColumns ( mat, 0UL );
         checkCapacity( mat, 0UL );
         checkNonZeros( mat, 0UL );

         if( isPositiveDefinite( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isPositiveDefinite evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-square matrix
      {
         blaze::DynamicMatrix<double,blaze::rowMajor> mat( 2UL, 3UL, 0 );

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 6UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );

         if( isPositiveDefinite( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isPositiveDefinite evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Positive definite matrix
      {
         blaze::DynamicMatrix<double,blaze::rowMajor> mat
            { {  2.0, -1.0,  0.0 }
            , { -1.0,  2.0, -1.0 }
            , {  0.0, -1.0,  2.0 } };

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 7UL );
         checkNonZeros( mat, 0UL, 2UL );
         checkNonZeros( mat, 1UL, 3UL );
         checkNonZeros( mat, 2UL, 2UL );

         if( isPositiveDefinite( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isPositiveDefinite evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-positive definite matrix
      {
         blaze::DynamicMatrix<double,blaze::rowMajor> mat
            { { 1.0, 2.0, 0.0 }
            , { 2.0, 1.0, 2.0 }
            , { 0.0, 2.0, 1.0 } };

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 7UL );
         checkNonZeros( mat, 0UL, 2UL );
         checkNonZeros( mat, 1UL, 3UL );
         checkNonZeros( mat, 2UL, 2UL );

         if( isPositiveDefinite( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isPositiveDefinite evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major isPositiveDefinite()";

      // 0x0 matrix
      {
         blaze::DynamicMatrix<double,blaze::columnMajor> mat( 0UL, 0UL );

         checkRows    ( mat, 0UL );
         checkColumns ( mat, 0UL );
         checkCapacity( mat, 0UL );
         checkNonZeros( mat, 0UL );

         if( isPositiveDefinite( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isPositiveDefinite evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-square matrix
      {
         blaze::DynamicMatrix<double,blaze::columnMajor> mat( 2UL, 3UL, 0 );

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 6UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );

         if( isPositiveDefinite( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isPositiveDefinite evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Positive definite matrix
      {
         blaze::DynamicMatrix<double,blaze::columnMajor> mat
            { {  2.0, -1.0,  0.0 }
            , { -1.0,  2.0, -1.0 }
            , {  0.0, -1.0,  2.0 } };

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 7UL );
         checkNonZeros( mat, 0UL, 2UL );
         checkNonZeros( mat, 1UL, 3UL );
         checkNonZeros( mat, 2UL, 2UL );

         if( isPositiveDefinite( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isPositiveDefinite evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-positive definite matrix
      {
         blaze::DynamicMatrix<double,blaze::columnMajor> mat
            { { 1.0, 2.0, 0.0 }
            , { 2.0, 1.0, 2.0 }
            , { 0.0, 2.0, 1.0 } };

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 7UL );
         checkNonZeros( mat, 0UL, 2UL );
         checkNonZeros( mat, 1UL, 3UL );
         checkNonZeros( mat, 2UL, 2UL );

         if( isPositiveDefinite( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isPositiveDefinite evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c min() function for dense matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c min() function for dense matrices. In case an error
// is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testMinimum()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major min()";

      // Attempt to find the minimum at the beginning in a fully filled matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 2UL, 0 );
         mat(0,0) = -1;
         mat(0,1) =  2;
         mat(1,0) =  3;
         mat(1,1) =  4;
         mat(2,0) =  5;
         mat(2,1) =  6;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 2UL );
         checkNonZeros( mat, 6UL );

         const int minimum = min( mat );

         if( minimum != -1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: First computation failed\n"
                << " Details:\n"
                << "   Result: " << minimum << "\n"
                << "   Expected result: -1\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Attempt to find the minimum at the end in a fully filled matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 3UL, 0 );
         mat(0,0) =  1;
         mat(0,1) =  2;
         mat(0,2) =  3;
         mat(1,0) =  4;
         mat(1,1) =  5;
         mat(1,2) = -6;

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
         checkNonZeros( mat, 6UL );

         const int minimum = min( mat );

         if( minimum != -6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Second computation failed\n"
                << " Details:\n"
                << "   Result: " << minimum << "\n"
                << "   Expected result: -6\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Attempt to find the minimum at the beginning in a partially filled matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 5UL, 3UL, 0 );
         mat(0,0) = -1;
         mat(0,2) =  2;
         mat(2,1) =  3;
         mat(4,0) =  4;
         mat(4,2) =  5;

         checkRows    ( mat, 5UL );
         checkColumns ( mat, 3UL );
         checkNonZeros( mat, 5UL );

         const int minimum = min( mat );

         if( minimum != -1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Third computation failed\n"
                << " Details:\n"
                << "   Result: " << minimum << "\n"
                << "   Expected result: -1\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Attempt to find the minimum at the end in a partially filled matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 5UL, 0 );
         mat(0,0) =  1;
         mat(0,4) =  2;
         mat(1,2) =  3;
         mat(2,0) =  4;
         mat(2,4) = -5;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 5UL );
         checkNonZeros( mat, 5UL );

         const int minimum = min( mat );

         if( minimum != -5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Fourth computation failed\n"
                << " Details:\n"
                << "   Result: " << minimum << "\n"
                << "   Expected result: -5\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Attempt to detect 0 as the minimum value
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(0,2) = 2;
         mat(1,1) = 3;
         mat(2,0) = 4;
         mat(2,2) = 5;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkNonZeros( mat, 5UL );

         const int minimum = min( mat );

         if( minimum != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Fifth computation failed\n"
                << " Details:\n"
                << "   Result: " << minimum << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major min()";

      // Attempt to find the minimum at the beginning in a partially filled matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 5UL, 3UL, 0 );
         mat(0,0) = -1;
         mat(0,2) =  2;
         mat(2,1) =  3;
         mat(4,0) =  4;
         mat(4,2) =  5;

         checkRows    ( mat, 5UL );
         checkColumns ( mat, 3UL );
         checkNonZeros( mat, 5UL );

         const int minimum = min( mat );

         if( minimum != -1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: First computation failed\n"
                << " Details:\n"
                << "   Result: " << minimum << "\n"
                << "   Expected result: -1\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Attempt to find the minimum at the end in a partially filled matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 5UL, 0 );
         mat(0,0) =  1;
         mat(0,4) =  2;
         mat(1,2) =  3;
         mat(2,0) =  4;
         mat(2,4) = -5;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 5UL );
         checkNonZeros( mat, 5UL );

         const int minimum = min( mat );

         if( minimum != -5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Second computation failed\n"
                << " Details:\n"
                << "   Result: " << minimum << "\n"
                << "   Expected result: -5\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Attempt to find the minimum at the beginning in a partially filled matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 5UL, 3UL, 0 );
         mat(0,0) = -1;
         mat(0,2) =  2;
         mat(2,1) =  3;
         mat(4,0) =  4;
         mat(4,2) =  5;

         checkRows    ( mat, 5UL );
         checkColumns ( mat, 3UL );
         checkNonZeros( mat, 5UL );

         const int minimum = min( mat );

         if( minimum != -1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Third computation failed\n"
                << " Details:\n"
                << "   Result: " << minimum << "\n"
                << "   Expected result: -1\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Attempt to find the minimum at the end in a partially filled matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 5UL, 0 );
         mat(0,0) =  1;
         mat(0,4) =  2;
         mat(1,2) =  3;
         mat(2,0) =  4;
         mat(2,4) = -5;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 5UL );
         checkNonZeros( mat, 5UL );

         const int minimum = min( mat );

         if( minimum != -5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Fourth computation failed\n"
                << " Details:\n"
                << "   Result: " << minimum << "\n"
                << "   Expected result: -5\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Attempt to detect 0 as the minimum value
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = 1;
         mat(0,2) = 2;
         mat(1,1) = 3;
         mat(2,0) = 4;
         mat(2,2) = 5;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkNonZeros( mat, 5UL );

         const int minimum = min( mat );

         if( minimum != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Fifth computation failed\n"
                << " Details:\n"
                << "   Result: " << minimum << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c max() function for dense matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c max() function for dense matrices. In case an error
// is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testMaximum()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major max()";

      // Attempt to find the maximum at the beginning in a fully filled matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 2UL, 0 );
         mat(0,0) =  1;
         mat(0,1) = -2;
         mat(1,0) = -3;
         mat(1,1) = -4;
         mat(2,0) = -5;
         mat(2,1) = -6;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 2UL );
         checkNonZeros( mat, 6UL );

         const int maximum = max( mat );

         if( maximum != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: First computation failed\n"
                << " Details:\n"
                << "   Result: " << maximum << "\n"
                << "   Expected result: 1\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Attempt to find the maximum at the end in a fully filled matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 3UL, 0 );
         mat(0,0) = -1;
         mat(0,1) = -2;
         mat(0,2) = -3;
         mat(1,0) = -4;
         mat(1,1) = -5;
         mat(1,2) =  6;

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
         checkNonZeros( mat, 6UL );

         const int maximum = max( mat );

         if( maximum != 6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Second computation failed\n"
                << " Details:\n"
                << "   Result: " << maximum << "\n"
                << "   Expected result: 6\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Attempt to find the maximum at the beginning in a partially filled matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 5UL, 3UL, 0 );
         mat(0,0) =  1;
         mat(0,2) = -2;
         mat(2,1) = -3;
         mat(4,0) = -4;
         mat(4,2) = -5;

         checkRows    ( mat, 5UL );
         checkColumns ( mat, 3UL );
         checkNonZeros( mat, 5UL );

         const int maximum = max( mat );

         if( maximum != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Third computation failed\n"
                << " Details:\n"
                << "   Result: " << maximum << "\n"
                << "   Expected result: 1\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Attempt to find the maximum at the end in a partially filled matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 5UL, 0 );
         mat(0,0) = -1;
         mat(0,4) = -2;
         mat(1,2) = -3;
         mat(2,0) = -4;
         mat(2,4) =  5;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 5UL );
         checkNonZeros( mat, 5UL );

         const int maximum = max( mat );

         if( maximum != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Fourth computation failed\n"
                << " Details:\n"
                << "   Result: " << maximum << "\n"
                << "   Expected result: 5\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Attempt to detect 0 as the maximum value
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = -1;
         mat(0,2) = -2;
         mat(1,1) = -3;
         mat(2,0) = -4;
         mat(2,2) = -5;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkNonZeros( mat, 5UL );

         const int maximum = max( mat );

         if( maximum != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Fifth computation failed\n"
                << " Details:\n"
                << "   Result: " << maximum << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major max()";

      // Attempt to find the maximum at the beginning in a fully filled matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 2UL, 0 );
         mat(0,0) =  1;
         mat(0,1) = -2;
         mat(1,0) = -3;
         mat(1,1) = -4;
         mat(2,0) = -5;
         mat(2,1) = -6;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 2UL );
         checkNonZeros( mat, 6UL );

         const int maximum = max( mat );

         if( maximum != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: First computation failed\n"
                << " Details:\n"
                << "   Result: " << maximum << "\n"
                << "   Expected result: 1\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Attempt to find the maximum at the end in a fully filled matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 2UL, 3UL, 0 );
         mat(0,0) = -1;
         mat(0,1) = -2;
         mat(0,2) = -3;
         mat(1,0) = -4;
         mat(1,1) = -5;
         mat(1,2) =  6;

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
         checkNonZeros( mat, 6UL );

         const int maximum = max( mat );

         if( maximum != 6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Second computation failed\n"
                << " Details:\n"
                << "   Result: " << maximum << "\n"
                << "   Expected result: 6\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Attempt to find the maximum at the beginning in a partially filled matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 5UL, 3UL, 0 );
         mat(0,0) =  1;
         mat(0,2) = -2;
         mat(2,1) = -3;
         mat(4,0) = -4;
         mat(4,2) = -5;

         checkRows    ( mat, 5UL );
         checkColumns ( mat, 3UL );
         checkNonZeros( mat, 5UL );

         const int maximum = max( mat );

         if( maximum != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Third computation failed\n"
                << " Details:\n"
                << "   Result: " << maximum << "\n"
                << "   Expected result: 1\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Attempt to find the maximum at the end in a partially filled matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 5UL, 0 );
         mat(0,0) = -1;
         mat(0,4) = -2;
         mat(1,2) = -3;
         mat(2,0) = -4;
         mat(2,4) =  5;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 5UL );
         checkNonZeros( mat, 5UL );

         const int maximum = max( mat );

         if( maximum != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Fourth computation failed\n"
                << " Details:\n"
                << "   Result: " << maximum << "\n"
                << "   Expected result: 5\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Attempt to detect 0 as the maximum value
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
         mat(0,0) = -1;
         mat(0,2) = -2;
         mat(1,1) = -3;
         mat(2,0) = -4;
         mat(2,2) = -5;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkNonZeros( mat, 5UL );

         const int maximum = max( mat );

         if( maximum != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Fifth computation failed\n"
                << " Details:\n"
                << "   Result: " << maximum << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c trace() function for dense matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c trace() function for dense matrices. In case an error
// is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testTrace()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major trace()";

      // Determining the trace of a 0x0 matrix
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat;

         checkRows   ( mat, 0UL );
         checkColumns( mat, 0UL );

         const int trace = blaze::trace( mat );

         if( trace != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: First computation failed\n"
                << " Details:\n"
                << "   Result: " << trace << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the trace of a 3x3 matrix
      {
         const blaze::DynamicMatrix<int,blaze::rowMajor> mat{ { -1,  2, -3 }
                                                            , { -4, -5,  6 }
                                                            , {  7, -8, -9 } };

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkNonZeros( mat, 9UL );

         const int trace = blaze::trace( mat );

         if( trace != -15 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Second computation failed\n"
                << " Details:\n"
                << "   Result: " << trace << "\n"
                << "   Expected result: -15\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the trace of a non-square matrix
      try
      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 3UL );

         checkRows   ( mat, 2UL );
         checkColumns( mat, 3UL );

         const int trace = blaze::trace( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Trace computation on a non-square matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << trace << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major trace()";

      // Determining the trace of a 0x0 matrix
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat;

         checkRows   ( mat, 0UL );
         checkColumns( mat, 0UL );

         const int trace = blaze::trace( mat );

         if( trace != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: First computation failed\n"
                << " Details:\n"
                << "   Result: " << trace << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the trace of a 3x3 matrix
      {
         const blaze::DynamicMatrix<int,blaze::columnMajor> mat{ { -1,  2, -3 }
                                                               , { -4, -5,  6 }
                                                               , {  7, -8, -9 } };

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkNonZeros( mat, 9UL );

         const int trace = blaze::trace( mat );

         if( trace != -15 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Second computation failed\n"
                << " Details:\n"
                << "   Result: " << trace << "\n"
                << "   Expected result: -15\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the trace of a non-square matrix
      try
      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 2UL, 3UL );

         checkRows   ( mat, 2UL );
         checkColumns( mat, 3UL );

         const int trace = blaze::trace( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Trace computation on a non-square matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << trace << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c rank() function for dense matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c rank() function for dense matrices. In case an error
// is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testRank()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major rank()";

      // Determining the rank of a 0x0 matrix
      {
         blaze::DynamicMatrix<double,blaze::rowMajor> mat;

         checkRows   ( mat, 0UL );
         checkColumns( mat, 0UL );

         const size_t rank = blaze::rank( mat );

         if( rank != 0UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Rank computation failed\n"
                << " Details:\n"
                << "   Result: " << rank << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the rank of a 3x3 matrix (rank deficient)
      {
         blaze::DynamicMatrix<double,blaze::rowMajor> mat{ { 1.0, 2.0, 3.0 },
                                                           { 0.0, 0.0, 1.0 },
                                                           { 0.0, 0.0, 1.0 } };

         checkRows   ( mat, 3UL );
         checkColumns( mat, 3UL );

         const size_t rank = blaze::rank( mat );

         if( rank != 2UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Rank computation failed\n"
                << " Details:\n"
                << "   Result: " << rank << "\n"
                << "   Expected result: 2\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the rank of a 3x3 matrix (full rank)
      {
         blaze::DynamicMatrix<double,blaze::rowMajor> mat{ { 1.0, 2.0, 3.0 },
                                                           { 0.0, 1.0, 2.0 },
                                                           { 0.0, 0.0, 1.0 } };

         checkRows   ( mat, 3UL );
         checkColumns( mat, 3UL );

         const size_t rank = blaze::rank( mat );

         if( rank != 3UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Rank computation failed\n"
                << " Details:\n"
                << "   Result: " << rank << "\n"
                << "   Expected result: 3\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major rank()";

      // Determining the rank of a 0x0 matrix
      {
         blaze::DynamicMatrix<double,blaze::columnMajor> mat;

         checkRows   ( mat, 0UL );
         checkColumns( mat, 0UL );

         const size_t rank = blaze::rank( mat );

         if( rank != 0UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Rank computation failed\n"
                << " Details:\n"
                << "   Result: " << rank << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the rank of a 3x3 matrix (rank deficient)
      {
         blaze::DynamicMatrix<double,blaze::columnMajor> mat{ { 1.0, 2.0, 3.0 },
                                                              { 0.0, 0.0, 1.0 },
                                                              { 0.0, 0.0, 1.0 } };

         checkRows   ( mat, 3UL );
         checkColumns( mat, 3UL );

         const size_t rank = blaze::rank( mat );

         if( rank != 2UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Rank computation failed\n"
                << " Details:\n"
                << "   Result: " << rank << "\n"
                << "   Expected result: 2\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the rank of a 3x3 matrix (full rank)
      {
         blaze::DynamicMatrix<double,blaze::columnMajor> mat{ { 1.0, 2.0, 3.0 },
                                                              { 0.0, 1.0, 2.0 },
                                                              { 0.0, 0.0, 1.0 } };

         checkRows   ( mat, 3UL );
         checkColumns( mat, 3UL );

         const size_t rank = blaze::rank( mat );

         if( rank != 3UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Rank computation failed\n"
                << " Details:\n"
                << "   Result: " << rank << "\n"
                << "   Expected result: 3\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c l1Norm() function for dense matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c l1Norm() function for dense matrices. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testL1Norm()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "l1Norm() function";

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat;

         const int norm = blaze::l1Norm( mat );

         if( !isEqual( norm, 0 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: L1 norm computation failed\n"
                << " Details:\n"
                << "   Result: " << norm << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 7UL, 0 );

         const int norm = blaze::l1Norm( mat );

         if( !isEqual( norm, 0 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: L1 norm computation failed\n"
                << " Details:\n"
                << "   Result: " << norm << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat{ { 0,  0,  1,  0,  1,  0,  0 },
                                                        { 0, -2,  0,  0,  0, -1,  0 },
                                                        { 0,  0,  0,  2,  0,  0,  0 } };

         const int norm = blaze::l1Norm( mat );

         if( !isEqual( norm, 7 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: L1 norm computation failed\n"
                << " Details:\n"
                << "   Result: " << norm << "\n"
                << "   Expected result: 7\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "l1Norm() function";

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat;

         const int norm = blaze::l1Norm( mat );

         if( !isEqual( norm, 0 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: L1 norm computation failed\n"
                << " Details:\n"
                << "   Result: " << norm << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 7UL, 0 );

         const int norm = blaze::l1Norm( mat );

         if( !isEqual( norm, 0 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: L1 norm computation failed\n"
                << " Details:\n"
                << "   Result: " << norm << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat{ { 0,  0,  0 },
                                                           { 0, -2,  0 },
                                                           { 1,  0,  0 },
                                                           { 0,  0,  2 },
                                                           { 1,  0,  0 },
                                                           { 0, -1,  0 },
                                                           { 0,  0,  0 } };

         const int norm = blaze::l1Norm( mat );

         if( !isEqual( norm, 7 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: L1 norm computation failed\n"
                << " Details:\n"
                << "   Result: " << norm << "\n"
                << "   Expected result: 7\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c l2Norm() function for dense matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c l2Norm() function for dense matrices. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testL2Norm()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "l2Norm() function";

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat;

         const double norm = blaze::l2Norm( mat );

         if( !isEqual( norm, 0.0 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: L2 norm computation failed\n"
                << " Details:\n"
                << "   Result: " << norm << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 7UL, 0 );

         const double norm = blaze::l2Norm( mat );

         if( !isEqual( norm, 0.0 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: L2 norm computation failed\n"
                << " Details:\n"
                << "   Result: " << norm << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat{ { 0,  0,  1,  0,  1, -2,  0 },
                                                        { 0, -2,  0,  0,  0, -1,  0 },
                                                        { 0,  1,  0,  2,  0,  0,  0 } };

         const double norm = blaze::l2Norm( mat );

         if( !isEqual( norm, 4.0 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: L2 norm computation failed\n"
                << " Details:\n"
                << "   Result: " << norm << "\n"
                << "   Expected result: 4\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "l2Norm() function";

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat;

         const double norm = blaze::l2Norm( mat );

         if( !isEqual( norm, 0.0 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: L2 norm computation failed\n"
                << " Details:\n"
                << "   Result: " << norm << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 7UL, 0 );

         const double norm = blaze::l2Norm( mat );

         if( !isEqual( norm, 0.0 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: L2 norm computation failed\n"
                << " Details:\n"
                << "   Result: " << norm << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat{ {  0,  0,  0 },
                                                           {  0, -2,  1 },
                                                           {  1,  0,  0 },
                                                           {  0,  0,  2 },
                                                           {  1,  0,  0 },
                                                           { -2, -1,  0 },
                                                           {  0,  0,  0 } };

         const double norm = blaze::l2Norm( mat );

         if( !isEqual( norm, 4.0 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: L2 norm computation failed\n"
                << " Details:\n"
                << "   Result: " << norm << "\n"
                << "   Expected result: 4\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c l3Norm() function for dense matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c l3Norm() function for dense matrices. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testL3Norm()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "l3Norm() function";

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat;

         const double norm = blaze::l3Norm( mat );

         if( !isEqual( norm, 0.0 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: L3 norm computation failed\n"
                << " Details:\n"
                << "   Result: " << norm << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 7UL, 0 );

         const double norm = blaze::l3Norm( mat );

         if( !isEqual( norm, 0.0 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: L3 norm computation failed\n"
                << " Details:\n"
                << "   Result: " << norm << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat{ { 0,  0,  1,  0,  1, -2,  0 },
                                                        { 0, -2,  0,  0,  0, -1,  0 },
                                                        { 0,  0,  0,  2,  0,  0,  0 } };

         const double norm = blaze::l3Norm( mat );

         if( !isEqual( norm, 3.0 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: L3 norm computation failed\n"
                << " Details:\n"
                << "   Result: " << norm << "\n"
                << "   Expected result: 3\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "l3Norm() function";

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat;

         const double norm = blaze::l3Norm( mat );

         if( !isEqual( norm, 0.0 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: L3 norm computation failed\n"
                << " Details:\n"
                << "   Result: " << norm << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 7UL, 0 );

         const double norm = blaze::l3Norm( mat );

         if( !isEqual( norm, 0.0 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: L3 norm computation failed\n"
                << " Details:\n"
                << "   Result: " << norm << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat{ {  0,  0,  0 },
                                                           {  0, -2,  0 },
                                                           {  1,  0,  0 },
                                                           {  0,  0,  2 },
                                                           {  1,  0,  0 },
                                                           { -2, -1,  0 },
                                                           {  0,  0,  0 } };

         const double norm = blaze::l3Norm( mat );

         if( !isEqual( norm, 3.0 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: L3 norm computation failed\n"
                << " Details:\n"
                << "   Result: " << norm << "\n"
                << "   Expected result: 3\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c l4Norm() function for dense matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c l4Norm() function for dense matrices. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testL4Norm()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "l4Norm() function";

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat;

         const double norm = blaze::l4Norm( mat );

         if( !isEqual( norm, 0.0 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: L4 norm computation failed\n"
                << " Details:\n"
                << "   Result: " << norm << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 7UL, 0 );

         const double norm = blaze::l4Norm( mat );

         if( !isEqual( norm, 0.0 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: L4 norm computation failed\n"
                << " Details:\n"
                << "   Result: " << norm << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat{ { 0,  0,  2,  0,  2, -2,  0 },
                                                        { 0, -2,  0,  0,  0, -1,  0 },
                                                        { 0,  0,  0,  2,  0,  0,  0 } };

         const double norm = blaze::l4Norm( mat );

         if( !isEqual( norm, 3.0 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: L4 norm computation failed\n"
                << " Details:\n"
                << "   Result: " << norm << "\n"
                << "   Expected result: 3\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "l4Norm() function";

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat;

         const double norm = blaze::l4Norm( mat );

         if( !isEqual( norm, 0.0 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: L4 norm computation failed\n"
                << " Details:\n"
                << "   Result: " << norm << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 7UL, 0 );

         const double norm = blaze::l4Norm( mat );

         if( !isEqual( norm, 0.0 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: L4 norm computation failed\n"
                << " Details:\n"
                << "   Result: " << norm << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat{ { 0,  0,  2,  0,  2, -2,  0 },
                                                           { 0, -2,  0,  0,  0, -1,  0 },
                                                           { 0,  0,  0,  2,  0,  0,  0 } };

         const double norm = blaze::l4Norm( mat );

         if( !isEqual( norm, 3.0 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: L4 norm computation failed\n"
                << " Details:\n"
                << "   Result: " << norm << "\n"
                << "   Expected result: 3\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c lpNorm() function for dense matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c lpNorm() function for dense matrices. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testLpNorm()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "lpNorm() function";

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat;

         const double norm1 = blaze::lpNorm( mat, 2 );
         const double norm2 = blaze::lpNorm<2UL>( mat );

         if( !isEqual( norm1, 0.0 ) || !isEqual( norm2, 0.0 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lp norm computation failed\n"
                << " Details:\n"
                << "   lpNorm<2>(): " << norm1 << "\n"
                << "   lpNorm(2): " << norm2 << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 7UL, 0 );

         const double norm1 = blaze::lpNorm( mat, 2 );
         const double norm2 = blaze::lpNorm<2UL>( mat );

         if( !isEqual( norm1, 0.0 ) || !isEqual( norm2, 0.0 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lp norm computation failed\n"
                << " Details:\n"
                << "   lpNorm<2>(): " << norm1 << "\n"
                << "   lpNorm(2): " << norm2 << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 5UL, 10UL );
         randomize( mat, -5, 5 );

         const int norm1( blaze::lpNorm( mat, 1 ) );
         const int norm2( blaze::lpNorm<1UL>( mat ) );
         const int norm3( blaze::l1Norm( mat ) );

         if( !isEqual( norm1, norm3 ) || !isEqual( norm2, norm3 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lp norm computation failed\n"
                << " Details:\n"
                << "   lpNorm<1>(): " << norm1 << "\n"
                << "   lpNorm(1): " << norm2 << "\n"
                << "   Expected result: " << norm3 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 5UL, 10UL );
         randomize( mat, -5, 5 );

         const double norm1( blaze::lpNorm( mat, 2 ) );
         const double norm2( blaze::lpNorm<2UL>( mat ) );
         const double norm3( blaze::l2Norm( mat ) );

         if( !isEqual( norm1, norm3 ) || !isEqual( norm2, norm3 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lp norm computation failed\n"
                << " Details:\n"
                << "   lpNorm<2>(): " << norm1 << "\n"
                << "   lpNorm(2): " << norm2 << "\n"
                << "   Expected result: " << norm3 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 5UL, 10UL );
         randomize( mat, -5, 5 );

         const double norm1( blaze::lpNorm( mat, 3 ) );
         const double norm2( blaze::lpNorm<3UL>( mat ) );
         const double norm3( blaze::l3Norm( mat ) );

         if( !isEqual( norm1, norm3 ) || !isEqual( norm2, norm3 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lp norm computation failed\n"
                << " Details:\n"
                << "   lpNorm<3>(): " << norm1 << "\n"
                << "   lpNorm(3): " << norm2 << "\n"
                << "   Expected result: " << norm3 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 5UL, 10UL );
         randomize( mat, -5, 5 );

         const double norm1( blaze::lpNorm( mat, 4 ) );
         const double norm2( blaze::lpNorm<4UL>( mat ) );
         const double norm3( blaze::l4Norm( mat ) );

         if( !isEqual( norm1, norm3 ) || !isEqual( norm2, norm3 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lp norm computation failed\n"
                << " Details:\n"
                << "   lpNorm<4>(): " << norm1 << "\n"
                << "   lpNorm(4): " << norm2 << "\n"
                << "   Expected result: " << norm3 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "lpNorm() function";

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat;

         const double norm1 = blaze::lpNorm( mat, 2 );
         const double norm2 = blaze::lpNorm<2UL>( mat );

         if( !isEqual( norm1, 0.0 ) || !isEqual( norm2, 0.0 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lp norm computation failed\n"
                << " Details:\n"
                << "   lpNorm<2>(): " << norm1 << "\n"
                << "   lpNorm(2): " << norm2 << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 7UL, 0 );

         const double norm1 = blaze::lpNorm( mat, 2 );
         const double norm2 = blaze::lpNorm<2UL>( mat );

         if( !isEqual( norm1, 0.0 ) || !isEqual( norm2, 0.0 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lp norm computation failed\n"
                << " Details:\n"
                << "   lpNorm<2>(): " << norm1 << "\n"
                << "   lpNorm(2): " << norm2 << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 5UL, 10UL );
         randomize( mat, -5, 5 );

         const int norm1( blaze::lpNorm( mat, 1 ) );
         const int norm2( blaze::lpNorm<1UL>( mat ) );
         const int norm3( blaze::l1Norm( mat ) );

         if( !isEqual( norm1, norm3 ) || !isEqual( norm2, norm3 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lp norm computation failed\n"
                << " Details:\n"
                << "   lpNorm<1>(): " << norm1 << "\n"
                << "   lpNorm(1): " << norm2 << "\n"
                << "   Expected result: " << norm3 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 5UL, 10UL );
         randomize( mat, -5, 5 );

         const double norm1( blaze::lpNorm( mat, 2 ) );
         const double norm2( blaze::lpNorm<2UL>( mat ) );
         const double norm3( blaze::l2Norm( mat ) );

         if( !isEqual( norm1, norm3 ) || !isEqual( norm2, norm3 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lp norm computation failed\n"
                << " Details:\n"
                << "   lpNorm<2>(): " << norm1 << "\n"
                << "   lpNorm(2): " << norm2 << "\n"
                << "   Expected result: " << norm3 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 5UL, 10UL );
         randomize( mat, -5, 5 );

         const double norm1( blaze::lpNorm( mat, 3 ) );
         const double norm2( blaze::lpNorm<3UL>( mat ) );
         const double norm3( blaze::l3Norm( mat ) );

         if( !isEqual( norm1, norm3 ) || !isEqual( norm2, norm3 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lp norm computation failed\n"
                << " Details:\n"
                << "   lpNorm<3>(): " << norm1 << "\n"
                << "   lpNorm(3): " << norm2 << "\n"
                << "   Expected result: " << norm3 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 5UL, 10UL );
         randomize( mat, -5, 5 );

         const double norm1( blaze::lpNorm( mat, 4 ) );
         const double norm2( blaze::lpNorm<4UL>( mat ) );
         const double norm3( blaze::l4Norm( mat ) );

         if( !isEqual( norm1, norm3 ) || !isEqual( norm2, norm3 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lp norm computation failed\n"
                << " Details:\n"
                << "   lpNorm<4>(): " << norm1 << "\n"
                << "   lpNorm(4): " << norm2 << "\n"
                << "   Expected result: " << norm3 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c linfNorm() function for dense matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c linfNorm() function for dense matrices. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testLinfNorm()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "linfNorm() function";

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat;

         const double norm = blaze::linfNorm( mat );

         if( !isEqual( norm, 0.0 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Infinity norm computation failed\n"
                << " Details:\n"
                << "   linfNorm(): " << norm << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 7UL, 0 );

         const double norm = blaze::linfNorm( mat );

         if( !isEqual( norm, 0.0 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Infinity norm computation failed\n"
                << " Details:\n"
                << "   linfNorm(): " << norm << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 5UL, 10UL );
         randomize( mat, -5, 5 );

         const int norm1( blaze::linfNorm( mat ) );
         const int norm2( blaze::max( blaze::abs( mat ) ) );

         if( !isEqual( norm1, norm2 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Infinity norm computation failed\n"
                << " Details:\n"
                << "   linfNorm(): " << norm1 << "\n"
                << "   Expected result: " << norm2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "linfNorm() function";

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat;

         const double norm = blaze::linfNorm( mat );

         if( !isEqual( norm, 0.0 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Infinity norm computation failed\n"
                << " Details:\n"
                << "   linfNorm(): " << norm << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 7UL, 0 );

         const double norm = blaze::linfNorm( mat );

         if( !isEqual( norm, 0.0 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Infinity norm computation failed\n"
                << " Details:\n"
                << "   linfNorm(): " << norm << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 5UL, 10UL );
         randomize( mat, -5, 5 );

         const int norm1( blaze::linfNorm( mat ) );
         const int norm2( blaze::max( blaze::abs( mat ) ) );

         if( !isEqual( norm1, norm2 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Infinity norm computation failed\n"
                << " Details:\n"
                << "   linfNorm(): " << norm1 << "\n"
                << "   Expected result: " << norm2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c mean() function for dense matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c mean() function for dense matrices. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testMean()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major mean()";

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );

         const double mean = blaze::mean( mat );

         if( !isEqual( mean, 0.0 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Mean computation failed\n"
                << " Details:\n"
                << "   Result: " << mean << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat{ { 1, 3, 2 },
                                                        { 2, 6, 4 },
                                                        { 9, 6, 3 } };

         const double mean = blaze::mean( mat );

         if( !isEqual( mean, 4.0 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Mean computation failed\n"
                << " Details:\n"
                << "   Result: " << mean << "\n"
                << "   Expected result: 4\n";
            throw std::runtime_error( oss.str() );
         }
      }

      try {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 0UL );

         const double mean = blaze::mean( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Mean computation of matrix with zero columns succeeded\n"
             << " Details:\n"
             << "   Result:\n" << mean << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      try {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 0UL, 3UL );

         const double mean = blaze::mean( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Mean computation of matrix with zero rows succeeded\n"
             << " Details:\n"
             << "   Result:\n" << mean << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   {
      test_ = "Row-major mean<rowwise>()";

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );

         blaze::DynamicVector<double,blaze::columnVector> mean;
         mean = blaze::mean<blaze::rowwise>( mat );

         if( !isEqual( mean[0], 0.0 ) || !isEqual( mean[1], 0.0 ) || !isEqual( mean[2], 0.0 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Mean computation failed\n"
                << " Details:\n"
                << "   Result: " << trans( mean ) << "\n"
                << "   Expected result: ( 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat{ { 1, 3, 2 },
                                                        { 2, 6, 4 },
                                                        { 9, 6, 3 } };

         blaze::DynamicVector<double,blaze::columnVector> mean;
         mean = blaze::mean<blaze::rowwise>( mat );

         if( !isEqual( mean[0], 2.0 ) || !isEqual( mean[1], 4.0 ) || !isEqual( mean[2], 6.0 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Mean computation failed\n"
                << " Details:\n"
                << "   Result: " << trans( mean ) << "\n"
                << "   Expected result: ( 2 4 6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      try {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 0UL );

         blaze::DynamicVector<double,blaze::columnVector> mean;
         mean = blaze::mean<blaze::rowwise>( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Mean computation of matrix with zero columns succeeded\n"
             << " Details:\n"
             << "   Result:\n" << trans( mean ) << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   {
      test_ = "Row-major mean<columnwise>()";

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );

         blaze::DynamicVector<double,blaze::rowVector> mean;
         mean = blaze::mean<blaze::columnwise>( mat );

         if( !isEqual( mean[0], 0.0 ) || !isEqual( mean[1], 0.0 ) || !isEqual( mean[2], 0.0 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Mean computation failed\n"
                << " Details:\n"
                << "   Result: " << mean << "\n"
                << "   Expected result: ( 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat{ { 1, 3, 2 },
                                                        { 2, 6, 4 },
                                                        { 9, 6, 3 } };

         blaze::DynamicVector<double,blaze::rowVector> mean;
         mean = blaze::mean<blaze::columnwise>( mat );

         if( !isEqual( mean[0], 4.0 ) || !isEqual( mean[1], 5.0 ) || !isEqual( mean[2], 3.0 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Mean computation failed\n"
                << " Details:\n"
                << "   Result: " << mean << "\n"
                << "   Expected result: ( 4 5 3 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      try {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 0UL, 3UL );

         blaze::DynamicVector<double,blaze::rowVector> mean;
         mean = blaze::mean<blaze::columnwise>( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Mean computation of matrix with zero rows succeeded\n"
             << " Details:\n"
             << "   Result:\n" << mean << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major mean()";

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );

         const double mean = blaze::mean( mat );

         if( !isEqual( mean, 0.0 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Mean computation failed\n"
                << " Details:\n"
                << "   Result: " << mean << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat{ { 1, 3, 2 },
                                                           { 2, 6, 4 },
                                                           { 9, 6, 3 } };

         const double mean = blaze::mean( mat );

         if( !isEqual( mean, 4.0 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Mean computation failed\n"
                << " Details:\n"
                << "   Result: " << mean << "\n"
                << "   Expected result: 4\n";
            throw std::runtime_error( oss.str() );
         }
      }

      try {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 0UL );

         const double mean = blaze::mean( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Mean computation of matrix with zero columns succeeded\n"
             << " Details:\n"
             << "   Result:\n" << mean << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      try {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 0UL, 3UL );

         const double mean = blaze::mean( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Mean computation of matrix with zero rows succeeded\n"
             << " Details:\n"
             << "   Result:\n" << mean << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   {
      test_ = "Column-major mean<rowwise>()";

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );

         blaze::DynamicVector<double,blaze::columnVector> mean;
         mean = blaze::mean<blaze::rowwise>( mat );

         if( !isEqual( mean[0], 0.0 ) || !isEqual( mean[1], 0.0 ) || !isEqual( mean[2], 0.0 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Mean computation failed\n"
                << " Details:\n"
                << "   Result: " << trans( mean ) << "\n"
                << "   Expected result: ( 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat{ { 1, 3, 2 },
                                                           { 2, 6, 4 },
                                                           { 9, 6, 3 } };

         blaze::DynamicVector<double,blaze::columnVector> mean;
         mean = blaze::mean<blaze::rowwise>( mat );

         if( !isEqual( mean[0], 2.0 ) || !isEqual( mean[1], 4.0 ) || !isEqual( mean[2], 6.0 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Mean computation failed\n"
                << " Details:\n"
                << "   Result: " << trans( mean ) << "\n"
                << "   Expected result: ( 2 4 6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      try {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 0UL );

         blaze::DynamicVector<double,blaze::columnVector> mean;
         mean = blaze::mean<blaze::rowwise>( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Mean computation of matrix with zero columns succeeded\n"
             << " Details:\n"
             << "   Result:\n" << trans( mean ) << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   {
      test_ = "Column-major mean<columnwise>()";

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );

         blaze::DynamicVector<double,blaze::rowVector> mean;
         mean = blaze::mean<blaze::columnwise>( mat );

         if( !isEqual( mean[0], 0.0 ) || !isEqual( mean[1], 0.0 ) || !isEqual( mean[2], 0.0 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Mean computation failed\n"
                << " Details:\n"
                << "   Result: " << mean << "\n"
                << "   Expected result: ( 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat{ { 1, 3, 2 },
                                                           { 2, 6, 4 },
                                                           { 9, 6, 3 } };

         blaze::DynamicVector<double,blaze::rowVector> mean;
         mean = blaze::mean<blaze::columnwise>( mat );

         if( !isEqual( mean[0], 4.0 ) || !isEqual( mean[1], 5.0 ) || !isEqual( mean[2], 3.0 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Mean computation failed\n"
                << " Details:\n"
                << "   Result: " << mean << "\n"
                << "   Expected result: ( 4 5 3 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      try {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 0UL, 3UL );

         blaze::DynamicVector<double,blaze::rowVector> mean;
         mean = blaze::mean<blaze::columnwise>( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Mean computation of matrix with zero rows succeeded\n"
             << " Details:\n"
             << "   Result:\n" << mean << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c var() function for dense matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c var() function for dense matrices. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testVar()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major var()";

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );

         const double var = blaze::var( mat );

         if( !isEqual( var, 0.0 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Variance computation failed\n"
                << " Details:\n"
                << "   Result: " << var << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat{ { 1, 3, 2 },
                                                        { 2, 6, 4 },
                                                        { 9, 6, 3 } };

         const double var = blaze::var( mat );

         if( !isEqual( var, 6.5 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Variance computation failed\n"
                << " Details:\n"
                << "   Result: " << var << "\n"
                << "   Expected result: 6.5\n";
            throw std::runtime_error( oss.str() );
         }
      }

      try {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 0UL );

         const double var = blaze::var( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Variance computation of matrix with zero columns succeeded\n"
             << " Details:\n"
             << "   Result:\n" << var << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      try {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 0UL, 3UL );

         const double var = blaze::var( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Variance computation of matrix with zero rows succeeded\n"
             << " Details:\n"
             << "   Result:\n" << var << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      try {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 1UL, 1UL );

         const double var = blaze::var( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Variance computation of 1x1 matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << var << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   {
      test_ = "Row-major var<rowwise>()";

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );

         blaze::DynamicVector<double,blaze::columnVector> var;
         var = blaze::var<blaze::rowwise>( mat );

         if( !isEqual( var[0], 0.0 ) || !isEqual( var[1], 0.0 ) || !isEqual( var[2], 0.0 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Variance computation failed\n"
                << " Details:\n"
                << "   Result: " << trans( var ) << "\n"
                << "   Expected result: ( 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat{ { 1, 3, 2 },
                                                        { 2, 6, 4 },
                                                        { 9, 6, 3 } };

         blaze::DynamicVector<double,blaze::columnVector> var;
         var = blaze::var<blaze::rowwise>( mat );

         if( !isEqual( var[0], 1.0 ) || !isEqual( var[1], 4.0 ) || !isEqual( var[2], 9.0 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Variance computation failed\n"
                << " Details:\n"
                << "   Result: " << trans( var ) << "\n"
                << "   Expected result: ( 1 4 9 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      try {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 0UL );

         blaze::DynamicVector<double,blaze::columnVector> var;
         var = blaze::var<blaze::rowwise>( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Variance computation of matrix with zero columns succeeded\n"
             << " Details:\n"
             << "   Result:\n" << trans( var ) << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      try {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 1UL );

         blaze::DynamicVector<double,blaze::columnVector> var;
         var = blaze::var<blaze::rowwise>( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Variance computation of matrix with one column succeeded\n"
             << " Details:\n"
             << "   Result:\n" << trans( var ) << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   {
      test_ = "Row-major var<columnwise>()";

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );

         blaze::DynamicVector<double,blaze::rowVector> var;
         var = blaze::var<blaze::columnwise>( mat );

         if( !isEqual( var[0], 0.0 ) || !isEqual( var[1], 0.0 ) || !isEqual( var[2], 0.0 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Variance computation failed\n"
                << " Details:\n"
                << "   Result: " << var << "\n"
                << "   Expected result: ( 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat{ { 1, 3, 2 },
                                                        { 2, 6, 4 },
                                                        { 9, 6, 3 } };

         blaze::DynamicVector<double,blaze::rowVector> var;
         var = blaze::var<blaze::columnwise>( mat );

         if( !isEqual( var[0], 19.0 ) || !isEqual( var[1], 3.0 ) || !isEqual( var[2], 1.0 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Variance computation failed\n"
                << " Details:\n"
                << "   Result: " << var << "\n"
                << "   Expected result: ( 19 3 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      try {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 0UL, 3UL );

         blaze::DynamicVector<double,blaze::rowVector> var;
         var = blaze::var<blaze::columnwise>( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Variance computation of matrix with zero rows succeeded\n"
             << " Details:\n"
             << "   Result:\n" << var << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      try {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 1UL, 3UL );

         blaze::DynamicVector<double,blaze::rowVector> var;
         var = blaze::var<blaze::columnwise>( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Variance computation of matrix with one row succeeded\n"
             << " Details:\n"
             << "   Result:\n" << var << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major var()";

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );

         const double var = blaze::var( mat );

         if( !isEqual( var, 0.0 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Variance computation failed\n"
                << " Details:\n"
                << "   Result: " << var << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat{ { 1, 3, 2 },
                                                           { 2, 6, 4 },
                                                           { 9, 6, 3 } };

         const double var = blaze::var( mat );

         if( !isEqual( var, 6.5 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Variance computation failed\n"
                << " Details:\n"
                << "   Result: " << var << "\n"
                << "   Expected result: 6.5\n";
            throw std::runtime_error( oss.str() );
         }
      }

      try {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 0UL );

         const double var = blaze::var( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Variance computation of matrix with zero columns succeeded\n"
             << " Details:\n"
             << "   Result:\n" << var << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      try {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 0UL, 3UL );

         const double var = blaze::var( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Variance computation of matrix with zero rows succeeded\n"
             << " Details:\n"
             << "   Result:\n" << var << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      try {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 1UL, 1UL );

         const double var = blaze::var( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Variance computation of 1x1 matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << var << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   {
      test_ = "Column-major var<rowwise>()";

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );

         blaze::DynamicVector<double,blaze::columnVector> var;
         var = blaze::var<blaze::rowwise>( mat );

         if( !isEqual( var[0], 0.0 ) || !isEqual( var[1], 0.0 ) || !isEqual( var[2], 0.0 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Variance computation failed\n"
                << " Details:\n"
                << "   Result: " << trans( var ) << "\n"
                << "   Expected result: ( 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat{ { 1, 3, 2 },
                                                           { 2, 6, 4 },
                                                           { 9, 6, 3 } };

         blaze::DynamicVector<double,blaze::columnVector> var;
         var = blaze::var<blaze::rowwise>( mat );

         if( !isEqual( var[0], 1.0 ) || !isEqual( var[1], 4.0 ) || !isEqual( var[2], 9.0 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Variance computation failed\n"
                << " Details:\n"
                << "   Result: " << trans( var ) << "\n"
                << "   Expected result: ( 1 4 9 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      try {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 0UL );

         blaze::DynamicVector<double,blaze::columnVector> var;
         var = blaze::var<blaze::rowwise>( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Variance computation of matrix with zero columns succeeded\n"
             << " Details:\n"
             << "   Result:\n" << trans( var ) << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      try {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 1UL );

         blaze::DynamicVector<double,blaze::columnVector> var;
         var = blaze::var<blaze::rowwise>( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Variance computation of matrix with one column succeeded\n"
             << " Details:\n"
             << "   Result:\n" << trans( var ) << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   {
      test_ = "Column-major var<columnwise>()";

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );

         blaze::DynamicVector<double,blaze::rowVector> var;
         var = blaze::var<blaze::columnwise>( mat );

         if( !isEqual( var[0], 0.0 ) || !isEqual( var[1], 0.0 ) || !isEqual( var[2], 0.0 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Variance computation failed\n"
                << " Details:\n"
                << "   Result: " << var << "\n"
                << "   Expected result: ( 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat{ { 1, 3, 2 },
                                                           { 2, 6, 4 },
                                                           { 9, 6, 3 } };

         blaze::DynamicVector<double,blaze::rowVector> var;
         var = blaze::var<blaze::columnwise>( mat );

         if( !isEqual( var[0], 19.0 ) || !isEqual( var[1], 3.0 ) || !isEqual( var[2], 1.0 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Variance computation failed\n"
                << " Details:\n"
                << "   Result: " << var << "\n"
                << "   Expected result: ( 19 3 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      try {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 0UL, 3UL );

         blaze::DynamicVector<double,blaze::rowVector> var;
         var = blaze::var<blaze::columnwise>( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Variance computation of matrix with zero rows succeeded\n"
             << " Details:\n"
             << "   Result:\n" << var << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      try {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 1UL, 3UL );

         blaze::DynamicVector<double,blaze::rowVector> var;
         var = blaze::var<blaze::columnwise>( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Variance computation of matrix with one row succeeded\n"
             << " Details:\n"
             << "   Result:\n" << var << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c stddev() function for dense matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c stddev() function for dense matrices. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testStdDev()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major stddev()";

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );

         const double stddev = blaze::stddev( mat );

         if( !isEqual( stddev, 0.0 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Standard deviation computation failed\n"
                << " Details:\n"
                << "   Result: " << stddev << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat{ { 1, 3, 2 },
                                                        { 2, 6, 4 },
                                                        { 9, 6, 3 } };

         const double stddev = blaze::stddev( mat );

         if( !isEqual( stddev, std::sqrt(6.5) ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Standard deviation computation failed\n"
                << " Details:\n"
                << "   Result: " << stddev << "\n"
                << "   Expected result: sqrt(6.5)\n";
            throw std::runtime_error( oss.str() );
         }
      }

      try {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 0UL );

         const double stddev = blaze::stddev( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Standard deviation computation of matrix with zero columns succeeded\n"
             << " Details:\n"
             << "   Result:\n" << stddev << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      try {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 0UL, 3UL );

         const double stddev = blaze::stddev( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Standard deviation computation of matrix with zero rows succeeded\n"
             << " Details:\n"
             << "   Result:\n" << stddev << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      try {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 1UL, 1UL );

         const double stddev = blaze::stddev( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Standard deviation computation of 1x1 matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << stddev << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   {
      test_ = "Row-major stddev<rowwise>()";

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );

         blaze::DynamicVector<double,blaze::columnVector> stddev;
         stddev = blaze::stddev<blaze::rowwise>( mat );

         if( !isEqual( stddev[0], 0.0 ) || !isEqual( stddev[1], 0.0 ) || !isEqual( stddev[2], 0.0 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Standard deviation computation failed\n"
                << " Details:\n"
                << "   Result: " << trans( stddev ) << "\n"
                << "   Expected result: ( 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat{ { 1, 3, 2 },
                                                        { 2, 6, 4 },
                                                        { 9, 6, 3 } };

         blaze::DynamicVector<double,blaze::columnVector> stddev;
         stddev = blaze::stddev<blaze::rowwise>( mat );

         if( !isEqual( stddev[0], 1.0 ) || !isEqual( stddev[1], 2.0 ) || !isEqual( stddev[2], 3.0 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Standard deviation computation failed\n"
                << " Details:\n"
                << "   Result: " << trans( stddev ) << "\n"
                << "   Expected result: ( 1 2 3 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      try {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 0UL );

         blaze::DynamicVector<double,blaze::columnVector> stddev;
         stddev = blaze::stddev<blaze::rowwise>( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Standard deviation computation of matrix with zero columns succeeded\n"
             << " Details:\n"
             << "   Result:\n" << trans( stddev ) << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      try {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 1UL );

         blaze::DynamicVector<double,blaze::columnVector> stddev;
         stddev = blaze::stddev<blaze::rowwise>( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Standard deviation computation of matrix with one column succeeded\n"
             << " Details:\n"
             << "   Result:\n" << trans( stddev ) << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   {
      test_ = "Row-major stddev<columnwise>()";

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );

         blaze::DynamicVector<double,blaze::rowVector> stddev;
         stddev = blaze::stddev<blaze::columnwise>( mat );

         if( !isEqual( stddev[0], 0.0 ) || !isEqual( stddev[1], 0.0 ) || !isEqual( stddev[2], 0.0 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Standard deviation computation failed\n"
                << " Details:\n"
                << "   Result: " << stddev << "\n"
                << "   Expected result: ( 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat{ { 1, 3, 2 },
                                                        { 2, 6, 4 },
                                                        { 9, 6, 3 } };

         blaze::DynamicVector<double,blaze::rowVector> stddev;
         stddev = blaze::stddev<blaze::columnwise>( mat );

         if( !isEqual( stddev[0], std::sqrt(19.0) ) || !isEqual( stddev[1], std::sqrt(3.0) ) ||
             !isEqual( stddev[2], 1.0 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Standard deviation computation failed\n"
                << " Details:\n"
                << "   Result: " << stddev << "\n"
                << "   Expected result: ( sqrt(19) sqrt(3) 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      try {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 0UL, 3UL );

         blaze::DynamicVector<double,blaze::rowVector> stddev;
         stddev = blaze::stddev<blaze::columnwise>( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Standard deviation computation of matrix with zero rows succeeded\n"
             << " Details:\n"
             << "   Result:\n" << stddev << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      try {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 1UL, 3UL );

         blaze::DynamicVector<double,blaze::rowVector> stddev;
         stddev = blaze::stddev<blaze::columnwise>( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Standard deviation computation of matrix with one row succeeded\n"
             << " Details:\n"
             << "   Result:\n" << stddev << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major stddev()";

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );

         const double stddev = blaze::stddev( mat );

         if( !isEqual( stddev, 0.0 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Standard deviation computation failed\n"
                << " Details:\n"
                << "   Result: " << stddev << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat{ { 1, 3, 2 },
                                                           { 2, 6, 4 },
                                                           { 9, 6, 3 } };

         const double stddev = blaze::stddev( mat );

         if( !isEqual( stddev, std::sqrt(6.5) ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Standard deviation computation failed\n"
                << " Details:\n"
                << "   Result: " << stddev << "\n"
                << "   Expected result: sqrt(6.5)\n";
            throw std::runtime_error( oss.str() );
         }
      }

      try {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 0UL );

         const double stddev = blaze::stddev( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Standard deviation computation of matrix with zero columns succeeded\n"
             << " Details:\n"
             << "   Result:\n" << stddev << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      try {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 0UL, 3UL );

         const double stddev = blaze::stddev( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Standard deviation computation of matrix with zero rows succeeded\n"
             << " Details:\n"
             << "   Result:\n" << stddev << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      try {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 1UL, 1UL );

         const double stddev = blaze::stddev( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Standard deviation computation of 1x1 matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << stddev << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   {
      test_ = "Column-major stddev<rowwise>()";

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );

         blaze::DynamicVector<double,blaze::columnVector> stddev;
         stddev = blaze::stddev<blaze::rowwise>( mat );

         if( !isEqual( stddev[0], 0.0 ) || !isEqual( stddev[1], 0.0 ) || !isEqual( stddev[2], 0.0 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Standard deviation computation failed\n"
                << " Details:\n"
                << "   Result: " << trans( stddev ) << "\n"
                << "   Expected result: ( 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat{ { 1, 3, 2 },
                                                           { 2, 6, 4 },
                                                           { 9, 6, 3 } };

         blaze::DynamicVector<double,blaze::columnVector> stddev;
         stddev = blaze::stddev<blaze::rowwise>( mat );

         if( !isEqual( stddev[0], 1.0 ) || !isEqual( stddev[1], 2.0 ) || !isEqual( stddev[2], 3.0 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Standard deviation computation failed\n"
                << " Details:\n"
                << "   Result: " << trans( stddev ) << "\n"
                << "   Expected result: ( 1 2 3 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      try {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 0UL );

         blaze::DynamicVector<double,blaze::columnVector> stddev;
         stddev = blaze::stddev<blaze::rowwise>( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Standard deviation computation of matrix with zero columns succeeded\n"
             << " Details:\n"
             << "   Result:\n" << trans( stddev ) << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      try {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 1UL );

         blaze::DynamicVector<double,blaze::columnVector> stddev;
         stddev = blaze::stddev<blaze::rowwise>( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Standard deviation computation of matrix with one column succeeded\n"
             << " Details:\n"
             << "   Result:\n" << trans( stddev ) << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   {
      test_ = "Column-major stddev<columnwise>()";

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );

         blaze::DynamicVector<double,blaze::rowVector> stddev;
         stddev = blaze::stddev<blaze::columnwise>( mat );

         if( !isEqual( stddev[0], 0.0 ) || !isEqual( stddev[1], 0.0 ) || !isEqual( stddev[2], 0.0 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Standard deviation computation failed\n"
                << " Details:\n"
                << "   Result: " << stddev << "\n"
                << "   Expected result: ( 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat{ { 1, 3, 2 },
                                                           { 2, 6, 4 },
                                                           { 9, 6, 3 } };

         blaze::DynamicVector<double,blaze::rowVector> stddev;
         stddev = blaze::stddev<blaze::columnwise>( mat );

         if( !isEqual( stddev[0], std::sqrt(19.0) ) || !isEqual( stddev[1], std::sqrt(3.0) ) ||
             !isEqual( stddev[2], 1.0 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Standard deviation computation failed\n"
                << " Details:\n"
                << "   Result: " << stddev << "\n"
                << "   Expected result: ( sqrt(19) sqrt(3) 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      try {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 0UL, 3UL );

         blaze::DynamicVector<double,blaze::rowVector> stddev;
         stddev = blaze::stddev<blaze::columnwise>( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Standard deviation computation of matrix with zero rows succeeded\n"
             << " Details:\n"
             << "   Result:\n" << stddev << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      try {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 1UL, 3UL );

         blaze::DynamicVector<double,blaze::rowVector> stddev;
         stddev = blaze::stddev<blaze::columnwise>( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Standard deviation computation of matrix with one row succeeded\n"
             << " Details:\n"
             << "   Result:\n" << stddev << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c softmax() function for dense matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c softmax() function for dense matrices. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testSoftmax()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major softmax()";

      blaze::DynamicMatrix<double,blaze::rowMajor> A( 2UL, 2UL );
      randomize( A, -5.0, 5.0 );

      const auto B = softmax( A );

      if( B(0,0) <= 0.0 || B(0,0) > 1.0 ||
          B(0,1) <= 0.0 || B(0,1) > 1.0 ||
          B(1,0) <= 0.0 || B(1,0) > 1.0 ||
          B(1,1) <= 0.0 || B(1,1) > 1.0 ||
          !isEqual( sum( B ), 1.0 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Softmax computation failed\n"
             << " Details:\n"
             << "   Result: " << sum( B ) << "\n"
             << "   Expected result: 1\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major softmax()";

      blaze::DynamicMatrix<double,blaze::columnMajor> A( 2UL, 2UL );
      randomize( A, -5.0, 5.0 );

      const auto B = softmax( A );

      if( B(0,0) <= 0.0 || B(0,0) > 1.0 ||
          B(0,1) <= 0.0 || B(0,1) > 1.0 ||
          B(1,0) <= 0.0 || B(1,0) > 1.0 ||
          B(1,1) <= 0.0 || B(1,1) > 1.0 ||
          !isEqual( sum( B ), 1.0 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Softmax computation failed\n"
             << " Details:\n"
             << "   Result: " << sum( B ) << "\n"
             << "   Expected result: 1\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the left-shift operator for dense matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the left-shift operator for dense matrices. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testLeftShift()
{
   //=====================================================================================
   // Row-major matrix/scalar left-shift tests
   //=====================================================================================

   {
      test_ = "Row-major matrix/scalar left-shift operator";

      // Matrix/scalar left-shift of an empty matrix
      {
         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> A;

         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> B( A << 2U );

         checkRows    ( B, 0UL );
         checkColumns ( B, 0UL );
         checkCapacity( B, 0UL );
         checkNonZeros( B, 0UL );
      }

      // Matrix/scalar left-shift of a general matrix
      {
         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> A{ { 1U, 2U,  4U,  8U, 16U },
                                                               { 2U, 4U,  8U, 16U, 32U },
                                                               { 4U, 8U, 16U, 32U, 64U } };

         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> B( A << 2U );

         checkRows    ( B,  3UL );
         checkColumns ( B,  5UL );
         checkCapacity( B, 15UL );
         checkNonZeros( B, 15UL );

         if( B(0,0) !=  4U || B(0,1) !=  8U || B(0,2) != 16U || B(0,3) !=  32U || B(0,4) !=  64U ||
             B(1,0) !=  8U || B(1,1) != 16U || B(1,2) != 32U || B(1,3) !=  64U || B(1,4) != 128U ||
             B(2,0) != 16U || B(2,1) != 32U || B(2,2) != 64U || B(2,3) != 128U || B(2,4) != 256U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/scalar left-shift operation failed\n"
                << " Details:\n"
                << "   Result:\n" << B << "\n"
                << "   Expected result:\n(  4  8 16  32  64 )\n"
                                        "(  8 16 32  64 128 )\n"
                                        "( 16 32 64 128 256 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Matrix/scalar left-shift assignment
      {
         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> A{ { 1U, 2U,  4U,  8U, 16U },
                                                               { 2U, 4U,  8U, 16U, 32U },
                                                               { 4U, 8U, 16U, 32U, 64U } };

         A <<= 2U;

         checkRows    ( A,  3UL );
         checkColumns ( A,  5UL );
         checkCapacity( A, 15UL );
         checkNonZeros( A, 15UL );

         if( A(0,0) !=  4U || A(0,1) !=  8U || A(0,2) != 16U || A(0,3) !=  32U || A(0,4) !=  64U ||
             A(1,0) !=  8U || A(1,1) != 16U || A(1,2) != 32U || A(1,3) !=  64U || A(1,4) != 128U ||
             A(2,0) != 16U || A(2,1) != 32U || A(2,2) != 64U || A(2,3) != 128U || A(2,4) != 256U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/scalar left-shift assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << A << "\n"
                << "   Expected result:\n(  4  8 16  32  64 )\n"
                                        "(  8 16 32  64 128 )\n"
                                        "( 16 32 64 128 256 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Row-major matrix/row-major matrix left-shift tests
   //=====================================================================================

   {
      test_ = "Row-major matrix/row-major matrix left-shift operator";

      // Matrix/matrix left-shift of an empty matrix
      {
         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> A;
         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> B;

         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> C( A << B );

         checkRows    ( C, 0UL );
         checkColumns ( C, 0UL );
         checkCapacity( C, 0UL );
         checkNonZeros( C, 0UL );
      }

      // Matrix/matrix left-shift of a general matrix
      {
         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> A{ { 1U, 2U,  4U,  8U, 16U },
                                                               { 2U, 4U,  8U, 16U, 32U },
                                                               { 4U, 8U, 16U, 32U, 64U } };

         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> B{ { 1U, 2U, 1U, 2U, 1U },
                                                               { 2U, 1U, 2U, 1U, 2U },
                                                               { 1U, 2U, 1U, 2U, 1U } };

         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> C( A << B );

         checkRows    ( C,  3UL );
         checkColumns ( C,  5UL );
         checkCapacity( C, 15UL );
         checkNonZeros( C, 15UL );

         if( C(0,0) != 2U || C(0,1) !=  8U || C(0,2) !=  8U || C(0,3) !=  32U || C(0,4) !=  32U ||
             C(1,0) != 8U || C(1,1) !=  8U || C(1,2) != 32U || C(1,3) !=  32U || C(1,4) != 128U ||
             C(2,0) != 8U || C(2,1) != 32U || C(2,2) != 32U || C(2,3) != 128U || C(2,4) != 128U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/matrix left-shift operation failed\n"
                << " Details:\n"
                << "   Result:\n" << C << "\n"
                << "   Expected result:\n( 2  8  8  32  32 )\n"
                                        "( 8  8 32  32 128 )\n"
                                        "( 8 32 32 128 128 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Matrix/matrix left-shift assignment
      {
         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> A{ { 1U, 2U,  4U,  8U, 16U },
                                                               { 2U, 4U,  8U, 16U, 32U },
                                                               { 4U, 8U, 16U, 32U, 64U } };

         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> B{ { 1U, 2U, 1U, 2U, 1U },
                                                               { 2U, 1U, 2U, 1U, 2U },
                                                               { 1U, 2U, 1U, 2U, 1U } };

         A <<= B;

         checkRows    ( A,  3UL );
         checkColumns ( A,  5UL );
         checkCapacity( A, 15UL );
         checkNonZeros( A, 15UL );

         if( A(0,0) != 2U || A(0,1) !=  8U || A(0,2) !=  8U || A(0,3) !=  32U || A(0,4) !=  32U ||
             A(1,0) != 8U || A(1,1) !=  8U || A(1,2) != 32U || A(1,3) !=  32U || A(1,4) != 128U ||
             A(2,0) != 8U || A(2,1) != 32U || A(2,2) != 32U || A(2,3) != 128U || A(2,4) != 128U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/matrix left-shift assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << A << "\n"
                << "   Expected result:\n( 2  8  8  32  32 )\n"
                                        "( 8  8 32  32 128 )\n"
                                        "( 8 32 32 128 128 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Row-major matrix/column-major matrix left-shift tests
   //=====================================================================================

   {
      test_ = "Row-major matrix/column-major matrix left-shift operator";

      // Matrix/matrix left-shift of an empty matrix
      {
         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> A;
         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> B;

         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> C( A << B );

         checkRows    ( C, 0UL );
         checkColumns ( C, 0UL );
         checkCapacity( C, 0UL );
         checkNonZeros( C, 0UL );
      }

      // Matrix/matrix left-shift of a general matrix
      {
         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> A{ { 1U, 2U,  4U,  8U, 16U },
                                                               { 2U, 4U,  8U, 16U, 32U },
                                                               { 4U, 8U, 16U, 32U, 64U } };

         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> B{ { 1U, 2U, 1U, 2U, 1U },
                                                                  { 2U, 1U, 2U, 1U, 2U },
                                                                  { 1U, 2U, 1U, 2U, 1U } };

         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> C( A << B );

         checkRows    ( C,  3UL );
         checkColumns ( C,  5UL );
         checkCapacity( C, 15UL );
         checkNonZeros( C, 15UL );

         if( C(0,0) != 2U || C(0,1) !=  8U || C(0,2) !=  8U || C(0,3) !=  32U || C(0,4) !=  32U ||
             C(1,0) != 8U || C(1,1) !=  8U || C(1,2) != 32U || C(1,3) !=  32U || C(1,4) != 128U ||
             C(2,0) != 8U || C(2,1) != 32U || C(2,2) != 32U || C(2,3) != 128U || C(2,4) != 128U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/matrix left-shift operation failed\n"
                << " Details:\n"
                << "   Result:\n" << C << "\n"
                << "   Expected result:\n( 2  8  8  32  32 )\n"
                                        "( 8  8 32  32 128 )\n"
                                        "( 8 32 32 128 128 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Matrix/matrix left-shift assignment
      {
         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> A{ { 1U, 2U,  4U,  8U, 16U },
                                                               { 2U, 4U,  8U, 16U, 32U },
                                                               { 4U, 8U, 16U, 32U, 64U } };

         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> B{ { 1U, 2U, 1U, 2U, 1U },
                                                                  { 2U, 1U, 2U, 1U, 2U },
                                                                  { 1U, 2U, 1U, 2U, 1U } };

         A <<= B;

         checkRows    ( A,  3UL );
         checkColumns ( A,  5UL );
         checkCapacity( A, 15UL );
         checkNonZeros( A, 15UL );

         if( A(0,0) != 2U || A(0,1) !=  8U || A(0,2) !=  8U || A(0,3) !=  32U || A(0,4) !=  32U ||
             A(1,0) != 8U || A(1,1) !=  8U || A(1,2) != 32U || A(1,3) !=  32U || A(1,4) != 128U ||
             A(2,0) != 8U || A(2,1) != 32U || A(2,2) != 32U || A(2,3) != 128U || A(2,4) != 128U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/matrix left-shift assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << A << "\n"
                << "   Expected result:\n( 2  8  8  32  32 )\n"
                                        "( 8  8 32  32 128 )\n"
                                        "( 8 32 32 128 128 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix/scalar left-shift tests
   //=====================================================================================

   {
      test_ = "Column-major matrix/scalar left-shift operator";

      // Matrix/scalar left-shift of an empty matrix
      {
         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> A;

         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> B( A << 2U );

         checkRows    ( B, 0UL );
         checkColumns ( B, 0UL );
         checkCapacity( B, 0UL );
         checkNonZeros( B, 0UL );
      }

      // Matrix/scalar left-shift of a general matrix
      {
         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> A{ { 1U, 2U,  4U,  8U, 16U },
                                                                  { 2U, 4U,  8U, 16U, 32U },
                                                                  { 4U, 8U, 16U, 32U, 64U } };

         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> B( A << 2U );

         checkRows    ( B,  3UL );
         checkColumns ( B,  5UL );
         checkCapacity( B, 15UL );
         checkNonZeros( B, 15UL );

         if( B(0,0) !=  4U || B(0,1) !=  8U || B(0,2) != 16U || B(0,3) !=  32U || B(0,4) !=  64U ||
             B(1,0) !=  8U || B(1,1) != 16U || B(1,2) != 32U || B(1,3) !=  64U || B(1,4) != 128U ||
             B(2,0) != 16U || B(2,1) != 32U || B(2,2) != 64U || B(2,3) != 128U || B(2,4) != 256U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/scalar left-shift operation failed\n"
                << " Details:\n"
                << "   Result:\n" << B << "\n"
                << "   Expected result:\n(  4  8 16  32  64 )\n"
                                        "(  8 16 32  64 128 )\n"
                                        "( 16 32 64 128 256 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Matrix/scalar left-shift assignment
      {
         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> A{ { 1U, 2U,  4U,  8U, 16U },
                                                                  { 2U, 4U,  8U, 16U, 32U },
                                                                  { 4U, 8U, 16U, 32U, 64U } };

         A <<= 2U;

         checkRows    ( A,  3UL );
         checkColumns ( A,  5UL );
         checkCapacity( A, 15UL );
         checkNonZeros( A, 15UL );

         if( A(0,0) !=  4U || A(0,1) !=  8U || A(0,2) != 16U || A(0,3) !=  32U || A(0,4) !=  64U ||
             A(1,0) !=  8U || A(1,1) != 16U || A(1,2) != 32U || A(1,3) !=  64U || A(1,4) != 128U ||
             A(2,0) != 16U || A(2,1) != 32U || A(2,2) != 64U || A(2,3) != 128U || A(2,4) != 256U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/scalar left-shift assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << A << "\n"
                << "   Expected result:\n(  4  8 16  32  64 )\n"
                                        "(  8 16 32  64 128 )\n"
                                        "( 16 32 64 128 256 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix/row-major matrix left-shift tests
   //=====================================================================================

   {
      test_ = "Column-major matrix/row-major matrix left-shift operator";

      // Matrix/matrix left-shift of an empty matrix
      {
         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> A;
         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> B;

         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> C( A << B );

         checkRows    ( C, 0UL );
         checkColumns ( C, 0UL );
         checkCapacity( C, 0UL );
         checkNonZeros( C, 0UL );
      }

      // Matrix/matrix left-shift of a general matrix
      {
         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> A{ { 1U, 2U,  4U,  8U, 16U },
                                                                  { 2U, 4U,  8U, 16U, 32U },
                                                                  { 4U, 8U, 16U, 32U, 64U } };

         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> B{ { 1U, 2U, 1U, 2U, 1U },
                                                               { 2U, 1U, 2U, 1U, 2U },
                                                               { 1U, 2U, 1U, 2U, 1U } };

         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> C( A << B );

         checkRows    ( C,  3UL );
         checkColumns ( C,  5UL );
         checkCapacity( C, 15UL );
         checkNonZeros( C, 15UL );

         if( C(0,0) != 2U || C(0,1) !=  8U || C(0,2) !=  8U || C(0,3) !=  32U || C(0,4) !=  32U ||
             C(1,0) != 8U || C(1,1) !=  8U || C(1,2) != 32U || C(1,3) !=  32U || C(1,4) != 128U ||
             C(2,0) != 8U || C(2,1) != 32U || C(2,2) != 32U || C(2,3) != 128U || C(2,4) != 128U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/matrix left-shift operation failed\n"
                << " Details:\n"
                << "   Result:\n" << C << "\n"
                << "   Expected result:\n( 2  8  8  32  32 )\n"
                                        "( 8  8 32  32 128 )\n"
                                        "( 8 32 32 128 128 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Matrix/matrix left-shift assignment
      {
         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> A{ { 1U, 2U,  4U,  8U, 16U },
                                                                  { 2U, 4U,  8U, 16U, 32U },
                                                                  { 4U, 8U, 16U, 32U, 64U } };

         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> B{ { 1U, 2U, 1U, 2U, 1U },
                                                               { 2U, 1U, 2U, 1U, 2U },
                                                               { 1U, 2U, 1U, 2U, 1U } };

         A <<= B;

         checkRows    ( A,  3UL );
         checkColumns ( A,  5UL );
         checkCapacity( A, 15UL );
         checkNonZeros( A, 15UL );

         if( A(0,0) != 2U || A(0,1) !=  8U || A(0,2) !=  8U || A(0,3) !=  32U || A(0,4) !=  32U ||
             A(1,0) != 8U || A(1,1) !=  8U || A(1,2) != 32U || A(1,3) !=  32U || A(1,4) != 128U ||
             A(2,0) != 8U || A(2,1) != 32U || A(2,2) != 32U || A(2,3) != 128U || A(2,4) != 128U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/matrix left-shift assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << A << "\n"
                << "   Expected result:\n( 2  8  8  32  32 )\n"
                                        "( 8  8 32  32 128 )\n"
                                        "( 8 32 32 128 128 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix/column-major matrix left-shift tests
   //=====================================================================================

   {
      test_ = "Column-major matrix/column-major matrix left-shift operator";

      // Matrix/matrix left-shift of an empty matrix
      {
         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> A;
         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> B;

         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> C( A << B );

         checkRows    ( C, 0UL );
         checkColumns ( C, 0UL );
         checkCapacity( C, 0UL );
         checkNonZeros( C, 0UL );
      }

      // Matrix/matrix left-shift of a general matrix
      {
         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> A{ { 1U, 2U,  4U,  8U, 16U },
                                                                  { 2U, 4U,  8U, 16U, 32U },
                                                                  { 4U, 8U, 16U, 32U, 64U } };

         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> B{ { 1U, 2U, 1U, 2U, 1U },
                                                                  { 2U, 1U, 2U, 1U, 2U },
                                                                  { 1U, 2U, 1U, 2U, 1U } };

         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> C( A << B );

         checkRows    ( C,  3UL );
         checkColumns ( C,  5UL );
         checkCapacity( C, 15UL );
         checkNonZeros( C, 15UL );

         if( C(0,0) != 2U || C(0,1) !=  8U || C(0,2) !=  8U || C(0,3) !=  32U || C(0,4) !=  32U ||
             C(1,0) != 8U || C(1,1) !=  8U || C(1,2) != 32U || C(1,3) !=  32U || C(1,4) != 128U ||
             C(2,0) != 8U || C(2,1) != 32U || C(2,2) != 32U || C(2,3) != 128U || C(2,4) != 128U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/matrix left-shift operation failed\n"
                << " Details:\n"
                << "   Result:\n" << C << "\n"
                << "   Expected result:\n( 2  8  8  32  32 )\n"
                                        "( 8  8 32  32 128 )\n"
                                        "( 8 32 32 128 128 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Matrix/matrix left-shift assignment
      {
         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> A{ { 1U, 2U,  4U,  8U, 16U },
                                                                  { 2U, 4U,  8U, 16U, 32U },
                                                                  { 4U, 8U, 16U, 32U, 64U } };

         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> B{ { 1U, 2U, 1U, 2U, 1U },
                                                                  { 2U, 1U, 2U, 1U, 2U },
                                                                  { 1U, 2U, 1U, 2U, 1U } };

         A <<= B;

         checkRows    ( A,  3UL );
         checkColumns ( A,  5UL );
         checkCapacity( A, 15UL );
         checkNonZeros( A, 15UL );

         if( A(0,0) != 2U || A(0,1) !=  8U || A(0,2) !=  8U || A(0,3) !=  32U || A(0,4) !=  32U ||
             A(1,0) != 8U || A(1,1) !=  8U || A(1,2) != 32U || A(1,3) !=  32U || A(1,4) != 128U ||
             A(2,0) != 8U || A(2,1) != 32U || A(2,2) != 32U || A(2,3) != 128U || A(2,4) != 128U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/matrix left-shift assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << A << "\n"
                << "   Expected result:\n( 2  8  8  32  32 )\n"
                                        "( 8  8 32  32 128 )\n"
                                        "( 8 32 32 128 128 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the right-shift operator for dense matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the right-shift operator for dense matrices. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testRightShift()
{
   //=====================================================================================
   // Row-major matrix/scalar right-shift tests
   //=====================================================================================

   {
      test_ = "Row-major matrix/scalar right-shift operator";

      // Matrix/scalar right-shift of an empty matrix
      {
         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> A;

         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> B( A >> 2U );

         checkRows    ( B, 0UL );
         checkColumns ( B, 0UL );
         checkCapacity( B, 0UL );
         checkNonZeros( B, 0UL );
      }

      // Matrix/scalar right-shift of a general matrix
      {
         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> A{ {  4U,  8U, 16U,  32U,  64U },
                                                               {  8U, 16U, 32U,  64U, 128U },
                                                               { 16U, 32U, 64U, 128U, 256U } };

         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> B( A >> 2U );

         checkRows    ( B,  3UL );
         checkColumns ( B,  5UL );
         checkCapacity( B, 15UL );
         checkNonZeros( B, 15UL );

         if( B(0,0) != 1U || B(0,1) != 2U || B(0,2) !=  4U || B(0,3) !=  8U || B(0,4) != 16U ||
             B(1,0) != 2U || B(1,1) != 4U || B(1,2) !=  8U || B(1,3) != 16U || B(1,4) != 32U ||
             B(2,0) != 4U || B(2,1) != 8U || B(2,2) != 16U || B(2,3) != 32U || B(2,4) != 64U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/scalar right-shift operation failed\n"
                << " Details:\n"
                << "   Result:\n" << B << "\n"
                << "   Expected result:\n( 1  2  4  8 16 )\n"
                                        "( 2  4  8 16 32 )\n"
                                        "( 4  8 16 32 64 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Matrix/scalar right-shift assignment
      {
         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> A{ {  4U,  8U, 16U,  32U,  64U },
                                                               {  8U, 16U, 32U,  64U, 128U },
                                                               { 16U, 32U, 64U, 128U, 256U } };

         A >>= 2U;

         checkRows    ( A,  3UL );
         checkColumns ( A,  5UL );
         checkCapacity( A, 15UL );
         checkNonZeros( A, 15UL );

         if( A(0,0) != 1U || A(0,1) != 2U || A(0,2) !=  4U || A(0,3) !=  8U || A(0,4) != 16U ||
             A(1,0) != 2U || A(1,1) != 4U || A(1,2) !=  8U || A(1,3) != 16U || A(1,4) != 32U ||
             A(2,0) != 4U || A(2,1) != 8U || A(2,2) != 16U || A(2,3) != 32U || A(2,4) != 64U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/scalar right-shift assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << A << "\n"
                << "   Expected result:\n( 1  2  4  8 16 )\n"
                                        "( 2  4  8 16 32 )\n"
                                        "( 4  8 16 32 64 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Row-major matrix/row-major matrix right-shift tests
   //=====================================================================================

   {
      test_ = "Row-major matrix/row-major matrix right-shift operator";

      // Matrix/matrix right-shift of an empty matrix
      {
         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> A;
         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> B;

         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> C( A << B );

         checkRows    ( C, 0UL );
         checkColumns ( C, 0UL );
         checkCapacity( C, 0UL );
         checkNonZeros( C, 0UL );
      }

      // Matrix/matrix right-shift of a general matrix
      {
         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> A{ {  4U,  8U, 16U,  32U,  64U },
                                                               {  8U, 16U, 32U,  64U, 128U },
                                                               { 16U, 32U, 64U, 128U, 256U } };

         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> B{ { 1U, 2U, 1U, 2U, 1U },
                                                               { 2U, 1U, 2U, 1U, 2U },
                                                               { 1U, 2U, 1U, 2U, 1U } };

         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> C( A >> B );

         checkRows    ( C,  3UL );
         checkColumns ( C,  5UL );
         checkCapacity( C, 15UL );
         checkNonZeros( C, 15UL );

         if( C(0,0) != 2U || C(0,1) != 2U || C(0,2) !=  8U || C(0,3) !=  8U || C(0,4) !=  32U ||
             C(1,0) != 2U || C(1,1) != 8U || C(1,2) !=  8U || C(1,3) != 32U || C(1,4) !=  32U ||
             C(2,0) != 8U || C(2,1) != 8U || C(2,2) != 32U || C(2,3) != 32U || C(2,4) != 128U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/matrix right-shift operation failed\n"
                << " Details:\n"
                << "   Result:\n" << C << "\n"
                << "   Expected result:\n( 2 2  8  8  32 )\n"
                                        "( 2 8  8 32  32 )\n"
                                        "( 8 8 32 32 128 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Matrix/matrix right-shift assignment
      {
         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> A{ {  4U,  8U, 16U,  32U,  64U },
                                                               {  8U, 16U, 32U,  64U, 128U },
                                                               { 16U, 32U, 64U, 128U, 256U } };

         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> B{ { 1U, 2U, 1U, 2U, 1U },
                                                               { 2U, 1U, 2U, 1U, 2U },
                                                               { 1U, 2U, 1U, 2U, 1U } };

         A >>= B;

         checkRows    ( A,  3UL );
         checkColumns ( A,  5UL );
         checkCapacity( A, 15UL );
         checkNonZeros( A, 15UL );

         if( A(0,0) != 2U || A(0,1) != 2U || A(0,2) !=  8U || A(0,3) !=  8U || A(0,4) !=  32U ||
             A(1,0) != 2U || A(1,1) != 8U || A(1,2) !=  8U || A(1,3) != 32U || A(1,4) !=  32U ||
             A(2,0) != 8U || A(2,1) != 8U || A(2,2) != 32U || A(2,3) != 32U || A(2,4) != 128U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/matrix right-shift assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << A << "\n"
                << "   Expected result:\n( 2 2  8  8  32 )\n"
                                        "( 2 8  8 32  32 )\n"
                                        "( 8 8 32 32 128 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Row-major matrix/column-major matrix right-shift tests
   //=====================================================================================

   {
      test_ = "Row-major matrix/column-major matrix right-shift operator";

      // Matrix/matrix right-shift of an empty matrix
      {
         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> A;
         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> B;

         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> C( A << B );

         checkRows    ( C, 0UL );
         checkColumns ( C, 0UL );
         checkCapacity( C, 0UL );
         checkNonZeros( C, 0UL );
      }

      // Matrix/matrix right-shift of a general matrix
      {
         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> A{ {  4U,  8U, 16U,  32U,  64U },
                                                               {  8U, 16U, 32U,  64U, 128U },
                                                               { 16U, 32U, 64U, 128U, 256U } };

         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> B{ { 1U, 2U, 1U, 2U, 1U },
                                                                  { 2U, 1U, 2U, 1U, 2U },
                                                                  { 1U, 2U, 1U, 2U, 1U } };

         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> C( A >> B );

         checkRows    ( C,  3UL );
         checkColumns ( C,  5UL );
         checkCapacity( C, 15UL );
         checkNonZeros( C, 15UL );

         if( C(0,0) != 2U || C(0,1) != 2U || C(0,2) !=  8U || C(0,3) !=  8U || C(0,4) !=  32U ||
             C(1,0) != 2U || C(1,1) != 8U || C(1,2) !=  8U || C(1,3) != 32U || C(1,4) !=  32U ||
             C(2,0) != 8U || C(2,1) != 8U || C(2,2) != 32U || C(2,3) != 32U || C(2,4) != 128U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/matrix right-shift operation failed\n"
                << " Details:\n"
                << "   Result:\n" << C << "\n"
                << "   Expected result:\n( 2 2  8  8  32 )\n"
                                        "( 2 8  8 32  32 )\n"
                                        "( 8 8 32 32 128 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Matrix/matrix right-shift assignment
      {
         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> A{ {  4U,  8U, 16U,  32U,  64U },
                                                               {  8U, 16U, 32U,  64U, 128U },
                                                               { 16U, 32U, 64U, 128U, 256U } };

         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> B{ { 1U, 2U, 1U, 2U, 1U },
                                                                  { 2U, 1U, 2U, 1U, 2U },
                                                                  { 1U, 2U, 1U, 2U, 1U } };

         A >>= B;

         checkRows    ( A,  3UL );
         checkColumns ( A,  5UL );
         checkCapacity( A, 15UL );
         checkNonZeros( A, 15UL );

         if( A(0,0) != 2U || A(0,1) != 2U || A(0,2) !=  8U || A(0,3) !=  8U || A(0,4) !=  32U ||
             A(1,0) != 2U || A(1,1) != 8U || A(1,2) !=  8U || A(1,3) != 32U || A(1,4) !=  32U ||
             A(2,0) != 8U || A(2,1) != 8U || A(2,2) != 32U || A(2,3) != 32U || A(2,4) != 128U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/matrix right-shift assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << A << "\n"
                << "   Expected result:\n( 2 2  8  8  32 )\n"
                                        "( 2 8  8 32  32 )\n"
                                        "( 8 8 32 32 128 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix/scalar right-shift tests
   //=====================================================================================

   {
      test_ = "Column-major matrix/scalar right-shift operator";

      // Matrix/scalar right-shift of an empty matrix
      {
         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> A;

         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> B( A >> 2U );

         checkRows    ( B, 0UL );
         checkColumns ( B, 0UL );
         checkCapacity( B, 0UL );
         checkNonZeros( B, 0UL );
      }

      // Matrix/scalar right-shift of a general matrix
      {
         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> A{ {  4U,  8U, 16U,  32U,  64U },
                                                                  {  8U, 16U, 32U,  64U, 128U },
                                                                  { 16U, 32U, 64U, 128U, 256U } };

         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> B( A >> 2U );

         checkRows    ( B,  3UL );
         checkColumns ( B,  5UL );
         checkCapacity( B, 15UL );
         checkNonZeros( B, 15UL );

         if( B(0,0) != 1U || B(0,1) != 2U || B(0,2) !=  4U || B(0,3) !=  8U || B(0,4) != 16U ||
             B(1,0) != 2U || B(1,1) != 4U || B(1,2) !=  8U || B(1,3) != 16U || B(1,4) != 32U ||
             B(2,0) != 4U || B(2,1) != 8U || B(2,2) != 16U || B(2,3) != 32U || B(2,4) != 64U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/scalar right-shift operation failed\n"
                << " Details:\n"
                << "   Result:\n" << B << "\n"
                << "   Expected result:\n( 1  2  4  8 16 )\n"
                                        "( 2  4  8 16 32 )\n"
                                        "( 4  8 16 32 64 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Matrix/scalar right-shift assignment
      {
         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> A{ {  4U,  8U, 16U,  32U,  64U },
                                                                  {  8U, 16U, 32U,  64U, 128U },
                                                                  { 16U, 32U, 64U, 128U, 256U } };

         A >>= 2U;

         checkRows    ( A,  3UL );
         checkColumns ( A,  5UL );
         checkCapacity( A, 15UL );
         checkNonZeros( A, 15UL );

         if( A(0,0) != 1U || A(0,1) != 2U || A(0,2) !=  4U || A(0,3) !=  8U || A(0,4) != 16U ||
             A(1,0) != 2U || A(1,1) != 4U || A(1,2) !=  8U || A(1,3) != 16U || A(1,4) != 32U ||
             A(2,0) != 4U || A(2,1) != 8U || A(2,2) != 16U || A(2,3) != 32U || A(2,4) != 64U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/scalar right-shift assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << A << "\n"
                << "   Expected result:\n( 1  2  4  8 16 )\n"
                                        "( 2  4  8 16 32 )\n"
                                        "( 4  8 16 32 64 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix/row-major matrix right-shift tests
   //=====================================================================================

   {
      test_ = "Column-major matrix/row-major matrix right-shift operator";

      // Matrix/matrix right-shift of an empty matrix
      {
         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> A;
         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> B;

         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> C( A << B );

         checkRows    ( C, 0UL );
         checkColumns ( C, 0UL );
         checkCapacity( C, 0UL );
         checkNonZeros( C, 0UL );
      }

      // Matrix/matrix right-shift of a general matrix
      {
         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> A{ {  4U,  8U, 16U,  32U,  64U },
                                                                  {  8U, 16U, 32U,  64U, 128U },
                                                                  { 16U, 32U, 64U, 128U, 256U } };

         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> B{ { 1U, 2U, 1U, 2U, 1U },
                                                               { 2U, 1U, 2U, 1U, 2U },
                                                               { 1U, 2U, 1U, 2U, 1U } };

         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> C( A >> B );

         checkRows    ( C,  3UL );
         checkColumns ( C,  5UL );
         checkCapacity( C, 15UL );
         checkNonZeros( C, 15UL );

         if( C(0,0) != 2U || C(0,1) != 2U || C(0,2) !=  8U || C(0,3) !=  8U || C(0,4) !=  32U ||
             C(1,0) != 2U || C(1,1) != 8U || C(1,2) !=  8U || C(1,3) != 32U || C(1,4) !=  32U ||
             C(2,0) != 8U || C(2,1) != 8U || C(2,2) != 32U || C(2,3) != 32U || C(2,4) != 128U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/matrix right-shift operation failed\n"
                << " Details:\n"
                << "   Result:\n" << C << "\n"
                << "   Expected result:\n( 2 2  8  8  32 )\n"
                                        "( 2 8  8 32  32 )\n"
                                        "( 8 8 32 32 128 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Matrix/matrix right-shift assignment
      {
         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> A{ {  4U,  8U, 16U,  32U,  64U },
                                                                  {  8U, 16U, 32U,  64U, 128U },
                                                                  { 16U, 32U, 64U, 128U, 256U } };

         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> B{ { 1U, 2U, 1U, 2U, 1U },
                                                               { 2U, 1U, 2U, 1U, 2U },
                                                               { 1U, 2U, 1U, 2U, 1U } };

         A >>= B;

         checkRows    ( A,  3UL );
         checkColumns ( A,  5UL );
         checkCapacity( A, 15UL );
         checkNonZeros( A, 15UL );

         if( A(0,0) != 2U || A(0,1) != 2U || A(0,2) !=  8U || A(0,3) !=  8U || A(0,4) !=  32U ||
             A(1,0) != 2U || A(1,1) != 8U || A(1,2) !=  8U || A(1,3) != 32U || A(1,4) !=  32U ||
             A(2,0) != 8U || A(2,1) != 8U || A(2,2) != 32U || A(2,3) != 32U || A(2,4) != 128U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/matrix right-shift assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << A << "\n"
                << "   Expected result:\n( 2 2  8  8  32 )\n"
                                        "( 2 8  8 32  32 )\n"
                                        "( 8 8 32 32 128 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix/column-major matrix right-shift tests
   //=====================================================================================

   {
      test_ = "Column-major matrix/column-major matrix right-shift operator";

      // Matrix/matrix right-shift of an empty matrix
      {
         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> A;
         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> B;

         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> C( A << B );

         checkRows    ( C, 0UL );
         checkColumns ( C, 0UL );
         checkCapacity( C, 0UL );
         checkNonZeros( C, 0UL );
      }

      // Matrix/matrix right-shift of a general matrix
      {
         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> A{ {  4U,  8U, 16U,  32U,  64U },
                                                                  {  8U, 16U, 32U,  64U, 128U },
                                                                  { 16U, 32U, 64U, 128U, 256U } };

         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> B{ { 1U, 2U, 1U, 2U, 1U },
                                                                  { 2U, 1U, 2U, 1U, 2U },
                                                                  { 1U, 2U, 1U, 2U, 1U } };

         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> C( A >> B );

         checkRows    ( C,  3UL );
         checkColumns ( C,  5UL );
         checkCapacity( C, 15UL );
         checkNonZeros( C, 15UL );

         if( C(0,0) != 2U || C(0,1) != 2U || C(0,2) !=  8U || C(0,3) !=  8U || C(0,4) !=  32U ||
             C(1,0) != 2U || C(1,1) != 8U || C(1,2) !=  8U || C(1,3) != 32U || C(1,4) !=  32U ||
             C(2,0) != 8U || C(2,1) != 8U || C(2,2) != 32U || C(2,3) != 32U || C(2,4) != 128U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/matrix right-shift operation failed\n"
                << " Details:\n"
                << "   Result:\n" << C << "\n"
                << "   Expected result:\n( 2 2  8  8  32 )\n"
                                        "( 2 8  8 32  32 )\n"
                                        "( 8 8 32 32 128 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Matrix/matrix right-shift assignment
      {
         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> A{ {  4U,  8U, 16U,  32U,  64U },
                                                                  {  8U, 16U, 32U,  64U, 128U },
                                                                  { 16U, 32U, 64U, 128U, 256U } };

         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> B{ { 1U, 2U, 1U, 2U, 1U },
                                                                  { 2U, 1U, 2U, 1U, 2U },
                                                                  { 1U, 2U, 1U, 2U, 1U } };

         A >>= B;

         checkRows    ( A,  3UL );
         checkColumns ( A,  5UL );
         checkCapacity( A, 15UL );
         checkNonZeros( A, 15UL );

         if( A(0,0) != 2U || A(0,1) != 2U || A(0,2) !=  8U || A(0,3) !=  8U || A(0,4) !=  32U ||
             A(1,0) != 2U || A(1,1) != 8U || A(1,2) !=  8U || A(1,3) != 32U || A(1,4) !=  32U ||
             A(2,0) != 8U || A(2,1) != 8U || A(2,2) != 32U || A(2,3) != 32U || A(2,4) != 128U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/matrix right-shift assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << A << "\n"
                << "   Expected result:\n( 2 2  8  8  32 )\n"
                                        "( 2 8  8 32  32 )\n"
                                        "( 8 8 32 32 128 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the bitwise AND operator for dense matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the bitwise AND operator for dense matrices. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testBitand()
{
   //=====================================================================================
   // Row-major matrix/scalar bitwise AND tests
   //=====================================================================================

   {
      test_ = "Row-major matrix/scalar bitwise AND operator";

      // Matrix/scalar bitwise AND of an empty matrix
      {
         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> A;

         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> B( A & 7U );

         checkRows    ( B, 0UL );
         checkColumns ( B, 0UL );
         checkCapacity( B, 0UL );
         checkNonZeros( B, 0UL );
      }

      // Matrix/scalar bitwise AND of a general matrix
      {
         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> A{ {  8U,  9U, 10U, 11U, 12U },
                                                               { 13U, 14U, 15U, 16U, 17U },
                                                               { 18U, 19U, 20U, 21U, 22U } };

         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> B( A & 7U );

         checkRows    ( B,  3UL );
         checkColumns ( B,  5UL );
         checkCapacity( B, 15UL );
         checkNonZeros( B, 13UL );

         if( B(0,0) != 0U || B(0,1) != 1U || B(0,2) != 2U || B(0,3) != 3U || B(0,4) != 4U ||
             B(1,0) != 5U || B(1,1) != 6U || B(1,2) != 7U || B(1,3) != 0U || B(1,4) != 1U ||
             B(2,0) != 2U || B(2,1) != 3U || B(2,2) != 4U || B(2,3) != 5U || B(2,4) != 6U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/scalar bitwise AND operation failed\n"
                << " Details:\n"
                << "   Result:\n" << B << "\n"
                << "   Expected result:\n( 0 1 2 3 4 )\n"
                                        "( 5 6 7 0 1 )\n"
                                        "( 2 3 4 5 6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Matrix/scalar bitwise AND assignment
      {
         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> A{ {  8U,  9U, 10U, 11U, 12U },
                                                               { 13U, 14U, 15U, 16U, 17U },
                                                               { 18U, 19U, 20U, 21U, 22U } };

         A &= 7U;

         checkRows    ( A,  3UL );
         checkColumns ( A,  5UL );
         checkCapacity( A, 15UL );
         checkNonZeros( A, 13UL );

         if( A(0,0) != 0U || A(0,1) != 1U || A(0,2) != 2U || A(0,3) != 3U || A(0,4) != 4U ||
             A(1,0) != 5U || A(1,1) != 6U || A(1,2) != 7U || A(1,3) != 0U || A(1,4) != 1U ||
             A(2,0) != 2U || A(2,1) != 3U || A(2,2) != 4U || A(2,3) != 5U || A(2,4) != 6U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/scalar bitwise AND assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << A << "\n"
                << "   Expected result:\n( 0 1 2 3 4 )\n"
                                        "( 5 6 7 0 1 )\n"
                                        "( 2 3 4 5 6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Row-major matrix/row-major matrix bitwise AND tests
   //=====================================================================================

   {
      test_ = "Row-major matrix/row-major matrix bitwise AND operator";

      // Matrix/matrix bitwise AND of an empty matrix
      {
         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> A;
         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> B;

         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> C( A & B );

         checkRows    ( C, 0UL );
         checkColumns ( C, 0UL );
         checkCapacity( C, 0UL );
         checkNonZeros( C, 0UL );
      }

      // Matrix/matrix bitwise AND of a general matrix
      {
         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> A{ {  8U,  9U, 10U, 11U, 12U },
                                                               { 13U, 14U, 15U, 16U, 17U },
                                                               { 18U, 19U, 20U, 21U, 22U } };

         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> B{ { 7U, 5U, 7U, 5U, 7U },
                                                               { 5U, 7U, 5U, 7U, 5U },
                                                               { 7U, 5U, 7U, 5U, 7U } };

         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> C( A & B );

         checkRows    ( C,  3UL );
         checkColumns ( C,  5UL );
         checkCapacity( C, 15UL );
         checkNonZeros( C, 13UL );

         if( C(0,0) != 0U || C(0,1) != 1U || C(0,2) != 2U || C(0,3) != 1U || C(0,4) != 4U ||
             C(1,0) != 5U || C(1,1) != 6U || C(1,2) != 5U || C(1,3) != 0U || C(1,4) != 1U ||
             C(2,0) != 2U || C(2,1) != 1U || C(2,2) != 4U || C(2,3) != 5U || C(2,4) != 6U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/matrix bitwise AND operation failed\n"
                << " Details:\n"
                << "   Result:\n" << C << "\n"
                << "   Expected result:\n( 0 1 2 1 4 )\n"
                                        "( 5 6 5 0 1 )\n"
                                        "( 2 1 4 5 6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Matrix/matrix bitwise AND assignment
      {
         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> A{ {  8U,  9U, 10U, 11U, 12U },
                                                               { 13U, 14U, 15U, 16U, 17U },
                                                               { 18U, 19U, 20U, 21U, 22U } };

         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> B{ { 7U, 5U, 7U, 5U, 7U },
                                                               { 5U, 7U, 5U, 7U, 5U },
                                                               { 7U, 5U, 7U, 5U, 7U } };

         A &= B;

         checkRows    ( A,  3UL );
         checkColumns ( A,  5UL );
         checkCapacity( A, 15UL );
         checkNonZeros( A, 13UL );

         if( A(0,0) != 0U || A(0,1) != 1U || A(0,2) != 2U || A(0,3) != 1U || A(0,4) != 4U ||
             A(1,0) != 5U || A(1,1) != 6U || A(1,2) != 5U || A(1,3) != 0U || A(1,4) != 1U ||
             A(2,0) != 2U || A(2,1) != 1U || A(2,2) != 4U || A(2,3) != 5U || A(2,4) != 6U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/matrix bitwise AND assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << A << "\n"
                << "   Expected result:\n( 0 1 2 1 4 )\n"
                                        "( 5 6 5 0 1 )\n"
                                        "( 2 1 4 5 6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Row-major matrix/column-major matrix bitwise AND tests
   //=====================================================================================

   {
      test_ = "Row-major matrix/column-major matrix bitwise AND operator";

      // Matrix/matrix bitwise AND of an empty matrix
      {
         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> A;
         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> B;

         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> C( A & B );

         checkRows    ( C, 0UL );
         checkColumns ( C, 0UL );
         checkCapacity( C, 0UL );
         checkNonZeros( C, 0UL );
      }

      // Matrix/matrix bitwise AND of a general matrix
      {
         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> A{ {  8U,  9U, 10U, 11U, 12U },
                                                               { 13U, 14U, 15U, 16U, 17U },
                                                               { 18U, 19U, 20U, 21U, 22U } };

         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> B{ { 7U, 5U, 7U, 5U, 7U },
                                                                  { 5U, 7U, 5U, 7U, 5U },
                                                                  { 7U, 5U, 7U, 5U, 7U } };

         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> C( A & B );

         checkRows    ( C,  3UL );
         checkColumns ( C,  5UL );
         checkCapacity( C, 15UL );
         checkNonZeros( C, 13UL );

         if( C(0,0) != 0U || C(0,1) != 1U || C(0,2) != 2U || C(0,3) != 1U || C(0,4) != 4U ||
             C(1,0) != 5U || C(1,1) != 6U || C(1,2) != 5U || C(1,3) != 0U || C(1,4) != 1U ||
             C(2,0) != 2U || C(2,1) != 1U || C(2,2) != 4U || C(2,3) != 5U || C(2,4) != 6U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/matrix bitwise AND operation failed\n"
                << " Details:\n"
                << "   Result:\n" << C << "\n"
                << "   Expected result:\n( 0 1 2 1 4 )\n"
                                        "( 5 6 5 0 1 )\n"
                                        "( 2 1 4 5 6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Matrix/matrix bitwise AND assignment
      {
         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> A{ {  8U,  9U, 10U, 11U, 12U },
                                                               { 13U, 14U, 15U, 16U, 17U },
                                                               { 18U, 19U, 20U, 21U, 22U } };

         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> B{ { 7U, 5U, 7U, 5U, 7U },
                                                                  { 5U, 7U, 5U, 7U, 5U },
                                                                  { 7U, 5U, 7U, 5U, 7U } };

         A &= B;

         checkRows    ( A,  3UL );
         checkColumns ( A,  5UL );
         checkCapacity( A, 15UL );
         checkNonZeros( A, 13UL );

         if( A(0,0) != 0U || A(0,1) != 1U || A(0,2) != 2U || A(0,3) != 1U || A(0,4) != 4U ||
             A(1,0) != 5U || A(1,1) != 6U || A(1,2) != 5U || A(1,3) != 0U || A(1,4) != 1U ||
             A(2,0) != 2U || A(2,1) != 1U || A(2,2) != 4U || A(2,3) != 5U || A(2,4) != 6U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/matrix bitwise AND assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << A << "\n"
                << "   Expected result:\n( 0 1 2 1 4 )\n"
                                        "( 5 6 5 0 1 )\n"
                                        "( 2 1 4 5 6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix/scalar bitwise AND tests
   //=====================================================================================

   {
      test_ = "Row-major matrix/scalar bitwise AND operator";

      // Matrix/scalar bitwise AND of an empty matrix
      {
         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> A;

         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> B( A & 7U );

         checkRows    ( B, 0UL );
         checkColumns ( B, 0UL );
         checkCapacity( B, 0UL );
         checkNonZeros( B, 0UL );
      }

      // Matrix/scalar bitwise AND of a general matrix
      {
         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> A{ {  8U,  9U, 10U, 11U, 12U },
                                                                  { 13U, 14U, 15U, 16U, 17U },
                                                                  { 18U, 19U, 20U, 21U, 22U } };

         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> B( A & 7U );

         checkRows    ( B,  3UL );
         checkColumns ( B,  5UL );
         checkCapacity( B, 15UL );
         checkNonZeros( B, 13UL );

         if( B(0,0) != 0U || B(0,1) != 1U || B(0,2) != 2U || B(0,3) != 3U || B(0,4) != 4U ||
             B(1,0) != 5U || B(1,1) != 6U || B(1,2) != 7U || B(1,3) != 0U || B(1,4) != 1U ||
             B(2,0) != 2U || B(2,1) != 3U || B(2,2) != 4U || B(2,3) != 5U || B(2,4) != 6U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/scalar bitwise AND operation failed\n"
                << " Details:\n"
                << "   Result:\n" << B << "\n"
                << "   Expected result:\n( 0 1 2 3 4 )\n"
                                        "( 5 6 7 0 1 )\n"
                                        "( 2 3 4 5 6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Matrix/scalar bitwise AND assignment
      {
         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> A{ {  8U,  9U, 10U, 11U, 12U },
                                                                  { 13U, 14U, 15U, 16U, 17U },
                                                                  { 18U, 19U, 20U, 21U, 22U } };

         A &= 7U;

         checkRows    ( A,  3UL );
         checkColumns ( A,  5UL );
         checkCapacity( A, 15UL );
         checkNonZeros( A, 13UL );

         if( A(0,0) != 0U || A(0,1) != 1U || A(0,2) != 2U || A(0,3) != 3U || A(0,4) != 4U ||
             A(1,0) != 5U || A(1,1) != 6U || A(1,2) != 7U || A(1,3) != 0U || A(1,4) != 1U ||
             A(2,0) != 2U || A(2,1) != 3U || A(2,2) != 4U || A(2,3) != 5U || A(2,4) != 6U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/scalar bitwise AND assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << A << "\n"
                << "   Expected result:\n( 0 1 2 3 4 )\n"
                                        "( 5 6 7 0 1 )\n"
                                        "( 2 3 4 5 6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix/row-major matrix bitwise AND tests
   //=====================================================================================

   {
      test_ = "Column-major matrix/row-major matrix bitwise AND operator";

      // Matrix/matrix bitwise AND of an empty matrix
      {
         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> A;
         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> B;

         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> C( A & B );

         checkRows    ( C, 0UL );
         checkColumns ( C, 0UL );
         checkCapacity( C, 0UL );
         checkNonZeros( C, 0UL );
      }

      // Matrix/matrix bitwise AND of a general matrix
      {
         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> A{ {  8U,  9U, 10U, 11U, 12U },
                                                                  { 13U, 14U, 15U, 16U, 17U },
                                                                  { 18U, 19U, 20U, 21U, 22U } };

         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> B{ { 7U, 5U, 7U, 5U, 7U },
                                                               { 5U, 7U, 5U, 7U, 5U },
                                                               { 7U, 5U, 7U, 5U, 7U } };

         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> C( A & B );

         checkRows    ( C,  3UL );
         checkColumns ( C,  5UL );
         checkCapacity( C, 15UL );
         checkNonZeros( C, 13UL );

         if( C(0,0) != 0U || C(0,1) != 1U || C(0,2) != 2U || C(0,3) != 1U || C(0,4) != 4U ||
             C(1,0) != 5U || C(1,1) != 6U || C(1,2) != 5U || C(1,3) != 0U || C(1,4) != 1U ||
             C(2,0) != 2U || C(2,1) != 1U || C(2,2) != 4U || C(2,3) != 5U || C(2,4) != 6U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/matrix bitwise AND operation failed\n"
                << " Details:\n"
                << "   Result:\n" << C << "\n"
                << "   Expected result:\n( 0 1 2 1 4 )\n"
                                        "( 5 6 5 0 1 )\n"
                                        "( 2 1 4 5 6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Matrix/matrix bitwise AND assignment
      {
         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> A{ {  8U,  9U, 10U, 11U, 12U },
                                                                  { 13U, 14U, 15U, 16U, 17U },
                                                                  { 18U, 19U, 20U, 21U, 22U } };

         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> B{ { 7U, 5U, 7U, 5U, 7U },
                                                               { 5U, 7U, 5U, 7U, 5U },
                                                               { 7U, 5U, 7U, 5U, 7U } };

         A &= B;

         checkRows    ( A,  3UL );
         checkColumns ( A,  5UL );
         checkCapacity( A, 15UL );
         checkNonZeros( A, 13UL );

         if( A(0,0) != 0U || A(0,1) != 1U || A(0,2) != 2U || A(0,3) != 1U || A(0,4) != 4U ||
             A(1,0) != 5U || A(1,1) != 6U || A(1,2) != 5U || A(1,3) != 0U || A(1,4) != 1U ||
             A(2,0) != 2U || A(2,1) != 1U || A(2,2) != 4U || A(2,3) != 5U || A(2,4) != 6U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/matrix bitwise AND assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << A << "\n"
                << "   Expected result:\n( 0 1 2 1 4 )\n"
                                        "( 5 6 5 0 1 )\n"
                                        "( 2 1 4 5 6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix/column-major matrix bitwise AND tests
   //=====================================================================================

   {
      test_ = "Column-major matrix/column-major matrix bitwise AND operator";

      // Matrix/matrix bitwise AND of an empty matrix
      {
         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> A;
         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> B;

         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> C( A & B );

         checkRows    ( C, 0UL );
         checkColumns ( C, 0UL );
         checkCapacity( C, 0UL );
         checkNonZeros( C, 0UL );
      }

      // Matrix/matrix bitwise AND of a general matrix
      {
         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> A{ {  8U,  9U, 10U, 11U, 12U },
                                                                  { 13U, 14U, 15U, 16U, 17U },
                                                                  { 18U, 19U, 20U, 21U, 22U } };

         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> B{ { 7U, 5U, 7U, 5U, 7U },
                                                                  { 5U, 7U, 5U, 7U, 5U },
                                                                  { 7U, 5U, 7U, 5U, 7U } };

         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> C( A & B );

         checkRows    ( C,  3UL );
         checkColumns ( C,  5UL );
         checkCapacity( C, 15UL );
         checkNonZeros( C, 13UL );

         if( C(0,0) != 0U || C(0,1) != 1U || C(0,2) != 2U || C(0,3) != 1U || C(0,4) != 4U ||
             C(1,0) != 5U || C(1,1) != 6U || C(1,2) != 5U || C(1,3) != 0U || C(1,4) != 1U ||
             C(2,0) != 2U || C(2,1) != 1U || C(2,2) != 4U || C(2,3) != 5U || C(2,4) != 6U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/matrix bitwise AND operation failed\n"
                << " Details:\n"
                << "   Result:\n" << C << "\n"
                << "   Expected result:\n( 0 1 2 1 4 )\n"
                                        "( 5 6 5 0 1 )\n"
                                        "( 2 1 4 5 6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Matrix/matrix bitwise AND assignment
      {
         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> A{ {  8U,  9U, 10U, 11U, 12U },
                                                                  { 13U, 14U, 15U, 16U, 17U },
                                                                  { 18U, 19U, 20U, 21U, 22U } };

         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> B{ { 7U, 5U, 7U, 5U, 7U },
                                                                  { 5U, 7U, 5U, 7U, 5U },
                                                                  { 7U, 5U, 7U, 5U, 7U } };

         A &= B;

         checkRows    ( A,  3UL );
         checkColumns ( A,  5UL );
         checkCapacity( A, 15UL );
         checkNonZeros( A, 13UL );

         if( A(0,0) != 0U || A(0,1) != 1U || A(0,2) != 2U || A(0,3) != 1U || A(0,4) != 4U ||
             A(1,0) != 5U || A(1,1) != 6U || A(1,2) != 5U || A(1,3) != 0U || A(1,4) != 1U ||
             A(2,0) != 2U || A(2,1) != 1U || A(2,2) != 4U || A(2,3) != 5U || A(2,4) != 6U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/matrix bitwise AND assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << A << "\n"
                << "   Expected result:\n( 0 1 2 1 4 )\n"
                                        "( 5 6 5 0 1 )\n"
                                        "( 2 1 4 5 6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the bitwise OR operator for dense matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the bitwise OR operator for dense matrices. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testBitor()
{
   //=====================================================================================
   // Row-major matrix/scalar bitwise OR tests
   //=====================================================================================

   {
      test_ = "Row-major matrix/scalar bitwise OR operator";

      // Matrix/scalar bitwise OR of an empty matrix
      {
         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> A;

         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> B( A | 7U );

         checkRows    ( B, 0UL );
         checkColumns ( B, 0UL );
         checkCapacity( B, 0UL );
         checkNonZeros( B, 0UL );
      }

      // Matrix/scalar bitwise OR of a general matrix
      {
         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> A{ {  8U,  9U, 10U, 11U, 12U },
                                                               { 13U, 14U, 15U, 16U, 17U },
                                                               { 18U, 19U, 20U, 21U, 22U } };

         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> B( A | 7U );

         checkRows    ( B,  3UL );
         checkColumns ( B,  5UL );
         checkCapacity( B, 15UL );
         checkNonZeros( B, 15UL );

         if( B(0,0) != 15U || B(0,1) != 15U || B(0,2) != 15U || B(0,3) != 15U || B(0,4) != 15U ||
             B(1,0) != 15U || B(1,1) != 15U || B(1,2) != 15U || B(1,3) != 23U || B(1,4) != 23U ||
             B(2,0) != 23U || B(2,1) != 23U || B(2,2) != 23U || B(2,3) != 23U || B(2,4) != 23U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/scalar bitwise OR operation failed\n"
                << " Details:\n"
                << "   Result:\n" << B << "\n"
                << "   Expected result:\n( 15 15 15 15 15 )\n"
                                        "( 15 15 15 23 23 )\n"
                                        "( 23 23 23 23 23 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Matrix/scalar bitwise OR assignment
      {
         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> A{ {  8U,  9U, 10U, 11U, 12U },
                                                               { 13U, 14U, 15U, 16U, 17U },
                                                               { 18U, 19U, 20U, 21U, 22U } };

         A |= 7U;

         checkRows    ( A,  3UL );
         checkColumns ( A,  5UL );
         checkCapacity( A, 15UL );
         checkNonZeros( A, 15UL );

         if( A(0,0) != 15U || A(0,1) != 15U || A(0,2) != 15U || A(0,3) != 15U || A(0,4) != 15U ||
             A(1,0) != 15U || A(1,1) != 15U || A(1,2) != 15U || A(1,3) != 23U || A(1,4) != 23U ||
             A(2,0) != 23U || A(2,1) != 23U || A(2,2) != 23U || A(2,3) != 23U || A(2,4) != 23U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/scalar bitwise OR assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << A << "\n"
                << "   Expected result:\n( 15 15 15 15 15 )\n"
                                        "( 15 15 15 23 23 )\n"
                                        "( 23 23 23 23 23 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Row-major matrix/row-major matrix bitwise OR tests
   //=====================================================================================

   {
      test_ = "Row-major matrix/row-major matrix bitwise OR operator";

      // Matrix/matrix bitwise OR of an empty matrix
      {
         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> A;
         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> B;

         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> C( A | B );

         checkRows    ( C, 0UL );
         checkColumns ( C, 0UL );
         checkCapacity( C, 0UL );
         checkNonZeros( C, 0UL );
      }

      // Matrix/matrix bitwise OR of a general matrix
      {
         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> A{ {  8U,  9U, 10U, 11U, 12U },
                                                               { 13U, 14U, 15U, 16U, 17U },
                                                               { 18U, 19U, 20U, 21U, 22U } };

         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> B{ { 7U, 5U, 7U, 5U, 7U },
                                                               { 5U, 7U, 5U, 7U, 5U },
                                                               { 7U, 5U, 7U, 5U, 7U } };

         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> C( A | B );

         checkRows    ( C,  3UL );
         checkColumns ( C,  5UL );
         checkCapacity( C, 15UL );
         checkNonZeros( C, 15UL );

         if( C(0,0) != 15U || C(0,1) != 13U || C(0,2) != 15U || C(0,3) != 15U || C(0,4) != 15U ||
             C(1,0) != 13U || C(1,1) != 15U || C(1,2) != 15U || C(1,3) != 23U || C(1,4) != 21U ||
             C(2,0) != 23U || C(2,1) != 23U || C(2,2) != 23U || C(2,3) != 21U || C(2,4) != 23U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/matrix bitwise OR operation failed\n"
                << " Details:\n"
                << "   Result:\n" << C << "\n"
                << "   Expected result:\n( 15 13 15 15 15 )\n"
                                        "( 13 15 15 23 21 )\n"
                                        "( 23 23 23 21 23 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Matrix/matrix bitwise OR assignment
      {
         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> A{ {  8U,  9U, 10U, 11U, 12U },
                                                               { 13U, 14U, 15U, 16U, 17U },
                                                               { 18U, 19U, 20U, 21U, 22U } };

         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> B{ { 7U, 5U, 7U, 5U, 7U },
                                                               { 5U, 7U, 5U, 7U, 5U },
                                                               { 7U, 5U, 7U, 5U, 7U } };

         A |= B;

         checkRows    ( A,  3UL );
         checkColumns ( A,  5UL );
         checkCapacity( A, 15UL );
         checkNonZeros( A, 15UL );

         if( A(0,0) != 15U || A(0,1) != 13U || A(0,2) != 15U || A(0,3) != 15U || A(0,4) != 15U ||
             A(1,0) != 13U || A(1,1) != 15U || A(1,2) != 15U || A(1,3) != 23U || A(1,4) != 21U ||
             A(2,0) != 23U || A(2,1) != 23U || A(2,2) != 23U || A(2,3) != 21U || A(2,4) != 23U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/matrix bitwise OR assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << A << "\n"
                << "   Expected result:\n( 15 13 15 15 15 )\n"
                                        "( 13 15 15 23 21 )\n"
                                        "( 23 23 23 21 23 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Row-major matrix/column-major matrix bitwise OR tests
   //=====================================================================================

   {
      test_ = "Row-major matrix/column-major matrix bitwise OR operator";

      // Matrix/matrix bitwise OR of an empty matrix
      {
         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> A;
         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> B;

         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> C( A | B );

         checkRows    ( C, 0UL );
         checkColumns ( C, 0UL );
         checkCapacity( C, 0UL );
         checkNonZeros( C, 0UL );
      }

      // Matrix/matrix bitwise OR of a general matrix
      {
         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> A{ {  8U,  9U, 10U, 11U, 12U },
                                                               { 13U, 14U, 15U, 16U, 17U },
                                                               { 18U, 19U, 20U, 21U, 22U } };

         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> B{ { 7U, 5U, 7U, 5U, 7U },
                                                                  { 5U, 7U, 5U, 7U, 5U },
                                                                  { 7U, 5U, 7U, 5U, 7U } };

         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> C( A | B );

         checkRows    ( C,  3UL );
         checkColumns ( C,  5UL );
         checkCapacity( C, 15UL );
         checkNonZeros( C, 15UL );

         if( C(0,0) != 15U || C(0,1) != 13U || C(0,2) != 15U || C(0,3) != 15U || C(0,4) != 15U ||
             C(1,0) != 13U || C(1,1) != 15U || C(1,2) != 15U || C(1,3) != 23U || C(1,4) != 21U ||
             C(2,0) != 23U || C(2,1) != 23U || C(2,2) != 23U || C(2,3) != 21U || C(2,4) != 23U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/matrix bitwise OR operation failed\n"
                << " Details:\n"
                << "   Result:\n" << C << "\n"
                << "   Expected result:\n( 15 13 15 15 15 )\n"
                                        "( 13 15 15 23 21 )\n"
                                        "( 23 23 23 21 23 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Matrix/matrix bitwise OR assignment
      {
         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> A{ {  8U,  9U, 10U, 11U, 12U },
                                                               { 13U, 14U, 15U, 16U, 17U },
                                                               { 18U, 19U, 20U, 21U, 22U } };

         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> B{ { 7U, 5U, 7U, 5U, 7U },
                                                                  { 5U, 7U, 5U, 7U, 5U },
                                                                  { 7U, 5U, 7U, 5U, 7U } };

         A |= B;

         checkRows    ( A,  3UL );
         checkColumns ( A,  5UL );
         checkCapacity( A, 15UL );
         checkNonZeros( A, 15UL );

         if( A(0,0) != 15U || A(0,1) != 13U || A(0,2) != 15U || A(0,3) != 15U || A(0,4) != 15U ||
             A(1,0) != 13U || A(1,1) != 15U || A(1,2) != 15U || A(1,3) != 23U || A(1,4) != 21U ||
             A(2,0) != 23U || A(2,1) != 23U || A(2,2) != 23U || A(2,3) != 21U || A(2,4) != 23U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/matrix bitwise OR assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << A << "\n"
                << "   Expected result:\n( 15 13 15 15 15 )\n"
                                        "( 13 15 15 23 21 )\n"
                                        "( 23 23 23 21 23 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix/scalar bitwise OR tests
   //=====================================================================================

   {
      test_ = "Column-major matrix/scalar bitwise OR operator";

      // Matrix/scalar bitwise OR of an empty matrix
      {
         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> A;

         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> B( A | 7U );

         checkRows    ( B, 0UL );
         checkColumns ( B, 0UL );
         checkCapacity( B, 0UL );
         checkNonZeros( B, 0UL );
      }

      // Matrix/scalar bitwise OR of a general matrix
      {
         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> A{ {  8U,  9U, 10U, 11U, 12U },
                                                                  { 13U, 14U, 15U, 16U, 17U },
                                                                  { 18U, 19U, 20U, 21U, 22U } };

         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> B( A | 7U );

         checkRows    ( B,  3UL );
         checkColumns ( B,  5UL );
         checkCapacity( B, 15UL );
         checkNonZeros( B, 15UL );

         if( B(0,0) != 15U || B(0,1) != 15U || B(0,2) != 15U || B(0,3) != 15U || B(0,4) != 15U ||
             B(1,0) != 15U || B(1,1) != 15U || B(1,2) != 15U || B(1,3) != 23U || B(1,4) != 23U ||
             B(2,0) != 23U || B(2,1) != 23U || B(2,2) != 23U || B(2,3) != 23U || B(2,4) != 23U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/scalar bitwise OR operation failed\n"
                << " Details:\n"
                << "   Result:\n" << B << "\n"
                << "   Expected result:\n( 15 15 15 15 15 )\n"
                                        "( 15 15 15 23 23 )\n"
                                        "( 23 23 23 23 23 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Matrix/scalar bitwise OR assignment
      {
         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> A{ {  8U,  9U, 10U, 11U, 12U },
                                                                  { 13U, 14U, 15U, 16U, 17U },
                                                                  { 18U, 19U, 20U, 21U, 22U } };

         A |= 7U;

         checkRows    ( A,  3UL );
         checkColumns ( A,  5UL );
         checkCapacity( A, 15UL );
         checkNonZeros( A, 15UL );

         if( A(0,0) != 15U || A(0,1) != 15U || A(0,2) != 15U || A(0,3) != 15U || A(0,4) != 15U ||
             A(1,0) != 15U || A(1,1) != 15U || A(1,2) != 15U || A(1,3) != 23U || A(1,4) != 23U ||
             A(2,0) != 23U || A(2,1) != 23U || A(2,2) != 23U || A(2,3) != 23U || A(2,4) != 23U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/scalar bitwise OR assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << A << "\n"
                << "   Expected result:\n( 15 15 15 15 15 )\n"
                                        "( 15 15 15 23 23 )\n"
                                        "( 23 23 23 23 23 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix/row-major matrix bitwise OR tests
   //=====================================================================================


   {
      test_ = "Column-major matrix/row-major matrix bitwise OR operator";

      // Matrix/matrix bitwise OR of an empty matrix
      {
         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> A;
         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> B;

         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> C( A | B );

         checkRows    ( C, 0UL );
         checkColumns ( C, 0UL );
         checkCapacity( C, 0UL );
         checkNonZeros( C, 0UL );
      }

      // Matrix/matrix bitwise OR of a general matrix
      {
         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> A{ {  8U,  9U, 10U, 11U, 12U },
                                                                  { 13U, 14U, 15U, 16U, 17U },
                                                                  { 18U, 19U, 20U, 21U, 22U } };

         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> B{ { 7U, 5U, 7U, 5U, 7U },
                                                               { 5U, 7U, 5U, 7U, 5U },
                                                               { 7U, 5U, 7U, 5U, 7U } };

         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> C( A | B );

         checkRows    ( C,  3UL );
         checkColumns ( C,  5UL );
         checkCapacity( C, 15UL );
         checkNonZeros( C, 15UL );

         if( C(0,0) != 15U || C(0,1) != 13U || C(0,2) != 15U || C(0,3) != 15U || C(0,4) != 15U ||
             C(1,0) != 13U || C(1,1) != 15U || C(1,2) != 15U || C(1,3) != 23U || C(1,4) != 21U ||
             C(2,0) != 23U || C(2,1) != 23U || C(2,2) != 23U || C(2,3) != 21U || C(2,4) != 23U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/matrix bitwise OR operation failed\n"
                << " Details:\n"
                << "   Result:\n" << C << "\n"
                << "   Expected result:\n( 15 13 15 15 15 )\n"
                                        "( 13 15 15 23 21 )\n"
                                        "( 23 23 23 21 23 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Matrix/matrix bitwise OR assignment
      {
         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> A{ {  8U,  9U, 10U, 11U, 12U },
                                                                  { 13U, 14U, 15U, 16U, 17U },
                                                                  { 18U, 19U, 20U, 21U, 22U } };

         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> B{ { 7U, 5U, 7U, 5U, 7U },
                                                               { 5U, 7U, 5U, 7U, 5U },
                                                               { 7U, 5U, 7U, 5U, 7U } };

         A |= B;

         checkRows    ( A,  3UL );
         checkColumns ( A,  5UL );
         checkCapacity( A, 15UL );
         checkNonZeros( A, 15UL );

         if( A(0,0) != 15U || A(0,1) != 13U || A(0,2) != 15U || A(0,3) != 15U || A(0,4) != 15U ||
             A(1,0) != 13U || A(1,1) != 15U || A(1,2) != 15U || A(1,3) != 23U || A(1,4) != 21U ||
             A(2,0) != 23U || A(2,1) != 23U || A(2,2) != 23U || A(2,3) != 21U || A(2,4) != 23U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/matrix bitwise OR assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << A << "\n"
                << "   Expected result:\n( 15 13 15 15 15 )\n"
                                        "( 13 15 15 23 21 )\n"
                                        "( 23 23 23 21 23 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix/column-major matrix bitwise OR tests
   //=====================================================================================

   {
      test_ = "Column-major matrix/column-major matrix bitwise OR operator";

      // Matrix/matrix bitwise OR of an empty matrix
      {
         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> A;
         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> B;

         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> C( A | B );

         checkRows    ( C, 0UL );
         checkColumns ( C, 0UL );
         checkCapacity( C, 0UL );
         checkNonZeros( C, 0UL );
      }

      // Matrix/matrix bitwise OR of a general matrix
      {
         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> A{ {  8U,  9U, 10U, 11U, 12U },
                                                                  { 13U, 14U, 15U, 16U, 17U },
                                                                  { 18U, 19U, 20U, 21U, 22U } };

         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> B{ { 7U, 5U, 7U, 5U, 7U },
                                                                  { 5U, 7U, 5U, 7U, 5U },
                                                                  { 7U, 5U, 7U, 5U, 7U } };

         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> C( A | B );

         checkRows    ( C,  3UL );
         checkColumns ( C,  5UL );
         checkCapacity( C, 15UL );
         checkNonZeros( C, 15UL );

         if( C(0,0) != 15U || C(0,1) != 13U || C(0,2) != 15U || C(0,3) != 15U || C(0,4) != 15U ||
             C(1,0) != 13U || C(1,1) != 15U || C(1,2) != 15U || C(1,3) != 23U || C(1,4) != 21U ||
             C(2,0) != 23U || C(2,1) != 23U || C(2,2) != 23U || C(2,3) != 21U || C(2,4) != 23U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/matrix bitwise OR operation failed\n"
                << " Details:\n"
                << "   Result:\n" << C << "\n"
                << "   Expected result:\n( 15 13 15 15 15 )\n"
                                        "( 13 15 15 23 21 )\n"
                                        "( 23 23 23 21 23 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Matrix/matrix bitwise OR assignment
      {
         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> A{ {  8U,  9U, 10U, 11U, 12U },
                                                                  { 13U, 14U, 15U, 16U, 17U },
                                                                  { 18U, 19U, 20U, 21U, 22U } };

         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> B{ { 7U, 5U, 7U, 5U, 7U },
                                                                  { 5U, 7U, 5U, 7U, 5U },
                                                                  { 7U, 5U, 7U, 5U, 7U } };

         A |= B;

         checkRows    ( A,  3UL );
         checkColumns ( A,  5UL );
         checkCapacity( A, 15UL );
         checkNonZeros( A, 15UL );

         if( A(0,0) != 15U || A(0,1) != 13U || A(0,2) != 15U || A(0,3) != 15U || A(0,4) != 15U ||
             A(1,0) != 13U || A(1,1) != 15U || A(1,2) != 15U || A(1,3) != 23U || A(1,4) != 21U ||
             A(2,0) != 23U || A(2,1) != 23U || A(2,2) != 23U || A(2,3) != 21U || A(2,4) != 23U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/matrix bitwise OR assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << A << "\n"
                << "   Expected result:\n( 15 13 15 15 15 )\n"
                                        "( 13 15 15 23 21 )\n"
                                        "( 23 23 23 21 23 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the bitwise XOR operator for dense matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the bitwise XOR operator for dense matrices. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testBitxor()
{
   //=====================================================================================
   // Row-major matrix/scalar bitwise XOR tests
   //=====================================================================================

   {
      test_ = "Row-major matrix/scalar bitwise XOR operator";

      // Matrix/scalar bitwise XOR of an empty matrix
      {
         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> A;

         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> B( A ^ 7U );

         checkRows    ( B, 0UL );
         checkColumns ( B, 0UL );
         checkCapacity( B, 0UL );
         checkNonZeros( B, 0UL );
      }

      // Matrix/scalar bitwise XOR of a general matrix
      {
         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> A{ {  8U,  9U, 10U, 11U, 12U },
                                                               { 13U, 14U, 15U, 16U, 17U },
                                                               { 18U, 19U, 20U, 21U, 22U } };

         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> B( A ^ 7U );

         checkRows    ( B,  3UL );
         checkColumns ( B,  5UL );
         checkCapacity( B, 15UL );
         checkNonZeros( B, 15UL );

         if( B(0,0) != 15U || B(0,1) != 14U || B(0,2) != 13U || B(0,3) != 12U || B(0,4) != 11U ||
             B(1,0) != 10U || B(1,1) !=  9U || B(1,2) !=  8U || B(1,3) != 23U || B(1,4) != 22U ||
             B(2,0) != 21U || B(2,1) != 20U || B(2,2) != 19U || B(2,3) != 18U || B(2,4) != 17U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/scalar bitwise XOR operation failed\n"
                << " Details:\n"
                << "   Result:\n" << B << "\n"
                << "   Expected result:\n( 15 14 13 12 11 )\n"
                                        "( 10  9  8 23 22 )\n"
                                        "( 21 20 19 18 17 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Matrix/scalar bitwise XOR assignment
      {
         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> A{ {  8U,  9U, 10U, 11U, 12U },
                                                               { 13U, 14U, 15U, 16U, 17U },
                                                               { 18U, 19U, 20U, 21U, 22U } };

         A ^= 7U;

         checkRows    ( A,  3UL );
         checkColumns ( A,  5UL );
         checkCapacity( A, 15UL );
         checkNonZeros( A, 15UL );

         if( A(0,0) != 15U || A(0,1) != 14U || A(0,2) != 13U || A(0,3) != 12U || A(0,4) != 11U ||
             A(1,0) != 10U || A(1,1) !=  9U || A(1,2) !=  8U || A(1,3) != 23U || A(1,4) != 22U ||
             A(2,0) != 21U || A(2,1) != 20U || A(2,2) != 19U || A(2,3) != 18U || A(2,4) != 17U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/scalar bitwise XOR assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << A << "\n"
                << "   Expected result:\n( 15 14 13 12 11 )\n"
                                        "( 10  9  8 23 22 )\n"
                                        "( 21 20 19 18 17 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Row-major matrix/row-major matrix bitwise XOR tests
   //=====================================================================================

   {
      test_ = "Row-major matrix/row-major matrix bitwise XOR operator";

      // Matrix/matrix bitwise XOR of an empty matrix
      {
         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> A;
         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> B;

         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> C( A ^ B );

         checkRows    ( C, 0UL );
         checkColumns ( C, 0UL );
         checkCapacity( C, 0UL );
         checkNonZeros( C, 0UL );
      }

      // Matrix/matrix bitwise XOR of a general matrix
      {
         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> A{ {  8U,  9U, 10U, 11U, 12U },
                                                               { 13U, 14U, 15U, 16U, 17U },
                                                               { 18U, 19U, 20U, 21U, 22U } };

         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> B{ { 7U, 5U, 7U, 5U, 7U },
                                                               { 5U, 7U, 5U, 7U, 5U },
                                                               { 7U, 5U, 7U, 5U, 7U } };

         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> C( A ^ B );

         checkRows    ( C,  3UL );
         checkColumns ( C,  5UL );
         checkCapacity( C, 15UL );
         checkNonZeros( C, 15UL );

         if( C(0,0) != 15U || C(0,1) != 12U || C(0,2) != 13U || C(0,3) != 14U || C(0,4) != 11U ||
             C(1,0) !=  8U || C(1,1) !=  9U || C(1,2) != 10U || C(1,3) != 23U || C(1,4) != 20U ||
             C(2,0) != 21U || C(2,1) != 22U || C(2,2) != 19U || C(2,3) != 16U || C(2,4) != 17U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/matrix bitwise XOR operation failed\n"
                << " Details:\n"
                << "   Result:\n" << C << "\n"
                << "   Expected result:\n( 15 12 13 14 11 )\n"
                                        "(  8  9 10 23 20 )\n"
                                        "( 21 22 19 16 17 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Matrix/matrix bitwise XOR assignment
      {
         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> A{ {  8U,  9U, 10U, 11U, 12U },
                                                               { 13U, 14U, 15U, 16U, 17U },
                                                               { 18U, 19U, 20U, 21U, 22U } };

         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> B{ { 7U, 5U, 7U, 5U, 7U },
                                                               { 5U, 7U, 5U, 7U, 5U },
                                                               { 7U, 5U, 7U, 5U, 7U } };

         A ^= B;

         checkRows    ( A,  3UL );
         checkColumns ( A,  5UL );
         checkCapacity( A, 15UL );
         checkNonZeros( A, 15UL );

         if( A(0,0) != 15U || A(0,1) != 12U || A(0,2) != 13U || A(0,3) != 14U || A(0,4) != 11U ||
             A(1,0) !=  8U || A(1,1) !=  9U || A(1,2) != 10U || A(1,3) != 23U || A(1,4) != 20U ||
             A(2,0) != 21U || A(2,1) != 22U || A(2,2) != 19U || A(2,3) != 16U || A(2,4) != 17U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/matrix bitwise XOR assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << A << "\n"
                << "   Expected result:\n( 15 12 13 14 11 )\n"
                                        "(  8  9 10 23 20 )\n"
                                        "( 21 22 19 16 17 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Row-major matrix/column-major matrix bitwise XOR tests
   //=====================================================================================

   {
      test_ = "Row-major matrix/column-major matrix bitwise XOR operator";

      // Matrix/matrix bitwise XOR of an empty matrix
      {
         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> A;
         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> B;

         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> C( A ^ B );

         checkRows    ( C, 0UL );
         checkColumns ( C, 0UL );
         checkCapacity( C, 0UL );
         checkNonZeros( C, 0UL );
      }

      // Matrix/matrix bitwise XOR of a general matrix
      {
         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> A{ {  8U,  9U, 10U, 11U, 12U },
                                                               { 13U, 14U, 15U, 16U, 17U },
                                                               { 18U, 19U, 20U, 21U, 22U } };

         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> B{ { 7U, 5U, 7U, 5U, 7U },
                                                                  { 5U, 7U, 5U, 7U, 5U },
                                                                  { 7U, 5U, 7U, 5U, 7U } };

         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> C( A ^ B );

         checkRows    ( C,  3UL );
         checkColumns ( C,  5UL );
         checkCapacity( C, 15UL );
         checkNonZeros( C, 15UL );

         if( C(0,0) != 15U || C(0,1) != 12U || C(0,2) != 13U || C(0,3) != 14U || C(0,4) != 11U ||
             C(1,0) !=  8U || C(1,1) !=  9U || C(1,2) != 10U || C(1,3) != 23U || C(1,4) != 20U ||
             C(2,0) != 21U || C(2,1) != 22U || C(2,2) != 19U || C(2,3) != 16U || C(2,4) != 17U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/matrix bitwise XOR operation failed\n"
                << " Details:\n"
                << "   Result:\n" << C << "\n"
                << "   Expected result:\n( 15 12 13 14 11 )\n"
                                        "(  8  9 10 23 20 )\n"
                                        "( 21 22 19 16 17 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Matrix/matrix bitwise XOR assignment
      {
         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> A{ {  8U,  9U, 10U, 11U, 12U },
                                                               { 13U, 14U, 15U, 16U, 17U },
                                                               { 18U, 19U, 20U, 21U, 22U } };

         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> B{ { 7U, 5U, 7U, 5U, 7U },
                                                                  { 5U, 7U, 5U, 7U, 5U },
                                                                  { 7U, 5U, 7U, 5U, 7U } };

         A ^= B;

         checkRows    ( A,  3UL );
         checkColumns ( A,  5UL );
         checkCapacity( A, 15UL );
         checkNonZeros( A, 15UL );

         if( A(0,0) != 15U || A(0,1) != 12U || A(0,2) != 13U || A(0,3) != 14U || A(0,4) != 11U ||
             A(1,0) !=  8U || A(1,1) !=  9U || A(1,2) != 10U || A(1,3) != 23U || A(1,4) != 20U ||
             A(2,0) != 21U || A(2,1) != 22U || A(2,2) != 19U || A(2,3) != 16U || A(2,4) != 17U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/matrix bitwise XOR assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << A << "\n"
                << "   Expected result:\n( 15 12 13 14 11 )\n"
                                        "(  8  9 10 23 20 )\n"
                                        "( 21 22 19 16 17 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix/scalar bitwise XOR tests
   //=====================================================================================

   {
      test_ = "Column-major matrix/scalar bitwise XOR operator";

      // Matrix/scalar bitwise XOR of an empty matrix
      {
         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> A;

         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> B( A ^ 7U );

         checkRows    ( B, 0UL );
         checkColumns ( B, 0UL );
         checkCapacity( B, 0UL );
         checkNonZeros( B, 0UL );
      }

      // Matrix/scalar bitwise XOR of a general matrix
      {
         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> A{ {  8U,  9U, 10U, 11U, 12U },
                                                                  { 13U, 14U, 15U, 16U, 17U },
                                                                  { 18U, 19U, 20U, 21U, 22U } };

         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> B( A ^ 7U );

         checkRows    ( B,  3UL );
         checkColumns ( B,  5UL );
         checkCapacity( B, 15UL );
         checkNonZeros( B, 15UL );

         if( B(0,0) != 15U || B(0,1) != 14U || B(0,2) != 13U || B(0,3) != 12U || B(0,4) != 11U ||
             B(1,0) != 10U || B(1,1) !=  9U || B(1,2) !=  8U || B(1,3) != 23U || B(1,4) != 22U ||
             B(2,0) != 21U || B(2,1) != 20U || B(2,2) != 19U || B(2,3) != 18U || B(2,4) != 17U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/scalar bitwise XOR operation failed\n"
                << " Details:\n"
                << "   Result:\n" << B << "\n"
                << "   Expected result:\n( 15 14 13 12 11 )\n"
                                        "( 10  9  8 23 22 )\n"
                                        "( 21 20 19 18 17 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Matrix/scalar bitwise XOR assignment
      {
         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> A{ {  8U,  9U, 10U, 11U, 12U },
                                                                  { 13U, 14U, 15U, 16U, 17U },
                                                                  { 18U, 19U, 20U, 21U, 22U } };

         A ^= 7U;

         checkRows    ( A,  3UL );
         checkColumns ( A,  5UL );
         checkCapacity( A, 15UL );
         checkNonZeros( A, 15UL );

         if( A(0,0) != 15U || A(0,1) != 14U || A(0,2) != 13U || A(0,3) != 12U || A(0,4) != 11U ||
             A(1,0) != 10U || A(1,1) !=  9U || A(1,2) !=  8U || A(1,3) != 23U || A(1,4) != 22U ||
             A(2,0) != 21U || A(2,1) != 20U || A(2,2) != 19U || A(2,3) != 18U || A(2,4) != 17U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/scalar bitwise XOR assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << A << "\n"
                << "   Expected result:\n( 15 14 13 12 11 )\n"
                                        "( 10  9  8 23 22 )\n"
                                        "( 21 20 19 18 17 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix/row-major matrix bitwise XOR tests
   //=====================================================================================


   {
      test_ = "Column-major matrix/row-major matrix bitwise XOR operator";

      // Matrix/matrix bitwise XOR of an empty matrix
      {
         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> A;
         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> B;

         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> C( A ^ B );

         checkRows    ( C, 0UL );
         checkColumns ( C, 0UL );
         checkCapacity( C, 0UL );
         checkNonZeros( C, 0UL );
      }

      // Matrix/matrix bitwise XOR of a general matrix
      {
         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> A{ {  8U,  9U, 10U, 11U, 12U },
                                                                  { 13U, 14U, 15U, 16U, 17U },
                                                                  { 18U, 19U, 20U, 21U, 22U } };

         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> B{ { 7U, 5U, 7U, 5U, 7U },
                                                               { 5U, 7U, 5U, 7U, 5U },
                                                               { 7U, 5U, 7U, 5U, 7U } };

         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> C( A ^ B );

         checkRows    ( C,  3UL );
         checkColumns ( C,  5UL );
         checkCapacity( C, 15UL );
         checkNonZeros( C, 15UL );

         if( C(0,0) != 15U || C(0,1) != 12U || C(0,2) != 13U || C(0,3) != 14U || C(0,4) != 11U ||
             C(1,0) !=  8U || C(1,1) !=  9U || C(1,2) != 10U || C(1,3) != 23U || C(1,4) != 20U ||
             C(2,0) != 21U || C(2,1) != 22U || C(2,2) != 19U || C(2,3) != 16U || C(2,4) != 17U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/matrix bitwise XOR operation failed\n"
                << " Details:\n"
                << "   Result:\n" << C << "\n"
                << "   Expected result:\n( 15 12 13 14 11 )\n"
                                        "(  8  9 10 23 20 )\n"
                                        "( 21 22 19 16 17 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Matrix/matrix bitwise XOR assignment
      {
         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> A{ {  8U,  9U, 10U, 11U, 12U },
                                                                  { 13U, 14U, 15U, 16U, 17U },
                                                                  { 18U, 19U, 20U, 21U, 22U } };

         blaze::DynamicMatrix<unsigned int,blaze::rowMajor> B{ { 7U, 5U, 7U, 5U, 7U },
                                                               { 5U, 7U, 5U, 7U, 5U },
                                                               { 7U, 5U, 7U, 5U, 7U } };

         A ^= B;

         checkRows    ( A,  3UL );
         checkColumns ( A,  5UL );
         checkCapacity( A, 15UL );
         checkNonZeros( A, 15UL );

         if( A(0,0) != 15U || A(0,1) != 12U || A(0,2) != 13U || A(0,3) != 14U || A(0,4) != 11U ||
             A(1,0) !=  8U || A(1,1) !=  9U || A(1,2) != 10U || A(1,3) != 23U || A(1,4) != 20U ||
             A(2,0) != 21U || A(2,1) != 22U || A(2,2) != 19U || A(2,3) != 16U || A(2,4) != 17U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/matrix bitwise XOR assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << A << "\n"
                << "   Expected result:\n( 15 12 13 14 11 )\n"
                                        "(  8  9 10 23 20 )\n"
                                        "( 21 22 19 16 17 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix/column-major matrix bitwise XOR tests
   //=====================================================================================

   {
      test_ = "Column-major matrix/column-major matrix bitwise XOR operator";

      // Matrix/matrix bitwise XOR of an empty matrix
      {
         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> A;
         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> B;

         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> C( A ^ B );

         checkRows    ( C, 0UL );
         checkColumns ( C, 0UL );
         checkCapacity( C, 0UL );
         checkNonZeros( C, 0UL );
      }

      // Matrix/matrix bitwise XOR of a general matrix
      {
         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> A{ {  8U,  9U, 10U, 11U, 12U },
                                                                  { 13U, 14U, 15U, 16U, 17U },
                                                                  { 18U, 19U, 20U, 21U, 22U } };

         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> B{ { 7U, 5U, 7U, 5U, 7U },
                                                                  { 5U, 7U, 5U, 7U, 5U },
                                                                  { 7U, 5U, 7U, 5U, 7U } };

         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> C( A ^ B );

         checkRows    ( C,  3UL );
         checkColumns ( C,  5UL );
         checkCapacity( C, 15UL );
         checkNonZeros( C, 15UL );

         if( C(0,0) != 15U || C(0,1) != 12U || C(0,2) != 13U || C(0,3) != 14U || C(0,4) != 11U ||
             C(1,0) !=  8U || C(1,1) !=  9U || C(1,2) != 10U || C(1,3) != 23U || C(1,4) != 20U ||
             C(2,0) != 21U || C(2,1) != 22U || C(2,2) != 19U || C(2,3) != 16U || C(2,4) != 17U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/matrix bitwise XOR operation failed\n"
                << " Details:\n"
                << "   Result:\n" << C << "\n"
                << "   Expected result:\n( 15 12 13 14 11 )\n"
                                        "(  8  9 10 23 20 )\n"
                                        "( 21 22 19 16 17 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Matrix/matrix bitwise XOR assignment
      {
         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> A{ {  8U,  9U, 10U, 11U, 12U },
                                                                  { 13U, 14U, 15U, 16U, 17U },
                                                                  { 18U, 19U, 20U, 21U, 22U } };

         blaze::DynamicMatrix<unsigned int,blaze::columnMajor> B{ { 7U, 5U, 7U, 5U, 7U },
                                                                  { 5U, 7U, 5U, 7U, 5U },
                                                                  { 7U, 5U, 7U, 5U, 7U } };

         A ^= B;

         checkRows    ( A,  3UL );
         checkColumns ( A,  5UL );
         checkCapacity( A, 15UL );
         checkNonZeros( A, 15UL );

         if( A(0,0) != 15U || A(0,1) != 12U || A(0,2) != 13U || A(0,3) != 14U || A(0,4) != 11U ||
             A(1,0) !=  8U || A(1,1) !=  9U || A(1,2) != 10U || A(1,3) != 23U || A(1,4) != 20U ||
             A(2,0) != 21U || A(2,1) != 22U || A(2,2) != 19U || A(2,3) != 16U || A(2,4) != 17U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/matrix bitwise XOR assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << A << "\n"
                << "   Expected result:\n( 15 12 13 14 11 )\n"
                                        "(  8  9 10 23 20 )\n"
                                        "( 21 22 19 16 17 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the logical NOT operator for dense matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the logical NOT operator for dense matrices. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testNot()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major logical NOT operator";

      // Matrix logical NOT of an empty matrix
      {
         blaze::DynamicMatrix<bool,blaze::rowMajor> A;
         blaze::DynamicMatrix<bool,blaze::rowMajor> B( !A );

         checkRows    ( B, 0UL );
         checkColumns ( B, 0UL );
         checkCapacity( B, 0UL );
         checkNonZeros( B, 0UL );
      }

      // Matrix logical NOT of a general matrix
      {
         blaze::DynamicMatrix<bool,blaze::rowMajor> A{ { true , false, true , false },
                                                       { false, true , false, true  },
                                                       { true , false, true , false } };

         blaze::DynamicMatrix<bool,blaze::rowMajor> B( !A );

         checkRows    ( B,  3UL );
         checkColumns ( B,  4UL );
         checkCapacity( B, 12UL );
         checkNonZeros( B,  6UL );

         if( B(0,0) != false || B(0,1) != true  || B(0,2) != false || B(0,3) != true  ||
             B(1,0) != true  || B(1,1) != false || B(1,2) != true  || B(1,3) != false ||
             B(2,0) != false || B(2,1) != true  || B(2,2) != false || B(2,3) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix logical NOT operation failed\n"
                << " Details:\n"
                << "   Result:\n" << B << "\n"
                << "   Expected result:\n( 0 1 0 1 )\n"
                                        "( 1 0 1 0 )\n"
                                        "( 0 1 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major logical NOT operator";

      // Matrix logical NOT of an empty matrix
      {
         blaze::DynamicMatrix<bool,blaze::columnMajor> A;
         blaze::DynamicMatrix<bool,blaze::columnMajor> B( !A );

         checkRows    ( B, 0UL );
         checkColumns ( B, 0UL );
         checkCapacity( B, 0UL );
         checkNonZeros( B, 0UL );
      }

      // Matrix logical NOT of a general matrix
      {
         blaze::DynamicMatrix<bool,blaze::columnMajor> A{ { true , false, true , false },
                                                          { false, true , false, true  },
                                                          { true , false, true , false } };

         blaze::DynamicMatrix<bool,blaze::columnMajor> B( !A );

         checkRows    ( B,  3UL );
         checkColumns ( B,  4UL );
         checkCapacity( B, 12UL );
         checkNonZeros( B,  6UL );

         if( B(0,0) != false || B(0,1) != true  || B(0,2) != false || B(0,3) != true  ||
             B(1,0) != true  || B(1,1) != false || B(1,2) != true  || B(1,3) != false ||
             B(2,0) != false || B(2,1) != true  || B(2,2) != false || B(2,3) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix logical NOT operation failed\n"
                << " Details:\n"
                << "   Result:\n" << B << "\n"
                << "   Expected result:\n( 0 1 0 1 )\n"
                                        "( 1 0 1 0 )\n"
                                        "( 0 1 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the logical AND operator for dense matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the logical AND operator for dense matrices. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testAnd()
{
   //=====================================================================================
   // Row-major matrix/row-major matrix logical AND tests
   //=====================================================================================

   {
      test_ = "Row-major matrix/row-major matrix logical AND operator";

      // Matrix/matrix logical AND of an empty matrix
      {
         blaze::DynamicMatrix<bool,blaze::rowMajor> A;
         blaze::DynamicMatrix<bool,blaze::rowMajor> B;

         blaze::DynamicMatrix<bool,blaze::rowMajor> C( A && B );

         checkRows    ( C, 0UL );
         checkColumns ( C, 0UL );
         checkCapacity( C, 0UL );
         checkNonZeros( C, 0UL );
      }

      // Matrix/matrix logical AND of a general matrix
      {
         blaze::DynamicMatrix<bool,blaze::rowMajor> A{ { true, false, true, false },
                                                       { false, true, false, true },
                                                       { true, false, true, false } };

         blaze::DynamicMatrix<bool,blaze::rowMajor> B{ { true, true, false, false },
                                                       { false, false, true, true },
                                                       { true, true, false, false } };

         blaze::DynamicMatrix<bool,blaze::rowMajor> C( A && B );

         checkRows    ( C,  3UL );
         checkColumns ( C,  4UL );
         checkCapacity( C, 12UL );
         checkNonZeros( C,  3UL );

         if( C(0,0) != true  || C(0,1) != false || C(0,2) != false || C(0,3) != false ||
             C(1,0) != false || C(1,1) != false || C(1,2) != false || C(1,3) != true  ||
             C(2,0) != true  || C(2,1) != false || C(2,2) != false || C(2,3) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/matrix logical AND operation failed\n"
                << " Details:\n"
                << "   Result:\n" << C << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n"
                                        "( 0 0 0 1 )\n"
                                        "( 1 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Row-major matrix/column-major matrix logical AND tests
   //=====================================================================================

   {
      test_ = "Row-major matrix/column-major matrix logical AND operator";

      // Matrix/matrix logical AND of an empty matrix
      {
         blaze::DynamicMatrix<bool,blaze::rowMajor> A;
         blaze::DynamicMatrix<bool,blaze::columnMajor> B;

         blaze::DynamicMatrix<bool,blaze::rowMajor> C( A && B );

         checkRows    ( C, 0UL );
         checkColumns ( C, 0UL );
         checkCapacity( C, 0UL );
         checkNonZeros( C, 0UL );
      }

      // Matrix/matrix logical AND of a general matrix
      {
         blaze::DynamicMatrix<bool,blaze::rowMajor> A{ { true, false, true, false },
                                                       { false, true, false, true },
                                                       { true, false, true, false } };

         blaze::DynamicMatrix<bool,blaze::columnMajor> B{ { true, true, false, false },
                                                          { false, false, true, true },
                                                          { true, true, false, false } };

         blaze::DynamicMatrix<bool,blaze::rowMajor> C( A && B );

         checkRows    ( C,  3UL );
         checkColumns ( C,  4UL );
         checkCapacity( C, 12UL );
         checkNonZeros( C,  3UL );

         if( C(0,0) != true  || C(0,1) != false || C(0,2) != false || C(0,3) != false ||
             C(1,0) != false || C(1,1) != false || C(1,2) != false || C(1,3) != true  ||
             C(2,0) != true  || C(2,1) != false || C(2,2) != false || C(2,3) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/matrix logical AND operation failed\n"
                << " Details:\n"
                << "   Result:\n" << C << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n"
                                        "( 0 0 0 1 )\n"
                                        "( 1 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix/row-major matrix logical AND tests
   //=====================================================================================

   {
      test_ = "Column-major matrix/row-major matrix logical AND operator";

      // Matrix/matrix logical AND of an empty matrix
      {
         blaze::DynamicMatrix<bool,blaze::columnMajor> A;
         blaze::DynamicMatrix<bool,blaze::rowMajor> B;

         blaze::DynamicMatrix<bool,blaze::columnMajor> C( A && B );

         checkRows    ( C, 0UL );
         checkColumns ( C, 0UL );
         checkCapacity( C, 0UL );
         checkNonZeros( C, 0UL );
      }

      // Matrix/matrix logical AND of a general matrix
      {
         blaze::DynamicMatrix<bool,blaze::columnMajor> A{ { true, false, true, false },
                                                          { false, true, false, true },
                                                          { true, false, true, false } };

         blaze::DynamicMatrix<bool,blaze::rowMajor> B{ { true, true, false, false },
                                                       { false, false, true, true },
                                                       { true, true, false, false } };

         blaze::DynamicMatrix<bool,blaze::columnMajor> C( A && B );

         checkRows    ( C,  3UL );
         checkColumns ( C,  4UL );
         checkCapacity( C, 12UL );
         checkNonZeros( C,  3UL );

         if( C(0,0) != true  || C(0,1) != false || C(0,2) != false || C(0,3) != false ||
             C(1,0) != false || C(1,1) != false || C(1,2) != false || C(1,3) != true  ||
             C(2,0) != true  || C(2,1) != false || C(2,2) != false || C(2,3) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/matrix logical AND operation failed\n"
                << " Details:\n"
                << "   Result:\n" << C << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n"
                                        "( 0 0 0 1 )\n"
                                        "( 1 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix/column-major matrix logical AND tests
   //=====================================================================================

   {
      test_ = "Column-major matrix/column-major matrix logical AND operator";

      // Matrix/matrix logical AND of an empty matrix
      {
         blaze::DynamicMatrix<bool,blaze::columnMajor> A;
         blaze::DynamicMatrix<bool,blaze::columnMajor> B;

         blaze::DynamicMatrix<bool,blaze::columnMajor> C( A && B );

         checkRows    ( C, 0UL );
         checkColumns ( C, 0UL );
         checkCapacity( C, 0UL );
         checkNonZeros( C, 0UL );
      }

      // Matrix/matrix logical AND of a general matrix
      {
         blaze::DynamicMatrix<bool,blaze::columnMajor> A{ { true, false, true, false },
                                                          { false, true, false, true },
                                                          { true, false, true, false } };

         blaze::DynamicMatrix<bool,blaze::columnMajor> B{ { true, true, false, false },
                                                          { false, false, true, true },
                                                          { true, true, false, false } };

         blaze::DynamicMatrix<bool,blaze::columnMajor> C( A && B );

         checkRows    ( C,  3UL );
         checkColumns ( C,  4UL );
         checkCapacity( C, 12UL );
         checkNonZeros( C,  3UL );

         if( C(0,0) != true  || C(0,1) != false || C(0,2) != false || C(0,3) != false ||
             C(1,0) != false || C(1,1) != false || C(1,2) != false || C(1,3) != true  ||
             C(2,0) != true  || C(2,1) != false || C(2,2) != false || C(2,3) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/matrix logical AND operation failed\n"
                << " Details:\n"
                << "   Result:\n" << C << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n"
                                        "( 0 0 0 1 )\n"
                                        "( 1 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the logical OR operator for dense matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the logical OR operator for dense matrices. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testOr()
{
   //=====================================================================================
   // Row-major matrix/row-major matrix logical OR tests
   //=====================================================================================

   {
      test_ = "Row-major matrix/row-major matrix logical OR operator";

      // Matrix/matrix logical OR of an empty matrix
      {
         blaze::DynamicMatrix<bool,blaze::rowMajor> A;
         blaze::DynamicMatrix<bool,blaze::rowMajor> B;

         blaze::DynamicMatrix<bool,blaze::rowMajor> C( A || B );

         checkRows    ( C, 0UL );
         checkColumns ( C, 0UL );
         checkCapacity( C, 0UL );
         checkNonZeros( C, 0UL );
      }

      // Matrix/matrix logical OR of a general matrix
      {
         blaze::DynamicMatrix<bool,blaze::rowMajor> A{ { true, false, true, false },
                                                       { false, true, false, true },
                                                       { true, false, true, false } };

         blaze::DynamicMatrix<bool,blaze::rowMajor> B{ { true, true, false, false },
                                                       { false, false, true, true },
                                                       { true, true, false, false } };

         blaze::DynamicMatrix<bool,blaze::rowMajor> C( A || B );

         checkRows    ( C,  3UL );
         checkColumns ( C,  4UL );
         checkCapacity( C, 12UL );
         checkNonZeros( C,  9UL );

         if( C(0,0) != true  || C(0,1) != true || C(0,2) != true || C(0,3) != false ||
             C(1,0) != false || C(1,1) != true || C(1,2) != true || C(1,3) != true  ||
             C(2,0) != true  || C(2,1) != true || C(2,2) != true || C(2,3) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/matrix logical OR operation failed\n"
                << " Details:\n"
                << "   Result:\n" << C << "\n"
                << "   Expected result:\n( 1 1 1 0 )\n"
                                        "( 0 1 1 1 )\n"
                                        "( 1 1 1 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Row-major matrix/column-major matrix logical OR tests
   //=====================================================================================

   {
      test_ = "Row-major matrix/column-major matrix logical OR operator";

      // Matrix/matrix logical OR of an empty matrix
      {
         blaze::DynamicMatrix<bool,blaze::rowMajor> A;
         blaze::DynamicMatrix<bool,blaze::columnMajor> B;

         blaze::DynamicMatrix<bool,blaze::rowMajor> C( A || B );

         checkRows    ( C, 0UL );
         checkColumns ( C, 0UL );
         checkCapacity( C, 0UL );
         checkNonZeros( C, 0UL );
      }

      // Matrix/matrix logical OR of a general matrix
      {
         blaze::DynamicMatrix<bool,blaze::rowMajor> A{ { true, false, true, false },
                                                       { false, true, false, true },
                                                       { true, false, true, false } };

         blaze::DynamicMatrix<bool,blaze::columnMajor> B{ { true, true, false, false },
                                                          { false, false, true, true },
                                                          { true, true, false, false } };

         blaze::DynamicMatrix<bool,blaze::rowMajor> C( A || B );

         checkRows    ( C,  3UL );
         checkColumns ( C,  4UL );
         checkCapacity( C, 12UL );
         checkNonZeros( C,  9UL );

         if( C(0,0) != true  || C(0,1) != true || C(0,2) != true || C(0,3) != false ||
             C(1,0) != false || C(1,1) != true || C(1,2) != true || C(1,3) != true  ||
             C(2,0) != true  || C(2,1) != true || C(2,2) != true || C(2,3) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/matrix logical OR operation failed\n"
                << " Details:\n"
                << "   Result:\n" << C << "\n"
                << "   Expected result:\n( 1 1 1 0 )\n"
                                        "( 0 1 1 1 )\n"
                                        "( 1 1 1 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix/row-major matrix logical OR tests
   //=====================================================================================

   {
      test_ = "Column-major matrix/row-major matrix logical OR operator";

      // Matrix/matrix logical OR of an empty matrix
      {
         blaze::DynamicMatrix<bool,blaze::columnMajor> A;
         blaze::DynamicMatrix<bool,blaze::rowMajor> B;

         blaze::DynamicMatrix<bool,blaze::columnMajor> C( A || B );

         checkRows    ( C, 0UL );
         checkColumns ( C, 0UL );
         checkCapacity( C, 0UL );
         checkNonZeros( C, 0UL );
      }

      // Matrix/matrix logical OR of a general matrix
      {
         blaze::DynamicMatrix<bool,blaze::columnMajor> A{ { true, false, true, false },
                                                          { false, true, false, true },
                                                          { true, false, true, false } };

         blaze::DynamicMatrix<bool,blaze::rowMajor> B{ { true, true, false, false },
                                                       { false, false, true, true },
                                                       { true, true, false, false } };

         blaze::DynamicMatrix<bool,blaze::columnMajor> C( A || B );

         checkRows    ( C,  3UL );
         checkColumns ( C,  4UL );
         checkCapacity( C, 12UL );
         checkNonZeros( C,  9UL );

         if( C(0,0) != true  || C(0,1) != true || C(0,2) != true || C(0,3) != false ||
             C(1,0) != false || C(1,1) != true || C(1,2) != true || C(1,3) != true  ||
             C(2,0) != true  || C(2,1) != true || C(2,2) != true || C(2,3) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/matrix logical OR operation failed\n"
                << " Details:\n"
                << "   Result:\n" << C << "\n"
                << "   Expected result:\n( 1 1 1 0 )\n"
                                        "( 0 1 1 1 )\n"
                                        "( 1 1 1 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix/column-major matrix logical OR tests
   //=====================================================================================

   {
      test_ = "Column-major matrix/column-major matrix logical OR operator";

      // Matrix/matrix logical OR of an empty matrix
      {
         blaze::DynamicMatrix<bool,blaze::columnMajor> A;
         blaze::DynamicMatrix<bool,blaze::columnMajor> B;

         blaze::DynamicMatrix<bool,blaze::columnMajor> C( A || B );

         checkRows    ( C, 0UL );
         checkColumns ( C, 0UL );
         checkCapacity( C, 0UL );
         checkNonZeros( C, 0UL );
      }

      // Matrix/matrix logical OR of a general matrix
      {
         blaze::DynamicMatrix<bool,blaze::columnMajor> A{ { true, false, true, false },
                                                          { false, true, false, true },
                                                          { true, false, true, false } };

         blaze::DynamicMatrix<bool,blaze::columnMajor> B{ { true, true, false, false },
                                                          { false, false, true, true },
                                                          { true, true, false, false } };

         blaze::DynamicMatrix<bool,blaze::columnMajor> C( A || B );

         checkRows    ( C,  3UL );
         checkColumns ( C,  4UL );
         checkCapacity( C, 12UL );
         checkNonZeros( C,  9UL );

         if( C(0,0) != true  || C(0,1) != true || C(0,2) != true || C(0,3) != false ||
             C(1,0) != false || C(1,1) != true || C(1,2) != true || C(1,3) != true  ||
             C(2,0) != true  || C(2,1) != true || C(2,2) != true || C(2,3) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix/matrix logical OR operation failed\n"
                << " Details:\n"
                << "   Result:\n" << C << "\n"
                << "   Expected result:\n( 1 1 1 0 )\n"
                                        "( 0 1 1 1 )\n"
                                        "( 1 1 1 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c generate() functions for dense matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c generate() functions for dense matrices. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testGenerate()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   // Empty integer matrix
   {
      const blaze::DynamicMatrix<int,blaze::rowMajor> mat(
         blaze::generate( 0UL, 0UL, []( size_t, size_t ){ return 2; } ) );

      const blaze::DynamicMatrix<int> ref{};

      if( mat != ref ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Generating empty integer matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n" << ref << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Single element integer matrix ( ( 2 ) )
   {
      const blaze::DynamicMatrix<int,blaze::rowMajor> mat(
         blaze::generate( 1UL, 1UL, []( size_t, size_t ){ return 2; } ) );

      const blaze::DynamicMatrix<int,blaze::rowMajor> ref{ { 2 } };

      if( mat != ref ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Generating single element integer matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n" << ref << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Uniform integer matrix ( ( 2, 2, 2 ), ( 2, 2, 2 ) )
   {
      const blaze::DynamicMatrix<int,blaze::rowMajor> mat(
         blaze::generate( 2UL, 3UL, []( size_t, size_t ){ return 2; } ) );

      const blaze::DynamicMatrix<int,blaze::rowMajor> ref{ { 2, 2, 2 }, { 2, 2, 2 } };

      if( mat != ref ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Generating uniform integer matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n" << ref << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Linearly spaced float matrix ( ( 2.1, 3.2, 4.3 ), ( 5.4, 6.5, 7.6 ) )
   {
      const blaze::DynamicMatrix<float,blaze::rowMajor> mat(
         blaze::generate( 2UL, 3UL, []( size_t i, size_t j ){ return 2.1F + 1.1F*(i*3UL+j); } ) );

      const blaze::DynamicMatrix<float,blaze::rowMajor> ref{ { 2.1F, 3.2F, 4.3F },
                                                             { 5.4F, 6.5F, 7.6F } };

      if( mat != ref ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Generating linearly spaced float matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n" << ref << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Logarithmically spaced double matrix ( ( 10.0, 100.0 ), ( 1000.0, 10000.0 ) )
   {
      const blaze::DynamicMatrix<double,blaze::rowMajor> mat(
         blaze::generate( 2UL, 2UL, []( size_t i, size_t j ){ return blaze::exp10( 1.0 + 1.0*(i*2UL+j) ); } ) );

      const blaze::DynamicMatrix<double,blaze::rowMajor> ref{ {   10.0,   100.0 },
                                                              { 1000.0, 10000.0 } };

      if( mat != ref ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Generating logarithmically spaced double matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n" << ref << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Vector of vectors
   {
      using VT = blaze::StaticVector<int,2UL>;

      const blaze::DynamicMatrix<VT,blaze::rowMajor> mat(
         blaze::generate( 2UL, 2UL, []( size_t i, size_t j ) { return evaluate( VT{ 1, 2 } + i*2UL + j ); } ) );

      const blaze::DynamicMatrix<VT,blaze::rowMajor> ref{ { { 1, 2 }, { 2, 3 } },
                                                          { { 3, 4 }, { 4, 5 } } };

      if( mat != ref ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Generating matrix of vectors failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n" << ref << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   // Empty integer matrix
   {
      const blaze::DynamicMatrix<int,blaze::columnMajor> mat(
         blaze::generate( 0UL, 0UL, []( size_t, size_t ){ return 2; } ) );

      const blaze::DynamicMatrix<int> ref{};

      if( mat != ref ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Generating empty integer matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n" << ref << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Single element integer matrix ( ( 2 ) )
   {
      const blaze::DynamicMatrix<int,blaze::columnMajor> mat(
         blaze::generate( 1UL, 1UL, []( size_t, size_t ){ return 2; } ) );

      const blaze::DynamicMatrix<int,blaze::columnMajor> ref{ { 2 } };

      if( mat != ref ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Generating single element integer matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n" << ref << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Uniform integer matrix ( ( 2, 2, 2 ), ( 2, 2, 2 ) )
   {
      const blaze::DynamicMatrix<int,blaze::columnMajor> mat(
         blaze::generate( 2UL, 3UL, []( size_t, size_t ){ return 2; } ) );

      const blaze::DynamicMatrix<int,blaze::columnMajor> ref{ { 2, 2, 2 }, { 2, 2, 2 } };

      if( mat != ref ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Generating uniform integer matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n" << ref << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Linearly spaced float matrix ( ( 2.1, 3.2, 4.3 ), ( 5.4, 6.5, 7.6 ) )
   {
      const blaze::DynamicMatrix<float,blaze::columnMajor> mat(
         blaze::generate( 2UL, 3UL, []( size_t i, size_t j ){ return 2.1F + 1.1F*(i*3UL+j); } ) );

      const blaze::DynamicMatrix<float,blaze::columnMajor> ref{ { 2.1F, 3.2F, 4.3F },
                                                                { 5.4F, 6.5F, 7.6F } };

      if( mat != ref ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Generating linearly spaced float matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n" << ref << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Logarithmically spaced double matrix ( ( 10.0, 100.0 ), ( 1000.0, 10000.0 ) )
   {
      const blaze::DynamicMatrix<double,blaze::columnMajor> mat(
         blaze::generate( 2UL, 2UL, []( size_t i, size_t j ){ return blaze::exp10( 1.0 + 1.0*(i*2UL+j) ); } ) );

      const blaze::DynamicMatrix<double,blaze::columnMajor> ref{ {   10.0,   100.0 },
                                                                 { 1000.0, 10000.0 } };

      if( mat != ref ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Generating logarithmically spaced double matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n" << ref << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Vector of vectors
   {
      using VT = blaze::StaticVector<int,2UL>;

      const blaze::DynamicMatrix<VT,blaze::columnMajor> mat(
         blaze::generate( 2UL, 2UL, []( size_t i, size_t j ) { return evaluate( VT{ 1, 2 } + i*2UL + j ); } ) );

      const blaze::DynamicMatrix<VT,blaze::columnMajor> ref{ { { 1, 2 }, { 2, 3 } },
                                                             { { 3, 4 }, { 4, 5 } } };

      if( mat != ref ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Generating matrix of vectors failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n" << ref << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c uniform() functions for dense matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c uniform() functions for dense matrices. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testUniform()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   // Empty integer matrix
   {
      const blaze::DynamicMatrix<int,blaze::rowMajor> mat( blaze::uniform( 0UL, 0UL, 5 ) );
      const blaze::DynamicMatrix<int,blaze::rowMajor> ref{};

      if( mat != ref ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Generating empty integer matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n" << ref << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Single element integer matrix ( ( 5 ) )
   {
      const blaze::DynamicMatrix<int,blaze::rowMajor> mat( blaze::uniform( 1UL, 1UL, 5 ) );
      const blaze::DynamicMatrix<int,blaze::rowMajor> ref{ { 5 } };

      if( mat != ref ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Generating single element integer matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n" << ref << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Uniform float matrix ( ( 2.1, 2.1, 2.1 ), ( 2.1, 2.1, 2.1 ) )
   {
      const blaze::DynamicMatrix<float,blaze::rowMajor> mat( blaze::uniform( 2UL, 3UL, 2.1f ) );
      const blaze::DynamicMatrix<float,blaze::rowMajor> ref{ { 2.1F, 2.1F, 2.1F },
                                                             { 2.1F, 2.1F, 2.1F } };

      if( mat != ref ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Generating uniform float matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n" << ref << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Uniform matrix of vectors
   {
      using VT = blaze::StaticVector<int,2UL>;

      const blaze::DynamicMatrix<VT,blaze::rowMajor> mat( blaze::uniform( 2UL, 3UL, VT{ 1, 2 } ) );
      const blaze::DynamicMatrix<VT,blaze::rowMajor> ref{ { { 1, 2 }, { 1, 2 }, { 1, 2 } },
                                                          { { 1, 2 }, { 1, 2 }, { 1, 2 } } };

      if( mat != ref ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Generating matrix of vectors failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n" << ref << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   // Empty integer matrix
   {
      const blaze::DynamicMatrix<int,blaze::columnMajor> mat( blaze::uniform( 0UL, 0UL, 5 ) );
      const blaze::DynamicMatrix<int,blaze::columnMajor> ref{};

      if( mat != ref ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Generating empty integer matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n" << ref << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Single element integer matrix ( ( 5 ) )
   {
      const blaze::DynamicMatrix<int,blaze::columnMajor> mat( blaze::uniform( 1UL, 1UL, 5 ) );
      const blaze::DynamicMatrix<int,blaze::columnMajor> ref{ { 5 } };

      if( mat != ref ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Generating single element integer matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n" << ref << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Uniform float matrix ( ( 2.1, 2.1, 2.1 ), ( 2.1, 2.1, 2.1 ) )
   {
      const blaze::DynamicMatrix<float,blaze::columnMajor> mat( blaze::uniform( 2UL, 3UL, 2.1f ) );
      const blaze::DynamicMatrix<float,blaze::columnMajor> ref{ { 2.1F, 2.1F, 2.1F },
                                                                { 2.1F, 2.1F, 2.1F } };

      if( mat != ref ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Generating uniform float matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n" << ref << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Uniform matrix of vectors
   {
      using VT = blaze::StaticVector<int,2UL>;

      const blaze::DynamicMatrix<VT,blaze::columnMajor> mat( blaze::uniform( 2UL, 3UL, VT{ 1, 2 } ) );
      const blaze::DynamicMatrix<VT,blaze::columnMajor> ref{ { { 1, 2 }, { 1, 2 }, { 1, 2 } },
                                                             { { 1, 2 }, { 1, 2 }, { 1, 2 } } };

      if( mat != ref ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Generating matrix of vectors failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n" << ref << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c zero() functions for dense matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c zero() functions for dense matrices. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testZero()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   // Empty integer matrix
   {
      const blaze::DynamicMatrix<int,blaze::rowMajor> mat( blaze::zero<int>( 0UL, 0UL ) );
      const blaze::DynamicMatrix<int,blaze::rowMajor> ref{};

      if( mat != ref ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Generating empty integer matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n" << ref << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Single element integer matrix ( ( 0 ) )
   {
      const blaze::DynamicMatrix<int,blaze::rowMajor> mat( blaze::zero<int>( 1UL, 1UL ) );
      const blaze::DynamicMatrix<int,blaze::rowMajor> ref{ { 0 } };

      if( mat != ref ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Generating single element integer matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n" << ref << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Uniform float matrix ( ( 0.0, 0.0, 0.0 ), ( 0.0, 0.0, 0.0 ) )
   {
      const blaze::DynamicMatrix<float,blaze::rowMajor> mat( blaze::zero<float>( 2UL, 3UL ) );
      const blaze::DynamicMatrix<float,blaze::rowMajor> ref{ { 0.0F, 0.0F, 0.0F },
                                                             { 0.0F, 0.0F, 0.0F } };

      if( mat != ref ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Generating zero float matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n" << ref << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Uniform matrix of vectors
   {
      using VT = blaze::StaticVector<int,2UL>;

      const blaze::DynamicMatrix<VT,blaze::rowMajor> mat( blaze::uniform( 2UL, 3UL, VT{} ) );
      const blaze::DynamicMatrix<VT,blaze::rowMajor> ref{ { VT{}, VT{}, VT{} },
                                                          { VT{}, VT{}, VT{} } };

      if( mat != ref ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Generating matrix of vectors failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n" << ref << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   // Empty integer matrix
   {
      const blaze::DynamicMatrix<int,blaze::columnMajor> mat( blaze::zero<int>( 0UL, 0UL ) );
      const blaze::DynamicMatrix<int,blaze::columnMajor> ref{};

      if( mat != ref ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Generating empty integer matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n" << ref << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Single element integer matrix ( ( 0 ) )
   {
      const blaze::DynamicMatrix<int,blaze::columnMajor> mat( blaze::zero<int>( 1UL, 1UL ) );
      const blaze::DynamicMatrix<int,blaze::columnMajor> ref{ { 0 } };

      if( mat != ref ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Generating single element integer matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n" << ref << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Uniform float matrix ( ( 0.0, 0.0, 0.0 ), ( 0.0, 0.0, 0.0 ) )
   {
      const blaze::DynamicMatrix<float,blaze::columnMajor> mat( blaze::zero<float>( 2UL, 3UL ) );
      const blaze::DynamicMatrix<float,blaze::columnMajor> ref{ { 0.0F, 0.0F, 0.0F },
                                                                { 0.0F, 0.0F, 0.0F } };

      if( mat != ref ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Generating zero float matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n" << ref << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Uniform matrix of vectors
   {
      using VT = blaze::StaticVector<int,2UL>;

      const blaze::DynamicMatrix<VT,blaze::columnMajor> mat( blaze::uniform( 2UL, 3UL, VT{} ) );
      const blaze::DynamicMatrix<VT,blaze::columnMajor> ref{ { VT{}, VT{}, VT{} },
                                                             { VT{}, VT{}, VT{} } };

      if( mat != ref ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Generating matrix of vectors failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n" << ref << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************

} // namespace densematrix

} // namespace matrices

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
   std::cout << "   Running general DenseMatrix operation test..." << std::endl;

   try
   {
      RUN_DENSEMATRIX_GENERAL_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during general DenseMatrix operation test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
