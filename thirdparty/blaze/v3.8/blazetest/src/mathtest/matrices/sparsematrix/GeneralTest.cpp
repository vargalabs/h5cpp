//=================================================================================================
/*!
//  \file src/mathtest/matrices/sparsematrix/GeneralTest.cpp
//  \brief Source file for the general SparseMatrix operation test
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
#include <blaze/math/sparse/SparseMatrix.h>
#include <blaze/math/CompressedMatrix.h>
#include <blaze/math/DynamicVector.h>
#include <blazetest/mathtest/IsEqual.h>
#include <blazetest/mathtest/matrices/sparsematrix/GeneralTest.h>

#ifdef BLAZE_USE_HPX_THREADS
#  include <hpx/hpx_main.hpp>
#endif


namespace blazetest {

namespace mathtest {

namespace matrices {

namespace sparsematrix {

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
   testMinimum();
   testMaximum();
   testL1Norm();
   testL2Norm();
   testL3Norm();
   testL4Norm();
   testLpNorm();
   testLinfNorm();
   testTrace();
   testMean();
   testVar();
   testStdDev();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the \c isnan() function for sparse matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isnan() function for sparse matrices. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testIsNan()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major isnan()";

      // isnan with 0x0 matrix
      {
         blaze::CompressedMatrix<float,blaze::rowMajor> mat;

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

      // isnan with empty 3x5 matrix (non-NAN)
      {
         blaze::CompressedMatrix<float,blaze::rowMajor> mat( 3UL, 5UL );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 5UL );
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

      // isnan with 4x2 matrix (non-NAN)
      {
         blaze::CompressedMatrix<float,blaze::rowMajor> mat( 4UL, 2UL );
         mat(1,1) =  1.0F;
         mat(2,0) = -2.0F;
         mat(2,1) =  3.0F;
         mat(3,0) =  4.0F;

         checkRows    ( mat, 4UL );
         checkColumns ( mat, 2UL );
         checkNonZeros( mat, 4UL );

         if( blaze::isnan( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isnan evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isnan with 4x2 matrix (NAN)
      {
         blaze::CompressedMatrix<float,blaze::rowMajor> mat( 4UL, 2UL );
         mat(1,1) =  1.0F;
         mat(2,0) = -2.0F;
         mat(2,1) =   NAN;
         mat(3,0) =  4.0F;

         checkRows    ( mat, 4UL );
         checkColumns ( mat, 2UL );
         checkNonZeros( mat, 4UL );

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

      // isnan with 0x0 matrix
      {
         blaze::CompressedMatrix<float,blaze::columnMajor> mat;

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

      // isnan with empty 3x5 matrix (non-NAN)
      {
         blaze::CompressedMatrix<float,blaze::columnMajor> mat( 3UL, 5UL );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 5UL );
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

      // isnan with 4x2 matrix (non-NAN)
      {
         blaze::CompressedMatrix<float,blaze::columnMajor> mat( 4UL, 2UL );
         mat(1,1) =  1.0F;
         mat(2,0) = -2.0F;
         mat(2,1) =  3.0F;
         mat(3,0) =  4.0F;

         checkRows    ( mat, 4UL );
         checkColumns ( mat, 2UL );
         checkNonZeros( mat, 4UL );

         if( blaze::isnan( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isnan evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isnan with 4x2 matrix (NAN)
      {
         blaze::CompressedMatrix<float,blaze::columnMajor> mat( 4UL, 2UL );
         mat(1,1) =  1.0F;
         mat(2,0) = -2.0F;
         mat(2,1) =   NAN;
         mat(3,0) =  4.0F;

         checkRows    ( mat, 4UL );
         checkColumns ( mat, 2UL );
         checkNonZeros( mat, 4UL );

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
/*!\brief Test of the \c isinf() function for sparse matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isinf() function for sparse matrices. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testIsInf()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major isinf()";

      // isinf with 0x0 matrix
      {
         blaze::CompressedMatrix<float,blaze::rowMajor> mat;

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

      // isinf with empty 3x5 matrix (non-inf)
      {
         blaze::CompressedMatrix<float,blaze::rowMajor> mat( 3UL, 5UL );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 5UL );
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

      // isinf with 4x2 matrix (non-inf)
      {
         blaze::CompressedMatrix<float,blaze::rowMajor> mat( 4UL, 2UL );
         mat(1,1) =  1.0F;
         mat(2,0) = -2.0F;
         mat(2,1) =  3.0F;
         mat(3,0) =  4.0F;

         checkRows    ( mat, 4UL );
         checkColumns ( mat, 2UL );
         checkNonZeros( mat, 4UL );

         if( blaze::isinf( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isinf evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isinf with 4x2 matrix (HUGE_VALF)
      {
         blaze::CompressedMatrix<float,blaze::rowMajor> mat( 4UL, 2UL );
         mat(1,1) =      1.0F;
         mat(2,0) =     -2.0F;
         mat(2,1) = HUGE_VALF;
         mat(3,0) =      4.0F;

         checkRows    ( mat, 4UL );
         checkColumns ( mat, 2UL );
         checkNonZeros( mat, 4UL );

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

      // isinf with 0x0 matrix
      {
         blaze::CompressedMatrix<float,blaze::columnMajor> mat;

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

      // isinf with empty 3x5 matrix (non-inf)
      {
         blaze::CompressedMatrix<float,blaze::columnMajor> mat( 3UL, 5UL );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 5UL );
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

      // isinf with 4x2 matrix (non-inf)
      {
         blaze::CompressedMatrix<float,blaze::columnMajor> mat( 4UL, 2UL );
         mat(1,1) =  1.0F;
         mat(2,0) = -2.0F;
         mat(2,1) =  3.0F;
         mat(3,0) =  4.0F;

         checkRows    ( mat, 4UL );
         checkColumns ( mat, 2UL );
         checkNonZeros( mat, 4UL );

         if( blaze::isinf( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isinf evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isinf with 4x2 matrix (HUGE_VALF)
      {
         blaze::CompressedMatrix<float,blaze::columnMajor> mat( 4UL, 2UL );
         mat(1,1) =      1.0F;
         mat(2,0) =     -2.0F;
         mat(2,1) = HUGE_VALF;
         mat(3,0) =      4.0F;

         checkRows    ( mat, 4UL );
         checkColumns ( mat, 2UL );
         checkNonZeros( mat, 4UL );

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
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c isfinite() function for sparse matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isfinite() function for sparse matrices. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testIsFinite()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major isfinite()";

      // isfinite with 0x0 matrix
      {
         blaze::CompressedMatrix<float,blaze::rowMajor> mat;

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

      // isfinite with empty 3x5 matrix (finite)
      {
         blaze::CompressedMatrix<float,blaze::rowMajor> mat( 3UL, 5UL );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 5UL );
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

      // isfinite with 4x2 matrix (finite)
      {
         blaze::CompressedMatrix<float,blaze::rowMajor> mat( 4UL, 2UL );
         mat(1,1) =  1.0F;
         mat(2,0) = -2.0F;
         mat(2,1) =  3.0F;
         mat(3,0) =  4.0F;

         checkRows    ( mat, 4UL );
         checkColumns ( mat, 2UL );
         checkNonZeros( mat, 4UL );

         if( blaze::isfinite( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isfinite evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isfinite with 4x2 matrix (NAN)
      {
         blaze::CompressedMatrix<float,blaze::rowMajor> mat( 4UL, 2UL );
         mat(1,1) =  1.0F;
         mat(2,0) = -2.0F;
         mat(2,1) =   NAN;
         mat(3,0) =  4.0F;

         checkRows    ( mat, 4UL );
         checkColumns ( mat, 2UL );
         checkNonZeros( mat, 4UL );

         if( blaze::isfinite( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isinf evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isfinite with 4x2 matrix (HUGE_VALF)
      {
         blaze::CompressedMatrix<float,blaze::rowMajor> mat( 4UL, 2UL );
         mat(1,1) =      1.0F;
         mat(2,0) =     -2.0F;
         mat(2,1) = HUGE_VALF;
         mat(3,0) =      4.0F;

         checkRows    ( mat, 4UL );
         checkColumns ( mat, 2UL );
         checkNonZeros( mat, 4UL );

         if( blaze::isfinite( mat ) != false ) {
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
      test_ = "Column-major isfinite()";

      // isfinite with 0x0 matrix
      {
         blaze::CompressedMatrix<float,blaze::columnMajor> mat;

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

      // isfinite with empty 3x5 matrix (finite)
      {
         blaze::CompressedMatrix<float,blaze::columnMajor> mat( 3UL, 5UL );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 5UL );
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

      // isfinite with 4x2 matrix (finite)
      {
         blaze::CompressedMatrix<float,blaze::columnMajor> mat( 4UL, 2UL );
         mat(1,1) =  1.0F;
         mat(2,0) = -2.0F;
         mat(2,1) =  3.0F;
         mat(3,0) =  4.0F;

         checkRows    ( mat, 4UL );
         checkColumns ( mat, 2UL );
         checkNonZeros( mat, 4UL );

         if( blaze::isfinite( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isfinite evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isfinite with 4x2 matrix (NAN)
      {
         blaze::CompressedMatrix<float,blaze::columnMajor> mat( 4UL, 2UL );
         mat(1,1) =  1.0F;
         mat(2,0) = -2.0F;
         mat(2,1) =   NAN;
         mat(3,0) =  4.0F;

         checkRows    ( mat, 4UL );
         checkColumns ( mat, 2UL );
         checkNonZeros( mat, 4UL );

         if( blaze::isfinite( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isinf evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isfinite with 4x2 matrix (HUGE_VALF)
      {
         blaze::CompressedMatrix<float,blaze::columnMajor> mat( 4UL, 2UL );
         mat(1,1) =      1.0F;
         mat(2,0) =     -2.0F;
         mat(2,1) = HUGE_VALF;
         mat(3,0) =      4.0F;

         checkRows    ( mat, 4UL );
         checkColumns ( mat, 2UL );
         checkNonZeros( mat, 4UL );

         if( blaze::isfinite( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isinf evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c isSquare() function for sparse matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isSquare() function for sparse matrices. In case an
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL );

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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 3UL, 0 );

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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );

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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 2UL, 0 );

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
/*!\brief Test of the \c isSymmetric() function for sparse matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isSymmetric() function for sparse matrices. In case an
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 3UL );

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 5UL );
         mat(0,0) = 1;
         mat(1,1) = 2;
         mat(2,2) = 3;
         mat.insert( 1UL, 2UL, 0 );
         mat.insert( 2UL, 0UL, 0 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 5UL );
         checkNonZeros( mat, 5UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 2UL );
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

      // Non-symmetric matrix (additional element in the lower part)
      {
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 3UL );
         mat(0,0) = 1;
         mat(1,1) = 2;
         mat(2,0) = 4;
         mat(2,2) = 3;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 4UL );
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 3UL );
         mat(0,0) = 1;
         mat(0,2) = 4;
         mat(1,1) = 2;
         mat(2,2) = 3;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 4UL );
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 5UL );
         mat(0,0) = 1;
         mat(0,2) = 4;
         mat(1,1) = 2;
         mat(2,0) = 4;
         mat(2,2) = 3;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 5UL );
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 2UL, 3UL );

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 5UL );
         mat(0,0) = 1;
         mat(1,1) = 2;
         mat(2,2) = 3;
         mat.insert( 1UL, 2UL, 0 );
         mat.insert( 2UL, 0UL, 0 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 5UL );
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

      // Non-symmetric matrix (additional element in the lower part)
      {
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 4UL );
         mat(0,0) = 1;
         mat(1,1) = 2;
         mat(2,0) = 4;
         mat(2,2) = 3;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 4UL );
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 4UL );
         mat(0,0) = 1;
         mat(0,2) = 4;
         mat(1,1) = 2;
         mat(2,2) = 3;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 4UL );
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 5UL );
         mat(0,0) = 1;
         mat(0,2) = 4;
         mat(1,1) = 2;
         mat(2,0) = 4;
         mat(2,2) = 3;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 5UL );
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
/*!\brief Test of the \c isHermitian() function for sparse matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isHermitian() function for sparse matrices. In case an
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
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 2UL, 3UL );

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
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
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
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
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL, 1UL );
         mat(1,1).imag( 1 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
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
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL, 4UL );
         mat(0,0).real( 1 );
         mat(1,1).real( 2 );
         mat(2,0).real( 4 );
         mat(2,2).real( 3 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
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
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL, 4UL );
         mat(0,0).real( 1 );
         mat(0,2).real( 4 );
         mat(1,1).real( 2 );
         mat(2,2).real( 3 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
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
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL, 5UL );
         mat(0,0).real( 1 );
         mat(0,2).imag( 4 );
         mat(1,1).real( 2 );
         mat(2,0).imag( 4 );
         mat(2,2).real( 3 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
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
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL, 7UL );
         mat(0,0).real(  1 );
         mat(0,2).imag(  4 );
         mat(1,1).real(  2 );
         mat(2,0).imag( -4 );
         mat(2,2).real(  3 );
         mat.insert( 0UL, 1UL, 0 );
         mat.insert( 1UL, 0UL, 0 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkNonZeros( mat, 7UL );
         checkNonZeros( mat, 0UL, 3UL );
         checkNonZeros( mat, 1UL, 2UL );
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
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 2UL, 3UL );

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
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
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
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
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL, 1UL );
         mat(1,1).imag( 1 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
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
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL, 4UL );
         mat(0,0).real( 1 );
         mat(1,1).real( 2 );
         mat(2,0).real( 4 );
         mat(2,2).real( 3 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
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
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL, 4UL );
         mat(0,0).real( 1 );
         mat(0,2).real( 4 );
         mat(1,1).real( 2 );
         mat(2,2).real( 3 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
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
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL, 5UL );
         mat(0,0).real( 1 );
         mat(0,2).imag( 4 );
         mat(1,1).real( 2 );
         mat(2,0).imag( 4 );
         mat(2,2).real( 3 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
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
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL, 5UL );
         mat(0,0).real(  1 );
         mat(0,2).imag(  4 );
         mat(1,1).real(  2 );
         mat(2,0).imag( -4 );
         mat(2,2).real(  3 );
         mat.insert( 0UL, 1UL, 0 );
         mat.insert( 1UL, 0UL, 0 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkNonZeros( mat, 7UL );
         checkNonZeros( mat, 0UL, 3UL );
         checkNonZeros( mat, 1UL, 2UL );
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
/*!\brief Test of the \c isUniform() function for sparse matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isUniform() function for sparse matrices. In case an
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 0UL, 3UL );

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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 0UL );

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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 1UL, 3UL, 3UL );
         mat(0,0) = 5;
         mat(0,1) = 5;
         mat(0,2) = 5;

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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 1UL, 3UL );
         mat(0,0) = 5;
         mat(1,0) = 5;
         mat(2,0) = 5;

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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 5UL, 3UL );
         mat.insert( 0UL, 1UL, 0 );
         mat.insert( 1UL, 2UL, 0 );
         mat.insert( 2UL, 0UL, 0 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 5UL );
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

      // Uniform matrix (5x3)
      {
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 5UL, 3UL, 5UL );
         mat.insert( 0UL, 1UL, 0 );
         mat.insert( 1UL, 2UL, 0 );
         mat.insert( 2UL, 0UL, 0 );
         mat.insert( 3UL, 1UL, 0 );
         mat.insert( 4UL, 2UL, 0 );

         checkRows    ( mat, 5UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 5UL );
         checkNonZeros( mat, 5UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 1UL );
         checkNonZeros( mat, 3UL, 1UL );
         checkNonZeros( mat, 4UL, 1UL );

         if( isUniform( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUniform evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-uniform matrix (3x3, 3 non-zero elements)
      {
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 3UL );
         mat.insert( 0UL, 1UL, 0 );
         mat.insert( 1UL, 0UL, 0 );
         mat.insert( 2UL, 2UL, 3 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 3UL );
         checkNonZeros( mat, 3UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( isUniform( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUniform evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-uniform matrix (3x3, 9 non-zero elements)
      {
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 9UL );
         for( size_t i=0UL; i<3UL; ++i )
            for( size_t j=0UL; j<3UL; ++j )
               mat.insert( i, j, 0UL );
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 0UL, 3UL );

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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 0UL );

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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 1UL, 3UL, 3UL );
         mat(0,0) = 5;
         mat(0,1) = 5;
         mat(0,2) = 5;

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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 1UL, 3UL );
         mat(0,0) = 5;
         mat(1,0) = 5;
         mat(2,0) = 5;

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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 5UL, 5UL );
         mat.insert( 0UL, 0UL, 0 );
         mat.insert( 2UL, 1UL, 0 );
         mat.insert( 1UL, 2UL, 0 );
         mat.insert( 2UL, 3UL, 0 );
         mat.insert( 0UL, 4UL, 0 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 5UL );
         checkCapacity( mat, 5UL );
         checkNonZeros( mat, 5UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 1UL );
         checkNonZeros( mat, 3UL, 1UL );
         checkNonZeros( mat, 4UL, 1UL );

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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 5UL, 3UL, 3UL );
         mat.insert( 1UL, 0UL, 0 );
         mat.insert( 2UL, 1UL, 0 );
         mat.insert( 0UL, 2UL, 0 );

         checkRows    ( mat, 5UL );
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

      // Non-uniform matrix (3x3, 3 non-zero elements)
      {
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 3UL );
         mat.insert( 1UL, 0UL, 0 );
         mat.insert( 0UL, 1UL, 0 );
         mat.insert( 2UL, 2UL, 3 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 3UL );
         checkNonZeros( mat, 3UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( isUniform( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUniform evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-uniform matrix (3x3, 9 non-zero elements)
      {
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 9UL );
         for( size_t i=0UL; i<3UL; ++i )
            for( size_t j=0UL; j<3UL; ++j )
               mat.insert( i, j, 0UL );
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
/*!\brief Test of the \c isZero() function for sparse matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isZero() function for sparse matrices. In case an
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 0UL, 3UL );

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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 0UL );

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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 1UL, 3UL, 1UL );
         mat.insert( 0UL, 0UL, 0 );
         mat.insert( 0UL, 1UL, 0 );
         mat.insert( 0UL, 2UL, 0 );

         checkRows    ( mat, 1UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 3UL );
         checkNonZeros( mat, 3UL );
         checkNonZeros( mat, 0UL, 3UL );

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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 1UL, 1UL );
         mat.insert( 0UL, 0UL, 0 );
         mat.insert( 1UL, 0UL, 0 );
         mat.insert( 2UL, 0UL, 0 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 1UL );
         checkCapacity( mat, 3UL );
         checkNonZeros( mat, 3UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 1UL );

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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 5UL, 3UL );
         mat.insert( 0UL, 1UL, 0 );
         mat.insert( 1UL, 2UL, 0 );
         mat.insert( 2UL, 0UL, 0 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 5UL );
         checkCapacity( mat, 3UL );
         checkNonZeros( mat, 3UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 1UL );

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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 5UL, 3UL, 5UL );
         mat.insert( 0UL, 1UL, 0 );
         mat.insert( 1UL, 2UL, 0 );
         mat.insert( 2UL, 0UL, 0 );
         mat.insert( 3UL, 1UL, 0 );
         mat.insert( 4UL, 2UL, 0 );

         checkRows    ( mat, 5UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 5UL );
         checkNonZeros( mat, 5UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 1UL );
         checkNonZeros( mat, 3UL, 1UL );
         checkNonZeros( mat, 4UL, 1UL );

         if( isZero( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isZero evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-zero matrix (3x3, 3 non-zero elements)
      {
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 3UL );
         mat.insert( 0UL, 1UL, 0 );
         mat.insert( 1UL, 0UL, 0 );
         mat.insert( 2UL, 2UL, 3 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 3UL );
         checkNonZeros( mat, 3UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 1UL );
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

      // Non-zero matrix (3x3, 9 non-zero elements)
      {
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 9UL );
         for( size_t i=0UL; i<3UL; ++i )
            for( size_t j=0UL; j<3UL; ++j )
               mat.insert( i, j, 0UL );
         mat(2,2) = 3;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 9UL );
         checkNonZeros( mat, 0UL, 3UL );
         checkNonZeros( mat, 1UL, 3UL );
         checkNonZeros( mat, 2UL, 3UL );

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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 0UL, 3UL );

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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 0UL );

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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 1UL, 3UL, 3UL );
         mat.insert( 0UL, 0UL, 0 );
         mat.insert( 0UL, 1UL, 0 );
         mat.insert( 0UL, 2UL, 0 );

         checkRows    ( mat, 1UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 3UL );
         checkNonZeros( mat, 3UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 1UL );

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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 1UL, 3UL );
         mat.insert( 0UL, 0UL, 0 );
         mat.insert( 1UL, 0UL, 0 );
         mat.insert( 2UL, 0UL, 0 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 1UL );
         checkCapacity( mat, 3UL );
         checkNonZeros( mat, 3UL );
         checkNonZeros( mat, 0UL, 3UL );

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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 5UL, 5UL );
         mat.insert( 0UL, 0UL, 0 );
         mat.insert( 2UL, 1UL, 0 );
         mat.insert( 1UL, 2UL, 0 );
         mat.insert( 2UL, 3UL, 0 );
         mat.insert( 0UL, 4UL, 0 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 5UL );
         checkCapacity( mat, 5UL );
         checkNonZeros( mat, 5UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 1UL );
         checkNonZeros( mat, 3UL, 1UL );
         checkNonZeros( mat, 4UL, 1UL );

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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 5UL, 3UL, 3UL );
         mat.insert( 1UL, 0UL, 0 );
         mat.insert( 2UL, 1UL, 0 );
         mat.insert( 0UL, 2UL, 0 );

         checkRows    ( mat, 5UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 3UL );
         checkNonZeros( mat, 3UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( isZero( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isZero evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-zero matrix (3x3, 3 non-zero elements)
      {
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 3UL );
         mat.insert( 1UL, 0UL, 0 );
         mat.insert( 0UL, 1UL, 0 );
         mat.insert( 2UL, 2UL, 3 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 3UL );
         checkNonZeros( mat, 3UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 1UL );
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

      // Non-zero matrix (3x3, 9 non-zero elements)
      {
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 9UL );
         for( size_t i=0UL; i<3UL; ++i )
            for( size_t j=0UL; j<3UL; ++j )
               mat.insert( i, j, 0UL );
         mat(2,2) = 3;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 9UL );
         checkNonZeros( mat, 9UL );
         checkNonZeros( mat, 0UL, 3UL );
         checkNonZeros( mat, 1UL, 3UL );
         checkNonZeros( mat, 2UL, 3UL );

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
/*!\brief Test of the \c isLower() function for sparse matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isLower() function for sparse matrices. In case an
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 3UL );

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 5UL );
         mat(0,0) = 1;
         mat(1,1) = 2;
         mat(2,2) = 3;
         mat.insert( 1UL, 2UL, 0 );
         mat.insert( 2UL, 0UL, 0 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 5UL );
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

      // Non-lower triangular matrix
      {
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 6UL );
         mat(0,0) = 1;
         mat(0,2) = 2;
         mat(1,0) = 3;
         mat(1,1) = 4;
         mat(2,2) = 5;
         mat(2,0) = 6;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 6UL );
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 5UL );
         mat(0,0) = 1;
         mat(1,0) = 2;
         mat(1,1) = 3;
         mat(2,2) = 4;
         mat(2,0) = 5;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 5UL );
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 2UL );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 2UL );
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 5UL );
         mat(0,0) = 1;
         mat(1,1) = 2;
         mat(2,2) = 3;
         mat.insert( 1UL, 2UL, 0 );
         mat.insert( 2UL, 0UL, 0 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 5UL );
         checkNonZeros( mat, 5UL );
         checkNonZeros( mat, 0UL, 2UL );
         checkNonZeros( mat, 1UL, 1UL );
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

      // Non-lower triangle matrix
      {
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 6UL );
         mat(0,0) = 1;
         mat(0,2) = 2;
         mat(1,0) = 3;
         mat(1,1) = 4;
         mat(2,2) = 5;
         mat(2,0) = 6;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 6UL );
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 5UL );
         mat(0,0) = 1;
         mat(1,0) = 2;
         mat(1,1) = 3;
         mat(2,2) = 4;
         mat(2,0) = 5;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 5UL );
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
/*!\brief Test of the \c isUniLower() function for sparse matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isUniLower() function for sparse matrices. In case an
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 3UL );

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 3UL );
         mat(0,0) = 1;
         mat(1,1) = 1;
         mat(2,2) = 1;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 3UL );
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 5UL );
         mat(0,0) = 1;
         mat(1,1) = 2;
         mat(2,2) = 3;
         mat.insert( 1UL, 2UL, 0 );
         mat.insert( 2UL, 0UL, 0 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 5UL );
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

      // Lower unitriangular matrix
      {
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 5UL );
         mat(0,0) = 1;
         mat(1,0) = 2;
         mat(1,1) = 1;
         mat(2,2) = 1;
         mat(2,0) = 3;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 5UL );
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 5UL );
         mat(0,0) = 1;
         mat(1,0) = 2;
         mat(1,1) = 3;
         mat(2,2) = 4;
         mat(2,0) = 5;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 5UL );
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 6UL );
         mat(0,0) = 1;
         mat(0,2) = 2;
         mat(1,0) = 3;
         mat(1,1) = 1;
         mat(2,2) = 1;
         mat(2,0) = 4;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 6UL );
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 2UL, 3UL );

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 3UL );
         mat(0,0) = 1;
         mat(1,1) = 1;
         mat(2,2) = 1;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 3UL );
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 5UL );
         mat(0,0) = 1;
         mat(1,1) = 2;
         mat(2,2) = 3;
         mat.insert( 1UL, 2UL, 0 );
         mat.insert( 2UL, 0UL, 0 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 5UL );
         checkNonZeros( mat, 5UL );
         checkNonZeros( mat, 0UL, 2UL );
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

      // Lower unitriangular matrix
      {
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 5UL );
         mat(0,0) = 1;
         mat(1,0) = 2;
         mat(1,1) = 1;
         mat(2,2) = 1;
         mat(2,0) = 3;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 5UL );
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 5UL );
         mat(0,0) = 1;
         mat(1,0) = 2;
         mat(1,1) = 3;
         mat(2,2) = 4;
         mat(2,0) = 5;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 5UL );
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 6UL );
         mat(0,0) = 1;
         mat(0,2) = 2;
         mat(1,0) = 3;
         mat(1,1) = 1;
         mat(2,2) = 1;
         mat(2,0) = 4;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 6UL );
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
/*!\brief Test of the \c isStrictlyLower() function for sparse matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isStrictlyLower() function for sparse matrices. In
// case an error is detected, a \a std::runtime_error exception is thrown.
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 3UL );

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 5UL );
         mat(0,0) = 1;
         mat(1,1) = 2;
         mat(2,2) = 3;
         mat.insert( 1UL, 2UL, 0 );
         mat.insert( 2UL, 0UL, 0 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 5UL );
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

      // Strictly lower triangular matrix
      {
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 2UL );
         mat(1,0) = 2;
         mat(2,0) = 5;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 2UL );
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 5UL );
         mat(0,0) = 1;
         mat(1,0) = 2;
         mat(1,1) = 3;
         mat(2,2) = 4;
         mat(2,0) = 5;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 5UL );
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 3UL );
         mat(0,2) = 2;
         mat(1,0) = 3;
         mat(2,0) = 4;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 3UL );
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 2UL, 3UL );

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 5UL );
         mat(0,0) = 1;
         mat(1,1) = 2;
         mat(2,2) = 3;
         mat.insert( 1UL, 2UL, 0 );
         mat.insert( 2UL, 0UL, 0 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 5UL );
         checkNonZeros( mat, 5UL );
         checkNonZeros( mat, 0UL, 2UL );
         checkNonZeros( mat, 1UL, 1UL );
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

      // Strictly lower triangular matrix
      {
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 2UL );
         mat(1,0) = 2;
         mat(2,0) = 5;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 2UL );
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 5UL );
         mat(0,0) = 1;
         mat(1,0) = 2;
         mat(1,1) = 3;
         mat(2,2) = 4;
         mat(2,0) = 5;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 5UL );
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 3UL );
         mat(0,2) = 2;
         mat(1,0) = 3;
         mat(2,0) = 4;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 3UL );
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
/*!\brief Test of the \c isUpper() function for sparse matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isUpper() function for sparse matrices. In case an
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 3UL );

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 5UL );
         mat(0,0) = 1;
         mat(1,1) = 2;
         mat(2,2) = 3;
         mat.insert( 1UL, 2UL, 0 );
         mat.insert( 2UL, 0UL, 0 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 5UL );
         checkNonZeros( mat, 5UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 2UL );
         checkNonZeros( mat, 2UL, 2UL );

         if( isUpper( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUpper evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-upper triangle matrix
      {
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 6UL );
         mat(0,0) = 1;
         mat(0,2) = 2;
         mat(1,1) = 3;
         mat(1,2) = 4;
         mat(2,0) = 5;
         mat(2,2) = 6;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 6UL );
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 5UL );
         mat(0,0) = 1;
         mat(0,2) = 2;
         mat(1,1) = 3;
         mat(1,2) = 4;
         mat(2,2) = 5;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 5UL );
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 2UL );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 2UL );
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 5UL );
         mat(0,0) = 1;
         mat(1,1) = 2;
         mat(2,2) = 3;
         mat.insert( 1UL, 2UL, 0 );
         mat.insert( 2UL, 0UL, 0 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 5UL );
         checkNonZeros( mat, 5UL );
         checkNonZeros( mat, 0UL, 2UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 2UL );

         if( isUpper( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUpper evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-upper triangle matrix
      {
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 6UL );
         mat(0,0) = 1;
         mat(0,2) = 2;
         mat(1,1) = 3;
         mat(1,2) = 4;
         mat(2,0) = 5;
         mat(2,2) = 6;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 6UL );
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 5UL );
         mat(0,0) = 1;
         mat(0,2) = 2;
         mat(1,1) = 3;
         mat(1,2) = 4;
         mat(2,2) = 5;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 5UL );
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
/*!\brief Test of the \c isUniUpper() function for sparse matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isUniUpper() function for sparse matrices. In case an
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 3UL );

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 3UL );
         mat(0,0) = 1;
         mat(1,1) = 1;
         mat(2,2) = 1;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 3UL );
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 5UL );
         mat(0,0) = 1;
         mat(1,1) = 2;
         mat(2,2) = 3;
         mat.insert( 1UL, 2UL, 0 );
         mat.insert( 2UL, 0UL, 0 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 5UL );
         checkNonZeros( mat, 5UL );
         checkNonZeros( mat, 0UL, 1UL );
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

      // Upper unitriangular matrix
      {
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 5UL );
         mat(0,0) = 1;
         mat(0,2) = 2;
         mat(1,1) = 1;
         mat(1,2) = 3;
         mat(2,2) = 1;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 5UL );
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 5UL );
         mat(0,0) = 1;
         mat(0,2) = 2;
         mat(1,1) = 3;
         mat(1,2) = 4;
         mat(2,2) = 5;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 5UL );
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 6UL );
         mat(0,0) = 1;
         mat(0,2) = 2;
         mat(1,1) = 1;
         mat(1,2) = 3;
         mat(2,0) = 4;
         mat(2,2) = 1;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 6UL );
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 2UL, 3UL );

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 3UL );
         mat(0,0) = 1;
         mat(1,1) = 1;
         mat(2,2) = 1;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 3UL );
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 5UL );
         mat(0,0) = 1;
         mat(1,1) = 2;
         mat(2,2) = 3;
         mat.insert( 1UL, 2UL, 0 );
         mat.insert( 2UL, 0UL, 0 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 5UL );
         checkNonZeros( mat, 5UL );
         checkNonZeros( mat, 0UL, 2UL );
         checkNonZeros( mat, 1UL, 1UL );
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

      // Upper unitriangular matrix
      {
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 5UL );
         mat(0,0) = 1;
         mat(0,2) = 2;
         mat(1,1) = 1;
         mat(1,2) = 3;
         mat(2,2) = 1;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 5UL );
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 5UL );
         mat(0,0) = 1;
         mat(0,2) = 2;
         mat(1,1) = 3;
         mat(1,2) = 4;
         mat(2,2) = 5;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 5UL );
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 6UL );
         mat(0,0) = 1;
         mat(0,2) = 2;
         mat(1,1) = 1;
         mat(1,2) = 3;
         mat(2,0) = 4;
         mat(2,2) = 1;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 6UL );
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
/*!\brief Test of the \c isStrictlyUpper() function for sparse matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isStrictlyUpper() function for sparse matrices. In
// case an error is detected, a \a std::runtime_error exception is thrown.
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 3UL );

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 5UL );
         mat(0,0) = 1;
         mat(1,1) = 2;
         mat(2,2) = 3;
         mat.insert( 1UL, 2UL, 0 );
         mat.insert( 2UL, 0UL, 0 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 5UL );
         checkNonZeros( mat, 5UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 2UL );
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

      // Strictly upper triangular matrix
      {
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 2UL );
         mat(0,2) = 2;
         mat(1,2) = 4;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 2UL );
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 5UL );
         mat(0,0) = 1;
         mat(0,2) = 2;
         mat(1,1) = 3;
         mat(1,2) = 4;
         mat(2,2) = 5;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 5UL );
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 3UL );
         mat(0,2) = 2;
         mat(1,2) = 3;
         mat(2,0) = 4;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 3UL );
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 2UL, 3UL );

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 5UL );
         mat(0,0) = 1;
         mat(1,1) = 2;
         mat(2,2) = 3;
         mat.insert( 1UL, 2UL, 0 );
         mat.insert( 2UL, 0UL, 0 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 5UL );
         checkNonZeros( mat, 5UL );
         checkNonZeros( mat, 0UL, 2UL );
         checkNonZeros( mat, 1UL, 1UL );
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

      // Strictly upper triangular matrix
      {
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 2UL );
         mat(0,2) = 2;
         mat(1,2) = 4;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 2UL );
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 5UL );
         mat(0,0) = 1;
         mat(0,2) = 2;
         mat(1,1) = 3;
         mat(1,2) = 4;
         mat(2,2) = 5;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 5UL );
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 3UL );
         mat(0,2) = 2;
         mat(1,2) = 3;
         mat(2,0) = 4;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 3UL );
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
/*!\brief Test of the \c isDiagonal() function for sparse matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isDiagonal() function for sparse matrices. In case
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 3UL );

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 5UL );
         mat(0,0) = 1;
         mat(1,1) = 2;
         mat(2,2) = 3;
         mat.insert( 1UL, 2UL, 0 );
         mat.insert( 2UL, 0UL, 0 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 5UL );
         checkNonZeros( mat, 5UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 2UL );
         checkNonZeros( mat, 2UL, 2UL );

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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 4UL );
         mat(0,0) = 1;
         mat(1,1) = 2;
         mat(2,0) = 4;
         mat(2,2) = 3;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 4UL );
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 4UL );
         mat(0,0) = 1;
         mat(0,2) = 4;
         mat(1,1) = 2;
         mat(2,2) = 3;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 4UL );
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 2UL, 3UL );

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 5UL );
         mat(0,0) = 1;
         mat(1,1) = 2;
         mat(2,2) = 3;
         mat.insert( 1UL, 2UL, 0 );
         mat.insert( 2UL, 0UL, 0 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 5UL );
         checkNonZeros( mat, 5UL );
         checkNonZeros( mat, 0UL, 2UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 2UL );

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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 4UL );
         mat(0,0) = 1;
         mat(1,1) = 2;
         mat(2,0) = 4;
         mat(2,2) = 3;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 4UL );
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 4UL );
         mat(0,0) = 1;
         mat(0,2) = 4;
         mat(1,1) = 2;
         mat(2,2) = 3;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 4UL );
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
/*!\brief Test of the \c isIdentity() function for sparse matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isIdentity() function for sparse matrices. In case
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 3UL );

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 5UL );
         mat(0,0) = 1;
         mat(1,1) = 1;
         mat(2,2) = 1;
         mat.insert( 1UL, 2UL, 0 );
         mat.insert( 2UL, 0UL, 0 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 5UL );
         checkNonZeros( mat, 5UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 2UL );
         checkNonZeros( mat, 2UL, 2UL );

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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 2UL );
         mat(0,0) = 1;
         mat(2,2) = 1;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 2UL );
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 3UL );
         mat(0,0) = 1;
         mat(1,1) = 2;
         mat(2,2) = 3;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 3UL );
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 4UL );
         mat(0,0) = 1;
         mat(1,1) = 1;
         mat(2,0) = 2;
         mat(2,2) = 1;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 4UL );
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 4UL );
         mat(0,0) = 1;
         mat(0,2) = 2;
         mat(1,1) = 1;
         mat(2,2) = 1;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 4UL );
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 2UL, 3UL );

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 5UL );
         mat(0,0) = 1;
         mat(1,1) = 1;
         mat(2,2) = 1;
         mat.insert( 1UL, 2UL, 0 );
         mat.insert( 2UL, 0UL, 0 );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 5UL );
         checkNonZeros( mat, 5UL );
         checkNonZeros( mat, 0UL, 2UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 2UL );

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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 2UL );
         mat(0,0) = 1;
         mat(2,2) = 1;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 2UL );
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 3UL );
         mat(0,0) = 1;
         mat(1,1) = 2;
         mat(2,2) = 3;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 3UL );
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 4UL );
         mat(0,0) = 1;
         mat(1,1) = 1;
         mat(2,0) = 2;
         mat(2,2) = 1;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 4UL );
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 4UL );
         mat(0,0) = 1;
         mat(0,2) = 2;
         mat(1,1) = 1;
         mat(2,2) = 1;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 4UL );
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
/*!\brief Test of the \c min() function for sparse matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c min() function for sparse matrices template. In case an
// error is detected, a \a std::runtime_error exception is thrown.
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 2UL, 5UL );
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 3UL, 5UL );
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 5UL, 3UL, 5UL );
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 5UL, 5UL );
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 5UL );
         mat(0,0) = 1;
         mat(0,2) = 2;
         mat(1,1) = 3;
         mat(2,0) = 4;
         mat(2,2) = 5;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkNonZeros( mat, 5UL );

         const int minimum = min( mat );

         if( minimum != 1 ) {
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 5UL, 3UL, 5UL );
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 5UL, 5UL );
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 5UL, 3UL, 5UL );
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 5UL, 5UL );
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 5UL );
         mat(0,0) = 1;
         mat(0,2) = 2;
         mat(1,1) = 3;
         mat(2,0) = 4;
         mat(2,2) = 5;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkNonZeros( mat, 5UL );

         const int minimum = min( mat );

         if( minimum != 1 ) {
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
/*!\brief Test of the \c max() function for sparse matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c max() function for sparse matrices template. In case an
// error is detected, a \a std::runtime_error exception is thrown.
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 2UL, 5UL );
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 3UL, 5UL );
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 5UL, 3UL, 5UL );
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 5UL, 5UL );
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 5UL );
         mat(0,0) = -1;
         mat(0,2) = -2;
         mat(1,1) = -3;
         mat(2,0) = -4;
         mat(2,2) = -5;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkNonZeros( mat, 5UL );

         const int maximum = max( mat );

         if( maximum != -1 ) {
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 2UL, 5UL );
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 2UL, 3UL, 5UL );
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 5UL, 3UL, 5UL );
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 5UL, 5UL );
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 5UL );
         mat(0,0) = -1;
         mat(0,2) = -2;
         mat(1,1) = -3;
         mat(2,0) = -4;
         mat(2,2) = -5;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkNonZeros( mat, 5UL );

         const int maximum = max( mat );

         if( maximum != -1 ) {
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
/*!\brief Test of the \c trace() function for sparse matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c trace() function for sparse matrices template. In case
// an error is detected, a \a std::runtime_error exception is thrown.
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat;

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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 6UL );
         mat(0,0) = -1;
         mat(0,2) = -3;
         mat(1,1) = -5;
         mat(1,2) =  6;
         mat(2,0) =  7;
         mat(2,2) = -9;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkNonZeros( mat, 6UL );

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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 3UL );

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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat;

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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 6UL );
         mat(0,0) = -1;
         mat(0,2) = -3;
         mat(1,1) = -5;
         mat(1,2) =  6;
         mat(2,0) =  7;
         mat(2,2) = -9;

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkNonZeros( mat, 6UL );

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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 2UL, 3UL );

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );

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
/*!\brief Test of the \c l1Norm() function for sparse matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c l1Norm() function for sparse matrices template. In
// case an error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testL1Norm()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "l1Norm() function";

      {
         blaze::CompressedMatrix<int,blaze::rowMajor> mat;

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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 7UL );

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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat{ { 0,  0,  1,  0,  1,  0,  0 },
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat;

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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 7UL, 0 );

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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat{ { 0,  0,  0 },
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
/*!\brief Test of the \c l2Norm() function for sparse matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c l2Norm() function for sparse matrices template. In
// case an error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testL2Norm()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "l2Norm() function";

      {
         blaze::CompressedMatrix<int,blaze::rowMajor> mat;

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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 7UL );

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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat{ { 0,  0,  1,  0,  1, -2,  0 },
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat;

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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 7UL );

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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat{ {  0,  0,  0 },
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
/*!\brief Test of the \c l3Norm() function for sparse matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c l3Norm() function for sparse matrices template. In
// case an error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testL3Norm()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "l3Norm() function";

      {
         blaze::CompressedMatrix<int,blaze::rowMajor> mat;

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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 7UL );

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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat{ { 0,  0,  1,  0,  1, -2,  0 },
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat;

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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 7UL );

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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat{ {  0,  0,  0 },
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
/*!\brief Test of the \c l4Norm() function for sparse matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c l4Norm() function for sparse matrices template. In
// case an error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testL4Norm()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "l4Norm() function";

      {
         blaze::CompressedMatrix<int,blaze::rowMajor> mat;

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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 7UL );

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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat{ { 0,  0,  2,  0,  2, -2,  0 },
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat;

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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 7UL );

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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat{ { 0,  0,  2,  0,  2, -2,  0 },
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
/*!\brief Test of the \c lpNorm() function for sparse matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c lpNorm() function for sparse matrices template. In
// case an error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testLpNorm()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "lpNorm() function";

      {
         blaze::CompressedMatrix<int,blaze::rowMajor> mat;

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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 7UL );

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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 5UL, 10UL );
         randomize( mat, 20UL, -5, 5 );

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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 5UL, 10UL );
         randomize( mat, 20UL, -5, 5 );

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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 5UL, 10UL );
         randomize( mat, 20UL, -5, 5 );

         const double norm1( blaze::lpNorm( mat, 3 ) );
         const double norm2( blaze::lpNorm<3UL>( mat ) );
         const double norm3( blaze::l3Norm( mat ) );

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

      {
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 5UL, 10UL );
         randomize( mat, 20UL, -5, 5 );

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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat;

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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 7UL );

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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 5UL, 10UL );
         randomize( mat, 20UL, -5, 5 );

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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 5UL, 10UL );
         randomize( mat, 20UL, -5, 5 );

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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 5UL, 10UL );
         randomize( mat, 20UL, -5, 5 );

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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 5UL, 10UL );
         randomize( mat, 20UL, -5, 5 );

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
/*!\brief Test of the \c linfNorm() function for sparse matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c linfNorm() function for sparse matrices template. In
// case an error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testLinfNorm()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "linfNorm() function";

      {
         blaze::CompressedMatrix<int,blaze::rowMajor> mat;

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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 7UL );

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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 5UL, 10UL );
         randomize( mat, 20UL, -5, 5 );

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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat;

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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 7UL );

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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 5UL, 10UL );
         randomize( mat, 20UL, -5, 5 );

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
/*!\brief Test of the \c mean() function for sparse matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c mean() function for sparse matrices. In case an
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL );

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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat{ { 1, 3, 2 },
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 0UL );

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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 0UL, 3UL );

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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL );

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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat{ { 1, 3, 2 },
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 0UL );

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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL );

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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat{ { 1, 3, 2 },
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 0UL, 3UL );

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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL );

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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat{ { 1, 3, 2 },
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 0UL );

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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 0UL, 3UL );

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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL );

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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat{ { 1, 3, 2 },
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 0UL );

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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL );

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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat{ { 1, 3, 2 },
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 0UL, 3UL );

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
/*!\brief Test of the \c var() function for sparse matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c var() function for sparse matrices. In case an
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL );

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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat{ { 1, 3, 2 },
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 0UL );

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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 0UL, 3UL );

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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 1UL, 1UL );

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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL );

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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat{ { 1, 3, 2 },
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 1UL );

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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL );

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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat{ { 1, 3, 2 },
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 0UL, 3UL );

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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 1UL, 3UL );

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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL );

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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat{ { 1, 3, 2 },
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 0UL );

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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 0UL, 3UL );

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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 1UL, 1UL );

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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL );

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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat{ { 1, 3, 2 },
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 0UL );

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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 1UL );

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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL );

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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat{ { 1, 3, 2 },
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 0UL, 3UL );

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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 1UL, 3UL );

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
/*!\brief Test of the \c stddev() function for sparse matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c stddev() function for sparse matrices. In case an
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL );

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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat{ { 1, 3, 2 },
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 0UL );

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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 0UL, 3UL );

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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 1UL, 1UL );

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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL );

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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat{ { 1, 3, 2 },
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 0UL );

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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 1UL );

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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL );

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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat{ { 1, 3, 2 },
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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 0UL, 3UL );

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
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 1UL, 3UL );

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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL );

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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat{ { 1, 3, 2 },
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 0UL );

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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 0UL, 3UL );

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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 1UL, 1UL );

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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL );

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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat{ { 1, 3, 2 },
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 0UL );

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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 1UL );

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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL );

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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat{ { 1, 3, 2 },
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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 0UL, 3UL );

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
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 1UL, 3UL );

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

} // namespace sparsematrix

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
   std::cout << "   Running general SparseMatrix operation test..." << std::endl;

   try
   {
      RUN_SPARSEMATRIX_GENERAL_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during general SparseMatrix operation test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
