//=================================================================================================
/*!
//  \file src/mathtest/matrices/sparsematrix/ZeroTest.cpp
//  \brief Source file for the zero SparseMatrix operation test
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
#include <blaze/math/UniformVector.h>
#include <blaze/math/ZeroMatrix.h>
#include <blazetest/mathtest/matrices/sparsematrix/ZeroTest.h>

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
/*!\brief Constructor for the ZeroTest class test.
//
// \exception std::runtime_error Operation error detected.
*/
ZeroTest::ZeroTest()
{
   testIsNan();
   testIsInf();
   testIsFinite();
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
void ZeroTest::testIsNan()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major isnan()";

      // isnan with 0x0 zero matrix
      {
         blaze::ZeroMatrix<float,blaze::rowMajor> mat;

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

      // isnan with 3x5 zero matrix
      {
         blaze::ZeroMatrix<float,blaze::rowMajor> mat( 3UL, 5UL );

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
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major isnan()";

      // isnan with 0x0 zero matrix
      {
         blaze::ZeroMatrix<float,blaze::columnMajor> mat;

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

      // isnan with 3x3 zero matrix
      {
         blaze::ZeroMatrix<float,blaze::columnMajor> mat( 3UL, 5UL );

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
void ZeroTest::testIsInf()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major isinf()";

      // isinf with 0x0 zero matrix
      {
         blaze::ZeroMatrix<float,blaze::rowMajor> mat;

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

      // isinf with 3x3 zero matrix
      {
         blaze::ZeroMatrix<float,blaze::rowMajor> mat( 3UL, 5UL );

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
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major isinf()";

      // isinf with 0x0 zero matrix
      {
         blaze::ZeroMatrix<float,blaze::columnMajor> mat;

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

      // isinf with 3x3 zero matrix
      {
         blaze::ZeroMatrix<float,blaze::columnMajor> mat( 3UL, 5UL );

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
void ZeroTest::testIsFinite()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major isfinite()";

      // isfinite with 0x0 zero matrix
      {
         blaze::ZeroMatrix<float,blaze::rowMajor> mat;

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

      // isfinite with 3x3 zero matrix
      {
         blaze::ZeroMatrix<float,blaze::rowMajor> mat( 3UL, 5UL );

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
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major isfinite()";

      // isfinite with 0x0 zero matrix
      {
         blaze::ZeroMatrix<float,blaze::columnMajor> mat;

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

      // isfinite with 3x3 zero matrix
      {
         blaze::ZeroMatrix<float,blaze::columnMajor> mat( 3UL, 5UL );

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
void ZeroTest::testIsSymmetric()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major isSymmetric()";

      // Default zero matrix
      {
         blaze::ZeroMatrix<int,blaze::rowMajor> mat;

         checkRows    ( mat, 0UL );
         checkColumns ( mat, 0UL );
         checkNonZeros( mat, 0UL );

         if( isSymmetric( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSymmetric evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-square zero matrix
      {
         blaze::ZeroMatrix<int,blaze::rowMajor> mat( 3UL, 5UL );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 5UL );
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

      // Square zero matrix
      {
         blaze::ZeroMatrix<int,blaze::rowMajor> mat( 3UL, 3UL );

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
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major isSymmetric()";

      // Default zero matrix
      {
         blaze::ZeroMatrix<int,blaze::columnMajor> mat;

         checkRows    ( mat, 0UL );
         checkColumns ( mat, 0UL );
         checkNonZeros( mat, 0UL );

         if( isSymmetric( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSymmetric evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-square zero matrix
      {
         blaze::ZeroMatrix<int,blaze::columnMajor> mat( 5UL, 3UL );

         checkRows    ( mat, 5UL );
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

      // Square zero matrix
      {
         blaze::ZeroMatrix<int,blaze::columnMajor> mat( 3UL, 3UL );

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
void ZeroTest::testIsHermitian()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major isHermitian()";

      // Default zero matrix
      {
         blaze::ZeroMatrix<cplx,blaze::rowMajor> mat;

         checkRows    ( mat, 0UL );
         checkColumns ( mat, 0UL );
         checkNonZeros( mat, 0UL );

         if( isHermitian( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isHermitian evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-square zero matrix
      {
         blaze::ZeroMatrix<cplx,blaze::rowMajor> mat( 3UL, 5UL );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 5UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
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

      // Square zero matrix
      {
         blaze::ZeroMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL );

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
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major isHermitian()";

      // Default zero matrix
      {
         blaze::ZeroMatrix<cplx,blaze::columnMajor> mat;

         checkRows    ( mat, 0UL );
         checkColumns ( mat, 0UL );
         checkNonZeros( mat, 0UL );

         if( isHermitian( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isHermitian evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-square zero matrix
      {
         blaze::ZeroMatrix<cplx,blaze::columnMajor> mat( 5UL, 3UL );

         checkRows    ( mat, 5UL );
         checkColumns ( mat, 3UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
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

      // Square zero matrix
      {
         blaze::ZeroMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL );

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
void ZeroTest::testIsUniform()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major isUniform()";

      // Default zero matrix
      {
         blaze::ZeroMatrix<int,blaze::rowMajor> mat;

         checkRows    ( mat, 0UL );
         checkColumns ( mat, 0UL );
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

      // Non-square zero matrix
      {
         blaze::ZeroMatrix<int,blaze::rowMajor> mat( 3UL, 5UL );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 5UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 0UL );

         if( isUniform( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUniform evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Square zero matrix
      {
         blaze::ZeroMatrix<int,blaze::rowMajor> mat( 3UL, 3UL );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 0UL );

         if( isUniform( mat ) != true ) {
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

      // Default zero matrix
      {
         blaze::ZeroMatrix<int,blaze::columnMajor> mat;

         checkRows    ( mat, 0UL );
         checkColumns ( mat, 0UL );
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

      // Non-square zero matrix
      {
         blaze::ZeroMatrix<int,blaze::columnMajor> mat( 5UL, 3UL );

         checkRows    ( mat, 5UL );
         checkColumns ( mat, 3UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 0UL );

         if( isUniform( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUniform evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Square zero matrix
      {
         blaze::ZeroMatrix<int,blaze::columnMajor> mat( 3UL, 3UL );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 0UL );

         if( isUniform( mat ) != true ) {
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
void ZeroTest::testIsZero()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major isZero()";

      // Default zero matrix
      {
         blaze::ZeroMatrix<int,blaze::rowMajor> mat;

         checkRows    ( mat, 0UL );
         checkColumns ( mat, 0UL );
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

      // Non-square zero matrix
      {
         blaze::ZeroMatrix<int,blaze::rowMajor> mat( 3UL, 5UL );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 5UL );
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

      // Square zero matrix
      {
         blaze::ZeroMatrix<int,blaze::rowMajor> mat( 3UL, 3UL );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
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
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major isZero()";

      // Default zero matrix
      {
         blaze::ZeroMatrix<int,blaze::columnMajor> mat;

         checkRows    ( mat, 0UL );
         checkColumns ( mat, 0UL );
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

      // Non-square zero matrix
      {
         blaze::ZeroMatrix<int,blaze::columnMajor> mat( 5UL, 3UL );

         checkRows    ( mat, 5UL );
         checkColumns ( mat, 3UL );
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

      // Square zero matrix
      {
         blaze::ZeroMatrix<int,blaze::columnMajor> mat( 3UL, 3UL );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 3UL );
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
void ZeroTest::testIsLower()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major isLower()";

      // Default zero matrix
      {
         blaze::ZeroMatrix<int,blaze::rowMajor> mat;

         checkRows    ( mat, 0UL );
         checkColumns ( mat, 0UL );
         checkNonZeros( mat, 0UL );

         if( isLower( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isLower evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-square zero matrix
      {
         blaze::ZeroMatrix<int,blaze::rowMajor> mat( 3UL, 5UL );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 5UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 0UL );

         if( isLower( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isLower evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Square zero matrix
      {
         blaze::ZeroMatrix<int,blaze::rowMajor> mat( 3UL, 3UL );

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
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major isLower()";

      // Default zero matrix
      {
         blaze::ZeroMatrix<int,blaze::columnMajor> mat;

         checkRows    ( mat, 0UL );
         checkColumns ( mat, 0UL );
         checkNonZeros( mat, 0UL );

         if( isLower( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isLower evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-square zero matrix
      {
         blaze::ZeroMatrix<int,blaze::columnMajor> mat( 5UL, 3UL );

         checkRows    ( mat, 5UL );
         checkColumns ( mat, 3UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 0UL );

         if( isLower( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isLower evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Square zero matrix
      {
         blaze::ZeroMatrix<int,blaze::columnMajor> mat( 3UL, 3UL );

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
void ZeroTest::testIsUniLower()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major isUniLower()";

      // Default zero matrix
      {
         blaze::ZeroMatrix<int,blaze::rowMajor> mat;

         checkRows    ( mat, 0UL );
         checkColumns ( mat, 0UL );
         checkNonZeros( mat, 0UL );

         if( isUniLower( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUniLower evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-square zero matrix
      {
         blaze::ZeroMatrix<int,blaze::rowMajor> mat( 3UL, 5UL );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 5UL );
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

      // Square zero matrix
      {
         blaze::ZeroMatrix<int,blaze::rowMajor> mat( 3UL, 3UL );

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
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major isUniLower()";

      // Default zero matrix
      {
         blaze::ZeroMatrix<int,blaze::columnMajor> mat;

         checkRows    ( mat, 0UL );
         checkColumns ( mat, 0UL );
         checkNonZeros( mat, 0UL );

         if( isUniLower( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUniLower evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-square zero matrix
      {
         blaze::ZeroMatrix<int,blaze::columnMajor> mat( 5UL, 3UL );

         checkRows    ( mat, 5UL );
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

      // Square zero matrix
      {
         blaze::ZeroMatrix<int,blaze::columnMajor> mat( 3UL, 3UL );

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
void ZeroTest::testIsStrictlyLower()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major isStrictlyLower()";

      // Default zero matrix
      {
         blaze::ZeroMatrix<int,blaze::rowMajor> mat;

         checkRows    ( mat, 0UL );
         checkColumns ( mat, 0UL );
         checkNonZeros( mat, 0UL );

         if( isStrictlyLower( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isStrictlyLower evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-square zero matrix
      {
         blaze::ZeroMatrix<int,blaze::rowMajor> mat( 3UL, 5UL );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 5UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 0UL );

         if( isStrictlyLower( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isStrictlyLower evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Square zero matrix
      {
         blaze::ZeroMatrix<int,blaze::rowMajor> mat( 3UL, 3UL );

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
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major isStrictlyLower()";

      // Default zero matrix
      {
         blaze::ZeroMatrix<int,blaze::columnMajor> mat;

         checkRows    ( mat, 0UL );
         checkColumns ( mat, 0UL );
         checkNonZeros( mat, 0UL );

         if( isStrictlyLower( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isStrictlyLower evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-square zero matrix
      {
         blaze::ZeroMatrix<int,blaze::columnMajor> mat( 5UL, 3UL );

         checkRows    ( mat, 5UL );
         checkColumns ( mat, 3UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 0UL );

         if( isStrictlyLower( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isStrictlyLower evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Square zero matrix
      {
         blaze::ZeroMatrix<int,blaze::columnMajor> mat( 3UL, 3UL );

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
void ZeroTest::testIsUpper()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major isUpper()";

      // Default zero matrix
      {
         blaze::ZeroMatrix<int,blaze::rowMajor> mat;

         checkRows    ( mat, 0UL );
         checkColumns ( mat, 0UL );
         checkNonZeros( mat, 0UL );

         if( isUpper( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUpper evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-square zero matrix
      {
         blaze::ZeroMatrix<int,blaze::rowMajor> mat( 3UL, 5UL );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 5UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 0UL );

         if( isUpper( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUpper evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Square zero matrix
      {
         blaze::ZeroMatrix<int,blaze::rowMajor> mat( 3UL, 3UL );

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
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major isUpper()";

      // Default zero matrix
      {
         blaze::ZeroMatrix<int,blaze::columnMajor> mat;

         checkRows    ( mat, 0UL );
         checkColumns ( mat, 0UL );
         checkNonZeros( mat, 0UL );

         if( isUpper( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUpper evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-square zero matrix
      {
         blaze::ZeroMatrix<int,blaze::columnMajor> mat( 5UL, 3UL );

         checkRows    ( mat, 5UL );
         checkColumns ( mat, 3UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 0UL );

         if( isUpper( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUpper evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Square zero matrix
      {
         blaze::ZeroMatrix<int,blaze::columnMajor> mat( 3UL, 3UL );

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
void ZeroTest::testIsUniUpper()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major isUniUpper()";

      // Default zero matrix
      {
         blaze::ZeroMatrix<int,blaze::rowMajor> mat;

         checkRows    ( mat, 0UL );
         checkColumns ( mat, 0UL );
         checkNonZeros( mat, 0UL );

         if( isUniUpper( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUniUpper evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-square zero matrix
      {
         blaze::ZeroMatrix<int,blaze::rowMajor> mat( 3UL, 5UL );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 5UL );
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

      // Square zero matrix
      {
         blaze::ZeroMatrix<int,blaze::rowMajor> mat( 3UL, 3UL );

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
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major isUniUpper()";

      // Default zero matrix
      {
         blaze::ZeroMatrix<int,blaze::columnMajor> mat;

         checkRows    ( mat, 0UL );
         checkColumns ( mat, 0UL );
         checkNonZeros( mat, 0UL );

         if( isUniUpper( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isUniUpper evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-square zero matrix
      {
         blaze::ZeroMatrix<int,blaze::columnMajor> mat( 5UL, 3UL );

         checkRows    ( mat, 5UL );
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

      // Square zero matrix
      {
         blaze::ZeroMatrix<int,blaze::columnMajor> mat( 3UL, 3UL );

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
void ZeroTest::testIsStrictlyUpper()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major isStrictlyUpper()";

      // Default zero matrix
      {
         blaze::ZeroMatrix<int,blaze::rowMajor> mat;

         checkRows    ( mat, 0UL );
         checkColumns ( mat, 0UL );
         checkNonZeros( mat, 0UL );

         if( isStrictlyUpper( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isStrictlyUpper evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-square zero matrix
      {
         blaze::ZeroMatrix<int,blaze::rowMajor> mat( 3UL, 5UL );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 5UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 0UL );

         if( isStrictlyUpper( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isStrictlyUpper evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Square zero matrix
      {
         blaze::ZeroMatrix<int,blaze::rowMajor> mat( 3UL, 3UL );

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
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major isStrictlyUpper()";

      // Default zero matrix
      {
         blaze::ZeroMatrix<int,blaze::columnMajor> mat;

         checkRows    ( mat, 0UL );
         checkColumns ( mat, 0UL );
         checkNonZeros( mat, 0UL );

         if( isStrictlyUpper( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isStrictlyUpper evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-square zero matrix
      {
         blaze::ZeroMatrix<int,blaze::columnMajor> mat( 5UL, 3UL );

         checkRows    ( mat, 5UL );
         checkColumns ( mat, 3UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 0UL );

         if( isStrictlyUpper( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isStrictlyUpper evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Square zero matrix
      {
         blaze::ZeroMatrix<int,blaze::columnMajor> mat( 3UL, 3UL );

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
void ZeroTest::testIsDiagonal()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major isDiagonal()";

      // Default zero matrix
      {
         blaze::ZeroMatrix<int,blaze::rowMajor> mat;

         checkRows    ( mat, 0UL );
         checkColumns ( mat, 0UL );
         checkNonZeros( mat, 0UL );

         if( isDiagonal( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDiagonal evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-square zero matrix
      {
         blaze::ZeroMatrix<int,blaze::rowMajor> mat( 3UL, 5UL );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 5UL );
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

      // Square zero matrix
      {
         blaze::ZeroMatrix<int,blaze::rowMajor> mat( 3UL, 3UL );

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
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major isDiagonal()";

      // Default zero matrix
      {
         blaze::ZeroMatrix<int,blaze::columnMajor> mat;

         checkRows    ( mat, 0UL );
         checkColumns ( mat, 0UL );
         checkNonZeros( mat, 0UL );

         if( isDiagonal( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDiagonal evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-square zero matrix
      {
         blaze::ZeroMatrix<int,blaze::columnMajor> mat( 5UL, 3UL );

         checkRows    ( mat, 5UL );
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

      // Square zero matrix
      {
         blaze::ZeroMatrix<int,blaze::columnMajor> mat( 3UL, 3UL );

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
void ZeroTest::testIsIdentity()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major isIdentity()";

      // Default zero matrix
      {
         blaze::ZeroMatrix<int,blaze::rowMajor> mat;

         checkRows    ( mat, 0UL );
         checkColumns ( mat, 0UL );
         checkNonZeros( mat, 0UL );

         if( isIdentity( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isIdentity evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-square zero matrix
      {
         blaze::ZeroMatrix<int,blaze::rowMajor> mat( 3UL, 5UL );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 5UL );
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

      // Square zero matrix
      {
         blaze::ZeroMatrix<int,blaze::rowMajor> mat( 3UL, 3UL );

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
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major isIdentity()";

      // Default zero matrix
      {
         blaze::ZeroMatrix<int,blaze::columnMajor> mat;

         checkRows    ( mat, 0UL );
         checkColumns ( mat, 0UL );
         checkNonZeros( mat, 0UL );

         if( isIdentity( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isIdentity evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Non-square zero matrix
      {
         blaze::ZeroMatrix<int,blaze::columnMajor> mat( 5UL, 3UL );

         checkRows    ( mat, 5UL );
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

      // Square zero matrix
      {
         blaze::ZeroMatrix<int,blaze::columnMajor> mat( 3UL, 3UL );

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
void ZeroTest::testMean()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major mean()";

      {
         blaze::ZeroMatrix<int,blaze::rowMajor> mat( 3UL, 3UL );

         const double mean = blaze::mean( mat );

         if( mean != 0.0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Mean computation failed\n"
                << " Details:\n"
                << "   Result: " << mean << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      try {
         blaze::ZeroMatrix<int,blaze::rowMajor> mat( 3UL, 0UL );

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
         blaze::ZeroMatrix<int,blaze::rowMajor> mat( 0UL, 3UL );

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
         blaze::ZeroMatrix<int,blaze::rowMajor> mat( 3UL, 3UL );

         blaze::UniformVector<double,blaze::columnVector> mean;
         mean = blaze::mean<blaze::rowwise>( mat );

         if( mean[0] != 0.0 || mean[1] != 0.0 || mean[2] != 0.0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Mean computation failed\n"
                << " Details:\n"
                << "   Result: " << trans( mean ) << "\n"
                << "   Expected result: ( 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      try {
         blaze::ZeroMatrix<int,blaze::rowMajor> mat( 3UL, 0UL );

         blaze::UniformVector<double,blaze::columnVector> mean;
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
         blaze::ZeroMatrix<int,blaze::rowMajor> mat( 3UL, 3UL );

         blaze::UniformVector<double,blaze::rowVector> mean;
         mean = blaze::mean<blaze::columnwise>( mat );

         if( mean[0] != 0.0 || mean[1] != 0.0 || mean[2] != 0.0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Mean computation failed\n"
                << " Details:\n"
                << "   Result: " << mean << "\n"
                << "   Expected result: ( 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      try {
         blaze::ZeroMatrix<int,blaze::rowMajor> mat( 0UL, 3UL );

         blaze::UniformVector<double,blaze::rowVector> mean;
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
         blaze::ZeroMatrix<int,blaze::columnMajor> mat( 3UL, 3UL );

         const double mean = blaze::mean( mat );

         if( mean != 0.0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Mean computation failed\n"
                << " Details:\n"
                << "   Result: " << mean << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      try {
         blaze::ZeroMatrix<int,blaze::columnMajor> mat( 3UL, 0UL );

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
         blaze::ZeroMatrix<int,blaze::columnMajor> mat( 0UL, 3UL );

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
         blaze::ZeroMatrix<int,blaze::columnMajor> mat( 3UL, 3UL );

         blaze::UniformVector<double,blaze::columnVector> mean;
         mean = blaze::mean<blaze::rowwise>( mat );

         if( mean[0] != 0.0 || mean[1] != 0.0 || mean[2] != 0.0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Mean computation failed\n"
                << " Details:\n"
                << "   Result: " << trans( mean ) << "\n"
                << "   Expected result: ( 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      try {
         blaze::ZeroMatrix<int,blaze::columnMajor> mat( 3UL, 0UL );

         blaze::UniformVector<double,blaze::columnVector> mean;
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
         blaze::ZeroMatrix<int,blaze::columnMajor> mat( 3UL, 3UL );

         blaze::UniformVector<double,blaze::rowVector> mean;
         mean = blaze::mean<blaze::columnwise>( mat );

         if( mean[0] != 0.0 || mean[1] != 0.0 || mean[2] != 0.0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Mean computation failed\n"
                << " Details:\n"
                << "   Result: " << mean << "\n"
                << "   Expected result: ( 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      try {
         blaze::ZeroMatrix<int,blaze::columnMajor> mat( 0UL, 3UL );

         blaze::UniformVector<double,blaze::rowVector> mean;
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
void ZeroTest::testVar()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major var()";

      {
         blaze::ZeroMatrix<int,blaze::rowMajor> mat( 3UL, 3UL );

         const double var = blaze::var( mat );

         if( var != 0.0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Variance computation failed\n"
                << " Details:\n"
                << "   Result: " << var << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      try {
         blaze::ZeroMatrix<int,blaze::rowMajor> mat( 3UL, 0UL );

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
         blaze::ZeroMatrix<int,blaze::rowMajor> mat( 0UL, 3UL );

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
         blaze::ZeroMatrix<int,blaze::rowMajor> mat( 1UL, 1UL );

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
         blaze::ZeroMatrix<int,blaze::rowMajor> mat( 3UL, 3UL );

         blaze::UniformVector<double,blaze::columnVector> var;
         var = blaze::var<blaze::rowwise>( mat );

         if( var[0] != 0.0 || var[1] != 0.0 || var[2] != 0.0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Variance computation failed\n"
                << " Details:\n"
                << "   Result: " << trans( var ) << "\n"
                << "   Expected result: ( 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      try {
         blaze::ZeroMatrix<int,blaze::rowMajor> mat( 3UL, 0UL );

         blaze::UniformVector<double,blaze::columnVector> var;
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
         blaze::ZeroMatrix<int,blaze::rowMajor> mat( 3UL, 1UL );

         blaze::UniformVector<double,blaze::columnVector> var;
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
         blaze::ZeroMatrix<int,blaze::rowMajor> mat( 3UL, 3UL );

         blaze::UniformVector<double,blaze::rowVector> var;
         var = blaze::var<blaze::columnwise>( mat );

         if( var[0] != 0.0 || var[1] != 0.0 || var[2] != 0.0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Variance computation failed\n"
                << " Details:\n"
                << "   Result: " << var << "\n"
                << "   Expected result: ( 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      try {
         blaze::ZeroMatrix<int,blaze::rowMajor> mat( 0UL, 3UL );

         blaze::UniformVector<double,blaze::rowVector> var;
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
         blaze::ZeroMatrix<int,blaze::rowMajor> mat( 1UL, 3UL );

         blaze::UniformVector<double,blaze::rowVector> var;
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
         blaze::ZeroMatrix<int,blaze::columnMajor> mat( 3UL, 3UL );

         const double var = blaze::var( mat );

         if( var != 0.0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Variance computation failed\n"
                << " Details:\n"
                << "   Result: " << var << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      try {
         blaze::ZeroMatrix<int,blaze::columnMajor> mat( 3UL, 0UL );

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
         blaze::ZeroMatrix<int,blaze::columnMajor> mat( 0UL, 3UL );

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
         blaze::ZeroMatrix<int,blaze::columnMajor> mat( 1UL, 1UL );

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
         blaze::ZeroMatrix<int,blaze::columnMajor> mat( 3UL, 3UL );

         blaze::UniformVector<double,blaze::columnVector> var;
         var = blaze::var<blaze::rowwise>( mat );

         if( var[0] != 0.0 || var[1] != 0.0 || var[2] != 0.0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Variance computation failed\n"
                << " Details:\n"
                << "   Result: " << trans( var ) << "\n"
                << "   Expected result: ( 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      try {
         blaze::ZeroMatrix<int,blaze::columnMajor> mat( 3UL, 0UL );

         blaze::UniformVector<double,blaze::columnVector> var;
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
         blaze::ZeroMatrix<int,blaze::columnMajor> mat( 3UL, 1UL );

         blaze::UniformVector<double,blaze::columnVector> var;
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
         blaze::ZeroMatrix<int,blaze::columnMajor> mat( 3UL, 3UL );

         blaze::UniformVector<double,blaze::rowVector> var;
         var = blaze::var<blaze::columnwise>( mat );

         if( var[0] != 0.0 || var[1] != 0.0 || var[2] != 0.0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Variance computation failed\n"
                << " Details:\n"
                << "   Result: " << var << "\n"
                << "   Expected result: ( 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      try {
         blaze::ZeroMatrix<int,blaze::columnMajor> mat( 0UL, 3UL );

         blaze::UniformVector<double,blaze::rowVector> var;
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
         blaze::ZeroMatrix<int,blaze::columnMajor> mat( 1UL, 3UL );

         blaze::UniformVector<double,blaze::rowVector> var;
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
void ZeroTest::testStdDev()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major stddev()";

      {
         blaze::ZeroMatrix<int,blaze::rowMajor> mat( 3UL, 3UL );

         const double stddev = blaze::stddev( mat );

         if( stddev != 0.0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Standard deviation computation failed\n"
                << " Details:\n"
                << "   Result: " << stddev << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      try {
         blaze::ZeroMatrix<int,blaze::rowMajor> mat( 3UL, 0UL );

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
         blaze::ZeroMatrix<int,blaze::rowMajor> mat( 0UL, 3UL );

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
         blaze::ZeroMatrix<int,blaze::rowMajor> mat( 1UL, 1UL );

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
         blaze::ZeroMatrix<int,blaze::rowMajor> mat( 3UL, 3UL );

         blaze::UniformVector<double,blaze::columnVector> stddev;
         stddev = blaze::stddev<blaze::rowwise>( mat );

         if( stddev[0] != 0.0 || stddev[1] != 0.0 || stddev[2] != 0.0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Standard deviation computation failed\n"
                << " Details:\n"
                << "   Result: " << trans( stddev ) << "\n"
                << "   Expected result: ( 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      try {
         blaze::ZeroMatrix<int,blaze::rowMajor> mat( 3UL, 0UL );

         blaze::UniformVector<double,blaze::columnVector> stddev;
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
         blaze::ZeroMatrix<int,blaze::rowMajor> mat( 3UL, 1UL );

         blaze::UniformVector<double,blaze::columnVector> stddev;
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
         blaze::ZeroMatrix<int,blaze::rowMajor> mat( 3UL, 3UL );

         blaze::UniformVector<double,blaze::rowVector> stddev;
         stddev = blaze::stddev<blaze::columnwise>( mat );

         if( stddev[0] != 0.0 || stddev[1] != 0.0 || stddev[2] != 0.0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Standard deviation computation failed\n"
                << " Details:\n"
                << "   Result: " << stddev << "\n"
                << "   Expected result: ( 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      try {
         blaze::ZeroMatrix<int,blaze::rowMajor> mat( 0UL, 3UL );

         blaze::UniformVector<double,blaze::rowVector> stddev;
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
         blaze::ZeroMatrix<int,blaze::rowMajor> mat( 1UL, 3UL );

         blaze::UniformVector<double,blaze::rowVector> stddev;
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
         blaze::ZeroMatrix<int,blaze::columnMajor> mat( 3UL, 3UL );

         const double stddev = blaze::stddev( mat );

         if( stddev != 0.0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Standard deviation computation failed\n"
                << " Details:\n"
                << "   Result: " << stddev << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      try {
         blaze::ZeroMatrix<int,blaze::columnMajor> mat( 3UL, 0UL );

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
         blaze::ZeroMatrix<int,blaze::columnMajor> mat( 0UL, 3UL );

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
         blaze::ZeroMatrix<int,blaze::columnMajor> mat( 1UL, 1UL );

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
         blaze::ZeroMatrix<int,blaze::columnMajor> mat( 3UL, 3UL );

         blaze::UniformVector<double,blaze::columnVector> stddev;
         stddev = blaze::stddev<blaze::rowwise>( mat );

         if( stddev[0] != 0.0 || stddev[1] != 0.0 || stddev[2] != 0.0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Standard deviation computation failed\n"
                << " Details:\n"
                << "   Result: " << trans( stddev ) << "\n"
                << "   Expected result: ( 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      try {
         blaze::ZeroMatrix<int,blaze::columnMajor> mat( 3UL, 0UL );

         blaze::UniformVector<double,blaze::columnVector> stddev;
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
         blaze::ZeroMatrix<int,blaze::columnMajor> mat( 3UL, 1UL );

         blaze::UniformVector<double,blaze::columnVector> stddev;
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
         blaze::ZeroMatrix<int,blaze::columnMajor> mat( 3UL, 3UL );

         blaze::UniformVector<double,blaze::rowVector> stddev;
         stddev = blaze::stddev<blaze::columnwise>( mat );

         if( stddev[0] != 0.0 || stddev[1] != 0.0 || stddev[2] != 0.0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Standard deviation computation failed\n"
                << " Details:\n"
                << "   Result: " << stddev << "\n"
                << "   Expected result: ( 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      try {
         blaze::ZeroMatrix<int,blaze::columnMajor> mat( 0UL, 3UL );

         blaze::UniformVector<double,blaze::rowVector> stddev;
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
         blaze::ZeroMatrix<int,blaze::columnMajor> mat( 1UL, 3UL );

         blaze::UniformVector<double,blaze::rowVector> stddev;
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
   std::cout << "   Running zero SparseMatrix operation test..." << std::endl;

   try
   {
      RUN_SPARSEMATRIX_ZERO_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during zero SparseMatrix operation test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
