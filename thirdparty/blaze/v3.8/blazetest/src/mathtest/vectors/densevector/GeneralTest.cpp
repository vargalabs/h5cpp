//=================================================================================================
/*!
//  \file src/mathtest/vectors/densevector/GeneralTest.cpp
//  \brief Source file for the general DenseVector operation test
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
#include <blaze/math/dense/DenseVector.h>
#include <blaze/math/DynamicVector.h>
#include <blaze/math/StaticVector.h>
#include <blazetest/mathtest/IsEqual.h>
#include <blazetest/mathtest/vectors/densevector/GeneralTest.h>

#ifdef BLAZE_USE_HPX_THREADS
#  include <hpx/hpx_main.hpp>
#endif


namespace blazetest {

namespace mathtest {

namespace vectors {

namespace densevector {

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
   testIsUniform();
   testIsZero();
   testNormalize();
   testMinimum();
   testMaximum();
   testArgmin();
   testArgmax();
   testL1Norm();
   testL2Norm();
   testL3Norm();
   testL4Norm();
   testLpNorm();
   testLinfNorm();
   testLength();
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
   testLinspace();
   testLogspace();
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
/*!\brief Test of the \c isnan() function for dense vectors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isnan() function for dense vectors. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testIsNan()
{
   test_ = "isnan() function";

   // isnan with 0-dimensional vector (non-NAN)
   {
      blaze::DynamicVector<float,blaze::rowVector> vec;

      if( blaze::isnan( vec ) != false ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isnan evaluation\n"
             << " Details:\n"
             << "   Vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // isnan with 6-dimensional vector (non-NAN)
   {
      blaze::DynamicVector<float,blaze::rowVector> vec{ 1.2F, -2.1F, 3.2F, -4.1F, 5.2F, -6.1F };

      if( blaze::isnan( vec ) != false ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isnan evaluation\n"
             << " Details:\n"
             << "   Vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // isnan with 6-dimensional vector (NAN)
   {
      blaze::DynamicVector<float,blaze::rowVector> vec{ 1.2F, -2.1F, 3.2F, -4.1F, NAN, -6.1F };

      if( blaze::isnan( vec ) != true ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isnan evaluation\n"
             << " Details:\n"
             << "   Vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c isinf() function for dense vectors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isinf() function for dense vectors. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testIsInf()
{
   test_ = "isinf() function";

   // isinf with 0-dimensional vector (non-inf)
   {
      blaze::DynamicVector<float,blaze::rowVector> vec;

      if( blaze::isinf( vec ) != false ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isinf evaluation\n"
             << " Details:\n"
             << "   Vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // isinf with 6-dimensional vector (non-inf)
   {
      blaze::DynamicVector<float,blaze::rowVector> vec{ 1.2F, -2.1F, 3.2F, -4.1F, 5.2F, -6.1F };

      if( blaze::isinf( vec ) != false ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isinf evaluation\n"
             << " Details:\n"
             << "   Vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // isinf with 6-dimensional vector (HUGE_VALF)
   {
      blaze::DynamicVector<float,blaze::rowVector> vec{ 1.2F, -2.1F, 3.2F, -4.1F, HUGE_VALF, -6.1F };

      if( blaze::isinf( vec ) != true ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isinf evaluation\n"
             << " Details:\n"
             << "   Vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c isfinite() function for dense vectors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isfinite() function for dense vectors. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testIsFinite()
{
   test_ = "isfinite() function";

   // isfinite with 0-dimensional vector (finite)
   {
      blaze::DynamicVector<float,blaze::rowVector> vec;

      if( blaze::isfinite( vec ) != true ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isfinite evaluation\n"
             << " Details:\n"
             << "   Vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // isfinite with 6-dimensional vector (finite)
   {
      blaze::DynamicVector<float,blaze::rowVector> vec{ 1.2F, -2.1F, 3.2F, -4.1F, 5.2F, -6.1F };

      if( blaze::isfinite( vec ) != true ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isfinite evaluation\n"
             << " Details:\n"
             << "   Vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // isfinite with 6-dimensional vector (NAN)
   {
      blaze::DynamicVector<float,blaze::rowVector> vec{ 1.2F, -2.1F, 3.2F, -4.1F, NAN, -6.1F };

      if( blaze::isfinite( vec ) != false ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isfinite evaluation\n"
             << " Details:\n"
             << "   Vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // isfinite with 6-dimensional vector (HUGE_VALF)
   {
      blaze::DynamicVector<float,blaze::rowVector> vec{ 1.2F, -2.1F, 3.2F, -4.1F, HUGE_VALF, -6.1F };

      if( blaze::isfinite( vec ) != false ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isfinite evaluation\n"
             << " Details:\n"
             << "   Vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c isUniform() function for dense vectors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isUniform() function for dense vectors. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testIsUniform()
{
   test_ = "isUniform() function";

   // Uniform vector (0-dimensional)
   {
      blaze::DynamicVector<int,blaze::rowVector> vec;

      if( blaze::isUniform( vec ) != true ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isUniform evaluation\n"
             << " Details:\n"
             << "   Vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Uniform vector (1-dimensional)
   {
      blaze::DynamicVector<int,blaze::rowVector> vec{ 5 };

      if( blaze::isUniform( vec ) != true ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isUniform evaluation\n"
             << " Details:\n"
             << "   Vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Uniform vector (5-dimensional)
   {
      blaze::DynamicVector<int,blaze::rowVector> vec{ 5, 5, 5, 5, 5 };

      if( blaze::isUniform( vec ) != true ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isUniform evaluation\n"
             << " Details:\n"
             << "   Vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Non-uniform vector (5-dimensional)
   {
      blaze::DynamicVector<int,blaze::rowVector> vec{ 5, 5, 5, 5, 3 };

      if( blaze::isUniform( vec ) != false ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isUniform evaluation\n"
             << " Details:\n"
             << "   Vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c isZero() function for dense vectors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isZero() function for dense vectors. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testIsZero()
{
   test_ = "isZero() function";

   // Zero vector (0-dimensional)
   {
      blaze::DynamicVector<int,blaze::rowVector> vec;

      if( blaze::isZero( vec ) != true ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isZero evaluation\n"
             << " Details:\n"
             << "   Vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Zero vector (1-dimensional)
   {
      blaze::DynamicVector<int,blaze::rowVector> vec{ 0 };

      if( blaze::isZero( vec ) != true ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isZero evaluation\n"
             << " Details:\n"
             << "   Vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Zero vector (5-dimensional)
   {
      blaze::DynamicVector<int,blaze::rowVector> vec{ 0, 0, 0, 0, 0 };

      if( blaze::isZero( vec ) != true ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isZero evaluation\n"
             << " Details:\n"
             << "   Vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Non-zero vector (5-dimensional)
   {
      blaze::DynamicVector<int,blaze::rowVector> vec{ 0, 0, 0, 0, 3 };

      if( blaze::isZero( vec ) != false ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isZero evaluation\n"
             << " Details:\n"
             << "   Vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c normalize() function for dense vectors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c normalize() function for dense vectors. In case
// an error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testNormalize()
{
   test_ = "normalize() function";

   // Initialization check
   blaze::DynamicVector<double,blaze::rowVector> vec( 4UL );
   vec[0] = 1.0;
   vec[1] = 2.0;
   vec[2] = 3.0;
   vec[3] = 4.0;

   checkSize    ( vec, 4UL );
   checkCapacity( vec, 4UL );
   checkNonZeros( vec, 4UL );

   if( vec[0] != 1.0 || vec[1] != 2.0 || vec[2] != 3.0 || vec[3] != 4.0 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Initialization failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 1 2 3 4 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Acquiring normalized vector
   const blaze::DynamicVector<double,blaze::rowVector> normalized( normalize( vec ) );

   if( !blaze::equal( length( normalized ), 1.0 ) ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Normalization failed\n"
          << " Details:\n"
          << "   Result: " << length( normalized ) << "\n"
          << "   Expected result: 1\n";
      throw std::runtime_error( oss.str() );
   }

   // Normalizing the vector
   vec = normalize( vec );

   if( !blaze::equal( length( vec ), 1.0 ) ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Self-normalization failed\n"
          << " Details:\n"
          << "   Result: " << length( vec ) << "\n"
          << "   Expected result: 1\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c min() function for dense vectors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c min() function for dense vectors. In case an error
// is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testMinimum()
{
   test_ = "min() function";

   {
      // Initialization check
      blaze::DynamicVector<int,blaze::rowVector> vec{ 1, -2, 3, -4 };

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );
      checkNonZeros( vec, 4UL );

      if( vec[0] != 1 || vec[1] != -2 || vec[2] != 3 || vec[3] != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 -2 3 -4 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Testing the min function
      const int minimum = min( vec );

      if( minimum != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: First computation failed\n"
             << " Details:\n"
             << "   Result: " << minimum << "\n"
             << "   Expected result: -4\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      // Initialization check
      blaze::DynamicVector<int,blaze::rowVector> vec{ -1, 2, 3, 4 };
      vec[0] = -1;
      vec[1] =  2;
      vec[2] =  3;
      vec[3] =  4;

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );
      checkNonZeros( vec, 4UL );

      if( vec[0] != -1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( -1 2 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Testing the min function
      const int minimum = min( vec );

      if( minimum != -1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Second computation failed\n"
             << " Details:\n"
             << "   Result: " << minimum << "\n"
             << "   Expected result: -1\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c max() function for dense vectors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c max() function for dense vectors. In case an error
// is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testMaximum()
{
   test_ = "max() function";

   {
      // Initialization check
      blaze::DynamicVector<int,blaze::rowVector> vec{ 1, -2, -3, -4 };

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );
      checkNonZeros( vec, 4UL );

      if( vec[0] != 1 || vec[1] != -2 || vec[2] != -3 || vec[3] != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 -2 -3 -4 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Testing the max function
      const int maximum = max( vec );

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

   {
      // Initialization check
      blaze::DynamicVector<int,blaze::rowVector> vec{ -1, 2, 3, 4 };

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );
      checkNonZeros( vec, 4UL );

      if( vec[0] != -1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( -1 2 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Testing the max function
      const int maximum = max( vec );

      if( maximum != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Second computation failed\n"
             << " Details:\n"
             << "   Result: " << maximum << "\n"
             << "   Expected result: 4\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c argmin() function for dense vectors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c argmin() function for dense vectors. In case an error
// is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testArgmin()
{
   test_ = "argmin() function";

   {
      blaze::DynamicVector<int,blaze::rowVector> vec;

      const size_t minimum = argmin( vec );

      checkSize    ( vec, 0UL );
      checkCapacity( vec, 0UL );
      checkNonZeros( vec, 0UL );

      if( minimum != 0UL ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Argmin evaluation failed\n"
             << " Details:\n"
             << "   Result: " << minimum << "\n"
             << "   Expected result: 0\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      blaze::DynamicVector<int,blaze::rowVector> vec{ 99 };

      const size_t minimum = argmin( vec );

      checkSize    ( vec, 1UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      if( minimum != 0UL ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Argmin evaluation failed\n"
             << " Details:\n"
             << "   Result: " << minimum << "\n"
             << "   Expected result: 0\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      blaze::DynamicVector<int,blaze::rowVector> vec{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };

      const size_t minimum = argmin( vec );

      checkSize    ( vec, 9UL );
      checkCapacity( vec, 9UL );
      checkNonZeros( vec, 9UL );

      if( minimum != 0UL ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Argmin evaluation failed\n"
             << " Details:\n"
             << "   Result: " << minimum << "\n"
             << "   Expected result: 0\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      blaze::DynamicVector<int,blaze::rowVector> vec{ 9, 8, 7, 6, 5, 4, 3, 2, 1 };

      const size_t minimum = argmin( vec );

      checkSize    ( vec, 9UL );
      checkCapacity( vec, 9UL );
      checkNonZeros( vec, 9UL );

      if( minimum != 8UL ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Argmin evaluation failed\n"
             << " Details:\n"
             << "   Result: " << minimum << "\n"
             << "   Expected result: 8\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      blaze::DynamicVector<int,blaze::rowVector> vec{ 2, 3, 4, 5, 1, 6, 7, 8, 9 };

      const size_t minimum = argmin( vec );

      checkSize    ( vec, 9UL );
      checkCapacity( vec, 9UL );
      checkNonZeros( vec, 9UL );

      if( minimum != 4UL ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Argmin evaluation failed\n"
             << " Details:\n"
             << "   Result: " << minimum << "\n"
             << "   Expected result: 4\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c argmax() function for dense vectors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c argmax() function for dense vectors. In case an error
// is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testArgmax()
{
   test_ = "argmax() function";

   {
      blaze::DynamicVector<int,blaze::rowVector> vec;

      const size_t maximum = argmax( vec );

      checkSize    ( vec, 0UL );
      checkCapacity( vec, 0UL );
      checkNonZeros( vec, 0UL );

      if( maximum != 0UL ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Argmax evaluation failed\n"
             << " Details:\n"
             << "   Result: " << maximum << "\n"
             << "   Expected result: 0\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      blaze::DynamicVector<int,blaze::rowVector> vec{ 99 };

      const size_t maximum = argmax( vec );

      checkSize    ( vec, 1UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      if( maximum != 0UL ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Argmax evaluation failed\n"
             << " Details:\n"
             << "   Result: " << maximum << "\n"
             << "   Expected result: 0\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      blaze::DynamicVector<int,blaze::rowVector> vec{ 9, 8, 7, 6, 5, 4, 3, 2, 1 };

      const size_t maximum = argmax( vec );

      checkSize    ( vec, 9UL );
      checkCapacity( vec, 9UL );
      checkNonZeros( vec, 9UL );

      if( maximum != 0UL ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Argmax evaluation failed\n"
             << " Details:\n"
             << "   Result: " << maximum << "\n"
             << "   Expected result: 0\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      blaze::DynamicVector<int,blaze::rowVector> vec{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };

      const size_t maximum = argmax( vec );

      checkSize    ( vec, 9UL );
      checkCapacity( vec, 9UL );
      checkNonZeros( vec, 9UL );

      if( maximum != 8UL ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Argmax evaluation failed\n"
             << " Details:\n"
             << "   Result: " << maximum << "\n"
             << "   Expected result: 8\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      blaze::DynamicVector<int,blaze::rowVector> vec{ 8, 7, 6, 5, 9, 4, 3, 2, 1 };

      const size_t maximum = argmax( vec );

      checkSize    ( vec, 9UL );
      checkCapacity( vec, 9UL );
      checkNonZeros( vec, 9UL );

      if( maximum != 4UL ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Argmax evaluation failed\n"
             << " Details:\n"
             << "   Result: " << maximum << "\n"
             << "   Expected result: 4\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c l1Norm() function for dense vectors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c l1Norm() function for dense vectors. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testL1Norm()
{
   test_ = "l1Norm() function";

   {
      blaze::DynamicVector<int,blaze::rowVector> vec;

      const int norm = blaze::l1Norm( vec );

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
      blaze::DynamicVector<int,blaze::rowVector> vec( 7UL, 0 );

      const int norm = blaze::l1Norm( vec );

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
      blaze::DynamicVector<int,blaze::rowVector> vec{ 0, -1, 2, -2, 0, 0, -1, 0, 1, 0 };

      const int norm = blaze::l1Norm( vec );

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
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c l2Norm() function for dense vectors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c l2Norm() function for dense vectors. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testL2Norm()
{
   test_ = "l2Norm() function";

   {
      blaze::DynamicVector<int,blaze::rowVector> vec;

      const double norm = blaze::l2Norm( vec );

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
      blaze::DynamicVector<int,blaze::rowVector> vec( 7UL, 0 );

      const double norm = blaze::l2Norm( vec );

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
      blaze::DynamicVector<int,blaze::rowVector> vec{ 0, -1, 2, -2, 2, 1, -1, 0, 1, 0 };

      const double norm = blaze::l2Norm( vec );

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
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c l3Norm() function for dense vectors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c l3Norm() function for dense vectors. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testL3Norm()
{
   test_ = "l3Norm() function";

   {
      blaze::DynamicVector<int,blaze::rowVector> vec;

      const double norm = blaze::l3Norm( vec );

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
      blaze::DynamicVector<int,blaze::rowVector> vec( 7UL, 0 );

      const double norm = l3Norm( vec );

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
      blaze::DynamicVector<int,blaze::rowVector> vec{ 0, -1, 2, -2, 2, 0, -1, 0, 1, 0 };

      const double norm = blaze::l3Norm( vec );

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
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c l4Norm() function for dense vectors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c l4Norm() function for dense vectors. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testL4Norm()
{
   test_ = "l4Norm() function";

   {
      blaze::DynamicVector<int,blaze::rowVector> vec;

      const double norm = blaze::l4Norm( vec );

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
      blaze::DynamicVector<int,blaze::rowVector> vec( 7UL, 0 );

      const double norm = blaze::l4Norm( vec );

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
      blaze::DynamicVector<int,blaze::rowVector> vec{ 0, 2, 0, -2, 2, -1, 0, -2, 0, 2 };

      const double norm = blaze::l4Norm( vec );

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
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c lpNorm() function for dense vectors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c lpNorm() function for dense vectors. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testLpNorm()
{
   test_ = "lpNorm() function";

   {
      blaze::DynamicVector<int,blaze::rowVector> vec;

      const double norm1 = blaze::lpNorm( vec, 2 );
      const double norm2 = blaze::lpNorm<2UL>( vec );

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
      blaze::DynamicVector<int,blaze::rowVector> vec( 7UL, 0 );

      const double norm1 = blaze::lpNorm( vec, 2 );
      const double norm2 = blaze::lpNorm<2UL>( vec );

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
      blaze::DynamicVector<int,blaze::rowVector> vec( 10UL );
      randomize( vec, -5, 5 );

      const int norm1( blaze::lpNorm( vec, 1 ) );
      const int norm2( blaze::lpNorm<1UL>( vec ) );
      const int norm3( blaze::l1Norm( vec ) );

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
      blaze::DynamicVector<int,blaze::rowVector> vec( 10UL );
      randomize( vec, -5, 5 );

      const double norm1( blaze::lpNorm( vec, 2 ) );
      const double norm2( blaze::lpNorm<2UL>( vec ) );
      const double norm3( blaze::l2Norm( vec ) );

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
      blaze::DynamicVector<int,blaze::rowVector> vec( 10UL );
      randomize( vec, -5, 5 );

      const double norm1( blaze::lpNorm( vec, 3 ) );
      const double norm2( blaze::lpNorm<3UL>( vec ) );
      const double norm3( blaze::l3Norm( vec ) );

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
      blaze::DynamicVector<int,blaze::rowVector> vec( 10UL );
      randomize( vec, -5, 5 );

      const double norm1( blaze::lpNorm( vec, 4 ) );
      const double norm2( blaze::lpNorm<4UL>( vec ) );
      const double norm3( blaze::l4Norm( vec ) );

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
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c linfNorm() function for dense vectors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c linfNorm() function for dense vectors. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testLinfNorm()
{
   test_ = "linfNorm() function";

   {
      blaze::DynamicVector<int,blaze::rowVector> vec;

      const double norm = blaze::linfNorm( vec );

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
      blaze::DynamicVector<int,blaze::rowVector> vec( 7UL, 0 );

      const double norm = blaze::linfNorm( vec );

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
      blaze::DynamicVector<int,blaze::rowVector> vec( 10UL );
      randomize( vec, -5, 5 );

      const int norm1( blaze::linfNorm( vec ) );
      const int norm2( blaze::max( blaze::abs( vec ) ) );

      if( !isEqual( norm1, norm2 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Lp norm computation failed\n"
             << " Details:\n"
             << "   linfNorm(): " << norm1 << "\n"
             << "   Expected result: " << norm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c length() and \c sqrLength() functions for dense vectors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c length() and \c sqrLength() functions for dense
// vectors. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testLength()
{
   test_ = "length() and sqrLength() functions";

   {
      blaze::DynamicVector<double,blaze::rowVector> vec;

      // Computing the vector length
      const double len( length( vec ) );

      if( !blaze::equal( len, 0.0 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Length computation failed\n"
             << " Details:\n"
             << "   Result: " << len << "\n"
             << "   Expected result: 0\n";
         throw std::runtime_error( oss.str() );
      }

      // Computing the vector square length
      const double sqrlen( sqrLength( vec ) );

      if( !blaze::equal( sqrlen, 0.0 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Square length computation failed\n"
             << " Details:\n"
             << "   Result: " << sqrlen << "\n"
             << "   Expected result: 0\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      blaze::DynamicVector<double,blaze::rowVector> vec( 2UL );
      vec[0] = 0.0;
      vec[1] = 0.0;

      // Computing the vector length
      const double len( length( vec ) );

      if( !blaze::equal( len, 0.0 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Length computation failed\n"
             << " Details:\n"
             << "   Result: " << len << "\n"
             << "   Expected result: 0\n";
         throw std::runtime_error( oss.str() );
      }

      // Computing the vector square length
      const double sqrlen( sqrLength( vec ) );

      if( !blaze::equal( sqrlen, 0.0 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Square length computation failed\n"
             << " Details:\n"
             << "   Result: " << sqrlen << "\n"
             << "   Expected result: 0\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      blaze::DynamicVector<double,blaze::rowVector> vec( 2UL );
      vec[0] = 3.0;
      vec[1] = 4.0;

      // Computing the vector length
      const double len( length( vec ) );

      if( !blaze::equal( len, 5.0 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Length computation failed\n"
             << " Details:\n"
             << "   Result: " << len << "\n"
             << "   Expected result: 5\n";
         throw std::runtime_error( oss.str() );
      }

      // Computing the vector square length
      const double sqrlen( sqrLength( vec ) );

      if( !blaze::equal( sqrlen, 25.0 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Square length computation failed\n"
             << " Details:\n"
             << "   Result: " << sqrlen << "\n"
             << "   Expected result: 25\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c mean() function for dense vectors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c mean() function for dense vectors. In case an error
// is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testMean()
{
   test_ = "mean() function";

   {
      blaze::DynamicVector<int,blaze::rowVector> vec( 5UL, 0 );

      const double mean = blaze::mean( vec );

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
      blaze::DynamicVector<int,blaze::rowVector> vec{ 1, 4, 3, 6, 7 };

      const double mean = blaze::mean( vec );

      if( !isEqual( mean, 4.2 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Mean computation failed\n"
             << " Details:\n"
             << "   Result: " << mean << "\n"
             << "   Expected result: 4.2\n";
         throw std::runtime_error( oss.str() );
      }
   }

   try {
      blaze::DynamicVector<int,blaze::rowVector> vec;

      const double mean = blaze::mean( vec );

      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Mean computation of empty vector succeeded\n"
          << " Details:\n"
          << "   Result:\n" << mean << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::invalid_argument& ) {}
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c var() function for dense vectors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c var() function for dense vectors. In case an error
// is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testVar()
{
   test_ = "var() function";

   {
      blaze::DynamicVector<int,blaze::rowVector> vec( 5UL, 0 );

      const double var = blaze::var( vec );

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
      blaze::DynamicVector<int,blaze::rowVector> vec{ 1, 4, 3, 6, 7 };

      const double var = blaze::var( vec );

      if( !isEqual( var, 5.7 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Variance computation failed\n"
             << " Details:\n"
             << "   Result: " << var << "\n"
             << "   Expected result: 5.7\n";
         throw std::runtime_error( oss.str() );
      }
   }

   try {
      blaze::DynamicVector<int,blaze::rowVector> vec;

      const double var = blaze::var( vec );

      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Variance computation of empty vector succeeded\n"
          << " Details:\n"
          << "   Result:\n" << var << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::invalid_argument& ) {}

   try {
      blaze::DynamicVector<int,blaze::rowVector> vec( 1UL );

      const double var = blaze::var( vec );

      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Variance computation of 1D vector succeeded\n"
          << " Details:\n"
          << "   Result:\n" << var << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::invalid_argument& ) {}
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c stddev() function for dense vectors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c stddev() function for dense vectors. In case an error
// is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testStdDev()
{
   test_ = "stddev() function";

   {
      blaze::DynamicVector<int,blaze::rowVector> vec( 5UL, 0 );

      const double stddev = blaze::stddev( vec );

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
      blaze::DynamicVector<int,blaze::rowVector> vec{ 1, 4, 3, 6, 7 };

      const double stddev = blaze::stddev( vec );

      if( !isEqual( stddev, std::sqrt(5.7) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Standard deviation computation failed\n"
             << " Details:\n"
             << "   Result: " << stddev << "\n"
             << "   Expected result: sqrt(5.7)\n";
         throw std::runtime_error( oss.str() );
      }
   }

   try {
      blaze::DynamicVector<int,blaze::rowVector> vec;

      const double stddev = blaze::stddev( vec );

      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Standard deviation computation of empty vector succeeded\n"
          << " Details:\n"
          << "   Result:\n" << stddev << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::invalid_argument& ) {}

   try {
      blaze::DynamicVector<int,blaze::rowVector> vec( 1UL );

      const double stddev = blaze::stddev( vec );

      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Standard deviation computation of 1D vector succeeded\n"
          << " Details:\n"
          << "   Result:\n" << stddev << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::invalid_argument& ) {}
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c softmax() function for dense vectors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c softmax() function for dense vectors. In case an error
// is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testSoftmax()
{
   test_ = "softmax() function";

   blaze::DynamicVector<double,blaze::rowVector> a( 4UL );
   randomize( a, -5.0, 5.0 );

   const auto b = softmax( a );

   if( b[0] <= 0.0 || b[0] > 1.0 ||
       b[1] <= 0.0 || b[1] > 1.0 ||
       b[2] <= 0.0 || b[2] > 1.0 ||
       b[3] <= 0.0 || b[3] > 1.0 ||
       !isEqual( sum( b ), 1.0 ) ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Softmax computation failed\n"
          << " Details:\n"
          << "   Result: " << sum( b ) << "\n"
          << "   Expected result: 1\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the left-shift operator for dense vectors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the left-shift operator for dense vectors. In case an error
// is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testLeftShift()
{
   test_ = "Left-shift operator";


   //=====================================================================================
   // Vector/scalar left-shift tests
   //=====================================================================================

   // Vector/scalar left-shift of an empty vector
   {
      blaze::DynamicVector<unsigned int> a;

      blaze::DynamicVector<unsigned int> b( a << 2U );

      checkSize    ( b, 0UL );
      checkCapacity( b, 0UL );
      checkNonZeros( b, 0UL );
   }

   // Vector/scalar left-shift of a general vector
   {
      blaze::DynamicVector<unsigned int> a{ 1U, 2U, 4U, 8U, 16U, 32U, 64U, 128U, 256U };

      blaze::DynamicVector<unsigned int> b( a << 2U );

      checkSize    ( b, 9UL );
      checkCapacity( b, 9UL );
      checkNonZeros( b, 9UL );

      if( b[0] !=   4U || b[1] !=   8U || b[2] !=  16U || b[3] !=   32U || b[4] != 64U ||
          b[5] != 128U || b[6] != 256U || b[7] != 512U || b[8] != 1024U ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Vector/scalar left-shift operation failed\n"
             << " Details:\n"
             << "   Result:\n" << b << "\n"
             << "   Expected result:\n( 4 8 16 32 64 128 256 512 1024 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Vector/scalar left-shift assignment
   {
      blaze::DynamicVector<unsigned int> a{ 1U, 2U, 4U, 8U, 16U, 32U, 64U, 128U, 256U };

      a <<= 2U;

      checkSize    ( a, 9UL );
      checkCapacity( a, 9UL );
      checkNonZeros( a, 9UL );

      if( a[0] !=   4U || a[1] !=   8U || a[2] !=  16U || a[3] !=   32U || a[4] != 64U ||
          a[5] != 128U || a[6] != 256U || a[7] != 512U || a[8] != 1024U ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Vector/scalar left-shift assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << a << "\n"
             << "   Expected result:\n( 4 8 16 32 64 128 256 512 1024 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Vector/vector left-shift tests
   //=====================================================================================

   // Vector/vector left-shift of an empty vector
   {
      blaze::DynamicVector<unsigned int> a;
      blaze::DynamicVector<unsigned int> b;

      blaze::DynamicVector<unsigned int> c( a << b );

      checkSize    ( b, 0UL );
      checkCapacity( b, 0UL );
      checkNonZeros( b, 0UL );
   }

   // Vector/vector left-shift of a general vector
   {
      blaze::DynamicVector<unsigned int> a{ 1U, 2U, 4U, 8U, 16U, 32U, 64U, 128U, 256U };
      blaze::DynamicVector<unsigned int> b{ 1U, 2U, 1U, 2U, 1U, 2U, 1U, 2U, 1U };

      blaze::DynamicVector<unsigned int> c( a << b );

      checkSize    ( c, 9UL );
      checkCapacity( c, 9UL );
      checkNonZeros( c, 9UL );

      if( c[0] !=   2U || c[1] !=   8U || c[2] !=   8U || c[3] !=  32U || c[4] != 32U ||
          c[5] != 128U || c[6] != 128U || c[7] != 512U || c[8] != 512U ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Vector/vector left-shift operation failed\n"
             << " Details:\n"
             << "   Result:\n" << c << "\n"
             << "   Expected result:\n( 2 8 8 32 32 128 128 512 512 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Vector/vector left-shift assignment
   {
      blaze::DynamicVector<unsigned int> a{ 1U, 2U, 4U, 8U, 16U, 32U, 64U, 128U, 256U };
      blaze::DynamicVector<unsigned int> b{ 1U, 2U, 1U, 2U, 1U, 2U, 1U, 2U, 1U };

      a <<= b;

      checkSize    ( a, 9UL );
      checkCapacity( a, 9UL );
      checkNonZeros( a, 9UL );

      if( a[0] !=   2U || a[1] !=   8U || a[2] !=   8U || a[3] !=  32U || a[4] != 32U ||
          a[5] != 128U || a[6] != 128U || a[7] != 512U || a[8] != 512U ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Vector/vector left-shift assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << a << "\n"
             << "   Expected result:\n( 2 8 8 32 32 128 128 512 512 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the right-shift operator for dense vectors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the right-shift operator for dense vectors. In case an error
// is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testRightShift()
{
   //=====================================================================================
   // Vector/scalar right-shift tests
   //=====================================================================================

   {
      test_ = "Vector/scalar right-shift operator";

      // Vector/scalar right-shift of an empty vector
      {
         blaze::DynamicVector<unsigned int> a;

         blaze::DynamicVector<unsigned int> b( a >> 2U );

         checkSize    ( b, 0UL );
         checkCapacity( b, 0UL );
         checkNonZeros( b, 0UL );
      }

      // Vector/scalar right-shift of a general vector
      {
         blaze::DynamicVector<unsigned int> a{ 4U, 8U, 16U, 32U, 64U, 128U, 256U, 512U, 1024U };

         blaze::DynamicVector<unsigned int> b( a >> 2U );

         checkSize    ( b, 9UL );
         checkCapacity( b, 9UL );
         checkNonZeros( b, 9UL );

         if( b[0] !=  1U || b[1] !=  2U || b[2] !=   4U || b[3] !=   8U || b[4] != 16U ||
             b[5] != 32U || b[6] != 64U || b[7] != 128U || b[8] != 256U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Vector/scalar right-shift operation failed\n"
                << " Details:\n"
                << "   Result:\n" << b << "\n"
                << "   Expected result:\n( 1 2 4 8 16 32 64 128 256 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Vector/scalar right-shift assignment
      {
         blaze::DynamicVector<unsigned int> a{ 4U, 8U, 16U, 32U, 64U, 128U, 256U, 512U, 1024U };

         a >>= 2U;

         checkSize    ( a, 9UL );
         checkCapacity( a, 9UL );
         checkNonZeros( a, 9UL );

         if( a[0] !=  1U || a[1] !=  2U || a[2] !=   4U || a[3] !=   8U || a[4] != 16U ||
             a[5] != 32U || a[6] != 64U || a[7] != 128U || a[8] != 256U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Vector/scalar right-shift assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << a << "\n"
                << "   Expected result:\n( 1 2 4 8 16 32 64 128 256 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Vector/vector right-shift tests
   //=====================================================================================

   {
      test_ = "Vector/vector right-shift operator";

      // Vector/vector right-shift of an empty vector
      {
         blaze::DynamicVector<unsigned int> a;
         blaze::DynamicVector<unsigned int> b;

         blaze::DynamicVector<unsigned int> c( a >> b );

         checkSize    ( b, 0UL );
         checkCapacity( b, 0UL );
         checkNonZeros( b, 0UL );
      }

      // Vector/vector right-shift of a general vector
      {
         blaze::DynamicVector<unsigned int> a{ 4U, 8U, 16U, 32U, 64U, 128U, 256U, 512U, 1024U };
         blaze::DynamicVector<unsigned int> b{ 1U, 2U,  1U, 2U, 1U, 2U, 1U, 2U, 1U };

         blaze::DynamicVector<unsigned int> c( a >> b );

         checkSize    ( c, 9UL );
         checkCapacity( c, 9UL );
         checkNonZeros( c, 9UL );

         if( c[0] !=  2U || c[1] !=   2U || c[2] !=   8U || c[3] !=   8U || c[4] != 32U ||
             c[5] != 32U || c[6] != 128U || c[7] != 128U || c[8] != 512U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Vector/vector right-shift operation failed\n"
                << " Details:\n"
                << "   Result:\n" << c << "\n"
                << "   Expected result:\n( 2 2 8 8 32 32 128 128 512 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Vector/vector right-shift assignment
      {
         blaze::DynamicVector<unsigned int> a{ 4U, 8U, 16U, 32U, 64U, 128U, 256U, 512U, 1024U };
         blaze::DynamicVector<unsigned int> b{ 1U, 2U,  1U, 2U, 1U, 2U, 1U, 2U, 1U };

         a >>= b;

         checkSize    ( a, 9UL );
         checkCapacity( a, 9UL );
         checkNonZeros( a, 9UL );

         if( a[0] !=  2U || a[1] !=   2U || a[2] !=   8U || a[3] !=   8U || a[4] != 32U ||
             a[5] != 32U || a[6] != 128U || a[7] != 128U || a[8] != 512U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Vector/vector right-shift assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << a << "\n"
                << "   Expected result:\n( 2 2 8 8 32 32 128 128 512 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the bitwise AND operator for dense vectors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the bitwise AND operator for dense vectors. In case an error
// is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testBitand()
{
   //=====================================================================================
   // Vector/scalar bitwise AND tests
   //=====================================================================================

   {
      test_ = "Vector/scalar bitwise AND operator";

      // Vector/scalar bitwise AND of an empty vector
      {
         blaze::DynamicVector<unsigned int> a;

         blaze::DynamicVector<unsigned int> b( a & 7U );

         checkSize    ( b, 0UL );
         checkCapacity( b, 0UL );
         checkNonZeros( b, 0UL );
      }

      // Vector/scalar bitwise AND of a general vector
      {
         blaze::DynamicVector<unsigned int> a{ 8U, 9U, 10U, 11U, 12U, 13U, 14U, 15U, 16U };

         blaze::DynamicVector<unsigned int> b( a & 7U );

         checkSize    ( b, 9UL );
         checkCapacity( b, 9UL );
         checkNonZeros( b, 7UL );

         if( b[0] != 0U || b[1] != 1U || b[2] != 2U || b[3] != 3U || b[4] != 4U ||
             b[5] != 5U || b[6] != 6U || b[7] != 7U || b[8] != 0U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Vector/scalar bitwise AND operation failed\n"
                << " Details:\n"
                << "   Result:\n" << b << "\n"
                << "   Expected result:\n( 0 1 2 3 4 5 6 7 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Vector/scalar bitwise AND assignment
      {
         blaze::DynamicVector<unsigned int> a{ 8U, 9U, 10U, 11U, 12U, 13U, 14U, 15U, 16U };

         a &= 7U;

         checkSize    ( a, 9UL );
         checkCapacity( a, 9UL );
         checkNonZeros( a, 7UL );

         if( a[0] != 0U || a[1] != 1U || a[2] != 2U || a[3] != 3U || a[4] != 4U ||
             a[5] != 5U || a[6] != 6U || a[7] != 7U || a[8] != 0U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Vector/scalar bitwise AND assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << a << "\n"
                << "   Expected result:\n( 0 1 2 3 4 5 6 7 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Vector/vector bitwise AND tests
   //=====================================================================================

   {
      test_ = "Vector/vector bitwise AND operator";

      // Vector/vector bitwise AND of an empty vector
      {
         blaze::DynamicVector<unsigned int> a;
         blaze::DynamicVector<unsigned int> b;

         blaze::DynamicVector<unsigned int> c( a & b );

         checkSize    ( c, 0UL );
         checkCapacity( c, 0UL );
         checkNonZeros( c, 0UL );
      }

      // Vector/vector bitwise AND of a general vector
      {
         blaze::DynamicVector<unsigned int> a{ 8U, 9U, 10U, 11U, 12U, 13U, 14U, 15U, 16U };
         blaze::DynamicVector<unsigned int> b{ 7U, 5U,  7U, 5U, 7U, 5U, 7U, 5U, 7U };

         blaze::DynamicVector<unsigned int> c( a & b );

         checkSize    ( c, 9UL );
         checkCapacity( c, 9UL );
         checkNonZeros( c, 7UL );

         if( c[0] != 0U || c[1] != 1U || c[2] != 2U || c[3] != 1U || c[4] != 4U ||
             c[5] != 5U || c[6] != 6U || c[7] != 5U || c[8] != 0U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Vector/vector bitwise AND operation failed\n"
                << " Details:\n"
                << "   Result:\n" << c << "\n"
                << "   Expected result:\n( 0 1 2 1 4 5 6 5 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Vector/vector bitwise AND assignment
      {
         blaze::DynamicVector<unsigned int> a{ 8U, 9U, 10U, 11U, 12U, 13U, 14U, 15U, 16U };
         blaze::DynamicVector<unsigned int> b{ 7U, 5U,  7U, 5U, 7U, 5U, 7U, 5U, 7U };

         a &= b;

         checkSize    ( a, 9UL );
         checkCapacity( a, 9UL );
         checkNonZeros( a, 7UL );

         if( a[0] != 0U || a[1] != 1U || a[2] != 2U || a[3] != 1U || a[4] != 4U ||
             a[5] != 5U || a[6] != 6U || a[7] != 5U || a[8] != 0U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Vector/vector bitwise AND assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << a << "\n"
                << "   Expected result:\n( 0 1 2 1 4 5 6 5 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the bitwise OR operator for dense vectors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the bitwise OR operator for dense vectors. In case an error
// is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testBitor()
{
   //=====================================================================================
   // Vector/scalar bitwise OR tests
   //=====================================================================================

   {
      test_ = "Vector/scalar bitwise OR operator";

      // Vector/scalar bitwise OR of an empty vector
      {
         blaze::DynamicVector<unsigned int> a;

         blaze::DynamicVector<unsigned int> b( a | 7U );

         checkSize    ( b, 0UL );
         checkCapacity( b, 0UL );
         checkNonZeros( b, 0UL );
      }

      // Vector/scalar bitwise OR of a general vector
      {
         blaze::DynamicVector<unsigned int> a{ 8U, 9U, 10U, 11U, 12U, 13U, 14U, 15U, 16U };

         blaze::DynamicVector<unsigned int> b( a | 7U );

         checkSize    ( b, 9UL );
         checkCapacity( b, 9UL );
         checkNonZeros( b, 9UL );

         if( b[0] != 15U || b[1] != 15U || b[2] != 15U || b[3] != 15U || b[4] != 15U ||
             b[5] != 15U || b[6] != 15U || b[7] != 15U || b[8] != 23U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Vector/scalar bitwise OR operation failed\n"
                << " Details:\n"
                << "   Result:\n" << b << "\n"
                << "   Expected result:\n( 15 15 15 15 15 15 15 15 23 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Vector/scalar bitwise OR assignment
      {
         blaze::DynamicVector<unsigned int> a{ 8U, 9U, 10U, 11U, 12U, 13U, 14U, 15U, 16U };

         a |= 7U;

         checkSize    ( a, 9UL );
         checkCapacity( a, 9UL );
         checkNonZeros( a, 9UL );

         if( a[0] != 15U || a[1] != 15U || a[2] != 15U || a[3] != 15U || a[4] != 15U ||
             a[5] != 15U || a[6] != 15U || a[7] != 15U || a[8] != 23U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Vector/scalar bitwise OR assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << a << "\n"
                << "   Expected result:\n( 15 15 15 15 15 15 15 15 23 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Vector/vector bitwise OR tests
   //=====================================================================================

   {
      test_ = "Vector/vector bitwise OR operator";

      // Vector/vector bitwise OR of an empty vector
      {
         blaze::DynamicVector<unsigned int> a;
         blaze::DynamicVector<unsigned int> b;

         blaze::DynamicVector<unsigned int> c( a | b );

         checkSize    ( c, 0UL );
         checkCapacity( c, 0UL );
         checkNonZeros( c, 0UL );
      }

      // Vector/vector bitwise OR of a general vector
      {
         blaze::DynamicVector<unsigned int> a{ 8U, 9U, 10U, 11U, 12U, 13U, 14U, 15U, 16U };
         blaze::DynamicVector<unsigned int> b{ 7U, 5U,  7U, 5U, 7U, 5U, 7U, 5U, 7U };

         blaze::DynamicVector<unsigned int> c( a | b );

         checkSize    ( c, 9UL );
         checkCapacity( c, 9UL );
         checkNonZeros( c, 9UL );

         if( c[0] != 15U || c[1] != 13U || c[2] != 15U || c[3] != 15U || c[4] != 15U ||
             c[5] != 13U || c[6] != 15U || c[7] != 15U || c[8] != 23U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Vector/vector bitwise OR operation failed\n"
                << " Details:\n"
                << "   Result:\n" << c << "\n"
                << "   Expected result:\n( 15 13 15 15 15 13 15 15 23 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Vector/vector bitwise OR assignment
      {
         blaze::DynamicVector<unsigned int> a{ 8U, 9U, 10U, 11U, 12U, 13U, 14U, 15U, 16U };
         blaze::DynamicVector<unsigned int> b{ 7U, 5U,  7U, 5U, 7U, 5U, 7U, 5U, 7U };

         a |= b;

         checkSize    ( a, 9UL );
         checkCapacity( a, 9UL );
         checkNonZeros( a, 9UL );

         if( a[0] != 15U || a[1] != 13U || a[2] != 15U || a[3] != 15U || a[4] != 15U ||
             a[5] != 13U || a[6] != 15U || a[7] != 15U || a[8] != 23U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Vector/vector bitwise OR assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << a << "\n"
                << "   Expected result:\n( 15 13 15 15 15 13 15 15 23 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the bitwise XOR operator for dense vectors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the bitwise XOR operator for dense vectors. In case an error
// is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testBitxor()
{
   //=====================================================================================
   // Vector/scalar bitwise XOR tests
   //=====================================================================================

   {
      test_ = "Vector/scalar bitwise XOR operator";

      // Vector/scalar bitwise XOR of an empty vector
      {
         blaze::DynamicVector<unsigned int> a;

         blaze::DynamicVector<unsigned int> b( a ^ 7U );

         checkSize    ( b, 0UL );
         checkCapacity( b, 0UL );
         checkNonZeros( b, 0UL );
      }

      // Vector/scalar bitwise XOR of a general vector
      {
         blaze::DynamicVector<unsigned int> a{ 8U, 9U, 10U, 11U, 12U, 13U, 14U, 15U, 16U };

         blaze::DynamicVector<unsigned int> b( a ^ 7U );

         checkSize    ( b, 9UL );
         checkCapacity( b, 9UL );
         checkNonZeros( b, 9UL );

         if( b[0] != 15U || b[1] != 14U || b[2] != 13U || b[3] != 12U || b[4] != 11U ||
             b[5] != 10U || b[6] !=  9U || b[7] !=  8U || b[8] != 23U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Vector/scalar bitwise XOR operation failed\n"
                << " Details:\n"
                << "   Result:\n" << b << "\n"
                << "   Expected result:\n( 15 14 13 12 11 10 9 8 23 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Vector/scalar bitwise XOR assignment
      {
         blaze::DynamicVector<unsigned int> a{ 8U, 9U, 10U, 11U, 12U, 13U, 14U, 15U, 16U };

         a ^= 7U;

         checkSize    ( a, 9UL );
         checkCapacity( a, 9UL );
         checkNonZeros( a, 9UL );

         if( a[0] != 15U || a[1] != 14U || a[2] != 13U || a[3] != 12U || a[4] != 11U ||
             a[5] != 10U || a[6] !=  9U || a[7] !=  8U || a[8] != 23U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Vector/scalar bitwise XOR assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << a << "\n"
                << "   Expected result:\n( 15 14 13 12 11 10 9 8 23 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Vector/vector bitwise XOR tests
   //=====================================================================================

   {
      test_ = "Vector/vector bitwise XOR operator";

      // Vector/vector bitwise XOR of an empty vector
      {
         blaze::DynamicVector<unsigned int> a;
         blaze::DynamicVector<unsigned int> b;

         blaze::DynamicVector<unsigned int> c( a ^ b );

         checkSize    ( c, 0UL );
         checkCapacity( c, 0UL );
         checkNonZeros( c, 0UL );
      }

      // Vector/vector bitwise XOR of a general vector
      {
         blaze::DynamicVector<unsigned int> a{ 8U, 9U, 10U, 11U, 12U, 13U, 14U, 15U, 16U };
         blaze::DynamicVector<unsigned int> b{ 7U, 5U,  7U, 5U, 7U, 5U, 7U, 5U, 7U };

         blaze::DynamicVector<unsigned int> c( a ^ b );

         checkSize    ( c, 9UL );
         checkCapacity( c, 9UL );
         checkNonZeros( c, 9UL );

         if( c[0] != 15U || c[1] != 12U || c[2] != 13U || c[3] != 14U || c[4] != 11U ||
             c[5] !=  8U || c[6] !=  9U || c[7] != 10U || c[8] != 23U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Vector/vector bitwise XOR operation failed\n"
                << " Details:\n"
                << "   Result:\n" << c << "\n"
                << "   Expected result:\n( 15 12 13 14 11 8 9 10 23 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Vector/vector bitwise XOR assignment
      {
         blaze::DynamicVector<unsigned int> a{ 8U, 9U, 10U, 11U, 12U, 13U, 14U, 15U, 16U };
         blaze::DynamicVector<unsigned int> b{ 7U, 5U,  7U, 5U, 7U, 5U, 7U, 5U, 7U };

         a ^= b;

         checkSize    ( a, 9UL );
         checkCapacity( a, 9UL );
         checkNonZeros( a, 9UL );

         if( a[0] != 15U || a[1] != 12U || a[2] != 13U || a[3] != 14U || a[4] != 11U ||
             a[5] !=  8U || a[6] !=  9U || a[7] != 10U || a[8] != 23U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Vector/vector bitwise XOR assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << a << "\n"
                << "   Expected result:\n( 15 12 13 14 11 8 9 10 23 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the logical NOT operator for dense vectors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the logical NOT operator for dense vectors. In case an error
// is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testNot()
{
   test_ = "Vector/vector logical NOT operator";

   // Vector/vector logical NOT of an empty vector
   {
      blaze::DynamicVector<bool> a;
      blaze::DynamicVector<bool> b( !a );

      checkSize    ( b, 0UL );
      checkCapacity( b, 0UL );
      checkNonZeros( b, 0UL );
   }

   // Vector/vector logical NOT of a general vector
   {
      blaze::DynamicVector<bool> a{ false, true, false, true, false };
      blaze::DynamicVector<bool> b{ !a };

      checkSize    ( b, 5UL );
      checkCapacity( b, 5UL );
      checkNonZeros( b, 3UL );

      if( b[0] != true || b[1] != false || b[2] != true || b[3] != false || b[4] != true ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Vector logical NOT operation failed\n"
             << " Details:\n"
             << "   Result:\n" << b << "\n"
             << "   Expected result:\n( 1 0 1 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the logical AND operator for dense vectors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the logical AND operator for dense vectors. In case an error
// is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testAnd()
{
   test_ = "Vector/vector logical AND operator";

   // Vector/vector logical AND of an empty vector
   {
      blaze::DynamicVector<bool> a;
      blaze::DynamicVector<bool> b;

      blaze::DynamicVector<bool> c( a && b );

      checkSize    ( c, 0UL );
      checkCapacity( c, 0UL );
      checkNonZeros( c, 0UL );
   }

   // Vector/vector logical AND of a general vector
   {
      blaze::DynamicVector<bool> a{ true, false, true, false, true };
      blaze::DynamicVector<bool> b{ true, true, false, false, true };

      blaze::DynamicVector<bool> c( a && b );

      checkSize    ( c, 5UL );
      checkCapacity( c, 5UL );
      checkNonZeros( c, 2UL );

      if( c[0] != true || c[1] != false || c[2] != false || c[3] != false || c[4] != true ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Vector/vector logical AND operation failed\n"
             << " Details:\n"
             << "   Result:\n" << c << "\n"
             << "   Expected result:\n( 1 0 0 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the logical OR operator for dense vectors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the logical OR operator for dense vectors. In case an error
// is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testOr()
{
   test_ = "Vector/vector logical OR operator";

   // Vector/vector logical OR of an empty vector
   {
      blaze::DynamicVector<bool> a;
      blaze::DynamicVector<bool> b;

      blaze::DynamicVector<bool> c( a || b );

      checkSize    ( c, 0UL );
      checkCapacity( c, 0UL );
      checkNonZeros( c, 0UL );
   }

   // Vector/vector logical OR of a general vector
   {
      blaze::DynamicVector<bool> a{ true, false, true, false, true };
      blaze::DynamicVector<bool> b{ true, true, false, false, true };

      blaze::DynamicVector<bool> c( a || b );

      checkSize    ( c, 5UL );
      checkCapacity( c, 5UL );
      checkNonZeros( c, 4UL );

      if( c[0] != true || c[1] != true || c[2] != true || c[3] != false || c[4] != true ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Vector/vector logical OR operation failed\n"
             << " Details:\n"
             << "   Result:\n" << c << "\n"
             << "   Expected result:\n( 1 1 1 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c generate() functions for dense vectors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c generate() functions for dense vectors. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testGenerate()
{
   test_ = "generate() function";

   // Empty integer vector
   {
      const blaze::DynamicVector<int> vec(
         blaze::generate( 0UL, []( size_t ){ return 2; } ) );

      const blaze::DynamicVector<int> ref{};

      if( vec != ref ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Generating empty integer vector failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n" << ref << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Single element integer vector ( 2 )
   {
      const blaze::DynamicVector<int> vec(
         blaze::generate( 1UL, []( size_t ){ return 2; } ) );

      const blaze::DynamicVector<int> ref{ 2 };

      if( vec != ref ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Generating single element integer vector failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n" << ref << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Uniform integer vector ( 2, 2, 2, 2, 2 )
   {
      const blaze::DynamicVector<int> vec(
         blaze::generate( 5UL, []( size_t ){ return 2; } ) );

      const blaze::DynamicVector<int> ref{ 2, 2, 2, 2, 2 };

      if( vec != ref ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Generating uniform integer vector failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n" << ref << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Linearly spaced float vector ( 2.1, 3.2, 4.3, 5.4 )
   {
      const blaze::DynamicVector<float> vec(
         blaze::generate( 4UL, []( size_t index ){ return 2.1F + 1.1F*index; } ) );

      const blaze::DynamicVector<float> ref{ 2.1F, 3.2F, 4.3F, 5.4F };

      if( vec != ref ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Generating linearly spaced float vector failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n" << ref << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Logarithmically spaced double vector ( 10.0, 100.0, 1000.0, 10000.0 )
   {
      const blaze::DynamicVector<double> vec(
         blaze::generate( 4UL, []( size_t index ){ return blaze::exp10( 1.0 + 1.0*index ); } ) );

      const blaze::DynamicVector<double> ref{ 10.0, 100.0, 1000.0, 10000.0 };

      if( vec != ref ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Generating logarithmically spaced double vector failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n" << ref << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Vector of vectors
   {
      using VT = blaze::StaticVector<int,2UL>;

      const blaze::DynamicVector<VT> vec(
         blaze::generate( 4UL, []( size_t index ) { return evaluate( VT{ 1, 2 } + index ); } ) );

      const blaze::DynamicVector<VT> ref{ { 1, 2 }, { 2, 3 }, { 3, 4 }, { 4, 5 } };

      if( vec != ref ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Generating vector of vectors failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n" << ref << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c linspace() functions for dense vectors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c linspace() functions for dense vectors. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testLinspace()
{
   test_ = "linspace() function";

   // Empty integer vector
   {
      const blaze::DynamicVector<int> vec( blaze::linspace( 0UL, 2, 5 ) );
      const blaze::DynamicVector<int> ref{};

      if( vec != ref ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Generating empty integer vector failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n" << ref << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Single element integer vector ( 5 )
   {
      const blaze::DynamicVector<int> vec( blaze::linspace( 1UL, 2, 5 ) );
      const blaze::DynamicVector<int> ref{ 5 };

      if( vec != ref ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Generating single element integer vector failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n" << ref << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Uniform integer vector ( 2, 2, 2, 2, 2 )
   {
      const blaze::DynamicVector<int> vec( blaze::linspace( 4UL, 2, 2 ) );
      const blaze::DynamicVector<int> ref{ 2, 2, 2, 2 };

      if( vec != ref ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Generating uniform integer vector failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n" << ref << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Linearly spaced float vector ( 2.1, 3.2, 4.3, 5.4 )
   {
      const blaze::DynamicVector<float> vec( blaze::linspace( 4UL, 2.1F, 5.4F ) );
      const blaze::DynamicVector<float> ref{ 2.1F, 3.2F, 4.3F, 5.4F };

      if( vec != ref ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Generating linearly spaced float vector failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n" << ref << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Linearly spaced float vector ( 5.4, 4.3, 3.2, 2.1 )
   {
      const blaze::DynamicVector<float> vec( blaze::linspace( 4UL, 5.4F, 2.1F ) );
      const blaze::DynamicVector<float> ref{ 5.4F, 4.3F, 3.2F, 2.1F };

      if( vec != ref ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Generating linearly spaced float vector failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n" << ref << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Vector of vectors
   {
      using VT = blaze::StaticVector<int,2UL>;

      const blaze::DynamicVector<VT> vec( blaze::linspace( 4UL, VT{ 1, 2 }, VT{ 4, 5 } ) );
      const blaze::DynamicVector<VT> ref{ { 1, 2 }, { 2, 3 }, { 3, 4 }, { 4, 5 } };

      if( vec != ref ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Generating vector of vectors failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n" << ref << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c logspace() functions for dense vectors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c logspace() functions for dense vectors. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testLogspace()
{
   test_ = "logspace() function";

   // Empty integer vector
   {
      const blaze::DynamicVector<int> vec( blaze::logspace( 0UL, 0, 3 ) );
      const blaze::DynamicVector<int> ref{};

      if( vec != ref ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Generating empty integer vector failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n" << ref << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Single element integer vector ( 1000 )
   {
      const blaze::DynamicVector<int> vec( blaze::logspace( 1UL, 0, 3 ) );
      const blaze::DynamicVector<int> ref{ 1000 };

      if( vec != ref ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Generating single element integer vector failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n" << ref << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Uniform integer vector ( 10, 10, 10, 10 )
   {
      const blaze::DynamicVector<int> vec( blaze::logspace( 4UL, 1, 1 ) );
      const blaze::DynamicVector<int> ref{ 10, 10, 10, 10 };

      if( vec != ref ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Generating uniform integer vector failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n" << ref << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Logarithmically spaced float vector ( 1.0, 10.0, 100.0, 1000.0 )
   {
      const blaze::DynamicVector<float> vec( blaze::logspace( 4UL, 0.0F, 3.0F ) );
      const blaze::DynamicVector<float> ref{ 1.0F, 10.0F, 100.0F, 1000.0F };

      if( vec != ref ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Generating logarithmically spaced float vector failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n" << ref << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Logarithmically spaced float vector ( 1000.0, 100.0, 10.0, 1.0 )
   {
      const blaze::DynamicVector<float> vec( blaze::logspace( 4UL, 3.0F, 0.0F ) );
      const blaze::DynamicVector<float> ref{ 1000.0F, 100.0F, 10.0F, 1.0F };

      if( vec != ref ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Generating logarithmically spaced float vector failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n" << ref << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Vector of vectors
   {
      using VT = blaze::StaticVector<int,2UL>;

      const blaze::DynamicVector<VT> vec( blaze::logspace( 4UL, VT{ 0, 1 }, VT{ 3, 4 } ) );
      const blaze::DynamicVector<VT> ref{ { 1, 10 }, { 10, 100 }, { 100, 1000 }, { 1000, 10000 } };

      if( vec != ref ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Generating vector of vectors failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n" << ref << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c uniform() functions for dense vectors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c uniform() functions for dense vectors. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testUniform()
{
   test_ = "uniform() function";

   // Empty integer vector
   {
      const blaze::DynamicVector<int> vec( blaze::uniform( 0UL, 5 ) );
      const blaze::DynamicVector<int> ref{};

      if( vec != ref ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Generating empty integer vector failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n" << ref << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Single element integer vector ( 5 )
   {
      const blaze::DynamicVector<int> vec( blaze::uniform( 1UL, 5 ) );
      const blaze::DynamicVector<int> ref{ 5 };

      if( vec != ref ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Generating single element integer vector failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n" << ref << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Uniform float vector ( 2.1, 2.1, 2.1, 2.1 )
   {
      const blaze::DynamicVector<float> vec( blaze::uniform( 4UL, 2.1f ) );
      const blaze::DynamicVector<float> ref{ 2.1F, 2.1F, 2.1F, 2.1F };

      if( vec != ref ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Generating uniform float vector failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n" << ref << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Uniform vector of vectors
   {
      using VT = blaze::StaticVector<int,2UL>;

      const blaze::DynamicVector<VT> vec( blaze::uniform( 4UL, VT{ 1, 2 } ) );
      const blaze::DynamicVector<VT> ref{ { 1, 2 }, { 1, 2 }, { 1, 2 }, { 1, 2 } };

      if( vec != ref ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Generating vector of vectors failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n" << ref << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c zero() functions for dense vectors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c zero() functions for dense vectors. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void GeneralTest::testZero()
{
   test_ = "zero() function";

   // Empty integer vector
   {
      const blaze::DynamicVector<int> vec( blaze::zero<int>( 0UL ) );
      const blaze::DynamicVector<int> ref{};

      if( vec != ref ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Generating empty integer vector failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n" << ref << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Single element integer vector ( 0 )
   {
      const blaze::DynamicVector<int> vec( blaze::zero<int>( 1UL ) );
      const blaze::DynamicVector<int> ref{ 0 };

      if( vec != ref ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Generating single element integer vector failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n" << ref << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Uniform float vector ( 0.0, 0.0, 0.0, 0.0 )
   {
      const blaze::DynamicVector<float> vec( blaze::zero<float>( 4UL ) );
      const blaze::DynamicVector<float> ref{ 0.0F, 0.0F, 0.0F, 0.0F };

      if( vec != ref ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Generating zero float vector failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n" << ref << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Uniform vector of vectors
   {
      using VT = blaze::StaticVector<int,2UL>;

      const blaze::DynamicVector<VT> vec( blaze::uniform( 4UL, VT{} ) );
      const blaze::DynamicVector<VT> ref{ VT{}, VT{}, VT{}, VT{} };

      if( vec != ref ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Generating vector of vectors failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n" << ref << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************

} // namespace densevector

} // namespace vectors

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
   std::cout << "   Running general DenseVector operation test..." << std::endl;

   try
   {
      RUN_DENSEVECTOR_GENERAL_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during general DenseVector operation test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
