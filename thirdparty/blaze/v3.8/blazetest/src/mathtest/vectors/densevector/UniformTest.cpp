//=================================================================================================
/*!
//  \file src/mathtest/vectors/densevector/UniformTest.cpp
//  \brief Source file for the uniform DenseVector operation test
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
#include <blaze/math/UniformVector.h>
#include <blazetest/mathtest/vectors/densevector/UniformTest.h>

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
/*!\brief Constructor for the UniformTest class test.
//
// \exception std::runtime_error Operation error detected.
*/
UniformTest::UniformTest()
{
   testIsNan();
   testIsInf();
   testIsFinite();
   testIsUniform();
   testIsZero();
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
/*!\brief Test of the \c isnan() function for dense vectors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isnan() function for dense vectors. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void UniformTest::testIsNan()
{
   test_ = "isnan() function";

   // isnan with 0-dimensional uniform vector (non-NAN)
   {
      blaze::UniformVector<float,blaze::rowVector> vec;

      if( blaze::isnan( vec ) != false ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isnan evaluation\n"
             << " Details:\n"
             << "   Vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // isnan with 6-dimensional uniform vector (non-NAN)
   {
      blaze::UniformVector<float,blaze::rowVector> vec( 6UL, 1.2 );

      if( blaze::isnan( vec ) != false ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isnan evaluation\n"
             << " Details:\n"
             << "   Vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // isnan with 6-dimensional uniform vector (NAN)
   {
      blaze::UniformVector<float,blaze::rowVector> vec( 6UL, NAN );

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
void UniformTest::testIsInf()
{
   test_ = "isinf() function";

   // isinf with 0-dimensional uniform vector (non-inf)
   {
      blaze::UniformVector<float,blaze::rowVector> vec;

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
      blaze::UniformVector<float,blaze::rowVector> vec( 6UL, 1.2F );

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
      blaze::UniformVector<float,blaze::rowVector> vec( 6UL, HUGE_VALF );

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
void UniformTest::testIsFinite()
{
   test_ = "isfinite() function";

   // isfinite with 0-dimensional vector (finite)
   {
      blaze::UniformVector<float,blaze::rowVector> vec;

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
      blaze::UniformVector<float,blaze::rowVector> vec( 6UL, 1.2F );

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
      blaze::UniformVector<float,blaze::rowVector> vec( 6UL, NAN );

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
      blaze::UniformVector<float,blaze::rowVector> vec( 6UL, HUGE_VALF );

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
void UniformTest::testIsUniform()
{
   test_ = "isUniform() function";

   // Uniform vector (0-dimensional)
   {
      blaze::UniformVector<int,blaze::rowVector> vec;

      if( blaze::isUniform( vec ) != true ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isUniform evaluation\n"
             << " Details:\n"
             << "   Vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Uniform vector (1-dimensional, default value)
   {
      blaze::UniformVector<int,blaze::rowVector> vec( 1UL );

      if( blaze::isUniform( vec ) != true ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isUniform evaluation\n"
             << " Details:\n"
             << "   Vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Uniform vector (1-dimensional, non-default value)
   {
      blaze::UniformVector<int,blaze::rowVector> vec( 1UL, 5 );

      if( blaze::isUniform( vec ) != true ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isUniform evaluation\n"
             << " Details:\n"
             << "   Vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Uniform vector (5-dimensional, default values)
   {
      blaze::UniformVector<int,blaze::rowVector> vec{ 5UL };

      if( blaze::isUniform( vec ) != true ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isUniform evaluation\n"
             << " Details:\n"
             << "   Vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Uniform vector (5-dimensional, non-default values)
   {
      blaze::UniformVector<int,blaze::rowVector> vec{ 5UL, 5 };

      if( blaze::isUniform( vec ) != true ) {
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
void UniformTest::testIsZero()
{
   test_ = "isZero() function";

   // Zero vector (0-dimensional)
   {
      blaze::UniformVector<int,blaze::rowVector> vec;

      if( blaze::isZero( vec ) != true ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isZero evaluation\n"
             << " Details:\n"
             << "   Vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Zero vector (1-dimensional, default values)
   {
      blaze::UniformVector<int,blaze::rowVector> vec( 1UL );

      if( blaze::isZero( vec ) != true ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isZero evaluation\n"
             << " Details:\n"
             << "   Vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Zero vector (1-dimensional, non-default values)
   {
      blaze::UniformVector<int,blaze::rowVector> vec( 1UL, 5 );

      if( blaze::isZero( vec ) != false ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isZero evaluation\n"
             << " Details:\n"
             << "   Vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Zero vector (5-dimensional, default values)
   {
      blaze::UniformVector<int,blaze::rowVector> vec( 5UL );

      if( blaze::isZero( vec ) != true ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isZero evaluation\n"
             << " Details:\n"
             << "   Vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Non-zero vector (5-dimensional, non-default values)
   {
      blaze::UniformVector<int,blaze::rowVector> vec( 5UL, 5 );

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
/*!\brief Test of the \c mean() function for dense vectors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c mean() function for dense vectors. In case an error
// is detected, a \a std::runtime_error exception is thrown.
*/

void UniformTest::testMean()
{
   test_ = "mean() function";

   {
      blaze::UniformVector<double,blaze::rowVector> vec( 5UL, 0.0 );

      const double mean = blaze::mean( vec );

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

   {
      blaze::UniformVector<double,blaze::rowVector> vec( 5UL, 4.0 );

      const double mean = blaze::mean( vec );

      if( mean != 4.0 ) {
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
      blaze::UniformVector<double,blaze::rowVector> vec;

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
void UniformTest::testVar()
{
   test_ = "var() function";

   {
      blaze::UniformVector<double,blaze::rowVector> vec( 5UL, 0.0 );

      const double var = blaze::var( vec );

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

   {
      blaze::UniformVector<double,blaze::rowVector> vec( 5UL, 4.0 );

      const double var = blaze::var( vec );

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
      blaze::UniformVector<double,blaze::rowVector> vec;

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
      blaze::UniformVector<double,blaze::rowVector> vec( 1UL );

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
void UniformTest::testStdDev()
{
   test_ = "stddev() function";

   {
      blaze::UniformVector<double,blaze::rowVector> vec( 5UL, 0.0 );

      const double stddev = blaze::stddev( vec );

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

   {
      blaze::UniformVector<double,blaze::rowVector> vec( 5UL, 4.0 );

      const double stddev = blaze::stddev( vec );

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
      blaze::UniformVector<double,blaze::rowVector> vec;

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
      blaze::UniformVector<double,blaze::rowVector> vec( 1UL );

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
   std::cout << "   Running uniform DenseVector operation test..." << std::endl;

   try
   {
      RUN_DENSEVECTOR_UNIFORM_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during uniform DenseVector operation test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
