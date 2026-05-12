//=================================================================================================
/*!
//  \file src/mathtest/vectors/sparsevector/ZeroTest.cpp
//  \brief Source file for the zero SparseVector operation test
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
#include <blaze/math/sparse/SparseVector.h>
#include <blaze/math/ZeroVector.h>
#include <blazetest/mathtest/vectors/sparsevector/ZeroTest.h>

#ifdef BLAZE_USE_HPX_THREADS
#  include <hpx/hpx_main.hpp>
#endif


namespace blazetest {

namespace mathtest {

namespace vectors {

namespace sparsevector {

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
/*!\brief Test of the \c isnan() function for sparse vectors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isnan() function for sparse vectors. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void ZeroTest::testIsNan()
{
   test_ = "isnan() function";

   // isnan with 0-dimensional vector
   {
      blaze::ZeroVector<float,blaze::rowVector> vec;

      if( blaze::isnan( vec ) != false ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isnan evaluation\n"
             << " Details:\n"
             << "   Vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // isnan with empty 9-dimensional vector (non-NAN)
   {
      blaze::ZeroVector<float,blaze::rowVector> vec( 9UL );

      if( blaze::isnan( vec ) != false ) {
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
/*!\brief Test of the \c isinf() function for sparse vectors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isinf() function for sparse vectors. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void ZeroTest::testIsInf()
{
   test_ = "isinf() function";

   // isinf with 0-dimensional vector
   {
      blaze::ZeroVector<float,blaze::rowVector> vec;

      if( blaze::isinf( vec ) != false ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isinf evaluation\n"
             << " Details:\n"
             << "   Vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // isinf with empty 9-dimensional vector (non-inf)
   {
      blaze::ZeroVector<float,blaze::rowVector> vec( 9UL );

      if( blaze::isinf( vec ) != false ) {
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
/*!\brief Test of the \c isfinite() function for sparse vectors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isfinite() function for sparse vectors. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void ZeroTest::testIsFinite()
{
   test_ = "isfinite() function";

   // isfinite with 0-dimensional vector
   {
      blaze::ZeroVector<float,blaze::rowVector> vec;

      if( blaze::isfinite( vec ) != true ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isfinite evaluation\n"
             << " Details:\n"
             << "   Vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // isfinite with empty 9-dimensional vector (non-finite)
   {
      blaze::ZeroVector<float,blaze::rowVector> vec( 9UL );

      if( blaze::isfinite( vec ) != true ) {
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
/*!\brief Test of the \c isUniform() function for sparse vectors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isUniform() function for sparse vectors. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void ZeroTest::testIsUniform()
{
   test_ = "isUniform() function";

   // Uniform vector (0-dimensional)
   {
      blaze::ZeroVector<int,blaze::rowVector> vec;

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
      blaze::ZeroVector<int,blaze::rowVector> vec( 1UL );

      if( blaze::isUniform( vec ) != true ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isUniform evaluation\n"
             << " Details:\n"
             << "   Vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Uniform 5-dimensional vector
   {
      blaze::ZeroVector<int,blaze::rowVector> vec( 5UL );

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
/*!\brief Test of the \c isZero() function for sparse vectors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isZero() function for sparse vectors. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void ZeroTest::testIsZero()
{
   test_ = "isZero() function";

   // Zero vector (0-dimensional)
   {
      blaze::ZeroVector<int,blaze::rowVector> vec;

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
      blaze::ZeroVector<int,blaze::rowVector> vec( 1UL );

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
      blaze::ZeroVector<int,blaze::rowVector> vec( 5UL );

      if( blaze::isZero( vec ) != true ) {
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
/*!\brief Test of the \c mean() function for sparse vectors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c mean() function for sparse vectors. In case an error
// is detected, a \a std::runtime_error exception is thrown.
*/
void ZeroTest::testMean()
{
   test_ = "mean() function";

   {
      blaze::ZeroVector<double,blaze::rowVector> vec( 5UL );

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

   try {
      blaze::ZeroVector<int,blaze::rowVector> vec;

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
/*!\brief Test of the \c var() function for sparse vectors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c var() function for sparse vectors. In case an error
// is detected, a \a std::runtime_error exception is thrown.
*/
void ZeroTest::testVar()
{
   test_ = "var() function";

   {
      blaze::ZeroVector<int,blaze::rowVector> vec( 5UL );

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
      blaze::ZeroVector<int,blaze::rowVector> vec;

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
      blaze::ZeroVector<int,blaze::rowVector> vec( 1UL );

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
/*!\brief Test of the \c stddev() function for sparse vectors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c stddev() function for sparse vectors. In case an error
// is detected, a \a std::runtime_error exception is thrown.
*/
void ZeroTest::testStdDev()
{
   test_ = "stddev() function";

   {
      blaze::ZeroVector<int,blaze::rowVector> vec( 5UL );

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
      blaze::ZeroVector<int,blaze::rowVector> vec;

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
      blaze::ZeroVector<int,blaze::rowVector> vec( 1UL );

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

} // namespace sparsevector

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
   std::cout << "   Running zero SparseVector operation test..." << std::endl;

   try
   {
      RUN_SPARSEVECTOR_ZERO_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during zero SparseVector operation test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
