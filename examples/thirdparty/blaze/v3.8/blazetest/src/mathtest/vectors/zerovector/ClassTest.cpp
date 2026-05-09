//=================================================================================================
/*!
//  \file src/mathtest/vectors/zerovector/ClassTest.cpp
//  \brief Source file for the ZeroVector class test
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
#include <blaze/math/DynamicVector.h>
#include <blaze/math/shims/Equal.h>
#include <blaze/util/Complex.h>
#include <blaze/util/Random.h>
#include <blazetest/mathtest/RandomMaximum.h>
#include <blazetest/mathtest/RandomMinimum.h>
#include <blazetest/mathtest/vectors/zerovector/ClassTest.h>

#ifdef BLAZE_USE_HPX_THREADS
#  include <hpx/hpx_main.hpp>
#endif


namespace blazetest {

namespace mathtest {

namespace vectors {

namespace zerovector {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the ZeroVector class test.
//
// \exception std::runtime_error Operation error detected.
*/
ClassTest::ClassTest()
{
   testConstructors();
   testAssignment();
   testSubscript();
   testAt();
   testIterator();
   testNonZeros();
   testReset();
   testClear();
   testResize();
   testSwap();
   testFind();
   testLowerBound();
   testUpperBound();
   testIsDefault();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the ZeroVector constructors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all constructors of the ZeroVector class template. In case
// an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testConstructors()
{
   //=====================================================================================
   // Default constructor
   //=====================================================================================

   {
      test_ = "ZeroVector default constructor";

      blaze::ZeroVector<int,blaze::rowVector> z;

      checkSize    ( z, 0UL );
      checkNonZeros( z, 0UL );
   }


   //=====================================================================================
   // Size constructor
   //=====================================================================================

   {
      test_ = "ZeroVector size constructor (size 0)";

      blaze::ZeroVector<int,blaze::rowVector> z( 0UL );

      checkSize    ( z, 0UL );
      checkNonZeros( z, 0UL );
   }

   {
      test_ = "ZeroVector size constructor (size 5)";

      blaze::ZeroVector<int,blaze::rowVector> z( 5UL );

      checkSize    ( z, 5UL );
      checkNonZeros( z, 0UL );
   }


   //=====================================================================================
   // Copy constructor
   //=====================================================================================

   {
      test_ = "ZeroVector copy constructor (size 0)";

      blaze::ZeroVector<int,blaze::rowVector> z1( 0UL );
      blaze::ZeroVector<int,blaze::rowVector> z2( z1 );

      checkSize    ( z2, 0UL );
      checkNonZeros( z2, 0UL );
   }

   {
      test_ = "ZeroVector copy constructor (size 7)";

      blaze::ZeroVector<int,blaze::rowVector> z1( 7UL );
      blaze::ZeroVector<int,blaze::rowVector> z2( z1 );

      checkSize    ( z2, 7UL );
      checkCapacity( z2, 0UL );
      checkNonZeros( z2, 0UL );

      if( z2[0] != 0 || z2[1] != 0 || z2[2] != 0 || z2[3] != 0 ||
          z2[4] != 0 || z2[5] != 0 || z2[6] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << z2 << "\n"
             << "   Expected result:\n( 0 0 0 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Move constructor
   //=====================================================================================

   {
      test_ = "ZeroVector move constructor (size 0)";

      blaze::ZeroVector<int,blaze::rowVector> z1( 0UL );
      blaze::ZeroVector<int,blaze::rowVector> z2( std::move( z1 ) );

      checkSize    ( z2, 0UL );
      checkNonZeros( z2, 0UL );
   }

   {
      test_ = "ZeroVector move constructor (size 7)";

      blaze::ZeroVector<int,blaze::rowVector> z1( 7UL );
      blaze::ZeroVector<int,blaze::rowVector> z2( std::move( z1 ) );

      checkSize    ( z2, 7UL );
      checkCapacity( z2, 0UL );
      checkNonZeros( z2, 0UL );

      if( z2[0] != 0 || z2[1] != 0 || z2[2] != 0 || z2[3] != 0 ||
          z2[4] != 0 || z2[5] != 0 || z2[6] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << z2 << "\n"
             << "   Expected result:\n( 0 0 0 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Dense vector constructor
   //=====================================================================================

   {
      test_ = "ZeroVector dense vector constructor";

      blaze::DynamicVector<int,blaze::rowVector> z1{ 0, 0, 0, 0, 0 };
      blaze::ZeroVector<int,blaze::rowVector> z2( z1 );

      checkSize    ( z2, 5UL );
      checkCapacity( z2, 0UL );
      checkNonZeros( z2, 0UL );

      if( z2[0] != 0 || z2[1] != 0 || z2[2] != 0 || z2[3] != 0 || z2[4] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << z2 << "\n"
             << "   Expected result:\n( 0 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "ZeroVector dense vector constructor (non-zero)";

      blaze::DynamicVector<int,blaze::rowVector> z1{ 0, 0, 1, 0, 0 };

      try {
         blaze::ZeroVector<int,blaze::rowVector> z2( z1 );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of non-zero ZeroVector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << z2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Sparse vector constructor
   //=====================================================================================

   {
      test_ = "ZeroVector sparse vector assignment";

      blaze::CompressedVector<int,blaze::columnVector> z1{ 0, 0, 0, 0, 0, 0, 0 };
      blaze::ZeroVector<int,blaze::rowVector> z2( trans( z1 ) );

      checkSize    ( z2, 7UL );
      checkNonZeros( z2, 0UL );

      if( z2[0] != 0 || z2[1] != 0 || z2[2] != 0 || z2[3] != 0 ||
          z2[4] != 0 || z2[5] != 0 || z2[6] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << z2 << "\n"
             << "   Expected result:\n( 0 0 0 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "ZeroVector dense vector constructor (non-zero)";

      blaze::CompressedVector<int,blaze::columnVector> z1{ 0, 0, 0, 1, 0, 0, 0 };

      try {
         blaze::ZeroVector<int,blaze::rowVector> z2( trans( z1 ) );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of non-zero ZeroVector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << z2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the ZeroVector assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all assignment operators of the ZeroVector class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testAssignment()
{
   //=====================================================================================
   // Copy assignment
   //=====================================================================================

   {
      test_ = "ZeroVector copy assignment";

      blaze::ZeroVector<int,blaze::rowVector> z1( 7UL );
      blaze::ZeroVector<int,blaze::rowVector> z2;
      z2 = z1;

      checkSize    ( z2, 7UL );
      checkCapacity( z2, 0UL );
      checkNonZeros( z2, 0UL );

      if( z2[0] != 0 || z2[1] != 0 || z2[2] != 0 || z2[3] != 0 ||
          z2[4] != 0 || z2[5] != 0 || z2[6] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << z2 << "\n"
             << "   Expected result:\n( 0 0 0 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "ZeroVector copy assignment stress test";

      using RandomVectorType = blaze::ZeroVector<int,blaze::rowVector>;

      blaze::ZeroVector<int,blaze::rowVector> z1;

      for( size_t i=0UL; i<100UL; ++i )
      {
         const size_t size( blaze::rand<size_t>( 0UL, 20UL ) );
         const RandomVectorType z2( blaze::rand<RandomVectorType>( size ) );

         z1 = z2;

         if( z1 != z2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << z1 << "\n"
                << "   Expected result:\n" << z2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Move assignment
   //=====================================================================================

   {
      test_ = "ZeroVector move assignment";

      blaze::ZeroVector<int,blaze::rowVector> z1( 7UL );
      blaze::ZeroVector<int,blaze::rowVector> z2( 4UL );

      z2 = std::move( z1 );

      checkSize    ( z2, 7UL );
      checkCapacity( z2, 0UL );
      checkNonZeros( z2, 0UL );

      if( z2[0] != 0 || z2[1] != 0 || z2[2] != 0 || z2[3] != 0 ||
          z2[4] != 0 || z2[5] != 0 || z2[6] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << z2 << "\n"
             << "   Expected result:\n( 0 0 0 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Dense vector assignment
   //=====================================================================================

   {
      test_ = "ZeroVector dense vector assignment";

      blaze::DynamicVector<int,blaze::rowVector> z1{ 0, 0, 0, 0, 0 };
      blaze::ZeroVector<int,blaze::rowVector> z2;
      z2 = z1;

      checkSize    ( z2, 5UL );
      checkNonZeros( z2, 0UL );

      if( z2[0] != 0 || z2[1] != 0 || z2[2] != 0 || z2[3] != 0 || z2[4] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << z2 << "\n"
             << "   Expected result:\n( 0 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "ZeroVector dense vector assignment (non-zero)";

      blaze::DynamicVector<int,blaze::rowVector> z1{ 0, 0, 1, 0, 0 };
      blaze::ZeroVector<int,blaze::rowVector> z2;

      try {
         z2 = z1;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-zero vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << z2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Sparse vector assignment
   //=====================================================================================

   {
      test_ = "ZeroVector sparse vector assignment";

      blaze::CompressedVector<int,blaze::columnVector> z1{ 0, 0, 0, 0, 0, 0, 0 };
      blaze::ZeroVector<int,blaze::rowVector> z2;

      z2 = trans( z1 );

      checkSize    ( z2, 7UL );
      checkNonZeros( z2, 0UL );

      if( z2[0] != 0 || z2[1] != 0 || z2[2] != 0 || z2[3] != 0 ||
          z2[4] != 0 || z2[5] != 0 || z2[6] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << z2 << "\n"
             << "   Expected result:\n( 0 0 0 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "ZeroVector sparse vector assignment (non-zero)";

      blaze::CompressedVector<int,blaze::columnVector> z1{ 0, 0, 0, 1, 0, 0, 0 };
      blaze::ZeroVector<int,blaze::rowVector> z2;

      try {
         z2 = trans( z1 );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-zero vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << z2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the ZeroVector subscript operator.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of adding and accessing elements via the subscript operator
// of the ZeroVector class template. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
void ClassTest::testSubscript()
{
   test_ = "ZeroVector::operator[]";

   blaze::ZeroVector<int,blaze::rowVector> z( 7UL );

   checkSize    ( z, 7UL );
   checkCapacity( z, 0UL );
   checkNonZeros( z, 0UL );

   if( z[0] != 0 || z[1] != 0 || z[2] != 0 || z[3] != 0 ||
       z[4] != 0 || z[5] != 0 || z[6] != 0 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << z << "\n"
          << "   Expected result:\n( 0 0 0 0 0 0 0 )\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c at() member function of the ZeroVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of adding and accessing elements via the \c at() member function
// of the ZeroVector class template. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
void ClassTest::testAt()
{
   test_ = "ZeroVector::at()";

   blaze::ZeroVector<int,blaze::rowVector> z( 7UL );

   checkSize    ( z, 7UL );
   checkCapacity( z, 0UL );
   checkNonZeros( z, 0UL );

   if( z.at(0) != 0 || z.at(1) != 0 || z.at(2) != 0 || z.at(3) != 0 ||
       z.at(4) != 0 || z.at(5) != 0 || z.at(6) != 0 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Access via at() function failed\n"
          << " Details:\n"
          << "   Result:\n" << z << "\n"
          << "   Expected result:\n( 0 0 0 0 0 0 0 )\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the ZeroVector iterator implementation.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the iterator implementation of the ZeroVector class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testIterator()
{
   using VectorType    = blaze::ZeroVector<int>;
   using ConstIterator = VectorType::ConstIterator;

   VectorType z( 4UL );

   // Testing the ConstIterator default constructor
   {
      test_ = "ConstIterator default constructor";

      ConstIterator it{};

      if( it != ConstIterator() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed iterator default constructor\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Counting the number of elements via ConstIterator (end-begin)
   {
      test_ = "ConstIterator subtraction (end-begin)";

      const ptrdiff_t number( cend( z ) - cbegin( z ) );

      if( number != 0L ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid number of elements detected\n"
             << " Details:\n"
             << "   Number of elements         : " << number << "\n"
             << "   Expected number of elements: 2\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing ConstIterator comparison
   {
      test_ = "ConstIterator comparison";

      ConstIterator it ( cbegin( z ) );
      ConstIterator end( cend( z ) );

      if( it != end ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator comparison failed\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c nonZeros() member function of the ZeroVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c nonZeros() member function of the ZeroVector class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testNonZeros()
{
   test_ = "ZeroZero::nonZeros()";

   blaze::ZeroVector<int,blaze::rowVector> z( 7UL );

   checkSize    ( z, 7UL );
   checkCapacity( z, 0UL );
   checkNonZeros( z, 0UL );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c reset() member function of the ZeroVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c reset() member function of the ZeroVector class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testReset()
{
   test_ = "ZeroVector::reset()";

   // Resetting a default constructed vector
   {
      blaze::ZeroVector<int,blaze::rowVector> vec;

      reset( vec );

      checkSize    ( vec, 0UL );
      checkNonZeros( vec, 0UL );
   }

   // Resetting an initialized vector
   {
      // Initialization check
      blaze::ZeroVector<int,blaze::rowVector> z( 9UL );

      checkSize    ( z, 9UL );
      checkCapacity( z, 0UL );
      checkNonZeros( z, 0UL );

      if( z[0] != 0 || z[1] != 0 || z[2] != 0 || z[3] != 0 || z[4] != 0 ||
          z[5] != 0 || z[6] != 0 || z[7] != 0 || z[8] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << z << "\n"
             << "   Expected result:\n( 0 0 0 0 0 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting the entire vector
      reset( z );

      checkSize    ( z, 9UL );
      checkCapacity( z, 0UL );
      checkNonZeros( z, 0UL );

      if( z[0] != 0 || z[1] != 0 || z[2] != 0 || z[3] != 0 || z[4] != 0 ||
          z[5] != 0 || z[6] != 0 || z[7] != 0 || z[8] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << z << "\n"
             << "   Expected result:\n( 0 0 0 0 0 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c clear() member function of the ZeroVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c clear() member function of the ZeroVector class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testClear()
{
   test_ = "ZeroVector::clear()";

   // Clearing a default constructed vector
   {
      blaze::ZeroVector<int,blaze::rowVector> z;

      clear( z );

      checkSize    ( z, 0UL );
      checkNonZeros( z, 0UL );
   }

   // Clearing an initialized vector
   {
      // Initialization check
      blaze::ZeroVector<int,blaze::rowVector> z( 9UL );

      checkSize    ( z, 9UL );
      checkCapacity( z, 0UL );
      checkNonZeros( z, 0UL );

      if( z[0] != 0 || z[1] != 0 || z[2] != 0 || z[3] != 0 || z[4] != 0 ||
          z[5] != 0 || z[6] != 0 || z[7] != 0 || z[8] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << z << "\n"
             << "   Expected result:\n( 0 0 0 0 0 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing the entire vector
      clear( z );

      checkSize    ( z, 0UL );
      checkNonZeros( z, 0UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c resize() member function of the ZeroVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c resize() member function of the ZeroVector class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testResize()
{
   test_ = "ZeroVector::resize()";

   // Initialization check
   blaze::ZeroVector<int,blaze::rowVector> z;

   checkSize    ( z, 0UL );
   checkNonZeros( z, 0UL );

   // Resizing to 0
   z.resize( 0UL );

   checkSize    ( z, 0UL );
   checkNonZeros( z, 0UL );

   // Resizing to 5
   z.resize( 5UL );

   checkSize    ( z, 5UL );
   checkNonZeros( z, 0UL );

   // Resizing to 2
   z.resize( 2UL );

   checkSize    ( z, 2UL );
   checkNonZeros( z, 0UL );

   // Resizing to 4
   z.resize( 4UL );

   checkSize    ( z, 4UL );
   checkNonZeros( z, 0UL );

   // Resizing to 1
   z.resize( 1UL );

   checkSize    ( z, 1UL );
   checkNonZeros( z, 0UL );

   // Resizing to 0
   z.resize( 0UL );

   checkSize    ( z, 0UL );
   checkNonZeros( z, 0UL );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c swap() functionality of the ZeroVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c swap() function of the ZeroVector class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testSwap()
{
   test_ = "ZeroVector swap";

   blaze::ZeroVector<int,blaze::rowVector> z1( 9UL );
   blaze::ZeroVector<int,blaze::rowVector> z2( 5UL );

   swap( z1, z2 );

   checkSize    ( z1, 5UL );
   checkCapacity( z1, 0UL );
   checkNonZeros( z1, 0UL );

   if( z1[0] != 0 || z1[1] != 0 || z1[2] != 0 || z1[3] != 0 || z1[4] != 0 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Swapping the first vector failed\n"
          << " Details:\n"
          << "   Result:\n" << z1 << "\n"
          << "   Expected result:\n( 0 0 0 0 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   checkSize    ( z2, 9UL );
   checkCapacity( z2, 0UL );
   checkNonZeros( z2, 0UL );

   if( z2[0] != 0 || z2[1] != 0 || z2[2] != 0 || z2[3] != 0 || z2[4] != 0 ||
       z2[5] != 0 || z2[6] != 0 || z2[7] != 0 || z2[8] != 0 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Swapping the second vector failed\n"
          << " Details:\n"
          << "   Result:\n" << z2 << "\n"
          << "   Expected result:\n( 0 0 0 0 0 0 0 0 0 )\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c find() member function of the ZeroVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c find() member function of the ZeroVector class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testFind()
{
   test_ = "ZeroVector::find()";

   using ConstIterator = blaze::ZeroVector<int,blaze::rowVector>::ConstIterator;

   // Initialization check
   blaze::ZeroVector<int,blaze::rowVector> z( 8UL );

   checkSize    ( z, 8UL );
   checkCapacity( z, 0UL );
   checkNonZeros( z, 0UL );

   // Searching for the first non-existing element
   {
      ConstIterator pos( z.find( 0UL ) );

      if( pos != z.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Non-existing element could be found\n"
             << " Details:\n"
             << "   Required index = 0\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = 0\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current vector:\n" << z << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Searching for the second non-existing element
   {
      ConstIterator pos( z.find( 4UL ) );

      if( pos != z.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Non-existing element could be found\n"
             << " Details:\n"
             << "   Required index = 4\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = 0\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current vector:\n" << z << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Searching for the third non-existing element
   {
      ConstIterator pos( z.find( 7UL ) );

      if( pos != z.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Non-existing element could be found\n"
             << " Details:\n"
             << "   Required index = 7\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = 0\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current vector:\n" << z << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c lowerBound() member function of the ZeroVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c lowerBound() member function of the ZeroVector class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testLowerBound()
{
   test_ = "ZeroVector::lowerBound()";

   using ConstIterator = blaze::ZeroVector<int,blaze::rowVector>::ConstIterator;

   // Initialization check
   blaze::ZeroVector<int,blaze::rowVector> z( 8UL );

   checkSize    ( z, 8UL );
   checkCapacity( z, 0UL );
   checkNonZeros( z, 0UL );

   // Determining the lower bound for index 0
   {
      ConstIterator pos( z.lowerBound( 0UL ) );

      if( pos != z.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Lower bound could not be determined\n"
             << " Details:\n"
             << "   Required index = 0\n"
             << "   Current vector:\n" << z << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Determining the lower bound for index 1
   {
      ConstIterator pos( z.lowerBound( 1UL ) );

      if( pos != z.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Lower bound could not be determined\n"
             << " Details:\n"
             << "   Required index = 1\n"
             << "   Current vector:\n" << z << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Determining the lower bound for index 4
   {
      ConstIterator pos( z.lowerBound( 4UL ) );

      if( pos != z.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Lower bound could not be determined\n"
             << " Details:\n"
             << "   Required index = 4\n"
             << "   Current vector:\n" << z << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Determining the lower bound for index 7
   {
      ConstIterator pos( z.lowerBound( 7UL ) );

      if( pos != z.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Lower bound could not be determined\n"
             << " Details:\n"
             << "   Required index = 7\n"
             << "   Current vector:\n" << z << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c upperBound() member function of the ZeroVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c upperBound() member function of the ZeroVector class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testUpperBound()
{
   test_ = "ZeroVector::upperBound()";

   using ConstIterator = blaze::ZeroVector<int,blaze::rowVector>::ConstIterator;

   // Initialization check
   blaze::ZeroVector<int,blaze::rowVector> z( 8UL );

   checkSize    ( z, 8UL );
   checkCapacity( z, 0UL );
   checkNonZeros( z, 0UL );

   // Determining the upper bound for index 0
   {
      ConstIterator pos( z.upperBound( 0UL ) );

      if( pos != z.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Upper bound could not be determined\n"
             << " Details:\n"
             << "   Required index = 0\n"
             << "   Current vector:\n" << z << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Determining the upper bound for index 1
   {
      ConstIterator pos( z.upperBound( 1UL ) );

      if( pos != z.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Upper bound could not be determined\n"
             << " Details:\n"
             << "   Required index = 1\n"
             << "   Current vector:\n" << z << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Determining the upper bound for index 4
   {
      ConstIterator pos( z.upperBound( 4UL ) );

      if( pos != z.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Upper bound could not be determined\n"
             << " Details:\n"
             << "   Required index = 4\n"
             << "   Current vector:\n" << z << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Determining the upper bound for index 7
   {
      ConstIterator pos( z.upperBound( 7UL ) );

      if( pos != z.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Upper bound could not be determined\n"
             << " Details:\n"
             << "   Required index = 7\n"
             << "   Current vector:\n" << z << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c isDefault() function with the ZeroVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isDefault() function with the ZeroVector class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testIsDefault()
{
   test_ = "isDefault() function";

   // isDefault with vector of size 0 (default)
   {
      blaze::ZeroVector<int,blaze::rowVector> z;

      if( blaze::isDefault( z ) != true ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isDefault evaluation\n"
             << " Details:\n"
             << "   Vector:\n" << z << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // isDefault with vector of size 5 (non-default)
   {
      blaze::ZeroVector<int,blaze::rowVector> z( 5UL );

      if( blaze::isDefault( z[1] ) != true ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isDefault evaluation\n"
             << " Details:\n"
             << "   Vector element: " << z[1] << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( blaze::isDefault( z ) != false ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isDefault evaluation\n"
             << " Details:\n"
             << "   Vector:\n" << z << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************

} // namespace zerovector

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
   std::cout << "   Running ZeroVector class test..." << std::endl;

   try
   {
      RUN_ZEROVECTOR_CLASS_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during ZeroVector class test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
