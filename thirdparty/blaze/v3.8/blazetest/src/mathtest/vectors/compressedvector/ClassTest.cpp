//=================================================================================================
/*!
//  \file src/mathtest/vectors/compressedvector/ClassTest.cpp
//  \brief Source file for the CompressedVector class test
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
#include <blazetest/mathtest/vectors/compressedvector/ClassTest.h>

#ifdef BLAZE_USE_HPX_THREADS
#  include <hpx/hpx_main.hpp>
#endif


namespace blazetest {

namespace mathtest {

namespace vectors {

namespace compressedvector {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the CompressedVector class test.
//
// \exception std::runtime_error Operation error detected.
*/
ClassTest::ClassTest()
{
   testConstructors();
   testAssignment();
   testAddAssign();
   testSubAssign();
   testMultAssign();
   testDivAssign();
   testCrossAssign();
   testScaling();
   testSubscript();
   testAt();
   testIterator();
   testNonZeros();
   testReset();
   testClear();
   testResize();
   testReserve();
   testShrinkToFit();
   testSwap();
   testSet();
   testInsert();
   testAppend();
   testErase();
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
/*!\brief Test of the CompressedVector constructors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all constructors of the CompressedVector class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testConstructors()
{
   //=====================================================================================
   // Default constructor
   //=====================================================================================

   {
      test_ = "CompressedVector default constructor";

      blaze::CompressedVector<int,blaze::rowVector> vec;

      checkSize    ( vec, 0UL );
      checkNonZeros( vec, 0UL );
   }


   //=====================================================================================
   // Size constructor
   //=====================================================================================

   {
      test_ = "CompressedVector size constructor (size 0)";

      blaze::CompressedVector<int,blaze::rowVector> vec( 0UL );

      checkSize    ( vec, 0UL );
      checkNonZeros( vec, 0UL );
   }

   {
      test_ = "CompressedVector size constructor (size 5)";

      blaze::CompressedVector<int,blaze::rowVector> vec( 5UL );

      checkSize    ( vec, 5UL );
      checkNonZeros( vec, 0UL );
   }


   //=====================================================================================
   // Size/Non-zeros constructor
   //=====================================================================================

   {
      test_ = "CompressedVector size/non-zeros constructor (size 0)";

      blaze::CompressedVector<int,blaze::rowVector> vec( 0UL, 3UL );

      checkSize    ( vec, 0UL );
      checkCapacity( vec, 3UL );
      checkNonZeros( vec, 0UL );
   }

   {
      test_ = "CompressedVector size/non-zeros constructor (size 7)";

      blaze::CompressedVector<int,blaze::rowVector> vec( 7UL, 3UL );

      checkSize    ( vec, 7UL );
      checkCapacity( vec, 3UL );
      checkNonZeros( vec, 0UL );
   }


   //=====================================================================================
   // List initialization
   //=====================================================================================

   {
      test_ = "CompressedVector initializer list constructor (size 5)";

      blaze::CompressedVector<int,blaze::rowVector> vec{ 0, 2, 0, 4, 0 };

      checkSize    ( vec, 5UL );
      checkCapacity( vec, 2UL );
      checkNonZeros( vec, 2UL );

      if( vec[0] != 0 || vec[1] != 2 || vec[2] != 0 || vec[3] != 4 || vec[4] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 2 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Copy constructor
   //=====================================================================================

   {
      test_ = "CompressedVector copy constructor (size 0)";

      blaze::CompressedVector<int,blaze::rowVector> vec1( 0UL, 3UL );
      blaze::CompressedVector<int,blaze::rowVector> vec2( vec1 );

      checkSize    ( vec2, 0UL );
      checkNonZeros( vec2, 0UL );
   }

   {
      test_ = "CompressedVector copy constructor (size 7)";

      blaze::CompressedVector<int,blaze::rowVector> vec1{ 1, 2, 0, 4, 0, 0, 0 };
      blaze::CompressedVector<int,blaze::rowVector> vec2( vec1 );

      checkSize    ( vec2, 7UL );
      checkCapacity( vec2, 3UL );
      checkNonZeros( vec2, 3UL );

      if( vec2[0] != 1 || vec2[1] != 2 || vec2[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 1 2 0 4 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Move constructor
   //=====================================================================================

   {
      test_ = "CompressedVector move constructor (size 0)";

      blaze::CompressedVector<int,blaze::rowVector> vec1( 0UL, 3UL );
      blaze::CompressedVector<int,blaze::rowVector> vec2( std::move( vec1 ) );

      checkSize    ( vec2, 0UL );
      checkNonZeros( vec2, 0UL );
   }

   {
      test_ = "CompressedVector move constructor (size 7)";

      blaze::CompressedVector<int,blaze::rowVector> vec1{ 1, 2, 0, 4, 0, 0, 0 };
      blaze::CompressedVector<int,blaze::rowVector> vec2( std::move( vec1 ) );

      checkSize    ( vec2, 7UL );
      checkCapacity( vec2, 3UL );
      checkNonZeros( vec2, 3UL );

      if( vec2[0] != 1 || vec2[1] != 2 || vec2[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 1 2 0 4 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Dense vector constructor
   //=====================================================================================

   {
      test_ = "CompressedVector dense vector constructor";

      blaze::DynamicVector<int,blaze::rowVector> vec1{ 10, 11, 12, 0, 13 };
      blaze::CompressedVector<int,blaze::rowVector> vec2( vec1 );

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 4UL );
      checkNonZeros( vec2, 4UL );

      if( vec2[0] != 10 || vec2[1] != 11 || vec2[2] != 12 || vec2[3] != 0 || vec2[4] != 13 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 10 11 12 0 13 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Sparse vector constructor
   //=====================================================================================

   {
      test_ = "CompressedVector sparse vector assignment";

      blaze::CompressedVector<int,blaze::columnVector> vec1{ 1, 2, 0, 4, 0, 0, 0 };
      blaze::CompressedVector<int,blaze::rowVector> vec2( trans( vec1 ) );

      checkSize    ( vec2, 7UL );
      checkNonZeros( vec2, 3UL );

      if( vec2[0] != 1 || vec2[1] != 2 || vec2[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 1 2 0 4 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the CompressedVector assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all assignment operators of the CompressedVector class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testAssignment()
{
   //=====================================================================================
   // List assignment
   //=====================================================================================

   {
      test_ = "CompressedVector initializer list assignment";

      blaze::CompressedVector<int,blaze::rowVector> vec;
      vec = { 0, 2, 0, 4, 0 };

      checkSize    ( vec, 5UL );
      checkCapacity( vec, 2UL );
      checkNonZeros( vec, 2UL );

      if( vec[0] != 0 || vec[1] != 2 || vec[2] != 0 || vec[3] != 4 || vec[4] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 2 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Copy assignment
   //=====================================================================================

   {
      test_ = "CompressedVector copy assignment";

      blaze::CompressedVector<int,blaze::rowVector> vec1{ 1, 2, 0, 4, 0, 0, 0 };
      blaze::CompressedVector<int,blaze::rowVector> vec2;
      vec2 = vec1;

      checkSize    ( vec2, 7UL );
      checkCapacity( vec2, 3UL );
      checkNonZeros( vec2, 3UL );

      if( vec2[0] != 1 || vec2[1] != 2 || vec2[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 1 2 0 4 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "CompressedVector copy assignment stress test";

      using RandomVectorType = blaze::CompressedVector<int,blaze::rowVector>;

      blaze::CompressedVector<int,blaze::rowVector> vec1;
      const int min( randmin );
      const int max( randmax );

      for( size_t i=0UL; i<100UL; ++i )
      {
         const size_t size( blaze::rand<size_t>( 0UL, 20UL ) );
         const RandomVectorType vec2( blaze::rand<RandomVectorType>( size, min, max ) );

         vec1 = vec2;

         if( vec1 != vec2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << vec1 << "\n"
                << "   Expected result:\n" << vec2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Move assignment
   //=====================================================================================

   {
      test_ = "CompressedVector move assignment";

      blaze::CompressedVector<int,blaze::rowVector> vec1{ 1, 2, 0, 4, 0, 0, 0 };
      blaze::CompressedVector<int,blaze::rowVector> vec2{ 0, 0, 11, 0 };

      vec2 = std::move( vec1 );

      checkSize    ( vec2, 7UL );
      checkCapacity( vec2, 3UL );
      checkNonZeros( vec2, 3UL );

      if( vec2[0] != 1 || vec2[1] != 2 || vec2[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 1 2 0 4 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Dense vector assignment
   //=====================================================================================

   {
      test_ = "CompressedVector dense vector assignment";

      blaze::DynamicVector<int,blaze::rowVector> vec1{ 10, 11, 12, 0, 13 };
      blaze::CompressedVector<int,blaze::rowVector> vec2;
      vec2 = vec1;

      checkSize    ( vec2, 5UL );
      checkNonZeros( vec2, 4UL );

      if( vec2[0] != 10 || vec2[1] != 11 || vec2[2] != 12 || vec2[3] != 0 || vec2[4] != 13 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 10 11 12 0 13 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "CompressedVector dense vector assignment stress test";

      using RandomVectorType = blaze::DynamicVector<int,blaze::rowVector>;

      blaze::CompressedVector<int,blaze::rowVector> vec1;
      const int min( randmin );
      const int max( randmax );

      for( size_t i=0UL; i<100UL; ++i )
      {
         const size_t size( blaze::rand<size_t>( 0UL, 20UL ) );
         const RandomVectorType vec2( blaze::rand<RandomVectorType>( size, min, max ) );

         vec1 = vec2;

         if( vec1 != vec2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << vec1 << "\n"
                << "   Expected result:\n" << vec2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Sparse vector assignment
   //=====================================================================================

   {
      test_ = "CompressedVector sparse vector assignment";

      blaze::CompressedVector<int,blaze::columnVector> vec1{ 1, 2, 0, 4, 0, 0, 0 };
      blaze::CompressedVector<int,blaze::rowVector> vec2;
      vec2 = trans( vec1 );

      checkSize    ( vec2, 7UL );
      checkNonZeros( vec2, 3UL );

      if( vec2[0] != 1 || vec2[1] != 2 || vec2[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 1 2 0 4 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "CompressedVector sparse vector assignment stress test";

      using RandomVectorType = blaze::CompressedVector<short,blaze::rowVector>;

      blaze::CompressedVector<int,blaze::rowVector> vec1;
      const short min( randmin );
      const short max( randmax );

      for( size_t i=0UL; i<100UL; ++i )
      {
         const size_t size( blaze::rand<size_t>( 0UL, 20UL ) );
         const RandomVectorType vec2( blaze::rand<RandomVectorType>( size, min, max ) );

         vec1 = vec2;

         if( vec1 != vec2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << vec1 << "\n"
                << "   Expected result:\n" << vec2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the CompressedVector addition assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the addition assignment operators of the CompressedVector
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testAddAssign()
{
   //=====================================================================================
   // Dense vector addition assignment
   //=====================================================================================

   {
      test_ = "CompressedVector dense vector addition assignment";

      blaze::DynamicVector<int,blaze::rowVector> vec1{ 10, 11, 12, 0, 13 };
      blaze::CompressedVector<int,blaze::rowVector> vec2{ 1, 2, 0, 4, 0 };

      vec2 += vec1;

      checkSize    ( vec2, 5UL );
      checkNonZeros( vec2, 5UL );

      if( vec2[0] != 11 || vec2[1] != 13 || vec2[2] != 12 || vec2[3] != 4 || vec2[4] != 13 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 11 13 12 4 13 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Sparse vector addition assignment
   //=====================================================================================

   {
      test_ = "CompressedVector sparse vector addition assignment";

      blaze::CompressedVector<int,blaze::columnVector> vec1{ 1, 2, 0, 4, 0 };
      blaze::CompressedVector<int,blaze::rowVector> vec2{ 0, 5, 6, 0, 0 };

      vec2 += trans( vec1 );

      checkSize    ( vec2, 5UL );
      checkNonZeros( vec2, 4UL );

      if( vec2[0] != 1 || vec2[1] != 7 || vec2[2] != 6 || vec2[3] != 4 || vec2[4] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 1 7 6 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the CompressedVector subtraction assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the subtraction assignment operators of the CompressedVector
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testSubAssign()
{
   //=====================================================================================
   // Dense vector subtraction assignment
   //=====================================================================================

   {
      test_ = "CompressedVector dense vector subtraction assignment";

      blaze::DynamicVector<int,blaze::rowVector> vec1{ 10, 11, 12, 0, 13 };
      blaze::CompressedVector<int,blaze::rowVector> vec2{ 1, 2, 0, 4, 0 };

      vec2 -= vec1;

      checkSize    ( vec2, 5UL );
      checkNonZeros( vec2, 5UL );

      if( vec2[0] != -9 || vec2[1] != -9 || vec2[2] != -12 || vec2[3] != 4 || vec2[4] != -13 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( -9 -9 -12 4 -13 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Sparse vector subtraction assignment
   //=====================================================================================

   {
      test_ = "CompressedVector sparse vector subtraction assignment";

      blaze::CompressedVector<int,blaze::columnVector> vec1{ 1, 2, 0, 4, 0 };
      blaze::CompressedVector<int,blaze::rowVector> vec2{ 0, 5, 6, 0, 0 };

      vec2 -= trans( vec1 );

      checkSize    ( vec2, 5UL );
      checkNonZeros( vec2, 4UL );

      if( vec2[0] != -1 || vec2[1] != 3 || vec2[2] != 6 || vec2[3] != -4 || vec2[4] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( -1 3 6 -4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the CompressedVector multiplication assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the multiplication assignment operators of the CompressedVector
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testMultAssign()
{
   //=====================================================================================
   // Dense vector multiplication assignment
   //=====================================================================================

   {
      test_ = "CompressedVector dense vector multiplication assignment";

      blaze::DynamicVector<int,blaze::rowVector> vec1{ 10, 11, 12, 0, 13 };
      blaze::CompressedVector<int,blaze::rowVector> vec2{ 1, 2, 0, 4, 0 };

      vec2 *= vec1;

      checkSize    ( vec2, 5UL );
      checkNonZeros( vec2, 3UL );

      if( vec2[0] != 10 || vec2[1] != 22 || vec2[2] != 0 || vec2[3] != 0 || vec2[4] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 10 22 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Sparse vector multiplication assignment
   //=====================================================================================

   {
      test_ = "CompressedVector sparse vector multiplication assignment";

      blaze::CompressedVector<int,blaze::columnVector> vec1{ 1, 2, 0, 4, 0 };
      blaze::CompressedVector<int,blaze::rowVector> vec2{ 0, 5, 6, 0, 0 };

      vec2 *= trans( vec1 );

      checkSize    ( vec2, 5UL );
      checkNonZeros( vec2, 1UL );

      if( vec2[0] != 0 || vec2[1] != 10 || vec2[2] != 0 || vec2[3] != 0 || vec2[4] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 0 10 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the CompressedVector division assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the division assignment operators of the CompressedVector
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testDivAssign()
{
   //=====================================================================================
   // Dense vector division assignment
   //=====================================================================================

   {
      test_ = "CompressedVector dense vector division assignment";

      blaze::DynamicVector<int,blaze::rowVector> vec1{ 1, 2, -3, 4, 1 };
      blaze::CompressedVector<int,blaze::rowVector> vec2{ 2, 0, -3, 8, 0 };

      vec2 /= vec1;

      checkSize    ( vec2, 5UL );
      checkNonZeros( vec2, 3UL );

      if( vec2[0] != 2 || vec2[1] != 0 || vec2[2] != 1 || vec2[3] != 2 || vec2[4] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 2 0 1 2 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the CompressedVector cross product assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the cross product assignment operators of the CompressedVector
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testCrossAssign()
{
   //=====================================================================================
   // Dense vector cross product assignment
   //=====================================================================================

   {
      test_ = "CompressedVector dense vector cross product assignment";

      blaze::DynamicVector<int,blaze::rowVector> vec1{ 1, 0, -2 };
      blaze::CompressedVector<int,blaze::rowVector> vec2{ 2, 0, -1 };

      vec2 %= vec1;

      checkSize    ( vec2, 3UL );
      checkNonZeros( vec2, 1UL );

      if( vec2[0] != 0 || vec2[1] != 3 || vec2[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cross product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 0 3 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Sparse vector cross product assignment
   //=====================================================================================

   {
      test_ = "CompressedVector sparse vector cross product assignment";

      blaze::CompressedVector<int,blaze::columnVector> vec1{ 1, 0, -2 };
      blaze::CompressedVector<int,blaze::rowVector> vec2{ 2, 0, -1 };

      vec2 %= trans( vec1 );

      checkSize    ( vec2, 3UL );
      checkNonZeros( vec2, 1UL );

      if( vec2[0] != 0 || vec2[1] != 3 || vec2[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cross product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 0 3 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of all CompressedVector (self-)scaling operations.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all available ways to scale an instance of the CompressedVector
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testScaling()
{
   //=====================================================================================
   // Self-scaling (v*=s)
   //=====================================================================================

   {
      test_ = "CompressedVector self-scaling (v*=s)";

      blaze::CompressedVector<int,blaze::columnVector> vec1{ 1, 2, 0, 4, 0 };

      vec1 *= 2;

      checkSize    ( vec1, 5UL );
      checkNonZeros( vec1, 3UL );

      if( vec1[0] != 2 || vec1[1] != 4 || vec1[2] != 0 || vec1[3] != 8 || vec1[4] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << vec1 << "\n"
             << "   Expected result:\n( 2 4 0 8 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Self-scaling (v=v*s)
   //=====================================================================================

   {
      test_ = "CompressedVector self-scaling (v=v*s)";

      blaze::CompressedVector<int,blaze::columnVector> vec1{ 1, 2, 0, 4, 0 };

      vec1 = vec1 * 2;

      checkSize    ( vec1, 5UL );
      checkNonZeros( vec1, 3UL );

      if( vec1[0] != 2 || vec1[1] != 4 || vec1[2] != 0 || vec1[3] != 8 || vec1[4] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << vec1 << "\n"
             << "   Expected result:\n( 2 4 0 8 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Self-scaling (v=s*v)
   //=====================================================================================

   {
      test_ = "CompressedVector self-scaling (v=v*s)";

      blaze::CompressedVector<int,blaze::columnVector> vec1{ 1, 2, 0, 4, 0 };

      vec1 = 2 * vec1;

      checkSize    ( vec1, 5UL );
      checkNonZeros( vec1, 3UL );

      if( vec1[0] != 2 || vec1[1] != 4 || vec1[2] != 0 || vec1[3] != 8 || vec1[4] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << vec1 << "\n"
             << "   Expected result:\n( 2 4 0 8 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Self-scaling (v/=s)
   //=====================================================================================

   {
      test_ = "CompressedVector self-scaling (v/=s)";

      blaze::CompressedVector<int,blaze::columnVector> vec1{ 2, 4, 0, 8, 0 };

      vec1 /= 2;

      checkSize    ( vec1, 5UL );
      checkNonZeros( vec1, 3UL );

      if( vec1[0] != 1 || vec1[1] != 2 || vec1[2] != 0 || vec1[3] != 4 || vec1[4] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << vec1 << "\n"
             << "   Expected result:\n( 1 2 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Self-scaling (v=v/s)
   //=====================================================================================

   {
      test_ = "CompressedVector self-scaling (v=v/s)";

      blaze::CompressedVector<int,blaze::columnVector> vec1{ 2, 4, 0, 8, 0 };

      vec1 = vec1 / 2;

      checkSize    ( vec1, 5UL );
      checkNonZeros( vec1, 3UL );

      if( vec1[0] != 1 || vec1[1] != 2 || vec1[2] != 0 || vec1[3] != 4 || vec1[4] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << vec1 << "\n"
             << "   Expected result:\n( 1 2 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // CompressedVector::scale()
   //=====================================================================================

   {
      test_ = "CompressedVector::scale() (int)";

      // Initialization check
      blaze::CompressedVector<int,blaze::rowVector> vec{ 0, 1, 0, 2, 0, 3 };

      checkSize    ( vec, 6UL );
      checkCapacity( vec, 3UL );
      checkNonZeros( vec, 3UL );

      if( vec[1] != 1 || vec[3] != 2 || vec[5] != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 1 0 2 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Integral scaling of the vector
      vec.scale( 2 );

      checkSize    ( vec, 6UL );
      checkCapacity( vec, 3UL );
      checkNonZeros( vec, 3UL );

      if( vec[1] != 2 || vec[3] != 4 || vec[5] != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 2 0 4 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Floating point scaling of the vector
      vec.scale( 0.5 );

      checkSize    ( vec, 6UL );
      checkCapacity( vec, 3UL );
      checkNonZeros( vec, 3UL );

      if( vec[1] != 1 || vec[3] != 2 || vec[5] != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 1 0 2 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "CompressedVector::scale() (complex)";

      using blaze::complex;

      blaze::CompressedVector<complex<float>,blaze::rowVector> vec( 2UL, 2UL );
      vec[0] = complex<float>( 1.0F, 0.0F );
      vec[1] = complex<float>( 2.0F, 0.0F );
      vec.scale( complex<float>( 3.0F, 0.0F ) );

      checkSize    ( vec, 2UL );
      checkCapacity( vec, 2UL );
      checkNonZeros( vec, 2UL );

      if( vec[0] != complex<float>( 3.0F, 0.0F ) || vec[1] != complex<float>( 6.0F, 0.0F ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( (3,0) (6,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the CompressedVector subscript operator.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of adding and accessing elements via the subscript operator
// of the CompressedVector class template. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
void ClassTest::testSubscript()
{
   test_ = "CompressedVector::operator[]";

   // Assignment to the element at index 2
   blaze::CompressedVector<int,blaze::rowVector> vec( 7UL, 3UL );
   vec[2] = 1;

   checkSize    ( vec, 7UL );
   checkCapacity( vec, 3UL );
   checkNonZeros( vec, 1UL );

   if( vec[2] != 1 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 0 0 1 0 0 0 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Assignment to the element at index 5
   vec[5] = 2;

   checkSize    ( vec, 7UL );
   checkCapacity( vec, 3UL );
   checkNonZeros( vec, 2UL );

   if( vec[2] != 1 || vec[5] != 2 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 0 0 1 0 0 2 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Assignment to the element at index 3
   vec[3] = 3;

   checkSize    ( vec, 7UL );
   checkCapacity( vec, 3UL );
   checkNonZeros( vec, 3UL );

   if( vec[2] != 1 || vec[3] != 3 || vec[5] != 2 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 0 0 1 3 0 2 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Assignment to the element at index 0
   vec[0] = 4;

   checkSize    ( vec, 7UL );
   checkCapacity( vec, 4UL );
   checkNonZeros( vec, 4UL );

   if( vec[0] != 4 || vec[2] != 1 || vec[3] != 3 || vec[5] != 2 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 4 0 1 3 0 2 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Addition assignment to the element at index 2
   vec[2] += vec[3];

   checkSize    ( vec, 7UL );
   checkCapacity( vec, 4UL );
   checkNonZeros( vec, 4UL );

   if( vec[0] != 4 || vec[2] != 4 || vec[3] != 3 || vec[5] != 2 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 4 0 4 3 0 2 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Subtraction assignment to the element at index 1
   vec[1] -= vec[5];

   checkSize    ( vec, 7UL );
   checkCapacity( vec, 5UL );
   checkNonZeros( vec, 5UL );

   if( vec[0] != 4 || vec[1] != -2 || vec[2] != 4 || vec[3] != 3 || vec[5] != 2 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 4 -2 4 3 0 2 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Multiplication assignment to the element at index 3
   vec[3] *= -3;

   checkSize    ( vec, 7UL );
   checkCapacity( vec, 5UL );
   checkNonZeros( vec, 5UL );

   if( vec[0] != 4 || vec[1] != -2 || vec[2] != 4 || vec[3] != -9 || vec[5] != 2 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 4 -2 4 -9 0 2 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Division assignment to the element at index 2
   vec[2] /= 2;

   checkSize    ( vec, 7UL );
   checkCapacity( vec, 5UL );
   checkNonZeros( vec, 5UL );

   if( vec[0] != 4 || vec[1] != -2 || vec[2] != 2 || vec[3] != -9 || vec[5] != 2 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 4 -2 2 -9 0 2 0 )\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c at() member function of the CompressedVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of adding and accessing elements via the \c at() member function
// of the CompressedVector class template. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
void ClassTest::testAt()
{
   test_ = "CompressedVector::at()";

   // Assignment to the element at index 2
   blaze::CompressedVector<int,blaze::rowVector> vec( 7UL, 3UL );
   vec.at(2) = 1;

   checkSize    ( vec, 7UL );
   checkCapacity( vec, 3UL );
   checkNonZeros( vec, 1UL );

   if( vec.at(2) != 1 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Access via at() function failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 0 0 1 0 0 0 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Assignment to the element at index 5
   vec.at(5) = 2;

   checkSize    ( vec, 7UL );
   checkCapacity( vec, 3UL );
   checkNonZeros( vec, 2UL );

   if( vec.at(2) != 1 || vec.at(5) != 2 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Access via at() function failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 0 0 1 0 0 2 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Assignment to the element at index 3
   vec.at(3) = 3;

   checkSize    ( vec, 7UL );
   checkCapacity( vec, 3UL );
   checkNonZeros( vec, 3UL );

   if( vec.at(2) != 1 || vec.at(3) != 3 || vec.at(5) != 2 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Access via at() function failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 0 0 1 3 0 2 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Assignment to the element at index 0
   vec.at(0) = 4;

   checkSize    ( vec, 7UL );
   checkCapacity( vec, 4UL );
   checkNonZeros( vec, 4UL );

   if( vec.at(0) != 4 || vec.at(2) != 1 || vec.at(3) != 3 || vec.at(5) != 2 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Access via at() function failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 4 0 1 3 0 2 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Addition assignment to the element at index 2
   vec.at(2) += vec.at(3);

   checkSize    ( vec, 7UL );
   checkCapacity( vec, 4UL );
   checkNonZeros( vec, 4UL );

   if( vec.at(0) != 4 || vec.at(2) != 4 || vec.at(3) != 3 || vec.at(5) != 2 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Access via at() function failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 4 0 4 3 0 2 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Subtraction assignment to the element at index 1
   vec.at(1) -= vec.at(5);

   checkSize    ( vec, 7UL );
   checkCapacity( vec, 5UL );
   checkNonZeros( vec, 5UL );

   if( vec.at(0) != 4 || vec.at(1) != -2 || vec.at(2) != 4 || vec.at(3) != 3 || vec.at(5) != 2 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Access via at() function failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 4 -2 4 3 0 2 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Multiplication assignment to the element at index 3
   vec.at(3) *= -3;

   checkSize    ( vec, 7UL );
   checkCapacity( vec, 5UL );
   checkNonZeros( vec, 5UL );

   if( vec.at(0) != 4 || vec.at(1) != -2 || vec.at(2) != 4 || vec.at(3) != -9 || vec.at(5) != 2 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Access via at() function failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 4 -2 4 -9 0 2 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Division assignment to the element at index 2
   vec.at(2) /= 2;

   checkSize    ( vec, 7UL );
   checkCapacity( vec, 5UL );
   checkNonZeros( vec, 5UL );

   if( vec.at(0) != 4 || vec.at(1) != -2 || vec.at(2) != 2 || vec.at(3) != -9 || vec.at(5) != 2 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Access via at() function failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 4 -2 2 -9 0 2 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Attempt to assign to the element at index 7
   try {
      vec.at(7) = 2;

      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Out-of-bound access succeeded\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 4 -2 2 -9 0 2 0 )\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::out_of_range& ) {}
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the CompressedVector iterator implementation.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the iterator implementation of the CompressedVector class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testIterator()
{
   using VectorType    = blaze::CompressedVector<int>;
   using Iterator      = VectorType::Iterator;
   using ConstIterator = VectorType::ConstIterator;

   VectorType vec{ 0, -2, -3, 0 };

   // Testing the Iterator default constructor
   {
      test_ = "Iterator default constructor";

      Iterator it{};

      if( it != Iterator() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed iterator default constructor\n";
         throw std::runtime_error( oss.str() );
      }
   }

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

   // Testing conversion from Iterator to ConstIterator
   {
      test_ = "Iterator/ConstIterator conversion";

      ConstIterator it( begin( vec ) );

      if( it == end( vec ) || it->value() != -2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed iterator conversion detected\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Counting the number of elements via Iterator (end-begin)
   {
      test_ = "Iterator subtraction (end-begin)";

      const ptrdiff_t number( end( vec ) - begin( vec ) );

      if( number != 2L ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid number of elements detected\n"
             << " Details:\n"
             << "   Number of elements         : " << number << "\n"
             << "   Expected number of elements: 2\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Counting the number of elements via ConstIterator (end-begin)
   {
      test_ = "ConstIterator subtraction (end-begin)";

      const ptrdiff_t number( cend( vec ) - cbegin( vec ) );

      if( number != 2L ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid number of elements detected\n"
             << " Details:\n"
             << "   Number of elements         : " << number << "\n"
             << "   Expected number of elements: 2\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing read-only access via ConstIterator
   {
      test_ = "Read-only access via ConstIterator";

      ConstIterator it ( cbegin( vec ) );
      ConstIterator end( cend( vec ) );

      if( it == end || it->value() != -2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid initial iterator detected\n";
         throw std::runtime_error( oss.str() );
      }

      ++it;

      if( it == end || it->value() != -3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator pre-increment failed\n";
         throw std::runtime_error( oss.str() );
      }

      it++;

      if( it != end ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator post-increment failed\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing assignment via Iterator
   {
      test_ = "Assignment via Iterator";

      int value = 6;

      for( Iterator it=begin( vec ); it!=end( vec ); ++it ) {
         *it = value++;
      }

      if( vec[0] != 0 || vec[1] != 6 || vec[2] != 7 || vec[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 6 7 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing addition assignment via Iterator
   {
      test_ = "Addition assignment via Iterator";

      int value = 2;

      for( Iterator it=begin( vec ); it!=end( vec ); ++it ) {
         *it += value++;
      }

      if( vec[0] != 0 || vec[1] != 8 || vec[2] != 10 || vec[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 8 10 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing subtraction assignment via Iterator
   {
      test_ = "Subtraction assignment via Iterator";

      int value = 2;

      for( Iterator it=begin( vec ); it!=end( vec ); ++it ) {
         *it -= value++;
      }

      if( vec[0] != 0 || vec[1] != 6 || vec[2] != 7 || vec[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 6 7 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing multiplication assignment via Iterator
   {
      test_ = "Multiplication assignment via Iterator";

      int value = 1;

      for( Iterator it=begin( vec ); it!=end( vec ); ++it ) {
         *it *= value++;
      }

      if( vec[0] != 0 || vec[1] != 6 || vec[2] != 14 || vec[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 6 14 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing division assignment via Iterator
   {
      test_ = "Division assignment via Iterator";

      for( Iterator it=begin( vec ); it!=end( vec ); ++it ) {
         *it /= 2;
      }

      if( vec[0] != 0 || vec[1] != 3 || vec[2] != 7 || vec[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 3 7 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c nonZeros() member function of the CompressedVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c nonZeros() member function of the CompressedVector
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testNonZeros()
{
   test_ = "CompressedVector::nonZeros()";

   // Initialization check
   blaze::CompressedVector<int,blaze::rowVector> vec( 7UL, 3UL );

   checkSize    ( vec, 7UL );
   checkCapacity( vec, 3UL );
   checkNonZeros( vec, 0UL );

   // Adding two non-zero elements
   vec[2] = 1;
   vec[5] = 2;

   checkSize    ( vec, 7UL );
   checkCapacity( vec, 3UL );
   checkNonZeros( vec, 2UL );

   // Adding a third non-zero elements
   vec[3] = 0;

   checkSize    ( vec, 7UL );
   checkCapacity( vec, 3UL );
   checkNonZeros( vec, 2UL );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c reset() member function of the CompressedVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c reset() member function of the CompressedVector
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testReset()
{
   test_ = "CompressedVector::reset()";

   // Resetting a default constructed vector
   {
      blaze::CompressedVector<int,blaze::rowVector> vec;

      reset( vec );

      checkSize    ( vec, 0UL );
      checkNonZeros( vec, 0UL );
   }

   // Resetting an initialized vector
   {
      // Initialization check
      blaze::CompressedVector<int,blaze::rowVector> vec{ 0, 1, 0, 2, 0, 0, 0, 3, 0, 4, 0 };

      checkSize    ( vec, 11UL );
      checkCapacity( vec,  4UL );
      checkNonZeros( vec,  4UL );

      if( vec[1] != 1 || vec[3] != 2 || vec[7] != 3 || vec[9] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 1 0 2 0 0 0 3 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting a single element
      reset( vec[7] );

      checkSize    ( vec, 11UL );
      checkCapacity( vec,  4UL );
      checkNonZeros( vec,  3UL );

      if( vec[1] != 1 || vec[3] != 2 || vec[7] != 0 || vec[9] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 1 0 2 0 0 0 0 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting the vector
      reset( vec );

      checkSize    ( vec, 11UL );
      checkCapacity( vec,  4UL );
      checkNonZeros( vec,  0UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c clear() member function of the CompressedVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c clear() member function of the CompressedVector
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testClear()
{
   test_ = "CompressedVector::clear()";

   // Clearing a default constructed vector
   {
      blaze::CompressedVector<int,blaze::rowVector> vec;

      clear( vec );

      checkSize    ( vec, 0UL );
      checkNonZeros( vec, 0UL );
   }

   // Clearing an initialized vector
   {
      // Initialization check
      blaze::CompressedVector<int,blaze::rowVector> vec{ 1, 0, 0, 0, 0, 0, 0, 2, 3 };

      checkSize    ( vec, 9UL );
      checkCapacity( vec, 3UL );
      checkNonZeros( vec, 3UL );

      if( vec[0] != 1 || vec[7] != 2 || vec[8] != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 0 0 0 0 0 0 2 3 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing a single element
      clear( vec[7] );

      checkSize    ( vec, 9UL );
      checkCapacity( vec, 3UL );
      checkNonZeros( vec, 2UL );

      if( vec[0] != 1 || vec[7] != 0 || vec[8] != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 0 0 0 0 0 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing the vector
      clear( vec );

      checkSize    ( vec, 0UL );
      checkNonZeros( vec, 0UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c resize() member function of the CompressedVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c resize() member function of the CompressedVector
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testResize()
{
   test_ = "CompressedVector::resize()";

   // Initialization check
   blaze::CompressedVector<int,blaze::rowVector> vec;

   checkSize    ( vec, 0UL );
   checkNonZeros( vec, 0UL );

   // Resizing to 0
   vec.resize( 0UL );

   checkSize    ( vec, 0UL );
   checkNonZeros( vec, 0UL );

   // Resizing to 5
   vec.resize( 5UL );

   checkSize    ( vec, 5UL );
   checkNonZeros( vec, 0UL );

   // Resizing to 9 and preserving the elements
   vec[0] = 1;
   vec[2] = 2;
   vec[4] = 3;
   vec.resize( 9UL, true );

   checkSize    ( vec, 9UL );
   checkCapacity( vec, 3UL );
   checkNonZeros( vec, 3UL );

   if( vec[0] != 1 || vec[2] != 2 || vec[4] != 3 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Resizing the vector failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 1 0 2 0 3 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Resizing to 2 and preserving the elements
   vec.resize( 2UL, true );

   checkSize    ( vec, 2UL );
   checkCapacity( vec, 1UL );
   checkNonZeros( vec, 1UL );

   if( vec[0] != 1 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Resizing the vector failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 1 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Resizing to 1
   vec.resize( 1UL );

   checkSize( vec, 1UL );

   // Resizing to 0
   vec.resize( 0UL );

   checkSize    ( vec, 0UL );
   checkNonZeros( vec, 0UL );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c reserve() member function of the CompressedVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c reserve() member function of the CompressedVector
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testReserve()
{
   test_ = "CompressedVector::reserve()";

   // Initialization check
   blaze::CompressedVector<int,blaze::rowVector> vec;

   checkSize    ( vec, 0UL );
   checkNonZeros( vec, 0UL );

   // Increasing the capacity of the vector
   vec.reserve( 10UL );

   checkSize    ( vec,  0UL );
   checkCapacity( vec, 10UL );
   checkNonZeros( vec,  0UL );

   // Further increasing the capacity of the vector
   vec.reserve( 20UL );

   checkSize    ( vec,  0UL );
   checkCapacity( vec, 20UL );
   checkNonZeros( vec,  0UL );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c shrinkToFit() member function of the CompressedVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c shrinkToFit() member function of the CompressedVector
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testShrinkToFit()
{
   test_ = "CompressedVector::shrinkToFit()";

   // Shrinking a vector without excessive capacity
   {
      blaze::CompressedVector<int,blaze::rowVector> vec( 5UL, 3UL );
      vec[0] = 1;
      vec[2] = 3;
      vec[4] = 5;

      vec.shrinkToFit();

      checkSize    ( vec, 5UL );
      checkCapacity( vec, 3UL );
      checkNonZeros( vec, 3UL );

      if( vec.capacity() != vec.nonZeros() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Shrinking the vector failed\n"
             << " Details:\n"
             << "   Capacity         : " << vec.capacity() << "\n"
             << "   Expected capacity: " << vec.nonZeros() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec[0] != 1 || vec[2] != 3 || vec[4] != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Shrinking the vector failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 0 3 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Shrinking a vector with excessive capacity
   {
      blaze::CompressedVector<int,blaze::rowVector> vec( 5UL, 100UL );
      vec[0] = 1;
      vec[2] = 3;
      vec[4] = 5;

      vec.shrinkToFit();

      checkSize    ( vec, 5UL );
      checkCapacity( vec, 3UL );
      checkNonZeros( vec, 3UL );

      if( vec.capacity() != vec.nonZeros() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Shrinking the vector failed\n"
             << " Details:\n"
             << "   Capacity         : " << vec.capacity() << "\n"
             << "   Expected capacity: " << vec.nonZeros() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec[0] != 1 || vec[2] != 3 || vec[4] != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Shrinking the vector failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 0 3 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c swap() functionality of the CompressedVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c swap() function of the CompressedVector class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testSwap()
{
   test_ = "CompressedVector swap";

   blaze::CompressedVector<int,blaze::rowVector> vec1{ 0, 1, 0, 0, 2, 0, 0, 3, 0, 0, 4, 0 };
   blaze::CompressedVector<int,blaze::rowVector> vec2{ 4, 0, 0, 0, 2 };

   swap( vec1, vec2 );

   checkSize    ( vec1, 5UL );
   checkCapacity( vec1, 2UL );
   checkNonZeros( vec1, 2UL );

   if( vec1[0] != 4 || vec1[4] != 2 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Swapping the first vector failed\n"
          << " Details:\n"
          << "   Result:\n" << vec1 << "\n"
          << "   Expected result:\n( 4 0 0 0 2 )\n";
      throw std::runtime_error( oss.str() );
   }

   checkSize    ( vec2, 12UL );
   checkCapacity( vec2,  4UL );
   checkNonZeros( vec2,  4UL );

   if( vec2[1] != 1 || vec2[4] != 2 || vec2[7] != 3 || vec2[10] != 4 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Swapping the second vector failed\n"
          << " Details:\n"
          << "   Result:\n" << vec1 << "\n"
          << "   Expected result:\n( 0 1 0 0 2 0 0 3 0 0 4 0 )\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c set() member function of the CompressedVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c set() member function of the CompressedVector
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testSet()
{
   test_ = "CompressedVector::set()";

   using Iterator = blaze::CompressedVector<int,blaze::rowVector>::Iterator;

   // Initialization check
   blaze::CompressedVector<int,blaze::rowVector> vec( 7UL );

   checkSize    ( vec, 7UL );
   checkNonZeros( vec, 0UL );

   // Setting a non-zero element
   {
      Iterator pos = vec.set( 4UL, 1 );

      checkSize    ( vec, 7UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      if( pos->value() != 1 || pos->index() != 4UL ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid iterator returned\n"
             << " Details:\n"
             << "   Value: " << pos->value() << "\n"
             << "   Index: " << pos->index() << "\n"
             << "   Expected value: 1\n"
             << "   Expected index: 4\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec[4] != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setting an element failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 0 0 0 1 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Setting a second non-zero element
   {
      Iterator pos = vec.set( 6UL, 2 );

      checkSize    ( vec, 7UL );
      checkCapacity( vec, 2UL );
      checkNonZeros( vec, 2UL );

      if( pos->value() != 2 || pos->index() != 6UL ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid iterator returned\n"
             << " Details:\n"
             << "   Value: " << pos->value() << "\n"
             << "   Index: " << pos->index() << "\n"
             << "   Expected value: 2\n"
             << "   Expected index: 6\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec[4] != 1 || vec[6] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setting an element failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 0 0 0 1 0 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Setting a third non-zero element
   {
      Iterator pos = vec.set( 2UL, 3 );

      checkSize    ( vec, 7UL );
      checkCapacity( vec, 3UL );
      checkNonZeros( vec, 3UL );

      if( pos->value() != 3 || pos->index() != 2UL ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid iterator returned\n"
             << " Details:\n"
             << "   Value: " << pos->value() << "\n"
             << "   Index: " << pos->index() << "\n"
             << "   Expected value: 3\n"
             << "   Expected index: 2\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec[2] != 3 || vec[4] != 1 || vec[6] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setting an element failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 0 3 0 1 0 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Setting a fourth non-zero element
   {
      Iterator pos = vec.set( 3UL, 4 );

      checkSize    ( vec, 7UL );
      checkCapacity( vec, 4UL );
      checkNonZeros( vec, 4UL );

      if( pos->value() != 4 || pos->index() != 3UL ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid iterator returned\n"
             << " Details:\n"
             << "   Value: " << pos->value() << "\n"
             << "   Index: " << pos->index() << "\n"
             << "   Expected value: 4\n"
             << "   Expected index: 3\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec[2] != 3 || vec[3] != 4 || vec[4] != 1 || vec[6] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setting an element failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 0 3 4 1 0 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Setting an already existing element
   {
      Iterator pos = vec.set( 3UL, 5 );

      checkSize    ( vec, 7UL );
      checkCapacity( vec, 4UL );
      checkNonZeros( vec, 4UL );

      if( pos->value() != 5 || pos->index() != 3UL ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid iterator returned\n"
             << " Details:\n"
             << "   Value: " << pos->value() << "\n"
             << "   Index: " << pos->index() << "\n"
             << "   Expected value: 5\n"
             << "   Expected index: 3\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec[2] != 3 || vec[3] != 5 || vec[4] != 1 || vec[6] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setting an element failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 0 3 5 1 0 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c insert() member function of the CompressedVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c insert() member function of the CompressedVector
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testInsert()
{
   test_ = "CompressedVector::insert()";

   using Iterator = blaze::CompressedVector<int,blaze::rowVector>::Iterator;

   // Initialization check
   blaze::CompressedVector<int,blaze::rowVector> vec( 7UL );

   checkSize    ( vec, 7UL );
   checkNonZeros( vec, 0UL );

   // Inserting a non-zero element
   {
      Iterator pos = vec.insert( 4UL, 1 );

      checkSize    ( vec, 7UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      if( pos->value() != 1 || pos->index() != 4UL ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid iterator returned\n"
             << " Details:\n"
             << "   Value: " << pos->value() << "\n"
             << "   Index: " << pos->index() << "\n"
             << "   Expected value: 1\n"
             << "   Expected index: 4\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec[4] != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Inserting an element failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 0 0 0 1 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Inserting a second non-zero element
   {
      Iterator pos = vec.insert( 6UL, 2 );

      checkSize    ( vec, 7UL );
      checkCapacity( vec, 2UL );
      checkNonZeros( vec, 2UL );

      if( pos->value() != 2 || pos->index() != 6UL ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid iterator returned\n"
             << " Details:\n"
             << "   Value: " << pos->value() << "\n"
             << "   Index: " << pos->index() << "\n"
             << "   Expected value: 2\n"
             << "   Expected index: 6\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec[4] != 1 || vec[6] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Inserting an element failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 0 0 0 1 0 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Inserting a third non-zero element
   {
      Iterator pos = vec.insert( 2UL, 3 );

      checkSize    ( vec, 7UL );
      checkCapacity( vec, 3UL );
      checkNonZeros( vec, 3UL );

      if( pos->value() != 3 || pos->index() != 2UL ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid iterator returned\n"
             << " Details:\n"
             << "   Value: " << pos->value() << "\n"
             << "   Index: " << pos->index() << "\n"
             << "   Expected value: 3\n"
             << "   Expected index: 2\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec[2] != 3 || vec[4] != 1 || vec[6] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Inserting an element failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 0 3 0 1 0 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Inserting a fourth non-zero element
   {
      Iterator pos = vec.insert( 3UL, 4 );

      checkSize    ( vec, 7UL );
      checkCapacity( vec, 4UL );
      checkNonZeros( vec, 4UL );

      if( pos->value() != 4 || pos->index() != 3UL ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid iterator returned\n"
             << " Details:\n"
             << "   Value: " << pos->value() << "\n"
             << "   Index: " << pos->index() << "\n"
             << "   Expected value: 4\n"
             << "   Expected index: 3\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec[2] != 3 || vec[3] != 4 || vec[4] != 1 || vec[6] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Inserting an element failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 0 3 4 1 0 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Trying to insert an already existing element
   try {
      vec.insert( 3UL, 5 );

      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Inserting an existing element succeeded\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 0 0 3 4 1 0 2 )\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::invalid_argument& ) {}
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c append() member function of the CompressedVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c append() member function of the CompressedVector
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testAppend()
{
   test_ = "CompressedVector::append()";

   // Initialization check
   blaze::CompressedVector<int,blaze::rowVector> vec( 9UL, 4UL );

   checkSize    ( vec, 9UL );
   checkCapacity( vec, 4UL );
   checkNonZeros( vec, 0UL );

   // Appending one non-zero element
   vec.append( 1UL, 1 );

   checkSize    ( vec, 9UL );
   checkCapacity( vec, 4UL );
   checkNonZeros( vec, 1UL );

   if( vec[1] != 1 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Append operation failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 0 1 0 0 0 0 0 0 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Appending three more non-zero elements
   vec.append( 3UL, 2 );
   vec.append( 4UL, 3 );
   vec.append( 8UL, 4 );

   checkSize    ( vec, 9UL );
   checkCapacity( vec, 4UL );
   checkNonZeros( vec, 4UL );

   if( vec[1] != 1 || vec[3] != 2 || vec[4] != 3 || vec[8] != 4 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Append operation failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 0 1 0 2 3 0 0 0 4 )\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c erase() member function of the CompressedVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c erase() member function of the CompressedVector
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testErase()
{
   //=====================================================================================
   // Index-based erase function
   //=====================================================================================

   {
      test_ = "CompressedVector::erase( size_t )";

      // Initialization check
      blaze::CompressedVector<int,blaze::rowVector> vec{ 1, 0, 2, 0, 0, 3, 0, 4, 5 };

      checkSize    ( vec, 9UL );
      checkCapacity( vec, 5UL );
      checkNonZeros( vec, 5UL );

      if( vec[0] != 1 || vec[2] != 2 || vec[5] != 3 || vec[7] != 4 || vec[8] != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 0 2 0 0 3 0 4 5 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element at index 0
      vec.erase( size_t(0) );

      checkSize    ( vec, 9UL );
      checkCapacity( vec, 5UL );
      checkNonZeros( vec, 4UL );

      if( vec[2] != 2 || vec[5] != 3 || vec[7] != 4 || vec[8] != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 0 2 0 0 3 0 4 5 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element at index 8
      vec.erase( 8UL );

      checkSize    ( vec, 9UL );
      checkCapacity( vec, 5UL );
      checkNonZeros( vec, 3UL );

      if( vec[2] != 2 || vec[5] != 3 || vec[7] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 0 2 0 0 3 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element at index 5
      vec.erase( 5UL );

      checkSize    ( vec, 9UL );
      checkCapacity( vec, 5UL );
      checkNonZeros( vec, 2UL );

      if( vec[2] != 2 || vec[7] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 0 2 0 0 0 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to erase a zero element
      vec.erase( 1UL );

      checkSize    ( vec, 9UL );
      checkCapacity( vec, 5UL );
      checkNonZeros( vec, 2UL );

      if( vec[2] != 2 || vec[7] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 0 2 0 0 0 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Iterator-based erase() function
   //=====================================================================================

   {
      test_ = "CompressedVector::erase( Iterator )";

      using VectorType = blaze::CompressedVector<int,blaze::rowVector>;
      using Iterator   = VectorType::Iterator;

      // Initialization check
      VectorType vec{ 1, 0, 2, 0, 0, 3, 0, 4, 5 };

      checkSize    ( vec, 9UL );
      checkCapacity( vec, 5UL );
      checkNonZeros( vec, 5UL );

      if( vec[0] != 1 || vec[2] != 2 || vec[5] != 3 || vec[7] != 4 || vec[8] != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 0 2 0 0 3 0 4 5 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element at index 0
      {
         Iterator pos = vec.erase( vec.find( 0UL ) );

         checkSize    ( vec, 9UL );
         checkCapacity( vec, 5UL );
         checkNonZeros( vec, 4UL );

         if( vec[2] != 2 || vec[5] != 3 || vec[7] != 4 || vec[8] != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << vec << "\n"
                << "   Expected result:\n( 0 0 2 0 0 3 0 4 5 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos->value() != 2 || pos->index() != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 2\n"
                << "   Expected index: 2\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the element at index 8
      {
         Iterator pos = vec.erase( vec.find( 8UL ) );

         checkSize    ( vec, 9UL );
         checkCapacity( vec, 5UL );
         checkNonZeros( vec, 3UL );

         if( vec[2] != 2 || vec[5] != 3 || vec[7] != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << vec << "\n"
                << "   Expected result:\n( 0 0 2 0 0 3 0 4 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos != vec.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the element at index 5
      {
         Iterator pos = vec.erase( vec.find( 5UL ) );

         checkSize    ( vec, 9UL );
         checkCapacity( vec, 5UL );
         checkNonZeros( vec, 2UL );

         if( vec[2] != 2 || vec[7] != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << vec << "\n"
                << "   Expected result:\n( 0 0 2 0 0 0 0 4 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos->value() != 4 || pos->index() != 7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 4\n"
                << "   Expected index: 7\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to erase a zero element
      {
         Iterator pos = vec.erase( vec.find( 1UL ) );

         checkSize    ( vec, 9UL );
         checkCapacity( vec, 5UL );
         checkNonZeros( vec, 2UL );

         if( vec[2] != 2 || vec[7] != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << vec << "\n"
                << "   Expected result:\n( 0 0 2 0 0 0 0 4 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos != vec.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Iterator-range-based erase() function
   //=====================================================================================

   {
      test_ = "CompressedVector::erase( Iterator, Iterator )";

      using VectorType = blaze::CompressedVector<int,blaze::rowVector>;
      using Iterator   = VectorType::Iterator;

      // Initialization check
      VectorType vec{ 1, 0, 2, 0, 0, 3, 0, 4, 5 };

      checkSize    ( vec, 9UL );
      checkCapacity( vec, 5UL );
      checkNonZeros( vec, 5UL );

      if( vec[0] != 1 || vec[2] != 2 || vec[5] != 3 || vec[7] != 4 || vec[8] != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 0 2 0 0 3 0 4 5 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the range from index 0 to index 2
      {
         Iterator pos = vec.erase( vec.find( 0UL ), vec.find( 2UL ) );

         checkSize    ( vec, 9UL );
         checkCapacity( vec, 5UL );
         checkNonZeros( vec, 4UL );

         if( vec[2] != 2 || vec[5] != 3 || vec[7] != 4 || vec[8] != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a single-element range failed\n"
                << " Details:\n"
                << "   Result:\n" << vec << "\n"
                << "   Expected result:\n( 0 0 2 0 0 3 0 4 5 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos->value() != 2 || pos->index() != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 2\n"
                << "   Expected index: 2\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the range from index 8 to the end
      {
         Iterator pos = vec.erase( vec.find( 8UL ), vec.end() );

         checkSize    ( vec, 9UL );
         checkCapacity( vec, 5UL );
         checkNonZeros( vec, 3UL );

         if( vec[2] != 2 || vec[5] != 3 || vec[7] != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a single-element range failed\n"
                << " Details:\n"
                << "   Result:\n" << vec << "\n"
                << "   Expected result:\n( 0 0 2 0 0 3 0 4 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos != vec.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the range from index 5 to index 7
      {
         Iterator pos = vec.erase( vec.find( 5UL ), vec.find( 7UL ) );

         checkSize    ( vec, 9UL );
         checkCapacity( vec, 5UL );
         checkNonZeros( vec, 2UL );

         if( vec[2] != 2 || vec[7] != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a single-element range failed\n"
                << " Details:\n"
                << "   Result:\n" << vec << "\n"
                << "   Expected result:\n( 0 0 2 0 0 0 0 4 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos->value() != 4 || pos->index() != 7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 4\n"
                << "   Expected index: 7\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to erase an empty range
      {
         Iterator pos = vec.erase( vec.find( 2UL ), vec.find( 2UL ) );

         checkSize    ( vec, 9UL );
         checkCapacity( vec, 5UL );
         checkNonZeros( vec, 2UL );

         if( vec[2] != 2 || vec[7] != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing an empty range failed\n"
                << " Details:\n"
                << "   Result:\n" << vec << "\n"
                << "   Expected result:\n( 0 0 2 0 0 0 0 4 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos != vec.find( 2UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   //  erase() function with predicate
   //=====================================================================================

   {
      test_ = "CompressedVector::erase( Predicate )";

      // Initialization check
      blaze::CompressedVector<int,blaze::rowVector> vec{ 1, 0, 2, 0, 0, 3, 0, 4, 5 };

      checkSize    ( vec, 9UL );
      checkCapacity( vec, 5UL );
      checkNonZeros( vec, 5UL );

      if( vec[0] != 1 || vec[2] != 2 || vec[5] != 3 || vec[7] != 4 || vec[8] != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 0 2 0 0 3 0 4 5 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing a selection of elements
      vec.erase( []( int value ){ return value == 1 || value == 3 || value == 5; } );

      checkSize    ( vec, 9UL );
      checkCapacity( vec, 5UL );
      checkNonZeros( vec, 2UL );

      if( vec[2] != 2 || vec[7] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a selection of elements failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 0 2 0 0 0 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to erase all elements with value 1
      vec.erase( []( int value ){ return value == 1; } );

      checkSize    ( vec, 9UL );
      checkCapacity( vec, 5UL );
      checkNonZeros( vec, 2UL );

      if( vec[2] != 2 || vec[7] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing all elements with value 1 failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 0 2 0 0 0 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Iterator-range-based erase() function with predicate
   //=====================================================================================

   {
      test_ = "CompressedVector::erase( Iterator, Iterator, Predicate )";

      using VectorType = blaze::CompressedVector<int,blaze::rowVector>;

      // Initialization check
      VectorType vec{ 1, 0, 2, 0, 0, 3, 0, 4, 5 };

      checkSize    ( vec, 9UL );
      checkCapacity( vec, 5UL );
      checkNonZeros( vec, 5UL );

      if( vec[0] != 1 || vec[2] != 2 || vec[5] != 3 || vec[7] != 4 || vec[8] != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 0 2 0 0 3 0 4 5 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing a selection of elements
      {
         vec.erase( vec.find( 2UL ), vec.find( 8UL ),
                    []( int value ){ return value == 2 || value == 4; } );

         checkSize    ( vec, 9UL );
         checkCapacity( vec, 5UL );
         checkNonZeros( vec, 3UL );

         if( vec[0] != 1 || vec[5] != 3 || vec[8] != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a selection of elements failed\n"
                << " Details:\n"
                << "   Result:\n" << vec << "\n"
                << "   Expected result:\n( 1 0 0 0 0 3 0 0 5 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to erase within an empty range
      {
         vec.erase( vec.find( 5UL ), vec.find( 5UL ), []( int ){ return true; } );

         checkSize    ( vec, 9UL );
         checkCapacity( vec, 5UL );
         checkNonZeros( vec, 3UL );

         if( vec[0] != 1 || vec[5] != 3 || vec[8] != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing an empty range failed\n"
                << " Details:\n"
                << "   Result:\n" << vec << "\n"
                << "   Expected result:\n( 1 0 0 0 0 3 0 0 5 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c find() member function of the CompressedVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c find() member function of the CompressedVector
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testFind()
{
   test_ = "CompressedVector::find()";

   using ConstIterator = blaze::CompressedVector<int,blaze::rowVector>::ConstIterator;

   // Initialization check
   blaze::CompressedVector<int,blaze::rowVector> vec{ 1, 0, 2, 0, 0, 0, 0, 3 };

   checkSize    ( vec, 8UL );
   checkCapacity( vec, 3UL );
   checkNonZeros( vec, 3UL );

   // Searching for the first element
   {
      ConstIterator pos( vec.find( 0UL ) );

      if( pos == vec.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Element could not be found\n"
             << " Details:\n"
             << "   Required index = 0\n"
             << "   Current vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
      else if( pos->index() != 0 || pos->value() != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Wrong element found\n"
             << " Details:\n"
             << "   Required index = 0\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = 1\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Searching for the second element
   {
      ConstIterator pos( vec.find( 2UL ) );

      if( pos == vec.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Element could not be found\n"
             << " Details:\n"
             << "   Required index = 2\n"
             << "   Current vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
      else if( pos->index() != 2 || pos->value() != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Wrong element found\n"
             << " Details:\n"
             << "   Required index = 2\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = 2\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Searching for the third element
   {
      ConstIterator pos( vec.find( 7UL ) );

      if( pos == vec.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Element could not be found\n"
             << " Details:\n"
             << "   Required index = 7\n"
             << "   Current vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
      else if( pos->index() != 7 || pos->value() != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Wrong element found\n"
             << " Details:\n"
             << "   Required index = 7\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = 3\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Searching for a non-existing non-zero element
   {
      ConstIterator pos( vec.find( 5UL ) );

      if( pos != vec.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Non-existing element could be found\n"
             << " Details:\n"
             << "   Required index = 5\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = 0\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c lowerBound() member function of the CompressedVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c lowerBound() member function of the CompressedVector
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testLowerBound()
{
   test_ = "CompressedVector::lowerBound()";

   using ConstIterator = blaze::CompressedVector<int,blaze::rowVector>::ConstIterator;

   // Initialization check
   blaze::CompressedVector<int,blaze::rowVector> vec{ 1, 0, 2, 0, 0, 0, 0, 3 };

   checkSize    ( vec, 8UL );
   checkCapacity( vec, 3UL );
   checkNonZeros( vec, 3UL );

   // Determining the lower bound for index 0
   {
      ConstIterator pos( vec.lowerBound( 0UL ) );

      if( pos == vec.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Lower bound could not be determined\n"
             << " Details:\n"
             << "   Required index = 0\n"
             << "   Current vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
      else if( pos->index() != 0 || pos->value() != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Wrong element found\n"
             << " Details:\n"
             << "   Required index = 0\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = 1\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Determining the lower bound for index 1
   {
      ConstIterator pos( vec.lowerBound( 1UL ) );

      if( pos == vec.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Lower bound could not be determined\n"
             << " Details:\n"
             << "   Required index = 1\n"
             << "   Current vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
      else if( pos->index() != 2 || pos->value() != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Wrong element found\n"
             << " Details:\n"
             << "   Required index = 2\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = 2\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Determining the lower bound for index 2
   {
      ConstIterator pos( vec.lowerBound( 2UL ) );

      if( pos == vec.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Lower bound could not be determined\n"
             << " Details:\n"
             << "   Required index = 2\n"
             << "   Current vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
      else if( pos->index() != 2 || pos->value() != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Wrong element found\n"
             << " Details:\n"
             << "   Required index = 2\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = 2\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Determining the lower bound for index 3
   {
      ConstIterator pos( vec.lowerBound( 3UL ) );

      if( pos == vec.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Lower bound could not be determined\n"
             << " Details:\n"
             << "   Required index = 3\n"
             << "   Current vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
      else if( pos->index() != 7 || pos->value() != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Wrong element found\n"
             << " Details:\n"
             << "   Required index = 7\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = 3\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Determining the lower bound for index 7
   {
      ConstIterator pos( vec.lowerBound( 7UL ) );

      if( pos == vec.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Lower bound could not be determined\n"
             << " Details:\n"
             << "   Required index = 7\n"
             << "   Current vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
      else if( pos->index() != 7 || pos->value() != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Wrong element found\n"
             << " Details:\n"
             << "   Required index = 7\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = 3\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c upperBound() member function of the CompressedVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c upperBound() member function of the CompressedVector
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testUpperBound()
{
   test_ = "CompressedVector::upperBound()";

   using ConstIterator = blaze::CompressedVector<int,blaze::rowVector>::ConstIterator;

   // Initialization check
   blaze::CompressedVector<int,blaze::rowVector> vec{ 1, 0, 2, 0, 0, 0, 0, 3 };

   checkSize    ( vec, 8UL );
   checkCapacity( vec, 3UL );
   checkNonZeros( vec, 3UL );

   // Determining the upper bound for index 0
   {
      ConstIterator pos( vec.upperBound( 0UL ) );

      if( pos == vec.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Upper bound could not be determined\n"
             << " Details:\n"
             << "   Required index = 0\n"
             << "   Current vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
      else if( pos->index() != 2 || pos->value() != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Wrong element found\n"
             << " Details:\n"
             << "   Required index = 2\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = 2\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Determining the upper bound for index 1
   {
      ConstIterator pos( vec.upperBound( 1UL ) );

      if( pos == vec.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Upper bound could not be determined\n"
             << " Details:\n"
             << "   Required index = 1\n"
             << "   Current vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
      else if( pos->index() != 2 || pos->value() != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Wrong element found\n"
             << " Details:\n"
             << "   Required index = 2\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = 2\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Determining the upper bound for index 2
   {
      ConstIterator pos( vec.upperBound( 2UL ) );

      if( pos == vec.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Upper bound could not be determined\n"
             << " Details:\n"
             << "   Required index = 2\n"
             << "   Current vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
      else if( pos->index() != 7 || pos->value() != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Wrong element found\n"
             << " Details:\n"
             << "   Required index = 7\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = 3\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Determining the upper bound for index 3
   {
      ConstIterator pos( vec.upperBound( 3UL ) );

      if( pos == vec.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Upper bound could not be determined\n"
             << " Details:\n"
             << "   Required index = 3\n"
             << "   Current vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
      else if( pos->index() != 7 || pos->value() != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Wrong element found\n"
             << " Details:\n"
             << "   Required index = 7\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = 3\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Determining the upper bound for index 7
   {
      ConstIterator pos( vec.upperBound( 7UL ) );

      if( pos != vec.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Upper bound could not be determined\n"
             << " Details:\n"
             << "   Required index = 7\n"
             << "   Current vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c isDefault() function with the CompressedVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isDefault() function with the CompressedVector class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testIsDefault()
{
   test_ = "isDefault() function";

   // isDefault with vector of size 0
   {
      blaze::CompressedVector<int,blaze::rowVector> vec;

      if( isDefault( vec ) != true ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isDefault evaluation\n"
             << " Details:\n"
             << "   Vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // isDefault with default vector
   {
      blaze::CompressedVector<int,blaze::rowVector> vec( 3UL );

      if( isDefault( vec[1] ) != true ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isDefault evaluation\n"
             << " Details:\n"
             << "   Vector element: " << vec[1] << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( isDefault( vec ) != false ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isDefault evaluation\n"
             << " Details:\n"
             << "   Vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // isDefault with non-default vector
   {
      blaze::CompressedVector<int,blaze::rowVector> vec{ 0, 1, 0 };

      if( isDefault( vec[1] ) != false ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isDefault evaluation\n"
             << " Details:\n"
             << "   Vector element: " << vec[1] << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( isDefault( vec ) != false ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isDefault evaluation\n"
             << " Details:\n"
             << "   Vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************

} // namespace compressedvector

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
   std::cout << "   Running CompressedVector class test..." << std::endl;

   try
   {
      RUN_COMPRESSEDVECTOR_CLASS_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during CompressedVector class test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
