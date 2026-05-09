//=================================================================================================
/*!
//  \file src/mathtest/vectors/uniformvector/ClassTest.cpp
//  \brief Source file for the UniformVector class test
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
#include <memory>
#include <blaze/math/CompressedVector.h>
#include <blaze/math/CustomVector.h>
#include <blaze/math/DynamicVector.h>
#include <blaze/math/shims/Equal.h>
#include <blaze/util/Complex.h>
#include <blaze/util/policies/Deallocate.h>
#include <blaze/util/Random.h>
#include <blazetest/mathtest/RandomMaximum.h>
#include <blazetest/mathtest/RandomMinimum.h>
#include <blazetest/mathtest/vectors/uniformvector/ClassTest.h>

#ifdef BLAZE_USE_HPX_THREADS
#  include <hpx/hpx_main.hpp>
#endif


namespace blazetest {

namespace mathtest {

namespace vectors {

namespace uniformvector {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the UniformVector class test.
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
   testScaling();
   testSubscript();
   testAt();
   testIterator();
   testNonZeros();
   testReset();
   testClear();
   testResize();
   testExtend();
   testSwap();
   testIsDefault();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the UniformVector constructors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all constructors of the UniformVector class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testConstructors()
{
   //=====================================================================================
   // Default constructor
   //=====================================================================================

   {
      test_ = "UniformVector default constructor";

      blaze::UniformVector<int,blaze::rowVector> vec;

      checkSize    ( vec, 0UL );
      checkNonZeros( vec, 0UL );
   }


   //=====================================================================================
   // Size constructor
   //=====================================================================================

   {
      test_ = "UniformVector size constructor (size 0)";

      blaze::UniformVector<int,blaze::rowVector> vec( 0UL );

      checkSize    ( vec, 0UL );
      checkNonZeros( vec, 0UL );
   }

   {
      test_ = "UniformVector size constructor (size 3)";

      blaze::UniformVector<int,blaze::rowVector> vec( 3UL );

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 3UL );
      checkNonZeros( vec, 0UL );

      if( vec[0] != 0 || vec[1] != 0 || vec[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Homogeneous initialization
   //=====================================================================================

   {
      test_ = "UniformVector homogeneous initialization constructor (size 0)";

      blaze::UniformVector<int,blaze::rowVector> vec( 0UL, 2 );

      checkSize    ( vec, 0UL );
      checkNonZeros( vec, 0UL );
   }

   {
      test_ = "UniformVector homogeneous initialization constructor (size 3)";

      blaze::UniformVector<int,blaze::rowVector> vec( 3UL, 2 );

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 3UL );
      checkNonZeros( vec, 3UL );

      if( vec[0] != 2 || vec[1] != 2 || vec[2] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Copy constructor
   //=====================================================================================

   {
      test_ = "UniformVector copy constructor (size 0)";

      blaze::UniformVector<int,blaze::rowVector> vec1( 0UL );
      blaze::UniformVector<int,blaze::rowVector> vec2( vec1 );

      checkSize    ( vec2, 0UL );
      checkNonZeros( vec2, 0UL );
   }

   {
      test_ = "UniformVector copy constructor (size 5)";

      blaze::UniformVector<int,blaze::rowVector> vec1( 5, 2 );
      blaze::UniformVector<int,blaze::rowVector> vec2( vec1 );

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 5UL );

      if( vec2[0] != 2 || vec2[1] != 2 || vec2[2] != 2 || vec2[3] != 2 || vec2[4] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 2 2 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Move constructor
   //=====================================================================================

   {
      test_ = "UniformVector move constructor (size 0)";

      blaze::UniformVector<int,blaze::rowVector> vec1( 0UL );
      blaze::UniformVector<int,blaze::rowVector> vec2( std::move( vec1 ) );

      checkSize    ( vec2, 0UL );
      checkNonZeros( vec2, 0UL );
   }

   {
      test_ = "UniformVector move constructor (size 5)";

      blaze::UniformVector<int,blaze::rowVector> vec1( 5, 2 );
      blaze::UniformVector<int,blaze::rowVector> vec2( std::move( vec1 ) );

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 5UL );

      if( vec2[0] != 2 || vec2[1] != 2 || vec2[2] != 2 || vec2[3] != 2 || vec2[4] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 2 2 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Dense vector constructor
   //=====================================================================================

   {
      test_ = "UniformVector dense vector constructor (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::rowVector;

      using AlignedPadded = blaze::CustomVector<int,aligned,padded,rowVector>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 16UL ) );
      AlignedPadded vec1( memory.get(), 5UL, 16UL );
      vec1 = 2;

      blaze::UniformVector<int,blaze::rowVector> vec2( vec1 );

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 5UL );

      if( vec2[0] != 2 || vec2[1] != 2 || vec2[2] != 2 || vec2[3] != 2 || vec2[4] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 2 2 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "UniformVector dense vector constructor (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowVector;

      using UnalignedUnpadded = blaze::CustomVector<int,unaligned,unpadded,rowVector>;
      std::unique_ptr<int[]> memory( new int[6] );
      UnalignedUnpadded vec1( memory.get()+1UL, 5UL );
      vec1 = 2;

      blaze::UniformVector<int,blaze::rowVector> vec2( vec1 );

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 5UL );

      if( vec2[0] != 2 || vec2[1] != 2 || vec2[2] != 2 || vec2[3] != 2 || vec2[4] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 2 2 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "UniformVector dense vector constructor (non-uniform)";

      blaze::DynamicVector<int,blaze::rowVector> vec1{ 2, 2, 2, 0, 2 };

      try {
         blaze::UniformVector<int,blaze::rowVector> vec2( vec1 );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of non-uniform UniformVector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Sparse vector constructor
   //=====================================================================================

   {
      test_ = "UniformVector sparse vector constructor (uniform)";

      blaze::CompressedVector<int,blaze::rowVector> vec1{ 2, 2, 2, 2, 2 };
      blaze::UniformVector<int,blaze::rowVector> vec2( vec1 );

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 5UL );

      if( vec2[0] != 2 || vec2[1] != 2 || vec2[2] != 2 || vec2[3] != 2 || vec2[4] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 2 2 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "UniformVector sparse vector constructor (non-uniform)";

      blaze::CompressedVector<int,blaze::rowVector> vec1{ 2, 2, 2, 0, 2 };

      try {
         blaze::UniformVector<int,blaze::rowVector> vec2( vec1 );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of non-uniform UniformVector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the UniformVector assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all assignment operators of the UniformVector class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testAssignment()
{
   //=====================================================================================
   // Homogeneous assignment
   //=====================================================================================

   {
      test_ = "UniformVector homogeneous assignment";

      blaze::UniformVector<int,blaze::rowVector> vec( 3UL );
      vec = 2;

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 3UL );
      checkNonZeros( vec, 3UL );

      if( vec[0] != 2 || vec[1] != 2 || vec[2] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Copy assignment
   //=====================================================================================

   {
      test_ = "UniformVector copy assignment";

      blaze::UniformVector<int,blaze::rowVector> vec1( 5, 2 );
      blaze::UniformVector<int,blaze::rowVector> vec2;
      vec2 = vec1;

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 5UL );

      if( vec2[0] != 2 || vec2[1] != 2 || vec2[2] != 2 || vec2[3] != 2 || vec2[4] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 2 2 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "UniformVector copy assignment stress test";

      using RandomVectorType = blaze::UniformVector<int,blaze::rowVector>;

      blaze::UniformVector<int,blaze::rowVector> vec1;
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
      test_ = "UniformVector move assignment";

      blaze::UniformVector<int,blaze::rowVector> vec1( 5, 2 );
      blaze::UniformVector<int,blaze::rowVector> vec2( 3, 4 );

      vec2 = std::move( vec1 );

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 5UL );

      if( vec2[0] != 2 || vec2[1] != 2 || vec2[2] != 2 || vec2[3] != 2 || vec2[4] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 2 2 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Dense vector assignment
   //=====================================================================================

   {
      test_ = "UniformVector dense vector assignment (mixed type)";

      blaze::UniformVector<short,blaze::rowVector> vec1( 5, 2 );
      blaze::UniformVector<int,blaze::rowVector> vec2;
      vec2 = vec1;

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 5UL );

      if( vec2[0] != 2 || vec2[1] != 2 || vec2[2] != 2 || vec2[3] != 2 || vec2[4] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 2 2 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "UniformVector dense vector assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::rowVector;

      using AlignedPadded = blaze::CustomVector<int,aligned,padded,rowVector>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 16UL ) );
      AlignedPadded vec1( memory.get(), 5UL, 16UL );
      vec1 = 2;

      blaze::UniformVector<int,blaze::rowVector> vec2;
      vec2 = vec1;

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 5UL );

      if( vec2[0] != 2 || vec2[1] != 2 || vec2[2] != 2 || vec2[3] != 2 || vec2[4] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 2 2 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "UniformVector dense vector assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowVector;

      using UnalignedUnpadded = blaze::CustomVector<int,unaligned,unpadded,rowVector>;
      std::unique_ptr<int[]> memory( new int[6] );
      UnalignedUnpadded vec1( memory.get()+1UL, 5UL );
      vec1 = 2;

      blaze::UniformVector<int,blaze::rowVector> vec2;
      vec2 = vec1;

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 5UL );

      if( vec2[0] != 2 || vec2[1] != 2 || vec2[2] != 2 || vec2[3] != 2 || vec2[4] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 2 2 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "UniformVector dense vector assignment (non-uniform)";

      blaze::DynamicVector<int,blaze::rowVector> vec1{ 2, 2, 2, 0, 2 };

      try {
         blaze::UniformVector<int,blaze::rowVector> vec2;
         vec2 = vec1;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-uniform dense vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Sparse vector assignment
   //=====================================================================================

   {
      test_ = "UniformVector sparse vector assignment (uniform)";

      blaze::CompressedVector<int,blaze::rowVector> vec1{ 2, 2, 2, 2, 2 };
      blaze::UniformVector<int,blaze::rowVector> vec2;
      vec2 = vec1;

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 5UL );

      if( vec2[0] != 2 || vec2[1] != 2 || vec2[2] != 2 || vec2[3] != 2 || vec2[4] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 2 2 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "UniformVector sparse vector assignment (non-uniform)";

      blaze::CompressedVector<int,blaze::rowVector> vec1{ 2, 2, 2, 0, 2 };

      try {
         blaze::UniformVector<int,blaze::rowVector> vec2;
         vec2 = vec1;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-uniform sparse vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the UniformVector addition assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the addition assignment operators of the UniformVector class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testAddAssign()
{
   //=====================================================================================
   // Dense vector addition assignment
   //=====================================================================================

   {
      test_ = "UniformVector dense vector addition assignment (mixed type)";

      blaze::UniformVector<short,blaze::rowVector> vec1( 5, 2 );
      blaze::UniformVector<int,blaze::rowVector> vec2( 5, 1 );

      vec2 += vec1;

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 5UL );

      if( vec2[0] != 3 || vec2[1] != 3 || vec2[2] != 3 || vec2[3] != 3 || vec2[4] != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 3 3 3 3 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "UniformVector dense vector addition assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::rowVector;

      using AlignedPadded = blaze::CustomVector<int,aligned,padded,rowVector>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 16UL ) );
      AlignedPadded vec1( memory.get(), 5UL, 16UL );
      vec1 = 2;

      blaze::UniformVector<int,blaze::rowVector> vec2( 5UL, 1 );

      vec2 += vec1;

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 5UL );

      if( vec2[0] != 3 || vec2[1] != 3 || vec2[2] != 3 || vec2[3] != 3 || vec2[4] != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 3 3 3 3 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "UniformVector dense vector addition assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowVector;

      using UnalignedUnpadded = blaze::CustomVector<int,unaligned,unpadded,rowVector>;
      std::unique_ptr<int[]> memory( new int[6] );
      UnalignedUnpadded vec1( memory.get()+1UL, 5UL );
      vec1 = 2;

      blaze::UniformVector<int,blaze::rowVector> vec2( 5UL, 1 );

      vec2 += vec1;

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 5UL );

      if( vec2[0] != 3 || vec2[1] != 3 || vec2[2] != 3 || vec2[3] != 3 || vec2[4] != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 3 3 3 3 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "UniformVector dense vector addition assignment (non-uniform)";

      blaze::DynamicVector<int,blaze::rowVector> vec1{ 2, 2, 2, 0, 2 };

      try {
         blaze::UniformVector<int,blaze::rowVector> vec2( 5UL, 1 );
         vec2 += vec1;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-uniform dense vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Sparse vector addition assignment
   //=====================================================================================

   {
      test_ = "UniformVector sparse vector addition assignment";

      blaze::CompressedVector<int,blaze::rowVector> vec1{ 2, 2, 2, 2, 2 };
      blaze::UniformVector<int,blaze::rowVector> vec2( 5UL, 1 );

      vec2 += vec1;

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 5UL );

      if( vec2[0] != 3 || vec2[1] != 3 || vec2[2] != 3 || vec2[3] != 3 || vec2[4] != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 3 3 3 3 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "UniformVector sparse vector addition assignment (non-uniform)";

      blaze::CompressedVector<int,blaze::rowVector> vec1{ 2, 2, 2, 0, 2 };

      try {
         blaze::UniformVector<int,blaze::rowVector> vec2( 5UL, 1 );
         vec2 += vec1;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-uniform sparse vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the UniformVector subtraction assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the subtraction assignment operators of the UniformVector class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testSubAssign()
{
   //=====================================================================================
   // Dense vector subtraction assignment
   //=====================================================================================

   {
      test_ = "UniformVector dense vector subtraction assignment (mixed type)";

      blaze::UniformVector<short,blaze::rowVector> vec1( 5, 2 );
      blaze::UniformVector<int,blaze::rowVector> vec2( 5, 1 );

      vec2 -= vec1;

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 5UL );

      if( vec2[0] != -1 || vec2[1] != -1 || vec2[2] != -1 || vec2[3] != -1 || vec2[4] != -1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( -1 -1 -1 -1 -1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "UniformVector dense vector subtraction assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::rowVector;

      using AlignedPadded = blaze::CustomVector<int,aligned,padded,rowVector>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 16UL ) );
      AlignedPadded vec1( memory.get(), 5UL, 16UL );
      vec1 = 2;

      blaze::UniformVector<int,blaze::rowVector> vec2( 5UL, 1 );

      vec2 -= vec1;

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 5UL );

      if( vec2[0] != -1 || vec2[1] != -1 || vec2[2] != -1 || vec2[3] != -1 || vec2[4] != -1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( -1 -1 -1 -1 -1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "UniformVector dense vector subtraction assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowVector;

      using UnalignedUnpadded = blaze::CustomVector<int,unaligned,unpadded,rowVector>;
      std::unique_ptr<int[]> memory( new int[6] );
      UnalignedUnpadded vec1( memory.get()+1UL, 5UL );
      vec1 = 2;

      blaze::UniformVector<int,blaze::rowVector> vec2( 5UL, 1 );

      vec2 -= vec1;

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 5UL );

      if( vec2[0] != -1 || vec2[1] != -1 || vec2[2] != -1 || vec2[3] != -1 || vec2[4] != -1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( -1 -1 -1 -1 -1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "UniformVector dense vector subtraction assignment (non-uniform)";

      blaze::DynamicVector<int,blaze::rowVector> vec1{ 2, 2, 2, 0, 2 };

      try {
         blaze::UniformVector<int,blaze::rowVector> vec2( 5UL, 1 );
         vec2 -= vec1;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-uniform dense vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Sparse vector subtraction assignment
   //=====================================================================================

   {
      test_ = "UniformVector sparse vector subtraction assignment";

      blaze::CompressedVector<int,blaze::rowVector> vec1{ 2, 2, 2, 2, 2 };
      blaze::UniformVector<int,blaze::rowVector> vec2( 5UL, 1 );

      vec2 -= vec1;

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 5UL );

      if( vec2[0] != -1 || vec2[1] != -1 || vec2[2] != -1 || vec2[3] != -1 || vec2[4] != -1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( -1 -1 -1 -1 -1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "UniformVector sparse vector subtraction assignment (non-uniform)";

      blaze::CompressedVector<int,blaze::rowVector> vec1{ 2, 2, 2, 0, 2 };

      try {
         blaze::UniformVector<int,blaze::rowVector> vec2( 5UL, 1 );
         vec2 -= vec1;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-uniform sparse vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the UniformVector multiplication assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the multiplication assignment operators of the UniformVector
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testMultAssign()
{
   //=====================================================================================
   // Dense vector multiplication assignment
   //=====================================================================================

   {
      test_ = "UniformVector dense vector multiplication assignment (mixed type)";

      blaze::UniformVector<short,blaze::rowVector> vec1( 5, 2 );
      blaze::UniformVector<int,blaze::rowVector> vec2( 5, 1 );

      vec2 *= vec1;

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 5UL );

      if( vec2[0] != 2 || vec2[1] != 2 || vec2[2] != 2 || vec2[3] != 2 || vec2[4] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 2 2 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "UniformVector dense vector multiplication assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::rowVector;

      using AlignedPadded = blaze::CustomVector<int,aligned,padded,rowVector>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 16UL ) );
      AlignedPadded vec1( memory.get(), 5UL, 16UL );
      vec1 = 2;

      blaze::UniformVector<int,blaze::rowVector> vec2( 5UL, 1 );

      vec2 *= vec1;

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 5UL );

      if( vec2[0] != 2 || vec2[1] != 2 || vec2[2] != 2 || vec2[3] != 2 || vec2[4] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 2 2 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "UniformVector dense vector multiplication assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowVector;

      using UnalignedUnpadded = blaze::CustomVector<int,unaligned,unpadded,rowVector>;
      std::unique_ptr<int[]> memory( new int[6] );
      UnalignedUnpadded vec1( memory.get()+1UL, 5UL );
      vec1 = 2;

      blaze::UniformVector<int,blaze::rowVector> vec2( 5UL, 1 );

      vec2 *= vec1;

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 5UL );

      if( vec2[0] != 2 || vec2[1] != 2 || vec2[2] != 2 || vec2[3] != 2 || vec2[4] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 2 2 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "UniformVector dense vector multiplication assignment (non-uniform)";

      blaze::DynamicVector<int,blaze::rowVector> vec1{ 2, 2, 2, 0, 2 };

      try {
         blaze::UniformVector<int,blaze::rowVector> vec2( 5UL, 1 );
         vec2 *= vec1;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-uniform dense vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Sparse vector multiplication assignment
   //=====================================================================================

   {
      test_ = "UniformVector sparse vector multiplication assignment";

      blaze::CompressedVector<int,blaze::rowVector> vec1{ 2, 2, 2, 2, 2 };
      blaze::UniformVector<int,blaze::rowVector> vec2( 5UL, 1 );

      vec2 *= vec1;

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 5UL );

      if( vec2[0] != 2 || vec2[1] != 2 || vec2[2] != 2 || vec2[3] != 2 || vec2[4] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 2 2 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "UniformVector sparse vector multiplication assignment (non-uniform)";

      blaze::CompressedVector<int,blaze::rowVector> vec1{ 2, 2, 2, 0, 2 };

      try {
         blaze::UniformVector<int,blaze::rowVector> vec2( 5UL, 1 );
         vec2 *= vec1;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-uniform sparse vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the UniformVector division assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the division assignment operators of the UniformVector
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testDivAssign()
{
   //=====================================================================================
   // Dense vector division assignment
   //=====================================================================================

   {
      test_ = "UniformVector dense vector division assignment (mixed type)";

      blaze::UniformVector<short,blaze::rowVector> vec1( 5, 3 );
      blaze::UniformVector<int,blaze::rowVector> vec2( 5, 6 );

      vec2 /= vec1;

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 5UL );

      if( vec2[0] != 2 || vec2[1] != 2 || vec2[2] != 2 || vec2[3] != 2 || vec2[4] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 2 2 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "UniformVector dense vector division assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::rowVector;

      using AlignedPadded = blaze::CustomVector<int,aligned,padded,rowVector>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 16UL ) );
      AlignedPadded vec1( memory.get(), 5UL, 16UL );
      vec1 = 3;

      blaze::UniformVector<int,blaze::rowVector> vec2( 5UL, 6 );

      vec2 /= vec1;

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 5UL );

      if( vec2[0] != 2 || vec2[1] != 2 || vec2[2] != 2 || vec2[3] != 2 || vec2[4] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 2 2 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "UniformVector dense vector division assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowVector;

      using UnalignedUnpadded = blaze::CustomVector<int,unaligned,unpadded,rowVector>;
      std::unique_ptr<int[]> memory( new int[6] );
      UnalignedUnpadded vec1( memory.get()+1UL, 5UL );
      vec1 = 3;

      blaze::UniformVector<int,blaze::rowVector> vec2( 5UL, 6 );

      vec2 /= vec1;

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 5UL );

      if( vec2[0] != 2 || vec2[1] != 2 || vec2[2] != 2 || vec2[3] != 2 || vec2[4] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 2 2 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "UniformVector dense vector multiplication assignment (non-uniform)";

      blaze::DynamicVector<int,blaze::rowVector> vec1{ 3, 3, 3, 0, 3 };

      try {
         blaze::UniformVector<int,blaze::rowVector> vec2( 5UL, 1 );
         vec2 /= vec1;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-uniform dense vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of all UniformVector (self-)scaling operations.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all available ways to scale an instance of the UniformVector
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testScaling()
{
   //=====================================================================================
   // Self-scaling (v*=s)
   //=====================================================================================

   {
      test_ = "UniformVector self-scaling (v*=s)";

      blaze::UniformVector<int,blaze::rowVector> vec( 5UL, 2 );

      vec *= 2;

      checkSize    ( vec, 5UL );
      checkCapacity( vec, 5UL );
      checkNonZeros( vec, 5UL );

      if( vec[0] != 4 || vec[1] != 4 || vec[2] != 4 || vec[3] != 4 || vec[4] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 4 4 4 4 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Self-scaling (v=v*s)
   //=====================================================================================

   {
      test_ = "UniformVector self-scaling (v=v*s)";

      blaze::UniformVector<int,blaze::rowVector> vec( 5UL, 2 );

      vec = vec * 2;

      checkSize    ( vec, 5UL );
      checkCapacity( vec, 5UL );
      checkNonZeros( vec, 5UL );

      if( vec[0] != 4 || vec[1] != 4 || vec[2] != 4 || vec[3] != 4 || vec[4] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 4 4 4 4 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Self-scaling (v=s*v)
   //=====================================================================================

   {
      test_ = "UniformVector self-scaling (v=s*v)";

      blaze::UniformVector<int,blaze::rowVector> vec( 5UL, 2 );

      vec = 2 * vec;

      checkSize    ( vec, 5UL );
      checkCapacity( vec, 5UL );
      checkNonZeros( vec, 5UL );

      if( vec[0] != 4 || vec[1] != 4 || vec[2] != 4 || vec[3] != 4 || vec[4] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 4 4 4 4 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Self-scaling (v/=s)
   //=====================================================================================

   {
      test_ = "UniformVector self-scaling (v/=s)";

      blaze::UniformVector<int,blaze::rowVector> vec( 5UL, 4 );

      vec /= 2;

      checkSize    ( vec, 5UL );
      checkCapacity( vec, 5UL );
      checkNonZeros( vec, 5UL );

      if( vec[0] != 2 || vec[1] != 2 || vec[2] != 2 || vec[3] != 2 || vec[4] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 2 2 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Self-scaling (v=v/s)
   //=====================================================================================

   {
      test_ = "UniformVector self-scaling (v=v/s)";

      blaze::UniformVector<int,blaze::rowVector> vec( 5UL, 4 );

      vec = vec / 2;

      checkSize    ( vec, 5UL );
      checkCapacity( vec, 5UL );
      checkNonZeros( vec, 5UL );

      if( vec[0] != 2 || vec[1] != 2 || vec[2] != 2 || vec[3] != 2 || vec[4] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 2 2 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // UniformVector::scale()
   //=====================================================================================

   {
      test_ = "UniformVector::scale() (int)";

      // Initialization check
      blaze::UniformVector<int,blaze::rowVector> vec( 4UL, 2 );

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );
      checkNonZeros( vec, 4UL );

      if( vec[0] != 2 || vec[1] != 2 || vec[2] != 2 || vec[3] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 2 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Integral scaling of the vector
      vec.scale( 2 );

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );
      checkNonZeros( vec, 4UL );

      if( vec[0] != 4 || vec[1] != 4 || vec[2] != 4 || vec[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 4 4 4 4 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Floating point scaling of the vector
      vec.scale( 0.5 );

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );
      checkNonZeros( vec, 4UL );

      if( vec[0] != 2 || vec[1] != 2 || vec[2] != 2 || vec[3] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 2 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "UniformVector::scale() (complex)";

      using blaze::complex;

      blaze::UniformVector<complex<float>,blaze::rowVector> vec( 2UL, complex<float>( 2.0F, 0.0F ) );
      vec.scale( complex<float>( 3.0F, 0.0F ) );

      checkSize    ( vec, 2UL );
      checkCapacity( vec, 2UL );
      checkNonZeros( vec, 2UL );

      if( vec[0] != complex<float>( 6.0F, 0.0F ) || vec[1] != complex<float>( 6.0F, 0.0F ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( (6,0) (6,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the UniformVector subscript operator.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of accessing elements via the subscript operator of the
// UniformVector class template. In case an error is detected, a \a std::runtime_error exception
// is thrown.
*/
void ClassTest::testSubscript()
{
   test_ = "UniformVector::operator[]";

   blaze::UniformVector<int,blaze::rowVector> vec( 5, 2 );

   checkSize    ( vec, 5UL );
   checkCapacity( vec, 5UL );
   checkNonZeros( vec, 5UL );

   // Accessing all elements
   if( vec[0] != 2 || vec[1] != 2 || vec[2] != 2 || vec[3] != 2 || vec[4] != 2 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 2 2 2 2 2 )\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c at() member function of the UniformVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of adding and accessing elements via the \c at() member function
// of the UniformVector class template. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
void ClassTest::testAt()
{
   test_ = "UniformVector::at()";

   blaze::UniformVector<int,blaze::rowVector> vec( 5, 2 );

   checkSize    ( vec, 5UL );
   checkCapacity( vec, 5UL );
   checkNonZeros( vec, 5UL );

   // Accessing the elements at index 0 through 4
   if( vec.at(0) != 2 || vec.at(1) != 2 || vec.at(2) != 2 || vec.at(3) != 2 || vec.at(4) != 2 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Access via at() function failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 2 2 2 2 2 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Attempt to assign to the element at index 5
   try {
      vec.at(5);

      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Out-of-bound access succeeded\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 2 2 2 2 2 )\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::out_of_range& ) {}
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the UniformVector iterator implementation.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the iterator implementation of the UniformVector class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testIterator()
{
   using VectorType    = blaze::UniformVector<int>;
   using ConstIterator = VectorType::ConstIterator;

   VectorType vec( 4, 2 );

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

      const ptrdiff_t number( cend( vec ) - cbegin( vec ) );

      if( number != 4L ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid number of elements detected\n"
             << " Details:\n"
             << "   Number of elements         : " << number << "\n"
             << "   Expected number of elements: 4\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Counting the number of elements via ConstIterator (begin-end)
   {
      test_ = "ConstIterator subtraction (begin-end)";

      const ptrdiff_t number( cbegin( vec ) - cend( vec ) );

      if( number != -4L ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid number of elements detected\n"
             << " Details:\n"
             << "   Number of elements         : " << number << "\n"
             << "   Expected number of elements: -4\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing read-only access via ConstIterator
   {
      test_ = "Read-only access via ConstIterator";

      ConstIterator it ( cbegin( vec ) );
      ConstIterator end( cend( vec ) );

      if( it == end || *it != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid initial iterator detected\n";
         throw std::runtime_error( oss.str() );
      }

      ++it;

      if( it == end || *it != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator pre-increment failed\n";
         throw std::runtime_error( oss.str() );
      }

      --it;

      if( it == end || *it != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator pre-decrement failed\n";
         throw std::runtime_error( oss.str() );
      }

      it++;

      if( it == end || *it != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator post-increment failed\n";
         throw std::runtime_error( oss.str() );
      }

      it--;

      if( it == end || *it != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator post-decrement failed\n";
         throw std::runtime_error( oss.str() );
      }

      it += 2UL;

      if( it == end || *it != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator addition assignment failed\n";
         throw std::runtime_error( oss.str() );
      }

      it -= 2UL;

      if( it == end || *it != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator subtraction assignment failed\n";
         throw std::runtime_error( oss.str() );
      }

      it = it + 3UL;

      if( it == end || *it != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator/scalar addition failed\n";
         throw std::runtime_error( oss.str() );
      }

      it = it - 3UL;

      if( it == end || *it != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator/scalar subtraction failed\n";
         throw std::runtime_error( oss.str() );
      }

      it = 4UL + it;

      if( it != end ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scalar/iterator addition failed\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c nonZeros() member function of the UniformVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c nonZeros() member function of the UniformVector class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testNonZeros()
{
   test_ = "UniformVector::nonZeros()";

   {
      blaze::UniformVector<int,blaze::rowVector> vec( 4UL, 0 );

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );
      checkNonZeros( vec, 0UL );

      if( vec[0] != 0 || vec[1] != 0 || vec[2] != 0 || vec[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      blaze::UniformVector<int,blaze::rowVector> vec( 4, 2 );

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );
      checkNonZeros( vec, 4UL );

      if( vec[0] != 2 || vec[1] != 2 || vec[2] != 2 || vec[3] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 2 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c reset() member function of the UniformVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c reset() member function of the UniformVector class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testReset()
{
   using blaze::reset;

   test_ = "UniformVector::reset()";

   // Resetting a default constructed vector
   {
      blaze::UniformVector<int,blaze::rowVector> vec;

      reset( vec );

      checkSize    ( vec, 0UL );
      checkNonZeros( vec, 0UL );
   }

   // Resetting an initialized vector
   {
      blaze::UniformVector<int,blaze::rowVector> vec( 4, 2 );

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );
      checkNonZeros( vec, 4UL );

      if( vec[0] != 2 || vec[1] != 2 || vec[2] != 2 || vec[3] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 2 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }

      reset( vec );

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );
      checkNonZeros( vec, 0UL );

      if( vec[0] != 0 || vec[1] != 0 || vec[2] != 0 || vec[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c clear() member function of the UniformVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c clear() member function of the UniformVector class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testClear()
{
   using blaze::clear;

   test_ = "UniformVector::clear()";

   // Clearing a default constructed vector
   {
      blaze::UniformVector<int,blaze::rowVector> vec;

      clear( vec );

      checkSize    ( vec, 0UL );
      checkNonZeros( vec, 0UL );
   }

   // Clearing an initialized vector
   {
      blaze::UniformVector<int,blaze::rowVector> vec( 4, 2 );

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );
      checkNonZeros( vec, 4UL );

      if( vec[0] != 2 || vec[1] != 2 || vec[2] != 2 || vec[3] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 2 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }

      clear( vec );

      checkSize    ( vec, 0UL );
      checkNonZeros( vec, 0UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c resize() member function of the UniformVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c resize() member function of the UniformVector class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testResize()
{
   test_ = "UniformVector::resize()";

   // Initialization check
   blaze::UniformVector<int,blaze::rowVector> vec;

   checkSize    ( vec, 0UL );
   checkNonZeros( vec, 0UL );

   // Resizing to 0
   vec.resize( 0UL );

   checkSize    ( vec, 0UL );
   checkNonZeros( vec, 0UL );

   // Resizing to 3
   vec.resize( 3UL );

   checkSize    ( vec, 3UL );
   checkCapacity( vec, 3UL );

   if( vec[0] != 0 || vec[1] != 0 || vec[2] != 0 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Resizing the vector failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 0 0 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Resizing to 5 and preserving the elements
   vec = 5;
   vec.resize( 5UL, true );

   checkSize    ( vec, 5UL );
   checkCapacity( vec, 5UL );

   if( vec[0] != 5 || vec[1] != 5 || vec[2] != 5 || vec[3] != 5 || vec[4] != 5 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Resizing the vector failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 5 5 5 5 5 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Resizing to 2 and preserving the elements
   vec.resize( 2UL, true );

   checkSize    ( vec, 2UL );
   checkCapacity( vec, 2UL );
   checkNonZeros( vec, 2UL );

   if( vec[0] != 5 || vec[1] != 5 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Resizing the vector failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 5 5 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Resizing to 0
   vec.resize( 0 );

   checkSize    ( vec, 0UL );
   checkNonZeros( vec, 0UL );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c extend() member function of the UniformVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c extend() member function of the UniformVector class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testExtend()
{
   test_ = "UniformVector::extend()";

   // Initialization check
   blaze::UniformVector<int,blaze::rowVector> vec;

   checkSize    ( vec, 0UL );
   checkNonZeros( vec, 0UL );

   // Increasing the size of the vector
   vec.extend( 3UL );

   checkSize    ( vec, 3UL );
   checkCapacity( vec, 3UL );
   checkNonZeros( vec, 0UL );

   if( vec[0] != 0 || vec[1] != 0 || vec[2] != 0 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Extending the vector failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 0 0 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Further increasing the size of the vector
   vec = 4;
   vec.extend( 2UL );

   checkSize    ( vec, 5UL );
   checkCapacity( vec, 5UL );

   if( vec[0] != 4 || vec[1] != 4 || vec[2] != 4 || vec[3] != 4 || vec[4] != 4 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Extending the vector failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 4 4 4 4 4 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Further increasing the size of the vector
   vec.extend( 10UL );

   checkSize    ( vec, 15UL );
   checkCapacity( vec, 15UL );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c swap() functionality of the UniformVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c swap() function of the UniformVector class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testSwap()
{
   test_ = "UniformVector swap";

   blaze::UniformVector<int,blaze::rowVector> vec1( 3, 2 );
   blaze::UniformVector<int,blaze::rowVector> vec2( 4, 5 );

   swap( vec1, vec2 );

   checkSize    ( vec1, 4UL );
   checkCapacity( vec1, 4UL );
   checkNonZeros( vec1, 4UL );

   if( vec1[0] != 5 || vec1[1] != 5 || vec1[2] != 5 || vec1[3] != 5 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Swapping the first vector failed\n"
          << " Details:\n"
          << "   Result:\n" << vec1 << "\n"
          << "   Expected result:\n( 5 5 5 5 )\n";
      throw std::runtime_error( oss.str() );
   }

   checkSize    ( vec2, 3UL );
   checkCapacity( vec2, 3UL );
   checkNonZeros( vec2, 3UL );

   if( vec2[0] != 2 || vec2[1] != 2 || vec2[2] != 2 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Swapping the second vector failed\n"
          << " Details:\n"
          << "   Result:\n" << vec1 << "\n"
          << "   Expected result:\n( 2 2 2 )\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c isDefault() function with the UniformVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isDefault() function with the UniformVector class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testIsDefault()
{
   using blaze::isDefault;

   test_ = "isDefault() function";

   // isDefault with vector of size 0
   {
      blaze::UniformVector<int,blaze::rowVector> vec;

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
      blaze::UniformVector<int,blaze::rowVector> vec( 3, 0 );

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
      blaze::UniformVector<int,blaze::rowVector> vec( 3, 1 );

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

} // namespace uniformvector

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
   std::cout << "   Running UniformVector class test..." << std::endl;

   try
   {
      RUN_UNIFORMVECTOR_CLASS_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during UniformVector class test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
