//=================================================================================================
/*!
//  \file src/mathtest/views/subvector/DenseUnalignedTest.cpp
//  \brief Source file for the Subvector dense unaligned test
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
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/Views.h>
#include <blaze/util/Memory.h>
#include <blaze/util/policies/Deallocate.h>
#include <blazetest/mathtest/views/subvector/DenseUnalignedTest.h>

#ifdef BLAZE_USE_HPX_THREADS
#  include <hpx/hpx_main.hpp>
#endif


namespace blazetest {

namespace mathtest {

namespace views {

namespace subvector {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the Subvector dense unaligned test.
//
// \exception std::runtime_error Operation error detected.
*/
DenseUnalignedTest::DenseUnalignedTest()
   : vec_( 8UL )
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
   testIterator();
   testNonZeros();
   testReset();
   testClear();
   testIsDefault();
   testIsSame();
   testSubvector();
   testElements();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the Subvector constructors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all constructors of the Subvector specialization. In case
// an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseUnalignedTest::testConstructors()
{
   test_ = "Subvector constructor";

   initialize();

   for( size_t start=0UL; start<vec_.size(); ++start ) {
      for( size_t size=1UL; start+size<vec_.size(); ++size )
      {
         SVT sv = blaze::subvector( vec_, start, size );

         for( size_t i=0UL; i<size; ++i )
         {
            if( sv[i] != vec_[start+i] ) {
               std::ostringstream oss;
               oss << " Test: " << test_ << "\n"
                   << " Error: Setup of dense subvector failed\n"
                   << " Details:\n"
                   << "   Start = " << start << "\n"
                   << "   Size  = " << size  << "\n"
                   << "   Subvector:\n" << sv << "\n"
                   << "   Vector:\n" << vec_ << "\n";
               throw std::runtime_error( oss.str() );
            }
         }
      }
   }

   try {
      SVT sv = blaze::subvector( vec_, 2UL, 7UL );

      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Setup of out-of-bounds subvector succeeded\n"
          << " Details:\n"
          << "   Result:\n" << sv << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::invalid_argument& ) {}

   try {
      SVT sv = blaze::subvector( vec_, 9UL, 0UL );

      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Setup of out-of-bounds subvector succeeded\n"
          << " Details:\n"
          << "   Result:\n" << sv << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::invalid_argument& ) {}
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Subvector assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all assignment operators of the Subvector specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseUnalignedTest::testAssignment()
{
   using blaze::aligned;
   using blaze::unaligned;
   using blaze::padded;
   using blaze::unpadded;
   using blaze::rowVector;


   //=====================================================================================
   // Homogeneous assignment
   //=====================================================================================

   {
      test_ = "Subvector homogeneous assignment";

      initialize();

      SVT sv = blaze::subvector( vec_, 2UL, 4UL );
      sv = 12;

      checkSize    ( sv  ,  4UL );
      checkNonZeros( sv  ,  4UL );
      checkSize    ( vec_,  8UL );
      checkNonZeros( vec_,  6UL );

      if( sv[0] != 12 || sv[1] != 12 || sv[2] != 12 || sv[3] != 12 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 12 12 12 12 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] !=  1 || vec_[2] != 12 || vec_[3] != 12 ||
          vec_[4] != 12 || vec_[5] != 12 || vec_[6] !=  4 || vec_[7] !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 1 12 12 12 12 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // List assignment
   //=====================================================================================

   {
      test_ = "Subvector initializer list assignment (complete list)";

      initialize();

      SVT sv = blaze::subvector( vec_, 2UL, 4UL );
      sv = { 1, 2, 3, 4 };

      checkSize    ( sv  ,  4UL );
      checkNonZeros( sv  ,  4UL );
      checkSize    ( vec_,  8UL );
      checkNonZeros( vec_,  6UL );

      if( sv[0] != 1 || sv[1] != 2 || sv[2] != 3 || sv[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 1 2 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] != 0 || vec_[1] != 1 || vec_[2] != 1 || vec_[3] != 2 ||
          vec_[4] != 3 || vec_[5] != 4 || vec_[6] != 4 || vec_[7] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 1 1 2 3 4 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Subvector initializer list assignment (incomplete list)";

      initialize();

      SVT sv = blaze::subvector( vec_, 2UL, 4UL );
      sv = { 1, 2 };

      checkSize    ( sv  ,  4UL );
      checkNonZeros( sv  ,  2UL );
      checkSize    ( vec_,  8UL );
      checkNonZeros( vec_,  4UL );

      if( sv[0] != 1 || sv[1] != 2 || sv[2] != 0 || sv[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 1 2 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] != 0 || vec_[1] != 1 || vec_[2] != 1 || vec_[3] != 2 ||
          vec_[4] != 0 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 1 1 2 0 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Copy assignment
   //=====================================================================================

   {
      test_ = "Subvector copy assignment (no aliasing)";

      initialize();

      VT vec( 10UL, 0 );
      vec[5] =  6;
      vec[6] = -8;

      SVT sv = blaze::subvector( vec, 5UL, 3UL );
      sv = blaze::subvector( vec_, 4UL, 3UL );

      checkSize    ( sv  ,  3UL );
      checkNonZeros( sv  ,  2UL );
      checkSize    ( vec_,  8UL );
      checkNonZeros( vec_,  4UL );
      checkSize    ( vec , 10UL );
      checkNonZeros( vec ,  2UL );

      if( sv[0] != -3 || sv[1] != 0 || sv[2] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( -3 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec[0] !=  0 || vec[1] != 0 || vec[2] != 0 || vec[3] != 0 || vec[4] != 0 ||
          vec[5] != -3 || vec[6] != 0 || vec[7] != 4 || vec[8] != 0 || vec[9] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 0 0 0 0 -3 0 4 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Subvector copy assignment (aliasing)";

      initialize();

      SVT sv = blaze::subvector( vec_, 1UL, 3UL );
      sv = blaze::subvector( vec_, 4UL, 3UL );

      checkSize    ( sv  , 3UL );
      checkNonZeros( sv  , 2UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 4UL );

      if( sv[0] != -3 || sv[1] != 0 || sv[2] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( -3 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != -3 || vec_[2] != 0 || vec_[3] != 4 ||
          vec_[4] != -3 || vec_[5] !=  0 || vec_[6] != 4 || vec_[7] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 -3 0 4 -3 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Dense vector assignment
   //=====================================================================================

   {
      test_ = "Dense vector assignment (mixed type)";

      initialize();

      SVT sv = blaze::subvector( vec_, 3UL, 4UL );

      const blaze::DynamicVector<short,rowVector> vec{ 0, 8, 0, 9 };

      sv = vec;

      checkSize    ( sv  , 4UL );
      checkNonZeros( sv  , 2UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 3UL );

      if( sv != vec ||
          sv[0] != 0 || sv[1] != 8 || sv[2] != 0 || sv[3] != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 0 8 0 9 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] != 0 || vec_[1] != 1 || vec_[2] != 0 || vec_[3] != 0 ||
          vec_[4] != 8 || vec_[5] != 0 || vec_[6] != 9 || vec_[7] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 1 0 0 8 0 9 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Dense vector assignment (aligned/padded)";

      initialize();

      SVT sv = blaze::subvector( vec_, 3UL, 4UL );

      using AlignedPadded = blaze::CustomVector<int,aligned,padded,rowVector>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 16UL ) );
      AlignedPadded vec( memory.get(), 4UL, 16UL );
      vec[0] = 0;
      vec[1] = 8;
      vec[2] = 0;
      vec[3] = 9;

      sv = vec;

      checkSize    ( sv  , 4UL );
      checkNonZeros( sv  , 2UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 3UL );

      if( sv != vec ||
          sv[0] != 0 || sv[1] != 8 || sv[2] != 0 || sv[3] != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 0 8 0 9 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] != 0 || vec_[1] != 1 || vec_[2] != 0 || vec_[3] != 0 ||
          vec_[4] != 8 || vec_[5] != 0 || vec_[6] != 9 || vec_[7] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 1 0 0 8 0 9 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Dense vector assignment (unaligned/unpadded)";

      initialize();

      SVT sv = blaze::subvector( vec_, 3UL, 4UL );

      using UnalignedUnpadded = blaze::CustomVector<int,unaligned,unpadded,rowVector>;
      std::unique_ptr<int[]> memory( new int[5] );
      UnalignedUnpadded vec( memory.get()+1UL, 4UL );
      vec[0] = 0;
      vec[1] = 8;
      vec[2] = 0;
      vec[3] = 9;

      sv = vec;

      checkSize    ( sv  , 4UL );
      checkNonZeros( sv  , 2UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 3UL );

      if( sv != vec ||
          sv[0] != 0 || sv[1] != 8 || sv[2] != 0 || sv[3] != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 0 8 0 9 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] != 0 || vec_[1] != 1 || vec_[2] != 0 || vec_[3] != 0 ||
          vec_[4] != 8 || vec_[5] != 0 || vec_[6] != 9 || vec_[7] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 1 0 0 8 0 9 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Sparse vector assignment
   //=====================================================================================

   {
      test_ = "Sparse vector assignment";

      initialize();

      SVT sv = blaze::subvector( vec_, 3UL, 4UL );

      blaze::CompressedVector<int,rowVector> vec( 4UL, 1UL );
      vec[3] = 9;

      sv = vec;

      checkSize    ( sv  , 4UL );
      checkNonZeros( sv  , 1UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 2UL );

      if( sv != vec ||
          sv[0] != 0 || sv[1] != 0 || sv[2] != 0 || sv[3] != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 0 0 0 9 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] != 0 || vec_[1] != 1 || vec_[2] != 0 || vec_[3] != 0 ||
          vec_[4] != 0 || vec_[5] != 0 || vec_[6] != 9 || vec_[7] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 1 0 0 0 0 9 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Subvector addition assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the addition assignment operators of the Subvector
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseUnalignedTest::testAddAssign()
{
   using blaze::aligned;
   using blaze::unaligned;
   using blaze::padded;
   using blaze::unpadded;
   using blaze::rowVector;


   //=====================================================================================
   // Subvector addition assignment
   //=====================================================================================

   {
      test_ = "Subvector addition assignment (no aliasing)";

      initialize();

      VT vec( 10UL, 0 );
      vec[5] =  6;
      vec[6] = -8;

      SVT sv = blaze::subvector( vec, 5UL, 3UL );
      sv += blaze::subvector( vec_, 4UL, 3UL );

      checkSize    ( sv  ,  3UL );
      checkNonZeros( sv  ,  3UL );
      checkSize    ( vec_,  8UL );
      checkNonZeros( vec_,  4UL );
      checkSize    ( vec , 10UL );
      checkNonZeros( vec ,  3UL );

      if( sv[0] != 3 || sv[1] != -8 || sv[2] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 3 -8 4 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec[0] != 0 || vec[1] !=  0 || vec[2] != 0 || vec[3] != 0 || vec[4] != 0 ||
          vec[5] != 3 || vec[6] != -8 || vec[7] != 4 || vec[8] != 0 || vec[9] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 0 0 0 0 3 -8 4 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Subvector addition assignment (aliasing)";

      initialize();

      SVT sv = blaze::subvector( vec_, 1UL, 3UL );
      sv += blaze::subvector( vec_, 3UL, 3UL );

      checkSize    ( sv  , 3UL );
      checkNonZeros( sv  , 3UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 5UL );

      if( sv[0] != -1 || sv[1] != -3 || sv[2] != -2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( -1 -3 -2 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != -1 || vec_[2] != -3 || vec_[3] != -2 ||
          vec_[4] != -3 || vec_[5] !=  0 || vec_[6] !=  4 || vec_[7] !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 -1 -3 -2 -3 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Dense vector addition assignment
   //=====================================================================================

   {
      test_ = "Dense vector addition assignment (mixed type)";

      initialize();

      SVT sv = blaze::subvector( vec_, 1UL, 3UL );

      const blaze::DynamicVector<short,rowVector> vec{ 2, -4, 0 };

      sv += vec;

      checkSize    ( sv  , 3UL );
      checkNonZeros( sv  , 3UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 5UL );

      if( sv[0] != 3 || sv[1] != -4 || sv[2] != -2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 3 -4 -2 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != 3 || vec_[2] != -4 || vec_[3] != -2 ||
          vec_[4] != -3 || vec_[5] != 0 || vec_[6] !=  4 || vec_[7] !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 3 -4 -2 -3 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Dense vector addition assignment (aligned/padded)";

      initialize();

      SVT sv = blaze::subvector( vec_, 1UL, 3UL );

      using AlignedPadded = blaze::CustomVector<int,aligned,padded,rowVector>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 16UL ) );
      AlignedPadded vec( memory.get(), 3UL, 16UL );
      vec[0] =  2;
      vec[1] = -4;
      vec[2] =  0;

      sv += vec;

      checkSize    ( sv  , 3UL );
      checkNonZeros( sv  , 3UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 5UL );

      if( sv[0] != 3 || sv[1] != -4 || sv[2] != -2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 3 -4 -2 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != 3 || vec_[2] != -4 || vec_[3] != -2 ||
          vec_[4] != -3 || vec_[5] != 0 || vec_[6] !=  4 || vec_[7] !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 3 -4 -2 -3 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Dense vector addition assignment (unaligned/unpadded)";

      initialize();

      SVT sv = blaze::subvector( vec_, 1UL, 3UL );

      using UnalignedUnpadded = blaze::CustomVector<int,unaligned,unpadded,rowVector>;
      std::unique_ptr<int[]> memory( new int[4] );
      UnalignedUnpadded vec( memory.get()+1UL, 3UL );
      vec[0] =  2;
      vec[1] = -4;
      vec[2] =  0;

      sv += vec;

      checkSize    ( sv  , 3UL );
      checkNonZeros( sv  , 3UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 5UL );

      if( sv[0] != 3 || sv[1] != -4 || sv[2] != -2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 3 -4 -2 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != 3 || vec_[2] != -4 || vec_[3] != -2 ||
          vec_[4] != -3 || vec_[5] != 0 || vec_[6] !=  4 || vec_[7] !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 3 -4 -2 -3 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Sparse vector addition assignment
   //=====================================================================================

   {
      test_ = "Sparse vector addition assignment";

      initialize();

      SVT sv = blaze::subvector( vec_, 1UL, 3UL );

      blaze::CompressedVector<int,rowVector> vec( 3UL, 2UL );
      vec[0] =  2;
      vec[1] = -4;

      sv += vec;

      checkSize    ( sv  , 3UL );
      checkNonZeros( sv  , 3UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 5UL );

      if( sv[0] != 3 || sv[1] != -4 || sv[2] != -2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 3 -4 -2 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != 3 || vec_[2] != -4 || vec_[3] != -2 ||
          vec_[4] != -3 || vec_[5] != 0 || vec_[6] !=  4 || vec_[7] !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 3 -4 -2 -3 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Subvector subtraction assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the subtraction assignment operators of the Subvector
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseUnalignedTest::testSubAssign()
{
   using blaze::aligned;
   using blaze::unaligned;
   using blaze::padded;
   using blaze::unpadded;
   using blaze::rowVector;


   //=====================================================================================
   // Subvector subtraction assignment
   //=====================================================================================

   {
      test_ = "Subvector subtraction assignment (no aliasing)";

      initialize();

      VT vec( 10UL, 0 );
      vec[5] =  6;
      vec[6] = -8;

      SVT sv = blaze::subvector( vec, 5UL, 3UL );
      sv -= blaze::subvector( vec_, 4UL, 3UL );

      checkSize    ( sv  ,  3UL );
      checkNonZeros( sv  ,  3UL );
      checkSize    ( vec_,  8UL );
      checkNonZeros( vec_,  4UL );
      checkSize    ( vec , 10UL );
      checkNonZeros( vec ,  3UL );

      if( sv[0] != 9 || sv[1] != -8 || sv[2] != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 9 -8 -4 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec[0] != 0 || vec[1] !=  0 || vec[2] !=  0 || vec[3] != 0 || vec[4] != 0 ||
          vec[5] != 9 || vec[6] != -8 || vec[7] != -4 || vec[8] != 0 || vec[9] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 0 0 0 0 9 -8 -4 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Subvector subtraction assignment (aliasing)";

      initialize();

      SVT sv = blaze::subvector( vec_, 1UL, 3UL );
      sv -= blaze::subvector( vec_, 3UL, 3UL );

      checkSize    ( sv  , 3UL );
      checkNonZeros( sv  , 3UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 5UL );

      if( sv[0] != 3 || sv[1] != 3 || sv[2] != -2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 3 3 -2 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != 3 || vec_[2] != 3 || vec_[3] != -2 ||
          vec_[4] != -3 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 3 3 -2 -3 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Dense vector subtraction assignment
   //=====================================================================================

   {
      test_ = "Dense vector subtraction assignment (mixed type)";

      initialize();

      SVT sv = blaze::subvector( vec_, 1UL, 3UL );

      const blaze::DynamicVector<short,rowVector> vec{ 2, -4, 0 };

      sv -= vec;

      checkSize    ( sv  , 3UL );
      checkNonZeros( sv  , 3UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 5UL );

      if( sv[0] != -1 || sv[1] != 4 || sv[2] != -2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( -1 4 -2 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != -1 || vec_[2] != 4 || vec_[3] != -2 ||
          vec_[4] != -3 || vec_[5] !=  0 || vec_[6] != 4 || vec_[7] !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 -1 4 -2 -3 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Dense vector subtraction assignment (aligned/padded)";

      initialize();

      SVT sv = blaze::subvector( vec_, 1UL, 3UL );

      using AlignedPadded = blaze::CustomVector<int,aligned,padded,rowVector>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 16UL ) );
      AlignedPadded vec( memory.get(), 3UL, 16UL );
      vec[0] =  2;
      vec[1] = -4;
      vec[2] =  0;

      sv -= vec;

      checkSize    ( sv  , 3UL );
      checkNonZeros( sv  , 3UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 5UL );

      if( sv[0] != -1 || sv[1] != 4 || sv[2] != -2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( -1 4 -2 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != -1 || vec_[2] != 4 || vec_[3] != -2 ||
          vec_[4] != -3 || vec_[5] !=  0 || vec_[6] != 4 || vec_[7] !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 -1 4 -2 -3 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Dense vector subtraction assignment (unaligned/unpadded)";

      initialize();

      SVT sv = blaze::subvector( vec_, 1UL, 3UL );

      using UnalignedUnpadded = blaze::CustomVector<int,unaligned,unpadded,rowVector>;
      std::unique_ptr<int[]> memory( new int[4] );
      UnalignedUnpadded vec( memory.get()+1UL, 3UL );
      vec[0] =  2;
      vec[1] = -4;
      vec[2] =  0;

      sv -= vec;

      checkSize    ( sv  , 3UL );
      checkNonZeros( sv  , 3UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 5UL );

      if( sv[0] != -1 || sv[1] != 4 || sv[2] != -2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( -1 4 -2 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != -1 || vec_[2] != 4 || vec_[3] != -2 ||
          vec_[4] != -3 || vec_[5] !=  0 || vec_[6] != 4 || vec_[7] !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 -1 4 -2 -3 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Sparse vector subtraction assignment
   //=====================================================================================

   {
      test_ = "Sparse vector subtraction assignment";

      initialize();

      SVT sv = blaze::subvector( vec_, 1UL, 3UL );

      blaze::CompressedVector<int,rowVector> vec( 3UL, 2UL );
      vec[0] =  2;
      vec[1] = -4;

      sv -= vec;

      checkSize    ( sv  , 3UL );
      checkNonZeros( sv  , 3UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 5UL );

      if( sv[0] != -1 || sv[1] != 4 || sv[2] != -2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( -1 4 -2 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != -1 || vec_[2] != 4 || vec_[3] != -2 ||
          vec_[4] != -3 || vec_[5] !=  0 || vec_[6] != 4 || vec_[7] !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 -1 4 -2 -3 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Subvector multiplication assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the multiplication assignment operators of the Subvector
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseUnalignedTest::testMultAssign()
{
   using blaze::aligned;
   using blaze::unaligned;
   using blaze::padded;
   using blaze::unpadded;
   using blaze::rowVector;


   //=====================================================================================
   // Subvector multiplication assignment
   //=====================================================================================

   {
      test_ = "Subvector multiplication assignment (no aliasing)";

      initialize();

      VT vec( 10UL, 0 );
      vec[5] =  6;
      vec[6] = -8;

      SVT sv = blaze::subvector( vec, 5UL, 3UL );
      sv *= blaze::subvector( vec_, 4UL, 3UL );

      checkSize    ( sv  ,  3UL );
      checkNonZeros( sv  ,  1UL );
      checkSize    ( vec_,  8UL );
      checkNonZeros( vec_,  4UL );
      checkSize    ( vec , 10UL );
      checkNonZeros( vec ,  1UL );

      if( sv[0] != -18 || sv[1] != 0 || sv[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( -18 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec[0] !=   0 || vec[1] != 0 || vec[2] != 0 || vec[3] != 0 || vec[4] != 0 ||
          vec[5] != -18 || vec[6] != 0 || vec[7] != 0 || vec[8] != 0 || vec[9] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 0 0 0 0 -18 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Subvector multiplication assignment (aliasing)";

      initialize();

      SVT sv = blaze::subvector( vec_, 1UL, 3UL );
      sv *= blaze::subvector( vec_, 3UL, 3UL );

      checkSize    ( sv  , 3UL );
      checkNonZeros( sv  , 1UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 3UL );

      if( sv[0] != -2 || sv[1] != 0 || sv[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( -2 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != -2 || vec_[2] != 0 || vec_[3] != 0 ||
          vec_[4] != -3 || vec_[5] !=  0 || vec_[6] != 4 || vec_[7] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 -2 0 0 -3 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Dense vector multiplication assignment
   //=====================================================================================

   {
      test_ = "Dense vector multiplication assignment (aligned/padded)";

      initialize();

      SVT sv = blaze::subvector( vec_, 1UL, 3UL );

      const blaze::DynamicVector<short,rowVector> vec{ 2, -4, 0 };

      sv *= vec;

      checkSize    ( sv  , 3UL );
      checkNonZeros( sv  , 1UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 3UL );

      if( sv[0] != 2 || sv[1] != 0 || sv[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 2 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != 2 || vec_[2] != 0 || vec_[3] != 0 ||
          vec_[4] != -3 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 2 0 0 -3 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Dense vector multiplication assignment (aligned/padded)";

      initialize();

      SVT sv = blaze::subvector( vec_, 1UL, 3UL );

      using AlignedPadded = blaze::CustomVector<int,aligned,padded,rowVector>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 16UL ) );
      AlignedPadded vec( memory.get(), 3UL, 16UL );
      vec[0] =  2;
      vec[1] = -4;
      vec[2] =  0;

      sv *= vec;

      checkSize    ( sv  , 3UL );
      checkNonZeros( sv  , 1UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 3UL );

      if( sv[0] != 2 || sv[1] != 0 || sv[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 2 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != 2 || vec_[2] != 0 || vec_[3] != 0 ||
          vec_[4] != -3 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 2 0 0 -3 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Dense vector multiplication assignment (unaligned/unpadded)";

      initialize();

      SVT sv = blaze::subvector( vec_, 1UL, 3UL );

      using UnalignedUnpadded = blaze::CustomVector<int,unaligned,unpadded,rowVector>;
      std::unique_ptr<int[]> memory( new int[4] );
      UnalignedUnpadded vec( memory.get()+1UL, 3UL );
      vec[0] =  2;
      vec[1] = -4;
      vec[2] =  0;

      sv *= vec;

      checkSize    ( sv  , 3UL );
      checkNonZeros( sv  , 1UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 3UL );

      if( sv[0] != 2 || sv[1] != 0 || sv[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 2 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != 2 || vec_[2] != 0 || vec_[3] != 0 ||
          vec_[4] != -3 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 2 0 0 -3 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Sparse vector multiplication assignment
   //=====================================================================================

   {
      test_ = "Sparse vector multiplication assignment";

      initialize();

      SVT sv = blaze::subvector( vec_, 1UL, 3UL );

      blaze::CompressedVector<int,rowVector> vec( 3UL, 2UL );
      vec[0] =  2;
      vec[1] = -4;

      sv *= vec;

      checkSize    ( sv  , 3UL );
      checkNonZeros( sv  , 1UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 3UL );

      if( sv[0] != 2 || sv[1] != 0 || sv[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 2 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != 2 || vec_[2] != 0 || vec_[3] != 0 ||
          vec_[4] != -3 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 2 0 0 -3 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Subvector division assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the division assignment operators of the Subvector
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseUnalignedTest::testDivAssign()
{
   using blaze::aligned;
   using blaze::unaligned;
   using blaze::padded;
   using blaze::unpadded;
   using blaze::rowVector;


   //=====================================================================================
   // Subvector division assignment
   //=====================================================================================

   {
      test_ = "Subvector division assignment (no aliasing)";

      initialize();

      VT vec( 10UL, 0 );
      vec[5] =  4;
      vec[6] = -6;

      SVT sv = blaze::subvector( vec, 5UL, 2UL );
      sv /= blaze::subvector( vec_, 3UL, 2UL );

      checkSize    ( sv  ,  2UL );
      checkNonZeros( sv  ,  2UL );
      checkSize    ( vec_,  8UL );
      checkNonZeros( vec_,  4UL );
      checkSize    ( vec , 10UL );
      checkNonZeros( vec ,  2UL );

      if( sv[0] != -2 || sv[1] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( -2 2 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec[0] !=  0 || vec[1] != 0 || vec[2] != 0 || vec[3] != 0 || vec[4] != 0 ||
          vec[5] != -2 || vec[6] != 2 || vec[7] != 0 || vec[8] != 0 || vec[9] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 0 0 0 0 -2 2 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Subvector division assignment (aliasing)";

      initialize();

      SVT sv = blaze::subvector( vec_, 6UL, 2UL );
      sv /= blaze::subvector( vec_, 3UL, 2UL );

      checkSize    ( sv  , 2UL );
      checkNonZeros( sv  , 1UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 4UL );

      if( sv[0] != -2 || sv[1] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( -2 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != 1 || vec_[2] !=  0 || vec_[3] != -2 ||
          vec_[4] != -3 || vec_[5] != 0 || vec_[6] != -2 || vec_[7] !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 1 0 -2 -3 0 -2 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Dense vector division assignment
   //=====================================================================================

   {
      test_ = "Dense vector division assignment (aligned/padded)";

      initialize();

      SVT sv = blaze::subvector( vec_, 1UL, 3UL );

      const blaze::DynamicVector<short,rowVector> vec{ 1, -4, 2 };

      sv /= vec;

      checkSize    ( sv  , 3UL );
      checkNonZeros( sv  , 2UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 4UL );

      if( sv[0] != 1 || sv[1] != 0 || sv[2] != -1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 1 0 -1 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != 1 || vec_[2] != 0 || vec_[3] != -1 ||
          vec_[4] != -3 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 1 0 -1 -3 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Dense vector division assignment (aligned/padded)";

      initialize();

      SVT sv = blaze::subvector( vec_, 1UL, 3UL );

      using AlignedPadded = blaze::CustomVector<int,aligned,padded,rowVector>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 16UL ) );
      AlignedPadded vec( memory.get(), 3UL, 16UL );
      vec[0] =  1;
      vec[1] = -4;
      vec[2] =  2;

      sv /= vec;

      checkSize    ( sv  , 3UL );
      checkNonZeros( sv  , 2UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 4UL );

      if( sv[0] != 1 || sv[1] != 0 || sv[2] != -1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 1 0 -1 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != 1 || vec_[2] != 0 || vec_[3] != -1 ||
          vec_[4] != -3 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 1 0 -1 -3 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Dense vector division assignment (unaligned/unpadded)";

      initialize();

      SVT sv = blaze::subvector( vec_, 1UL, 3UL );

      using UnalignedUnpadded = blaze::CustomVector<int,unaligned,unpadded,rowVector>;
      std::unique_ptr<int[]> memory( new int[4] );
      UnalignedUnpadded vec( memory.get()+1UL, 3UL );
      vec[0] =  1;
      vec[1] = -4;
      vec[2] =  2;

      sv /= vec;

      checkSize    ( sv  , 3UL );
      checkNonZeros( sv  , 2UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 4UL );

      if( sv[0] != 1 || sv[1] != 0 || sv[2] != -1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 1 0 -1 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != 1 || vec_[2] != 0 || vec_[3] != -1 ||
          vec_[4] != -3 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 1 0 -1 -3 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Subvector cross product assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the cross product assignment operators of the Subvector
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseUnalignedTest::testCrossAssign()
{
   using blaze::aligned;
   using blaze::unaligned;
   using blaze::padded;
   using blaze::unpadded;
   using blaze::rowVector;


   //=====================================================================================
   // Subvector cross product assignment
   //=====================================================================================

   {
      test_ = "Subvector cross product assignment (no aliasing)";

      initialize();

      VT vec( 10UL, 0 );
      vec[4] =  2;
      vec[6] = -1;
      vec[7] =  4;

      SVT sv = blaze::subvector( vec, 4UL, 3UL );
      sv %= blaze::subvector( vec_, 1UL, 3UL );

      checkSize    ( sv  ,  3UL );
      checkNonZeros( sv  ,  1UL );
      checkSize    ( vec_,  8UL );
      checkNonZeros( vec_,  4UL );
      checkSize    ( vec , 10UL );
      checkNonZeros( vec ,  2UL );

      if( sv[0] != 0 || sv[1] != 3 || sv[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cross product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 0 3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec[0] != 0 || vec[1] != 0 || vec[2] != 0 || vec[3] != 0 || vec[4] != 0 ||
          vec[5] != 3 || vec[6] != 0 || vec[7] != 4 || vec[8] != 0 || vec[9] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cross product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 0 0 0 0 3 0 4 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Subvector cross product assignment (aliasing)";

      initialize();

      SVT sv = blaze::subvector( vec_, 1UL, 3UL );
      sv %= blaze::subvector( vec_, 3UL, 3UL );

      checkSize    ( sv  , 3UL );
      checkNonZeros( sv  , 3UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 5UL );

      if( sv[0] != -6 || sv[1] != 4 || sv[2] != -3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cross product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( -6 4 -3 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != -6 || vec_[2] != 4 || vec_[3] != -3 ||
          vec_[4] != -3 || vec_[5] !=  0 || vec_[6] != 4 || vec_[7] !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cross product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 -6 4 -3 -3 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Dense vector cross product assignment
   //=====================================================================================

   {
      test_ = "Dense vector cross product assignment (aligned/padded)";

      initialize();

      SVT sv = blaze::subvector( vec_, 1UL, 3UL );

      const blaze::DynamicVector<short,rowVector> vec{ -2, 0, 1 };

      sv %= vec;

      checkSize    ( sv  , 3UL );
      checkNonZeros( sv  , 1UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 3UL );

      if( sv[0] != 0 || sv[1] != 3 || sv[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cross product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 0 3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != 0 || vec_[2] != 3 || vec_[3] != 0 ||
          vec_[4] != -3 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cross product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 0 3 0 -3 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Dense vector cross product assignment (aligned/padded)";

      initialize();

      SVT sv = blaze::subvector( vec_, 1UL, 3UL );

      using AlignedPadded = blaze::CustomVector<int,aligned,padded,rowVector>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 16UL ) );
      AlignedPadded vec( memory.get(), 3UL, 16UL );
      vec[0] = -2;
      vec[1] =  0;
      vec[2] =  1;

      sv %= vec;

      checkSize    ( sv  , 3UL );
      checkNonZeros( sv  , 1UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 3UL );

      if( sv[0] != 0 || sv[1] != 3 || sv[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cross product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 0 3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != 0 || vec_[2] != 3 || vec_[3] != 0 ||
          vec_[4] != -3 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cross product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 0 3 0 -3 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Dense vector cross product assignment (unaligned/unpadded)";

      initialize();

      SVT sv = blaze::subvector( vec_, 1UL, 3UL );

      using UnalignedUnpadded = blaze::CustomVector<int,unaligned,unpadded,rowVector>;
      std::unique_ptr<int[]> memory( new int[4] );
      UnalignedUnpadded vec( memory.get()+1UL, 3UL );
      vec[0] = -2;
      vec[1] =  0;
      vec[2] =  1;

      sv %= vec;

      checkSize    ( sv  , 3UL );
      checkNonZeros( sv  , 1UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 3UL );

      if( sv[0] != 0 || sv[1] != 3 || sv[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cross product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 0 3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != 0 || vec_[2] != 3 || vec_[3] != 0 ||
          vec_[4] != -3 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cross product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 0 3 0 -3 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Sparse vector cross product assignment
   //=====================================================================================

   {
      test_ = "Sparse vector cross product assignment";

      initialize();

      SVT sv = blaze::subvector( vec_, 1UL, 3UL );

      blaze::CompressedVector<int,rowVector> vec( 3UL, 2UL );
      vec[0] = -2;
      vec[2] =  1;

      sv %= vec;

      checkSize    ( sv  , 3UL );
      checkNonZeros( sv  , 1UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 3UL );

      if( sv[0] != 0 || sv[1] != 3 || sv[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cross product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 0 3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != 0 || vec_[2] != 3 || vec_[3] != 0 ||
          vec_[4] != -3 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cross product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 0 3 0 -3 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of all Subvector (self-)scaling operations.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all available ways to scale an instance of the Subvector
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseUnalignedTest::testScaling()
{
   //=====================================================================================
   // Self-scaling (v*=s)
   //=====================================================================================

   {
      test_ = "Subvector self-scaling (v*=s)";

      initialize();

      SVT sv = blaze::subvector( vec_, 1UL, 3UL );

      sv *= 3;

      checkSize    ( sv  , 3UL );
      checkNonZeros( sv  , 2UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 4UL );

      if( sv[0] != 3 || sv[1] != 0 || sv[2] != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 3 0 -6 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != 3 || vec_[2] != 0 || vec_[3] != -6 ||
          vec_[4] != -3 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 3 0 -6 -3 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Self-scaling (v=v*s)
   //=====================================================================================

   {
      test_ = "Subvector self-scaling (v=v*s)";

      initialize();

      SVT sv = blaze::subvector( vec_, 1UL, 3UL );

      sv = sv * 3;

      checkSize    ( sv  , 3UL );
      checkNonZeros( sv  , 2UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 4UL );

      if( sv[0] != 3 || sv[1] != 0 || sv[2] != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 3 0 -6 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != 3 || vec_[2] != 0 || vec_[3] != -6 ||
          vec_[4] != -3 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 3 0 -6 -3 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Self-scaling (v=s*v)
   //=====================================================================================

   {
      test_ = "Subvector self-scaling (v=v*s)";

      initialize();

      SVT sv = blaze::subvector( vec_, 1UL, 3UL );

      sv = 3 * sv;

      checkSize    ( sv  , 3UL );
      checkNonZeros( sv  , 2UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 4UL );

      if( sv[0] != 3 || sv[1] != 0 || sv[2] != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 3 0 -6 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != 3 || vec_[2] != 0 || vec_[3] != -6 ||
          vec_[4] != -3 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 3 0 -6 -3 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Self-scaling (v/=s)
   //=====================================================================================

   {
      test_ = "Subvector self-scaling (v/=s)";

      initialize();

      SVT sv = blaze::subvector( vec_, 1UL, 3UL );

      sv /= 0.5;

      checkSize    ( sv  , 3UL );
      checkNonZeros( sv  , 2UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 4UL );

      if( sv[0] != 2 || sv[1] != 0 || sv[2] != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 2 0 -4 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != 2 || vec_[2] != 0 || vec_[3] != -4 ||
          vec_[4] != -3 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 2 0 -4 -3 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Self-scaling (v=v/s)
   //=====================================================================================

   {
      test_ = "Subvector self-scaling (v/=s)";

      initialize();

      SVT sv = blaze::subvector( vec_, 1UL, 3UL );

      sv = sv / 0.5;

      checkSize    ( sv  , 3UL );
      checkNonZeros( sv  , 2UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 4UL );

      if( sv[0] != 2 || sv[1] != 0 || sv[2] != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 2 0 -4 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != 2 || vec_[2] != 0 || vec_[3] != -4 ||
          vec_[4] != -3 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 2 0 -4 -3 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Subvector::scale()
   //=====================================================================================

   {
      test_ = "Subvector::scale()";

      initialize();

      SVT sv = blaze::subvector( vec_, 1UL, 4UL );

      // Integral scaling of the subvector in the range [1,4]
      sv.scale( 3 );

      checkSize    ( sv  , 4UL );
      checkNonZeros( sv  , 3UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 4UL );

      if( sv[0] != 3 || sv[1] != 0 || sv[2] != -6 || sv[3] != -9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Integral scale operation of range [1,4] failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 3 0 -6 -9 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != 3 || vec_[2] != 0 || vec_[3] != -6 ||
          vec_[4] != -9 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Integral scale operation of range [1,4] failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 3 0 -6 -9 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Floating point scaling of the subvector in the range [1,4]
      sv.scale( 0.5 );

      checkSize    ( sv  , 4UL );
      checkNonZeros( sv  , 3UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 4UL );

      if( sv[0] != 1 || sv[1] != 0 || sv[2] != -3 || sv[3] != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Floating point scale operation of range [1,4] failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 1 0 -3 -4 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != 1 || vec_[2] != 0 || vec_[3] != -3 ||
          vec_[4] != -4 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Floating point scale operation of range [1,4] failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 1 0 -3 -4 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Subvector subscript operator.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of adding and accessing elements via the subscript operator
// of the Subvector specialization. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
void DenseUnalignedTest::testSubscript()
{
   test_ = "Subvector::operator[]";

   initialize();

   SVT sv = blaze::subvector( vec_, 1UL, 4UL );

   // Assignment to the element at index 1
   sv[1] = 9;

   checkSize    ( sv  , 4UL );
   checkNonZeros( sv  , 4UL );
   checkSize    ( vec_, 8UL );
   checkNonZeros( vec_, 5UL );

   if( sv[0] != 1 || sv[1] != 9 || sv[2] != -2 || sv[3] != -3 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << sv << "\n"
          << "   Expected result:\n( 1 9 -2 -3 )\n";
      throw std::runtime_error( oss.str() );
   }

   if( vec_[0] !=  0 || vec_[1] != 1 || vec_[2] != 9 || vec_[3] != -2 ||
       vec_[4] != -3 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] !=  0 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << vec_ << "\n"
          << "   Expected result:\n( 0 1 9 -2 -3 0 4 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Assignment to the element at index 2
   sv[2] = 0;

   checkSize    ( sv  , 4UL );
   checkNonZeros( sv  , 3UL );
   checkSize    ( vec_, 8UL );
   checkNonZeros( vec_, 4UL );

   if( sv[0] != 1 || sv[1] != 9 || sv[2] != 0 || sv[3] != -3 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << sv << "\n"
          << "   Expected result:\n( 1 9 0 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   if( vec_[0] !=  0 || vec_[1] != 1 || vec_[2] != 9 || vec_[3] != 0 ||
       vec_[4] != -3 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] != 0 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << vec_ << "\n"
          << "   Expected result:\n( 0 1 9 0 -3 0 4 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Assignment to the element at index 8
   sv[3] = -8;

   checkSize    ( sv  , 4UL );
   checkNonZeros( sv  , 3UL );
   checkSize    ( vec_, 8UL );
   checkNonZeros( vec_, 4UL );

   if( sv[0] != 1 || sv[1] != 9 || sv[2] != 0 || sv[3] != -8 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << sv << "\n"
          << "   Expected result:\n( 1 9 0 -8 )\n";
      throw std::runtime_error( oss.str() );
   }

   if( vec_[0] !=  0 || vec_[1] != 1 || vec_[2] != 9 || vec_[3] != 0 ||
       vec_[4] != -8 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] != 0 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << vec_ << "\n"
          << "   Expected result:\n( 0 1 9 0 -8 0 4 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Addition assignment to the element at index 0
   sv[0] += -3;

   checkSize    ( sv  , 4UL );
   checkNonZeros( sv  , 3UL );
   checkSize    ( vec_, 8UL );
   checkNonZeros( vec_, 4UL );

   if( sv[0] != -2 || sv[1] != 9 || sv[2] != 0 || sv[3] != -8 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << sv << "\n"
          << "   Expected result:\n( -2 9 0 -8 )\n";
      throw std::runtime_error( oss.str() );
   }

   if( vec_[0] !=  0 || vec_[1] != -2 || vec_[2] != 9 || vec_[3] != 0 ||
       vec_[4] != -8 || vec_[5] !=  0 || vec_[6] != 4 || vec_[7] != 0 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << vec_ << "\n"
          << "   Expected result:\n( 0 -2 9 0 -8 0 4 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Subtraction assignment to the element at index 1
   sv[1] -= 6;

   checkSize    ( sv  , 4UL );
   checkNonZeros( sv  , 3UL );
   checkSize    ( vec_, 8UL );
   checkNonZeros( vec_, 4UL );

   if( sv[0] != -2 || sv[1] != 3 || sv[2] != 0 || sv[3] != -8 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << sv << "\n"
          << "   Expected result:\n( -2 3 0 -8 )\n";
      throw std::runtime_error( oss.str() );
   }

   if( vec_[0] !=  0 || vec_[1] != -2 || vec_[2] != 3 || vec_[3] != 0 ||
       vec_[4] != -8 || vec_[5] !=  0 || vec_[6] != 4 || vec_[7] != 0 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << vec_ << "\n"
          << "   Expected result:\n( 0 -2 3 0 -8 0 4 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Multiplication assignment to the element at index 1
   sv[1] *= -3;

   checkSize    ( sv  , 4UL );
   checkNonZeros( sv  , 3UL );
   checkSize    ( vec_, 8UL );
   checkNonZeros( vec_, 4UL );

   if( sv[0] != -2 || sv[1] != -9 || sv[2] != 0 || sv[3] != -8 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << sv << "\n"
          << "   Expected result:\n( -2 -9 0 -8 )\n";
      throw std::runtime_error( oss.str() );
   }

   if( vec_[0] !=  0 || vec_[1] != -2 || vec_[2] != -9 || vec_[3] != 0 ||
       vec_[4] != -8 || vec_[5] !=  0 || vec_[6] !=  4 || vec_[7] != 0 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << vec_ << "\n"
          << "   Expected result:\n( 0 -2 -9 0 -8 0 4 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Multiplication assignment to the element at index 3
   sv[3] /= 2;

   checkSize    ( sv  , 4UL );
   checkNonZeros( sv  , 3UL );
   checkSize    ( vec_, 8UL );
   checkNonZeros( vec_, 4UL );

   if( sv[0] != -2 || sv[1] != -9 || sv[2] != 0 || sv[3] != -4 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << sv << "\n"
          << "   Expected result:\n( -2 -9 0 -4 )\n";
      throw std::runtime_error( oss.str() );
   }

   if( vec_[0] !=  0 || vec_[1] != -2 || vec_[2] != -9 || vec_[3] != 0 ||
       vec_[4] != -4 || vec_[5] !=  0 || vec_[6] !=  4 || vec_[7] != 0 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << vec_ << "\n"
          << "   Expected result:\n( 0 -2 -9 0 -4 0 4 0 )\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Subvector iterator implementation.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the iterator implementation of the Subvector specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseUnalignedTest::testIterator()
{
   initialize();

   // Testing the Iterator default constructor
   {
      test_ = "Iterator default constructor";

      SVT::Iterator it{};

      if( it != SVT::Iterator() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed iterator default constructor\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing the ConstIterator default constructor
   {
      test_ = "ConstIterator default constructor";

      SVT::ConstIterator it{};

      if( it != SVT::ConstIterator() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed iterator default constructor\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing conversion from Iterator to ConstIterator
   {
      test_ = "Iterator/ConstIterator conversion";

      SVT sv = blaze::subvector( vec_, 1UL, 4UL );
      SVT::ConstIterator it( begin( sv ) );

      if( it == end( sv ) || *it != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed iterator conversion detected\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Counting the number of elements in first half of the vector via Iterator (end-begin)
   {
      test_ = "Iterator subtraction (end-begin)";

      SVT sv = blaze::subvector( vec_, 0UL, 5UL );
      const ptrdiff_t number( end( sv ) - begin( sv ) );

      if( number != 5L ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid number of elements detected\n"
             << " Details:\n"
             << "   Number of elements         : " << number << "\n"
             << "   Expected number of elements: 5\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Counting the number of elements in first half of the vector via Iterator (begin-end)
   {
      test_ = "Iterator subtraction (begin-end)";

      SVT sv = blaze::subvector( vec_, 0UL, 5UL );
      const ptrdiff_t number( begin( sv ) - end( sv ) );

      if( number != -5L ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid number of elements detected\n"
             << " Details:\n"
             << "   Number of elements         : " << number << "\n"
             << "   Expected number of elements: -5\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Counting the number of elements in second half of the vector via ConstIterator (end-begin)
   {
      test_ = "ConstIterator subtraction (end-begin)";

      SVT sv = blaze::subvector( vec_, 5UL, 3UL );
      const ptrdiff_t number( cend( sv ) - cbegin( sv ) );

      if( number != 3L ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid number of elements detected\n"
             << " Details:\n"
             << "   Number of elements         : " << number << "\n"
             << "   Expected number of elements: 3\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Counting the number of elements in second half of the vector via ConstIterator (begin-end)
   {
      test_ = "ConstIterator subtraction (begin-end)";

      SVT sv = blaze::subvector( vec_, 5UL, 3UL );
      const ptrdiff_t number( cbegin( sv ) - cend( sv ) );

      if( number != -3L ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid number of elements detected\n"
             << " Details:\n"
             << "   Number of elements         : " << number << "\n"
             << "   Expected number of elements: -3\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing read-only access via ConstIterator
   {
      test_ = "Read-only access via ConstIterator";

      SVT sv = blaze::subvector( vec_, 1UL, 4UL );
      SVT::ConstIterator it ( cbegin( sv ) );
      SVT::ConstIterator end( cend( sv ) );

      if( it == end || *it != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid initial iterator detected\n";
         throw std::runtime_error( oss.str() );
      }

      ++it;

      if( it == end || *it != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator pre-increment failed\n";
         throw std::runtime_error( oss.str() );
      }

      --it;

      if( it == end || *it != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator pre-decrement failed\n";
         throw std::runtime_error( oss.str() );
      }

      it++;

      if( it == end || *it != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator post-increment failed\n";
         throw std::runtime_error( oss.str() );
      }

      it--;

      if( it == end || *it != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator post-decrement failed\n";
         throw std::runtime_error( oss.str() );
      }

      it += 2UL;

      if( it == end || *it != -2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator addition assignment failed\n";
         throw std::runtime_error( oss.str() );
      }

      it -= 2UL;

      if( it == end || *it != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator subtraction assignment failed\n";
         throw std::runtime_error( oss.str() );
      }

      it = it + 3UL;

      if( it == end || *it != -3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator/scalar addition failed\n";
         throw std::runtime_error( oss.str() );
      }

      it = it - 3UL;

      if( it == end || *it != 1 ) {
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

   // Testing assignment via Iterator
   {
      test_ = "Assignment via Iterator";

      SVT sv = blaze::subvector( vec_, 2UL, 4UL );
      int value = 6;

      for( SVT::Iterator it=begin( sv ); it!=end( sv ); ++it ) {
         *it = value++;
      }

      if( sv[0] != 6 || sv[1] != 7 || sv[2] != 8 || sv[3] != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 6 7 8 9 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] != 0 || vec_[1] != 1 || vec_[2] != 6 || vec_[3] != 7 ||
          vec_[4] != 8 || vec_[5] != 9 || vec_[6] != 4 || vec_[7] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 1 6 7 8 9 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing addition assignment via Iterator
   {
      test_ = "Addition assignment via Iterator";

      SVT sv = blaze::subvector( vec_, 2UL, 4UL );
      int value = 2;

      for( SVT::Iterator it=begin( sv ); it!=end( sv ); ++it ) {
         *it += value++;
      }

      if( sv[0] != 8 || sv[1] != 10 || sv[2] != 12 || sv[3] != 14 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 8 10 12 14 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] !=  1 || vec_[2] != 8 || vec_[3] != 10 ||
          vec_[4] != 12 || vec_[5] != 14 || vec_[6] != 4 || vec_[7] !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 1 8 10 12 14 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing subtraction assignment via Iterator
   {
      test_ = "Subtraction assignment via Iterator";

      SVT sv = blaze::subvector( vec_, 2UL, 4UL );
      int value = 2;

      for( SVT::Iterator it=begin( sv ); it!=end( sv ); ++it ) {
         *it -= value++;
      }

      if( sv[0] != 6 || sv[1] != 7 || sv[2] != 8 || sv[3] != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 6 7 8 9 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] != 0 || vec_[1] != 1 || vec_[2] != 6 || vec_[3] != 7 ||
          vec_[4] != 8 || vec_[5] != 9 || vec_[6] != 4 || vec_[7] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 1 6 7 8 9 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing multiplication assignment via Iterator
   {
      test_ = "Multiplication assignment via Iterator";

      SVT sv = blaze::subvector( vec_, 2UL, 4UL );
      int value = 1;

      for( SVT::Iterator it=begin( sv ); it!=end( sv ); ++it ) {
         *it *= value++;
      }

      if( sv[0] != 6 || sv[1] != 14 || sv[2] != 24 || sv[3] != 36 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 6 14 24 36 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] !=  1 || vec_[2] != 6 || vec_[3] != 14 ||
          vec_[4] != 24 || vec_[5] != 36 || vec_[6] != 4 || vec_[7] !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 1 6 14 24 36 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing division assignment via Iterator
   {
      test_ = "Division assignment via Iterator";

      SVT sv = blaze::subvector( vec_, 2UL, 4UL );

      for( SVT::Iterator it=begin( sv ); it!=end( sv ); ++it ) {
         *it /= 2;
      }

      if( sv[0] != 3 || sv[1] != 7 || sv[2] != 12 || sv[3] != 18 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 3 7 12 18 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] !=  1 || vec_[2] != 3 || vec_[3] != 7 ||
          vec_[4] != 12 || vec_[5] != 18 || vec_[6] != 4 || vec_[7] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 1 0 3 7 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c nonZeros() member function of the Subvector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c nonZeros() member function of the Subvector
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseUnalignedTest::testNonZeros()
{
   test_ = "Subvector::nonZeros()";

   initialize();

   // Initialization check
   SVT sv = blaze::subvector( vec_, 0UL, 4UL );

   checkSize    ( sv  , 4UL );
   checkNonZeros( sv  , 2UL );
   checkSize    ( vec_, 8UL );
   checkNonZeros( vec_, 4UL );

   if( sv[0] != 0 || sv[1] != 1 || sv[2] != 0 || sv[3] != -2 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Initialization failed\n"
          << " Details:\n"
          << "   Result:\n" << sv << "\n"
          << "   Expected result:\n( 0 1 0 -2 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Changing the number of non-zeros via the dense subvector
   sv[3] = 0;

   checkSize    ( sv  , 4UL );
   checkNonZeros( sv  , 1UL );
   checkSize    ( vec_, 8UL );
   checkNonZeros( vec_, 3UL );

   if( sv[0] != 0 || sv[1] != 1 || sv[2] != 0 || sv[3] != 0 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << sv << "\n"
          << "   Expected result:\n( 0 1 0 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Changing the number of non-zeros via the dense vector
   vec_[2UL] = 5;

   checkSize    ( sv  , 4UL );
   checkNonZeros( sv  , 2UL );
   checkSize    ( vec_, 8UL );
   checkNonZeros( vec_, 4UL );

   if( sv[0] != 0 || sv[1] != 1 || sv[2] != 5 || sv[3] != 0 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << sv << "\n"
          << "   Expected result:\n( 0 1 5 0 )\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c reset() member function of the Subvector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c reset() member function of the Subvector
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseUnalignedTest::testReset()
{
   test_ = "Subvector::reset()";

   using blaze::reset;

   // Resetting a single element of the range [1,6]
   {
      initialize();

      SVT sv = blaze::subvector( vec_, 1UL, 6UL );
      reset( sv[2] );

      checkSize    ( sv  , 6UL );
      checkNonZeros( sv  , 3UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 3UL );

      if( sv[0] != 1 || sv[1] != 0 || sv[2] != 0 || sv[3] != -3 || sv[4] != 0 || sv[5] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 1 0 0 -3 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Resetting the range [0,3] (lvalue)
   {
      initialize();

      SVT sv = blaze::subvector( vec_, 0UL, 4UL );
      reset( sv );

      checkSize    ( sv  , 4UL );
      checkNonZeros( sv  , 0UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 2UL );

      if( sv[0] != 0 || sv[1] != 0 || sv[2] != 0 || sv[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation of range [0,3] failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != 0 || vec_[2] != 0 || vec_[3] != 0 ||
          vec_[4] != -3 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation of range [0,3] failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 0 0 0 -3 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Resetting the range [4,7] (rvalue)
   {
      initialize();

      reset( blaze::subvector( vec_, 4UL, 4UL ) );

      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 2UL );

      if( vec_[0] != 0 || vec_[1] != 1 || vec_[2] != 0 || vec_[3] != -2 ||
          vec_[4] != 0 || vec_[5] != 0 || vec_[6] != 0 || vec_[7] !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation of range [4,7] failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 1 0 -2 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c clear() function with the Subvector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c clear() function of the Subvector specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseUnalignedTest::testClear()
{
   test_ = "clear() function";

   using blaze::clear;

   // Clearing a single element of the range [1,6]
   {
      initialize();

      SVT sv = blaze::subvector( vec_, 1UL, 6UL );
      clear( sv[2] );

      checkSize    ( sv  , 6UL );
      checkNonZeros( sv  , 3UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 3UL );

      if( sv[0] != 1 || sv[1] != 0 || sv[2] != 0 || sv[3] != -3 || sv[4] != 0 || sv[5] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 1 0 0 -3 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Clearing the range [0,3] (lvalue)
   {
      initialize();

      SVT sv = blaze::subvector( vec_, 0UL, 4UL );
      clear( sv );

      checkSize    ( sv  , 4UL );
      checkNonZeros( sv  , 0UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 2UL );

      if( sv[0] != 0 || sv[1] != 0 || sv[2] != 0 || sv[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation of range [0,3] failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != 0 || vec_[2] != 0 || vec_[3] != 0 ||
          vec_[4] != -3 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation of range [0,3] failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 0 0 0 -3 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Clearing the range [4,7] (rvalue)
   {
      initialize();

      clear( blaze::subvector( vec_, 4UL, 4UL ) );

      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 2UL );

      if( vec_[0] != 0 || vec_[1] != 1 || vec_[2] != 0 || vec_[3] != -2 ||
          vec_[4] != 0 || vec_[5] != 0 || vec_[6] != 0 || vec_[7] !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation of range [4,7] failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 1 0 -2 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c isDefault() function with the Subvector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isDefault() function with the Subvector specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseUnalignedTest::testIsDefault()
{
   test_ = "isDefault() function";

   using blaze::isDefault;

   initialize();

   // isDefault with default vector
   {
      VT vec( 8UL, 0 );
      SVT sv = blaze::subvector( vec, 2UL, 5UL );

      if( isDefault( sv[1] ) != true ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isDefault evaluation\n"
             << " Details:\n"
             << "   Subvector element: " << sv[1] << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( isDefault( sv ) != true ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isDefault evaluation\n"
             << " Details:\n"
             << "   Subvector:\n" << sv << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // isDefault with non-default vector
   {
      SVT sv = blaze::subvector( vec_, 2UL, 5UL );

      if( isDefault( sv[1] ) != false ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isDefault evaluation\n"
             << " Details:\n"
             << "   Subvector element: " << sv[1] << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( isDefault( sv ) != false ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isDefault evaluation\n"
             << " Details:\n"
             << "   Subvector:\n" << sv << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c isSame() function with the Subvector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isSame() function with the Subvector specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseUnalignedTest::testIsSame()
{
   //=====================================================================================
   // Vector-based tests
   //=====================================================================================

   {
      test_ = "isSame() function (vector-based)";

      // isSame with vector and matching subvector
      {
         SVT sv = blaze::subvector( vec_, 0UL, 8UL );

         if( blaze::isSame( sv, vec_ ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Vector:\n" << vec_ << "\n"
                << "   Subvector:\n" << sv << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( vec_, sv ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Vector:\n" << vec_ << "\n"
                << "   Subvector:\n" << sv << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with vector and non-matching subvector (different size)
      {
         SVT sv = blaze::subvector( vec_, 0UL, 6UL );

         if( blaze::isSame( sv, vec_ ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Vector:\n" << vec_ << "\n"
                << "   Subvector:\n" << sv << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( vec_, sv ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Vector:\n" << vec_ << "\n"
                << "   Subvector:\n" << sv << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with vector and non-matching subvector (different offset)
      {
         SVT sv = blaze::subvector( vec_, 1UL, 7UL );

         if( blaze::isSame( sv, vec_ ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Vector:\n" << vec_ << "\n"
                << "   Subvector:\n" << sv << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( vec_, sv ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Vector:\n" << vec_ << "\n"
                << "   Subvector:\n" << sv << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with matching subvectors
      {
         SVT sv1 = blaze::subvector( vec_, 3UL, 4UL );
         SVT sv2 = blaze::subvector( vec_, 3UL, 4UL );

         if( blaze::isSame( sv1, sv2 ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First subvector:\n" << sv1 << "\n"
                << "   Second subvector:\n" << sv2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with non-matching subvectors (different size)
      {
         SVT sv1 = blaze::subvector( vec_, 3UL, 4UL );
         SVT sv2 = blaze::subvector( vec_, 3UL, 3UL );

         if( blaze::isSame( sv1, sv2 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First subvector:\n" << sv1 << "\n"
                << "   Second subvector:\n" << sv2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with non-matching subvectors (different offset)
      {
         SVT sv1 = blaze::subvector( vec_, 3UL, 4UL );
         SVT sv2 = blaze::subvector( vec_, 2UL, 4UL );

         if( blaze::isSame( sv1, sv2 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First subvector:\n" << sv1 << "\n"
                << "   Second subvector:\n" << sv2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Row-based tests
   //=====================================================================================

   {
      test_ = "isSame() function (row-based)";

      const blaze::DynamicMatrix<int,blaze::rowMajor> mat{ { 1, 2, 3 },
                                                           { 4, 5, 6 },
                                                           { 7, 8, 9 } };

      // isSame with row and matching subvector
      {
         auto r  = blaze::row( mat, 1UL );
         auto sv = blaze::subvector( r, 0UL, 3UL );

         if( blaze::isSame( sv, r ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Row:\n" << r << "\n"
                << "   Subvector:\n" << sv << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( r, sv ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Row:\n" << r << "\n"
                << "   Subvector:\n" << sv << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with row and non-matching subvector (different size)
      {
         auto r  = blaze::row( mat, 1UL );
         auto sv = blaze::subvector( r, 0UL, 2UL );

         if( blaze::isSame( sv, r ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Row:\n" << r << "\n"
                << "   Subvector:\n" << sv << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( r, sv ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Row:\n" << r << "\n"
                << "   Subvector:\n" << sv << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with row and non-matching subvector (different offset)
      {
         auto r  = blaze::row( mat, 1UL );
         auto sv = blaze::subvector( r, 1UL, 2UL );

         if( blaze::isSame( sv, r ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Row:\n" << r << "\n"
                << "   Subvector:\n" << sv << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( r, sv ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Row:\n" << r << "\n"
                << "   Subvector:\n" << sv << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with matching subvectors
      {
         auto r   = blaze::row( mat, 1UL );
         auto sv1 = blaze::subvector( r, 0UL, 2UL );
         auto sv2 = blaze::subvector( r, 0UL, 2UL );

         if( blaze::isSame( sv1, sv2 ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First subvector:\n" << sv1 << "\n"
                << "   Second subvector:\n" << sv2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with non-matching subvectors (different size)
      {
         auto r   = blaze::row( mat, 1UL );
         auto sv1 = blaze::subvector( r, 0UL, 2UL );
         auto sv2 = blaze::subvector( r, 0UL, 3UL );

         if( blaze::isSame( sv1, sv2 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First subvector:\n" << sv1 << "\n"
                << "   Second subvector:\n" << sv2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with non-matching subvectors (different offset)
      {
         auto r   = blaze::row( mat, 1UL );
         auto sv1 = blaze::subvector( r, 0UL, 2UL );
         auto sv2 = blaze::subvector( r, 1UL, 2UL );

         if( blaze::isSame( sv1, sv2 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First subvector:\n" << sv1 << "\n"
                << "   Second subvector:\n" << sv2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-based tests
   //=====================================================================================

   {
      test_ = "isSame() function (column-based)";

      const blaze::DynamicMatrix<int,blaze::columnMajor> mat{ { 1, 2, 3 },
                                                              { 4, 5, 6 },
                                                              { 7, 8, 9 } };

      // isSame with column and matching subvector
      {
         auto c  = blaze::column( mat, 1UL );
         auto sv = blaze::subvector( c, 0UL, 3UL );

         if( blaze::isSame( sv, c ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Column:\n" << c << "\n"
                << "   Subvector:\n" << sv << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( c, sv ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Column:\n" << c << "\n"
                << "   Subvector:\n" << sv << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with column and non-matching subvector (different size)
      {
         auto c  = blaze::column( mat, 1UL );
         auto sv = blaze::subvector( c, 0UL, 2UL );

         if( blaze::isSame( sv, c ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Column:\n" << c << "\n"
                << "   Subvector:\n" << sv << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( c, sv ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Column:\n" << c << "\n"
                << "   Subvector:\n" << sv << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with column and non-matching subvector (different offset)
      {
         auto c  = blaze::column( mat, 1UL );
         auto sv = blaze::subvector( c, 1UL, 2UL );

         if( blaze::isSame( sv, c ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Column:\n" << c << "\n"
                << "   Subvector:\n" << sv << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( c, sv ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Column:\n" << c << "\n"
                << "   Subvector:\n" << sv << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with matching subvectors
      {
         auto c   = blaze::column( mat, 1UL );
         auto sv1 = blaze::subvector( c, 0UL, 2UL );
         auto sv2 = blaze::subvector( c, 0UL, 2UL );

         if( blaze::isSame( sv1, sv2 ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First subvector:\n" << sv1 << "\n"
                << "   Second subvector:\n" << sv2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with non-matching subvectors (different size)
      {
         auto c   = blaze::column( mat, 1UL );
         auto sv1 = blaze::subvector( c, 0UL, 2UL );
         auto sv2 = blaze::subvector( c, 0UL, 3UL );

         if( blaze::isSame( sv1, sv2 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First subvector:\n" << sv1 << "\n"
                << "   Second subvector:\n" << sv2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with non-matching subvectors (different offset)
      {
         auto c   = blaze::column( mat, 1UL );
         auto sv1 = blaze::subvector( c, 0UL, 2UL );
         auto sv2 = blaze::subvector( c, 1UL, 2UL );

         if( blaze::isSame( sv1, sv2 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First subvector:\n" << sv1 << "\n"
                << "   Second subvector:\n" << sv2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c subvector() function with the Subvector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c subvector() function used with the Subvector
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseUnalignedTest::testSubvector()
{
   test_ = "subvector() function";

   initialize();

   {
      SVT sv1 = blaze::subvector( vec_, 1UL, 6UL );
      SVT sv2 = blaze::subvector( sv1 , 1UL, 4UL );

      if( sv2[0] != 0 || sv2[1] != -2 || sv2[2] != -3 || sv2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator access failed\n"
             << " Details:\n"
             << "   Result:\n" << sv2 << "\n"
             << "   Expected result:\n( 0 -2 -3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( *sv2.begin() != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << *sv2.begin() << "\n"
             << "   Expected result: 0\n";
         throw std::runtime_error( oss.str() );
      }
   }

   try {
      SVT sv1 = blaze::subvector( vec_, 1UL, 6UL );
      SVT sv2 = blaze::subvector( sv1 , 6UL, 2UL );

      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Setup of out-of-bounds subvector succeeded\n"
          << " Details:\n"
          << "   Result:\n" << sv2 << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::invalid_argument& ) {}

   try {
      SVT sv1 = blaze::subvector( vec_, 1UL, 6UL );
      SVT sv2 = blaze::subvector( sv1 , 2UL, 5UL );

      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Setup of out-of-bounds subvector succeeded\n"
          << " Details:\n"
          << "   Result:\n" << sv2 << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::invalid_argument& ) {}
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c elements() function with the Subvector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c elements() function used with the Subvector
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseUnalignedTest::testElements()
{
   //=====================================================================================
   // Setup via index_sequence
   //=====================================================================================

   {
      test_ = "elements() function (index_sequence)";

      initialize();

      {
         SVT sv = blaze::subvector( vec_, 1UL, 6UL );
         auto e = blaze::elements( sv, { 4UL, 3UL, 2UL, 1UL } );

         if( e[0] != 0 || e[1] != -3 || e[2] != -2 || e[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subscript operator access failed\n"
                << " Details:\n"
                << "   Result:\n" << e << "\n"
                << "   Expected result:\n( 0 -3 -2 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( *e.begin() != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator access failed\n"
                << " Details:\n"
                << "   Result: " << *e.begin() << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      try {
         SVT sv = blaze::subvector( vec_, 1UL, 6UL );
         auto e = blaze::elements( sv, { 6UL } );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of out-of-bounds element selection succeeded\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Setup via std::array
   //=====================================================================================

   {
      test_ = "elements() function (std::array)";

      initialize();

      {
         std::array<int,4UL> indices{ 4UL, 3UL, 2UL, 1UL };

         SVT sv = blaze::subvector( vec_, 1UL, 6UL );
         auto e = blaze::elements( sv, indices );

         if( e[0] != 0 || e[1] != -3 || e[2] != -2 || e[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subscript operator access failed\n"
                << " Details:\n"
                << "   Result:\n" << e << "\n"
                << "   Expected result:\n( 0 -3 -2 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( *e.begin() != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator access failed\n"
                << " Details:\n"
                << "   Result: " << *e.begin() << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      try {
         std::array<int,1UL> indices{ 6UL };

         SVT sv = blaze::subvector( vec_, 1UL, 6UL );
         auto e = blaze::elements( sv, indices );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of out-of-bounds element selection succeeded\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Setup via lambda expression
   //=====================================================================================

   {
      test_ = "elements() function (lambda expression)";

      initialize();

      {
         SVT sv = blaze::subvector( vec_, 1UL, 6UL );
         auto e = blaze::elements( sv, []( size_t i ){ return 4UL-i; }, 4UL );

         if( e[0] != 0 || e[1] != -3 || e[2] != -2 || e[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subscript operator access failed\n"
                << " Details:\n"
                << "   Result:\n" << e << "\n"
                << "   Expected result:\n( 0 -3 -2 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( *e.begin() != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator access failed\n"
                << " Details:\n"
                << "   Result: " << *e.begin() << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      try {
         SVT sv = blaze::subvector( vec_, 1UL, 6UL );
         auto e = blaze::elements( sv, []( size_t i ){ return i+6UL; }, 1UL );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of out-of-bounds element selection succeeded\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Initialization of all member vectors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function initializes all member vectors to specific predetermined values.
*/
void DenseUnalignedTest::initialize()
{
   // Initializing the dynamic row vector
   vec_[0] =  0;
   vec_[1] =  1;
   vec_[2] =  0;
   vec_[3] = -2;
   vec_[4] = -3;
   vec_[5] =  0;
   vec_[6] =  4;
   vec_[7] =  0;
}
//*************************************************************************************************

} // namespace subvector

} // namespace views

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
   std::cout << "   Running Subvector dense unaligned test..." << std::endl;

   try
   {
      RUN_SUBVECTOR_DENSEUNALIGNED_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during Subvector dense unaligned test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
