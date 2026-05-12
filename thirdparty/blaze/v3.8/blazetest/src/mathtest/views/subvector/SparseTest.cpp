//=================================================================================================
/*!
//  \file src/mathtest/views/subvector/SparseTest.cpp
//  \brief Source file for the Subvector sparse test
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
#include <blaze/math/Views.h>
#include <blazetest/mathtest/views/subvector/SparseTest.h>

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
/*!\brief Constructor for the Subvector sparse test.
//
// \exception std::runtime_error Operation error detected.
*/
SparseTest::SparseTest()
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
   testReserve();
   testSet();
   testInsert();
   testAppend();
   testErase();
   testFind();
   testLowerBound();
   testUpperBound();
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
void SparseTest::testConstructors()
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
                   << " Error: Setup of sparse subvector failed\n"
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
void SparseTest::testAssignment()
{
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
             << "   Expected result:\n( 1 2 3 4 )\n";
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

      VT vec( 10UL );
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
      test_ = "Dense vector assignment";

      initialize();

      SVT sv = blaze::subvector( vec_, 3UL, 4UL );

      blaze::DynamicVector<int,blaze::rowVector> vec{ 0, 8, 0, 9 };

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

      blaze::CompressedVector<int,blaze::rowVector> vec( 4UL );
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
void SparseTest::testAddAssign()
{
   //=====================================================================================
   // Subvector addition assignment
   //=====================================================================================

   {
      test_ = "Subvector addition assignment (no aliasing)";

      initialize();

      VT vec( 10UL );
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
      test_ = "Dense vector addition assignment";

      initialize();

      SVT sv = blaze::subvector( vec_, 1UL, 3UL );

      blaze::DynamicVector<int,blaze::rowVector> vec{ 2, -4, 0 };

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

      blaze::CompressedVector<int,blaze::rowVector> vec( 3UL );
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
void SparseTest::testSubAssign()
{
   //=====================================================================================
   // Subvector subtraction assignment
   //=====================================================================================

   {
      test_ = "Subvector subtraction assignment (no aliasing)";

      initialize();

      VT vec( 10UL );
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
      test_ = "Dense vector subtraction assignment";

      initialize();

      SVT sv = blaze::subvector( vec_, 1UL, 3UL );

      blaze::DynamicVector<int,blaze::rowVector> vec{ 2, -4, 0 };

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

      blaze::CompressedVector<int,blaze::rowVector> vec( 3UL );
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
void SparseTest::testMultAssign()
{
   //=====================================================================================
   // Subvector multiplication assignment
   //=====================================================================================

   {
      test_ = "Subvector multiplication assignment (no aliasing)";

      initialize();

      VT vec( 10UL );
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
      test_ = "Dense vector multiplication assignment";

      initialize();

      SVT sv = blaze::subvector( vec_, 1UL, 3UL );

      blaze::DynamicVector<int,blaze::rowVector> vec{ 2, -4, 0 };

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

      blaze::CompressedVector<int,blaze::rowVector> vec( 3UL, 0 );
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
void SparseTest::testDivAssign()
{
   //=====================================================================================
   // Dense vector division assignment
   //=====================================================================================

   {
      test_ = "Dense vector division assignment";

      initialize();

      SVT sv = blaze::subvector( vec_, 1UL, 3UL );

      blaze::DynamicVector<int,blaze::rowVector> vec{ 1, -4, 2 };

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
void SparseTest::testCrossAssign()
{
   //=====================================================================================
   // Subvector cross product assignment
   //=====================================================================================

   {
      test_ = "Subvector cross product assignment (no aliasing)";

      initialize();

      VT vec( 10UL, 3UL );
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
      test_ = "Dense vector cross product assignment";

      initialize();

      SVT sv = blaze::subvector( vec_, 1UL, 3UL );

      blaze::DynamicVector<int,blaze::rowVector> vec{ -2, 0, 1 };

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

      blaze::CompressedVector<int,blaze::rowVector> vec( 3UL, 0 );
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
void SparseTest::testScaling()
{
   //=====================================================================================
   // Self-scaling (v*=s)
   //=====================================================================================

   {
      test_ = "Self-scaling (v*=s)";

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
      test_ = "Self-scaling (v=v*s)";

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
      test_ = "Self-scaling (v=s*v)";

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
      test_ = "Self-scaling (v/=s)";

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
      test_ = "Self-scaling (v=v/s)";

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

      // Integral scaling the subvector in the range [1,4]
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

      // Floating point scaling the subvector in the range [1,4]
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
void SparseTest::testSubscript()
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

   // Assignment to the element at index 3
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
void SparseTest::testIterator()
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

      if( it == end( sv ) || it->value() != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed iterator conversion detected\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Counting the number of elements in first half of the vector via Iterator (end-begin)
   {
      test_ = "Iterator subtraction (end-begin)";

      SVT sv = blaze::subvector( vec_, 0UL, 4UL );
      const ptrdiff_t number( end( sv ) - begin( sv ) );

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

   // Counting the number of elements in second half of the vector via ConstIterator (end-begin)
   {
      test_ = "ConstIterator subtraction (end-begin)";

      SVT sv = blaze::subvector( vec_, 4UL, 4UL );
      const ptrdiff_t number( cend( sv ) - cbegin( sv ) );

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

      SVT sv = blaze::subvector( vec_, 1UL, 3UL );
      SVT::ConstIterator it ( cbegin( sv ) );
      SVT::ConstIterator end( cend( sv ) );

      if( it == end || it->value() != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid initial iterator detected\n";
         throw std::runtime_error( oss.str() );
      }

      ++it;

      if( it == end || it->value() != -2 ) {
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

      SVT sv = blaze::subvector( vec_, 2UL, 4UL );
      int value = 6;

      for( SVT::Iterator it=begin( sv ); it!=end( sv ); ++it ) {
         *it = value++;
      }

      if( sv[0] != 0 || sv[1] != 6 || sv[2] != 7 || sv[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 0 6 7 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] != 0 || vec_[1] != 1 || vec_[2] != 0 || vec_[3] != 6 ||
          vec_[4] != 7 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 1 0 6 7 0 4 0 )\n";
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

      if( sv[0] != 0 || sv[1] != 8 || sv[2] != 10 || sv[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 0 8 10 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != 1 || vec_[2] != 0 || vec_[3] != 8 ||
          vec_[4] != 10 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 1 0 8 10 0 4 0 )\n";
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

      if( sv[0] != 0 || sv[1] != 6 || sv[2] != 7 || sv[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 0 6 7 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] != 0 || vec_[1] != 1 || vec_[2] != 0 || vec_[3] != 6 ||
          vec_[4] != 7 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 1 0 6 7 0 4 0 )\n";
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

      if( sv[0] != 0 || sv[1] != 6 || sv[2] != 14 || sv[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 0 6 14 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != 1 || vec_[2] != 0 || vec_[3] != 6 ||
          vec_[4] != 14 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 1 0 6 14 0 4 0 )\n";
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

      if( sv[0] != 0 || sv[1] != 3 || sv[2] != 7 || sv[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 0 3 7 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] != 0 || vec_[1] != 1 || vec_[2] != 0 || vec_[3] != 3 ||
          vec_[4] != 7 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] != 0 ) {
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
void SparseTest::testNonZeros()
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

   // Changing the number of non-zeros via the sparse subvector
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

   // Changing the number of non-zeros via the sparse vector
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
void SparseTest::testReset()
{
   using blaze::reset;

   test_ = "Subvector::reset()";

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
void SparseTest::testClear()
{
   using blaze::clear;

   test_ = "clear() function";

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
/*!\brief Test of the \c reserve() member function of the Subvector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c reserve() member function of the Subvector
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testReserve()
{
   test_ = "Subvector::reserve()";

   VT vec( 10UL );

   SVT sv = blaze::subvector( vec, 2UL, 4UL );

   // Increasing the capacity of the vector
   sv.reserve( 10UL );

   checkSize    ( sv,  4UL );
   checkCapacity( sv, 10UL );
   checkNonZeros( sv,  0UL );

   // Further increasing the capacity of the vector
   sv.reserve( 20UL );

   checkSize    ( sv,  4UL );
   checkCapacity( sv, 20UL );
   checkNonZeros( sv,  0UL );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c set() member function of the Subvector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c set() member function of the Subvector specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testSet()
{
   test_ = "Subvector::set()";

   initialize();

   SVT sv = blaze::subvector( vec_, 0UL, 8UL );

   // Setting a non-zero element at the end of the subvector
   {
      SVT::Iterator pos = sv.set( 7UL, 9 );

      checkSize    ( sv  , 8UL );
      checkNonZeros( sv  , 5UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 5UL );

      if( pos->value() != 9 || pos->index() != 7UL ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid iterator returned\n"
             << " Details:\n"
             << "   Value: " << pos->value() << "\n"
             << "   Index: " << pos->index() << "\n"
             << "   Expected value: 9\n"
             << "   Expected index: 7\n";
         throw std::runtime_error( oss.str() );
      }

      if( sv[0] !=  0 || sv[1] != 1 || sv[2] != 0 || sv[3] != -2 ||
          sv[4] != -3 || sv[5] != 0 || sv[6] != 4 || sv[7] !=  9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setting a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 0 2 0 -2 -3 0 4 9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Setting a non-zero element at the beginning of the subvector
   {
      SVT::Iterator pos = sv.set( 0UL, 9 );

      checkSize    ( sv  , 8UL );
      checkNonZeros( sv  , 6UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 6UL );

      if( pos->value() != 9 || pos->index() != 0UL ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid iterator returned\n"
             << " Details:\n"
             << "   Value: " << pos->value() << "\n"
             << "   Index: " << pos->index() << "\n"
             << "   Expected value: 9\n"
             << "   Expected index: 0\n";
         throw std::runtime_error( oss.str() );
      }

      if( sv[0] !=  9 || sv[1] != 1 || sv[2] != 0 || sv[3] != -2 ||
          sv[4] != -3 || sv[5] != 0 || sv[6] != 4 || sv[7] !=  9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setting a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 9 2 0 -2 -3 0 4 9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Setting a non-zero element at the center of the subvector
   {
      SVT::Iterator pos = sv.set( 2UL, 9 );

      checkSize    ( sv  , 8UL );
      checkNonZeros( sv  , 7UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 7UL );

      if( pos->value() != 9 || pos->index() != 2UL ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid iterator returned\n"
             << " Details:\n"
             << "   Value: " << pos->value() << "\n"
             << "   Index: " << pos->index() << "\n"
             << "   Expected value: 9\n"
             << "   Expected index: 2\n";
         throw std::runtime_error( oss.str() );
      }

      if( sv[0] !=  9 || sv[1] != 1 || sv[2] != 9 || sv[3] != -2 ||
          sv[4] != -3 || sv[5] != 0 || sv[6] != 4 || sv[7] !=  9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setting a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 9 2 9 -2 -3 0 4 9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Setting an already existing element
   {
      SVT::Iterator pos = sv.set( 3UL, 9 );

      checkSize    ( sv  , 8UL );
      checkNonZeros( sv  , 7UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 7UL );

      if( pos->value() != 9 || pos->index() != 3UL ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid iterator returned\n"
             << " Details:\n"
             << "   Value: " << pos->value() << "\n"
             << "   Index: " << pos->index() << "\n"
             << "   Expected value: 9\n"
             << "   Expected index: 3\n";
         throw std::runtime_error( oss.str() );
      }

      if( sv[0] !=  9 || sv[1] != 1 || sv[2] != 9 || sv[3] != 9 ||
          sv[4] != -3 || sv[5] != 0 || sv[6] != 4 || sv[7] != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setting a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 9 2 9 9 -3 0 4 9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c insert() member function of the Subvector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c insert() member function of the Subvector
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testInsert()
{
   test_ = "Subvector::insert()";

   initialize();

   SVT sv = blaze::subvector( vec_, 0UL, 8UL );

   // Inserting a non-zero element at the end of the subvector
   {
      SVT::Iterator pos = sv.insert( 7UL, 9 );

      checkSize    ( sv  , 8UL );
      checkNonZeros( sv  , 5UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 5UL );

      if( pos->value() != 9 || pos->index() != 7UL ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid iterator returned\n"
             << " Details:\n"
             << "   Value: " << pos->value() << "\n"
             << "   Index: " << pos->index() << "\n"
             << "   Expected value: 9\n"
             << "   Expected index: 7\n";
         throw std::runtime_error( oss.str() );
      }

      if( sv[0] !=  0 || sv[1] != 1 || sv[2] != 0 || sv[3] != -2 ||
          sv[4] != -3 || sv[5] != 0 || sv[6] != 4 || sv[7] !=  9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Inserting a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 0 2 0 -2 -3 0 4 9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Inserting a non-zero element at the beginning of the subvector
   {
      SVT::Iterator pos = sv.insert( 0UL, 9 );

      checkSize    ( sv  , 8UL );
      checkNonZeros( sv  , 6UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 6UL );

      if( pos->value() != 9 || pos->index() != 0UL ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid iterator returned\n"
             << " Details:\n"
             << "   Value: " << pos->value() << "\n"
             << "   Index: " << pos->index() << "\n"
             << "   Expected value: 9\n"
             << "   Expected index: 0\n";
         throw std::runtime_error( oss.str() );
      }

      if( sv[0] !=  9 || sv[1] != 1 || sv[2] != 0 || sv[3] != -2 ||
          sv[4] != -3 || sv[5] != 0 || sv[6] != 4 || sv[7] !=  9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Inserting a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 9 2 0 -2 -3 0 4 9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Inserting a non-zero element at the center of the subvector
   {
      SVT::Iterator pos = sv.insert( 2UL, 9 );

      checkSize    ( sv  , 8UL );
      checkNonZeros( sv  , 7UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 7UL );

      if( pos->value() != 9 || pos->index() != 2UL ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid iterator returned\n"
             << " Details:\n"
             << "   Value: " << pos->value() << "\n"
             << "   Index: " << pos->index() << "\n"
             << "   Expected value: 9\n"
             << "   Expected index: 2\n";
         throw std::runtime_error( oss.str() );
      }

      if( sv[0] !=  9 || sv[1] != 1 || sv[2] != 9 || sv[3] != -2 ||
          sv[4] != -3 || sv[5] != 0 || sv[6] != 4 || sv[7] !=  9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Inserting a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 9 2 9 -2 -3 0 4 9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Trying to insert an already existing element
   try {
      sv.insert( 3UL, 9 );

      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Inserting an existing element succeeded\n"
          << " Details:\n"
          << "   Result:\n" << sv << "\n"
          << "   Expected result:\n( 9 2 0 9 -3 0 4 9 )\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::invalid_argument& ) {}
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c append() member function of the Subvector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c append() member function of the Subvector
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testAppend()
{
   test_ = "Subvector::append()";

   VT vec( 10UL );

   SVT sv = blaze::subvector( vec, 2UL, 4UL );
   sv.reserve( 4UL );

   // Appending one non-zero element
   sv.append( 0UL, 1 );

   checkSize    ( sv , 4UL );
   checkCapacity( sv , 4UL );
   checkNonZeros( sv , 1UL );
   checkNonZeros( vec, 1UL );

   if( sv[0] != 1 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Append operation failed\n"
          << " Details:\n"
          << "   Result:\n" << sv << "\n"
          << "   Expected result:\n( 0 0 1 0 0 0 0 0 0 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Appending three more non-zero elements
   sv.append( 1UL, 2 );
   sv.append( 2UL, 3 );
   sv.append( 3UL, 4 );

   checkSize    ( sv , 4UL );
   checkCapacity( sv , 4UL );
   checkNonZeros( sv , 4UL );
   checkNonZeros( vec, 4UL );

   if( sv[0] != 1 || sv[1] != 2 || sv[2] != 3 || sv[3] != 4 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Append operation failed\n"
          << " Details:\n"
          << "   Result:\n" << sv << "\n"
          << "   Expected result:\n( 0 0 1 2 3 4 0 0 0 )\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c erase() member function of the Subvector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c erase() member function of the Subvector
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testErase()
{
   //=====================================================================================
   // Index-based erase function
   //=====================================================================================

   {
      test_ = "Subvector::erase( size_t )";

      initialize();

      SVT sv = blaze::subvector( vec_, 1UL, 6UL );

      // Erasing the non-zero element at the end of the subvector
      sv.erase( 5UL );

      checkSize    ( sv  , 6UL );
      checkNonZeros( sv  , 3UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 3UL );

      if( sv[0] !=  1 || sv[1] != 0 || sv[2] != -2 ||
          sv[3] != -3 || sv[4] != 0 || sv[5] !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 1 0 -2 -3 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the non-zero element at the beginning of the subvector
      sv.erase( 0UL );

      checkSize    ( sv  , 6UL );
      checkNonZeros( sv  , 2UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 2UL );

      if( sv[0] !=  0 || sv[1] != 0 || sv[2] != -2 ||
          sv[3] != -3 || sv[4] != 0 || sv[5] !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 0 0 -2 -3 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the non-zero element at the beginning of the subvector
      sv.erase( 2UL );

      checkSize    ( sv  , 6UL );
      checkNonZeros( sv  , 1UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 1UL );

      if( sv[0] !=  0 || sv[1] != 0 || sv[2] != 0 ||
          sv[3] != -3 || sv[4] != 0 || sv[5] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 0 0 0 -3 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to erase an already erased element
      sv.erase( 2UL );

      checkSize    ( sv  , 6UL );
      checkNonZeros( sv  , 1UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 1UL );

      if( sv[0] !=  0 || sv[1] != 0 || sv[2] != 0 ||
          sv[3] != -3 || sv[4] != 0 || sv[5] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 0 0 0 -3 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Iterator-based erase function
   //=====================================================================================

   {
      test_ = "Subvector::erase( Iterator )";

      initialize();

      SVT sv = blaze::subvector( vec_, 1UL, 6UL );

      // Erasing the non-zero element at the end of the subvector
      {
         SVT::Iterator pos = sv.erase( sv.find( 5UL ) );

         checkSize    ( sv  , 6UL );
         checkNonZeros( sv  , 3UL );
         checkSize    ( vec_, 8UL );
         checkNonZeros( vec_, 3UL );

         if( pos != sv.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( sv[0] !=  1 || sv[1] != 0 || sv[2] != -2 ||
             sv[3] != -3 || sv[4] != 0 || sv[5] !=  0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << sv << "\n"
                << "   Expected result:\n( 1 0 -2 -3 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the non-zero element at the beginning of the subvector
      {
         SVT::Iterator pos = sv.erase( sv.find( 0UL ) );

         checkSize    ( sv  , 6UL );
         checkNonZeros( sv  , 2UL );
         checkSize    ( vec_, 8UL );
         checkNonZeros( vec_, 2UL );

         if( pos->value() != -2 || pos->index() != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: -2\n"
                << "   Expected index:  2\n";
            throw std::runtime_error( oss.str() );
         }

         if( sv[0] !=  0 || sv[1] != 0 || sv[2] != -2 ||
             sv[3] != -3 || sv[4] != 0 || sv[5] !=  0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << sv << "\n"
                << "   Expected result:\n( 0 0 -2 -3 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the non-zero element at the beginning of the subvector
      {
         SVT::Iterator pos = sv.erase( sv.find( 2UL ) );

         checkSize    ( sv  , 6UL );
         checkNonZeros( sv  , 1UL );
         checkSize    ( vec_, 8UL );
         checkNonZeros( vec_, 1UL );

         if( pos->value() != -3 || pos->index() != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: -3\n"
                << "   Expected index:  3\n";
            throw std::runtime_error( oss.str() );
         }

         if( sv[0] !=  0 || sv[1] != 0 || sv[2] != 0 ||
             sv[3] != -3 || sv[4] != 0 || sv[5] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << sv << "\n"
                << "   Expected result:\n( 0 0 0 -3 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to erase an already erased element
      {
         SVT::Iterator pos = sv.erase( sv.find( 2UL ) );

         checkSize    ( sv  , 6UL );
         checkNonZeros( sv  , 1UL );
         checkSize    ( vec_, 8UL );
         checkNonZeros( vec_, 1UL );

         if( pos != sv.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( sv[0] !=  0 || sv[1] != 0 || sv[2] != 0 ||
             sv[3] != -3 || sv[4] != 0 || sv[5] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << sv << "\n"
                << "   Expected result:\n( 0 0 0 -3 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Iterator-range-based erase function
   //=====================================================================================

   {
      test_ = "Subvector::erase( Iterator, Iterator )";

      // Erasing the entire vector
      {
         initialize();

         SVT sv = blaze::subvector( vec_, 0UL, 8UL );

         SVT::Iterator pos = sv.erase( sv.begin(), sv.end() );

         checkSize    ( sv  , 8UL );
         checkNonZeros( sv  , 0UL );
         checkSize    ( vec_, 8UL );
         checkNonZeros( vec_, 0UL );

         if( pos != sv.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( sv[0] != 0 || sv[1] != 0 || sv[2] != 0 || sv[3] != 0 ||
             sv[4] != 0 || sv[5] != 0 || sv[6] != 0 || sv[7] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing the subvector failed\n"
                << " Details:\n"
                << "   Result:\n" << sv << "\n"
                << "   Expected result:\n( 0 0 0 0 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the first half of the vector
      {
         initialize();

         SVT sv = blaze::subvector( vec_, 0UL, 8UL );

         SVT::Iterator pos = sv.erase( sv.begin(), sv.find( 4UL ) );

         checkSize    ( sv  , 8UL );
         checkNonZeros( sv  , 2UL );
         checkSize    ( vec_, 8UL );
         checkNonZeros( vec_, 2UL );

         if( pos->value() != -3 || pos->index() != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: -3\n"
                << "   Expected index:  4\n";
            throw std::runtime_error( oss.str() );
         }

         if( sv[0] !=  0 || sv[1] != 0 || sv[2] != 0 || sv[3] != 0 ||
             sv[4] != -3 || sv[5] != 0 || sv[6] != 4 || sv[7] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a partial subvector failed\n"
                << " Details:\n"
                << "   Result:\n" << sv << "\n"
                << "   Expected result:\n( 0 0 0 0 -3 0 4 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the second half of the vector
      {
         initialize();

         SVT sv = blaze::subvector( vec_, 0UL, 8UL );

         SVT::Iterator pos = sv.erase( sv.find( 4UL ), sv.end() );

         checkSize    ( sv  , 8UL );
         checkNonZeros( sv  , 2UL );
         checkSize    ( vec_, 8UL );
         checkNonZeros( vec_, 2UL );

         if( pos != sv.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( sv[0] != 0 || sv[1] != 1 || sv[2] != 0 || sv[3] != -2 ||
             sv[4] != 0 || sv[5] != 0 || sv[6] != 0 || sv[7] !=  0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a partial subvector failed\n"
                << " Details:\n"
                << "   Result:\n" << sv << "\n"
                << "   Expected result:\n( 0 1 0 -2 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to erase an empty range
      {
         initialize();

         SVT sv = blaze::subvector( vec_, 0UL, 8UL );

         SVT::Iterator pos = sv.erase( sv.find( 1UL ), sv.find( 1UL ) );

         checkSize    ( sv  , 8UL );
         checkNonZeros( sv  , 4UL );
         checkSize    ( vec_, 8UL );
         checkNonZeros( vec_, 4UL );

         if( pos != sv.find( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( sv[0] !=  0 || sv[1] != 1 || sv[2] != 0 || sv[3] != -2 ||
             sv[4] != -3 || sv[5] != 0 || sv[6] != 4 || sv[7] !=  0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing an empty range failed\n"
                << " Details:\n"
                << "   Result:\n" << sv << "\n"
                << "   Expected result:\n( 0 1 0 -2 -3 0 4 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   //  erase() function with predicate
   //=====================================================================================

   {
      test_ = "Subvector::erase( Predicate )";

      initialize();

      SVT sv = blaze::subvector( vec_, 1UL, 6UL );

      // Erasing a selection of elements
      sv.erase( []( int value ){ return value == 1 || value == 4; } );

      checkSize    ( sv  , 6UL );
      checkNonZeros( sv  , 2UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 2UL );

      if( sv[0] !=  0 || sv[1] != 0 || sv[2] != -2 ||
          sv[3] != -3 || sv[4] != 0 || sv[5] !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a selection of elements failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 0 0 -2 -3 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to erase all elements with value 1
      sv.erase( []( int value ){ return value == 1; } );

      checkSize    ( sv  , 6UL );
      checkNonZeros( sv  , 2UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 2UL );

      if( sv[0] !=  0 || sv[1] != 0 || sv[2] != -2 ||
          sv[3] != -3 || sv[4] != 0 || sv[5] !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing all elements with value 1 failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 0 0 -2 -3 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Iterator-range-based erase() function with predicate
   //=====================================================================================

   {
      test_ = "Subvector::erase( Iterator, Iterator, Predicate )";

      initialize();

      SVT sv = blaze::subvector( vec_, 1UL, 6UL );

      // Erasing a selection of elements
      sv.erase( sv.begin(), sv.find( 3UL ), []( int value ){ return value == 1; } );

      checkSize    ( sv  , 6UL );
      checkNonZeros( sv  , 3UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 3UL );

      if( sv[0] !=  0 || sv[1] != 0 || sv[2] != -2 ||
          sv[3] != -3 || sv[4] != 0 || sv[5] !=  4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a selection of elements failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 0 0 -2 -3 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to erase from an empty range
      sv.erase( sv.begin(), sv.begin(), []( int value ){ return value == 1; } );

      checkSize    ( sv  , 6UL );
      checkNonZeros( sv  , 3UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 3UL );

      if( sv[0] !=  0 || sv[1] != 0 || sv[2] != -2 ||
          sv[3] != -3 || sv[4] != 0 || sv[5] !=  4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing from an empty range failed\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n"
             << "   Expected result:\n( 0 0 -2 -3 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c find() member function of the Subvector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c find() member function of the Subvector
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testFind()
{
   test_ = "Subvector::find()";

   initialize();

   SVT sv = blaze::subvector( vec_, 1UL, 5UL );

   // Searching for the first element
   {
      SVT::Iterator pos = sv.find( 0UL );

      if( pos == sv.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Element could not be found\n"
             << " Details:\n"
             << "   Required index = 0\n"
             << "   Current subvector:\n" << sv << "\n";
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
             << "   Current subvector:\n" << sv << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Searching for the second element
   {
      SVT::Iterator pos = sv.find( 2UL );

      if( pos == sv.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Element could not be found\n"
             << " Details:\n"
             << "   Required index = 2\n"
             << "   Current subvector:\n" << sv << "\n";
         throw std::runtime_error( oss.str() );
      }
      else if( pos->index() != 2 || pos->value() != -2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Wrong element found\n"
             << " Details:\n"
             << "   Required index = 2\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = -2\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current subvector:\n" << sv << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Searching for a non-existing non-zero element
   {
      SVT::Iterator pos = sv.find( 1UL );

      if( pos != sv.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Non-existing element could be found\n"
             << " Details:\n"
             << "   Required index = 1\n"
             << "   Current subvector:\n" << sv << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c lowerBound() member function of the Subvector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c lowerBound() member function of the Subvector
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testLowerBound()
{
   test_ = "Subvector::lowerBound()";

   initialize();

   SVT sv = blaze::subvector( vec_, 0UL, 3UL );

   // Determining the lower bound for index 0
   {
      SVT::Iterator pos = sv.lowerBound( 0UL );

      if( pos == sv.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Lower bound could not be determined\n"
             << " Details:\n"
             << "   Required index = 0\n"
             << "   Current subvector:\n" << sv << "\n";
         throw std::runtime_error( oss.str() );
      }
      else if( pos->index() != 1 || pos->value() != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Wrong element found\n"
             << " Details:\n"
             << "   Required index = 1\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = 1\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current subvector:\n" << sv << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Determining the lower bound for index 1
   {
      SVT::Iterator pos = sv.lowerBound( 1UL );

      if( pos == sv.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Lower bound could not be determined\n"
             << " Details:\n"
             << "   Required index = 1\n"
             << "   Current subvector:\n" << sv << "\n";
         throw std::runtime_error( oss.str() );
      }
      else if( pos->index() != 1 || pos->value() != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Wrong element found\n"
             << " Details:\n"
             << "   Required index = 1\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = 1\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current subvector:\n" << sv << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Determining the lower bound for index 2
   {
      SVT::Iterator pos = sv.lowerBound( 2UL );

      if( pos != sv.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Lower bound could not be determined\n"
             << " Details:\n"
             << "   Required index = 2\n"
             << "   Current subvector:\n" << sv << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c upperBound() member function of the Subvector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c upperBound() member function of the Subvector
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testUpperBound()
{
   test_ = "Subvector::upperBound()";

   initialize();

   SVT sv = blaze::subvector( vec_, 0UL, 3UL );

   // Determining the upper bound for index 0
   {
      SVT::Iterator pos = sv.upperBound( 0UL );

      if( pos == sv.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Upper bound could not be determined\n"
             << " Details:\n"
             << "   Required index = 0\n"
             << "   Current subvector:\n" << sv << "\n";
         throw std::runtime_error( oss.str() );
      }
      else if( pos->index() != 1 || pos->value() != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Wrong element found\n"
             << " Details:\n"
             << "   Required index = 1\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = 1\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current subvector:\n" << sv << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Determining the upper bound for index 1
   {
      SVT::Iterator pos = sv.upperBound( 1UL );

      if( pos != sv.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Upper bound could not be determined\n"
             << " Details:\n"
             << "   Required index = 1\n"
             << "   Current subvector:\n" << sv << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Determining the upper bound for index 2
   {
      SVT::Iterator pos = sv.upperBound( 2UL );

      if( pos != sv.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Upper bound could not be determined\n"
             << " Details:\n"
             << "   Required index = 2\n"
             << "   Current subvector:\n" << sv << "\n";
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
void SparseTest::testIsDefault()
{
   using blaze::isDefault;

   test_ = "isDefault() function";

   initialize();

   // isDefault with default vector
   {
      VT vec( 8UL );
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
void SparseTest::testIsSame()
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
void SparseTest::testSubvector()
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

      if( sv2.begin()->value() != -2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << sv2.begin()->value() << "\n"
             << "   Expected result: -2\n";
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
void SparseTest::testElements()
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

         if( e.begin()->value() != -3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator access failed\n"
                << " Details:\n"
                << "   Result: " << e.begin()->value() << "\n"
                << "   Expected result: -3\n";
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

         if( e.begin()->value() != -3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator access failed\n"
                << " Details:\n"
                << "   Result: " << e.begin()->value() << "\n"
                << "   Expected result: -3\n";
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

         if( e.begin()->value() != -3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator access failed\n"
                << " Details:\n"
                << "   Result: " << e.begin()->value() << "\n"
                << "   Expected result: -3\n";
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
void SparseTest::initialize()
{
   // Initializing the compressed row vector
   vec_.reset();
   vec_[1] =  1;
   vec_[3] = -2;
   vec_[4] = -3;
   vec_[6] =  4;
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
   std::cout << "   Running Subvector sparse test..." << std::endl;

   try
   {
      RUN_SUBVECTOR_SPARSE_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during Subvector sparse test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
