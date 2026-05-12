//=================================================================================================
/*!
//  \file src/mathtest/views/column/SparseGeneralTest.cpp
//  \brief Source file for the Column sparse general test
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
#include <blaze/math/CompressedVector.h>
#include <blaze/math/DynamicVector.h>
#include <blaze/math/Views.h>
#include <blazetest/mathtest/views/column/SparseGeneralTest.h>

#ifdef BLAZE_USE_HPX_THREADS
#  include <hpx/hpx_main.hpp>
#endif


namespace blazetest {

namespace mathtest {

namespace views {

namespace column {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the Column sparse general test.
//
// \exception std::runtime_error Operation error detected.
*/
SparseGeneralTest::SparseGeneralTest()
   : mat_ ( 4UL, 5UL )
   , tmat_( 4UL, 5UL )
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
/*!\brief Test of the Column constructors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all constructors of the Column specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseGeneralTest::testConstructors()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major Column constructor (0x0)";

      MT mat;

      // 0th matrix column
      try {
         blaze::column( mat, 0UL );
      }
      catch( std::invalid_argument& ) {}
   }

   {
      test_ = "Row-major Column constructor (0x2)";

      MT mat( 0UL, 2UL );

      // 0th matrix column
      {
         CT col0 = blaze::column( mat, 0UL );

         checkSize    ( col0, 0UL );
         checkCapacity( col0, 0UL );
         checkNonZeros( col0, 0UL );
      }

      // 1st matrix column
      {
         CT col1 = blaze::column( mat, 1UL );

         checkSize    ( col1, 0UL );
         checkCapacity( col1, 0UL );
         checkNonZeros( col1, 0UL );
      }

      // 2nd matrix column
      try {
         blaze::column( mat, 2UL );
      }
      catch( std::invalid_argument& ) {}
   }

   {
      test_ = "Row-major Column constructor (4x5)";

      initialize();

      // 0th matrix column
      {
         CT col0 = blaze::column( mat_, 0UL );

         checkSize    ( col0, 4UL );
         checkNonZeros( col0, 0UL );

         if( col0[0] != 0 || col0[1] != 0 || col0[2] != 0 || col0[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of 0th sparse column failed\n"
                << " Details:\n"
                << "   Result:\n" << col0 << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // 1st matrix column
      {
         CT col1 = blaze::column( mat_, 1UL );

         checkSize    ( col1, 4UL );
         checkNonZeros( col1, 1UL );

         if( col1[0] != 0 || col1[1] != 1 || col1[2] != 0 || col1[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of 1st sparse column failed\n"
                << " Details:\n"
                << "   Result:\n" << col1 << "\n"
                << "   Expected result:\n( 0 1 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // 2nd matrix column
      {
         CT col2 = blaze::column( mat_, 2UL );

         checkSize    ( col2, 4UL );
         checkNonZeros( col2, 2UL );

         if( col2[0] != -2 || col2[1] != 0 || col2[2] != -3 || col2[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of 2nd sparse column failed\n"
                << " Details:\n"
                << "   Result:\n" << col2 << "\n"
                << "   Expected result:\n( -2 0 -3 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // 3rd matrix column
      {
         CT col3 = blaze::column( mat_, 3UL );

         checkSize    ( col3, 4UL );
         checkNonZeros( col3, 3UL );

         if( col3[0] != 0 || col3[1] != 4 || col3[2] != 5 || col3[3] != -6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of 3rd sparse column failed\n"
                << " Details:\n"
                << "   Result:\n" << col3 << "\n"
                << "   Expected result:\n( 0 4 5 -6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // 4th matrix column
      {
         CT col4 = blaze::column( mat_, 4UL );

         checkSize    ( col4, 4UL );
         checkNonZeros( col4, 4UL );

         if( col4[0] != 7 || col4[1] != -8 || col4[2] != 9 || col4[3] != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of 4th sparse column failed\n"
                << " Details:\n"
                << "   Result:\n" << col4 << "\n"
                << "   Expected result:\n( 7 -8 9 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // 5th matrix column
      try {
         blaze::column( mat_, 5UL );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major Column constructor (0x0)";

      MT tmat;

      // 0th matrix column
      try {
         blaze::column( tmat, 0UL );
      }
      catch( std::invalid_argument& ) {}
   }

   {
      test_ = "Column-major Column constructor (0x2)";

      MT tmat( 0UL, 2UL );

      // 0th matrix column
      {
         CT col0 = blaze::column( tmat, 0UL );

         checkSize    ( col0, 0UL );
         checkCapacity( col0, 0UL );
         checkNonZeros( col0, 0UL );
      }

      // 1st matrix column
      {
         CT col1 = blaze::column( tmat, 1UL );

         checkSize    ( col1, 0UL );
         checkCapacity( col1, 0UL );
         checkNonZeros( col1, 0UL );
      }

      // 2nd matrix column
      try {
         blaze::column( tmat, 2UL );
      }
      catch( std::invalid_argument& ) {}
   }

   {
      test_ = "Column-major Column constructor (4x5)";

      initialize();

      // 0th matrix column
      {
         OCT col0 = blaze::column( tmat_, 0UL );

         checkSize    ( col0, 4UL );
         checkNonZeros( col0, 0UL );

         if( col0[0] != 0 || col0[1] != 0 || col0[2] != 0 || col0[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of 0th sparse column failed\n"
                << " Details:\n"
                << "   Result:\n" << col0 << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // 1st matrix column
      {
         OCT col1 = blaze::column( tmat_, 1UL );

         checkSize    ( col1, 4UL );
         checkNonZeros( col1, 1UL );

         if( col1[0] != 0 || col1[1] != 1 || col1[2] != 0 || col1[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of 1st sparse column failed\n"
                << " Details:\n"
                << "   Result:\n" << col1 << "\n"
                << "   Expected result:\n( 0 1 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // 2nd matrix column
      {
         OCT col2 = blaze::column( tmat_, 2UL );

         checkSize    ( col2, 4UL );
         checkNonZeros( col2, 2UL );

         if( col2[0] != -2 || col2[1] != 0 || col2[2] != -3 || col2[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of 2nd sparse column failed\n"
                << " Details:\n"
                << "   Result:\n" << col2 << "\n"
                << "   Expected result:\n( -2 0 -3 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // 3rd matrix column
      {
         OCT col3 = blaze::column( tmat_, 3UL );

         checkSize    ( col3, 4UL );
         checkNonZeros( col3, 3UL );

         if( col3[0] != 0 || col3[1] != 4 || col3[2] != 5 || col3[3] != -6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of 3rd sparse column failed\n"
                << " Details:\n"
                << "   Result:\n" << col3 << "\n"
                << "   Expected result:\n( 0 4 5 -6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // 4th matrix column
      {
         OCT col4 = blaze::column( tmat_, 4UL );

         checkSize    ( col4, 4UL );
         checkNonZeros( col4, 4UL );

         if( col4[0] != 7 || col4[1] != -8 || col4[2] != 9 || col4[3] != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of 4th sparse column failed\n"
                << " Details:\n"
                << "   Result:\n" << col4 << "\n"
                << "   Expected result:\n( 7 -8 9 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // 5th matrix column
      try {
         blaze::column( tmat_, 5UL );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Column assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all assignment operators of the Column specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseGeneralTest::testAssignment()
{
   //=====================================================================================
   // Row-major list assignment
   //=====================================================================================

   {
      test_ = "Row-major initializer list assignment (complete list)";

      initialize();

      CT col3 = blaze::column( mat_, 3UL );
      col3 = { 1, 2, 3, 4 };

      checkSize    ( col3,  4UL );
      checkNonZeros( col3,  4UL );
      checkRows    ( mat_,  4UL );
      checkColumns ( mat_,  5UL );
      checkNonZeros( mat_, 11UL );

      if( col3[0] != 1 || col3[1] != 2 || col3[2] != 3 || col3[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << col3 << "\n"
             << "   Expected result:\n( 1 2 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) != -2 || mat_(0,3) != 1 || mat_(0,4) !=  7 ||
          mat_(1,0) != 0 || mat_(1,1) != 1 || mat_(1,2) !=  0 || mat_(1,3) != 2 || mat_(1,4) != -8 ||
          mat_(2,0) != 0 || mat_(2,1) != 0 || mat_(2,2) != -3 || mat_(2,3) != 3 || mat_(2,4) !=  9 ||
          mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) != 4 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0 -2  1  7 )\n"
                                     "( 0  1  0  2 -8 )\n"
                                     "( 0  0 -3  3  9 )\n"
                                     "( 0  0  0  4 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major initializer list assignment (incomplete list)";

      initialize();

      CT col3 = blaze::column( mat_, 3UL );
      col3 = { 1, 2 };

      checkSize    ( col3, 4UL );
      checkNonZeros( col3, 2UL );
      checkRows    ( mat_, 4UL );
      checkColumns ( mat_, 5UL );
      checkNonZeros( mat_, 9UL );

      if( col3[0] != 1 || col3[1] != 2 || col3[2] != 0 || col3[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << col3 << "\n"
             << "   Expected result:\n( 1 2 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) != -2 || mat_(0,3) != 1 || mat_(0,4) !=  7 ||
          mat_(1,0) != 0 || mat_(1,1) != 1 || mat_(1,2) !=  0 || mat_(1,3) != 2 || mat_(1,4) != -8 ||
          mat_(2,0) != 0 || mat_(2,1) != 0 || mat_(2,2) != -3 || mat_(2,3) != 0 || mat_(2,4) !=  9 ||
          mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) != 0 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0 -2  1  7 )\n"
                                     "( 0  1  0  2 -8 )\n"
                                     "( 0  0 -3  0  9 )\n"
                                     "( 0  0  0  0 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major copy assignment
   //=====================================================================================

   {
      test_ = "Row-major Column copy assignment";

      initialize();

      CT col1 = blaze::column( mat_, 1UL );
      col1 = blaze::column( mat_, 2UL );

      checkSize    ( col1,  4UL );
      checkNonZeros( col1,  2UL );
      checkRows    ( mat_,  4UL );
      checkColumns ( mat_,  5UL );
      checkNonZeros( mat_, 11UL );

      if( col1[0] != -2 || col1[1] != 0 || col1[2] != -3 || col1[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( -2 0 -3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) != -2 || mat_(0,2) != -2 || mat_(0,3) !=  0 || mat_(0,4) !=  7 ||
          mat_(1,0) != 0 || mat_(1,1) !=  0 || mat_(1,2) !=  0 || mat_(1,3) !=  4 || mat_(1,4) != -8 ||
          mat_(2,0) != 0 || mat_(2,1) != -3 || mat_(2,2) != -3 || mat_(2,3) !=  5 || mat_(2,4) !=  9 ||
          mat_(3,0) != 0 || mat_(3,1) !=  0 || mat_(3,2) !=  0 || mat_(3,3) != -6 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0 -2 -2  0  7 )\n"
                                     "( 0  0  0  4 -8 )\n"
                                     "( 0 -3 -3  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major dense vector assignment
   //=====================================================================================

   {
      test_ = "Row-major dense vector assignment";

      initialize();

      CT col1 = blaze::column( mat_, 1UL );

      blaze::DynamicVector<int,blaze::columnVector> vec1{ 0, 8, 0, 9 };

      col1 = vec1;

      checkSize    ( col1,  4UL );
      checkNonZeros( col1,  2UL );
      checkRows    ( mat_,  4UL );
      checkColumns ( mat_,  5UL );
      checkNonZeros( mat_, 11UL );

      if( col1[0] != 0 || col1[1] != 8 || col1[2] != 0 || col1[3] != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( 0 8 0 9 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) != -2 || mat_(0,3) !=  0 || mat_(0,4) !=  7 ||
          mat_(1,0) != 0 || mat_(1,1) != 8 || mat_(1,2) !=  0 || mat_(1,3) !=  4 || mat_(1,4) != -8 ||
          mat_(2,0) != 0 || mat_(2,1) != 0 || mat_(2,2) != -3 || mat_(2,3) !=  5 || mat_(2,4) !=  9 ||
          mat_(3,0) != 0 || mat_(3,1) != 9 || mat_(3,2) !=  0 || mat_(3,3) != -6 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0 -2  0  7 )\n"
                                     "( 0  8  0  4 -8 )\n"
                                     "( 0  0 -3  5  9 )\n"
                                     "( 0  9  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse vector assignment
   //=====================================================================================

   {
      test_ = "Row-major sparse vector assignment";

      initialize();

      CT col4 = blaze::column( mat_, 4UL );

      blaze::CompressedVector<int,blaze::columnVector> vec1( 4UL );
      vec1[3] = 9;

      col4 = vec1;

      checkSize    ( col4, 4UL );
      checkNonZeros( col4, 1UL );
      checkRows    ( mat_, 4UL );
      checkColumns ( mat_, 5UL );
      checkNonZeros( mat_, 7UL );

      if( col4[0] != 0 || col4[1] != 0 || col4[2] != 0 || col4[3] != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << col4 << "\n"
             << "   Expected result:\n( 0 0 0 9 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) != -2 || mat_(0,3) !=  0 || mat_(0,4) != 0 ||
          mat_(1,0) != 0 || mat_(1,1) != 1 || mat_(1,2) !=  0 || mat_(1,3) !=  4 || mat_(1,4) != 0 ||
          mat_(2,0) != 0 || mat_(2,1) != 0 || mat_(2,2) != -3 || mat_(2,3) !=  5 || mat_(2,4) != 0 ||
          mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) != -6 || mat_(3,4) != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0 -2  0  0 )\n"
                                     "( 0  1  0  4  0 )\n"
                                     "( 0  0 -3  5  0 )\n"
                                     "( 0  0  0 -6  9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major list assignment
   //=====================================================================================

   {
      test_ = "Column-major initializer list assignment (complete list)";

      initialize();

      OCT col3 = blaze::column( tmat_, 3UL );
      col3 = { 1, 2, 3, 4 };

      checkSize    ( col3 ,  4UL );
      checkNonZeros( col3 ,  4UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 11UL );

      if( col3[0] != 1 || col3[1] != 2 || col3[2] != 3 || col3[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << col3 << "\n"
             << "   Expected result:\n( 1 2 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -2 || tmat_(0,3) != 1 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) !=  0 || tmat_(1,3) != 2 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != -3 || tmat_(2,3) != 3 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != 4 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -2  1  7 )\n"
                                     "( 0  1  0  2 -8 )\n"
                                     "( 0  0 -3  3  9 )\n"
                                     "( 0  0  0  4 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major initializer list assignment (incomplete list)";

      initialize();

      OCT col3 = blaze::column( tmat_, 3UL );
      col3 = { 1, 2 };

      checkSize    ( col3 , 4UL );
      checkNonZeros( col3 , 2UL );
      checkRows    ( tmat_, 4UL );
      checkColumns ( tmat_, 5UL );
      checkNonZeros( tmat_, 9UL );

      if( col3[0] != 1 || col3[1] != 2 || col3[2] != 0 || col3[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << col3 << "\n"
             << "   Expected result:\n( 1 2 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -2 || tmat_(0,3) != 1 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) !=  0 || tmat_(1,3) != 2 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != -3 || tmat_(2,3) != 0 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != 0 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -2  1  7 )\n"
                                     "( 0  1  0  2 -8 )\n"
                                     "( 0  0 -3  0  9 )\n"
                                     "( 0  0  0  0 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major copy assignment
   //=====================================================================================

   {
      test_ = "Column-major Column copy assignment";

      initialize();

      OCT col1 = blaze::column( tmat_, 1UL );
      col1 = blaze::column( tmat_, 2UL );

      checkSize    ( col1 ,  4UL );
      checkNonZeros( col1 ,  2UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 11UL );

      if( col1[0] != -2 || col1[1] != 0 || col1[2] != -3 || col1[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( -2 0 -3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != -2 || tmat_(0,2) != -2 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) !=  0 || tmat_(1,2) !=  0 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != -3 || tmat_(2,2) != -3 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0 -2 -2  0  7 )\n"
                                     "( 0  0  0  4 -8 )\n"
                                     "( 0 -3 -3  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense vector assignment
   //=====================================================================================

   {
      test_ = "Column-major dense vector assignment";

      initialize();

      OCT col1 = blaze::column( tmat_, 1UL );

      blaze::DynamicVector<int,blaze::columnVector> vec1{ 0, 8, 0, 9 };

      col1 = vec1;

      checkSize    ( col1 ,  4UL );
      checkNonZeros( col1 ,  2UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 11UL );

      if( col1[0] != 0 || col1[1] != 8 || col1[2] != 0 || col1[3] != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( 0 8 0 9 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -2 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 8 || tmat_(1,2) !=  0 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != -3 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 9 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -2  0  7 )\n"
                                     "( 0  8  0  4 -8 )\n"
                                     "( 0  0 -3  5  9 )\n"
                                     "( 0  9  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse vector assignment
   //=====================================================================================

   {
      test_ = "Column-major sparse vector assignment";

      initialize();

      OCT col4 = blaze::column( tmat_, 4UL );

      blaze::CompressedVector<int,blaze::columnVector> vec1( 4UL );
      vec1[3] = 9;

      col4 = vec1;

      checkSize    ( col4 , 4UL );
      checkNonZeros( col4 , 1UL );
      checkRows    ( tmat_, 4UL );
      checkColumns ( tmat_, 5UL );
      checkNonZeros( tmat_, 7UL );

      if( col4[0] != 0 || col4[1] != 0 || col4[2] != 0 || col4[3] != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << col4 << "\n"
             << "   Expected result:\n( 0 0 0 9 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -2 || tmat_(0,3) !=  0 || tmat_(0,4) != 0 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) !=  0 || tmat_(1,3) !=  4 || tmat_(1,4) != 0 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != -3 || tmat_(2,3) !=  5 || tmat_(2,4) != 0 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -2  0  0 )\n"
                                     "( 0  1  0  4  0 )\n"
                                     "( 0  0 -3  5  0 )\n"
                                     "( 0  0  0 -6  9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Column addition assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the addition assignment operators of the Column
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseGeneralTest::testAddAssign()
{
   //=====================================================================================
   // Row-major Column addition assignment
   //=====================================================================================

   {
      test_ = "Row-major Column addition assignment";

      initialize();

      CT col2 = blaze::column( mat_, 2UL );
      col2 += blaze::column( mat_, 3UL );

      checkSize    ( col2,  4UL );
      checkNonZeros( col2,  4UL );
      checkRows    ( mat_,  4UL );
      checkColumns ( mat_,  5UL );
      checkNonZeros( mat_, 12UL );

      if( col2[0] != -2 || col2[1] != 4 || col2[2] != 2 || col2[3] != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( -2 4 2 -6 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) != -2 || mat_(0,3) !=  0 || mat_(0,4) !=  7 ||
          mat_(1,0) != 0 || mat_(1,1) != 1 || mat_(1,2) !=  4 || mat_(1,3) !=  4 || mat_(1,4) != -8 ||
          mat_(2,0) != 0 || mat_(2,1) != 0 || mat_(2,2) !=  2 || mat_(2,3) !=  5 || mat_(2,4) !=  9 ||
          mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) != -6 || mat_(3,3) != -6 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0 -2  0  7 )\n"
                                     "( 0  1  4  4 -8 )\n"
                                     "( 0  0  2  5  9 )\n"
                                     "( 0  0 -6 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major dense vector addition assignment
   //=====================================================================================

   {
      test_ = "Row-major dense vector addition assignment";

      initialize();

      CT col2 = blaze::column( mat_, 2UL );

      blaze::DynamicVector<int,blaze::columnVector> vec{ 2, -4, 0, 0 };

      col2 += vec;

      checkSize    ( col2,  4UL );
      checkNonZeros( col2,  2UL );
      checkRows    ( mat_,  4UL );
      checkColumns ( mat_,  5UL );
      checkNonZeros( mat_, 10UL );

      if( col2[0] != 0 || col2[1] != -4 || col2[2] != -3 || col2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( 0 -4 -3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 || mat_(0,4) !=  7 ||
          mat_(1,0) != 0 || mat_(1,1) != 1 || mat_(1,2) != -4 || mat_(1,3) !=  4 || mat_(1,4) != -8 ||
          mat_(2,0) != 0 || mat_(2,1) != 0 || mat_(2,2) != -3 || mat_(2,3) !=  5 || mat_(2,4) !=  9 ||
          mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) != -6 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0  0  0  7 )\n"
                                     "( 0  1 -4  4 -8 )\n"
                                     "( 0  0 -3  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse vector addition assignment
   //=====================================================================================

   {
      test_ = "Row-major sparse vector addition assignment";

      initialize();

      CT col2 = blaze::column( mat_, 2UL );

      blaze::CompressedVector<int,blaze::columnVector> vec( 4UL );
      vec[0] =  2;
      vec[1] = -4;

      col2 += vec;

      checkSize    ( col2,  4UL );
      checkNonZeros( col2,  2UL );
      checkRows    ( mat_,  4UL );
      checkColumns ( mat_,  5UL );
      checkNonZeros( mat_, 10UL );

      if( col2[0] != 0 || col2[1] != -4 || col2[2] != -3 || col2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( 0 -4 -3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) !=  0 || mat_(0,3) !=  0 || mat_(0,4) !=  7 ||
          mat_(1,0) != 0 || mat_(1,1) != 1 || mat_(1,2) != -4 || mat_(1,3) !=  4 || mat_(1,4) != -8 ||
          mat_(2,0) != 0 || mat_(2,1) != 0 || mat_(2,2) != -3 || mat_(2,3) !=  5 || mat_(2,4) !=  9 ||
          mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) != -6 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0  0  0  7 )\n"
                                     "( 0  1 -4  4 -8 )\n"
                                     "( 0  0 -3  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major Column addition assignment
   //=====================================================================================

   {
      test_ = "Column-major Column addition assignment";

      initialize();

      OCT col2 = blaze::column( tmat_, 2UL );
      col2 += blaze::column( tmat_, 3UL );

      checkSize    ( col2 ,  4UL );
      checkNonZeros( col2 ,  4UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 12UL );

      if( col2[0] != -2 || col2[1] != 4 || col2[2] != 2 || col2[3] != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( -2 4 2 -6 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -2 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) !=  4 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) !=  2 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) != -6 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -2  0  7 )\n"
                                     "( 0  1  4  4 -8 )\n"
                                     "( 0  0  2  5  9 )\n"
                                     "( 0  0 -6 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense vector addition assignment
   //=====================================================================================

   {
      test_ = "Column-major dense vector addition assignment";

      initialize();

      OCT col2 = blaze::column( tmat_, 2UL );

      blaze::DynamicVector<int,blaze::columnVector> vec{ 2, -4, 0, 0 };

      col2 += vec;

      checkSize    ( col2 ,  4UL );
      checkNonZeros( col2 ,  2UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 10UL );

      if( col2[0] != 0 || col2[1] != -4 || col2[2] != -3 || col2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( 0 -4 -3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) != -4 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != -3 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0  0  0  7 )\n"
                                     "( 0  1 -4  4 -8 )\n"
                                     "( 0  0 -3  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse vector addition assignment
   //=====================================================================================

   {
      test_ = "Column-major sparse vector addition assignment";

      initialize();

      OCT col2 = blaze::column( tmat_, 2UL );

      blaze::CompressedVector<int,blaze::columnVector> vec( 4UL );
      vec[0] =  2;
      vec[1] = -4;

      col2 += vec;

      checkSize    ( col2 ,  4UL );
      checkNonZeros( col2 ,  2UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 10UL );

      if( col2[0] != 0 || col2[1] != -4 || col2[2] != -3 || col2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( 0 -4 -3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) != -4 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != -3 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0  0  0  7 )\n"
                                     "( 0  1 -4  4 -8 )\n"
                                     "( 0  0 -3  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Column subtraction assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the subtraction assignment operators of the Column
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseGeneralTest::testSubAssign()
{
   //=====================================================================================
   // Row-major Column subtraction assignment
   //=====================================================================================

   {
      test_ = "Row-major Column subtraction assignment";

      initialize();

      CT col2 = blaze::column( mat_, 2UL );
      col2 -= blaze::column( mat_, 3UL );

      checkSize    ( col2,  4UL );
      checkNonZeros( col2,  4UL );
      checkRows    ( mat_,  4UL );
      checkColumns ( mat_,  5UL );
      checkNonZeros( mat_, 12UL );

      if( col2[0] != -2 || col2[1] != -4 || col2[2] != -8 || col2[3] != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( -2 -4 -8 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) != -2 || mat_(0,3) !=  0 || mat_(0,4) !=  7 ||
          mat_(1,0) != 0 || mat_(1,1) != 1 || mat_(1,2) != -4 || mat_(1,3) !=  4 || mat_(1,4) != -8 ||
          mat_(2,0) != 0 || mat_(2,1) != 0 || mat_(2,2) != -8 || mat_(2,3) !=  5 || mat_(2,4) !=  9 ||
          mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) !=  6 || mat_(3,3) != -6 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0 -2  0  7 )\n"
                                     "( 0  1 -4  4 -8 )\n"
                                     "( 0  0 -8  5  9 )\n"
                                     "( 0  0  6 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major dense vector subtraction assignment
   //=====================================================================================

   {
      test_ = "Row-major dense vector subtraction assignment";

      initialize();

      CT col2 = blaze::column( mat_, 2UL );

      blaze::DynamicVector<int,blaze::columnVector> vec{ 2, -4, 0, 0 };

      col2 -= vec;

      checkSize    ( col2,  4UL );
      checkNonZeros( col2,  3UL );
      checkRows    ( mat_,  4UL );
      checkColumns ( mat_,  5UL );
      checkNonZeros( mat_, 11UL );

      if( col2[0] != -4 || col2[1] != 4 || col2[2] != -3 || col2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( -4 4 -3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) != -4 || mat_(0,3) !=  0 || mat_(0,4) !=  7 ||
          mat_(1,0) != 0 || mat_(1,1) != 1 || mat_(1,2) !=  4 || mat_(1,3) !=  4 || mat_(1,4) != -8 ||
          mat_(2,0) != 0 || mat_(2,1) != 0 || mat_(2,2) != -3 || mat_(2,3) !=  5 || mat_(2,4) !=  9 ||
          mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) != -6 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0 -4  0  7 )\n"
                                     "( 0  1  4  4 -8 )\n"
                                     "( 0  0 -3  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse vector subtraction assignment
   //=====================================================================================

   {
      test_ = "Row-major sparse vector subtraction assignment";

      initialize();

      CT col2 = blaze::column( mat_, 2UL );

      blaze::CompressedVector<int,blaze::columnVector> vec( 4UL );
      vec[0] =  2;
      vec[1] = -4;

      col2 -= vec;

      checkSize    ( col2,  4UL );
      checkNonZeros( col2,  3UL );
      checkRows    ( mat_,  4UL );
      checkColumns ( mat_,  5UL );
      checkNonZeros( mat_, 11UL );

      if( col2[0] != -4 || col2[1] != 4 || col2[2] != -3 || col2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( -4 4 -3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) != -4 || mat_(0,3) !=  0 || mat_(0,4) !=  7 ||
          mat_(1,0) != 0 || mat_(1,1) != 1 || mat_(1,2) !=  4 || mat_(1,3) !=  4 || mat_(1,4) != -8 ||
          mat_(2,0) != 0 || mat_(2,1) != 0 || mat_(2,2) != -3 || mat_(2,3) !=  5 || mat_(2,4) !=  9 ||
          mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) != -6 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0 -4  0  7 )\n"
                                     "( 0  1  4  4 -8 )\n"
                                     "( 0  0 -3  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major Column subtraction assignment
   //=====================================================================================

   {
      test_ = "Column-major Column subtraction assignment";

      initialize();

      OCT col2 = blaze::column( tmat_, 2UL );
      col2 -= blaze::column( tmat_, 3UL );

      checkSize    ( col2 ,  4UL );
      checkNonZeros( col2 ,  4UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 12UL );

      if( col2[0] != -2 || col2[1] != -4 || col2[2] != -8 || col2[3] != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( -2 -4 -8 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -2 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) != -4 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != -8 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  6 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -2  0  7 )\n"
                                     "( 0  1 -4  4 -8 )\n"
                                     "( 0  0 -8  5  9 )\n"
                                     "( 0  0  6 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense vector subtraction assignment
   //=====================================================================================

   {
      test_ = "Column-major dense vector subtraction assignment";

      initialize();

      OCT col2 = blaze::column( tmat_, 2UL );

      blaze::DynamicVector<int,blaze::columnVector> vec{ 2, -4, 0, 0 };

      col2 -= vec;

      checkSize    ( col2 ,  4UL );
      checkNonZeros( col2 ,  3UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 11UL );

      if( col2[0] != -4 || col2[1] != 4 || col2[2] != -3 || col2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( -4 4 -3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -4 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) !=  4 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != -3 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -4  0  7 )\n"
                                     "( 0  1  4  4 -8 )\n"
                                     "( 0  0 -3  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse vector subtraction assignment
   //=====================================================================================

   {
      test_ = "Column-major sparse vector subtraction assignment";

      initialize();

      OCT col2 = blaze::column( tmat_, 2UL );

      blaze::CompressedVector<int,blaze::columnVector> vec( 4UL );
      vec[0] =  2;
      vec[1] = -4;

      col2 -= vec;

      checkSize    ( col2 ,  4UL );
      checkNonZeros( col2 ,  3UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 11UL );

      if( col2[0] != -4 || col2[1] != 4 || col2[2] != -3 || col2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( -4 4 -3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -4 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) !=  4 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != -3 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -4  0  7 )\n"
                                     "( 0  1  4  4 -8 )\n"
                                     "( 0  0 -3  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Column multiplication assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the multiplication assignment operators of the Column
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseGeneralTest::testMultAssign()
{
   //=====================================================================================
   // Row-major Column multiplication assignment
   //=====================================================================================

   {
      test_ = "Row-major Column multiplication assignment";

      initialize();

      CT col2 = blaze::column( mat_, 2UL );
      col2 *= blaze::column( mat_, 3UL );

      checkSize    ( col2, 4UL );
      checkNonZeros( col2, 1UL );
      checkRows    ( mat_, 4UL );
      checkColumns ( mat_, 5UL );
      checkNonZeros( mat_, 9UL );

      if( col2[0] != 0 || col2[1] != 0 || col2[2] != -15 || col2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( 0 0 -15 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) !=   0 || mat_(0,3) !=  0 || mat_(0,4) !=  7 ||
          mat_(1,0) != 0 || mat_(1,1) != 1 || mat_(1,2) !=   0 || mat_(1,3) !=  4 || mat_(1,4) != -8 ||
          mat_(2,0) != 0 || mat_(2,1) != 0 || mat_(2,2) != -15 || mat_(2,3) !=  5 || mat_(2,4) !=  9 ||
          mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) !=   0 || mat_(3,3) != -6 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0   0  0  7 )\n"
                                     "( 0  1   0  4 -8 )\n"
                                     "( 0  0 -15  5  9 )\n"
                                     "( 0  0   0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major dense vector multiplication assignment
   //=====================================================================================

   {
      test_ = "Row-major dense vector multiplication assignment";

      initialize();

      CT col2 = blaze::column( mat_, 2UL );

      blaze::DynamicVector<int,blaze::columnVector> vec{ 2, -4, 0, 0 };

      col2 *= vec;

      checkSize    ( col2, 4UL );
      checkNonZeros( col2, 1UL );
      checkRows    ( mat_, 4UL );
      checkColumns ( mat_, 5UL );
      checkNonZeros( mat_, 9UL );

      if( col2[0] != -4 || col2[1] != 0 || col2[2] != 0 || col2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( -4 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) != -4 || mat_(0,3) !=  0 || mat_(0,4) !=  7 ||
          mat_(1,0) != 0 || mat_(1,1) != 1 || mat_(1,2) !=  0 || mat_(1,3) !=  4 || mat_(1,4) != -8 ||
          mat_(2,0) != 0 || mat_(2,1) != 0 || mat_(2,2) !=  0 || mat_(2,3) !=  5 || mat_(2,4) !=  9 ||
          mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) != -6 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0 -4  0  7 )\n"
                                     "( 0  1  0  4 -8 )\n"
                                     "( 0  0  0  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse vector multiplication assignment
   //=====================================================================================

   {
      test_ = "Row-major sparse vector multiplication assignment";

      initialize();

      CT col2 = blaze::column( mat_, 2UL );

      blaze::CompressedVector<int,blaze::columnVector> vec( 4UL );
      vec[0] =  2;
      vec[1] = -4;

      col2 *= vec;

      checkSize    ( col2, 4UL );
      checkNonZeros( col2, 1UL );
      checkRows    ( mat_, 4UL );
      checkColumns ( mat_, 5UL );
      checkNonZeros( mat_, 9UL );

      if( col2[0] != -4 || col2[1] != 0 || col2[2] != 0 || col2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( -4 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) != -4 || mat_(0,3) !=  0 || mat_(0,4) !=  7 ||
          mat_(1,0) != 0 || mat_(1,1) != 1 || mat_(1,2) !=  0 || mat_(1,3) !=  4 || mat_(1,4) != -8 ||
          mat_(2,0) != 0 || mat_(2,1) != 0 || mat_(2,2) !=  0 || mat_(2,3) !=  5 || mat_(2,4) !=  9 ||
          mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) != -6 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0 -4  0  7 )\n"
                                     "( 0  1  0  4 -8 )\n"
                                     "( 0  0  0  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major Column multiplication assignment
   //=====================================================================================

   {
      test_ = "Column-major Column multiplication assignment";

      initialize();

      OCT col2 = blaze::column( tmat_, 2UL );
      col2 *= blaze::column( tmat_, 3UL );

      checkSize    ( col2 , 4UL );
      checkNonZeros( col2 , 1UL );
      checkRows    ( tmat_, 4UL );
      checkColumns ( tmat_, 5UL );
      checkNonZeros( tmat_, 9UL );

      if( col2[0] != 0 || col2[1] != 0 || col2[2] != -15 || col2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( 0 0 -15 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) !=   0 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) !=   0 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != -15 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=   0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0   0  0  7 )\n"
                                     "( 0  1   0  4 -8 )\n"
                                     "( 0  0 -15  5  9 )\n"
                                     "( 0  0   0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense vector multiplication assignment
   //=====================================================================================

   {
      test_ = "Column-major dense vector multiplication assignment";

      initialize();

      OCT col2 = blaze::column( tmat_, 2UL );

      blaze::DynamicVector<int,blaze::columnVector> vec{ 2, -4, 0, 0 };

      col2 *= vec;

      checkSize    ( col2 , 4UL );
      checkNonZeros( col2 , 1UL );
      checkRows    ( tmat_, 4UL );
      checkColumns ( tmat_, 5UL );
      checkNonZeros( tmat_, 9UL );

      if( col2[0] != -4 || col2[1] != 0 || col2[2] != 0 || col2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( -4 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -4 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) !=  0 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) !=  0 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -4  0  7 )\n"
                                     "( 0  1  0  4 -8 )\n"
                                     "( 0  0  0  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse vector multiplication assignment
   //=====================================================================================

   {
      test_ = "Column-major sparse vector multiplication assignment";

      initialize();

      OCT col2 = blaze::column( tmat_, 2UL );

      blaze::CompressedVector<int,blaze::columnVector> vec( 4UL );
      vec[0] =  2;
      vec[1] = -4;

      col2 *= vec;

      checkSize    ( col2 , 4UL );
      checkNonZeros( col2 , 1UL );
      checkRows    ( tmat_, 4UL );
      checkColumns ( tmat_, 5UL );
      checkNonZeros( tmat_, 9UL );

      if( col2[0] != -4 || col2[1] != 0 || col2[2] != 0 || col2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( -4 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -4 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) !=  0 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) !=  0 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -4  0  7 )\n"
                                     "( 0  1  0  4 -8 )\n"
                                     "( 0  0  0  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Column division assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the division assignment operators of the Column
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseGeneralTest::testDivAssign()
{
   //=====================================================================================
   // Row-major dense vector division assignment
   //=====================================================================================

   {
      test_ = "Row-major dense vector division assignment";

      initialize();

      CT col2 = blaze::column( mat_, 2UL );

      blaze::DynamicVector<int,blaze::columnVector> vec{ -1, 2, 3, 4 };

      col2 /= vec;

      checkSize    ( col2,  4UL );
      checkNonZeros( col2,  2UL );
      checkRows    ( mat_,  4UL );
      checkColumns ( mat_,  5UL );
      checkNonZeros( mat_, 10UL );

      if( col2[0] != 2 || col2[1] != 0 || col2[2] != -1 || col2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( 2 0 -1 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) !=  2 || mat_(0,3) !=  0 || mat_(0,4) !=  7 ||
          mat_(1,0) != 0 || mat_(1,1) != 1 || mat_(1,2) !=  0 || mat_(1,3) !=  4 || mat_(1,4) != -8 ||
          mat_(2,0) != 0 || mat_(2,1) != 0 || mat_(2,2) != -1 || mat_(2,3) !=  5 || mat_(2,4) !=  9 ||
          mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) != -6 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0  2  0  7 )\n"
                                     "( 0  1  0  4 -8 )\n"
                                     "( 0  0 -1  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense vector division assignment
   //=====================================================================================

   {
      test_ = "Column-major dense vector division assignment";

      initialize();

      OCT col2 = blaze::column( tmat_, 2UL );

      blaze::DynamicVector<int,blaze::columnVector> vec{ -1, 2, 3, 4 };

      col2 /= vec;

      checkSize    ( col2 ,  4UL );
      checkNonZeros( col2 ,  2UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 10UL );

      if( col2[0] != 2 || col2[1] != 0 || col2[2] != -1 || col2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( 2 0 -1 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) !=  2 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) !=  0 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != -1 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0  2  0  7 )\n"
                                     "( 0  1  0  4 -8 )\n"
                                     "( 0  0 -1  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Column cross product assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the cross product assignment operators of the Column
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseGeneralTest::testCrossAssign()
{
   //=====================================================================================
   // Row-major Column cross product assignment
   //=====================================================================================

   {
      test_ = "Row-major Column cross product assignment";

      MT mat( 3UL, 2UL, 4UL );
      mat(0,0) =  2;
      mat(0,1) =  1;
      mat(2,0) = -1;
      mat(2,1) = -2;

      CT col0 = blaze::column( mat, 0UL );
      col0 %= blaze::column( mat, 1UL );

      checkSize    ( col0, 3UL );
      checkNonZeros( col0, 1UL );
      checkRows    ( mat , 3UL );
      checkColumns ( mat , 2UL );
      checkNonZeros( mat , 3UL );

      if( col0[0] != 0 || col0[1] != 3 || col0[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cross product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << col0 << "\n"
             << "   Expected result:\n( 0 3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat(0,0) != 0 || mat(0,1) !=  1 ||
          mat(1,0) != 3 || mat(1,1) !=  0 ||
          mat(2,0) != 0 || mat(2,1) != -2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cross product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0  1 )\n"
                                     "( 3  0 )\n"
                                     "( 0 -2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major dense vector cross product assignment
   //=====================================================================================

   {
      test_ = "Row-major dense vector cross product assignment";

      MT mat( 3UL, 2UL, 4UL );
      mat(0,0) =  2;
      mat(0,1) =  1;
      mat(2,0) = -1;
      mat(2,1) = -2;

      CT col0 = blaze::column( mat, 0UL );

      const blaze::DynamicVector<int,blaze::columnVector> vec{ 1, 0, -2 };

      col0 %= vec;

      checkSize    ( col0, 3UL );
      checkNonZeros( col0, 1UL );
      checkRows    ( mat , 3UL );
      checkColumns ( mat , 2UL );
      checkNonZeros( mat , 3UL );

      if( col0[0] != 0 || col0[1] != 3 || col0[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cross product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << col0 << "\n"
             << "   Expected result:\n( 0 3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat(0,0) != 0 || mat(0,1) !=  1 ||
          mat(1,0) != 3 || mat(1,1) !=  0 ||
          mat(2,0) != 0 || mat(2,1) != -2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cross product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0  1 )\n"
                                     "( 3  0 )\n"
                                     "( 0 -2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse vector cross product assignment
   //=====================================================================================

   {
      test_ = "Row-major sparse vector cross product assignment";

      MT mat( 3UL, 2UL, 4UL );
      mat(0,0) =  2;
      mat(0,1) =  1;
      mat(2,0) = -1;
      mat(2,1) = -2;

      CT col0 = blaze::column( mat, 0UL );

      blaze::CompressedVector<int,blaze::columnVector> vec( 3UL );
      vec[0] =  1;
      vec[2] = -2;

      col0 %= vec;

      checkSize    ( col0, 3UL );
      checkNonZeros( col0, 1UL );
      checkRows    ( mat , 3UL );
      checkColumns ( mat , 2UL );
      checkNonZeros( mat , 3UL );

      if( col0[0] != 0 || col0[1] != 3 || col0[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cross product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << col0 << "\n"
             << "   Expected result:\n( 0 3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat(0,0) != 0 || mat(0,1) !=  1 ||
          mat(1,0) != 3 || mat(1,1) !=  0 ||
          mat(2,0) != 0 || mat(2,1) != -2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cross product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0  1 )\n"
                                     "( 3  0 )\n"
                                     "( 0 -2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major Column cross product assignment
   //=====================================================================================

   {
      test_ = "Column-major Column cross product assignment";

      OMT mat( 3UL, 2UL, 4UL );
      mat(0,0) =  2;
      mat(0,1) =  1;
      mat(2,0) = -1;
      mat(2,1) = -2;

      OCT col0 = blaze::column( mat, 0UL );
      col0 %= blaze::column( mat, 1UL );

      checkSize    ( col0, 3UL );
      checkNonZeros( col0, 1UL );
      checkRows    ( mat , 3UL );
      checkColumns ( mat , 2UL );
      checkNonZeros( mat , 3UL );

      if( col0[0] != 0 || col0[1] != 3 || col0[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cross product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << col0 << "\n"
             << "   Expected result:\n( 0 3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat(0,0) != 0 || mat(0,1) !=  1 ||
          mat(1,0) != 3 || mat(1,1) !=  0 ||
          mat(2,0) != 0 || mat(2,1) != -2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cross product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0  1 )\n"
                                     "( 3  0 )\n"
                                     "( 0 -2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense vector cross product assignment
   //=====================================================================================

   {
      test_ = "Column-major dense vector cross product assignment";

      OMT mat( 3UL, 2UL, 4UL );
      mat(0,0) =  2;
      mat(0,1) =  1;
      mat(2,0) = -1;
      mat(2,1) = -2;

      OCT col0 = blaze::column( mat, 0UL );

      const blaze::DynamicVector<int,blaze::columnVector> vec{ 1, 0, -2 };

      col0 %= vec;

      checkSize    ( col0, 3UL );
      checkNonZeros( col0, 1UL );
      checkRows    ( mat , 3UL );
      checkColumns ( mat , 2UL );
      checkNonZeros( mat , 3UL );

      if( col0[0] != 0 || col0[1] != 3 || col0[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cross product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << col0 << "\n"
             << "   Expected result:\n( 0 3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat(0,0) != 0 || mat(0,1) !=  1 ||
          mat(1,0) != 3 || mat(1,1) !=  0 ||
          mat(2,0) != 0 || mat(2,1) != -2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cross product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0  1 )\n"
                                     "( 3  0 )\n"
                                     "( 0 -2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse vector cross product assignment
   //=====================================================================================

   {
      test_ = "Column-major sparse vector cross product assignment";

      OMT mat( 3UL, 2UL, 4UL );
      mat(0,0) =  2;
      mat(0,1) =  1;
      mat(2,0) = -1;
      mat(2,1) = -2;

      OCT col0 = blaze::column( mat, 0UL );

      blaze::CompressedVector<int,blaze::columnVector> vec( 3UL );
      vec[0] =  1;
      vec[2] = -2;

      col0 %= vec;

      checkSize    ( col0, 3UL );
      checkNonZeros( col0, 1UL );
      checkRows    ( mat , 3UL );
      checkColumns ( mat , 2UL );
      checkNonZeros( mat , 3UL );

      if( col0[0] != 0 || col0[1] != 3 || col0[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cross product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << col0 << "\n"
             << "   Expected result:\n( 0 3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat(0,0) != 0 || mat(0,1) !=  1 ||
          mat(1,0) != 3 || mat(1,1) !=  0 ||
          mat(2,0) != 0 || mat(2,1) != -2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cross product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0  1 )\n"
                                     "( 3  0 )\n"
                                     "( 0 -2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of all Column (self-)scaling operations.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all available ways to scale an instance of the Column
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseGeneralTest::testScaling()
{
   //=====================================================================================
   // Row-major self-scaling (v*=s)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (v*=s)";

      initialize();

      CT col2 = blaze::column( mat_, 2UL );

      col2 *= 3;

      checkSize    ( col2,  4UL );
      checkNonZeros( col2,  2UL );
      checkRows    ( mat_,  4UL );
      checkColumns ( mat_,  5UL );
      checkNonZeros( mat_, 10UL );

      if( col2[0] != -6 || col2[1] != 0 || col2[2] != -9 || col2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( -6 0 -9 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) != -6 || mat_(0,3) !=  0 || mat_(0,4) !=  7 ||
          mat_(1,0) != 0 || mat_(1,1) != 1 || mat_(1,2) !=  0 || mat_(1,3) !=  4 || mat_(1,4) != -8 ||
          mat_(2,0) != 0 || mat_(2,1) != 0 || mat_(2,2) != -9 || mat_(2,3) !=  5 || mat_(2,4) !=  9 ||
          mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) != -6 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0 -6  0  7 )\n"
                                     "( 0  1  0  4 -8 )\n"
                                     "( 0  0 -9  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major self-scaling (v=v*s)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (v=v*s)";

      initialize();

      CT col2 = blaze::column( mat_, 2UL );

      col2 = col2 * 3;

      checkSize    ( col2,  4UL );
      checkNonZeros( col2,  2UL );
      checkRows    ( mat_,  4UL );
      checkColumns ( mat_,  5UL );
      checkNonZeros( mat_, 10UL );

      if( col2[0] != -6 || col2[1] != 0 || col2[2] != -9 || col2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( -6 0 -9 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) != -6 || mat_(0,3) !=  0 || mat_(0,4) !=  7 ||
          mat_(1,0) != 0 || mat_(1,1) != 1 || mat_(1,2) !=  0 || mat_(1,3) !=  4 || mat_(1,4) != -8 ||
          mat_(2,0) != 0 || mat_(2,1) != 0 || mat_(2,2) != -9 || mat_(2,3) !=  5 || mat_(2,4) !=  9 ||
          mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) != -6 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0 -6  0  7 )\n"
                                     "( 0  1  0  4 -8 )\n"
                                     "( 0  0 -9  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major self-scaling (v=s*v)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (v=s*v)";

      initialize();

      CT col2 = blaze::column( mat_, 2UL );

      col2 = 3 * col2;

      checkSize    ( col2,  4UL );
      checkNonZeros( col2,  2UL );
      checkRows    ( mat_,  4UL );
      checkColumns ( mat_,  5UL );
      checkNonZeros( mat_, 10UL );

      if( col2[0] != -6 || col2[1] != 0 || col2[2] != -9 || col2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( -6 0 -9 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) != -6 || mat_(0,3) !=  0 || mat_(0,4) !=  7 ||
          mat_(1,0) != 0 || mat_(1,1) != 1 || mat_(1,2) !=  0 || mat_(1,3) !=  4 || mat_(1,4) != -8 ||
          mat_(2,0) != 0 || mat_(2,1) != 0 || mat_(2,2) != -9 || mat_(2,3) !=  5 || mat_(2,4) !=  9 ||
          mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) != -6 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0 -6  0  7 )\n"
                                     "( 0  1  0  4 -8 )\n"
                                     "( 0  0 -9  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major self-scaling (v/=s)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (v/=s)";

      initialize();

      CT col2 = blaze::column( mat_, 2UL );

      col2 /= 0.5;

      checkSize    ( col2,  4UL );
      checkNonZeros( col2,  2UL );
      checkRows    ( mat_,  4UL );
      checkColumns ( mat_,  5UL );
      checkNonZeros( mat_, 10UL );

      if( col2[0] != -4 || col2[1] != 0 || col2[2] != -6 || col2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( -4 0 -6 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) != -4 || mat_(0,3) !=  0 || mat_(0,4) !=  7 ||
          mat_(1,0) != 0 || mat_(1,1) != 1 || mat_(1,2) !=  0 || mat_(1,3) !=  4 || mat_(1,4) != -8 ||
          mat_(2,0) != 0 || mat_(2,1) != 0 || mat_(2,2) != -6 || mat_(2,3) !=  5 || mat_(2,4) !=  9 ||
          mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) != -6 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0 -4  0  7 )\n"
                                     "( 0  1  0  4 -8 )\n"
                                     "( 0  0 -6  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major self-scaling (v=v/s)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (v=v/s)";

      initialize();

      CT col2 = blaze::column( mat_, 2UL );

      col2 = col2 / 0.5;

      checkSize    ( col2,  4UL );
      checkNonZeros( col2,  2UL );
      checkRows    ( mat_,  4UL );
      checkColumns ( mat_,  5UL );
      checkNonZeros( mat_, 10UL );

      if( col2[0] != -4 || col2[1] != 0 || col2[2] != -6 || col2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( -4 0 -6 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) != -4 || mat_(0,3) !=  0 || mat_(0,4) !=  7 ||
          mat_(1,0) != 0 || mat_(1,1) != 1 || mat_(1,2) !=  0 || mat_(1,3) !=  4 || mat_(1,4) != -8 ||
          mat_(2,0) != 0 || mat_(2,1) != 0 || mat_(2,2) != -6 || mat_(2,3) !=  5 || mat_(2,4) !=  9 ||
          mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) != -6 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0 -4  0  7 )\n"
                                     "( 0  1  0  4 -8 )\n"
                                     "( 0  0 -6  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major Column::scale()
   //=====================================================================================

   {
      test_ = "Row-major Column::scale()";

      initialize();

      // Integral scaling the 3rd column
      {
         CT col3 = blaze::column( mat_, 3UL );
         col3.scale( 3 );

         checkSize    ( col3,  4UL );
         checkNonZeros( col3,  3UL );
         checkRows    ( mat_,  4UL );
         checkColumns ( mat_,  5UL );
         checkNonZeros( mat_, 10UL );

         if( col3[0] != 0 || col3[1] != 12 || col3[2] != 15 || col3[3] != -18 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Integral scale operation of 3rd column failed\n"
                << " Details:\n"
                << "   Result:\n" << col3 << "\n"
                << "   Expected result:\n( 0 12 15 -18 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) != -2 || mat_(0,3) !=   0 || mat_(0,4) !=  7 ||
             mat_(1,0) != 0 || mat_(1,1) != 1 || mat_(1,2) !=  0 || mat_(1,3) !=  12 || mat_(1,4) != -8 ||
             mat_(2,0) != 0 || mat_(2,1) != 0 || mat_(2,2) != -3 || mat_(2,3) !=  15 || mat_(2,4) !=  9 ||
             mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) != -18 || mat_(3,4) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Integral scale operation of 3rd column failed\n"
                << " Details:\n"
                << "   Result:\n" << mat_ << "\n"
                << "   Expected result:\n( 0  0 -2   0  7 )\n"
                                        "( 0  1  0  12 -8 )\n"
                                        "( 0  0 -3  15  9 )\n"
                                        "( 0  0  0 -18 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Floating point scaling the 3rd column
      {
         CT col3 = blaze::column( mat_, 3UL );
         col3.scale( 0.5 );

         checkSize    ( col3,  4UL );
         checkNonZeros( col3,  3UL );
         checkRows    ( mat_,  4UL );
         checkColumns ( mat_,  5UL );
         checkNonZeros( mat_, 10UL );

         if( col3[0] != 0 || col3[1] != 6 || col3[2] != 7 || col3[3] != -9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Floating point scale operation of 3rd column failed\n"
                << " Details:\n"
                << "   Result:\n" << col3 << "\n"
                << "   Expected result:\n( 0 6 7 -9 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) != -2 || mat_(0,3) !=  0 || mat_(0,4) !=  7 ||
             mat_(1,0) != 0 || mat_(1,1) != 1 || mat_(1,2) !=  0 || mat_(1,3) !=  6 || mat_(1,4) != -8 ||
             mat_(2,0) != 0 || mat_(2,1) != 0 || mat_(2,2) != -3 || mat_(2,3) !=  7 || mat_(2,4) !=  9 ||
             mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) != -9 || mat_(3,4) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Floating point scale operation of 3rd column failed\n"
                << " Details:\n"
                << "   Result:\n" << mat_ << "\n"
                << "   Expected result:\n( 0  0 -2  0  7 )\n"
                                        "( 0  1  0  6 -8 )\n"
                                        "( 0  0 -3  7  9 )\n"
                                        "( 0  0  0 -9 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major self-scaling (v*=s)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (v*=s)";

      initialize();

      OCT col2 = blaze::column( tmat_, 2UL );

      col2 *= 3;

      checkSize    ( col2 ,  4UL );
      checkNonZeros( col2 ,  2UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 10UL );

      if( col2[0] != -6 || col2[1] != 0 || col2[2] != -9 || col2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( -6 0 -9 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -6 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) !=  0 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != -9 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -6  0  7 )\n"
                                     "( 0  1  0  4 -8 )\n"
                                     "( 0  0 -9  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (v=v*s)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (v=v*s)";

      initialize();

      OCT col2 = blaze::column( tmat_, 2UL );

      col2 = col2 * 3;

      checkSize    ( col2 ,  4UL );
      checkNonZeros( col2 ,  2UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 10UL );

      if( col2[0] != -6 || col2[1] != 0 || col2[2] != -9 || col2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( -6 0 -9 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -6 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) !=  0 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != -9 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -6  0  7 )\n"
                                     "( 0  1  0  4 -8 )\n"
                                     "( 0  0 -9  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (v=s*v)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (v=s*v)";

      initialize();

      OCT col2 = blaze::column( tmat_, 2UL );

      col2 = 3 * col2;

      checkSize    ( col2 ,  4UL );
      checkNonZeros( col2 ,  2UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 10UL );

      if( col2[0] != -6 || col2[1] != 0 || col2[2] != -9 || col2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( -6 0 -9 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -6 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) !=  0 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != -9 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -6  0  7 )\n"
                                     "( 0  1  0  4 -8 )\n"
                                     "( 0  0 -9  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (v/=s)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (v/=s)";

      initialize();

      OCT col2 = blaze::column( tmat_, 2UL );

      col2 /= 0.5;

      checkSize    ( col2 ,  4UL );
      checkNonZeros( col2 ,  2UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 10UL );

      if( col2[0] != -4 || col2[1] != 0 || col2[2] != -6 || col2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( -4 0 -6 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -4 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) !=  0 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != -6 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -4  0  7 )\n"
                                     "( 0  1  0  4 -8 )\n"
                                     "( 0  0 -6  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (v=v/s)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (v=v/s)";

      initialize();

      OCT col2 = blaze::column( tmat_, 2UL );

      col2 = col2 / 0.5;

      checkSize    ( col2 ,  4UL );
      checkNonZeros( col2 ,  2UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 10UL );

      if( col2[0] != -4 || col2[1] != 0 || col2[2] != -6 || col2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( -4 0 -6 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -4 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) !=  0 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != -6 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -4  0  7 )\n"
                                     "( 0  1  0  4 -8 )\n"
                                     "( 0  0 -6  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major Column::scale()
   //=====================================================================================

   {
      test_ = "Column-major Column::scale()";

      initialize();

      // Integral scaling the 3rd column
      {
         OCT col3 = blaze::column( tmat_, 3UL );
         col3.scale( 3 );

         checkSize    ( col3 ,  4UL );
         checkNonZeros( col3 ,  3UL );
         checkRows    ( tmat_,  4UL );
         checkColumns ( tmat_,  5UL );
         checkNonZeros( tmat_, 10UL );

         if( col3[0] != 0 || col3[1] != 12 || col3[2] != 15 || col3[3] != -18 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Integral scale operation of 3rd column failed\n"
                << " Details:\n"
                << "   Result:\n" << col3 << "\n"
                << "   Expected result:\n( 0 12 15 -18 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -2 || tmat_(0,3) !=   0 || tmat_(0,4) !=  7 ||
             tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) !=  0 || tmat_(1,3) !=  12 || tmat_(1,4) != -8 ||
             tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != -3 || tmat_(2,3) !=  15 || tmat_(2,4) !=  9 ||
             tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != -18 || tmat_(3,4) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Integral scale operation of 3rd column failed\n"
                << " Details:\n"
                << "   Result:\n" << tmat_ << "\n"
                << "   Expected result:\n( 0  0 -2   0  7 )\n"
                                        "( 0  1  0  12 -8 )\n"
                                        "( 0  0 -3  15  9 )\n"
                                        "( 0  0  0 -18 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Floating point scaling the 3rd column
      {
         OCT col3 = blaze::column( tmat_, 3UL );
         col3.scale( 0.5 );

         checkSize    ( col3 ,  4UL );
         checkNonZeros( col3 ,  3UL );
         checkRows    ( tmat_,  4UL );
         checkColumns ( tmat_,  5UL );
         checkNonZeros( tmat_, 10UL );

         if( col3[0] != 0 || col3[1] != 6 || col3[2] != 7 || col3[3] != -9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Floating point scale operation of 3rd column failed\n"
                << " Details:\n"
                << "   Result:\n" << col3 << "\n"
                << "   Expected result:\n( 0 6 7 -9 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -2 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
             tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) !=  0 || tmat_(1,3) !=  6 || tmat_(1,4) != -8 ||
             tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != -3 || tmat_(2,3) !=  7 || tmat_(2,4) !=  9 ||
             tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != -9 || tmat_(3,4) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Floating point scale operation of 3rd column failed\n"
                << " Details:\n"
                << "   Result:\n" << tmat_ << "\n"
                << "   Expected result:\n( 0  0 -2  0  7 )\n"
                                        "( 0  1  0  6 -8 )\n"
                                        "( 0  0 -3  7  9 )\n"
                                        "( 0  0  0 -9 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Column subscript operator.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of adding and accessing elements via the subscript operator
// of the Column specialization. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
void SparseGeneralTest::testSubscript()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major Column::operator[]";

      initialize();

      CT col2 = blaze::column( mat_, 2UL );

      // Assignment to the element at index 1
      col2[1] = 9;

      checkSize    ( col2, 4UL );
      checkNonZeros( col2, 3UL );

      if( col2[0] != -2 || col2[1] != 9 || col2[2] != -3 || col2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( -2 9 -3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) != -2 || mat_(0,3) !=  0 || mat_(0,4) !=  7 ||
          mat_(1,0) != 0 || mat_(1,1) != 1 || mat_(1,2) !=  9 || mat_(1,3) !=  4 || mat_(1,4) != -8 ||
          mat_(2,0) != 0 || mat_(2,1) != 0 || mat_(2,2) != -3 || mat_(2,3) !=  5 || mat_(2,4) !=  9 ||
          mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) != -6 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0 -2  0  7 )\n"
                                     "( 0  1  9  4 -8 )\n"
                                     "( 0  0 -3  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Assignment to the element at index 2
      col2[2] = 0;

      checkSize    ( col2, 4UL );
      checkNonZeros( col2, 2UL );

      if( col2[0] != -2 || col2[1] != 9 || col2[2] != 0 || col2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( -2 9 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) != -2 || mat_(0,3) !=  0 || mat_(0,4) !=  7 ||
          mat_(1,0) != 0 || mat_(1,1) != 1 || mat_(1,2) !=  9 || mat_(1,3) !=  4 || mat_(1,4) != -8 ||
          mat_(2,0) != 0 || mat_(2,1) != 0 || mat_(2,2) !=  0 || mat_(2,3) !=  5 || mat_(2,4) !=  9 ||
          mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) != -6 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0 -2  0  7 )\n"
                                     "( 0  1  9  4 -8 )\n"
                                     "( 0  0  0  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Assignment to the element at index 3
      col2[3] = -8;

      checkSize    ( col2, 4UL );
      checkNonZeros( col2, 3UL );

      if( col2[0] != -2 || col2[1] != 9 || col2[2] != 0 || col2[3] != -8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( -2 9 0 -8 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) != -2 || mat_(0,3) !=  0 || mat_(0,4) !=  7 ||
          mat_(1,0) != 0 || mat_(1,1) != 1 || mat_(1,2) !=  9 || mat_(1,3) !=  4 || mat_(1,4) != -8 ||
          mat_(2,0) != 0 || mat_(2,1) != 0 || mat_(2,2) !=  0 || mat_(2,3) !=  5 || mat_(2,4) !=  9 ||
          mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) != -8 || mat_(3,3) != -6 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0 -2  0  7 )\n"
                                     "( 0  1  9  4 -8 )\n"
                                     "( 0  0  0  5  9 )\n"
                                     "( 0  0 -8 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Addition assignment to the element at index 0
      col2[0] += -3;

      checkSize    ( col2, 4UL );
      checkCapacity( col2, 4UL );
      checkNonZeros( col2, 3UL );

      if( col2[0] != -5 || col2[1] != 9 || col2[2] != 0 || col2[3] != -8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( -5 9 0 -8 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) != -5 || mat_(0,3) !=  0 || mat_(0,4) !=  7 ||
          mat_(1,0) != 0 || mat_(1,1) != 1 || mat_(1,2) !=  9 || mat_(1,3) !=  4 || mat_(1,4) != -8 ||
          mat_(2,0) != 0 || mat_(2,1) != 0 || mat_(2,2) !=  0 || mat_(2,3) !=  5 || mat_(2,4) !=  9 ||
          mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) != -8 || mat_(3,3) != -6 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0 -5  0  7 )\n"
                                     "( 0  1  9  4 -8 )\n"
                                     "( 0  0  0  5  9 )\n"
                                     "( 0  0 -8 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Subtraction assignment to the element at index 1
      col2[1] -= 6;

      checkSize    ( col2, 4UL );
      checkCapacity( col2, 4UL );
      checkNonZeros( col2, 3UL );

      if( col2[0] != -5 || col2[1] != 3 || col2[2] != 0 || col2[3] != -8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( -5 3 0 -8 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) != -5 || mat_(0,3) !=  0 || mat_(0,4) !=  7 ||
          mat_(1,0) != 0 || mat_(1,1) != 1 || mat_(1,2) !=  3 || mat_(1,3) !=  4 || mat_(1,4) != -8 ||
          mat_(2,0) != 0 || mat_(2,1) != 0 || mat_(2,2) !=  0 || mat_(2,3) !=  5 || mat_(2,4) !=  9 ||
          mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) != -8 || mat_(3,3) != -6 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0 -5  0  7 )\n"
                                     "( 0  1  3  4 -8 )\n"
                                     "( 0  0  0  5  9 )\n"
                                     "( 0  0 -8 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Multiplication assignment to the element at index 1
      col2[1] *= -3;

      checkSize    ( col2, 4UL );
      checkCapacity( col2, 4UL );
      checkNonZeros( col2, 3UL );

      if( col2[0] != -5 || col2[1] != -9 || col2[2] != 0 || col2[3] != -8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( -5 -9 0 -8 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) != -5 || mat_(0,3) !=  0 || mat_(0,4) !=  7 ||
          mat_(1,0) != 0 || mat_(1,1) != 1 || mat_(1,2) != -9 || mat_(1,3) !=  4 || mat_(1,4) != -8 ||
          mat_(2,0) != 0 || mat_(2,1) != 0 || mat_(2,2) !=  0 || mat_(2,3) !=  5 || mat_(2,4) !=  9 ||
          mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) != -8 || mat_(3,3) != -6 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0 -5  0  7 )\n"
                                     "( 0  1 -9  4 -8 )\n"
                                     "( 0  0  0  5  9 )\n"
                                     "( 0  0 -8 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Division assignment to the element at index 3
      col2[3] /= 2;

      checkSize    ( col2, 4UL );
      checkCapacity( col2, 4UL );
      checkNonZeros( col2, 3UL );

      if( col2[0] != -5 || col2[1] != -9 || col2[2] != 0 || col2[3] != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( -5 -9 0 -4 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) != -5 || mat_(0,3) !=  0 || mat_(0,4) !=  7 ||
          mat_(1,0) != 0 || mat_(1,1) != 1 || mat_(1,2) != -9 || mat_(1,3) !=  4 || mat_(1,4) != -8 ||
          mat_(2,0) != 0 || mat_(2,1) != 0 || mat_(2,2) !=  0 || mat_(2,3) !=  5 || mat_(2,4) !=  9 ||
          mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) != -4 || mat_(3,3) != -6 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0 -5  0  7 )\n"
                                     "( 0  1 -9  4 -8 )\n"
                                     "( 0  0  0  5  9 )\n"
                                     "( 0  0 -4 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major Column::operator[]";

      initialize();

      OCT col2 = blaze::column( tmat_, 2UL );

      // Assignment to the element at index 1
      col2[1] = 9;

      checkSize    ( col2, 4UL );
      checkNonZeros( col2, 3UL );

      if( col2[0] != -2 || col2[1] != 9 || col2[2] != -3 || col2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( -2 9 -3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -2 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) !=  9 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != -3 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -2  0  7 )\n"
                                     "( 0  1  9  4 -8 )\n"
                                     "( 0  0 -3  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Assignment to the element at index 2
      col2[2] = 0;

      checkSize    ( col2, 4UL );
      checkNonZeros( col2, 2UL );

      if( col2[0] != -2 || col2[1] != 9 || col2[2] != 0 || col2[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( -2 9 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -2 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) !=  9 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) !=  0 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -2  0  7 )\n"
                                     "( 0  1  9  4 -8 )\n"
                                     "( 0  0  0  5  9 )\n"
                                     "( 0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Assignment to the element at index 3
      col2[3] = -8;

      checkSize    ( col2, 4UL );
      checkNonZeros( col2, 3UL );

      if( col2[0] != -2 || col2[1] != 9 || col2[2] != 0 || col2[3] != -8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( -2 9 0 -8 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -2 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) !=  9 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) !=  0 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) != -8 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -2  0  7 )\n"
                                     "( 0  1  9  4 -8 )\n"
                                     "( 0  0  0  5  9 )\n"
                                     "( 0  0 -8 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Addition assignment to the element at index 0
      col2[0] += -3;

      checkSize    ( col2, 4UL );
      checkCapacity( col2, 4UL );
      checkNonZeros( col2, 3UL );

      if( col2[0] != -5 || col2[1] != 9 || col2[2] != 0 || col2[3] != -8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( -5 9 0 -8 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -5 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) !=  9 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) !=  0 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) != -8 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -5  0  7 )\n"
                                     "( 0  1  9  4 -8 )\n"
                                     "( 0  0  0  5  9 )\n"
                                     "( 0  0 -8 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Subtraction assignment to the element at index 1
      col2[1] -= 6;

      checkSize    ( col2, 4UL );
      checkCapacity( col2, 4UL );
      checkNonZeros( col2, 3UL );

      if( col2[0] != -5 || col2[1] != 3 || col2[2] != 0 || col2[3] != -8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( -5 3 0 -8 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -5 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) !=  3 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) !=  0 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) != -8 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -5  0  7 )\n"
                                     "( 0  1  3  4 -8 )\n"
                                     "( 0  0  0  5  9 )\n"
                                     "( 0  0 -8 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Multiplication assignment to the element at index 1
      col2[1] *= -3;

      checkSize    ( col2, 4UL );
      checkCapacity( col2, 4UL );
      checkNonZeros( col2, 3UL );

      if( col2[0] != -5 || col2[1] != -9 || col2[2] != 0 || col2[3] != -8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( -5 -9 0 -8 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -5 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) != -9 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) !=  0 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) != -8 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -5  0  7 )\n"
                                     "( 0  1 -9  4 -8 )\n"
                                     "( 0  0  0  5  9 )\n"
                                     "( 0  0 -8 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Division assignment to the element at index 3
      col2[3] /= 2;

      checkSize    ( col2, 4UL );
      checkCapacity( col2, 4UL );
      checkNonZeros( col2, 3UL );

      if( col2[0] != -5 || col2[1] != -9 || col2[2] != 0 || col2[3] != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << col2 << "\n"
             << "   Expected result:\n( -5 -9 0 -4 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -5 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) != -9 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) !=  0 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) != -4 || tmat_(3,3) != -6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -5  0  7 )\n"
                                     "( 0  1 -9  4 -8 )\n"
                                     "( 0  0  0  5  9 )\n"
                                     "( 0  0 -4 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Column iterator implementation.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the iterator implementation of the Column specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseGeneralTest::testIterator()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      initialize();

      // Testing the Iterator default constructor
      {
         test_ = "Row-major Iterator default constructor";

         CT::Iterator it{};

         if( it != CT::Iterator() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Failed iterator default constructor\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing the ConstIterator default constructor
      {
         test_ = "Row-major ConstIterator default constructor";

         CT::ConstIterator it{};

         if( it != CT::ConstIterator() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Failed iterator default constructor\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing conversion from Iterator to ConstIterator
      {
         test_ = "Row-major Iterator/ConstIterator conversion";

         CT col2 = blaze::column( mat_, 2UL );
         CT::ConstIterator it( begin( col2 ) );

         if( it == end( col2 ) || it->value() != -2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Failed iterator conversion detected\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Counting the number of elements in 1st column via Iterator (end-begin)
      {
         test_ = "Row-major Iterator subtraction (end-begin)";

         CT col1 = blaze::column( mat_, 1UL );
         const ptrdiff_t number( end( col1 ) - begin( col1 ) );

         if( number != 1L ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid number of elements detected\n"
                << " Details:\n"
                << "   Number of elements         : " << number << "\n"
                << "   Expected number of elements: 1\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Counting the number of elements in 2nd column via ConstIterator (end-begin)
      {
         test_ = "Row-major ConstIterator subtraction (end-begin)";

         CT col2 = blaze::column( mat_, 2UL );
         const ptrdiff_t number( cend( col2 ) - cbegin( col2 ) );

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
         test_ = "Row-major read-only access via ConstIterator";

         CT col2 = blaze::column( mat_, 2UL );
         CT::ConstIterator it ( cbegin( col2 ) );
         CT::ConstIterator end( cend( col2 ) );

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
         test_ = "Row-major assignment via Iterator";

         CT col4 = blaze::column( mat_, 4UL );
         int value = 6;

         for( CT::Iterator it=begin( col4 ); it!=end( col4 ); ++it ) {
            *it = value++;
         }

         if( col4[0] != 6 || col4[1] != 7 || col4[2] != 8 || col4[3] != 9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << col4 << "\n"
                << "   Expected result:\n( 6 7 8 9 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) != -2 || mat_(0,3) !=  0 || mat_(0,4) != 6 ||
             mat_(1,0) != 0 || mat_(1,1) != 1 || mat_(1,2) !=  0 || mat_(1,3) !=  4 || mat_(1,4) != 7 ||
             mat_(2,0) != 0 || mat_(2,1) != 0 || mat_(2,2) != -3 || mat_(2,3) !=  5 || mat_(2,4) != 8 ||
             mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) != -6 || mat_(3,4) != 9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << mat_ << "\n"
                << "   Expected result:\n( 0  0 -2  0  6 )\n"
                                        "( 0  1  0  4  7 )\n"
                                        "( 0  0 -3  5  8 )\n"
                                        "( 0  0  0 -6  9 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing addition assignment via Iterator
      {
         test_ = "Row-major addition assignment via Iterator";

         CT col4 = blaze::column( mat_, 4UL );
         int value = 2;

         for( CT::Iterator it=begin( col4 ); it!=end( col4 ); ++it ) {
            *it += value++;
         }

         if( col4[0] != 8 || col4[1] != 10 || col4[2] != 12 || col4[3] != 14 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Addition assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << col4 << "\n"
                << "   Expected result:\n( 8 10 12 14 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) != -2 || mat_(0,3) !=  0 || mat_(0,4) !=  8 ||
             mat_(1,0) != 0 || mat_(1,1) != 1 || mat_(1,2) !=  0 || mat_(1,3) !=  4 || mat_(1,4) != 10 ||
             mat_(2,0) != 0 || mat_(2,1) != 0 || mat_(2,2) != -3 || mat_(2,3) !=  5 || mat_(2,4) != 12 ||
             mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) != -6 || mat_(3,4) != 14 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Addition assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << mat_ << "\n"
                << "   Expected result:\n( 0  0 -2  0  8 )\n"
                                        "( 0  1  0  4 10 )\n"
                                        "( 0  0 -3  5 12 )\n"
                                        "( 0  0  0 -6 14 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing subtraction assignment via Iterator
      {
         test_ = "Row-major subtraction assignment via Iterator";

         CT col4 = blaze::column( mat_, 4UL );
         int value = 2;

         for( CT::Iterator it=begin( col4 ); it!=end( col4 ); ++it ) {
            *it -= value++;
         }

         if( col4[0] != 6 || col4[1] != 7 || col4[2] != 8 || col4[3] != 9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subtraction assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << col4 << "\n"
                << "   Expected result:\n( 6 7 8 9 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) != -2 || mat_(0,3) !=  0 || mat_(0,4) != 6 ||
             mat_(1,0) != 0 || mat_(1,1) != 1 || mat_(1,2) !=  0 || mat_(1,3) !=  4 || mat_(1,4) != 7 ||
             mat_(2,0) != 0 || mat_(2,1) != 0 || mat_(2,2) != -3 || mat_(2,3) !=  5 || mat_(2,4) != 8 ||
             mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) != -6 || mat_(3,4) != 9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << mat_ << "\n"
                << "   Expected result:\n( 0  0 -2  0  6 )\n"
                                        "( 0  1  0  4  7 )\n"
                                        "( 0  0 -3  5  8 )\n"
                                        "( 0  0  0 -6  9 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing multiplication assignment via Iterator
      {
         test_ = "Row-major multiplication assignment via Iterator";

         CT col4 = blaze::column( mat_, 4UL );
         int value = 1;

         for( CT::Iterator it=begin( col4 ); it!=end( col4 ); ++it ) {
            *it *= value++;
         }

         if( col4[0] != 6 || col4[1] != 14 || col4[2] != 24 || col4[3] != 36 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Multiplication assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << col4 << "\n"
                << "   Expected result:\n( 6 14 24 36 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) != -2 || mat_(0,3) !=  0 || mat_(0,4) !=  6 ||
             mat_(1,0) != 0 || mat_(1,1) != 1 || mat_(1,2) !=  0 || mat_(1,3) !=  4 || mat_(1,4) != 14 ||
             mat_(2,0) != 0 || mat_(2,1) != 0 || mat_(2,2) != -3 || mat_(2,3) !=  5 || mat_(2,4) != 24 ||
             mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) != -6 || mat_(3,4) != 36 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << mat_ << "\n"
                << "   Expected result:\n( 0  0 -2  0  6 )\n"
                                        "( 0  1  0  4 14 )\n"
                                        "( 0  0 -3  5 24 )\n"
                                        "( 0  0  0 -6 36 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing division assignment via Iterator
      {
         test_ = "Row-major division assignment via Iterator";

         CT col4 = blaze::column( mat_, 4UL );

         for( CT::Iterator it=begin( col4 ); it!=end( col4 ); ++it ) {
            *it /= 2;
         }

         if( col4[0] != 3 || col4[1] != 7 || col4[2] != 12 || col4[3] != 18 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Division assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << col4 << "\n"
                << "   Expected result:\n( 3 7 12 18 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) != -2 || mat_(0,3) !=  0 || mat_(0,4) !=  3 ||
             mat_(1,0) != 0 || mat_(1,1) != 1 || mat_(1,2) !=  0 || mat_(1,3) !=  4 || mat_(1,4) !=  7 ||
             mat_(2,0) != 0 || mat_(2,1) != 0 || mat_(2,2) != -3 || mat_(2,3) !=  5 || mat_(2,4) != 12 ||
             mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) != -6 || mat_(3,4) != 18 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << mat_ << "\n"
                << "   Expected result:\n( 0  0 -2  0  3 )\n"
                                        "( 0  1  0  4  7 )\n"
                                        "( 0  0 -3  5 12 )\n"
                                        "( 0  0  0 -6 18 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      initialize();

      // Testing the Iterator default constructor
      {
         test_ = "Column-major Iterator default constructor";

         OCT::Iterator it{};

         if( it != OCT::Iterator() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Failed iterator default constructor\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing the ConstIterator default constructor
      {
         test_ = "Column-major ConstIterator default constructor";

         OCT::ConstIterator it{};

         if( it != OCT::ConstIterator() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Failed iterator default constructor\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing conversion from Iterator to ConstIterator
      {
         test_ = "Column-major Iterator/ConstIterator conversion";

         OCT col2 = blaze::column( tmat_, 2UL );
         OCT::ConstIterator it( begin( col2 ) );

         if( it == end( col2 ) || it->value() != -2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Failed iterator conversion detected\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Counting the number of elements in 1st column via Iterator (end-begin)
      {
         test_ = "Column-major Iterator subtraction (end-begin)";

         OCT col1 = blaze::column( tmat_, 1UL );
         const ptrdiff_t number( end( col1 ) - begin( col1 ) );

         if( number != 1L ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid number of elements detected\n"
                << " Details:\n"
                << "   Number of elements         : " << number << "\n"
                << "   Expected number of elements: 1\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Counting the number of elements in 2nd column via ConstIterator (end-begin)
      {
         test_ = "Column-major ConstIterator subtraction (end-begin)";

         OCT col2 = blaze::column( tmat_, 2UL );
         const ptrdiff_t number( cend( col2 ) - cbegin( col2 ) );

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
         test_ = "Column-major read-only access via ConstIterator";

         OCT col2 = blaze::column( tmat_, 2UL );
         OCT::ConstIterator it ( cbegin( col2 ) );
         OCT::ConstIterator end( cend( col2 ) );

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
         test_ = "Column-major assignment via Iterator";

         OCT col4 = blaze::column( tmat_, 4UL );
         int value = 6;

         for( OCT::Iterator it=begin( col4 ); it!=end( col4 ); ++it ) {
            *it = value++;
         }

         if( col4[0] != 6 || col4[1] != 7 || col4[2] != 8 || col4[3] != 9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << col4 << "\n"
                << "   Expected result:\n( 6 7 8 9 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -2 || tmat_(0,3) !=  0 || tmat_(0,4) != 6 ||
             tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) !=  0 || tmat_(1,3) !=  4 || tmat_(1,4) != 7 ||
             tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != -3 || tmat_(2,3) !=  5 || tmat_(2,4) != 8 ||
             tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << tmat_ << "\n"
                << "   Expected result:\n( 0  0 -2  0  6 )\n"
                                        "( 0  1  0  4  7 )\n"
                                        "( 0  0 -3  5  8 )\n"
                                        "( 0  0  0 -6  9 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing addition assignment via Iterator
      {
         test_ = "Column-major addition assignment via Iterator";

         OCT col4 = blaze::column( tmat_, 4UL );
         int value = 2;

         for( OCT::Iterator it=begin( col4 ); it!=end( col4 ); ++it ) {
            *it += value++;
         }

         if( col4[0] != 8 || col4[1] != 10 || col4[2] != 12 || col4[3] != 14 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Addition assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << col4 << "\n"
                << "   Expected result:\n( 8 10 12 14 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -2 || tmat_(0,3) !=  0 || tmat_(0,4) !=  8 ||
             tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) !=  0 || tmat_(1,3) !=  4 || tmat_(1,4) != 10 ||
             tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != -3 || tmat_(2,3) !=  5 || tmat_(2,4) != 12 ||
             tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 14 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << tmat_ << "\n"
                << "   Expected result:\n( 0  0 -2  0  8 )\n"
                                        "( 0  1  0  4 10 )\n"
                                        "( 0  0 -3  5 12 )\n"
                                        "( 0  0  0 -6 14 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing subtraction assignment via Iterator
      {
         test_ = "Column-major subtraction assignment via Iterator";

         OCT col4 = blaze::column( tmat_, 4UL );
         int value = 2;

         for( OCT::Iterator it=begin( col4 ); it!=end( col4 ); ++it ) {
            *it -= value++;
         }

         if( col4[0] != 6 || col4[1] != 7 || col4[2] != 8 || col4[3] != 9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subtraction assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << col4 << "\n"
                << "   Expected result:\n( 6 7 8 9 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -2 || tmat_(0,3) !=  0 || tmat_(0,4) != 6 ||
             tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) !=  0 || tmat_(1,3) !=  4 || tmat_(1,4) != 7 ||
             tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != -3 || tmat_(2,3) !=  5 || tmat_(2,4) != 8 ||
             tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << tmat_ << "\n"
                << "   Expected result:\n( 0  0 -2  0  6 )\n"
                                        "( 0  1  0  4  7 )\n"
                                        "( 0  0 -3  5  8 )\n"
                                        "( 0  0  0 -6  9 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing multiplication assignment via Iterator
      {
         test_ = "Column-major multiplication assignment via Iterator";

         OCT col4 = blaze::column( tmat_, 4UL );
         int value = 1;

         for( OCT::Iterator it=begin( col4 ); it!=end( col4 ); ++it ) {
            *it *= value++;
         }

         if( col4[0] != 6 || col4[1] != 14 || col4[2] != 24 || col4[3] != 36 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Multiplication assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << col4 << "\n"
                << "   Expected result:\n( 6 14 24 36 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -2 || tmat_(0,3) !=  0 || tmat_(0,4) !=  6 ||
             tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) !=  0 || tmat_(1,3) !=  4 || tmat_(1,4) != 14 ||
             tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != -3 || tmat_(2,3) !=  5 || tmat_(2,4) != 24 ||
             tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 36 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << tmat_ << "\n"
                << "   Expected result:\n( 0  0 -2  0  6 )\n"
                                        "( 0  1  0  4 14 )\n"
                                        "( 0  0 -3  5 24 )\n"
                                        "( 0  0  0 -6 36 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing division assignment via Iterator
      {
         test_ = "Column-major division assignment via Iterator";

         OCT col4 = blaze::column( tmat_, 4UL );

         for( OCT::Iterator it=begin( col4 ); it!=end( col4 ); ++it ) {
            *it /= 2;
         }

         if( col4[0] != 3 || col4[1] != 7 || col4[2] != 12 || col4[3] != 18 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Division assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << col4 << "\n"
                << "   Expected result:\n( 3 7 12 18 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -2 || tmat_(0,3) !=  0 || tmat_(0,4) !=  3 ||
             tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) !=  0 || tmat_(1,3) !=  4 || tmat_(1,4) !=  7 ||
             tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != -3 || tmat_(2,3) !=  5 || tmat_(2,4) != 12 ||
             tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) != 18 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << tmat_ << "\n"
                << "   Expected result:\n( 0  0 -2  0  3 )\n"
                                        "( 0  1  0  4  7 )\n"
                                        "( 0  0 -3  5 12 )\n"
                                        "( 0  0  0 -6 18 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c nonZeros() member function of the Column specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c nonZeros() member function of the Column
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseGeneralTest::testNonZeros()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major Column::nonZeros()";

      initialize();

      // Initialization check
      CT col3 = blaze::column( mat_, 3UL );

      checkSize    ( col3, 4UL );
      checkNonZeros( col3, 3UL );

      if( col3[0] != 0 || col3[1] != 4 || col3[2] != 5 || col3[3] != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << col3 << "\n"
             << "   Expected result:\n( 0 4 5 -6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Changing the number of non-zeros via the sparse column
      col3[2] = 0;

      checkSize    ( col3, 4UL );
      checkNonZeros( col3, 2UL );

      if( col3[0] != 0 || col3[1] != 4 || col3[2] != 0 || col3[3] != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << col3 << "\n"
             << "   Expected result:\n( 0 4 0 -6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Changing the number of non-zeros via the sparse matrix
      mat_(0,3) = 5;

      checkSize    ( col3, 4UL );
      checkNonZeros( col3, 3UL );

      if( col3[0] != 5 || col3[1] != 4 || col3[2] != 0 || col3[3] != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << col3 << "\n"
             << "   Expected result:\n( 5 4 0 -6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major Column::nonZeros()";

      initialize();

      // Initialization check
      OCT col3 = blaze::column( tmat_, 3UL );

      checkSize    ( col3, 4UL );
      checkNonZeros( col3, 3UL );

      if( col3[0] != 0 || col3[1] != 4 || col3[2] != 5 || col3[3] != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << col3 << "\n"
             << "   Expected result:\n( 0 4 5 -6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Changing the number of non-zeros via the sparse column
      col3[2] = 0;

      checkSize    ( col3, 4UL );
      checkNonZeros( col3, 2UL );

      if( col3[0] != 0 || col3[1] != 4 || col3[2] != 0 || col3[3] != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << col3 << "\n"
             << "   Expected result:\n( 0 4 0 -6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Changing the number of non-zeros via the sparse matrix
      tmat_(0,3) = 5;

      checkSize    ( col3, 4UL );
      checkNonZeros( col3, 3UL );

      if( col3[0] != 5 || col3[1] != 4 || col3[2] != 0 || col3[3] != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << col3 << "\n"
             << "   Expected result:\n( 5 4 0 -6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c reset() member function of the Column specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c reset() member function of the Column specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseGeneralTest::testReset()
{
   using blaze::reset;


   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major Column::reset()";

      // Resetting a single element in column 3
      {
         initialize();

         CT col3 = blaze::column( mat_, 3UL );
         reset( col3[1] );

         checkSize    ( col3, 4UL );
         checkCapacity( col3, 4UL );
         checkNonZeros( col3, 2UL );
         checkRows    ( mat_, 4UL );
         checkColumns ( mat_, 5UL );
         checkNonZeros( mat_, 9UL );

         if( col3[0] != 0 || col3[1] != 0 || col3[2] != 5 || col3[3] != -6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Reset operation failed\n"
                << " Details:\n"
                << "   Result:\n" << col3 << "\n"
                << "   Expected result:\n( 0 0 5 -6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Resetting the 3rd column (lvalue)
      {
         initialize();

         CT col3 = blaze::column( mat_, 3UL );
         reset( col3 );

         checkSize    ( col3, 4UL );
         checkNonZeros( col3, 0UL );
         checkRows    ( mat_, 4UL );
         checkColumns ( mat_, 5UL );
         checkNonZeros( mat_, 7UL );

         if( col3[0] != 0 || col3[1] != 0 || col3[2] != 0 || col3[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Reset operation of 3rd column failed\n"
                << " Details:\n"
                << "   Result:\n" << col3 << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Resetting the 4th column (rvalue)
      {
         initialize();

         reset( blaze::column( mat_, 4UL ) );

         checkRows    ( mat_, 4UL );
         checkColumns ( mat_, 5UL );
         checkNonZeros( mat_, 6UL );

         if( mat_(0,4) != 0 || mat_(1,4) != 0 || mat_(2,4) != 0 || mat_(3,4) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Reset operation of 4th column failed\n"
                << " Details:\n"
                << "   Result:\n" << mat_ << "\n"
                << "   Expected result:\n( 0  0 -2  0  0 )\n"
                                        "( 0  1  0  4  0 )\n"
                                        "( 0  0 -3  5  0 )\n"
                                        "( 0  0  0 -6  0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major Column::reset()";

      // Resetting a single element in column 3
      {
         initialize();

         OCT col3 = blaze::column( tmat_, 3UL );
         reset( col3[1] );

         checkSize    ( col3 , 4UL );
         checkCapacity( col3 , 4UL );
         checkNonZeros( col3 , 2UL );
         checkRows    ( tmat_, 4UL );
         checkColumns ( tmat_, 5UL );
         checkNonZeros( tmat_, 9UL );

         if( col3[0] != 0 || col3[1] != 0 || col3[2] != 5 || col3[3] != -6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Reset operation failed\n"
                << " Details:\n"
                << "   Result:\n" << col3 << "\n"
                << "   Expected result:\n( 0 0 5 -6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Resetting the 3rd column (lvalue)
      {
         initialize();

         OCT col3 = blaze::column( tmat_, 3UL );
         reset( col3 );

         checkSize    ( col3 , 4UL );
         checkNonZeros( col3 , 0UL );
         checkRows    ( tmat_, 4UL );
         checkColumns ( tmat_, 5UL );
         checkNonZeros( tmat_, 7UL );

         if( col3[0] != 0 || col3[1] != 0 || col3[2] != 0 || col3[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Reset operation of 3rd column failed\n"
                << " Details:\n"
                << "   Result:\n" << col3 << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Resetting the 4th column (rvalue)
      {
         initialize();

         reset( blaze::column( tmat_, 4UL ) );

         checkRows    ( tmat_, 4UL );
         checkColumns ( tmat_, 5UL );
         checkNonZeros( tmat_, 6UL );

         if( tmat_(0,4) != 0 || tmat_(1,4) != 0 || tmat_(2,4) != 0 || tmat_(3,4) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Reset operation of 4th column failed\n"
                << " Details:\n"
                << "   Result:\n" << tmat_ << "\n"
                << "   Expected result:\n( 0  0 -2  0  0 )\n"
                                        "( 0  1  0  4  0 )\n"
                                        "( 0  0 -3  5  0 )\n"
                                        "( 0  0  0 -6  0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c clear() function with the Column specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c clear() function with the Column specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseGeneralTest::testClear()
{
   using blaze::clear;


   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major clear() function";

      // Clearing a single element in column 3
      {
         initialize();

         CT col3 = blaze::column( mat_, 3UL );
         clear( col3[1] );

         checkSize    ( col3, 4UL );
         checkCapacity( col3, 4UL );
         checkNonZeros( col3, 2UL );
         checkRows    ( mat_, 4UL );
         checkColumns ( mat_, 5UL );
         checkNonZeros( mat_, 9UL );

         if( col3[0] != 0 || col3[1] != 0 || col3[2] != 5 || col3[3] != -6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Clear operation failed\n"
                << " Details:\n"
                << "   Result:\n" << col3 << "\n"
                << "   Expected result:\n( 0 0 5 -6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Clearing the 3rd column (lvalue)
      {
         initialize();

         CT col3 = blaze::column( mat_, 3UL );
         clear( col3 );

         checkSize    ( col3, 4UL );
         checkNonZeros( col3, 0UL );
         checkRows    ( mat_, 4UL );
         checkColumns ( mat_, 5UL );
         checkNonZeros( mat_, 7UL );

         if( col3[0] != 0 || col3[1] != 0 || col3[2] != 0 || col3[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Clear operation of 3rd column failed\n"
                << " Details:\n"
                << "   Result:\n" << col3 << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Clearing the 4th column (rvalue)
      {
         initialize();

         clear( blaze::column( mat_, 4UL ) );

         checkRows    ( mat_, 4UL );
         checkColumns ( mat_, 5UL );
         checkNonZeros( mat_, 6UL );

         if( mat_(0,4) != 0 || mat_(1,4) != 0 || mat_(2,4) != 0 || mat_(3,4) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Clear operation of 4th column failed\n"
                << " Details:\n"
                << "   Result:\n" << mat_ << "\n"
                << "   Expected result:\n( 0  0 -2  0  0 )\n"
                                        "( 0  1  0  4  0 )\n"
                                        "( 0  0 -3  5  0 )\n"
                                        "( 0  0  0 -6  0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major clear() function";

      // Clearing a single element in column 3
      {
         initialize();

         OCT col3 = blaze::column( tmat_, 3UL );
         clear( col3[1] );

         checkSize    ( col3 , 4UL );
         checkCapacity( col3 , 4UL );
         checkNonZeros( col3 , 2UL );
         checkRows    ( tmat_, 4UL );
         checkColumns ( tmat_, 5UL );
         checkNonZeros( tmat_, 9UL );

         if( col3[0] != 0 || col3[1] != 0 || col3[2] != 5 || col3[3] != -6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Clear operation failed\n"
                << " Details:\n"
                << "   Result:\n" << col3 << "\n"
                << "   Expected result:\n( 0 0 5 -6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Clearing the 3rd column (lvalue)
      {
         initialize();

         OCT col3 = blaze::column( tmat_, 3UL );
         clear( col3 );

         checkSize    ( col3 , 4UL );
         checkNonZeros( col3 , 0UL );
         checkRows    ( tmat_, 4UL );
         checkColumns ( tmat_, 5UL );
         checkNonZeros( tmat_, 7UL );

         if( col3[0] != 0 || col3[1] != 0 || col3[2] != 0 || col3[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Clear operation of 3rd column failed\n"
                << " Details:\n"
                << "   Result:\n" << col3 << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Clearing the 4th column (rvalue)
      {
         initialize();

         clear( blaze::column( tmat_, 4UL ) );

         checkRows    ( tmat_, 4UL );
         checkColumns ( tmat_, 5UL );
         checkNonZeros( tmat_, 6UL );

         if( tmat_(0,4) != 0 || tmat_(1,4) != 0 || tmat_(2,4) != 0 || tmat_(3,4) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Clear operation of 4th column failed\n"
                << " Details:\n"
                << "   Result:\n" << tmat_ << "\n"
                << "   Expected result:\n( 0  0 -2  0  0 )\n"
                                        "( 0  1  0  4  0 )\n"
                                        "( 0  0 -3  5  0 )\n"
                                        "( 0  0  0 -6  0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c reserve() member function of the Column specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c reserve() member function of the Column specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseGeneralTest::testReserve()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major Column::reserve()";

      MT mat( 20UL, 3UL );

      CT col0 = blaze::column( mat, 0UL );

      // Increasing the capacity of the column
      col0.reserve( 10UL );

      checkSize    ( col0, 20UL );
      checkCapacity( col0, 10UL );
      checkNonZeros( col0,  0UL );

      // Further increasing the capacity of the column
      col0.reserve( 15UL );

      checkSize    ( col0, 20UL );
      checkCapacity( col0, 15UL );
      checkNonZeros( col0,  0UL );
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major Column::reserve()";

      OMT mat( 20UL, 3UL );

      OCT col0 = blaze::column( mat, 0UL );

      // Increasing the capacity of the column
      col0.reserve( 10UL );

      checkSize    ( col0, 20UL );
      checkCapacity( col0, 10UL );
      checkNonZeros( col0,  0UL );

      // Further increasing the capacity of the column
      col0.reserve( 15UL );

      checkSize    ( col0, 20UL );
      checkCapacity( col0, 15UL );
      checkNonZeros( col0,  0UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c set() member function of the Column specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c set() member function of the Column specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseGeneralTest::testSet()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major Column::set()";

      initialize();

      CT col0 = blaze::column( mat_, 0UL );

      // Setting a non-zero element at the end of the column
      {
         CT::Iterator pos = col0.set( 3UL, 1 );

         checkSize    ( col0,  4UL );
         checkNonZeros( col0,  1UL );
         checkRows    ( mat_,  4UL );
         checkColumns ( mat_,  5UL );
         checkNonZeros( mat_, 11UL );

         if( pos->value() != 1 || pos->index() != 3UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 1\n"
                << "   Expected index: 3\n";
            throw std::runtime_error( oss.str() );
         }

         if( col0[0] != 0 || col0[1] != 0 || col0[2] != 0 || col0[3] != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << col0 << "\n"
                << "   Expected result:\n( 0 0 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setting a non-zero element at the beginning of the column
      {
         CT::Iterator pos = col0.set( 0UL, 2 );

         checkSize    ( col0,  4UL );
         checkNonZeros( col0,  2UL );
         checkRows    ( mat_,  4UL );
         checkColumns ( mat_,  5UL );
         checkNonZeros( mat_, 12UL );

         if( pos->value() != 2 || pos->index() != 0UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 2\n"
                << "   Expected index: 0\n";
            throw std::runtime_error( oss.str() );
         }

         if( col0[0] != 2 || col0[1] != 0 || col0[2] != 0 || col0[3] != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << col0 << "\n"
                << "   Expected result:\n( 2 0 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setting a non-zero element at the center of the column
      {
         CT::Iterator pos = col0.set( 2UL, 3 );

         checkSize    ( col0,  4UL );
         checkNonZeros( col0,  3UL );
         checkRows    ( mat_,  4UL );
         checkColumns ( mat_,  5UL );
         checkNonZeros( mat_, 13UL );

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

         if( col0[0] != 2 || col0[1] != 0 || col0[2] != 3 || col0[3] != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << col0 << "\n"
                << "   Expected result:\n( 2 0 3 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setting an already existing element
      {
         CT::Iterator pos = col0.set( 3UL, 4 );

         checkSize    ( col0,  4UL );
         checkNonZeros( col0,  3UL );
         checkRows    ( mat_,  4UL );
         checkColumns ( mat_,  5UL );
         checkNonZeros( mat_, 13UL );

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

         if( col0[0] != 2 || col0[1] != 0 || col0[2] != 3 || col0[3] != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << col0 << "\n"
                << "   Expected result:\n( 2 0 3 4 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major Column::set()";

      initialize();

      OCT col0 = blaze::column( tmat_, 0UL );

      // Setting a non-zero element at the end of the column
      {
         OCT::Iterator pos = col0.set( 3UL, 1 );

         checkSize    ( col0 ,  4UL );
         checkNonZeros( col0 ,  1UL );
         checkRows    ( tmat_,  4UL );
         checkColumns ( tmat_,  5UL );
         checkNonZeros( tmat_, 11UL );

         if( pos->value() != 1 || pos->index() != 3UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 1\n"
                << "   Expected index: 3\n";
            throw std::runtime_error( oss.str() );
         }

         if( col0[0] != 0 || col0[1] != 0 || col0[2] != 0 || col0[3] != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << col0 << "\n"
                << "   Expected result:\n( 0 0 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setting a non-zero element at the beginning of the column
      {
         OCT::Iterator pos = col0.set( 0UL, 2 );

         checkSize    ( col0 ,  4UL );
         checkNonZeros( col0 ,  2UL );
         checkRows    ( tmat_,  4UL );
         checkColumns ( tmat_,  5UL );
         checkNonZeros( tmat_, 12UL );

         if( pos->value() != 2 || pos->index() != 0UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 2\n"
                << "   Expected index: 0\n";
            throw std::runtime_error( oss.str() );
         }

         if( col0[0] != 2 || col0[1] != 0 || col0[2] != 0 || col0[3] != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << col0 << "\n"
                << "   Expected result:\n( 2 0 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setting a non-zero element at the center of the column
      {
         OCT::Iterator pos = col0.set( 2UL, 3 );

         checkSize    ( col0 ,  4UL );
         checkNonZeros( col0 ,  3UL );
         checkRows    ( tmat_,  4UL );
         checkColumns ( tmat_,  5UL );
         checkNonZeros( tmat_, 13UL );

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

         if( col0[0] != 2 || col0[1] != 0 || col0[2] != 3 || col0[3] != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << col0 << "\n"
                << "   Expected result:\n( 2 0 3 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setting an already existing element
      {
         OCT::Iterator pos = col0.set( 3UL, 4 );

         checkSize    ( col0 ,  4UL );
         checkNonZeros( col0 ,  3UL );
         checkRows    ( tmat_,  4UL );
         checkColumns ( tmat_,  5UL );
         checkNonZeros( tmat_, 13UL );

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

         if( col0[0] != 2 || col0[1] != 0 || col0[2] != 3 || col0[3] != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << col0 << "\n"
                << "   Expected result:\n( 2 0 3 4 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c insert() member function of the Column specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c insert() member function of the Column specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseGeneralTest::testInsert()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major Column::insert()";

      initialize();

      CT col0 = blaze::column( mat_, 0UL );

      // Inserting a non-zero element at the end of the column
      {
         CT::Iterator pos = col0.insert( 3UL, 1 );

         checkSize    ( col0,  4UL );
         checkNonZeros( col0,  1UL );
         checkRows    ( mat_,  4UL );
         checkColumns ( mat_,  5UL );
         checkNonZeros( mat_, 11UL );

         if( pos->value() != 1 || pos->index() != 3UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 1\n"
                << "   Expected index: 3\n";
            throw std::runtime_error( oss.str() );
         }

         if( col0[0] != 0 || col0[1] != 0 || col0[2] != 0 || col0[3] != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << col0 << "\n"
                << "   Expected result:\n( 0 0 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting a non-zero element at the beginning of the column
      {
         CT::Iterator pos = col0.insert( 0UL, 2 );

         checkSize    ( col0,  4UL );
         checkNonZeros( col0,  2UL );
         checkRows    ( mat_,  4UL );
         checkColumns ( mat_,  5UL );
         checkNonZeros( mat_, 12UL );

         if( pos->value() != 2 || pos->index() != 0UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 2\n"
                << "   Expected index: 0\n";
            throw std::runtime_error( oss.str() );
         }

         if( col0[0] != 2 || col0[1] != 0 || col0[2] != 0 || col0[3] != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << col0 << "\n"
                << "   Expected result:\n( 2 0 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting a non-zero element at the center of the column
      {
         CT::Iterator pos = col0.insert( 2UL, 3 );

         checkSize    ( col0,  4UL );
         checkNonZeros( col0,  3UL );
         checkRows    ( mat_,  4UL );
         checkColumns ( mat_,  5UL );
         checkNonZeros( mat_, 13UL );

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

         if( col0[0] != 2 || col0[1] != 0 || col0[2] != 3 || col0[3] != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << col0 << "\n"
                << "   Expected result:\n( 2 0 3 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to insert an already existing element
      try {
         col0.insert( 3UL, 4 );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Inserting an existing element succeeded\n"
             << " Details:\n"
             << "   Result:\n" << col0 << "\n"
             << "   Expected result:\n( 2 0 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major Column::insert()";

      initialize();

      OCT col0 = blaze::column( tmat_, 0UL );

      // Inserting a non-zero element at the end of the column
      {
         OCT::Iterator pos = col0.insert( 3UL, 1 );

         checkSize    ( col0 ,  4UL );
         checkNonZeros( col0 ,  1UL );
         checkRows    ( tmat_,  4UL );
         checkColumns ( tmat_,  5UL );
         checkNonZeros( tmat_, 11UL );

         if( pos->value() != 1 || pos->index() != 3UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 1\n"
                << "   Expected index: 3\n";
            throw std::runtime_error( oss.str() );
         }

         if( col0[0] != 0 || col0[1] != 0 || col0[2] != 0 || col0[3] != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << col0 << "\n"
                << "   Expected result:\n( 0 0 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting a non-zero element at the beginning of the column
      {
         OCT::Iterator pos = col0.insert( 0UL, 2 );

         checkSize    ( col0 ,  4UL );
         checkNonZeros( col0 ,  2UL );
         checkRows    ( tmat_,  4UL );
         checkColumns ( tmat_,  5UL );
         checkNonZeros( tmat_, 12UL );

         if( pos->value() != 2 || pos->index() != 0UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 2\n"
                << "   Expected index: 0\n";
            throw std::runtime_error( oss.str() );
         }

         if( col0[0] != 2 || col0[1] != 0 || col0[2] != 0 || col0[3] != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << col0 << "\n"
                << "   Expected result:\n( 2 0 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting a non-zero element at the center of the column
      {
         OCT::Iterator pos = col0.insert( 2UL, 3 );

         checkSize    ( col0 ,  4UL );
         checkNonZeros( col0 ,  3UL );
         checkRows    ( tmat_,  4UL );
         checkColumns ( tmat_,  5UL );
         checkNonZeros( tmat_, 13UL );

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

         if( col0[0] != 2 || col0[1] != 0 || col0[2] != 3 || col0[3] != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << col0 << "\n"
                << "   Expected result:\n( 2 0 3 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to insert an already existing element
      try {
         col0.insert( 3UL, 4 );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Inserting an existing element succeeded\n"
             << " Details:\n"
             << "   Result:\n" << col0 << "\n"
             << "   Expected result:\n( 2 0 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c append() member function of the Column specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c append() member function of the Column specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseGeneralTest::testAppend()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major Column::append()";

      MT mat( 9UL, 3UL );

      CT col1 = blaze::column( mat, 1UL );
      col1.reserve( 4UL );

      // Appending one non-zero element
      col1.append( 1UL, 1 );

      checkSize    ( col1, 9UL );
      checkCapacity( col1, 4UL );
      checkNonZeros( col1, 1UL );

      if( col1[1] != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Append operation failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( 0 1 0 0 0 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Appending three more non-zero elements
      col1.append( 3UL, 2 );
      col1.append( 4UL, 3 );
      col1.append( 8UL, 4 );

      checkSize    ( col1, 9UL );
      checkCapacity( col1, 4UL );
      checkNonZeros( col1, 4UL );

      if( col1[1] != 1 || col1[3] != 2 || col1[4] != 3 || col1[8] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Append operation failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( 0 1 0 2 3 0 0 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major Column::append()";

      OMT mat( 9UL, 3UL );

      OCT col1 = blaze::column( mat, 1UL );
      col1.reserve( 4UL );

      // Appending one non-zero element
      col1.append( 1UL, 1 );

      checkSize    ( col1, 9UL );
      checkCapacity( col1, 4UL );
      checkNonZeros( col1, 1UL );

      if( col1[1] != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Append operation failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( 0 1 0 0 0 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Appending three more non-zero elements
      col1.append( 3UL, 2 );
      col1.append( 4UL, 3 );
      col1.append( 8UL, 4 );

      checkSize    ( col1, 9UL );
      checkCapacity( col1, 4UL );
      checkNonZeros( col1, 4UL );

      if( col1[1] != 1 || col1[3] != 2 || col1[4] != 3 || col1[8] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Append operation failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( 0 1 0 2 3 0 0 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c erase() member function of the Column specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c erase() member function of the Column specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseGeneralTest::testErase()
{
   //=====================================================================================
   // Row-major index-based erase function
   //=====================================================================================

   {
      test_ = "Row-major Column::erase( size_t )";

      initialize();

      CT col4 = blaze::column( mat_, 4UL );

      // Erasing the non-zero element at the end of the column
      col4.erase( 3UL );

      checkSize    ( col4, 4UL );
      checkNonZeros( col4, 3UL );
      checkRows    ( mat_, 4UL );
      checkColumns ( mat_, 5UL );
      checkNonZeros( mat_, 9UL );

      if( col4[0] != 7 || col4[1] != -8 || col4[2] != 9 || col4[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << col4 << "\n"
             << "   Expected result:\n( 7 -8 9 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the non-zero element at the beginning of the column
      col4.erase( 0UL );

      checkSize    ( col4, 4UL );
      checkNonZeros( col4, 2UL );
      checkRows    ( mat_, 4UL );
      checkColumns ( mat_, 5UL );
      checkNonZeros( mat_, 8UL );

      if( col4[0] != 0 || col4[1] != -8 || col4[2] != 9 || col4[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << col4 << "\n"
             << "   Expected result:\n( 0 -8 9 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the non-zero element at the beginning of the column
      col4.erase( 2UL );

      checkSize    ( col4, 4UL );
      checkNonZeros( col4, 1UL );
      checkRows    ( mat_, 4UL );
      checkColumns ( mat_, 5UL );
      checkNonZeros( mat_, 7UL );

      if( col4[0] != 0 || col4[1] != -8 || col4[2] != 0 || col4[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << col4 << "\n"
             << "   Expected result:\n( 0 -8 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to erase an already erased element
      col4.erase( 3UL );

      checkSize    ( col4, 4UL );
      checkNonZeros( col4, 1UL );
      checkRows    ( mat_, 4UL );
      checkColumns ( mat_, 5UL );
      checkNonZeros( mat_, 7UL );

      if( col4[0] != 0 || col4[1] != -8 || col4[2] != 0 || col4[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << col4 << "\n"
             << "   Expected result:\n( 0 -8 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major iterator-based erase function
   //=====================================================================================

   {
      test_ = "Row-major Column::erase( Iterator )";

      initialize();

      CT col4 = blaze::column( mat_, 4UL );

      // Erasing the non-zero element at the end of the column
      {
         CT::Iterator pos = col4.erase( col4.find( 3UL ) );

         checkSize    ( col4, 4UL );
         checkNonZeros( col4, 3UL );
         checkRows    ( mat_, 4UL );
         checkColumns ( mat_, 5UL );
         checkNonZeros( mat_, 9UL );

         if( pos != col4.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( col4[0] != 7 || col4[1] != -8 || col4[2] != 9 || col4[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << col4 << "\n"
                << "   Expected result:\n( 7 -8 9 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the non-zero element at the beginning of the column
      {
         CT::Iterator pos = col4.erase( col4.find( 0UL ) );

         checkSize    ( col4, 4UL );
         checkNonZeros( col4, 2UL );
         checkRows    ( mat_, 4UL );
         checkColumns ( mat_, 5UL );
         checkNonZeros( mat_, 8UL );

         if( pos->value() != -8 || pos->index() != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: -8\n"
                << "   Expected index:  1\n";
            throw std::runtime_error( oss.str() );
         }

         if( col4[0] != 0 || col4[1] != -8 || col4[2] != 9 || col4[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << col4 << "\n"
                << "   Expected result:\n( 0 -8 9 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the non-zero element at the beginning of the column
      {
         CT::Iterator pos = col4.erase( col4.find( 2UL ) );

         checkSize    ( col4, 4UL );
         checkNonZeros( col4, 1UL );
         checkRows    ( mat_, 4UL );
         checkColumns ( mat_, 5UL );
         checkNonZeros( mat_, 7UL );

         if( pos != col4.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( col4[0] != 0 || col4[1] != -8 || col4[2] != 0 || col4[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << col4 << "\n"
                << "   Expected result:\n( 0 -8 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to erase an already erased element
      {
         CT::Iterator pos = col4.erase( col4.find( 3UL ) );

         checkSize    ( col4, 4UL );
         checkNonZeros( col4, 1UL );
         checkRows    ( mat_, 4UL );
         checkColumns ( mat_, 5UL );
         checkNonZeros( mat_, 7UL );

         if( pos != col4.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( col4[0] != 0 || col4[1] != -8 || col4[2] != 0 || col4[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << col4 << "\n"
                << "   Expected result:\n( 0 -8 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Row-major iterator-range-based erase function
   //=====================================================================================

   {
      test_ = "Row-major Column::erase( Iterator, Iterator )";

      initialize();

      // Erasing the 2nd column
      {
         CT col2 = blaze::column( mat_, 2UL );

         CT::Iterator pos = col2.erase( col2.begin(), col2.end() );

         checkSize    ( col2, 4UL );
         checkNonZeros( col2, 0UL );
         checkRows    ( mat_, 4UL );
         checkColumns ( mat_, 5UL );
         checkNonZeros( mat_, 8UL );

         if( pos != col2.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( col2[0] != 0 || col2[1] != 0 || col2[2] != 0 || col2[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing the column failed\n"
                << " Details:\n"
                << "   Result:\n" << col2 << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the first half of the 4th column
      {
         CT col4 = blaze::column( mat_, 4UL );

         CT::Iterator pos = col4.erase( col4.begin(), col4.find( 2UL ) );

         checkSize    ( col4, 4UL );
         checkNonZeros( col4, 2UL );
         checkRows    ( mat_, 4UL );
         checkColumns ( mat_, 5UL );
         checkNonZeros( mat_, 6UL );

         if( pos->value() != 9 || pos->index() != 2 ) {
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

         if( col4[0] != 0 || col4[1] != 0 || col4[2] != 9 || col4[3] != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a partial column failed\n"
                << " Details:\n"
                << "   Result:\n" << col4 << "\n"
                << "   Expected result:\n( 0 0 9 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the second half of the 4th column
      {
         CT col4 = blaze::column( mat_, 4UL );

         CT::Iterator pos = col4.erase( col4.find( 2UL ), col4.end() );

         checkSize    ( col4, 4UL );
         checkNonZeros( col4, 0UL );
         checkRows    ( mat_, 4UL );
         checkColumns ( mat_, 5UL );
         checkNonZeros( mat_, 4UL );

         if( pos != col4.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( col4[0] != 0 || col4[1] != 0 || col4[2] != 0 || col4[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a partial column failed\n"
                << " Details:\n"
                << "   Result:\n" << col4 << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to erase an empty range
      {
         CT col3 = blaze::column( mat_, 3UL );

         CT::Iterator pos = col3.erase( col3.find( 1UL ), col3.find( 1UL ) );

         checkSize    ( col3, 4UL );
         checkNonZeros( col3, 3UL );
         checkRows    ( mat_, 4UL );
         checkColumns ( mat_, 5UL );
         checkNonZeros( mat_, 4UL );

         if( pos != col3.find( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the given end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( col3[0] != 0 || col3[1] != 4 || col3[2] != 5 || col3[3] != -6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing an empty range failed\n"
                << " Details:\n"
                << "   Result:\n" << col3 << "\n"
                << "   Expected result:\n( 0 4 5 -6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Row-major erase function with predicate
   //=====================================================================================

   {
      test_ = "Row-major Column::erase( Predicate )";

      initialize();

      CT col4 = blaze::column( mat_, 4UL );

      // Erasing a selection of elements
      col4.erase( []( int value ) { return value == 7 || value == 10; } );

      checkSize    ( col4, 4UL );
      checkNonZeros( col4, 2UL );
      checkRows    ( mat_, 4UL );
      checkColumns ( mat_, 5UL );
      checkNonZeros( mat_, 8UL );

      if( col4[0] != 0 || col4[1] != -8 || col4[2] != 9 || col4[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a selection of elements failed\n"
             << " Details:\n"
             << "   Result:\n" << col4 << "\n"
             << "   Expected result:\n( 0 -8 9 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to erase all elements with value 1
      col4.erase( []( int value ){ return value == 1; } );

      checkSize    ( col4, 4UL );
      checkNonZeros( col4, 2UL );
      checkRows    ( mat_, 4UL );
      checkColumns ( mat_, 5UL );
      checkNonZeros( mat_, 8UL );

      if( col4[0] != 0 || col4[1] != -8 || col4[2] != 9 || col4[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing all elements with value 1 failed\n"
             << " Details:\n"
             << "   Result:\n" << col4 << "\n"
             << "   Expected result:\n( 0 -8 9 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major iterator-range-based erase function with predicate
   //=====================================================================================

   {
      test_ = "Row-major Column::erase( Iterator, Iterator, Predicate )";

      initialize();

      CT col4 = blaze::column( mat_, 4UL );

      // Erasing a selection of elements
      col4.erase( col4.find( 1UL ), col4.end(),
                  []( int value ) { return value == -8 || value == 10; } );

      checkSize    ( col4, 4UL );
      checkNonZeros( col4, 2UL );
      checkRows    ( mat_, 4UL );
      checkColumns ( mat_, 5UL );
      checkNonZeros( mat_, 8UL );

      if( col4[0] != 7 || col4[1] != 0 || col4[2] != 9 || col4[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a selection of elements failed\n"
             << " Details:\n"
             << "   Result:\n" << col4 << "\n"
             << "   Expected result:\n( 7 0 9 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to erase from an empty range
      col4.erase( col4.begin(), col4.begin(), []( int ){ return true; } );

      checkSize    ( col4, 4UL );
      checkNonZeros( col4, 2UL );
      checkRows    ( mat_, 4UL );
      checkColumns ( mat_, 5UL );
      checkNonZeros( mat_, 8UL );

      if( col4[0] != 7 || col4[1] != 0 || col4[2] != 9 || col4[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing from an empty range failed\n"
             << " Details:\n"
             << "   Result:\n" << col4 << "\n"
             << "   Expected result:\n( 7 0 9 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major index-based erase function
   //=====================================================================================

   {
      test_ = "Column-major Column::erase( size_t )";

      initialize();

      OCT col4 = blaze::column( tmat_, 4UL );

      // Erasing the non-zero element at the end of the column
      col4.erase( 3UL );

      checkSize    ( col4 , 4UL );
      checkNonZeros( col4 , 3UL );
      checkRows    ( tmat_, 4UL );
      checkColumns ( tmat_, 5UL );
      checkNonZeros( tmat_, 9UL );

      if( col4[0] != 7 || col4[1] != -8 || col4[2] != 9 || col4[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << col4 << "\n"
             << "   Expected result:\n( 7 -8 9 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the non-zero element at the beginning of the column
      col4.erase( size_t(0) );

      checkSize    ( col4 , 4UL );
      checkNonZeros( col4 , 2UL );
      checkRows    ( tmat_, 4UL );
      checkColumns ( tmat_, 5UL );
      checkNonZeros( tmat_, 8UL );

      if( col4[0] != 0 || col4[1] != -8 || col4[2] != 9 || col4[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << col4 << "\n"
             << "   Expected result:\n( 0 -8 9 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the non-zero element at the beginning of the column
      col4.erase( 2UL );

      checkSize    ( col4 , 4UL );
      checkNonZeros( col4 , 1UL );
      checkRows    ( tmat_, 4UL );
      checkColumns ( tmat_, 5UL );
      checkNonZeros( tmat_, 7UL );

      if( col4[0] != 0 || col4[1] != -8 || col4[2] != 0 || col4[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << col4 << "\n"
             << "   Expected result:\n( 0 -8 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to erase an already erased element
      col4.erase( 3UL );

      checkSize    ( col4 , 4UL );
      checkNonZeros( col4 , 1UL );
      checkRows    ( tmat_, 4UL );
      checkColumns ( tmat_, 5UL );
      checkNonZeros( tmat_, 7UL );

      if( col4[0] != 0 || col4[1] != -8 || col4[2] != 0 || col4[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << col4 << "\n"
             << "   Expected result:\n( 0 -8 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major iterator-based erase function
   //=====================================================================================

   {
      test_ = "Column-major Column::erase( size_t )";

      initialize();

      OCT col4 = blaze::column( tmat_, 4UL );

      // Erasing the non-zero element at the end of the column
      {
         OCT::Iterator pos = col4.erase( col4.find( 3UL ) );

         checkSize    ( col4 , 4UL );
         checkNonZeros( col4 , 3UL );
         checkRows    ( tmat_, 4UL );
         checkColumns ( tmat_, 5UL );
         checkNonZeros( tmat_, 9UL );

         if( pos != col4.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( col4[0] != 7 || col4[1] != -8 || col4[2] != 9 || col4[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << col4 << "\n"
                << "   Expected result:\n( 7 -8 9 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the non-zero element at the beginning of the column
      {
         OCT::Iterator pos = col4.erase( col4.find( 0UL ) );

         checkSize    ( col4 , 4UL );
         checkNonZeros( col4 , 2UL );
         checkRows    ( tmat_, 4UL );
         checkColumns ( tmat_, 5UL );
         checkNonZeros( tmat_, 8UL );

         if( pos->value() != -8 || pos->index() != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: -8\n"
                << "   Expected index:  1\n";
            throw std::runtime_error( oss.str() );
         }

         if( col4[0] != 0 || col4[1] != -8 || col4[2] != 9 || col4[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << col4 << "\n"
                << "   Expected result:\n( 0 -8 9 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the non-zero element at the beginning of the column
      {
         OCT::Iterator pos = col4.erase( col4.find( 2UL ) );

         checkSize    ( col4 , 4UL );
         checkNonZeros( col4 , 1UL );
         checkRows    ( tmat_, 4UL );
         checkColumns ( tmat_, 5UL );
         checkNonZeros( tmat_, 7UL );

         if( pos != col4.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( col4[0] != 0 || col4[1] != -8 || col4[2] != 0 || col4[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << col4 << "\n"
                << "   Expected result:\n( 0 -8 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to erase an already erased element
      {
         OCT::Iterator pos = col4.erase( col4.find( 3UL ) );

         checkSize    ( col4 , 4UL );
         checkNonZeros( col4 , 1UL );
         checkRows    ( tmat_, 4UL );
         checkColumns ( tmat_, 5UL );
         checkNonZeros( tmat_, 7UL );

         if( pos != col4.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( col4[0] != 0 || col4[1] != -8 || col4[2] != 0 || col4[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << col4 << "\n"
                << "   Expected result:\n( 0 -8 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major iterator-range-based erase function
   //=====================================================================================

   {
      test_ = "Column-major Column::erase( Iterator, Iterator )";

      initialize();

      // Erasing the 2nd column
      {
         OCT col2 = blaze::column( tmat_, 2UL );

         OCT::Iterator pos = col2.erase( col2.begin(), col2.end() );

         checkSize    ( col2 , 4UL );
         checkNonZeros( col2 , 0UL );
         checkRows    ( tmat_, 4UL );
         checkColumns ( tmat_, 5UL );
         checkNonZeros( tmat_, 8UL );

         if( pos != col2.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( col2[0] != 0 || col2[1] != 0 || col2[2] != 0 || col2[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing the column failed\n"
                << " Details:\n"
                << "   Result:\n" << col2 << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the first half of the 4th column
      {
         OCT col4 = blaze::column( tmat_, 4UL );

         OCT::Iterator pos = col4.erase( col4.begin(), col4.find( 2UL ) );

         checkSize    ( col4 , 4UL );
         checkNonZeros( col4 , 2UL );
         checkRows    ( tmat_, 4UL );
         checkColumns ( tmat_, 5UL );
         checkNonZeros( tmat_, 6UL );

         if( pos->value() != 9 || pos->index() != 2 ) {
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

         if( col4[0] != 0 || col4[1] != 0 || col4[2] != 9 || col4[3] != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a partial column failed\n"
                << " Details:\n"
                << "   Result:\n" << col4 << "\n"
                << "   Expected result:\n( 0 0 9 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the second half of the 4th column
      {
         OCT col4 = blaze::column( tmat_, 4UL );

         OCT::Iterator pos = col4.erase( col4.find( 2UL ), col4.end() );

         checkSize    ( col4 , 4UL );
         checkNonZeros( col4 , 0UL );
         checkRows    ( tmat_, 4UL );
         checkColumns ( tmat_, 5UL );
         checkNonZeros( tmat_, 4UL );

         if( pos != col4.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( col4[0] != 0 || col4[1] != 0 || col4[2] != 0 || col4[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a partial column failed\n"
                << " Details:\n"
                << "   Result:\n" << col4 << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to erase an empty range
      {
         OCT col3 = blaze::column( tmat_, 3UL );

         OCT::Iterator pos = col3.erase( col3.find( 1UL ), col3.find( 1UL ) );

         checkSize    ( col3 , 4UL );
         checkNonZeros( col3 , 3UL );
         checkRows    ( tmat_, 4UL );
         checkColumns ( tmat_, 5UL );
         checkNonZeros( tmat_, 4UL );

         if( pos != col3.find( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the given end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( col3[0] != 0 || col3[1] != 4 || col3[2] != 5 || col3[3] != -6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing an empty range failed\n"
                << " Details:\n"
                << "   Result:\n" << col3 << "\n"
                << "   Expected result:\n( 0 4 5 -6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major erase function with predicate
   //=====================================================================================

   {
      test_ = "Column-major Column::erase( Predicate )";

      initialize();

      OCT col4 = blaze::column( tmat_, 4UL );

      // Erasing a selection of elements
      col4.erase( []( int value ) { return value == 7 || value == 10; } );

      checkSize    ( col4 , 4UL );
      checkNonZeros( col4 , 2UL );
      checkRows    ( tmat_, 4UL );
      checkColumns ( tmat_, 5UL );
      checkNonZeros( tmat_, 8UL );

      if( col4[0] != 0 || col4[1] != -8 || col4[2] != 9 || col4[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a selection of elements failed\n"
             << " Details:\n"
             << "   Result:\n" << col4 << "\n"
             << "   Expected result:\n( 0 -8 9 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to erase all elements with value 1
      col4.erase( []( int value ){ return value == 1; } );

      checkSize    ( col4 , 4UL );
      checkNonZeros( col4 , 2UL );
      checkRows    ( tmat_, 4UL );
      checkColumns ( tmat_, 5UL );
      checkNonZeros( tmat_, 8UL );

      if( col4[0] != 0 || col4[1] != -8 || col4[2] != 9 || col4[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing all elements with value 1 failed\n"
             << " Details:\n"
             << "   Result:\n" << col4 << "\n"
             << "   Expected result:\n( 0 -8 9 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major iterator-range-based erase function with predicate
   //=====================================================================================

   {
      test_ = "Column-major Column::erase( Iterator, Iterator, Predicate )";

      initialize();

      OCT col4 = blaze::column( tmat_, 4UL );

      // Erasing a selection of elements
      col4.erase( col4.find( 1UL ), col4.end(),
                  []( int value ) { return value == -8 || value == 10; } );

      checkSize    ( col4 , 4UL );
      checkNonZeros( col4 , 2UL );
      checkRows    ( tmat_, 4UL );
      checkColumns ( tmat_, 5UL );
      checkNonZeros( tmat_, 8UL );

      if( col4[0] != 7 || col4[1] != 0 || col4[2] != 9 || col4[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a selection of elements failed\n"
             << " Details:\n"
             << "   Result:\n" << col4 << "\n"
             << "   Expected result:\n( 0 -8 9 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to erase from an empty range
      col4.erase( col4.begin(), col4.begin(), []( int ){ return true; } );

      checkSize    ( col4 , 4UL );
      checkNonZeros( col4 , 2UL );
      checkRows    ( tmat_, 4UL );
      checkColumns ( tmat_, 5UL );
      checkNonZeros( tmat_, 8UL );

      if( col4[0] != 7 || col4[1] != 0 || col4[2] != 9 || col4[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing from an empty range failed\n"
             << " Details:\n"
             << "   Result:\n" << col4 << "\n"
             << "   Expected result:\n( 0 -8 9 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************



//*************************************************************************************************
/*!\brief Test of the \c find() member function of the Column specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c find() member function of the Column specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseGeneralTest::testFind()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major Column::find()";

      initialize();

      CT col2 = blaze::column( mat_, 2UL );

      // Searching for the first element
      {
         CT::Iterator pos = col2.find( 0UL );

         if( pos == col2.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required index = 0\n"
                << "   Current column:\n" << col2 << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 0 || pos->value() != -2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 0\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = -2\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current column:\n" << col2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for the second element
      {
         CT::Iterator pos = col2.find( 2UL );

         if( pos == col2.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required index = 2\n"
                << "   Current column:\n" << col2 << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 2 || pos->value() != -3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 2\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = -3\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current column:\n" << col2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for a non-existing non-zero element
      {
         CT::Iterator pos = col2.find( 1UL );

         if( pos != col2.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Non-existing element could be found\n"
                << " Details:\n"
                << "   Required index = 1\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 0\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current column:\n" << col2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major Column::find()";

      initialize();

      OCT col2 = blaze::column( tmat_, 2UL );

      // Searching for the first element
      {
         OCT::Iterator pos = col2.find( 0UL );

         if( pos == col2.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required index = 0\n"
                << "   Current column:\n" << col2 << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 0 || pos->value() != -2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 0\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = -2\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current column:\n" << col2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for the second element
      {
         OCT::Iterator pos = col2.find( 2UL );

         if( pos == col2.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required index = 2\n"
                << "   Current column:\n" << col2 << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 2 || pos->value() != -3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 2\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = -3\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current column:\n" << col2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for a non-existing non-zero element
      {
         OCT::Iterator pos = col2.find( 1UL );

         if( pos != col2.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Non-existing element could be found\n"
                << " Details:\n"
                << "   Required index = 1\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 0\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current column:\n" << col2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c lowerBound() member function of the Column specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c lowerBound() member function of the Column
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseGeneralTest::testLowerBound()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major Column::lowerBound()";

      initialize();

      CT col1 = blaze::column( mat_, 1UL );

      // Determining the lower bound for index 0
      {
         CT::Iterator pos = col1.lowerBound( 0UL );

         if( pos == col1.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required index = 0\n"
                << "   Current column:\n" << col1 << "\n";
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
                << "   Current column:\n" << col1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for index 1
      {
         CT::Iterator pos = col1.lowerBound( 1UL );

         if( pos == col1.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required index = 1\n"
                << "   Current column:\n" << col1 << "\n";
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
                << "   Current column:\n" << col1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for index 2
      {
         CT::Iterator pos = col1.lowerBound( 2UL );

         if( pos != col1.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required index = 2\n"
                << "   Current column:\n" << col1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major Column::lowerBound()";

      initialize();

      OCT col1 = blaze::column( tmat_, 1UL );

      // Determining the lower bound for index 0
      {
         OCT::Iterator pos = col1.lowerBound( 0UL );

         if( pos == col1.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required index = 0\n"
                << "   Current column:\n" << col1 << "\n";
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
                << "   Current column:\n" << col1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for index 1
      {
         OCT::Iterator pos = col1.lowerBound( 1UL );

         if( pos == col1.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required index = 1\n"
                << "   Current column:\n" << col1 << "\n";
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
                << "   Current column:\n" << col1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for index 2
      {
         OCT::Iterator pos = col1.lowerBound( 2UL );

         if( pos != col1.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required index = 2\n"
                << "   Current column:\n" << col1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c upperBound() member function of the Column specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c upperBound() member function of the Column
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseGeneralTest::testUpperBound()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major Column::upperBound()";

      initialize();

      CT col1 = blaze::column( mat_, 1UL );

      // Determining the upper bound for index 0
      {
         CT::Iterator pos = col1.upperBound( 0UL );

         if( pos == col1.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required index = 0\n"
                << "   Current column:\n" << col1 << "\n";
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
                << "   Current column:\n" << col1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for index 1
      {
         CT::Iterator pos = col1.upperBound( 1UL );

         if( pos != col1.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required index = 1\n"
                << "   Current column:\n" << col1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for index 2
      {
         CT::Iterator pos = col1.upperBound( 2UL );

         if( pos != col1.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required index = 2\n"
                << "   Current column:\n" << col1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major Column::upperBound()";

      initialize();

      OCT col1 = blaze::column( tmat_, 1UL );

      // Determining the upper bound for index 0
      {
         OCT::Iterator pos = col1.upperBound( 0UL );

         if( pos == col1.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required index = 0\n"
                << "   Current column:\n" << col1 << "\n";
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
                << "   Current column:\n" << col1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for index 1
      {
         OCT::Iterator pos = col1.upperBound( 1UL );

         if( pos != col1.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required index = 1\n"
                << "   Current column:\n" << col1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for index 2
      {
         OCT::Iterator pos = col1.upperBound( 2UL );

         if( pos != col1.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required index = 2\n"
                << "   Current column:\n" << col1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c isDefault() function with the Column specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isDefault() function with the Column specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseGeneralTest::testIsDefault()
{
   using blaze::isDefault;


   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major isDefault() function";

      initialize();

      // isDefault with default column
      {
         CT col0 = blaze::column( mat_, 0UL );

         if( isDefault( col0[1] ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Column element: " << col0[1] << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( col0 ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Column:\n" << col0 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with non-default column
      {
         CT col1 = blaze::column( mat_, 1UL );

         if( isDefault( col1[1] ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Column element: " << col1[1] << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( col1 ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Column:\n" << col1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major isDefault() function";

      initialize();

      // isDefault with default column
      {
         OCT col0 = blaze::column( tmat_, 0UL );

         if( isDefault( col0[1] ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Column element: " << col0[1] << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( col0 ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Column:\n" << col0 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with non-default column
      {
         OCT col1 = blaze::column( tmat_, 1UL );

         if( isDefault( col1[1] ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Column element: " << col1[1] << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( col1 ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Column:\n" << col1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c isSame() function with the Column specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isSame() function with the Column specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseGeneralTest::testIsSame()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major isSame() function";

      // isSame with matching columns
      {
         CT col1 = blaze::column( mat_, 1UL );
         CT col2 = blaze::column( mat_, 1UL );

         if( blaze::isSame( col1, col2 ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First column:\n" << col1 << "\n"
                << "   Second column:\n" << col2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with non-matching columns
      {
         CT col1 = blaze::column( mat_, 1UL );
         CT col2 = blaze::column( mat_, 2UL );

         if( blaze::isSame( col1, col2 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First column:\n" << col1 << "\n"
                << "   Second column:\n" << col2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with column and matching subvector
      {
         CT   col1 = blaze::column( mat_, 1UL );
         auto sv   = blaze::subvector( col1, 0UL, 4UL );

         if( blaze::isSame( col1, sv ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Sparse column:\n" << col1 << "\n"
                << "   Sparse subvector:\n" << sv << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( sv, col1 ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Sparse column:\n" << col1 << "\n"
                << "   Sparse subvector:\n" << sv << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with column and non-matching subvector (different size)
      {
         CT   col1 = blaze::column( mat_, 1UL );
         auto sv   = blaze::subvector( col1, 0UL, 3UL );

         if( blaze::isSame( col1, sv ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Sparse column:\n" << col1 << "\n"
                << "   Sparse subvector:\n" << sv << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( sv, col1 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Sparse column:\n" << col1 << "\n"
                << "   Sparse subvector:\n" << sv << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with column and non-matching subvector (different offset)
      {
         CT   col1 = blaze::column( mat_, 1UL );
         auto sv   = blaze::subvector( col1, 1UL, 3UL );

         if( blaze::isSame( col1, sv ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Sparse column:\n" << col1 << "\n"
                << "   Sparse subvector:\n" << sv << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( sv, col1 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Sparse column:\n" << col1 << "\n"
                << "   Sparse subvector:\n" << sv << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with matching columns on a common submatrix
      {
         auto sm   = blaze::submatrix( mat_, 1UL, 1UL, 3UL, 2UL );
         auto col1 = blaze::column( sm, 1UL );
         auto col2 = blaze::column( sm, 1UL );

         if( blaze::isSame( col1, col2 ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First column:\n" << col1 << "\n"
                << "   Second column:\n" << col2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with non-matching columns on a common submatrix
      {
         auto sm   = blaze::submatrix( mat_, 1UL, 1UL, 3UL, 2UL );
         auto col1 = blaze::column( sm, 0UL );
         auto col2 = blaze::column( sm, 1UL );

         if( blaze::isSame( col1, col2 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First column:\n" << col1 << "\n"
                << "   Second column:\n" << col2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with matching columns on matrix and submatrix
      {
         auto sm   = blaze::submatrix( mat_, 0UL, 1UL, 4UL, 3UL );
         auto col1 = blaze::column( mat_, 2UL );
         auto col2 = blaze::column( sm  , 1UL );

         if( blaze::isSame( col1, col2 ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First column:\n" << col1 << "\n"
                << "   Second column:\n" << col2 << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( col2, col1 ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First column:\n" << col1 << "\n"
                << "   Second column:\n" << col2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with non-matching columns on matrix and submatrix (different column)
      {
         auto sm   = blaze::submatrix( mat_, 0UL, 1UL, 4UL, 3UL );
         auto col1 = blaze::column( mat_, 1UL );
         auto col2 = blaze::column( sm  , 1UL );

         if( blaze::isSame( col1, col2 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First column:\n" << col1 << "\n"
                << "   Second column:\n" << col2 << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( col2, col1 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First column:\n" << col1 << "\n"
                << "   Second column:\n" << col2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with non-matching columns on matrix and submatrix (different size)
      {
         auto sm   = blaze::submatrix( mat_, 0UL, 1UL, 3UL, 3UL );
         auto col1 = blaze::column( mat_, 2UL );
         auto col2 = blaze::column( sm  , 1UL );

         if( blaze::isSame( col1, col2 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First column:\n" << col1 << "\n"
                << "   Second column:\n" << col2 << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( col2, col1 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First column:\n" << col1 << "\n"
                << "   Second column:\n" << col2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with matching columns on two submatrices
      {
         auto sm1  = blaze::submatrix( mat_, 0UL, 1UL, 4UL, 3UL );
         auto sm2  = blaze::submatrix( mat_, 0UL, 2UL, 4UL, 3UL );
         auto col1 = blaze::column( sm1, 1UL );
         auto col2 = blaze::column( sm2, 0UL );

         if( blaze::isSame( col1, col2 ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First column:\n" << col1 << "\n"
                << "   Second column:\n" << col2 << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( col2, col1 ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First column:\n" << col1 << "\n"
                << "   Second column:\n" << col2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with non-matching columns on two submatrices (different column)
      {
         auto sm1  = blaze::submatrix( mat_, 0UL, 1UL, 4UL, 3UL );
         auto sm2  = blaze::submatrix( mat_, 0UL, 2UL, 4UL, 3UL );
         auto col1 = blaze::column( sm1, 1UL );
         auto col2 = blaze::column( sm2, 1UL );

         if( blaze::isSame( col1, col2 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First column:\n" << col1 << "\n"
                << "   Second column:\n" << col2 << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( col2, col1 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First column:\n" << col1 << "\n"
                << "   Second column:\n" << col2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with non-matching columns on two submatrices (different size)
      {
         auto sm1  = blaze::submatrix( mat_, 0UL, 1UL, 4UL, 3UL );
         auto sm2  = blaze::submatrix( mat_, 0UL, 2UL, 3UL, 3UL );
         auto col1 = blaze::column( sm1, 1UL );
         auto col2 = blaze::column( sm2, 0UL );

         if( blaze::isSame( col1, col2 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First column:\n" << col1 << "\n"
                << "   Second column:\n" << col2 << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( col2, col1 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First column:\n" << col1 << "\n"
                << "   Second column:\n" << col2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with non-matching columns on two submatrices (different offset)
      {
         auto sm1  = blaze::submatrix( mat_, 0UL, 1UL, 3UL, 3UL );
         auto sm2  = blaze::submatrix( mat_, 1UL, 2UL, 3UL, 3UL );
         auto col1 = blaze::column( sm1, 1UL );
         auto col2 = blaze::column( sm2, 0UL );

         if( blaze::isSame( col1, col2 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First column:\n" << col1 << "\n"
                << "   Second column:\n" << col2 << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( col2, col1 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First column:\n" << col1 << "\n"
                << "   Second column:\n" << col2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with matching column subvectors on a common submatrix
      {
         auto sm   = blaze::submatrix( mat_, 1UL, 1UL, 3UL, 2UL );
         auto col1 = blaze::column( sm, 1UL );
         auto sv1  = blaze::subvector( col1, 0UL, 2UL );
         auto sv2  = blaze::subvector( col1, 0UL, 2UL );

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

      // isSame with non-matching column subvectors on submatrices (different size)
      {
         auto sm   = blaze::submatrix( mat_, 1UL, 1UL, 3UL, 2UL );
         auto col1 = blaze::column( sm, 1UL );
         auto sv1  = blaze::subvector( col1, 0UL, 2UL );
         auto sv2  = blaze::subvector( col1, 0UL, 3UL );

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

      // isSame with non-matching column subvectors on submatrices (different offset)
      {
         auto sm   = blaze::submatrix( mat_, 1UL, 1UL, 3UL, 2UL );
         auto col1 = blaze::column( sm, 1UL );
         auto sv1  = blaze::subvector( col1, 0UL, 2UL );
         auto sv2  = blaze::subvector( col1, 1UL, 2UL );

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

      // isSame with matching column subvectors on two submatrices
      {
         auto sm1  = blaze::submatrix( mat_, 0UL, 1UL, 4UL, 3UL );
         auto sm2  = blaze::submatrix( mat_, 0UL, 2UL, 4UL, 3UL );
         auto col1 = blaze::column( sm1, 1UL );
         auto col2 = blaze::column( sm2, 0UL );
         auto sv1  = blaze::subvector( col1, 0UL, 2UL );
         auto sv2  = blaze::subvector( col2, 0UL, 2UL );

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

      // isSame with non-matching column subvectors on two submatrices (different size)
      {
         auto sm1  = blaze::submatrix( mat_, 0UL, 1UL, 4UL, 3UL );
         auto sm2  = blaze::submatrix( mat_, 0UL, 2UL, 4UL, 3UL );
         auto col1 = blaze::column( sm1, 1UL );
         auto col2 = blaze::column( sm2, 0UL );
         auto sv1  = blaze::subvector( col1, 0UL, 2UL );
         auto sv2  = blaze::subvector( col2, 0UL, 3UL );

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

      // isSame with non-matching column subvectors on two submatrices (different offset)
      {
         auto sm1  = blaze::submatrix( mat_, 0UL, 1UL, 4UL, 3UL );
         auto sm2  = blaze::submatrix( mat_, 0UL, 2UL, 4UL, 3UL );
         auto col1 = blaze::column( sm1, 1UL );
         auto col2 = blaze::column( sm2, 0UL );
         auto sv1  = blaze::subvector( col1, 0UL, 2UL );
         auto sv2  = blaze::subvector( col2, 1UL, 2UL );

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
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major isSame() function";

      // isSame with matching columns
      {
         OCT col1 = blaze::column( tmat_, 1UL );
         OCT col2 = blaze::column( tmat_, 1UL );

         if( blaze::isSame( col1, col2 ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First column:\n" << col1 << "\n"
                << "   Second column:\n" << col2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with non-matching columns
      {
         OCT col1 = blaze::column( tmat_, 1UL );
         OCT col2 = blaze::column( tmat_, 2UL );

         if( blaze::isSame( col1, col2 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First column:\n" << col1 << "\n"
                << "   Second column:\n" << col2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with column and matching subvector
      {
         OCT  col1 = blaze::column( tmat_, 1UL );
         auto sv   = blaze::subvector( col1, 0UL, 4UL );

         if( blaze::isSame( col1, sv ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Sparse column:\n" << col1 << "\n"
                << "   Sparse subvector:\n" << sv << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( sv, col1 ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Sparse column:\n" << col1 << "\n"
                << "   Sparse subvector:\n" << sv << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with column and non-matching subvector (different size)
      {
         OCT  col1 = blaze::column( tmat_, 1UL );
         auto sv   = blaze::subvector( col1, 0UL, 3UL );

         if( blaze::isSame( col1, sv ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Sparse column:\n" << col1 << "\n"
                << "   Sparse subvector:\n" << sv << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( sv, col1 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Sparse column:\n" << col1 << "\n"
                << "   Sparse subvector:\n" << sv << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with column and non-matching subvector (different offset)
      {
         OCT  col1 = blaze::column( tmat_, 1UL );
         auto sv   = blaze::subvector( col1, 1UL, 3UL );

         if( blaze::isSame( col1, sv ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Sparse column:\n" << col1 << "\n"
                << "   Sparse subvector:\n" << sv << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( sv, col1 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Sparse column:\n" << col1 << "\n"
                << "   Sparse subvector:\n" << sv << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with matching columns on a common submatrix
      {
         auto sm   = blaze::submatrix( tmat_, 1UL, 1UL, 3UL, 2UL );
         auto col1 = blaze::column( sm, 1UL );
         auto col2 = blaze::column( sm, 1UL );

         if( blaze::isSame( col1, col2 ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First column:\n" << col1 << "\n"
                << "   Second column:\n" << col2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with non-matching columns on a common submatrix
      {
         auto sm   = blaze::submatrix( tmat_, 1UL, 1UL, 3UL, 2UL );
         auto col1 = blaze::column( sm, 0UL );
         auto col2 = blaze::column( sm, 1UL );

         if( blaze::isSame( col1, col2 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First column:\n" << col1 << "\n"
                << "   Second column:\n" << col2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with matching columns on matrix and submatrix
      {
         auto sm   = blaze::submatrix( tmat_, 0UL, 1UL, 4UL, 3UL );
         auto col1 = blaze::column( tmat_, 2UL );
         auto col2 = blaze::column( sm   , 1UL );

         if( blaze::isSame( col1, col2 ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First column:\n" << col1 << "\n"
                << "   Second column:\n" << col2 << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( col2, col1 ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First column:\n" << col1 << "\n"
                << "   Second column:\n" << col2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with non-matching columns on matrix and submatrix (different column)
      {
         auto sm   = blaze::submatrix( tmat_, 0UL, 1UL, 4UL, 3UL );
         auto col1 = blaze::column( tmat_, 1UL );
         auto col2 = blaze::column( sm   , 1UL );

         if( blaze::isSame( col1, col2 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First column:\n" << col1 << "\n"
                << "   Second column:\n" << col2 << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( col2, col1 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First column:\n" << col1 << "\n"
                << "   Second column:\n" << col2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with non-matching columns on matrix and submatrix (different size)
      {
         auto sm   = blaze::submatrix( tmat_, 0UL, 1UL, 3UL, 3UL );
         auto col1 = blaze::column( tmat_, 2UL );
         auto col2 = blaze::column( sm   , 1UL );

         if( blaze::isSame( col1, col2 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First column:\n" << col1 << "\n"
                << "   Second column:\n" << col2 << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( col2, col1 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First column:\n" << col1 << "\n"
                << "   Second column:\n" << col2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with matching columns on two submatrices
      {
         auto sm1  = blaze::submatrix( tmat_, 0UL, 1UL, 4UL, 3UL );
         auto sm2  = blaze::submatrix( tmat_, 0UL, 2UL, 4UL, 3UL );
         auto col1 = blaze::column( sm1, 1UL );
         auto col2 = blaze::column( sm2, 0UL );

         if( blaze::isSame( col1, col2 ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First column:\n" << col1 << "\n"
                << "   Second column:\n" << col2 << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( col2, col1 ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First column:\n" << col1 << "\n"
                << "   Second column:\n" << col2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with non-matching columns on two submatrices (different column)
      {
         auto sm1  = blaze::submatrix( tmat_, 0UL, 1UL, 4UL, 3UL );
         auto sm2  = blaze::submatrix( tmat_, 0UL, 2UL, 4UL, 3UL );
         auto col1 = blaze::column( sm1, 1UL );
         auto col2 = blaze::column( sm2, 1UL );

         if( blaze::isSame( col1, col2 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First column:\n" << col1 << "\n"
                << "   Second column:\n" << col2 << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( col2, col1 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First column:\n" << col1 << "\n"
                << "   Second column:\n" << col2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with non-matching columns on two submatrices (different size)
      {
         auto sm1  = blaze::submatrix( tmat_, 0UL, 1UL, 4UL, 3UL );
         auto sm2  = blaze::submatrix( tmat_, 0UL, 2UL, 3UL, 3UL );
         auto col1 = blaze::column( sm1, 1UL );
         auto col2 = blaze::column( sm2, 0UL );

         if( blaze::isSame( col1, col2 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First column:\n" << col1 << "\n"
                << "   Second column:\n" << col2 << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( col2, col1 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First column:\n" << col1 << "\n"
                << "   Second column:\n" << col2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with non-matching columns on two submatrices (different offset)
      {
         auto sm1  = blaze::submatrix( tmat_, 0UL, 1UL, 3UL, 3UL );
         auto sm2  = blaze::submatrix( tmat_, 1UL, 2UL, 3UL, 3UL );
         auto col1 = blaze::column( sm1, 1UL );
         auto col2 = blaze::column( sm2, 0UL );

         if( blaze::isSame( col1, col2 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First column:\n" << col1 << "\n"
                << "   Second column:\n" << col2 << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( col2, col1 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First column:\n" << col1 << "\n"
                << "   Second column:\n" << col2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with matching column subvectors on a common submatrix
      {
         auto sm   = blaze::submatrix( tmat_, 1UL, 1UL, 3UL, 2UL );
         auto col1 = blaze::column( sm, 1UL );
         auto sv1  = blaze::subvector( col1, 0UL, 2UL );
         auto sv2  = blaze::subvector( col1, 0UL, 2UL );

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

      // isSame with non-matching column subvectors on submatrices (different size)
      {
         auto sm   = blaze::submatrix( tmat_, 1UL, 1UL, 3UL, 2UL );
         auto col1 = blaze::column( sm, 1UL );
         auto sv1  = blaze::subvector( col1, 0UL, 2UL );
         auto sv2  = blaze::subvector( col1, 0UL, 3UL );

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

      // isSame with non-matching column subvectors on submatrices (different offset)
      {
         auto sm   = blaze::submatrix( tmat_, 1UL, 1UL, 3UL, 2UL );
         auto col1 = blaze::column( sm, 1UL );
         auto sv1  = blaze::subvector( col1, 0UL, 2UL );
         auto sv2  = blaze::subvector( col1, 1UL, 2UL );

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

      // isSame with matching column subvectors on two submatrices
      {
         auto sm1  = blaze::submatrix( tmat_, 0UL, 1UL, 4UL, 3UL );
         auto sm2  = blaze::submatrix( tmat_, 0UL, 2UL, 4UL, 3UL );
         auto col1 = blaze::column( sm1, 1UL );
         auto col2 = blaze::column( sm2, 0UL );
         auto sv1  = blaze::subvector( col1, 0UL, 2UL );
         auto sv2  = blaze::subvector( col2, 0UL, 2UL );

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

      // isSame with non-matching column subvectors on two submatrices (different size)
      {
         auto sm1  = blaze::submatrix( tmat_, 0UL, 1UL, 4UL, 3UL );
         auto sm2  = blaze::submatrix( tmat_, 0UL, 2UL, 4UL, 3UL );
         auto col1 = blaze::column( sm1, 1UL );
         auto col2 = blaze::column( sm2, 0UL );
         auto sv1  = blaze::subvector( col1, 0UL, 2UL );
         auto sv2  = blaze::subvector( col2, 0UL, 3UL );

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

      // isSame with non-matching column subvectors on two submatrices (different offset)
      {
         auto sm1  = blaze::submatrix( tmat_, 0UL, 1UL, 4UL, 3UL );
         auto sm2  = blaze::submatrix( tmat_, 0UL, 2UL, 4UL, 3UL );
         auto col1 = blaze::column( sm1, 1UL );
         auto col2 = blaze::column( sm2, 0UL );
         auto sv1  = blaze::subvector( col1, 0UL, 2UL );
         auto sv2  = blaze::subvector( col2, 1UL, 2UL );

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
/*!\brief Test of the \c subvector() function with the Column specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c subvector() function used with the Column
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseGeneralTest::testSubvector()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major subvector() function";

      initialize();

      {
         CT   col1 = blaze::column( mat_, 1UL );
         auto sv   = blaze::subvector( col1, 0UL, 4UL );

         if( sv[1] != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subscript operator access failed\n"
                << " Details:\n"
                << "   Result: " << sv[1] << "\n"
                << "   Expected result: 1\n";
            throw std::runtime_error( oss.str() );
         }

         if( sv.begin()->value() != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator access failed\n"
                << " Details:\n"
                << "   Result: " << sv.begin()->value() << "\n"
                << "   Expected result: 1\n";
            throw std::runtime_error( oss.str() );
         }
      }

      try {
         CT   col1 = blaze::column( mat_, 1UL );
         auto sv   = blaze::subvector( col1, 4UL, 4UL );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of out-of-bounds subvector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      try {
         CT   col1 = blaze::column( mat_, 1UL );
         auto sv   = blaze::subvector( col1, 0UL, 5UL );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of out-of-bounds subvector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major subvector() function";

      initialize();

      {
         OCT  col1 = blaze::column( tmat_, 1UL );
         auto sv   = blaze::subvector( col1, 0UL, 4UL );

         if( sv[1] != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subscript operator access failed\n"
                << " Details:\n"
                << "   Result: " << sv[1] << "\n"
                << "   Expected result: 1\n";
            throw std::runtime_error( oss.str() );
         }

         if( sv.begin()->value() != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator access failed\n"
                << " Details:\n"
                << "   Result: " << sv.begin()->value() << "\n"
                << "   Expected result: 1\n";
            throw std::runtime_error( oss.str() );
         }
      }

      try {
         OCT  col1 = blaze::column( tmat_, 1UL );
         auto sv   = blaze::subvector( col1, 4UL, 4UL );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of out-of-bounds subvector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      try {
         OCT  col1 = blaze::column( tmat_, 1UL );
         auto sv   = blaze::subvector( col1, 0UL, 5UL );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of out-of-bounds subvector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c elements() function with the Column specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c elements() function used with the Column
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseGeneralTest::testElements()
{
   //=====================================================================================
   // Row-major matrix tests (initializer_list)
   //=====================================================================================

   {
      test_ = "Row-major elements() function (initializer_list)";

      initialize();

      {
         CT   col2 = blaze::column( mat_, 2UL );
         auto e    = blaze::elements( col2, { 2UL, 0UL } );

         if( e[1] != -2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subscript operator access failed\n"
                << " Details:\n"
                << "   Result: " << e[1] << "\n"
                << "   Expected result: -2\n";
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
         CT   col2 = blaze::column( mat_, 2UL );
         auto e    = blaze::elements( col2, { 4UL } );

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
   // Row-major matrix tests (std::array)
   //=====================================================================================

   {
      test_ = "Row-major elements() function (std::array)";

      initialize();

      {
         std::array<int,2UL> indices{ 2UL, 0UL };

         CT   col2 = blaze::column( mat_, 2UL );
         auto e    = blaze::elements( col2, indices );

         if( e[1] != -2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subscript operator access failed\n"
                << " Details:\n"
                << "   Result: " << e[1] << "\n"
                << "   Expected result: -2\n";
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
         std::array<int,1UL> indices{ 4UL };

         CT   col2 = blaze::column( mat_, 2UL );
         auto e    = blaze::elements( col2, indices );

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
   // Row-major matrix tests (lambda expression)
   //=====================================================================================

   {
      test_ = "Row-major elements() function (lambda expression)";

      initialize();

      {
         CT   col2 = blaze::column( mat_, 2UL );
         auto e    = blaze::elements( col2, []( size_t i ){ return 2UL-2UL*i; }, 2UL );

         if( e[1] != -2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subscript operator access failed\n"
                << " Details:\n"
                << "   Result: " << e[1] << "\n"
                << "   Expected result: -2\n";
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
         CT   col2 = blaze::column( mat_, 2UL );
         auto e    = blaze::elements( col2, []( size_t ){ return 4UL; }, 1UL );

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
   // Column-major matrix tests (initializer_list)
   //=====================================================================================

   {
      test_ = "Column-major elements() function (initializer_list)";

      initialize();

      {
         OCT  col2 = blaze::column( tmat_, 2UL );
         auto e    = blaze::elements( col2, { 2UL, 0UL } );

         if( e[1] != -2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subscript operator access failed\n"
                << " Details:\n"
                << "   Result: " << e[1] << "\n"
                << "   Expected result: -2\n";
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
         OCT  col2 = blaze::column( tmat_, 2UL );
         auto e    = blaze::elements( col2, { 4UL } );

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
   // Column-major matrix tests (std::array)
   //=====================================================================================

   {
      test_ = "Column-major elements() function (std::array)";

      initialize();

      {
         std::array<int,2UL> indices{ 2UL, 0UL };

         OCT  col2 = blaze::column( tmat_, 2UL );
         auto e    = blaze::elements( col2, indices );

         if( e[1] != -2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subscript operator access failed\n"
                << " Details:\n"
                << "   Result: " << e[1] << "\n"
                << "   Expected result: -2\n";
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
         std::array<int,1UL> indices{ 4UL };

         OCT  col2 = blaze::column( tmat_, 2UL );
         auto e    = blaze::elements( col2, indices );

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
   // Column-major matrix tests (lambda expression)
   //=====================================================================================

   {
      test_ = "Column-major elements() function (lambda expression)";

      initialize();

      {
         OCT  col2 = blaze::column( tmat_, 2UL );
         auto e    = blaze::elements( col2, []( size_t i ){ return 2UL-2UL*i; }, 2UL );

         if( e[1] != -2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subscript operator access failed\n"
                << " Details:\n"
                << "   Result: " << e[1] << "\n"
                << "   Expected result: -2\n";
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
         OCT  col2 = blaze::column( tmat_, 2UL );
         auto e    = blaze::elements( col2, []( size_t ){ return 4UL; }, 1UL );

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
/*!\brief Initialization of all member matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function initializes all member matrices to specific predetermined values.
*/
void SparseGeneralTest::initialize()
{
   // Initializing the row-major compressed matrix
   mat_.reset();
   mat_(1,1) =  1;
   mat_(0,2) = -2;
   mat_(2,2) = -3;
   mat_(1,3) =  4;
   mat_(2,3) =  5;
   mat_(3,3) = -6;
   mat_(0,4) =  7;
   mat_(1,4) = -8;
   mat_(2,4) =  9;
   mat_(3,4) = 10;

   // Initializing the column-major compressed matrix
   tmat_.reset();
   tmat_(1,1) =  1;
   tmat_(0,2) = -2;
   tmat_(2,2) = -3;
   tmat_(1,3) =  4;
   tmat_(2,3) =  5;
   tmat_(3,3) = -6;
   tmat_(0,4) =  7;
   tmat_(1,4) = -8;
   tmat_(2,4) =  9;
   tmat_(3,4) = 10;
}
//*************************************************************************************************

} // namespace column

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
   std::cout << "   Running Column sparse general test..." << std::endl;

   try
   {
      RUN_COLUMN_SPARSEGENERAL_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during Column sparse general test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
