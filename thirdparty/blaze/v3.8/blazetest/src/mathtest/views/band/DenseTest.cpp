//=================================================================================================
/*!
//  \file src/mathtest/views/band/DenseTest.cpp
//  \brief Source file for the Band dense test
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
#include <blaze/math/Views.h>
#include <blaze/util/policies/Deallocate.h>
#include <blazetest/mathtest/views/band/DenseTest.h>

#ifdef BLAZE_USE_HPX_THREADS
#  include <hpx/hpx_main.hpp>
#endif


namespace blazetest {

namespace mathtest {

namespace views {

namespace band {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the Band dense test.
//
// \exception std::runtime_error Operation error detected.
*/
DenseTest::DenseTest()
   : mat_ ( 4UL, 6UL )
   , tmat_( 6UL, 4UL )
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
/*!\brief Test of the Band constructors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all constructors of the Band specialization. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testConstructors()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major Band constructor (0x0)";

      MT mat;

      // 1st lower matrix band
      try {
         blaze::band( mat, -1L );
      }
      catch( std::invalid_argument& ) {}

      // 0th matrix band (diagonal)
      {
         BT band0 = blaze::band( mat, 0L );

         checkSize    ( band0, 0UL );
         checkCapacity( band0, 0UL );
         checkNonZeros( band0, 0UL );
      }

      // 1st upper matrix band
      try {
         blaze::band( mat, 1L );
      }
      catch( std::invalid_argument& ) {}
   }

   {
      test_ = "Row-major Band constructor (0x2)";

      MT mat( 0UL, 2UL );

      // 1st lower matrix band
      try {
         blaze::band( mat, -1L );
      }
      catch( std::invalid_argument& ) {}

      // 0th matrix band (diagonal)
      {
         BT band0 = blaze::band( mat, 0L );

         checkSize    ( band0, 0UL );
         checkCapacity( band0, 0UL );
         checkNonZeros( band0, 0UL );
      }

      // 1st upper matrix
      {
         BT band1 = blaze::band( mat, 1L );

         checkSize    ( band1, 0UL );
         checkCapacity( band1, 0UL );
         checkNonZeros( band1, 0UL );
      }

      // 2nd upper matrix band
      try {
         blaze::band( mat, 2L );
      }
      catch( std::invalid_argument& ) {}
   }

   {
      test_ = "Row-major Band constructor (2x0)";

      MT mat( 2UL, 0UL );

      // 2nd lower matrix band
      try {
         blaze::band( mat, -2L );
      }
      catch( std::invalid_argument& ) {}

      // 1st lower matrix band
      {
         BT band1 = blaze::band( mat, -1L );

         checkSize    ( band1, 0UL );
         checkCapacity( band1, 0UL );
         checkNonZeros( band1, 0UL );
      }

      // 0th matrix band (diagonal)
      {
         BT band0 = blaze::band( mat, 0L );

         checkSize    ( band0, 0UL );
         checkCapacity( band0, 0UL );
         checkNonZeros( band0, 0UL );
      }

      // 1st upper matrix band
      try {
         blaze::band( mat, 1L );
      }
      catch( std::invalid_argument& ) {}
   }

   {
      test_ = "Row-major Band constructor (4x6)";

      initialize();

      // 4th lower matrix band
      try {
         blaze::band( mat_, -4L );
      }
      catch( std::invalid_argument& ) {}

      // 3rd lower matrix band
      {
         BT band3 = blaze::band( mat_, -3L );

         checkSize    ( band3, 1UL );
         checkCapacity( band3, 1UL );
         checkNonZeros( band3, 0UL );

         if( band3[0] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of 3rd lower band failed\n"
                << " Details:\n"
                << "   Result:\n" << band3 << "\n"
                << "   Expected result:\n( 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // 2nd lower matrix band
      {
         BT band2 = blaze::band( mat_, -2L );

         checkSize    ( band2, 2UL );
         checkCapacity( band2, 2UL );
         checkNonZeros( band2, 0UL );

         if( band2[0] != 0 || band2[1] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of 2nd lower band failed\n"
                << " Details:\n"
                << "   Result:\n" << band2 << "\n"
                << "   Expected result:\n( 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // 1st lower matrix band
      {
         BT band1 = blaze::band( mat_, -1L );

         checkSize    ( band1, 3UL );
         checkCapacity( band1, 3UL );
         checkNonZeros( band1, 1UL );

         if( band1[0] != 0 || band1[1] != 1 || band1[2] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of 1st lower band failed\n"
                << " Details:\n"
                << "   Result:\n" << band1 << "\n"
                << "   Expected result:\n( 0 1 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // 0th matrix band (diagonal)
      {
         BT band0 = blaze::band( mat_, 0L );

         checkSize    ( band0, 4UL );
         checkCapacity( band0, 4UL );
         checkNonZeros( band0, 2UL );

         if( band0[0] != -2 || band0[1] != 0 || band0[2] != -3 || band0[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of 0th band (diagonal) failed\n"
                << " Details:\n"
                << "   Result:\n" << band0 << "\n"
                << "   Expected result:\n( -2 0 -3 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // 1st upper matrix band
      {
         BT band1 = blaze::band( mat_, 1L );

         checkSize    ( band1, 4UL );
         checkCapacity( band1, 4UL );
         checkNonZeros( band1, 3UL );

         if( band1[0] != 0 || band1[1] != 4 || band1[2] != 5 || band1[3] != -6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of 1st upper band failed\n"
                << " Details:\n"
                << "   Result:\n" << band1 << "\n"
                << "   Expected result:\n( 0 4 5 -6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // 2nd upper matrix band
      {
         BT band2 = blaze::band( mat_, 2L );

         checkSize    ( band2, 4UL );
         checkCapacity( band2, 4UL );
         checkNonZeros( band2, 4UL );

         if( band2[0] != 7 || band2[1] != -8 || band2[2] != 9 || band2[3] != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of 2nd upper band failed\n"
                << " Details:\n"
                << "   Result:\n" << band2 << "\n"
                << "   Expected result:\n( 7 -8 9 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // 3rd upper matrix band
      {
         BT band3 = blaze::band( mat_, 3L );

         checkSize    ( band3, 3UL );
         checkCapacity( band3, 3UL );
         checkNonZeros( band3, 0UL );

         if( band3[0] != 0 || band3[1] != 0 || band3[2] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of 3rd upper band failed\n"
                << " Details:\n"
                << "   Result:\n" << band3 << "\n"
                << "   Expected result:\n( 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // 4th upper matrix band
      {
         BT band4 = blaze::band( mat_, 4L );

         checkSize    ( band4, 2UL );
         checkCapacity( band4, 2UL );
         checkNonZeros( band4, 0UL );

         if( band4[0] != 0 || band4[1] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of 4th upper band failed\n"
                << " Details:\n"
                << "   Result:\n" << band4 << "\n"
                << "   Expected result:\n( 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // 5th upper matrix band
      {
         BT band5 = blaze::band( mat_, 5L );

         checkSize    ( band5, 1UL );
         checkCapacity( band5, 1UL );
         checkNonZeros( band5, 0UL );

         if( band5[0] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of 5th upper band failed\n"
                << " Details:\n"
                << "   Result:\n" << band5 << "\n"
                << "   Expected result:\n( 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // 6th upper matrix band
      try {
         blaze::band( mat_, 6L );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major Band constructor (0x0)";

      OMT tmat;

      // 1st lower matrix band
      try {
         blaze::band( tmat, -1L );
      }
      catch( std::invalid_argument& ) {}

      // 0th matrix band (diagonal)
      {
         OBT band0 = blaze::band( tmat, 0L );

         checkSize    ( band0, 0UL );
         checkCapacity( band0, 0UL );
         checkNonZeros( band0, 0UL );
      }

      // 1st upper matrix band
      try {
         blaze::band( tmat, 1L );
      }
      catch( std::invalid_argument& ) {}
   }

   {
      test_ = "Column-major Band constructor (0x2)";

      OMT tmat( 0UL, 2UL );

      // 1st lower matrix band
      try {
         blaze::band( tmat, -1L );
      }
      catch( std::invalid_argument& ) {}

      // 0th matrix band (diagonal)
      {
         OBT band0 = blaze::band( tmat, 0L );

         checkSize    ( band0, 0UL );
         checkCapacity( band0, 0UL );
         checkNonZeros( band0, 0UL );
      }

      // 1st upper matrix
      {
         OBT band1 = blaze::band( tmat, 1L );

         checkSize    ( band1, 0UL );
         checkCapacity( band1, 0UL );
         checkNonZeros( band1, 0UL );
      }

      // 2nd upper matrix band
      try {
         blaze::band( tmat, 2L );
      }
      catch( std::invalid_argument& ) {}
   }

   {
      test_ = "Column-major Band constructor (2x0)";

      OMT tmat( 2UL, 0UL );

      // 2nd lower matrix band
      try {
         blaze::band( tmat, -2L );
      }
      catch( std::invalid_argument& ) {}

      // 1st lower matrix band
      {
         OBT band1 = blaze::band( tmat, -1L );

         checkSize    ( band1, 0UL );
         checkCapacity( band1, 0UL );
         checkNonZeros( band1, 0UL );
      }

      // 0th matrix band (diagonal)
      {
         OBT band0 = blaze::band( tmat, 0L );

         checkSize    ( band0, 0UL );
         checkCapacity( band0, 0UL );
         checkNonZeros( band0, 0UL );
      }

      // 1st upper matrix band
      try {
         blaze::band( tmat, 1L );
      }
      catch( std::invalid_argument& ) {}
   }

   {
      test_ = "Column-major Band constructor (6x4)";

      initialize();

      // 6th lower matrix band
      try {
         blaze::band( tmat_, -6L );
      }
      catch( std::invalid_argument& ) {}

      // 5th lower matrix band
      {
         OBT band5 = blaze::band( tmat_, -5L );

         checkSize    ( band5, 1UL );
         checkCapacity( band5, 1UL );
         checkNonZeros( band5, 0UL );

         if( band5[0] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of 5th lower band failed\n"
                << " Details:\n"
                << "   Result:\n" << band5 << "\n"
                << "   Expected result:\n( 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // 4th lower matrix band
      {
         OBT band4 = blaze::band( tmat_, -4L );

         checkSize    ( band4, 2UL );
         checkCapacity( band4, 2UL );
         checkNonZeros( band4, 0UL );

         if( band4[0] != 0 || band4[1] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of 4th lower band failed\n"
                << " Details:\n"
                << "   Result:\n" << band4 << "\n"
                << "   Expected result:\n( 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // 3rd lower matrix band
      {
         OBT band3 = blaze::band( tmat_, -3L );

         checkSize    ( band3, 3UL );
         checkCapacity( band3, 3UL );
         checkNonZeros( band3, 0UL );

         if( band3[0] != 0 || band3[1] != 0 || band3[2] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of 3rd lower band failed\n"
                << " Details:\n"
                << "   Result:\n" << band3 << "\n"
                << "   Expected result:\n( 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // 2nd lower matrix band
      {
         OBT band2 = blaze::band( tmat_, -2L );

         checkSize    ( band2, 4UL );
         checkCapacity( band2, 4UL );
         checkNonZeros( band2, 4UL );

         if( band2[0] != 7 || band2[1] != -8 || band2[2] != 9 || band2[3] != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of 2nd lower band failed\n"
                << " Details:\n"
                << "   Result:\n" << band2 << "\n"
                << "   Expected result:\n( 7 -8 9 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // 1st lower matrix band
      {
         OBT band1 = blaze::band( tmat_, -1L );

         checkSize    ( band1, 4UL );
         checkCapacity( band1, 4UL );
         checkNonZeros( band1, 3UL );

         if( band1[0] != 0 || band1[1] != 4 || band1[2] != 5 || band1[3] != -6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of 1st lower band failed\n"
                << " Details:\n"
                << "   Result:\n" << band1 << "\n"
                << "   Expected result:\n( 0 4 5 -6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // 0th matrix band (diagonal)
      {
         OBT band0 = blaze::band( tmat_, 0L );

         checkSize    ( band0, 4UL );
         checkCapacity( band0, 4UL );
         checkNonZeros( band0, 2UL );

         if( band0[0] != -2 || band0[1] != 0 || band0[2] != -3 || band0[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of 0th band (diagonal) failed\n"
                << " Details:\n"
                << "   Result:\n" << band0 << "\n"
                << "   Expected result:\n( -2 0 -3 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // 1st upper matrix band
      {
         OBT band1 = blaze::band( tmat_, 1L );

         checkSize    ( band1, 3UL );
         checkCapacity( band1, 3UL );
         checkNonZeros( band1, 1UL );

         if( band1[0] != 0 || band1[1] != 1 || band1[2] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of 1st upper band failed\n"
                << " Details:\n"
                << "   Result:\n" << band1 << "\n"
                << "   Expected result:\n( 0 1 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // 2nd upper matrix band
      {
         OBT band2 = blaze::band( tmat_, 2L );

         checkSize    ( band2, 2UL );
         checkCapacity( band2, 2UL );
         checkNonZeros( band2, 0UL );

         if( band2[0] != 0 || band2[1] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of 2nd upper band failed\n"
                << " Details:\n"
                << "   Result:\n" << band2 << "\n"
                << "   Expected result:\n( 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // 3rd upper matrix band
      {
         OBT band3 = blaze::band( tmat_, 3L );

         checkSize    ( band3, 1UL );
         checkCapacity( band3, 1UL );
         checkNonZeros( band3, 0UL );

         if( band3[0] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of 3rd upper band failed\n"
                << " Details:\n"
                << "   Result:\n" << band3 << "\n"
                << "   Expected result:\n( 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // 4th upper matrix band
      try {
         blaze::band( tmat_, 4L );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Band assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all assignment operators of the Band specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testAssignment()
{
   //=====================================================================================
   // Row-major homogeneous assignment
   //=====================================================================================

   {
      test_ = "Row-major Band homogeneous assignment";

      initialize();

      BT band1 = blaze::band( mat_, -1L );
      band1 = 8;

      checkSize    ( band1,  3UL );
      checkCapacity( band1,  3UL );
      checkNonZeros( band1,  3UL );
      checkRows    ( mat_ ,  4UL );
      checkColumns ( mat_ ,  6UL );
      checkNonZeros( mat_ , 12UL );

      if( band1[0] != 8 || band1[1] != 8 || band1[2] != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << band1 << "\n"
             << "   Expected result:\n( 8 8 8 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != -2 || mat_(0,1) != 0 || mat_(0,2) !=  7 || mat_(0,3) !=  0 || mat_(0,4) !=  0 || mat_(0,5) !=  0 ||
          mat_(1,0) !=  8 || mat_(1,1) != 0 || mat_(1,2) !=  4 || mat_(1,3) != -8 || mat_(1,4) !=  0 || mat_(1,5) !=  0 ||
          mat_(2,0) !=  0 || mat_(2,1) != 8 || mat_(2,2) != -3 || mat_(2,3) !=  5 || mat_(2,4) !=  9 || mat_(2,5) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) != 0 || mat_(3,2) !=  8 || mat_(3,3) !=  0 || mat_(3,4) != -6 || mat_(3,5) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( -2  0  7  0  0  0 )\n"
                                     "(  8  0  4 -8  0  0 )\n"
                                     "(  0  8 -3  5  9  0 )\n"
                                     "(  0  0  8  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major list assignment
   //=====================================================================================

   {
      test_ = "Row-major initializer list assignment (complete list)";

      initialize();

      BT band1 = blaze::band( mat_, 1L );
      band1 = { 1, 2, 3, 4 };

      checkSize    ( band1,  4UL );
      checkCapacity( band1,  4UL );
      checkNonZeros( band1,  4UL );
      checkRows    ( mat_ ,  4UL );
      checkColumns ( mat_ ,  6UL );
      checkNonZeros( mat_ , 11UL );

      if( band1[0] != 1 || band1[1] != 2 || band1[2] != 3 || band1[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << band1 << "\n"
             << "   Expected result:\n( 1 2 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != -2 || mat_(0,1) != 1 || mat_(0,2) !=  7 || mat_(0,3) !=  0 || mat_(0,4) != 0 || mat_(0,5) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 0 || mat_(1,2) !=  2 || mat_(1,3) != -8 || mat_(1,4) != 0 || mat_(1,5) !=  0 ||
          mat_(2,0) !=  0 || mat_(2,1) != 1 || mat_(2,2) != -3 || mat_(2,3) !=  3 || mat_(2,4) != 9 || mat_(2,5) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) !=  0 || mat_(3,4) != 4 || mat_(3,5) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( -2  1  7  0  0  0 )\n"
                                     "(  0  0  2 -8  0  0 )\n"
                                     "(  0  1 -3  3  9  0 )\n"
                                     "(  0  0  0  0  4 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major initializer list assignment (incomplete list)";

      initialize();

      BT band1 = blaze::band( mat_, 1L );
      band1 = { 1, 2 };

      checkSize    ( band1, 4UL );
      checkCapacity( band1, 4UL );
      checkNonZeros( band1, 2UL );
      checkRows    ( mat_ , 4UL );
      checkColumns ( mat_ , 6UL );
      checkNonZeros( mat_ , 9UL );

      if( band1[0] != 1 || band1[1] != 2 || band1[2] != 0 || band1[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << band1 << "\n"
             << "   Expected result:\n( 1 2 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != -2 || mat_(0,1) != 1 || mat_(0,2) !=  7 || mat_(0,3) !=  0 || mat_(0,4) != 0 || mat_(0,5) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 0 || mat_(1,2) !=  2 || mat_(1,3) != -8 || mat_(1,4) != 0 || mat_(1,5) !=  0 ||
          mat_(2,0) !=  0 || mat_(2,1) != 1 || mat_(2,2) != -3 || mat_(2,3) !=  0 || mat_(2,4) != 9 || mat_(2,5) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) !=  0 || mat_(3,4) != 0 || mat_(3,5) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( -2  1  7  0  0  0 )\n"
                                     "(  0  0  2 -8  0  0 )\n"
                                     "(  0  1 -3  0  9  0 )\n"
                                     "(  0  0  0  0  0 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major copy assignment
   //=====================================================================================

   {
      test_ = "Row-major Band copy assignment";

      initialize();

      BT band0 = blaze::band( mat_, 0L );
      band0 = blaze::band( mat_, 1L );

      checkSize    ( band0,  4UL );
      checkCapacity( band0,  4UL );
      checkNonZeros( band0,  3UL );
      checkRows    ( mat_ ,  4UL );
      checkColumns ( mat_ ,  6UL );
      checkNonZeros( mat_ , 11UL );

      if( band0[0] != 0 || band0[1] != 4 || band0[2] != 5 || band0[3] != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( 0  4  5 -6 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) != 7 || mat_(0,3) !=  0 || mat_(0,4) !=  0 || mat_(0,5) !=  0 ||
          mat_(1,0) != 0 || mat_(1,1) != 4 || mat_(1,2) != 4 || mat_(1,3) != -8 || mat_(1,4) !=  0 || mat_(1,5) !=  0 ||
          mat_(2,0) != 0 || mat_(2,1) != 1 || mat_(2,2) != 5 || mat_(2,3) !=  5 || mat_(2,4) !=  9 || mat_(2,5) !=  0 ||
          mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) != 0 || mat_(3,3) != -6 || mat_(3,4) != -6 || mat_(3,5) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0  7  0  0  0 )\n"
                                     "( 0  4  4 -8  0  0 )\n"
                                     "( 0  1  5  5  9  0 )\n"
                                     "( 0  0  0 -6 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major dense vector assignment
   //=====================================================================================

   {
      test_ = "Row-major dense vector assignment (mixed type)";

      initialize();

      BT band1 = blaze::band( mat_, -1L );

      const blaze::DynamicVector<short,blaze::columnVector> vec1{ 8, 0, 9 };

      band1 = vec1;

      checkSize    ( band1,  3UL );
      checkCapacity( band1,  3UL );
      checkNonZeros( band1,  2UL );
      checkRows    ( mat_ ,  4UL );
      checkColumns ( mat_ ,  6UL );
      checkNonZeros( mat_ , 11UL );

      if( band1[0] != 8 || band1[1] != 0 || band1[2] != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" <<  band1 << "\n"
             << "   Expected result:\n( 8 0 9 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != -2 || mat_(0,1) != 0 || mat_(0,2) !=  7 || mat_(0,3) !=  0 || mat_(0,4) !=  0 || mat_(0,5) !=  0 ||
          mat_(1,0) !=  8 || mat_(1,1) != 0 || mat_(1,2) !=  4 || mat_(1,3) != -8 || mat_(1,4) !=  0 || mat_(1,5) !=  0 ||
          mat_(2,0) !=  0 || mat_(2,1) != 0 || mat_(2,2) != -3 || mat_(2,3) !=  5 || mat_(2,4) !=  9 || mat_(2,5) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) != 0 || mat_(3,2) !=  9 || mat_(3,3) !=  0 || mat_(3,4) != -6 || mat_(3,5) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( -2  0  7  0  0  0 )\n"
                                     "(  8  0  4 -8  0  0 )\n"
                                     "(  0  0 -3  5  9  0 )\n"
                                     "(  0  0  9  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major dense vector assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::columnVector;

      initialize();

      BT band1 = blaze::band( mat_, -1L );

      using AlignedPadded = blaze::CustomVector<int,aligned,padded,columnVector>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 16UL ) );
      AlignedPadded vec1( memory.get(), 3UL, 16UL );
      vec1[0] = 8;
      vec1[1] = 0;
      vec1[2] = 9;

      band1 = vec1;

      checkSize    ( band1,  3UL );
      checkCapacity( band1,  3UL );
      checkNonZeros( band1,  2UL );
      checkRows    ( mat_ ,  4UL );
      checkColumns ( mat_ ,  6UL );
      checkNonZeros( mat_ , 11UL );

      if( band1[0] != 8 || band1[1] != 0 || band1[2] != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << band1 << "\n"
             << "   Expected result:\n( 8 0 9 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != -2 || mat_(0,1) != 0 || mat_(0,2) !=  7 || mat_(0,3) !=  0 || mat_(0,4) !=  0 || mat_(0,5) !=  0 ||
          mat_(1,0) !=  8 || mat_(1,1) != 0 || mat_(1,2) !=  4 || mat_(1,3) != -8 || mat_(1,4) !=  0 || mat_(1,5) !=  0 ||
          mat_(2,0) !=  0 || mat_(2,1) != 0 || mat_(2,2) != -3 || mat_(2,3) !=  5 || mat_(2,4) !=  9 || mat_(2,5) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) != 0 || mat_(3,2) !=  9 || mat_(3,3) !=  0 || mat_(3,4) != -6 || mat_(3,5) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( -2  0  7  0  0  0 )\n"
                                     "(  8  0  4 -8  0  0 )\n"
                                     "(  0  0 -3  5  9  0 )\n"
                                     "(  0  0  9  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major dense vector assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::columnVector;

      initialize();

      BT band1 = blaze::band( mat_, -1L );

      using UnalignedUnpadded = blaze::CustomVector<int,unaligned,unpadded,columnVector>;
      std::unique_ptr<int[]> memory( new int[4] );
      UnalignedUnpadded vec1( memory.get()+1UL, 3UL );
      vec1[0] = 8;
      vec1[1] = 0;
      vec1[2] = 9;

      band1 = vec1;

      checkSize    ( band1,  3UL );
      checkCapacity( band1,  3UL );
      checkNonZeros( band1,  2UL );
      checkRows    ( mat_ ,  4UL );
      checkColumns ( mat_ ,  6UL );
      checkNonZeros( mat_ , 11UL );

      if( band1[0] != 8 || band1[1] != 0 || band1[2] != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << band1 << "\n"
             << "   Expected result:\n( 8 0 9 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != -2 || mat_(0,1) != 0 || mat_(0,2) !=  7 || mat_(0,3) !=  0 || mat_(0,4) !=  0 || mat_(0,5) !=  0 ||
          mat_(1,0) !=  8 || mat_(1,1) != 0 || mat_(1,2) !=  4 || mat_(1,3) != -8 || mat_(1,4) !=  0 || mat_(1,5) !=  0 ||
          mat_(2,0) !=  0 || mat_(2,1) != 0 || mat_(2,2) != -3 || mat_(2,3) !=  5 || mat_(2,4) !=  9 || mat_(2,5) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) != 0 || mat_(3,2) !=  9 || mat_(3,3) !=  0 || mat_(3,4) != -6 || mat_(3,5) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( -2  0  7  0  0  0 )\n"
                                     "(  8  0  4 -8  0  0 )\n"
                                     "(  0  0 -3  5  9  0 )\n"
                                     "(  0  0  9  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse vector assignment
   //=====================================================================================

   {
      test_ = "Row-major sparse vector assignment";

      initialize();

      BT band2 = blaze::band( mat_, 2L );

      blaze::CompressedVector<int,blaze::columnVector> vec1( 4UL );
      vec1[3] = 9;

      band2 = vec1;

      checkSize    ( band2, 4UL );
      checkCapacity( band2, 4UL );
      checkNonZeros( band2, 1UL );
      checkRows    ( mat_ , 4UL );
      checkColumns ( mat_ , 6UL );
      checkNonZeros( mat_ , 7UL );

      if( band2[0] != 0 || band2[1] != 0 || band2[2] != 0 || band2[3] != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << band2 << "\n"
             << "   Expected result:\n( 0 0 0 9 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != -2 || mat_(0,1) != 0 || mat_(0,2) !=  0 || mat_(0,3) != 0 || mat_(0,4) !=  0 || mat_(0,5) != 0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 0 || mat_(1,2) !=  4 || mat_(1,3) != 0 || mat_(1,4) !=  0 || mat_(1,5) != 0 ||
          mat_(2,0) !=  0 || mat_(2,1) != 1 || mat_(2,2) != -3 || mat_(2,3) != 5 || mat_(2,4) !=  0 || mat_(2,5) != 0 ||
          mat_(3,0) !=  0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) != 0 || mat_(3,4) != -6 || mat_(3,5) != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( -2  0  0  0  0  0 )\n"
                                     "(  0  0  4  0  0  0 )\n"
                                     "(  0  1 -3  5  0  0 )\n"
                                     "(  0  0  0  0 -6  9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major homogeneous assignment
   //=====================================================================================

   {
      test_ = "Column-major Band homogeneous assignment";

      initialize();

      OBT band1 = blaze::band( tmat_, 1L );
      band1 = 8;

      checkSize    ( band1,  3UL );
      checkCapacity( band1,  3UL );
      checkNonZeros( band1,  3UL );
      checkRows    ( tmat_,  6UL );
      checkColumns ( tmat_,  4UL );
      checkNonZeros( tmat_, 12UL );

      if( band1[0] != 8 || band1[1] != 8 || band1[2] != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << band1 << "\n"
             << "   Expected result:\n( 8 8 8 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != -2 || tmat_(0,1) !=  8 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) !=  0 || tmat_(1,2) !=  8 || tmat_(1,3) !=  0 ||
          tmat_(2,0) !=  7 || tmat_(2,1) !=  4 || tmat_(2,2) != -3 || tmat_(2,3) !=  8 ||
          tmat_(3,0) !=  0 || tmat_(3,1) != -8 || tmat_(3,2) !=  5 || tmat_(3,3) !=  0 ||
          tmat_(4,0) !=  0 || tmat_(4,1) !=  0 || tmat_(4,2) !=  9 || tmat_(4,3) != -6 ||
          tmat_(5,0) !=  0 || tmat_(5,1) !=  0 || tmat_(5,2) !=  0 || tmat_(5,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( -2  8  0  0 )\n"
                                     "(  0  0  8  0 )\n"
                                     "(  7  4 -3  8 )\n"
                                     "(  0 -8  5  0 )\n"
                                     "(  0  0  9 -6 )\n"
                                     "(  0  0  0 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major list assignment
   //=====================================================================================

   {
      test_ = "Column-major initializer list assignment (complete list)";;

      initialize();

      OBT band1 = blaze::band( tmat_, -1L );
      band1 = { 1, 2, 3, 4 };

      checkSize    ( band1,  4UL );
      checkCapacity( band1,  4UL );
      checkNonZeros( band1,  4UL );
      checkRows    ( tmat_,  6UL );
      checkColumns ( tmat_,  4UL );
      checkNonZeros( tmat_, 11UL );

      if( band1[0] != 1 || band1[1] != 2 || band1[2] != 3 || band1[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << band1 << "\n"
             << "   Expected result:\n( 1 2 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != -2 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  1 || tmat_(1,1) !=  0 || tmat_(1,2) !=  1 || tmat_(1,3) !=  0 ||
          tmat_(2,0) !=  7 || tmat_(2,1) !=  2 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 ||
          tmat_(3,0) !=  0 || tmat_(3,1) != -8 || tmat_(3,2) !=  3 || tmat_(3,3) !=  0 ||
          tmat_(4,0) !=  0 || tmat_(4,1) !=  0 || tmat_(4,2) !=  9 || tmat_(4,3) !=  4 ||
          tmat_(5,0) !=  0 || tmat_(5,1) !=  0 || tmat_(5,2) !=  0 || tmat_(5,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( -2  0  0  0 )\n"
                                     "(  1  0  1  0 )\n"
                                     "(  7  2 -3  0 )\n"
                                     "(  0 -8  3  0 )\n"
                                     "(  0  0  9  4 )\n"
                                     "(  0  0  0 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major initializer list assignment (incomplete list)";

      initialize();

      OBT band1 = blaze::band( tmat_, -1L );
      band1 = { 1, 2 };

      checkSize    ( band1, 4UL );
      checkCapacity( band1, 4UL );
      checkNonZeros( band1, 2UL );
      checkRows    ( tmat_, 6UL );
      checkColumns ( tmat_, 4UL );
      checkNonZeros( tmat_, 9UL );

      if( band1[0] != 1 || band1[1] != 2 || band1[2] != 0 || band1[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << band1 << "\n"
             << "   Expected result:\n( 1 2 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != -2 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  1 || tmat_(1,1) !=  0 || tmat_(1,2) !=  1 || tmat_(1,3) !=  0 ||
          tmat_(2,0) !=  7 || tmat_(2,1) !=  2 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 ||
          tmat_(3,0) !=  0 || tmat_(3,1) != -8 || tmat_(3,2) !=  0 || tmat_(3,3) !=  0 ||
          tmat_(4,0) !=  0 || tmat_(4,1) !=  0 || tmat_(4,2) !=  9 || tmat_(4,3) !=  0 ||
          tmat_(5,0) !=  0 || tmat_(5,1) !=  0 || tmat_(5,2) !=  0 || tmat_(5,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( -2  0  0  0 )\n"
                                     "(  1  0  1  0 )\n"
                                     "(  7  2 -3  0 )\n"
                                     "(  0 -8  0  0 )\n"
                                     "(  0  0  9  0 )\n"
                                     "(  0  0  0 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major copy assignment
   //=====================================================================================

   {
      test_ = "Column-major Band copy assignment";

      initialize();

      OBT band0 = blaze::band( tmat_, 0L );
      band0 = blaze::band( tmat_, -1L );

      checkSize    ( band0,  4UL );
      checkCapacity( band0,  4UL );
      checkNonZeros( band0,  3UL );
      checkRows    ( tmat_,  6UL );
      checkColumns ( tmat_,  4UL );
      checkNonZeros( tmat_, 11UL );

      if( band0[0] != 0 || band0[1] != 4 || band0[2] != 5 || band0[3] != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( 0  4  5 -6 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) != 0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) != 0 || tmat_(1,1) !=  4 || tmat_(1,2) != 1 || tmat_(1,3) !=  0 ||
          tmat_(2,0) != 7 || tmat_(2,1) !=  4 || tmat_(2,2) != 5 || tmat_(2,3) !=  0 ||
          tmat_(3,0) != 0 || tmat_(3,1) != -8 || tmat_(3,2) != 5 || tmat_(3,3) != -6 ||
          tmat_(4,0) != 0 || tmat_(4,1) !=  0 || tmat_(4,2) != 9 || tmat_(4,3) != -6 ||
          tmat_(5,0) != 0 || tmat_(5,1) !=  0 || tmat_(5,2) != 0 || tmat_(5,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0  0  0 )\n"
                                     "( 0  4  1  0 )\n"
                                     "( 7  4  5  0 )\n"
                                     "( 0 -8  5 -6 )\n"
                                     "( 0  0  9 -6 )\n"
                                     "( 0  0  0 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense vector assignment
   //=====================================================================================

   {
      test_ = "Column-major dense vector assignment (mixed type)";

      initialize();

      OBT band1 = blaze::band( tmat_, 1L );

      const blaze::DynamicVector<short,blaze::columnVector> vec1{ 8, 0, 9 };

      band1 = vec1;

      checkSize    ( band1,  3UL );
      checkCapacity( band1,  3UL );
      checkNonZeros( band1,  2UL );
      checkRows    ( tmat_,  6UL );
      checkColumns ( tmat_,  4UL );
      checkNonZeros( tmat_, 11UL );

      if( band1[0] != 8 || band1[1] != 0 || band1[2] != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << band1 << "\n"
             << "   Expected result:\n( 8 0 9 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != -2 || tmat_(0,1) !=  8 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) !=  0 || tmat_(1,2) !=  0 || tmat_(1,3) !=  0 ||
          tmat_(2,0) !=  7 || tmat_(2,1) !=  4 || tmat_(2,2) != -3 || tmat_(2,3) !=  9 ||
          tmat_(3,0) !=  0 || tmat_(3,1) != -8 || tmat_(3,2) !=  5 || tmat_(3,3) !=  0 ||
          tmat_(4,0) !=  0 || tmat_(4,1) !=  0 || tmat_(4,2) !=  9 || tmat_(4,3) != -6 ||
          tmat_(5,0) !=  0 || tmat_(5,1) !=  0 || tmat_(5,2) !=  0 || tmat_(5,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( -2  8  0  0 )\n"
                                     "(  0  0  0  0 )\n"
                                     "(  7  4 -3  9 )\n"
                                     "(  0 -8  5  0 )\n"
                                     "(  0  0  9 -6 )\n"
                                     "(  0  0  0 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major dense vector assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::columnVector;

      initialize();

      OBT band1 = blaze::band( tmat_, 1L );

      using AlignedPadded = blaze::CustomVector<int,aligned,padded,columnVector>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 16UL ) );
      AlignedPadded vec1( memory.get(), 3UL, 16UL );
      vec1[0] = 8;
      vec1[1] = 0;
      vec1[2] = 9;

      band1 = vec1;

      checkSize    ( band1,  3UL );
      checkCapacity( band1,  3UL );
      checkNonZeros( band1,  2UL );
      checkRows    ( tmat_,  6UL );
      checkColumns ( tmat_,  4UL );
      checkNonZeros( tmat_, 11UL );

      if( band1[0] != 8 || band1[1] != 0 || band1[2] != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << band1 << "\n"
             << "   Expected result:\n( 8 0 9 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != -2 || tmat_(0,1) !=  8 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) !=  0 || tmat_(1,2) !=  0 || tmat_(1,3) !=  0 ||
          tmat_(2,0) !=  7 || tmat_(2,1) !=  4 || tmat_(2,2) != -3 || tmat_(2,3) !=  9 ||
          tmat_(3,0) !=  0 || tmat_(3,1) != -8 || tmat_(3,2) !=  5 || tmat_(3,3) !=  0 ||
          tmat_(4,0) !=  0 || tmat_(4,1) !=  0 || tmat_(4,2) !=  9 || tmat_(4,3) != -6 ||
          tmat_(5,0) !=  0 || tmat_(5,1) !=  0 || tmat_(5,2) !=  0 || tmat_(5,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( -2  8  0  0 )\n"
                                     "(  0  0  0  0 )\n"
                                     "(  7  4 -3  9 )\n"
                                     "(  0 -8  5  0 )\n"
                                     "(  0  0  9 -6 )\n"
                                     "(  0  0  0 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major dense vector assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::columnVector;

      initialize();

      OBT band1 = blaze::band( tmat_, 1L );

      using UnalignedUnpadded = blaze::CustomVector<int,unaligned,unpadded,columnVector>;
      std::unique_ptr<int[]> memory( new int[4] );
      UnalignedUnpadded vec1( memory.get()+1UL, 3UL );
      vec1[0] = 8;
      vec1[1] = 0;
      vec1[2] = 9;

      band1 = vec1;

      checkSize    ( band1,  3UL );
      checkCapacity( band1,  3UL );
      checkNonZeros( band1,  2UL );
      checkRows    ( tmat_,  6UL );
      checkColumns ( tmat_,  4UL );
      checkNonZeros( tmat_, 11UL );

      if( band1[0] != 8 || band1[1] != 0 || band1[2] != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << band1 << "\n"
             << "   Expected result:\n( 8 0 9 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != -2 || tmat_(0,1) !=  8 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) !=  0 || tmat_(1,2) !=  0 || tmat_(1,3) !=  0 ||
          tmat_(2,0) !=  7 || tmat_(2,1) !=  4 || tmat_(2,2) != -3 || tmat_(2,3) !=  9 ||
          tmat_(3,0) !=  0 || tmat_(3,1) != -8 || tmat_(3,2) !=  5 || tmat_(3,3) !=  0 ||
          tmat_(4,0) !=  0 || tmat_(4,1) !=  0 || tmat_(4,2) !=  9 || tmat_(4,3) != -6 ||
          tmat_(5,0) !=  0 || tmat_(5,1) !=  0 || tmat_(5,2) !=  0 || tmat_(5,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( -2  8  0  0 )\n"
                                     "(  0  0  0  0 )\n"
                                     "(  7  4 -3  9 )\n"
                                     "(  0 -8  5  0 )\n"
                                     "(  0  0  9 -6 )\n"
                                     "(  0  0  0 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse vector assignment
   //=====================================================================================

   {
      test_ = "Column-major sparse vector assignment";

      initialize();

      OBT band2 = blaze::band( tmat_, -2L );

      blaze::CompressedVector<int,blaze::columnVector> vec1( 4UL );
      vec1[3] = 9;

      band2 = vec1;

      checkSize    ( band2, 4UL );
      checkCapacity( band2, 4UL );
      checkNonZeros( band2, 1UL );
      checkRows    ( tmat_, 6UL );
      checkColumns ( tmat_, 4UL );
      checkNonZeros( tmat_, 7UL );

      if( band2[0] != 0 || band2[1] != 0 || band2[2] != 0 || band2[3] != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << band2 << "\n"
             << "   Expected result:\n( 0 0 0 9 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != -2 || tmat_(0,1) != 0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) != 0 || tmat_(1,2) !=  1 || tmat_(1,3) !=  0 ||
          tmat_(2,0) !=  0 || tmat_(2,1) != 4 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 ||
          tmat_(3,0) !=  0 || tmat_(3,1) != 0 || tmat_(3,2) !=  5 || tmat_(3,3) !=  0 ||
          tmat_(4,0) !=  0 || tmat_(4,1) != 0 || tmat_(4,2) !=  0 || tmat_(4,3) != -6 ||
          tmat_(5,0) !=  0 || tmat_(5,1) != 0 || tmat_(5,2) !=  0 || tmat_(5,3) !=  9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( -2  0  0  0 )\n"
                                     "(  0  0  1  0 )\n"
                                     "(  0  4 -3  0 )\n"
                                     "(  0  0  5  0 )\n"
                                     "(  0  0  0 -6 )\n"
                                     "(  0  0  0  9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Band addition assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the addition assignment operators of the Band specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testAddAssign()
{
   //=====================================================================================
   // Row-major Band addition assignment
   //=====================================================================================

   {
      test_ = "Row-major Band addition assignment";

      initialize();

      BT band0 = blaze::band( mat_, 0L );
      band0 += blaze::band( mat_, 1L );

      checkSize    ( band0,  4UL );
      checkCapacity( band0,  4UL );
      checkNonZeros( band0,  4UL );
      checkRows    ( mat_ ,  4UL );
      checkColumns ( mat_ ,  6UL );
      checkNonZeros( mat_ , 12UL );

      if( band0[0] != -2 || band0[1] != 4 || band0[2] != 2 || band0[3] != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( -2 4 2 -6 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != -2 || mat_(0,1) != 0 || mat_(0,2) != 7 || mat_(0,3) !=  0 || mat_(0,4) !=  0 || mat_(0,5) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 4 || mat_(1,2) != 4 || mat_(1,3) != -8 || mat_(1,4) !=  0 || mat_(1,5) !=  0 ||
          mat_(2,0) !=  0 || mat_(2,1) != 1 || mat_(2,2) != 2 || mat_(2,3) !=  5 || mat_(2,4) !=  9 || mat_(2,5) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) != 0 || mat_(3,2) != 0 || mat_(3,3) != -6 || mat_(3,4) != -6 || mat_(3,5) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( -2  0  7  0  0  0 )\n"
                                     "(  0  4  4 -8  0  0 )\n"
                                     "(  0  1  2  5  9  0 )\n"
                                     "(  0  0  0 -6 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major dense vector addition assignment
   //=====================================================================================

   {
      test_ = "Row-major dense vector addition assignment (mixed type)";

      initialize();

      BT band0 = blaze::band( mat_, 0L );

      const blaze::DynamicVector<short,blaze::columnVector> vec{ 2, -4, 0, 0 };

      band0 += vec;

      checkSize    ( band0,  4UL );
      checkCapacity( band0,  4UL );
      checkNonZeros( band0,  2UL );
      checkRows    ( mat_ ,  4UL );
      checkColumns ( mat_ ,  6UL );
      checkNonZeros( mat_ , 10UL );

      if( band0[0] != 0 || band0[1] != -4 || band0[2] != -3 || band0[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( 0 -4 -3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) !=  0 || mat_(0,2) !=  7 || mat_(0,3) !=  0 || mat_(0,4) !=  0 || mat_(0,5) !=  0 ||
          mat_(1,0) != 0 || mat_(1,1) != -4 || mat_(1,2) !=  4 || mat_(1,3) != -8 || mat_(1,4) !=  0 || mat_(1,5) !=  0 ||
          mat_(2,0) != 0 || mat_(2,1) !=  1 || mat_(2,2) != -3 || mat_(2,3) !=  5 || mat_(2,4) !=  9 || mat_(2,5) !=  0 ||
          mat_(3,0) != 0 || mat_(3,1) !=  0 || mat_(3,2) !=  0 || mat_(3,3) !=  0 || mat_(3,4) != -6 || mat_(3,5) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0  7  0  0  0 )\n"
                                     "( 0 -4  4 -8  0  0 )\n"
                                     "( 0  1 -3  5  9  0 )\n"
                                     "( 0  0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major dense vector addition assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::columnVector;

      initialize();

      BT band0 = blaze::band( mat_, 0L );

      using AlignedPadded = blaze::CustomVector<int,aligned,padded,columnVector>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 16UL ) );
      AlignedPadded vec( memory.get(), 4UL, 16UL );
      vec[0] =  2;
      vec[1] = -4;
      vec[2] =  0;
      vec[3] =  0;

      band0 += vec;

      checkSize    ( band0,  4UL );
      checkCapacity( band0,  4UL );
      checkNonZeros( band0,  2UL );
      checkRows    ( mat_ ,  4UL );
      checkColumns ( mat_ ,  6UL );
      checkNonZeros( mat_ , 10UL );

      if( band0[0] != 0 || band0[1] != -4 || band0[2] != -3 || band0[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( 0 -4 -3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) !=  0 || mat_(0,2) !=  7 || mat_(0,3) !=  0 || mat_(0,4) !=  0 || mat_(0,5) !=  0 ||
          mat_(1,0) != 0 || mat_(1,1) != -4 || mat_(1,2) !=  4 || mat_(1,3) != -8 || mat_(1,4) !=  0 || mat_(1,5) !=  0 ||
          mat_(2,0) != 0 || mat_(2,1) !=  1 || mat_(2,2) != -3 || mat_(2,3) !=  5 || mat_(2,4) !=  9 || mat_(2,5) !=  0 ||
          mat_(3,0) != 0 || mat_(3,1) !=  0 || mat_(3,2) !=  0 || mat_(3,3) !=  0 || mat_(3,4) != -6 || mat_(3,5) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0  7  0  0  0 )\n"
                                     "( 0 -4  4 -8  0  0 )\n"
                                     "( 0  1 -3  5  9  0 )\n"
                                     "( 0  0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major dense vector addition assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::columnVector;

      initialize();

      BT band0 = blaze::band( mat_, 0L );

      using UnalignedUnpadded = blaze::CustomVector<int,unaligned,unpadded,columnVector>;
      std::unique_ptr<int[]> memory( new int[5] );
      UnalignedUnpadded vec( memory.get()+1UL, 4UL );
      vec[0] =  2;
      vec[1] = -4;
      vec[2] =  0;
      vec[3] =  0;

      band0 += vec;

      checkSize    ( band0,  4UL );
      checkCapacity( band0,  4UL );
      checkNonZeros( band0,  2UL );
      checkRows    ( mat_ ,  4UL );
      checkColumns ( mat_ ,  6UL );
      checkNonZeros( mat_ , 10UL );

      if( band0[0] != 0 || band0[1] != -4 || band0[2] != -3 || band0[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( 0 -4 -3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) !=  0 || mat_(0,2) !=  7 || mat_(0,3) !=  0 || mat_(0,4) !=  0 || mat_(0,5) !=  0 ||
          mat_(1,0) != 0 || mat_(1,1) != -4 || mat_(1,2) !=  4 || mat_(1,3) != -8 || mat_(1,4) !=  0 || mat_(1,5) !=  0 ||
          mat_(2,0) != 0 || mat_(2,1) !=  1 || mat_(2,2) != -3 || mat_(2,3) !=  5 || mat_(2,4) !=  9 || mat_(2,5) !=  0 ||
          mat_(3,0) != 0 || mat_(3,1) !=  0 || mat_(3,2) !=  0 || mat_(3,3) !=  0 || mat_(3,4) != -6 || mat_(3,5) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0  7  0  0  0 )\n"
                                     "( 0 -4  4 -8  0  0 )\n"
                                     "( 0  1 -3  5  9  0 )\n"
                                     "( 0  0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse vector addition assignment
   //=====================================================================================

   {
      test_ = "Row-major sparse vector addition assignment";

      initialize();

      BT band0 = blaze::band( mat_, 0L );

      blaze::CompressedVector<int,blaze::columnVector> vec( 4UL );
      vec[0] =  2;
      vec[1] = -4;

      band0 += vec;

      checkSize    ( band0,  4UL );
      checkCapacity( band0,  4UL );
      checkNonZeros( band0,  2UL );
      checkRows    ( mat_ ,  4UL );
      checkColumns ( mat_ ,  6UL );
      checkNonZeros( mat_ , 10UL );

      if( band0[0] != 0 || band0[1] != -4 || band0[2] != -3 || band0[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( 0 -4 -3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) !=  0 || mat_(0,2) !=  7 || mat_(0,3) !=  0 || mat_(0,4) !=  0 || mat_(0,5) !=  0 ||
          mat_(1,0) != 0 || mat_(1,1) != -4 || mat_(1,2) !=  4 || mat_(1,3) != -8 || mat_(1,4) !=  0 || mat_(1,5) !=  0 ||
          mat_(2,0) != 0 || mat_(2,1) !=  1 || mat_(2,2) != -3 || mat_(2,3) !=  5 || mat_(2,4) !=  9 || mat_(2,5) !=  0 ||
          mat_(3,0) != 0 || mat_(3,1) !=  0 || mat_(3,2) !=  0 || mat_(3,3) !=  0 || mat_(3,4) != -6 || mat_(3,5) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0  7  0  0  0 )\n"
                                     "( 0 -4  4 -8  0  0 )\n"
                                     "( 0  1 -3  5  9  0 )\n"
                                     "( 0  0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major Band addition assignment
   //=====================================================================================

   {
      test_ = "Column-major Band addition assignment";

      initialize();

      OBT band0 = blaze::band( tmat_, 0L );
      band0 += blaze::band( tmat_, -1L );

      checkSize    ( band0,  4UL );
      checkCapacity( band0,  4UL );
      checkNonZeros( band0,  4UL );
      checkRows    ( tmat_,  6UL );
      checkColumns ( tmat_,  4UL );
      checkNonZeros( tmat_, 12UL );

      if( band0[0] != -2 || band0[1] != 4 || band0[2] != 2 || band0[3] != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( -2 4 2 -6 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != -2 || tmat_(0,1) !=  0 || tmat_(0,2) != 0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) !=  4 || tmat_(1,2) != 1 || tmat_(1,3) !=  0 ||
          tmat_(2,0) !=  7 || tmat_(2,1) !=  4 || tmat_(2,2) != 2 || tmat_(2,3) !=  0 ||
          tmat_(3,0) !=  0 || tmat_(3,1) != -8 || tmat_(3,2) != 5 || tmat_(3,3) != -6 ||
          tmat_(4,0) !=  0 || tmat_(4,1) !=  0 || tmat_(4,2) != 9 || tmat_(4,3) != -6 ||
          tmat_(5,0) !=  0 || tmat_(5,1) !=  0 || tmat_(5,2) != 0 || tmat_(5,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( -2  0  0  0 )\n"
                                     "(  0  4  1  0 )\n"
                                     "(  7  4  2  0 )\n"
                                     "(  0 -8  5 -6 )\n"
                                     "(  0  0  9 -6 )\n"
                                     "(  0  0  0 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense vector addition assignment
   //=====================================================================================

   {
      test_ = "Column-major dense vector addition assignment (mixed type)";

      initialize();

      OBT band0 = blaze::band( tmat_, 0L );

      const blaze::DynamicVector<short,blaze::columnVector> vec{ 2, -4, 0, 0 };

      band0 += vec;

      checkSize    ( band0,  4UL );
      checkCapacity( band0,  4UL );
      checkNonZeros( band0,  2UL );
      checkRows    ( tmat_,  6UL );
      checkColumns ( tmat_,  4UL );
      checkNonZeros( tmat_, 10UL );

      if( band0[0] != 0 || band0[1] != -4 || band0[2] != -3 || band0[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( 0 -4 -3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) != 0 || tmat_(1,1) != -4 || tmat_(1,2) !=  1 || tmat_(1,3) !=  0 ||
          tmat_(2,0) != 7 || tmat_(2,1) !=  4 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 ||
          tmat_(3,0) != 0 || tmat_(3,1) != -8 || tmat_(3,2) !=  5 || tmat_(3,3) !=  0 ||
          tmat_(4,0) != 0 || tmat_(4,1) !=  0 || tmat_(4,2) !=  9 || tmat_(4,3) != -6 ||
          tmat_(5,0) != 0 || tmat_(5,1) !=  0 || tmat_(5,2) !=  0 || tmat_(5,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0  0  0 )\n"
                                     "( 0 -4  1  0 )\n"
                                     "( 7  4 -3  0 )\n"
                                     "( 0 -8  5  0 )\n"
                                     "( 0  0  9 -6 )\n"
                                     "( 0  0  0 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major dense vector addition assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::columnVector;

      initialize();

      OBT band0 = blaze::band( tmat_, 0L );

      using AlignedPadded = blaze::CustomVector<int,aligned,padded,columnVector>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 16UL ) );
      AlignedPadded vec( memory.get(), 4UL, 16UL );
      vec[0] =  2;
      vec[1] = -4;
      vec[2] =  0;
      vec[3] =  0;

      band0 += vec;

      checkSize    ( band0,  4UL );
      checkCapacity( band0,  4UL );
      checkNonZeros( band0,  2UL );
      checkRows    ( tmat_,  6UL );
      checkColumns ( tmat_,  4UL );
      checkNonZeros( tmat_, 10UL );

      if( band0[0] != 0 || band0[1] != -4 || band0[2] != -3 || band0[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( 0 -4 -3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) != 0 || tmat_(1,1) != -4 || tmat_(1,2) !=  1 || tmat_(1,3) !=  0 ||
          tmat_(2,0) != 7 || tmat_(2,1) !=  4 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 ||
          tmat_(3,0) != 0 || tmat_(3,1) != -8 || tmat_(3,2) !=  5 || tmat_(3,3) !=  0 ||
          tmat_(4,0) != 0 || tmat_(4,1) !=  0 || tmat_(4,2) !=  9 || tmat_(4,3) != -6 ||
          tmat_(5,0) != 0 || tmat_(5,1) !=  0 || tmat_(5,2) !=  0 || tmat_(5,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0  0  0 )\n"
                                     "( 0 -4  1  0 )\n"
                                     "( 7  4 -3  0 )\n"
                                     "( 0 -8  5  0 )\n"
                                     "( 0  0  9 -6 )\n"
                                     "( 0  0  0 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major dense vector addition assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::columnVector;

      initialize();

      OBT band0 = blaze::band( tmat_, 0L );

      using UnalignedUnpadded = blaze::CustomVector<int,unaligned,unpadded,columnVector>;
      std::unique_ptr<int[]> memory( new int[5] );
      UnalignedUnpadded vec( memory.get()+1UL, 4UL );
      vec[0] =  2;
      vec[1] = -4;
      vec[2] =  0;
      vec[3] =  0;

      band0 += vec;

      checkSize    ( band0,  4UL );
      checkCapacity( band0,  4UL );
      checkNonZeros( band0,  2UL );
      checkRows    ( tmat_,  6UL );
      checkColumns ( tmat_,  4UL );
      checkNonZeros( tmat_, 10UL );

      if( band0[0] != 0 || band0[1] != -4 || band0[2] != -3 || band0[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( 0 -4 -3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) != 0 || tmat_(1,1) != -4 || tmat_(1,2) !=  1 || tmat_(1,3) !=  0 ||
          tmat_(2,0) != 7 || tmat_(2,1) !=  4 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 ||
          tmat_(3,0) != 0 || tmat_(3,1) != -8 || tmat_(3,2) !=  5 || tmat_(3,3) !=  0 ||
          tmat_(4,0) != 0 || tmat_(4,1) !=  0 || tmat_(4,2) !=  9 || tmat_(4,3) != -6 ||
          tmat_(5,0) != 0 || tmat_(5,1) !=  0 || tmat_(5,2) !=  0 || tmat_(5,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0  0  0 )\n"
                                     "( 0 -4  1  0 )\n"
                                     "( 7  4 -3  0 )\n"
                                     "( 0 -8  5  0 )\n"
                                     "( 0  0  9 -6 )\n"
                                     "( 0  0  0 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse vector addition assignment
   //=====================================================================================

   {
      test_ = "Column-major sparse vector addition assignment";

      initialize();

      OBT band0 = blaze::band( tmat_, 0L );

      blaze::CompressedVector<int,blaze::columnVector> vec( 4UL );
      vec[0] =  2;
      vec[1] = -4;

      band0 += vec;

      checkSize    ( band0,  4UL );
      checkCapacity( band0,  4UL );
      checkNonZeros( band0,  2UL );
      checkRows    ( tmat_,  6UL );
      checkColumns ( tmat_,  4UL );
      checkNonZeros( tmat_, 10UL );

      if( band0[0] != 0 || band0[1] != -4 || band0[2] != -3 || band0[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( 0 -4 -3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) != 0 || tmat_(1,1) != -4 || tmat_(1,2) !=  1 || tmat_(1,3) !=  0 ||
          tmat_(2,0) != 7 || tmat_(2,1) !=  4 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 ||
          tmat_(3,0) != 0 || tmat_(3,1) != -8 || tmat_(3,2) !=  5 || tmat_(3,3) !=  0 ||
          tmat_(4,0) != 0 || tmat_(4,1) !=  0 || tmat_(4,2) !=  9 || tmat_(4,3) != -6 ||
          tmat_(5,0) != 0 || tmat_(5,1) !=  0 || tmat_(5,2) !=  0 || tmat_(5,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0  0  0 )\n"
                                     "( 0 -4  1  0 )\n"
                                     "( 7  4 -3  0 )\n"
                                     "( 0 -8  5  0 )\n"
                                     "( 0  0  9 -6 )\n"
                                     "( 0  0  0 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Band subtraction assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the subtraction assignment operators of the Band
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testSubAssign()
{
   //=====================================================================================
   // Row-major Band subtraction assignment
   //=====================================================================================

   {
      test_ = "Row-major Band subtraction assignment";

      initialize();

      BT band0 = blaze::band( mat_, 0L );
      band0 -= blaze::band( mat_, 1L );

      checkSize    ( band0,  4UL );
      checkCapacity( band0,  4UL );
      checkNonZeros( band0,  4UL );
      checkRows    ( mat_ ,  4UL );
      checkColumns ( mat_ ,  6UL );
      checkNonZeros( mat_ , 12UL );

      if( band0[0] != -2 || band0[1] != -4 || band0[2] != -8 || band0[3] != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( -2 -4 -8 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != -2 || mat_(0,1) !=  0 || mat_(0,2) !=  7 || mat_(0,3) !=  0 || mat_(0,4) !=  0 || mat_(0,5) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != -4 || mat_(1,2) !=  4 || mat_(1,3) != -8 || mat_(1,4) !=  0 || mat_(1,5) !=  0 ||
          mat_(2,0) !=  0 || mat_(2,1) !=  1 || mat_(2,2) != -8 || mat_(2,3) !=  5 || mat_(2,4) !=  9 || mat_(2,5) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) !=  0 || mat_(3,2) !=  0 || mat_(3,3) !=  6 || mat_(3,4) != -6 || mat_(3,5) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( -2  0  7  0  0  0 )\n"
                                     "(  0 -4  4 -8  0  0 )\n"
                                     "(  0  1 -8  5  9  0 )\n"
                                     "(  0  0  0  6 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major dense vector subtraction assignment
   //=====================================================================================

   {
      test_ = "Row-major dense vector subtraction assignment (mixed type)";

      initialize();

      BT band0 = blaze::band( mat_, 0L );

      const blaze::DynamicVector<short,blaze::columnVector> vec{ 2, -4, 0, 0 };

      band0 -= vec;

      checkSize    ( band0,  4UL );
      checkCapacity( band0,  4UL );
      checkNonZeros( band0,  3UL );
      checkRows    ( mat_ ,  4UL );
      checkColumns ( mat_ ,  6UL );
      checkNonZeros( mat_ , 11UL );

      if( band0[0] != -4 || band0[1] != 4 || band0[2] != -3 || band0[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( -4 4 -3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != -4 || mat_(0,1) != 0 || mat_(0,2) !=  7 || mat_(0,3) !=  0 || mat_(0,4) !=  0 || mat_(0,5) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 4 || mat_(1,2) !=  4 || mat_(1,3) != -8 || mat_(1,4) !=  0 || mat_(1,5) !=  0 ||
          mat_(2,0) !=  0 || mat_(2,1) != 1 || mat_(2,2) != -3 || mat_(2,3) !=  5 || mat_(2,4) !=  9 || mat_(2,5) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) !=  0 || mat_(3,4) != -6 || mat_(3,5) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( -4  0  7  0  0  0 )\n"
                                     "(  0  4  4 -8  0  0 )\n"
                                     "(  0  1 -3  5  9  0 )\n"
                                     "(  0  0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major dense vector subtraction assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::columnVector;

      initialize();

      BT band0 = blaze::band( mat_, 0L );

      using AlignedPadded = blaze::CustomVector<int,aligned,padded,columnVector>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 16UL ) );
      AlignedPadded vec( memory.get(), 4UL, 16UL );
      vec[0] =  2;
      vec[1] = -4;
      vec[2] =  0;
      vec[3] =  0;

      band0 -= vec;

      checkSize    ( band0,  4UL );
      checkCapacity( band0,  4UL );
      checkNonZeros( band0,  3UL );
      checkRows    ( mat_ ,  4UL );
      checkColumns ( mat_ ,  6UL );
      checkNonZeros( mat_ , 11UL );

      if( band0[0] != -4 || band0[1] != 4 || band0[2] != -3 || band0[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( -4 4 -3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != -4 || mat_(0,1) != 0 || mat_(0,2) !=  7 || mat_(0,3) !=  0 || mat_(0,4) !=  0 || mat_(0,5) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 4 || mat_(1,2) !=  4 || mat_(1,3) != -8 || mat_(1,4) !=  0 || mat_(1,5) !=  0 ||
          mat_(2,0) !=  0 || mat_(2,1) != 1 || mat_(2,2) != -3 || mat_(2,3) !=  5 || mat_(2,4) !=  9 || mat_(2,5) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) !=  0 || mat_(3,4) != -6 || mat_(3,5) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( -4  0  7  0  0  0 )\n"
                                     "(  0  4  4 -8  0  0 )\n"
                                     "(  0  1 -3  5  9  0 )\n"
                                     "(  0  0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major dense vector subtraction assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::columnVector;

      initialize();

      BT band0 = blaze::band( mat_, 0L );

      using UnalignedUnpadded = blaze::CustomVector<int,unaligned,unpadded,columnVector>;
      std::unique_ptr<int[]> memory( new int[5] );
      UnalignedUnpadded vec( memory.get()+1UL, 4UL );
      vec[0] =  2;
      vec[1] = -4;
      vec[2] =  0;
      vec[3] =  0;

      band0 -= vec;

      checkSize    ( band0,  4UL );
      checkCapacity( band0,  4UL );
      checkNonZeros( band0,  3UL );
      checkRows    ( mat_ ,  4UL );
      checkColumns ( mat_ ,  6UL );
      checkNonZeros( mat_ , 11UL );

      if( band0[0] != -4 || band0[1] != 4 || band0[2] != -3 || band0[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( -4 4 -3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != -4 || mat_(0,1) != 0 || mat_(0,2) !=  7 || mat_(0,3) !=  0 || mat_(0,4) !=  0 || mat_(0,5) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 4 || mat_(1,2) !=  4 || mat_(1,3) != -8 || mat_(1,4) !=  0 || mat_(1,5) !=  0 ||
          mat_(2,0) !=  0 || mat_(2,1) != 1 || mat_(2,2) != -3 || mat_(2,3) !=  5 || mat_(2,4) !=  9 || mat_(2,5) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) !=  0 || mat_(3,4) != -6 || mat_(3,5) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( -4  0  7  0  0  0 )\n"
                                     "(  0  4  4 -8  0  0 )\n"
                                     "(  0  1 -3  5  9  0 )\n"
                                     "(  0  0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse vector subtraction assignment
   //=====================================================================================

   {
      test_ = "Row-major sparse vector subtraction assignment";

      initialize();

      BT band0 = blaze::band( mat_, 0L );

      blaze::CompressedVector<int,blaze::columnVector> vec( 4UL );
      vec[0] =  2;
      vec[1] = -4;

      band0 -= vec;

      checkSize    ( band0,  4UL );
      checkCapacity( band0,  4UL );
      checkNonZeros( band0,  3UL );
      checkRows    ( mat_ ,  4UL );
      checkColumns ( mat_ ,  6UL );
      checkNonZeros( mat_ , 11UL );

      if( band0[0] != -4 || band0[1] != 4 || band0[2] != -3 || band0[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( -4 4 -3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != -4 || mat_(0,1) != 0 || mat_(0,2) !=  7 || mat_(0,3) !=  0 || mat_(0,4) !=  0 || mat_(0,5) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 4 || mat_(1,2) !=  4 || mat_(1,3) != -8 || mat_(1,4) !=  0 || mat_(1,5) !=  0 ||
          mat_(2,0) !=  0 || mat_(2,1) != 1 || mat_(2,2) != -3 || mat_(2,3) !=  5 || mat_(2,4) !=  9 || mat_(2,5) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) !=  0 || mat_(3,4) != -6 || mat_(3,5) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( -4  0  7  0  0  0 )\n"
                                     "(  0  4  4 -8  0  0 )\n"
                                     "(  0  1 -3  5  9  0 )\n"
                                     "(  0  0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major Band subtraction assignment
   //=====================================================================================

   {
      test_ = "Column-major Band subtraction assignment";

      initialize();

      OBT band0 = blaze::band( tmat_, 0L );
      band0 -= blaze::band( tmat_, -1L );

      checkSize    ( band0,  4UL );
      checkCapacity( band0,  4UL );
      checkNonZeros( band0,  4UL );
      checkRows    ( tmat_,  6UL );
      checkColumns ( tmat_,  4UL );
      checkNonZeros( tmat_, 12UL );

      if( band0[0] != -2 || band0[1] != -4 || band0[2] != -8 || band0[3] != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( -2 -4 -8 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != -2 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) != -4 || tmat_(1,2) !=  1 || tmat_(1,3) !=  0 ||
          tmat_(2,0) !=  7 || tmat_(2,1) !=  4 || tmat_(2,2) != -8 || tmat_(2,3) !=  0 ||
          tmat_(3,0) !=  0 || tmat_(3,1) != -8 || tmat_(3,2) !=  5 || tmat_(3,3) !=  6 ||
          tmat_(4,0) !=  0 || tmat_(4,1) !=  0 || tmat_(4,2) !=  9 || tmat_(4,3) != -6 ||
          tmat_(5,0) !=  0 || tmat_(5,1) !=  0 || tmat_(5,2) !=  0 || tmat_(5,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( -2  0  0  0 )\n"
                                     "(  0 -4  1  0 )\n"
                                     "(  7  4 -8  0 )\n"
                                     "(  0 -8  5  6 )\n"
                                     "(  0  0  9 -6 )\n"
                                     "(  0  0  0 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense vector subtraction assignment
   //=====================================================================================

   {
      test_ = "Column-major dense vector subtraction assignment (mixed type)";

      initialize();

      OBT band0 = blaze::band( tmat_, 0L );

      const blaze::DynamicVector<short,blaze::columnVector> vec{ 2, -4, 0, 0 };

      band0 -= vec;

      checkSize    ( band0,  4UL );
      checkCapacity( band0,  4UL );
      checkNonZeros( band0,  3UL );
      checkRows    ( tmat_,  6UL );
      checkColumns ( tmat_,  4UL );
      checkNonZeros( tmat_, 11UL );

      if( band0[0] != -4 || band0[1] != 4 || band0[2] != -3 || band0[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( -4 4 -3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != -4 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) !=  4 || tmat_(1,2) !=  1 || tmat_(1,3) !=  0 ||
          tmat_(2,0) !=  7 || tmat_(2,1) !=  4 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 ||
          tmat_(3,0) !=  0 || tmat_(3,1) != -8 || tmat_(3,2) !=  5 || tmat_(3,3) !=  0 ||
          tmat_(4,0) !=  0 || tmat_(4,1) !=  0 || tmat_(4,2) !=  9 || tmat_(4,3) != -6 ||
          tmat_(5,0) !=  0 || tmat_(5,1) !=  0 || tmat_(5,2) !=  0 || tmat_(5,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( -4  0  0  0 )\n"
                                     "(  0  4  1  0 )\n"
                                     "(  7  4 -3  0 )\n"
                                     "(  0 -8  5  0 )\n"
                                     "(  0  0  9 -6 )\n"
                                     "(  0  0  0 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major dense vector subtraction assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::columnVector;

      initialize();

      OBT band0 = blaze::band( tmat_, 0L );

      using AlignedPadded = blaze::CustomVector<int,aligned,padded,columnVector>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 16UL ) );
      AlignedPadded vec( memory.get(), 4UL, 16UL );
      vec[0] =  2;
      vec[1] = -4;
      vec[2] =  0;
      vec[3] =  0;

      band0 -= vec;

      checkSize    ( band0,  4UL );
      checkCapacity( band0,  4UL );
      checkNonZeros( band0,  3UL );
      checkRows    ( tmat_,  6UL );
      checkColumns ( tmat_,  4UL );
      checkNonZeros( tmat_, 11UL );

      if( band0[0] != -4 || band0[1] != 4 || band0[2] != -3 || band0[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( -4 4 -3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != -4 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) !=  4 || tmat_(1,2) !=  1 || tmat_(1,3) !=  0 ||
          tmat_(2,0) !=  7 || tmat_(2,1) !=  4 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 ||
          tmat_(3,0) !=  0 || tmat_(3,1) != -8 || tmat_(3,2) !=  5 || tmat_(3,3) !=  0 ||
          tmat_(4,0) !=  0 || tmat_(4,1) !=  0 || tmat_(4,2) !=  9 || tmat_(4,3) != -6 ||
          tmat_(5,0) !=  0 || tmat_(5,1) !=  0 || tmat_(5,2) !=  0 || tmat_(5,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( -4  0  0  0 )\n"
                                     "(  0  4  1  0 )\n"
                                     "(  7  4 -3  0 )\n"
                                     "(  0 -8  5  0 )\n"
                                     "(  0  0  9 -6 )\n"
                                     "(  0  0  0 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major dense vector subtraction assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::columnVector;

      initialize();

      OBT band0 = blaze::band( tmat_, 0L );

      using UnalignedUnpadded = blaze::CustomVector<int,unaligned,unpadded,columnVector>;
      std::unique_ptr<int[]> memory( new int[5] );
      UnalignedUnpadded vec( memory.get()+1UL, 4UL );
      vec[0] =  2;
      vec[1] = -4;
      vec[2] =  0;
      vec[3] =  0;

      band0 -= vec;

      checkSize    ( band0,  4UL );
      checkCapacity( band0,  4UL );
      checkNonZeros( band0,  3UL );
      checkRows    ( tmat_,  6UL );
      checkColumns ( tmat_,  4UL );
      checkNonZeros( tmat_, 11UL );

      if( band0[0] != -4 || band0[1] != 4 || band0[2] != -3 || band0[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( -4 4 -3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != -4 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) !=  4 || tmat_(1,2) !=  1 || tmat_(1,3) !=  0 ||
          tmat_(2,0) !=  7 || tmat_(2,1) !=  4 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 ||
          tmat_(3,0) !=  0 || tmat_(3,1) != -8 || tmat_(3,2) !=  5 || tmat_(3,3) !=  0 ||
          tmat_(4,0) !=  0 || tmat_(4,1) !=  0 || tmat_(4,2) !=  9 || tmat_(4,3) != -6 ||
          tmat_(5,0) !=  0 || tmat_(5,1) !=  0 || tmat_(5,2) !=  0 || tmat_(5,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( -4  0  0  0 )\n"
                                     "(  0  4  1  0 )\n"
                                     "(  7  4 -3  0 )\n"
                                     "(  0 -8  5  0 )\n"
                                     "(  0  0  9 -6 )\n"
                                     "(  0  0  0 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse vector subtraction assignment
   //=====================================================================================

   {
      test_ = "Column-major sparse vector subtraction assignment";

      initialize();

      OBT band0 = blaze::band( tmat_, 0L );

      blaze::CompressedVector<int,blaze::columnVector> vec( 4UL );
      vec[0] =  2;
      vec[1] = -4;

      band0 -= vec;

      checkSize    ( band0,  4UL );
      checkCapacity( band0,  4UL );
      checkNonZeros( band0,  3UL );
      checkRows    ( tmat_,  6UL );
      checkColumns ( tmat_,  4UL );
      checkNonZeros( tmat_, 11UL );

      if( band0[0] != -4 || band0[1] != 4 || band0[2] != -3 || band0[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( -4 4 -3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != -4 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) !=  4 || tmat_(1,2) !=  1 || tmat_(1,3) !=  0 ||
          tmat_(2,0) !=  7 || tmat_(2,1) !=  4 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 ||
          tmat_(3,0) !=  0 || tmat_(3,1) != -8 || tmat_(3,2) !=  5 || tmat_(3,3) !=  0 ||
          tmat_(4,0) !=  0 || tmat_(4,1) !=  0 || tmat_(4,2) !=  9 || tmat_(4,3) != -6 ||
          tmat_(5,0) !=  0 || tmat_(5,1) !=  0 || tmat_(5,2) !=  0 || tmat_(5,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( -4  0  0  0 )\n"
                                     "(  0  4  1  0 )\n"
                                     "(  7  4 -3  0 )\n"
                                     "(  0 -8  5  0 )\n"
                                     "(  0  0  9 -6 )\n"
                                     "(  0  0  0 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Band multiplication assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the multiplication assignment operators of the Band
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testMultAssign()
{
   //=====================================================================================
   // Row-major Band multiplication assignment
   //=====================================================================================

   {
      test_ = "Row-major Band multiplication assignment";

      initialize();

      BT band0 = blaze::band( mat_, 0L );
      band0 *= blaze::band( mat_, 1L );

      checkSize    ( band0, 4UL );
      checkCapacity( band0, 4UL );
      checkNonZeros( band0, 1UL );
      checkRows    ( mat_ , 4UL );
      checkColumns ( mat_ , 6UL );
      checkNonZeros( mat_ , 9UL );

      if( band0[0] != 0 || band0[1] != 0 || band0[2] != -15 || band0[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( 0 0 -15 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) !=   7 || mat_(0,3) !=  0 || mat_(0,4) !=  0 || mat_(0,5) !=  0 ||
          mat_(1,0) != 0 || mat_(1,1) != 0 || mat_(1,2) !=   4 || mat_(1,3) != -8 || mat_(1,4) !=  0 || mat_(1,5) !=  0 ||
          mat_(2,0) != 0 || mat_(2,1) != 1 || mat_(2,2) != -15 || mat_(2,3) !=  5 || mat_(2,4) !=  9 || mat_(2,5) !=  0 ||
          mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) !=   0 || mat_(3,3) !=  0 || mat_(3,4) != -6 || mat_(3,5) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0   7  0  0  0 )\n"
                                     "( 0  0   4 -8  0  0 )\n"
                                     "( 0  1 -15  5  9  0 )\n"
                                     "( 0  0   0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major dense vector multiplication assignment
   //=====================================================================================

   {
      test_ = "Row-major dense vector multiplication assignment (mixed type)";

      initialize();

      BT band0 = blaze::band( mat_, 0L );

      const blaze::DynamicVector<short,blaze::columnVector> vec{ 2, -4, 0, 0 };

      band0 *= vec;

      checkSize    ( band0, 4UL );
      checkCapacity( band0, 4UL );
      checkNonZeros( band0, 1UL );
      checkRows    ( mat_ , 4UL );
      checkColumns ( mat_ , 6UL );
      checkNonZeros( mat_ , 9UL );

      if( band0[0] != -4 || band0[1] != 0 || band0[2] != 0 || band0[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( -4 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != -4 || mat_(0,1) != 0 || mat_(0,2) != 7 || mat_(0,3) !=  0 || mat_(0,4) !=  0 || mat_(0,5) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 0 || mat_(1,2) != 4 || mat_(1,3) != -8 || mat_(1,4) !=  0 || mat_(1,5) !=  0 ||
          mat_(2,0) !=  0 || mat_(2,1) != 1 || mat_(2,2) != 0 || mat_(2,3) !=  5 || mat_(2,4) !=  9 || mat_(2,5) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) != 0 || mat_(3,2) != 0 || mat_(3,3) !=  0 || mat_(3,4) != -6 || mat_(3,5) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( -4  0  7  0  0  0 )\n"
                                     "(  0  0  4 -8  0  0 )\n"
                                     "(  0  1  0  5  9  0 )\n"
                                     "(  0  0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major dense vector multiplication assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::columnVector;

      initialize();

      BT band0 = blaze::band( mat_, 0L );

      using AlignedPadded = blaze::CustomVector<int,aligned,padded,columnVector>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 16UL ) );
      AlignedPadded vec( memory.get(), 4UL, 16UL );
      vec[0] =  2;
      vec[1] = -4;
      vec[2] =  0;
      vec[3] =  0;

      band0 *= vec;

      checkSize    ( band0, 4UL );
      checkCapacity( band0, 4UL );
      checkNonZeros( band0, 1UL );
      checkRows    ( mat_ , 4UL );
      checkColumns ( mat_ , 6UL );
      checkNonZeros( mat_ , 9UL );

      if( band0[0] != -4 || band0[1] != 0 || band0[2] != 0 || band0[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( -4 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != -4 || mat_(0,1) != 0 || mat_(0,2) != 7 || mat_(0,3) !=  0 || mat_(0,4) !=  0 || mat_(0,5) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 0 || mat_(1,2) != 4 || mat_(1,3) != -8 || mat_(1,4) !=  0 || mat_(1,5) !=  0 ||
          mat_(2,0) !=  0 || mat_(2,1) != 1 || mat_(2,2) != 0 || mat_(2,3) !=  5 || mat_(2,4) !=  9 || mat_(2,5) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) != 0 || mat_(3,2) != 0 || mat_(3,3) !=  0 || mat_(3,4) != -6 || mat_(3,5) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( -4  0  7  0  0  0 )\n"
                                     "(  0  0  4 -8  0  0 )\n"
                                     "(  0  1  0  5  9  0 )\n"
                                     "(  0  0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major dense vector multiplication assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::columnVector;

      initialize();

      BT band0 = blaze::band( mat_, 0L );

      using UnalignedUnpadded = blaze::CustomVector<int,unaligned,unpadded,columnVector>;
      std::unique_ptr<int[]> memory( new int[5] );
      UnalignedUnpadded vec( memory.get()+1UL, 4UL );
      vec[0] =  2;
      vec[1] = -4;
      vec[2] =  0;
      vec[3] =  0;

      band0 *= vec;

      checkSize    ( band0, 4UL );
      checkCapacity( band0, 4UL );
      checkNonZeros( band0, 1UL );
      checkRows    ( mat_ , 4UL );
      checkColumns ( mat_ , 6UL );
      checkNonZeros( mat_ , 9UL );

      if( band0[0] != -4 || band0[1] != 0 || band0[2] != 0 || band0[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( -4 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != -4 || mat_(0,1) != 0 || mat_(0,2) != 7 || mat_(0,3) !=  0 || mat_(0,4) !=  0 || mat_(0,5) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 0 || mat_(1,2) != 4 || mat_(1,3) != -8 || mat_(1,4) !=  0 || mat_(1,5) !=  0 ||
          mat_(2,0) !=  0 || mat_(2,1) != 1 || mat_(2,2) != 0 || mat_(2,3) !=  5 || mat_(2,4) !=  9 || mat_(2,5) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) != 0 || mat_(3,2) != 0 || mat_(3,3) !=  0 || mat_(3,4) != -6 || mat_(3,5) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( -4  0  7  0  0  0 )\n"
                                     "(  0  0  4 -8  0  0 )\n"
                                     "(  0  1  0  5  9  0 )\n"
                                     "(  0  0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse vector multiplication assignment
   //=====================================================================================

   {
      test_ = "Row-major sparse vector multiplication assignment";

      initialize();

      BT band0 = blaze::band( mat_, 0L );

      blaze::CompressedVector<int,blaze::columnVector> vec( 4UL );
      vec[0] =  2;
      vec[1] = -4;

      band0 *= vec;

      checkSize    ( band0, 4UL );
      checkCapacity( band0, 4UL );
      checkNonZeros( band0, 1UL );
      checkRows    ( mat_ , 4UL );
      checkColumns ( mat_ , 6UL );
      checkNonZeros( mat_ , 9UL );

      if( band0[0] != -4 || band0[1] != 0 || band0[2] != 0 || band0[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( -4 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != -4 || mat_(0,1) != 0 || mat_(0,2) != 7 || mat_(0,3) !=  0 || mat_(0,4) !=  0 || mat_(0,5) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 0 || mat_(1,2) != 4 || mat_(1,3) != -8 || mat_(1,4) !=  0 || mat_(1,5) !=  0 ||
          mat_(2,0) !=  0 || mat_(2,1) != 1 || mat_(2,2) != 0 || mat_(2,3) !=  5 || mat_(2,4) !=  9 || mat_(2,5) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) != 0 || mat_(3,2) != 0 || mat_(3,3) !=  0 || mat_(3,4) != -6 || mat_(3,5) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( -4  0  7  0  0  0 )\n"
                                     "(  0  0  4 -8  0  0 )\n"
                                     "(  0  1  0  5  9  0 )\n"
                                     "(  0  0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major Band multiplication assignment
   //=====================================================================================

   {
      test_ = "Column-major Band multiplication assignment";

      initialize();

      OBT band0 = blaze::band( tmat_, 0L );
      band0 *= blaze::band( tmat_, -1L );

      checkSize    ( band0, 4UL );
      checkCapacity( band0, 4UL );
      checkNonZeros( band0, 1UL );
      checkRows    ( tmat_, 6UL );
      checkColumns ( tmat_, 4UL );
      checkNonZeros( tmat_, 9UL );

      if( band0[0] != 0 || band0[1] != 0 || band0[2] != -15 || band0[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( 0 0 -15 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) !=   0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) != 0 || tmat_(1,1) !=  0 || tmat_(1,2) !=   1 || tmat_(1,3) !=  0 ||
          tmat_(2,0) != 7 || tmat_(2,1) !=  4 || tmat_(2,2) != -15 || tmat_(2,3) !=  0 ||
          tmat_(3,0) != 0 || tmat_(3,1) != -8 || tmat_(3,2) !=   5 || tmat_(3,3) !=  0 ||
          tmat_(4,0) != 0 || tmat_(4,1) !=  0 || tmat_(4,2) !=   9 || tmat_(4,3) != -6 ||
          tmat_(5,0) != 0 || tmat_(5,1) !=  0 || tmat_(5,2) !=   0 || tmat_(5,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0   0  0 )\n"
                                     "( 0  0   1  0 )\n"
                                     "( 7  4 -15  0 )\n"
                                     "( 0 -8   5  0 )\n"
                                     "( 0  0   9 -6 )\n"
                                     "( 0  0   0 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense vector multiplication assignment
   //=====================================================================================

   {
      test_ = "Column-major dense vector multiplication assignment (mixed type)";

      initialize();

      OBT band0 = blaze::band( tmat_, 0L );

      const blaze::DynamicVector<short,blaze::columnVector> vec{ 2, -4, 0, 0 };

      band0 *= vec;

      checkSize    ( band0, 4UL );
      checkCapacity( band0, 4UL );
      checkNonZeros( band0, 1UL );
      checkRows    ( tmat_, 6UL );
      checkColumns ( tmat_, 4UL );
      checkNonZeros( tmat_, 9UL );

      if( band0[0] != -4 || band0[1] != 0 || band0[2] != 0 || band0[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( -4 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != -4 || tmat_(0,1) !=  0 || tmat_(0,2) != 0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) !=  0 || tmat_(1,2) != 1 || tmat_(1,3) !=  0 ||
          tmat_(2,0) !=  7 || tmat_(2,1) !=  4 || tmat_(2,2) != 0 || tmat_(2,3) !=  0 ||
          tmat_(3,0) !=  0 || tmat_(3,1) != -8 || tmat_(3,2) != 5 || tmat_(3,3) !=  0 ||
          tmat_(4,0) !=  0 || tmat_(4,1) !=  0 || tmat_(4,2) != 9 || tmat_(4,3) != -6 ||
          tmat_(5,0) !=  0 || tmat_(5,1) !=  0 || tmat_(5,2) != 0 || tmat_(5,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( -4  0  0  0 )\n"
                                     "(  0  0  1  0 )\n"
                                     "(  7  4  0  0 )\n"
                                     "(  0 -8  5  0 )\n"
                                     "(  0  0  9 -6 )\n"
                                     "(  0  0  0 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major dense vector multiplication assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::columnVector;

      initialize();

      OBT band0 = blaze::band( tmat_, 0L );

      using AlignedPadded = blaze::CustomVector<int,aligned,padded,columnVector>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 16UL ) );
      AlignedPadded vec( memory.get(), 4UL, 16UL );
      vec[0] =  2;
      vec[1] = -4;
      vec[2] =  0;
      vec[3] =  0;

      band0 *= vec;

      checkSize    ( band0, 4UL );
      checkCapacity( band0, 4UL );
      checkNonZeros( band0, 1UL );
      checkRows    ( tmat_, 6UL );
      checkColumns ( tmat_, 4UL );
      checkNonZeros( tmat_, 9UL );

      if( band0[0] != -4 || band0[1] != 0 || band0[2] != 0 || band0[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( -4 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != -4 || tmat_(0,1) !=  0 || tmat_(0,2) != 0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) !=  0 || tmat_(1,2) != 1 || tmat_(1,3) !=  0 ||
          tmat_(2,0) !=  7 || tmat_(2,1) !=  4 || tmat_(2,2) != 0 || tmat_(2,3) !=  0 ||
          tmat_(3,0) !=  0 || tmat_(3,1) != -8 || tmat_(3,2) != 5 || tmat_(3,3) !=  0 ||
          tmat_(4,0) !=  0 || tmat_(4,1) !=  0 || tmat_(4,2) != 9 || tmat_(4,3) != -6 ||
          tmat_(5,0) !=  0 || tmat_(5,1) !=  0 || tmat_(5,2) != 0 || tmat_(5,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( -4  0  0  0 )\n"
                                     "(  0  0  1  0 )\n"
                                     "(  7  4  0  0 )\n"
                                     "(  0 -8  5  0 )\n"
                                     "(  0  0  9 -6 )\n"
                                     "(  0  0  0 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major dense vector multiplication assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::columnVector;

      initialize();

      OBT band0 = blaze::band( tmat_, 0L );

      using UnalignedUnpadded = blaze::CustomVector<int,unaligned,unpadded,columnVector>;
      std::unique_ptr<int[]> memory( new int[5] );
      UnalignedUnpadded vec( memory.get()+1UL, 4UL );
      vec[0] =  2;
      vec[1] = -4;
      vec[2] =  0;
      vec[3] =  0;

      band0 *= vec;

      checkSize    ( band0, 4UL );
      checkCapacity( band0, 4UL );
      checkNonZeros( band0, 1UL );
      checkRows    ( tmat_, 6UL );
      checkColumns ( tmat_, 4UL );
      checkNonZeros( tmat_, 9UL );

      if( band0[0] != -4 || band0[1] != 0 || band0[2] != 0 || band0[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( -4 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != -4 || tmat_(0,1) !=  0 || tmat_(0,2) != 0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) !=  0 || tmat_(1,2) != 1 || tmat_(1,3) !=  0 ||
          tmat_(2,0) !=  7 || tmat_(2,1) !=  4 || tmat_(2,2) != 0 || tmat_(2,3) !=  0 ||
          tmat_(3,0) !=  0 || tmat_(3,1) != -8 || tmat_(3,2) != 5 || tmat_(3,3) !=  0 ||
          tmat_(4,0) !=  0 || tmat_(4,1) !=  0 || tmat_(4,2) != 9 || tmat_(4,3) != -6 ||
          tmat_(5,0) !=  0 || tmat_(5,1) !=  0 || tmat_(5,2) != 0 || tmat_(5,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( -4  0  0  0 )\n"
                                     "(  0  0  1  0 )\n"
                                     "(  7  4  0  0 )\n"
                                     "(  0 -8  5  0 )\n"
                                     "(  0  0  9 -6 )\n"
                                     "(  0  0  0 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse vector multiplication assignment
   //=====================================================================================

   {
      test_ = "Column-major sparse vector multiplication assignment";

      initialize();

      OBT band0 = blaze::band( tmat_, 0L );

      blaze::CompressedVector<int,blaze::columnVector> vec( 4UL );
      vec[0] =  2;
      vec[1] = -4;

      band0 *= vec;

      checkSize    ( band0, 4UL );
      checkCapacity( band0, 4UL );
      checkNonZeros( band0, 1UL );
      checkRows    ( tmat_, 6UL );
      checkColumns ( tmat_, 4UL );
      checkNonZeros( tmat_, 9UL );

      if( band0[0] != -4 || band0[1] != 0 || band0[2] != 0 || band0[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( -4 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != -4 || tmat_(0,1) !=  0 || tmat_(0,2) != 0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) !=  0 || tmat_(1,2) != 1 || tmat_(1,3) !=  0 ||
          tmat_(2,0) !=  7 || tmat_(2,1) !=  4 || tmat_(2,2) != 0 || tmat_(2,3) !=  0 ||
          tmat_(3,0) !=  0 || tmat_(3,1) != -8 || tmat_(3,2) != 5 || tmat_(3,3) !=  0 ||
          tmat_(4,0) !=  0 || tmat_(4,1) !=  0 || tmat_(4,2) != 9 || tmat_(4,3) != -6 ||
          tmat_(5,0) !=  0 || tmat_(5,1) !=  0 || tmat_(5,2) != 0 || tmat_(5,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( -4  0  0  0 )\n"
                                     "(  0  0  1  0 )\n"
                                     "(  7  4  0  0 )\n"
                                     "(  0 -8  5  0 )\n"
                                     "(  0  0  9 -6 )\n"
                                     "(  0  0  0 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Band division assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the division assignment operators of the Band specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testDivAssign()
{
   //=====================================================================================
   // Row-major Band division assignment
   //=====================================================================================

   {
      test_ = "Row-major Band division assignment";

      initialize();

      BT band0 = blaze::band( mat_, 0L );
      band0 /= blaze::band( mat_, 2L );

      checkSize    ( band0, 4UL );
      checkCapacity( band0, 4UL );
      checkNonZeros( band0, 0UL );
      checkRows    ( mat_ , 4UL );
      checkColumns ( mat_ , 6UL );
      checkNonZeros( mat_ , 8UL );

      if( band0[0] != 0 || band0[1] != 0 || band0[2] != 0 || band0[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) != 7 || mat_(0,3) !=  0 || mat_(0,4) !=  0 || mat_(0,5) !=  0 ||
          mat_(1,0) != 0 || mat_(1,1) != 0 || mat_(1,2) != 4 || mat_(1,3) != -8 || mat_(1,4) !=  0 || mat_(1,5) !=  0 ||
          mat_(2,0) != 0 || mat_(2,1) != 1 || mat_(2,2) != 0 || mat_(2,3) !=  5 || mat_(2,4) !=  9 || mat_(2,5) !=  0 ||
          mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) != 0 || mat_(3,3) !=  0 || mat_(3,4) != -6 || mat_(3,5) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0  7  0  0  0 )\n"
                                     "( 0  0  4 -8  0  0 )\n"
                                     "( 0  1  0  5  9  0 )\n"
                                     "( 0  0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major dense vector division assignment
   //=====================================================================================

   {
      test_ = "Row-major dense vector division assignment (mixed type)";

      initialize();

      BT band0 = blaze::band( mat_, 0L );

      const blaze::DynamicVector<short,blaze::columnVector> vec{ -1, 2, 3, 4 };

      band0 /= vec;

      checkSize    ( band0,  4UL );
      checkCapacity( band0,  4UL );
      checkNonZeros( band0,  2UL );
      checkRows    ( mat_ ,  4UL );
      checkColumns ( mat_ ,  6UL );
      checkNonZeros( mat_ , 10UL );

      if( band0[0] != 2 || band0[1] != 0 || band0[2] != -1 || band0[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( 2 0 -1 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 2 || mat_(0,1) != 0 || mat_(0,2) !=  7 || mat_(0,3) !=  0 || mat_(0,4) !=  0 || mat_(0,5) !=  0 ||
          mat_(1,0) != 0 || mat_(1,1) != 0 || mat_(1,2) !=  4 || mat_(1,3) != -8 || mat_(1,4) !=  0 || mat_(1,5) !=  0 ||
          mat_(2,0) != 0 || mat_(2,1) != 1 || mat_(2,2) != -1 || mat_(2,3) !=  5 || mat_(2,4) !=  9 || mat_(2,5) !=  0 ||
          mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) !=  0 || mat_(3,4) != -6 || mat_(3,5) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 2  0  7  0  0  0 )\n"
                                     "( 0  0  4 -8  0  0 )\n"
                                     "( 0  1 -1  5  9  0 )\n"
                                     "( 0  0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major dense vector division assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::columnVector;

      initialize();

      BT band0 = blaze::band( mat_, 0L );

      using AlignedPadded = blaze::CustomVector<int,aligned,padded,columnVector>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 16UL ) );
      AlignedPadded vec( memory.get(), 4UL, 16UL );
      vec[0] = -1;
      vec[1] =  2;
      vec[2] =  3;
      vec[3] =  4;

      band0 /= vec;

      checkSize    ( band0,  4UL );
      checkCapacity( band0,  4UL );
      checkNonZeros( band0,  2UL );
      checkRows    ( mat_ ,  4UL );
      checkColumns ( mat_ ,  6UL );
      checkNonZeros( mat_ , 10UL );

      if( band0[0] != 2 || band0[1] != 0 || band0[2] != -1 || band0[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( 2 0 -1 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 2 || mat_(0,1) != 0 || mat_(0,2) !=  7 || mat_(0,3) !=  0 || mat_(0,4) !=  0 || mat_(0,5) !=  0 ||
          mat_(1,0) != 0 || mat_(1,1) != 0 || mat_(1,2) !=  4 || mat_(1,3) != -8 || mat_(1,4) !=  0 || mat_(1,5) !=  0 ||
          mat_(2,0) != 0 || mat_(2,1) != 1 || mat_(2,2) != -1 || mat_(2,3) !=  5 || mat_(2,4) !=  9 || mat_(2,5) !=  0 ||
          mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) !=  0 || mat_(3,4) != -6 || mat_(3,5) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0  0  0 )\n"
                                     "( 0  1  0  0 )\n"
                                     "( 2  0 -1  0 )\n"
                                     "( 0  4  5 -6 )\n"
                                     "( 7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major dense vector division assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::columnVector;

      initialize();

      BT band0 = blaze::band( mat_, 0L );

      using UnalignedUnpadded = blaze::CustomVector<int,unaligned,unpadded,columnVector>;
      std::unique_ptr<int[]> memory( new int[5] );
      UnalignedUnpadded vec( memory.get()+1UL, 4UL );
      vec[0] = -1;
      vec[1] =  2;
      vec[2] =  3;
      vec[3] =  4;

      band0 /= vec;

      checkSize    ( band0,  4UL );
      checkCapacity( band0,  4UL );
      checkNonZeros( band0,  2UL );
      checkRows    ( mat_ ,  4UL );
      checkColumns ( mat_ ,  6UL );
      checkNonZeros( mat_ , 10UL );

      if( band0[0] != 2 || band0[1] != 0 || band0[2] != -1 || band0[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( 2 0 -1 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 2 || mat_(0,1) != 0 || mat_(0,2) !=  7 || mat_(0,3) !=  0 || mat_(0,4) !=  0 || mat_(0,5) !=  0 ||
          mat_(1,0) != 0 || mat_(1,1) != 0 || mat_(1,2) !=  4 || mat_(1,3) != -8 || mat_(1,4) !=  0 || mat_(1,5) !=  0 ||
          mat_(2,0) != 0 || mat_(2,1) != 1 || mat_(2,2) != -1 || mat_(2,3) !=  5 || mat_(2,4) !=  9 || mat_(2,5) !=  0 ||
          mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) !=  0 || mat_(3,4) != -6 || mat_(3,5) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0  0  0 )\n"
                                     "( 0  1  0  0 )\n"
                                     "( 2  0 -1  0 )\n"
                                     "( 0  4  5 -6 )\n"
                                     "( 7 -8  9 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major Band division assignment
   //=====================================================================================

   {
      test_ = "Column-major Band division assignment";

      initialize();

      OBT band0 = blaze::band( tmat_, 0L );
      band0 /= blaze::band( tmat_, -2L );

      checkSize    ( band0, 4UL );
      checkCapacity( band0, 4UL );
      checkNonZeros( band0, 0UL );
      checkRows    ( tmat_, 6UL );
      checkColumns ( tmat_, 4UL );
      checkNonZeros( tmat_, 8UL );

      if( band0[0] != 0 || band0[1] != 0 || band0[2] != 0 || band0[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) != 0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) != 0 || tmat_(1,1) !=  0 || tmat_(1,2) != 1 || tmat_(1,3) !=  0 ||
          tmat_(2,0) != 7 || tmat_(2,1) !=  4 || tmat_(2,2) != 0 || tmat_(2,3) !=  0 ||
          tmat_(3,0) != 0 || tmat_(3,1) != -8 || tmat_(3,2) != 5 || tmat_(3,3) !=  0 ||
          tmat_(4,0) != 0 || tmat_(4,1) !=  0 || tmat_(4,2) != 9 || tmat_(4,3) != -6 ||
          tmat_(5,0) != 0 || tmat_(5,1) !=  0 || tmat_(5,2) != 0 || tmat_(5,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0  0  0 )\n"
                                     "( 0  0  1  0 )\n"
                                     "( 7  4  0  0 )\n"
                                     "( 0 -8  5  0 )\n"
                                     "( 0  0  9 -6 )\n"
                                     "( 0  0  0 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense vector division assignment
   //=====================================================================================

   {
      test_ = "Column-major dense vector division assignment (mixed type)";

      initialize();

      OBT band0 = blaze::band( tmat_, 0L );

      const blaze::DynamicVector<short,blaze::columnVector> vec{ -1, 2, 3, 4 };

      band0 /= vec;

      checkSize    ( band0,  4UL );
      checkCapacity( band0,  4UL );
      checkNonZeros( band0,  2UL );
      checkRows    ( tmat_,  6UL );
      checkColumns ( tmat_,  4UL );
      checkNonZeros( tmat_, 10UL );

      if( band0[0] != 2 || band0[1] != 0 || band0[2] != -1 || band0[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( 2 0 -1 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 2 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) != 0 || tmat_(1,1) !=  0 || tmat_(1,2) !=  1 || tmat_(1,3) !=  0 ||
          tmat_(2,0) != 7 || tmat_(2,1) !=  4 || tmat_(2,2) != -1 || tmat_(2,3) !=  0 ||
          tmat_(3,0) != 0 || tmat_(3,1) != -8 || tmat_(3,2) !=  5 || tmat_(3,3) !=  0 ||
          tmat_(4,0) != 0 || tmat_(4,1) !=  0 || tmat_(4,2) !=  9 || tmat_(4,3) != -6 ||
          tmat_(5,0) != 0 || tmat_(5,1) !=  0 || tmat_(5,2) !=  0 || tmat_(5,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 2  0  0  0 )\n"
                                     "( 0  0  1  0 )\n"
                                     "( 7  4 -1  0 )\n"
                                     "( 0 -8  5  0 )\n"
                                     "( 0  0  9 -6 )\n"
                                     "( 0  0  0 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major dense vector division assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::columnVector;

      initialize();

      OBT band0 = blaze::band( tmat_, 0L );

      using AlignedPadded = blaze::CustomVector<int,aligned,padded,columnVector>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 16UL ) );
      AlignedPadded vec( memory.get(), 4UL, 16UL );
      vec[0] = -1;
      vec[1] =  2;
      vec[2] =  3;
      vec[3] =  4;

      band0 /= vec;

      checkSize    ( band0,  4UL );
      checkCapacity( band0,  4UL );
      checkNonZeros( band0,  2UL );
      checkRows    ( tmat_,  6UL );
      checkColumns ( tmat_,  4UL );
      checkNonZeros( tmat_, 10UL );

      if( band0[0] != 2 || band0[1] != 0 || band0[2] != -1 || band0[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( 2 0 -1 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 2 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) != 0 || tmat_(1,1) !=  0 || tmat_(1,2) !=  1 || tmat_(1,3) !=  0 ||
          tmat_(2,0) != 7 || tmat_(2,1) !=  4 || tmat_(2,2) != -1 || tmat_(2,3) !=  0 ||
          tmat_(3,0) != 0 || tmat_(3,1) != -8 || tmat_(3,2) !=  5 || tmat_(3,3) !=  0 ||
          tmat_(4,0) != 0 || tmat_(4,1) !=  0 || tmat_(4,2) !=  9 || tmat_(4,3) != -6 ||
          tmat_(5,0) != 0 || tmat_(5,1) !=  0 || tmat_(5,2) !=  0 || tmat_(5,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 2  0  0  0 )\n"
                                     "( 0  0  1  0 )\n"
                                     "( 7  4 -1  0 )\n"
                                     "( 0 -8  5  0 )\n"
                                     "( 0  0  9 -6 )\n"
                                     "( 0  0  0 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major dense vector division assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::columnVector;

      initialize();

      OBT band0 = blaze::band( tmat_, 0L );

      using UnalignedUnpadded = blaze::CustomVector<int,unaligned,unpadded,columnVector>;
      std::unique_ptr<int[]> memory( new int[5] );
      UnalignedUnpadded vec( memory.get()+1UL, 4UL );
      vec[0] = -1;
      vec[1] =  2;
      vec[2] =  3;
      vec[3] =  4;

      band0 /= vec;

      checkSize    ( band0,  4UL );
      checkCapacity( band0,  4UL );
      checkNonZeros( band0,  2UL );
      checkRows    ( tmat_,  6UL );
      checkColumns ( tmat_,  4UL );
      checkNonZeros( tmat_, 10UL );

      if( band0[0] != 2 || band0[1] != 0 || band0[2] != -1 || band0[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( 2 0 -1 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 2 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) != 0 || tmat_(1,1) !=  0 || tmat_(1,2) !=  1 || tmat_(1,3) !=  0 ||
          tmat_(2,0) != 7 || tmat_(2,1) !=  4 || tmat_(2,2) != -1 || tmat_(2,3) !=  0 ||
          tmat_(3,0) != 0 || tmat_(3,1) != -8 || tmat_(3,2) !=  5 || tmat_(3,3) !=  0 ||
          tmat_(4,0) != 0 || tmat_(4,1) !=  0 || tmat_(4,2) !=  9 || tmat_(4,3) != -6 ||
          tmat_(5,0) != 0 || tmat_(5,1) !=  0 || tmat_(5,2) !=  0 || tmat_(5,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 2  0  0  0 )\n"
                                     "( 0  0  1  0 )\n"
                                     "( 7  4 -1  0 )\n"
                                     "( 0 -8  5  0 )\n"
                                     "( 0  0  9 -6 )\n"
                                     "( 0  0  0 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Band cross product assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the cross product assignment operators of the Band
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testCrossAssign()
{
   //=====================================================================================
   // Row-major Band cross product assignment
   //=====================================================================================

   {
      test_ = "Row-major Band cross product assignment";

      MT mat{ { 2, 1, 0, 0 }, { 0, 0, 0, 0 }, { 0, 0, -1, -2 } };

      BT band0 = blaze::band( mat, 0L );
      band0 %= blaze::band( mat, 1L );

      checkSize    ( band0, 3UL );
      checkCapacity( band0, 3UL );
      checkNonZeros( band0, 1UL );
      checkRows    ( mat  , 3UL );
      checkColumns ( mat  , 4UL );
      checkNonZeros( mat  , 3UL );

      if( band0[0] != 0 || band0[1] != 3 || band0[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cross product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( 0 3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat(0,0) != 0 || mat(0,1) != 1 || mat(0,2) != 0 || mat(0,3) !=  0 ||
          mat(1,0) != 0 || mat(1,1) != 3 || mat(1,2) != 0 || mat(1,3) !=  0 ||
          mat(2,0) != 0 || mat(2,1) != 0 || mat(2,2) != 0 || mat(2,3) != -2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cross product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0  1  0  0 )\n"
                                     "( 0  3  0  0 )\n"
                                     "( 0  0  0 -2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major dense vector cross product assignment
   //=====================================================================================

   {
      test_ = "Row-major dense vector cross product assignment (mixed type)";

      MT mat{ { 2, 1, 0, 0 }, { 0, 0, 0, 0 }, { 0, 0, -1, -2 } };

      BT band0 = blaze::band( mat, 0L );

      const blaze::DynamicVector<short,blaze::columnVector> vec{ 1, 0, -2 };

      band0 %= vec;

      checkSize    ( band0, 3UL );
      checkCapacity( band0, 3UL );
      checkNonZeros( band0, 1UL );
      checkRows    ( mat  , 3UL );
      checkColumns ( mat  , 4UL );
      checkNonZeros( mat  , 3UL );

      if( band0[0] != 0 || band0[1] != 3 || band0[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cross product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( 0 3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat(0,0) != 0 || mat(0,1) != 1 || mat(0,2) != 0 || mat(0,3) !=  0 ||
          mat(1,0) != 0 || mat(1,1) != 3 || mat(1,2) != 0 || mat(1,3) !=  0 ||
          mat(2,0) != 0 || mat(2,1) != 0 || mat(2,2) != 0 || mat(2,3) != -2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cross product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0  1  0  0 )\n"
                                     "( 0  3  0  0 )\n"
                                     "( 0  0  0 -2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major dense vector cross product assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::columnVector;

      MT mat{ { 2, 1, 0, 0 }, { 0, 0, 0, 0 }, { 0, 0, -1, -2 } };

      BT band0 = blaze::band( mat, 0L );

      using AlignedPadded = blaze::CustomVector<int,aligned,padded,columnVector>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 16UL ) );
      AlignedPadded vec( memory.get(), 3UL, 16UL );
      vec[0] =  1;
      vec[1] =  0;
      vec[2] = -2;

      band0 %= vec;

      checkSize    ( band0, 3UL );
      checkCapacity( band0, 3UL );
      checkNonZeros( band0, 1UL );
      checkRows    ( mat  , 3UL );
      checkColumns ( mat  , 4UL );
      checkNonZeros( mat  , 3UL );

      if( band0[0] != 0 || band0[1] != 3 || band0[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cross product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( 0 3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat(0,0) != 0 || mat(0,1) != 1 || mat(0,2) != 0 || mat(0,3) !=  0 ||
          mat(1,0) != 0 || mat(1,1) != 3 || mat(1,2) != 0 || mat(1,3) !=  0 ||
          mat(2,0) != 0 || mat(2,1) != 0 || mat(2,2) != 0 || mat(2,3) != -2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cross product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0  1  0  0 )\n"
                                     "( 0  3  0  0 )\n"
                                     "( 0  0  0 -2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major dense vector cross product assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::columnVector;

      MT mat{ { 2, 1, 0, 0 }, { 0, 0, 0, 0 }, { 0, 0, -1, -2 } };

      BT band0 = blaze::band( mat, 0L );

      using UnalignedUnpadded = blaze::CustomVector<int,unaligned,unpadded,columnVector>;
      std::unique_ptr<int[]> memory( new int[4] );
      UnalignedUnpadded vec( memory.get()+1UL, 3UL );
      vec[0] =  1;
      vec[1] =  0;
      vec[2] = -2;

      band0 %= vec;

      checkSize    ( band0, 3UL );
      checkCapacity( band0, 3UL );
      checkNonZeros( band0, 1UL );
      checkRows    ( mat  , 3UL );
      checkColumns ( mat  , 4UL );
      checkNonZeros( mat  , 3UL );

      if( band0[0] != 0 || band0[1] != 3 || band0[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cross product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( 0 3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat(0,0) != 0 || mat(0,1) != 1 || mat(0,2) != 0 || mat(0,3) !=  0 ||
          mat(1,0) != 0 || mat(1,1) != 3 || mat(1,2) != 0 || mat(1,3) !=  0 ||
          mat(2,0) != 0 || mat(2,1) != 0 || mat(2,2) != 0 || mat(2,3) != -2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cross product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0  1  0  0 )\n"
                                     "( 0  3  0  0 )\n"
                                     "( 0  0  0 -2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse vector cross product assignment
   //=====================================================================================

   {
      test_ = "Row-major sparse vector cross product assignment";

      MT mat{ { 2, 1, 0, 0 }, { 0, 0, 0, 0 }, { 0, 0, -1, -2 } };

      BT band0 = blaze::band( mat, 0L );

      blaze::CompressedVector<int,blaze::columnVector> vec( 3UL );
      vec[0] =  1;
      vec[2] = -2;

      band0 %= vec;

      checkSize    ( band0, 3UL );
      checkCapacity( band0, 3UL );
      checkNonZeros( band0, 1UL );
      checkRows    ( mat  , 3UL );
      checkColumns ( mat  , 4UL );
      checkNonZeros( mat  , 3UL );

      if( band0[0] != 0 || band0[1] != 3 || band0[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cross product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( 0 3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat(0,0) != 0 || mat(0,1) != 1 || mat(0,2) != 0 || mat(0,3) !=  0 ||
          mat(1,0) != 0 || mat(1,1) != 3 || mat(1,2) != 0 || mat(1,3) !=  0 ||
          mat(2,0) != 0 || mat(2,1) != 0 || mat(2,2) != 0 || mat(2,3) != -2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cross product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0  1  0  0 )\n"
                                     "( 0  3  0  0 )\n"
                                     "( 0  0  0 -2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major Band cross product assignment
   //=====================================================================================

   {
      test_ = "Column-major Band cross product assignment";

      OMT mat{ { 2, 0, 0 }, { 1, 0, 0 }, { 0, 0, -1 }, { 0, 0, -2 } };

      OBT band0 = blaze::band( mat, 0L );
      band0 %= blaze::band( mat, -1L );

      checkSize    ( band0, 3UL );
      checkCapacity( band0, 3UL );
      checkNonZeros( band0, 1UL );
      checkRows    ( mat  , 4UL );
      checkColumns ( mat  , 3UL );
      checkNonZeros( mat  , 3UL );

      if( band0[0] != 0 || band0[1] != 3 || band0[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cross product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( 0 3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat(0,0) != 0 || mat(0,1) != 0 || mat(0,2) !=  0 ||
          mat(1,0) != 1 || mat(1,1) != 3 || mat(1,2) !=  0 ||
          mat(2,0) != 0 || mat(2,1) != 0 || mat(2,2) !=  0 ||
          mat(3,0) != 0 || mat(3,1) != 0 || mat(3,2) != -2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cross product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0  0  0 )\n"
                                     "( 1  3  0 )\n"
                                     "( 0  0  0 )\n"
                                     "( 0  0 -2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense vector cross product assignment
   //=====================================================================================

   {
      test_ = "Column-major dense vector cross product assignment (mixed type)";

      OMT mat{ { 2, 0, 0 }, { 1, 0, 0 }, { 0, 0, -1 }, { 0, 0, -2 } };

      OBT band0 = blaze::band( mat, 0L );

      const blaze::DynamicVector<short,blaze::columnVector> vec{ 1, 0, -2 };

      band0 %= vec;

      checkSize    ( band0, 3UL );
      checkCapacity( band0, 3UL );
      checkNonZeros( band0, 1UL );
      checkRows    ( mat  , 4UL );
      checkColumns ( mat  , 3UL );
      checkNonZeros( mat  , 3UL );

      if( band0[0] != 0 || band0[1] != 3 || band0[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cross product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( 0 3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat(0,0) != 0 || mat(0,1) != 0 || mat(0,2) !=  0 ||
          mat(1,0) != 1 || mat(1,1) != 3 || mat(1,2) !=  0 ||
          mat(2,0) != 0 || mat(2,1) != 0 || mat(2,2) !=  0 ||
          mat(3,0) != 0 || mat(3,1) != 0 || mat(3,2) != -2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cross product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0  0  0 )\n"
                                     "( 1  3  0 )\n"
                                     "( 0  0  0 )\n"
                                     "( 0  0 -2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major dense vector cross product assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::columnVector;

      OMT mat{ { 2, 0, 0 }, { 1, 0, 0 }, { 0, 0, -1 }, { 0, 0, -2 } };

      OBT band0 = blaze::band( mat, 0L );

      using AlignedPadded = blaze::CustomVector<int,aligned,padded,columnVector>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 16UL ) );
      AlignedPadded vec( memory.get(), 3UL, 16UL );
      vec[0] =  1;
      vec[1] =  0;
      vec[2] = -2;

      band0 %= vec;

      checkSize    ( band0, 3UL );
      checkCapacity( band0, 3UL );
      checkNonZeros( band0, 1UL );
      checkRows    ( mat  , 4UL );
      checkColumns ( mat  , 3UL );
      checkNonZeros( mat  , 3UL );

      if( band0[0] != 0 || band0[1] != 3 || band0[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cross product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( 0 3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat(0,0) != 0 || mat(0,1) != 0 || mat(0,2) !=  0 ||
          mat(1,0) != 1 || mat(1,1) != 3 || mat(1,2) !=  0 ||
          mat(2,0) != 0 || mat(2,1) != 0 || mat(2,2) !=  0 ||
          mat(3,0) != 0 || mat(3,1) != 0 || mat(3,2) != -2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cross product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0  0  0 )\n"
                                     "( 1  3  0 )\n"
                                     "( 0  0  0 )\n"
                                     "( 0  0 -2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major dense vector cross product assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::columnVector;

      OMT mat{ { 2, 0, 0 }, { 1, 0, 0 }, { 0, 0, -1 }, { 0, 0, -2 } };

      OBT band0 = blaze::band( mat, 0L );

      using UnalignedUnpadded = blaze::CustomVector<int,unaligned,unpadded,columnVector>;
      std::unique_ptr<int[]> memory( new int[4] );
      UnalignedUnpadded vec( memory.get()+1UL, 3UL );
      vec[0] =  1;
      vec[1] =  0;
      vec[2] = -2;

      band0 %= vec;

      checkSize    ( band0, 3UL );
      checkCapacity( band0, 3UL );
      checkNonZeros( band0, 1UL );
      checkRows    ( mat  , 4UL );
      checkColumns ( mat  , 3UL );
      checkNonZeros( mat  , 3UL );

      if( band0[0] != 0 || band0[1] != 3 || band0[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cross product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( 0 3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat(0,0) != 0 || mat(0,1) != 0 || mat(0,2) !=  0 ||
          mat(1,0) != 1 || mat(1,1) != 3 || mat(1,2) !=  0 ||
          mat(2,0) != 0 || mat(2,1) != 0 || mat(2,2) !=  0 ||
          mat(3,0) != 0 || mat(3,1) != 0 || mat(3,2) != -2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cross product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0  0  0 )\n"
                                     "( 1  3  0 )\n"
                                     "( 0  0  0 )\n"
                                     "( 0  0 -2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse vector cross product assignment
   //=====================================================================================

   {
      test_ = "Column-major sparse vector cross product assignment";

      OMT mat{ { 2, 0, 0 }, { 1, 0, 0 }, { 0, 0, -1 }, { 0, 0, -2 } };

      OBT band0 = blaze::band( mat, 0L );

      blaze::CompressedVector<int,blaze::columnVector> vec( 3UL );
      vec[0] =  1;
      vec[2] = -2;

      band0 %= vec;

      checkSize    ( band0, 3UL );
      checkCapacity( band0, 3UL );
      checkNonZeros( band0, 1UL );
      checkRows    ( mat  , 4UL );
      checkColumns ( mat  , 3UL );
      checkNonZeros( mat  , 3UL );

      if( band0[0] != 0 || band0[1] != 3 || band0[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cross product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( 0 3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat(0,0) != 0 || mat(0,1) != 0 || mat(0,2) !=  0 ||
          mat(1,0) != 1 || mat(1,1) != 3 || mat(1,2) !=  0 ||
          mat(2,0) != 0 || mat(2,1) != 0 || mat(2,2) !=  0 ||
          mat(3,0) != 0 || mat(3,1) != 0 || mat(3,2) != -2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cross product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0  0  0 )\n"
                                     "( 1  3  0 )\n"
                                     "( 0  0  0 )\n"
                                     "( 0  0 -2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of all Band (self-)scaling operations.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all available ways to scale an instance of the Band
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testScaling()
{
   //=====================================================================================
   // Row-major self-scaling (v*=s)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (v*=s)";

      initialize();

      BT band0 = blaze::band( mat_, 0L );

      band0 *= 3;

      checkSize    ( band0,  4UL );
      checkCapacity( band0,  4UL );
      checkNonZeros( band0,  2UL );
      checkRows    ( mat_ ,  4UL );
      checkColumns ( mat_ ,  6UL );
      checkNonZeros( mat_ , 10UL );

      if( band0[0] != -6 || band0[1] != 0 || band0[2] != -9 || band0[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( -6 0 -9 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != -6 || mat_(0,1) != 0 || mat_(0,2) !=  7 || mat_(0,3) !=  0 || mat_(0,4) !=  0 || mat_(0,5) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 0 || mat_(1,2) !=  4 || mat_(1,3) != -8 || mat_(1,4) !=  0 || mat_(1,5) !=  0 ||
          mat_(2,0) !=  0 || mat_(2,1) != 1 || mat_(2,2) != -9 || mat_(2,3) !=  5 || mat_(2,4) !=  9 || mat_(2,5) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) !=  0 || mat_(3,4) != -6 || mat_(3,5) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( -6  0  7  0  0  0 )\n"
                                     "(  0  0  4 -8  0  0 )\n"
                                     "(  0  1 -9  5  9  0 )\n"
                                     "(  0  0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major self-scaling (v=v*s)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (v=v*s)";

      initialize();

      BT band0 = blaze::band( mat_, 0L );

      band0 = band0 * 3;

      checkSize    ( band0,  4UL );
      checkCapacity( band0,  4UL );
      checkNonZeros( band0,  2UL );
      checkRows    ( mat_ ,  4UL );
      checkColumns ( mat_ ,  6UL );
      checkNonZeros( mat_ , 10UL );

      if( band0[0] != -6 || band0[1] != 0 || band0[2] != -9 || band0[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( -6 0 -9 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != -6 || mat_(0,1) != 0 || mat_(0,2) !=  7 || mat_(0,3) !=  0 || mat_(0,4) !=  0 || mat_(0,5) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 0 || mat_(1,2) !=  4 || mat_(1,3) != -8 || mat_(1,4) !=  0 || mat_(1,5) !=  0 ||
          mat_(2,0) !=  0 || mat_(2,1) != 1 || mat_(2,2) != -9 || mat_(2,3) !=  5 || mat_(2,4) !=  9 || mat_(2,5) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) !=  0 || mat_(3,4) != -6 || mat_(3,5) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( -6  0  7  0  0  0 )\n"
                                     "(  0  0  4 -8  0  0 )\n"
                                     "(  0  1 -9  5  9  0 )\n"
                                     "(  0  0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major self-scaling (v=s*v)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (v=s*v)";

      initialize();

      BT band0 = blaze::band( mat_, 0L );

      band0 = 3 * band0;

      checkSize    ( band0,  4UL );
      checkCapacity( band0,  4UL );
      checkNonZeros( band0,  2UL );
      checkRows    ( mat_ ,  4UL );
      checkColumns ( mat_ ,  6UL );
      checkNonZeros( mat_ , 10UL );

      if( band0[0] != -6 || band0[1] != 0 || band0[2] != -9 || band0[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( -6 0 -9 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != -6 || mat_(0,1) != 0 || mat_(0,2) !=  7 || mat_(0,3) !=  0 || mat_(0,4) !=  0 || mat_(0,5) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 0 || mat_(1,2) !=  4 || mat_(1,3) != -8 || mat_(1,4) !=  0 || mat_(1,5) !=  0 ||
          mat_(2,0) !=  0 || mat_(2,1) != 1 || mat_(2,2) != -9 || mat_(2,3) !=  5 || mat_(2,4) !=  9 || mat_(2,5) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) !=  0 || mat_(3,4) != -6 || mat_(3,5) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( -6  0  7  0  0  0 )\n"
                                     "(  0  0  4 -8  0  0 )\n"
                                     "(  0  1 -9  5  9  0 )\n"
                                     "(  0  0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major self-scaling (v/=s)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (v/=s)";

      initialize();

      BT band0 = blaze::band( mat_, 0L );

      band0 /= 0.5;

      checkSize    ( band0,  4UL );
      checkCapacity( band0,  4UL );
      checkNonZeros( band0,  2UL );
      checkRows    ( mat_ ,  4UL );
      checkColumns ( mat_ ,  6UL );
      checkNonZeros( mat_ , 10UL );

      if( band0[0] != -4 || band0[1] != 0 || band0[2] != -6 || band0[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( -4 0 -6 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != -4 || mat_(0,1) != 0 || mat_(0,2) !=  7 || mat_(0,3) !=  0 || mat_(0,4) !=  0 || mat_(0,5) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 0 || mat_(1,2) !=  4 || mat_(1,3) != -8 || mat_(1,4) !=  0 || mat_(1,5) !=  0 ||
          mat_(2,0) !=  0 || mat_(2,1) != 1 || mat_(2,2) != -6 || mat_(2,3) !=  5 || mat_(2,4) !=  9 || mat_(2,5) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) !=  0 || mat_(3,4) != -6 || mat_(3,5) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( -4  0  7  0  0  0 )\n"
                                     "(  0  0  4 -8  0  0 )\n"
                                     "(  0  1 -6  5  9  0 )\n"
                                     "(  0  0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major self-scaling (v=v/s)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (v=v/s)";

      initialize();

      BT band0 = blaze::band( mat_, 0L );

      band0 = band0 / 0.5;

      checkSize    ( band0,  4UL );
      checkCapacity( band0,  4UL );
      checkNonZeros( band0,  2UL );
      checkRows    ( mat_ ,  4UL );
      checkColumns ( mat_ ,  6UL );
      checkNonZeros( mat_ , 10UL );

      if( band0[0] != -4 || band0[1] != 0 || band0[2] != -6 || band0[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( -4 0 -6 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != -4 || mat_(0,1) != 0 || mat_(0,2) !=  7 || mat_(0,3) !=  0 || mat_(0,4) !=  0 || mat_(0,5) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 0 || mat_(1,2) !=  4 || mat_(1,3) != -8 || mat_(1,4) !=  0 || mat_(1,5) !=  0 ||
          mat_(2,0) !=  0 || mat_(2,1) != 1 || mat_(2,2) != -6 || mat_(2,3) !=  5 || mat_(2,4) !=  9 || mat_(2,5) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) !=  0 || mat_(3,4) != -6 || mat_(3,5) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( -4  0  7  0  0  0 )\n"
                                     "(  0  0  4 -8  0  0 )\n"
                                     "(  0  1 -6  5  9  0 )\n"
                                     "(  0  0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major Band::scale()
   //=====================================================================================

   {
      test_ = "Row-major Band::scale()";

      initialize();

      // Integral scaling the 1st upper band
      {
         BT band1 = blaze::band( mat_, 1L );
         band1.scale( 3 );

         checkSize    ( band1,  4UL );
         checkCapacity( band1,  4UL );
         checkNonZeros( band1,  3UL );
         checkRows    ( mat_ ,  4UL );
         checkColumns ( mat_ ,  6UL );
         checkNonZeros( mat_ , 10UL );

         if( band1[0] != 0 || band1[1] != 12 || band1[2] != 15 || band1[3] != -18 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Integral scale operation of 1st upper band failed\n"
                << " Details:\n"
                << "   Result:\n" << band1 << "\n"
                << "   Expected result:\n( 0 12 15 -18 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( mat_(0,0) != -2 || mat_(0,1) != 0 || mat_(0,2) !=  7 || mat_(0,3) !=  0 || mat_(0,4) !=   0 || mat_(0,5) !=  0 ||
             mat_(1,0) !=  0 || mat_(1,1) != 0 || mat_(1,2) != 12 || mat_(1,3) != -8 || mat_(1,4) !=   0 || mat_(1,5) !=  0 ||
             mat_(2,0) !=  0 || mat_(2,1) != 1 || mat_(2,2) != -3 || mat_(2,3) != 15 || mat_(2,4) !=   9 || mat_(2,5) !=  0 ||
             mat_(3,0) !=  0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) !=  0 || mat_(3,4) != -18 || mat_(3,5) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Integral scale operation of 3rd row failed\n"
                << " Details:\n"
                << "   Result:\n" << mat_ << "\n"
                << "   Expected result:\n( -2  0  7  0   0  0 )\n"
                                        "(  0  0 12 -8   0  0 )\n"
                                        "(  0  1 -3 15   9  0 )\n"
                                        "(  0  0  0  0 -18 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Floating point scaling the 1st upper band
      {
         BT band1 = blaze::band( mat_, 1L );
         band1.scale( 0.5 );

         checkSize    ( band1,  4UL );
         checkCapacity( band1,  4UL );
         checkNonZeros( band1,  3UL );
         checkRows    ( mat_ ,  4UL );
         checkColumns ( mat_ ,  6UL );
         checkNonZeros( mat_ , 10UL );

         if( band1[0] != 0 || band1[1] != 6 || band1[2] != 7 || band1[3] != -9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Floating point scale operation of 3rd row failed\n"
                << " Details:\n"
                << "   Result:\n" << band1 << "\n"
                << "   Expected result:\n( 0 6 7 -9 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( mat_(0,0) != -2 || mat_(0,1) != 0 || mat_(0,2) !=  7 || mat_(0,3) !=  0 || mat_(0,4) !=  0 || mat_(0,5) !=  0 ||
             mat_(1,0) !=  0 || mat_(1,1) != 0 || mat_(1,2) !=  6 || mat_(1,3) != -8 || mat_(1,4) !=  0 || mat_(1,5) !=  0 ||
             mat_(2,0) !=  0 || mat_(2,1) != 1 || mat_(2,2) != -3 || mat_(2,3) !=  7 || mat_(2,4) !=  9 || mat_(2,5) !=  0 ||
             mat_(3,0) !=  0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) !=  0 || mat_(3,4) != -9 || mat_(3,5) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Floating point scale operation of 3rd row failed\n"
                << " Details:\n"
                << "   Result:\n" << mat_ << "\n"
                << "   Expected result:\n( -2  0  7  0  0  0 )\n"
                                        "(  0  0  6 -8  0  0 )\n"
                                        "(  0  1 -3  7  9  0 )\n"
                                        "(  0  0  0  0 -9 10 )\n";
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

      OBT band0 = blaze::band( tmat_, 0L );

      band0 *= 3;

      checkSize    ( band0,  4UL );
      checkCapacity( band0,  4UL );
      checkNonZeros( band0,  2UL );
      checkRows    ( tmat_,  6UL );
      checkColumns ( tmat_,  4UL );
      checkNonZeros( tmat_, 10UL );

      if( band0[0] != -6 || band0[1] != 0 || band0[2] != -9 || band0[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( -6 0 -9 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != -6 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) !=  0 || tmat_(1,2) !=  1 || tmat_(1,3) !=  0 ||
          tmat_(2,0) !=  7 || tmat_(2,1) !=  4 || tmat_(2,2) != -9 || tmat_(2,3) !=  0 ||
          tmat_(3,0) !=  0 || tmat_(3,1) != -8 || tmat_(3,2) !=  5 || tmat_(3,3) !=  0 ||
          tmat_(4,0) !=  0 || tmat_(4,1) !=  0 || tmat_(4,2) !=  9 || tmat_(4,3) != -6 ||
          tmat_(5,0) !=  0 || tmat_(5,1) !=  0 || tmat_(5,2) !=  0 || tmat_(5,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( -6  0  0  0 )\n"
                                     "(  0  0  1  0 )\n"
                                     "(  7  4 -9  0 )\n"
                                     "(  0 -8  5  0 )\n"
                                     "(  0  0  9 -6 )\n"
                                     "(  0  0  0 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (v=v*s)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (v=v*s)";

      initialize();

      OBT band0 = blaze::band( tmat_, 0L );

      band0 = band0 * 3;

      checkSize    ( band0,  4UL );
      checkCapacity( band0,  4UL );
      checkNonZeros( band0,  2UL );
      checkRows    ( tmat_,  6UL );
      checkColumns ( tmat_,  4UL );
      checkNonZeros( tmat_, 10UL );

      if( band0[0] != -6 || band0[1] != 0 || band0[2] != -9 || band0[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( -6 0 -9 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != -6 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) !=  0 || tmat_(1,2) !=  1 || tmat_(1,3) !=  0 ||
          tmat_(2,0) !=  7 || tmat_(2,1) !=  4 || tmat_(2,2) != -9 || tmat_(2,3) !=  0 ||
          tmat_(3,0) !=  0 || tmat_(3,1) != -8 || tmat_(3,2) !=  5 || tmat_(3,3) !=  0 ||
          tmat_(4,0) !=  0 || tmat_(4,1) !=  0 || tmat_(4,2) !=  9 || tmat_(4,3) != -6 ||
          tmat_(5,0) !=  0 || tmat_(5,1) !=  0 || tmat_(5,2) !=  0 || tmat_(5,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( -6  0  0  0 )\n"
                                     "(  0  0  1  0 )\n"
                                     "(  7  4 -9  0 )\n"
                                     "(  0 -8  5  0 )\n"
                                     "(  0  0  9 -6 )\n"
                                     "(  0  0  0 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (v=s*v)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (v=s*v)";

      initialize();

      OBT band0 = blaze::band( tmat_, 0L );

      band0 = 3 * band0;

      checkSize    ( band0,  4UL );
      checkCapacity( band0,  4UL );
      checkNonZeros( band0,  2UL );
      checkRows    ( tmat_,  6UL );
      checkColumns ( tmat_,  4UL );
      checkNonZeros( tmat_, 10UL );

      if( band0[0] != -6 || band0[1] != 0 || band0[2] != -9 || band0[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( -6 0 -9 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != -6 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) !=  0 || tmat_(1,2) !=  1 || tmat_(1,3) !=  0 ||
          tmat_(2,0) !=  7 || tmat_(2,1) !=  4 || tmat_(2,2) != -9 || tmat_(2,3) !=  0 ||
          tmat_(3,0) !=  0 || tmat_(3,1) != -8 || tmat_(3,2) !=  5 || tmat_(3,3) !=  0 ||
          tmat_(4,0) !=  0 || tmat_(4,1) !=  0 || tmat_(4,2) !=  9 || tmat_(4,3) != -6 ||
          tmat_(5,0) !=  0 || tmat_(5,1) !=  0 || tmat_(5,2) !=  0 || tmat_(5,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( -6  0  0  0 )\n"
                                     "(  0  0  1  0 )\n"
                                     "(  7  4 -9  0 )\n"
                                     "(  0 -8  5  0 )\n"
                                     "(  0  0  9 -6 )\n"
                                     "(  0  0  0 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (v/=s)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (v/=s)";

      initialize();

      OBT band0 = blaze::band( tmat_, 0L );

      band0 /= 0.5;

      checkSize    ( band0,  4UL );
      checkCapacity( band0,  4UL );
      checkNonZeros( band0,  2UL );
      checkRows    ( tmat_,  6UL );
      checkColumns ( tmat_,  4UL );
      checkNonZeros( tmat_, 10UL );

      if( band0[0] != -4 || band0[1] != 0 || band0[2] != -6 || band0[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( -4 0 -6 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != -4 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) !=  0 || tmat_(1,2) !=  1 || tmat_(1,3) !=  0 ||
          tmat_(2,0) !=  7 || tmat_(2,1) !=  4 || tmat_(2,2) != -6 || tmat_(2,3) !=  0 ||
          tmat_(3,0) !=  0 || tmat_(3,1) != -8 || tmat_(3,2) !=  5 || tmat_(3,3) !=  0 ||
          tmat_(4,0) !=  0 || tmat_(4,1) !=  0 || tmat_(4,2) !=  9 || tmat_(4,3) != -6 ||
          tmat_(5,0) !=  0 || tmat_(5,1) !=  0 || tmat_(5,2) !=  0 || tmat_(5,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( -4  0  0  0 )\n"
                                     "(  0  0  1  0 )\n"
                                     "(  7  4 -6  0 )\n"
                                     "(  0 -8  5  0 )\n"
                                     "(  0  0  9 -6 )\n"
                                     "(  0  0  0 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (v=v/s)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (v=v/s)";

      initialize();

      OBT band0 = blaze::band( tmat_, 0L );

      band0 = band0 / 0.5;

      checkSize    ( band0,  4UL );
      checkCapacity( band0,  4UL );
      checkNonZeros( band0,  2UL );
      checkRows    ( tmat_,  6UL );
      checkColumns ( tmat_,  4UL );
      checkNonZeros( tmat_, 10UL );

      if( band0[0] != -4 || band0[1] != 0 || band0[2] != -6 || band0[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( -4 0 -6 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != -4 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) !=  0 || tmat_(1,2) !=  1 || tmat_(1,3) !=  0 ||
          tmat_(2,0) !=  7 || tmat_(2,1) !=  4 || tmat_(2,2) != -6 || tmat_(2,3) !=  0 ||
          tmat_(3,0) !=  0 || tmat_(3,1) != -8 || tmat_(3,2) !=  5 || tmat_(3,3) !=  0 ||
          tmat_(4,0) !=  0 || tmat_(4,1) !=  0 || tmat_(4,2) !=  9 || tmat_(4,3) != -6 ||
          tmat_(5,0) !=  0 || tmat_(5,1) !=  0 || tmat_(5,2) !=  0 || tmat_(5,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( -4  0  0  0 )\n"
                                     "(  0  0  1  0 )\n"
                                     "(  7  4 -6  0 )\n"
                                     "(  0 -8  5  0 )\n"
                                     "(  0  0  9 -6 )\n"
                                     "(  0  0  0 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major Band::scale()
   //=====================================================================================

   {
      test_ = "Column-major Band::scale()";

      initialize();

      // Integral scaling the 1st lower band
      {
         OBT band1 = blaze::band( tmat_, -1L );
         band1.scale( 3 );

         checkSize    ( band1,  4UL );
         checkCapacity( band1,  4UL );
         checkNonZeros( band1,  3UL );
         checkRows    ( tmat_,  6UL );
         checkColumns ( tmat_,  4UL );
         checkNonZeros( tmat_, 10UL );

         if( band1[0] != 0 || band1[1] != 12 || band1[2] != 15 || band1[3] != -18 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Integral scale operation of 1st lower band failed\n"
                << " Details:\n"
                << "   Result:\n" << band1 << "\n"
                << "   Expected result:\n( 0 12 15 -18 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( tmat_(0,0) != -2 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=   0 ||
             tmat_(1,0) !=  0 || tmat_(1,1) !=  0 || tmat_(1,2) !=  1 || tmat_(1,3) !=   0 ||
             tmat_(2,0) !=  7 || tmat_(2,1) != 12 || tmat_(2,2) != -3 || tmat_(2,3) !=   0 ||
             tmat_(3,0) !=  0 || tmat_(3,1) != -8 || tmat_(3,2) != 15 || tmat_(3,3) !=   0 ||
             tmat_(4,0) !=  0 || tmat_(4,1) !=  0 || tmat_(4,2) !=  9 || tmat_(4,3) != -18 ||
             tmat_(5,0) !=  0 || tmat_(5,1) !=  0 || tmat_(5,2) !=  0 || tmat_(5,3) !=  10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Integral scale operation of 1st lower band failed\n"
                << " Details:\n"
                << "   Result:\n" << tmat_ << "\n"
                << "   Expected result:\n( -2  0  0   0 )\n"
                                        "(  0  0  1   0 )\n"
                                        "(  7 12 -3   0 )\n"
                                        "(  0 -8 15   0 )\n"
                                        "(  0  0  9 -18 )\n"
                                        "(  0  0  0  10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Floating point scaling the first lower band
      {
         OBT band1 = blaze::band( tmat_, -1L );
         band1.scale( 0.5 );

         checkSize    ( band1,  4UL );
         checkCapacity( band1,  4UL );
         checkNonZeros( band1,  3UL );
         checkRows    ( tmat_,  6UL );
         checkColumns ( tmat_,  4UL );
         checkNonZeros( tmat_, 10UL );

         if( band1[0] != 0 || band1[1] != 6 || band1[2] != 7 || band1[3] != -9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Floating point scale operation of 1st lower band failed\n"
                << " Details:\n"
                << "   Result:\n" << band1 << "\n"
                << "   Expected result:\n( 0 6 7 -9 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( tmat_(0,0) != -2 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
             tmat_(1,0) !=  0 || tmat_(1,1) !=  0 || tmat_(1,2) !=  1 || tmat_(1,3) !=  0 ||
             tmat_(2,0) !=  7 || tmat_(2,1) !=  6 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 ||
             tmat_(3,0) !=  0 || tmat_(3,1) != -8 || tmat_(3,2) !=  7 || tmat_(3,3) !=  0 ||
             tmat_(4,0) !=  0 || tmat_(4,1) !=  0 || tmat_(4,2) !=  9 || tmat_(4,3) != -9 ||
             tmat_(5,0) !=  0 || tmat_(5,1) !=  0 || tmat_(5,2) !=  0 || tmat_(5,3) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Integral scale operation of 1st lower band failed\n"
                << " Details:\n"
                << "   Result:\n" << tmat_ << "\n"
                << "   Expected result:\n( -2  0  0  0 )\n"
                                        "(  0  0  1  0 )\n"
                                        "(  7  6 -3  0 )\n"
                                        "(  0 -8  7  0 )\n"
                                        "(  0  0  9 -9 )\n"
                                        "(  0  0  0 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Band subscript operator.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of adding and accessing elements via the subscript operator
// of the Band specialization. In case an error is detected, a \a std::runtime_error exception
// is thrown.
*/
void DenseTest::testSubscript()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major Band::operator[]";

      initialize();

      BT band0 = blaze::band( mat_, 0L );

      // Assignment to the element at index 1
      band0[1] = 9;

      checkSize    ( band0, 4UL );
      checkCapacity( band0, 4UL );
      checkNonZeros( band0, 3UL );

      if( band0[0] != -2 || band0[1] != 9 || band0[2] != -3 || band0[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( -2 9 -3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != -2 || mat_(0,1) != 0 || mat_(0,2) !=  7 || mat_(0,3) !=  0 || mat_(0,4) !=  0 || mat_(0,5) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 9 || mat_(1,2) !=  4 || mat_(1,3) != -8 || mat_(1,4) !=  0 || mat_(1,5) !=  0 ||
          mat_(2,0) !=  0 || mat_(2,1) != 1 || mat_(2,2) != -3 || mat_(2,3) !=  5 || mat_(2,4) !=  9 || mat_(2,5) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) !=  0 || mat_(3,4) != -6 || mat_(3,5) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( -2  0  7  0  0  0 )\n"
                                     "(  0  9  4 -8  0  0 )\n"
                                     "(  0  1 -3  5  9  0 )\n"
                                     "(  0  0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Assignment to the element at index 2
      band0[2] = 0;

      checkSize    ( band0, 4UL );
      checkCapacity( band0, 4UL );
      checkNonZeros( band0, 2UL );

      if( band0[0] != -2 || band0[1] != 9 || band0[2] != 0 || band0[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( -2 9 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != -2 || mat_(0,1) != 0 || mat_(0,2) != 7 || mat_(0,3) !=  0 || mat_(0,4) !=  0 || mat_(0,5) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 9 || mat_(1,2) != 4 || mat_(1,3) != -8 || mat_(1,4) !=  0 || mat_(1,5) !=  0 ||
          mat_(2,0) !=  0 || mat_(2,1) != 1 || mat_(2,2) != 0 || mat_(2,3) !=  5 || mat_(2,4) !=  9 || mat_(2,5) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) != 0 || mat_(3,2) != 0 || mat_(3,3) !=  0 || mat_(3,4) != -6 || mat_(3,5) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( -2  0  7  0  0  0 )\n"
                                     "(  0  9  4 -8  0  0 )\n"
                                     "(  0  1  0  5  9  0 )\n"
                                     "(  0  0  0  0 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Assignment to the element at index 3
      band0[3] = -8;

      checkSize    ( band0, 4UL );
      checkCapacity( band0, 4UL );
      checkNonZeros( band0, 3UL );

      if( band0[0] != -2 || band0[1] != 9 || band0[2] != 0 || band0[3] != -8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( -2 9 0 -8 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != -2 || mat_(0,1) != 0 || mat_(0,2) != 7 || mat_(0,3) !=  0 || mat_(0,4) !=  0 || mat_(0,5) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 9 || mat_(1,2) != 4 || mat_(1,3) != -8 || mat_(1,4) !=  0 || mat_(1,5) !=  0 ||
          mat_(2,0) !=  0 || mat_(2,1) != 1 || mat_(2,2) != 0 || mat_(2,3) !=  5 || mat_(2,4) !=  9 || mat_(2,5) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) != 0 || mat_(3,2) != 0 || mat_(3,3) != -8 || mat_(3,4) != -6 || mat_(3,5) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( -2  0  7  0  0  0 )\n"
                                     "(  0  9  4 -8  0  0 )\n"
                                     "(  0  1  0  5  9  0 )\n"
                                     "(  0  0  0 -8 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Addition assignment to the element at index 0
      band0[0] += -3;

      checkSize    ( band0, 4UL );
      checkCapacity( band0, 4UL );
      checkNonZeros( band0, 3UL );

      if( band0[0] != -5 || band0[1] != 9 || band0[2] != 0 || band0[3] != -8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( -5 9 0 -8 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != -5 || mat_(0,1) != 0 || mat_(0,2) != 7 || mat_(0,3) !=  0 || mat_(0,4) !=  0 || mat_(0,5) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 9 || mat_(1,2) != 4 || mat_(1,3) != -8 || mat_(1,4) !=  0 || mat_(1,5) !=  0 ||
          mat_(2,0) !=  0 || mat_(2,1) != 1 || mat_(2,2) != 0 || mat_(2,3) !=  5 || mat_(2,4) !=  9 || mat_(2,5) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) != 0 || mat_(3,2) != 0 || mat_(3,3) != -8 || mat_(3,4) != -6 || mat_(3,5) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( -5  0  7  0  0  0 )\n"
                                     "(  0  9  4 -8  0  0 )\n"
                                     "(  0  1  0  5  9  0 )\n"
                                     "(  0  0  0 -8 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Subtraction assignment to the element at index 1
      band0[1] -= 6;

      checkSize    ( band0, 4UL );
      checkCapacity( band0, 4UL );
      checkNonZeros( band0, 3UL );

      if( band0[0] != -5 || band0[1] != 3 || band0[2] != 0 || band0[3] != -8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( -5 3 0 -8 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != -5 || mat_(0,1) != 0 || mat_(0,2) != 7 || mat_(0,3) !=  0 || mat_(0,4) !=  0 || mat_(0,5) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != 3 || mat_(1,2) != 4 || mat_(1,3) != -8 || mat_(1,4) !=  0 || mat_(1,5) !=  0 ||
          mat_(2,0) !=  0 || mat_(2,1) != 1 || mat_(2,2) != 0 || mat_(2,3) !=  5 || mat_(2,4) !=  9 || mat_(2,5) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) != 0 || mat_(3,2) != 0 || mat_(3,3) != -8 || mat_(3,4) != -6 || mat_(3,5) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( -5  0  7  0  0  0 )\n"
                                     "(  0  3  4 -8  0  0 )\n"
                                     "(  0  1  0  5  9  0 )\n"
                                     "(  0  0  0 -8 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Multiplication assignment to the element at index 1
      band0[1] *= -3;

      checkSize    ( band0, 4UL );
      checkCapacity( band0, 4UL );
      checkNonZeros( band0, 3UL );

      if( band0[0] != -5 || band0[1] != -9 || band0[2] != 0 || band0[3] != -8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( -5 -9 0 -8 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != -5 || mat_(0,1) !=  0 || mat_(0,2) != 7 || mat_(0,3) !=  0 || mat_(0,4) !=  0 || mat_(0,5) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != -9 || mat_(1,2) != 4 || mat_(1,3) != -8 || mat_(1,4) !=  0 || mat_(1,5) !=  0 ||
          mat_(2,0) !=  0 || mat_(2,1) !=  1 || mat_(2,2) != 0 || mat_(2,3) !=  5 || mat_(2,4) !=  9 || mat_(2,5) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) !=  0 || mat_(3,2) != 0 || mat_(3,3) != -8 || mat_(3,4) != -6 || mat_(3,5) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( -5  0  7  0  0  0 )\n"
                                     "(  0 -9  4 -8  0  0 )\n"
                                     "(  0  1  0  5  9  0 )\n"
                                     "(  0  0  0 -8 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Division assignment to the element at index 3
      band0[3] /= 2;

      checkSize    ( band0, 4UL );
      checkCapacity( band0, 4UL );
      checkNonZeros( band0, 3UL );

      if( band0[0] != -5 || band0[1] != -9 || band0[2] != 0 || band0[3] != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( -5 -9 0 -4 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != -5 || mat_(0,1) !=  0 || mat_(0,2) != 7 || mat_(0,3) !=  0 || mat_(0,4) !=  0 || mat_(0,5) !=  0 ||
          mat_(1,0) !=  0 || mat_(1,1) != -9 || mat_(1,2) != 4 || mat_(1,3) != -8 || mat_(1,4) !=  0 || mat_(1,5) !=  0 ||
          mat_(2,0) !=  0 || mat_(2,1) !=  1 || mat_(2,2) != 0 || mat_(2,3) !=  5 || mat_(2,4) !=  9 || mat_(2,5) !=  0 ||
          mat_(3,0) !=  0 || mat_(3,1) !=  0 || mat_(3,2) != 0 || mat_(3,3) != -4 || mat_(3,4) != -6 || mat_(3,5) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( -5  0  7  0  0  0 )\n"
                                     "(  0 -9  4 -8  0  0 )\n"
                                     "(  0  1  0  5  9  0 )\n"
                                     "(  0  0  0 -4 -6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major Band::operator[]";

      initialize();

      OBT band0 = blaze::band( tmat_, 0L );

      // Assignment to the element at index 1
      band0[1] = 9;

      checkSize    ( band0, 4UL );
      checkCapacity( band0, 4UL );
      checkNonZeros( band0, 3UL );

      if( band0[0] != -2 || band0[1] != 9 || band0[2] != -3 || band0[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( -2 9 -3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != -2 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) !=  9 || tmat_(1,2) !=  1 || tmat_(1,3) !=  0 ||
          tmat_(2,0) !=  7 || tmat_(2,1) !=  4 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 ||
          tmat_(3,0) !=  0 || tmat_(3,1) != -8 || tmat_(3,2) !=  5 || tmat_(3,3) !=  0 ||
          tmat_(4,0) !=  0 || tmat_(4,1) !=  0 || tmat_(4,2) !=  9 || tmat_(4,3) != -6 ||
          tmat_(5,0) !=  0 || tmat_(5,1) !=  0 || tmat_(5,2) !=  0 || tmat_(5,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( -2  0  0  0 )\n"
                                     "(  0  9  1  0 )\n"
                                     "(  7  4 -3  0 )\n"
                                     "(  0 -8  5  0 )\n"
                                     "(  0  0  9 -6 )\n"
                                     "(  0  0  0 10 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Assignment to the element at index 2
      band0[2] = 0;

      checkSize    ( band0, 4UL );
      checkCapacity( band0, 4UL );
      checkNonZeros( band0, 2UL );

      if( band0[0] != -2 || band0[1] != 9 || band0[2] != 0 || band0[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( -2 9 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != -2 || tmat_(0,1) !=  0 || tmat_(0,2) != 0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) !=  9 || tmat_(1,2) != 1 || tmat_(1,3) !=  0 ||
          tmat_(2,0) !=  7 || tmat_(2,1) !=  4 || tmat_(2,2) != 0 || tmat_(2,3) !=  0 ||
          tmat_(3,0) !=  0 || tmat_(3,1) != -8 || tmat_(3,2) != 5 || tmat_(3,3) !=  0 ||
          tmat_(4,0) !=  0 || tmat_(4,1) !=  0 || tmat_(4,2) != 9 || tmat_(4,3) != -6 ||
          tmat_(5,0) !=  0 || tmat_(5,1) !=  0 || tmat_(5,2) != 0 || tmat_(5,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( -2  0  0  0 )\n"
                                     "(  0  9  1  0 )\n"
                                     "(  7  4  0  0 )\n"
                                     "(  0 -8  5  0 )\n"
                                     "(  0  0  9 -6 )\n"
                                     "(  0  0  0 10 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Assignment to the element at index 3
      band0[3] = -8;

      checkSize    ( band0, 4UL );
      checkCapacity( band0, 4UL );
      checkNonZeros( band0, 3UL );

      if( band0[0] != -2 || band0[1] != 9 || band0[2] != 0 || band0[3] != -8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( -2 9 0 -8 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != -2 || tmat_(0,1) !=  0 || tmat_(0,2) != 0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) !=  9 || tmat_(1,2) != 1 || tmat_(1,3) !=  0 ||
          tmat_(2,0) !=  7 || tmat_(2,1) !=  4 || tmat_(2,2) != 0 || tmat_(2,3) !=  0 ||
          tmat_(3,0) !=  0 || tmat_(3,1) != -8 || tmat_(3,2) != 5 || tmat_(3,3) != -8 ||
          tmat_(4,0) !=  0 || tmat_(4,1) !=  0 || tmat_(4,2) != 9 || tmat_(4,3) != -6 ||
          tmat_(5,0) !=  0 || tmat_(5,1) !=  0 || tmat_(5,2) != 0 || tmat_(5,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( -2  0  0  0 )\n"
                                     "(  0  9  1  0 )\n"
                                     "(  7  4  0  0 )\n"
                                     "(  0 -8  5 -8 )\n"
                                     "(  0  0  9 -6 )\n"
                                     "(  0  0  0 10 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Addition assignment to the element at index 0
      band0[0] += -3;

      checkSize    ( band0, 4UL );
      checkCapacity( band0, 4UL );
      checkNonZeros( band0, 3UL );

      if( band0[0] != -5 || band0[1] != 9 || band0[2] != 0 || band0[3] != -8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( -5 9 0 -8 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != -5 || tmat_(0,1) !=  0 || tmat_(0,2) != 0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) !=  9 || tmat_(1,2) != 1 || tmat_(1,3) !=  0 ||
          tmat_(2,0) !=  7 || tmat_(2,1) !=  4 || tmat_(2,2) != 0 || tmat_(2,3) !=  0 ||
          tmat_(3,0) !=  0 || tmat_(3,1) != -8 || tmat_(3,2) != 5 || tmat_(3,3) != -8 ||
          tmat_(4,0) !=  0 || tmat_(4,1) !=  0 || tmat_(4,2) != 9 || tmat_(4,3) != -6 ||
          tmat_(5,0) !=  0 || tmat_(5,1) !=  0 || tmat_(5,2) != 0 || tmat_(5,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( -5  0  0  0 )\n"
                                     "(  0  9  1  0 )\n"
                                     "(  7  4  0  0 )\n"
                                     "(  0 -8  5 -8 )\n"
                                     "(  0  0  9 -6 )\n"
                                     "(  0  0  0 10 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Subtraction assignment to the element at index 1
      band0[1] -= 6;

      checkSize    ( band0, 4UL );
      checkCapacity( band0, 4UL );
      checkNonZeros( band0, 3UL );

      if( band0[0] != -5 || band0[1] != 3 || band0[2] != 0 || band0[3] != -8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( -5 3 0 -8 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != -5 || tmat_(0,1) !=  0 || tmat_(0,2) != 0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) !=  3 || tmat_(1,2) != 1 || tmat_(1,3) !=  0 ||
          tmat_(2,0) !=  7 || tmat_(2,1) !=  4 || tmat_(2,2) != 0 || tmat_(2,3) !=  0 ||
          tmat_(3,0) !=  0 || tmat_(3,1) != -8 || tmat_(3,2) != 5 || tmat_(3,3) != -8 ||
          tmat_(4,0) !=  0 || tmat_(4,1) !=  0 || tmat_(4,2) != 9 || tmat_(4,3) != -6 ||
          tmat_(5,0) !=  0 || tmat_(5,1) !=  0 || tmat_(5,2) != 0 || tmat_(5,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( -5  0  0  0 )\n"
                                     "(  0  3  1  0 )\n"
                                     "(  7  4  0  0 )\n"
                                     "(  0 -8  5 -8 )\n"
                                     "(  0  0  9 -6 )\n"
                                     "(  0  0  0 10 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Multiplication assignment to the element at index 1
      band0[1] *= -3;

      checkSize    ( band0, 4UL );
      checkCapacity( band0, 4UL );
      checkNonZeros( band0, 3UL );

      if( band0[0] != -5 || band0[1] != -9 || band0[2] != 0 || band0[3] != -8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( -5 -9 0 -8 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != -5 || tmat_(0,1) !=  0 || tmat_(0,2) != 0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) != -9 || tmat_(1,2) != 1 || tmat_(1,3) !=  0 ||
          tmat_(2,0) !=  7 || tmat_(2,1) !=  4 || tmat_(2,2) != 0 || tmat_(2,3) !=  0 ||
          tmat_(3,0) !=  0 || tmat_(3,1) != -8 || tmat_(3,2) != 5 || tmat_(3,3) != -8 ||
          tmat_(4,0) !=  0 || tmat_(4,1) !=  0 || tmat_(4,2) != 9 || tmat_(4,3) != -6 ||
          tmat_(5,0) !=  0 || tmat_(5,1) !=  0 || tmat_(5,2) != 0 || tmat_(5,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( -5  0  0  0 )\n"
                                     "(  0 -9  1  0 )\n"
                                     "(  7  4  0  0 )\n"
                                     "(  0 -8  5 -8 )\n"
                                     "(  0  0  9 -6 )\n"
                                     "(  0  0  0 10 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Division assignment to the element at index 3
      band0[3] /= 2;

      checkSize    ( band0, 4UL );
      checkCapacity( band0, 4UL );
      checkNonZeros( band0, 3UL );

      if( band0[0] != -5 || band0[1] != -9 || band0[2] != 0 || band0[3] != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << band0 << "\n"
             << "   Expected result:\n( -5 -9 0 -4 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != -5 || tmat_(0,1) !=  0 || tmat_(0,2) != 0 || tmat_(0,3) !=  0 ||
          tmat_(1,0) !=  0 || tmat_(1,1) != -9 || tmat_(1,2) != 1 || tmat_(1,3) !=  0 ||
          tmat_(2,0) !=  7 || tmat_(2,1) !=  4 || tmat_(2,2) != 0 || tmat_(2,3) !=  0 ||
          tmat_(3,0) !=  0 || tmat_(3,1) != -8 || tmat_(3,2) != 5 || tmat_(3,3) != -4 ||
          tmat_(4,0) !=  0 || tmat_(4,1) !=  0 || tmat_(4,2) != 9 || tmat_(4,3) != -6 ||
          tmat_(5,0) !=  0 || tmat_(5,1) !=  0 || tmat_(5,2) != 0 || tmat_(5,3) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( -5  0  0  0 )\n"
                                     "(  0 -9  1  0 )\n"
                                     "(  7  4  0  0 )\n"
                                     "(  0 -8  5 -4 )\n"
                                     "(  0  0  9 -6 )\n"
                                     "(  0  0  0 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Band iterator implementation.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the iterator implementation of the Band specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testIterator()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      initialize();

      // Testing the Iterator default constructor
      {
         test_ = "Row-major Iterator default constructor";

         BT::Iterator it{};

         if( it != BT::Iterator() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Failed iterator default constructor\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing the ConstIterator default constructor
      {
         test_ = "Row-major ConstIterator default constructor";

         BT::ConstIterator it{};

         if( it != BT::ConstIterator() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Failed iterator default constructor\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing conversion from Iterator to ConstIterator
      {
         test_ = "Row-major Iterator/ConstIterator conversion";

         BT band0 = blaze::band( mat_, 0L );
         BT::ConstIterator it( begin( band0 ) );

         if( it == end( band0 ) || *it != -2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Failed iterator conversion detected\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Counting the number of elements in 1st lower band via Iterator (end-begin)
      {
         test_ = "Row-major Iterator subtraction (end-begin)";

         BT band1 = blaze::band( mat_, -1L );
         const ptrdiff_t number( end( band1 ) - begin( band1 ) );

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

      // Counting the number of elements in 1st lower band via Iterator (begin-end)
      {
         test_ = "Row-major Iterator subtraction (begin-end)";

         BT band1 = blaze::band( mat_, -1L );
         const ptrdiff_t number( begin( band1 ) - end( band1 ) );

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

      // Counting the number of elements on the digaonal via ConstIterator (end-begin)
      {
         test_ = "Row-major ConstIterator subtraction (end-begin)";

         BT band0 = blaze::band( mat_, 0L );
         const ptrdiff_t number( cend( band0 ) - cbegin( band0 ) );

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

      // Counting the number of elements on the digaonal via ConstIterator (begin-end)
      {
         test_ = "Row-major ConstIterator subtraction (begin-end)";

         BT band0 = blaze::band( mat_, 0L );
         const ptrdiff_t number( cbegin( band0 ) - cend( band0 ) );

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
         test_ = "Row-major read-only access via ConstIterator";

         BT band1 = blaze::band( mat_, 1L );
         BT::ConstIterator it ( cbegin( band1 ) );
         BT::ConstIterator end( cend( band1 ) );

         if( it == end || *it != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid initial iterator detected\n";
            throw std::runtime_error( oss.str() );
         }

         ++it;

         if( it == end || *it != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator pre-increment failed\n";
            throw std::runtime_error( oss.str() );
         }

         --it;

         if( it == end || *it != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator pre-decrement failed\n";
            throw std::runtime_error( oss.str() );
         }

         it++;

         if( it == end || *it != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator post-increment failed\n";
            throw std::runtime_error( oss.str() );
         }

         it--;

         if( it == end || *it != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator post-decrement failed\n";
            throw std::runtime_error( oss.str() );
         }

         it += 2UL;

         if( it == end || *it != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator addition assignment failed\n";
            throw std::runtime_error( oss.str() );
         }

         it -= 2UL;

         if( it == end || *it != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator subtraction assignment failed\n";
            throw std::runtime_error( oss.str() );
         }

         it = it + 3UL;

         if( it == end || *it != -6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator/scalar addition failed\n";
            throw std::runtime_error( oss.str() );
         }

         it = it - 3UL;

         if( it == end || *it != 0 ) {
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
         test_ = "Row-major assignment via Iterator";

         BT band0 = blaze::band( mat_, 0L );
         int value = 6;

         for( BT::Iterator it=begin( band0 ); it!=end( band0 ); ++it ) {
            *it = value++;
         }

         if( band0[0] != 6 || band0[1] != 7 || band0[2] != 8 || band0[3] != 9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << band0 << "\n"
                << "   Expected result:\n( 6 7 8 9 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( mat_(0,0) != 6 || mat_(0,1) != 0 || mat_(0,2) != 7 || mat_(0,3) !=  0 || mat_(0,4) !=  0 || mat_(0,5) !=  0 ||
             mat_(1,0) != 0 || mat_(1,1) != 7 || mat_(1,2) != 4 || mat_(1,3) != -8 || mat_(1,4) !=  0 || mat_(1,5) !=  0 ||
             mat_(2,0) != 0 || mat_(2,1) != 1 || mat_(2,2) != 8 || mat_(2,3) !=  5 || mat_(2,4) !=  9 || mat_(2,5) !=  0 ||
             mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) != 0 || mat_(3,3) !=  9 || mat_(3,4) != -6 || mat_(3,5) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << mat_ << "\n"
                << "   Expected result:\n( 6  0  7  0  0  0 )\n"
                                        "( 0  7  4 -8  0  0 )\n"
                                        "( 0  1  8  5  9  0 )\n"
                                        "( 0  0  0  9 -6 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing addition assignment via Iterator
      {
         test_ = "Row-major addition assignment via Iterator";

         BT band0 = blaze::band( mat_, 0L );
         int value = 2;

         for( BT::Iterator it=begin( band0 ); it!=end( band0 ); ++it ) {
            *it += value++;
         }

         if( band0[0] != 8 || band0[1] != 10 || band0[2] != 12 || band0[3] != 14 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Addition assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << band0 << "\n"
                << "   Expected result:\n( 8 10 12 14 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( mat_(0,0) != 8 || mat_(0,1) !=  0 || mat_(0,2) !=  7 || mat_(0,3) !=  0 || mat_(0,4) !=  0 || mat_(0,5) !=  0 ||
             mat_(1,0) != 0 || mat_(1,1) != 10 || mat_(1,2) !=  4 || mat_(1,3) != -8 || mat_(1,4) !=  0 || mat_(1,5) !=  0 ||
             mat_(2,0) != 0 || mat_(2,1) !=  1 || mat_(2,2) != 12 || mat_(2,3) !=  5 || mat_(2,4) !=  9 || mat_(2,5) !=  0 ||
             mat_(3,0) != 0 || mat_(3,1) !=  0 || mat_(3,2) !=  0 || mat_(3,3) != 14 || mat_(3,4) != -6 || mat_(3,5) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << mat_ << "\n"
                << "   Expected result:\n( 8  0  7  0  0  0 )\n"
                                        "( 0 10  4 -8  0  0 )\n"
                                        "( 0  1 12  5  9  0 )\n"
                                        "( 0  0  0 14 -6 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing subtraction assignment via Iterator
      {
         test_ = "Row-major subtraction assignment via Iterator";

         BT band0 = blaze::band( mat_, 0L );
         int value = 2;

         for( BT::Iterator it=begin( band0 ); it!=end( band0 ); ++it ) {
            *it -= value++;
         }

         if( band0[0] != 6 || band0[1] != 7 || band0[2] != 8 || band0[3] != 9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subtraction assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << band0 << "\n"
                << "   Expected result:\n( 6 7 8 9 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( mat_(0,0) != 6 || mat_(0,1) != 0 || mat_(0,2) != 7 || mat_(0,3) !=  0 || mat_(0,4) !=  0 || mat_(0,5) !=  0 ||
             mat_(1,0) != 0 || mat_(1,1) != 7 || mat_(1,2) != 4 || mat_(1,3) != -8 || mat_(1,4) !=  0 || mat_(1,5) !=  0 ||
             mat_(2,0) != 0 || mat_(2,1) != 1 || mat_(2,2) != 8 || mat_(2,3) !=  5 || mat_(2,4) !=  9 || mat_(2,5) !=  0 ||
             mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) != 0 || mat_(3,3) !=  9 || mat_(3,4) != -6 || mat_(3,5) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << mat_ << "\n"
                << "   Expected result:\n( 6  0  7  0  0  0 )\n"
                                        "( 0  7  4 -8  0  0 )\n"
                                        "( 0  1  8  5  9  0 )\n"
                                        "( 0  0  0  9 -6 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing multiplication assignment via Iterator
      {
         test_ = "Row-major multiplication assignment via Iterator";

         BT band0 = blaze::band( mat_, 0L );
         int value = 1;

         for( BT::Iterator it=begin( band0 ); it!=end( band0 ); ++it ) {
            *it *= value++;
         }

         if( band0[0] != 6 || band0[1] != 14 || band0[2] != 24 || band0[3] != 36 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Multiplication assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << band0 << "\n"
                << "   Expected result:\n( 6 14 24 36 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( mat_(0,0) != 6 || mat_(0,1) !=  0 || mat_(0,2) !=  7 || mat_(0,3) !=  0 || mat_(0,4) !=  0 || mat_(0,5) !=  0 ||
             mat_(1,0) != 0 || mat_(1,1) != 14 || mat_(1,2) !=  4 || mat_(1,3) != -8 || mat_(1,4) !=  0 || mat_(1,5) !=  0 ||
             mat_(2,0) != 0 || mat_(2,1) !=  1 || mat_(2,2) != 24 || mat_(2,3) !=  5 || mat_(2,4) !=  9 || mat_(2,5) !=  0 ||
             mat_(3,0) != 0 || mat_(3,1) !=  0 || mat_(3,2) !=  0 || mat_(3,3) != 36 || mat_(3,4) != -6 || mat_(3,5) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << mat_ << "\n"
                << "   Expected result:\n( 6  0  7  0  0  0 )\n"
                                        "( 0 14  4 -8  0  0 )\n"
                                        "( 0  1 24  5  9  0 )\n"
                                        "( 0  0  0 36 -6 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing division assignment via Iterator
      {
         test_ = "Row-major division assignment via Iterator";

         BT band0 = blaze::band( mat_, 0L );

         for( BT::Iterator it=begin( band0 ); it!=end( band0 ); ++it ) {
            *it /= 2;
         }

         if( band0[0] != 3 || band0[1] != 7 || band0[2] != 12 || band0[3] != 18 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Division assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << band0 << "\n"
                << "   Expected result:\n( 3 7 12 18 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( mat_(0,0) != 3 || mat_(0,1) != 0 || mat_(0,2) !=  7 || mat_(0,3) !=  0 || mat_(0,4) !=  0 || mat_(0,5) !=  0 ||
             mat_(1,0) != 0 || mat_(1,1) != 7 || mat_(1,2) !=  4 || mat_(1,3) != -8 || mat_(1,4) !=  0 || mat_(1,5) !=  0 ||
             mat_(2,0) != 0 || mat_(2,1) != 1 || mat_(2,2) != 12 || mat_(2,3) !=  5 || mat_(2,4) !=  9 || mat_(2,5) !=  0 ||
             mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) != 18 || mat_(3,4) != -6 || mat_(3,5) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << mat_ << "\n"
                << "   Expected result:\n( 3  0  7  0  0  0 )\n"
                                        "( 0  7  4 -8  0  0 )\n"
                                        "( 0  1 12  5  9  0 )\n"
                                        "( 0  0  0 18 -6 10 )\n";
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

         OBT::Iterator it{};

         if( it != OBT::Iterator() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Failed iterator default constructor\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing the ConstIterator default constructor
      {
         test_ = "Column-major ConstIterator default constructor";

         OBT::ConstIterator it{};

         if( it != OBT::ConstIterator() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Failed iterator default constructor\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing conversion from Iterator to ConstIterator
      {
         test_ = "Column-major Iterator/ConstIterator conversion";

         OBT band0 = blaze::band( tmat_, 0L );
         OBT::ConstIterator it( begin( band0 ) );

         if( it == end( band0 ) || *it != -2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Failed iterator conversion detected\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Counting the number of elements in 1st upper band via Iterator (end-begin)
      {
         test_ = "Column-major Iterator subtraction (end-begin)";

         OBT band1 = blaze::band( tmat_, 1L );
         const ptrdiff_t number( end( band1 ) - begin( band1 ) );

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

      // Counting the number of elements in 1st upper band via Iterator (begin-end)
      {
         test_ = "Column-major Iterator subtraction (begin-end)";

         OBT band1 = blaze::band( tmat_, 1L );
         const ptrdiff_t number( begin( band1 ) - end( band1 ) );

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

      // Counting the number of elements on the diagonal via ConstIterator (end-begin)
      {
         test_ = "Column-major ConstIterator subtraction (end-begin)";

         OBT band0 = blaze::band( tmat_, 0L );
         const ptrdiff_t number( cend( band0 ) - cbegin( band0 ) );

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

      // Counting the number of elements on the diagonal via ConstIterator (begin-end)
      {
         test_ = "Column-major ConstIterator subtraction (begin-end)";

         OBT band0 = blaze::band( tmat_, 0L );
         const ptrdiff_t number( cbegin( band0 ) - cend( band0 ) );

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
         test_ = "Column-major read-only access via ConstIterator";

         OBT band1 = blaze::band( tmat_, -1L );
         OBT::ConstIterator it ( cbegin( band1 ) );
         OBT::ConstIterator end( cend( band1 ) );

         if( it == end || *it != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid initial iterator detected\n";
            throw std::runtime_error( oss.str() );
         }

         ++it;

         if( it == end || *it != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator pre-increment failed\n";
            throw std::runtime_error( oss.str() );
         }

         --it;

         if( it == end || *it != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator pre-decrement failed\n";
            throw std::runtime_error( oss.str() );
         }

         it++;

         if( it == end || *it != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator post-increment failed\n";
            throw std::runtime_error( oss.str() );
         }

         it--;

         if( it == end || *it != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator post-decrement failed\n";
            throw std::runtime_error( oss.str() );
         }

         it += 2UL;

         if( it == end || *it != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator addition assignment failed\n";
            throw std::runtime_error( oss.str() );
         }

         it -= 2UL;

         if( it == end || *it != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator subtraction assignment failed\n";
            throw std::runtime_error( oss.str() );
         }

         it = it + 3UL;

         if( it == end || *it != -6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator/scalar addition failed\n";
            throw std::runtime_error( oss.str() );
         }

         it = it - 3UL;

         if( it == end || *it != 0 ) {
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
         test_ = "Column-major assignment via Iterator";

         OBT band0 = blaze::band( tmat_, 0L );
         int value = 6;

         for( OBT::Iterator it=begin( band0 ); it!=end( band0 ); ++it ) {
            *it = value++;
         }

         if( band0[0] != 6 || band0[1] != 7 || band0[2] != 8 || band0[3] != 9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << band0 << "\n"
                << "   Expected result:\n( 6 7 8 9 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( tmat_(0,0) != 6 || tmat_(0,1) !=  0 || tmat_(0,2) != 0 || tmat_(0,3) !=  0 ||
             tmat_(1,0) != 0 || tmat_(1,1) !=  7 || tmat_(1,2) != 1 || tmat_(1,3) !=  0 ||
             tmat_(2,0) != 7 || tmat_(2,1) !=  4 || tmat_(2,2) != 8 || tmat_(2,3) !=  0 ||
             tmat_(3,0) != 0 || tmat_(3,1) != -8 || tmat_(3,2) != 5 || tmat_(3,3) !=  9 ||
             tmat_(4,0) != 0 || tmat_(4,1) !=  0 || tmat_(4,2) != 9 || tmat_(4,3) != -6 ||
             tmat_(5,0) != 0 || tmat_(5,1) !=  0 || tmat_(5,2) != 0 || tmat_(5,3) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << tmat_ << "\n"
                << "   Expected result:\n( 6  0  0  0 )\n"
                                        "( 0  7  1  0 )\n"
                                        "( 7  4  8  0 )\n"
                                        "( 0 -8  5  9 )\n"
                                        "( 0  0  9 -6 )\n"
                                        "( 0  0  0 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing addition assignment via Iterator
      {
         test_ = "Column-major addition assignment via Iterator";

         OBT band0 = blaze::band( tmat_, 0L );
         int value = 2;

         for( OBT::Iterator it=begin( band0 ); it!=end( band0 ); ++it ) {
            *it += value++;
         }

         if( band0[0] != 8 || band0[1] != 10 || band0[2] != 12 || band0[3] != 14 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Addition assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << band0 << "\n"
                << "   Expected result:\n( 8 10 12 14 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( tmat_(0,0) != 8 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
             tmat_(1,0) != 0 || tmat_(1,1) != 10 || tmat_(1,2) !=  1 || tmat_(1,3) !=  0 ||
             tmat_(2,0) != 7 || tmat_(2,1) !=  4 || tmat_(2,2) != 12 || tmat_(2,3) !=  0 ||
             tmat_(3,0) != 0 || tmat_(3,1) != -8 || tmat_(3,2) !=  5 || tmat_(3,3) != 14 ||
             tmat_(4,0) != 0 || tmat_(4,1) !=  0 || tmat_(4,2) !=  9 || tmat_(4,3) != -6 ||
             tmat_(5,0) != 0 || tmat_(5,1) !=  0 || tmat_(5,2) !=  0 || tmat_(5,3) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << tmat_ << "\n"
                << "   Expected result:\n( 8  0  0  0 )\n"
                                        "( 0 10  1  0 )\n"
                                        "( 7  4 12  0 )\n"
                                        "( 0 -8  5 14 )\n"
                                        "( 0  0  9 -6 )\n"
                                        "( 0  0  0 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing subtraction assignment via Iterator
      {
         test_ = "Column-major subtraction assignment via Iterator";

         OBT band0 = blaze::band( tmat_, 0L );
         int value = 2;

         for( OBT::Iterator it=begin( band0 ); it!=end( band0 ); ++it ) {
            *it -= value++;
         }

         if( band0[0] != 6 || band0[1] != 7 || band0[2] != 8 || band0[3] != 9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subtraction assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << band0 << "\n"
                << "   Expected result:\n( 6 7 8 9 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( tmat_(0,0) != 6 || tmat_(0,1) !=  0 || tmat_(0,2) != 0 || tmat_(0,3) !=  0 ||
             tmat_(1,0) != 0 || tmat_(1,1) !=  7 || tmat_(1,2) != 1 || tmat_(1,3) !=  0 ||
             tmat_(2,0) != 7 || tmat_(2,1) !=  4 || tmat_(2,2) != 8 || tmat_(2,3) !=  0 ||
             tmat_(3,0) != 0 || tmat_(3,1) != -8 || tmat_(3,2) != 5 || tmat_(3,3) !=  9 ||
             tmat_(4,0) != 0 || tmat_(4,1) !=  0 || tmat_(4,2) != 9 || tmat_(4,3) != -6 ||
             tmat_(5,0) != 0 || tmat_(5,1) !=  0 || tmat_(5,2) != 0 || tmat_(5,3) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << tmat_ << "\n"
                << "   Expected result:\n( 6  0  0  0 )\n"
                                        "( 0  7  1  0 )\n"
                                        "( 7  4  8  0 )\n"
                                        "( 0 -8  5  9 )\n"
                                        "( 0  0  9 -6 )\n"
                                        "( 0  0  0 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing multiplication assignment via Iterator
      {
         test_ = "Column-major multiplication assignment via Iterator";

         OBT band0 = blaze::band( tmat_, 0L );
         int value = 1;

         for( OBT::Iterator it=begin( band0 ); it!=end( band0 ); ++it ) {
            *it *= value++;
         }

         if( band0[0] != 6 || band0[1] != 14 || band0[2] != 24 || band0[3] != 36 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Multiplication assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << band0 << "\n"
                << "   Expected result:\n( 6 14 24 36 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( tmat_(0,0) != 6 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
             tmat_(1,0) != 0 || tmat_(1,1) != 14 || tmat_(1,2) !=  1 || tmat_(1,3) !=  0 ||
             tmat_(2,0) != 7 || tmat_(2,1) !=  4 || tmat_(2,2) != 24 || tmat_(2,3) !=  0 ||
             tmat_(3,0) != 0 || tmat_(3,1) != -8 || tmat_(3,2) !=  5 || tmat_(3,3) != 36 ||
             tmat_(4,0) != 0 || tmat_(4,1) !=  0 || tmat_(4,2) !=  9 || tmat_(4,3) != -6 ||
             tmat_(5,0) != 0 || tmat_(5,1) !=  0 || tmat_(5,2) !=  0 || tmat_(5,3) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << tmat_ << "\n"
                << "   Expected result:\n( 6  0  0  0 )\n"
                                        "( 0 14  1  0 )\n"
                                        "( 7  4 24  0 )\n"
                                        "( 0 -8  5 36 )\n"
                                        "( 0  0  9 -6 )\n"
                                        "( 0  0  0 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing division assignment via Iterator
      {
         test_ = "Column-major division assignment via Iterator";

         OBT band0 = blaze::band( tmat_, 0L );

         for( OBT::Iterator it=begin( band0 ); it!=end( band0 ); ++it ) {
            *it /= 2;
         }

         if( band0[0] != 3 || band0[1] != 7 || band0[2] != 12 || band0[3] != 18 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Division assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << band0 << "\n"
                << "   Expected result:\n( 3 7 12 18 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( tmat_(0,0) != 3 || tmat_(0,1) !=  0 || tmat_(0,2) !=  0 || tmat_(0,3) !=  0 ||
             tmat_(1,0) != 0 || tmat_(1,1) !=  7 || tmat_(1,2) !=  1 || tmat_(1,3) !=  0 ||
             tmat_(2,0) != 7 || tmat_(2,1) !=  4 || tmat_(2,2) != 12 || tmat_(2,3) !=  0 ||
             tmat_(3,0) != 0 || tmat_(3,1) != -8 || tmat_(3,2) !=  5 || tmat_(3,3) != 18 ||
             tmat_(4,0) != 0 || tmat_(4,1) !=  0 || tmat_(4,2) !=  9 || tmat_(4,3) != -6 ||
             tmat_(5,0) != 0 || tmat_(5,1) !=  0 || tmat_(5,2) !=  0 || tmat_(5,3) != 10 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << tmat_ << "\n"
                << "   Expected result:\n( 3  0  0  0 )\n"
                                        "( 0  7  1  0 )\n"
                                        "( 7  4 12  0 )\n"
                                        "( 0 -8  5 18 )\n"
                                        "( 0  0  9 -6 )\n"
                                        "( 0  0  0 10 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c nonZeros() member function of the Band specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c nonZeros() member function of the Band specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testNonZeros()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major Band::nonZeros()";

      initialize();

      // Initialization check
      BT band1 = blaze::band( mat_, 1L );

      checkSize    ( band1, 4UL );
      checkCapacity( band1, 4UL );
      checkNonZeros( band1, 3UL );

      if( band1[0] != 0 || band1[1] != 4 || band1[2] != 5 || band1[3] != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << band1 << "\n"
             << "   Expected result:\n( 0 4 5 -6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Changing the number of non-zeros via the dense band
      band1[2] = 0;

      checkSize    ( band1, 4UL );
      checkCapacity( band1, 4UL );
      checkNonZeros( band1, 2UL );

      if( band1[0] != 0 || band1[1] != 4 || band1[2] != 0 || band1[3] != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << band1 << "\n"
             << "   Expected result:\n( 0 4 0 -6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Changing the number of non-zeros via the dense matrix
      mat_(0,1) = 5;

      checkSize    ( band1, 4UL );
      checkCapacity( band1, 4UL );
      checkNonZeros( band1, 3UL );

      if( band1[0] != 5 || band1[1] != 4 || band1[2] != 0 || band1[3] != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << band1 << "\n"
             << "   Expected result:\n( 5 4 0 -6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major Band::nonZeros()";

      initialize();

      // Initialization check
      OBT band1 = blaze::band( tmat_, -1L );

      checkSize    ( band1, 4UL );
      checkCapacity( band1, 4UL );
      checkNonZeros( band1, 3UL );

      if( band1[0] != 0 || band1[1] != 4 || band1[2] != 5 || band1[3] != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << band1 << "\n"
             << "   Expected result:\n( 0 4 5 -6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Changing the number of non-zeros via the dense band
      band1[2] = 0;

      checkSize    ( band1, 4UL );
      checkCapacity( band1, 4UL );
      checkNonZeros( band1, 2UL );

      if( band1[0] != 0 || band1[1] != 4 || band1[2] != 0 || band1[3] != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << band1 << "\n"
             << "   Expected result:\n( 0 4 0 -6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Changing the number of non-zeros via the dense matrix
      tmat_(1,0) = 5;

      checkSize    ( band1, 4UL );
      checkCapacity( band1, 4UL );
      checkNonZeros( band1, 3UL );

      if( band1[0] != 5 || band1[1] != 4 || band1[2] != 0 || band1[3] != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << band1 << "\n"
             << "   Expected result:\n( 5 4 0 -6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c reset() member function of the Band specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c reset() member function of the Band specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testReset()
{
   using blaze::reset;


   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major Band::reset()";

      // Resetting a single element in the 1st upper band
      {
         initialize();

         BT band1 = blaze::band( mat_, 1L );
         reset( band1[1] );

         checkSize    ( band1, 4UL );
         checkCapacity( band1, 4UL );
         checkNonZeros( band1, 2UL );
         checkRows    ( mat_ , 4UL );
         checkColumns ( mat_ , 6UL );
         checkNonZeros( mat_ , 9UL );

         if( band1[0] != 0 || band1[1] != 0 || band1[2] != 5 || band1[3] != -6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Reset operation failed\n"
                << " Details:\n"
                << "   Result:\n" << band1 << "\n"
                << "   Expected result:\n( 0 0 5 -6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Resetting the 1st upper band (lvalue)
      {
         initialize();

         BT band1 = blaze::band( mat_, 1L );
         reset( band1 );

         checkSize    ( band1, 4UL );
         checkCapacity( band1, 4UL );
         checkNonZeros( band1, 0UL );
         checkRows    ( mat_ , 4UL );
         checkColumns ( mat_ , 6UL );
         checkNonZeros( mat_ , 7UL );

         if( band1[0] != 0 || band1[1] != 0 || band1[2] != 0 || band1[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Reset operation of 1st upper band failed\n"
                << " Details:\n"
                << "   Result:\n" << band1 << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Resetting the 2nd upper band (rvalue)
      {
         initialize();

         reset( blaze::band( mat_, 2L ) );

         checkRows    ( mat_, 4UL );
         checkColumns ( mat_, 6UL );
         checkNonZeros( mat_, 6UL );

         if( mat_(0,2) != 0 || mat_(1,3) != 0 || mat_(2,4) != 0 || mat_(3,5) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Reset operation of 2nd upper band failed\n"
                << " Details:\n"
                << "   Result:\n" << mat_ << "\n"
                << "   Expected result:\n( -2  0  0  0  0  0 )\n"
                                        "(  0  0  4  0  0  0 )\n"
                                        "(  0  1 -3  5  0  0 )\n"
                                        "(  0  0  0  0 -6  0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major Band::reset()";

      // Resetting a single element in the 1st lower band
      {
         initialize();

         OBT band1 = blaze::band( tmat_, -1L );
         reset( band1[1] );

         checkSize    ( band1, 4UL );
         checkCapacity( band1, 4UL );
         checkNonZeros( band1, 2UL );
         checkRows    ( tmat_, 6UL );
         checkColumns ( tmat_, 4UL );
         checkNonZeros( tmat_, 9UL );

         if( band1[0] != 0 || band1[1] != 0 || band1[2] != 5 || band1[3] != -6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Reset operation failed\n"
                << " Details:\n"
                << "   Result:\n" << band1 << "\n"
                << "   Expected result:\n( 0 0 5 -6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Resetting the 1st lower band (lvalue)
      {
         initialize();

         OBT band1 = blaze::band( tmat_, -1L );
         reset( band1 );

         checkSize    ( band1, 4UL );
         checkCapacity( band1, 4UL );
         checkNonZeros( band1, 0UL );
         checkRows    ( tmat_, 6UL );
         checkColumns ( tmat_, 4UL );
         checkNonZeros( tmat_, 7UL );

         if( band1[0] != 0 || band1[1] != 0 || band1[2] != 0 || band1[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Reset operation of 1st lower band failed\n"
                << " Details:\n"
                << "   Result:\n" << band1 << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Resetting the 2nd lower band (rvalue)
      {
         initialize();

         reset( blaze::band( tmat_, -2L ) );

         checkRows    ( tmat_, 6UL );
         checkColumns ( tmat_, 4UL );
         checkNonZeros( tmat_, 6UL );

         if( tmat_(2,0) != 0 || tmat_(3,1) != 0 || tmat_(4,2) != 0 || tmat_(5,3) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Reset operation of 2nd lower band failed\n"
                << " Details:\n"
                << "   Result:\n" << tmat_ << "\n"
                << "   Expected result:\n( -2  0  0  0 )\n"
                                        "(  0  0  1  0 )\n"
                                        "(  0  4 -3  0 )\n"
                                        "(  0  0  5  0 )\n"
                                        "(  0  0  0 -6 )\n"
                                        "(  0  0  0  0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c clear() function with the Band specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c clear() function with the Band specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testClear()
{
   using blaze::clear;


   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major clear() function";

      // Clearing a single element in the 1st upper band
      {
         initialize();

         BT band1 = blaze::band( mat_, 1L );
         clear( band1[1] );

         checkSize    ( band1, 4UL );
         checkCapacity( band1, 4UL );
         checkNonZeros( band1, 2UL );
         checkRows    ( mat_ , 4UL );
         checkColumns ( mat_ , 6UL );
         checkNonZeros( mat_ , 9UL );

         if( band1[0] != 0 || band1[1] != 0 || band1[2] != 5 || band1[3] != -6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Clear operation failed\n"
                << " Details:\n"
                << "   Result:\n" << band1 << "\n"
                << "   Expected result:\n( 0 0 5 -6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Clearing the 1st upper band (lvalue)
      {
         initialize();

         BT band1 = blaze::band( mat_, 1L );
         clear( band1 );

         checkSize    ( band1, 4UL );
         checkCapacity( band1, 4UL );
         checkNonZeros( band1, 0UL );
         checkRows    ( mat_ , 4UL );
         checkColumns ( mat_ , 6UL );
         checkNonZeros( mat_ , 7UL );

         if( band1[0] != 0 || band1[1] != 0 || band1[2] != 0 || band1[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Clear operation of 1st upper band failed\n"
                << " Details:\n"
                << "   Result:\n" << band1 << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Clearing the 2nd upper band (rvalue)
      {
         initialize();

         clear( blaze::band( mat_, 2L ) );

         checkRows    ( mat_, 4UL );
         checkColumns ( mat_, 6UL );
         checkNonZeros( mat_, 6UL );

         if( mat_(0,2) != 0 || mat_(1,3) != 0 || mat_(2,4) != 0 || mat_(3,5) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Clear operation of 2nd upper band failed\n"
                << " Details:\n"
                << "   Result:\n" << mat_ << "\n"
                << "   Expected result:\n( -2  0  0  0  0  0 )\n"
                                        "(  0  0  4  0  0  0 )\n"
                                        "(  0  1 -3  5  0  0 )\n"
                                        "(  0  0  0  0 -6  0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major clear() function";

      // Clearing a single element in the 1st lower band
      {
         initialize();

         OBT band1 = blaze::band( tmat_, -1L );
         clear( band1[1] );

         checkSize    ( band1, 4UL );
         checkCapacity( band1, 4UL );
         checkNonZeros( band1, 2UL );
         checkRows    ( tmat_, 6UL );
         checkColumns ( tmat_, 4UL );
         checkNonZeros( tmat_, 9UL );

         if( band1[0] != 0 || band1[1] != 0 || band1[2] != 5 || band1[3] != -6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Clear operation failed\n"
                << " Details:\n"
                << "   Result:\n" << band1 << "\n"
                << "   Expected result:\n( 0 0 5 -6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Clearing the 1st lower band (lvalue)
      {
         initialize();

         OBT band1 = blaze::band( tmat_, -1L );
         clear( band1 );

         checkSize    ( band1, 4UL );
         checkCapacity( band1, 4UL );
         checkNonZeros( band1, 0UL );
         checkRows    ( tmat_, 6UL );
         checkColumns ( tmat_, 4UL );
         checkNonZeros( tmat_, 7UL );

         if( band1[0] != 0 || band1[1] != 0 || band1[2] != 0 || band1[3] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Clear operation of 1st lower band failed\n"
                << " Details:\n"
                << "   Result:\n" << band1 << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Clearing the 2nd lower band (rvalue)
      {
         initialize();

         clear( blaze::band( tmat_, -2L ) );

         checkRows    ( tmat_, 6UL );
         checkColumns ( tmat_, 4UL );
         checkNonZeros( tmat_, 6UL );

         if( tmat_(2,0) != 0 || tmat_(3,1) != 0 || tmat_(4,2) != 0 || tmat_(5,3) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Clear operation of 2nd lower band failed\n"
                << " Details:\n"
                << "   Result:\n" << tmat_ << "\n"
                << "   Expected result:\n( -2  0  0  0 )\n"
                                        "(  0  0  1  0 )\n"
                                        "(  0  4 -3  0 )\n"
                                        "(  0  0  5  0 )\n"
                                        "(  0  0  0 -6 )\n"
                                        "(  0  0  0  0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c isDefault() function with the Band specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isDefault() function with the Band specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testIsDefault()
{
   using blaze::isDefault;


   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major isDefault() function";

      initialize();

      // isDefault with default band
      {
         BT band3 = blaze::band( mat_, 3L );

         if( isDefault( band3[1] ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Band element: " << band3[1] << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( band3 ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Band:\n" << band3 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with non-default band
      {
         BT band2 = blaze::band( mat_, 2L );

         if( isDefault( band2[1] ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Band element: " << band2[1] << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( band2 ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Band:\n" << band2 << "\n";
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

      // isDefault with default band
      {
         OBT band3 = blaze::band( tmat_, -3L );

         if( isDefault( band3[1] ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Band element: " << band3[1] << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( band3 ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Band:\n" << band3 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with non-default band
      {
         OBT band2 = blaze::band( tmat_, -2L );

         if( isDefault( band2[1] ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Band element: " << band2[1] << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( band2 ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Band:\n" << band2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c isSame() function with the Band specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isSame() function with the Band specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testIsSame()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major isSame() function";

      // isSame with matching bands
      {
         BT band1 = blaze::band( mat_, 1L );
         BT band2 = blaze::band( mat_, 1L );

         if( blaze::isSame( band1, band2 ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First band:\n" << band1 << "\n"
                << "   Second band:\n" << band2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with non-matching bands
      {
         BT band1 = blaze::band( mat_, 0L );
         BT band2 = blaze::band( mat_, 1L );

         if( blaze::isSame( band1, band2 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First band:\n" << band1 << "\n"
                << "   Second band:\n" << band2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with band and matching subvector
      {
         BT   band1 = blaze::band( mat_, 1L );
         auto sv    = blaze::subvector( band1, 0UL, 4UL );

         if( blaze::isSame( band1, sv ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Dense band:\n" << band1 << "\n"
                << "   Dense subvector:\n" << sv << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( sv, band1 ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Dense band:\n" << band1 << "\n"
                << "   Dense subvector:\n" << sv << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with band and non-matching subvector (different size)
      {
         BT   band1 = blaze::band( mat_, 1L );
         auto sv    = blaze::subvector( band1, 0UL, 3UL );

         if( blaze::isSame( band1, sv ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Dense band:\n" << band1 << "\n"
                << "   Dense subvector:\n" << sv << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( sv, band1 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Dense band:\n" << band1 << "\n"
                << "   Dense subvector:\n" << sv << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with band and non-matching subvector (different offset)
      {
         BT   band1 = blaze::band( mat_, 1L );
         auto sv    = blaze::subvector( band1, 1UL, 3UL );

         if( blaze::isSame( band1, sv ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Dense band:\n" << band1 << "\n"
                << "   Dense subvector:\n" << sv << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( sv, band1 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Dense band:\n" << band1 << "\n"
                << "   Dense subvector:\n" << sv << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with matching bands on a common submatrix
      {
         auto sm    = blaze::submatrix( mat_, 1UL, 1UL, 3UL, 4UL );
         auto band1 = blaze::band( sm, 1L );
         auto band2 = blaze::band( sm, 1L );

         if( blaze::isSame( band1, band2 ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First band:\n" << band1 << "\n"
                << "   Second band:\n" << band2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with non-matching bands on a common submatrix
      {
         auto sm    = blaze::submatrix( mat_, 1UL, 1UL, 3UL, 4UL );
         auto band1 = blaze::band( sm, 0L );
         auto band2 = blaze::band( sm, 1L );

         if( blaze::isSame( band1, band2 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First band:\n" << band1 << "\n"
                << "   Second band:\n" << band2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with matching rows on matrix and submatrix
      {
         auto sm    = blaze::submatrix( mat_, 0UL, 1UL, 4UL, 4UL );
         auto band1 = blaze::band( mat_, 1L );
         auto band2 = blaze::band( sm  , 0L );

         if( blaze::isSame( band1, band2 ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First band:\n" << band1 << "\n"
                << "   Second band:\n" << band2 << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( band2, band1 ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First band:\n" << band1 << "\n"
                << "   Second band:\n" << band2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with non-matching rows on matrix and submatrix (different band)
      {
         auto sm    = blaze::submatrix( mat_, 0UL, 1UL, 4UL, 4UL );
         auto band1 = blaze::band( mat_, 2L );
         auto band2 = blaze::band( sm  , 0L );

         if( blaze::isSame( band1, band2 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First band:\n" << band1 << "\n"
                << "   Second band:\n" << band2 << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( band2, band1 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First band:\n" << band1 << "\n"
                << "   Second band:\n" << band2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with non-matching rows on matrix and submatrix (different size)
      {
         auto sm    = blaze::submatrix( mat_, 0UL, 1UL, 3UL, 4UL );
         auto band1 = blaze::band( mat_, 1L );
         auto band2 = blaze::band( sm  , 0L );

         if( blaze::isSame( band1, band2 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First band:\n" << band1 << "\n"
                << "   Second band:\n" << band2 << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( band2, band1 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First band:\n" << band1 << "\n"
                << "   Second band:\n" << band2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with matching rows on two submatrices
      {
         auto sm1   = blaze::submatrix( mat_, 1UL, 0UL, 3UL, 4UL );
         auto sm2   = blaze::submatrix( mat_, 1UL, 1UL, 3UL, 5UL );
         auto band1 = blaze::band( sm1, 1L );
         auto band2 = blaze::band( sm2, 0L );

         if( blaze::isSame( band1, band2 ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First band:\n" << band1 << "\n"
                << "   Second band:\n" << band2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with non-matching rows on two submatrices (different band)
      {
         auto sm1   = blaze::submatrix( mat_, 1UL, 0UL, 3UL, 4UL );
         auto sm2   = blaze::submatrix( mat_, 1UL, 1UL, 3UL, 5UL );
         auto band1 = blaze::band( sm1, 1L );
         auto band2 = blaze::band( sm2, 1L );

         if( blaze::isSame( band1, band2 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First band:\n" << band1 << "\n"
                << "   Second band:\n" << band2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with non-matching rows on two submatrices (different size)
      {
         auto sm1   = blaze::submatrix( mat_, 1UL, 0UL, 3UL, 4UL );
         auto sm2   = blaze::submatrix( mat_, 1UL, 1UL, 2UL, 5UL );
         auto band1 = blaze::band( sm1, 1L );
         auto band2 = blaze::band( sm2, 0L );

         if( blaze::isSame( band1, band2 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First band:\n" << band1 << "\n"
                << "   Second band:\n" << band2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with non-matching rows on two submatrices (different offset)
      {
         auto sm1   = blaze::submatrix( mat_, 0UL, 0UL, 3UL, 4UL );
         auto sm2   = blaze::submatrix( mat_, 1UL, 1UL, 2UL, 5UL );
         auto band1 = blaze::band( sm1, 0L );
         auto band2 = blaze::band( sm2, 0L );

         if( blaze::isSame( band1, band2 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First band:\n" << band1 << "\n"
                << "   Second band:\n" << band2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with matching band subvectors on submatrices
      {
         auto sm    = blaze::submatrix( mat_, 1UL, 1UL, 3UL, 4UL );
         auto band1 = blaze::band( sm, 1L );
         auto sv1   = blaze::subvector( band1, 0UL, 2UL );
         auto sv2   = blaze::subvector( band1, 0UL, 2UL );

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

      // isSame with non-matching band subvectors on submatrices (different size)
      {
         auto sm    = blaze::submatrix( mat_, 1UL, 1UL, 3UL, 4UL );
         auto band1 = blaze::band( sm, 1L );
         auto sv1   = blaze::subvector( band1, 0UL, 3UL );
         auto sv2   = blaze::subvector( band1, 0UL, 2UL );

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

      // isSame with non-matching band subvectors on submatrices (different offset)
      {
         auto sm    = blaze::submatrix( mat_, 1UL, 1UL, 3UL, 4UL );
         auto band1 = blaze::band( sm, 1L );
         auto sv1   = blaze::subvector( band1, 0UL, 2UL );
         auto sv2   = blaze::subvector( band1, 1UL, 2UL );

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

      // isSame with matching row subvectors on two submatrices
      {
         auto sm1   = blaze::submatrix( mat_, 1UL, 0UL, 3UL, 4UL );
         auto sm2   = blaze::submatrix( mat_, 1UL, 1UL, 3UL, 5UL );
         auto band1 = blaze::band( sm1, 1L );
         auto band2 = blaze::band( sm2, 0L );
         auto sv1   = blaze::subvector( band1, 0UL, 2UL );
         auto sv2   = blaze::subvector( band2, 0UL, 2UL );

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

      // isSame with non-matching row subvectors on two submatrices (different size)
      {
         auto sm1   = blaze::submatrix( mat_, 1UL, 0UL, 3UL, 4UL );
         auto sm2   = blaze::submatrix( mat_, 1UL, 1UL, 3UL, 5UL );
         auto band1 = blaze::band( sm1, 1L );
         auto band2 = blaze::band( sm2, 0L );
         auto sv1   = blaze::subvector( band1, 0UL, 3UL );
         auto sv2   = blaze::subvector( band2, 0UL, 2UL );

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

      // isSame with non-matching row subvectors on two submatrices (different offset)
      {
         auto sm1   = blaze::submatrix( mat_, 1UL, 0UL, 3UL, 4UL );
         auto sm2   = blaze::submatrix( mat_, 1UL, 1UL, 3UL, 5UL );
         auto band1 = blaze::band( sm1, 1L );
         auto band2 = blaze::band( sm2, 0L );
         auto sv1   = blaze::subvector( band1, 0UL, 2UL );
         auto sv2   = blaze::subvector( band2, 1UL, 2UL );

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

      // isSame with matching bands
      {
         OBT band1 = blaze::band( tmat_, -1L );
         OBT band2 = blaze::band( tmat_, -1L );

         if( blaze::isSame( band1, band2 ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First band:\n" << band1 << "\n"
                << "   Second band:\n" << band2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with non-matching bands
      {
         OBT band1 = blaze::band( tmat_, -1L );
         OBT band2 = blaze::band( tmat_,  0L );

         if( blaze::isSame( band1, band2 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First band:\n" << band1 << "\n"
                << "   Second band:\n" << band2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with band and matching subvector
      {
         OBT  band1 = blaze::band( tmat_, -1L );
         auto sv    = blaze::subvector( band1, 0UL, 4UL );

         if( blaze::isSame( band1, sv ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Dense band:\n" << band1 << "\n"
                << "   Dense subvector:\n" << sv << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( sv, band1 ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Dense band:\n" << band1 << "\n"
                << "   Dense subvector:\n" << sv << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with band and non-matching subvector (different size)
      {
         OBT  band1 = blaze::band( tmat_, -1L );
         auto sv    = blaze::subvector( band1, 0UL, 3UL );

         if( blaze::isSame( band1, sv ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Dense band:\n" << band1 << "\n"
                << "   Dense subvector:\n" << sv << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( sv, band1 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Dense band:\n" << band1 << "\n"
                << "   Dense subvector:\n" << sv << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with band and non-matching subvector (different offset)
      {
         OBT  band1 = blaze::band( tmat_, -1L );
         auto sv    = blaze::subvector( band1, 1UL, 3UL );

         if( blaze::isSame( band1, sv ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Dense band:\n" << band1 << "\n"
                << "   Dense subvector:\n" << sv << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( sv, band1 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Dense band:\n" << band1 << "\n"
                << "   Dense subvector:\n" << sv << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with matching bands on a common submatrix
      {
         auto sm    = blaze::submatrix( tmat_, 1UL, 1UL, 4UL, 3UL );
         auto band1 = blaze::band( sm, -1L );
         auto band2 = blaze::band( sm, -1L );

         if( blaze::isSame( band1, band2 ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First band:\n" << band1 << "\n"
                << "   Second band:\n" << band2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with non-matching bands on a common submatrix
      {
         auto sm    = blaze::submatrix( tmat_, 1UL, 1UL, 4UL, 3UL );
         auto band1 = blaze::band( sm, -1L );
         auto band2 = blaze::band( sm,  0L );

         if( blaze::isSame( band1, band2 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First band:\n" << band1 << "\n"
                << "   Second band:\n" << band2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with matching rows on matrix and submatrix
      {
         auto sm    = blaze::submatrix( tmat_, 1UL, 0UL, 4UL, 4UL );
         auto band1 = blaze::band( tmat_, -1L );
         auto band2 = blaze::band( sm   ,  0L );

         if( blaze::isSame( band1, band2 ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First band:\n" << band1 << "\n"
                << "   Second band:\n" << band2 << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( band2, band1 ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First band:\n" << band1 << "\n"
                << "   Second band:\n" << band2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with non-matching rows on matrix and submatrix (different band)
      {
         auto sm    = blaze::submatrix( tmat_, 1UL, 0UL, 4UL, 4UL );
         auto band1 = blaze::band( tmat_, -2L );
         auto band2 = blaze::band( sm   ,  0L );

         if( blaze::isSame( band1, band2 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First band:\n" << band1 << "\n"
                << "   Second band:\n" << band2 << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( band2, band1 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First band:\n" << band1 << "\n"
                << "   Second band:\n" << band2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with non-matching rows on matrix and submatrix (different size)
      {
         auto sm    = blaze::submatrix( tmat_, 1UL, 0UL, 4UL, 3UL );
         auto band1 = blaze::band( tmat_, -1L );
         auto band2 = blaze::band( sm   ,  0L );

         if( blaze::isSame( band1, band2 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First band:\n" << band1 << "\n"
                << "   Second band:\n" << band2 << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( band2, band1 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First band:\n" << band1 << "\n"
                << "   Second band:\n" << band2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with matching rows on two submatrices
      {
         auto sm1   = blaze::submatrix( tmat_, 0UL, 1UL, 4UL, 3UL );
         auto sm2   = blaze::submatrix( tmat_, 1UL, 1UL, 5UL, 3UL );
         auto band1 = blaze::band( sm1, -1L );
         auto band2 = blaze::band( sm2,  0L );

         if( blaze::isSame( band1, band2 ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First band:\n" << band1 << "\n"
                << "   Second band:\n" << band2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with non-matching rows on two submatrices (different band)
      {
         auto sm1   = blaze::submatrix( tmat_, 0UL, 1UL, 4UL, 3UL );
         auto sm2   = blaze::submatrix( tmat_, 1UL, 1UL, 5UL, 3UL );
         auto band1 = blaze::band( sm1, -1L );
         auto band2 = blaze::band( sm2, -1L );

         if( blaze::isSame( band1, band2 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First band:\n" << band1 << "\n"
                << "   Second band:\n" << band2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with non-matching rows on two submatrices (different size)
      {
         auto sm1   = blaze::submatrix( tmat_, 0UL, 1UL, 4UL, 3UL );
         auto sm2   = blaze::submatrix( tmat_, 1UL, 1UL, 5UL, 2UL );
         auto band1 = blaze::band( sm1, -1L );
         auto band2 = blaze::band( sm2,  0L );

         if( blaze::isSame( band1, band2 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First band:\n" << band1 << "\n"
                << "   Second band:\n" << band2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with non-matching rows on two submatrices (different offset)
      {
         auto sm1   = blaze::submatrix( tmat_, 0UL, 0UL, 4UL, 3UL );
         auto sm2   = blaze::submatrix( tmat_, 1UL, 1UL, 5UL, 2UL );
         auto band1 = blaze::band( sm1, 0L );
         auto band2 = blaze::band( sm2, 0L );

         if( blaze::isSame( band1, band2 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First band:\n" << band1 << "\n"
                << "   Second band:\n" << band2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with matching band subvectors on submatrices
      {
         auto sm    = blaze::submatrix( tmat_, 1UL, 1UL, 4UL, 3UL );
         auto band1 = blaze::band( sm, -1L );
         auto sv1   = blaze::subvector( band1, 0UL, 2UL );
         auto sv2   = blaze::subvector( band1, 0UL, 2UL );

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

      // isSame with non-matching band subvectors on submatrices (different size)
      {
         auto sm    = blaze::submatrix( tmat_, 1UL, 1UL, 4UL, 3UL );
         auto band1 = blaze::band( sm, -1L );
         auto sv1   = blaze::subvector( band1, 0UL, 3UL );
         auto sv2   = blaze::subvector( band1, 0UL, 2UL );

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

      // isSame with non-matching band subvectors on submatrices (different offset)
      {
         auto sm    = blaze::submatrix( tmat_, 1UL, 1UL, 4UL, 3UL );
         auto band1 = blaze::band( sm, -1L );
         auto sv1   = blaze::subvector( band1, 0UL, 2UL );
         auto sv2   = blaze::subvector( band1, 1UL, 2UL );

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

      // isSame with matching row subvectors on two submatrices
      {
         auto sm1   = blaze::submatrix( tmat_, 0UL, 1UL, 4UL, 3UL );
         auto sm2   = blaze::submatrix( tmat_, 1UL, 1UL, 5UL, 3UL );
         auto band1 = blaze::band( sm1, -1L );
         auto band2 = blaze::band( sm2,  0L );
         auto sv1   = blaze::subvector( band1, 0UL, 2UL );
         auto sv2   = blaze::subvector( band2, 0UL, 2UL );

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

      // isSame with non-matching row subvectors on two submatrices (different size)
      {
         auto sm1   = blaze::submatrix( tmat_, 0UL, 1UL, 4UL, 3UL );
         auto sm2   = blaze::submatrix( tmat_, 1UL, 1UL, 5UL, 3UL );
         auto band1 = blaze::band( sm1, -1L );
         auto band2 = blaze::band( sm2,  0L );
         auto sv1   = blaze::subvector( band1, 0UL, 3UL );
         auto sv2   = blaze::subvector( band2, 0UL, 2UL );

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

      // isSame with non-matching row subvectors on two submatrices (different offset)
      {
         auto sm1   = blaze::submatrix( tmat_, 0UL, 1UL, 4UL, 3UL );
         auto sm2   = blaze::submatrix( tmat_, 1UL, 1UL, 5UL, 3UL );
         auto band1 = blaze::band( sm1, -1L );
         auto band2 = blaze::band( sm2,  0L );
         auto sv1   = blaze::subvector( band1, 0UL, 2UL );
         auto sv2   = blaze::subvector( band2, 1UL, 2UL );

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
/*!\brief Test of the \c subvector() function with the Band specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c subvector() function used with the Band specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testSubvector()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major subvector() function";

      initialize();

      {
         BT   band1 = blaze::band( mat_, 1L );
         auto sv    = blaze::subvector( band1, 0UL, 4UL );

         if( sv[0] != 0 || sv[1] != 4 || sv[2] != 5 || sv[3] != -6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subscript operator access failed\n"
                << " Details:\n"
                << "   Result:\n" << sv << "\n"
                << "   Expected result:\n( 0 4 5 -6 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( *sv.begin() != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator access failed\n"
                << " Details:\n"
                << "   Result: " << *sv.begin() << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      try {
         BT   band1 = blaze::band( mat_, 1L );
         auto sv    = blaze::subvector( band1, 4UL, 4UL );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of out-of-bounds subvector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      try {
         BT   band1 = blaze::band( mat_, 1L );
         auto sv    = blaze::subvector( band1, 0UL, 5UL );

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
         OBT  band1 = blaze::band( tmat_, -1L );
         auto sv    = blaze::subvector( band1, 0UL, 4UL );

         if( sv[0] != 0 || sv[1] != 4 || sv[2] != 5 || sv[3] != -6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subscript operator access failed\n"
                << " Details:\n"
                << "   Result:\n" << sv << "\n"
                << "   Expected result:\n( 0 4 5 -6 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( *sv.begin() != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator access failed\n"
                << " Details:\n"
                << "   Result: " << *sv.begin() << "\n"
                << "   Expected result: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      try {
         OBT  band1 = blaze::band( tmat_, -1L );
         auto sv    = blaze::subvector( band1, 4UL, 4UL );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of out-of-bounds subvector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      try {
         OBT  band1 = blaze::band( tmat_, -1L );
         auto sv    = blaze::subvector( band1, 0UL, 5UL );

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
/*!\brief Test of the \c elements() function with the Band specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c elements() function used with the Band specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testElements()
{
   //=====================================================================================
   // Row-major matrix tests (initializer_list)
   //=====================================================================================

   {
      test_ = "Row-major elements() function (initializer_list)";

      initialize();

      {
         BT   band1 = blaze::band( mat_, 1L );
         auto e     = blaze::elements( band1, { 3UL, 2UL } );

         if( e[0] != -6 || e[1] != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subscript operator access failed\n"
                << " Details:\n"
                << "   Result:\n" << e << "\n"
                << "   Expected result:\n( -6 5 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( *e.begin() != -6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator access failed\n"
                << " Details:\n"
                << "   Result: " << *e.begin() << "\n"
                << "   Expected result: -6\n";
            throw std::runtime_error( oss.str() );
         }
      }

      try {
         BT   band1 = blaze::band( mat_, 1L );
         auto e     = blaze::elements( band1, { 4UL } );

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
         std::array<int,2UL> indices{ 3UL, 2UL };

         BT   band1 = blaze::band( mat_, 1L );
         auto e     = blaze::elements( band1, indices );

         if( e[0] != -6 || e[1] != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subscript operator access failed\n"
                << " Details:\n"
                << "   Result:\n" << e << "\n"
                << "   Expected result:\n( -6 5 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( *e.begin() != -6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator access failed\n"
                << " Details:\n"
                << "   Result: " << *e.begin() << "\n"
                << "   Expected result: -6\n";
            throw std::runtime_error( oss.str() );
         }
      }

      try {
         std::array<int,2UL> indices{ 4UL };

         BT   band1 = blaze::band( mat_, 1L );
         auto e     = blaze::elements( band1, indices );

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
         BT   band1 = blaze::band( mat_, 1L );
         auto e     = blaze::elements( band1, []( size_t i ){ return 3UL-i; }, 2UL );

         if( e[0] != -6 || e[1] != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subscript operator access failed\n"
                << " Details:\n"
                << "   Result:\n" << e << "\n"
                << "   Expected result:\n( -6 5 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( *e.begin() != -6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator access failed\n"
                << " Details:\n"
                << "   Result: " << *e.begin() << "\n"
                << "   Expected result: -6\n";
            throw std::runtime_error( oss.str() );
         }
      }

      try {
         BT   band1 = blaze::band( mat_, 1L );
         auto e     = blaze::elements( band1, []( size_t ){ return 4UL; }, 1UL );

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
         OBT  band1 = blaze::band( tmat_, -1L );
         auto e     = blaze::elements( band1, { 3UL, 2UL } );

         if( e[0] != -6 || e[1] != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subscript operator access failed\n"
                << " Details:\n"
                << "   Result:\n" << e << "\n"
                << "   Expected result:\n( -6 5 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( *e.begin() != -6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator access failed\n"
                << " Details:\n"
                << "   Result: " << *e.begin() << "\n"
                << "   Expected result: -6\n";
            throw std::runtime_error( oss.str() );
         }
      }

      try {
         OBT  band1 = blaze::band( tmat_, -1L );
         auto e     = blaze::elements( band1, { 4UL } );

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
         std::array<int,2UL> indices{ 3UL, 2UL };

         OBT  band1 = blaze::band( tmat_, -1L );
         auto e     = blaze::elements( band1, indices );

         if( e[0] != -6 || e[1] != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subscript operator access failed\n"
                << " Details:\n"
                << "   Result:\n" << e << "\n"
                << "   Expected result:\n( -6 5 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( *e.begin() != -6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator access failed\n"
                << " Details:\n"
                << "   Result: " << *e.begin() << "\n"
                << "   Expected result: -6\n";
            throw std::runtime_error( oss.str() );
         }
      }

      try {
         std::array<int,1UL> indices{ 4UL };

         OBT  band1 = blaze::band( tmat_, -1L );
         auto e     = blaze::elements( band1, indices );

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
         OBT  band1 = blaze::band( tmat_, -1L );
         auto e     = blaze::elements( band1, []( size_t i ){ return 3UL-i; }, 2UL );

         if( e[0] != -6 || e[1] != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subscript operator access failed\n"
                << " Details:\n"
                << "   Result:\n" << e << "\n"
                << "   Expected result:\n( -6 5 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( *e.begin() != -6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator access failed\n"
                << " Details:\n"
                << "   Result: " << *e.begin() << "\n"
                << "   Expected result: -6\n";
            throw std::runtime_error( oss.str() );
         }
      }

      try {
         OBT  band1 = blaze::band( tmat_, -1L );
         auto e     = blaze::elements( band1, []( size_t ){ return 4UL; }, 1UL );

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
void DenseTest::initialize()
{
   // Initializing the row-major dynamic matrix
   mat_.reset();
   mat_(0,0) = -2;
   mat_(0,2) =  7;
   mat_(1,2) =  4;
   mat_(1,3) = -8;
   mat_(2,1) =  1;
   mat_(2,2) = -3;
   mat_(2,3) =  5;
   mat_(2,4) =  9;
   mat_(3,4) = -6;
   mat_(3,5) = 10;

   // Initializing the column-major dynamic matrix
   tmat_ = trans( mat_ );
}
//*************************************************************************************************

} // namespace band

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
   std::cout << "   Running Band dense test..." << std::endl;

   try
   {
      RUN_BAND_DENSE_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during Band dense test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
