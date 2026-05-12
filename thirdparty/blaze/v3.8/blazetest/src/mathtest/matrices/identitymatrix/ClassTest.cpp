//=================================================================================================
/*!
//  \file src/mathtest/matrices/identitymatrix/ClassTest.cpp
//  \brief Source file for the IdentityMatrix class test
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
#include <blaze/math/CompressedMatrix.h>
#include <blaze/math/DynamicMatrix.h>
#include <blaze/util/Complex.h>
#include <blaze/util/Random.h>
#include <blazetest/mathtest/matrices/identitymatrix/ClassTest.h>

#ifdef BLAZE_USE_HPX_THREADS
#  include <hpx/hpx_main.hpp>
#endif


namespace blazetest {

namespace mathtest {

namespace matrices {

namespace identitymatrix {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the IdentityMatrix class test.
//
// \exception std::runtime_error Operation error detected.
*/
ClassTest::ClassTest()
{
   testConstructors();
   testAssignment();
   testFunctionCall();
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
   testTranspose();
   testCTranspose();
   testIsDefault();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the IdentityMatrix constructors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all constructors of the IdentityMatrix class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testConstructors()
{
   //=====================================================================================
   // Row-major default constructor
   //=====================================================================================

   {
      test_ = "Row-major IdentityMatrix default constructor";

      blaze::IdentityMatrix<int,blaze::rowMajor> I;

      checkRows    ( I, 0UL );
      checkColumns ( I, 0UL );
      checkNonZeros( I, 0UL );
   }


   //=====================================================================================
   // Row-major size constructor
   //=====================================================================================

   {
      test_ = "Row-major IdentityMatrix size constructor (0x0)";

      blaze::IdentityMatrix<int,blaze::rowMajor> I( 0UL );

      checkRows    ( I, 0UL );
      checkColumns ( I, 0UL );
      checkNonZeros( I, 0UL );
   }

   {
      test_ = "Row-major IdentityMatrix size constructor (3x3)";

      blaze::IdentityMatrix<int,blaze::rowMajor> I( 3UL );

      checkRows    ( I, 3UL );
      checkColumns ( I, 3UL );
      checkNonZeros( I, 3UL );
      checkNonZeros( I, 0UL, 1UL );
      checkNonZeros( I, 1UL, 1UL );
      checkNonZeros( I, 2UL, 1UL );

      if( I(0,0) != 1 || I(0,1) != 0 || I(0,2) != 0 ||
          I(1,0) != 0 || I(1,1) != 1 || I(1,2) != 0 ||
          I(2,0) != 0 || I(2,1) != 0 || I(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << I << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 1 0 )\n( 0 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major copy constructor
   //=====================================================================================

   {
      test_ = "Row-major IdentityMatrix copy constructor (0x0)";

      blaze::IdentityMatrix<int,blaze::rowMajor> I1( 0UL );
      blaze::IdentityMatrix<int,blaze::rowMajor> I2( I1 );

      checkRows    ( I2, 0UL );
      checkColumns ( I2, 0UL );
      checkNonZeros( I2, 0UL );
   }

   {
      test_ = "Row-major IdentityMatrix copy constructor (3x3)";

      blaze::IdentityMatrix<int,blaze::rowMajor> I1( 3UL );
      blaze::IdentityMatrix<int,blaze::rowMajor> I2( I1 );

      checkRows    ( I2, 3UL );
      checkColumns ( I2, 3UL );
      checkNonZeros( I2, 3UL );
      checkNonZeros( I2, 0UL, 1UL );
      checkNonZeros( I2, 1UL, 1UL );
      checkNonZeros( I2, 2UL, 1UL );

      if( I2(0,0) != 1 || I2(0,1) != 0 || I2(0,2) != 0 ||
          I2(1,0) != 0 || I2(1,1) != 1 || I2(1,2) != 0 ||
          I2(2,0) != 0 || I2(2,1) != 0 || I2(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << I2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 1 0 )\n( 0 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major move constructor
   //=====================================================================================

   {
      test_ = "Row-major IdentityMatrix move constructor (0x0)";

      blaze::IdentityMatrix<int,blaze::rowMajor> I1( 0UL );
      blaze::IdentityMatrix<int,blaze::rowMajor> I2( std::move( I1 ) );

      checkRows    ( I2, 0UL );
      checkColumns ( I2, 0UL );
      checkNonZeros( I2, 0UL );
   }

   {
      test_ = "Row-major IdentityMatrix move constructor (3x3)";

      blaze::IdentityMatrix<int,blaze::rowMajor> I1( 3UL );
      blaze::IdentityMatrix<int,blaze::rowMajor> I2( std::move( I1 ) );

      checkRows    ( I2, 3UL );
      checkColumns ( I2, 3UL );
      checkNonZeros( I2, 3UL );
      checkNonZeros( I2, 0UL, 1UL );
      checkNonZeros( I2, 1UL, 1UL );
      checkNonZeros( I2, 2UL, 1UL );

      if( I2(0,0) != 1 || I2(0,1) != 0 || I2(0,2) != 0 ||
          I2(1,0) != 0 || I2(1,1) != 1 || I2(1,2) != 0 ||
          I2(2,0) != 0 || I2(2,1) != 0 || I2(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << I2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 1 0 )\n( 0 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major dense matrix constructor
   //=====================================================================================

   {
      test_ = "Row-major/row-major IdentityMatrix dense matrix constructor";

      blaze::DynamicMatrix<int,blaze::rowMajor> I1{ { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
      blaze::IdentityMatrix<int,blaze::rowMajor> I2( I1 );

      checkRows    ( I2, 3UL );
      checkColumns ( I2, 3UL );
      checkNonZeros( I2, 3UL );
      checkNonZeros( I2, 0UL, 1UL );
      checkNonZeros( I2, 1UL, 1UL );
      checkNonZeros( I2, 2UL, 1UL );

      if( I2(0,0) != 1 || I2(0,1) != 0 || I2(0,2) != 0 ||
          I2(1,0) != 0 || I2(1,1) != 1 || I2(1,2) != 0 ||
          I2(2,0) != 0 || I2(2,1) != 0 || I2(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << I2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 1 0 )\n( 0 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major IdentityMatrix dense matrix constructor";

      blaze::DynamicMatrix<int,blaze::columnMajor> I1{ { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
      blaze::IdentityMatrix<int,blaze::rowMajor> I2( I1 );

      checkRows    ( I2, 3UL );
      checkColumns ( I2, 3UL );
      checkNonZeros( I2, 3UL );
      checkNonZeros( I2, 0UL, 1UL );
      checkNonZeros( I2, 1UL, 1UL );
      checkNonZeros( I2, 2UL, 1UL );

      if( I2(0,0) != 1 || I2(0,1) != 0 || I2(0,2) != 0 ||
          I2(1,0) != 0 || I2(1,1) != 1 || I2(1,2) != 0 ||
          I2(2,0) != 0 || I2(2,1) != 0 || I2(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << I2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 1 0 )\n( 0 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major IdentityMatrix dense matrix constructor (non-square)";

      blaze::DynamicMatrix<int,blaze::rowMajor> I1{ { 1, 0, 0 }, { 0, 1, 0 } };

      try {
         blaze::IdentityMatrix<int,blaze::rowMajor> I2( I1 );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of non-identity IdentityMatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << I2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   {
      test_ = "Row-major IdentityMatrix dense matrix constructor (non-identity)";

      blaze::DynamicMatrix<int,blaze::rowMajor> I1( 3UL, 3UL, 0 );

      try {
         blaze::IdentityMatrix<int,blaze::rowMajor> I2( I1 );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of non-identity IdentityMatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << I2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Row-major sparse matrix constructor
   //=====================================================================================

   {
      test_ = "Row-major/row-major IdentityMatrix sparse matrix constructor";

      blaze::CompressedMatrix<int,blaze::rowMajor> I1( 3UL, 3UL, 3UL );
      I1(0,0) = 1;
      I1(1,1) = 1;
      I1(2,2) = 1;

      blaze::IdentityMatrix<int,blaze::rowMajor> I2( I1 );

      checkRows    ( I2, 3UL );
      checkColumns ( I2, 3UL );
      checkNonZeros( I2, 3UL );
      checkNonZeros( I2, 0UL, 1UL );
      checkNonZeros( I2, 1UL, 1UL );
      checkNonZeros( I2, 2UL, 1UL );

      if( I2(0,0) != 1 || I2(0,1) != 0 || I2(0,2) != 0 ||
          I2(1,0) != 0 || I2(1,1) != 1 || I2(1,2) != 0 ||
          I2(2,0) != 0 || I2(2,1) != 0 || I2(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << I2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 1 0 )\n( 0 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major IdentityMatrix sparse matrix constructor";

      blaze::CompressedMatrix<int,blaze::columnMajor> I1( 3UL, 3UL, 3UL );
      I1(0,0) = 1;
      I1(1,1) = 1;
      I1(2,2) = 1;

      blaze::IdentityMatrix<int,blaze::rowMajor> I2( I1 );

      checkRows    ( I2, 3UL );
      checkColumns ( I2, 3UL );
      checkNonZeros( I2, 3UL );
      checkNonZeros( I2, 0UL, 1UL );
      checkNonZeros( I2, 1UL, 1UL );
      checkNonZeros( I2, 2UL, 1UL );

      if( I2(0,0) != 1 || I2(0,1) != 0 || I2(0,2) != 0 ||
          I2(1,0) != 0 || I2(1,1) != 1 || I2(1,2) != 0 ||
          I2(2,0) != 0 || I2(2,1) != 0 || I2(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << I2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 1 0 )\n( 0 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major IdentityMatrix sparse matrix constructor (non-square)";

      blaze::CompressedMatrix<int,blaze::rowMajor> I1( 2UL, 3UL, 2UL );
      I1(0,0) = 1;
      I1(1,1) = 1;

      try {
         blaze::IdentityMatrix<int,blaze::rowMajor> I2( I1 );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of non-identity IdentityMatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << I2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   {
      test_ = "Row-major IdentityMatrix sparse matrix constructor (non-identity)";

      blaze::CompressedMatrix<int,blaze::rowMajor> I1( 3UL, 3UL );

      try {
         blaze::IdentityMatrix<int,blaze::rowMajor> I2( I1 );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of non-identity IdentityMatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << I2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Column-major default constructor
   //=====================================================================================

   {
      test_ = "Column-major IdentityMatrix default constructor";

      blaze::IdentityMatrix<int,blaze::columnMajor> I;

      checkRows    ( I, 0UL );
      checkColumns ( I, 0UL );
      checkNonZeros( I, 0UL );
   }


   //=====================================================================================
   // Column-major size constructor
   //=====================================================================================

   {
      test_ = "Column-major IdentityMatrix size constructor (0x0)";

      blaze::IdentityMatrix<int,blaze::columnMajor> I( 0UL );

      checkRows    ( I, 0UL );
      checkColumns ( I, 0UL );
      checkNonZeros( I, 0UL );
   }

   {
      test_ = "Column-major IdentityMatrix size constructor (3x3)";

      blaze::IdentityMatrix<int,blaze::columnMajor> I( 3UL );

      checkRows    ( I, 3UL );
      checkColumns ( I, 3UL );
      checkNonZeros( I, 3UL );
      checkNonZeros( I, 0UL, 1UL );
      checkNonZeros( I, 1UL, 1UL );
      checkNonZeros( I, 2UL, 1UL );

      if( I(0,0) != 1 || I(0,1) != 0 || I(0,2) != 0 ||
          I(1,0) != 0 || I(1,1) != 1 || I(1,2) != 0 ||
          I(2,0) != 0 || I(2,1) != 0 || I(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << I << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 1 0 )\n( 0 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major copy constructor
   //=====================================================================================

   {
      test_ = "Column-major IdentityMatrix copy constructor (0x0)";

      blaze::IdentityMatrix<int,blaze::columnMajor> I1( 0UL );
      blaze::IdentityMatrix<int,blaze::columnMajor> I2( I1 );

      checkRows    ( I2, 0UL );
      checkColumns ( I2, 0UL );
      checkNonZeros( I2, 0UL );
   }

   {
      test_ = "Column-major IdentityMatrix copy constructor (3x3)";

      blaze::IdentityMatrix<int,blaze::columnMajor> I1( 3UL );
      blaze::IdentityMatrix<int,blaze::columnMajor> I2( I1 );

      checkRows    ( I2, 3UL );
      checkColumns ( I2, 3UL );
      checkNonZeros( I2, 3UL );
      checkNonZeros( I2, 0UL, 1UL );
      checkNonZeros( I2, 1UL, 1UL );
      checkNonZeros( I2, 2UL, 1UL );

      if( I2(0,0) != 1 || I2(0,1) != 0 || I2(0,2) != 0 ||
          I2(1,0) != 0 || I2(1,1) != 1 || I2(1,2) != 0 ||
          I2(2,0) != 0 || I2(2,1) != 0 || I2(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << I2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 1 0 )\n( 0 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major move constructor
   //=====================================================================================

   {
      test_ = "Column-major IdentityMatrix move constructor (0x0)";

      blaze::IdentityMatrix<int,blaze::columnMajor> I1( 0UL );
      blaze::IdentityMatrix<int,blaze::columnMajor> I2( std::move( I1 ) );

      checkRows    ( I2, 0UL );
      checkColumns ( I2, 0UL );
      checkNonZeros( I2, 0UL );
   }

   {
      test_ = "Column-major IdentityMatrix move constructor (3x3)";

      blaze::IdentityMatrix<int,blaze::columnMajor> I1( 3UL );
      blaze::IdentityMatrix<int,blaze::columnMajor> I2( std::move( I1 ) );

      checkRows    ( I2, 3UL );
      checkColumns ( I2, 3UL );
      checkNonZeros( I2, 3UL );
      checkNonZeros( I2, 0UL, 1UL );
      checkNonZeros( I2, 1UL, 1UL );
      checkNonZeros( I2, 2UL, 1UL );

      if( I2(0,0) != 1 || I2(0,1) != 0 || I2(0,2) != 0 ||
          I2(1,0) != 0 || I2(1,1) != 1 || I2(1,2) != 0 ||
          I2(2,0) != 0 || I2(2,1) != 0 || I2(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << I2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 1 0 )\n( 0 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix constructor
   //=====================================================================================

   {
      test_ = "Column-major/row-major IdentityMatrix dense matrix constructor";

      blaze::DynamicMatrix<int,blaze::rowMajor> I1{ { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
      blaze::IdentityMatrix<int,blaze::columnMajor> I2( I1 );

      checkRows    ( I2, 3UL );
      checkColumns ( I2, 3UL );
      checkNonZeros( I2, 3UL );
      checkNonZeros( I2, 0UL, 1UL );
      checkNonZeros( I2, 1UL, 1UL );
      checkNonZeros( I2, 2UL, 1UL );

      if( I2(0,0) != 1 || I2(0,1) != 0 || I2(0,2) != 0 ||
          I2(1,0) != 0 || I2(1,1) != 1 || I2(1,2) != 0 ||
          I2(2,0) != 0 || I2(2,1) != 0 || I2(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << I2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 1 0 )\n( 0 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major IdentityMatrix dense matrix constructor";

      blaze::DynamicMatrix<int,blaze::columnMajor> I1{ { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
      blaze::IdentityMatrix<int,blaze::columnMajor> I2( I1 );

      checkRows    ( I2, 3UL );
      checkColumns ( I2, 3UL );
      checkNonZeros( I2, 3UL );
      checkNonZeros( I2, 0UL, 1UL );
      checkNonZeros( I2, 1UL, 1UL );
      checkNonZeros( I2, 2UL, 1UL );

      if( I2(0,0) != 1 || I2(0,1) != 0 || I2(0,2) != 0 ||
          I2(1,0) != 0 || I2(1,1) != 1 || I2(1,2) != 0 ||
          I2(2,0) != 0 || I2(2,1) != 0 || I2(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << I2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 1 0 )\n( 0 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major IdentityMatrix dense matrix constructor (non-square)";

      blaze::DynamicMatrix<int,blaze::columnMajor> I1{ { 1, 0, 0 }, { 0, 1, 0 } };

      try {
         blaze::IdentityMatrix<int,blaze::columnMajor> I2( I1 );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of non-identity IdentityMatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << I2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   {
      test_ = "Column-major IdentityMatrix dense matrix constructor (non-identity)";

      blaze::DynamicMatrix<int,blaze::columnMajor> I1( 3UL, 3UL, 0 );

      try {
         blaze::IdentityMatrix<int,blaze::columnMajor> I2( I1 );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of non-identity IdentityMatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << I2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Column-major sparse matrix constructor
   //=====================================================================================

   {
      test_ = "Column-major/row-major IdentityMatrix sparse matrix constructor";

      blaze::CompressedMatrix<int,blaze::rowMajor> I1( 3UL, 3UL, 3UL );
      I1(0,0) = 1;
      I1(1,1) = 1;
      I1(2,2) = 1;

      blaze::IdentityMatrix<int,blaze::columnMajor> I2( I1 );

      checkRows    ( I2, 3UL );
      checkColumns ( I2, 3UL );
      checkNonZeros( I2, 3UL );
      checkNonZeros( I2, 0UL, 1UL );
      checkNonZeros( I2, 1UL, 1UL );
      checkNonZeros( I2, 2UL, 1UL );

      if( I2(0,0) != 1 || I2(0,1) != 0 || I2(0,2) != 0 ||
          I2(1,0) != 0 || I2(1,1) != 1 || I2(1,2) != 0 ||
          I2(2,0) != 0 || I2(2,1) != 0 || I2(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << I2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 1 0 )\n( 0 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major IdentityMatrix sparse matrix constructor";

      blaze::CompressedMatrix<int,blaze::columnMajor> I1( 3UL, 3UL, 3UL );
      I1(0,0) = 1;
      I1(1,1) = 1;
      I1(2,2) = 1;

      blaze::IdentityMatrix<int,blaze::columnMajor> I2( I1 );

      checkRows    ( I2, 3UL );
      checkColumns ( I2, 3UL );
      checkNonZeros( I2, 3UL );
      checkNonZeros( I2, 0UL, 1UL );
      checkNonZeros( I2, 1UL, 1UL );
      checkNonZeros( I2, 2UL, 1UL );

      if( I2(0,0) != 1 || I2(0,1) != 0 || I2(0,2) != 0 ||
          I2(1,0) != 0 || I2(1,1) != 1 || I2(1,2) != 0 ||
          I2(2,0) != 0 || I2(2,1) != 0 || I2(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << I2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 1 0 )\n( 0 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major IdentityMatrix sparse matrix constructor (non-square)";

      blaze::CompressedMatrix<int,blaze::columnMajor> I1( 2UL, 3UL, 2UL );
      I1(0,0) = 1;
      I1(1,1) = 1;

      try {
         blaze::IdentityMatrix<int,blaze::columnMajor> I2( I1 );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of non-identity IdentityMatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << I2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   {
      test_ = "Column-major IdentityMatrix sparse matrix constructor (non-identity)";

      blaze::CompressedMatrix<int,blaze::columnMajor> I1( 3UL, 3UL );

      try {
         blaze::IdentityMatrix<int,blaze::columnMajor> I2( I1 );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of non-identity IdentityMatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << I2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the IdentityMatrix assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all assignment operators of the IdentityMatrix class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testAssignment()
{
   //=====================================================================================
   // Row-major copy assignment
   //=====================================================================================

   {
      test_ = "Row-major IdentityMatrix copy assignment";

      blaze::IdentityMatrix<int,blaze::rowMajor> I1( 3UL );
      blaze::IdentityMatrix<int,blaze::rowMajor> I2;
      I2 = I1;

      checkRows    ( I2, 3UL );
      checkColumns ( I2, 3UL );
      checkNonZeros( I2, 3UL );
      checkNonZeros( I2, 0UL, 1UL );
      checkNonZeros( I2, 1UL, 1UL );
      checkNonZeros( I2, 2UL, 1UL );

      if( I2(0,0) != 1 || I2(0,1) != 0 || I2(0,2) != 0 ||
          I2(1,0) != 0 || I2(1,1) != 1 || I2(1,2) != 0 ||
          I2(2,0) != 0 || I2(2,1) != 0 || I2(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << I2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 1 0 )\n( 0 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major IdentityMatrix copy assignment stress test";

      blaze::IdentityMatrix<int,blaze::rowMajor> I1;

      for( size_t i=0UL; i<100UL; ++i )
      {
         const size_t n( blaze::rand<size_t>( 0UL, 10UL ) );
         const blaze::IdentityMatrix<int,blaze::rowMajor> I2( n );

         I1 = I2;

         if( I1 != I2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << I1 << "\n"
                << "   Expected result:\n" << I2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Row-major move assignment
   //=====================================================================================

   {
      test_ = "Row-major IdentityMatrix move assignment";

      blaze::IdentityMatrix<int,blaze::rowMajor> I1( 3UL );
      blaze::IdentityMatrix<int,blaze::rowMajor> I2;

      I2 = std::move( I1 );

      checkRows    ( I2, 3UL );
      checkColumns ( I2, 3UL );
      checkNonZeros( I2, 3UL );
      checkNonZeros( I2, 0UL, 1UL );
      checkNonZeros( I2, 1UL, 1UL );
      checkNonZeros( I2, 2UL, 1UL );

      if( I2(0,0) != 1 || I2(0,1) != 0 || I2(0,2) != 0 ||
          I2(1,0) != 0 || I2(1,1) != 1 || I2(1,2) != 0 ||
          I2(2,0) != 0 || I2(2,1) != 0 || I2(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << I2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 1 0 )\n( 0 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major copy assignment
   //=====================================================================================

   {
      test_ = "Column-major IdentityMatrix copy assignment";

      blaze::IdentityMatrix<int,blaze::columnMajor> I1( 3UL );
      blaze::IdentityMatrix<int,blaze::columnMajor> I2;
      I2 = I1;

      checkRows    ( I2, 3UL );
      checkColumns ( I2, 3UL );
      checkNonZeros( I2, 3UL );
      checkNonZeros( I2, 0UL, 1UL );
      checkNonZeros( I2, 1UL, 1UL );
      checkNonZeros( I2, 2UL, 1UL );

      if( I2(0,0) != 1 || I2(0,1) != 0 || I2(0,2) != 0 ||
          I2(1,0) != 0 || I2(1,1) != 1 || I2(1,2) != 0 ||
          I2(2,0) != 0 || I2(2,1) != 0 || I2(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << I2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 1 0 )\n( 0 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major IdentityMatrix copy assignment stress test";

      blaze::IdentityMatrix<int,blaze::columnMajor> I1;

      for( size_t i=0UL; i<100UL; ++i )
      {
         const size_t n( blaze::rand<size_t>( 0UL, 10UL ) );
         const blaze::IdentityMatrix<int,blaze::columnMajor> I2( n );

         I1 = I2;

         if( I1 != I2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << I1 << "\n"
                << "   Expected result:\n" << I2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major move assignment
   //=====================================================================================

   {
      test_ = "Column-major IdentityMatrix move assignment";

      blaze::IdentityMatrix<int,blaze::columnMajor> I1( 3UL );
      blaze::IdentityMatrix<int,blaze::columnMajor> I2;

      I2 = std::move( I1 );

      checkRows    ( I2, 3UL );
      checkColumns ( I2, 3UL );
      checkNonZeros( I2, 3UL );
      checkNonZeros( I2, 0UL, 1UL );
      checkNonZeros( I2, 1UL, 1UL );
      checkNonZeros( I2, 2UL, 1UL );

      if( I2(0,0) != 1 || I2(0,1) != 0 || I2(0,2) != 0 ||
          I2(1,0) != 0 || I2(1,1) != 1 || I2(1,2) != 0 ||
          I2(2,0) != 0 || I2(2,1) != 0 || I2(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << I2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 1 0 )\n( 0 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the IdentityMatrix function call operator.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of adding and accessing elements via the function call operator
// of the IdentityMatrix class template. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
void ClassTest::testFunctionCall()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major IdentityMatrix::operator()";

      blaze::IdentityMatrix<int,blaze::rowMajor> I( 3UL );

      checkRows    ( I, 3UL );
      checkColumns ( I, 3UL );
      checkNonZeros( I, 3UL );
      checkNonZeros( I, 0UL, 1UL );
      checkNonZeros( I, 1UL, 1UL );
      checkNonZeros( I, 2UL, 1UL );

      if( I(0,0) != 1 || I(0,1) != 0 || I(0,2) != 0 ||
          I(1,0) != 0 || I(1,1) != 1 || I(1,2) != 0 ||
          I(2,0) != 0 || I(2,1) != 0 || I(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << I << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 1 0 )\n( 0 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major IdentityMatrix::operator()";

      blaze::IdentityMatrix<int,blaze::columnMajor> I( 3UL );

      checkRows    ( I, 3UL );
      checkColumns ( I, 3UL );
      checkNonZeros( I, 3UL );
      checkNonZeros( I, 0UL, 1UL );
      checkNonZeros( I, 1UL, 1UL );
      checkNonZeros( I, 2UL, 1UL );

      if( I(0,0) != 1 || I(0,1) != 0 || I(0,2) != 0 ||
          I(1,0) != 0 || I(1,1) != 1 || I(1,2) != 0 ||
          I(2,0) != 0 || I(2,1) != 0 || I(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << I << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 1 0 )\n( 0 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c at() member function of the IdentityMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of adding and accessing elements via the \c at() member function
// of the IdentityMatrix class template. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
void ClassTest::testAt()
{
   //=====================================================================================
   // Row-major matrix tests
   //==========≈===========================================================================

   {
      test_ = "Row-major IdentityMatrix::at()";

      blaze::IdentityMatrix<int,blaze::rowMajor> I( 3UL );

      checkRows    ( I, 3UL );
      checkColumns ( I, 3UL );
      checkNonZeros( I, 3UL );
      checkNonZeros( I, 0UL, 1UL );
      checkNonZeros( I, 1UL, 1UL );
      checkNonZeros( I, 2UL, 1UL );

      if( I.at(0,0) != 1 || I.at(0,1) != 0 || I.at(0,2) != 0 ||
          I.at(1,0) != 0 || I.at(1,1) != 1 || I.at(1,2) != 0 ||
          I.at(2,0) != 0 || I.at(2,1) != 0 || I.at(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Access via at() function failed\n"
             << " Details:\n"
             << "   Result:\n" << I << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 1 0 )\n( 0 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major IdentityMatrix::at()";

      blaze::IdentityMatrix<int,blaze::columnMajor> I( 3UL );

      checkRows    ( I, 3UL );
      checkColumns ( I, 3UL );
      checkNonZeros( I, 3UL );
      checkNonZeros( I, 0UL, 1UL );
      checkNonZeros( I, 1UL, 1UL );
      checkNonZeros( I, 2UL, 1UL );

      if( I.at(0,0) != 1 || I.at(0,1) != 0 || I.at(0,2) != 0 ||
          I.at(1,0) != 0 || I.at(1,1) != 1 || I.at(1,2) != 0 ||
          I.at(2,0) != 0 || I.at(2,1) != 0 || I.at(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Access via at() function failed\n"
             << " Details:\n"
             << "   Result:\n" << I << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 1 0 )\n( 0 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the IdentityMatrix iterator implementation.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the iterator implementation of the IdentityMatrix class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testIterator()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      using MatrixType    = blaze::IdentityMatrix<int,blaze::rowMajor>;
      using ConstIterator = MatrixType::ConstIterator;

      MatrixType I( 3UL );

      // Testing the ConstIterator default constructor
      {
         test_ = "Row-major ConstIterator default constructor";

         ConstIterator it{};

         if( it != ConstIterator() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Failed iterator default constructor\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Counting the number of elements in 1st row via ConstIterator (end-begin)
      {
         test_ = "Row-major ConstIterator subtraction";

         const ptrdiff_t number( cend( I, 1UL ) - cbegin( I, 1UL ) );

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

      // Testing read-only access via ConstIterator (pre-increment)
      {
         test_ = "Row-major read-only access via ConstIterator";

         ConstIterator it ( cbegin( I, 1UL ) );
         ConstIterator end( cend( I, 1UL ) );

         if( it == end || it->value() != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid initial iterator detected\n";
            throw std::runtime_error( oss.str() );
         }

         ++it;

         if( it != cend( I, 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator pre-increment failed\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing read-only access via ConstIterator (post-increment)
      {
         test_ = "Row-major read-only access via ConstIterator";

         ConstIterator it ( cbegin( I, 1UL ) );
         ConstIterator end( cend( I, 1UL ) );

         if( it == end || it->value() != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid initial iterator detected\n";
            throw std::runtime_error( oss.str() );
         }

         it++;

         if( it != cend( I, 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator post-increment failed\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      using MatrixType    = blaze::IdentityMatrix<int,blaze::columnMajor>;
      using ConstIterator = MatrixType::ConstIterator;

      MatrixType I( 3UL );

      // Testing the ConstIterator default constructor
      {
         test_ = "Column-major ConstIterator default constructor";

         ConstIterator it{};

         if( it != ConstIterator() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Failed iterator default constructor\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Counting the number of elements in 1st row via ConstIterator (end-begin)
      {
         test_ = "Column-major ConstIterator subtraction";

         const ptrdiff_t number( cend( I, 1UL ) - cbegin( I, 1UL ) );

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

      // Testing read-only access via ConstIterator (pre-increment)
      {
         test_ = "Column-major read-only access via ConstIterator";

         ConstIterator it ( cbegin( I, 1UL ) );
         ConstIterator end( cend( I, 1UL ) );

         if( it == end || it->value() != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid initial iterator detected\n";
            throw std::runtime_error( oss.str() );
         }

         ++it;

         if( it != cend( I, 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator pre-increment failed\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing read-only access via ConstIterator (post-increment)
      {
         test_ = "Column-major read-only access via ConstIterator";

         ConstIterator it ( cbegin( I, 1UL ) );
         ConstIterator end( cend( I, 1UL ) );

         if( it == end || it->value() != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid initial iterator detected\n";
            throw std::runtime_error( oss.str() );
         }

         it++;

         if( it != cend( I, 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator post-increment failed\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c nonZeros() member function of the IdentityMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c nonZeros() member function of the IdentityMatrix
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testNonZeros()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major IdentityMatrix::nonZeros()";

      blaze::IdentityMatrix<int,blaze::rowMajor> I( 6UL );

      checkRows    ( I, 6UL );
      checkColumns ( I, 6UL );
      checkNonZeros( I, 6UL );
      checkNonZeros( I, 0UL, 1UL );
      checkNonZeros( I, 1UL, 1UL );
      checkNonZeros( I, 2UL, 1UL );
      checkNonZeros( I, 3UL, 1UL );
      checkNonZeros( I, 4UL, 1UL );
      checkNonZeros( I, 5UL, 1UL );
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major IdentityMatrix::nonZeros()";

      blaze::IdentityMatrix<int,blaze::columnMajor> I( 6UL );

      checkRows    ( I, 6UL );
      checkColumns ( I, 6UL );
      checkNonZeros( I, 6UL );
      checkNonZeros( I, 0UL, 1UL );
      checkNonZeros( I, 1UL, 1UL );
      checkNonZeros( I, 2UL, 1UL );
      checkNonZeros( I, 3UL, 1UL );
      checkNonZeros( I, 4UL, 1UL );
      checkNonZeros( I, 5UL, 1UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c reset() member function of the IdentityMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c reset() member function of the IdentityMatrix class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testReset()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major IdentityMatrix::reset()";

      // Resetting a default constructed matrix
      {
         blaze::IdentityMatrix<int,blaze::rowMajor> I;

         reset( I );

         checkRows    ( I, 0UL );
         checkColumns ( I, 0UL );
         checkNonZeros( I, 0UL );
      }

      // Resetting an initialized matrix
      {
         // Initialization check
         blaze::IdentityMatrix<int,blaze::rowMajor> I( 4UL );

         checkRows    ( I, 4UL );
         checkColumns ( I, 4UL );
         checkNonZeros( I, 4UL );
         checkNonZeros( I, 0UL, 1UL );
         checkNonZeros( I, 1UL, 1UL );
         checkNonZeros( I, 2UL, 1UL );
         checkNonZeros( I, 3UL, 1UL );

         if( I(0,0) != 1 || I(0,1) != 0 || I(0,2) != 0 || I(0,3) != 0 ||
             I(1,0) != 0 || I(1,1) != 1 || I(1,2) != 0 || I(1,3) != 0 ||
             I(2,0) != 0 || I(2,1) != 0 || I(2,2) != 1 || I(2,3) != 0 ||
             I(3,0) != 0 || I(3,1) != 0 || I(3,2) != 0 || I(3,3) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << I << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 0 1 0 0 )\n( 0 0 1 0 )\n( 0 0 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Resetting row 1
         reset( I, 1UL );

         checkRows    ( I, 4UL );
         checkColumns ( I, 4UL );
         checkNonZeros( I, 4UL );
         checkNonZeros( I, 0UL, 1UL );
         checkNonZeros( I, 1UL, 1UL );
         checkNonZeros( I, 2UL, 1UL );
         checkNonZeros( I, 3UL, 1UL );

         if( I(0,0) != 1 || I(0,1) != 0 || I(0,2) != 0 || I(0,3) != 0 ||
             I(1,0) != 0 || I(1,1) != 1 || I(1,2) != 0 || I(1,3) != 0 ||
             I(2,0) != 0 || I(2,1) != 0 || I(2,2) != 1 || I(2,3) != 0 ||
             I(3,0) != 0 || I(3,1) != 0 || I(3,2) != 0 || I(3,3) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Reset operation failed\n"
                << " Details:\n"
                << "   Result:\n" << I << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 0 1 0 0 )\n( 0 0 1 0 )\n( 0 0 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Resetting the entire matrix
         reset( I );

         checkRows    ( I, 4UL );
         checkColumns ( I, 4UL );
         checkNonZeros( I, 4UL );
         checkNonZeros( I, 0UL, 1UL );
         checkNonZeros( I, 1UL, 1UL );
         checkNonZeros( I, 2UL, 1UL );
         checkNonZeros( I, 3UL, 1UL );

         if( I(0,0) != 1 || I(0,1) != 0 || I(0,2) != 0 || I(0,3) != 0 ||
             I(1,0) != 0 || I(1,1) != 1 || I(1,2) != 0 || I(1,3) != 0 ||
             I(2,0) != 0 || I(2,1) != 0 || I(2,2) != 1 || I(2,3) != 0 ||
             I(3,0) != 0 || I(3,1) != 0 || I(3,2) != 0 || I(3,3) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Reset operation failed\n"
                << " Details:\n"
                << "   Result:\n" << I << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 0 1 0 0 )\n( 0 0 1 0 )\n( 0 0 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major IdentityMatrix::reset()";

      // Resetting a default constructed matrix
      {
         blaze::IdentityMatrix<int,blaze::columnMajor> I;

         reset( I );

         checkRows    ( I, 0UL );
         checkColumns ( I, 0UL );
         checkNonZeros( I, 0UL );
      }

      // Resetting an initialized matrix
      {
         // Initialization check
         blaze::IdentityMatrix<int,blaze::columnMajor> I( 4UL );

         checkRows    ( I, 4UL );
         checkColumns ( I, 4UL );
         checkNonZeros( I, 4UL );
         checkNonZeros( I, 0UL, 1UL );
         checkNonZeros( I, 1UL, 1UL );
         checkNonZeros( I, 2UL, 1UL );
         checkNonZeros( I, 3UL, 1UL );

         if( I(0,0) != 1 || I(0,1) != 0 || I(0,2) != 0 || I(0,3) != 0 ||
             I(1,0) != 0 || I(1,1) != 1 || I(1,2) != 0 || I(1,3) != 0 ||
             I(2,0) != 0 || I(2,1) != 0 || I(2,2) != 1 || I(2,3) != 0 ||
             I(3,0) != 0 || I(3,1) != 0 || I(3,2) != 0 || I(3,3) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << I << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 0 1 0 0 )\n( 0 0 1 0 )\n( 0 0 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Resetting column 1
         reset( I, 1UL );

         checkRows    ( I, 4UL );
         checkColumns ( I, 4UL );
         checkNonZeros( I, 4UL );
         checkNonZeros( I, 0UL, 1UL );
         checkNonZeros( I, 1UL, 1UL );
         checkNonZeros( I, 2UL, 1UL );
         checkNonZeros( I, 3UL, 1UL );

         if( I(0,0) != 1 || I(0,1) != 0 || I(0,2) != 0 || I(0,3) != 0 ||
             I(1,0) != 0 || I(1,1) != 1 || I(1,2) != 0 || I(1,3) != 0 ||
             I(2,0) != 0 || I(2,1) != 0 || I(2,2) != 1 || I(2,3) != 0 ||
             I(3,0) != 0 || I(3,1) != 0 || I(3,2) != 0 || I(3,3) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Reset operation failed\n"
                << " Details:\n"
                << "   Result:\n" << I << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 0 1 0 0 )\n( 0 0 1 0 )\n( 0 0 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Resetting the entire matrix
         reset( I );

         checkRows    ( I, 4UL );
         checkColumns ( I, 4UL );
         checkNonZeros( I, 4UL );
         checkNonZeros( I, 0UL, 1UL );
         checkNonZeros( I, 1UL, 1UL );
         checkNonZeros( I, 2UL, 1UL );
         checkNonZeros( I, 3UL, 1UL );

         if( I(0,0) != 1 || I(0,1) != 0 || I(0,2) != 0 || I(0,3) != 0 ||
             I(1,0) != 0 || I(1,1) != 1 || I(1,2) != 0 || I(1,3) != 0 ||
             I(2,0) != 0 || I(2,1) != 0 || I(2,2) != 1 || I(2,3) != 0 ||
             I(3,0) != 0 || I(3,1) != 0 || I(3,2) != 0 || I(3,3) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Reset operation failed\n"
                << " Details:\n"
                << "   Result:\n" << I << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 0 1 0 0 )\n( 0 0 1 0 )\n( 0 0 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c clear() member function of the IdentityMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c clear() member function of the IdentityMatrix class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testClear()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major IdentityMatrix::clear()";

      // Clearing a default constructed matrix
      {
         blaze::IdentityMatrix<int,blaze::rowMajor> I;

         clear( I );

         checkRows    ( I, 0UL );
         checkColumns ( I, 0UL );
         checkNonZeros( I, 0UL );
      }

      // Clearing an initialized matrix
      {
         // Initialization check
         blaze::IdentityMatrix<int,blaze::rowMajor> I( 4UL );

         checkRows    ( I, 4UL );
         checkColumns ( I, 4UL );
         checkNonZeros( I, 4UL );
         checkNonZeros( I, 0UL, 1UL );
         checkNonZeros( I, 1UL, 1UL );
         checkNonZeros( I, 2UL, 1UL );
         checkNonZeros( I, 3UL, 1UL );

         if( I(0,0) != 1 || I(0,1) != 0 || I(0,2) != 0 || I(0,3) != 0 ||
             I(1,0) != 0 || I(1,1) != 1 || I(1,2) != 0 || I(1,3) != 0 ||
             I(2,0) != 0 || I(2,1) != 0 || I(2,2) != 1 || I(2,3) != 0 ||
             I(3,0) != 0 || I(3,1) != 0 || I(3,2) != 0 || I(3,3) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << I << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 0 1 0 0 )\n( 0 0 1 0 )\n( 0 0 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Clearing the matrix
         clear( I );

         checkRows    ( I, 0UL );
         checkColumns ( I, 0UL );
         checkNonZeros( I, 0UL );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major IdentityMatrix::clear()";

      // Clearing a default constructed matrix
      {
         blaze::IdentityMatrix<int,blaze::columnMajor> I;

         clear( I );

         checkRows    ( I, 0UL );
         checkColumns ( I, 0UL );
         checkNonZeros( I, 0UL );
      }

      // Clearing an initialized matrix
      {
         // Initialization check
         blaze::IdentityMatrix<int,blaze::columnMajor> I( 4UL );

         checkRows    ( I, 4UL );
         checkColumns ( I, 4UL );
         checkNonZeros( I, 4UL );
         checkNonZeros( I, 0UL, 1UL );
         checkNonZeros( I, 1UL, 1UL );
         checkNonZeros( I, 2UL, 1UL );
         checkNonZeros( I, 3UL, 1UL );

         if( I(0,0) != 1 || I(0,1) != 0 || I(0,2) != 0 || I(0,3) != 0 ||
             I(1,0) != 0 || I(1,1) != 1 || I(1,2) != 0 || I(1,3) != 0 ||
             I(2,0) != 0 || I(2,1) != 0 || I(2,2) != 1 || I(2,3) != 0 ||
             I(3,0) != 0 || I(3,1) != 0 || I(3,2) != 0 || I(3,3) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Clear operation failed\n"
                << " Details:\n"
                << "   Result:\n" << I << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 0 1 0 0 )\n( 0 0 1 0 )\n( 0 0 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Clearing the matrix
         clear( I );

         checkRows    ( I, 0UL );
         checkColumns ( I, 0UL );
         checkNonZeros( I, 0UL );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c resize() member function of the IdentityMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c resize() member function of the IdentityMatrix class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testResize()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major IdentityMatrix::resize()";

      // Initialization check
      blaze::IdentityMatrix<int,blaze::rowMajor> I;

      checkRows    ( I, 0UL );
      checkColumns ( I, 0UL );
      checkNonZeros( I, 0UL );

      // Resizing to 4x4
      I.resize( 4UL );

      checkRows    ( I, 4UL );
      checkColumns ( I, 4UL );
      checkNonZeros( I, 4UL );
      checkNonZeros( I, 0UL, 1UL );
      checkNonZeros( I, 1UL, 1UL );
      checkNonZeros( I, 2UL, 1UL );
      checkNonZeros( I, 3UL, 1UL );

      // Resizing to 2x2
      I.resize( 2UL );

      checkRows    ( I, 2UL );
      checkColumns ( I, 2UL );
      checkNonZeros( I, 2UL );
      checkNonZeros( I, 0UL, 1UL );
      checkNonZeros( I, 1UL, 1UL );

      // Resizing to 0x0
      I.resize( 0UL );

      checkRows    ( I, 0UL );
      checkColumns ( I, 0UL );
      checkNonZeros( I, 0UL );
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major IdentityMatrix::resize()";

      // Initialization check
      blaze::IdentityMatrix<int,blaze::columnMajor> I;

      checkRows    ( I, 0UL );
      checkColumns ( I, 0UL );
      checkNonZeros( I, 0UL );

      // Resizing to 4x4
      I.resize( 4UL );

      checkRows    ( I, 4UL );
      checkColumns ( I, 4UL );
      checkNonZeros( I, 4UL );
      checkNonZeros( I, 0UL, 1UL );
      checkNonZeros( I, 1UL, 1UL );
      checkNonZeros( I, 2UL, 1UL );
      checkNonZeros( I, 3UL, 1UL );

      // Resizing to 2x2
      I.resize( 2UL );

      checkRows    ( I, 2UL );
      checkColumns ( I, 2UL );
      checkNonZeros( I, 2UL );
      checkNonZeros( I, 0UL, 1UL );
      checkNonZeros( I, 1UL, 1UL );

      // Resizing to 0x0
      I.resize( 0UL );

      checkRows    ( I, 0UL );
      checkColumns ( I, 0UL );
      checkNonZeros( I, 0UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c swap() functionality of the IdentityMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c swap() function of the IdentityMatrix class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testSwap()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major IdentityMatrix swap";

      blaze::IdentityMatrix<int,blaze::rowMajor> I1( 4UL );
      blaze::IdentityMatrix<int,blaze::rowMajor> I2( 2UL );

      swap( I1, I2 );

      checkRows    ( I1, 2UL );
      checkColumns ( I1, 2UL );
      checkNonZeros( I1, 2UL );
      checkNonZeros( I1, 0UL, 1UL );
      checkNonZeros( I1, 1UL, 1UL );

      if( I1(0,0) != 1 || I1(0,1) != 0 ||
          I1(1,0) != 0 || I1(1,1) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the first matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << I1 << "\n"
             << "   Expected result:\n( 1 0 )\n( 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      checkRows    ( I2, 4UL );
      checkColumns ( I2, 4UL );
      checkNonZeros( I2, 4UL );
      checkNonZeros( I2, 0UL, 1UL );
      checkNonZeros( I2, 1UL, 1UL );
      checkNonZeros( I2, 2UL, 1UL );
      checkNonZeros( I2, 3UL, 1UL );

      if( I2(0,0) != 1 || I2(0,1) != 0 || I2(0,2) != 0 || I2(0,3) != 0 ||
          I2(1,0) != 0 || I2(1,1) != 1 || I2(1,2) != 0 || I2(1,3) != 0 ||
          I2(2,0) != 0 || I2(2,1) != 0 || I2(2,2) != 1 || I2(2,3) != 0 ||
          I2(3,0) != 0 || I2(3,1) != 0 || I2(3,2) != 0 || I2(3,3) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the second matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << I2 << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 0 1 0 0 )\n( 0 0 1 0 )\n( 0 0 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major IdentityMatrix swap";

      blaze::IdentityMatrix<int,blaze::columnMajor> I1( 4UL );
      blaze::IdentityMatrix<int,blaze::columnMajor> I2( 2UL );

      swap( I1, I2 );

      checkRows    ( I1, 2UL );
      checkColumns ( I1, 2UL );
      checkNonZeros( I1, 2UL );
      checkNonZeros( I1, 0UL, 1UL );
      checkNonZeros( I1, 1UL, 1UL );

      if( I1(0,0) != 1 || I1(0,1) != 0 ||
          I1(1,0) != 0 || I1(1,1) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the first matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << I1 << "\n"
             << "   Expected result:\n( 1 0 )\n( 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      checkRows    ( I2, 4UL );
      checkColumns ( I2, 4UL );
      checkNonZeros( I2, 4UL );
      checkNonZeros( I2, 0UL, 1UL );
      checkNonZeros( I2, 1UL, 1UL );
      checkNonZeros( I2, 2UL, 1UL );
      checkNonZeros( I2, 3UL, 1UL );

      if( I2(0,0) != 1 || I2(0,1) != 0 || I2(0,2) != 0 || I2(0,3) != 0 ||
          I2(1,0) != 0 || I2(1,1) != 1 || I2(1,2) != 0 || I2(1,3) != 0 ||
          I2(2,0) != 0 || I2(2,1) != 0 || I2(2,2) != 1 || I2(2,3) != 0 ||
          I2(3,0) != 0 || I2(3,1) != 0 || I2(3,2) != 0 || I2(3,3) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the second matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << I2 << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 0 1 0 0 )\n( 0 0 1 0 )\n( 0 0 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c find() member function of the IdentityMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c find() member function of the IdentityMatrix class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testFind()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major IdentityMatrix::find()";

      using ConstIterator = blaze::IdentityMatrix<int,blaze::rowMajor>::ConstIterator;

      // Initialization check
      blaze::IdentityMatrix<int,blaze::rowMajor> I( 8UL );

      checkRows    ( I, 8UL );
      checkColumns ( I, 8UL );
      checkNonZeros( I, 8UL );
      checkNonZeros( I, 0UL, 1UL );
      checkNonZeros( I, 1UL, 1UL );
      checkNonZeros( I, 2UL, 1UL );
      checkNonZeros( I, 3UL, 1UL );
      checkNonZeros( I, 4UL, 1UL );
      checkNonZeros( I, 5UL, 1UL );
      checkNonZeros( I, 6UL, 1UL );
      checkNonZeros( I, 7UL, 1UL );

      // Searching for the first element
      {
         ConstIterator pos( I.find( 0UL, 0UL ) );

         if( pos == I.end( 0UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required position = (1,1)\n"
                << "   Current matrix:\n" << I << "\n";
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
                << "   Current matrix:\n" << I << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for the second element
      {
         ConstIterator pos( I.find( 4UL, 4UL ) );

         if( pos == I.end( 4UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required position = (4,4)\n"
                << "   Current matrix:\n" << I << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 4 || pos->value() != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 4\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 1\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << I << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for the third element
      {
         ConstIterator pos( I.find( 7UL, 7UL ) );

         if( pos == I.end( 7UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required position = (7,7)\n"
                << "   Current matrix:\n" << I << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 7 || pos->value() != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 7\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 1\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << I << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for a non-existing non-zero element
      {
         ConstIterator pos( I.find( 4UL, 0UL ) );

         if( pos != I.end( 4UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Non-existing element could be found\n"
                << " Details:\n"
                << "   Required index = 4\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 0\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << I << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major IdentityMatrix::find()";

      using ConstIterator = blaze::IdentityMatrix<int,blaze::columnMajor>::ConstIterator;

      // Initialization check
      blaze::IdentityMatrix<int,blaze::columnMajor> I( 8UL );

      checkRows    ( I, 8UL );
      checkColumns ( I, 8UL );
      checkNonZeros( I, 8UL );
      checkNonZeros( I, 0UL, 1UL );
      checkNonZeros( I, 1UL, 1UL );
      checkNonZeros( I, 2UL, 1UL );
      checkNonZeros( I, 3UL, 1UL );
      checkNonZeros( I, 4UL, 1UL );
      checkNonZeros( I, 5UL, 1UL );
      checkNonZeros( I, 6UL, 1UL );
      checkNonZeros( I, 7UL, 1UL );

      // Searching for the first element
      {
         ConstIterator pos( I.find( 0UL, 0UL ) );

         if( pos == I.end( 0UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required position = (1,1)\n"
                << "   Current matrix:\n" << I << "\n";
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
                << "   Current matrix:\n" << I << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for the second element
      {
         ConstIterator pos( I.find( 4UL, 4UL ) );

         if( pos == I.end( 4UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required position = (4,4)\n"
                << "   Current matrix:\n" << I << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 4 || pos->value() != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 4\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 1\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << I << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for the third element
      {
         ConstIterator pos( I.find( 7UL, 7UL ) );

         if( pos == I.end( 7UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required position = (7,7)\n"
                << "   Current matrix:\n" << I << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 7 || pos->value() != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 7\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 1\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << I << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for a non-existing non-zero element
      {
         ConstIterator pos( I.find( 4UL, 0UL ) );

         if( pos != I.end( 0UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Non-existing element could be found\n"
                << " Details:\n"
                << "   Required index = 4\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 0\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << I << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c lowerBound() member function of the IdentityMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c lowerBound() member function of the IdentityMatrix
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testLowerBound()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major IdentityMatrix::lowerBound()";

      using ConstIterator = blaze::IdentityMatrix<int,blaze::rowMajor>::ConstIterator;

      // Initialization check
      blaze::IdentityMatrix<int,blaze::rowMajor> I( 3UL );

      checkRows    ( I, 3UL );
      checkColumns ( I, 3UL );
      checkNonZeros( I, 3UL );
      checkNonZeros( I, 0UL, 1UL );
      checkNonZeros( I, 1UL, 1UL );
      checkNonZeros( I, 2UL, 1UL );

      // Determining the lower bound for position (1,0)
      {
         ConstIterator pos( I.lowerBound( 1UL, 0UL ) );

         if( pos == I.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,0)\n"
                << "   Current matrix:\n" << I << "\n";
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
                << "   Current matrix:\n" << I << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for position (1,1)
      {
         ConstIterator pos( I.lowerBound( 1UL, 1UL ) );

         if( pos == I.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,1)\n"
                << "   Current matrix:\n" << I << "\n";
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
                << "   Current matrix:\n" << I << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for position (1,2)
      {
         ConstIterator pos( I.lowerBound( 1UL, 2UL ) );

         if( pos != I.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,2)\n"
                << "   Current matrix:\n" << I << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major IdentityMatrix::lowerBound()";

      using ConstIterator = blaze::IdentityMatrix<int,blaze::columnMajor>::ConstIterator;

      // Initialization check
      blaze::IdentityMatrix<int,blaze::columnMajor> I( 3UL );

      checkRows    ( I, 3UL );
      checkColumns ( I, 3UL );
      checkNonZeros( I, 3UL );
      checkNonZeros( I, 0UL, 1UL );
      checkNonZeros( I, 1UL, 1UL );
      checkNonZeros( I, 2UL, 1UL );

      // Determining the lower bound for position (0,1)
      {
         ConstIterator pos( I.lowerBound( 0UL, 1UL ) );

         if( pos == I.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (0,1)\n"
                << "   Current matrix:\n" << I << "\n";
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
                << "   Current matrix:\n" << I << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for position (1,1)
      {
         ConstIterator pos( I.lowerBound( 1UL, 1UL ) );

         if( pos == I.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,1)\n"
                << "   Current matrix:\n" << I << "\n";
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
                << "   Current matrix:\n" << I << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for position (2,1)
      {
         ConstIterator pos( I.lowerBound( 2UL, 1UL ) );

         if( pos != I.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (2,1)\n"
                << "   Current matrix:\n" << I << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c upperBound() member function of the IdentityMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c upperBound() member function of the IdentityMatrix
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testUpperBound()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major IdentityMatrix::upperBound()";

      using ConstIterator = blaze::IdentityMatrix<int,blaze::rowMajor>::ConstIterator;

      // Initialization check
      blaze::IdentityMatrix<int,blaze::rowMajor> I( 3UL );

      checkRows    ( I, 3UL );
      checkColumns ( I, 3UL );
      checkNonZeros( I, 3UL );
      checkNonZeros( I, 0UL, 1UL );
      checkNonZeros( I, 1UL, 1UL );
      checkNonZeros( I, 2UL, 1UL );

      // Determining the upper bound for position (1,0)
      {
         ConstIterator pos( I.upperBound( 1UL, 0UL ) );

         if( pos == I.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,0)\n"
                << "   Current matrix:\n" << I << "\n";
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
                << "   Current matrix:\n" << I << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for position (1,1)
      {
         ConstIterator pos( I.upperBound( 1UL, 1UL ) );

         if( pos != I.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,1)\n"
                << "   Current matrix:\n" << I << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for position (1,2)
      {
         ConstIterator pos( I.upperBound( 1UL, 2UL ) );

         if( pos != I.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,2)\n"
                << "   Current matrix:\n" << I << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major IdentityMatrix::upperBound()";

      using ConstIterator = blaze::IdentityMatrix<int,blaze::columnMajor>::ConstIterator;

      // Initialization check
      blaze::IdentityMatrix<int,blaze::columnMajor> I( 3UL );

      checkRows    ( I, 3UL );
      checkColumns ( I, 3UL );
      checkNonZeros( I, 3UL );
      checkNonZeros( I, 0UL, 1UL );
      checkNonZeros( I, 1UL, 1UL );
      checkNonZeros( I, 2UL, 1UL );

      // Determining the upper bound for position (0,1)
      {
         ConstIterator pos( I.upperBound( 0UL, 1UL ) );

         if( pos == I.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (0,1)\n"
                << "   Current matrix:\n" << I << "\n";
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
                << "   Current matrix:\n" << I << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for position (1,1)
      {
         ConstIterator pos( I.upperBound( 1UL, 1UL ) );

         if( pos != I.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,1)\n"
                << "   Current matrix:\n" << I << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for position (2,1)
      {
         ConstIterator pos( I.upperBound( 2UL, 1UL ) );

         if( pos != I.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (2,1)\n"
                << "   Current matrix:\n" << I << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c transpose() member function of the IdentityMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c transpose() member function of the IdentityMatrix
// class template. Additionally, it performs a test of self-transpose via the \c trans()
// function. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testTranspose()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major self-transpose via transpose()";

      blaze::IdentityMatrix<int,blaze::rowMajor> I( 4UL );

      transpose( I );

      checkRows    ( I, 4UL );
      checkColumns ( I, 4UL );
      checkNonZeros( I, 4UL );
      checkNonZeros( I, 0UL, 1UL );
      checkNonZeros( I, 1UL, 1UL );
      checkNonZeros( I, 2UL, 1UL );
      checkNonZeros( I, 3UL, 1UL );

      if( I(0,0) != 1 || I(0,1) != 0 || I(0,2) != 0 || I(0,3) != 0 ||
          I(1,0) != 0 || I(1,1) != 1 || I(1,2) != 0 || I(1,3) != 0 ||
          I(2,0) != 0 || I(2,1) != 0 || I(2,2) != 1 || I(2,3) != 0 ||
          I(3,0) != 0 || I(3,1) != 0 || I(3,2) != 0 || I(3,3) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Transposition failed\n"
             << " Details:\n"
             << "   Result:\n" << I << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 0 1 0 0 )\n( 0 0 1 0 )\n( 0 0 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major self-transpose via trans()";

      blaze::IdentityMatrix<int,blaze::rowMajor> I( 4UL );

      I = trans( I );

      checkRows    ( I, 4UL );
      checkColumns ( I, 4UL );
      checkNonZeros( I, 4UL );
      checkNonZeros( I, 0UL, 1UL );
      checkNonZeros( I, 1UL, 1UL );
      checkNonZeros( I, 2UL, 1UL );
      checkNonZeros( I, 3UL, 1UL );

      if( I(0,0) != 1 || I(0,1) != 0 || I(0,2) != 0 || I(0,3) != 0 ||
          I(1,0) != 0 || I(1,1) != 1 || I(1,2) != 0 || I(1,3) != 0 ||
          I(2,0) != 0 || I(2,1) != 0 || I(2,2) != 1 || I(2,3) != 0 ||
          I(3,0) != 0 || I(3,1) != 0 || I(3,2) != 0 || I(3,3) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Transposition failed\n"
             << " Details:\n"
             << "   Result:\n" << I << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 0 1 0 0 )\n( 0 0 1 0 )\n( 0 0 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major self-transpose via transpose()";

      blaze::IdentityMatrix<int,blaze::columnMajor> I( 4UL );

      transpose( I );

      checkRows    ( I, 4UL );
      checkColumns ( I, 4UL );
      checkNonZeros( I, 4UL );
      checkNonZeros( I, 0UL, 1UL );
      checkNonZeros( I, 1UL, 1UL );
      checkNonZeros( I, 2UL, 1UL );
      checkNonZeros( I, 3UL, 1UL );

      if( I(0,0) != 1 || I(0,1) != 0 || I(0,2) != 0 || I(0,3) != 0 ||
          I(1,0) != 0 || I(1,1) != 1 || I(1,2) != 0 || I(1,3) != 0 ||
          I(2,0) != 0 || I(2,1) != 0 || I(2,2) != 1 || I(2,3) != 0 ||
          I(3,0) != 0 || I(3,1) != 0 || I(3,2) != 0 || I(3,3) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Transposition failed\n"
             << " Details:\n"
             << "   Result:\n" << I << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 0 1 0 0 )\n( 0 0 1 0 )\n( 0 0 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major self-transpose via trans()";

      blaze::IdentityMatrix<int,blaze::columnMajor> I( 4UL );

      I = trans( I );

      checkRows    ( I, 4UL );
      checkColumns ( I, 4UL );
      checkNonZeros( I, 4UL );
      checkNonZeros( I, 0UL, 1UL );
      checkNonZeros( I, 1UL, 1UL );
      checkNonZeros( I, 2UL, 1UL );
      checkNonZeros( I, 3UL, 1UL );

      if( I(0,0) != 1 || I(0,1) != 0 || I(0,2) != 0 || I(0,3) != 0 ||
          I(1,0) != 0 || I(1,1) != 1 || I(1,2) != 0 || I(1,3) != 0 ||
          I(2,0) != 0 || I(2,1) != 0 || I(2,2) != 1 || I(2,3) != 0 ||
          I(3,0) != 0 || I(3,1) != 0 || I(3,2) != 0 || I(3,3) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Transposition failed\n"
             << " Details:\n"
             << "   Result:\n" << I << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 0 1 0 0 )\n( 0 0 1 0 )\n( 0 0 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c ctranspose() member function of the IdentityMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c ctranspose() member function of the IdentityMatrix
// class template. Additionally, it performs a test of self-transpose via the \c ctrans()
// function. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testCTranspose()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major self-transpose via ctranspose()";

      using cplx = blaze::complex<int>;

      blaze::IdentityMatrix<cplx,blaze::rowMajor> I( 4UL );

      ctranspose( I );

      checkRows    ( I, 4UL );
      checkColumns ( I, 4UL );
      checkNonZeros( I, 4UL );
      checkNonZeros( I, 0UL, 1UL );
      checkNonZeros( I, 1UL, 1UL );
      checkNonZeros( I, 2UL, 1UL );
      checkNonZeros( I, 3UL, 1UL );

      if( I(0,0) != cplx(1,0) || I(0,1) != cplx(0,0) || I(0,2) != cplx(0,0) || I(0,3) != cplx(0,0) ||
          I(1,0) != cplx(0,0) || I(1,1) != cplx(1,0) || I(1,2) != cplx(0,0) || I(1,3) != cplx(0,0) ||
          I(2,0) != cplx(0,0) || I(2,1) != cplx(0,0) || I(2,2) != cplx(1,0) || I(2,3) != cplx(0,0) ||
          I(3,0) != cplx(0,0) || I(3,1) != cplx(0,0) || I(3,2) != cplx(0,0) || I(3,3) != cplx(1,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Transposition failed\n"
             << " Details:\n"
             << "   Result:\n" << I << "\n"
             << "   Expected result:\n( (1,0) (0,0) (0,0) (0,0) )\n"
                                     "( (0,0) (1,0) (0,0) (0,0) )\n"
                                     "( (0,0) (0,0) (1,0) (0,0) )\n"
                                     "( (0,0) (0,0) (0,0) (1,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major self-transpose via ctranspose()";

      using cplx = blaze::complex<int>;

      blaze::IdentityMatrix<cplx,blaze::rowMajor> I( 4UL );

      I = ctrans( I );

      checkRows    ( I, 4UL );
      checkColumns ( I, 4UL );
      checkNonZeros( I, 4UL );
      checkNonZeros( I, 0UL, 1UL );
      checkNonZeros( I, 1UL, 1UL );
      checkNonZeros( I, 2UL, 1UL );
      checkNonZeros( I, 3UL, 1UL );

      if( I(0,0) != cplx(1,0) || I(0,1) != cplx(0,0) || I(0,2) != cplx(0,0) || I(0,3) != cplx(0,0) ||
          I(1,0) != cplx(0,0) || I(1,1) != cplx(1,0) || I(1,2) != cplx(0,0) || I(1,3) != cplx(0,0) ||
          I(2,0) != cplx(0,0) || I(2,1) != cplx(0,0) || I(2,2) != cplx(1,0) || I(2,3) != cplx(0,0) ||
          I(3,0) != cplx(0,0) || I(3,1) != cplx(0,0) || I(3,2) != cplx(0,0) || I(3,3) != cplx(1,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Transposition failed\n"
             << " Details:\n"
             << "   Result:\n" << I << "\n"
             << "   Expected result:\n( (1,0) (0,0) (0,0) (0,0) )\n"
                                     "( (0,0) (1,0) (0,0) (0,0) )\n"
                                     "( (0,0) (0,0) (1,0) (0,0) )\n"
                                     "( (0,0) (0,0) (0,0) (1,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major self-transpose via ctranspose()";

      using cplx = blaze::complex<int>;

      blaze::IdentityMatrix<cplx,blaze::columnMajor> I( 4UL );

      ctranspose( I );

      checkRows    ( I, 4UL );
      checkColumns ( I, 4UL );
      checkNonZeros( I, 4UL );
      checkNonZeros( I, 0UL, 1UL );
      checkNonZeros( I, 1UL, 1UL );
      checkNonZeros( I, 2UL, 1UL );
      checkNonZeros( I, 2UL, 1UL );

      if( I(0,0) != cplx(1,0) || I(0,1) != cplx(0,0) || I(0,2) != cplx(0,0) || I(0,3) != cplx(0,0) ||
          I(1,0) != cplx(0,0) || I(1,1) != cplx(1,0) || I(1,2) != cplx(0,0) || I(1,3) != cplx(0,0) ||
          I(2,0) != cplx(0,0) || I(2,1) != cplx(0,0) || I(2,2) != cplx(1,0) || I(2,3) != cplx(0,0) ||
          I(3,0) != cplx(0,0) || I(3,1) != cplx(0,0) || I(3,2) != cplx(0,0) || I(3,3) != cplx(1,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Transposition failed\n"
             << " Details:\n"
             << "   Result:\n" << I << "\n"
             << "   Expected result:\n( (1,0) (0,0) (0,0) (0,0) )\n"
                                     "( (0,0) (1,0) (0,0) (0,0) )\n"
                                     "( (0,0) (0,0) (1,0) (0,0) )\n"
                                     "( (0,0) (0,0) (0,0) (1,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major self-transpose via ctranspose()";

      using cplx = blaze::complex<int>;

      blaze::IdentityMatrix<cplx,blaze::columnMajor> I( 4UL );

      I = ctrans( I );

      checkRows    ( I, 4UL );
      checkColumns ( I, 4UL );
      checkNonZeros( I, 4UL );
      checkNonZeros( I, 0UL, 1UL );
      checkNonZeros( I, 1UL, 1UL );
      checkNonZeros( I, 2UL, 1UL );
      checkNonZeros( I, 2UL, 1UL );

      if( I(0,0) != cplx(1,0) || I(0,1) != cplx(0,0) || I(0,2) != cplx(0,0) || I(0,3) != cplx(0,0) ||
          I(1,0) != cplx(0,0) || I(1,1) != cplx(1,0) || I(1,2) != cplx(0,0) || I(1,3) != cplx(0,0) ||
          I(2,0) != cplx(0,0) || I(2,1) != cplx(0,0) || I(2,2) != cplx(1,0) || I(2,3) != cplx(0,0) ||
          I(3,0) != cplx(0,0) || I(3,1) != cplx(0,0) || I(3,2) != cplx(0,0) || I(3,3) != cplx(1,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Transposition failed\n"
             << " Details:\n"
             << "   Result:\n" << I << "\n"
             << "   Expected result:\n( (1,0) (0,0) (0,0) (0,0) )\n"
                                     "( (0,0) (1,0) (0,0) (0,0) )\n"
                                     "( (0,0) (0,0) (1,0) (0,0) )\n"
                                     "( (0,0) (0,0) (0,0) (1,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c isDefault() function with the IdentityMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isDefault() function with the IdentityMatrix class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testIsDefault()
{
   using blaze::isDefault;


   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major isDefault() function";

      // isDefault with 0x0 matrix (default)
      {
         blaze::IdentityMatrix<int,blaze::rowMajor> I;

         if( isDefault( I ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << I << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with 3x3 matrix (non-default)
      {
         blaze::IdentityMatrix<int,blaze::rowMajor> I( 3UL );

         if( isDefault( I(0,1) ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix element: " << I(0,1) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( I(1,1) ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix element: " << I(1,1) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( I ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << I << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major isDefault() function";

      // isDefault with 0x0 matrix (default)
      {
         blaze::IdentityMatrix<int,blaze::columnMajor> I;

         if( isDefault( I ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << I << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with 3x3 matrix (non-default)
      {
         blaze::IdentityMatrix<int,blaze::columnMajor> I( 3UL );

         if( isDefault( I(1,0) ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix element: " << I(1,0) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( I(1,1) ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix element: " << I(1,1) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( I ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << I << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************

} // namespace identitymatrix

} // namespace matrices

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
   std::cout << "   Running IdentityMatrix class test..." << std::endl;

   try
   {
      RUN_IDENTITYMATRIX_CLASS_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during IdentityMatrix class test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
