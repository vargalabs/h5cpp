//=================================================================================================
/*!
//  \file src/mathtest/matrices/zeromatrix/ClassTest.cpp
//  \brief Source file for the ZeroMatrix class test
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
#include <blazetest/mathtest/matrices/zeromatrix/ClassTest.h>

#ifdef BLAZE_USE_HPX_THREADS
#  include <hpx/hpx_main.hpp>
#endif


namespace blazetest {

namespace mathtest {

namespace matrices {

namespace zeromatrix {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the ZeroMatrix class test.
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
/*!\brief Test of the ZeroMatrix constructors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all constructors of the ZeroMatrix class template. In case
// an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testConstructors()
{
   //=====================================================================================
   // Row-major default constructor
   //=====================================================================================

   {
      test_ = "Row-major ZeroMatrix default constructor";

      blaze::ZeroMatrix<int,blaze::rowMajor> Z;

      checkRows    ( Z, 0UL );
      checkColumns ( Z, 0UL );
      checkNonZeros( Z, 0UL );
   }


   //=====================================================================================
   // Row-major size constructor
   //=====================================================================================

   {
      test_ = "Row-major ZeroMatrix size constructor (0x0)";

      blaze::ZeroMatrix<int,blaze::rowMajor> Z( 0UL, 0UL );

      checkRows    ( Z, 0UL );
      checkColumns ( Z, 0UL );
      checkNonZeros( Z, 0UL );
   }

   {
      test_ = "Row-major ZeroMatrix size constructor (3x4)";

      blaze::ZeroMatrix<int,blaze::rowMajor> Z( 3UL, 4UL );

      checkRows    ( Z, 3UL );
      checkColumns ( Z, 4UL );
      checkNonZeros( Z, 0UL );
      checkNonZeros( Z, 0UL, 0UL );
      checkNonZeros( Z, 1UL, 0UL );
      checkNonZeros( Z, 2UL, 0UL );

      if( Z(0,0) != 0 || Z(0,1) != 0 || Z(0,2) != 0 || Z(0,3) != 0 ||
          Z(1,0) != 0 || Z(1,1) != 0 || Z(1,2) != 0 || Z(1,3) != 0 ||
          Z(2,0) != 0 || Z(2,1) != 0 || Z(2,2) != 0 || Z(2,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << Z << "\n"
             << "   Expected result:\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major copy constructor
   //=====================================================================================

   {
      test_ = "Row-major ZeroMatrix copy constructor (0x0)";

      blaze::ZeroMatrix<int,blaze::rowMajor> Z1( 0UL, 0UL );
      blaze::ZeroMatrix<int,blaze::rowMajor> Z2( Z1 );

      checkRows    ( Z2, 0UL );
      checkColumns ( Z2, 0UL );
      checkNonZeros( Z2, 0UL );
   }

   {
      test_ = "Row-major ZeroMatrix copy constructor (3x4)";

      blaze::ZeroMatrix<int,blaze::rowMajor> Z1( 3UL, 4UL );
      blaze::ZeroMatrix<int,blaze::rowMajor> Z2( Z1 );

      checkRows    ( Z2, 3UL );
      checkColumns ( Z2, 4UL );
      checkNonZeros( Z2, 0UL );
      checkNonZeros( Z2, 0UL, 0UL );
      checkNonZeros( Z2, 1UL, 0UL );
      checkNonZeros( Z2, 2UL, 0UL );

      if( Z2(0,0) != 0 || Z2(0,1) != 0 || Z2(0,2) != 0 ||
          Z2(1,0) != 0 || Z2(1,1) != 0 || Z2(1,2) != 0 ||
          Z2(2,0) != 0 || Z2(2,1) != 0 || Z2(2,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << Z2 << "\n"
             << "   Expected result:\n( 0 0 0 )\n( 0 0 0 )\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major move constructor
   //=====================================================================================

   {
      test_ = "Row-major ZeroMatrix move constructor (0x0)";

      blaze::ZeroMatrix<int,blaze::rowMajor> Z1( 0UL, 0UL );
      blaze::ZeroMatrix<int,blaze::rowMajor> Z2( std::move( Z1 ) );

      checkRows    ( Z2, 0UL );
      checkColumns ( Z2, 0UL );
      checkNonZeros( Z2, 0UL );
   }

   {
      test_ = "Row-major ZeroMatrix move constructor (3x4)";

      blaze::ZeroMatrix<int,blaze::rowMajor> Z1( 3UL, 4UL );
      blaze::ZeroMatrix<int,blaze::rowMajor> Z2( std::move( Z1 ) );

      checkRows    ( Z2, 3UL );
      checkColumns ( Z2, 4UL );
      checkNonZeros( Z2, 0UL );
      checkNonZeros( Z2, 0UL, 0UL );
      checkNonZeros( Z2, 1UL, 0UL );
      checkNonZeros( Z2, 2UL, 0UL );

      if( Z2(0,0) != 0 || Z2(0,1) != 0 || Z2(0,2) != 0 ||
          Z2(1,0) != 0 || Z2(1,1) != 0 || Z2(1,2) != 0 ||
          Z2(2,0) != 0 || Z2(2,1) != 0 || Z2(2,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << Z2 << "\n"
             << "   Expected result:\n( 0 0 0 )\n( 0 0 0 )\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major dense matrix constructor
   //=====================================================================================

   {
      test_ = "Row-major/row-major ZeroMatrix dense matrix constructor";

      blaze::DynamicMatrix<int,blaze::rowMajor> Z1{ { 0, 0, 0 }, { 0, 0, 0 } };
      blaze::ZeroMatrix<int,blaze::rowMajor> Z2( Z1 );

      checkRows    ( Z2, 2UL );
      checkColumns ( Z2, 3UL );
      checkNonZeros( Z2, 0UL );
      checkNonZeros( Z2, 0UL, 0UL );
      checkNonZeros( Z2, 1UL, 0UL );

      if( Z2(0,0) != 0 || Z2(0,1) != 0 || Z2(0,2) != 0 ||
          Z2(1,0) != 0 || Z2(1,1) != 0 || Z2(1,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << Z2 << "\n"
             << "   Expected result:\n( 0 0 0 )\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major ZeroMatrix dense matrix constructor";

      blaze::DynamicMatrix<int,blaze::columnMajor> Z1{ { 0, 0, 0 }, { 0, 0, 0 } };
      blaze::ZeroMatrix<int,blaze::rowMajor> Z2( Z1 );

      checkRows    ( Z2, 2UL );
      checkColumns ( Z2, 3UL );
      checkNonZeros( Z2, 0UL );
      checkNonZeros( Z2, 0UL, 0UL );
      checkNonZeros( Z2, 1UL, 0UL );

      if( Z2(0,0) != 0 || Z2(0,1) != 0 || Z2(0,2) != 0 ||
          Z2(1,0) != 0 || Z2(1,1) != 0 || Z2(1,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << Z2 << "\n"
             << "   Expected result:\n( 0 0 0 )\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major ZeroMatrix dense matrix constructor (non-zero)";

      blaze::DynamicMatrix<int,blaze::rowMajor> Z1{ { 0, 0, 0 }, { 0, 1, 0 } };

      try {
         blaze::ZeroMatrix<int,blaze::rowMajor> Z2( Z1 );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of non-zero ZeroMatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << Z2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Row-major sparse matrix constructor
   //=====================================================================================

   {
      test_ = "Row-major/row-major ZeroMatrix sparse matrix constructor";

      blaze::CompressedMatrix<int,blaze::rowMajor> Z1( 2UL, 3UL, 2UL );
      Z1.insert( 0UL, 1UL, 0 );
      Z1.insert( 1UL, 2UL, 0 );

      blaze::ZeroMatrix<int,blaze::rowMajor> Z2( Z1 );

      checkRows    ( Z2, 2UL );
      checkColumns ( Z2, 3UL );
      checkNonZeros( Z2, 0UL );
      checkNonZeros( Z2, 0UL, 0UL );
      checkNonZeros( Z2, 1UL, 0UL );

      if( Z2(0,0) != 0 || Z2(0,1) != 0 || Z2(0,2) != 0 ||
          Z2(1,0) != 0 || Z2(1,1) != 0 || Z2(1,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << Z2 << "\n"
             << "   Expected result:\n( 0 0 0 )\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major ZeroMatrix sparse matrix constructor";

      blaze::CompressedMatrix<int,blaze::columnMajor> Z1( 2UL, 3UL, 2UL );
      Z1.insert( 0UL, 1UL, 0 );
      Z1.insert( 1UL, 2UL, 0 );

      blaze::ZeroMatrix<int,blaze::rowMajor> Z2( Z1 );

      checkRows    ( Z2, 2UL );
      checkColumns ( Z2, 3UL );
      checkNonZeros( Z2, 0UL );
      checkNonZeros( Z2, 0UL, 0UL );
      checkNonZeros( Z2, 1UL, 0UL );

      if( Z2(0,0) != 0 || Z2(0,1) != 0 || Z2(0,2) != 0 ||
          Z2(1,0) != 0 || Z2(1,1) != 0 || Z2(1,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << Z2 << "\n"
             << "   Expected result:\n( 0 0 0 )\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major ZeroMatrix sparse matrix constructor (non-zero)";

      blaze::CompressedMatrix<int,blaze::rowMajor> Z1{ { 0, 0, 0 }, { 0, 1, 0 } };

      try {
         blaze::ZeroMatrix<int,blaze::rowMajor> Z2( Z1 );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of non-zero ZeroMatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << Z2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Column-major default constructor
   //=====================================================================================

   {
      test_ = "Column-major ZeroMatrix default constructor";

      blaze::ZeroMatrix<int,blaze::columnMajor> Z;

      checkRows    ( Z, 0UL );
      checkColumns ( Z, 0UL );
      checkNonZeros( Z, 0UL );
   }


   //=====================================================================================
   // Column-major size constructor
   //=====================================================================================

   {
      test_ = "Column-major ZeroMatrix size constructor (0x0)";

      blaze::ZeroMatrix<int,blaze::columnMajor> Z( 0UL, 0UL );

      checkRows    ( Z, 0UL );
      checkColumns ( Z, 0UL );
      checkNonZeros( Z, 0UL );
   }

   {
      test_ = "Column-major ZeroMatrix size constructor (4x3)";

      blaze::ZeroMatrix<int,blaze::columnMajor> Z( 4UL, 3UL );

      checkRows    ( Z, 4UL );
      checkColumns ( Z, 3UL );
      checkNonZeros( Z, 0UL );
      checkNonZeros( Z, 0UL, 0UL );
      checkNonZeros( Z, 1UL, 0UL );
      checkNonZeros( Z, 2UL, 0UL );

      if( Z(0,0) != 0 || Z(0,1) != 0 || Z(0,2) != 0 ||
          Z(1,0) != 0 || Z(1,1) != 0 || Z(1,2) != 0 ||
          Z(2,0) != 0 || Z(2,1) != 0 || Z(2,2) != 0 ||
          Z(3,0) != 0 || Z(3,1) != 0 || Z(3,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << Z << "\n"
             << "   Expected result:\n( 0 0 0 )\n( 0 0 0 )\n( 0 0 0 )\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major copy constructor
   //=====================================================================================

   {
      test_ = "Column-major ZeroMatrix copy constructor (0x0)";

      blaze::ZeroMatrix<int,blaze::columnMajor> Z1( 0UL, 0UL );
      blaze::ZeroMatrix<int,blaze::columnMajor> Z2( Z1 );

      checkRows    ( Z2, 0UL );
      checkColumns ( Z2, 0UL );
      checkNonZeros( Z2, 0UL );
   }

   {
      test_ = "Column-major ZeroMatrix copy constructor (4x3)";

      blaze::ZeroMatrix<int,blaze::columnMajor> Z1( 4UL, 3UL );
      blaze::ZeroMatrix<int,blaze::columnMajor> Z2( Z1 );

      checkRows    ( Z2, 4UL );
      checkColumns ( Z2, 3UL );
      checkNonZeros( Z2, 0UL );
      checkNonZeros( Z2, 0UL, 0UL );
      checkNonZeros( Z2, 1UL, 0UL );
      checkNonZeros( Z2, 2UL, 0UL );

      if( Z2(0,0) != 0 || Z2(0,1) != 0 || Z2(0,2) != 0 ||
          Z2(1,0) != 0 || Z2(1,1) != 0 || Z2(1,2) != 0 ||
          Z2(2,0) != 0 || Z2(2,1) != 0 || Z2(2,2) != 0 ||
          Z2(3,0) != 0 || Z2(3,1) != 0 || Z2(3,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << Z2 << "\n"
             << "   Expected result:\n( 0 0 0 )\n( 0 0 0 )\n( 0 0 0 )\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major move constructor
   //=====================================================================================

   {
      test_ = "Column-major ZeroMatrix move constructor (0x0)";

      blaze::ZeroMatrix<int,blaze::columnMajor> Z1( 0UL, 0UL );
      blaze::ZeroMatrix<int,blaze::columnMajor> Z2( std::move( Z1 ) );

      checkRows    ( Z2, 0UL );
      checkColumns ( Z2, 0UL );
      checkNonZeros( Z2, 0UL );
   }

   {
      test_ = "Column-major ZeroMatrix move constructor (4x3)";

      blaze::ZeroMatrix<int,blaze::columnMajor> Z1( 4UL, 3UL );
      blaze::ZeroMatrix<int,blaze::columnMajor> Z2( std::move( Z1 ) );

      checkRows    ( Z2, 4UL );
      checkColumns ( Z2, 3UL );
      checkNonZeros( Z2, 0UL );
      checkNonZeros( Z2, 0UL, 0UL );
      checkNonZeros( Z2, 1UL, 0UL );
      checkNonZeros( Z2, 2UL, 0UL );

      if( Z2(0,0) != 0 || Z2(0,1) != 0 || Z2(0,2) != 0 ||
          Z2(1,0) != 0 || Z2(1,1) != 0 || Z2(1,2) != 0 ||
          Z2(2,0) != 0 || Z2(2,1) != 0 || Z2(2,2) != 0 ||
          Z2(3,0) != 0 || Z2(3,1) != 0 || Z2(3,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << Z2 << "\n"
             << "   Expected result:\n( 0 0 0 )\n( 0 0 0 )\n( 0 0 0 )\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix constructor
   //=====================================================================================

   {
      test_ = "Column-major/row-major ZeroMatrix dense matrix constructor";

      blaze::DynamicMatrix<int,blaze::rowMajor> Z1{ { 0, 0 }, { 0, 0 }, { 0, 0 } };
      blaze::ZeroMatrix<int,blaze::columnMajor> Z2( Z1 );

      checkRows    ( Z2, 3UL );
      checkColumns ( Z2, 2UL );
      checkNonZeros( Z2, 0UL );
      checkNonZeros( Z2, 0UL, 0UL );
      checkNonZeros( Z2, 1UL, 0UL );

      if( Z2(0,0) != 0 || Z2(0,1) != 0 ||
          Z2(1,0) != 0 || Z2(1,1) != 0 ||
          Z2(2,0) != 0 || Z2(2,1) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << Z2 << "\n"
             << "   Expected result:\n( 0 0 )\n( 0 0 )\n( 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major ZeroMatrix dense matrix constructor";

      blaze::DynamicMatrix<int,blaze::columnMajor> Z1{ { 0, 0 }, { 0, 0 }, { 0, 0 } };
      blaze::ZeroMatrix<int,blaze::columnMajor> Z2( Z1 );

      checkRows    ( Z2, 3UL );
      checkColumns ( Z2, 2UL );
      checkNonZeros( Z2, 0UL );
      checkNonZeros( Z2, 0UL, 0UL );
      checkNonZeros( Z2, 1UL, 0UL );

      if( Z2(0,0) != 0 || Z2(0,1) != 0 ||
          Z2(1,0) != 0 || Z2(1,1) != 0 ||
          Z2(2,0) != 0 || Z2(2,1) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << Z2 << "\n"
             << "   Expected result:\n( 0 0 )\n( 0 0 )\n( 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major ZeroMatrix dense matrix constructor (non-zero)";

      blaze::DynamicMatrix<int,blaze::columnMajor> Z1{ { 0, 0 }, { 0, 1 }, { 0, 0 } };

      try {
         blaze::ZeroMatrix<int,blaze::columnMajor> Z2( Z1 );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of non-zero ZeroMatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << Z2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Column-major sparse matrix constructor
   //=====================================================================================

   {
      test_ = "Column-major/row-major ZeroMatrix sparse matrix constructor";

      blaze::CompressedMatrix<int,blaze::rowMajor> Z1( 3UL, 2UL, 2UL );
      Z1.insert( 1UL, 0UL, 0 );
      Z1.insert( 2UL, 1UL, 0 );

      blaze::ZeroMatrix<int,blaze::columnMajor> Z2( Z1 );

      checkRows    ( Z2, 3UL );
      checkColumns ( Z2, 2UL );
      checkNonZeros( Z2, 0UL );
      checkNonZeros( Z2, 0UL, 0UL );
      checkNonZeros( Z2, 1UL, 0UL );

      if( Z2(0,0) != 0 || Z2(0,1) != 0 ||
          Z2(1,0) != 0 || Z2(1,1) != 0 ||
          Z2(2,0) != 0 || Z2(2,1) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << Z2 << "\n"
             << "   Expected result:\n( 0 0 )\n( 0 0 )\n( 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major ZeroMatrix sparse matrix constructor";

      blaze::CompressedMatrix<int,blaze::columnMajor> Z1( 3UL, 2UL, 2UL );
      Z1.insert( 1UL, 0UL, 0 );
      Z1.insert( 2UL, 1UL, 0 );

      blaze::ZeroMatrix<int,blaze::columnMajor> Z2( Z1 );

      checkRows    ( Z2, 3UL );
      checkColumns ( Z2, 2UL );
      checkNonZeros( Z2, 0UL );
      checkNonZeros( Z2, 0UL, 0UL );
      checkNonZeros( Z2, 1UL, 0UL );

      if( Z2(0,0) != 0 || Z2(0,1) != 0 ||
          Z2(1,0) != 0 || Z2(1,1) != 0 ||
          Z2(2,0) != 0 || Z2(2,1) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << Z2 << "\n"
             << "   Expected result:\n( 0 0 )\n( 0 0 )\n( 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major ZeroMatrix sparse matrix constructor (non-zero)";

      blaze::CompressedMatrix<int,blaze::columnMajor> Z1{ { 0, 0 }, { 0, 1 }, { 0, 0 } };

      try {
         blaze::ZeroMatrix<int,blaze::columnMajor> Z2( Z1 );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of non-zero ZeroMatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << Z2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the ZeroMatrix assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all assignment operators of the ZeroMatrix class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testAssignment()
{
   //=====================================================================================
   // Row-major copy assignment
   //=====================================================================================

   {
      test_ = "Row-major ZeroMatrix copy assignment";

      blaze::ZeroMatrix<int,blaze::rowMajor> Z1( 3UL, 4UL );
      blaze::ZeroMatrix<int,blaze::rowMajor> Z2;
      Z2 = Z1;

      checkRows    ( Z2, 3UL );
      checkColumns ( Z2, 4UL );
      checkNonZeros( Z2, 0UL );
      checkNonZeros( Z2, 0UL, 0UL );
      checkNonZeros( Z2, 1UL, 0UL );
      checkNonZeros( Z2, 2UL, 0UL );

      if( Z2(0,0) != 0 || Z2(0,1) != 0 || Z2(0,2) != 0 || Z2(0,3) != 0 ||
          Z2(1,0) != 0 || Z2(1,1) != 0 || Z2(1,2) != 0 || Z2(1,3) != 0 ||
          Z2(2,0) != 0 || Z2(2,1) != 0 || Z2(2,2) != 0 || Z2(2,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << Z2 << "\n"
             << "   Expected result:\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major ZeroMatrix copy assignment stress test";

      blaze::ZeroMatrix<int,blaze::rowMajor> Z1;

      for( size_t i=0UL; i<100UL; ++i )
      {
         const size_t m( blaze::rand<size_t>( 0UL, 10UL ) );
         const size_t n( blaze::rand<size_t>( 0UL, 10UL ) );
         const blaze::ZeroMatrix<int,blaze::rowMajor> Z2( m, n );

         Z1 = Z2;

         if( Z1 != Z2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << Z1 << "\n"
                << "   Expected result:\n" << Z2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Row-major move assignment
   //=====================================================================================

   {
      test_ = "Row-major ZeroMatrix move assignment";

      blaze::ZeroMatrix<int,blaze::rowMajor> Z1( 3UL, 4UL );
      blaze::ZeroMatrix<int,blaze::rowMajor> Z2;

      Z2 = std::move( Z1 );

      checkRows    ( Z2, 3UL );
      checkColumns ( Z2, 4UL );
      checkNonZeros( Z2, 0UL );
      checkNonZeros( Z2, 0UL, 0UL );
      checkNonZeros( Z2, 1UL, 0UL );
      checkNonZeros( Z2, 2UL, 0UL );

      if( Z2(0,0) != 0 || Z2(0,1) != 0 || Z2(0,2) != 0 || Z2(0,3) != 0 ||
          Z2(1,0) != 0 || Z2(1,1) != 0 || Z2(1,2) != 0 || Z2(1,3) != 0 ||
          Z2(2,0) != 0 || Z2(2,1) != 0 || Z2(2,2) != 0 || Z2(2,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << Z2 << "\n"
             << "   Expected result:\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major copy assignment
   //=====================================================================================

   {
      test_ = "Column-major ZeroMatrix copy assignment";

      blaze::ZeroMatrix<int,blaze::columnMajor> Z1( 4UL, 3UL );
      blaze::ZeroMatrix<int,blaze::columnMajor> Z2;
      Z2 = Z1;

      checkRows    ( Z2, 4UL );
      checkColumns ( Z2, 3UL );
      checkNonZeros( Z2, 0UL );
      checkNonZeros( Z2, 0UL, 0UL );
      checkNonZeros( Z2, 1UL, 0UL );
      checkNonZeros( Z2, 2UL, 0UL );

      if( Z2(0,0) != 0 || Z2(0,1) != 0 || Z2(0,2) != 0 ||
          Z2(1,0) != 0 || Z2(1,1) != 0 || Z2(1,2) != 0 ||
          Z2(2,0) != 0 || Z2(2,1) != 0 || Z2(2,2) != 0 ||
          Z2(3,0) != 0 || Z2(3,1) != 0 || Z2(3,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << Z2 << "\n"
             << "   Expected result:\n( 0 0 0 )\n( 0 0 0 )\n( 0 0 0 )\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major ZeroMatrix copy assignment stress test";

      blaze::ZeroMatrix<int,blaze::columnMajor> Z1;

      for( size_t i=0UL; i<100UL; ++i )
      {
         const size_t m( blaze::rand<size_t>( 0UL, 10UL ) );
         const size_t n( blaze::rand<size_t>( 0UL, 10UL ) );
         const blaze::ZeroMatrix<int,blaze::columnMajor> Z2( m, n );

         Z1 = Z2;

         if( Z1 != Z2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << Z1 << "\n"
                << "   Expected result:\n" << Z2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major move assignment
   //=====================================================================================

   {
      test_ = "Column-major ZeroMatrix move assignment";

      blaze::ZeroMatrix<int,blaze::columnMajor> Z1( 4UL, 3UL );
      blaze::ZeroMatrix<int,blaze::columnMajor> Z2;

      Z2 = std::move( Z1 );

      checkRows    ( Z2, 4UL );
      checkColumns ( Z2, 3UL );
      checkNonZeros( Z2, 0UL );
      checkNonZeros( Z2, 0UL, 0UL );
      checkNonZeros( Z2, 1UL, 0UL );
      checkNonZeros( Z2, 2UL, 0UL );

      if( Z2(0,0) != 0 || Z2(0,1) != 0 || Z2(0,2) != 0 ||
          Z2(1,0) != 0 || Z2(1,1) != 0 || Z2(1,2) != 0 ||
          Z2(2,0) != 0 || Z2(2,1) != 0 || Z2(2,2) != 0 ||
          Z2(3,0) != 0 || Z2(3,1) != 0 || Z2(3,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << Z2 << "\n"
             << "   Expected result:\n( 0 0 0 )\n( 0 0 0 )\n( 0 0 0 )\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the ZeroMatrix function call operator.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of adding and accessing elements via the function call operator
// of the ZeroMatrix class template. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
void ClassTest::testFunctionCall()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major ZeroMatrix::operator()";

      blaze::ZeroMatrix<int,blaze::rowMajor> Z( 3UL, 4UL );

      checkRows    ( Z, 3UL );
      checkColumns ( Z, 4UL );
      checkNonZeros( Z, 0UL );
      checkNonZeros( Z, 0UL, 0UL );
      checkNonZeros( Z, 1UL, 0UL );
      checkNonZeros( Z, 2UL, 0UL );

      if( Z(0,0) != 0 || Z(0,1) != 0 || Z(0,2) != 0 || Z(0,3) != 0 ||
          Z(1,0) != 0 || Z(1,1) != 0 || Z(1,2) != 0 || Z(1,3) != 0 ||
          Z(2,0) != 0 || Z(2,1) != 0 || Z(2,2) != 0 || Z(2,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << Z << "\n"
             << "   Expected result:\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major ZeroMatrix::operator()";

      blaze::ZeroMatrix<int,blaze::columnMajor> Z( 4UL, 3UL );

      checkRows    ( Z, 4UL );
      checkColumns ( Z, 3UL );
      checkNonZeros( Z, 0UL );
      checkNonZeros( Z, 0UL, 0UL );
      checkNonZeros( Z, 1UL, 0UL );
      checkNonZeros( Z, 2UL, 0UL );

      if( Z(0,0) != 0 || Z(0,1) != 0 || Z(0,2) != 0 ||
          Z(1,0) != 0 || Z(1,1) != 0 || Z(1,2) != 0 ||
          Z(2,0) != 0 || Z(2,1) != 0 || Z(2,2) != 0 ||
          Z(3,0) != 0 || Z(3,1) != 0 || Z(3,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << Z << "\n"
             << "   Expected result:\n( 0 0 0 )\n( 0 0 0 )\n( 0 0 0 )\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c at() member function of the ZeroMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of adding and accessing elements via the \c at() member function
// of the ZeroMatrix class template. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
void ClassTest::testAt()
{
   //=====================================================================================
   // Row-major matrix tests
   //==========≈===========================================================================

   {
      test_ = "Row-major ZeroMatrix::at()";

      blaze::ZeroMatrix<int,blaze::rowMajor> Z( 3UL, 4UL );

      checkRows    ( Z, 3UL );
      checkColumns ( Z, 4UL );
      checkNonZeros( Z, 0UL );
      checkNonZeros( Z, 0UL, 0UL );
      checkNonZeros( Z, 1UL, 0UL );
      checkNonZeros( Z, 2UL, 0UL );

      if( Z.at(0,0) != 0 || Z.at(0,1) != 0 || Z.at(0,2) != 0 || Z.at(0,3) != 0 ||
          Z.at(1,0) != 0 || Z.at(1,1) != 0 || Z.at(1,2) != 0 || Z.at(1,3) != 0 ||
          Z.at(2,0) != 0 || Z.at(2,1) != 0 || Z.at(2,2) != 0 || Z.at(2,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Access via at() function failed\n"
             << " Details:\n"
             << "   Result:\n" << Z << "\n"
             << "   Expected result:\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major ZeroMatrix::at()";

      blaze::ZeroMatrix<int,blaze::columnMajor> Z( 4UL, 3UL );

      checkRows    ( Z, 4UL );
      checkColumns ( Z, 3UL );
      checkNonZeros( Z, 0UL );
      checkNonZeros( Z, 0UL, 0UL );
      checkNonZeros( Z, 1UL, 0UL );
      checkNonZeros( Z, 2UL, 0UL );

      if( Z.at(0,0) != 0 || Z.at(0,1) != 0 || Z.at(0,2) != 0 ||
          Z.at(1,0) != 0 || Z.at(1,1) != 0 || Z.at(1,2) != 0 ||
          Z.at(2,0) != 0 || Z.at(2,1) != 0 || Z.at(2,2) != 0 ||
          Z.at(3,0) != 0 || Z.at(3,1) != 0 || Z.at(3,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Access via at() function failed\n"
             << " Details:\n"
             << "   Result:\n" << Z << "\n"
             << "   Expected result:\n( 0 0 0 )\n( 0 0 0 )\n( 0 0 0 )\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the ZeroMatrix iterator implementation.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the iterator implementation of the ZeroMatrix class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testIterator()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      using MatrixType    = blaze::ZeroMatrix<int,blaze::rowMajor>;
      using ConstIterator = MatrixType::ConstIterator;

      MatrixType Z( 3UL, 4UL );

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

         const ptrdiff_t number( cend( Z, 1UL ) - cbegin( Z, 1UL ) );

         if( number != 0L ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid number of elements detected\n"
                << " Details:\n"
                << "   Number of elements         : " << number << "\n"
                << "   Expected number of elements: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing ConstIterator comparison
      {
         test_ = "Row-major ConstIterator comparison";

         ConstIterator it ( cbegin( Z, 1UL ) );
         ConstIterator end( cend( Z, 1UL ) );

         if( it != end ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator comparison failed\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      using MatrixType    = blaze::ZeroMatrix<int,blaze::columnMajor>;
      using ConstIterator = MatrixType::ConstIterator;

      MatrixType Z( 4UL, 3UL );

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

         const ptrdiff_t number( cend( Z, 1UL ) - cbegin( Z, 1UL ) );

         if( number != 0L ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid number of elements detected\n"
                << " Details:\n"
                << "   Number of elements         : " << number << "\n"
                << "   Expected number of elements: 0\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing ConstIterator comparison
      {
         test_ = "Column-major ConstIterator comparison";

         ConstIterator it ( cbegin( Z, 1UL ) );
         ConstIterator end( cend( Z, 1UL ) );

         if( it != end ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator comparison failed\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c nonZeros() member function of the ZeroMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c nonZeros() member function of the ZeroMatrix class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testNonZeros()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major ZeroMatrix::nonZeros()";

      blaze::ZeroMatrix<int,blaze::rowMajor> Z( 6UL, 8UL );

      checkRows    ( Z, 6UL );
      checkColumns ( Z, 8UL );
      checkNonZeros( Z, 0UL );
      checkNonZeros( Z, 0UL, 0UL );
      checkNonZeros( Z, 1UL, 0UL );
      checkNonZeros( Z, 2UL, 0UL );
      checkNonZeros( Z, 3UL, 0UL );
      checkNonZeros( Z, 4UL, 0UL );
      checkNonZeros( Z, 5UL, 0UL );
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major ZeroMatrix::nonZeros()";

      blaze::ZeroMatrix<int,blaze::columnMajor> Z( 8UL, 6UL );

      checkRows    ( Z, 8UL );
      checkColumns ( Z, 6UL );
      checkNonZeros( Z, 0UL );
      checkNonZeros( Z, 0UL, 0UL );
      checkNonZeros( Z, 1UL, 0UL );
      checkNonZeros( Z, 2UL, 0UL );
      checkNonZeros( Z, 3UL, 0UL );
      checkNonZeros( Z, 4UL, 0UL );
      checkNonZeros( Z, 5UL, 0UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c reset() member function of the ZeroMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c reset() member function of the ZeroMatrix class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testReset()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major ZeroMatrix::reset()";

      // Resetting a default constructed matrix
      {
         blaze::ZeroMatrix<int,blaze::rowMajor> Z;

         reset( Z );

         checkRows    ( Z, 0UL );
         checkColumns ( Z, 0UL );
         checkNonZeros( Z, 0UL );
      }

      // Resetting an initialized matrix
      {
         // Initialization check
         blaze::ZeroMatrix<int,blaze::rowMajor> Z( 3UL, 4UL );

         checkRows    ( Z, 3UL );
         checkColumns ( Z, 4UL );
         checkNonZeros( Z, 0UL );
         checkNonZeros( Z, 0UL, 0UL );
         checkNonZeros( Z, 1UL, 0UL );
         checkNonZeros( Z, 2UL, 0UL );

         if( Z(0,0) != 0 || Z(0,1) != 0 || Z(0,2) != 0 || Z(0,3) != 0 ||
             Z(1,0) != 0 || Z(1,1) != 0 || Z(1,2) != 0 || Z(1,3) != 0 ||
             Z(2,0) != 0 || Z(2,1) != 0 || Z(2,2) != 0 || Z(2,3) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << Z << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Resetting row 1
         reset( Z, 1UL );

         checkRows    ( Z, 3UL );
         checkColumns ( Z, 4UL );
         checkNonZeros( Z, 0UL );
         checkNonZeros( Z, 0UL, 0UL );
         checkNonZeros( Z, 1UL, 0UL );
         checkNonZeros( Z, 2UL, 0UL );

         if( Z(0,0) != 0 || Z(0,1) != 0 || Z(0,2) != 0 || Z(0,3) != 0 ||
             Z(1,0) != 0 || Z(1,1) != 0 || Z(1,2) != 0 || Z(1,3) != 0 ||
             Z(2,0) != 0 || Z(2,1) != 0 || Z(2,2) != 0 || Z(2,3) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Reset operation failed\n"
                << " Details:\n"
                << "   Result:\n" << Z << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Resetting the entire matrix
         reset( Z );

         checkRows    ( Z, 3UL );
         checkColumns ( Z, 4UL );
         checkNonZeros( Z, 0UL );
         checkNonZeros( Z, 0UL, 0UL );
         checkNonZeros( Z, 1UL, 0UL );
         checkNonZeros( Z, 2UL, 0UL );

         if( Z(0,0) != 0 || Z(0,1) != 0 || Z(0,2) != 0 || Z(0,3) != 0 ||
             Z(1,0) != 0 || Z(1,1) != 0 || Z(1,2) != 0 || Z(1,3) != 0 ||
             Z(2,0) != 0 || Z(2,1) != 0 || Z(2,2) != 0 || Z(2,3) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Reset operation failed\n"
                << " Details:\n"
                << "   Result:\n" << Z << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major ZeroMatrix::reset()";

      // Resetting a default constructed matrix
      {
         blaze::ZeroMatrix<int,blaze::columnMajor> Z;

         reset( Z );

         checkRows    ( Z, 0UL );
         checkColumns ( Z, 0UL );
         checkNonZeros( Z, 0UL );
      }

      // Resetting an initialized matrix
      {
         // Initialization check
         blaze::ZeroMatrix<int,blaze::columnMajor> Z( 4UL, 3UL );

         checkRows    ( Z, 4UL );
         checkColumns ( Z, 3UL );
         checkNonZeros( Z, 0UL );
         checkNonZeros( Z, 0UL, 0UL );
         checkNonZeros( Z, 1UL, 0UL );
         checkNonZeros( Z, 2UL, 0UL );

         if( Z(0,0) != 0 || Z(0,1) != 0 || Z(0,2) != 0 ||
             Z(1,0) != 0 || Z(1,1) != 0 || Z(1,2) != 0 ||
             Z(2,0) != 0 || Z(2,1) != 0 || Z(2,2) != 0 ||
             Z(3,0) != 0 || Z(3,1) != 0 || Z(3,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << Z << "\n"
                << "   Expected result:\n( 0 0 0 )\n( 0 0 0 )\n( 0 0 0 )\n( 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Resetting column 1
         reset( Z, 1UL );

         checkRows    ( Z, 4UL );
         checkColumns ( Z, 3UL );
         checkNonZeros( Z, 0UL );
         checkNonZeros( Z, 0UL, 0UL );
         checkNonZeros( Z, 1UL, 0UL );
         checkNonZeros( Z, 2UL, 0UL );

         if( Z(0,0) != 0 || Z(0,1) != 0 || Z(0,2) != 0 ||
             Z(1,0) != 0 || Z(1,1) != 0 || Z(1,2) != 0 ||
             Z(2,0) != 0 || Z(2,1) != 0 || Z(2,2) != 0 ||
             Z(3,0) != 0 || Z(3,1) != 0 || Z(3,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Reset operation failed\n"
                << " Details:\n"
                << "   Result:\n" << Z << "\n"
                << "   Expected result:\n( 0 0 0 )\n( 0 0 0 )\n( 0 0 0 )\n( 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Resetting the entire matrix
         reset( Z );

         checkRows    ( Z, 4UL );
         checkColumns ( Z, 3UL );
         checkNonZeros( Z, 0UL );
         checkNonZeros( Z, 0UL, 0UL );
         checkNonZeros( Z, 1UL, 0UL );
         checkNonZeros( Z, 2UL, 0UL );

         if( Z(0,0) != 0 || Z(0,1) != 0 || Z(0,2) != 0 ||
             Z(1,0) != 0 || Z(1,1) != 0 || Z(1,2) != 0 ||
             Z(2,0) != 0 || Z(2,1) != 0 || Z(2,2) != 0 ||
             Z(3,0) != 0 || Z(3,1) != 0 || Z(3,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Reset operation failed\n"
                << " Details:\n"
                << "   Result:\n" << Z << "\n"
                << "   Expected result:\n( 0 0 0 )\n( 0 0 0 )\n( 0 0 0 )\n( 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c clear() member function of the ZeroMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c clear() member function of the ZeroMatrix class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testClear()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major ZeroMatrix::clear()";

      // Clearing a default constructed matrix
      {
         blaze::ZeroMatrix<int,blaze::rowMajor> Z;

         clear( Z );

         checkRows    ( Z, 0UL );
         checkColumns ( Z, 0UL );
         checkNonZeros( Z, 0UL );
      }

      // Clearing an initialized matrix
      {
         // Initialization check
         blaze::ZeroMatrix<int,blaze::rowMajor> Z( 3UL, 4UL );

         checkRows    ( Z, 3UL );
         checkColumns ( Z, 4UL );
         checkNonZeros( Z, 0UL );
         checkNonZeros( Z, 0UL, 0UL );
         checkNonZeros( Z, 1UL, 0UL );
         checkNonZeros( Z, 2UL, 0UL );

         if( Z(0,0) != 0 || Z(0,1) != 0 || Z(0,2) != 0 || Z(0,3) != 0 ||
             Z(1,0) != 0 || Z(1,1) != 0 || Z(1,2) != 0 || Z(1,3) != 0 ||
             Z(2,0) != 0 || Z(2,1) != 0 || Z(2,2) != 0 || Z(2,3) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << Z << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Clearing the matrix
         clear( Z );

         checkRows    ( Z, 0UL );
         checkColumns ( Z, 0UL );
         checkNonZeros( Z, 0UL );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major ZeroMatrix::clear()";

      // Clearing a default constructed matrix
      {
         blaze::ZeroMatrix<int,blaze::columnMajor> Z;

         clear( Z );

         checkRows    ( Z, 0UL );
         checkColumns ( Z, 0UL );
         checkNonZeros( Z, 0UL );
      }

      // Clearing an initialized matrix
      {
         // Initialization check
         blaze::ZeroMatrix<int,blaze::columnMajor> Z( 4UL, 3UL );

         checkRows    ( Z, 4UL );
         checkColumns ( Z, 3UL );
         checkNonZeros( Z, 0UL );
         checkNonZeros( Z, 0UL, 0UL );
         checkNonZeros( Z, 1UL, 0UL );
         checkNonZeros( Z, 2UL, 0UL );

         if( Z(0,0) != 0 || Z(0,1) != 0 || Z(0,2) != 0 ||
             Z(1,0) != 0 || Z(1,1) != 0 || Z(1,2) != 0 ||
             Z(2,0) != 0 || Z(2,1) != 0 || Z(2,2) != 0 ||
             Z(3,0) != 0 || Z(3,1) != 0 || Z(3,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Clear operation failed\n"
                << " Details:\n"
                << "   Result:\n" << Z << "\n"
                << "   Expected result:\n( 0 0 0 )\n( 0 0 0 )\n( 0 0 0 )\n( 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Clearing the matrix
         clear( Z );

         checkRows    ( Z, 0UL );
         checkColumns ( Z, 0UL );
         checkNonZeros( Z, 0UL );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c resize() member function of the ZeroMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c resize() member function of the ZeroMatrix class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testResize()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major ZeroMatrix::resize()";

      // Initialization check
      blaze::ZeroMatrix<int,blaze::rowMajor> Z;

      checkRows    ( Z, 0UL );
      checkColumns ( Z, 0UL );
      checkNonZeros( Z, 0UL );

      // Resizing to 0x3
      Z.resize( 0UL, 3UL );

      checkRows    ( Z, 0UL );
      checkColumns ( Z, 3UL );
      checkNonZeros( Z, 0UL );

      // Resizing to 5x0
      Z.resize( 5UL, 0UL );

      checkRows    ( Z, 5UL );
      checkColumns ( Z, 0UL );
      checkNonZeros( Z, 0UL );

      // Resizing to 3x4
      Z.resize( 3UL, 4UL );

      checkRows    ( Z, 3UL );
      checkColumns ( Z, 4UL );
      checkNonZeros( Z, 0UL );
      checkNonZeros( Z, 0UL, 0UL );
      checkNonZeros( Z, 1UL, 0UL );
      checkNonZeros( Z, 2UL, 0UL );

      // Resizing to 2x1
      Z.resize( 2UL, 1UL );

      checkRows    ( Z, 2UL );
      checkColumns ( Z, 1UL );
      checkNonZeros( Z, 0UL );
      checkNonZeros( Z, 0UL, 0UL );
      checkNonZeros( Z, 1UL, 0UL );

      // Resizing to 3x2
      Z.resize( 3UL, 2UL );

      checkRows    ( Z, 3UL );
      checkColumns ( Z, 2UL );
      checkNonZeros( Z, 0UL );
      checkNonZeros( Z, 0UL, 0UL );
      checkNonZeros( Z, 1UL, 0UL );
      checkNonZeros( Z, 2UL, 0UL );

      // Resizing to 2x2
      Z.resize( 2UL, 2UL );

      checkRows    ( Z, 2UL );
      checkColumns ( Z, 2UL );
      checkNonZeros( Z, 0UL );
      checkNonZeros( Z, 0UL, 0UL );
      checkNonZeros( Z, 1UL, 0UL );

      // Resizing to 0x0
      Z.resize( 0UL, 0UL );

      checkRows    ( Z, 0UL );
      checkColumns ( Z, 0UL );
      checkNonZeros( Z, 0UL );
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major ZeroMatrix::resize()";

      // Initialization check
      blaze::ZeroMatrix<int,blaze::columnMajor> Z;

      checkRows    ( Z, 0UL );
      checkColumns ( Z, 0UL );
      checkNonZeros( Z, 0UL );

      // Resizing to 0x3
      Z.resize( 0UL, 3UL );

      checkRows    ( Z, 0UL );
      checkColumns ( Z, 3UL );
      checkNonZeros( Z, 0UL );

      // Resizing to 5x0
      Z.resize( 5UL, 0UL );

      checkRows    ( Z, 5UL );
      checkColumns ( Z, 0UL );
      checkNonZeros( Z, 0UL );

      // Resizing to 4x3
      Z.resize( 4UL, 3UL );

      checkRows    ( Z, 4UL );
      checkColumns ( Z, 3UL );
      checkNonZeros( Z, 0UL );
      checkNonZeros( Z, 0UL, 0UL );
      checkNonZeros( Z, 1UL, 0UL );
      checkNonZeros( Z, 2UL, 0UL );

      // Resizing to 1x2
      Z.resize( 1UL, 2UL );

      checkRows    ( Z, 1UL );
      checkColumns ( Z, 2UL );
      checkNonZeros( Z, 0UL );
      checkNonZeros( Z, 0UL, 0UL );
      checkNonZeros( Z, 1UL, 0UL );

      // Resizing to 2x3
      Z.resize( 2UL, 3UL );

      checkRows    ( Z, 2UL );
      checkColumns ( Z, 3UL );
      checkNonZeros( Z, 0UL );
      checkNonZeros( Z, 0UL, 0UL );
      checkNonZeros( Z, 1UL, 0UL );
      checkNonZeros( Z, 2UL, 0UL );

      // Resizing to 2x2
      Z.resize( 2UL, 2UL );

      checkRows    ( Z, 2UL );
      checkColumns ( Z, 2UL );
      checkNonZeros( Z, 0UL );
      checkNonZeros( Z, 0UL, 0UL );
      checkNonZeros( Z, 1UL, 0UL );

      // Resizing to 0x0
      Z.resize( 0UL, 0UL );

      checkRows    ( Z, 0UL );
      checkColumns ( Z, 0UL );
      checkNonZeros( Z, 0UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c swap() functionality of the ZeroMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c swap() function of the ZeroMatrix class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testSwap()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major ZeroMatrix swap";

      blaze::ZeroMatrix<int,blaze::rowMajor> Z1( 2UL, 3UL );
      blaze::ZeroMatrix<int,blaze::rowMajor> Z2( 3UL, 2UL );

      swap( Z1, Z2 );

      checkRows    ( Z1, 3UL );
      checkColumns ( Z1, 2UL );
      checkNonZeros( Z1, 0UL );
      checkNonZeros( Z1, 0UL, 0UL );
      checkNonZeros( Z1, 1UL, 0UL );
      checkNonZeros( Z1, 2UL, 0UL );

      if( Z1(0,0) != 0 || Z1(0,1) != 0 ||
          Z1(1,0) != 0 || Z1(1,1) != 0 ||
          Z1(2,0) != 0 || Z1(2,1) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the first matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << Z1 << "\n"
             << "   Expected result:\n( 0 0 )\n( 0 0 )\n( 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      checkRows    ( Z2, 2UL );
      checkColumns ( Z2, 3UL );
      checkNonZeros( Z2, 0UL );
      checkNonZeros( Z2, 0UL, 0UL );
      checkNonZeros( Z2, 1UL, 0UL );

      if( Z2(0,0) != 0 || Z2(0,1) != 0 || Z2(0,2) != 0 ||
          Z2(1,0) != 0 || Z2(1,1) != 0 || Z2(1,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the second matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << Z2 << "\n"
             << "   Expected result:\n( 0 0 0 )\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major ZeroMatrix swap";

      blaze::ZeroMatrix<int,blaze::columnMajor> Z1( 2UL, 3UL );
      blaze::ZeroMatrix<int,blaze::columnMajor> Z2( 3UL, 2UL );

      swap( Z1, Z2 );

      checkRows    ( Z1, 3UL );
      checkColumns ( Z1, 2UL );
      checkNonZeros( Z1, 0UL );
      checkNonZeros( Z1, 0UL, 0UL );
      checkNonZeros( Z1, 1UL, 0UL );

      if( Z1(0,0) != 0 || Z1(0,1) != 0 ||
          Z1(1,0) != 0 || Z1(1,1) != 0 ||
          Z1(2,0) != 0 || Z1(2,1) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the first matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << Z1 << "\n"
             << "   Expected result:\n( 0 0 )\n( 0 0 )\n( 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      checkRows    ( Z2, 2UL );
      checkColumns ( Z2, 3UL );
      checkNonZeros( Z2, 0UL );
      checkNonZeros( Z2, 0UL, 0UL );
      checkNonZeros( Z2, 1UL, 0UL );
      checkNonZeros( Z2, 2UL, 0UL );

      if( Z2(0,0) != 0 || Z2(0,1) != 0 || Z2(0,2) != 0 ||
          Z2(1,0) != 0 || Z2(1,1) != 0 || Z2(1,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the second matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << Z2 << "\n"
             << "   Expected result:\n( 0 0 0 )\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c find() member function of the ZeroMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c find() member function of the ZeroMatrix class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testFind()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major ZeroMatrix::find()";

      using ConstIterator = blaze::ZeroMatrix<int,blaze::rowMajor>::ConstIterator;

      // Initialization check
      blaze::ZeroMatrix<int,blaze::rowMajor> Z( 6UL, 8UL );

      checkRows    ( Z, 6UL );
      checkColumns ( Z, 8UL );
      checkNonZeros( Z, 0UL );
      checkNonZeros( Z, 0UL, 0UL );
      checkNonZeros( Z, 1UL, 0UL );
      checkNonZeros( Z, 2UL, 0UL );
      checkNonZeros( Z, 3UL, 0UL );
      checkNonZeros( Z, 4UL, 0UL );
      checkNonZeros( Z, 5UL, 0UL );

      // Searching for the first non-existing element
      {
         ConstIterator pos( Z.find( 0UL, 0UL ) );

         if( pos != Z.end( 0UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Non-existing element could be found\n"
                << " Details:\n"
                << "   Required index = 0\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 0\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << Z << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for the second non-existing element
      {
         ConstIterator pos( Z.find( 2UL, 4UL ) );

         if( pos != Z.end( 2UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Non-existing element could be found\n"
                << " Details:\n"
                << "   Required index = 4\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 0\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << Z << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for the third non-existing element
      {
         ConstIterator pos( Z.find( 5UL, 7UL ) );

         if( pos != Z.end( 5UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Non-existing element could be found\n"
                << " Details:\n"
                << "   Required index = 7\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 0\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << Z << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major ZeroMatrix::find()";

      using ConstIterator = blaze::ZeroMatrix<int,blaze::columnMajor>::ConstIterator;

      // Initialization check
      blaze::ZeroMatrix<int,blaze::columnMajor> Z( 8UL, 6UL );

      checkRows    ( Z, 8UL );
      checkColumns ( Z, 6UL );
      checkNonZeros( Z, 0UL );
      checkNonZeros( Z, 0UL, 0UL );
      checkNonZeros( Z, 1UL, 0UL );
      checkNonZeros( Z, 2UL, 0UL );
      checkNonZeros( Z, 3UL, 0UL );
      checkNonZeros( Z, 4UL, 0UL );
      checkNonZeros( Z, 5UL, 0UL );

      // Searching for the first non-existing element
      {
         ConstIterator pos( Z.find( 0UL, 0UL ) );

         if( pos != Z.end( 0UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Non-existing element could be found\n"
                << " Details:\n"
                << "   Required index = 0\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 0\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << Z << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for the second non-existing element
      {
         ConstIterator pos( Z.find( 4UL, 2UL ) );

         if( pos != Z.end( 2UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Non-existing element could be found\n"
                << " Details:\n"
                << "   Required index = 4\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 0\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << Z << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for the third non-existing element
      {
         ConstIterator pos( Z.find( 7UL, 5UL ) );

         if( pos != Z.end( 5UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Non-existing element could be found\n"
                << " Details:\n"
                << "   Required index = 7\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 0\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << Z << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c lowerBound() member function of the ZeroMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c lowerBound() member function of the ZeroMatrix class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testLowerBound()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major ZeroMatrix::lowerBound()";

      using ConstIterator = blaze::ZeroMatrix<int,blaze::rowMajor>::ConstIterator;

      // Initialization check
      blaze::ZeroMatrix<int,blaze::rowMajor> Z( 3UL, 4UL );

      checkRows    ( Z, 3UL );
      checkColumns ( Z, 4UL );
      checkNonZeros( Z, 0UL );
      checkNonZeros( Z, 0UL, 0UL );
      checkNonZeros( Z, 1UL, 0UL );
      checkNonZeros( Z, 2UL, 0UL );

      // Determining the lower bound for position (1,0)
      {
         ConstIterator pos( Z.lowerBound( 1UL, 0UL ) );

         if( pos != Z.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,0)\n"
                << "   Current matrix:\n" << Z << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for position (1,1)
      {
         ConstIterator pos( Z.lowerBound( 1UL, 1UL ) );

         if( pos != Z.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,1)\n"
                << "   Current matrix:\n" << Z << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for position (1,2)
      {
         ConstIterator pos( Z.lowerBound( 1UL, 2UL ) );

         if( pos != Z.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,2)\n"
                << "   Current matrix:\n" << Z << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major ZeroMatrix::lowerBound()";

      using ConstIterator = blaze::ZeroMatrix<int,blaze::columnMajor>::ConstIterator;

      // Initialization check
      blaze::ZeroMatrix<int,blaze::columnMajor> Z( 4UL, 3UL );

      checkRows    ( Z, 4UL );
      checkColumns ( Z, 3UL );
      checkNonZeros( Z, 0UL );
      checkNonZeros( Z, 0UL, 0UL );
      checkNonZeros( Z, 1UL, 0UL );
      checkNonZeros( Z, 2UL, 0UL );

      // Determining the lower bound for position (0,1)
      {
         ConstIterator pos( Z.lowerBound( 0UL, 1UL ) );

         if( pos != Z.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,0)\n"
                << "   Current matrix:\n" << Z << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for position (1,1)
      {
         ConstIterator pos( Z.lowerBound( 1UL, 1UL ) );

         if( pos != Z.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,1)\n"
                << "   Current matrix:\n" << Z << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for position (2,1)
      {
         ConstIterator pos( Z.lowerBound( 2UL, 1UL ) );

         if( pos != Z.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (2,1)\n"
                << "   Current matrix:\n" << Z << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c upperBound() member function of the ZeroMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c upperBound() member function of the ZeroMatrix class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testUpperBound()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major ZeroMatrix::upperBound()";

      using ConstIterator = blaze::ZeroMatrix<int,blaze::rowMajor>::ConstIterator;

      // Initialization check
      blaze::ZeroMatrix<int,blaze::rowMajor> Z( 3UL, 4UL );

      checkRows    ( Z, 3UL );
      checkColumns ( Z, 4UL );
      checkNonZeros( Z, 0UL );
      checkNonZeros( Z, 0UL, 0UL );
      checkNonZeros( Z, 1UL, 0UL );
      checkNonZeros( Z, 2UL, 0UL );

      // Determining the upper bound for position (1,0)
      {
         ConstIterator pos( Z.upperBound( 1UL, 0UL ) );

         if( pos != Z.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,0)\n"
                << "   Current matrix:\n" << Z << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for position (1,1)
      {
         ConstIterator pos( Z.upperBound( 1UL, 1UL ) );

         if( pos != Z.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,1)\n"
                << "   Current matrix:\n" << Z << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for position (1,2)
      {
         ConstIterator pos( Z.upperBound( 1UL, 2UL ) );

         if( pos != Z.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,2)\n"
                << "   Current matrix:\n" << Z << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major ZeroMatrix::upperBound()";

      using ConstIterator = blaze::ZeroMatrix<int,blaze::columnMajor>::ConstIterator;

      // Initialization check
      blaze::ZeroMatrix<int,blaze::columnMajor> Z( 4UL, 3UL );

      checkRows    ( Z, 4UL );
      checkColumns ( Z, 3UL );
      checkNonZeros( Z, 0UL );
      checkNonZeros( Z, 0UL, 0UL );
      checkNonZeros( Z, 1UL, 0UL );
      checkNonZeros( Z, 2UL, 0UL );

      // Determining the upper bound for position (0,1)
      {
         ConstIterator pos( Z.upperBound( 0UL, 1UL ) );

         if( pos != Z.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (0,1)\n"
                << "   Current matrix:\n" << Z << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for position (1,1)
      {
         ConstIterator pos( Z.upperBound( 1UL, 1UL ) );

         if( pos != Z.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,1)\n"
                << "   Current matrix:\n" << Z << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for position (2,1)
      {
         ConstIterator pos( Z.upperBound( 2UL, 1UL ) );

         if( pos != Z.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (2,1)\n"
                << "   Current matrix:\n" << Z << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c transpose() member function of the ZeroMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c transpose() member function of the ZeroMatrix
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

      blaze::ZeroMatrix<int,blaze::rowMajor> Z( 3UL, 4UL );

      transpose( Z );

      checkRows    ( Z, 4UL );
      checkColumns ( Z, 3UL );
      checkNonZeros( Z, 0UL );
      checkNonZeros( Z, 0UL, 0UL );
      checkNonZeros( Z, 1UL, 0UL );
      checkNonZeros( Z, 2UL, 0UL );
      checkNonZeros( Z, 3UL, 0UL );

      if( Z(0,0) != 0 || Z(0,1) != 0 || Z(0,2) != 0 ||
          Z(1,0) != 0 || Z(1,1) != 0 || Z(1,2) != 0 ||
          Z(2,0) != 0 || Z(2,1) != 0 || Z(2,2) != 0 ||
          Z(3,0) != 0 || Z(3,1) != 0 || Z(3,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Transposition failed\n"
             << " Details:\n"
             << "   Result:\n" << Z << "\n"
             << "   Expected result:\n( 0 0 0 )\n( 0 0 0 )\n( 0 0 0 )\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major self-transpose via trans()";

      blaze::ZeroMatrix<int,blaze::rowMajor> Z( 3UL, 4UL );

      Z = trans( Z );

      checkRows    ( Z, 4UL );
      checkColumns ( Z, 3UL );
      checkNonZeros( Z, 0UL );
      checkNonZeros( Z, 0UL, 0UL );
      checkNonZeros( Z, 1UL, 0UL );
      checkNonZeros( Z, 2UL, 0UL );
      checkNonZeros( Z, 3UL, 0UL );

      if( Z(0,0) != 0 || Z(0,1) != 0 || Z(0,2) != 0 ||
          Z(1,0) != 0 || Z(1,1) != 0 || Z(1,2) != 0 ||
          Z(2,0) != 0 || Z(2,1) != 0 || Z(2,2) != 0 ||
          Z(3,0) != 0 || Z(3,1) != 0 || Z(3,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Transposition failed\n"
             << " Details:\n"
             << "   Result:\n" << Z << "\n"
             << "   Expected result:\n( 0 0 0 )\n( 0 0 0 )\n( 0 0 0 )\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major self-transpose via transpose()";

      blaze::ZeroMatrix<int,blaze::columnMajor> Z( 4UL, 3UL );

      transpose( Z );

      checkRows    ( Z, 3UL );
      checkColumns ( Z, 4UL );
      checkNonZeros( Z, 0UL );
      checkNonZeros( Z, 0UL, 0UL );
      checkNonZeros( Z, 1UL, 0UL );
      checkNonZeros( Z, 2UL, 0UL );
      checkNonZeros( Z, 3UL, 0UL );

      if( Z(0,0) != 0 || Z(0,1) != 0 || Z(0,2) != 0 || Z(0,3) != 0 ||
          Z(1,0) != 0 || Z(1,1) != 0 || Z(1,2) != 0 || Z(1,3) != 0 ||
          Z(2,0) != 0 || Z(2,1) != 0 || Z(2,2) != 0 || Z(2,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Transposition failed\n"
             << " Details:\n"
             << "   Result:\n" << Z << "\n"
             << "   Expected result:\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major self-transpose via trans()";

      blaze::ZeroMatrix<int,blaze::columnMajor> Z( 4UL, 3UL );

      Z = trans( Z );

      checkRows    ( Z, 3UL );
      checkColumns ( Z, 4UL );
      checkNonZeros( Z, 0UL );
      checkNonZeros( Z, 0UL, 0UL );
      checkNonZeros( Z, 1UL, 0UL );
      checkNonZeros( Z, 2UL, 0UL );
      checkNonZeros( Z, 3UL, 0UL );

      if( Z(0,0) != 0 || Z(0,1) != 0 || Z(0,2) != 0 || Z(0,3) != 0 ||
          Z(1,0) != 0 || Z(1,1) != 0 || Z(1,2) != 0 || Z(1,3) != 0 ||
          Z(2,0) != 0 || Z(2,1) != 0 || Z(2,2) != 0 || Z(2,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Transposition failed\n"
             << " Details:\n"
             << "   Result:\n" << Z << "\n"
             << "   Expected result:\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c ctranspose() member function of the ZeroMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c ctranspose() member function of the ZeroMatrix
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

      blaze::ZeroMatrix<cplx,blaze::rowMajor> Z( 3UL, 4UL );

      ctranspose( Z );

      checkRows    ( Z, 4UL );
      checkColumns ( Z, 3UL );
      checkNonZeros( Z, 0UL );
      checkNonZeros( Z, 0UL, 0UL );
      checkNonZeros( Z, 1UL, 0UL );
      checkNonZeros( Z, 2UL, 0UL );
      checkNonZeros( Z, 3UL, 0UL );

      if( Z(0,0) != cplx(0,0) || Z(0,1) != cplx(0,0) || Z(0,2) != cplx(0,0) ||
          Z(1,0) != cplx(0,0) || Z(1,1) != cplx(0,0) || Z(1,2) != cplx(0,0) ||
          Z(2,0) != cplx(0,0) || Z(2,1) != cplx(0,0) || Z(2,2) != cplx(0,0) ||
          Z(3,0) != cplx(0,0) || Z(3,1) != cplx(0,0) || Z(3,2) != cplx(0,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Transposition failed\n"
             << " Details:\n"
             << "   Result:\n" << Z << "\n"
             << "   Expected result:\n( (0,0) (0,0) (0,0) )\n"
                                     "( (0,0) (0,0) (0,0) )\n"
                                     "( (0,0) (0,0) (0,0) )\n"
                                     "( (0,0) (0,0) (0,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major self-transpose via ctranspose()";

      using cplx = blaze::complex<int>;

      blaze::ZeroMatrix<cplx,blaze::rowMajor> Z( 3UL, 4UL );

      Z = ctrans( Z );

      checkRows    ( Z, 4UL );
      checkColumns ( Z, 3UL );
      checkNonZeros( Z, 0UL );
      checkNonZeros( Z, 0UL, 0UL );
      checkNonZeros( Z, 1UL, 0UL );
      checkNonZeros( Z, 2UL, 0UL );
      checkNonZeros( Z, 3UL, 0UL );

      if( Z(0,0) != cplx(0,0) || Z(0,1) != cplx(0,0) || Z(0,2) != cplx(0,0) ||
          Z(1,0) != cplx(0,0) || Z(1,1) != cplx(0,0) || Z(1,2) != cplx(0,0) ||
          Z(2,0) != cplx(0,0) || Z(2,1) != cplx(0,0) || Z(2,2) != cplx(0,0) ||
          Z(3,0) != cplx(0,0) || Z(3,1) != cplx(0,0) || Z(3,2) != cplx(0,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Transposition failed\n"
             << " Details:\n"
             << "   Result:\n" << Z << "\n"
             << "   Expected result:\n( (1,0) (0,0) (0,0) )\n"
                                     "( (0,0) (1,0) (0,0) )\n"
                                     "( (0,0) (0,0) (1,0) )\n"
                                     "( (0,0) (0,0) (0,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major self-transpose via ctranspose()";

      using cplx = blaze::complex<int>;

      blaze::ZeroMatrix<cplx,blaze::columnMajor> Z( 4UL, 3UL );

      ctranspose( Z );

      checkRows    ( Z, 3UL );
      checkColumns ( Z, 4UL );
      checkNonZeros( Z, 0UL );
      checkNonZeros( Z, 0UL, 0UL );
      checkNonZeros( Z, 1UL, 0UL );
      checkNonZeros( Z, 2UL, 0UL );
      checkNonZeros( Z, 2UL, 0UL );

      if( Z(0,0) != cplx(0,0) || Z(0,1) != cplx(0,0) || Z(0,2) != cplx(0,0) || Z(0,3) != cplx(0,0) ||
          Z(1,0) != cplx(0,0) || Z(1,1) != cplx(0,0) || Z(1,2) != cplx(0,0) || Z(1,3) != cplx(0,0) ||
          Z(2,0) != cplx(0,0) || Z(2,1) != cplx(0,0) || Z(2,2) != cplx(0,0) || Z(2,3) != cplx(0,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Transposition failed\n"
             << " Details:\n"
             << "   Result:\n" << Z << "\n"
             << "   Expected result:\n( (0,0) (0,0) (0,0) (0,0) )\n"
                                     "( (0,0) (0,0) (0,0) (0,0) )\n"
                                     "( (0,0) (0,0) (0,0) (0,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major self-transpose via ctranspose()";

      using cplx = blaze::complex<int>;

      blaze::ZeroMatrix<cplx,blaze::columnMajor> Z( 4UL, 3UL );

      Z = ctrans( Z );

      checkRows    ( Z, 3UL );
      checkColumns ( Z, 4UL );
      checkNonZeros( Z, 0UL );
      checkNonZeros( Z, 0UL, 0UL );
      checkNonZeros( Z, 1UL, 0UL );
      checkNonZeros( Z, 2UL, 0UL );
      checkNonZeros( Z, 2UL, 0UL );

      if( Z(0,0) != cplx(0,0) || Z(0,1) != cplx(0,0) || Z(0,2) != cplx(0,0) || Z(0,3) != cplx(0,0) ||
          Z(1,0) != cplx(0,0) || Z(1,1) != cplx(0,0) || Z(1,2) != cplx(0,0) || Z(1,3) != cplx(0,0) ||
          Z(2,0) != cplx(0,0) || Z(2,1) != cplx(0,0) || Z(2,2) != cplx(0,0) || Z(2,3) != cplx(0,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Transposition failed\n"
             << " Details:\n"
             << "   Result:\n" << Z << "\n"
             << "   Expected result:\n( (0,0) (0,0) (0,0) (0,0) )\n"
                                     "( (0,0) (0,0) (0,0) (0,0) )\n"
                                     "( (0,0) (0,0) (0,0) (0,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c isDefault() function with the ZeroMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isDefault() function with the ZeroMatrix class
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
         blaze::ZeroMatrix<int,blaze::rowMajor> Z;

         if( isDefault( Z ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << Z << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with 0x4 matrix (non-default)
      {
         blaze::ZeroMatrix<int,blaze::rowMajor> Z( 0UL, 4UL );

         if( isDefault( Z ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << Z << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with 3x0 matrix (non-default)
      {
         blaze::ZeroMatrix<int,blaze::rowMajor> Z( 3UL, 0UL );

         if( isDefault( Z ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << Z << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with 3x4 matrix (non-default)
      {
         blaze::ZeroMatrix<int,blaze::rowMajor> Z( 3UL, 4UL );

         if( isDefault( Z(0,1) ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix element: " << Z(0,1) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( Z ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << Z << "\n";
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
         blaze::ZeroMatrix<int,blaze::columnMajor> Z;

         if( isDefault( Z ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << Z << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with 0x3 matrix (non-default)
      {
         blaze::ZeroMatrix<int,blaze::rowMajor> Z( 0UL, 3UL );

         if( isDefault( Z ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << Z << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with 4x0 matrix (non-default)
      {
         blaze::ZeroMatrix<int,blaze::rowMajor> Z( 4UL, 0UL );

         if( isDefault( Z ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << Z << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with 4x3 matrix (non-default)
      {
         blaze::ZeroMatrix<int,blaze::columnMajor> Z( 4UL, 3UL );

         if( isDefault( Z(1,0) ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix element: " << Z(1,0) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( Z ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << Z << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************

} // namespace zeromatrix

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
   std::cout << "   Running ZeroMatrix class test..." << std::endl;

   try
   {
      RUN_ZEROMATRIX_CLASS_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during ZeroMatrix class test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
