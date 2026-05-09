//=================================================================================================
/*!
//  \file src/mathtest/matrices/compressedmatrix/ProxyTest.cpp
//  \brief Source file for the CompressedMatrix proxy test
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
#include <blazetest/mathtest/matrices/compressedmatrix/ProxyTest.h>
#include <blazetest/system/LAPACK.h>

#ifdef BLAZE_USE_HPX_THREADS
#  include <hpx/hpx_main.hpp>
#endif


namespace blazetest {

namespace mathtest {

namespace matrices {

namespace compressedmatrix {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the CompressedMatrix proxy test.
//
// \exception std::runtime_error Operation error detected.
*/
ProxyTest::ProxyTest()
{
   testAssignment();
   testAddAssign();
   testSubAssign();
   testMultAssign();
   testDivAssign();
   testModAssign();
   testScaling();
   testSubscript();
   testFunctionCall();
   testIterator();
   testNonZeros();
   testReset();
   testClear();
   testResize();
   testExtend();
   testReserve();
   testTrim();
   testSwap();
   testSet();
   testInsert();
   testAppend();
   testErase();
   testFind();
   testLowerBound();
   testUpperBound();
   testTranspose();
   testCTranspose();
   testInvert();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the MatrixAccessProxy assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all assignment operators of the MatrixAccessProxy class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ProxyTest::testAssignment()
{
   //=====================================================================================
   // Row-major homogeneous assignment
   //=====================================================================================

   {
      test_ = "Row-major MatrixAccessProxy homogeneous assignment";

      DVM mat( 2UL, 2UL, 1UL );
      mat(0,1) = DV( 3UL, 2 );

      mat(0,1) = 4;

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 0UL );
      checkSize    ( mat(0,1), 3UL );
      checkCapacity( mat(0,1), 3UL );
      checkNonZeros( mat(0,1), 3UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 0UL );

      if( mat(0,1)[0] != 4 || mat(0,1)[1] != 4 || mat(0,1)[2] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat(0,1) << "\n"
             << "   Expected result:\n( 4 4 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major list assignment
   //=====================================================================================

   {
      test_ = "Row-major MatrixAccessProxy 1D initializer list assignment";

      DVM mat( 2UL, 2UL, 1UL );
      mat(0,1) = DV( 3UL, 2 );

      mat(0,1) = { 1, -2, 3 };

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 0UL );
      checkSize    ( mat(0,1), 3UL );
      checkCapacity( mat(0,1), 3UL );
      checkNonZeros( mat(0,1), 3UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 0UL );

      if( mat(0,1)[0] != 1 || mat(0,1)[1] != -2 || mat(0,1)[2] != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat(0,1) << "\n"
             << "   Expected result:\n( 1 -2 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major MatrixAccessProxy 2D initializer list assignment";

      using blaze::initializer_list;

      DMM mat( 2UL, 2UL, 1UL );
      mat(0,1) = DM( 3UL, 3UL, 2 );

      initializer_list< initializer_list<int> > list = { { 1, -2, 3 }, { -2, 4, -6 }, { 3, -6, 9 } };
      mat(0,1) = list;

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0), 0UL );
      checkColumns ( mat(0,0), 0UL );
      checkRows    ( mat(0,1), 3UL );
      checkColumns ( mat(0,1), 3UL );
      checkCapacity( mat(0,1), 9UL );
      checkNonZeros( mat(0,1), 9UL );
      checkRows    ( mat(1,0), 0UL );
      checkColumns ( mat(1,0), 0UL );
      checkRows    ( mat(1,1), 0UL );
      checkColumns ( mat(1,1), 0UL );

      if( mat(0,1)(0,0) !=  1 || mat(0,1)(0,1) != -2 || mat(0,1)(0,2) !=  3 ||
          mat(0,1)(1,0) != -2 || mat(0,1)(1,1) !=  4 || mat(0,1)(1,2) != -6 ||
          mat(0,1)(2,0) !=  3 || mat(0,1)(2,1) != -6 || mat(0,1)(2,2) !=  9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat(0,1) << "\n"
             << "   Expected result:\n(  1 -2  3 )\n( -2  4 -6 )\n(  3 -6  9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major array assignment
   //=====================================================================================

   {
      test_ = "Row-major MatrixAccessProxy array assignment";

      const int array[3] = { 1, 2, 3 };
      DVM mat( 2UL, 2UL, 1UL );
      mat(0,1) = DV( 1UL, 2 );

      mat(0,1) = array;

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 0UL );
      checkSize    ( mat(0,1), 3UL );
      checkCapacity( mat(0,1), 3UL );
      checkNonZeros( mat(0,1), 3UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 0UL );

      if( mat(0,1)[0] != 1 || mat(0,1)[1] != 2 || mat(0,1)[2] != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat(0,1) << "\n"
             << "   Expected result:\n( 1 2 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major copy assignment
   //=====================================================================================

   {
      test_ = "Row-major MatrixAccessProxy copy assignment";

      DVM mat( 2UL, 2UL, 1UL );
      mat(0,1) = DV( 3UL, 2 );

      mat(1,0) = mat(0,1);

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 2UL );
      checkNonZeros( mat, 2UL );

      checkSize    ( mat(0,0), 0UL );
      checkSize    ( mat(0,1), 3UL );
      checkCapacity( mat(0,1), 3UL );
      checkNonZeros( mat(0,1), 3UL );
      checkSize    ( mat(1,0), 3UL );
      checkCapacity( mat(1,0), 3UL );
      checkNonZeros( mat(1,0), 3UL );
      checkSize    ( mat(1,1), 0UL );

      if( mat(0,1)[0] != 2 || mat(0,1)[1] != 2 || mat(0,1)[2] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat(0,1) << "\n"
             << "   Expected result:\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major dense vector assignment
   //=====================================================================================

   {
      test_ = "Row-major MatrixAccessProxy dense vector assignment";

      DV tmp( 3UL );
      tmp[0] = 1;
      tmp[1] = 2;
      tmp[2] = 3;
      DVM mat( 2UL, 2UL, 1UL );

      mat(0,1) = tmp;

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 0UL );
      checkSize    ( mat(0,1), 3UL );
      checkCapacity( mat(0,1), 3UL );
      checkNonZeros( mat(0,1), 3UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 0UL );

      if( mat(0,1)[0] != 1 || mat(0,1)[1] != 2 || mat(0,1)[2] != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat(0,1) << "\n"
             << "   Expected result:\n( 1 2 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse vector assignment
   //=====================================================================================

   {
      test_ = "Row-major MatrixAccessProxy sparse vector assignment";

      SV tmp( 3UL );
      tmp[1] = 2;
      DVM mat( 2UL, 2UL, 1UL );

      mat(0,1) = tmp;

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 0UL );
      checkSize    ( mat(0,1), 3UL );
      checkCapacity( mat(0,1), 3UL );
      checkNonZeros( mat(0,1), 1UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 0UL );

      if( mat(0,1)[0] != 0 || mat(0,1)[1] != 2 || mat(0,1)[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat(0,1) << "\n"
             << "   Expected result:\n( 0 2 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major homogeneous assignment
   //=====================================================================================

   {
      test_ = "Column-major MatrixAccessProxy homogeneous assignment";

      ODVM mat( 2UL, 2UL, 1UL );
      mat(0,1) = DV( 3UL, 2 );

      mat(0,1) = 4;

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 0UL );
      checkSize    ( mat(0,1), 3UL );
      checkCapacity( mat(0,1), 3UL );
      checkNonZeros( mat(0,1), 3UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 0UL );

      if( mat(0,1)[0] != 4 || mat(0,1)[1] != 4 || mat(0,1)[2] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat(0,1) << "\n"
             << "   Expected result:\n( 4 4 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major list assignment
   //=====================================================================================

   {
      test_ = "Column-major MatrixAccessProxy 1D initializer list assignment";

      ODVM mat( 2UL, 2UL, 1UL );
      mat(0,1) = DV( 3UL, 2 );

      mat(0,1) = { 1, -2, 3 };

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 0UL );
      checkSize    ( mat(0,1), 3UL );
      checkCapacity( mat(0,1), 3UL );
      checkNonZeros( mat(0,1), 3UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 0UL );

      if( mat(0,1)[0] != 1 || mat(0,1)[1] != -2 || mat(0,1)[2] != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat(0,1) << "\n"
             << "   Expected result:\n( 1 -2 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major MatrixAccessProxy 2D initializer list assignment";

      using blaze::initializer_list;

      ODMM mat( 2UL, 2UL, 1UL );
      mat(0,1) = DM( 3UL, 3UL, 2 );

      initializer_list< initializer_list<int> > list = { { 1, -2, 3 }, { -2, 4, -6 }, { 3, -6, 9 } };
      mat(0,1) = list;

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0), 0UL );
      checkColumns ( mat(0,0), 0UL );
      checkRows    ( mat(0,1), 3UL );
      checkColumns ( mat(0,1), 3UL );
      checkCapacity( mat(0,1), 9UL );
      checkNonZeros( mat(0,1), 9UL );
      checkRows    ( mat(1,0), 0UL );
      checkColumns ( mat(1,0), 0UL );
      checkRows    ( mat(1,1), 0UL );
      checkColumns ( mat(1,1), 0UL );

      if( mat(0,1)(0,0) !=  1 || mat(0,1)(0,1) != -2 || mat(0,1)(0,2) !=  3 ||
          mat(0,1)(1,0) != -2 || mat(0,1)(1,1) !=  4 || mat(0,1)(1,2) != -6 ||
          mat(0,1)(2,0) !=  3 || mat(0,1)(2,1) != -6 || mat(0,1)(2,2) !=  9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat(0,1) << "\n"
             << "   Expected result:\n(  1 -2  3 )\n( -2  4 -6 )\n(  3 -6  9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major array assignment
   //=====================================================================================

   {
      test_ = "Column-major MatrixAccessProxy array assignment";

      const int array[3] = { 1, 2, 3 };
      ODVM mat( 2UL, 2UL, 1UL );
      mat(0,1) = DV( 1UL, 2 );

      mat(0,1) = array;

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 0UL );
      checkSize    ( mat(0,1), 3UL );
      checkCapacity( mat(0,1), 3UL );
      checkNonZeros( mat(0,1), 3UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 0UL );

      if( mat(0,1)[0] != 1 || mat(0,1)[1] != 2 || mat(0,1)[2] != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat(0,1) << "\n"
             << "   Expected result:\n( 1 2 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major copy assignment
   //=====================================================================================

   {
      test_ = "Column-major MatrixAccessProxy copy assignment";

      ODVM mat( 2UL, 2UL, 1UL );
      mat(0,1) = DV( 3UL, 2 );

      mat(1,0) = mat(0,1);

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 2UL );
      checkNonZeros( mat, 2UL );

      checkSize    ( mat(0,0), 0UL );
      checkSize    ( mat(0,1), 3UL );
      checkCapacity( mat(0,1), 3UL );
      checkNonZeros( mat(0,1), 3UL );
      checkSize    ( mat(1,0), 3UL );
      checkCapacity( mat(1,0), 3UL );
      checkNonZeros( mat(1,0), 3UL );
      checkSize    ( mat(1,1), 0UL );

      if( mat(0,1)[0] != 2 || mat(0,1)[1] != 2 || mat(0,1)[2] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat(0,1) << "\n"
             << "   Expected result:\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense vector assignment
   //=====================================================================================

   {
      test_ = "Column-major MatrixAccessProxy dense vector assignment";

      DV tmp( 3UL );
      tmp[0] = 1;
      tmp[1] = 2;
      tmp[2] = 3;
      ODVM mat( 2UL, 2UL, 1UL );

      mat(0,1) = tmp;

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 0UL );
      checkSize    ( mat(0,1), 3UL );
      checkCapacity( mat(0,1), 3UL );
      checkNonZeros( mat(0,1), 3UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 0UL );

      if( mat(0,1)[0] != 1 || mat(0,1)[1] != 2 || mat(0,1)[2] != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat(0,1) << "\n"
             << "   Expected result:\n( 1 2 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse vector assignment
   //=====================================================================================

   {
      test_ = "Column-major MatrixAccessProxy sparse vector assignment";

      SV tmp( 3UL );
      tmp[1] = 2;
      DVM mat( 2UL, 2UL, 1UL );

      mat(0,1) = tmp;

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 0UL );
      checkSize    ( mat(0,1), 3UL );
      checkCapacity( mat(0,1), 3UL );
      checkNonZeros( mat(0,1), 1UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 0UL );

      if( mat(0,1)[0] != 0 || mat(0,1)[1] != 2 || mat(0,1)[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat(0,1) << "\n"
             << "   Expected result:\n( 0 2 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the MatrixAccessProxy addition assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the addition assignment operators of the MatrixAccessProxy
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ProxyTest::testAddAssign()
{
   //=====================================================================================
   // Row-major dense vector addition assignment
   //=====================================================================================

   {
      test_ = "Row-major MatrixAccessProxy dense vector addition assignment";

      DV tmp( 3UL );
      tmp[0] = 1;
      tmp[1] = 2;
      tmp[2] = 3;
      DVM mat( 2UL, 2UL, 1UL );
      mat(0,1) = tmp;

      mat(0,1) += tmp;

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 0UL );
      checkSize    ( mat(0,1), 3UL );
      checkCapacity( mat(0,1), 3UL );
      checkNonZeros( mat(0,1), 3UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 0UL );

      if( mat(0,1)[0] != 2 || mat(0,1)[1] != 4 || mat(0,1)[2] != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat(0,1) << "\n"
             << "   Expected result:\n( 2 4 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse vector addition assignment
   //=====================================================================================

   {
      test_ = "Row-major MatrixAccessProxy sparse vector addition assignment";

      SV tmp( 3UL );
      tmp[1] = 2;
      DVM mat( 2UL, 2UL, 1UL );
      mat(0,1) = tmp;

      mat(0,1) += tmp;

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 0UL );
      checkSize    ( mat(0,1), 3UL );
      checkCapacity( mat(0,1), 3UL );
      checkNonZeros( mat(0,1), 1UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 0UL );

      if( mat(0,1)[0] != 0 || mat(0,1)[1] != 4 || mat(0,1)[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat(0,1) << "\n"
             << "   Expected result:\n( 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense vector addition assignment
   //=====================================================================================

   {
      test_ = "Column-major MatrixAccessProxy dense vector addition assignment";

      DV tmp( 3UL );
      tmp[0] = 1;
      tmp[1] = 2;
      tmp[2] = 3;
      DVM mat( 2UL, 2UL, 1UL );
      mat(0,1) = tmp;

      mat(0,1) += tmp;

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 0UL );
      checkSize    ( mat(0,1), 3UL );
      checkCapacity( mat(0,1), 3UL );
      checkNonZeros( mat(0,1), 3UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 0UL );

      if( mat(0,1)[0] != 2 || mat(0,1)[1] != 4 || mat(0,1)[2] != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat(0,1) << "\n"
             << "   Expected result:\n( 2 4 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse vector addition assignment
   //=====================================================================================

   {
      test_ = "Column-major MatrixAccessProxy sparse vector addition assignment";

      SV tmp( 3UL );
      tmp[1] = 2;
      DVM mat( 2UL, 2UL, 1UL );
      mat(0,1) = tmp;

      mat(0,1) += tmp;

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 0UL );
      checkSize    ( mat(0,1), 3UL );
      checkCapacity( mat(0,1), 3UL );
      checkNonZeros( mat(0,1), 1UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 0UL );

      if( mat(0,1)[0] != 0 || mat(0,1)[1] != 4 || mat(0,1)[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat(0,1) << "\n"
             << "   Expected result:\n( 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the MatrixAccessProxy subtraction assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the subtraction assignment operators of the MatrixAccessProxy
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ProxyTest::testSubAssign()
{
   //=====================================================================================
   // Row-major dense vector subtraction assignment
   //=====================================================================================

   {
      test_ = "Row-major MatrixAccessProxy dense vector subtraction assignment";

      DV tmp( 3UL );
      tmp[0] = 1;
      tmp[1] = 2;
      tmp[2] = 3;
      DVM mat( 2UL, 2UL, 1UL );
      mat(0,1) = tmp;

      mat(0,1) -= tmp;

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 0UL );
      checkSize    ( mat(0,1), 3UL );
      checkCapacity( mat(0,1), 3UL );
      checkNonZeros( mat(0,1), 0UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 0UL );

      if( mat(0,1)[0] != 0 || mat(0,1)[1] != 0 || mat(0,1)[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat(0,1) << "\n"
             << "   Expected result:\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse vector subtraction assignment
   //=====================================================================================

   {
      test_ = "Row-major MatrixAccessProxy sparse vector subtraction assignment";

      SV tmp( 3UL );
      tmp[1] = 2;
      DVM mat( 2UL, 2UL, 1UL );
      mat(0,1) = tmp;

      mat(0,1) -= tmp;

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 0UL );
      checkSize    ( mat(0,1), 3UL );
      checkCapacity( mat(0,1), 3UL );
      checkNonZeros( mat(0,1), 0UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 0UL );

      if( mat(0,1)[0] != 0 || mat(0,1)[1] != 0 || mat(0,1)[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat(0,1) << "\n"
             << "   Expected result:\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense vector subtraction assignment
   //=====================================================================================

   {
      test_ = "Column-major MatrixAccessProxy dense vector subtraction assignment";

      DV tmp( 3UL );
      tmp[0] = 1;
      tmp[1] = 2;
      tmp[2] = 3;
      DVM mat( 2UL, 2UL, 1UL );
      mat(0,1) = tmp;

      mat(0,1) -= tmp;

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 0UL );
      checkSize    ( mat(0,1), 3UL );
      checkCapacity( mat(0,1), 3UL );
      checkNonZeros( mat(0,1), 0UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 0UL );

      if( mat(0,1)[0] != 0 || mat(0,1)[1] != 0 || mat(0,1)[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat(0,1) << "\n"
             << "   Expected result:\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse vector subtraction assignment
   //=====================================================================================

   {
      test_ = "Column-major MatrixAccessProxy sparse vector subtraction assignment";

      SV tmp( 3UL );
      tmp[1] = 2;
      DVM mat( 2UL, 2UL, 1UL );
      mat(0,1) = tmp;

      mat(0,1) -= tmp;

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 0UL );
      checkSize    ( mat(0,1), 3UL );
      checkCapacity( mat(0,1), 3UL );
      checkNonZeros( mat(0,1), 0UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 0UL );

      if( mat(0,1)[0] != 0 || mat(0,1)[1] != 0 || mat(0,1)[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat(0,1) << "\n"
             << "   Expected result:\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the MatrixAccessProxy multiplication assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the multiplication assignment operators of the MatrixAccessProxy
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ProxyTest::testMultAssign()
{
   //=====================================================================================
   // Row-major dense vector multiplication assignment
   //=====================================================================================

   {
      test_ = "Row-major MatrixAccessProxy dense vector multiplication assignment";

      DV tmp( 3UL );
      tmp[0] = 1;
      tmp[1] = 2;
      tmp[2] = 3;
      DVM mat( 2UL, 2UL, 1UL );
      mat(0,1) = tmp;

      mat(0,1) *= tmp;

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 0UL );
      checkSize    ( mat(0,1), 3UL );
      checkCapacity( mat(0,1), 3UL );
      checkNonZeros( mat(0,1), 3UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 0UL );

      if( mat(0,1)[0] != 1 || mat(0,1)[1] != 4 || mat(0,1)[2] != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat(0,1) << "\n"
             << "   Expected result:\n( 1 4 9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse vector multiplication assignment
   //=====================================================================================

   {
      test_ = "Row-major MatrixAccessProxy sparse vector multiplication assignment";

      SV tmp( 3UL );
      tmp[1] = 2;
      DVM mat( 2UL, 2UL, 1UL );
      mat(0,1) = tmp;

      mat(0,1) *= tmp;

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 0UL );
      checkSize    ( mat(0,1), 3UL );
      checkCapacity( mat(0,1), 3UL );
      checkNonZeros( mat(0,1), 1UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 0UL );

      if( mat(0,1)[0] != 0 || mat(0,1)[1] != 4 || mat(0,1)[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat(0,1) << "\n"
             << "   Expected result:\n( 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense vector multiplication assignment
   //=====================================================================================

   {
      test_ = "Column-major MatrixAccessProxy dense vector multiplication assignment";

      DV tmp( 3UL );
      tmp[0] = 1;
      tmp[1] = 2;
      tmp[2] = 3;
      DVM mat( 2UL, 2UL, 1UL );
      mat(0,1) = tmp;

      mat(0,1) *= tmp;

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 0UL );
      checkSize    ( mat(0,1), 3UL );
      checkCapacity( mat(0,1), 3UL );
      checkNonZeros( mat(0,1), 3UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 0UL );

      if( mat(0,1)[0] != 1 || mat(0,1)[1] != 4 || mat(0,1)[2] != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat(0,1) << "\n"
             << "   Expected result:\n( 1 4 9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse vector multiplication assignment
   //=====================================================================================

   {
      test_ = "Column-major MatrixAccessProxy sparse vector multiplication assignment";

      SV tmp( 3UL );
      tmp[1] = 2;
      DVM mat( 2UL, 2UL, 1UL );
      mat(0,1) = tmp;

      mat(0,1) *= tmp;

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 0UL );
      checkSize    ( mat(0,1), 3UL );
      checkCapacity( mat(0,1), 3UL );
      checkNonZeros( mat(0,1), 1UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 0UL );

      if( mat(0,1)[0] != 0 || mat(0,1)[1] != 4 || mat(0,1)[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat(0,1) << "\n"
             << "   Expected result:\n( 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the MatrixAccessProxy division assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the division assignment operators of the MatrixAccessProxy
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ProxyTest::testDivAssign()
{
   //=====================================================================================
   // Row-major dense vector division assignment
   //=====================================================================================

   {
      test_ = "Row-major MatrixAccessProxy dense vector division assignment";

      DV tmp( 3UL );
      tmp[0] = 1;
      tmp[1] = 2;
      tmp[2] = 3;
      DVM mat( 2UL, 2UL, 1UL );
      mat(0,1) = tmp;

      mat(0,1) /= tmp;

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 0UL );
      checkSize    ( mat(0,1), 3UL );
      checkCapacity( mat(0,1), 3UL );
      checkNonZeros( mat(0,1), 3UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 0UL );

      if( mat(0,1)[0] != 1 || mat(0,1)[1] != 1 || mat(0,1)[2] != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat(0,1) << "\n"
             << "   Expected result:\n( 1 1 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense vector division assignment
   //=====================================================================================

   {
      test_ = "Column-major MatrixAccessProxy dense vector division assignment";

      DV tmp( 3UL );
      tmp[0] = 1;
      tmp[1] = 2;
      tmp[2] = 3;
      DVM mat( 2UL, 2UL, 1UL );
      mat(0,1) = tmp;

      mat(0,1) /= tmp;

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 0UL );
      checkSize    ( mat(0,1), 3UL );
      checkCapacity( mat(0,1), 3UL );
      checkNonZeros( mat(0,1), 3UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 0UL );

      if( mat(0,1)[0] != 1 || mat(0,1)[1] != 1 || mat(0,1)[2] != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat(0,1) << "\n"
             << "   Expected result:\n( 1 1 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the MatrixAccessProxy modulo assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the modulo assignment operators of the MatrixAccessProxy
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ProxyTest::testModAssign()
{
   //=====================================================================================
   // Row-major dense vector cross product assignment
   //=====================================================================================

   {
      test_ = "Row-major MatrixAccessProxy dense vector cross product assignment";

      DV tmp( 3UL );
      tmp[0] = 1;
      tmp[1] = 2;
      tmp[2] = 3;
      DVM mat( 2UL, 2UL, 1UL );
      mat(0,1) = tmp;

      mat(0,1) %= tmp;

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 0UL );
      checkSize    ( mat(0,1), 3UL );
      checkCapacity( mat(0,1), 3UL );
      checkNonZeros( mat(0,1), 0UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 0UL );

      if( mat(0,1)[0] != 0 || mat(0,1)[1] != 0 || mat(0,1)[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat(0,1) << "\n"
             << "   Expected result:\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse vector cross product assignment
   //=====================================================================================

   {
      test_ = "Row-major MatrixAccessProxy sparse vector cross product assignment";

      SV tmp( 3UL );
      tmp[1] = 2;
      DVM mat( 2UL, 2UL, 1UL );
      mat(0,1) = tmp;

      mat(0,1) %= tmp;

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 0UL );
      checkSize    ( mat(0,1), 3UL );
      checkCapacity( mat(0,1), 1UL );
      checkNonZeros( mat(0,1), 0UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 0UL );

      if( mat(0,1)[0] != 0 || mat(0,1)[1] != 0 || mat(0,1)[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat(0,1) << "\n"
             << "   Expected result:\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense vector cross product assignment
   //=====================================================================================

   {
      test_ = "Column-major MatrixAccessProxy dense vector cross product assignment";

      DV tmp( 3UL );
      tmp[0] = 1;
      tmp[1] = 2;
      tmp[2] = 3;
      DVM mat( 2UL, 2UL, 1UL );
      mat(0,1) = tmp;

      mat(0,1) %= tmp;

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 0UL );
      checkSize    ( mat(0,1), 3UL );
      checkCapacity( mat(0,1), 3UL );
      checkNonZeros( mat(0,1), 0UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 0UL );

      if( mat(0,1)[0] != 0 || mat(0,1)[1] != 0 || mat(0,1)[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat(0,1) << "\n"
             << "   Expected result:\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse vector cross product assignment
   //=====================================================================================

   {
      test_ = "Column-major MatrixAccessProxy sparse vector cross product assignment";

      SV tmp( 3UL );
      tmp[1] = 2;
      DVM mat( 2UL, 2UL, 1UL );
      mat(0,1) = tmp;

      mat(0,1) %= tmp;

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 0UL );
      checkSize    ( mat(0,1), 3UL );
      checkCapacity( mat(0,1), 1UL );
      checkNonZeros( mat(0,1), 0UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 0UL );

      if( mat(0,1)[0] != 0 || mat(0,1)[1] != 0 || mat(0,1)[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat(0,1) << "\n"
             << "   Expected result:\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of all MatrixAccessProxy (self-)scaling operations.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all available ways to scale an instance of the MatrixAccessProxy
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ProxyTest::testScaling()
{
   //=====================================================================================
   // Row-major self-scaling (v*=s)
   //=====================================================================================

   {
      test_ = "Row-major MatrixAccessProxy self-scaling (v*=s)";

      DVM mat( 2UL, 2UL, 1UL );
      mat(1,1) = DV( 1UL, 2 );

      mat(1,1) *= 2;

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 0UL );
      checkSize    ( mat(0,1), 0UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 1UL );
      checkCapacity( mat(1,1), 1UL );
      checkNonZeros( mat(1,1), 1UL );

      if( mat(1,1)[0] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << mat(1,1) << "\n"
             << "   Expected result:\n( 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major self-scaling (v/=s)
   //=====================================================================================

   {
      test_ = "Row-major MatrixAccessProxy self-scaling (v*=s)";

      DVM mat( 2UL, 2UL, 1UL );
      mat(1,1) = DV( 1UL, 2 );

      mat(1,1) /= 2;

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 0UL );
      checkSize    ( mat(0,1), 0UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 1UL );
      checkCapacity( mat(1,1), 1UL );
      checkNonZeros( mat(1,1), 1UL );

      if( mat(1,1)[0] != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << mat(1,1) << "\n"
             << "   Expected result:\n( 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major MatrixAccessProxy::scale()
   //=====================================================================================

   {
      test_ = "Row-major MatrixAccessProxy::scale()";

      DVM mat( 2UL, 2UL, 1UL );
      mat(1,1) = DV( 1UL, 2 );

      mat(1,1).scale( 2 );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 0UL );
      checkSize    ( mat(0,1), 0UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 1UL );
      checkCapacity( mat(1,1), 1UL );
      checkNonZeros( mat(1,1), 1UL );

      if( mat(1,1)[0] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << mat(1,1) << "\n"
             << "   Expected result:\n( 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (v*=s)
   //=====================================================================================

   {
      test_ = "Column-major MatrixAccessProxy self-scaling (v*=s)";

      ODVM mat( 2UL, 2UL, 1UL );
      mat(1,1) = DV( 1UL, 2 );

      mat(1,1) *= 2;

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 0UL );
      checkSize    ( mat(0,1), 0UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 1UL );
      checkCapacity( mat(1,1), 1UL );
      checkNonZeros( mat(1,1), 1UL );

      if( mat(1,1)[0] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << mat(1,1) << "\n"
             << "   Expected result:\n( 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (v/=s)
   //=====================================================================================

   {
      test_ = "Column-major MatrixAccessProxy self-scaling (v*=s)";

      ODVM mat( 2UL, 2UL, 1UL );
      mat(1,1) = DV( 1UL, 2 );

      mat(1,1) /= 2;

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 0UL );
      checkSize    ( mat(0,1), 0UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 1UL );
      checkCapacity( mat(1,1), 1UL );
      checkNonZeros( mat(1,1), 1UL );

      if( mat(1,1)[0] != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << mat(1,1) << "\n"
             << "   Expected result:\n( 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major MatrixAccessProxy::scale()
   //=====================================================================================

   {
      test_ = "Column-major MatrixAccessProxy::scale()";

      ODVM mat( 2UL, 2UL, 1UL );
      mat(1,1) = DV( 1UL, 2 );

      mat(1,1).scale( 2 );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 0UL );
      checkSize    ( mat(0,1), 0UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 1UL );
      checkCapacity( mat(1,1), 1UL );
      checkNonZeros( mat(1,1), 1UL );

      if( mat(1,1)[0] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << mat(1,1) << "\n"
             << "   Expected result:\n( 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the MatrixAccessProxy subscript operator.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of adding and accessing elements via the subscript operator
// of the MatrixAccessProxy class template. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
void ProxyTest::testSubscript()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major MatrixAccessProxy::operator[]";

      DVM mat( 2UL, 2UL, 1UL );
      mat(1,1) = DV( 1UL, 2 );
      mat(1,1)[0] = 3;

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 0UL );
      checkSize    ( mat(0,1), 0UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 1UL );
      checkCapacity( mat(1,1), 1UL );
      checkNonZeros( mat(1,1), 1UL );

      if( mat(1,1)[0] != DV( 1UL, 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << mat(1,1) << "\n"
             << "   Expected result:\n( 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major MatrixAccessProxy::operator[]";

      ODVM mat( 2UL, 2UL, 1UL );
      mat(1,1) = DV( 1UL, 2 );
      mat(1,1)[0] = 3;

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 0UL );
      checkSize    ( mat(0,1), 0UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 1UL );
      checkCapacity( mat(1,1), 1UL );
      checkNonZeros( mat(1,1), 1UL );

      if( mat(1,1)[0] != DV( 1UL, 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << mat(1,1) << "\n"
             << "   Expected result:\n( 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the MatrixAccessProxy function call operator.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of adding and accessing elements via the function call operator
// of the MatrixAccessProxy class template. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
void ProxyTest::testFunctionCall()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major MatrixAccessProxy::operator()";

      DMM mat( 2UL, 2UL, 1UL );
      mat(1,1) = DM( 1UL, 1UL, 2 );
      mat(1,1)(0,0) = 3;

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0), 0UL );
      checkColumns ( mat(0,0), 0UL );
      checkRows    ( mat(0,1), 0UL );
      checkColumns ( mat(0,1), 0UL );
      checkRows    ( mat(1,0), 0UL );
      checkColumns ( mat(1,0), 0UL );
      checkRows    ( mat(1,1), 1UL );
      checkColumns ( mat(1,1), 1UL );
      checkCapacity( mat(1,1), 1UL );
      checkNonZeros( mat(1,1), 1UL );

      if( mat(1,1)(0,0) != DM( 1UL, 1UL, 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << mat(1,1) << "\n"
             << "   Expected result:\n( 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major MatrixAccessProxy::operator()";

      ODMM mat( 2UL, 2UL, 1UL );
      mat(1,1) = DM( 1UL, 1UL, 2 );
      mat(1,1)(0,0) = 3;

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0), 0UL );
      checkColumns ( mat(0,0), 0UL );
      checkRows    ( mat(0,1), 0UL );
      checkColumns ( mat(0,1), 0UL );
      checkRows    ( mat(1,0), 0UL );
      checkColumns ( mat(1,0), 0UL );
      checkRows    ( mat(1,1), 1UL );
      checkColumns ( mat(1,1), 1UL );
      checkCapacity( mat(1,1), 1UL );
      checkNonZeros( mat(1,1), 1UL );

      if( mat(1,1)(0,0) != DM( 1UL, 1UL, 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << mat(1,1) << "\n"
             << "   Expected result:\n( 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the MatrixAccessProxy iterator implementation.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the iterator implementation of the MatrixAccessProxy class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ProxyTest::testIterator()
{
   //=====================================================================================
   // Row-major matrix tests with vector elements
   //=====================================================================================

   {
      DVM mat( 2UL, 2UL, 1UL );
      mat(0,1) = DV( 4UL, 4 );

      // Counting the number of elements via Iterator (end-begin)
      {
         test_ = "Row-major MatrixAccessProxy::begin() and MatrixAccessProxy::end()";

         const ptrdiff_t number( end( mat(0,1) ) - begin( mat(0,1) ) );

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

      // Counting the number of elements via ConstIterator (end-begin)
      {
         test_ = "Row-major MatrixAccessProxy::cbegin() and MatrixAccessProxy::cend()";

         const ptrdiff_t number( cend( mat(0,1) ) - cbegin( mat(0,1) ) );

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
   }


   //=====================================================================================
   // Row-major matrix tests with matrix elements
   //=====================================================================================

   {
      DMM mat( 2UL, 2UL, 1UL );
      mat(0,1) = DM( 4UL, 4UL, 4 );

      // Counting the number of elements via Iterator (end-begin)
      {
         test_ = "Row-major MatrixAccessProxy::begin( size_t ) and MatrixAccessProxy::end( size_t )";

         const ptrdiff_t number( end( mat(0,1), 1UL ) - begin( mat(0,1), 1UL ) );

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

      // Counting the number of elements via ConstIterator (end-begin)
      {
         test_ = "Row-major MatrixAccessProxy::cbegin( size_t ) and MatrixAccessProxy::cend( size_t )";

         const ptrdiff_t number( cend( mat(0,1), 1UL ) - cbegin( mat(0,1), 1UL ) );

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
   }


   //=====================================================================================
   // Column-major matrix tests with vector elements
   //=====================================================================================

   {
      ODVM mat( 2UL, 2UL, 1UL );
      mat(0,1) = DV( 4UL, 4 );

      // Counting the number of elements via Iterator (end-begin)
      {
         test_ = "Column-major MatrixAccessProxy::begin() and MatrixAccessProxy::end()";

         const ptrdiff_t number( end( mat(0,1) ) - begin( mat(0,1) ) );

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

      // Counting the number of elements via ConstIterator (end-begin)
      {
         test_ = "Column-major MatrixAccessProxy::cbegin() and MatrixAccessProxy::cend()";

         const ptrdiff_t number( cend( mat(0,1) ) - cbegin( mat(0,1) ) );

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
   }


   //=====================================================================================
   // Column-major matrix tests with matrix elements
   //=====================================================================================

   {
      ODMM mat( 2UL, 2UL, 1UL );
      mat(0,1) = DM( 4UL, 4UL, 4 );

      // Counting the number of elements via Iterator (end-begin)
      {
         test_ = "Column-major MatrixAccessProxy::begin( size_t ) and MatrixAccessProxy::end( size_t )";

         const ptrdiff_t number( end( mat(0,1), 1UL ) - begin( mat(0,1), 1UL ) );

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

      // Counting the number of elements via ConstIterator (end-begin)
      {
         test_ = "Column-major MatrixAccessProxy::cbegin( size_t ) and MatrixAccessProxy::cend( size_t )";

         const ptrdiff_t number( cend( mat(0,1), 1UL ) - cbegin( mat(0,1), 1UL ) );

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
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c nonZeros() member functions of the MatrixAccessProxy class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c nonZeros() member functions of the MatrixAccessProxy
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ProxyTest::testNonZeros()
{
   //=====================================================================================
   // Row-major matrix tests with vector elements
   //=====================================================================================

   {
      test_ = "Row-major MatrixAccessProxy::nonZeros()";

      DVM mat( 2UL, 2UL, 1UL );
      mat(1,1) = DV( 8UL, 8 );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 0UL );
      checkSize    ( mat(0,1), 0UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 8UL );
      checkCapacity( mat(1,1), 8UL );
      checkNonZeros( mat(1,1), 8UL );
   }


   //=====================================================================================
   // Row-major matrix tests with matrix elements
   //=====================================================================================

   {
      test_ = "Row-major MatrixAccessProxy::nonZeros()";

      DMM mat( 2UL, 2UL, 1UL );
      mat(1,1) = DM( 3UL, 3UL, 3 );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0), 0UL );
      checkColumns ( mat(0,0), 0UL );
      checkRows    ( mat(0,1), 0UL );
      checkColumns ( mat(0,1), 0UL );
      checkRows    ( mat(1,0), 0UL );
      checkColumns ( mat(1,0), 0UL );
      checkRows    ( mat(1,1), 3UL );
      checkColumns ( mat(1,1), 3UL );
      checkCapacity( mat(1,1), 9UL );
      checkNonZeros( mat(1,1), 9UL );
      checkNonZeros( mat(1,1), 0UL, 3UL );
      checkNonZeros( mat(1,1), 1UL, 3UL );
      checkNonZeros( mat(1,1), 2UL, 3UL );
   }


   //=====================================================================================
   // Column-major matrix tests with vector elements
   //=====================================================================================

   {
      test_ = "Column-major MatrixAccessProxy::nonZeros()";

      ODVM mat( 2UL, 2UL, 1UL );
      mat(1,1) = DV( 8UL, 8 );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 0UL );
      checkSize    ( mat(0,1), 0UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 8UL );
      checkCapacity( mat(1,1), 8UL );
      checkNonZeros( mat(1,1), 8UL );
   }


   //=====================================================================================
   // Column-major matrix tests with matrix elements
   //=====================================================================================

   {
      test_ = "Row-major MatrixAccessProxy::nonZeros()";

      ODMM mat( 2UL, 2UL, 1UL );
      mat(1,1) = DM( 3UL, 3UL, 3 );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0), 0UL );
      checkColumns ( mat(0,0), 0UL );
      checkRows    ( mat(0,1), 0UL );
      checkColumns ( mat(0,1), 0UL );
      checkRows    ( mat(1,0), 0UL );
      checkColumns ( mat(1,0), 0UL );
      checkRows    ( mat(1,1), 3UL );
      checkColumns ( mat(1,1), 3UL );
      checkCapacity( mat(1,1), 9UL );
      checkNonZeros( mat(1,1), 9UL );
      checkNonZeros( mat(1,1), 0UL, 3UL );
      checkNonZeros( mat(1,1), 1UL, 3UL );
      checkNonZeros( mat(1,1), 2UL, 3UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c reset() member functions of the MatrixAccessProxy class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c reset() member functions of the MatrixAccessProxy
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ProxyTest::testReset()
{
   //=====================================================================================
   // Row-major matrix tests with vector elements
   //=====================================================================================

   {
      test_ = "Row-major MatrixAccessProxy::reset()";

      DVM mat( 2UL, 2UL, 1UL );
      mat(0,1) = DV( 8UL, 8 );
      mat(0,1).reset();

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 0UL );
      checkSize    ( mat(0,1), 8UL );
      checkCapacity( mat(0,1), 8UL );
      checkNonZeros( mat(0,1), 0UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 0UL );
   }


   //=====================================================================================
   // Row-major matrix tests with matrix elements
   //=====================================================================================

   {
      test_ = "Row-major MatrixAccessProxy::reset( size_t )";

      DMM mat( 2UL, 2UL, 1UL );
      mat(0,1) = DM( 3UL, 3UL, 3 );
      mat(0,1).reset( 1UL );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0), 0UL );
      checkColumns ( mat(0,0), 0UL );
      checkRows    ( mat(0,1), 3UL );
      checkColumns ( mat(0,1), 3UL );
      checkCapacity( mat(0,1), 9UL );
      checkNonZeros( mat(0,1), 6UL );
      checkNonZeros( mat(0,1), 0UL, 3UL );
      checkNonZeros( mat(0,1), 1UL, 0UL );
      checkNonZeros( mat(0,1), 2UL, 3UL );
      checkRows    ( mat(1,0), 0UL );
      checkColumns ( mat(1,0), 0UL );
      checkRows    ( mat(1,1), 0UL );
      checkColumns ( mat(1,1), 0UL );
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major MatrixAccessProxy::reset()";

      ODVM mat( 2UL, 2UL, 1UL );
      mat(0,1) = DV( 8UL, 8 );
      mat(0,1).reset();

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 0UL );
      checkSize    ( mat(0,1), 8UL );
      checkCapacity( mat(0,1), 8UL );
      checkNonZeros( mat(0,1), 0UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 0UL );
   }


   //=====================================================================================
   // Column-major matrix tests with matrix elements
   //=====================================================================================

   {
      test_ = "Column-major MatrixAccessProxy::reset( size_t )";

      ODMM mat( 2UL, 2UL, 1UL );
      mat(0,1) = DM( 3UL, 3UL, 3 );
      mat(0,1).reset( 1UL );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0), 0UL );
      checkColumns ( mat(0,0), 0UL );
      checkRows    ( mat(0,1), 3UL );
      checkColumns ( mat(0,1), 3UL );
      checkCapacity( mat(0,1), 9UL );
      checkNonZeros( mat(0,1), 6UL );
      checkNonZeros( mat(0,1), 0UL, 3UL );
      checkNonZeros( mat(0,1), 1UL, 0UL );
      checkNonZeros( mat(0,1), 2UL, 3UL );
      checkRows    ( mat(1,0), 0UL );
      checkColumns ( mat(1,0), 0UL );
      checkRows    ( mat(1,1), 0UL );
      checkColumns ( mat(1,1), 0UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c clear() member function of the MatrixAccessProxy class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c clear() member function of the MatrixAccessProxy
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ProxyTest::testClear()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major MatrixAccessProxy::clear()";

      DVM mat( 2UL, 2UL, 1UL );
      mat(0,1) = DV( 8UL, 8 );
      mat(0,1).clear();

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 0UL );

      checkSize( mat(0,0), 0UL );
      checkSize( mat(0,1), 0UL );
      checkSize( mat(1,0), 0UL );
      checkSize( mat(1,1), 0UL );
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major MatrixAccessProxy::clear()";

      ODVM mat( 2UL, 2UL, 1UL );
      mat(0,1) = DV( 8UL, 8 );
      mat(0,1).clear();

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 0UL );

      checkSize( mat(0,0), 0UL );
      checkSize( mat(0,1), 0UL );
      checkSize( mat(1,0), 0UL );
      checkSize( mat(1,1), 0UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c resize() member functions of the MatrixAccessProxy class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c resize() member functions of the MatrixAccessProxy
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ProxyTest::testResize()
{
   //=====================================================================================
   // Row-major matrix tests with vector elements
   //=====================================================================================

   {
      test_ = "Row-major MatrixAccessProxy::resize( size_t )";

      DVM mat( 2UL, 2UL, 1UL );
      mat(0,0).resize( 5UL );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 5UL );
      checkCapacity( mat(0,0), 5UL );
      checkSize    ( mat(0,1), 0UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 0UL );
   }

   {
      test_ = "Row-major resize( MatrixAccessProxy, size_t )";

      DVM mat( 2UL, 2UL, 1UL );
      resize( mat(0,0), 5UL );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 5UL );
      checkCapacity( mat(0,0), 5UL );
      checkSize    ( mat(0,1), 0UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 0UL );
   }


   //=====================================================================================
   // Row-major matrix tests with matrix elements
   //=====================================================================================

   {
      test_ = "Row-major MatrixAccessProxy::resize( size_t, size_t )";

      DMM mat( 2UL, 2UL, 1UL );
      mat(0,0).resize( 5UL, 5UL );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0),  5UL );
      checkColumns ( mat(0,0),  5UL );
      checkCapacity( mat(0,0), 25UL );
      checkRows    ( mat(0,1),  0UL );
      checkColumns ( mat(0,1),  0UL );
      checkRows    ( mat(1,0),  0UL );
      checkColumns ( mat(1,0),  0UL );
      checkRows    ( mat(1,1),  0UL );
      checkColumns ( mat(1,1),  0UL );
   }

   {
      test_ = "Row-major resize( MatrixAccessProxy, size_t, size_t )";

      DMM mat( 2UL, 2UL, 1UL );
      resize( mat(0,0), 5UL, 5UL );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0),  5UL );
      checkColumns ( mat(0,0),  5UL );
      checkCapacity( mat(0,0), 25UL );
      checkRows    ( mat(0,1),  0UL );
      checkColumns ( mat(0,1),  0UL );
      checkRows    ( mat(1,0),  0UL );
      checkColumns ( mat(1,0),  0UL );
      checkRows    ( mat(1,1),  0UL );
      checkColumns ( mat(1,1),  0UL );
   }


   //=====================================================================================
   // Column-major matrix tests with vector elements
   //=====================================================================================

   {
      test_ = "Row-major MatrixAccessProxy::resize( size_t )";

      ODVM mat( 2UL, 2UL, 1UL );
      mat(0,0).resize( 5UL );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 5UL );
      checkCapacity( mat(0,0), 5UL );
      checkSize    ( mat(0,1), 0UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 0UL );
   }

   {
      test_ = "Row-major resize( MatrixAccessProxy, size_t )";

      ODVM mat( 2UL, 2UL, 1UL );
      resize( mat(0,0), 5UL );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 5UL );
      checkCapacity( mat(0,0), 5UL );
      checkSize    ( mat(0,1), 0UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 0UL );
   }


   //=====================================================================================
   // Column-major matrix tests with matrix elements
   //=====================================================================================

   {
      test_ = "Column-major MatrixAccessProxy::resize( size_t, size_t )";

      ODMM mat( 2UL, 2UL, 1UL );
      mat(0,0).resize( 5UL, 5UL );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0),  5UL );
      checkColumns ( mat(0,0),  5UL );
      checkCapacity( mat(0,0), 25UL );
      checkRows    ( mat(0,1),  0UL );
      checkColumns ( mat(0,1),  0UL );
      checkRows    ( mat(1,0),  0UL );
      checkColumns ( mat(1,0),  0UL );
      checkRows    ( mat(1,1),  0UL );
      checkColumns ( mat(1,1),  0UL );
   }

   {
      test_ = "Column-major MatrixAccessProxy::resize( size_t, size_t )";

      ODMM mat( 2UL, 2UL, 1UL );
      resize( mat(0,0), 5UL, 5UL );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0),  5UL );
      checkColumns ( mat(0,0),  5UL );
      checkCapacity( mat(0,0), 25UL );
      checkRows    ( mat(0,1),  0UL );
      checkColumns ( mat(0,1),  0UL );
      checkRows    ( mat(1,0),  0UL );
      checkColumns ( mat(1,0),  0UL );
      checkRows    ( mat(1,1),  0UL );
      checkColumns ( mat(1,1),  0UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c extend() member functions of the MatrixAccessProxy class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c extend() member functions of the MatrixAccessProxy
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ProxyTest::testExtend()
{
   //=====================================================================================
   // Row-major matrix tests with vector elements
   //=====================================================================================

   {
      test_ = "Row-major MatrixAccessProxy::extend( size_t )";

      DVM mat( 2UL, 2UL, 1UL );
      mat(0,0).extend( 5UL );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 5UL );
      checkCapacity( mat(0,0), 5UL );
      checkSize    ( mat(0,1), 0UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 0UL );
   }


   //=====================================================================================
   // Row-major matrix tests with matrix elements
   //=====================================================================================

   {
      test_ = "Row-major MatrixAccessProxy::extend( size_t, size_t )";

      DMM mat( 2UL, 2UL, 1UL );
      mat(0,0).extend( 5UL, 5UL );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0),  5UL );
      checkColumns ( mat(0,0),  5UL );
      checkCapacity( mat(0,0), 25UL );
      checkRows    ( mat(0,1),  0UL );
      checkColumns ( mat(0,1),  0UL );
      checkRows    ( mat(1,0),  0UL );
      checkColumns ( mat(1,0),  0UL );
      checkRows    ( mat(1,1),  0UL );
      checkColumns ( mat(1,1),  0UL );
   }


   //=====================================================================================
   // Column-major matrix tests with vector elements
   //=====================================================================================

   {
      test_ = "Row-major MatrixAccessProxy::extend( size_t )";

      ODVM mat( 2UL, 2UL, 1UL );
      mat(0,0).extend( 5UL );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 5UL );
      checkCapacity( mat(0,0), 5UL );
      checkSize    ( mat(0,1), 0UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 0UL );
   }


   //=====================================================================================
   // Column-major matrix tests with matrix elements
   //=====================================================================================

   {
      test_ = "Column-major MatrixAccessProxy::extend( size_t, size_t )";

      ODMM mat( 2UL, 2UL, 1UL );
      mat(0,0).extend( 5UL, 5UL );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0),  5UL );
      checkColumns ( mat(0,0),  5UL );
      checkCapacity( mat(0,0), 25UL );
      checkRows    ( mat(0,1),  0UL );
      checkColumns ( mat(0,1),  0UL );
      checkRows    ( mat(1,0),  0UL );
      checkColumns ( mat(1,0),  0UL );
      checkRows    ( mat(1,1),  0UL );
      checkColumns ( mat(1,1),  0UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c reserve() member functions of the MatrixAccessProxy class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c reserve() member functions of the MatrixAccessProxy
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ProxyTest::testReserve()
{
   //=====================================================================================
   // Row-major matrix tests with vector elements
   //=====================================================================================

   {
      test_ = "Row-major MatrixAccessProxy::reserve( size_t )";

      DVM mat( 2UL, 2UL, 1UL );
      mat(0,0).resize( 5UL );
      mat(0,0).reserve( 10UL );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0),  5UL );
      checkCapacity( mat(0,0), 10UL );
      checkSize    ( mat(0,1),  0UL );
      checkSize    ( mat(1,0),  0UL );
      checkSize    ( mat(1,1),  0UL );
   }


   //=====================================================================================
   // Row-major matrix tests with matrix elements
   //=====================================================================================

   {
      test_ = "Row-major MatrixAccessProxy::reserve( size_t, size_t )";

      SMM mat( 2UL, 2UL, 1UL );
      mat(0,0) = SM( 2UL, 2UL, 1UL );
      mat(0,0).reserve( 0UL, 1UL );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0), 2UL );
      checkColumns ( mat(0,0), 2UL );
      checkCapacity( mat(0,0), 1UL );
      checkCapacity( mat(0,0), 0UL, 1UL );
      checkCapacity( mat(0,0), 1UL, 0UL );
      checkRows    ( mat(0,1), 0UL );
      checkColumns ( mat(0,1), 0UL );
      checkRows    ( mat(1,0), 0UL );
      checkColumns ( mat(1,0), 0UL );
      checkRows    ( mat(1,1), 0UL );
      checkColumns ( mat(1,1), 0UL );
   }


   //=====================================================================================
   // Column-major matrix tests with vector elements
   //=====================================================================================

   {
      test_ = "Column-major MatrixAccessProxy::reserve( size_t )";

      ODVM mat( 2UL, 2UL, 1UL );
      mat(0,0).resize( 5UL );
      mat(0,0).reserve( 10UL );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0),  5UL );
      checkCapacity( mat(0,0), 10UL );
      checkSize    ( mat(0,1),  0UL );
      checkSize    ( mat(1,0),  0UL );
      checkSize    ( mat(1,1),  0UL );
   }


   //=====================================================================================
   // Column-major matrix tests with matrix elements
   //=====================================================================================

   {
      test_ = "Column-major MatrixAccessProxy::reserve( size_t, size_t )";

      OSMM mat( 2UL, 2UL, 1UL );
      mat(0,0) = SM( 2UL, 2UL, 1UL );
      mat(0,0).reserve( 0UL, 1UL );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0), 2UL );
      checkColumns ( mat(0,0), 2UL );
      checkCapacity( mat(0,0), 1UL );
      checkCapacity( mat(0,0), 0UL, 1UL );
      checkCapacity( mat(0,0), 1UL, 0UL );
      checkRows    ( mat(0,1), 0UL );
      checkColumns ( mat(0,1), 0UL );
      checkRows    ( mat(1,0), 0UL );
      checkColumns ( mat(1,0), 0UL );
      checkRows    ( mat(1,1), 0UL );
      checkColumns ( mat(1,1), 0UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c trim() member functions of the VectorAccessProxy class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c trim() member functions of the VectorAccessProxy
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ProxyTest::testTrim()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major MatrixAccessProxy::trim()";

      SMM mat( 2UL, 2UL, 4UL );
      mat(0,0).resize( 2UL, 2UL );
      mat(0,0).reserve( 10UL );
      mat(0,0).reserve( 0UL, 6UL );
      mat(0,0).reserve( 1UL, 4UL );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 4UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0),  2UL );
      checkColumns ( mat(0,0),  2UL );
      checkCapacity( mat(0,0), 10UL );
      checkCapacity( mat(0,0),  0UL, 6UL );
      checkCapacity( mat(0,0),  1UL, 4UL );
      checkRows    ( mat(0,1),  0UL );
      checkColumns ( mat(0,1),  0UL );
      checkRows    ( mat(1,0),  0UL );
      checkColumns ( mat(1,0),  0UL );
      checkRows    ( mat(1,1),  0UL );
      checkColumns ( mat(1,1),  0UL );

      mat(0,0).trim();

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 4UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0),  2UL );
      checkColumns ( mat(0,0),  2UL );
      checkCapacity( mat(0,0), 10UL );
      checkCapacity( mat(0,0),  0UL, 0UL );
      checkCapacity( mat(0,0),  1UL, 0UL );
      checkRows    ( mat(0,1),  0UL );
      checkColumns ( mat(0,1),  0UL );
      checkRows    ( mat(1,0),  0UL );
      checkColumns ( mat(1,0),  0UL );
      checkRows    ( mat(1,1),  0UL );
      checkColumns ( mat(1,1),  0UL );
   }

   {
      test_ = "Row-major MatrixAccessProxy::trim( size_t )";

      SMM mat( 2UL, 2UL, 4UL );
      mat(0,0).resize( 2UL, 2UL );
      mat(0,0).reserve( 10UL );
      mat(0,0).reserve( 0UL, 6UL );
      mat(0,0).reserve( 1UL, 4UL );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 4UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0),  2UL );
      checkColumns ( mat(0,0),  2UL );
      checkCapacity( mat(0,0), 10UL );
      checkCapacity( mat(0,0),  0UL, 6UL );
      checkCapacity( mat(0,0),  1UL, 4UL );
      checkRows    ( mat(0,1),  0UL );
      checkColumns ( mat(0,1),  0UL );
      checkRows    ( mat(1,0),  0UL );
      checkColumns ( mat(1,0),  0UL );
      checkRows    ( mat(1,1),  0UL );
      checkColumns ( mat(1,1),  0UL );

      mat(0,0).trim( 0UL );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 4UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0),  2UL );
      checkColumns ( mat(0,0),  2UL );
      checkCapacity( mat(0,0), 10UL );
      checkCapacity( mat(0,0),  0UL, 0UL );
      checkCapacity( mat(0,0),  1UL, 4UL );
      checkRows    ( mat(0,1),  0UL );
      checkColumns ( mat(0,1),  0UL );
      checkRows    ( mat(1,0),  0UL );
      checkColumns ( mat(1,0),  0UL );
      checkRows    ( mat(1,1),  0UL );
      checkColumns ( mat(1,1),  0UL );
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major MatrixAccessProxy::trim()";

      OSMM mat( 2UL, 2UL, 4UL );
      mat(0,0).resize( 2UL, 2UL );
      mat(0,0).reserve( 10UL );
      mat(0,0).reserve( 0UL, 6UL );
      mat(0,0).reserve( 1UL, 4UL );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 4UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0),  2UL );
      checkColumns ( mat(0,0),  2UL );
      checkCapacity( mat(0,0), 10UL );
      checkCapacity( mat(0,0),  0UL, 6UL );
      checkCapacity( mat(0,0),  1UL, 4UL );
      checkRows    ( mat(0,1),  0UL );
      checkColumns ( mat(0,1),  0UL );
      checkRows    ( mat(1,0),  0UL );
      checkColumns ( mat(1,0),  0UL );
      checkRows    ( mat(1,1),  0UL );
      checkColumns ( mat(1,1),  0UL );

      mat(0,0).trim();

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 4UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0),  2UL );
      checkColumns ( mat(0,0),  2UL );
      checkCapacity( mat(0,0), 10UL );
      checkCapacity( mat(0,0),  0UL, 0UL );
      checkCapacity( mat(0,0),  1UL, 0UL );
      checkRows    ( mat(0,1),  0UL );
      checkColumns ( mat(0,1),  0UL );
      checkRows    ( mat(1,0),  0UL );
      checkColumns ( mat(1,0),  0UL );
      checkRows    ( mat(1,1),  0UL );
      checkColumns ( mat(1,1),  0UL );
   }

   {
      test_ = "Column-major MatrixAccessProxy::trim( size_t )";

      OSMM mat( 2UL, 2UL, 4UL );
      mat(0,0).resize( 2UL, 2UL );
      mat(0,0).reserve( 10UL );
      mat(0,0).reserve( 0UL, 6UL );
      mat(0,0).reserve( 1UL, 4UL );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 4UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0),  2UL );
      checkColumns ( mat(0,0),  2UL );
      checkCapacity( mat(0,0), 10UL );
      checkCapacity( mat(0,0),  0UL, 6UL );
      checkCapacity( mat(0,0),  1UL, 4UL );
      checkRows    ( mat(0,1),  0UL );
      checkColumns ( mat(0,1),  0UL );
      checkRows    ( mat(1,0),  0UL );
      checkColumns ( mat(1,0),  0UL );
      checkRows    ( mat(1,1),  0UL );
      checkColumns ( mat(1,1),  0UL );

      mat(0,0).trim( 0UL );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 4UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0),  2UL );
      checkColumns ( mat(0,0),  2UL );
      checkCapacity( mat(0,0), 10UL );
      checkCapacity( mat(0,0),  0UL, 0UL );
      checkCapacity( mat(0,0),  1UL, 4UL );
      checkRows    ( mat(0,1),  0UL );
      checkColumns ( mat(0,1),  0UL );
      checkRows    ( mat(1,0),  0UL );
      checkColumns ( mat(1,0),  0UL );
      checkRows    ( mat(1,1),  0UL );
      checkColumns ( mat(1,1),  0UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c swap() functionality of the MatrixAccessProxy class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c swap() function of the MatrixAccessProxy class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ProxyTest::testSwap()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   test_ = "Row-major MatrixAccessProxy swap";

   {
      DVM mat( 2UL, 2UL, 2UL );
      mat(0,0) = DV( 2UL, 0 );
      mat(1,1) = DV( 6UL, 0 );

      swap( mat(0,0), mat(1,1) );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 2UL );
      checkNonZeros( mat, 2UL );

      checkSize    ( mat(0,0), 6UL );
      checkCapacity( mat(0,0), 6UL );
      checkNonZeros( mat(0,0), 0UL );
      checkSize    ( mat(0,1), 0UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 2UL );
      checkCapacity( mat(1,1), 2UL );
      checkNonZeros( mat(1,1), 0UL );
   }

   {
      DVM mat( 2UL, 2UL, 1UL );
      mat(0,1) = DV( 2UL, 2 );
      DV tmp( 6UL, 6 );

      swap( mat(0,1), tmp );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 0UL );
      checkSize    ( mat(0,1), 6UL );
      checkCapacity( mat(0,1), 6UL );
      checkNonZeros( mat(0,1), 6UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 0UL );
      checkSize    ( tmp     , 2UL );
      checkCapacity( tmp     , 2UL );
      checkNonZeros( tmp     , 2UL );
   }

   {
      DVM mat( 2UL, 2UL, 1UL );
      mat(0,1) = DV( 2UL, 2 );
      DV tmp( 6UL, 6 );

      swap( tmp, mat(0,1) );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 0UL );
      checkSize    ( mat(0,1), 6UL );
      checkCapacity( mat(0,1), 6UL );
      checkNonZeros( mat(0,1), 6UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 0UL );
      checkSize    ( tmp     , 2UL );
      checkCapacity( tmp     , 2UL );
      checkNonZeros( tmp     , 2UL );
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   test_ = "Column-major MatrixAccessProxy swap";

   {
      ODVM mat( 2UL, 2UL, 2UL );
      mat(0,0) = DV( 2UL, 0 );
      mat(1,1) = DV( 6UL, 0 );

      swap( mat(0,0), mat(1,1) );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 2UL );
      checkNonZeros( mat, 2UL );

      checkSize    ( mat(0,0), 6UL );
      checkCapacity( mat(0,0), 6UL );
      checkNonZeros( mat(0,0), 0UL );
      checkSize    ( mat(0,1), 0UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 2UL );
      checkCapacity( mat(1,1), 2UL );
      checkNonZeros( mat(1,1), 0UL );
   }

   {
      ODVM mat( 2UL, 2UL, 1UL );
      mat(0,1) = DV( 2UL, 2 );
      DV tmp( 6UL, 6 );

      swap( mat(0,1), tmp );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 0UL );
      checkSize    ( mat(0,1), 6UL );
      checkCapacity( mat(0,1), 6UL );
      checkNonZeros( mat(0,1), 6UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 0UL );
      checkSize    ( tmp     , 2UL );
      checkCapacity( tmp     , 2UL );
      checkNonZeros( tmp     , 2UL );
   }

   {
      ODVM mat( 2UL, 2UL, 1UL );
      mat(0,1) = DV( 2UL, 2 );
      DV tmp( 6UL, 6 );

      swap( tmp, mat(0,1) );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 0UL );
      checkSize    ( mat(0,1), 6UL );
      checkCapacity( mat(0,1), 6UL );
      checkNonZeros( mat(0,1), 6UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 0UL );
      checkSize    ( tmp     , 2UL );
      checkCapacity( tmp     , 2UL );
      checkNonZeros( tmp     , 2UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c set() member functions of the MatrixAccessProxy class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c set() member functions of the MatrixAccessProxy
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ProxyTest::testSet()
{
   //=====================================================================================
   // Row-major matrix tests with vector elements
   //=====================================================================================

   {
      test_ = "Row-major MatrixAccessProxy::set( size_t, ElementType )";

      SVM mat( 2UL, 2UL, 1UL );
      mat(1,1) = SV( 3UL, 1UL );
      mat(1,1).set( 1UL, 5 );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 0UL );
      checkSize    ( mat(0,1), 0UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 3UL );
      checkCapacity( mat(1,1), 1UL );
      checkNonZeros( mat(1,1), 1UL );

      if( mat(1,1)[1] != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setting an element failed\n"
             << " Details:\n"
             << "   Result:\n" << mat(1,1) << "\n"
             << "   Expected result:\n( 0 5 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major matrix tests with matrix elements
   //=====================================================================================

   {
      test_ = "Row-major MatrixAccessProxy::set( size_t, size_t, ElementType )";

      SMM mat( 2UL, 2UL, 1UL );
      mat(1,1) = SM( 2UL, 2UL, 1UL );
      mat(1,1).set( 0UL, 1UL, 5 );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0), 0UL );
      checkColumns ( mat(0,0), 0UL );
      checkRows    ( mat(0,1), 0UL );
      checkColumns ( mat(0,1), 0UL );
      checkRows    ( mat(1,0), 0UL );
      checkColumns ( mat(1,0), 0UL );
      checkRows    ( mat(1,1), 2UL );
      checkColumns ( mat(1,1), 2UL );
      checkCapacity( mat(1,1), 1UL );
      checkNonZeros( mat(1,1), 1UL );

      if( mat(1,1)(0,1) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setting an element failed\n"
             << " Details:\n"
             << "   Result:\n" << mat(1,1) << "\n"
             << "   Expected result:\n( 0 5 )\n( 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests with vector elements
   //=====================================================================================

   {
      test_ = "Column-major MatrixAccessProxy::set( size_t, ElementType )";

      OSVM mat( 2UL, 2UL, 1UL );
      mat(1,1) = SV( 3UL, 1UL );
      mat(1,1).set( 1UL, 5 );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 0UL );
      checkSize    ( mat(0,1), 0UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 3UL );
      checkCapacity( mat(1,1), 1UL );
      checkNonZeros( mat(1,1), 1UL );

      if( mat(1,1)[1] != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setting an element failed\n"
             << " Details:\n"
             << "   Result:\n" << mat(1,1) << "\n"
             << "   Expected result:\n( 0 5 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests with matrix elements
   //=====================================================================================

   {
      test_ = "Column-major MatrixAccessProxy::set( size_t, size_t, ElementType )";

      OSMM mat( 2UL, 2UL, 1UL );
      mat(1,1) = SM( 2UL, 2UL, 1UL );
      mat(1,1).set( 0UL, 1UL, 5 );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0), 0UL );
      checkColumns ( mat(0,0), 0UL );
      checkRows    ( mat(0,1), 0UL );
      checkColumns ( mat(0,1), 0UL );
      checkRows    ( mat(1,0), 0UL );
      checkColumns ( mat(1,0), 0UL );
      checkRows    ( mat(1,1), 2UL );
      checkColumns ( mat(1,1), 2UL );
      checkCapacity( mat(1,1), 1UL );
      checkNonZeros( mat(1,1), 1UL );

      if( mat(1,1)(0,1) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setting an element failed\n"
             << " Details:\n"
             << "   Result:\n" << mat(1,1) << "\n"
             << "   Expected result:\n( 0 5 )\n( 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c insert() member functions of the MatrixAccessProxy class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c insert() member functions of the MatrixAccessProxy
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ProxyTest::testInsert()
{
   //=====================================================================================
   // Row-major matrix tests with vector elements
   //=====================================================================================

   {
      test_ = "Row-major MatrixAccessProxy::insert( size_t, ElementType )";

      SVM mat( 2UL, 2UL, 1UL );
      mat(1,1) = SV( 3UL, 1UL );
      mat(1,1).insert( 1UL, 5 );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 0UL );
      checkSize    ( mat(0,1), 0UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 3UL );
      checkCapacity( mat(1,1), 1UL );
      checkNonZeros( mat(1,1), 1UL );

      if( mat(1,1)[1] != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Inserting an element failed\n"
             << " Details:\n"
             << "   Result:\n" << mat(1,1) << "\n"
             << "   Expected result:\n( 0 5 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major matrix tests with matrix elements
   //=====================================================================================

   {
      test_ = "Row-major MatrixAccessProxy::insert( size_t, size_t, ElementType )";

      SMM mat( 2UL, 2UL, 1UL );
      mat(1,1) = SM( 2UL, 2UL, 1UL );
      mat(1,1).insert( 0UL, 1UL, 5 );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0), 0UL );
      checkColumns ( mat(0,0), 0UL );
      checkRows    ( mat(0,1), 0UL );
      checkColumns ( mat(0,1), 0UL );
      checkRows    ( mat(1,0), 0UL );
      checkColumns ( mat(1,0), 0UL );
      checkRows    ( mat(1,1), 2UL );
      checkColumns ( mat(1,1), 2UL );
      checkCapacity( mat(1,1), 1UL );
      checkNonZeros( mat(1,1), 1UL );

      if( mat(1,1)(0,1) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Inserting an element failed\n"
             << " Details:\n"
             << "   Result:\n" << mat(1,1) << "\n"
             << "   Expected result:\n( 0 5 )\n( 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests with vector elements
   //=====================================================================================

   {
      test_ = "Column-major MatrixAccessProxy::insert( size_t, ElementType )";

      OSVM mat( 2UL, 2UL, 1UL );
      mat(1,1) = SV( 3UL, 1UL );
      mat(1,1).insert( 1UL, 5 );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 0UL );
      checkSize    ( mat(0,1), 0UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 3UL );
      checkCapacity( mat(1,1), 1UL );
      checkNonZeros( mat(1,1), 1UL );

      if( mat(1,1)[1] != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Inserting an element failed\n"
             << " Details:\n"
             << "   Result:\n" << mat(1,1) << "\n"
             << "   Expected result:\n( 0 5 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests with matrix elements
   //=====================================================================================

   {
      test_ = "Column-major MatrixAccessProxy::insert( size_t, size_t, ElementType )";

      OSMM mat( 2UL, 2UL, 1UL );
      mat(1,1) = SM( 2UL, 2UL, 1UL );
      mat(1,1).insert( 0UL, 1UL, 5 );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0), 0UL );
      checkColumns ( mat(0,0), 0UL );
      checkRows    ( mat(0,1), 0UL );
      checkColumns ( mat(0,1), 0UL );
      checkRows    ( mat(1,0), 0UL );
      checkColumns ( mat(1,0), 0UL );
      checkRows    ( mat(1,1), 2UL );
      checkColumns ( mat(1,1), 2UL );
      checkCapacity( mat(1,1), 1UL );
      checkNonZeros( mat(1,1), 1UL );

      if( mat(1,1)(0,1) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Inserting an element failed\n"
             << " Details:\n"
             << "   Result:\n" << mat(1,1) << "\n"
             << "   Expected result:\n( 0 5 )\n( 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c append() member functions of the MatrixAccessProxy class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c append() member functions of the MatrixAccessProxy
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ProxyTest::testAppend()
{
   //=====================================================================================
   // Row-major matrix tests with vector elements
   //=====================================================================================

   {
      test_ = "Row-major MatrixAccessProxy::append( size_t, ElementType )";

      SVM mat( 2UL, 2UL, 1UL );
      mat(1,1) = SV( 3UL );
      mat(1,1).reserve( 1UL );
      mat(1,1).append( 1UL, 5 );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 0UL );
      checkSize    ( mat(0,1), 0UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 3UL );
      checkCapacity( mat(1,1), 1UL );
      checkNonZeros( mat(1,1), 1UL );

      if( mat(1,1)[1] != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Inserting an element failed\n"
             << " Details:\n"
             << "   Result:\n" << mat(1,1) << "\n"
             << "   Expected result:\n( 0 5 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major matrix tests with matrix elements
   //=====================================================================================

   {
      test_ = "Row-major MatrixAccessProxy::append( size_t, size_t, ElementType )";

      SMM mat( 2UL, 2UL, 1UL );
      mat(1,1) = SM( 2UL, 2UL );
      mat(1,1).reserve( 0UL, 1UL );
      mat(1,1).append( 0UL, 1UL, 5 );
      mat(1,1).finalize( 0UL );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0), 0UL );
      checkColumns ( mat(0,0), 0UL );
      checkRows    ( mat(0,1), 0UL );
      checkColumns ( mat(0,1), 0UL );
      checkRows    ( mat(1,0), 0UL );
      checkColumns ( mat(1,0), 0UL );
      checkRows    ( mat(1,1), 2UL );
      checkColumns ( mat(1,1), 2UL );
      checkCapacity( mat(1,1), 1UL );
      checkNonZeros( mat(1,1), 1UL );

      if( mat(1,1)(0,1) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Inserting an element failed\n"
             << " Details:\n"
             << "   Result:\n" << mat(1,1) << "\n"
             << "   Expected result:\n( 0 5 )\n( 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests with vector elements
   //=====================================================================================

   {
      test_ = "Column-major MatrixAccessProxy::append( size_t, ElementType )";

      OSVM mat( 2UL, 2UL, 1UL );
      mat(1,1) = SV( 3UL );
      mat(1,1).reserve( 1UL );
      mat(1,1).append( 1UL, 5 );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 0UL );
      checkSize    ( mat(0,1), 0UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 3UL );
      checkCapacity( mat(1,1), 1UL );
      checkNonZeros( mat(1,1), 1UL );

      if( mat(1,1)[1] != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Inserting an element failed\n"
             << " Details:\n"
             << "   Result:\n" << mat(1,1) << "\n"
             << "   Expected result:\n( 0 5 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests with matrix elements
   //=====================================================================================

   {
      test_ = "Column-major MatrixAccessProxy::append( size_t, size_t, ElementType )";

      OSMM mat( 2UL, 2UL, 1UL );
      mat(1,1) = SM( 2UL, 2UL );
      mat(1,1).reserve( 0UL, 1UL );
      mat(1,1).append( 0UL, 1UL, 5 );
      mat(1,1).finalize( 0UL );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0), 0UL );
      checkColumns ( mat(0,0), 0UL );
      checkRows    ( mat(0,1), 0UL );
      checkColumns ( mat(0,1), 0UL );
      checkRows    ( mat(1,0), 0UL );
      checkColumns ( mat(1,0), 0UL );
      checkRows    ( mat(1,1), 2UL );
      checkColumns ( mat(1,1), 2UL );
      checkCapacity( mat(1,1), 1UL );
      checkNonZeros( mat(1,1), 1UL );

      if( mat(1,1)(0,1) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Inserting an element failed\n"
             << " Details:\n"
             << "   Result:\n" << mat(1,1) << "\n"
             << "   Expected result:\n( 0 5 )\n( 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c erase() member functions of the MatrixAccessProxy class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c erase() member functions of the MatrixAccessProxy
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ProxyTest::testErase()
{
   //=====================================================================================
   // Row-major matrix tests with vector elements
   //=====================================================================================

   {
      test_ = "Row-major MatrixAccessProxy::erase( size_t )";

      SVM mat( 2UL, 2UL, 1UL );
      mat(0,0) = SV( 3UL, 1UL );
      mat(0,0).insert( 1UL, 5 );
      mat(0,0).erase( 1UL );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 3UL );
      checkNonZeros( mat(0,0), 0UL );
      checkSize    ( mat(0,1), 0UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 0UL );
   }

   {
      test_ = "Row-major MatrixAccessProxy::erase( Iterator )";

      SVM mat( 2UL, 2UL, 1UL );
      mat(0,0) = SV( 3UL, 1UL );
      mat(0,0).erase( mat(0,0).insert( 1UL, 5 ) );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 3UL );
      checkNonZeros( mat(0,0), 0UL );
      checkSize    ( mat(0,1), 0UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 0UL );
   }

   {
      test_ = "Row-major MatrixAccessProxy::erase( Iterator, Iterator )";

      SVM mat( 2UL, 2UL, 1UL );
      mat(0,0) = SV( 3UL, 1UL );
      mat(0,0).insert( 0UL, 1 );
      mat(0,0).insert( 1UL, 2 );
      mat(0,0).insert( 2UL, 3 );
      mat(0,0).erase( begin( mat(0,0) ), end( mat(0,0) ) );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 3UL );
      checkNonZeros( mat(0,0), 0UL );
      checkSize    ( mat(0,1), 0UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 0UL );
   }

   {
      test_ = "Row-major MatrixAccessProxy::erase( Predicate )";

      SVM mat( 2UL, 2UL, 1UL );
      mat(0,0) = SV( 3UL, 1UL );
      mat(0,0).insert( 1UL, 5 );
      mat(0,0).erase( []( int value ){ return value == 5; } );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 3UL );
      checkNonZeros( mat(0,0), 0UL );
      checkSize    ( mat(0,1), 0UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 0UL );
   }

   {
      test_ = "Row-major MatrixAccessProxy::erase( Iterator, Iterator, Predicate )";

      SVM mat( 2UL, 2UL, 1UL );
      mat(0,0) = SV( 3UL, 1UL );
      mat(0,0).insert( 1UL, 5 );
      mat(0,0).erase( begin( mat(0,0) ), end( mat(0,0) ), []( int value ){ return value == 5; } );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 3UL );
      checkNonZeros( mat(0,0), 0UL );
      checkSize    ( mat(0,1), 0UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 0UL );
   }


   //=====================================================================================
   // Row-major matrix tests with matrix elements
   //=====================================================================================

   {
      test_ = "Row-major MatrixAccessProxy::erase( size_t, size_t )";

      SMM mat( 2UL, 2UL, 1UL );
      mat(0,0) = SM( 2UL, 2UL, 1UL );
      mat(0,0).insert( 0UL, 1UL, 5 );
      mat(0,0).erase( 0UL, 1UL );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0), 2UL );
      checkColumns ( mat(0,0), 2UL );
      checkNonZeros( mat(0,0), 0UL );
      checkRows    ( mat(0,1), 0UL );
      checkColumns ( mat(0,1), 0UL );
      checkRows    ( mat(1,0), 0UL );
      checkColumns ( mat(1,0), 0UL );
      checkRows    ( mat(1,1), 0UL );
      checkColumns ( mat(1,1), 0UL );
   }

   {
      test_ = "Row-major MatrixAccessProxy::erase( Iterator )";

      SMM mat( 2UL, 2UL, 1UL );
      mat(0,0) = SM( 2UL, 2UL, 1UL );
      mat(0,0).erase( 0UL, mat(0,0).insert( 0UL, 1UL, 5 ) );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0), 2UL );
      checkColumns ( mat(0,0), 2UL );
      checkNonZeros( mat(0,0), 0UL );
      checkRows    ( mat(0,1), 0UL );
      checkColumns ( mat(0,1), 0UL );
      checkRows    ( mat(1,0), 0UL );
      checkColumns ( mat(1,0), 0UL );
      checkRows    ( mat(1,1), 0UL );
      checkColumns ( mat(1,1), 0UL );
   }

   {
      test_ = "Row-major MatrixAccessProxy::erase( Iterator )";

      SMM mat( 2UL, 2UL, 1UL );
      mat(0,0) = SM( 2UL, 2UL, 1UL );
      mat(0,0).insert( 0UL, 0UL, 1 );
      mat(0,0).insert( 0UL, 1UL, 2 );
      mat(0,0).erase( 0UL, begin( mat(0,0), 0UL ), end( mat(0,0), 0UL ) );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0), 2UL );
      checkColumns ( mat(0,0), 2UL );
      checkNonZeros( mat(0,0), 0UL );
      checkRows    ( mat(0,1), 0UL );
      checkColumns ( mat(0,1), 0UL );
      checkRows    ( mat(1,0), 0UL );
      checkColumns ( mat(1,0), 0UL );
      checkRows    ( mat(1,1), 0UL );
      checkColumns ( mat(1,1), 0UL );
   }

   {
      test_ = "Row-major MatrixAccessProxy::erase( Predicate )";

      SMM mat( 2UL, 2UL, 1UL );
      mat(0,0) = SM( 2UL, 2UL, 1UL );
      mat(0,0).insert( 0UL, 1UL, 5 );
      mat(0,0).erase( []( int value ){ return value == 5; } );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0), 2UL );
      checkColumns ( mat(0,0), 2UL );
      checkNonZeros( mat(0,0), 0UL );
      checkRows    ( mat(0,1), 0UL );
      checkColumns ( mat(0,1), 0UL );
      checkRows    ( mat(1,0), 0UL );
      checkColumns ( mat(1,0), 0UL );
      checkRows    ( mat(1,1), 0UL );
      checkColumns ( mat(1,1), 0UL );
   }

   {
      test_ = "Row-major MatrixAccessProxy::erase( size_t, Iterator, Iterator, Predicate )";

      SMM mat( 2UL, 2UL, 1UL );
      mat(0,0) = SM( 2UL, 2UL, 1UL );
      mat(0,0).insert( 0UL, 1UL, 5 );
      mat(0,0).erase( 0UL, begin( mat(0,0), 0UL ), end( mat(0,0), 0UL ),
                      []( int value ){ return value == 5; } );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0), 2UL );
      checkColumns ( mat(0,0), 2UL );
      checkNonZeros( mat(0,0), 0UL );
      checkRows    ( mat(0,1), 0UL );
      checkColumns ( mat(0,1), 0UL );
      checkRows    ( mat(1,0), 0UL );
      checkColumns ( mat(1,0), 0UL );
      checkRows    ( mat(1,1), 0UL );
      checkColumns ( mat(1,1), 0UL );
   }


   //=====================================================================================
   // Column-major matrix tests with vector elements
   //=====================================================================================

   {
      test_ = "Column-major MatrixAccessProxy::erase( size_t )";

      OSVM mat( 2UL, 2UL, 1UL );
      mat(0,0) = SV( 3UL, 1UL );
      mat(0,0).insert( 1UL, 5 );
      mat(0,0).erase( 1UL );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 3UL );
      checkNonZeros( mat(0,0), 0UL );
      checkSize    ( mat(0,1), 0UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 0UL );
   }

   {
      test_ = "Column-major MatrixAccessProxy::erase( Iterator )";

      OSVM mat( 2UL, 2UL, 1UL );
      mat(0,0) = SV( 3UL, 1UL );
      mat(0,0).erase( mat(0,0).insert( 1UL, 5 ) );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 3UL );
      checkNonZeros( mat(0,0), 0UL );
      checkSize    ( mat(0,1), 0UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 0UL );
   }

   {
      test_ = "Column-major MatrixAccessProxy::erase( Iterator, Iterator )";

      OSVM mat( 2UL, 2UL, 1UL );
      mat(0,0) = SV( 3UL, 1UL );
      mat(0,0).insert( 0UL, 1 );
      mat(0,0).insert( 1UL, 2 );
      mat(0,0).insert( 2UL, 3 );
      mat(0,0).erase( begin( mat(0,0) ), end( mat(0,0) ) );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 3UL );
      checkNonZeros( mat(0,0), 0UL );
      checkSize    ( mat(0,1), 0UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 0UL );
   }

   {
      test_ = "Column-major MatrixAccessProxy::erase( Predicate )";

      OSVM mat( 2UL, 2UL, 1UL );
      mat(0,0) = SV( 3UL, 1UL );
      mat(0,0).insert( 1UL, 5 );
      mat(0,0).erase( []( int value ){ return value == 5; } );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 3UL );
      checkNonZeros( mat(0,0), 0UL );
      checkSize    ( mat(0,1), 0UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 0UL );
   }

   {
      test_ = "Column-major MatrixAccessProxy::erase( Iterator, Iterator, Predicate )";

      OSVM mat( 2UL, 2UL, 1UL );
      mat(0,0) = SV( 3UL, 1UL );
      mat(0,0).insert( 1UL, 5 );
      mat(0,0).erase( begin( mat(0,0) ), end( mat(0,0) ), []( int value ){ return value == 5; } );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 3UL );
      checkNonZeros( mat(0,0), 0UL );
      checkSize    ( mat(0,1), 0UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 0UL );
   }


   //=====================================================================================
   // Column-major matrix tests with matrix elements
   //=====================================================================================

   {
      test_ = "Column-major MatrixAccessProxy::erase( size_t, size_t )";

      OSMM mat( 2UL, 2UL, 1UL );
      mat(0,0) = SM( 2UL, 2UL, 1UL );
      mat(0,0).insert( 0UL, 1UL, 5 );
      mat(0,0).erase( 0UL, 1UL );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0), 2UL );
      checkColumns ( mat(0,0), 2UL );
      checkNonZeros( mat(0,0), 0UL );
      checkRows    ( mat(0,1), 0UL );
      checkColumns ( mat(0,1), 0UL );
      checkRows    ( mat(1,0), 0UL );
      checkColumns ( mat(1,0), 0UL );
      checkRows    ( mat(1,1), 0UL );
      checkColumns ( mat(1,1), 0UL );
   }

   {
      test_ = "Column-major MatrixAccessProxy::erase( Iterator )";

      OSMM mat( 2UL, 2UL, 1UL );
      mat(0,0) = SM( 2UL, 2UL, 1UL );
      mat(0,0).erase( 0UL, mat(0,0).insert( 0UL, 1UL, 5 ) );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0), 2UL );
      checkColumns ( mat(0,0), 2UL );
      checkNonZeros( mat(0,0), 0UL );
      checkRows    ( mat(0,1), 0UL );
      checkColumns ( mat(0,1), 0UL );
      checkRows    ( mat(1,0), 0UL );
      checkColumns ( mat(1,0), 0UL );
      checkRows    ( mat(1,1), 0UL );
      checkColumns ( mat(1,1), 0UL );
   }

   {
      test_ = "Column-major MatrixAccessProxy::erase( Iterator )";

      OSMM mat( 2UL, 2UL, 1UL );
      mat(0,0) = SM( 2UL, 2UL, 1UL );
      mat(0,0).insert( 0UL, 0UL, 1 );
      mat(0,0).insert( 0UL, 1UL, 2 );
      mat(0,0).erase( 0UL, begin( mat(0,0), 0UL ), end( mat(0,0), 0UL ) );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0), 2UL );
      checkColumns ( mat(0,0), 2UL );
      checkNonZeros( mat(0,0), 0UL );
      checkRows    ( mat(0,1), 0UL );
      checkColumns ( mat(0,1), 0UL );
      checkRows    ( mat(1,0), 0UL );
      checkColumns ( mat(1,0), 0UL );
      checkRows    ( mat(1,1), 0UL );
      checkColumns ( mat(1,1), 0UL );
   }

   {
      test_ = "Column-major MatrixAccessProxy::erase( Predicate )";

      OSMM mat( 2UL, 2UL, 1UL );
      mat(0,0) = SM( 2UL, 2UL, 1UL );
      mat(0,0).insert( 0UL, 1UL, 5 );
      mat(0,0).erase( []( int value ){ return value == 5; } );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0), 2UL );
      checkColumns ( mat(0,0), 2UL );
      checkNonZeros( mat(0,0), 0UL );
      checkRows    ( mat(0,1), 0UL );
      checkColumns ( mat(0,1), 0UL );
      checkRows    ( mat(1,0), 0UL );
      checkColumns ( mat(1,0), 0UL );
      checkRows    ( mat(1,1), 0UL );
      checkColumns ( mat(1,1), 0UL );
   }

   {
      test_ = "Column-major MatrixAccessProxy::erase( size_t, Iterator, Iterator, Predicate )";

      OSMM mat( 2UL, 2UL, 1UL );
      mat(0,0) = SM( 2UL, 2UL, 1UL );
      mat(0,0).insert( 0UL, 1UL, 5 );
      mat(0,0).erase( 0UL, begin( mat(0,0), 0UL ), end( mat(0,0), 0UL ),
                      []( int value ){ return value == 5; } );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0), 2UL );
      checkColumns ( mat(0,0), 2UL );
      checkNonZeros( mat(0,0), 0UL );
      checkRows    ( mat(0,1), 0UL );
      checkColumns ( mat(0,1), 0UL );
      checkRows    ( mat(1,0), 0UL );
      checkColumns ( mat(1,0), 0UL );
      checkRows    ( mat(1,1), 0UL );
      checkColumns ( mat(1,1), 0UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c find() member functions of the MatrixAccessProxy class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c find() member functions of the MatrixAccessProxy
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ProxyTest::testFind()
{
   //=====================================================================================
   // Row-major matrix tests with vector elements
   //=====================================================================================

   {
      test_ = "Row-major MatrixAccessProxy::find( size_t )";

      SVM mat( 2UL, 2UL, 4UL );
      mat(0,0) = SV( 5UL, 3UL );
      mat(0,0)[1] = 2;
      mat(0,0)[2] = 3;
      mat(0,0)[4] = 5;

      SV::Iterator pos( mat(0,0).find( 2UL ) );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 4UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 5UL );
      checkCapacity( mat(0,0), 3UL );
      checkNonZeros( mat(0,0), 3UL );
      checkSize    ( mat(0,1), 0UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 0UL );

      if( pos == mat(0,0).end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Element could not be found\n"
             << " Details:\n"
             << "   Required index = 2\n"
             << "   Current vector:\n" << mat(0,0) << "\n";
         throw std::runtime_error( oss.str() );
      }
      else if( pos->index() != 2 || pos->value() != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Wrong element found\n"
             << " Details:\n"
             << "   Required index = 2\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = 3\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current vector:\n" << mat(0,0) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major find( MatrixAccessProxy, size_t )";

      SVM mat( 2UL, 2UL, 4UL );
      mat(0,0) = SV( 5UL, 3UL );
      mat(0,0)[1] = 2;
      mat(0,0)[2] = 3;
      mat(0,0)[4] = 5;

      SV::Iterator pos( find( mat(0,0), 2UL ) );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 4UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 5UL );
      checkCapacity( mat(0,0), 3UL );
      checkNonZeros( mat(0,0), 3UL );
      checkSize    ( mat(0,1), 0UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 0UL );

      if( pos == mat(0,0).end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Element could not be found\n"
             << " Details:\n"
             << "   Required index = 2\n"
             << "   Current vector:\n" << mat(0,0) << "\n";
         throw std::runtime_error( oss.str() );
      }
      else if( pos->index() != 2 || pos->value() != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Wrong element found\n"
             << " Details:\n"
             << "   Required index = 2\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = 3\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current vector:\n" << mat(0,0) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major matrix tests with matrix elements
   //=====================================================================================

   {
      test_ = "Row-major MatrixAccessProxy::find( size_t, size_t )";

      SMM mat( 2UL, 2UL, 4UL );
      mat(0,0) = SM( 2UL, 5UL, 3UL );
      mat(0,0)(1,1) = 2;
      mat(0,0)(1,2) = 3;
      mat(0,0)(1,4) = 5;

      SM::Iterator pos( mat(0,0).find( 1UL, 2 ) );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 4UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0), 2UL );
      checkColumns ( mat(0,0), 5UL );
      checkCapacity( mat(0,0), 3UL );
      checkNonZeros( mat(0,0), 3UL );
      checkRows    ( mat(0,1), 0UL );
      checkColumns ( mat(0,1), 0UL );
      checkRows    ( mat(1,0), 0UL );
      checkColumns ( mat(1,0), 0UL );
      checkRows    ( mat(1,1), 0UL );
      checkColumns ( mat(1,1), 0UL );

      if( pos == mat(0,0).end( 1UL ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Element could not be found\n"
             << " Details:\n"
             << "   Required index = 2\n"
             << "   Current matrix:\n" << mat(0,0) << "\n";
         throw std::runtime_error( oss.str() );
      }
      else if( pos->index() != 2 || pos->value() != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Wrong element found\n"
             << " Details:\n"
             << "   Required index = 2\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = 3\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current matrix:\n" << mat(0,0) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major find( MatrixAccessProxy, size_t, size_t )";

      SMM mat( 2UL, 2UL, 4UL );
      mat(0,0) = SM( 2UL, 5UL, 3UL );
      mat(0,0)(1,1) = 2;
      mat(0,0)(1,2) = 3;
      mat(0,0)(1,4) = 5;

      SM::Iterator pos( find( mat(0,0), 1UL, 2UL ) );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 4UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0), 2UL );
      checkColumns ( mat(0,0), 5UL );
      checkCapacity( mat(0,0), 3UL );
      checkNonZeros( mat(0,0), 3UL );
      checkRows    ( mat(0,1), 0UL );
      checkColumns ( mat(0,1), 0UL );
      checkRows    ( mat(1,0), 0UL );
      checkColumns ( mat(1,0), 0UL );
      checkRows    ( mat(1,1), 0UL );
      checkColumns ( mat(1,1), 0UL );

      if( pos == mat(0,0).end( 1UL ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Element could not be found\n"
             << " Details:\n"
             << "   Required index = 2\n"
             << "   Current matrix:\n" << mat(0,0) << "\n";
         throw std::runtime_error( oss.str() );
      }
      else if( pos->index() != 2 || pos->value() != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Wrong element found\n"
             << " Details:\n"
             << "   Required index = 2\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = 3\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current matrix:\n" << mat(0,0) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests with vector elements
   //=====================================================================================

   {
      test_ = "Column-major MatrixAccessProxy::find( size_t )";

      OSVM mat( 2UL, 2UL, 4UL );
      mat(0,0) = SV( 5UL, 3UL );
      mat(0,0)[1] = 2;
      mat(0,0)[2] = 3;
      mat(0,0)[4] = 5;

      SV::Iterator pos( mat(0,0).find( 2UL ) );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 4UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 5UL );
      checkCapacity( mat(0,0), 3UL );
      checkNonZeros( mat(0,0), 3UL );
      checkSize    ( mat(0,1), 0UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 0UL );

      if( pos == mat(0,0).end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Element could not be found\n"
             << " Details:\n"
             << "   Required index = 2\n"
             << "   Current vector:\n" << mat(0,0) << "\n";
         throw std::runtime_error( oss.str() );
      }
      else if( pos->index() != 2 || pos->value() != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Wrong element found\n"
             << " Details:\n"
             << "   Required index = 2\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = 3\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current vector:\n" << mat(0,0) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major find( MatrixAccessProxy, size_t )";

      OSVM mat( 2UL, 2UL, 4UL );
      mat(0,0) = SV( 5UL, 3UL );
      mat(0,0)[1] = 2;
      mat(0,0)[2] = 3;
      mat(0,0)[4] = 5;

      SV::Iterator pos( find( mat(0,0), 2UL ) );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 4UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 5UL );
      checkCapacity( mat(0,0), 3UL );
      checkNonZeros( mat(0,0), 3UL );
      checkSize    ( mat(0,1), 0UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 0UL );

      if( pos == mat(0,0).end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Element could not be found\n"
             << " Details:\n"
             << "   Required index = 2\n"
             << "   Current vector:\n" << mat(0,0) << "\n";
         throw std::runtime_error( oss.str() );
      }
      else if( pos->index() != 2 || pos->value() != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Wrong element found\n"
             << " Details:\n"
             << "   Required index = 2\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = 3\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current vector:\n" << mat(0,0) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests with matrix elements
   //=====================================================================================

   {
      test_ = "Column-major MatrixAccessProxy::find( size_t, size_t )";

      OSMM mat( 2UL, 2UL, 4UL );
      mat(0,0) = SM( 2UL, 5UL, 3UL );
      mat(0,0)(1,1) = 2;
      mat(0,0)(1,2) = 3;
      mat(0,0)(1,4) = 5;

      SM::Iterator pos( mat(0,0).find( 1UL, 2 ) );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 4UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0), 2UL );
      checkColumns ( mat(0,0), 5UL );
      checkCapacity( mat(0,0), 3UL );
      checkNonZeros( mat(0,0), 3UL );
      checkRows    ( mat(0,1), 0UL );
      checkColumns ( mat(0,1), 0UL );
      checkRows    ( mat(1,0), 0UL );
      checkColumns ( mat(1,0), 0UL );
      checkRows    ( mat(1,1), 0UL );
      checkColumns ( mat(1,1), 0UL );

      if( pos == mat(0,0).end( 1UL ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Element could not be found\n"
             << " Details:\n"
             << "   Required index = 2\n"
             << "   Current matrix:\n" << mat(0,0) << "\n";
         throw std::runtime_error( oss.str() );
      }
      else if( pos->index() != 2 || pos->value() != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Wrong element found\n"
             << " Details:\n"
             << "   Required index = 2\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = 3\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current matrix:\n" << mat(0,0) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major find( MatrixAccessProxy, size_t, size_t )";

      OSMM mat( 2UL, 2UL, 4UL );
      mat(0,0) = SM( 2UL, 5UL, 3UL );
      mat(0,0)(1,1) = 2;
      mat(0,0)(1,2) = 3;
      mat(0,0)(1,4) = 5;

      SM::Iterator pos( find( mat(0,0), 1UL, 2UL ) );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 4UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0), 2UL );
      checkColumns ( mat(0,0), 5UL );
      checkCapacity( mat(0,0), 3UL );
      checkNonZeros( mat(0,0), 3UL );
      checkRows    ( mat(0,1), 0UL );
      checkColumns ( mat(0,1), 0UL );
      checkRows    ( mat(1,0), 0UL );
      checkColumns ( mat(1,0), 0UL );
      checkRows    ( mat(1,1), 0UL );
      checkColumns ( mat(1,1), 0UL );

      if( pos == mat(0,0).end( 1UL ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Element could not be found\n"
             << " Details:\n"
             << "   Required index = 2\n"
             << "   Current matrix:\n" << mat(0,0) << "\n";
         throw std::runtime_error( oss.str() );
      }
      else if( pos->index() != 2 || pos->value() != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Wrong element found\n"
             << " Details:\n"
             << "   Required index = 2\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = 3\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current matrix:\n" << mat(0,0) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c lowerBound() member functions of the MatrixAccessProxy class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c lowerBound() member functions of the MatrixAccessProxy
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ProxyTest::testLowerBound()
{
   //=====================================================================================
   // Row-major matrix tests with vector elements
   //=====================================================================================

   {
      test_ = "Row-major MatrixAccessProxy::lowerBound( size_t )";

      SVM mat( 2UL, 2UL, 4UL );
      mat(0,0) = SV( 5UL, 3UL );
      mat(0,0)[1] = 2;
      mat(0,0)[2] = 3;
      mat(0,0)[4] = 5;

      SV::Iterator pos( mat(0,0).lowerBound( 3UL ) );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 4UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 5UL );
      checkCapacity( mat(0,0), 3UL );
      checkNonZeros( mat(0,0), 3UL );
      checkSize    ( mat(0,1), 0UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 0UL );

      if( pos == mat(0,0).end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Element could not be found\n"
             << " Details:\n"
             << "   Required index = 3\n"
             << "   Current vector:\n" << mat(0,0) << "\n";
         throw std::runtime_error( oss.str() );
      }
      else if( pos->index() != 4 || pos->value() != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Wrong element found\n"
             << " Details:\n"
             << "   Required index = 4\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = 5\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current vector:\n" << mat(0,0) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major lowerBound( MatrixAccessProxy, size_t )";

      SVM mat( 2UL, 2UL, 4UL );
      mat(0,0) = SV( 5UL, 3UL );
      mat(0,0)[1] = 2;
      mat(0,0)[2] = 3;
      mat(0,0)[4] = 5;

      SV::Iterator pos( lowerBound( mat(0,0), 3UL ) );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 4UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 5UL );
      checkCapacity( mat(0,0), 3UL );
      checkNonZeros( mat(0,0), 3UL );
      checkSize    ( mat(0,1), 0UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 0UL );

      if( pos == mat(0,0).end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Element could not be found\n"
             << " Details:\n"
             << "   Required index = 3\n"
             << "   Current vector:\n" << mat(0,0) << "\n";
         throw std::runtime_error( oss.str() );
      }
      else if( pos->index() != 4 || pos->value() != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Wrong element found\n"
             << " Details:\n"
             << "   Required index = 4\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = 5\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current vector:\n" << mat(0,0) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major matrix tests with matrix elements
   //=====================================================================================

   {
      test_ = "Row-major MatrixAccessProxy::lowerBound( size_t, size_t )";

      SMM mat( 2UL, 2UL, 4UL );
      mat(0,0) = SM( 2UL, 5UL, 3UL );
      mat(0,0)(1,1) = 2;
      mat(0,0)(1,2) = 3;
      mat(0,0)(1,4) = 5;

      SM::Iterator pos( mat(0,0).lowerBound( 1UL, 3 ) );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 4UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0), 2UL );
      checkColumns ( mat(0,0), 5UL );
      checkCapacity( mat(0,0), 3UL );
      checkNonZeros( mat(0,0), 3UL );
      checkRows    ( mat(0,1), 0UL );
      checkColumns ( mat(0,1), 0UL );
      checkRows    ( mat(1,0), 0UL );
      checkColumns ( mat(1,0), 0UL );
      checkRows    ( mat(1,1), 0UL );
      checkColumns ( mat(1,1), 0UL );

      if( pos == mat(0,0).end( 1UL ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Element could not be found\n"
             << " Details:\n"
             << "   Required index = 3\n"
             << "   Current matrix:\n" << mat(0,0) << "\n";
         throw std::runtime_error( oss.str() );
      }
      else if( pos->index() != 4 || pos->value() != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Wrong element found\n"
             << " Details:\n"
             << "   Required index = 4\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = 5\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current matrix:\n" << mat(0,0) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major lowerBound( MatrixAccessProxy, size_t, size_t )";

      SMM mat( 2UL, 2UL, 4UL );
      mat(0,0) = SM( 2UL, 5UL, 3UL );
      mat(0,0)(1,1) = 2;
      mat(0,0)(1,2) = 3;
      mat(0,0)(1,4) = 5;

      SM::Iterator pos( lowerBound( mat(0,0), 1UL, 3UL ) );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 4UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0), 2UL );
      checkColumns ( mat(0,0), 5UL );
      checkCapacity( mat(0,0), 3UL );
      checkNonZeros( mat(0,0), 3UL );
      checkRows    ( mat(0,1), 0UL );
      checkColumns ( mat(0,1), 0UL );
      checkRows    ( mat(1,0), 0UL );
      checkColumns ( mat(1,0), 0UL );
      checkRows    ( mat(1,1), 0UL );
      checkColumns ( mat(1,1), 0UL );

      if( pos == mat(0,0).end( 1UL ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Element could not be found\n"
             << " Details:\n"
             << "   Required index = 3\n"
             << "   Current matrix:\n" << mat(0,0) << "\n";
         throw std::runtime_error( oss.str() );
      }
      else if( pos->index() != 4 || pos->value() != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Wrong element found\n"
             << " Details:\n"
             << "   Required index = 4\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = 5\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current matrix:\n" << mat(0,0) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests with vector elements
   //=====================================================================================

   {
      test_ = "Column-major MatrixAccessProxy::lowerBound( size_t )";

      SVM mat( 2UL, 2UL, 4UL );
      mat(0,0) = SV( 5UL, 3UL );
      mat(0,0)[1] = 2;
      mat(0,0)[2] = 3;
      mat(0,0)[4] = 5;

      SV::Iterator pos( mat(0,0).lowerBound( 3UL ) );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 4UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 5UL );
      checkCapacity( mat(0,0), 3UL );
      checkNonZeros( mat(0,0), 3UL );
      checkSize    ( mat(0,1), 0UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 0UL );

      if( pos == mat(0,0).end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Element could not be found\n"
             << " Details:\n"
             << "   Required index = 3\n"
             << "   Current vector:\n" << mat(0,0) << "\n";
         throw std::runtime_error( oss.str() );
      }
      else if( pos->index() != 4 || pos->value() != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Wrong element found\n"
             << " Details:\n"
             << "   Required index = 4\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = 5\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current vector:\n" << mat(0,0) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major lowerBound( MatrixAccessProxy, size_t )";

      SVM mat( 2UL, 2UL, 4UL );
      mat(0,0) = SV( 5UL, 3UL );
      mat(0,0)[1] = 2;
      mat(0,0)[2] = 3;
      mat(0,0)[4] = 5;

      SV::Iterator pos( lowerBound( mat(0,0), 3UL ) );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 4UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 5UL );
      checkCapacity( mat(0,0), 3UL );
      checkNonZeros( mat(0,0), 3UL );
      checkSize    ( mat(0,1), 0UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 0UL );

      if( pos == mat(0,0).end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Element could not be found\n"
             << " Details:\n"
             << "   Required index = 3\n"
             << "   Current vector:\n" << mat(0,0) << "\n";
         throw std::runtime_error( oss.str() );
      }
      else if( pos->index() != 4 || pos->value() != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Wrong element found\n"
             << " Details:\n"
             << "   Required index = 4\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = 5\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current vector:\n" << mat(0,0) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests with matrix elements
   //=====================================================================================

   {
      test_ = "Column-major MatrixAccessProxy::lowerBound( size_t, size_t )";

      SMM mat( 2UL, 2UL, 4UL );
      mat(0,0) = SM( 2UL, 5UL, 3UL );
      mat(0,0)(1,1) = 2;
      mat(0,0)(1,2) = 3;
      mat(0,0)(1,4) = 5;

      SM::Iterator pos( mat(0,0).lowerBound( 1UL, 3 ) );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 4UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0), 2UL );
      checkColumns ( mat(0,0), 5UL );
      checkCapacity( mat(0,0), 3UL );
      checkNonZeros( mat(0,0), 3UL );
      checkRows    ( mat(0,1), 0UL );
      checkColumns ( mat(0,1), 0UL );
      checkRows    ( mat(1,0), 0UL );
      checkColumns ( mat(1,0), 0UL );
      checkRows    ( mat(1,1), 0UL );
      checkColumns ( mat(1,1), 0UL );

      if( pos == mat(0,0).end( 1UL ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Element could not be found\n"
             << " Details:\n"
             << "   Required index = 3\n"
             << "   Current matrix:\n" << mat(0,0) << "\n";
         throw std::runtime_error( oss.str() );
      }
      else if( pos->index() != 4 || pos->value() != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Wrong element found\n"
             << " Details:\n"
             << "   Required index = 4\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = 5\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current matrix:\n" << mat(0,0) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major lowerBound( MatrixAccessProxy, size_t, size_t )";

      SMM mat( 2UL, 2UL, 4UL );
      mat(0,0) = SM( 2UL, 5UL, 3UL );
      mat(0,0)(1,1) = 2;
      mat(0,0)(1,2) = 3;
      mat(0,0)(1,4) = 5;

      SM::Iterator pos( lowerBound( mat(0,0), 1UL, 3UL ) );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 4UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0), 2UL );
      checkColumns ( mat(0,0), 5UL );
      checkCapacity( mat(0,0), 3UL );
      checkNonZeros( mat(0,0), 3UL );
      checkRows    ( mat(0,1), 0UL );
      checkColumns ( mat(0,1), 0UL );
      checkRows    ( mat(1,0), 0UL );
      checkColumns ( mat(1,0), 0UL );
      checkRows    ( mat(1,1), 0UL );
      checkColumns ( mat(1,1), 0UL );

      if( pos == mat(0,0).end( 1UL ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Element could not be found\n"
             << " Details:\n"
             << "   Required index = 3\n"
             << "   Current matrix:\n" << mat(0,0) << "\n";
         throw std::runtime_error( oss.str() );
      }
      else if( pos->index() != 4 || pos->value() != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Wrong element found\n"
             << " Details:\n"
             << "   Required index = 4\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = 5\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current matrix:\n" << mat(0,0) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c upperBound() member functions of the MatrixAccessProxy class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c upperBound() member functions of the MatrixAccessProxy
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ProxyTest::testUpperBound()
{
   //=====================================================================================
   // Row-major matrix tests with vector elements
   //=====================================================================================

   {
      test_ = "Row-major MatrixAccessProxy::upperBound( size_t )";

      SVM mat( 2UL, 2UL, 4UL );
      mat(0,0) = SV( 5UL, 3UL );
      mat(0,0)[1] = 2;
      mat(0,0)[2] = 3;
      mat(0,0)[4] = 5;

      SV::Iterator pos( mat(0,0).upperBound( 3UL ) );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 4UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 5UL );
      checkCapacity( mat(0,0), 3UL );
      checkNonZeros( mat(0,0), 3UL );
      checkSize    ( mat(0,1), 0UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 0UL );

      if( pos == mat(0,0).end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Element could not be found\n"
             << " Details:\n"
             << "   Required index = 3\n"
             << "   Current vector:\n" << mat(0,0) << "\n";
         throw std::runtime_error( oss.str() );
      }
      else if( pos->index() != 4 || pos->value() != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Wrong element found\n"
             << " Details:\n"
             << "   Required index = 4\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = 5\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current vector:\n" << mat(0,0) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major upperBound( MatrixAccessProxy, size_t )";

      SVM mat( 2UL, 2UL, 4UL );
      mat(0,0) = SV( 5UL, 3UL );
      mat(0,0)[1] = 2;
      mat(0,0)[2] = 3;
      mat(0,0)[4] = 5;

      SV::Iterator pos( upperBound( mat(0,0), 3UL ) );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 4UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 5UL );
      checkCapacity( mat(0,0), 3UL );
      checkNonZeros( mat(0,0), 3UL );
      checkSize    ( mat(0,1), 0UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 0UL );

      if( pos == mat(0,0).end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Element could not be found\n"
             << " Details:\n"
             << "   Required index = 3\n"
             << "   Current vector:\n" << mat(0,0) << "\n";
         throw std::runtime_error( oss.str() );
      }
      else if( pos->index() != 4 || pos->value() != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Wrong element found\n"
             << " Details:\n"
             << "   Required index = 4\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = 5\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current vector:\n" << mat(0,0) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major matrix tests with matrix elements
   //=====================================================================================

   {
      test_ = "Row-major MatrixAccessProxy::upperBound( size_t, size_t )";

      SMM mat( 2UL, 2UL, 4UL );
      mat(0,0) = SM( 2UL, 5UL, 3UL );
      mat(0,0)(1,1) = 2;
      mat(0,0)(1,2) = 3;
      mat(0,0)(1,4) = 5;

      SM::Iterator pos( mat(0,0).upperBound( 1UL, 3 ) );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 4UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0), 2UL );
      checkColumns ( mat(0,0), 5UL );
      checkCapacity( mat(0,0), 3UL );
      checkNonZeros( mat(0,0), 3UL );
      checkRows    ( mat(0,1), 0UL );
      checkColumns ( mat(0,1), 0UL );
      checkRows    ( mat(1,0), 0UL );
      checkColumns ( mat(1,0), 0UL );
      checkRows    ( mat(1,1), 0UL );
      checkColumns ( mat(1,1), 0UL );

      if( pos == mat(0,0).end( 1UL ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Element could not be found\n"
             << " Details:\n"
             << "   Required index = 3\n"
             << "   Current matrix:\n" << mat(0,0) << "\n";
         throw std::runtime_error( oss.str() );
      }
      else if( pos->index() != 4 || pos->value() != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Wrong element found\n"
             << " Details:\n"
             << "   Required index = 4\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = 5\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current matrix:\n" << mat(0,0) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major upperBound( MatrixAccessProxy, size_t, size_t )";

      SMM mat( 2UL, 2UL, 4UL );
      mat(0,0) = SM( 2UL, 5UL, 3UL );
      mat(0,0)(1,1) = 2;
      mat(0,0)(1,2) = 3;
      mat(0,0)(1,4) = 5;

      SM::Iterator pos( upperBound( mat(0,0), 1UL, 3UL ) );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 4UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0), 2UL );
      checkColumns ( mat(0,0), 5UL );
      checkCapacity( mat(0,0), 3UL );
      checkNonZeros( mat(0,0), 3UL );
      checkRows    ( mat(0,1), 0UL );
      checkColumns ( mat(0,1), 0UL );
      checkRows    ( mat(1,0), 0UL );
      checkColumns ( mat(1,0), 0UL );
      checkRows    ( mat(1,1), 0UL );
      checkColumns ( mat(1,1), 0UL );

      if( pos == mat(0,0).end( 1UL ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Element could not be found\n"
             << " Details:\n"
             << "   Required index = 3\n"
             << "   Current matrix:\n" << mat(0,0) << "\n";
         throw std::runtime_error( oss.str() );
      }
      else if( pos->index() != 4 || pos->value() != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Wrong element found\n"
             << " Details:\n"
             << "   Required index = 4\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = 5\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current matrix:\n" << mat(0,0) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests with vector elements
   //=====================================================================================

   {
      test_ = "Column-major MatrixAccessProxy::upperBound( size_t )";

      SVM mat( 2UL, 2UL, 4UL );
      mat(0,0) = SV( 5UL, 3UL );
      mat(0,0)[1] = 2;
      mat(0,0)[2] = 3;
      mat(0,0)[4] = 5;

      SV::Iterator pos( mat(0,0).upperBound( 3UL ) );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 4UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 5UL );
      checkCapacity( mat(0,0), 3UL );
      checkNonZeros( mat(0,0), 3UL );
      checkSize    ( mat(0,1), 0UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 0UL );

      if( pos == mat(0,0).end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Element could not be found\n"
             << " Details:\n"
             << "   Required index = 3\n"
             << "   Current vector:\n" << mat(0,0) << "\n";
         throw std::runtime_error( oss.str() );
      }
      else if( pos->index() != 4 || pos->value() != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Wrong element found\n"
             << " Details:\n"
             << "   Required index = 4\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = 5\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current vector:\n" << mat(0,0) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major MatrixAccessProxy::upperBound( size_t )";

      SVM mat( 2UL, 2UL, 4UL );
      mat(0,0) = SV( 5UL, 3UL );
      mat(0,0)[1] = 2;
      mat(0,0)[2] = 3;
      mat(0,0)[4] = 5;

      SV::Iterator pos( upperBound( mat(0,0), 3UL ) );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 4UL );
      checkNonZeros( mat, 1UL );

      checkSize    ( mat(0,0), 5UL );
      checkCapacity( mat(0,0), 3UL );
      checkNonZeros( mat(0,0), 3UL );
      checkSize    ( mat(0,1), 0UL );
      checkSize    ( mat(1,0), 0UL );
      checkSize    ( mat(1,1), 0UL );

      if( pos == mat(0,0).end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Element could not be found\n"
             << " Details:\n"
             << "   Required index = 3\n"
             << "   Current vector:\n" << mat(0,0) << "\n";
         throw std::runtime_error( oss.str() );
      }
      else if( pos->index() != 4 || pos->value() != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Wrong element found\n"
             << " Details:\n"
             << "   Required index = 4\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = 5\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current vector:\n" << mat(0,0) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests with matrix elements
   //=====================================================================================

   {
      test_ = "Column-major MatrixAccessProxy::upperBound( size_t, size_t )";

      SMM mat( 2UL, 2UL, 4UL );
      mat(0,0) = SM( 2UL, 5UL, 3UL );
      mat(0,0)(1,1) = 2;
      mat(0,0)(1,2) = 3;
      mat(0,0)(1,4) = 5;

      SM::Iterator pos( mat(0,0).upperBound( 1UL, 3 ) );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 4UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0), 2UL );
      checkColumns ( mat(0,0), 5UL );
      checkCapacity( mat(0,0), 3UL );
      checkNonZeros( mat(0,0), 3UL );
      checkRows    ( mat(0,1), 0UL );
      checkColumns ( mat(0,1), 0UL );
      checkRows    ( mat(1,0), 0UL );
      checkColumns ( mat(1,0), 0UL );
      checkRows    ( mat(1,1), 0UL );
      checkColumns ( mat(1,1), 0UL );

      if( pos == mat(0,0).end( 1UL ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Element could not be found\n"
             << " Details:\n"
             << "   Required index = 3\n"
             << "   Current matrix:\n" << mat(0,0) << "\n";
         throw std::runtime_error( oss.str() );
      }
      else if( pos->index() != 4 || pos->value() != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Wrong element found\n"
             << " Details:\n"
             << "   Required index = 4\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = 5\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current matrix:\n" << mat(0,0) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major upperBound( MatrixAccessProxy, size_t, size_t )";

      SMM mat( 2UL, 2UL, 4UL );
      mat(0,0) = SM( 2UL, 5UL, 3UL );
      mat(0,0)(1,1) = 2;
      mat(0,0)(1,2) = 3;
      mat(0,0)(1,4) = 5;

      SM::Iterator pos( upperBound( mat(0,0), 1UL, 3UL ) );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 4UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0), 2UL );
      checkColumns ( mat(0,0), 5UL );
      checkCapacity( mat(0,0), 3UL );
      checkNonZeros( mat(0,0), 3UL );
      checkRows    ( mat(0,1), 0UL );
      checkColumns ( mat(0,1), 0UL );
      checkRows    ( mat(1,0), 0UL );
      checkColumns ( mat(1,0), 0UL );
      checkRows    ( mat(1,1), 0UL );
      checkColumns ( mat(1,1), 0UL );

      if( pos == mat(0,0).end( 1UL ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Element could not be found\n"
             << " Details:\n"
             << "   Required index = 3\n"
             << "   Current matrix:\n" << mat(0,0) << "\n";
         throw std::runtime_error( oss.str() );
      }
      else if( pos->index() != 4 || pos->value() != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Wrong element found\n"
             << " Details:\n"
             << "   Required index = 4\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = 5\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current matrix:\n" << mat(0,0) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c transpose() functions of the MatrixAccessProxy class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c transpose() functions of the MatrixAccessProxy class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ProxyTest::testTranspose()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major MatrixAccessProxy::transpose()";

      DMM mat( 2UL, 2UL, 1UL );
      mat(1,1).resize( 5UL, 3UL );
      mat(1,1).transpose();

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0),  0UL );
      checkColumns ( mat(0,0),  0UL );
      checkRows    ( mat(0,1),  0UL );
      checkColumns ( mat(0,1),  0UL );
      checkRows    ( mat(1,0),  0UL );
      checkColumns ( mat(1,0),  0UL );
      checkRows    ( mat(1,1),  3UL );
      checkColumns ( mat(1,1),  5UL );
      checkCapacity( mat(1,1), 15UL );
   }

   {
      test_ = "Row-major transpose( MatrixAccessProxy )";

      DMM mat( 2UL, 2UL, 1UL );
      mat(1,1).resize( 5UL, 3UL );
      transpose( mat(1,1) );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0),  0UL );
      checkColumns ( mat(0,0),  0UL );
      checkRows    ( mat(0,1),  0UL );
      checkColumns ( mat(0,1),  0UL );
      checkRows    ( mat(1,0),  0UL );
      checkColumns ( mat(1,0),  0UL );
      checkRows    ( mat(1,1),  3UL );
      checkColumns ( mat(1,1),  5UL );
      checkCapacity( mat(1,1), 15UL );
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major MatrixAccessProxy::transpose()";

      ODMM mat( 2UL, 2UL, 1UL );
      mat(1,1).resize( 5UL, 3UL );
      mat(1,1).transpose();

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0),  0UL );
      checkColumns ( mat(0,0),  0UL );
      checkRows    ( mat(0,1),  0UL );
      checkColumns ( mat(0,1),  0UL );
      checkRows    ( mat(1,0),  0UL );
      checkColumns ( mat(1,0),  0UL );
      checkRows    ( mat(1,1),  3UL );
      checkColumns ( mat(1,1),  5UL );
      checkCapacity( mat(1,1), 15UL );
   }

   {
      test_ = "Column-major transpose( MatrixAccessProxy )";

      ODMM mat( 2UL, 2UL, 1UL );
      mat(1,1).resize( 5UL, 3UL );
      transpose( mat(1,1) );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0),  0UL );
      checkColumns ( mat(0,0),  0UL );
      checkRows    ( mat(0,1),  0UL );
      checkColumns ( mat(0,1),  0UL );
      checkRows    ( mat(1,0),  0UL );
      checkColumns ( mat(1,0),  0UL );
      checkRows    ( mat(1,1),  3UL );
      checkColumns ( mat(1,1),  5UL );
      checkCapacity( mat(1,1), 15UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c ctranspose() functions of the MatrixAccessProxy class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c ctranspose() functions of the MatrixAccessProxy class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ProxyTest::testCTranspose()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major MatrixAccessProxy::ctranspose()";

      DMM mat( 2UL, 2UL, 1UL );
      mat(1,1).resize( 5UL, 3UL );
      mat(1,1).ctranspose();

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0),  0UL );
      checkColumns ( mat(0,0),  0UL );
      checkRows    ( mat(0,1),  0UL );
      checkColumns ( mat(0,1),  0UL );
      checkRows    ( mat(1,0),  0UL );
      checkColumns ( mat(1,0),  0UL );
      checkRows    ( mat(1,1),  3UL );
      checkColumns ( mat(1,1),  5UL );
      checkCapacity( mat(1,1), 15UL );
   }

   {
      test_ = "Row-major ctranspose( MatrixAccessProxy )";

      DMM mat( 2UL, 2UL, 1UL );
      mat(1,1).resize( 5UL, 3UL );
      ctranspose( mat(1,1) );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0),  0UL );
      checkColumns ( mat(0,0),  0UL );
      checkRows    ( mat(0,1),  0UL );
      checkColumns ( mat(0,1),  0UL );
      checkRows    ( mat(1,0),  0UL );
      checkColumns ( mat(1,0),  0UL );
      checkRows    ( mat(1,1),  3UL );
      checkColumns ( mat(1,1),  5UL );
      checkCapacity( mat(1,1), 15UL );
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major MatrixAccessProxy ctranspose()";

      ODMM mat( 2UL, 2UL, 1UL );
      mat(1,1).resize( 5UL, 3UL );
      mat(1,1).ctranspose();

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0),  0UL );
      checkColumns ( mat(0,0),  0UL );
      checkRows    ( mat(0,1),  0UL );
      checkColumns ( mat(0,1),  0UL );
      checkRows    ( mat(1,0),  0UL );
      checkColumns ( mat(1,0),  0UL );
      checkRows    ( mat(1,1),  3UL );
      checkColumns ( mat(1,1),  5UL );
      checkCapacity( mat(1,1), 15UL );
   }

   {
      test_ = "Column-major ctranspose( MatrixAccessProxy )";

      ODMM mat( 2UL, 2UL, 1UL );
      mat(1,1).resize( 5UL, 3UL );
      ctranspose( mat(1,1) );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0),  0UL );
      checkColumns ( mat(0,0),  0UL );
      checkRows    ( mat(0,1),  0UL );
      checkColumns ( mat(0,1),  0UL );
      checkRows    ( mat(1,0),  0UL );
      checkColumns ( mat(1,0),  0UL );
      checkRows    ( mat(1,1),  3UL );
      checkColumns ( mat(1,1),  5UL );
      checkCapacity( mat(1,1), 15UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c invert() function with the MatrixAccessProxy class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c invert() function with the MatrixAccessProxy class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ProxyTest::testInvert()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   using blaze::invert;
   using blaze::byLU;
   using blaze::byLLH;


   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major invert( MatrixAccessProxy )";

      blaze::CompressedMatrix< blaze::DynamicMatrix<double>, blaze::rowMajor > mat( 2UL, 2UL, 1UL );
      mat(1,1).resize( 3UL, 3UL );
      mat(1,1) = 0.0;
      mat(1,1)(0,0) = 1.0;
      mat(1,1)(1,1) = 1.0;
      mat(1,1)(2,2) = 1.0;
      invert( mat(1,1) );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0), 0UL );
      checkColumns ( mat(0,0), 0UL );
      checkRows    ( mat(0,1), 0UL );
      checkColumns ( mat(0,1), 0UL );
      checkRows    ( mat(1,0), 0UL );
      checkColumns ( mat(1,0), 0UL );
      checkRows    ( mat(1,1), 3UL );
      checkColumns ( mat(1,1), 3UL );
      checkCapacity( mat(1,1), 9UL );
      checkNonZeros( mat(1,1), 3UL );
   }

   {
      test_ = "Row-major invert<byLU>( MatrixAccessProxy )";

      blaze::CompressedMatrix< blaze::DynamicMatrix<double>, blaze::rowMajor > mat( 2UL, 2UL, 1UL );
      mat(1,1).resize( 3UL, 3UL );
      mat(1,1) = 0.0;
      mat(1,1)(0,0) = 1.0;
      mat(1,1)(1,1) = 1.0;
      mat(1,1)(2,2) = 1.0;
      invert<byLU>( mat(1,1) );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0), 0UL );
      checkColumns ( mat(0,0), 0UL );
      checkRows    ( mat(0,1), 0UL );
      checkColumns ( mat(0,1), 0UL );
      checkRows    ( mat(1,0), 0UL );
      checkColumns ( mat(1,0), 0UL );
      checkRows    ( mat(1,1), 3UL );
      checkColumns ( mat(1,1), 3UL );
      checkCapacity( mat(1,1), 9UL );
      checkNonZeros( mat(1,1), 3UL );
   }

   {
      test_ = "Row-major invert<byLLH>( MatrixAccessProxy )";

      blaze::CompressedMatrix< blaze::DynamicMatrix<double>, blaze::rowMajor > mat( 2UL, 2UL, 1UL );
      mat(1,1).resize( 3UL, 3UL );
      mat(1,1) = 0.0;
      mat(1,1)(0,0) = 1.0;
      mat(1,1)(1,1) = 1.0;
      mat(1,1)(2,2) = 1.0;
      invert<byLLH>( mat(1,1) );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0), 0UL );
      checkColumns ( mat(0,0), 0UL );
      checkRows    ( mat(0,1), 0UL );
      checkColumns ( mat(0,1), 0UL );
      checkRows    ( mat(1,0), 0UL );
      checkColumns ( mat(1,0), 0UL );
      checkRows    ( mat(1,1), 3UL );
      checkColumns ( mat(1,1), 3UL );
      checkCapacity( mat(1,1), 9UL );
      checkNonZeros( mat(1,1), 3UL );
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major invert( MatrixAccessProxy )";

      blaze::CompressedMatrix< blaze::DynamicMatrix<double>, blaze::columnMajor > mat( 2UL, 2UL, 1UL );
      mat(1,1).resize( 3UL, 3UL );
      mat(1,1) = 0.0;
      mat(1,1)(0,0) = 1.0;
      mat(1,1)(1,1) = 1.0;
      mat(1,1)(2,2) = 1.0;
      invert( mat(1,1) );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0), 0UL );
      checkColumns ( mat(0,0), 0UL );
      checkRows    ( mat(0,1), 0UL );
      checkColumns ( mat(0,1), 0UL );
      checkRows    ( mat(1,0), 0UL );
      checkColumns ( mat(1,0), 0UL );
      checkRows    ( mat(1,1), 3UL );
      checkColumns ( mat(1,1), 3UL );
      checkCapacity( mat(1,1), 9UL );
      checkNonZeros( mat(1,1), 3UL );
   }

   {
      test_ = "Column-major invert<byLU>( MatrixAccessProxy )";

      blaze::CompressedMatrix< blaze::DynamicMatrix<double>, blaze::columnMajor > mat( 2UL, 2UL, 1UL );
      mat(1,1).resize( 3UL, 3UL );
      mat(1,1) = 0.0;
      mat(1,1)(0,0) = 1.0;
      mat(1,1)(1,1) = 1.0;
      mat(1,1)(2,2) = 1.0;
      invert<byLU>( mat(1,1) );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0), 0UL );
      checkColumns ( mat(0,0), 0UL );
      checkRows    ( mat(0,1), 0UL );
      checkColumns ( mat(0,1), 0UL );
      checkRows    ( mat(1,0), 0UL );
      checkColumns ( mat(1,0), 0UL );
      checkRows    ( mat(1,1), 3UL );
      checkColumns ( mat(1,1), 3UL );
      checkCapacity( mat(1,1), 9UL );
      checkNonZeros( mat(1,1), 3UL );
   }

   {
      test_ = "Column-major invert<byLLH>( MatrixAccessProxy )";

      blaze::CompressedMatrix< blaze::DynamicMatrix<double>, blaze::columnMajor > mat( 2UL, 2UL, 1UL );
      mat(1,1).resize( 3UL, 3UL );
      mat(1,1) = 0.0;
      mat(1,1)(0,0) = 1.0;
      mat(1,1)(1,1) = 1.0;
      mat(1,1)(2,2) = 1.0;
      invert<byLLH>( mat(1,1) );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 1UL );
      checkNonZeros( mat, 1UL );

      checkRows    ( mat(0,0), 0UL );
      checkColumns ( mat(0,0), 0UL );
      checkRows    ( mat(0,1), 0UL );
      checkColumns ( mat(0,1), 0UL );
      checkRows    ( mat(1,0), 0UL );
      checkColumns ( mat(1,0), 0UL );
      checkRows    ( mat(1,1), 3UL );
      checkColumns ( mat(1,1), 3UL );
      checkCapacity( mat(1,1), 9UL );
      checkNonZeros( mat(1,1), 3UL );
   }

#endif
}
//*************************************************************************************************

} // namespace compressedmatrix

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
   std::cout << "   Running CompressedMatrix proxy test..." << std::endl;

   try
   {
      RUN_COMPRESSEDMATRIX_PROXY_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during CompressedMatrix proxy test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
