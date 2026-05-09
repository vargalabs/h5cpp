//=================================================================================================
/*!
//  \file src/mathtest/vectors/compressedvector/ProxyTest.cpp
//  \brief Source file for the CompressedVector proxy test
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
#include <blazetest/mathtest/vectors/compressedvector/ProxyTest.h>
#include <blazetest/system/LAPACK.h>

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
/*!\brief Constructor for the CompressedVector proxy test.
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
/*!\brief Test of the VectorAccessProxy assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all assignment operators of the VectorAccessProxy class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ProxyTest::testAssignment()
{
   //=====================================================================================
   // Homogeneous assignment
   //=====================================================================================

   {
      test_ = "VectorAccessProxy homogeneous assignment";

      DVV vec( 3UL, 1UL );
      vec[1] = DV( 3UL, 2 );

      vec[1] = 4;

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkSize    ( vec[0], 0UL );
      checkSize    ( vec[1], 3UL );
      checkCapacity( vec[1], 3UL );
      checkNonZeros( vec[1], 3UL );
      checkSize    ( vec[2], 0UL );

      if( vec[1][0] != 4 || vec[1][1] != 4 || vec[1][2] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec[1] << "\n"
             << "   Expected result:\n( 4 4 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // List assignment
   //=====================================================================================

   {
      test_ = "VectorAccessProxy 1D initializer list assignment";

      DVV vec( 3UL, 1UL );
      vec[1] = DV( 3UL, 2 );

      vec[1] = { 1, -2, 3 };

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkSize    ( vec[0], 0UL );
      checkSize    ( vec[1], 3UL );
      checkCapacity( vec[1], 3UL );
      checkNonZeros( vec[1], 3UL );
      checkSize    ( vec[2], 0UL );

      if( vec[1][0] != 1 || vec[1][1] != -2 || vec[1][2] != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec[1] << "\n"
             << "   Expected result:\n( 1 -2 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "VectorAccessProxy 2D initializer list assignment";

      using blaze::initializer_list;

      DMV vec( 3UL, 1UL );
      vec[1] = DM( 3UL, 3UL, 2 );

      initializer_list< initializer_list<int> > list = { { 1, -2, 3 }, { -2, 4, -6 }, { 3, -6, 9 } };
      vec[1] = list;

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkRows    ( vec[0], 0UL );
      checkColumns ( vec[0], 0UL );
      checkRows    ( vec[1], 3UL );
      checkColumns ( vec[1], 3UL );
      checkCapacity( vec[1], 9UL );
      checkNonZeros( vec[1], 9UL );
      checkRows    ( vec[2], 0UL );
      checkColumns ( vec[2], 0UL );

      if( vec[1](0,0) !=  1 || vec[1](0,1) != -2 || vec[1](0,2) !=  3 ||
          vec[1](1,0) != -2 || vec[1](1,1) !=  4 || vec[1](1,2) != -6 ||
          vec[1](2,0) !=  3 || vec[1](2,1) != -6 || vec[1](2,2) !=  9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec[1] << "\n"
             << "   Expected result:\n(  1 -2  3 )\n( -2  4 -6 )\n(  3 -6  9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Array assignment
   //=====================================================================================

   {
      test_ = "VectorAccessProxy array assignment";

      const int array[3] = { 1, 2, 3 };
      DVV vec( 3UL, 1UL );
      vec[1] = DV( 3UL, 0 );

      vec[1] = array;

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkSize    ( vec[0], 0UL );
      checkSize    ( vec[1], 3UL );
      checkCapacity( vec[1], 3UL );
      checkNonZeros( vec[1], 3UL );
      checkSize    ( vec[2], 0UL );

      if( vec[1][0] != 1 || vec[1][1] != 2 || vec[1][2] != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec[1] << "\n"
             << "   Expected result:\n( 1 2 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Copy assignment
   //=====================================================================================

   {
      test_ = "VectorAccessProxy copy assignment";

      DVV vec( 3UL, 1UL );
      vec[0] = DV( 3UL, 5 );

      vec[1] = vec[0];

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 2UL );

      checkSize    ( vec[0], 3UL );
      checkCapacity( vec[0], 3UL );
      checkNonZeros( vec[0], 3UL );
      checkSize    ( vec[1], 3UL );
      checkCapacity( vec[1], 3UL );
      checkNonZeros( vec[1], 3UL );
      checkSize    ( vec[2], 0UL );

      if( vec[1][0] != 5 || vec[1][1] != 5 || vec[1][2] != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec[1] << "\n"
             << "   Expected result:\n( 5 5 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Dense vector assignment
   //=====================================================================================

   {
      test_ = "VectorAccessProxy dense vector assignment";

      DV tmp( 3UL );
      tmp[0] = 1;
      tmp[1] = 2;
      tmp[2] = 3;
      DVV vec( 3UL, 1UL );

      vec[1] = tmp;

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkSize    ( vec[0], 0UL );
      checkSize    ( vec[1], 3UL );
      checkCapacity( vec[1], 3UL );
      checkNonZeros( vec[1], 3UL );
      checkSize    ( vec[2], 0UL );

      if( vec[1][0] != 1 || vec[1][1] != 2 || vec[1][2] != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec[1] << "\n"
             << "   Expected result:\n( 1 2 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Sparse vector assignment
   //=====================================================================================

   {
      test_ = "VectorAccessProxy sparse vector assignment";

      SV tmp( 3UL );
      tmp[1] = 2;
      DVV vec( 3UL, 1UL );

      vec[1] = tmp;

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkSize    ( vec[0], 0UL );
      checkSize    ( vec[1], 3UL );
      checkCapacity( vec[1], 3UL );
      checkNonZeros( vec[1], 1UL );
      checkSize    ( vec[2], 0UL );

      if( vec[1][0] != 0 || vec[1][1] != 2 || vec[1][2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec[1] << "\n"
             << "   Expected result:\n( 0 2 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the VectorAccessProxy addition assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the addition assignment operators of the VectorAccessProxy
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ProxyTest::testAddAssign()
{
   //=====================================================================================
   // Dense vector addition assignment
   //=====================================================================================

   {
      test_ = "VectorAccessProxy dense vector addition assignment";

      DV tmp( 3UL );
      tmp[0] = 1;
      tmp[1] = 2;
      tmp[2] = 3;
      DVV vec( 3UL, 1UL );
      vec[1] = tmp;

      vec[1] += tmp;

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkSize    ( vec[0], 0UL );
      checkSize    ( vec[1], 3UL );
      checkCapacity( vec[1], 3UL );
      checkNonZeros( vec[1], 3UL );
      checkSize    ( vec[2], 0UL );

      if( vec[1][0] != 2 || vec[1][1] != 4 || vec[1][2] != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec[1] << "\n"
             << "   Expected result:\n( 2 4 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Sparse vector addition assignment
   //=====================================================================================

   {
      test_ = "VectorAccessProxy sparse vector addition assignment";

      SV tmp( 3UL );
      tmp[1] = 2;
      DVV vec( 3UL, 1UL );
      vec[1] = tmp;

      vec[1] += tmp;

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkSize    ( vec[0], 0UL );
      checkSize    ( vec[1], 3UL );
      checkCapacity( vec[1], 3UL );
      checkNonZeros( vec[1], 1UL );
      checkSize    ( vec[2], 0UL );

      if( vec[1][0] != 0 || vec[1][1] != 4 || vec[1][2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec[1] << "\n"
             << "   Expected result:\n( 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the VectorAccessProxy subtraction assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the subtraction assignment operators of the VectorAccessProxy
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ProxyTest::testSubAssign()
{
   //=====================================================================================
   // Dense vector subtraction assignment
   //=====================================================================================

   {
      test_ = "VectorAccessProxy dense vector subtraction assignment";

      DV tmp( 3UL );
      tmp[0] = 1;
      tmp[1] = 2;
      tmp[2] = 3;
      DVV vec( 3UL, 1UL );
      vec[1] = tmp;

      vec[1] -= tmp;

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkSize    ( vec[0], 0UL );
      checkSize    ( vec[1], 3UL );
      checkCapacity( vec[1], 3UL );
      checkNonZeros( vec[1], 0UL );
      checkSize    ( vec[2], 0UL );

      if( vec[1][0] != 0 || vec[1][1] != 0 || vec[1][2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec[1] << "\n"
             << "   Expected result:\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Sparse vector subtraction assignment
   //=====================================================================================

   {
      test_ = "VectorAccessProxy sparse vector subtraction assignment";

      SV tmp( 3UL );
      tmp[1] = 2;
      DVV vec( 3UL, 1UL );
      vec[1] = tmp;

      vec[1] -= tmp;

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkSize    ( vec[0], 0UL );
      checkSize    ( vec[1], 3UL );
      checkCapacity( vec[1], 1UL );
      checkNonZeros( vec[1], 0UL );
      checkSize    ( vec[2], 0UL );

      if( vec[1][0] != 0 || vec[1][1] != 0 || vec[1][2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec[1] << "\n"
             << "   Expected result:\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the VectorAccessProxy multiplication assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the multiplication assignment operators of the VectorAccessProxy
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ProxyTest::testMultAssign()
{
   //=====================================================================================
   // Dense vector multiplication assignment
   //=====================================================================================

   {
      test_ = "VectorAccessProxy dense vector multiplication assignment";

      DV tmp( 3UL );
      tmp[0] = 1;
      tmp[1] = 2;
      tmp[2] = 3;
      DVV vec( 3UL, 1UL );
      vec[1] = tmp;

      vec[1] *= tmp;

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkSize    ( vec[0], 0UL );
      checkSize    ( vec[1], 3UL );
      checkCapacity( vec[1], 3UL );
      checkNonZeros( vec[1], 3UL );
      checkSize    ( vec[2], 0UL );

      if( vec[1][0] != 1 || vec[1][1] != 4 || vec[1][2] != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec[1] << "\n"
             << "   Expected result:\n( 1 4 9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Sparse vector multiplication assignment
   //=====================================================================================

   {
      test_ = "VectorAccessProxy sparse vector multiplication assignment";

      SV tmp( 3UL );
      tmp[1] = 2;
      DVV vec( 3UL, 1UL );
      vec[1] = tmp;

      vec[1] *= tmp;

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkSize    ( vec[0], 0UL );
      checkSize    ( vec[1], 3UL );
      checkCapacity( vec[1], 3UL );
      checkNonZeros( vec[1], 1UL );
      checkSize    ( vec[2], 0UL );

      if( vec[1][0] != 0 || vec[1][1] != 4 || vec[1][2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec[1] << "\n"
             << "   Expected result:\n( 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the VectorAccessProxy division assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the division assignment operators of the VectorAccessProxy
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ProxyTest::testDivAssign()
{
   //=====================================================================================
   // Dense vector division assignment
   //=====================================================================================

   {
      test_ = "VectorAccessProxy dense vector division assignment";

      DV tmp( 3UL );
      tmp[0] = 1;
      tmp[1] = 2;
      tmp[2] = 3;
      DVV vec( 3UL, 1UL );
      vec[1] = tmp;

      vec[1] /= tmp;

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkSize    ( vec[0], 0UL );
      checkSize    ( vec[1], 3UL );
      checkCapacity( vec[1], 3UL );
      checkNonZeros( vec[1], 3UL );
      checkSize    ( vec[2], 0UL );

      if( vec[1][0] != 1 || vec[1][1] != 1 || vec[1][2] != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec[1] << "\n"
             << "   Expected result:\n( 1 1 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the VectorAccessProxy modulo assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the modulo assignment operators of the VectorAccessProxy
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ProxyTest::testModAssign()
{
   //=====================================================================================
   // Dense vector cross product assignment
   //=====================================================================================

   {
      test_ = "VectorAccessProxy dense vector cross product assignment";

      DV tmp( 3UL );
      tmp[0] = 1;
      tmp[1] = 2;
      tmp[2] = 3;
      DVV vec( 3UL, 1UL );
      vec[1] = tmp;

      vec[1] %= tmp;

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkSize    ( vec[0], 0UL );
      checkSize    ( vec[1], 3UL );
      checkCapacity( vec[1], 3UL );
      checkNonZeros( vec[1], 0UL );
      checkSize    ( vec[2], 0UL );

      if( vec[1][0] != 0 || vec[1][1] != 0 || vec[1][2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec[1] << "\n"
             << "   Expected result:\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Sparse vector cross product assignment
   //=====================================================================================

   {
      test_ = "VectorAccessProxy sparse vector cross product assignment";

      SV tmp( 3UL );
      tmp[1] = 2;
      DVV vec( 3UL, 1UL );
      vec[1] = tmp;

      vec[1] %= tmp;

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkSize    ( vec[0], 0UL );
      checkSize    ( vec[1], 3UL );
      checkCapacity( vec[1], 1UL );
      checkNonZeros( vec[1], 0UL );
      checkSize    ( vec[2], 0UL );

      if( vec[1][0] != 0 || vec[1][1] != 0 || vec[1][2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec[1] << "\n"
             << "   Expected result:\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of all VectorAccessProxy (self-)scaling operations.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all available ways to scale an instance of the VectorAccessProxy
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ProxyTest::testScaling()
{
   //=====================================================================================
   // Self-scaling (v*=s)
   //=====================================================================================

   {
      test_ = "VectorAccessProxy self-scaling (v*=s)";

      DVV vec( 3UL, 1UL );
      vec[1] = DV( 1UL, 2 );

      vec[1] *= 2;

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkSize    ( vec[0], 0UL );
      checkSize    ( vec[1], 1UL );
      checkCapacity( vec[1], 1UL );
      checkNonZeros( vec[1], 1UL );
      checkSize    ( vec[2], 0UL );

      if( vec[1][0] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << vec[1] << "\n"
             << "   Expected result:\n( 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Self-scaling (v/=s)
   //=====================================================================================

   {
      test_ = "VectorAccessProxy self-scaling (v/=s)";

      DVV vec( 3UL, 1UL );
      vec[1] = DV( 1UL, 2 );

      vec[1] /= 2;

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkSize    ( vec[0], 0UL );
      checkSize    ( vec[1], 1UL );
      checkCapacity( vec[1], 1UL );
      checkNonZeros( vec[1], 1UL );
      checkSize    ( vec[2], 0UL );

      if( vec[1][0] != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << vec[1] << "\n"
             << "   Expected result:\n( 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // VectorAccessProxy::scale()
   //=====================================================================================

   {
      test_ = "VectorAccessProxy::scale()";

      DVV vec( 3UL, 1UL );
      vec[1] = DV( 1UL, 2 );

      vec[1].scale( 2 );

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkSize    ( vec[0], 0UL );
      checkSize    ( vec[1], 1UL );
      checkCapacity( vec[1], 1UL );
      checkNonZeros( vec[1], 1UL );
      checkSize    ( vec[2], 0UL );

      if( vec[1][0] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << vec[1] << "\n"
             << "   Expected result:\n( 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the VectorAccessProxy subscript operator.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of adding and accessing elements via the subscript operator
// of the VectorAccessProxy class template. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
void ProxyTest::testSubscript()
{
   test_ = "VectorAccessProxy::operator[]";

   DVV vec( 3UL, 1UL );
   vec[1] = DV( 1UL, 2 );
   vec[1][0] = 3;

   checkSize    ( vec, 3UL );
   checkCapacity( vec, 1UL );
   checkNonZeros( vec, 1UL );

   checkSize    ( vec[0], 0UL );
   checkSize    ( vec[1], 1UL );
   checkCapacity( vec[1], 1UL );
   checkNonZeros( vec[1], 1UL );
   checkSize    ( vec[2], 0UL );

   if( vec[1][0] != DV( 1UL, 3 ) ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << vec[1] << "\n"
          << "   Expected result:\n( 3 )\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the VectorAccessProxy function call operator.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of adding and accessing elements via the function call operator
// of the VectorAccessProxy class template. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
void ProxyTest::testFunctionCall()
{
   test_ = "VectorAccessProxy::operator()";

   DMV vec( 3UL, 1UL );
   vec[1] = DM( 1UL, 1UL, 2 );
   vec[1](0,0) = 3;

   checkSize    ( vec, 3UL );
   checkCapacity( vec, 1UL );
   checkNonZeros( vec, 1UL );

   checkRows    ( vec[0], 0UL );
   checkColumns ( vec[0], 0UL );
   checkRows    ( vec[1], 1UL );
   checkColumns ( vec[1], 1UL );
   checkCapacity( vec[1], 1UL );
   checkNonZeros( vec[1], 1UL );
   checkRows    ( vec[2], 0UL );
   checkColumns ( vec[2], 0UL );

   if( vec[1](0,0) != DM( 1UL, 1UL, 3 ) ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Function call operator failed\n"
          << " Details:\n"
          << "   Result:\n" << vec[1] << "\n"
          << "   Expected result:\n( 3 )\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the VectorAccessProxy iterator implementation.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the iterator implementation of the VectorAccessProxy class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ProxyTest::testIterator()
{
   //=====================================================================================
   // Vector elements
   //=====================================================================================

   {
      DVV vec( 3UL, 1UL );
      vec[1] = DV( 4UL, 4 );

      // Counting the number of elements via Iterator (end-begin)
      {
         test_ = "VectorAccessProxy::begin() and VectorAccessProxy::end()";

         const ptrdiff_t number( end( vec[1] ) - begin( vec[1] ) );

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
         test_ = "VectorAccessProxy::cbegin() and VectorAccessProxy::cend()";

         const ptrdiff_t number( cend( vec[1] ) - cbegin( vec[1] ) );

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
   // Matrix elements
   //=====================================================================================

   {
      DMV vec( 3UL, 1UL );
      vec[1] = DM( 4UL, 4UL, 4 );

      // Counting the number of elements via Iterator (end-begin)
      {
         test_ = "VectorAccessProxy::begin( size_t ) and VectorAccessProxy::end( size_t )";

         const ptrdiff_t number( end( vec[1], 1UL ) - begin( vec[1], 1UL ) );

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
         test_ = "VectorAccessProxy::cbegin( size_t ) and VectorAccessProxy::cend( size_t )";

         const ptrdiff_t number( cend( vec[1], 1UL ) - cbegin( vec[1], 1UL ) );

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
/*!\brief Test of the \c nonZeros() member functions of the VectorAccessProxy class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c nonZeros() member functions of the VectorAccessProxy
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ProxyTest::testNonZeros()
{
   //=====================================================================================
   // Vector elements
   //=====================================================================================

   {
      test_ = "VectorAccessProxy::nonZeros()";

      DVV vec( 3UL, 1UL );
      vec[1] = DV( 8UL, 8 );

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkSize    ( vec[0], 0UL );
      checkSize    ( vec[1], 8UL );
      checkCapacity( vec[1], 8UL );
      checkNonZeros( vec[1], 8UL );
      checkSize    ( vec[2], 0UL );
   }


   //=====================================================================================
   // Matrix elements
   //=====================================================================================

   {
      test_ = "VectorAccessProxy::nonZeros( size_t )";

      DMV vec( 3UL, 1UL );
      vec[1] = DM( 3UL, 3UL, 3 );

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkRows    ( vec[0],  0UL );
      checkColumns ( vec[0],  0UL );
      checkRows    ( vec[1],  3UL );
      checkColumns ( vec[1],  3UL );
      checkCapacity( vec[1],  9UL );
      checkNonZeros( vec[1],  0UL, 3UL );
      checkNonZeros( vec[1],  1UL, 3UL );
      checkNonZeros( vec[1],  2UL, 3UL );
      checkRows    ( vec[2],  0UL );
      checkColumns ( vec[2],  0UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c reset() member functions of the VectorAccessProxy class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c reset() member functions of the VectorAccessProxy
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ProxyTest::testReset()
{
   //=====================================================================================
   // Vector elements
   //=====================================================================================

   {
      test_ = "VectorAccessProxy::reset()";

      DVV vec( 3UL, 1UL );
      vec[1] = DV( 8UL, 8 );
      vec[1].reset();

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkSize    ( vec[0], 0UL );
      checkSize    ( vec[1], 8UL );
      checkCapacity( vec[1], 8UL );
      checkNonZeros( vec[1], 0UL );
      checkSize    ( vec[2], 0UL );
   }


   //=====================================================================================
   // Matrix elements
   //=====================================================================================

   {
      test_ = "VectorAccessProxy::reset( size_t )";

      DMV vec( 3UL, 1UL );
      vec[1] = DM( 3UL, 3UL, 3 );
      vec[1].reset( 1UL );

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkRows    ( vec[0], 0UL );
      checkColumns ( vec[0], 0UL );
      checkRows    ( vec[1], 3UL );
      checkColumns ( vec[1], 3UL );
      checkCapacity( vec[1], 9UL );
      checkNonZeros( vec[1], 6UL );
      checkNonZeros( vec[1], 0UL, 3UL );
      checkNonZeros( vec[1], 1UL, 0UL );
      checkNonZeros( vec[1], 2UL, 3UL );
      checkRows    ( vec[2], 0UL );
      checkColumns ( vec[2], 0UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c clear() member function of the VectorAccessProxy class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c clear() member function of the VectorAccessProxy
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ProxyTest::testClear()
{
   test_ = "VectorAccessProxy::clear()";

   DVV vec( 3UL, 1UL );
   vec[1] = DV( 8UL, 8 );
   vec[1].clear();

   checkSize    ( vec, 3UL );
   checkCapacity( vec, 1UL );
   checkNonZeros( vec, 0UL );

   checkSize( vec[0], 0UL );
   checkSize( vec[1], 0UL );
   checkSize( vec[2], 0UL );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c resize() member functions of the VectorAccessProxy class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c resize() member functions of the VectorAccessProxy
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ProxyTest::testResize()
{
   //=====================================================================================
   // Vector elements
   //=====================================================================================

   {
      test_ = "VectorAccessProxy::resize( size_t )";

      DVV vec( 3UL, 1UL );
      vec[1].resize( 10UL );

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkSize    ( vec[0],  0UL );
      checkSize    ( vec[1], 10UL );
      checkCapacity( vec[1], 10UL );
      checkSize    ( vec[2],  0UL );
   }

   {
      test_ = "resize( VectorAccessProxy, size_t )";

      DVV vec( 3UL, 1UL );
      resize( vec[1], 10UL );

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkSize    ( vec[0],  0UL );
      checkSize    ( vec[1], 10UL );
      checkCapacity( vec[1], 10UL );
      checkSize    ( vec[2],  0UL );
   }


   //=====================================================================================
   // Matrix elements
   //=====================================================================================

   {
      test_ = "VectorAccessProxy::resize( size_t, size_t )";

      DMV vec( 3UL, 1UL );
      vec[1].resize( 5UL, 5UL );

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkRows    ( vec[0],  0UL );
      checkColumns ( vec[0],  0UL );
      checkRows    ( vec[1],  5UL );
      checkColumns ( vec[1],  5UL );
      checkCapacity( vec[1], 25UL );
      checkRows    ( vec[2],  0UL );
      checkColumns ( vec[2],  0UL );
   }

   {
      test_ = "resize( VectorAccessProxy, size_t, size_t )";

      DMV vec( 3UL, 1UL );
      resize( vec[1], 5UL, 5UL );

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkRows    ( vec[0],  0UL );
      checkColumns ( vec[0],  0UL );
      checkRows    ( vec[1],  5UL );
      checkColumns ( vec[1],  5UL );
      checkCapacity( vec[1], 25UL );
      checkRows    ( vec[2],  0UL );
      checkColumns ( vec[2],  0UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c extend() member functions of the VectorAccessProxy class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c extend() member functions of the VectorAccessProxy
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ProxyTest::testExtend()
{
   //=====================================================================================
   // Vector elements
   //=====================================================================================

   {
      test_ = "VectorAccessProxy::extend( size_t )";

      DVV vec( 3UL, 1UL );
      vec[1].extend( 10UL );

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkSize    ( vec[0],  0UL );
      checkSize    ( vec[1], 10UL );
      checkCapacity( vec[1], 10UL );
      checkSize    ( vec[2],  0UL );
   }


   //=====================================================================================
   // Matrix elements
   //=====================================================================================

   {
      test_ = "VectorAccessProxy::extend( size_t, size_t )";

      DMV vec( 3UL, 1UL );
      vec[1].extend( 5UL, 5UL );

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkRows    ( vec[0],  0UL );
      checkColumns ( vec[0],  0UL );
      checkRows    ( vec[1],  5UL );
      checkColumns ( vec[1],  5UL );
      checkCapacity( vec[1], 25UL );
      checkRows    ( vec[2],  0UL );
      checkColumns ( vec[2],  0UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c reserve() member functions of the VectorAccessProxy class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c reserve() member functions of the VectorAccessProxy
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ProxyTest::testReserve()
{
   //=====================================================================================
   // Vector elements
   //=====================================================================================

   {
      test_ = "VectorAccessProxy::reserve( size_t )";

      DVV vec( 3UL, 1UL );
      vec[0].resize( 5UL );
      vec[0].reserve( 10UL );

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkSize    ( vec[0],  5UL );
      checkCapacity( vec[0], 10UL );
      checkSize    ( vec[1],  0UL );
      checkSize    ( vec[2],  0UL );
   }


   //=====================================================================================
   // Matrix elements
   //=====================================================================================

   {
      test_ = "VectorAccessProxy::reserve( size_t, size_t )";

      SMV vec( 3UL, 1UL );
      vec[0] = SM( 2UL, 2UL, 1UL );
      vec[0].reserve( 0UL, 1UL );

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkRows    ( vec[0], 2UL );
      checkColumns ( vec[0], 2UL );
      checkCapacity( vec[0], 1UL );
      checkCapacity( vec[0], 0UL, 1UL );
      checkCapacity( vec[0], 0UL, 0UL );
      checkRows    ( vec[1], 0UL );
      checkColumns ( vec[1], 0UL );
      checkRows    ( vec[2], 0UL );
      checkColumns ( vec[2], 0UL );
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
   {
      test_ = "VectorAccessProxy::trim()";

      SMV vec( 3UL, 3UL );
      vec[0].resize( 2UL, 2UL );
      vec[0].reserve( 10UL );
      vec[0].reserve( 0UL, 6UL );
      vec[0].reserve( 1UL, 4UL );

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 3UL );
      checkNonZeros( vec, 1UL );

      checkRows    ( vec[0],  2UL );
      checkColumns ( vec[0],  2UL );
      checkCapacity( vec[0], 10UL );
      checkCapacity( vec[0],  0UL, 6UL );
      checkCapacity( vec[0],  1UL, 4UL );
      checkRows    ( vec[1],  0UL );
      checkColumns ( vec[1],  0UL );
      checkRows    ( vec[2],  0UL );
      checkColumns ( vec[2],  0UL );

      vec[0].trim();

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 3UL );
      checkNonZeros( vec, 1UL );

      checkRows    ( vec[0],  2UL );
      checkColumns ( vec[0],  2UL );
      checkCapacity( vec[0], 10UL );
      checkCapacity( vec[0],  0UL, 0UL );
      checkCapacity( vec[0],  1UL, 0UL );
      checkRows    ( vec[1],  0UL );
      checkColumns ( vec[1],  0UL );
      checkRows    ( vec[2],  0UL );
      checkColumns ( vec[2],  0UL );
   }

   {
      test_ = "VectorAccessProxy::trim( size_t )";

      SMV vec( 3UL, 3UL );
      vec[0].resize( 2UL, 2UL );
      vec[0].reserve( 10UL );
      vec[0].reserve( 0UL, 6UL );
      vec[0].reserve( 1UL, 4UL );

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkRows    ( vec[0],  2UL );
      checkColumns ( vec[0],  2UL );
      checkCapacity( vec[0], 10UL );
      checkCapacity( vec[0],  0UL, 6UL );
      checkCapacity( vec[0],  1UL, 4UL );
      checkRows    ( vec[1],  0UL );
      checkColumns ( vec[1],  0UL );
      checkRows    ( vec[2],  0UL );
      checkColumns ( vec[2],  0UL );

      vec[0].trim( 0UL );

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkRows    ( vec[0],  2UL );
      checkColumns ( vec[0],  2UL );
      checkCapacity( vec[0], 10UL );
      checkCapacity( vec[0],  0UL, 0UL );
      checkCapacity( vec[0],  1UL, 4UL );
      checkRows    ( vec[1],  0UL );
      checkColumns ( vec[1],  0UL );
      checkRows    ( vec[2],  0UL );
      checkColumns ( vec[2],  0UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c swap() functionality of the VectorAccessProxy class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c swap() function of the VectorAccessProxy class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ProxyTest::testSwap()
{
   test_ = "VectorAccessProxy swap";

   {
      DVV vec( 3UL, 2UL );
      vec[0].resize( 2UL );
      vec[0] = 0;
      vec[2].resize( 6UL );
      vec[2] = 0;

      swap( vec[0], vec[2] );

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 2UL );
      checkNonZeros( vec, 2UL );

      checkSize    ( vec[0], 6UL );
      checkCapacity( vec[0], 6UL );
      checkNonZeros( vec[0], 0UL );
      checkSize    ( vec[1], 0UL );
      checkSize    ( vec[2], 2UL );
      checkCapacity( vec[2], 2UL );
      checkNonZeros( vec[2], 0UL );
   }

   {
      DVV vec( 3UL, 1UL );
      vec[1] = DV( 2UL, 2 );
      DV tmp( 6UL, 6 );

      swap( vec[1], tmp );

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkSize    ( vec[0], 0UL );
      checkSize    ( vec[1], 6UL );
      checkCapacity( vec[1], 6UL );
      checkNonZeros( vec[1], 6UL );
      checkSize    ( vec[2], 0UL );
      checkSize    ( tmp   , 2UL );
      checkCapacity( tmp   , 2UL );
      checkNonZeros( tmp   , 2UL );
   }

   {
      DVV vec( 3UL, 1UL );
      vec[1] = DV( 2UL, 2 );
      DV tmp( 6UL, 6 );

      swap( tmp, vec[1] );

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkSize    ( vec[0], 0UL );
      checkSize    ( vec[1], 6UL );
      checkCapacity( vec[1], 6UL );
      checkNonZeros( vec[1], 6UL );
      checkSize    ( vec[2], 0UL );
      checkSize    ( tmp   , 2UL );
      checkCapacity( tmp   , 2UL );
      checkNonZeros( tmp   , 2UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c set() member functions of the VectorAccessProxy class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c set() member functions of the VectorAccessProxy
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ProxyTest::testSet()
{
   //=====================================================================================
   // Vector elements
   //=====================================================================================

   {
      test_ = "VectorAccessProxy::set( size_t, ElementType )";

      SVV vec( 3UL, 1UL );
      vec[1] = SV( 3UL, 1UL );
      vec[1].set( 1UL, 5 );

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkSize    ( vec[0], 0UL );
      checkSize    ( vec[1], 3UL );
      checkCapacity( vec[1], 1UL );
      checkNonZeros( vec[1], 1UL );
      checkSize    ( vec[2], 0UL );

      if( vec[1][1] != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setting an element failed\n"
             << " Details:\n"
             << "   Result:\n" << vec[1] << "\n"
             << "   Expected result:\n( 0 5 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Matrix elements
   //=====================================================================================

   {
      test_ = "VectorAccessProxy::set( size_t, size_t, ElementType )";

      SMV vec( 3UL, 1UL );
      vec[1] = SM( 2UL, 2UL, 1UL );
      vec[1].set( 0UL, 1UL, 5 );

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkRows    ( vec[0], 0UL );
      checkColumns ( vec[0], 0UL );
      checkRows    ( vec[1], 2UL );
      checkColumns ( vec[1], 2UL );
      checkCapacity( vec[1], 1UL );
      checkNonZeros( vec[1], 1UL );
      checkRows    ( vec[2], 0UL );
      checkColumns ( vec[2], 0UL );

      if( vec[1](0,1) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setting an element failed\n"
             << " Details:\n"
             << "   Result:\n" << vec[1] << "\n"
             << "   Expected result:\n( 0 5 )\n( 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c insert() member functions of the VectorAccessProxy class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c insert() member functions of the VectorAccessProxy
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ProxyTest::testInsert()
{
   //=====================================================================================
   // Vector elements
   //=====================================================================================

   {
      test_ = "VectorAccessProxy::insert( size_t, ElementType )";

      SVV vec( 3UL, 1UL );
      vec[1] = SV( 3UL, 1UL );
      vec[1].insert( 1UL, 5 );

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkSize    ( vec[0], 0UL );
      checkSize    ( vec[1], 3UL );
      checkCapacity( vec[1], 1UL );
      checkNonZeros( vec[1], 1UL );
      checkSize    ( vec[2], 0UL );

      if( vec[1][1] != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Inserting an element failed\n"
             << " Details:\n"
             << "   Result:\n" << vec[1] << "\n"
             << "   Expected result:\n( 0 5 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Matrix elements
   //=====================================================================================

   {
      test_ = "VectorAccessProxy::insert( size_t, size_t, ElementType )";

      SMV vec( 3UL, 1UL );
      vec[1] = SM( 2UL, 2UL, 1UL );
      vec[1].insert( 0UL, 1UL, 5 );

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkRows    ( vec[0], 0UL );
      checkColumns ( vec[0], 0UL );
      checkRows    ( vec[1], 2UL );
      checkColumns ( vec[1], 2UL );
      checkCapacity( vec[1], 1UL );
      checkNonZeros( vec[1], 1UL );
      checkRows    ( vec[2], 0UL );
      checkColumns ( vec[2], 0UL );

      if( vec[1](0,1) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Inserting an element failed\n"
             << " Details:\n"
             << "   Result:\n" << vec[1] << "\n"
             << "   Expected result:\n( 0 5 )\n( 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c append() member functions of the VectorAccessProxy class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c append() member functions of the VectorAccessProxy
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ProxyTest::testAppend()
{
   //=====================================================================================
   // Vector elements
   //=====================================================================================

   {
      test_ = "VectorAccessProxy::append( size_t, ElementType )";

      SVV vec( 3UL, 1UL );
      vec[1] = SV( 3UL );
      vec[1].reserve( 1UL );
      vec[1].append( 1UL, 5 );

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkSize    ( vec[0], 0UL );
      checkSize    ( vec[1], 3UL );
      checkCapacity( vec[1], 1UL );
      checkNonZeros( vec[1], 1UL );
      checkSize    ( vec[2], 0UL );

      if( vec[1][1] != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Append operation failed\n"
             << " Details:\n"
             << "   Result:\n" << vec[1] << "\n"
             << "   Expected result:\n( 0 5 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Matrix elements
   //=====================================================================================

   {
      test_ = "VectorAccessProxy::append( size_t, size_t, ElementType )";

      SMV vec( 3UL, 1UL );
      vec[1] = SM( 2UL, 2UL );
      vec[1].reserve( 0UL, 1UL );
      vec[1].append( 0UL, 1UL, 5 );
      vec[1].finalize( 0UL );

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkRows    ( vec[0], 0UL );
      checkColumns ( vec[0], 0UL );
      checkRows    ( vec[1], 2UL );
      checkColumns ( vec[1], 2UL );
      checkCapacity( vec[1], 1UL );
      checkNonZeros( vec[1], 1UL );
      checkRows    ( vec[2], 0UL );
      checkColumns ( vec[2], 0UL );

      if( vec[1](0,1) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Append operation failed\n"
             << " Details:\n"
             << "   Result:\n" << vec[1] << "\n"
             << "   Expected result:\n( 0 5 )\n( 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c erase() member functions of the VectorAccessProxy class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c erase() member functions of the VectorAccessProxy
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ProxyTest::testErase()
{
   //=====================================================================================
   // Vector elements
   //=====================================================================================

   {
      test_ = "VectorAccessProxy::erase( size_t )";

      SVV vec( 3UL, 1UL );
      vec[1] = SV( 3UL, 1UL );
      vec[1].insert( 1UL, 5 );
      vec[1].erase( 1UL );

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkSize    ( vec[0], 0UL );
      checkSize    ( vec[1], 3UL );
      checkNonZeros( vec[1], 0UL );
      checkSize    ( vec[2], 0UL );
   }

   {
      test_ = "VectorAccessProxy::erase( Iterator )";

      SVV vec( 3UL, 1UL );
      vec[1] = SV( 3UL, 1UL );
      vec[1].erase( vec[1].insert( 1UL, 5 ) );

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkSize    ( vec[0], 0UL );
      checkSize    ( vec[1], 3UL );
      checkNonZeros( vec[1], 0UL );
      checkSize    ( vec[2], 0UL );
   }

   {
      test_ = "VectorAccessProxy::erase( Iterator, Iterator )";

      SVV vec( 3UL, 1UL );
      vec[1] = SV( 3UL, 1UL );
      vec[1].insert( 0UL, 1 );
      vec[1].insert( 1UL, 2 );
      vec[1].insert( 2UL, 3 );
      vec[1].erase( begin( vec[1] ), end( vec[1] ) );

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkSize    ( vec[0], 0UL );
      checkSize    ( vec[1], 3UL );
      checkNonZeros( vec[1], 0UL );
      checkSize    ( vec[2], 0UL );
   }

   {
      test_ = "VectorAccessProxy::erase( Predicate )";

      SVV vec( 3UL, 1UL );
      vec[1] = SV( 3UL, 1UL );
      vec[1].insert( 1UL, 5 );
      vec[1].erase( []( int value ){ return value == 5; } );

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkSize    ( vec[0], 0UL );
      checkSize    ( vec[1], 3UL );
      checkNonZeros( vec[1], 0UL );
      checkSize    ( vec[2], 0UL );
   }

   {
      test_ = "VectorAccessProxy::erase( Iterator, Iterator, Predicate )";

      SVV vec( 3UL, 1UL );
      vec[1] = SV( 3UL, 1UL );
      vec[1].insert( 1UL, 5 );
      vec[1].erase( begin( vec[1] ), end( vec[1] ), []( int value ){ return value == 5; } );

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkSize    ( vec[0], 0UL );
      checkSize    ( vec[1], 3UL );
      checkNonZeros( vec[1], 0UL );
      checkSize    ( vec[2], 0UL );
   }


   //=====================================================================================
   // Matrix elements
   //=====================================================================================

   {
      test_ = "VectorAccessProxy::erase( size_t, size_t )";

      SMV vec( 3UL, 1UL );
      vec[1] = SM( 2UL, 2UL, 1UL );
      vec[1].insert( 0UL, 1UL, 5 );
      vec[1].erase( 0UL, 1UL );

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkRows    ( vec[0], 0UL );
      checkColumns ( vec[0], 0UL );
      checkRows    ( vec[1], 2UL );
      checkColumns ( vec[1], 2UL );
      checkNonZeros( vec[1], 0UL );
      checkRows    ( vec[2], 0UL );
      checkColumns ( vec[2], 0UL );
   }

   {
      test_ = "VectorAccessProxy::erase( size_t, Iterator )";

      SMV vec( 3UL, 1UL );
      vec[1] = SM( 2UL, 2UL, 1UL );
      vec[1].erase( 0UL, vec[1].insert( 0UL, 1UL, 5 ) );

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkRows    ( vec[0], 0UL );
      checkColumns ( vec[0], 0UL );
      checkRows    ( vec[1], 2UL );
      checkColumns ( vec[1], 2UL );
      checkNonZeros( vec[1], 0UL );
      checkRows    ( vec[2], 0UL );
      checkColumns ( vec[2], 0UL );
   }

   {
      test_ = "VectorAccessProxy::erase( size_t, Iterator, Iterator )";

      SMV vec( 3UL, 1UL );
      vec[1] = SM( 2UL, 2UL, 1UL );
      vec[1].insert( 0UL, 0UL, 1 );
      vec[1].insert( 0UL, 1UL, 2 );
      vec[1].erase( 0UL, begin( vec[1], 0UL ), end( vec[1], 0UL ) );

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkRows    ( vec[0], 0UL );
      checkColumns ( vec[0], 0UL );
      checkRows    ( vec[1], 2UL );
      checkColumns ( vec[1], 2UL );
      checkNonZeros( vec[1], 0UL );
      checkRows    ( vec[2], 0UL );
      checkColumns ( vec[2], 0UL );
   }

   {
      test_ = "VectorAccessProxy::erase( Predicate )";

      SMV vec( 3UL, 1UL );
      vec[1] = SM( 2UL, 2UL, 1UL );
      vec[1].insert( 0UL, 1UL, 5 );
      vec[1].erase( []( int value ){ return value == 5; } );

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkRows    ( vec[0], 0UL );
      checkColumns ( vec[0], 0UL );
      checkRows    ( vec[1], 2UL );
      checkColumns ( vec[1], 2UL );
      checkNonZeros( vec[1], 0UL );
      checkRows    ( vec[2], 0UL );
      checkColumns ( vec[2], 0UL );
   }

   {
      test_ = "VectorAccessProxy::erase( size_t, Iterator, Iterator, Predicate )";

      SMV vec( 3UL, 1UL );
      vec[1] = SM( 2UL, 2UL, 1UL );
      vec[1].insert( 0UL, 1UL, 5 );
      vec[1].erase( 0UL, begin( vec[1], 0UL ), end( vec[1], 0UL ),
                    []( int value ){ return value == 5; } );

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkRows    ( vec[0], 0UL );
      checkColumns ( vec[0], 0UL );
      checkRows    ( vec[1], 2UL );
      checkColumns ( vec[1], 2UL );
      checkNonZeros( vec[1], 0UL );
      checkRows    ( vec[2], 0UL );
      checkColumns ( vec[2], 0UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c find() member functions of the VectorAccessProxy class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c find() member functions of the VectorAccessProxy
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ProxyTest::testFind()
{
   //=====================================================================================
   // Vector elements
   //=====================================================================================

   {
      test_ = "VectorAccessProxy::find( size_t )";

      SVV vec( 3UL, 1UL );
      vec[0] = SV( 5UL, 3UL );
      vec[0][1] = 2;
      vec[0][2] = 3;
      vec[0][4] = 5;

      SV::Iterator pos( vec[0].find( 2UL ) );

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkSize    ( vec[0], 5UL );
      checkCapacity( vec[0], 3UL );
      checkNonZeros( vec[0], 3UL );
      checkSize    ( vec[1], 0UL );
      checkSize    ( vec[2], 0UL );

      if( pos == vec[0].end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Element could not be found\n"
             << " Details:\n"
             << "   Required index = 2\n"
             << "   Current vector:\n" << vec[0] << "\n";
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
             << "   Current vector:\n" << vec[0] << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "find( VectorAccessProxy, size_t )";

      SVV vec( 3UL, 1UL );
      vec[0] = SV( 5UL, 3UL );
      vec[0][1] = 2;
      vec[0][2] = 3;
      vec[0][4] = 5;

      SV::Iterator pos( find( vec[0], 2UL ) );

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkSize    ( vec[0], 5UL );
      checkCapacity( vec[0], 3UL );
      checkNonZeros( vec[0], 3UL );
      checkSize    ( vec[1], 0UL );
      checkSize    ( vec[2], 0UL );

      if( pos == vec[0].end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Element could not be found\n"
             << " Details:\n"
             << "   Required index = 2\n"
             << "   Current vector:\n" << vec[0] << "\n";
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
             << "   Current vector:\n" << vec[0] << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Matrix elements
   //=====================================================================================

   {
      test_ = "VectorAccessProxy::find( size_t, size_t )";

      SMV vec( 3UL, 1UL );
      vec[0] = SM( 2UL, 5UL, 3UL );
      vec[0](1,1) = 2;
      vec[0](1,2) = 3;
      vec[0](1,4) = 5;

      SM::Iterator pos( vec[0].find( 1UL, 2UL ) );

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkRows    ( vec[0], 2UL );
      checkColumns ( vec[0], 5UL );
      checkCapacity( vec[0], 3UL );
      checkNonZeros( vec[0], 3UL );
      checkRows    ( vec[1], 0UL );
      checkColumns ( vec[1], 0UL );
      checkRows    ( vec[2], 0UL );
      checkColumns ( vec[2], 0UL );

      if( pos == vec[0].end( 1UL ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Element could not be found\n"
             << " Details:\n"
             << "   Required index = 2\n"
             << "   Current matrix:\n" << vec[0] << "\n";
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
             << "   Current matrix:\n" << vec[0] << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "find( VectorAccessProxy, size_t, size_t )";

      SMV vec( 3UL, 1UL );
      vec[0] = SM( 2UL, 5UL, 3UL );
      vec[0](1,1) = 2;
      vec[0](1,2) = 3;
      vec[0](1,4) = 5;

      SM::Iterator pos( find( vec[0], 1UL, 2UL ) );

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkRows    ( vec[0], 2UL );
      checkColumns ( vec[0], 5UL );
      checkCapacity( vec[0], 3UL );
      checkNonZeros( vec[0], 3UL );
      checkRows    ( vec[1], 0UL );
      checkColumns ( vec[1], 0UL );
      checkRows    ( vec[2], 0UL );
      checkColumns ( vec[2], 0UL );

      if( pos == vec[0].end( 1UL ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Element could not be found\n"
             << " Details:\n"
             << "   Required index = 2\n"
             << "   Current matrix:\n" << vec[0] << "\n";
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
             << "   Current matrix:\n" << vec[0] << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c lowerBound() member functions of the VectorAccessProxy class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c lowerBound() member functions of the VectorAccessProxy
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ProxyTest::testLowerBound()
{
   //=====================================================================================
   // Vector elements
   //=====================================================================================

   {
      test_ = "VectorAccessProxy::lowerBound( size_t )";

      SVV vec( 3UL, 1UL );
      vec[0] = SV( 5UL, 3UL );
      vec[0][1] = 2;
      vec[0][2] = 3;
      vec[0][4] = 5;

      SV::Iterator pos( vec[0].lowerBound( 3UL ) );

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkSize    ( vec[0], 5UL );
      checkCapacity( vec[0], 3UL );
      checkNonZeros( vec[0], 3UL );
      checkSize    ( vec[1], 0UL );
      checkSize    ( vec[2], 0UL );

      if( pos == vec[0].end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Element could not be found\n"
             << " Details:\n"
             << "   Required index = 3\n"
             << "   Current vector:\n" << vec[0] << "\n";
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
             << "   Current vector:\n" << vec[0] << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "lowerBound( VectorAccessProxy, size_t )";

      SVV vec( 3UL, 1UL );
      vec[0] = SV( 5UL, 3UL );
      vec[0][1] = 2;
      vec[0][2] = 3;
      vec[0][4] = 5;

      SV::Iterator pos( lowerBound( vec[0], 3UL ) );

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkSize    ( vec[0], 5UL );
      checkCapacity( vec[0], 3UL );
      checkNonZeros( vec[0], 3UL );
      checkSize    ( vec[1], 0UL );
      checkSize    ( vec[2], 0UL );

      if( pos == vec[0].end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Element could not be found\n"
             << " Details:\n"
             << "   Required index = 3\n"
             << "   Current vector:\n" << vec[0] << "\n";
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
             << "   Current vector:\n" << vec[0] << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Matrix elements
   //=====================================================================================

   {
      test_ = "VectorAccessProxy::lowerBound( size_t, size_t )";

      SMV vec( 3UL, 1UL );
      vec[0] = SM( 2UL, 5UL, 3UL );
      vec[0](1,1) = 2;
      vec[0](1,2) = 3;
      vec[0](1,4) = 5;

      SM::Iterator pos( vec[0].lowerBound( 1UL, 3 ) );

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkRows    ( vec[0], 2UL );
      checkColumns ( vec[0], 5UL );
      checkCapacity( vec[0], 3UL );
      checkNonZeros( vec[0], 3UL );
      checkRows    ( vec[1], 0UL );
      checkColumns ( vec[1], 0UL );
      checkRows    ( vec[2], 0UL );
      checkColumns ( vec[2], 0UL );

      if( pos == vec[0].end( 1UL ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Element could not be found\n"
             << " Details:\n"
             << "   Required index = 3\n"
             << "   Current matrix:\n" << vec[0] << "\n";
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
             << "   Current matrix:\n" << vec[0] << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "lowerBound( VectorAccessProxy, size_t, size_t )";

      SMV vec( 3UL, 1UL );
      vec[0] = SM( 2UL, 5UL, 3UL );
      vec[0](1,1) = 2;
      vec[0](1,2) = 3;
      vec[0](1,4) = 5;

      SM::Iterator pos( lowerBound( vec[0], 1UL, 3UL ) );

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkRows    ( vec[0], 2UL );
      checkColumns ( vec[0], 5UL );
      checkCapacity( vec[0], 3UL );
      checkNonZeros( vec[0], 3UL );
      checkRows    ( vec[1], 0UL );
      checkColumns ( vec[1], 0UL );
      checkRows    ( vec[2], 0UL );
      checkColumns ( vec[2], 0UL );

      if( pos == vec[0].end( 1UL ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Element could not be found\n"
             << " Details:\n"
             << "   Required index = 3\n"
             << "   Current matrix:\n" << vec[0] << "\n";
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
             << "   Current matrix:\n" << vec[0] << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c upperBound() member functions of the VectorAccessProxy class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c upperBound() member functions of the VectorAccessProxy
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ProxyTest::testUpperBound()
{
   //=====================================================================================
   // Vector elements
   //=====================================================================================

   {
      test_ = "VectorAccessProxy::upperBound( size_t )";

      SVV vec( 3UL, 1UL );
      vec[0] = SV( 5UL, 3UL );
      vec[0][1] = 2;
      vec[0][2] = 3;
      vec[0][4] = 5;

      SV::Iterator pos( vec[0].upperBound( 3UL ) );

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkSize    ( vec[0], 5UL );
      checkCapacity( vec[0], 3UL );
      checkNonZeros( vec[0], 3UL );
      checkSize    ( vec[1], 0UL );
      checkSize    ( vec[2], 0UL );

      if( pos == vec[0].end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Element could not be found\n"
             << " Details:\n"
             << "   Required index = 3\n"
             << "   Current vector:\n" << vec[0] << "\n";
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
             << "   Current vector:\n" << vec[0] << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "upperBound( VectorAccessProxy, size_t )";

      SVV vec( 3UL, 1UL );
      vec[0] = SV( 5UL, 3UL );
      vec[0][1] = 2;
      vec[0][2] = 3;
      vec[0][4] = 5;

      SV::Iterator pos( upperBound( vec[0], 3UL ) );

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkSize    ( vec[0], 5UL );
      checkCapacity( vec[0], 3UL );
      checkNonZeros( vec[0], 3UL );
      checkSize    ( vec[1], 0UL );
      checkSize    ( vec[2], 0UL );

      if( pos == vec[0].end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Element could not be found\n"
             << " Details:\n"
             << "   Required index = 3\n"
             << "   Current vector:\n" << vec[0] << "\n";
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
             << "   Current vector:\n" << vec[0] << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Matrix elements
   //=====================================================================================

   {
      test_ = "VectorAccessProxy::upperBound( size_t, size_t )";

      SMV vec( 3UL, 1UL );
      vec[0] = SM( 2UL, 5UL, 3UL );
      vec[0](1,1) = 2;
      vec[0](1,2) = 3;
      vec[0](1,4) = 5;

      SM::Iterator pos( vec[0].upperBound( 1UL, 3 ) );

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkRows    ( vec[0], 2UL );
      checkColumns ( vec[0], 5UL );
      checkCapacity( vec[0], 3UL );
      checkNonZeros( vec[0], 3UL );
      checkRows    ( vec[1], 0UL );
      checkColumns ( vec[1], 0UL );
      checkRows    ( vec[2], 0UL );
      checkColumns ( vec[2], 0UL );

      if( pos == vec[0].end( 1UL ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Element could not be found\n"
             << " Details:\n"
             << "   Required index = 3\n"
             << "   Current matrix:\n" << vec[0] << "\n";
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
             << "   Current matrix:\n" << vec[0] << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "upperBound( VectorAccessProxy, size_t, size_t )";

      SMV vec( 3UL, 1UL );
      vec[0] = SM( 2UL, 5UL, 3UL );
      vec[0](1,1) = 2;
      vec[0](1,2) = 3;
      vec[0](1,4) = 5;

      SM::Iterator pos( upperBound( vec[0], 1UL, 3UL ) );

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkRows    ( vec[0], 2UL );
      checkColumns ( vec[0], 5UL );
      checkCapacity( vec[0], 3UL );
      checkNonZeros( vec[0], 3UL );
      checkRows    ( vec[1], 0UL );
      checkColumns ( vec[1], 0UL );
      checkRows    ( vec[2], 0UL );
      checkColumns ( vec[2], 0UL );

      if( pos == vec[0].end( 1UL ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Element could not be found\n"
             << " Details:\n"
             << "   Required index = 3\n"
             << "   Current matrix:\n" << vec[0] << "\n";
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
             << "   Current matrix:\n" << vec[0] << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c transpose() functions of the VectorAccessProxy class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c transpose() functions of the VectorAccessProxy class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ProxyTest::testTranspose()
{
   {
      test_ = "VectorAccessProxy::transpose()";

      DMV vec( 3UL, 1UL );
      vec[0].resize( 5UL, 3UL );
      vec[0].transpose();

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkRows    ( vec[0],  3UL );
      checkColumns ( vec[0],  5UL );
      checkCapacity( vec[0], 15UL );
      checkRows    ( vec[1],  0UL );
      checkColumns ( vec[1],  0UL );
      checkRows    ( vec[2],  0UL );
      checkColumns ( vec[2],  0UL );
   }

   {
      test_ = "transpose( VectorAccessProxy )";

      DMV vec( 3UL, 1UL );
      vec[0].resize( 5UL, 3UL );
      transpose( vec[0] );

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkRows    ( vec[0],  3UL );
      checkColumns ( vec[0],  5UL );
      checkCapacity( vec[0], 15UL );
      checkRows    ( vec[1],  0UL );
      checkColumns ( vec[1],  0UL );
      checkRows    ( vec[2],  0UL );
      checkColumns ( vec[2],  0UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c ctranspose() functions of the VectorAccessProxy class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c ctranspose() functions of the VectorAccessProxy class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ProxyTest::testCTranspose()
{
   {
      test_ = "VectorAccessProxy::ctranspose()";

      DMV vec( 3UL, 1UL );
      vec[0].resize( 5UL, 3UL );
      vec[0].ctranspose();

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkRows    ( vec[0],  3UL );
      checkColumns ( vec[0],  5UL );
      checkCapacity( vec[0], 15UL );
      checkRows    ( vec[1],  0UL );
      checkColumns ( vec[1],  0UL );
      checkRows    ( vec[2],  0UL );
      checkColumns ( vec[2],  0UL );
   }

   {
      test_ = "ctranspose( VectorAccessProxy )";

      DMV vec( 3UL, 1UL );
      vec[0].resize( 5UL, 3UL );
      ctranspose( vec[0] );

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkRows    ( vec[0],  3UL );
      checkColumns ( vec[0],  5UL );
      checkCapacity( vec[0], 15UL );
      checkRows    ( vec[1],  0UL );
      checkColumns ( vec[1],  0UL );
      checkRows    ( vec[2],  0UL );
      checkColumns ( vec[2],  0UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c invert() function with the VectorAccessProxy class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c invert() functions with the VectorAccessProxy class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ProxyTest::testInvert()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   using blaze::invert;
   using blaze::byLU;
   using blaze::byLLH;


   {
      test_ = "invert( VectorAccessProxy )";

      blaze::CompressedVector< blaze::DynamicMatrix<double> > vec( 3UL, 1UL );
      vec[0].resize( 3UL, 3UL );
      vec[0] = 0.0;
      vec[0](0,0) = 1.0;
      vec[0](1,1) = 1.0;
      vec[0](2,2) = 1.0;
      invert( vec[0] );

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkRows    ( vec[0], 3UL );
      checkColumns ( vec[0], 3UL );
      checkCapacity( vec[0], 9UL );
      checkNonZeros( vec[0], 3UL );
      checkRows    ( vec[1], 0UL );
      checkColumns ( vec[1], 0UL );
      checkRows    ( vec[2], 0UL );
      checkColumns ( vec[2], 0UL );
   }

   {
      test_ = "invert<byLU>( VectorAccessProxy )";

      blaze::CompressedVector< blaze::DynamicMatrix<double> > vec( 3UL, 1UL );
      vec[0].resize( 3UL, 3UL );
      vec[0] = 0.0;
      vec[0](0,0) = 1.0;
      vec[0](1,1) = 1.0;
      vec[0](2,2) = 1.0;
      invert<byLU>( vec[0] );

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkRows    ( vec[0], 3UL );
      checkColumns ( vec[0], 3UL );
      checkCapacity( vec[0], 9UL );
      checkNonZeros( vec[0], 3UL );
      checkRows    ( vec[1], 0UL );
      checkColumns ( vec[1], 0UL );
      checkRows    ( vec[2], 0UL );
      checkColumns ( vec[2], 0UL );
   }

   {
      test_ = "invert<byLLH>( VectorAccessProxy )";

      blaze::CompressedVector< blaze::DynamicMatrix<double> > vec( 3UL, 1UL );
      vec[0].resize( 3UL, 3UL );
      vec[0] = 0.0;
      vec[0](0,0) = 1.0;
      vec[0](1,1) = 1.0;
      vec[0](2,2) = 1.0;
      invert<byLLH>( vec[0] );

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 1UL );
      checkNonZeros( vec, 1UL );

      checkRows    ( vec[0], 3UL );
      checkColumns ( vec[0], 3UL );
      checkCapacity( vec[0], 9UL );
      checkNonZeros( vec[0], 3UL );
      checkRows    ( vec[1], 0UL );
      checkColumns ( vec[1], 0UL );
      checkRows    ( vec[2], 0UL );
      checkColumns ( vec[2], 0UL );
   }

#endif
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
   std::cout << "   Running CompressedVector proxy test..." << std::endl;

   try
   {
      RUN_COMPRESSEDVECTOR_PROXY_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during CompressedVector proxy test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
