//=================================================================================================
/*!
//  \file src/mathtest/matrices/compressedmatrix/ClassTest2.cpp
//  \brief Source file for the CompressedMatrix class test (part 2)
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
#include <blaze/util/Complex.h>
#include <blaze/util/Random.h>
#include <blazetest/mathtest/matrices/compressedmatrix/ClassTest.h>
#include <blazetest/mathtest/RandomMaximum.h>
#include <blazetest/mathtest/RandomMinimum.h>

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
/*!\brief Constructor for the CompressedMatrix class test.
//
// \exception std::runtime_error Operation error detected.
*/
ClassTest::ClassTest()
{
   testFunctionCall();
   testAt();
   testIterator();
   testNonZeros();
   testReset();
   testClear();
   testResize();
   testReserve();
   testTrim();
   testShrinkToFit();
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
   testIsDefault();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the CompressedMatrix function call operator.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of adding and accessing elements via the function call operator
// of the CompressedMatrix class template. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
void ClassTest::testFunctionCall()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major CompressedMatrix::operator()";

      // Assignment to the element (2,1)
      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 5UL, 3UL );
      mat(2,1) = 1;

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 5UL );
      checkCapacity( mat, 3UL );
      checkNonZeros( mat, 1UL );
      checkNonZeros( mat, 0UL, 0UL );
      checkNonZeros( mat, 1UL, 0UL );
      checkNonZeros( mat, 2UL, 1UL );

      if( mat(2,1) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 0 0 0 )\n( 0 0 0 0 0 )\n( 0 1 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Assignment to the element (1,4)
      mat(1,4) = 2;

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 5UL );
      checkCapacity( mat, 3UL );
      checkNonZeros( mat, 2UL );
      checkNonZeros( mat, 0UL, 0UL );
      checkNonZeros( mat, 1UL, 1UL );
      checkNonZeros( mat, 2UL, 1UL );

      if( mat(1,4) != 2 || mat(2,1) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 0 0 0 )\n( 0 0 0 0 2 )\n( 0 1 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Assignment to the element (0,3)
      mat(0,3) = 3;

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 5UL );
      checkCapacity( mat, 3UL );
      checkNonZeros( mat, 3UL );
      checkNonZeros( mat, 0UL, 1UL );
      checkNonZeros( mat, 1UL, 1UL );
      checkNonZeros( mat, 2UL, 1UL );

      if( mat(0,3) != 3 || mat(1,4) != 2 || mat(2,1) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 0 3 0 )\n( 0 0 0 0 2 )\n( 0 1 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Assignment to the element (2,2)
      mat(2,2) = 4;

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 5UL );
      checkCapacity( mat, 4UL );
      checkNonZeros( mat, 4UL );
      checkNonZeros( mat, 0UL, 1UL );
      checkNonZeros( mat, 1UL, 1UL );
      checkNonZeros( mat, 2UL, 2UL );

      if( mat(0,3) != 3 || mat(1,4) != 2 || mat(2,1) != 1 || mat(2,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 0 3 0 )\n( 0 0 0 0 2 )\n( 0 1 4 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Addition assignment to the element (2,1)
      mat(2,1) += mat(0,3);

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 5UL );
      checkCapacity( mat, 4UL );
      checkNonZeros( mat, 4UL );
      checkNonZeros( mat, 0UL, 1UL );
      checkNonZeros( mat, 1UL, 1UL );
      checkNonZeros( mat, 2UL, 2UL );

      if( mat(0,3) != 3 || mat(1,4) != 2 || mat(2,1) != 4 || mat(2,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 0 3 0 )\n( 0 0 0 0 2 )\n( 0 4 4 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Subtraction assignment to the element (1,0)
      mat(1,0) -= mat(1,4);

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 5UL );
      checkCapacity( mat, 4UL );
      checkNonZeros( mat, 5UL );
      checkNonZeros( mat, 0UL, 1UL );
      checkNonZeros( mat, 1UL, 2UL );
      checkNonZeros( mat, 2UL, 2UL );

      if( mat(0,3) != 3 || mat(1,0) != -2 || mat(1,4) != 2 || mat(2,1) != 4 || mat(2,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  0 0 0 3 0 )\n( -2 0 0 0 2 )\n(  0 4 4 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Multiplication assignment to the element (0,3)
      mat(0,3) *= -3;

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 5UL );
      checkCapacity( mat, 4UL );
      checkNonZeros( mat, 5UL );
      checkNonZeros( mat, 0UL, 1UL );
      checkNonZeros( mat, 1UL, 2UL );
      checkNonZeros( mat, 2UL, 2UL );

      if( mat(0,3) != -9 || mat(1,0) != -2 || mat(1,4) != 2 || mat(2,1) != 4 || mat(2,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  0 0 0 -3 0 )\n( -2 0 0  0 2 )\n(  0 4 4  0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Division assignment to the element (2,1)
      mat(2,1) /= 2;

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 5UL );
      checkCapacity( mat, 4UL );
      checkNonZeros( mat, 5UL );
      checkNonZeros( mat, 0UL, 1UL );
      checkNonZeros( mat, 1UL, 2UL );
      checkNonZeros( mat, 2UL, 2UL );

      if( mat(0,3) != -9 || mat(1,0) != -2 || mat(1,4) != 2 || mat(2,1) != 2 || mat(2,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  0 0 0 -3 0 )\n( -2 0 0  0 2 )\n(  0 2 4  0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major CompressedMatrix::operator()";

      // Assignment to the element (2,1)
      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 5UL, 3UL );
      mat(2,1) = 1;

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 5UL );
      checkCapacity( mat, 3UL );
      checkNonZeros( mat, 1UL );
      checkNonZeros( mat, 0UL, 0UL );
      checkNonZeros( mat, 1UL, 1UL );
      checkNonZeros( mat, 2UL, 0UL );
      checkNonZeros( mat, 3UL, 0UL );
      checkNonZeros( mat, 4UL, 0UL );

      if( mat(2,1) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 0 0 0 )\n( 0 0 0 0 0 )\n( 0 1 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Assignment to the element (1,4)
      mat(1,4) = 2;

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 5UL );
      checkCapacity( mat, 3UL );
      checkNonZeros( mat, 2UL );
      checkNonZeros( mat, 0UL, 0UL );
      checkNonZeros( mat, 1UL, 1UL );
      checkNonZeros( mat, 2UL, 0UL );
      checkNonZeros( mat, 3UL, 0UL );
      checkNonZeros( mat, 4UL, 1UL );

      if( mat(2,1) != 1 || mat(1,4) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 0 0 0 )\n( 0 0 0 0 2 )\n( 0 1 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Assignment to the element (0,3)
      mat(0,3) = 3;

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 5UL );
      checkCapacity( mat, 3UL );
      checkNonZeros( mat, 3UL );
      checkNonZeros( mat, 0UL, 0UL );
      checkNonZeros( mat, 1UL, 1UL );
      checkNonZeros( mat, 2UL, 0UL );
      checkNonZeros( mat, 3UL, 1UL );
      checkNonZeros( mat, 4UL, 1UL );

      if( mat(2,1) != 1 || mat(0,3) != 3 || mat(1,4) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 0 3 0 )\n( 0 0 0 0 2 )\n( 0 1 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Assignment to the element (2,2)
      mat(2,2) = 4;

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 5UL );
      checkCapacity( mat, 3UL );
      checkNonZeros( mat, 4UL );
      checkNonZeros( mat, 0UL, 0UL );
      checkNonZeros( mat, 1UL, 1UL );
      checkNonZeros( mat, 2UL, 1UL );
      checkNonZeros( mat, 3UL, 1UL );
      checkNonZeros( mat, 4UL, 1UL );

      if( mat(2,1) != 1 || mat(2,2) != 4 || mat(0,3) != 3 || mat(1,4) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 0 3 0 )\n( 0 0 0 0 2 )\n( 0 1 4 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Addition assignment to the element (2,1)
      mat(2,1) += mat(0,3);

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 5UL );
      checkCapacity( mat, 3UL );
      checkNonZeros( mat, 4UL );
      checkNonZeros( mat, 0UL, 0UL );
      checkNonZeros( mat, 1UL, 1UL );
      checkNonZeros( mat, 2UL, 1UL );
      checkNonZeros( mat, 3UL, 1UL );
      checkNonZeros( mat, 4UL, 1UL );

      if( mat(2,1) != 4 || mat(2,2) != 4 || mat(0,3) != 3 || mat(1,4) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 0 3 0 )\n( 0 0 0 0 2 )\n( 0 4 4 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Subtraction assignment to the element (1,0)
      mat(1,0) -= mat(1,4);

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 5UL );
      checkCapacity( mat, 3UL );
      checkNonZeros( mat, 5UL );
      checkNonZeros( mat, 0UL, 1UL );
      checkNonZeros( mat, 1UL, 1UL );
      checkNonZeros( mat, 2UL, 1UL );
      checkNonZeros( mat, 3UL, 1UL );
      checkNonZeros( mat, 4UL, 1UL );

      if( mat(1,0) != -2 || mat(2,1) != 4 || mat(2,2) != 4 || mat(0,3) != 3 || mat(1,4) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  0 0 0 3 0 )\n( -2 0 0 0 2 )\n(  0 4 4 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Multiplication assignment to the element (0,3)
      mat(0,3) *= -3;

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 5UL );
      checkCapacity( mat, 3UL );
      checkNonZeros( mat, 5UL );
      checkNonZeros( mat, 0UL, 1UL );
      checkNonZeros( mat, 1UL, 1UL );
      checkNonZeros( mat, 2UL, 1UL );
      checkNonZeros( mat, 3UL, 1UL );
      checkNonZeros( mat, 4UL, 1UL );

      if( mat(1,0) != -2 || mat(2,1) != 4 || mat(2,2) != 4 || mat(0,3) != -9 || mat(1,4) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  0 0 0 -9 0 )\n( -2 0 0  0 2 )\n(  0 4 4  0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Division assignment to the element (2,1)
      mat(2,1) /= 2;

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 5UL );
      checkCapacity( mat, 3UL );
      checkNonZeros( mat, 5UL );
      checkNonZeros( mat, 0UL, 1UL );
      checkNonZeros( mat, 1UL, 1UL );
      checkNonZeros( mat, 2UL, 1UL );
      checkNonZeros( mat, 3UL, 1UL );
      checkNonZeros( mat, 4UL, 1UL );

      if( mat(1,0) != -2 || mat(2,1) != 2 || mat(2,2) != 4 || mat(0,3) != -9 || mat(1,4) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  0 0 0 -9 0 )\n( -2 0 0  0 2 )\n(  0 2 4  0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c at() member function of the CompressedMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of adding and accessing elements via the \c at() member function
// of the CompressedMatrix class template. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
void ClassTest::testAt()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major CompressedMatrix::at()";

      // Assignment to the element (2,1)
      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 5UL, 3UL );
      mat.at(2,1) = 1;

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 5UL );
      checkCapacity( mat, 3UL );
      checkNonZeros( mat, 1UL );
      checkNonZeros( mat, 0UL, 0UL );
      checkNonZeros( mat, 1UL, 0UL );
      checkNonZeros( mat, 2UL, 1UL );

      if( mat.at(2,1) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Access via at() function failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 0 0 0 )\n( 0 0 0 0 0 )\n( 0 1 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Assignment to the element (1,4)
      mat.at(1,4) = 2;

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 5UL );
      checkCapacity( mat, 3UL );
      checkNonZeros( mat, 2UL );
      checkNonZeros( mat, 0UL, 0UL );
      checkNonZeros( mat, 1UL, 1UL );
      checkNonZeros( mat, 2UL, 1UL );

      if( mat.at(1,4) != 2 || mat.at(2,1) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Access via at() function failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 0 0 0 )\n( 0 0 0 0 2 )\n( 0 1 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Assignment to the element (0,3)
      mat.at(0,3) = 3;

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 5UL );
      checkCapacity( mat, 3UL );
      checkNonZeros( mat, 3UL );
      checkNonZeros( mat, 0UL, 1UL );
      checkNonZeros( mat, 1UL, 1UL );
      checkNonZeros( mat, 2UL, 1UL );

      if( mat.at(0,3) != 3 || mat.at(1,4) != 2 || mat.at(2,1) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Access via at() function failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 0 3 0 )\n( 0 0 0 0 2 )\n( 0 1 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Assignment to the element (2,2)
      mat.at(2,2) = 4;

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 5UL );
      checkCapacity( mat, 4UL );
      checkNonZeros( mat, 4UL );
      checkNonZeros( mat, 0UL, 1UL );
      checkNonZeros( mat, 1UL, 1UL );
      checkNonZeros( mat, 2UL, 2UL );

      if( mat.at(0,3) != 3 || mat.at(1,4) != 2 || mat.at(2,1) != 1 || mat.at(2,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Access via at() function failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 0 3 0 )\n( 0 0 0 0 2 )\n( 0 1 4 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Addition assignment to the element (2,1)
      mat.at(2,1) += mat.at(0,3);

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 5UL );
      checkCapacity( mat, 4UL );
      checkNonZeros( mat, 4UL );
      checkNonZeros( mat, 0UL, 1UL );
      checkNonZeros( mat, 1UL, 1UL );
      checkNonZeros( mat, 2UL, 2UL );

      if( mat.at(0,3) != 3 || mat.at(1,4) != 2 || mat.at(2,1) != 4 || mat.at(2,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Access via at() function failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 0 3 0 )\n( 0 0 0 0 2 )\n( 0 4 4 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Subtraction assignment to the element (1,0)
      mat.at(1,0) -= mat.at(1,4);

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 5UL );
      checkCapacity( mat, 4UL );
      checkNonZeros( mat, 5UL );
      checkNonZeros( mat, 0UL, 1UL );
      checkNonZeros( mat, 1UL, 2UL );
      checkNonZeros( mat, 2UL, 2UL );

      if( mat.at(0,3) != 3 || mat.at(1,0) != -2 || mat.at(1,4) != 2 || mat.at(2,1) != 4 || mat.at(2,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Access via at() function failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  0 0 0 3 0 )\n( -2 0 0 0 2 )\n(  0 4 4 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Multiplication assignment to the element (0,3)
      mat.at(0,3) *= -3;

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 5UL );
      checkCapacity( mat, 4UL );
      checkNonZeros( mat, 5UL );
      checkNonZeros( mat, 0UL, 1UL );
      checkNonZeros( mat, 1UL, 2UL );
      checkNonZeros( mat, 2UL, 2UL );

      if( mat.at(0,3) != -9 || mat.at(1,0) != -2 || mat.at(1,4) != 2 || mat.at(2,1) != 4 || mat.at(2,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Access via at() function failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  0 0 0 -3 0 )\n( -2 0 0  0 2 )\n(  0 4 4  0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Division assignment to the element (2,1)
      mat.at(2,1) /= 2;

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 5UL );
      checkCapacity( mat, 4UL );
      checkNonZeros( mat, 5UL );
      checkNonZeros( mat, 0UL, 1UL );
      checkNonZeros( mat, 1UL, 2UL );
      checkNonZeros( mat, 2UL, 2UL );

      if( mat.at(0,3) != -9 || mat.at(1,0) != -2 || mat.at(1,4) != 2 || mat.at(2,1) != 2 || mat.at(2,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Access via at() function failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  0 0 0 -3 0 )\n( -2 0 0  0 2 )\n(  0 2 4  0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Attempt to assign to the element (3,0)
      try {
         mat.at(3,0) = 2;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Out-of-bound access succeeded\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  0 0 0 -3 0 )\n( -2 0 0  0 2 )\n(  0 2 4  0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::out_of_range& ) {}

      // Attempt to assign to the element (0,5)
      try {
         mat.at(0,5) = 2;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Out-of-bound access succeeded\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  0 0 0 -3 0 )\n( -2 0 0  0 2 )\n(  0 2 4  0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::out_of_range& ) {}
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major CompressedMatrix::at()";

      // Assignment to the element (2,1)
      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 5UL, 3UL );
      mat.at(2,1) = 1;

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 5UL );
      checkCapacity( mat, 3UL );
      checkNonZeros( mat, 1UL );
      checkNonZeros( mat, 0UL, 0UL );
      checkNonZeros( mat, 1UL, 1UL );
      checkNonZeros( mat, 2UL, 0UL );
      checkNonZeros( mat, 3UL, 0UL );
      checkNonZeros( mat, 4UL, 0UL );

      if( mat.at(2,1) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Access via at() function failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 0 0 0 )\n( 0 0 0 0 0 )\n( 0 1 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Assignment to the element (1,4)
      mat.at(1,4) = 2;

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 5UL );
      checkCapacity( mat, 3UL );
      checkNonZeros( mat, 2UL );
      checkNonZeros( mat, 0UL, 0UL );
      checkNonZeros( mat, 1UL, 1UL );
      checkNonZeros( mat, 2UL, 0UL );
      checkNonZeros( mat, 3UL, 0UL );
      checkNonZeros( mat, 4UL, 1UL );

      if( mat.at(2,1) != 1 || mat.at(1,4) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Access via at() function failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 0 0 0 )\n( 0 0 0 0 2 )\n( 0 1 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Assignment to the element (0,3)
      mat.at(0,3) = 3;

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 5UL );
      checkCapacity( mat, 3UL );
      checkNonZeros( mat, 3UL );
      checkNonZeros( mat, 0UL, 0UL );
      checkNonZeros( mat, 1UL, 1UL );
      checkNonZeros( mat, 2UL, 0UL );
      checkNonZeros( mat, 3UL, 1UL );
      checkNonZeros( mat, 4UL, 1UL );

      if( mat.at(2,1) != 1 || mat.at(0,3) != 3 || mat.at(1,4) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Access via at() function failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 0 3 0 )\n( 0 0 0 0 2 )\n( 0 1 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Assignment to the element (2,2)
      mat.at(2,2) = 4;

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 5UL );
      checkCapacity( mat, 3UL );
      checkNonZeros( mat, 4UL );
      checkNonZeros( mat, 0UL, 0UL );
      checkNonZeros( mat, 1UL, 1UL );
      checkNonZeros( mat, 2UL, 1UL );
      checkNonZeros( mat, 3UL, 1UL );
      checkNonZeros( mat, 4UL, 1UL );

      if( mat.at(2,1) != 1 || mat.at(2,2) != 4 || mat.at(0,3) != 3 || mat.at(1,4) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Access via at() function failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 0 3 0 )\n( 0 0 0 0 2 )\n( 0 1 4 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Addition assignment to the element (2,1)
      mat.at(2,1) += mat.at(0,3);

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 5UL );
      checkCapacity( mat, 3UL );
      checkNonZeros( mat, 4UL );
      checkNonZeros( mat, 0UL, 0UL );
      checkNonZeros( mat, 1UL, 1UL );
      checkNonZeros( mat, 2UL, 1UL );
      checkNonZeros( mat, 3UL, 1UL );
      checkNonZeros( mat, 4UL, 1UL );

      if( mat.at(2,1) != 4 || mat.at(2,2) != 4 || mat.at(0,3) != 3 || mat.at(1,4) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Access via at() function failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 0 3 0 )\n( 0 0 0 0 2 )\n( 0 4 4 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Subtraction assignment to the element (1,0)
      mat.at(1,0) -= mat.at(1,4);

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 5UL );
      checkCapacity( mat, 3UL );
      checkNonZeros( mat, 5UL );
      checkNonZeros( mat, 0UL, 1UL );
      checkNonZeros( mat, 1UL, 1UL );
      checkNonZeros( mat, 2UL, 1UL );
      checkNonZeros( mat, 3UL, 1UL );
      checkNonZeros( mat, 4UL, 1UL );

      if( mat.at(1,0) != -2 || mat.at(2,1) != 4 || mat.at(2,2) != 4 || mat.at(0,3) != 3 || mat.at(1,4) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Access via at() function failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  0 0 0 3 0 )\n( -2 0 0 0 2 )\n(  0 4 4 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Multiplication assignment to the element (0,3)
      mat.at(0,3) *= -3;

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 5UL );
      checkCapacity( mat, 3UL );
      checkNonZeros( mat, 5UL );
      checkNonZeros( mat, 0UL, 1UL );
      checkNonZeros( mat, 1UL, 1UL );
      checkNonZeros( mat, 2UL, 1UL );
      checkNonZeros( mat, 3UL, 1UL );
      checkNonZeros( mat, 4UL, 1UL );

      if( mat.at(1,0) != -2 || mat.at(2,1) != 4 || mat.at(2,2) != 4 || mat.at(0,3) != -9 || mat.at(1,4) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Access via at() function failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  0 0 0 -9 0 )\n( -2 0 0  0 2 )\n(  0 4 4  0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Division assignment to the element (2,1)
      mat.at(2,1) /= 2;

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 5UL );
      checkCapacity( mat, 3UL );
      checkNonZeros( mat, 5UL );
      checkNonZeros( mat, 0UL, 1UL );
      checkNonZeros( mat, 1UL, 1UL );
      checkNonZeros( mat, 2UL, 1UL );
      checkNonZeros( mat, 3UL, 1UL );
      checkNonZeros( mat, 4UL, 1UL );

      if( mat.at(1,0) != -2 || mat.at(2,1) != 2 || mat.at(2,2) != 4 || mat.at(0,3) != -9 || mat.at(1,4) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Access via at() function failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  0 0 0 -9 0 )\n( -2 0 0  0 2 )\n(  0 2 4  0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Attempt to assign to the element (3,0)
      try {
         mat.at(3,0) = 2;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Out-of-bound access succeeded\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  0 0 0 -9 0 )\n( -2 0 0  0 2 )\n(  0 2 4  0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::out_of_range& ) {}

      // Attempt to assign to the element (0,5)
      try {
         mat.at(0,5) = 2;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Out-of-bound access succeeded\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  0 0 0 -9 0 )\n( -2 0 0  0 2 )\n(  0 2 4  0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::out_of_range& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the CompressedMatrix iterator implementation.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the iterator implementation of the CompressedMatrix class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testIterator()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      using MatrixType    = blaze::CompressedMatrix<int,blaze::rowMajor>;
      using Iterator      = MatrixType::Iterator;
      using ConstIterator = MatrixType::ConstIterator;

      MatrixType mat{ {  0, 1,  0 },
                      { -2, 0, -3 },
                      {  0, 4,  5 } };

      // Testing the Iterator default constructor
      {
         test_ = "Row-major Iterator default constructor";

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
         test_ = "Row-major ConstIterator default constructor";

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
         test_ = "Row-major Iterator/ConstIterator conversion";

         ConstIterator it( begin( mat, 1UL ) );

         if( it == end( mat, 1UL ) || it->value() != -2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Failed iterator conversion detected\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Counting the number of elements in 0th row via Iterator (end-begin)
      {
         test_ = "Row-major Iterator subtraction (end-begin)";

         const ptrdiff_t number( end( mat, 0UL ) - begin( mat, 0UL ) );

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

      // Counting the number of elements in 1st row via ConstIterator (end-begin)
      {
         test_ = "Row-major ConstIterator subtraction (end-begin)";

         const ptrdiff_t number( cend( mat, 1UL ) - cbegin( mat, 1UL ) );

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

         ConstIterator it ( cbegin( mat, 2UL ) );
         ConstIterator end( cend( mat, 2UL ) );

         if( it == end || it->value() != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid initial iterator detected\n";
            throw std::runtime_error( oss.str() );
         }

         ++it;

         if( it == end || it->value() != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator pre-increment failed\n";
            throw std::runtime_error( oss.str() );
         }

         it++;

         if( it != cend( mat, 2UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator post-increment failed\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing assignment via Iterator
      {
         test_ = "Row-major assignment via Iterator";

         int value = 8;

         for( Iterator it=begin( mat, 2UL ); it!=end( mat, 2UL ); ++it ) {
            *it = value++;
         }

         if( mat(0,0) !=  0 || mat(0,1) != 1 || mat(0,2) !=  0 ||
             mat(1,0) != -2 || mat(1,1) != 0 || mat(1,2) != -3 ||
             mat(2,0) !=  0 || mat(2,1) != 8 || mat(2,2) !=  9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n(  0  1  0 )\n( -2  0 -3 )\n(  0  8  9 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing addition assignment via Iterator
      {
         test_ = "Row-major addition assignment via Iterator";

         int value = 4;

         for( Iterator it=begin( mat, 1UL ); it!=end( mat, 1UL ); ++it ) {
            *it += value++;
         }

         if( mat(0,0) != 0 || mat(0,1) != 1 || mat(0,2) != 0 ||
             mat(1,0) != 2 || mat(1,1) != 0 || mat(1,2) != 2 ||
             mat(2,0) != 0 || mat(2,1) != 8 || mat(2,2) != 9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Addition assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 1 0 )\n( 2 0 2 )\n( 0 8 9 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing subtraction assignment via Iterator
      {
         test_ = "Row-major subtraction assignment via Iterator";

         int value = 4;

         for( Iterator it=begin( mat, 1UL ); it!=end( mat, 1UL ); ++it ) {
            *it -= value++;
         }

         if( mat(0,0) !=  0 || mat(0,1) != 1 || mat(0,2) !=  0 ||
             mat(1,0) != -2 || mat(1,1) != 0 || mat(1,2) != -3 ||
             mat(2,0) !=  0 || mat(2,1) != 8 || mat(2,2) !=  9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subtraction assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n(  0  1  0 )\n( -2  0 -3 )\n(  0  8  9 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing multiplication assignment via Iterator
      {
         test_ = "Row-major multiplication assignment via Iterator";

         int value = 1;

         for( Iterator it=begin( mat, 1UL ); it!=end( mat, 1UL ); ++it ) {
            *it *= value++;
         }

         if( mat(0,0) !=  0 || mat(0,1) != 1 || mat(0,2) !=  0 ||
             mat(1,0) != -2 || mat(1,1) != 0 || mat(1,2) != -6 ||
             mat(2,0) !=  0 || mat(2,1) != 8 || mat(2,2) !=  9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Multiplication assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n(  0  1  0 )\n( -2  0 -6 )\n(  0  8  9 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing division assignment via Iterator
      {
         test_ = "Row-major division assignment via Iterator";

         for( Iterator it=begin( mat, 1UL ); it!=end( mat, 1UL ); ++it ) {
            *it /= 2;
         }

         if( mat(0,0) !=  0 || mat(0,1) != 1 || mat(0,2) !=  0 ||
             mat(1,0) != -1 || mat(1,1) != 0 || mat(1,2) != -3 ||
             mat(2,0) !=  0 || mat(2,1) != 8 || mat(2,2) !=  9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Division assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n(  0  1  0 )\n( -1  0 -3 )\n(  0  8  9 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      using MatrixType    = blaze::CompressedMatrix<int,blaze::columnMajor>;
      using Iterator      = MatrixType::Iterator;
      using ConstIterator = MatrixType::ConstIterator;

      MatrixType mat{ { 0, -2, 0 },
                      { 1,  0, 4 },
                      { 0, -3, 5 } };

      // Testing the Iterator default constructor
      {
         test_ = "Column-major Iterator default constructor";

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
         test_ = "Column-major ConstIterator default constructor";

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
         test_ = "Column-major Iterator/ConstIterator conversion";

         ConstIterator it( begin( mat, 1UL ) );

         if( it == end( mat, 1UL ) || it->value() != -2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Failed iterator conversion detected\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Counting the number of elements in 0th column via Iterator (end-begin)
      {
         test_ = "Column-major Iterator subtraction (end-begin)";

         const ptrdiff_t number( end( mat, 0UL ) - begin( mat, 0UL ) );

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

      // Counting the number of elements in 1st row via ConstIterator (end-begin)
      {
         test_ = "Column-major ConstIterator subtraction (end-begin)";

         const ptrdiff_t number( cend( mat, 1UL ) - cbegin( mat, 1UL ) );

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

         ConstIterator it ( cbegin( mat, 2UL ) );
         ConstIterator end( cend( mat, 2UL ) );

         if( it == end || it->value() != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid initial iterator detected\n";
            throw std::runtime_error( oss.str() );
         }

         ++it;

         if( it == end || it->value() != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator pre-increment failed\n";
            throw std::runtime_error( oss.str() );
         }

         it++;

         if( it != cend( mat, 2UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator post-increment failed\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing assignment via Iterator
      {
         test_ = "Column-major assignment via Iterator";

         int value = 8;

         for( Iterator it=begin( mat, 2UL ); it!=end( mat, 2UL ); ++it ) {
            *it = value++;
         }

         if( mat(0,0) != 0 || mat(0,1) != -2 || mat(0,2) != 0 ||
             mat(1,0) != 1 || mat(1,1) !=  0 || mat(1,2) != 8 ||
             mat(2,0) != 0 || mat(2,1) != -3 || mat(2,2) != 9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 -2  0 )\n( 1  0  8 )\n( 0 -3  9 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing addition assignment via Iterator
      {
         test_ = "Column-major addition assignment via Iterator";

         int value = 4;

         for( Iterator it=begin( mat, 1UL ); it!=end( mat, 1UL ); ++it ) {
            *it += value++;
         }

         if( mat(0,0) != 0 || mat(0,1) != 2 || mat(0,2) != 0 ||
             mat(1,0) != 1 || mat(1,1) != 0 || mat(1,2) != 8 ||
             mat(2,0) != 0 || mat(2,1) != 2 || mat(2,2) != 9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Addition assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 2 0 )\n( 1 0 8 )\n( 0 2 9 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing subtraction assignment via Iterator
      {
         test_ = "Column-major subtraction assignment via Iterator";

         int value = 4;

         for( Iterator it=begin( mat, 1UL ); it!=end( mat, 1UL ); ++it ) {
            *it -= value++;
         }

         if( mat(0,0) != 0 || mat(0,1) != -2 || mat(0,2) != 0 ||
             mat(1,0) != 1 || mat(1,1) !=  0 || mat(1,2) != 8 ||
             mat(2,0) != 0 || mat(2,1) != -3 || mat(2,2) != 9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subtraction assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 -2  0 )\n( 1  0  8 )\n( 0 -3  9 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing multiplication assignment via Iterator
      {
         test_ = "Column-major multiplication assignment via Iterator";

         int value = 1;

         for( Iterator it=begin( mat, 1UL ); it!=end( mat, 1UL ); ++it ) {
            *it *= value++;
         }

         if( mat(0,0) != 0 || mat(0,1) != -2 || mat(0,2) != 0 ||
             mat(1,0) != 1 || mat(1,1) !=  0 || mat(1,2) != 8 ||
             mat(2,0) != 0 || mat(2,1) != -6 || mat(2,2) != 9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Multiplication assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 -2  0 )\n( 1  0  8 )\n( 0 -6  9 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing division assignment via Iterator
      {
         test_ = "Column-major division assignment via Iterator";

         for( Iterator it=begin( mat, 1UL ); it!=end( mat, 1UL ); ++it ) {
            *it /= 2;
         }

         if( mat(0,0) != 0 || mat(0,1) != -1 || mat(0,2) != 0 ||
             mat(1,0) != 1 || mat(1,1) !=  0 || mat(1,2) != 8 ||
             mat(2,0) != 0 || mat(2,1) != -3 || mat(2,2) != 9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Division assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 -1  0 )\n( 1  0  8 )\n( 0 -3  9 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c nonZeros() member function of the CompressedMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c nonZeros() member function of the CompressedMatrix
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testNonZeros()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major CompressedMatrix::nonZeros()";

      // Initial check
      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 6UL, 5UL, 2UL );

      checkRows    ( mat, 6UL );
      checkColumns ( mat, 5UL );
      checkCapacity( mat, 2UL );
      checkNonZeros( mat, 0UL );
      checkNonZeros( mat, 0UL, 0UL );
      checkNonZeros( mat, 1UL, 0UL );
      checkNonZeros( mat, 2UL, 0UL );
      checkNonZeros( mat, 3UL, 0UL );
      checkNonZeros( mat, 4UL, 0UL );
      checkNonZeros( mat, 5UL, 0UL );

      // Adding two non-zero elements
      mat(2,2) = 1;
      mat(4,0) = 2;

      checkRows    ( mat, 6UL );
      checkColumns ( mat, 5UL );
      checkCapacity( mat, 2UL );
      checkNonZeros( mat, 2UL );
      checkNonZeros( mat, 0UL, 0UL );
      checkNonZeros( mat, 1UL, 0UL );
      checkNonZeros( mat, 2UL, 1UL );
      checkNonZeros( mat, 3UL, 0UL );
      checkNonZeros( mat, 4UL, 1UL );
      checkNonZeros( mat, 5UL, 0UL );

      // Adding a third non-zero element
      mat(1,4) = 3;

      checkRows    ( mat, 6UL );
      checkColumns ( mat, 5UL );
      checkCapacity( mat, 3UL );
      checkNonZeros( mat, 3UL );
      checkNonZeros( mat, 0UL, 0UL );
      checkNonZeros( mat, 1UL, 1UL );
      checkNonZeros( mat, 2UL, 1UL );
      checkNonZeros( mat, 3UL, 0UL );
      checkNonZeros( mat, 4UL, 1UL );
      checkNonZeros( mat, 5UL, 0UL );
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major CompressedMatrix::nonZeros()";

      // Initial check
      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 6UL, 5UL, 2UL );

      checkRows    ( mat, 6UL );
      checkColumns ( mat, 5UL );
      checkCapacity( mat, 2UL );
      checkNonZeros( mat, 0UL );
      checkNonZeros( mat, 0UL, 0UL );
      checkNonZeros( mat, 1UL, 0UL );
      checkNonZeros( mat, 2UL, 0UL );
      checkNonZeros( mat, 3UL, 0UL );
      checkNonZeros( mat, 4UL, 0UL );

      // Adding two non-zero elements
      mat(2,2) = 1;
      mat(4,0) = 2;

      checkRows    ( mat, 6UL );
      checkColumns ( mat, 5UL );
      checkCapacity( mat, 2UL );
      checkNonZeros( mat, 2UL );
      checkNonZeros( mat, 0UL, 1UL );
      checkNonZeros( mat, 1UL, 0UL );
      checkNonZeros( mat, 2UL, 1UL );
      checkNonZeros( mat, 3UL, 0UL );
      checkNonZeros( mat, 4UL, 0UL );

      // Adding a third non-zero element
      mat(1,4) = 3;

      checkRows    ( mat, 6UL );
      checkColumns ( mat, 5UL );
      checkCapacity( mat, 3UL );
      checkNonZeros( mat, 3UL );
      checkNonZeros( mat, 0UL, 1UL );
      checkNonZeros( mat, 1UL, 0UL );
      checkNonZeros( mat, 2UL, 1UL );
      checkNonZeros( mat, 3UL, 0UL );
      checkNonZeros( mat, 4UL, 1UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c reset() member function of the CompressedMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c reset() member function of the CompressedMatrix
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testReset()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major CompressedMatrix::reset()";

      // Resetting a default constructed matrix
      {
         blaze::CompressedMatrix<int,blaze::rowMajor> mat;

         reset( mat );

         checkRows    ( mat, 0UL );
         checkColumns ( mat, 0UL );
         checkNonZeros( mat, 0UL );
      }

      // Resetting an initialized matrix
      {
         // Initialization check

         blaze::CompressedMatrix<int,blaze::rowMajor> mat{ { 1, 0, 0 },
                                                           { 0, 2, 3 },
                                                           { 0, 0, 0 },
                                                           { 0, 4, 5 } };

         checkRows    ( mat, 4UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 5UL );
         checkNonZeros( mat, 5UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 2UL );
         checkNonZeros( mat, 2UL, 0UL );
         checkNonZeros( mat, 3UL, 2UL );

         if( mat(0,0) != 1 || mat(1,1) != 2 || mat(1,2) != 3 || mat(3,1) != 4 || mat(3,2) != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 1 0 0 )\n( 0 2 3 )\n( 0 0 0 )\n( 0 4 5 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Resetting a single element
         reset( mat(3,1) );

         checkRows    ( mat, 4UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 5UL );
         checkNonZeros( mat, 4UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 2UL );
         checkNonZeros( mat, 2UL, 0UL );
         checkNonZeros( mat, 3UL, 1UL );

         if( mat(0,0) != 1 || mat(1,1) != 2 || mat(1,2) != 3 || mat(3,1) != 0 || mat(3,2) != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Reset operation failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 1 0 0 )\n( 0 2 3 )\n( 0 0 0 )\n( 0 0 5 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Resetting row 1
         reset( mat, 1UL );

         checkRows    ( mat, 4UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 5UL );
         checkNonZeros( mat, 2UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 0UL );
         checkNonZeros( mat, 3UL, 1UL );

         if( mat(0,0) != 1 || mat(1,1) != 0 || mat(1,2) != 0 || mat(3,1) != 0 || mat(3,2) != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Reset operation failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 0 )\n( 0 0 5 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Resetting the entire matrix
         reset( mat );

         checkRows    ( mat, 4UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 5UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 0UL );
         checkNonZeros( mat, 3UL, 0UL );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major CompressedMatrix::reset()";

      // Resetting a default constructed matrix
      {
         blaze::CompressedMatrix<int,blaze::columnMajor> mat;

         reset( mat );

         checkRows    ( mat, 0UL );
         checkColumns ( mat, 0UL );
         checkNonZeros( mat, 0UL );
      }

      // Resetting an initialized matrix
      {
         // Initialization check
         blaze::CompressedMatrix<int,blaze::columnMajor> mat{ { 1, 0, 0 },
                                                              { 0, 2, 3 },
                                                              { 0, 0, 0 },
                                                              { 0, 4, 5 } };

         checkRows    ( mat, 4UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 5UL );
         checkNonZeros( mat, 5UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 2UL );
         checkNonZeros( mat, 2UL, 2UL );

         if( mat(0,0) != 1 || mat(1,1) != 2 || mat(1,2) != 3 || mat(3,1) != 4 || mat(3,2) != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 1 0 0 )\n( 0 2 3 )\n( 0 0 0 )\n( 0 4 5 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Resetting a single element
         reset( mat(3,1) );

         checkRows    ( mat, 4UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 5UL );
         checkNonZeros( mat, 4UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 2UL );

         if( mat(0,0) != 1 || mat(1,1) != 2 || mat(1,2) != 3 || mat(3,1) != 0 || mat(3,2) != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Reset operation failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 1 0 0 )\n( 0 2 3 )\n( 0 0 0 )\n( 0 0 5 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Resetting column 1
         reset( mat, 1UL );

         checkRows    ( mat, 4UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 5UL );
         checkNonZeros( mat, 3UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 2UL );

         if( mat(0,0) != 1 || mat(1,1) != 0 || mat(1,2) != 3 || mat(3,1) != 0 || mat(3,2) != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Reset operation failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 1 0 0 )\n( 0 0 3 )\n( 0 0 0 )\n( 0 0 5 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Resetting the entire matrix
         reset( mat );

         checkRows    ( mat, 4UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 5UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 0UL );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c clear() member function of the CompressedMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c clear() member function of the CompressedMatrix
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testClear()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major CompressedMatrix::clear()";

      // Clearing a default constructed matrix
      {
         blaze::CompressedMatrix<int,blaze::rowMajor> mat;

         clear( mat );

         checkRows    ( mat, 0UL );
         checkColumns ( mat, 0UL );
         checkNonZeros( mat, 0UL );
      }

      // Clearing an initialized matrix
      {
         // Initialization check
         blaze::CompressedMatrix<int,blaze::rowMajor> mat{ { 1, 0, 0 },
                                                           { 0, 2, 3 },
                                                           { 0, 0, 0 },
                                                           { 0, 4, 0 } };

         checkRows    ( mat, 4UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 4UL );
         checkNonZeros( mat, 4UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 2UL );
         checkNonZeros( mat, 2UL, 0UL );
         checkNonZeros( mat, 3UL, 1UL );

         if( mat(0,0) != 1 || mat(1,1) != 2 || mat(1,2) != 3 || mat(3,1) != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 1 0 0 )\n( 0 2 3 )\n( 0 0 0 )\n( 0 4 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Clearing a single element
         clear( mat(1,2) );

         checkRows    ( mat, 4UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 4UL );
         checkNonZeros( mat, 3UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 0UL );
         checkNonZeros( mat, 3UL, 1UL );

         if( mat(0,0) != 1 || mat(1,1) != 2 || mat(1,2) != 0 || mat(3,1) != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Clear operation failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 0 )\n( 0 4 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Clearing the matrix
         clear( mat );

         checkRows    ( mat, 0UL );
         checkColumns ( mat, 0UL );
         checkNonZeros( mat, 0UL );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major CompressedMatrix::clear()";

      // Clearing a default constructed matrix
      {
         blaze::CompressedMatrix<int,blaze::columnMajor> mat;

         clear( mat );

         checkRows    ( mat, 0UL );
         checkColumns ( mat, 0UL );
         checkNonZeros( mat, 0UL );
      }

      // Clearing an initialized matrix
      {
         // Initialization check
         blaze::CompressedMatrix<int,blaze::columnMajor> mat{ { 1, 0, 0 },
                                                              { 0, 2, 3 },
                                                              { 0, 0, 0 },
                                                              { 0, 4, 0 } };

         checkRows    ( mat, 4UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 4UL );
         checkNonZeros( mat, 4UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 2UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( mat(0,0) != 1 || mat(1,1) != 2 || mat(1,2) != 3 || mat(3,1) != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Clear operation failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 1 0 0 )\n( 0 2 3 )\n( 0 0 0 )\n( 0 4 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Clearing a single element
         clear( mat(1,2) );

         checkRows    ( mat, 4UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 4UL );
         checkNonZeros( mat, 3UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 2UL );
         checkNonZeros( mat, 2UL, 0UL );

         if( mat(0,0) != 1 || mat(1,1) != 2 || mat(1,2) != 0 || mat(3,1) != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 0 )\n( 0 4 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Clearing the matrix
         clear( mat );

         checkRows    ( mat, 0UL );
         checkColumns ( mat, 0UL );
         checkNonZeros( mat, 0UL );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c resize() member function of the CompressedMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c resize() member function of the CompressedMatrix
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testResize()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major CompressedMatrix::resize()";

      // Initialization check
      blaze::CompressedMatrix<int,blaze::rowMajor> mat;

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );

      // Resizing to 0x3
      mat.resize( 0UL, 3UL );

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 3UL );
      checkNonZeros( mat, 0UL );

      // Resizing to 5x0
      mat.resize( 5UL, 0UL );

      checkRows    ( mat, 5UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );

      // Resizing to 3x4
      mat.resize( 3UL, 4UL );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 4UL );
      checkNonZeros( mat, 0UL );
      checkNonZeros( mat, 0UL, 0UL );
      checkNonZeros( mat, 1UL, 0UL );
      checkNonZeros( mat, 2UL, 0UL );

      // Resizing to 5x3 and preserving the elements
      mat(1,0) = 1;
      mat(2,2) = 2;
      mat.resize( 5UL, 3UL, true );

      checkRows    ( mat, 5UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 2UL );
      checkNonZeros( mat, 2UL );
      checkNonZeros( mat, 0UL, 0UL );
      checkNonZeros( mat, 1UL, 1UL );
      checkNonZeros( mat, 2UL, 1UL );
      checkNonZeros( mat, 3UL, 0UL );
      checkNonZeros( mat, 4UL, 0UL );

      if( mat(1,0) != 1 || mat(2,2) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 0 )\n( 1 0 0 )\n( 0 0 2 )\n( 0 0 0 )\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 4x4 and preserving the elements
      mat(0,1) = 3;
      mat.resize( 4UL, 4UL, true );

      checkRows    ( mat, 4UL );
      checkColumns ( mat, 4UL );
      checkCapacity( mat, 3UL );
      checkNonZeros( mat, 3UL );
      checkNonZeros( mat, 0UL, 1UL );
      checkNonZeros( mat, 1UL, 1UL );
      checkNonZeros( mat, 2UL, 1UL );
      checkNonZeros( mat, 3UL, 0UL );

      if( mat(1,0) != 1 || mat(2,2) != 2 || mat(0,1) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 3 0 0 )\n( 1 0 0 0 )\n( 0 0 2 0 )\n( 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 6x5 and preserving the elements
      mat(3,2) = 4;
      mat.resize( 6UL, 5UL, true );

      checkRows    ( mat, 6UL );
      checkColumns ( mat, 5UL );
      checkCapacity( mat, 4UL );
      checkNonZeros( mat, 4UL );
      checkNonZeros( mat, 0UL, 1UL );
      checkNonZeros( mat, 1UL, 1UL );
      checkNonZeros( mat, 2UL, 1UL );
      checkNonZeros( mat, 3UL, 1UL );
      checkNonZeros( mat, 4UL, 0UL );
      checkNonZeros( mat, 5UL, 0UL );

      if( mat(1,0) != 1 || mat(2,2) != 2 || mat(0,1) != 3 || mat(3,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 3 0 0 0 )\n( 1 0 0 0 0 )\n( 0 0 2 0 0 )\n"
                                     "( 0 0 4 0 0 )\n( 0 0 0 0 0 )\n( 0 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 4x3 and preserving the elements
      mat(0,4) = 5;
      mat(5,2) = 6;
      mat(5,4) = 7;
      mat.resize( 4UL, 3UL, true );

      checkRows    ( mat, 4UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 4UL );
      checkNonZeros( mat, 4UL );
      checkNonZeros( mat, 0UL, 1UL );
      checkNonZeros( mat, 1UL, 1UL );
      checkNonZeros( mat, 2UL, 1UL );
      checkNonZeros( mat, 3UL, 1UL );

      if( mat(1,0) != 1 || mat(2,2) != 2 || mat(0,1) != 3 || mat(3,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 3 0 )\n( 1 0 0 )\n( 0 0 2 )\n( 0 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 2x2
      mat.resize( 2UL, 2UL );

      checkRows   ( mat, 2UL );
      checkColumns( mat, 2UL );

      // Resizing to 0x0
      mat.resize( 0UL, 0UL );

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major CompressedMatrix::resize()";

      // Initialization check
      blaze::CompressedMatrix<int,blaze::columnMajor> mat;

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );

      // Resizing to 0x3
      mat.resize( 0UL, 3UL );

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 3UL );
      checkNonZeros( mat, 0UL );

      // Resizing to 5x0
      mat.resize( 5UL, 0UL );

      checkRows    ( mat, 5UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );

      // Resizing to 3x4
      mat.resize( 3UL, 4UL );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 4UL );
      checkNonZeros( mat, 0UL );
      checkNonZeros( mat, 0UL, 0UL );
      checkNonZeros( mat, 1UL, 0UL );
      checkNonZeros( mat, 2UL, 0UL );
      checkNonZeros( mat, 3UL, 0UL );

      // Resizing to 5x3 and preserving the elements
      mat(1,0) = 1;
      mat(2,2) = 2;
      mat.resize( 5UL, 3UL, true );

      checkRows    ( mat, 5UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 2UL );
      checkNonZeros( mat, 2UL );
      checkNonZeros( mat, 0UL, 1UL );
      checkNonZeros( mat, 1UL, 0UL );
      checkNonZeros( mat, 2UL, 1UL );

      if( mat(1,0) != 1 || mat(2,2) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 0 )\n( 1 0 0 )\n( 0 0 2 )\n( 0 0 0 )\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 4x4 and preserving the elements
      mat(0,1) = 3;
      mat.resize( 4UL, 4UL, true );

      checkRows    ( mat, 4UL );
      checkColumns ( mat, 4UL );
      checkCapacity( mat, 3UL );
      checkNonZeros( mat, 3UL );
      checkNonZeros( mat, 0UL, 1UL );
      checkNonZeros( mat, 1UL, 1UL );
      checkNonZeros( mat, 2UL, 1UL );
      checkNonZeros( mat, 3UL, 0UL );

      if( mat(1,0) != 1 || mat(2,2) != 2 || mat(0,1) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 3 0 0 )\n( 1 0 0 0 )\n( 0 0 2 0 )\n( 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 6x5 and preserving the elements
      mat(3,2) = 4;
      mat.resize( 6UL, 5UL, true );

      checkRows    ( mat, 6UL );
      checkColumns ( mat, 5UL );
      checkCapacity( mat, 4UL );
      checkNonZeros( mat, 4UL );
      checkNonZeros( mat, 0UL, 1UL );
      checkNonZeros( mat, 1UL, 1UL );
      checkNonZeros( mat, 2UL, 2UL );
      checkNonZeros( mat, 3UL, 0UL );
      checkNonZeros( mat, 4UL, 0UL );

      if( mat(1,0) != 1 || mat(2,2) != 2 || mat(0,1) != 3 || mat(3,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 3 0 0 0 )\n( 1 0 0 0 0 )\n( 0 0 2 0 0 )\n"
                                     "( 0 0 4 0 0 )\n( 0 0 0 0 0 )\n( 0 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 4x3 and preserving the elements
      mat(0,4) = 5;
      mat(5,2) = 6;
      mat(5,4) = 7;
      mat.resize( 4UL, 3UL, true );

      checkRows    ( mat, 4UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 4UL );
      checkNonZeros( mat, 4UL );
      checkNonZeros( mat, 0UL, 1UL );
      checkNonZeros( mat, 1UL, 1UL );
      checkNonZeros( mat, 2UL, 2UL );

      if( mat(1,0) != 1 || mat(2,2) != 2 || mat(0,1) != 3 || mat(3,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 3 0 )\n( 1 0 0 )\n( 0 0 2 )\n( 0 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 2x2
      mat.resize( 2UL, 2UL );

      checkRows   ( mat, 2UL );
      checkColumns( mat, 2UL );

      // Resizing to 0x0
      mat.resize( 0UL, 0UL );

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c reserve() member function of the CompressedMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c reserve() member function of the CompressedMatrix
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testReserve()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major CompressedMatrix::reserve()";

      // Initialization check
      blaze::CompressedMatrix<int,blaze::rowMajor> mat;

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );

      // Increasing the capacity of the matrix
      mat.reserve( 10UL );

      checkRows    ( mat,  0UL );
      checkColumns ( mat,  0UL );
      checkCapacity( mat, 10UL );
      checkNonZeros( mat,  0UL );

      // Further increasing the capacity of the matrix
      mat.reserve( 20UL );

      checkRows    ( mat,  0UL );
      checkColumns ( mat,  0UL );
      checkCapacity( mat, 20UL );
      checkNonZeros( mat,  0UL );
   }

   {
      test_ = "Row-major CompressedMatrix::reserve( size_t )";

      // Initialization check
      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 4UL );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 4UL );
      checkNonZeros( mat, 0UL );

      // Increasing the capacity of the 2nd row
      mat.reserve( 2UL, 10UL );

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  4UL );
      checkCapacity( mat, 10UL );
      checkCapacity( mat,  0UL,  0UL );
      checkCapacity( mat,  1UL,  0UL );
      checkCapacity( mat,  2UL, 10UL );

      // Increasing the capacity of the 0th row
      mat.reserve( 0UL, 20UL );

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  4UL );
      checkCapacity( mat, 30UL );
      checkCapacity( mat,  0UL, 20UL );
      checkCapacity( mat,  1UL,  0UL );
      checkCapacity( mat,  2UL, 10UL );

      // Increasing the capacity of the 1st row
      mat.reserve( 1UL, 15UL );

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  4UL );
      checkCapacity( mat, 45UL );
      checkCapacity( mat,  0UL, 20UL );
      checkCapacity( mat,  1UL, 15UL );
      checkCapacity( mat,  2UL, 10UL );
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major CompressedMatrix::reserve()";

      // Initialization check
      blaze::CompressedMatrix<int,blaze::columnMajor> mat;

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );

      // Increasing the capacity of the matrix
      mat.reserve( 10UL );

      checkRows    ( mat,  0UL );
      checkColumns ( mat,  0UL );
      checkCapacity( mat, 10UL );
      checkNonZeros( mat,  0UL );

      // Further increasing the capacity of the matrix
      mat.reserve( 20UL );

      checkRows    ( mat,  0UL );
      checkColumns ( mat,  0UL );
      checkCapacity( mat, 20UL );
      checkNonZeros( mat,  0UL );
   }

   {
      test_ = "Column-major CompressedMatrix::reserve( size_t )";

      // Initialization check
      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 4UL, 3UL );

      checkRows    ( mat, 4UL );
      checkColumns ( mat, 3UL );
      checkNonZeros( mat, 0UL );

      // Increasing the capacity of the 2nd column
      mat.reserve( 2UL, 10UL );

      checkRows    ( mat,  4UL );
      checkColumns ( mat,  3UL );
      checkCapacity( mat, 10UL );
      checkCapacity( mat,  0UL,  0UL );
      checkCapacity( mat,  1UL,  0UL );
      checkCapacity( mat,  2UL, 10UL );

      // Increasing the capacity of the 0th column
      mat.reserve( 0UL, 20UL );

      checkRows    ( mat,  4UL );
      checkColumns ( mat,  3UL );
      checkCapacity( mat, 30UL );
      checkCapacity( mat,  0UL, 20UL );
      checkCapacity( mat,  1UL,  0UL );
      checkCapacity( mat,  2UL, 10UL );

      // Increasing the capacity of the 1st column
      mat.reserve( 1UL, 15UL );

      checkRows    ( mat,  4UL );
      checkColumns ( mat,  3UL );
      checkCapacity( mat, 45UL );
      checkCapacity( mat,  0UL, 20UL );
      checkCapacity( mat,  1UL, 15UL );
      checkCapacity( mat,  2UL, 10UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c trim() member functions of the CompressedMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c trim() member functions of the CompressedMatrix
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testTrim()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major CompressedMatrix::trim()";

      // Initialization check
      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 4UL );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 4UL );
      checkNonZeros( mat, 0UL );

      // Increasing the row capacity of the matrix
      mat.reserve( 0UL, 10UL );
      mat.reserve( 1UL, 15UL );
      mat.reserve( 2UL, 20UL );

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  4UL );
      checkCapacity( mat, 45UL );
      checkCapacity( mat,  0UL, 10UL );
      checkCapacity( mat,  1UL, 15UL );
      checkCapacity( mat,  2UL, 20UL );

      // Trimming the matrix
      mat.trim();

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  4UL );
      checkCapacity( mat, 45UL );
      checkCapacity( mat,  0UL, 0UL );
      checkCapacity( mat,  1UL, 0UL );
      checkCapacity( mat,  2UL, 0UL );
   }

   {
      test_ = "Row-major CompressedMatrix::trim( size_t )";

      // Initialization check
      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 4UL );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 4UL );
      checkNonZeros( mat, 0UL );

      // Increasing the row capacity of the matrix
      mat.reserve( 0UL, 10UL );
      mat.reserve( 1UL, 15UL );
      mat.reserve( 2UL, 20UL );

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  4UL );
      checkCapacity( mat, 45UL );
      checkCapacity( mat,  0UL, 10UL );
      checkCapacity( mat,  1UL, 15UL );
      checkCapacity( mat,  2UL, 20UL );

      // Trimming the 0th row
      mat.trim( 0UL );

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  4UL );
      checkCapacity( mat, 45UL );
      checkCapacity( mat,  0UL,  0UL );
      checkCapacity( mat,  1UL, 25UL );
      checkCapacity( mat,  2UL, 20UL );

      // Trimming the 1st row
      mat.trim( 1UL );

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  4UL );
      checkCapacity( mat, 45UL );
      checkCapacity( mat,  0UL,  0UL );
      checkCapacity( mat,  1UL,  0UL );
      checkCapacity( mat,  2UL, 45UL );

      // Trimming the 2nd row
      mat.trim( 2UL );

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  4UL );
      checkCapacity( mat, 45UL );
      checkCapacity( mat,  0UL, 0UL );
      checkCapacity( mat,  1UL, 0UL );
      checkCapacity( mat,  2UL, 0UL );
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major CompressedMatrix::trim()";

      // Initialization check
      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 4UL, 3UL );

      checkRows    ( mat, 4UL );
      checkColumns ( mat, 3UL );
      checkNonZeros( mat, 0UL );

      // Increasing the column capacity of the matrix
      mat.reserve( 0UL, 10UL );
      mat.reserve( 1UL, 15UL );
      mat.reserve( 2UL, 20UL );

      checkRows    ( mat,  4UL );
      checkColumns ( mat,  3UL );
      checkCapacity( mat, 45UL );
      checkCapacity( mat,  0UL, 10UL );
      checkCapacity( mat,  1UL, 15UL );
      checkCapacity( mat,  2UL, 20UL );

      // Trimming the matrix
      mat.trim();

      checkRows    ( mat,  4UL );
      checkColumns ( mat,  3UL );
      checkCapacity( mat, 45UL );
      checkCapacity( mat,  0UL, 0UL );
      checkCapacity( mat,  1UL, 0UL );
      checkCapacity( mat,  2UL, 0UL );
   }

   {
      test_ = "Column-major CompressedMatrix::trim( size_t )";

      // Initialization check
      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 4UL, 3UL );

      checkRows    ( mat, 4UL );
      checkColumns ( mat, 3UL );
      checkNonZeros( mat, 0UL );

      // Increasing the column capacity of the matrix
      mat.reserve( 0UL, 10UL );
      mat.reserve( 1UL, 15UL );
      mat.reserve( 2UL, 20UL );

      checkRows    ( mat,  4UL );
      checkColumns ( mat,  3UL );
      checkCapacity( mat, 45UL );
      checkCapacity( mat,  0UL, 10UL );
      checkCapacity( mat,  1UL, 15UL );
      checkCapacity( mat,  2UL, 20UL );

      // Trimming the 0th column
      mat.trim( 0UL );

      checkRows    ( mat,  4UL );
      checkColumns ( mat,  3UL );
      checkCapacity( mat, 45UL );
      checkCapacity( mat,  0UL,  0UL );
      checkCapacity( mat,  1UL, 25UL );
      checkCapacity( mat,  2UL, 20UL );

      // Trimming the 1st column
      mat.trim( 1UL );

      checkRows    ( mat,  4UL );
      checkColumns ( mat,  3UL );
      checkCapacity( mat, 45UL );
      checkCapacity( mat,  0UL,  0UL );
      checkCapacity( mat,  1UL,  0UL );
      checkCapacity( mat,  2UL, 45UL );

      // Trimming the 2nd column
      mat.trim( 2UL );

      checkRows    ( mat,  4UL );
      checkColumns ( mat,  3UL );
      checkCapacity( mat, 45UL );
      checkCapacity( mat,  0UL, 0UL );
      checkCapacity( mat,  1UL, 0UL );
      checkCapacity( mat,  2UL, 0UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c shrinkToFit() member functions of the CompressedMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c shrinkToFit() member functions of the CompressedMatrix
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testShrinkToFit()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major DynamicMatrix::shrinkToFit()";

      // Shrinking a matrix without excessive capacity
      {
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 3UL, 3UL );
         mat(0,0) = 1;
         mat(0,2) = 3;
         mat(1,1) = 5;

         mat.shrinkToFit();

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 3UL );
         checkNonZeros( mat, 3UL );
         checkNonZeros( mat, 0UL, 2UL );
         checkNonZeros( mat, 1UL, 1UL );

         if( mat.capacity() != mat.nonZeros() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Shrinking the matrix failed\n"
                << " Details:\n"
                << "   Capacity         : " << mat.capacity() << "\n"
                << "   Expected capacity: " << mat.nonZeros() << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( mat(0,0) != 1 || mat(0,1) != 0 || mat(0,2) != 3 ||
             mat(1,0) != 0 || mat(1,1) != 5 || mat(1,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Shrinking the matrix failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 1 0 3 )\n( 0 5 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Shrinking a matrix with excessive capacity
      {
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 3UL, 100UL );
         mat(0,0) = 1;
         mat(0,2) = 3;
         mat(1,1) = 5;

         mat.shrinkToFit();

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 3UL );
         checkNonZeros( mat, 3UL );
         checkNonZeros( mat, 0UL, 2UL );
         checkNonZeros( mat, 1UL, 1UL );

         if( mat.capacity() != mat.nonZeros() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Shrinking the matrix failed\n"
                << " Details:\n"
                << "   Capacity         : " << mat.capacity() << "\n"
                << "   Expected capacity: " << mat.nonZeros() << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( mat(0,0) != 1 || mat(0,1) != 0 || mat(0,2) != 3 ||
             mat(1,0) != 0 || mat(1,1) != 5 || mat(1,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Shrinking the matrix failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 1 0 3 )\n( 0 5 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major DynamicMatrix::shrinkToFit()";

      // Shrinking a matrix without excessive capacity
      {
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 2UL, 3UL, 3UL );
         mat(0,0) = 1;
         mat(0,2) = 3;
         mat(1,1) = 5;

         mat.shrinkToFit();

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 3UL );
         checkNonZeros( mat, 3UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( mat.capacity() != mat.nonZeros() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Shrinking the matrix failed\n"
                << " Details:\n"
                << "   Capacity         : " << mat.capacity() << "\n"
                << "   Expected capacity: " << mat.nonZeros() << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( mat(0,0) != 1 || mat(0,1) != 0 || mat(0,2) != 3 ||
             mat(1,0) != 0 || mat(1,1) != 5 || mat(1,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Shrinking the matrix failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 1 0 3 )\n( 0 5 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Shrinking a matrix with excessive capacity
      {
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 2UL, 3UL, 100UL );
         mat(0,0) = 1;
         mat(0,2) = 3;
         mat(1,1) = 5;

         mat.shrinkToFit();

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 3UL );
         checkNonZeros( mat, 3UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( mat.capacity() != mat.nonZeros() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Shrinking the matrix failed\n"
                << " Details:\n"
                << "   Capacity         : " << mat.capacity() << "\n"
                << "   Expected capacity: " << mat.nonZeros() << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( mat(0,0) != 1 || mat(0,1) != 0 || mat(0,2) != 3 ||
             mat(1,0) != 0 || mat(1,1) != 5 || mat(1,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Shrinking the matrix failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 1 0 3 )\n( 0 5 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c swap() functionality of the CompressedMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c swap() function of the CompressedMatrix class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testSwap()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major CompressedMatrix swap";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1{ { 1, 0 },
                                                         { 0, 0 },
                                                         { 0, 0 },
                                                         { 0, 2 },
                                                         { 0, 0 } };

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2{ { 0, 3, 4, 0 },
                                                         { 0, 0, 0, 0 },
                                                         { 5, 0, 0, 0 } };

      swap( mat1, mat2 );

      checkRows    ( mat1, 3UL );
      checkColumns ( mat1, 4UL );
      checkCapacity( mat1, 3UL );
      checkNonZeros( mat1, 3UL );
      checkNonZeros( mat1, 0UL, 2UL );
      checkNonZeros( mat1, 1UL, 0UL );
      checkNonZeros( mat1, 2UL, 1UL );

      if( mat1(0,1) != 3 || mat1(0,2) != 4 || mat1(2,0) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the first matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n( 0 3 4 0 )\n( 0 0 0 0 )\n( 5 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      checkRows    ( mat2, 5UL );
      checkColumns ( mat2, 2UL );
      checkCapacity( mat2, 2UL );
      checkNonZeros( mat2, 2UL );
      checkNonZeros( mat2, 0UL, 1UL );
      checkNonZeros( mat2, 1UL, 0UL );
      checkNonZeros( mat2, 2UL, 0UL );
      checkNonZeros( mat2, 3UL, 1UL );
      checkNonZeros( mat2, 4UL, 0UL );

      if( mat2(0,0) != 1 || mat2(3,1) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the second matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 )\n( 0 0 )\n( 0 0 )\n( 0 2 )\n( 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major CompressedMatrix swap";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1{ { 1, 0 },
                                                            { 0, 0 },
                                                            { 0, 0 },
                                                            { 0, 2 },
                                                            { 0, 0 } };

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2{ { 0, 3, 4, 0 },
                                                            { 0, 0, 0, 0 },
                                                            { 5, 0, 0, 0 } };

      swap( mat1, mat2 );

      checkRows    ( mat1, 3UL );
      checkColumns ( mat1, 4UL );
      checkCapacity( mat1, 3UL );
      checkNonZeros( mat1, 3UL );
      checkNonZeros( mat1, 0UL, 1UL );
      checkNonZeros( mat1, 1UL, 1UL );
      checkNonZeros( mat1, 2UL, 1UL );
      checkNonZeros( mat1, 3UL, 0UL );

      if( mat1(0,1) != 3 || mat1(0,2) != 4 || mat1(2,0) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the first matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n( 0 3 4 0 )\n( 0 0 0 0 )\n( 5 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      checkRows    ( mat2, 5UL );
      checkColumns ( mat2, 2UL );
      checkCapacity( mat2, 2UL );
      checkNonZeros( mat2, 2UL );
      checkNonZeros( mat2, 0UL, 1UL );
      checkNonZeros( mat2, 1UL, 1UL );

      if( mat2(0,0) != 1 || mat2(3,1) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the second matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 )\n( 0 0 )\n( 0 0 )\n( 0 2 )\n( 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c set() member function of the CompressedMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c set() member function of the CompressedMatrix
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testSet()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major CompressedMatrix::set()";

      using Iterator = blaze::CompressedMatrix<int,blaze::rowMajor>::Iterator;

      // Initialization check
      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 4UL, 5UL );

      checkRows    ( mat, 4UL );
      checkColumns ( mat, 5UL );
      checkNonZeros( mat, 0UL );
      checkNonZeros( mat, 0UL, 0UL );
      checkNonZeros( mat, 1UL, 0UL );
      checkNonZeros( mat, 2UL, 0UL );
      checkNonZeros( mat, 3UL, 0UL );

      // Setting a non-zero element
      {
         Iterator pos = mat.set( 2UL, 3UL, 1 );

         checkRows    ( mat, 4UL );
         checkColumns ( mat, 5UL );
         checkCapacity( mat, 1UL );
         checkNonZeros( mat, 1UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 1UL );
         checkNonZeros( mat, 3UL, 0UL );

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

         if( mat(2,3) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 0 0 0 0 )\n( 0 0 0 0 0 )\n( 0 0 0 1 0 )\n( 0 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setting a second non-zero element
      {
         Iterator pos = mat.set( 2UL, 4UL, 2 );

         checkRows    ( mat, 4UL );
         checkColumns ( mat, 5UL );
         checkCapacity( mat, 2UL );
         checkNonZeros( mat, 2UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 2UL );
         checkNonZeros( mat, 3UL, 0UL );

         if( pos->value() != 2 || pos->index() != 4UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 2\n"
                << "   Expected index: 4\n";
            throw std::runtime_error( oss.str() );
         }

         if( mat(2,3) != 1 || mat(2,4) != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 0 0 0 0 )\n( 0 0 0 0 0 )\n( 0 0 0 1 2 )\n( 0 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setting a third non-zero element
      {
         Iterator pos = mat.set( 2UL, 2UL, 3 );

         checkRows    ( mat, 4UL );
         checkColumns ( mat, 5UL );
         checkCapacity( mat, 3UL );
         checkNonZeros( mat, 3UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 3UL );
         checkNonZeros( mat, 3UL, 0UL );

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

         if( mat(2,3) != 1 || mat(2,4) != 2 || mat(2,2) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 0 0 0 0 )\n( 0 0 0 0 0 )\n( 0 0 3 1 2 )\n( 0 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setting a fourth non-zero element
      {
         Iterator pos = mat.set( 0UL, 1UL, 4 );

         checkRows    ( mat, 4UL );
         checkColumns ( mat, 5UL );
         checkCapacity( mat, 4UL );
         checkNonZeros( mat, 4UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 3UL );
         checkNonZeros( mat, 3UL, 0UL );

         if( pos->value() != 4 || pos->index() != 1UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 4\n"
                << "   Expected index: 1\n";
            throw std::runtime_error( oss.str() );
         }

         if( mat(2,3) != 1 || mat(2,4) != 2 || mat(2,2) != 3 || mat(0,1) != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 4 0 0 0 )\n( 0 0 0 0 0 )\n( 0 0 3 1 2 )\n( 0 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setting a fifth non-zero element
      {
         Iterator pos = mat.set( 3UL, 2UL, 5 );

         checkRows    ( mat, 4UL );
         checkColumns ( mat, 5UL );
         checkCapacity( mat, 5UL );
         checkNonZeros( mat, 5UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 3UL );
         checkNonZeros( mat, 3UL, 1UL );

         if( pos->value() != 5 || pos->index() != 2UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 5\n"
                << "   Expected index: 2\n";
            throw std::runtime_error( oss.str() );
         }

         if( mat(2,3) != 1 || mat(2,4) != 2 || mat(2,2) != 3 || mat(0,1) != 4 || mat(3,2) != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 4 0 0 0 )\n( 0 0 0 0 0 )\n( 0 0 3 1 2 )\n( 0 0 5 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setting an already existing element
      {
         Iterator pos = mat.set( 3UL, 2UL, 6 );

         checkRows    ( mat, 4UL );
         checkColumns ( mat, 5UL );
         checkCapacity( mat, 5UL );
         checkNonZeros( mat, 5UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 3UL );
         checkNonZeros( mat, 3UL, 1UL );

         if( pos->value() != 6 || pos->index() != 2UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 5\n"
                << "   Expected index: 2\n";
            throw std::runtime_error( oss.str() );
         }

         if( mat(2,3) != 1 || mat(2,4) != 2 || mat(2,2) != 3 || mat(0,1) != 4 || mat(3,2) != 6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 4 0 0 0 )\n( 0 0 0 0 0 )\n( 0 0 3 1 2 )\n( 0 0 6 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major CompressedMatrix::set()";

      using Iterator = blaze::CompressedMatrix<int,blaze::columnMajor>::Iterator;

      // Initialization check
      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 5UL, 4UL );

      checkRows    ( mat, 5UL );
      checkColumns ( mat, 4UL );
      checkNonZeros( mat, 0UL );

      // Setting a non-zero element
      {
         Iterator pos = mat.set( 3UL, 2UL, 1 );

         checkRows    ( mat, 5UL );
         checkColumns ( mat, 4UL );
         checkCapacity( mat, 1UL );
         checkNonZeros( mat, 1UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 1UL );
         checkNonZeros( mat, 3UL, 0UL );

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

         if( mat(3,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 0 1 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setting a second non-zero element
      {
         Iterator pos = mat.set( 4UL, 2UL, 2 );

         checkRows    ( mat, 5UL );
         checkColumns ( mat, 4UL );
         checkCapacity( mat, 2UL );
         checkNonZeros( mat, 2UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 2UL );
         checkNonZeros( mat, 3UL, 0UL );

         if( pos->value() != 2 || pos->index() != 4UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 2\n"
                << "   Expected index: 4\n";
            throw std::runtime_error( oss.str() );
         }

         if( mat(3,2) != 1 || mat(4,2) != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 0 1 0 )\n( 0 0 2 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setting a third non-zero element
      {
         Iterator pos = mat.set( 2UL, 2UL, 3 );

         checkRows    ( mat, 5UL );
         checkColumns ( mat, 4UL );
         checkCapacity( mat, 3UL );
         checkNonZeros( mat, 3UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 3UL );
         checkNonZeros( mat, 3UL, 0UL );

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

         if( mat(3,2) != 1 || mat(4,2) != 2 || mat(2,2) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 0 3 0 )\n( 0 0 1 0 )\n( 0 0 2 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setting a fourth non-zero element
      {
         Iterator pos = mat.set( 1UL, 0UL, 4 );

         checkRows    ( mat, 5UL );
         checkColumns ( mat, 4UL );
         checkCapacity( mat, 4UL );
         checkNonZeros( mat, 4UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 3UL );
         checkNonZeros( mat, 3UL, 0UL );

         if( pos->value() != 4 || pos->index() != 1UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 4\n"
                << "   Expected index: 1\n";
            throw std::runtime_error( oss.str() );
         }

         if( mat(3,2) != 1 || mat(4,2) != 2 || mat(2,2) != 3 || mat(1,0) != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 4 0 0 0 )\n( 0 0 3 0 )\n( 0 0 1 0 )\n( 0 0 2 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setting a fifth non-zero element
      {
         Iterator pos = mat.set( 2UL, 3UL, 5 );

         checkRows    ( mat, 5UL );
         checkColumns ( mat, 4UL );
         checkCapacity( mat, 5UL );
         checkNonZeros( mat, 5UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 3UL );
         checkNonZeros( mat, 3UL, 1UL );

         if( pos->value() != 5 || pos->index() != 2UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 5\n"
                << "   Expected index: 2\n";
            throw std::runtime_error( oss.str() );
         }

         if( mat(3,2) != 1 || mat(4,2) != 2 || mat(2,2) != 3 || mat(1,0) != 4 || mat(2,3) != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 4 0 0 0 )\n( 0 0 3 5 )\n( 0 0 1 0 )\n( 0 0 2 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setting an already existing element
      {
         Iterator pos = mat.set( 2UL, 3UL, 6 );

         checkRows    ( mat, 5UL );
         checkColumns ( mat, 4UL );
         checkCapacity( mat, 5UL );
         checkNonZeros( mat, 5UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 3UL );
         checkNonZeros( mat, 3UL, 1UL );

         if( pos->value() != 6 || pos->index() != 2UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 6\n"
                << "   Expected index: 2\n";
            throw std::runtime_error( oss.str() );
         }

         if( mat(3,2) != 1 || mat(4,2) != 2 || mat(2,2) != 3 || mat(1,0) != 4 || mat(2,3) != 6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 4 0 0 0 )\n( 0 0 3 6 )\n( 0 0 1 0 )\n( 0 0 2 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c insert() member function of the CompressedMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c insert() member function of the CompressedMatrix
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testInsert()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major CompressedMatrix::insert()";

      using Iterator = blaze::CompressedMatrix<int,blaze::rowMajor>::Iterator;

      // Initialization check
      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 4UL, 5UL );

      checkRows    ( mat, 4UL );
      checkColumns ( mat, 5UL );
      checkNonZeros( mat, 0UL );
      checkNonZeros( mat, 0UL, 0UL );
      checkNonZeros( mat, 1UL, 0UL );
      checkNonZeros( mat, 2UL, 0UL );
      checkNonZeros( mat, 3UL, 0UL );

      // Inserting a non-zero element
      {
         Iterator pos = mat.insert( 2UL, 3UL, 1 );

         checkRows    ( mat, 4UL );
         checkColumns ( mat, 5UL );
         checkCapacity( mat, 1UL );
         checkNonZeros( mat, 1UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 1UL );
         checkNonZeros( mat, 3UL, 0UL );

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

         if( mat(2,3) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 0 0 0 0 )\n( 0 0 0 0 0 )\n( 0 0 0 1 0 )\n( 0 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting a second non-zero element
      {
         Iterator pos = mat.insert( 2UL, 4UL, 2 );

         checkRows    ( mat, 4UL );
         checkColumns ( mat, 5UL );
         checkCapacity( mat, 2UL );
         checkNonZeros( mat, 2UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 2UL );
         checkNonZeros( mat, 3UL, 0UL );

         if( pos->value() != 2 || pos->index() != 4UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 2\n"
                << "   Expected index: 4\n";
            throw std::runtime_error( oss.str() );
         }

         if( mat(2,3) != 1 || mat(2,4) != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 0 0 0 0 )\n( 0 0 0 0 0 )\n( 0 0 0 1 2 )\n( 0 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting a third non-zero element
      {
         Iterator pos = mat.insert( 2UL, 2UL, 3 );

         checkRows    ( mat, 4UL );
         checkColumns ( mat, 5UL );
         checkCapacity( mat, 3UL );
         checkNonZeros( mat, 3UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 3UL );
         checkNonZeros( mat, 3UL, 0UL );

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

         if( mat(2,3) != 1 || mat(2,4) != 2 || mat(2,2) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 0 0 0 0 )\n( 0 0 0 0 0 )\n( 0 0 3 1 2 )\n( 0 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting a fourth non-zero element
      {
         Iterator pos = mat.insert( 0UL, 1UL, 4 );

         checkRows    ( mat, 4UL );
         checkColumns ( mat, 5UL );
         checkCapacity( mat, 4UL );
         checkNonZeros( mat, 4UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 3UL );
         checkNonZeros( mat, 3UL, 0UL );

         if( pos->value() != 4 || pos->index() != 1UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 4\n"
                << "   Expected index: 1\n";
            throw std::runtime_error( oss.str() );
         }

         if( mat(2,3) != 1 || mat(2,4) != 2 || mat(2,2) != 3 || mat(0,1) != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 4 0 0 0 )\n( 0 0 0 0 0 )\n( 0 0 3 1 2 )\n( 0 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting a fifth non-zero element
      {
         Iterator pos = mat.insert( 3UL, 2UL, 5 );

         checkRows    ( mat, 4UL );
         checkColumns ( mat, 5UL );
         checkCapacity( mat, 5UL );
         checkNonZeros( mat, 5UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 3UL );
         checkNonZeros( mat, 3UL, 1UL );

         if( pos->value() != 5 || pos->index() != 2UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 5\n"
                << "   Expected index: 2\n";
            throw std::runtime_error( oss.str() );
         }

         if( mat(2,3) != 1 || mat(2,4) != 2 || mat(2,2) != 3 || mat(0,1) != 4 || mat(3,2) != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 4 0 0 0 )\n( 0 0 0 0 0 )\n( 0 0 3 1 2 )\n( 0 0 5 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to insert an already existing element
      try {
         mat.insert( 3UL, 2UL, 6 );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Inserting an existing element succeeded\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 4 0 0 0 )\n( 0 0 0 0 0 )\n( 0 0 3 1 2 )\n( 0 0 5 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major CompressedMatrix::insert()";

      using Iterator = blaze::CompressedMatrix<int,blaze::columnMajor>::Iterator;

      // Initialization check
      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 5UL, 4UL );

      checkRows    ( mat, 5UL );
      checkColumns ( mat, 4UL );
      checkNonZeros( mat, 0UL );

      // Inserting a non-zero element
      {
         Iterator pos = mat.insert( 3UL, 2UL, 1 );

         checkRows    ( mat, 5UL );
         checkColumns ( mat, 4UL );
         checkCapacity( mat, 1UL );
         checkNonZeros( mat, 1UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 1UL );
         checkNonZeros( mat, 3UL, 0UL );

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

         if( mat(3,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 0 1 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting a second non-zero element
      {
         Iterator pos = mat.insert( 4UL, 2UL, 2 );

         checkRows    ( mat, 5UL );
         checkColumns ( mat, 4UL );
         checkCapacity( mat, 2UL );
         checkNonZeros( mat, 2UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 2UL );
         checkNonZeros( mat, 3UL, 0UL );

         if( pos->value() != 2 || pos->index() != 4UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 2\n"
                << "   Expected index: 4\n";
            throw std::runtime_error( oss.str() );
         }

         if( mat(3,2) != 1 || mat(4,2) != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 0 1 0 )\n( 0 0 2 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting a third non-zero element
      {
         Iterator pos = mat.insert( 2UL, 2UL, 3 );

         checkRows    ( mat, 5UL );
         checkColumns ( mat, 4UL );
         checkCapacity( mat, 3UL );
         checkNonZeros( mat, 3UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 3UL );
         checkNonZeros( mat, 3UL, 0UL );

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

         if( mat(3,2) != 1 || mat(4,2) != 2 || mat(2,2) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 0 3 0 )\n( 0 0 1 0 )\n( 0 0 2 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting a fourth non-zero element
      {
         Iterator pos = mat.insert( 1UL, 0UL, 4 );

         checkRows    ( mat, 5UL );
         checkColumns ( mat, 4UL );
         checkCapacity( mat, 4UL );
         checkNonZeros( mat, 4UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 3UL );
         checkNonZeros( mat, 3UL, 0UL );

         if( pos->value() != 4 || pos->index() != 1UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 4\n"
                << "   Expected index: 1\n";
            throw std::runtime_error( oss.str() );
         }

         if( mat(3,2) != 1 || mat(4,2) != 2 || mat(2,2) != 3 || mat(1,0) != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 4 0 0 0 )\n( 0 0 3 0 )\n( 0 0 1 0 )\n( 0 0 2 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting a fifth non-zero element
      {
         Iterator pos = mat.insert( 2UL, 3UL, 5 );

         checkRows    ( mat, 5UL );
         checkColumns ( mat, 4UL );
         checkCapacity( mat, 5UL );
         checkNonZeros( mat, 5UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 3UL );
         checkNonZeros( mat, 3UL, 1UL );

         if( pos->value() != 5 || pos->index() != 2UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 5\n"
                << "   Expected index: 2\n";
            throw std::runtime_error( oss.str() );
         }

         if( mat(3,2) != 1 || mat(4,2) != 2 || mat(2,2) != 3 || mat(1,0) != 4 || mat(2,3) != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 4 0 0 0 )\n( 0 0 3 5 )\n( 0 0 1 0 )\n( 0 0 2 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to insert an already existing element
      try {
         mat.insert( 2UL, 3UL, 6 );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Inserting an existing element succeeded\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 0 0 )\n( 4 0 0 0 )\n( 0 0 3 5 )\n( 0 0 1 0 )\n( 0 0 2 0 )\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c append() member function of the CompressedMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c append() member function of the CompressedMatrix
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testAppend()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major CompressedMatrix::append()";

      // Appending with pre-allocation in each row
      {
         // Initialization check
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 4UL, 4UL, 5UL );
         mat.reserve( 0UL, 2UL );
         mat.reserve( 2UL, 1UL );
         mat.reserve( 3UL, 2UL );

         checkRows    ( mat, 4UL );
         checkColumns ( mat, 4UL );
         checkCapacity( mat, 5UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 0UL );
         checkNonZeros( mat, 3UL, 0UL );

         // Appending one non-zero element
         mat.append( 2UL, 1UL, 1 );

         checkRows    ( mat, 4UL );
         checkColumns ( mat, 4UL );
         checkCapacity( mat, 5UL );
         checkNonZeros( mat, 1UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 1UL );
         checkNonZeros( mat, 3UL, 0UL );

         if( mat(2,1) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 1 0 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Appending two more non-zero elements
         mat.append( 0UL, 0UL, 2 );
         mat.append( 0UL, 3UL, 3 );

         checkRows    ( mat, 4UL );
         checkColumns ( mat, 4UL );
         checkCapacity( mat, 5UL );
         checkNonZeros( mat, 3UL );
         checkNonZeros( mat, 0UL, 2UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 1UL );
         checkNonZeros( mat, 3UL, 0UL );

         if( mat(2,1) != 1 || mat(0,0) != 2 || mat(0,3) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 2 0 0 3 )\n( 0 0 0 0 )\n( 0 1 0 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Appending two more non-zero elements
         mat.append( 3UL, 1UL, 4 );
         mat.append( 3UL, 2UL, 5 );

         checkRows    ( mat, 4UL );
         checkColumns ( mat, 4UL );
         checkCapacity( mat, 5UL );
         checkNonZeros( mat, 5UL );
         checkNonZeros( mat, 0UL, 2UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 1UL );
         checkNonZeros( mat, 3UL, 2UL );

         if( mat(2,1) != 1 || mat(0,0) != 2 || mat(0,3) != 3 || mat(3,1) != 4 || mat(3,2) != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 2 0 0 3 )\n( 0 0 0 0 )\n( 0 1 0 0 )\n( 0 4 5 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Appending with row finalization
      {
         // Initialization check
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 4UL, 4UL, 5UL );
         mat.reserve( 0UL, 2UL );
         mat.reserve( 2UL, 1UL );
         mat.reserve( 3UL, 2UL );

         // Appending one non-zero element
         mat.append( 0UL, 1UL, 1 );
         mat.finalize( 0UL );

         checkRows    ( mat, 4UL );
         checkColumns ( mat, 4UL );
         checkCapacity( mat, 5UL );
         checkNonZeros( mat, 1UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 0UL );
         checkNonZeros( mat, 3UL, 0UL );

         if( mat(0,1) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 1 0 0 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Appending two more non-zero elements
         mat.append( 1UL, 1UL, 2 );
         mat.append( 1UL, 3UL, 3 );
         mat.finalize( 1UL );

         checkRows    ( mat, 4UL );
         checkColumns ( mat, 4UL );
         checkCapacity( mat, 5UL );
         checkNonZeros( mat, 3UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 2UL );
         checkNonZeros( mat, 2UL, 0UL );
         checkNonZeros( mat, 3UL, 0UL );

         if( mat(0,1) != 1 || mat(1,1) != 2 || mat(1,3) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 1 0 0 )\n( 0 2 0 3 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Appending two more non-zero elements
         mat.append( 3UL, 0UL, 4 );
         mat.append( 3UL, 1UL, 5 );
         mat.finalize( 3UL );

         checkRows    ( mat, 4UL );
         checkColumns ( mat, 4UL );
         checkCapacity( mat, 5UL );
         checkNonZeros( mat, 5UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 2UL );
         checkNonZeros( mat, 2UL, 0UL );
         checkNonZeros( mat, 3UL, 2UL );

         if( mat(0,1) != 1 || mat(1,1) != 2 || mat(1,3) != 3 || mat(3,0) != 4 || mat(3,1) != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 1 0 0 )\n( 0 2 0 3 )\n( 0 0 0 0 )\n( 4 5 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major CompressedMatrix::append()";

      // Appending with pre-allocation in each column
      {
         // Initialization check
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 4UL, 4UL, 5UL );
         mat.reserve( 0UL, 2UL );
         mat.reserve( 2UL, 1UL );
         mat.reserve( 3UL, 2UL );

         checkRows    ( mat, 4UL );
         checkColumns ( mat, 4UL );
         checkCapacity( mat, 5UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 0UL );
         checkNonZeros( mat, 3UL, 0UL );

         // Appending one non-zero element
         mat.append( 1UL, 2UL, 1 );

         checkRows    ( mat, 4UL );
         checkColumns ( mat, 4UL );
         checkCapacity( mat, 5UL );
         checkNonZeros( mat, 1UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 1UL );
         checkNonZeros( mat, 3UL, 0UL );

         if( mat(1,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 1 0 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Appending two more non-zero elements
         mat.append( 0UL, 0UL, 2 );
         mat.append( 3UL, 0UL, 3 );

         checkRows    ( mat, 4UL );
         checkColumns ( mat, 4UL );
         checkCapacity( mat, 5UL );
         checkNonZeros( mat, 3UL );
         checkNonZeros( mat, 0UL, 2UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 1UL );
         checkNonZeros( mat, 3UL, 0UL );

         if( mat(1,2) != 1 || mat(0,0) != 2 || mat(3,0) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 2 0 0 0 )\n( 0 0 1 0 )\n( 0 0 0 0 )\n( 3 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Appending two more non-zero elements
         mat.append( 1UL, 3UL, 4 );
         mat.append( 2UL, 3UL, 5 );

         checkRows    ( mat, 4UL );
         checkColumns ( mat, 4UL );
         checkCapacity( mat, 5UL );
         checkNonZeros( mat, 5UL );
         checkNonZeros( mat, 0UL, 2UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 1UL );
         checkNonZeros( mat, 3UL, 2UL );

         if( mat(1,2) != 1 || mat(0,0) != 2 || mat(3,0) != 3 || mat(1,3) != 4 || mat(2,3) != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 2 0 0 0 )\n( 0 0 1 4 )\n( 0 0 0 5 )\n( 3 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Appending with column finalization
      {
         // Initialization check
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 4UL, 4UL, 5UL );
         mat.reserve( 0UL, 2UL );
         mat.reserve( 2UL, 1UL );
         mat.reserve( 3UL, 2UL );

         // Appending one non-zero element
         mat.append( 1UL, 0UL, 1 );
         mat.finalize( 0UL );

         checkRows    ( mat, 4UL );
         checkColumns ( mat, 4UL );
         checkCapacity( mat, 5UL );
         checkNonZeros( mat, 1UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 0UL );
         checkNonZeros( mat, 3UL, 0UL );

         if( mat(1,0) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 1 0 0 0 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Appending two more non-zero elements
         mat.append( 1UL, 1UL, 2 );
         mat.append( 3UL, 1UL, 3 );
         mat.finalize( 1UL );

         checkRows    ( mat, 4UL );
         checkColumns ( mat, 4UL );
         checkCapacity( mat, 5UL );
         checkNonZeros( mat, 3UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 2UL );
         checkNonZeros( mat, 2UL, 0UL );
         checkNonZeros( mat, 3UL, 0UL );

         if( mat(1,0) != 1 || mat(1,1) != 2 || mat(3,1) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 1 2 0 0 )\n( 0 0 0 0 )\n( 0 3 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Appending two more non-zero elements
         mat.append( 0UL, 3UL, 4 );
         mat.append( 1UL, 3UL, 5 );
         mat.finalize( 3UL );

         checkRows    ( mat, 4UL );
         checkColumns ( mat, 4UL );
         checkCapacity( mat, 5UL );
         checkNonZeros( mat, 5UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 2UL );
         checkNonZeros( mat, 2UL, 0UL );
         checkNonZeros( mat, 3UL, 2UL );

         if( mat(1,0) != 1 || mat(1,1) != 2 || mat(3,1) != 3 || mat(0,3) != 4 || mat(1,3) != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 0 0 4 )\n( 1 2 0 5 )\n( 0 0 0 0 )\n( 0 3 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c erase() member function of the CompressedMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c erase() member function of the CompressedMatrix
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testErase()
{
   //=====================================================================================
   // Row-major index-based erase function
   //=====================================================================================

   {
      test_ = "Row-major CompressedMatrix::erase( size_t, size_t )";

      // Initialization check
      blaze::CompressedMatrix<int,blaze::rowMajor> mat{ { 1, 0, 2, 0, 0 },
                                                        { 0, 3, 4, 0, 5 },
                                                        { 0, 6, 0, 0, 7 } };

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 5UL );
      checkCapacity( mat, 7UL );
      checkNonZeros( mat, 7UL );
      checkNonZeros( mat, 0UL, 2UL );
      checkNonZeros( mat, 1UL, 3UL );
      checkNonZeros( mat, 2UL, 2UL );

      if( mat(0,0) != 1 || mat(0,2) != 2 ||
          mat(1,1) != 3 || mat(1,2) != 4 || mat(1,4) != 5 ||
          mat(2,1) != 6 || mat(2,4) != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 0 2 0 0 )\n( 0 3 4 0 5 )\n( 0 6 0 0 7 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element at (0,0)
      mat.erase( 0UL, size_t(0) );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 5UL );
      checkCapacity( mat, 7UL );
      checkNonZeros( mat, 6UL );
      checkNonZeros( mat, 0UL, 1UL );
      checkNonZeros( mat, 1UL, 3UL );
      checkNonZeros( mat, 2UL, 2UL );

      if( mat(0,2) != 2 ||
          mat(1,1) != 3 || mat(1,2) != 4 || mat(1,4) != 5 ||
          mat(2,1) != 6 || mat(2,4) != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 2 0 0 )\n( 0 3 4 0 5 )\n( 0 6 0 0 7 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element at (1,2)
      mat.erase( 1UL, 2UL );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 5UL );
      checkCapacity( mat, 7UL );
      checkNonZeros( mat, 5UL );
      checkNonZeros( mat, 0UL, 1UL );
      checkNonZeros( mat, 1UL, 2UL );
      checkNonZeros( mat, 2UL, 2UL );

      if( mat(0,2) != 2 ||
          mat(1,1) != 3 || mat(1,4) != 5 ||
          mat(2,1) != 6 || mat(2,4) != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 2 0 0 )\n( 0 3 0 0 5 )\n( 0 6 0 0 7 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element at (2,4)
      mat.erase( 2UL, 4UL );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 5UL );
      checkCapacity( mat, 7UL );
      checkNonZeros( mat, 4UL );
      checkNonZeros( mat, 0UL, 1UL );
      checkNonZeros( mat, 1UL, 2UL );
      checkNonZeros( mat, 2UL, 1UL );

      if( mat(0,2) != 2 ||
          mat(1,1) != 3 || mat(1,4) != 5 ||
          mat(2,1) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 2 0 0 )\n( 0 3 0 0 5 )\n( 0 6 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to erase a zero element
      mat.erase( 0UL, 1UL );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 5UL );
      checkCapacity( mat, 7UL );
      checkNonZeros( mat, 4UL );
      checkNonZeros( mat, 0UL, 1UL );
      checkNonZeros( mat, 1UL, 2UL );
      checkNonZeros( mat, 2UL, 1UL );

      if( mat(0,2) != 2 ||
          mat(1,1) != 3 || mat(1,4) != 5 ||
          mat(2,1) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 2 0 0 )\n( 0 3 0 0 5 )\n( 0 6 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major iterator-based erase function
   //=====================================================================================

   {
      test_ = "Row-major CompressedMatrix::erase( size_t, Iterator )";

      using MatrixType = blaze::CompressedMatrix<int,blaze::rowMajor>;
      using Iterator   = MatrixType::Iterator;

      // Initialization check
      MatrixType mat{ { 1, 0, 2, 0, 0 },
                      { 0, 3, 4, 0, 5 },
                      { 0, 6, 0, 0, 7 } };

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 5UL );
      checkCapacity( mat, 7UL );
      checkNonZeros( mat, 7UL );
      checkNonZeros( mat, 0UL, 2UL );
      checkNonZeros( mat, 1UL, 3UL );
      checkNonZeros( mat, 2UL, 2UL );

      if( mat(0,0) != 1 || mat(0,2) != 2 ||
          mat(1,1) != 3 || mat(1,2) != 4 || mat(1,4) != 5 ||
          mat(2,1) != 6 || mat(2,4) != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 0 2 0 0 )\n( 0 3 4 0 5 )\n( 0 6 0 0 7 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element at (0,0)
      {
         Iterator pos = mat.erase( 0UL, mat.find( 0UL, 0UL ) );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 5UL );
         checkCapacity( mat, 7UL );
         checkNonZeros( mat, 6UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 3UL );
         checkNonZeros( mat, 2UL, 2UL );

         if( mat(0,2) != 2 ||
             mat(1,1) != 3 || mat(1,2) != 4 || mat(1,4) != 5 ||
             mat(2,1) != 6 || mat(2,4) != 7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 0 2 0 0 )\n( 0 3 4 0 5 )\n( 0 6 0 0 7 )\n";
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

      // Erasing the element at (1,2)
      {
         Iterator pos = mat.erase( 1UL, mat.find( 1UL, 2UL ) );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 5UL );
         checkCapacity( mat, 7UL );
         checkNonZeros( mat, 5UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 2UL );
         checkNonZeros( mat, 2UL, 2UL );

         if( mat(0,2) != 2 ||
             mat(1,1) != 3 || mat(1,4) != 5 ||
             mat(2,1) != 6 || mat(2,4) != 7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 0 2 0 0 )\n( 0 3 0 0 5 )\n( 0 6 0 0 7 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos->value() != 5 || pos->index() != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 5\n"
                << "   Expected index: 4\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the element at (2,4)
      {
         Iterator pos = mat.erase( 2UL, mat.find( 2UL, 4UL ) );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 5UL );
         checkCapacity( mat, 7UL );
         checkNonZeros( mat, 4UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 2UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( mat(0,2) != 2 ||
             mat(1,1) != 3 || mat(1,4) != 5 ||
             mat(2,1) != 6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 0 2 0 0 )\n( 0 3 0 0 5 )\n( 0 6 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos != mat.end( 2UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to erase a zero element
      {
         Iterator pos = mat.erase( 0UL, mat.find( 0UL, 1UL ) );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 5UL );
         checkCapacity( mat, 7UL );
         checkNonZeros( mat, 4UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 2UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( mat(0,2) != 2 ||
             mat(1,1) != 3 || mat(1,4) != 5 ||
             mat(2,1) != 6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 0 2 0 0 )\n( 0 3 0 0 5 )\n( 0 6 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos != mat.end( 0UL ) ) {
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
   // Row-major iterator-range-based erase function
   //=====================================================================================

   {
      test_ = "Row-major CompressedMatrix::erase( size_t, Iterator, Iterator )";

      using MatrixType = blaze::CompressedMatrix<int,blaze::rowMajor>;
      using Iterator   = MatrixType::Iterator;

      // Initialization check
      MatrixType mat{ { 1, 0, 2, 0, 0 },
                      { 0, 3, 4, 0, 5 },
                      { 0, 6, 0, 0, 7 } };

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 5UL );
      checkCapacity( mat, 7UL );
      checkNonZeros( mat, 7UL );
      checkNonZeros( mat, 0UL, 2UL );
      checkNonZeros( mat, 1UL, 3UL );
      checkNonZeros( mat, 2UL, 2UL );

      if( mat(0,0) != 1 || mat(0,2) != 2 ||
          mat(1,1) != 3 || mat(1,2) != 4 || mat(1,4) != 5 ||
          mat(2,1) != 6 || mat(2,4) != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 0 2 0 0 )\n( 0 3 4 0 5 )\n( 0 6 0 0 7 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the elements from (0,0) to (0,2)
      {
         Iterator pos = mat.erase( 0UL, mat.find( 0UL, 0UL ), mat.find( 0UL, 2UL ) );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 5UL );
         checkCapacity( mat, 7UL );
         checkNonZeros( mat, 6UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 3UL );
         checkNonZeros( mat, 2UL, 2UL );

         if( mat(0,2) != 2 ||
             mat(1,1) != 3 || mat(1,2) != 4 || mat(1,4) != 5 ||
             mat(2,1) != 6 || mat(2,4) != 7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a single-element range failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 0 2 0 0 )\n( 0 3 4 0 5 )\n( 0 6 0 0 7 )\n";
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

      // Erasing the elements from (1,2) to (1,4)
      {
         Iterator pos = mat.erase( 1UL, mat.find( 1UL, 2UL ), mat.find( 1UL, 4UL ) );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 5UL );
         checkCapacity( mat, 7UL );
         checkNonZeros( mat, 5UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 2UL );
         checkNonZeros( mat, 2UL, 2UL );

         if( mat(0,2) != 2 ||
             mat(1,1) != 3 || mat(1,4) != 5 ||
             mat(2,1) != 6 || mat(2,4) != 7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a single-element range failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 0 2 0 0 )\n( 0 3 0 0 5 )\n( 0 6 0 0 7 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos->value() != 5 || pos->index() != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 5\n"
                << "   Expected index: 4\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the elements from (2,4) to the row end
      {
         Iterator pos = mat.erase( 2UL, mat.find( 2UL, 4UL ), mat.end( 2UL ) );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 5UL );
         checkCapacity( mat, 7UL );
         checkNonZeros( mat, 4UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 2UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( mat(0,2) != 2 ||
             mat(1,1) != 3 || mat(1,4) != 5 ||
             mat(2,1) != 6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a single-element range failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 0 2 0 0 )\n( 0 3 0 0 5 )\n( 0 6 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos != mat.end( 2UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to erase an empty range
      {
         Iterator pos = mat.erase( 0UL, mat.find( 0UL, 2UL ), mat.find( 0UL, 2UL ) );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 5UL );
         checkCapacity( mat, 7UL );
         checkNonZeros( mat, 4UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 2UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( mat(0,2) != 2 ||
             mat(1,1) != 3 || mat(1,4) != 5 ||
             mat(2,1) != 6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing an empty range failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 0 2 0 0 )\n( 0 3 0 0 5 )\n( 0 6 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos != mat.find( 0UL, 2UL ) ) {
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
   // Row-major erase function with predicate
   //=====================================================================================

   {
      test_ = "Row-major CompressedMatrix::erase( Predicate )";

      // Initialization check
      blaze::CompressedMatrix<int,blaze::rowMajor> mat{ { 1, 0, 2, 0, 0 },
                                                        { 0, 3, 4, 0, 5 },
                                                        { 0, 6, 0, 0, 7 } };

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 5UL );
      checkCapacity( mat, 7UL );
      checkNonZeros( mat, 7UL );
      checkNonZeros( mat, 0UL, 2UL );
      checkNonZeros( mat, 1UL, 3UL );
      checkNonZeros( mat, 2UL, 2UL );

      if( mat(0,0) != 1 || mat(0,2) != 2 ||
          mat(1,1) != 3 || mat(1,2) != 4 || mat(1,4) != 5 ||
          mat(2,1) != 6 || mat(2,4) != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 0 2 0 0 )\n( 0 3 4 0 5 )\n( 0 6 0 0 7 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing a selection of elements
      mat.erase( []( int value ){ return value == 1 || value == 4 || value == 7; } );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 5UL );
      checkCapacity( mat, 7UL );
      checkNonZeros( mat, 4UL );
      checkNonZeros( mat, 0UL, 1UL );
      checkNonZeros( mat, 1UL, 2UL );
      checkNonZeros( mat, 2UL, 1UL );

      if( mat(0,2) != 2 ||
          mat(1,1) != 3 || mat(1,4) != 5 ||
          mat(2,1) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a selection of elements failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 2 0 0 )\n( 0 3 0 0 5 )\n( 0 6 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to erase all elements with value 1
      mat.erase( []( int value ){ return value == 1; } );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 5UL );
      checkCapacity( mat, 7UL );
      checkNonZeros( mat, 4UL );
      checkNonZeros( mat, 0UL, 1UL );
      checkNonZeros( mat, 1UL, 2UL );
      checkNonZeros( mat, 2UL, 1UL );

      if( mat(0,2) != 2 ||
          mat(1,1) != 3 || mat(1,4) != 5 ||
          mat(2,1) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing all elements with value 1 failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 2 0 0 )\n( 0 3 0 0 5 )\n( 0 6 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major iterator-range-based erase function with predicate
   //=====================================================================================

   {
      test_ = "Row-major CompressedMatrix::erase( size_t, Iterator, Iterator, Predicate )";

      // Initialization check
      blaze::CompressedMatrix<int,blaze::rowMajor> mat{ { 1, 0, 2, 0, 0 },
                                                        { 0, 3, 4, 0, 5 },
                                                        { 0, 6, 0, 0, 7 } };

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 5UL );
      checkCapacity( mat, 7UL );
      checkNonZeros( mat, 7UL );
      checkNonZeros( mat, 0UL, 2UL );
      checkNonZeros( mat, 1UL, 3UL );
      checkNonZeros( mat, 2UL, 2UL );

      if( mat(0,0) != 1 || mat(0,2) != 2 ||
          mat(1,1) != 3 || mat(1,2) != 4 || mat(1,4) != 5 ||
          mat(2,1) != 6 || mat(2,4) != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 0 2 0 0 )\n( 0 3 4 0 5 )\n( 0 6 0 0 7 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing a selection of elements
      mat.erase( 1UL, mat.begin( 1UL ), mat.find( 1UL, 4UL ),
                 []( int value ){ return value == 4; } );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 5UL );
      checkCapacity( mat, 7UL );
      checkNonZeros( mat, 6UL );
      checkNonZeros( mat, 0UL, 2UL );
      checkNonZeros( mat, 1UL, 2UL );
      checkNonZeros( mat, 2UL, 2UL );

      if( mat(0,0) != 1 || mat(0,2) != 2 ||
          mat(1,1) != 3 || mat(1,4) != 5 ||
          mat(2,1) != 6 || mat(2,4) != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a selection of elements failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 0 2 0 0 )\n( 0 3 0 0 5 )\n( 0 6 0 0 7 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to erase from an empty range
      mat.erase( 1UL, mat.begin( 1UL ), mat.begin( 1UL ), []( int ){ return true; } );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 5UL );
      checkCapacity( mat, 7UL );
      checkNonZeros( mat, 6UL );
      checkNonZeros( mat, 0UL, 2UL );
      checkNonZeros( mat, 1UL, 2UL );
      checkNonZeros( mat, 2UL, 2UL );

      if( mat(0,0) != 1 || mat(0,2) != 2 ||
          mat(1,1) != 3 || mat(1,4) != 5 ||
          mat(2,1) != 6 || mat(2,4) != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing from an empty range failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 0 2 0 0 )\n( 0 3 0 0 5 )\n( 0 6 0 0 7 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major index-based erase function
   //=====================================================================================

   {
      test_ = "Column-major CompressedMatrix::erase( size_t, size_t )";

      // Initialization check
      blaze::CompressedMatrix<int,blaze::columnMajor> mat{ { 1, 0, 0 },
                                                           { 0, 3, 6 },
                                                           { 2, 4, 0 },
                                                           { 0, 0, 0 },
                                                           { 0, 5, 7 } };

      checkRows    ( mat, 5UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 7UL );
      checkNonZeros( mat, 7UL );
      checkNonZeros( mat, 0UL, 2UL );
      checkNonZeros( mat, 1UL, 3UL );
      checkNonZeros( mat, 2UL, 2UL );

      if( mat(0,0) != 1 || mat(2,0) != 2 ||
          mat(1,1) != 3 || mat(2,1) != 4 || mat(4,1) != 5 ||
          mat(1,2) != 6 || mat(4,2) != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 3 6 )\n( 2 4 0 )\n( 0 0 0 )\n( 0 5 7 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element at (0,0)
      mat.erase( 0UL, size_t(0) );

      checkRows    ( mat, 5UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 7UL );
      checkNonZeros( mat, 6UL );
      checkNonZeros( mat, 0UL, 1UL );
      checkNonZeros( mat, 1UL, 3UL );
      checkNonZeros( mat, 2UL, 2UL );

      if( mat(2,0) != 2 ||
          mat(1,1) != 3 || mat(2,1) != 4 || mat(4,1) != 5 ||
          mat(1,2) != 6 || mat(4,2) != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 0 )\n( 0 3 6 )\n( 2 4 0 )\n( 0 0 0 )\n( 0 5 7 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element at (2,1)
      mat.erase( 2UL, 1UL );

      checkRows    ( mat, 5UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 7UL );
      checkNonZeros( mat, 5UL );
      checkNonZeros( mat, 0UL, 1UL );
      checkNonZeros( mat, 1UL, 2UL );
      checkNonZeros( mat, 2UL, 2UL );

      if( mat(2,0) != 2 ||
          mat(1,1) != 3 || mat(4,1) != 5 ||
          mat(1,2) != 6 || mat(4,2) != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 0 )\n( 0 3 6 )\n( 2 0 0 )\n( 0 0 0 )\n( 0 5 7 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element at (4,2)
      mat.erase( 4UL, 2UL );

      checkRows    ( mat, 5UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 7UL );
      checkNonZeros( mat, 4UL );
      checkNonZeros( mat, 0UL, 1UL );
      checkNonZeros( mat, 1UL, 2UL );
      checkNonZeros( mat, 2UL, 1UL );

      if( mat(2,0) != 2 ||
          mat(1,1) != 3 || mat(4,1) != 5 ||
          mat(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 0 )\n( 0 3 6 )\n( 2 0 0 )\n( 0 0 0 )\n( 0 5 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to erase a zero element
      mat.erase( 0UL, 1UL );

      checkRows    ( mat, 5UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 7UL );
      checkNonZeros( mat, 4UL );
      checkNonZeros( mat, 0UL, 1UL );
      checkNonZeros( mat, 1UL, 2UL );
      checkNonZeros( mat, 2UL, 1UL );

      if( mat(2,0) != 2 ||
          mat(1,1) != 3 || mat(4,1) != 5 ||
          mat(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 0 )\n( 0 3 6 )\n( 2 0 0 )\n( 0 0 0 )\n( 0 5 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major iterator-based erase function
   //=====================================================================================

   {
      test_ = "Column-major CompressedMatrix::erase( size_t, Iterator )";

      using MatrixType = blaze::CompressedMatrix<int,blaze::columnMajor>;
      using Iterator   = MatrixType::Iterator;

      // Initialization check
      MatrixType mat{ { 1, 0, 0 },
                      { 0, 3, 6 },
                      { 2, 4, 0 },
                      { 0, 0, 0 },
                      { 0, 5, 7 } };

      checkRows    ( mat, 5UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 7UL );
      checkNonZeros( mat, 7UL );
      checkNonZeros( mat, 0UL, 2UL );
      checkNonZeros( mat, 1UL, 3UL );
      checkNonZeros( mat, 2UL, 2UL );

      if( mat(0,0) != 1 || mat(2,0) != 2 ||
          mat(1,1) != 3 || mat(2,1) != 4 || mat(4,1) != 5 ||
          mat(1,2) != 6 || mat(4,2) != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 3 6 )\n( 2 4 0 )\n( 0 0 0 )\n( 0 5 7 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element at (0,0)
      {
         Iterator pos = mat.erase( 0UL, mat.find( 0UL, 0UL ) );

         checkRows    ( mat, 5UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 7UL );
         checkNonZeros( mat, 6UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 3UL );
         checkNonZeros( mat, 2UL, 2UL );

         if( mat(2,0) != 2 ||
             mat(1,1) != 3 || mat(2,1) != 4 || mat(4,1) != 5 ||
             mat(1,2) != 6 || mat(4,2) != 7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 0 0 )\n( 0 3 6 )\n( 2 4 0 )\n( 0 0 0 )\n( 0 5 7 )\n";
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

      // Erasing the element at (2,1)
      {
         Iterator pos = mat.erase( 1UL, mat.find( 2UL, 1UL ) );

         checkRows    ( mat, 5UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 7UL );
         checkNonZeros( mat, 5UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 2UL );
         checkNonZeros( mat, 2UL, 2UL );

         if( mat(2,0) != 2 ||
             mat(1,1) != 3 || mat(4,1) != 5 ||
             mat(1,2) != 6 || mat(4,2) != 7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 0 0 )\n( 0 3 6 )\n( 2 0 0 )\n( 0 0 0 )\n( 0 5 7 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos->value() != 5 || pos->index() != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 5\n"
                << "   Expected index: 4\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the element at (4,2)
      {
         Iterator pos = mat.erase( 2UL, mat.find( 4UL, 2UL ) );

         checkRows    ( mat, 5UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 7UL );
         checkNonZeros( mat, 4UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 2UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( mat(2,0) != 2 ||
             mat(1,1) != 3 || mat(4,1) != 5 ||
             mat(1,2) != 6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 0 0 )\n( 0 3 6 )\n( 2 0 0 )\n( 0 0 0 )\n( 0 5 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos != mat.end( 2UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to erase a zero element
      {
         Iterator pos = mat.erase( 1UL, mat.find( 0UL, 1UL ) );

         checkRows    ( mat, 5UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 7UL );
         checkNonZeros( mat, 4UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 2UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( mat(2,0) != 2 ||
             mat(1,1) != 3 || mat(4,1) != 5 ||
             mat(1,2) != 6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 0 0 )\n( 0 3 6 )\n( 2 0 0 )\n( 0 0 0 )\n( 0 5 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos != mat.end( 1UL ) ) {
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
   // Column-major iterator-range-based erase function
   //=====================================================================================

   {
      test_ = "Column-major CompressedMatrix::erase( size_t, Iterator, Iterator )";

      using MatrixType = blaze::CompressedMatrix<int,blaze::columnMajor>;
      using Iterator   = MatrixType::Iterator;

      // Initialization check
      MatrixType mat{ { 1, 0, 0 },
                      { 0, 3, 6 },
                      { 2, 4, 0 },
                      { 0, 0, 0 },
                      { 0, 5, 7 } };

      checkRows    ( mat, 5UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 7UL );
      checkNonZeros( mat, 7UL );
      checkNonZeros( mat, 0UL, 2UL );
      checkNonZeros( mat, 1UL, 3UL );
      checkNonZeros( mat, 2UL, 2UL );

      if( mat(0,0) != 1 || mat(2,0) != 2 ||
          mat(1,1) != 3 || mat(2,1) != 4 || mat(4,1) != 5 ||
          mat(1,2) != 6 || mat(4,2) != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 3 6 )\n( 2 4 0 )\n( 0 0 0 )\n( 0 5 7 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the elements from (0,0) to (2,0)
      {
         Iterator pos = mat.erase( 0UL, mat.find( 0UL, 0UL ), mat.find( 2UL, 0UL ) );

         checkRows    ( mat, 5UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 7UL );
         checkNonZeros( mat, 6UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 3UL );
         checkNonZeros( mat, 2UL, 2UL );

         if( mat(2,0) != 2 ||
             mat(1,1) != 3 || mat(2,1) != 4 || mat(4,1) != 5 ||
             mat(1,2) != 6 || mat(4,2) != 7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a single-element range failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 0 0 )\n( 0 3 6 )\n( 2 4 0 )\n( 0 0 0 )\n( 0 5 7 )\n";
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

      // Erasing the elements from (2,1) to (4,1)
      {
         Iterator pos = mat.erase( 1UL, mat.find( 2UL, 1UL ), mat.find( 4UL, 1UL ) );

         checkRows    ( mat, 5UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 7UL );
         checkNonZeros( mat, 5UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 2UL );
         checkNonZeros( mat, 2UL, 2UL );

         if( mat(2,0) != 2 ||
             mat(1,1) != 3 || mat(4,1) != 5 ||
             mat(1,2) != 6 || mat(4,2) != 7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a single-element range failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 0 0 )\n( 0 3 6 )\n( 2 0 0 )\n( 0 0 0 )\n( 0 5 7 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos->value() != 5 || pos->index() != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 5\n"
                << "   Expected index: 4\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the elements from (4,2) to the column end
      {
         Iterator pos = mat.erase( 2UL, mat.find( 4UL, 2UL ), mat.end( 2UL ) );

         checkRows    ( mat, 5UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 7UL );
         checkNonZeros( mat, 4UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 2UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( mat(2,0) != 2 ||
             mat(1,1) != 3 || mat(4,1) != 5 ||
             mat(1,2) != 6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a single-element range failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 0 0 )\n( 0 3 6 )\n( 2 0 0 )\n( 0 0 0 )\n( 0 5 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos != mat.end( 2UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to erase an empty range
      {
         Iterator pos = mat.erase( 0UL, mat.find( 2UL, 0UL ), mat.find( 2UL, 0UL ) );

         checkRows    ( mat, 5UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 7UL );
         checkNonZeros( mat, 4UL );
         checkNonZeros( mat, 0UL, 1UL );
         checkNonZeros( mat, 1UL, 2UL );
         checkNonZeros( mat, 2UL, 1UL );

         if( mat(2,0) != 2 ||
             mat(1,1) != 3 || mat(4,1) != 5 ||
             mat(1,2) != 6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing an empty range failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 0 0 )\n( 0 3 6 )\n( 2 0 0 )\n( 0 0 0 )\n( 0 5 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos != mat.find( 2UL, 0UL ) ) {
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
   // Column-major predicate-based erase function
   //=====================================================================================

   {
      test_ = "Column-major CompressedMatrix::erase( Predicate )";

      // Initialization check
      blaze::CompressedMatrix<int,blaze::columnMajor> mat{ { 1, 0, 0 },
                                                           { 0, 3, 6 },
                                                           { 2, 4, 0 },
                                                           { 0, 0, 0 },
                                                           { 0, 5, 7 } };

      checkRows    ( mat, 5UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 7UL );
      checkNonZeros( mat, 7UL );
      checkNonZeros( mat, 0UL, 2UL );
      checkNonZeros( mat, 1UL, 3UL );
      checkNonZeros( mat, 2UL, 2UL );

      if( mat(0,0) != 1 || mat(2,0) != 2 ||
          mat(1,1) != 3 || mat(2,1) != 4 || mat(4,1) != 5 ||
          mat(1,2) != 6 || mat(4,2) != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 3 6 )\n( 2 4 0 )\n( 0 0 0 )\n( 0 5 7 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing a selection of elements
      mat.erase( []( int value ){ return value == 1 || value == 4 || value == 7; } );

      checkRows    ( mat, 5UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 7UL );
      checkNonZeros( mat, 4UL );
      checkNonZeros( mat, 0UL, 1UL );
      checkNonZeros( mat, 1UL, 2UL );
      checkNonZeros( mat, 2UL, 1UL );

      if( mat(2,0) != 2 ||
          mat(1,1) != 3 || mat(4,1) != 5 ||
          mat(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a selection of elements failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 0 )\n( 0 3 6 )\n( 2 0 0 )\n( 0 0 0 )\n( 0 5 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to erase all elements with value 1
      mat.erase( []( int value ){ return value == 1; } );

      checkRows    ( mat, 5UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 7UL );
      checkNonZeros( mat, 4UL );
      checkNonZeros( mat, 0UL, 1UL );
      checkNonZeros( mat, 1UL, 2UL );
      checkNonZeros( mat, 2UL, 1UL );

      if( mat(2,0) != 2 ||
          mat(1,1) != 3 || mat(4,1) != 5 ||
          mat(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing all elements with value 1 failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 0 )\n( 0 3 6 )\n( 2 0 0 )\n( 0 0 0 )\n( 0 5 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major iterator-range-based erase function with predicate
   //=====================================================================================

   {
      test_ = "Column-major CompressedMatrix::erase( size_t, Iterator, Iterator, Predicate )";

      // Initialization check
      blaze::CompressedMatrix<int,blaze::columnMajor> mat{ { 1, 0, 0 },
                                                           { 0, 3, 6 },
                                                           { 2, 4, 0 },
                                                           { 0, 0, 0 },
                                                           { 0, 5, 7 } };

      checkRows    ( mat, 5UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 7UL );
      checkNonZeros( mat, 7UL );
      checkNonZeros( mat, 0UL, 2UL );
      checkNonZeros( mat, 1UL, 3UL );
      checkNonZeros( mat, 2UL, 2UL );

      if( mat(0,0) != 1 || mat(2,0) != 2 ||
          mat(1,1) != 3 || mat(2,1) != 4 || mat(4,1) != 5 ||
          mat(1,2) != 6 || mat(4,2) != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 3 6 )\n( 2 4 0 )\n( 0 0 0 )\n( 0 5 7 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing a selection of elements
      mat.erase( 1UL, mat.begin( 1UL ), mat.find( 4UL, 1UL ),
                 []( int value ){ return value == 4; } );

      checkRows    ( mat, 5UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 7UL );
      checkNonZeros( mat, 6UL );
      checkNonZeros( mat, 0UL, 2UL );
      checkNonZeros( mat, 1UL, 2UL );
      checkNonZeros( mat, 2UL, 2UL );

      if( mat(0,0) != 1 || mat(2,0) != 2 ||
          mat(1,1) != 3 || mat(4,1) != 5 ||
          mat(1,2) != 6 || mat(4,2) != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a selection of elements failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 3 6 )\n( 2 0 0 )\n( 0 0 0 )\n( 0 5 7 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to erase from an empty range
      mat.erase( 1UL, mat.begin( 1UL ), mat.begin( 1UL ), []( int ){ return true; } );

      checkRows    ( mat, 5UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 7UL );
      checkNonZeros( mat, 6UL );
      checkNonZeros( mat, 0UL, 2UL );
      checkNonZeros( mat, 1UL, 2UL );
      checkNonZeros( mat, 2UL, 2UL );

      if( mat(0,0) != 1 || mat(2,0) != 2 ||
          mat(1,1) != 3 || mat(4,1) != 5 ||
          mat(1,2) != 6 || mat(4,2) != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing from an empty range failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 3 6 )\n( 2 0 0 )\n( 0 0 0 )\n( 0 5 7 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c find() member function of the CompressedMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c find() member function of the CompressedMatrix
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testFind()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major CompressedMatrix::find()";

      using ConstIterator = blaze::CompressedMatrix<int,blaze::rowMajor>::ConstIterator;

      // Initialization check
      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 8UL, 6UL, 3UL );
      mat(1,2) = 1;
      mat(2,3) = 2;
      mat(6,5) = 3;

      checkRows    ( mat, 8UL );
      checkColumns ( mat, 6UL );
      checkCapacity( mat, 3UL );
      checkNonZeros( mat, 3UL );
      checkNonZeros( mat, 0UL, 0UL );
      checkNonZeros( mat, 1UL, 1UL );
      checkNonZeros( mat, 2UL, 1UL );
      checkNonZeros( mat, 3UL, 0UL );
      checkNonZeros( mat, 4UL, 0UL );
      checkNonZeros( mat, 5UL, 0UL );
      checkNonZeros( mat, 6UL, 1UL );
      checkNonZeros( mat, 7UL, 0UL );

      // Searching for the first element
      {
         ConstIterator pos( mat.find( 1UL, 2UL ) );

         if( pos == mat.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required position = (1,2)\n"
                << "   Current matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 2 || pos->value() != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 2\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 1\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for the second element
      {
         ConstIterator pos( mat.find( 2UL, 3UL ) );

         if( pos == mat.end( 2UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required position = (2,3)\n"
                << "   Current matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 3 || pos->value() != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 3\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 2\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for the third element
      {
         ConstIterator pos( mat.find( 6UL, 5UL ) );

         if( pos == mat.end( 6UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required position = (6,5)\n"
                << "   Current matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 5 || pos->value() != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 5\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 3\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for a non-existing non-zero element
      {
         ConstIterator pos( mat.find( 4UL, 0UL ) );

         if( pos != mat.end( 4UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Non-existing element could be found\n"
                << " Details:\n"
                << "   Required index = 4\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 0\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major CompressedMatrix::find()";

      using ConstIterator = blaze::CompressedMatrix<int,blaze::columnMajor>::ConstIterator;

      // Initialization check
      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 8UL, 6UL, 3UL );
      mat(1,2) = 1;
      mat(2,3) = 2;
      mat(6,5) = 3;

      checkRows    ( mat, 8UL );
      checkColumns ( mat, 6UL );
      checkCapacity( mat, 3UL );
      checkNonZeros( mat, 3UL );
      checkNonZeros( mat, 0UL, 0UL );
      checkNonZeros( mat, 1UL, 0UL );
      checkNonZeros( mat, 2UL, 1UL );
      checkNonZeros( mat, 3UL, 1UL );
      checkNonZeros( mat, 4UL, 0UL );
      checkNonZeros( mat, 5UL, 1UL );

      // Searching for the first element
      {
         ConstIterator pos( mat.find( 1UL, 2UL ) );

         if( pos == mat.end( 2UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required position = (1,2)\n"
                << "   Current matrix:\n" << mat << "\n";
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
                << "   Current matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for the second element
      {
         ConstIterator pos( mat.find( 2UL, 3UL ) );

         if( pos == mat.end( 3UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required position = (2,3)\n"
                << "   Current matrix:\n" << mat << "\n";
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
                << "   Current matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for the third element
      {
         ConstIterator pos( mat.find( 6UL, 5UL ) );

         if( pos == mat.end( 5UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required position = (6,5)\n"
                << "   Current matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 6 || pos->value() != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 6\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 3\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for a non-existing non-zero element
      {
         ConstIterator pos( mat.find( 4UL, 0UL ) );

         if( pos != mat.end( 0UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Non-existing element could be found\n"
                << " Details:\n"
                << "   Required index = 0\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 0\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c lowerBound() member function of the CompressedMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c lowerBound() member function of the CompressedMatrix
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testLowerBound()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major CompressedMatrix::lowerBound()";

      using ConstIterator = blaze::CompressedMatrix<int,blaze::rowMajor>::ConstIterator;

      // Initialization check
      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 6UL, 3UL );
      mat(1,2) = 1;
      mat(1,4) = 2;

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 6UL );
      checkCapacity( mat, 2UL );
      checkNonZeros( mat, 2UL );
      checkNonZeros( mat, 0UL, 0UL );
      checkNonZeros( mat, 1UL, 2UL );
      checkNonZeros( mat, 2UL, 0UL );

      // Determining the lower bound for position (1,1)
      {
         ConstIterator pos( mat.lowerBound( 1UL, 1UL ) );

         if( pos == mat.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,1)\n"
                << "   Current matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 2 || pos->value() != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 2\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 1\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for position (1,2)
      {
         ConstIterator pos( mat.lowerBound( 1UL, 2UL ) );

         if( pos == mat.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,2)\n"
                << "   Current matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 2 || pos->value() != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 2\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 1\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for position (1,3)
      {
         ConstIterator pos( mat.lowerBound( 1UL, 3UL ) );

         if( pos == mat.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,3)\n"
                << "   Current matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 4 || pos->value() != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 4\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 2\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for position (1,4)
      {
         ConstIterator pos( mat.lowerBound( 1UL, 4UL ) );

         if( pos == mat.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,4)\n"
                << "   Current matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 4 || pos->value() != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 4\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 2\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for position (1,5)
      {
         ConstIterator pos( mat.lowerBound( 1UL, 5UL ) );

         if( pos != mat.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,5)\n"
                << "   Current matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major CompressedMatrix::lowerBound()";

      using ConstIterator = blaze::CompressedMatrix<int,blaze::columnMajor>::ConstIterator;

      // Initialization check
      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 6UL, 3UL, 3UL );
      mat(2,1) = 1;
      mat(4,1) = 2;

      checkRows    ( mat, 6UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 2UL );
      checkNonZeros( mat, 2UL );
      checkNonZeros( mat, 0UL, 0UL );
      checkNonZeros( mat, 1UL, 2UL );
      checkNonZeros( mat, 2UL, 0UL );

      // Determining the lower bound for position (1,1)
      {
         ConstIterator pos( mat.lowerBound( 1UL, 1UL ) );

         if( pos == mat.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,1)\n"
                << "   Current matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 2 || pos->value() != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 2\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 1\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for position (2,1)
      {
         ConstIterator pos( mat.lowerBound( 2UL, 1UL ) );

         if( pos == mat.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (2,1)\n"
                << "   Current matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 2 || pos->value() != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 2\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 1\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for position (3,1)
      {
         ConstIterator pos( mat.lowerBound( 3UL, 1UL ) );

         if( pos == mat.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (3,1)\n"
                << "   Current matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 4 || pos->value() != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 4\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 2\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for position (4,1)
      {
         ConstIterator pos( mat.lowerBound( 4UL, 1UL ) );

         if( pos == mat.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (4,1)\n"
                << "   Current matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 4 || pos->value() != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 4\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 2\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for position (5,1)
      {
         ConstIterator pos( mat.lowerBound( 5UL, 1UL ) );

         if( pos != mat.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (5,1)\n"
                << "   Current matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c upperBound() member function of the CompressedMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c upperBound() member function of the CompressedMatrix
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testUpperBound()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major CompressedMatrix::upperBound()";

      using ConstIterator = blaze::CompressedMatrix<int,blaze::rowMajor>::ConstIterator;

      // Initialization check
      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 6UL, 3UL );
      mat(1,2) = 1;
      mat(1,4) = 2;

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 6UL );
      checkCapacity( mat, 2UL );
      checkNonZeros( mat, 2UL );
      checkNonZeros( mat, 0UL, 0UL );
      checkNonZeros( mat, 1UL, 2UL );
      checkNonZeros( mat, 2UL, 0UL );

      // Determining the upper bound for position (1,1)
      {
         ConstIterator pos( mat.upperBound( 1UL, 1UL ) );

         if( pos == mat.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,1)\n"
                << "   Current matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 2 || pos->value() != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 2\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 1\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for position (1,2)
      {
         ConstIterator pos( mat.upperBound( 1UL, 2UL ) );

         if( pos == mat.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,2)\n"
                << "   Current matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 4 || pos->value() != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 4\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 2\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for position (1,3)
      {
         ConstIterator pos( mat.upperBound( 1UL, 3UL ) );

         if( pos == mat.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,3)\n"
                << "   Current matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 4 || pos->value() != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 4\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 2\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for position (1,4)
      {
         ConstIterator pos( mat.upperBound( 1UL, 4UL ) );

         if( pos != mat.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,4)\n"
                << "   Current matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for position (1,5)
      {
         ConstIterator pos( mat.upperBound( 1UL, 5UL ) );

         if( pos != mat.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,5)\n"
                << "   Current matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major CompressedMatrix::upperBound()";

      using ConstIterator = blaze::CompressedMatrix<int,blaze::columnMajor>::ConstIterator;

      // Initialization check
      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 6UL, 3UL, 3UL );
      mat(2,1) = 1;
      mat(4,1) = 2;

      checkRows    ( mat, 6UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 2UL );
      checkNonZeros( mat, 2UL );
      checkNonZeros( mat, 0UL, 0UL );
      checkNonZeros( mat, 1UL, 2UL );
      checkNonZeros( mat, 2UL, 0UL );

      // Determining the upper bound for position (1,1)
      {
         ConstIterator pos( mat.upperBound( 1UL, 1UL ) );

         if( pos == mat.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,1)\n"
                << "   Current matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 2 || pos->value() != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 2\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 1\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for position (2,1)
      {
         ConstIterator pos( mat.upperBound( 2UL, 1UL ) );

         if( pos == mat.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (2,1)\n"
                << "   Current matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 4 || pos->value() != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 4\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 2\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for position (3,1)
      {
         ConstIterator pos( mat.upperBound( 3UL, 1UL ) );

         if( pos == mat.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (3,1)\n"
                << "   Current matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 4 || pos->value() != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 4\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 2\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for position (4,1)
      {
         ConstIterator pos( mat.upperBound( 4UL, 1UL ) );

         if( pos != mat.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (4,1)\n"
                << "   Current matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for position (5,1)
      {
         ConstIterator pos( mat.upperBound( 5UL, 1UL ) );

         if( pos != mat.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (5,1)\n"
                << "   Current matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c transpose() member function of the CompressedMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c transpose() member function of the CompressedMatrix
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

      // Self-transpose of a 3x5 matrix
      {
         blaze::CompressedMatrix<int,blaze::rowMajor> mat{ { 1, 0, 2, 0, 3 },
                                                           { 0, 4, 0, 5, 0 },
                                                           { 6, 0, 7, 0, 8 } };

         transpose( mat );

         checkRows    ( mat, 5UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 8UL );
         checkNonZeros( mat, 8UL );
         checkNonZeros( mat, 0UL, 2UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 2UL );
         checkNonZeros( mat, 3UL, 1UL );
         checkNonZeros( mat, 4UL, 2UL );

         if( mat(0,0) != 1 || mat(2,0) != 2 || mat(4,0) != 3 || mat(1,1) != 4 ||
             mat(3,1) != 5 || mat(0,2) != 6 || mat(2,2) != 7 || mat(4,2) != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Transposition failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 1 0 6 )\n( 0 4 0 )\n( 2 0 7 )\n( 0 5 0 )\n( 3 0 8 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Self-transpose of a 5x3 matrix
      {
         blaze::CompressedMatrix<int,blaze::rowMajor> mat{ { 1, 0, 6 },
                                                           { 0, 4, 0 },
                                                           { 2, 0, 7 },
                                                           { 0, 5, 0 },
                                                           { 3, 0, 8 } };

         transpose( mat );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 5UL );
         checkCapacity( mat, 8UL );
         checkNonZeros( mat, 8UL );
         checkNonZeros( mat, 0UL, 3UL );
         checkNonZeros( mat, 1UL, 2UL );
         checkNonZeros( mat, 2UL, 3UL );

         if( mat(0,0) != 1 || mat(0,1) != 0 || mat(0,2) != 2 || mat(0,3) != 0 || mat(0,4) != 3 ||
             mat(1,0) != 0 || mat(1,1) != 4 || mat(1,2) != 0 || mat(1,3) != 5 || mat(1,4) != 0 ||
             mat(2,0) != 6 || mat(2,1) != 0 || mat(2,2) != 7 || mat(2,3) != 0 || mat(2,4) != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Transposition failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 1 0 2 0 3 )\n( 0 4 0 5 0 )\n( 6 0 7 0 8 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   {
      test_ = "Row-major self-transpose via trans()";

      // Self-transpose of a 3x5 matrix
      {
         blaze::CompressedMatrix<int,blaze::rowMajor> mat{ { 1, 0, 2, 0, 3 },
                                                           { 0, 4, 0, 5, 0 },
                                                           { 6, 0, 7, 0, 8 } };

         mat = trans( mat );

         checkRows    ( mat, 5UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 8UL );
         checkNonZeros( mat, 8UL );
         checkNonZeros( mat, 0UL, 2UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 2UL );
         checkNonZeros( mat, 3UL, 1UL );
         checkNonZeros( mat, 4UL, 2UL );

         if( mat(0,0) != 1 || mat(2,0) != 2 || mat(4,0) != 3 || mat(1,1) != 4 ||
             mat(3,1) != 5 || mat(0,2) != 6 || mat(2,2) != 7 || mat(4,2) != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Transposition failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 1 0 6 )\n( 0 4 0 )\n( 2 0 7 )\n( 0 5 0 )\n( 3 0 8 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Self-transpose of a 5x3 matrix
      {
         blaze::CompressedMatrix<int,blaze::rowMajor> mat{ { 1, 0, 6 },
                                                           { 0, 4, 0 },
                                                           { 2, 0, 7 },
                                                           { 0, 5, 0 },
                                                           { 3, 0, 8 } };

         mat = trans( mat );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 5UL );
         checkCapacity( mat, 8UL );
         checkNonZeros( mat, 8UL );
         checkNonZeros( mat, 0UL, 3UL );
         checkNonZeros( mat, 1UL, 2UL );
         checkNonZeros( mat, 2UL, 3UL );

         if( mat(0,0) != 1 || mat(0,1) != 0 || mat(0,2) != 2 || mat(0,3) != 0 || mat(0,4) != 3 ||
             mat(1,0) != 0 || mat(1,1) != 4 || mat(1,2) != 0 || mat(1,3) != 5 || mat(1,4) != 0 ||
             mat(2,0) != 6 || mat(2,1) != 0 || mat(2,2) != 7 || mat(2,3) != 0 || mat(2,4) != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Transposition failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 1 0 2 0 3 )\n( 0 4 0 5 0 )\n( 6 0 7 0 8 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major self-transpose via transpose()";

      // Self-transpose of a 3x5 matrix
      {
         blaze::CompressedMatrix<int,blaze::columnMajor> mat{ { 1, 0, 2, 0, 3 },
                                                              { 0, 4, 0, 5, 0 },
                                                              { 6, 0, 7, 0, 8 } };

         transpose( mat );

         checkRows    ( mat, 5UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 8UL );
         checkNonZeros( mat, 8UL );
         checkNonZeros( mat, 0UL, 3UL );
         checkNonZeros( mat, 1UL, 2UL );
         checkNonZeros( mat, 2UL, 3UL );

         if( mat(0,0) != 1 || mat(2,0) != 2 || mat(4,0) != 3 || mat(1,1) != 4 ||
             mat(3,1) != 5 || mat(0,2) != 6 || mat(2,2) != 7 || mat(4,2) != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Transposition failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 1 0 6 )\n( 0 4 0 )\n( 2 0 7 )\n( 0 5 0 )\n( 3 0 8 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Self-transpose of a 5x3 matrix
      {
         blaze::CompressedMatrix<int,blaze::columnMajor> mat{ { 1, 0, 6 },
                                                              { 0, 4, 0 },
                                                              { 2, 0, 7 },
                                                              { 0, 5, 0 },
                                                              { 3, 0, 8 } };

         transpose( mat );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 5UL );
         checkCapacity( mat, 8UL );
         checkNonZeros( mat, 8UL );
         checkNonZeros( mat, 0UL, 2UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 2UL );
         checkNonZeros( mat, 3UL, 1UL );
         checkNonZeros( mat, 4UL, 2UL );

         if( mat(0,0) != 1 || mat(0,1) != 0 || mat(0,2) != 2 || mat(0,3) != 0 || mat(0,4) != 3 ||
             mat(1,0) != 0 || mat(1,1) != 4 || mat(1,2) != 0 || mat(1,3) != 5 || mat(1,4) != 0 ||
             mat(2,0) != 6 || mat(2,1) != 0 || mat(2,2) != 7 || mat(2,3) != 0 || mat(2,4) != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Transposition failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 1 0 2 0 3 )\n( 0 4 0 5 0 )\n( 6 0 7 0 8 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   {
      test_ = "Column-major self-transpose via trans()";

      // Self-transpose of a 3x5 matrix
      {
         blaze::CompressedMatrix<int,blaze::columnMajor> mat{ { 1, 0, 2, 0, 3 },
                                                              { 0, 4, 0, 5, 0 },
                                                              { 6, 0, 7, 0, 8 } };

         mat = trans( mat );

         checkRows    ( mat, 5UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 8UL );
         checkNonZeros( mat, 8UL );
         checkNonZeros( mat, 0UL, 3UL );
         checkNonZeros( mat, 1UL, 2UL );
         checkNonZeros( mat, 2UL, 3UL );

         if( mat(0,0) != 1 || mat(2,0) != 2 || mat(4,0) != 3 || mat(1,1) != 4 ||
             mat(3,1) != 5 || mat(0,2) != 6 || mat(2,2) != 7 || mat(4,2) != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Transposition failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 1 0 6 )\n( 0 4 0 )\n( 2 0 7 )\n( 0 5 0 )\n( 3 0 8 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Self-transpose of a 5x3 matrix
      {
         blaze::CompressedMatrix<int,blaze::columnMajor> mat{ { 1, 0, 6 },
                                                              { 0, 4, 0 },
                                                              { 2, 0, 7 },
                                                              { 0, 5, 0 },
                                                              { 3, 0, 8 } };

         mat = trans( mat );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 5UL );
         checkCapacity( mat, 8UL );
         checkNonZeros( mat, 8UL );
         checkNonZeros( mat, 0UL, 2UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 2UL );
         checkNonZeros( mat, 3UL, 1UL );
         checkNonZeros( mat, 4UL, 2UL );

         if( mat(0,0) != 1 || mat(0,1) != 0 || mat(0,2) != 2 || mat(0,3) != 0 || mat(0,4) != 3 ||
             mat(1,0) != 0 || mat(1,1) != 4 || mat(1,2) != 0 || mat(1,3) != 5 || mat(1,4) != 0 ||
             mat(2,0) != 6 || mat(2,1) != 0 || mat(2,2) != 7 || mat(2,3) != 0 || mat(2,4) != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Transposition failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 1 0 2 0 3 )\n( 0 4 0 5 0 )\n( 6 0 7 0 8 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c ctranspose() member function of the CompressedMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c ctranspose() member function of the CompressedMatrix
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

      // Self-transpose of a 3x5 matrix
      {
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 3UL, 5UL, 8UL );
         mat(0,0) = cplx(1,-1);
         mat(0,2) = cplx(2,-2);
         mat(0,4) = cplx(3,-3);
         mat(1,1) = cplx(4,-4);
         mat(1,3) = cplx(5,-5);
         mat(2,0) = cplx(6,-6);
         mat(2,2) = cplx(7,-7);
         mat(2,4) = cplx(8,-8);

         ctranspose( mat );

         checkRows    ( mat, 5UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 8UL );
         checkNonZeros( mat, 8UL );
         checkNonZeros( mat, 0UL, 2UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 2UL );
         checkNonZeros( mat, 3UL, 1UL );
         checkNonZeros( mat, 4UL, 2UL );

         if( mat(0,0) != cplx(1,1) || mat(0,1) != cplx(0,0) || mat(0,2) != cplx(6,6) ||
             mat(1,0) != cplx(0,0) || mat(1,1) != cplx(4,4) || mat(1,2) != cplx(0,0) ||
             mat(2,0) != cplx(2,2) || mat(2,1) != cplx(0,0) || mat(2,2) != cplx(7,7) ||
             mat(3,0) != cplx(0,0) || mat(3,1) != cplx(5,5) || mat(3,2) != cplx(0,0) ||
             mat(4,0) != cplx(3,3) || mat(4,1) != cplx(0,0) || mat(4,2) != cplx(8,8) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Transposition failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( (1,1) (0,0) (6,6) )\n"
                                        "( (0,0) (4,4) (0,0) )\n"
                                        "( (2,2) (0,0) (7,7) )\n"
                                        "( (0,0) (5,5) (0,0) )\n"
                                        "( (3,3) (0,0) (8,8) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Self-transpose of a 5x3 matrix
      {
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 5UL, 3UL, 8UL );
         mat(0,0) = cplx(1,-1);
         mat(0,2) = cplx(6,-6);
         mat(1,1) = cplx(4,-4);
         mat(2,0) = cplx(2,-2);
         mat(2,2) = cplx(7,-7);
         mat(3,1) = cplx(5,-5);
         mat(4,0) = cplx(3,-3);
         mat(4,2) = cplx(8,-8);

         ctranspose( mat );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 5UL );
         checkCapacity( mat, 8UL );
         checkNonZeros( mat, 8UL );
         checkNonZeros( mat, 0UL, 3UL );
         checkNonZeros( mat, 1UL, 2UL );
         checkNonZeros( mat, 2UL, 3UL );

         if( mat(0,0) != cplx(1,1) || mat(0,1) != cplx(0,0) || mat(0,2) != cplx(2,2) || mat(0,3) != cplx(0,0) || mat(0,4) != cplx(3,3) ||
             mat(1,0) != cplx(0,0) || mat(1,1) != cplx(4,4) || mat(1,2) != cplx(0,0) || mat(1,3) != cplx(5,5) || mat(1,4) != cplx(0,0) ||
             mat(2,0) != cplx(6,6) || mat(2,1) != cplx(0,0) || mat(2,2) != cplx(7,7) || mat(2,3) != cplx(0,0) || mat(2,4) != cplx(8,8) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Transposition failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( (1,1) (0,0) (2,2) (0,0) (3,3) )\n"
                                        "( (0,0) (4,4) (0,0) (5,5) (0,0) )\n"
                                        "( (6,6) (0,0) (7,7) (0,0) (8,8) )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   {
      test_ = "Row-major self-transpose via ctrans()";

      using cplx = blaze::complex<int>;

      // Self-transpose of a 3x5 matrix
      {
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 3UL, 5UL, 8UL );
         mat(0,0) = cplx(1,-1);
         mat(0,2) = cplx(2,-2);
         mat(0,4) = cplx(3,-3);
         mat(1,1) = cplx(4,-4);
         mat(1,3) = cplx(5,-5);
         mat(2,0) = cplx(6,-6);
         mat(2,2) = cplx(7,-7);
         mat(2,4) = cplx(8,-8);

         mat = ctrans( mat );

         checkRows    ( mat, 5UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 8UL );
         checkNonZeros( mat, 8UL );
         checkNonZeros( mat, 0UL, 2UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 2UL );
         checkNonZeros( mat, 3UL, 1UL );
         checkNonZeros( mat, 4UL, 2UL );

         if( mat(0,0) != cplx(1,1) || mat(0,1) != cplx(0,0) || mat(0,2) != cplx(6,6) ||
             mat(1,0) != cplx(0,0) || mat(1,1) != cplx(4,4) || mat(1,2) != cplx(0,0) ||
             mat(2,0) != cplx(2,2) || mat(2,1) != cplx(0,0) || mat(2,2) != cplx(7,7) ||
             mat(3,0) != cplx(0,0) || mat(3,1) != cplx(5,5) || mat(3,2) != cplx(0,0) ||
             mat(4,0) != cplx(3,3) || mat(4,1) != cplx(0,0) || mat(4,2) != cplx(8,8) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Transposition failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( (1,1) (0,0) (6,6) )\n"
                                        "( (0,0) (4,4) (0,0) )\n"
                                        "( (2,2) (0,0) (7,7) )\n"
                                        "( (0,0) (5,5) (0,0) )\n"
                                        "( (3,3) (0,0) (8,8) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Self-transpose of a 5x3 matrix
      {
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 5UL, 3UL, 8UL );
         mat(0,0) = cplx(1,-1);
         mat(0,2) = cplx(6,-6);
         mat(1,1) = cplx(4,-4);
         mat(2,0) = cplx(2,-2);
         mat(2,2) = cplx(7,-7);
         mat(3,1) = cplx(5,-5);
         mat(4,0) = cplx(3,-3);
         mat(4,2) = cplx(8,-8);

         mat = ctrans( mat );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 5UL );
         checkCapacity( mat, 8UL );
         checkNonZeros( mat, 8UL );
         checkNonZeros( mat, 0UL, 3UL );
         checkNonZeros( mat, 1UL, 2UL );
         checkNonZeros( mat, 2UL, 3UL );

         if( mat(0,0) != cplx(1,1) || mat(0,1) != cplx(0,0) || mat(0,2) != cplx(2,2) || mat(0,3) != cplx(0,0) || mat(0,4) != cplx(3,3) ||
             mat(1,0) != cplx(0,0) || mat(1,1) != cplx(4,4) || mat(1,2) != cplx(0,0) || mat(1,3) != cplx(5,5) || mat(1,4) != cplx(0,0) ||
             mat(2,0) != cplx(6,6) || mat(2,1) != cplx(0,0) || mat(2,2) != cplx(7,7) || mat(2,3) != cplx(0,0) || mat(2,4) != cplx(8,8) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Transposition failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( (1,1) (0,0) (2,2) (0,0) (3,3) )\n"
                                        "( (0,0) (4,4) (0,0) (5,5) (0,0) )\n"
                                        "( (6,6) (0,0) (7,7) (0,0) (8,8) )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major self-transpose via ctranspose()";

      using cplx = blaze::complex<int>;

      // Self-transpose of a 3x5 matrix
      {
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 3UL, 5UL, 8UL );
         mat(0,0) = cplx(1,-1);
         mat(0,2) = cplx(2,-2);
         mat(0,4) = cplx(3,-3);
         mat(1,1) = cplx(4,-4);
         mat(1,3) = cplx(5,-5);
         mat(2,0) = cplx(6,-6);
         mat(2,2) = cplx(7,-7);
         mat(2,4) = cplx(8,-8);

         ctranspose( mat );

         checkRows    ( mat, 5UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 8UL );
         checkNonZeros( mat, 8UL );
         checkNonZeros( mat, 0UL, 3UL );
         checkNonZeros( mat, 1UL, 2UL );
         checkNonZeros( mat, 2UL, 3UL );

         if( mat(0,0) != cplx(1,1) || mat(0,1) != cplx(0,0) || mat(0,2) != cplx(6,6) ||
             mat(1,0) != cplx(0,0) || mat(1,1) != cplx(4,4) || mat(1,2) != cplx(0,0) ||
             mat(2,0) != cplx(2,2) || mat(2,1) != cplx(0,0) || mat(2,2) != cplx(7,7) ||
             mat(3,0) != cplx(0,0) || mat(3,1) != cplx(5,5) || mat(3,2) != cplx(0,0) ||
             mat(4,0) != cplx(3,3) || mat(4,1) != cplx(0,0) || mat(4,2) != cplx(8,8) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Transposition failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( (1,1) (0,0) (6,6) )\n"
                                        "( (0,0) (4,4) (0,0) )\n"
                                        "( (2,2) (0,0) (7,7) )\n"
                                        "( (0,0) (5,5) (0,0) )\n"
                                        "( (3,3) (0,0) (8,8) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Self-transpose of a 5x3 matrix
      {
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 5UL, 3UL, 8UL );
         mat(0,0) = cplx(1,-1);
         mat(0,2) = cplx(6,-6);
         mat(1,1) = cplx(4,-4);
         mat(2,0) = cplx(2,-2);
         mat(2,2) = cplx(7,-7);
         mat(3,1) = cplx(5,-5);
         mat(4,0) = cplx(3,-3);
         mat(4,2) = cplx(8,-8);

         ctranspose( mat );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 5UL );
         checkCapacity( mat, 8UL );
         checkNonZeros( mat, 8UL );
         checkNonZeros( mat, 0UL, 2UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 2UL );
         checkNonZeros( mat, 3UL, 1UL );
         checkNonZeros( mat, 4UL, 2UL );

         if( mat(0,0) != cplx(1,1) || mat(0,1) != cplx(0,0) || mat(0,2) != cplx(2,2) || mat(0,3) != cplx(0,0) || mat(0,4) != cplx(3,3) ||
             mat(1,0) != cplx(0,0) || mat(1,1) != cplx(4,4) || mat(1,2) != cplx(0,0) || mat(1,3) != cplx(5,5) || mat(1,4) != cplx(0,0) ||
             mat(2,0) != cplx(6,6) || mat(2,1) != cplx(0,0) || mat(2,2) != cplx(7,7) || mat(2,3) != cplx(0,0) || mat(2,4) != cplx(8,8) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Transposition failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( (1,1) (0,0) (2,2) (0,0) (3,3) )\n"
                                        "( (0,0) (4,4) (0,0) (5,5) (0,0) )\n"
                                        "( (6,6) (0,0) (7,7) (0,0) (8,8) )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   {
      test_ = "Column-major self-transpose via ctrans()";

      using cplx = blaze::complex<int>;

      // Self-transpose of a 3x5 matrix
      {
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 3UL, 5UL, 8UL );
         mat(0,0) = cplx(1,-1);
         mat(0,2) = cplx(2,-2);
         mat(0,4) = cplx(3,-3);
         mat(1,1) = cplx(4,-4);
         mat(1,3) = cplx(5,-5);
         mat(2,0) = cplx(6,-6);
         mat(2,2) = cplx(7,-7);
         mat(2,4) = cplx(8,-8);

         mat = ctrans( mat );

         checkRows    ( mat, 5UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 8UL );
         checkNonZeros( mat, 8UL );
         checkNonZeros( mat, 0UL, 3UL );
         checkNonZeros( mat, 1UL, 2UL );
         checkNonZeros( mat, 2UL, 3UL );

         if( mat(0,0) != cplx(1,1) || mat(0,1) != cplx(0,0) || mat(0,2) != cplx(6,6) ||
             mat(1,0) != cplx(0,0) || mat(1,1) != cplx(4,4) || mat(1,2) != cplx(0,0) ||
             mat(2,0) != cplx(2,2) || mat(2,1) != cplx(0,0) || mat(2,2) != cplx(7,7) ||
             mat(3,0) != cplx(0,0) || mat(3,1) != cplx(5,5) || mat(3,2) != cplx(0,0) ||
             mat(4,0) != cplx(3,3) || mat(4,1) != cplx(0,0) || mat(4,2) != cplx(8,8) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Transposition failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( (1,1) (0,0) (6,6) )\n"
                                        "( (0,0) (4,4) (0,0) )\n"
                                        "( (2,2) (0,0) (7,7) )\n"
                                        "( (0,0) (5,5) (0,0) )\n"
                                        "( (3,3) (0,0) (8,8) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Self-transpose of a 5x3 matrix
      {
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 5UL, 3UL, 8UL );
         mat(0,0) = cplx(1,-1);
         mat(0,2) = cplx(6,-6);
         mat(1,1) = cplx(4,-4);
         mat(2,0) = cplx(2,-2);
         mat(2,2) = cplx(7,-7);
         mat(3,1) = cplx(5,-5);
         mat(4,0) = cplx(3,-3);
         mat(4,2) = cplx(8,-8);

         mat = ctrans( mat );

         checkRows    ( mat, 3UL );
         checkColumns ( mat, 5UL );
         checkCapacity( mat, 8UL );
         checkNonZeros( mat, 8UL );
         checkNonZeros( mat, 0UL, 2UL );
         checkNonZeros( mat, 1UL, 1UL );
         checkNonZeros( mat, 2UL, 2UL );
         checkNonZeros( mat, 3UL, 1UL );
         checkNonZeros( mat, 4UL, 2UL );

         if( mat(0,0) != cplx(1,1) || mat(0,1) != cplx(0,0) || mat(0,2) != cplx(2,2) || mat(0,3) != cplx(0,0) || mat(0,4) != cplx(3,3) ||
             mat(1,0) != cplx(0,0) || mat(1,1) != cplx(4,4) || mat(1,2) != cplx(0,0) || mat(1,3) != cplx(5,5) || mat(1,4) != cplx(0,0) ||
             mat(2,0) != cplx(6,6) || mat(2,1) != cplx(0,0) || mat(2,2) != cplx(7,7) || mat(2,3) != cplx(0,0) || mat(2,4) != cplx(8,8) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Transposition failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( (1,1) (0,0) (2,2) (0,0) (3,3) )\n"
                                        "( (0,0) (4,4) (0,0) (5,5) (0,0) )\n"
                                        "( (6,6) (0,0) (7,7) (0,0) (8,8) )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c isDefault() function with the CompressedMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isDefault() function with the CompressedMatrix class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testIsDefault()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major isDefault() function";

      // isDefault with 0x0 matrix
      {
         blaze::CompressedMatrix<int,blaze::rowMajor> mat;

         if( isDefault( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with default matrix
      {
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 3UL );

         if( isDefault( mat(0,1) ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix element: " << mat(0,1) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with non-default matrix
      {
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 2UL );
         mat(0,1) = 1;

         if( isDefault( mat(0,1) ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix element: " << mat(0,1) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major isDefault() function";

      // isDefault with 0x0 matrix
      {
         blaze::CompressedMatrix<int,blaze::columnMajor> mat;

         if( isDefault( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with default matrix
      {
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 2UL, 3UL );

         if( isDefault( mat(1,0) ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix element: " << mat(1,0) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with non-default matrix
      {
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 2UL, 1UL );
         mat(1,0) = 1;

         if( isDefault( mat(1,0) ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix element: " << mat(1,0) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
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
   std::cout << "   Running CompressedMatrix class test (part 2)..." << std::endl;

   try
   {
      RUN_COMPRESSEDMATRIX_CLASS_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during CompressedMatrix class test (part 2):\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
