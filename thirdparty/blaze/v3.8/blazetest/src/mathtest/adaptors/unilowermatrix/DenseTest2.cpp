//=================================================================================================
/*!
//  \file src/mathtest/adaptors/unilowermatrix/DenseTest2.cpp
//  \brief Source file for the UniLowerMatrix dense test (part 2)
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
#include <blaze/math/Column.h>
#include <blaze/math/Row.h>
#include <blaze/math/StaticMatrix.h>
#include <blaze/math/Submatrix.h>
#include <blazetest/mathtest/adaptors/unilowermatrix/DenseTest.h>

#ifdef BLAZE_USE_HPX_THREADS
#  include <hpx/hpx_main.hpp>
#endif


namespace blazetest {

namespace mathtest {

namespace adaptors {

namespace unilowermatrix {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the UniLowerMatrix dense test.
//
// \exception std::runtime_error Operation error detected.
*/
DenseTest::DenseTest()
{
   testFunctionCall();
   testIterator();
   testNonZeros();
   testReset();
   testClear();
   testResize();
   testExtend();
   testReserve();
   testShrinkToFit();
   testSwap();
   testIsDefault();
   testSubmatrix();
   testRow();
   testColumn();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the UniLowerMatrix function call operator.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of adding and accessing elements via the function call operator
// of the UniLowerMatrix specialization. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
void DenseTest::testFunctionCall()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major UniLowerMatrix::operator()";

      // Good cases
      {
         LT lower( 3UL );

         // Writing the lower element (2,1)
         lower(2,1) = 2;

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 9UL );
         checkNonZeros( lower, 4UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 1UL );
         checkNonZeros( lower, 2UL, 2UL );

         if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
             lower(1,0) != 0 || lower(1,1) != 1 || lower(1,2) != 0 ||
             lower(2,0) != 0 || lower(2,1) != 2 || lower(2,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 )\n( 0 1 0 )\n( 0 2 1 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Writing the lower element (1,0)
         lower(1,0) = lower(2,1);

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 9UL );
         checkNonZeros( lower, 5UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 2UL );

         if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
             lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 ||
             lower(2,0) != 0 || lower(2,1) != 2 || lower(2,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 )\n( 2 1 0 )\n( 0 2 1 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Adding to the lower element (2,0)
         lower(2,0) += 3;

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 9UL );
         checkNonZeros( lower, 6UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 3UL );

         if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
             lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 ||
             lower(2,0) != 3 || lower(2,1) != 2 || lower(2,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 )\n( 2 1 0 )\n( 3 2 1 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Subtracting from the lower element (1,0)
         lower(1,0) -= 4;

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 9UL );
         checkNonZeros( lower, 6UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 3UL );

         if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
             lower(1,0) != -2 || lower(1,1) != 1 || lower(1,2) != 0 ||
             lower(2,0) !=  3 || lower(2,1) != 2 || lower(2,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n(  1 0 0 )\n( -2 1 0 )\n(  3 2 1 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Multiplying the lower element (2,1)
         lower(2,1) *= -3;

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 9UL );
         checkNonZeros( lower, 6UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 3UL );

         if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) != 0 ||
             lower(1,0) != -2 || lower(1,1) !=  1 || lower(1,2) != 0 ||
             lower(2,0) !=  3 || lower(2,1) != -6 || lower(2,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n(  1  0  0 )\n( -2  1  0 )\n(  3 -6  1 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Dividing the lower element (2,1)
         lower(2,1) /= 2;

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 9UL );
         checkNonZeros( lower, 6UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 3UL );

         if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) != 0 ||
             lower(1,0) != -2 || lower(1,1) !=  1 || lower(1,2) != 0 ||
             lower(2,0) !=  3 || lower(2,1) != -3 || lower(2,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n(  1  0  0 )\n( -2  1  0 )\n(  3 -3  1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Failure cases
      {
         LT lower( 3UL );

         // Trying to write the diagonal element (1,1)
         try {
            lower(1,1) = 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to write the upper element (1,2)
         try {
            lower(1,2) = 2;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to upper matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to write the diagonal element (2,2)
         try {
            lower(2,2) = lower(1,1);

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to write the upper element (0,1)
         try {
            lower(0,1) = lower(2,1);

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to upper matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to add to the diagonal element (1,1)
         try {
            lower(1,1) += 6;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Addition assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to add to the upper element (0,2)
         try {
            lower(0,2) += 3;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Addition assignment to upper matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to subtract from the diagonal element (1,1)
         try {
            lower(1,1) -= 8;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subtraction assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to subtract from the upper element (0,1)
         try {
            lower(0,1) -= 4;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subtraction assignment to upper matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to multiply the diagonal element (1,1)
         try {
            lower(1,1) *= -6;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Multiplication assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to multiply the upper element (1,2)
         try {
            lower(1,2) *= -3;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Multiplication assignment to upper matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to divide the diagonal element (1,1)
         try {
            lower(1,1) /= 4;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Division assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to divide the upper element (1,2)
         try {
            lower(1,2) /= 2;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Division assignment to upper matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major UniLowerMatrix::operator()";

      // Good cases
      {
         OLT lower( 3UL );

         // Writing the lower element (2,1)
         lower(2,1) = 2;

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 9UL );
         checkNonZeros( lower, 4UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 1UL );

         if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
             lower(1,0) != 0 || lower(1,1) != 1 || lower(1,2) != 0 ||
             lower(2,0) != 0 || lower(2,1) != 2 || lower(2,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 )\n( 0 1 0 )\n( 0 2 1 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Writing the lower element (1,0)
         lower(1,0) = lower(2,1);

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 9UL );
         checkNonZeros( lower, 5UL );
         checkNonZeros( lower, 0UL, 2UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 1UL );

         if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
             lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 ||
             lower(2,0) != 0 || lower(2,1) != 2 || lower(2,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 )\n( 2 1 0 )\n( 0 2 1 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Adding to the lower element (2,0)
         lower(2,0) += 3;

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 9UL );
         checkNonZeros( lower, 6UL );
         checkNonZeros( lower, 0UL, 3UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 1UL );

         if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
             lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 ||
             lower(2,0) != 3 || lower(2,1) != 2 || lower(2,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 )\n( 2 1 0 )\n( 3 2 1 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Subtracting from the lower element (1,0)
         lower(1,0) -= 4;

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 9UL );
         checkNonZeros( lower, 6UL );
         checkNonZeros( lower, 0UL, 3UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 1UL );

         if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
             lower(1,0) != -2 || lower(1,1) != 1 || lower(1,2) != 0 ||
             lower(2,0) !=  3 || lower(2,1) != 2 || lower(2,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n(  1 0 0 )\n( -2 1 0 )\n(  3 2 1 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Multiplying the lower element (2,1)
         lower(2,1) *= -3;

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 9UL );
         checkNonZeros( lower, 6UL );
         checkNonZeros( lower, 0UL, 3UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 1UL );

         if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) != 0 ||
             lower(1,0) != -2 || lower(1,1) !=  1 || lower(1,2) != 0 ||
             lower(2,0) !=  3 || lower(2,1) != -6 || lower(2,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n(  1  0  0 )\n( -2  1  0 )\n(  3 -6  1 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Dividing the lower element (2,1)
         lower(2,1) /= 2;

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 9UL );
         checkNonZeros( lower, 6UL );
         checkNonZeros( lower, 0UL, 3UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 1UL );

         if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) != 0 ||
             lower(1,0) != -2 || lower(1,1) !=  1 || lower(1,2) != 0 ||
             lower(2,0) !=  3 || lower(2,1) != -3 || lower(2,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n(  1  0  0 )\n( -2  1  0 )\n(  3 -3  1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Failure cases
      {
         OLT lower( 3UL );

         // Trying to write the diagonal element (1,1)
         try {
            lower(1,1) = 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to write the upper element (1,2)
         try {
            lower(1,2) = 2;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to upper matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to write the diagonal element (2,2)
         try {
            lower(2,2) = lower(1,1);

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to write the upper element (0,1)
         try {
            lower(0,1) = lower(2,1);

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to upper matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to add to the diagonal element (1,1)
         try {
            lower(1,1) += 6;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Addition assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to add to the upper element (0,2)
         try {
            lower(0,2) += 3;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Addition assignment to upper matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to subtract from the diagonal element (1,1)
         try {
            lower(1,1) -= 8;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subtraction assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to subtract from the upper element (0,1)
         try {
            lower(0,1) -= 4;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subtraction assignment to upper matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to multiply the diagonal element (1,1)
         try {
            lower(1,1) *= -6;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Multiplication assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to multiply the upper element (1,2)
         try {
            lower(1,2) *= -3;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Multiplication assignment to upper matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to divide the diagonal element (1,1)
         try {
            lower(1,1) /= 4;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Division assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to divide the upper element (1,2)
         try {
            lower(1,2) /= 2;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Division assignment to upper matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the UniLowerMatrix iterator implementation.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the iterator implementation of the UniLowerMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testIterator()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      using Iterator      = LT::Iterator;
      using ConstIterator = LT::ConstIterator;

      LT lower( 3UL );
      lower(2,0) = 3;

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

         ConstIterator it( begin( lower, 1UL ) );

         if( it == end( lower, 1UL ) || *it != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Failed iterator conversion detected\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Counting the number of elements in 0th row via Iterator (end-begin)
      {
         test_ = "Row-major Iterator subtraction (end-begin)";

         const ptrdiff_t number( end( lower, 0UL ) - begin( lower, 0UL ) );

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

      // Counting the number of elements in 0th row via Iterator (begin-end)
      {
         test_ = "Row-major Iterator subtraction (begin-end)";

         const ptrdiff_t number( begin( lower, 0UL ) - end( lower, 0UL ) );

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

      // Counting the number of elements in 1st row via ConstIterator (end-begin)
      {
         test_ = "Row-major ConstIterator subtraction (end-begin)";

         const ptrdiff_t number( cend( lower, 1UL ) - cbegin( lower, 1UL ) );

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

      // Counting the number of elements in 1st row via ConstIterator (begin-end)
      {
         test_ = "Row-major ConstIterator subtraction (begin-end)";

         const ptrdiff_t number( cbegin( lower, 1UL ) - cend( lower, 1UL ) );

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
         test_ = "Row-major read-only access via ConstIterator";

         ConstIterator it ( cbegin( lower, 2UL ) );
         ConstIterator end( cend( lower, 2UL ) );

         if( it == end || *it != 3 ) {
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

         if( it == end || *it != 3 ) {
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

         if( it == end || *it != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator post-decrement failed\n";
            throw std::runtime_error( oss.str() );
         }

         it += 2UL;

         if( it == end || *it != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator addition assignment failed\n";
            throw std::runtime_error( oss.str() );
         }

         it -= 2UL;

         if( it == end || *it != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator subtraction assignment failed\n";
            throw std::runtime_error( oss.str() );
         }

         it = it + 2UL;

         if( it == end || *it != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator/scalar addition failed\n";
            throw std::runtime_error( oss.str() );
         }

         it = it - 2UL;

         if( it == end || *it != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator/scalar subtraction failed\n";
            throw std::runtime_error( oss.str() );
         }

         it = 3UL + it;

         if( it != end ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Scalar/iterator addition failed\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing assignment to lower elements via Iterator
      {
         test_ = "Row-major assignment to lower elements via Iterator";

         Iterator it = begin( lower, 2UL );
         *it = 2;
         ++it;
         *it = -2;

         if( lower(0,0) != 1 || lower(0,1) !=  0 || lower(0,2) != 0 ||
             lower(1,0) != 0 || lower(1,1) !=  1 || lower(1,2) != 0 ||
             lower(2,0) != 2 || lower(2,1) != -2 || lower(2,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1  0  0 )\n( 0  1  0 )\n( 2 -2  1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing assignment to diagonal elements via Iterator
      {
         test_ = "Row-major assignment to diagonal elements via Iterator";

         try {
            const Iterator it = begin( lower, 0UL );
            *it = 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Testing assignment to upper elements via Iterator
      {
         test_ = "Row-major assignment to upper elements via Iterator";

         try {
            const Iterator it = begin( lower, 0UL ) + 1UL;
            *it = 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to upper matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Testing addition assignment to lower elements via Iterator
      {
         test_ = "Row-major addition assignment to lower elements via Iterator";

         Iterator it = begin( lower, 2UL );
         *it += 2;
         ++it;
         *it += -2;

         if( lower(0,0) != 1 || lower(0,1) !=  0 || lower(0,2) != 0 ||
             lower(1,0) != 0 || lower(1,1) !=  1 || lower(1,2) != 0 ||
             lower(2,0) != 4 || lower(2,1) != -4 || lower(2,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1  0  0 )\n( 0  1  0 )\n( 4 -4  1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing addition assignment to diagonal elements via Iterator
      {
         test_ = "Row-major addition assignment to diagonal elements via Iterator";

         try {
            const Iterator it = begin( lower, 1UL ) + 1UL;
            *it += 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Testing addition assignment to upper elements via Iterator
      {
         test_ = "Row-major addition assignment to upper elements via Iterator";

         try {
            const Iterator it = begin( lower, 1UL ) + 2UL;
            *it += 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to upper matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Testing subtraction assignment to lower elements via Iterator
      {
         test_ = "Row-major subtraction assignment to lower elements via Iterator";

         Iterator it = begin( lower, 2UL );
         *it -= 2;
         ++it;
         *it -= -2;

         if( lower(0,0) != 1 || lower(0,1) !=  0 || lower(0,2) != 0 ||
             lower(1,0) != 0 || lower(1,1) !=  1 || lower(1,2) != 0 ||
             lower(2,0) != 2 || lower(2,1) != -2 || lower(2,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1  0  0 )\n( 0  1  0 )\n( 2 -2  1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing subtraction assignment to diagonal elements via Iterator
      {
         test_ = "Row-major subtraction assignment to diagonal elements via Iterator";

         try {
            const Iterator it = begin( lower, 0UL );
            *it -= 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Testing subtraction assignment to upper elements via Iterator
      {
         test_ = "Row-major subtraction assignment to upper elements via Iterator";

         try {
            const Iterator it = begin( lower, 0UL ) + 1UL;
            *it -= 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to upper matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Testing multiplication assignment to lower elements via Iterator
      {
         test_ = "Row-major multiplication assignment to lower elements via Iterator";

         Iterator it = begin( lower, 2UL );
         *it *= 2;
         ++it;
         *it *= 2;

         if( lower(0,0) != 1 || lower(0,1) !=  0 || lower(0,2) != 0 ||
             lower(1,0) != 0 || lower(1,1) !=  1 || lower(1,2) != 0 ||
             lower(2,0) != 4 || lower(2,1) != -4 || lower(2,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1  0  0 )\n( 0  1  0 )\n( 4 -4  1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing multiplication assignment to diagonal elements via Iterator
      {
         test_ = "Row-major multiplication assignment to diagonal elements via Iterator";

         try {
            const Iterator it = begin( lower, 1UL ) + 1UL;
            *it *= 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Testing multiplication assignment to upper elements via Iterator
      {
         test_ = "Row-major multiplication assignment to upper elements via Iterator";

         try {
            const Iterator it = begin( lower, 1UL ) + 2UL;
            *it *= 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to upper matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Testing division assignment to lower elements via Iterator
      {
         test_ = "Row-major division assignment to lower elements via Iterator";

         Iterator it = begin( lower, 2UL );
         *it /= 2;
         ++it;
         *it /= 2;

         if( lower(0,0) != 1 || lower(0,1) !=  0 || lower(0,2) != 0 ||
             lower(1,0) != 0 || lower(1,1) !=  1 || lower(1,2) != 0 ||
             lower(2,0) != 2 || lower(2,1) != -2 || lower(2,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1  0  0 )\n( 0  1  0 )\n( 2 -2  1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing division assignment to diagonal elements via Iterator
      {
         test_ = "Row-major division assignment to diagonal elements via Iterator";

         try {
            const Iterator it = begin( lower, 0UL );
            *it /= 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Testing division assignment to upper elements via Iterator
      {
         test_ = "Row-major division assignment to upper elements via Iterator";

         try {
            const Iterator it = begin( lower, 0UL ) + 1UL;
            *it /= 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to upper matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      using Iterator      = OLT::Iterator;
      using ConstIterator = OLT::ConstIterator;

      OLT lower( 3UL );
      lower(2,0) = 3;

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

         ConstIterator it( begin( lower, 1UL ) );

         if( it == end( lower, 1UL ) || *it != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Failed iterator conversion detected\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Counting the number of elements in 0th column via Iterator (end-begin)
      {
         test_ = "Column-major Iterator subtraction (end-begin)";

         const ptrdiff_t number( end( lower, 0UL ) - begin( lower, 0UL ) );

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

      // Counting the number of elements in 0th column via Iterator (begin-end)
      {
         test_ = "Column-major Iterator subtraction (begin-end)";

         const ptrdiff_t number( begin( lower, 0UL ) - end( lower, 0UL ) );

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

      // Counting the number of elements in 1st column via ConstIterator (end-begin)
      {
         test_ = "Column-major ConstIterator subtraction (end-begin)";

         const ptrdiff_t number( cend( lower, 1UL ) - cbegin( lower, 1UL ) );

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

      // Counting the number of elements in 1st column via ConstIterator (begin-end)
      {
         test_ = "Column-major ConstIterator subtraction (begin-end)";

         const ptrdiff_t number( cbegin( lower, 1UL ) - cend( lower, 1UL ) );

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
         test_ = "Column-major read-only access via ConstIterator";

         ConstIterator it ( cbegin( lower, 0UL ) );
         ConstIterator end( cend( lower, 0UL ) );

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

         if( it == end || *it != 3 ) {
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

         it = it + 2UL;

         if( it == end || *it != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator/scalar addition failed\n";
            throw std::runtime_error( oss.str() );
         }

         it = it - 2UL;

         if( it == end || *it != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator/scalar subtraction failed\n";
            throw std::runtime_error( oss.str() );
         }

         it = 3UL + it;

         if( it != end ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Scalar/iterator addition failed\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing assignment to lower elements via Iterator
      {
         test_ = "Column-major assignment to lower elements via Iterator";

         Iterator it = begin( lower, 0UL );
         ++it;
         *it = -2;
         ++it;
         *it = 2;

         if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
             lower(1,0) != -2 || lower(1,1) != 1 || lower(1,2) != 0 ||
             lower(2,0) !=  2 || lower(2,1) != 0 || lower(2,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n(  1 0 0 )\n( -2 1 0 )\n(  2 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing assignment to diagonal elements via Iterator
      {
         test_ = "Column-major assignment to diagonal elements via Iterator";

         try {
            const Iterator it = begin( lower, 1UL ) + 1UL;
            *it = 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Testing assignment to upper elements via Iterator
      {
         test_ = "Column-major assignment to upper elements via Iterator";

         try {
            const Iterator it = begin( lower, 1UL );
            *it = 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to upper matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Testing addition assignment to lower elements via Iterator
      {
         test_ = "Column-major addition assignment to lower elements via Iterator";

         Iterator it = begin( lower, 0UL );
         ++it;
         *it += -2;
         ++it;
         *it += 2;

         if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
             lower(1,0) != -4 || lower(1,1) != 1 || lower(1,2) != 0 ||
             lower(2,0) !=  4 || lower(2,1) != 0 || lower(2,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n(  1 0 0 )\n( -4 1 0 )\n(  4 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing addition assignment to diagonal elements via Iterator
      {
         test_ = "Column-major addition assignment to diagonal elements via Iterator";

         try {
            const Iterator it = begin( lower, 1UL ) + 1UL;
            *it += 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Testing addition assignment to upper elements via Iterator
      {
         test_ = "Column-major addition assignment to upper elements via Iterator";

         try {
            const Iterator it = begin( lower, 1UL );
            *it += 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to upper matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Testing subtraction assignment to lower elements via Iterator
      {
         test_ = "Column-major subtraction assignment to lower elements via Iterator";

         Iterator it = begin( lower, 0UL );
         ++it;
         *it -= -2;
         ++it;
         *it -= 2;

         if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
             lower(1,0) != -2 || lower(1,1) != 1 || lower(1,2) != 0 ||
             lower(2,0) !=  2 || lower(2,1) != 0 || lower(2,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n(  1 0 0 )\n( -2 1 0 )\n(  2 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing subtraction assignment to diagonal elements via Iterator
      {
         test_ = "Column-major subtraction assignment to diagonal elements via Iterator";

         try {
            const Iterator it = begin( lower, 1UL ) + 1UL;
            *it -= 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Testing subtraction assignment to upper elements via Iterator
      {
         test_ = "Column-major subtraction assignment to upper elements via Iterator";

         try {
            const Iterator it = begin( lower, 1UL );
            *it -= 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to upper matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Testing multiplication assignment to lower elements via Iterator
      {
         test_ = "Column-major multiplication assignment to lower elements via Iterator";

         Iterator it = begin( lower, 0UL );
         ++it;
         *it *= 2;
         ++it;
         *it *= 2;

         if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
             lower(1,0) != -4 || lower(1,1) != 1 || lower(1,2) != 0 ||
             lower(2,0) !=  4 || lower(2,1) != 0 || lower(2,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n(  1 0 0 )\n( -4 1 0 )\n(  4 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing multiplication assignment to diagonal elements via Iterator
      {
         test_ = "Column-major multiplication assignment to diagonal elements via Iterator";

         try {
            const Iterator it = begin( lower, 1UL ) + 1UL;
            *it *= 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Testing multiplication assignment to upper elements via Iterator
      {
         test_ = "Column-major multiplication assignment to upper elements via Iterator";

         try {
            const Iterator it = begin( lower, 1UL );
            *it *= 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to upper matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Testing division assignment to lower elements via Iterator
      {
         test_ = "Column-major division assignment to lower elements via Iterator";

         Iterator it = begin( lower, 0UL );
         ++it;
         *it /= 2;
         ++it;
         *it /= 2;

         if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
             lower(1,0) != -2 || lower(1,1) != 1 || lower(1,2) != 0 ||
             lower(2,0) !=  2 || lower(2,1) != 0 || lower(2,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n(  1 0 0 )\n( -2 1 0 )\n(  2 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing division assignment to diagonal elements via Iterator
      {
         test_ = "Column-major division assignment to diagonal elements via Iterator";

         try {
            const Iterator it = begin( lower, 1UL ) + 1UL;
            *it /= 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Testing division assignment to upper elements via Iterator
      {
         test_ = "Column-major division assignment to upper elements via Iterator";

         try {
            const Iterator it = begin( lower, 1UL );
            *it /= 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to upper matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c nonZeros() member function of the UniLowerMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c nonZeros() member function of the UniLowerMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testNonZeros()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major UniLowerMatrix::nonZeros()";

      // Default matrix
      {
         LT lower( 3UL );

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 9UL );
         checkNonZeros( lower, 3UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 1UL );
         checkNonZeros( lower, 2UL, 1UL );

         if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
             lower(1,0) != 0 || lower(1,1) != 1 || lower(1,2) != 0 ||
             lower(2,0) != 0 || lower(2,1) != 0 || lower(2,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 )\n( 0 1 0 )\n( 0 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Fully filled matrix
      {
         LT lower( 3UL );
         lower(1,0) =  2;
         lower(2,0) = -4;
         lower(2,1) = -5;

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 9UL );
         checkNonZeros( lower, 6UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 3UL );

         if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) != 0 ||
             lower(1,0) !=  2 || lower(1,1) !=  1 || lower(1,2) != 0 ||
             lower(2,0) != -4 || lower(2,1) != -5 || lower(2,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n(  1  0  0 )\n(  2  1  0 )\n( -4 -5  1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major UniLowerMatrix::nonZeros()";

      // Default matrix
      {
         OLT lower( 3UL );

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 9UL );
         checkNonZeros( lower, 3UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 1UL );
         checkNonZeros( lower, 2UL, 1UL );

         if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
             lower(1,0) != 0 || lower(1,1) != 1 || lower(1,2) != 0 ||
             lower(2,0) != 0 || lower(2,1) != 0 || lower(2,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 )\n( 0 1 0 )\n( 0 0 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Fully filled matrix
      {
         OLT lower( 3UL );
         lower(1,0) =  2;
         lower(2,0) = -4;
         lower(2,1) = -5;

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 9UL );
         checkNonZeros( lower, 6UL );
         checkNonZeros( lower, 0UL, 3UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 1UL );

         if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) != 0 ||
             lower(1,0) !=  2 || lower(1,1) !=  1 || lower(1,2) != 0 ||
             lower(2,0) != -4 || lower(2,1) != -5 || lower(2,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n(  1  0  0 )\n(  2  1  0 )\n( -4 -5  1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c reset() member function of the UniLowerMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c reset() member function of the UniLowerMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testReset()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major UniLowerMatrix::reset()";

      // Initialization check
      LT lower( 3UL );
      lower(1,0) = 2;
      lower(2,0) = 4;
      lower(2,1) = 5;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 6UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 3UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) != 4 || lower(2,1) != 5 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 2 1 0 )\n( 4 5 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting a lower element
      reset( lower(2,0) );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 2UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) != 0 || lower(2,1) != 5 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 2 1 0 )\n( 0 5 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting a diagonal element
      reset( lower(1,1) );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 2UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) != 0 || lower(2,1) != 5 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 2 1 0 )\n( 0 5 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting an upper element
      reset( lower(0,2) );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 2UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) != 0 || lower(2,1) != 5 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 2 1 0 )\n( 0 5 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting row 1
      reset( lower, 1UL );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 4UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 2UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 0 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) != 0 || lower(2,1) != 5 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 1 0 )\n( 0 5 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting the entire matrix
      reset( lower );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 3UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 0 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) != 0 || lower(2,1) != 0 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 1 0 )\n( 0 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major UniLowerMatrix::reset()";

      // Initialization check
      OLT lower( 3UL );
      lower(1,0) = 2;
      lower(2,0) = 4;
      lower(2,1) = 5;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 6UL );
      checkNonZeros( lower, 0UL, 3UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) != 4 || lower(2,1) != 5 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 2 1 0 )\n( 4 5 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting a lower element
      reset( lower(2,0) );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 2UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) != 0 || lower(2,1) != 5 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 2 1 0 )\n( 0 5 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting a diagonal element
      reset( lower(1,1) );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 2UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) != 0 || lower(2,1) != 5 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 2 1 0 )\n( 0 5 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting an upper element
      reset( lower(0,2) );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 2UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) != 0 || lower(2,1) != 5 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 2 1 0 )\n( 0 5 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting column 1
      reset( lower, 1UL );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 4UL );
      checkNonZeros( lower, 0UL, 2UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) != 0 || lower(2,1) != 0 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 2 1 0 )\n( 0 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting the entire matrix
      reset( lower );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 3UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 0 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) != 0 || lower(2,1) != 0 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 1 0 )\n( 0 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c clear() member function of the UniLowerMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c clear() member function of the UniLowerMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testClear()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major UniLowerMatrix::clear() (StaticMatrix)";

      // Initialization check
      blaze::UniLowerMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> > lower;
      lower(1,0) = 2;
      lower(2,0) = 4;
      lower(2,1) = 5;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 6UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 3UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) != 4 || lower(2,1) != 5 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 2 1 0 )\n( 4 5 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing a lower element
      clear( lower(2,0) );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 2UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) != 0 || lower(2,1) != 5 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 2 1 0 )\n( 0 5 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing a diagonal element
      clear( lower(1,1) );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 2UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) != 0 || lower(2,1) != 5 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 2 1 0 )\n( 0 5 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing an upper element
      clear( lower(0,2) );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 2UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) != 0 || lower(2,1) != 5 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 2 1 0 )\n( 0 5 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing the matrix
      clear( lower );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 3UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 0 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) != 0 || lower(2,1) != 0 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 1 0 )\n( 0 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major UniLowerMatrix::clear() (DynamicMatrix)";

      // Initialization check
      LT lower( 3UL );
      lower(1,0) = 2;
      lower(2,0) = 4;
      lower(2,1) = 5;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 6UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 3UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) != 4 || lower(2,1) != 5 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 2 1 0 )\n( 4 5 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing a lower element
      clear( lower(2,0) );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 2UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) != 0 || lower(2,1) != 5 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 2 1 0 )\n( 0 5 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing a diagonal element
      clear( lower(1,1) );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 2UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) != 0 || lower(2,1) != 5 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 2 1 0 )\n( 0 5 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing an upper element
      clear( lower(0,2) );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 2UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) != 0 || lower(2,1) != 5 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 2 1 0 )\n( 0 5 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing the matrix
      clear( lower );

      checkRows    ( lower, 0UL );
      checkColumns ( lower, 0UL );
      checkNonZeros( lower, 0UL );
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major UniLowerMatrix::clear() (StaticMatrix)";

      // Initialization check
      blaze::UniLowerMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> > lower;
      lower(1,0) = 2;
      lower(2,0) = 4;
      lower(2,1) = 5;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 6UL );
      checkNonZeros( lower, 0UL, 3UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) != 4 || lower(2,1) != 5 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 2 1 0 )\n( 4 5 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing a lower element
      clear( lower(2,0) );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 2UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) != 0 || lower(2,1) != 5 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 2 1 0 )\n( 0 5 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing a diagonal element
      clear( lower(1,1) );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 2UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) != 0 || lower(2,1) != 5 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 2 1 0 )\n( 0 5 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing an upper element
      clear( lower(0,2) );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 2UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) != 0 || lower(2,1) != 5 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 2 1 0 )\n( 0 5 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing the matrix
      clear( lower );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 3UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 0 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) != 0 || lower(2,1) != 0 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 1 0 )\n( 0 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major UniLowerMatrix::clear() (DynamicMatrix)";

      // Initialization check
      OLT lower( 3UL );
      lower(1,0) = 2;
      lower(2,0) = 4;
      lower(2,1) = 5;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 6UL );
      checkNonZeros( lower, 0UL, 3UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) != 4 || lower(2,1) != 5 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 2 1 0 )\n( 4 5 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing a lower element
      clear( lower(2,0) );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 2UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) != 0 || lower(2,1) != 5 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 2 1 0 )\n( 0 5 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing a diagonal element
      clear( lower(1,1) );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 2UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) != 0 || lower(2,1) != 5 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 2 1 0 )\n( 0 5 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing an upper element
      clear( lower(0,2) );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 2UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) != 0 || lower(2,1) != 5 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 2 1 0 )\n( 0 5 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing the matrix
      clear( lower );

      checkRows    ( lower, 0UL );
      checkColumns ( lower, 0UL );
      checkNonZeros( lower, 0UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c resize() member function of the UniLowerMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c resize() member function of the UniLowerMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testResize()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major UniLowerMatrix::resize()";

      // Initialization check
      LT lower;

      checkRows    ( lower, 0UL );
      checkColumns ( lower, 0UL );
      checkNonZeros( lower, 0UL );

      // Resizing to 2x2
      lower.resize( 2UL );

      checkRows    ( lower, 2UL );
      checkColumns ( lower, 2UL );
      checkCapacity( lower, 4UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(1,1) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 )\n( x 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 4x4 and preserving the elements
      lower(1,0) = 2;
      lower.resize( 4UL, true );

      checkRows    ( lower,  4UL );
      checkColumns ( lower,  4UL );
      checkCapacity( lower, 16UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 || lower(0,3) != 0 ||
          lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 || lower(1,3) != 0 ||
                                                lower(2,2) != 1 || lower(2,3) != 0 ||
                                                                   lower(3,3) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 2 1 0 0 )\n( x x 1 0 )\n( x x x 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 2x2
      lower(2,1) = 4;
      lower.resize( 2UL );

      checkRows    ( lower, 2UL );
      checkColumns ( lower, 2UL );
      checkCapacity( lower, 4UL );
      checkNonZeros( lower, 3UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 ||
          lower(1,0) != 2 || lower(1,1) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 )\n( 2 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 0x0
      lower.resize( 0UL );

      checkRows    ( lower, 0UL );
      checkColumns ( lower, 0UL );
      checkNonZeros( lower, 0UL );
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major UniLowerMatrix::resize()";

      // Initialization check
      OLT lower;

      checkRows    ( lower, 0UL );
      checkColumns ( lower, 0UL );
      checkNonZeros( lower, 0UL );

      // Resizing to 2x2
      lower.resize( 2UL );

      checkRows    ( lower, 2UL );
      checkColumns ( lower, 2UL );
      checkCapacity( lower, 4UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(1,1) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 )\n( x 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 4x4 and preserving the elements
      lower(1,0) = 2;
      lower.resize( 4UL, true );

      checkRows    ( lower,  4UL );
      checkColumns ( lower,  4UL );
      checkCapacity( lower, 16UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 || lower(0,3) != 0 ||
          lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 || lower(1,3) != 0 ||
                                                lower(2,2) != 1 || lower(2,3) != 0 ||
                                                                   lower(3,3) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 2 1 0 0 )\n( x x 1 0 )\n( x x x 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 2x2
      lower(2,1) = 4;
      lower.resize( 2UL );

      checkRows    ( lower, 2UL );
      checkColumns ( lower, 2UL );
      checkCapacity( lower, 4UL );
      checkNonZeros( lower, 3UL );
      checkNonZeros( lower, 0UL, 2UL );
      checkNonZeros( lower, 1UL, 1UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 ||
          lower(1,0) != 2 || lower(1,1) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 )\n( 2 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 0x0
      lower.resize( 0UL );

      checkRows    ( lower, 0UL );
      checkColumns ( lower, 0UL );
      checkNonZeros( lower, 0UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c extend() member function of the UniLowerMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c extend() member function of the UniLowerMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testExtend()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major UniLowerMatrix::extend()";

      // Initialization check
      LT lower;

      checkRows    ( lower, 0UL );
      checkColumns ( lower, 0UL );
      checkNonZeros( lower, 0UL );

      // Extending the size of the matrix to 2x2
      lower.extend( 2UL );

      checkRows    ( lower, 2UL );
      checkColumns ( lower, 2UL );
      checkCapacity( lower, 4UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(1,1) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Extending the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 )\n( x 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Extending to 4x4 and preserving the elements
      lower(1,0) = 2;
      lower.extend( 2UL, true );

      checkRows    ( lower,  4UL );
      checkColumns ( lower,  4UL );
      checkCapacity( lower, 16UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 || lower(0,3) != 0 ||
          lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 || lower(1,3) != 0 ||
                                                lower(2,2) != 1 || lower(2,3) != 0 ||
                                                                   lower(3,3) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Extending the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 2 1 0 0 )\n( x x 1 0 )\n( x x x 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major UniLowerMatrix::extend()";

      // Initialization check
      OLT lower;

      checkRows    ( lower, 0UL );
      checkColumns ( lower, 0UL );
      checkNonZeros( lower, 0UL );

      // Extending the size of the matrix to 2x2
      lower.extend( 2UL );

      checkRows    ( lower, 2UL );
      checkColumns ( lower, 2UL );
      checkCapacity( lower, 4UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(1,1) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Extending the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 )\n( x 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Extending to 4x4 and preserving the elements
      lower(1,0) = 2;
      lower.extend( 2UL, true );

      checkRows    ( lower,  4UL );
      checkColumns ( lower,  4UL );
      checkCapacity( lower, 16UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 || lower(0,3) != 0 ||
          lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 || lower(1,3) != 0 ||
                                                lower(2,2) != 1 || lower(2,3) != 0 ||
                                                                   lower(3,3) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Extending the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 2 1 0 0 )\n( x x 1 0 )\n( x x x 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c reserve() member function of the UniLowerMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c reserve() member function of the UniLowerMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testReserve()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major UniLowerMatrix::reserve()";

      // Initialization check
      LT lower;

      checkRows    ( lower, 0UL );
      checkColumns ( lower, 0UL );
      checkNonZeros( lower, 0UL );

      // Increasing the capacity of the matrix
      lower.reserve( 10UL );

      checkRows    ( lower,  0UL );
      checkColumns ( lower,  0UL );
      checkCapacity( lower, 10UL );
      checkNonZeros( lower,  0UL );

      // Further increasing the capacity of the matrix
      lower.reserve( 20UL );

      checkRows    ( lower,  0UL );
      checkColumns ( lower,  0UL );
      checkCapacity( lower, 20UL );
      checkNonZeros( lower,  0UL );
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major UniLowerMatrix::reserve()";

      // Initialization check
      OLT lower;

      checkRows    ( lower, 0UL );
      checkColumns ( lower, 0UL );
      checkNonZeros( lower, 0UL );

      // Increasing the capacity of the matrix
      lower.reserve( 10UL );

      checkRows    ( lower,  0UL );
      checkColumns ( lower,  0UL );
      checkCapacity( lower, 10UL );
      checkNonZeros( lower,  0UL );

      // Further increasing the capacity of the matrix
      lower.reserve( 20UL );

      checkRows    ( lower,  0UL );
      checkColumns ( lower,  0UL );
      checkCapacity( lower, 20UL );
      checkNonZeros( lower,  0UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c shrinkToFit() member function of the UniLowerMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c shrinkToFit() member function of the UniLowerMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testShrinkToFit()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major UniLowerMatrix::shrinkToFit()";

      // Shrinking a matrix without excessive capacity
      {
         LT lower{ { 1, 0, 0 }, { 2, 1, 0 }, { 3, 4, 1 } };

         lower.shrinkToFit();

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 9UL );
         checkNonZeros( lower, 6UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 3UL );

         if( lower.capacity() != lower.rows() * lower.spacing() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Shrinking the matrix failed\n"
                << " Details:\n"
                << "   Capacity         : " << lower.capacity() << "\n"
                << "   Expected capacity: " << ( lower.rows() * lower.spacing() ) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
             lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 ||
             lower(2,0) != 3 || lower(2,1) != 4 || lower(2,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Shrinking the matrix failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 )\n( 2 1 0 )\n( 3 4 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Shrinking a matrix with excessive capacity
      {
         LT lower{ { 1, 0, 0 }, { 2, 1, 0 }, { 3, 4, 1 } };
         lower.reserve( 100UL );

         lower.shrinkToFit();

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 9UL );
         checkNonZeros( lower, 6UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 3UL );

         if( lower.capacity() != lower.rows() * lower.spacing() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Shrinking the matrix failed\n"
                << " Details:\n"
                << "   Capacity         : " << lower.capacity() << "\n"
                << "   Expected capacity: " << ( lower.rows() * lower.spacing() ) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
             lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 ||
             lower(2,0) != 3 || lower(2,1) != 4 || lower(2,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Shrinking the matrix failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 )\n( 2 1 0 )\n( 3 4 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major UniLowerMatrix::shrinkToFit()";

      // Shrinking a matrix without excessive capacity
      {
         OLT lower{ { 1, 0, 0 }, { 2, 1, 0 }, { 3, 4, 1 } };

         lower.shrinkToFit();

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 9UL );
         checkNonZeros( lower, 6UL );
         checkNonZeros( lower, 0UL, 3UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 1UL );

         if( lower.capacity() != lower.spacing() * lower.columns() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Shrinking the matrix failed\n"
                << " Details:\n"
                << "   Capacity         : " << lower.capacity() << "\n"
                << "   Expected capacity: " << ( lower.spacing() * lower.columns() ) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
             lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 ||
             lower(2,0) != 3 || lower(2,1) != 4 || lower(2,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Shrinking the matrix failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 )\n( 2 1 0 )\n( 3 4 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Shrinking a matrix with excessive capacity
      {
         OLT lower{ { 1, 0, 0 }, { 2, 1, 0 }, { 3, 4, 1 } };
         lower.reserve( 100UL );

         lower.shrinkToFit();

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 9UL );
         checkNonZeros( lower, 6UL );
         checkNonZeros( lower, 0UL, 3UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 1UL );

         if( lower.capacity() != lower.spacing() * lower.columns() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Shrinking the matrix failed\n"
                << " Details:\n"
                << "   Capacity         : " << lower.capacity() << "\n"
                << "   Expected capacity: " << ( lower.spacing() * lower.columns() ) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
             lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 ||
             lower(2,0) != 3 || lower(2,1) != 4 || lower(2,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Shrinking the matrix failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 )\n( 2 1 0 )\n( 3 4 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c swap() functionality of the UniLowerMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c swap() function of the UniLowerMatrix specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testSwap()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major UniLowerMatrix swap";

      LT lower1( 2UL );
      lower1(1,0) = 2;

      LT lower2( 3UL );
      lower2(1,0) = 3;
      lower2(2,0) = 4;
      lower2(2,1) = 5;

      swap( lower1, lower2 );

      checkRows    ( lower1, 3UL );
      checkColumns ( lower1, 3UL );
      checkCapacity( lower1, 9UL );
      checkNonZeros( lower1, 6UL );
      checkNonZeros( lower1, 0UL, 1UL );
      checkNonZeros( lower1, 1UL, 2UL );
      checkNonZeros( lower1, 2UL, 3UL );

      if( lower1(0,0) != 1 || lower1(0,1) != 0 || lower1(0,2) != 0 ||
          lower1(1,0) != 3 || lower1(1,1) != 1 || lower1(1,2) != 0 ||
          lower1(2,0) != 4 || lower1(2,1) != 5 || lower1(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the first matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower1 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 3 1 0 )\n( 4 5 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      checkRows    ( lower2, 2UL );
      checkColumns ( lower2, 2UL );
      checkCapacity( lower2, 4UL );
      checkNonZeros( lower2, 3UL );
      checkNonZeros( lower2, 0UL, 1UL );
      checkNonZeros( lower2, 1UL, 2UL );

      if( lower2(0,0) != 1 || lower2(0,1) != 0 || lower2(1,0) != 2 || lower2(1,1) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the second matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n( 1 0 )\n( 2 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major UniLowerMatrix swap";

      OLT lower1( 2UL );
      lower1(1,0) = 2;

      OLT lower2( 3UL );
      lower2(1,0) = 3;
      lower2(2,0) = 4;
      lower2(2,1) = 5;

      swap( lower1, lower2 );

      checkRows    ( lower1, 3UL );
      checkColumns ( lower1, 3UL );
      checkCapacity( lower1, 9UL );
      checkNonZeros( lower1, 6UL );
      checkNonZeros( lower1, 0UL, 3UL );
      checkNonZeros( lower1, 1UL, 2UL );
      checkNonZeros( lower1, 2UL, 1UL );

      if( lower1(0,0) != 1 || lower1(0,1) != 0 || lower1(0,2) != 0 ||
          lower1(1,0) != 3 || lower1(1,1) != 1 || lower1(1,2) != 0 ||
          lower1(2,0) != 4 || lower1(2,1) != 5 || lower1(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the first matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower1 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 3 1 0 )\n( 4 5 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      checkRows    ( lower2, 2UL );
      checkColumns ( lower2, 2UL );
      checkCapacity( lower2, 4UL );
      checkNonZeros( lower2, 3UL );
      checkNonZeros( lower2, 0UL, 2UL );
      checkNonZeros( lower2, 1UL, 1UL );

      if( lower2(0,0) != 1 || lower2(0,1) != 0 || lower2(1,0) != 2 || lower2(1,1) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the second matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n( 1 0 )\n( 2 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c isDefault() function with the UniLowerMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isDefault() function with the UniLowerMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testIsDefault()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major isDefault() function (StaticMatrix)";

      // isDefault with default matrix
      {
         blaze::UniLowerMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> > lower;

         if( isDefault( lower(1,0) ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix element:\n" << lower(1,0) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( lower(1,1) ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix element:\n" << lower(1,1) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( lower ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with non-default matrix
      {
         blaze::UniLowerMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> > lower;
         lower(1,0) = 1;

         if( isDefault( lower(1,0) ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix element:\n" << lower(1,0) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( lower(1,1) ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix element:\n" << lower(1,1) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( lower ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   {
      test_ = "Row-major isDefault() function (DynamicMatrix)";

      // isDefault with 0x0 matrix
      {
         LT lower;

         if( isDefault( lower ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with default matrix
      {
         LT lower( 3UL );

         if( isDefault( lower(1,0) ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix element:\n" << lower(1,0) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( lower(1,1) ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix element:\n" << lower(1,1) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( lower ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with non-default matrix
      {
         LT lower( 3UL );
         lower(1,0) = 1;

         if( isDefault( lower(1,0) ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix element:\n" << lower(1,0) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( lower(1,1) ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix element:\n" << lower(1,1) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( lower ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major isDefault() function (StaticMatrix)";

      // isDefault with default matrix
      {
         blaze::UniLowerMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> > lower;

         if( isDefault( lower(1,0) ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix element:\n" << lower(1,0) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( lower(1,1) ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix element:\n" << lower(1,1) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( lower ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with non-default matrix
      {
         blaze::UniLowerMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> > lower;
         lower(1,0) = 1;

         if( isDefault( lower(1,0) ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix element:\n" << lower(1,0) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( lower(1,1) ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix element:\n" << lower(1,1) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( lower ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   {
      test_ = "Column-major isDefault() function (DynamicMatrix)";

      // isDefault with 0x0 matrix
      {
         OLT lower;

         if( isDefault( lower ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with default matrix
      {
         OLT lower( 3UL );

         if( isDefault( lower(1,0) ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix element:\n" << lower(1,0) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( lower(1,1) ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix element:\n" << lower(1,1) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( lower ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with non-default matrix
      {
         OLT lower( 3UL );
         lower(1,0) = 1;

         if( isDefault( lower(1,0) ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix element:\n" << lower(1,0) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( lower(1,1) ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix element:\n" << lower(1,1) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( lower ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << lower << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the assignment to submatrices of the UniLowerMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the assignment to submatrices of the UniLowerMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testSubmatrix()
{
   //=====================================================================================
   // Row-major general tests
   //=====================================================================================

   {
      test_ = "Row-major submatrix() function";

      using SMT = blaze::Submatrix<LT>;

      LT lower( 3UL );
      lower(1,0) = -4;
      lower(2,0) =  7;

      SMT sm = submatrix( lower, 1UL, 1UL, 2UL, 2UL );

      if( sm(1,1) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator access failed\n"
             << " Details:\n"
             << "   Result: " << sm(1,1) << "\n"
             << "   Expected result: 1\n";
         throw std::runtime_error( oss.str() );
      }

      SMT::Iterator it = sm.begin(0UL);

      if( it == sm.end(0UL) || *it != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << *it << "\n"
             << "   Expected result: 1\n";
         throw std::runtime_error( oss.str() );
      }

      sm(1,0) = -5;

      if( sm(0,0) !=  1 || sm(0,1) != 0 ||
          sm(1,0) != -5 || sm(1,1) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix access failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  1  0 )\n( -5  1 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) != 0 ||
          lower(1,0) != -4 || lower(1,1) !=  1 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != -5 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix access failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0 )\n( -4  1  0 )\n(  7 -5  1 )\n";
         throw std::runtime_error( oss.str() );
      }

      reset( sm );

      if( sm(0,0) != 1 || sm(0,1) != 0 ||
          sm(1,0) != 0 || sm(1,1) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix reset failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 1 0 )\n( 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -4 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 0 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix reset failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 1 0 )\n(  7 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major scalar assignment
   //=====================================================================================

   // (  1  0  0  0 )      (  1  0  0  0 )
   // ( -4  1  0  0 )  =>  ( -4  1  0  0 )
   // (  7  0  1  0 )      (  7 12  1  0 )
   // ( -2  0  1  1 )      ( -2 12 12  1 )
   {
      test_ = "Row-major submatrix() function (scalar assignment test 1)";

      using SMT = blaze::Submatrix<LT>;

      LT lower( 4UL );
      lower(1,0) = -4;
      lower(2,0) =  7;
      lower(3,0) = -2;
      lower(3,2) =  1;

      SMT sm = submatrix( lower, 0UL, 1UL, 4UL, 2UL );
      sm = 12;

      checkRows    ( lower,  4UL );
      checkColumns ( lower,  4UL );
      checkNonZeros( lower, 10UL );
      checkNonZeros( lower,  0UL, 1UL );
      checkNonZeros( lower,  1UL, 2UL );
      checkNonZeros( lower,  2UL, 3UL );
      checkNonZeros( lower,  3UL, 4UL );

      if( sm(0,0) !=  0 || sm(0,1) !=  0 ||
          sm(1,0) !=  1 || sm(1,1) !=  0 ||
          sm(2,0) != 12 || sm(2,1) !=  1 ||
          sm(3,0) != 12 || sm(3,1) != 12 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0  0 )\n(  1  0 )\n( 12  1 )\n( 12 12 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) !=  0 || lower(0,3) != 0 ||
          lower(1,0) != -4 || lower(1,1) !=  1 || lower(1,2) !=  0 || lower(1,3) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 12 || lower(2,2) !=  1 || lower(2,3) != 0 ||
          lower(3,0) != -2 || lower(3,1) != 12 || lower(3,2) != 12 || lower(3,3) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0  0 )\n"
                                     "( -4  1  0  0 )\n"
                                     "(  7 12  1  0 )\n"
                                     "( -2 12 12  1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // (  1  0  0  0 )      (  1  0  0  0 )
   // ( -4  1  0  0 )  =>  ( 12  1  0  0 )
   // (  7  0  1  0 )      ( 12 12  1  0 )
   // ( -2  0  1  1 )      ( -2  0  1  1 )
   {
      test_ = "Row-major submatrix() function (scalar assignment test 2)";

      using SMT = blaze::Submatrix<LT>;

      LT lower( 4UL );
      lower(1,0) = -4;
      lower(2,0) =  7;
      lower(3,0) = -2;
      lower(3,2) =  1;

      SMT sm = submatrix( lower, 1UL, 0UL, 2UL, 4UL );
      sm = 12;

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkNonZeros( lower, 9UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 3UL );
      checkNonZeros( lower, 3UL, 3UL );

      if( sm(0,0) != 12 || sm(0,1) !=  1 || sm(0,2) != 0 || sm(0,3) != 0 ||
          sm(1,0) != 12 || sm(1,1) != 12 || sm(1,2) != 1 || sm(1,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 12  1  0  0 )\n( 12 12  1  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) !=  0 || lower(0,3) != 0 ||
          lower(1,0) != 12 || lower(1,1) !=  1 || lower(1,2) !=  0 || lower(1,3) != 0 ||
          lower(2,0) != 12 || lower(2,1) != 12 || lower(2,2) !=  1 || lower(2,3) != 0 ||
          lower(3,0) != -2 || lower(3,1) !=  0 || lower(3,2) !=  1 || lower(3,3) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0  0 )\n"
                                     "( 12  1  0  0 )\n"
                                     "( 12 12  1  0 )\n"
                                     "( -2  0  1  1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // (  1  0  0  0 )      (  1  0  0  0 )
   // ( -4  1  0  0 )  =>  ( -4  1  0  0 )
   // (  7  0  1  0 )      (  7  0  1  0 )
   // ( -2  0  1  1 )      ( -2  0  1  1 )
   {
      test_ = "Row-major submatrix() function (scalar assignment test 3)";

      using SMT = blaze::Submatrix<LT>;

      LT lower( 4UL );
      lower(1,0) = -4;
      lower(2,0) =  7;
      lower(3,0) = -2;
      lower(3,2) =  1;

      SMT sm = submatrix( lower, 0UL, 2UL, 2UL, 2UL );
      sm = 12;

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkNonZeros( lower, 8UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 2UL );
      checkNonZeros( lower, 3UL, 3UL );

      if( sm(0,0) != 0 || sm(0,1) != 0 ||
          sm(1,0) != 0 || sm(1,1) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0 0 )\n( 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 || lower(0,3) != 0 ||
          lower(1,0) != -4 || lower(1,1) != 1 || lower(1,2) != 0 || lower(1,3) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 0 || lower(2,2) != 1 || lower(2,3) != 0 ||
          lower(3,0) != -2 || lower(3,1) != 0 || lower(3,2) != 1 || lower(3,3) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0  0 )\n"
                                     "( -4  1  0  0 )\n"
                                     "(  7  0  1  0 )\n"
                                     "( -2  0  1  1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major general tests
   //=====================================================================================

   {
      test_ = "Column-major submatrix() function";

      using SMT = blaze::Submatrix<OLT>;

      OLT lower( 3UL );
      lower(1,0) = -4;
      lower(2,0) =  7;

      SMT sm = submatrix( lower, 1UL, 1UL, 2UL, 2UL );

      if( sm(1,1) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator access failed\n"
             << " Details:\n"
             << "   Result: " << sm(1,1) << "\n"
             << "   Expected result: 1\n";
         throw std::runtime_error( oss.str() );
      }

      SMT::Iterator it = sm.begin(0UL);

      if( it == sm.end(0UL) || *it != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << *it << "\n"
             << "   Expected result: 1\n";
         throw std::runtime_error( oss.str() );
      }

      sm(1,0) = -5;

      if( sm(0,0) !=  1 || sm(0,1) != 0 ||
          sm(1,0) != -5 || sm(1,1) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix access failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  1  0 )\n( -5  1 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) != 0 ||
          lower(1,0) != -4 || lower(1,1) !=  1 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != -5 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix access failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0 )\n( -4  1  0 )\n(  7 -5  1 )\n";
         throw std::runtime_error( oss.str() );
      }

      reset( sm );

      if( sm(0,0) != 1 || sm(0,1) != 0 ||
          sm(1,0) != 0 || sm(1,1) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix reset failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 1 0 )\n( 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -4 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 0 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix reset failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 1 0 )\n(  7 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major scalar assignment
   //=====================================================================================

   // (  1  0  0  0 )      (  1  0  0  0 )
   // ( -4  2  0  0 )  =>  ( -4  1  0  0 )
   // (  7  0  3  0 )      (  7 12  1  0 )
   // ( -2  0  1  5 )      ( -2 12 12  1 )
   {
      test_ = "Column-major submatrix() function (scalar assignment test 1)";

      using SMT = blaze::Submatrix<OLT>;

      OLT lower( 4UL );
      lower(1,0) = -4;
      lower(2,0) =  7;
      lower(3,0) = -2;
      lower(3,2) =  1;

      SMT sm = submatrix( lower, 0UL, 1UL, 4UL, 2UL );
      sm = 12;

      checkRows    ( lower,  4UL );
      checkColumns ( lower,  4UL );
      checkNonZeros( lower, 10UL );
      checkNonZeros( lower,  0UL, 4UL );
      checkNonZeros( lower,  1UL, 3UL );
      checkNonZeros( lower,  2UL, 2UL );
      checkNonZeros( lower,  3UL, 1UL );

      if( sm(0,0) !=  0 || sm(0,1) !=  0 ||
          sm(1,0) !=  1 || sm(1,1) !=  0 ||
          sm(2,0) != 12 || sm(2,1) !=  1 ||
          sm(3,0) != 12 || sm(3,1) != 12 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0  0 )\n(  1  0 )\n( 12  1 )\n( 12 12 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) !=  0 || lower(0,3) != 0 ||
          lower(1,0) != -4 || lower(1,1) !=  1 || lower(1,2) !=  0 || lower(1,3) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 12 || lower(2,2) !=  1 || lower(2,3) != 0 ||
          lower(3,0) != -2 || lower(3,1) != 12 || lower(3,2) != 12 || lower(3,3) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0  0 )\n"
                                     "( -4  1  0  0 )\n"
                                     "(  7 12  1  0 )\n"
                                     "( -2 12 12  1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // (  1  0  0  0 )      (  1  0  0  0 )
   // ( -4  1  0  0 )  =>  ( 12  1  0  0 )
   // (  7  0  1  0 )      ( 12 12  1  0 )
   // ( -2  0  1  1 )      ( -2  0  1  1 )
   {
      test_ = "Column-major submatrix() function (scalar assignment test 2)";

      using SMT = blaze::Submatrix<OLT>;

      OLT lower( 4UL );
      lower(1,0) = -4;
      lower(2,0) =  7;
      lower(3,0) = -2;
      lower(3,2) =  1;

      SMT sm = submatrix( lower, 1UL, 0UL, 2UL, 4UL );
      sm = 12;

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkNonZeros( lower, 9UL );
      checkNonZeros( lower, 0UL, 4UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 2UL );
      checkNonZeros( lower, 3UL, 1UL );

      if( sm(0,0) != 12 || sm(0,1) !=  1 || sm(0,2) != 0 || sm(0,3) != 0 ||
          sm(1,0) != 12 || sm(1,1) != 12 || sm(1,2) != 1 || sm(1,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 12  1  0  0 )\n( 12 12  1  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) != 0 || lower(0,3) != 0 ||
          lower(1,0) != 12 || lower(1,1) !=  1 || lower(1,2) != 0 || lower(1,3) != 0 ||
          lower(2,0) != 12 || lower(2,1) != 12 || lower(2,2) != 1 || lower(2,3) != 0 ||
          lower(3,0) != -2 || lower(3,1) !=  0 || lower(3,2) != 1 || lower(3,3) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0  0 )\n"
                                     "( 12  1  0  0 )\n"
                                     "( 12 12  1  0 )\n"
                                     "( -2  0  1  1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // (  1  0  0  0 )      (  1  0  0  0 )
   // ( -4  1  0  0 )  =>  ( -4  1  0  0 )
   // (  7  0  1  0 )      (  7  0  1  0 )
   // ( -2  0  1  1 )      ( -2  0  1  1 )
   {
      test_ = "Column-major submatrix() function (scalar assignment test 3)";

      using SMT = blaze::Submatrix<OLT>;

      OLT lower( 4UL );
      lower(1,0) = -4;
      lower(2,0) =  7;
      lower(3,0) = -2;
      lower(3,2) =  1;

      SMT sm = submatrix( lower, 0UL, 2UL, 2UL, 2UL );
      sm = 12;

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkNonZeros( lower, 8UL );
      checkNonZeros( lower, 0UL, 4UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 2UL );
      checkNonZeros( lower, 3UL, 1UL );

      if( sm(0,0) != 0 || sm(0,1) != 0 ||
          sm(1,0) != 0 || sm(1,1) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0 0 )\n( 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 || lower(0,3) != 0 ||
          lower(1,0) != -4 || lower(1,1) != 1 || lower(1,2) != 0 || lower(1,3) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 0 || lower(2,2) != 1 || lower(2,3) != 0 ||
          lower(3,0) != -2 || lower(3,1) != 0 || lower(3,2) != 1 || lower(3,3) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0  0 )\n"
                                     "( -4  1  0  0 )\n"
                                     "(  7  0  1  0 )\n"
                                     "( -2  0  1  1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the assignment to rows of the UniLowerMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the assignment to rows of the UniLowerMatrix specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testRow()
{
   //=====================================================================================
   // Row-major general tests
   //=====================================================================================

   {
      test_ = "Row-major row() function";

      using RT = blaze::Row<LT>;

      LT lower( 3UL );
      lower(1,0) = -4;
      lower(2,0) =  7;

      RT row1 = row( lower, 1UL );

      if( row1[0] != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator access failed\n"
             << " Details:\n"
             << "   Result: " << row1[0] << "\n"
             << "   Expected result: -4\n";
         throw std::runtime_error( oss.str() );
      }

      RT::Iterator it( row1.begin() );

      if( it == row1.end() || *it != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << *it << "\n"
             << "   Expected result: -4\n";
         throw std::runtime_error( oss.str() );
      }

      row1[0] = -5;

      if( row1[0] != -5 || row1[1] != 1 || row1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row access failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( -5 1 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -5 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 0 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row access failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -5 1 0 )\n(  7 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      reset( row1 );

      if( row1[0] != 0 || row1[1] != 1 || row1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row reset failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( 0 1 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 0 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) != 7 || lower(2,1) != 0 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row reset failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 1 0 )\n( 7 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major scalar assignment
   //=====================================================================================

   {
      test_ = "Row-major row() function (scalar assignment test)";

      using RT = blaze::Row<LT>;

      LT lower( 3UL );
      lower(1,0) = -4;
      lower(2,0) =  7;

      RT row1 = row( lower, 1UL );
      row1 = 8;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 2UL );

      if( row1[0] != 8 || row1[1] != 1 || row1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row access failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( 8 1 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 8 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) != 7 || lower(2,1) != 0 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row access failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 8 1 0 )\n( 7 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major general tests
   //=====================================================================================

   {
      test_ = "Column-major row() function";

      using RT = blaze::Row<OLT>;

      OLT lower( 3UL );
      lower(1,0) = -4;
      lower(2,0) =  7;

      RT row1 = row( lower, 1UL );

      if( row1[0] != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator access failed\n"
             << " Details:\n"
             << "   Result: " << row1[0] << "\n"
             << "   Expected result: -4\n";
         throw std::runtime_error( oss.str() );
      }

      RT::Iterator it( row1.begin() );

      if( it == row1.end() || *it != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << *it << "\n"
             << "   Expected result: -4\n";
         throw std::runtime_error( oss.str() );
      }

      row1[0] = -5;

      if( row1[0] != -5 || row1[1] != 1 || row1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row access failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( -5 1 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -5 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 0 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row access failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -5 1 0 )\n(  7 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      reset( row1 );

      if( row1[0] != 0 || row1[1] != 1 || row1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row reset failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( 0 1 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 0 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) != 7 || lower(2,1) != 0 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row reset failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 1 0 )\n( 7 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major scalar assignment
   //=====================================================================================

   {
      test_ = "Column-major row() function (scalar assignment test)";

      using RT = blaze::Row<OLT>;

      OLT lower( 3UL );
      lower(1,0) = -4;
      lower(2,0) =  7;

      RT row1 = row( lower, 1UL );
      row1 = 8;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 3UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( row1[0] != 8 || row1[1] != 1 || row1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row access failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( 8 1 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 8 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) != 7 || lower(2,1) != 0 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row access failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 8 1 0 )\n( 7 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the assignment to columns of the UniLowerMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the assignment to columns of the UniLowerMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testColumn()
{
   //=====================================================================================
   // Row-major general tests
   //=====================================================================================

   {
      test_ = "Row-major column() function";

      using CT = blaze::Column<LT>;

      LT lower( 3UL );
      lower(1,0) = -4;
      lower(2,0) =  7;

      CT col0 = column( lower, 0UL );

      if( col0[0] != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator access failed\n"
             << " Details:\n"
             << "   Result: " << col0[0] << "\n"
             << "   Expected result: 1\n";
         throw std::runtime_error( oss.str() );
      }

      CT::Iterator it( col0.begin() );

      if( it == col0.end() || *it != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << *it << "\n"
             << "   Expected result: 1\n";
         throw std::runtime_error( oss.str() );
      }

      col0[1] = -5;

      if( col0[0] != 1 || col0[1] != -5 || col0[2] != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column access failed\n"
             << " Details:\n"
             << "   Result:\n" << col0 << "\n"
             << "   Expected result:\n( 1 -5  7 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -5 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 0 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column access failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -5 1 0 )\n(  7 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      reset( col0 );

      if( col0[0] != 1 || col0[1] != 0 || col0[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column reset failed\n"
             << " Details:\n"
             << "   Result:\n" << col0 << "\n"
             << "   Expected result:\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 0 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) != 0 || lower(2,1) != 0 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column reset failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 1 0 )\n( 0 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major scalar assignment
   //=====================================================================================

   {
      test_ = "Row-major column() function (scalar assignment test)";

      using CT = blaze::Column<LT>;

      LT lower( 3UL );
      lower(1,0) = -4;
      lower(2,0) =  7;

      CT col1 = column( lower, 1UL );
      col1 = 8;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkNonZeros( lower, 6UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 3UL );

      if( col1[0] != 0 || col1[1] != 1 || col1[2] != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column access failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( 0 1 8 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -4 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 8 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column access failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 1 0 )\n(  7 8 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major general tests
   //=====================================================================================

   {
      test_ = "Column-major column() function";

      using CT = blaze::Column<OLT>;

      OLT lower( 3UL );
      lower(1,0) = -4;
      lower(2,0) =  7;

      CT col0 = column( lower, 0UL );

      if( col0[0] != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator access failed\n"
             << " Details:\n"
             << "   Result: " << col0[0] << "\n"
             << "   Expected result: 1\n";
         throw std::runtime_error( oss.str() );
      }

      CT::Iterator it( col0.begin() );

      if( it == col0.end() || *it != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << *it << "\n"
             << "   Expected result: 1\n";
         throw std::runtime_error( oss.str() );
      }

      col0[1] = -5;

      if( col0[0] != 1 || col0[1] != -5 || col0[2] != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column access failed\n"
             << " Details:\n"
             << "   Result:\n" << col0 << "\n"
             << "   Expected result:\n( 1 -5  7 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -5 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 0 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column access failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -5 1 0 )\n(  7 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      reset( col0 );

      if( col0[0] != 1 || col0[1] != 0 || col0[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column reset failed\n"
             << " Details:\n"
             << "   Result:\n" << col0 << "\n"
             << "   Expected result:\n( 1 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 0 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) != 0 || lower(2,1) != 0 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column reset failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 1 0 )\n( 0 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major scalar assignment
   //=====================================================================================

   {
      test_ = "Column-major column() function (scalar assignment test)";

      using CT = blaze::Column<OLT>;

      OLT lower( 3UL );
      lower(1,0) = -4;
      lower(2,0) =  7;

      CT col1 = column( lower, 1UL );
      col1 = 8;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkNonZeros( lower, 6UL );
      checkNonZeros( lower, 0UL, 3UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( col1[0] != 0 || col1[1] != 1 || col1[2] != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column access failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( 0 1 8 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -4 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 8 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column access failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 1 0 )\n(  7 8 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************

} // namespace unilowermatrix

} // namespace adaptors

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
   std::cout << "   Running UniLowerMatrix dense test (part 2)..." << std::endl;

   try
   {
      RUN_UNILOWERMATRIX_DENSE_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during UniLowerMatrix dense test (part 2):\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
