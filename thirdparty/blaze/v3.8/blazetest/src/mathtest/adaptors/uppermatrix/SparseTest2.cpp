//=================================================================================================
/*!
//  \file src/mathtest/adaptors/uppermatrix/SparseTest2.cpp
//  \brief Source file for the UpperMatrix sparse test (part 2)
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
#include <blaze/math/Column.h>
#include <blaze/math/Row.h>
#include <blaze/math/Submatrix.h>
#include <blaze/util/Complex.h>
#include <blazetest/mathtest/adaptors/uppermatrix/SparseTest.h>

#ifdef BLAZE_USE_HPX_THREADS
#  include <hpx/hpx_main.hpp>
#endif


namespace blazetest {

namespace mathtest {

namespace adaptors {

namespace uppermatrix {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the UpperMatrix sparse test.
//
// \exception std::runtime_error Operation error detected.
*/
SparseTest::SparseTest()
{
   testScaling();
   testFunctionCall();
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
/*!\brief Test of all UpperMatrix (self-)scaling operations.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all available ways to scale an instance of the UpperMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testScaling()
{
   //=====================================================================================
   // Row-major self-scaling (M*=s)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (M*=s)";

      UT upper( 3UL );
      upper(0,1) =  1;
      upper(0,2) = -2;
      upper(1,2) =  3;
      upper(2,2) = -4;

      upper *= 2;

      checkRows    ( upper, 3UL );
      checkColumns ( upper, 3UL );
      checkCapacity( upper, 4UL );
      checkNonZeros( upper, 4UL );
      checkNonZeros( upper, 0UL, 2UL );
      checkNonZeros( upper, 1UL, 1UL );
      checkNonZeros( upper, 2UL, 1UL );

      if( upper(0,0) != 0 || upper(0,1) != 2 || upper(0,2) != -4 ||
          upper(1,0) != 0 || upper(1,1) != 0 || upper(1,2) !=  6 ||
          upper(2,0) != 0 || upper(2,1) != 0 || upper(2,2) != -8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 0  2 -4 )\n( 0  0  6 )\n( 0  0 -8 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major self-scaling (M=M*s)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (M=M*s)";

      UT upper( 3UL );
      upper(0,1) =  1;
      upper(0,2) = -2;
      upper(1,2) =  3;
      upper(2,2) = -4;

      upper = upper * 2;

      checkRows    ( upper, 3UL );
      checkColumns ( upper, 3UL );
      checkCapacity( upper, 4UL );
      checkNonZeros( upper, 4UL );
      checkNonZeros( upper, 0UL, 2UL );
      checkNonZeros( upper, 1UL, 1UL );
      checkNonZeros( upper, 2UL, 1UL );

      if( upper(0,0) != 0 || upper(0,1) != 2 || upper(0,2) != -4 ||
          upper(1,0) != 0 || upper(1,1) != 0 || upper(1,2) !=  6 ||
          upper(2,0) != 0 || upper(2,1) != 0 || upper(2,2) != -8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 0  2 -4 )\n( 0  0  6 )\n( 0  0 -8 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major self-scaling (M=s*M)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (M=s*M)";

      UT upper( 3UL );
      upper(0,1) =  1;
      upper(0,2) = -2;
      upper(1,2) =  3;
      upper(2,2) = -4;

      upper = 2 * upper;

      checkRows    ( upper, 3UL );
      checkColumns ( upper, 3UL );
      checkCapacity( upper, 4UL );
      checkNonZeros( upper, 4UL );
      checkNonZeros( upper, 0UL, 2UL );
      checkNonZeros( upper, 1UL, 1UL );
      checkNonZeros( upper, 2UL, 1UL );

      if( upper(0,0) != 0 || upper(0,1) != 2 || upper(0,2) != -4 ||
          upper(1,0) != 0 || upper(1,1) != 0 || upper(1,2) !=  6 ||
          upper(2,0) != 0 || upper(2,1) != 0 || upper(2,2) != -8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 0  2 -4 )\n( 0  0  6 )\n( 0  0 -8 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major self-scaling (M/=s)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (M/=s)";

      UT upper( 3UL );
      upper(0,1) =  2;
      upper(0,2) = -4;
      upper(1,2) =  6;
      upper(2,2) = -8;

      upper /= 2;

      checkRows    ( upper, 3UL );
      checkColumns ( upper, 3UL );
      checkCapacity( upper, 4UL );
      checkNonZeros( upper, 4UL );
      checkNonZeros( upper, 0UL, 2UL );
      checkNonZeros( upper, 1UL, 1UL );
      checkNonZeros( upper, 2UL, 1UL );

      if( upper(0,0) != 0 || upper(0,1) != 1 || upper(0,2) != -2 ||
          upper(1,0) != 0 || upper(1,1) != 0 || upper(1,2) !=  3 ||
          upper(2,0) != 0 || upper(2,1) != 0 || upper(2,2) != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 0  1 -2 )\n( 0  0  3 )\n( 0  0 -4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major self-scaling (M=M/s)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (M=M/s)";

      UT upper( 3UL );
      upper(0,1) =  2;
      upper(0,2) = -4;
      upper(1,2) =  6;
      upper(2,2) = -8;

      upper = upper / 2;

      checkRows    ( upper, 3UL );
      checkColumns ( upper, 3UL );
      checkCapacity( upper, 4UL );
      checkNonZeros( upper, 4UL );
      checkNonZeros( upper, 0UL, 2UL );
      checkNonZeros( upper, 1UL, 1UL );
      checkNonZeros( upper, 2UL, 1UL );

      if( upper(0,0) != 0 || upper(0,1) != 1 || upper(0,2) != -2 ||
          upper(1,0) != 0 || upper(1,1) != 0 || upper(1,2) !=  3 ||
          upper(2,0) != 0 || upper(2,1) != 0 || upper(2,2) != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 0  1 -2 )\n( 0  0  3 )\n( 0  0 -4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major UpperMatrix::scale()
   //=====================================================================================

   {
      test_ = "Row-major UpperMatrix::scale()";

      // Initialization check
      UT upper( 3UL );
      upper(0,1) =  1;
      upper(0,2) = -2;
      upper(1,2) =  3;
      upper(2,2) = -4;

      checkRows    ( upper, 3UL );
      checkColumns ( upper, 3UL );
      checkCapacity( upper, 4UL );
      checkNonZeros( upper, 4UL );
      checkNonZeros( upper, 0UL, 2UL );
      checkNonZeros( upper, 1UL, 1UL );
      checkNonZeros( upper, 2UL, 1UL );

      if( upper(0,0) != 0 || upper(0,1) != 1 || upper(0,2) != -2 ||
          upper(1,0) != 0 || upper(1,1) != 0 || upper(1,2) !=  3 ||
          upper(2,0) != 0 || upper(2,1) != 0 || upper(2,2) != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 0  1 -2 )\n( 0  0  3 )\n( 0  0 -4 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Integral scaling of the matrix
      upper.scale( 2 );

      checkRows    ( upper, 3UL );
      checkColumns ( upper, 3UL );
      checkCapacity( upper, 4UL );
      checkNonZeros( upper, 4UL );
      checkNonZeros( upper, 0UL, 2UL );
      checkNonZeros( upper, 1UL, 1UL );
      checkNonZeros( upper, 2UL, 1UL );

      if( upper(0,0) != 0 || upper(0,1) != 2 || upper(0,2) != -4 ||
          upper(1,0) != 0 || upper(1,1) != 0 || upper(1,2) !=  6 ||
          upper(2,0) != 0 || upper(2,1) != 0 || upper(2,2) != -8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 0  2 -4 )\n( 0  0  6 )\n( 0  0 -8 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Floating point scaling of the matrix
      upper.scale( 0.5 );

      checkRows    ( upper, 3UL );
      checkColumns ( upper, 3UL );
      checkCapacity( upper, 4UL );
      checkNonZeros( upper, 4UL );
      checkNonZeros( upper, 0UL, 2UL );
      checkNonZeros( upper, 1UL, 1UL );
      checkNonZeros( upper, 2UL, 1UL );

      if( upper(0,0) != 0 || upper(0,1) != 1 || upper(0,2) != -2 ||
          upper(1,0) != 0 || upper(1,1) != 0 || upper(1,2) !=  3 ||
          upper(2,0) != 0 || upper(2,1) != 0 || upper(2,2) != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 0  1 -2 )\n( 0  0  3 )\n( 0  0 -4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major UpperMatrix::scale() (complex)";

      using blaze::complex;

      blaze::UpperMatrix< blaze::CompressedMatrix<complex<float>,blaze::rowMajor> > upper( 2UL );
      upper(0,0) = complex<float>( 1.0F, 0.0F );
      upper(0,1) = complex<float>( 2.0F, 0.0F );
      upper(1,1) = complex<float>( 4.0F, 0.0F );

      upper.scale( complex<float>( 3.0F, 0.0F ) );

      checkRows    ( upper, 2UL );
      checkColumns ( upper, 2UL );
      checkCapacity( upper, 3UL );
      checkNonZeros( upper, 3UL );
      checkNonZeros( upper, 0UL, 2UL );
      checkNonZeros( upper, 1UL, 1UL );

      if( upper(0,0) != complex<float>( 3.0F, 0.0F ) || upper(0,1) != complex<float>(  6.0F, 0.0F ) ||
          upper(1,0) != complex<float>( 0.0F, 0.0F ) || upper(1,1) != complex<float>( 12.0F, 0.0F ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( ( 3,0) ( 6,0)\n( 0,0) (12,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (M*=s)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (M*=s)";

      OUT upper( 3UL );
      upper(0,1) =  1;
      upper(0,2) = -2;
      upper(1,2) =  3;
      upper(2,2) = -4;

      upper *= 2;

      checkRows    ( upper, 3UL );
      checkColumns ( upper, 3UL );
      checkCapacity( upper, 4UL );
      checkNonZeros( upper, 4UL );
      checkNonZeros( upper, 0UL, 0UL );
      checkNonZeros( upper, 1UL, 1UL );
      checkNonZeros( upper, 2UL, 3UL );

      if( upper(0,0) != 0 || upper(0,1) != 2 || upper(0,2) != -4 ||
          upper(1,0) != 0 || upper(1,1) != 0 || upper(1,2) !=  6 ||
          upper(2,0) != 0 || upper(2,1) != 0 || upper(2,2) != -8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 0  2 -4 )\n( 0  0  6 )\n( 0  0 -8 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (M=M*s)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (M=M*s)";

      OUT upper( 3UL );
      upper(0,1) =  1;
      upper(0,2) = -2;
      upper(1,2) =  3;
      upper(2,2) = -4;

      upper = upper * 2;

      checkRows    ( upper, 3UL );
      checkColumns ( upper, 3UL );
      checkCapacity( upper, 4UL );
      checkNonZeros( upper, 4UL );
      checkNonZeros( upper, 0UL, 0UL );
      checkNonZeros( upper, 1UL, 1UL );
      checkNonZeros( upper, 2UL, 3UL );

      if( upper(0,0) != 0 || upper(0,1) != 2 || upper(0,2) != -4 ||
          upper(1,0) != 0 || upper(1,1) != 0 || upper(1,2) !=  6 ||
          upper(2,0) != 0 || upper(2,1) != 0 || upper(2,2) != -8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 0  2 -4 )\n( 0  0  6 )\n( 0  0 -8 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (M=s*M)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (M=s*M)";

      OUT upper( 3UL );
      upper(0,1) =  1;
      upper(0,2) = -2;
      upper(1,2) =  3;
      upper(2,2) = -4;

      upper = 2 * upper;

      checkRows    ( upper, 3UL );
      checkColumns ( upper, 3UL );
      checkCapacity( upper, 4UL );
      checkNonZeros( upper, 4UL );
      checkNonZeros( upper, 0UL, 0UL );
      checkNonZeros( upper, 1UL, 1UL );
      checkNonZeros( upper, 2UL, 3UL );

      if( upper(0,0) != 0 || upper(0,1) != 2 || upper(0,2) != -4 ||
          upper(1,0) != 0 || upper(1,1) != 0 || upper(1,2) !=  6 ||
          upper(2,0) != 0 || upper(2,1) != 0 || upper(2,2) != -8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 0  2 -4 )\n( 0  0  6 )\n( 0  0 -8 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (M/=s)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (M/=s)";

      OUT upper( 3UL );
      upper(0,1) =  2;
      upper(0,2) = -4;
      upper(1,2) =  6;
      upper(2,2) = -8;

      upper /= 2;

      checkRows    ( upper, 3UL );
      checkColumns ( upper, 3UL );
      checkCapacity( upper, 4UL );
      checkNonZeros( upper, 4UL );
      checkNonZeros( upper, 0UL, 0UL );
      checkNonZeros( upper, 1UL, 1UL );
      checkNonZeros( upper, 2UL, 3UL );

      if( upper(0,0) != 0 || upper(0,1) != 1 || upper(0,2) != -2 ||
          upper(1,0) != 0 || upper(1,1) != 0 || upper(1,2) !=  3 ||
          upper(2,0) != 0 || upper(2,1) != 0 || upper(2,2) != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 0  1 -2 )\n( 0  0  3 )\n( 0  0 -4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (M=M/s)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (M=M/s)";

      OUT upper( 3UL );
      upper(0,1) =  2;
      upper(0,2) = -4;
      upper(1,2) =  6;
      upper(2,2) = -8;

      upper = upper / 2;

      checkRows    ( upper, 3UL );
      checkColumns ( upper, 3UL );
      checkCapacity( upper, 4UL );
      checkNonZeros( upper, 4UL );
      checkNonZeros( upper, 0UL, 0UL );
      checkNonZeros( upper, 1UL, 1UL );
      checkNonZeros( upper, 2UL, 3UL );

      if( upper(0,0) != 0 || upper(0,1) != 1 || upper(0,2) != -2 ||
          upper(1,0) != 0 || upper(1,1) != 0 || upper(1,2) !=  3 ||
          upper(2,0) != 0 || upper(2,1) != 0 || upper(2,2) != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 0  1 -2 )\n( 0  0  3 )\n( 0  0 -4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major UpperMatrix::scale()
   //=====================================================================================

   {
      test_ = "Column-major UpperMatrix::scale()";

      // Initialization check
      OUT upper( 3UL );
      upper(0,1) =  1;
      upper(0,2) = -2;
      upper(1,2) =  3;
      upper(2,2) = -4;

      checkRows    ( upper, 3UL );
      checkColumns ( upper, 3UL );
      checkCapacity( upper, 4UL );
      checkNonZeros( upper, 4UL );
      checkNonZeros( upper, 0UL, 0UL );
      checkNonZeros( upper, 1UL, 1UL );
      checkNonZeros( upper, 2UL, 3UL );

      if( upper(0,0) != 0 || upper(0,1) != 1 || upper(0,2) != -2 ||
          upper(1,0) != 0 || upper(1,1) != 0 || upper(1,2) !=  3 ||
          upper(2,0) != 0 || upper(2,1) != 0 || upper(2,2) != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 0  1 -2 )\n( 0  0  3 )\n( 0  0 -4 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Integral scaling of the matrix
      upper.scale( 2 );

      checkRows    ( upper, 3UL );
      checkColumns ( upper, 3UL );
      checkCapacity( upper, 4UL );
      checkNonZeros( upper, 4UL );
      checkNonZeros( upper, 0UL, 0UL );
      checkNonZeros( upper, 1UL, 1UL );
      checkNonZeros( upper, 2UL, 3UL );

      if( upper(0,0) != 0 || upper(0,1) != 2 || upper(0,2) != -4 ||
          upper(1,0) != 0 || upper(1,1) != 0 || upper(1,2) !=  6 ||
          upper(2,0) != 0 || upper(2,1) != 0 || upper(2,2) != -8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 0  2 -4 )\n( 0  0  6 )\n( 0  0 -8 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Floating point scaling of the matrix
      upper.scale( 0.5 );

      checkRows    ( upper, 3UL );
      checkColumns ( upper, 3UL );
      checkCapacity( upper, 4UL );
      checkNonZeros( upper, 4UL );
      checkNonZeros( upper, 0UL, 0UL );
      checkNonZeros( upper, 1UL, 1UL );
      checkNonZeros( upper, 2UL, 3UL );

      if( upper(0,0) != 0 || upper(0,1) != 1 || upper(0,2) != -2 ||
          upper(1,0) != 0 || upper(1,1) != 0 || upper(1,2) !=  3 ||
          upper(2,0) != 0 || upper(2,1) != 0 || upper(2,2) != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 0  1 -2 )\n( 0  0  3 )\n( 0  0 -4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major UpperMatrix::scale() (complex)";

      using blaze::complex;

      blaze::UpperMatrix< blaze::CompressedMatrix<complex<float>,blaze::columnMajor> > upper( 2UL );
      upper(0,0) = complex<float>( 1.0F, 0.0F );
      upper(0,1) = complex<float>( 2.0F, 0.0F );
      upper(1,1) = complex<float>( 4.0F, 0.0F );

      upper.scale( complex<float>( 3.0F, 0.0F ) );

      checkRows    ( upper, 2UL );
      checkColumns ( upper, 2UL );
      checkCapacity( upper, 3UL );
      checkNonZeros( upper, 3UL );
      checkNonZeros( upper, 0UL, 1UL );
      checkNonZeros( upper, 1UL, 2UL );

      if( upper(0,0) != complex<float>( 3.0F, 0.0F ) || upper(0,1) != complex<float>(  6.0F, 0.0F ) ||
          upper(1,0) != complex<float>( 0.0F, 0.0F ) || upper(1,1) != complex<float>( 12.0F, 0.0F ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( ( 3,0) ( 6,0)\n( 0,0) (12,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the UpperMatrix function call operator.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of adding and accessing elements via the function call operator
// of the UpperMatrix specialization. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
void SparseTest::testFunctionCall()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major UpperMatrix::operator()";

      // Good cases
      {
         UT upper( 3UL );

         // Writing the diagonal element (1,1)
         upper(1,1) = 1;

         checkRows    ( upper, 3UL );
         checkColumns ( upper, 3UL );
         checkCapacity( upper, 1UL );
         checkNonZeros( upper, 1UL );
         checkNonZeros( upper, 0UL, 0UL );
         checkNonZeros( upper, 1UL, 1UL );
         checkNonZeros( upper, 2UL, 0UL );

         if( upper(0,0) != 0 || upper(0,1) != 0 || upper(0,2) != 0 ||
             upper(1,0) != 0 || upper(1,1) != 1 || upper(1,2) != 0 ||
             upper(2,0) != 0 || upper(2,1) != 0 || upper(2,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 0 0 0 )\n( 0 1 0 )\n( 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Writing the upper element (1,2)
         upper(1,2) = 2;

         checkRows    ( upper, 3UL );
         checkColumns ( upper, 3UL );
         checkCapacity( upper, 2UL );
         checkNonZeros( upper, 2UL );
         checkNonZeros( upper, 0UL, 0UL );
         checkNonZeros( upper, 1UL, 2UL );
         checkNonZeros( upper, 2UL, 0UL );

         if( upper(0,0) != 0 || upper(0,1) != 0 || upper(0,2) != 0 ||
             upper(1,0) != 0 || upper(1,1) != 1 || upper(1,2) != 2 ||
             upper(2,0) != 0 || upper(2,1) != 0 || upper(2,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 0 0 0 )\n( 0 1 2 )\n( 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Writing the element (0,1)
         upper(0,1) = upper(1,2);

         checkRows    ( upper, 3UL );
         checkColumns ( upper, 3UL );
         checkCapacity( upper, 3UL );
         checkNonZeros( upper, 3UL );
         checkNonZeros( upper, 0UL, 1UL );
         checkNonZeros( upper, 1UL, 2UL );
         checkNonZeros( upper, 2UL, 0UL );

         if( upper(0,0) != 0 || upper(0,1) != 2 || upper(0,2) != 0 ||
             upper(1,0) != 0 || upper(1,1) != 1 || upper(1,2) != 2 ||
             upper(2,0) != 0 || upper(2,1) != 0 || upper(2,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 0 2 0 )\n( 0 1 2 )\n( 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Adding to the upper element (0,2)
         upper(0,2) += 3;

         checkRows    ( upper, 3UL );
         checkColumns ( upper, 3UL );
         checkCapacity( upper, 4UL );
         checkNonZeros( upper, 4UL );
         checkNonZeros( upper, 0UL, 2UL );
         checkNonZeros( upper, 1UL, 2UL );
         checkNonZeros( upper, 2UL, 0UL );

         if( upper(0,0) != 0 || upper(0,1) != 2 || upper(0,2) != 3 ||
             upper(1,0) != 0 || upper(1,1) != 1 || upper(1,2) != 2 ||
             upper(2,0) != 0 || upper(2,1) != 0 || upper(2,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 0 0 0 )\n( 2 1 0 )\n( 3 2 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Subtracting from the upper element (0,1)
         upper(0,1) -= 4;

         checkRows    ( upper, 3UL );
         checkColumns ( upper, 3UL );
         checkCapacity( upper, 4UL );
         checkNonZeros( upper, 4UL );
         checkNonZeros( upper, 0UL, 2UL );
         checkNonZeros( upper, 1UL, 2UL );
         checkNonZeros( upper, 2UL, 0UL );

         if( upper(0,0) != 0 || upper(0,1) != -2 || upper(0,2) != 3 ||
             upper(1,0) != 0 || upper(1,1) !=  1 || upper(1,2) != 2 ||
             upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 0 -2  3 )\n( 0  1  2 )\n( 0  0  0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Multiplying the upper element (1,2)
         upper(1,2) *= -3;

         checkRows    ( upper, 3UL );
         checkColumns ( upper, 3UL );
         checkCapacity( upper, 4UL );
         checkNonZeros( upper, 4UL );
         checkNonZeros( upper, 0UL, 2UL );
         checkNonZeros( upper, 1UL, 2UL );
         checkNonZeros( upper, 2UL, 0UL );

         if( upper(0,0) != 0 || upper(0,1) != -2 || upper(0,2) !=  3 ||
             upper(1,0) != 0 || upper(1,1) !=  1 || upper(1,2) != -6 ||
             upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) !=  0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 0 -2  3 )\n( 0  1 -6 )\n( 0  0  0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Dividing the upper element (1,2)
         upper(1,2) /= 2;

         checkRows    ( upper, 3UL );
         checkColumns ( upper, 3UL );
         checkCapacity( upper, 4UL );
         checkNonZeros( upper, 4UL );
         checkNonZeros( upper, 0UL, 2UL );
         checkNonZeros( upper, 1UL, 2UL );
         checkNonZeros( upper, 2UL, 0UL );

         if( upper(0,0) != 0 || upper(0,1) != -2 || upper(0,2) !=  3 ||
             upper(1,0) != 0 || upper(1,1) !=  1 || upper(1,2) != -3 ||
             upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) !=  0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 0 -2  3 )\n( 0  1 -3 )\n( 0  0  0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Failure cases
      {
         UT upper( 3UL );

         // Trying to write the lower element (2,1)
         try {
            upper(2,1) = 2;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to lower matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to write the lower element (1,0)
         try {
            upper(1,0) = upper(1,2);

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to lower matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to add to the lower element (2,0)
         try {
            upper(2,0) += 3;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Addition assignment to lower matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to subtract from the lower element (1,0)
         try {
            upper(1,0) -= 4;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subtraction assignment to lower matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to multiply the lower element (2,1)
         try {
            upper(2,1) *= -3;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Multiplication assignment to lower matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to divide the lower element (2,1)
         try {
            upper(2,1) /= 2;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Division assignment to lower matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major UpperMatrix::operator()";

      // Good cases
      {
         OUT upper( 3UL );

         // Writing the diagonal element (1,1)
         upper(1,1) = 1;

         checkRows    ( upper, 3UL );
         checkColumns ( upper, 3UL );
         checkCapacity( upper, 1UL );
         checkNonZeros( upper, 1UL );
         checkNonZeros( upper, 0UL, 0UL );
         checkNonZeros( upper, 1UL, 1UL );
         checkNonZeros( upper, 2UL, 0UL );

         if( upper(0,0) != 0 || upper(0,1) != 0 || upper(0,2) != 0 ||
             upper(1,0) != 0 || upper(1,1) != 1 || upper(1,2) != 0 ||
             upper(2,0) != 0 || upper(2,1) != 0 || upper(2,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 0 0 0 )\n( 0 1 0 )\n( 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Writing the upper element (1,2)
         upper(1,2) = 2;

         checkRows    ( upper, 3UL );
         checkColumns ( upper, 3UL );
         checkCapacity( upper, 2UL );
         checkNonZeros( upper, 2UL );
         checkNonZeros( upper, 0UL, 0UL );
         checkNonZeros( upper, 1UL, 1UL );
         checkNonZeros( upper, 2UL, 1UL );

         if( upper(0,0) != 0 || upper(0,1) != 0 || upper(0,2) != 0 ||
             upper(1,0) != 0 || upper(1,1) != 1 || upper(1,2) != 2 ||
             upper(2,0) != 0 || upper(2,1) != 0 || upper(2,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 0 0 0 )\n( 0 1 2 )\n( 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Writing the upper element (0,1)
         upper(0,1) = upper(1,2);

         checkRows    ( upper, 3UL );
         checkColumns ( upper, 3UL );
         checkCapacity( upper, 3UL );
         checkNonZeros( upper, 3UL );
         checkNonZeros( upper, 0UL, 0UL );
         checkNonZeros( upper, 1UL, 2UL );
         checkNonZeros( upper, 2UL, 1UL );

         if( upper(0,0) != 0 || upper(0,1) != 2 || upper(0,2) != 0 ||
             upper(1,0) != 0 || upper(1,1) != 1 || upper(1,2) != 2 ||
             upper(2,0) != 0 || upper(2,1) != 0 || upper(2,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 0 2 0 )\n( 0 1 2 )\n( 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Adding to the upper element (0,2)
         upper(0,2) += 3;

         checkRows    ( upper, 3UL );
         checkColumns ( upper, 3UL );
         checkCapacity( upper, 4UL );
         checkNonZeros( upper, 4UL );
         checkNonZeros( upper, 0UL, 0UL );
         checkNonZeros( upper, 1UL, 2UL );
         checkNonZeros( upper, 2UL, 2UL );

         if( upper(0,0) != 0 || upper(0,1) != 2 || upper(0,2) != 3 ||
             upper(1,0) != 0 || upper(1,1) != 1 || upper(1,2) != 2 ||
             upper(2,0) != 0 || upper(2,1) != 0 || upper(2,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 0 2 3 )\n( 0 1 2 )\n( 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Subtracting from the upper element (0,1)
         upper(0,1) -= 4;

         checkRows    ( upper, 3UL );
         checkColumns ( upper, 3UL );
         checkCapacity( upper, 4UL );
         checkNonZeros( upper, 4UL );
         checkNonZeros( upper, 0UL, 0UL );
         checkNonZeros( upper, 1UL, 2UL );
         checkNonZeros( upper, 2UL, 2UL );

         if( upper(0,0) != 0 || upper(0,1) != -2 || upper(0,2) != 3 ||
             upper(1,0) != 0 || upper(1,1) !=  1 || upper(1,2) != 2 ||
             upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 0 -2  3 )\n( 0  1  2 )\n( 0  0  0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Multiplying the upper element (1,2)
         upper(1,2) *= -3;

         checkRows    ( upper, 3UL );
         checkColumns ( upper, 3UL );
         checkCapacity( upper, 4UL );
         checkNonZeros( upper, 4UL );
         checkNonZeros( upper, 0UL, 0UL );
         checkNonZeros( upper, 1UL, 2UL );
         checkNonZeros( upper, 2UL, 2UL );

         if( upper(0,0) != 0 || upper(0,1) != -2 || upper(0,2) !=  3 ||
             upper(1,0) != 0 || upper(1,1) !=  1 || upper(1,2) != -6 ||
             upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) !=  0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 0 -2  3 )\n( 0  1 -6 )\n( 0  0  0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Dividing the upper element (1,2)
         upper(1,2) /= 2;

         checkRows    ( upper, 3UL );
         checkColumns ( upper, 3UL );
         checkCapacity( upper, 4UL );
         checkNonZeros( upper, 4UL );
         checkNonZeros( upper, 0UL, 0UL );
         checkNonZeros( upper, 1UL, 2UL );
         checkNonZeros( upper, 2UL, 2UL );

         if( upper(0,0) != 0 || upper(0,1) != -2 || upper(0,2) !=  3 ||
             upper(1,0) != 0 || upper(1,1) !=  1 || upper(1,2) != -3 ||
             upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) !=  0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 0 -2  3 )\n( 0  1 -3 )\n( 0  0  0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Failure cases
      {
         OUT upper( 3UL );

         // Trying to write the lower element (2,1)
         try {
            upper(2,1) = 2;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to lower matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to write the lower element (1,0)
         try {
            upper(1,0) = upper(1,2);

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to lower matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to add to the lower element (2,0)
         try {
            upper(2,0) += 3;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Addition assignment to lower matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to subtract from the lower element (1,0)
         try {
            upper(1,0) -= 4;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subtraction assignment to lower matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to multiply the lower element (2,1)
         try {
            upper(2,1) *= -3;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Multiplication assignment to lower matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to divide the lower element (2,1)
         try {
            upper(2,1) /= 2;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Division assignment to lower matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the UpperMatrix iterator implementation.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the iterator implementation of the UpperMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testIterator()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      using Iterator      = UT::Iterator;
      using ConstIterator = UT::ConstIterator;

      UT upper( 3UL );
      upper(0,0) =  1;
      upper(0,2) =  3;
      upper(1,1) = -2;
      upper(2,2) =  4;

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

         ConstIterator it( begin( upper, 1UL ) );

         if( it == end( upper, 1UL ) || it->value() != -2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Failed iterator conversion detected\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Counting the number of elements in 0th row via Iterator (end-begin)
      {
         test_ = "Row-major Iterator subtraction (end-begin)";

         const ptrdiff_t number( end( upper, 0UL ) - begin( upper, 0UL ) );

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

      // Counting the number of elements in 1st row via ConstIterator (end-begin)
      {
         test_ = "Row-major ConstIterator subtraction (end-begin)";

         const ptrdiff_t number( cend( upper, 1UL ) - cbegin( upper, 1UL ) );

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

      // Testing read-only access via ConstIterator
      {
         test_ = "Row-major read-only access via ConstIterator";

         ConstIterator it ( cbegin( upper, 0UL ) );
         ConstIterator end( cend( upper, 0UL ) );

         if( it == end || it->value() != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid initial iterator detected\n";
            throw std::runtime_error( oss.str() );
         }

         ++it;

         if( it == end || it->value() != 3 ) {
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

      // Testing assignment to upper elements via Iterator
      {
         test_ = "Row-major assignment to upper elements via Iterator";

         int value = 7;

         for( Iterator it=begin( upper, 0UL ); it!=end( upper, 0UL ); ++it ) {
            *it = value++;
         }

         if( upper(0,0) != 7 || upper(0,1) !=  0 || upper(0,2) != 8 ||
             upper(1,0) != 0 || upper(1,1) != -2 || upper(1,2) != 0 ||
             upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 7  0  8 )\n( 0 -2  0 )\n( 0  0  4 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing addition assignment to upper elements via Iterator
      {
         test_ = "Row-major addition assignment to upper elements via Iterator";

         int value = 4;

         for( Iterator it=begin( upper, 0UL ); it!=end( upper, 0UL ); ++it ) {
            *it += value++;
         }

         if( upper(0,0) != 11 || upper(0,1) !=  0 || upper(0,2) != 13 ||
             upper(1,0) !=  0 || upper(1,1) != -2 || upper(1,2) !=  0 ||
             upper(2,0) !=  0 || upper(2,1) !=  0 || upper(2,2) !=  4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 11  0 13 )\n(  0 -2  0 )\n(  0  0  4 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing subtraction assignment to upper elements via Iterator
      {
         test_ = "Row-major subtraction assignment to upper elements via Iterator";

         int value = 4;

         for( Iterator it=begin( upper, 0UL ); it!=end( upper, 0UL ); ++it ) {
            *it -= value++;
         }

         if( upper(0,0) != 7 || upper(0,1) !=  0 || upper(0,2) != 8 ||
             upper(1,0) != 0 || upper(1,1) != -2 || upper(1,2) != 0 ||
             upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 7  0  8 )\n( 0 -2  0 )\n( 0  0  4 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing multiplication assignment to upper elements via Iterator
      {
         test_ = "Row-major multiplication assignment to upper elements via Iterator";

         for( Iterator it=begin( upper, 0UL ); it!=end( upper, 0UL ); ++it ) {
            *it *= 2;
         }

         if( upper(0,0) != 14 || upper(0,1) !=  0 || upper(0,2) != 16 ||
             upper(1,0) !=  0 || upper(1,1) != -2 || upper(1,2) !=  0 ||
             upper(2,0) !=  0 || upper(2,1) !=  0 || upper(2,2) !=  4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 14  0 16 )\n(  0 -2  0 )\n(  0  0  4 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing division assignment to upper elements via Iterator
      {
         test_ = "Row-major division assignment to upper elements via Iterator";

         for( Iterator it=begin( upper, 0UL ); it!=end( upper, 0UL ); ++it ) {
            *it /= 2;
         }

         if( upper(0,0) != 7 || upper(0,1) !=  0 || upper(0,2) != 8 ||
             upper(1,0) != 0 || upper(1,1) != -2 || upper(1,2) != 0 ||
             upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 7  0  8 )\n( 0 -2  0 )\n( 0  0  4 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      using Iterator      = OUT::Iterator;
      using ConstIterator = OUT::ConstIterator;

      OUT upper( 3UL );
      upper(0,0) =  1;
      upper(0,2) =  3;
      upper(1,1) = -2;
      upper(2,2) =  4;

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

         ConstIterator it( begin( upper, 1UL ) );

         if( it == end( upper, 1UL ) || it->value() != -2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Failed iterator conversion detected\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Counting the number of elements in 0th column via Iterator (end-begin)
      {
         test_ = "Column-major Iterator subtraction (end-begin)";

         const ptrdiff_t number( end( upper, 0UL ) - begin( upper, 0UL ) );

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

      // Counting the number of elements in 1st column via ConstIterator (end-begin)
      {
         test_ = "Column-major ConstIterator subtraction (end-begin)";

         const ptrdiff_t number( cend( upper, 1UL ) - cbegin( upper, 1UL ) );

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

      // Testing read-only access via ConstIterator
      {
         test_ = "Column-major read-only access via ConstIterator";

         ConstIterator it ( cbegin( upper, 2UL ) );
         ConstIterator end( cend( upper, 2UL ) );

         if( it == end || it->value() != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid initial iterator detected\n";
            throw std::runtime_error( oss.str() );
         }

         ++it;

         if( it == end || it->value() != 4 ) {
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

      // Testing assignment to upper elements via Iterator
      {
         test_ = "Column-major assignment to upper elements via Iterator";

         int value = 7;

         for( Iterator it=begin( upper, 2UL ); it!=end( upper, 2UL ); ++it ) {
            *it = value++;
         }

         if( upper(0,0) != 1 || upper(0,1) !=  0 || upper(0,2) != 7 ||
             upper(1,0) != 0 || upper(1,1) != -2 || upper(1,2) != 0 ||
             upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 1  0  7 )\n( 0 -2  0 )\n( 0  0  8 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing addition assignment to upper elements via Iterator
      {
         test_ = "Column-major addition assignment to upper elements via Iterator";

         int value = 4;

         for( Iterator it=begin( upper, 2UL ); it!=end( upper, 2UL ); ++it ) {
            *it += value++;
         }

         if( upper(0,0) != 1 || upper(0,1) !=  0 || upper(0,2) != 11 ||
             upper(1,0) != 0 || upper(1,1) != -2 || upper(1,2) !=  0 ||
             upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) != 13 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 1  0 11 )\n( 0 -2  0 )\n( 0  0 13 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing subtraction assignment to upper elements via Iterator
      {
         test_ = "Column-major subtraction assignment to upper elements via Iterator";

         int value = 4;

         for( Iterator it=begin( upper, 2UL ); it!=end( upper, 2UL ); ++it ) {
            *it -= value++;
         }

         if( upper(0,0) != 1 || upper(0,1) !=  0 || upper(0,2) != 7 ||
             upper(1,0) != 0 || upper(1,1) != -2 || upper(1,2) != 0 ||
             upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 1  0  7 )\n( 0 -2  0 )\n( 0  0  8 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing multiplication assignment to upper elements via Iterator
      {
         test_ = "Column-major multiplication assignment to upper elements via Iterator";

         for( Iterator it=begin( upper, 2UL ); it!=end( upper, 2UL ); ++it ) {
            *it *= 2;
         }

         if( upper(0,0) != 1 || upper(0,1) !=  0 || upper(0,2) != 14 ||
             upper(1,0) != 0 || upper(1,1) != -2 || upper(1,2) !=  0 ||
             upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) != 16 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 1  0 14 )\n( 0 -2  0 )\n( 0  0 16 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing division assignment to upper elements via Iterator
      {
         test_ = "Column-major division assignment to upper elements via Iterator";

         for( Iterator it=begin( upper, 2UL ); it!=end( upper, 2UL ); ++it ) {
            *it /= 2;
         }

         if( upper(0,0) != 1 || upper(0,1) !=  0 || upper(0,2) != 7 ||
             upper(1,0) != 0 || upper(1,1) != -2 || upper(1,2) != 0 ||
             upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 1  0  7 )\n( 0 -2  0 )\n( 0  0  8 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c nonZeros() member function of the UpperMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c nonZeros() member function of the UpperMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testNonZeros()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major UpperMatrix::nonZeros()";

      // Empty matrix
      {
         UT upper( 3UL );

         checkRows    ( upper, 3UL );
         checkColumns ( upper, 3UL );
         checkNonZeros( upper, 0UL );
         checkNonZeros( upper, 0UL, 0UL );
         checkNonZeros( upper, 1UL, 0UL );
         checkNonZeros( upper, 2UL, 0UL );

         if( upper(0,0) != 0 || upper(0,1) != 0 || upper(0,2) != 0 ||
             upper(1,0) != 0 || upper(1,1) != 0 || upper(1,2) != 0 ||
             upper(2,0) != 0 || upper(2,1) != 0 || upper(2,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 0 0 0 )\n( 0 0 0 )\n( 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Partially filled matrix
      {
         UT upper( 3UL );
         upper(0,0) =  1;
         upper(1,1) = -2;
         upper(1,2) =  3;
         upper(2,2) = -4;

         checkRows    ( upper, 3UL );
         checkColumns ( upper, 3UL );
         checkCapacity( upper, 4UL );
         checkNonZeros( upper, 4UL );
         checkNonZeros( upper, 0UL, 1UL );
         checkNonZeros( upper, 1UL, 2UL );
         checkNonZeros( upper, 2UL, 1UL );

         if( upper(0,0) != 1 || upper(0,1) !=  0 || upper(0,2) !=  0 ||
             upper(1,0) != 0 || upper(1,1) != -2 || upper(1,2) !=  3 ||
             upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) != -4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 1  0  0 )\n( 0 -2  3 )\n( 0  0 -4 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Fully filled matrix
      {
         UT upper( 3UL );
         upper(0,0) = -1;
         upper(0,1) =  2;
         upper(0,2) =  3;
         upper(1,1) = -4;
         upper(1,2) = -5;
         upper(2,2) =  6;

         checkRows    ( upper, 3UL );
         checkColumns ( upper, 3UL );
         checkCapacity( upper, 6UL );
         checkNonZeros( upper, 6UL );
         checkNonZeros( upper, 0UL, 3UL );
         checkNonZeros( upper, 1UL, 2UL );
         checkNonZeros( upper, 2UL, 1UL );

         if( upper(0,0) != -1 || upper(0,1) !=  2 || upper(0,2) !=  3 ||
             upper(1,0) !=  0 || upper(1,1) != -4 || upper(1,2) != -5 ||
             upper(2,0) !=  0 || upper(2,1) !=  0 || upper(2,2) !=  6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( -1  2  3 )\n(  0 -4 -5 )\n(  0  0  6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major UpperMatrix::nonZeros()";

      // Empty matrix
      {
         OUT upper( 3UL );

         checkRows    ( upper, 3UL );
         checkColumns ( upper, 3UL );
         checkNonZeros( upper, 0UL );
         checkNonZeros( upper, 0UL, 0UL );
         checkNonZeros( upper, 1UL, 0UL );
         checkNonZeros( upper, 2UL, 0UL );

         if( upper(0,0) != 0 || upper(0,1) != 0 || upper(0,2) != 0 ||
             upper(1,0) != 0 || upper(1,1) != 0 || upper(1,2) != 0 ||
             upper(2,0) != 0 || upper(2,1) != 0 || upper(2,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 0 0 0 )\n( 0 0 0 )\n( 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Partially filled matrix
      {
         OUT upper( 3UL );
         upper(0,0) =  1;
         upper(1,1) = -2;
         upper(1,2) =  3;
         upper(2,2) = -4;

         checkRows    ( upper, 3UL );
         checkColumns ( upper, 3UL );
         checkCapacity( upper, 4UL );
         checkNonZeros( upper, 4UL );
         checkNonZeros( upper, 0UL, 1UL );
         checkNonZeros( upper, 1UL, 1UL );
         checkNonZeros( upper, 2UL, 2UL );

         if( upper(0,0) != 1 || upper(0,1) !=  0 || upper(0,2) !=  0 ||
             upper(1,0) != 0 || upper(1,1) != -2 || upper(1,2) !=  3 ||
             upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) != -4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 1  0  0 )\n( 0 -2  3 )\n( 0  0 -4 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Fully filled matrix
      {
         OUT upper( 3UL );
         upper(0,0) = -1;
         upper(0,1) =  2;
         upper(0,2) =  3;
         upper(1,1) = -4;
         upper(1,2) = -5;
         upper(2,2) =  6;

         checkRows    ( upper, 3UL );
         checkColumns ( upper, 3UL );
         checkCapacity( upper, 6UL );
         checkNonZeros( upper, 6UL );
         checkNonZeros( upper, 0UL, 1UL );
         checkNonZeros( upper, 1UL, 2UL );
         checkNonZeros( upper, 2UL, 3UL );

         if( upper(0,0) != -1 || upper(0,1) !=  2 || upper(0,2) !=  3 ||
             upper(1,0) !=  0 || upper(1,1) != -4 || upper(1,2) != -5 ||
             upper(2,0) !=  0 || upper(2,1) !=  0 || upper(2,2) !=  6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( -1  2  3 )\n(  0 -4 -5 )\n(  0  0  6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c reset() member function of the UpperMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c reset() member function of the UpperMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testReset()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major UpperMatrix::reset()";

      // Initialization check
      UT upper( 3UL );
      upper(0,0) = 1;
      upper(0,1) = 2;
      upper(0,2) = 3;
      upper(1,1) = 4;
      upper(1,2) = 5;
      upper(2,2) = 6;

      checkRows    ( upper, 3UL );
      checkColumns ( upper, 3UL );
      checkCapacity( upper, 6UL );
      checkNonZeros( upper, 6UL );
      checkNonZeros( upper, 0UL, 3UL );
      checkNonZeros( upper, 1UL, 2UL );
      checkNonZeros( upper, 2UL, 1UL );

      if( upper(0,0) != 1 || upper(0,1) != 2 || upper(0,2) != 3 ||
          upper(1,0) != 0 || upper(1,1) != 4 || upper(1,2) != 5 ||
          upper(2,0) != 0 || upper(2,1) != 0 || upper(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 0 4 5 )\n( 0 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting an upper element
      reset( upper(0,1) );

      checkRows    ( upper, 3UL );
      checkColumns ( upper, 3UL );
      checkCapacity( upper, 6UL );
      checkNonZeros( upper, 5UL );
      checkNonZeros( upper, 0UL, 2UL );
      checkNonZeros( upper, 1UL, 2UL );
      checkNonZeros( upper, 2UL, 1UL );

      if( upper(0,0) != 1 || upper(0,1) != 0 || upper(0,2) != 3 ||
          upper(1,0) != 0 || upper(1,1) != 4 || upper(1,2) != 5 ||
          upper(2,0) != 0 || upper(2,1) != 0 || upper(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 0 3 )\n( 0 4 5 )\n( 0 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting a lower element
      reset( upper(1,0) );

      checkRows    ( upper, 3UL );
      checkColumns ( upper, 3UL );
      checkCapacity( upper, 6UL );
      checkNonZeros( upper, 5UL );
      checkNonZeros( upper, 0UL, 2UL );
      checkNonZeros( upper, 1UL, 2UL );
      checkNonZeros( upper, 2UL, 1UL );

      if( upper(0,0) != 1 || upper(0,1) != 0 || upper(0,2) != 3 ||
          upper(1,0) != 0 || upper(1,1) != 4 || upper(1,2) != 5 ||
          upper(2,0) != 0 || upper(2,1) != 0 || upper(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 0 3 )\n( 0 4 5 )\n( 0 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting row 1
      reset( upper, 1UL );

      checkRows    ( upper, 3UL );
      checkColumns ( upper, 3UL );
      checkCapacity( upper, 6UL );
      checkNonZeros( upper, 3UL );
      checkNonZeros( upper, 0UL, 2UL );
      checkNonZeros( upper, 1UL, 0UL );
      checkNonZeros( upper, 2UL, 1UL );

      if( upper(0,0) != 1 || upper(0,1) != 0 || upper(0,2) != 3 ||
          upper(1,0) != 0 || upper(1,1) != 0 || upper(1,2) != 0 ||
          upper(2,0) != 0 || upper(2,1) != 0 || upper(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 0 3 )\n( 0 0 0 )\n( 0 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting the entire matrix
      reset( upper );

      checkRows    ( upper, 3UL );
      checkColumns ( upper, 3UL );
      checkCapacity( upper, 6UL );
      checkNonZeros( upper, 0UL );
      checkNonZeros( upper, 0UL, 0UL );
      checkNonZeros( upper, 1UL, 0UL );
      checkNonZeros( upper, 2UL, 0UL );

      if( upper(0,0) != 0 || upper(0,1) != 0 || upper(0,2) != 0 ||
          upper(1,0) != 0 || upper(1,1) != 0 || upper(1,2) != 0 ||
          upper(2,0) != 0 || upper(2,1) != 0 || upper(2,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 0 0 0 )\n( 0 0 0 )\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major UpperMatrix::reset()";

      // Initialization check
      OUT upper( 3UL );
      upper(0,0) = 1;
      upper(0,1) = 2;
      upper(0,2) = 3;
      upper(1,1) = 4;
      upper(1,2) = 5;
      upper(2,2) = 6;

      checkRows    ( upper, 3UL );
      checkColumns ( upper, 3UL );
      checkCapacity( upper, 6UL );
      checkNonZeros( upper, 6UL );
      checkNonZeros( upper, 0UL, 1UL );
      checkNonZeros( upper, 1UL, 2UL );
      checkNonZeros( upper, 2UL, 3UL );

      if( upper(0,0) != 1 || upper(0,1) != 2 || upper(0,2) != 3 ||
          upper(1,0) != 0 || upper(1,1) != 4 || upper(1,2) != 5 ||
          upper(2,0) != 0 || upper(2,1) != 0 || upper(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 0 4 5 )\n( 0 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting an upper element
      reset( upper(0,1) );

      checkRows    ( upper, 3UL );
      checkColumns ( upper, 3UL );
      checkCapacity( upper, 6UL );
      checkNonZeros( upper, 5UL );
      checkNonZeros( upper, 0UL, 1UL );
      checkNonZeros( upper, 1UL, 1UL );
      checkNonZeros( upper, 2UL, 3UL );

      if( upper(0,0) != 1 || upper(0,1) != 0 || upper(0,2) != 3 ||
          upper(1,0) != 0 || upper(1,1) != 4 || upper(1,2) != 5 ||
          upper(2,0) != 0 || upper(2,1) != 0 || upper(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 0 3 )\n( 0 4 5 )\n( 0 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting a lower element
      reset( upper(1,0) );

      checkRows    ( upper, 3UL );
      checkColumns ( upper, 3UL );
      checkCapacity( upper, 6UL );
      checkNonZeros( upper, 5UL );
      checkNonZeros( upper, 0UL, 1UL );
      checkNonZeros( upper, 1UL, 1UL );
      checkNonZeros( upper, 2UL, 3UL );

      if( upper(0,0) != 1 || upper(0,1) != 0 || upper(0,2) != 3 ||
          upper(1,0) != 0 || upper(1,1) != 4 || upper(1,2) != 5 ||
          upper(2,0) != 0 || upper(2,1) != 0 || upper(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 0 3 )\n( 0 4 5 )\n( 0 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting column 1
      reset( upper, 1UL );

      checkRows    ( upper, 3UL );
      checkColumns ( upper, 3UL );
      checkCapacity( upper, 6UL );
      checkNonZeros( upper, 4UL );
      checkNonZeros( upper, 0UL, 1UL );
      checkNonZeros( upper, 1UL, 0UL );
      checkNonZeros( upper, 2UL, 3UL );

      if( upper(0,0) != 1 || upper(0,1) != 0 || upper(0,2) != 3 ||
          upper(1,0) != 0 || upper(1,1) != 0 || upper(1,2) != 5 ||
          upper(2,0) != 0 || upper(2,1) != 0 || upper(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 0 3 )\n( 0 0 5 )\n( 0 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting the entire matrix
      reset( upper );

      checkRows    ( upper, 3UL );
      checkColumns ( upper, 3UL );
      checkCapacity( upper, 6UL );
      checkNonZeros( upper, 0UL );
      checkNonZeros( upper, 0UL, 0UL );
      checkNonZeros( upper, 1UL, 0UL );
      checkNonZeros( upper, 2UL, 0UL );

      if( upper(0,0) != 0 || upper(0,1) != 0 || upper(0,2) != 0 ||
          upper(1,0) != 0 || upper(1,1) != 0 || upper(1,2) != 0 ||
          upper(2,0) != 0 || upper(2,1) != 0 || upper(2,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 0 0 0 )\n( 0 0 0 )\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c clear() member function of the UpperMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c clear() member function of the UpperMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testClear()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major UpperMatrix::clear()";

      // Initialization check
      UT upper( 3UL );
      upper(0,0) = 1;
      upper(0,1) = 2;
      upper(0,2) = 3;
      upper(1,1) = 4;
      upper(1,2) = 5;
      upper(2,2) = 6;

      checkRows    ( upper, 3UL );
      checkColumns ( upper, 3UL );
      checkCapacity( upper, 6UL );
      checkNonZeros( upper, 6UL );
      checkNonZeros( upper, 0UL, 3UL );
      checkNonZeros( upper, 1UL, 2UL );
      checkNonZeros( upper, 2UL, 1UL );

      if( upper(0,0) != 1 || upper(0,1) != 2 || upper(0,2) != 3 ||
          upper(1,0) != 0 || upper(1,1) != 4 || upper(1,2) != 5 ||
          upper(2,0) != 0 || upper(2,1) != 0 || upper(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 0 4 5 )\n( 0 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing an upper element
      clear( upper(0,1) );

      checkRows    ( upper, 3UL );
      checkColumns ( upper, 3UL );
      checkCapacity( upper, 6UL );
      checkNonZeros( upper, 5UL );
      checkNonZeros( upper, 0UL, 2UL );
      checkNonZeros( upper, 1UL, 2UL );
      checkNonZeros( upper, 2UL, 1UL );

      if( upper(0,0) != 1 || upper(0,1) != 0 || upper(0,2) != 3 ||
          upper(1,0) != 0 || upper(1,1) != 4 || upper(1,2) != 5 ||
          upper(2,0) != 0 || upper(2,1) != 0 || upper(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 0 3 )\n( 0 4 5 )\n( 0 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing a lower element
      clear( upper(1,0) );

      checkRows    ( upper, 3UL );
      checkColumns ( upper, 3UL );
      checkCapacity( upper, 6UL );
      checkNonZeros( upper, 5UL );
      checkNonZeros( upper, 0UL, 2UL );
      checkNonZeros( upper, 1UL, 2UL );
      checkNonZeros( upper, 2UL, 1UL );

      if( upper(0,0) != 1 || upper(0,1) != 0 || upper(0,2) != 3 ||
          upper(1,0) != 0 || upper(1,1) != 4 || upper(1,2) != 5 ||
          upper(2,0) != 0 || upper(2,1) != 0 || upper(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 0 3 )\n( 0 4 5 )\n( 0 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing the matrix
      clear( upper );

      checkRows    ( upper, 0UL );
      checkColumns ( upper, 0UL );
      checkNonZeros( upper, 0UL );
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major UpperMatrix::clear()";

      // Initialization check
      OUT upper( 3UL );
      upper(0,0) = 1;
      upper(0,1) = 2;
      upper(0,2) = 3;
      upper(1,1) = 4;
      upper(1,2) = 5;
      upper(2,2) = 6;

      checkRows    ( upper, 3UL );
      checkColumns ( upper, 3UL );
      checkCapacity( upper, 6UL );
      checkNonZeros( upper, 6UL );
      checkNonZeros( upper, 0UL, 1UL );
      checkNonZeros( upper, 1UL, 2UL );
      checkNonZeros( upper, 2UL, 3UL );

      if( upper(0,0) != 1 || upper(0,1) != 2 || upper(0,2) != 3 ||
          upper(1,0) != 0 || upper(1,1) != 4 || upper(1,2) != 5 ||
          upper(2,0) != 0 || upper(2,1) != 0 || upper(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 0 4 5 )\n( 0 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing an upper element
      clear( upper(0,1) );

      checkRows    ( upper, 3UL );
      checkColumns ( upper, 3UL );
      checkCapacity( upper, 6UL );
      checkNonZeros( upper, 5UL );
      checkNonZeros( upper, 0UL, 1UL );
      checkNonZeros( upper, 1UL, 1UL );
      checkNonZeros( upper, 2UL, 3UL );

      if( upper(0,0) != 1 || upper(0,1) != 0 || upper(0,2) != 3 ||
          upper(1,0) != 0 || upper(1,1) != 4 || upper(1,2) != 5 ||
          upper(2,0) != 0 || upper(2,1) != 0 || upper(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 0 3 )\n( 0 4 5 )\n( 0 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing a lower element
      clear( upper(1,0) );

      checkRows    ( upper, 3UL );
      checkColumns ( upper, 3UL );
      checkCapacity( upper, 6UL );
      checkNonZeros( upper, 5UL );
      checkNonZeros( upper, 0UL, 1UL );
      checkNonZeros( upper, 1UL, 1UL );
      checkNonZeros( upper, 2UL, 3UL );

      if( upper(0,0) != 1 || upper(0,1) != 0 || upper(0,2) != 3 ||
          upper(1,0) != 0 || upper(1,1) != 4 || upper(1,2) != 5 ||
          upper(2,0) != 0 || upper(2,1) != 0 || upper(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 0 3 )\n( 0 4 5 )\n( 0 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing the matrix
      clear( upper );

      checkRows    ( upper, 0UL );
      checkColumns ( upper, 0UL );
      checkNonZeros( upper, 0UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c resize() member function of the UpperMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c resize() member function of the UpperMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testResize()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major UpperMatrix::resize()";

      // Initialization check
      UT upper;

      checkRows    ( upper, 0UL );
      checkColumns ( upper, 0UL );
      checkNonZeros( upper, 0UL );

      // Resizing to 2x2
      upper.resize( 2UL );

      checkRows    ( upper, 2UL );
      checkColumns ( upper, 2UL );
      checkNonZeros( upper, 0UL );
      checkNonZeros( upper, 0UL, 0UL );
      checkNonZeros( upper, 1UL, 0UL );

      // Resizing to 4x4 and preserving the elements
      upper(0,0) = 1;
      upper(0,1) = 2;
      upper(1,1) = 3;
      upper.resize( 4UL, true );

      checkRows    ( upper, 4UL );
      checkColumns ( upper, 4UL );
      checkCapacity( upper, 3UL );
      checkNonZeros( upper, 3UL );
      checkNonZeros( upper, 0UL, 2UL );
      checkNonZeros( upper, 1UL, 1UL );
      checkNonZeros( upper, 2UL, 0UL );
      checkNonZeros( upper, 3UL, 0UL );

      // Resizing to 2x2
      upper(2,2) = 4;
      upper.resize( 2UL );

      checkRows    ( upper, 2UL );
      checkColumns ( upper, 2UL );
      checkCapacity( upper, 3UL );
      checkNonZeros( upper, 3UL );
      checkNonZeros( upper, 0UL, 2UL );
      checkNonZeros( upper, 1UL, 1UL );

      // Resizing to 0x0
      upper.resize( 0UL );

      checkRows    ( upper, 0UL );
      checkColumns ( upper, 0UL );
      checkNonZeros( upper, 0UL );
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major UpperMatrix::resize()";

      // Initialization check
      OUT upper;

      checkRows    ( upper, 0UL );
      checkColumns ( upper, 0UL );
      checkNonZeros( upper, 0UL );

      // Resizing to 2x2
      upper.resize( 2UL );

      checkRows    ( upper, 2UL );
      checkColumns ( upper, 2UL );
      checkNonZeros( upper, 0UL );
      checkNonZeros( upper, 0UL, 0UL );
      checkNonZeros( upper, 1UL, 0UL );

      // Resizing to 4x4 and preserving the elements
      upper(0,0) = 1;
      upper(0,1) = 2;
      upper(1,1) = 3;
      upper.resize( 4UL, true );

      checkRows    ( upper,  4UL );
      checkColumns ( upper,  4UL );
      checkCapacity( upper, 3UL );
      checkNonZeros( upper, 3UL );
      checkNonZeros( upper, 0UL, 1UL );
      checkNonZeros( upper, 1UL, 2UL );
      checkNonZeros( upper, 2UL, 0UL );
      checkNonZeros( upper, 3UL, 0UL );

      // Resizing to 2x2
      upper(2,2) = 4;
      upper.resize( 2UL );

      checkRows    ( upper, 2UL );
      checkColumns ( upper, 2UL );
      checkCapacity( upper, 4UL );
      checkNonZeros( upper, 3UL );
      checkNonZeros( upper, 0UL, 1UL );
      checkNonZeros( upper, 1UL, 2UL );

      // Resizing to 0x0
      upper.resize( 0UL );

      checkRows    ( upper, 0UL );
      checkColumns ( upper, 0UL );
      checkNonZeros( upper, 0UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c reserve() member function of the UpperMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c reserve() member function of the UpperMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testReserve()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major UpperMatrix::reserve()";

      // Initialization check
      UT upper;

      checkRows    ( upper, 0UL );
      checkColumns ( upper, 0UL );
      checkNonZeros( upper, 0UL );

      // Increasing the capacity of the matrix
      upper.reserve( 10UL );

      checkRows    ( upper,  0UL );
      checkColumns ( upper,  0UL );
      checkCapacity( upper, 10UL );
      checkNonZeros( upper,  0UL );

      // Further increasing the capacity of the matrix
      upper.reserve( 20UL );

      checkRows    ( upper,  0UL );
      checkColumns ( upper,  0UL );
      checkCapacity( upper, 20UL );
      checkNonZeros( upper,  0UL );
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major UpperMatrix::reserve()";

      // Initialization check
      OUT upper;

      checkRows    ( upper, 0UL );
      checkColumns ( upper, 0UL );
      checkNonZeros( upper, 0UL );

      // Increasing the capacity of the matrix
      upper.reserve( 10UL );

      checkRows    ( upper,  0UL );
      checkColumns ( upper,  0UL );
      checkCapacity( upper, 10UL );
      checkNonZeros( upper,  0UL );

      // Further increasing the capacity of the matrix
      upper.reserve( 20UL );

      checkRows    ( upper,  0UL );
      checkColumns ( upper,  0UL );
      checkCapacity( upper, 20UL );
      checkNonZeros( upper,  0UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c trim() member function of the UpperMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c trim() member function of the UpperMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testTrim()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major UpperMatrix::trim()";

      // Initialization check
      UT upper( 3UL );

      checkRows    ( upper, 3UL );
      checkColumns ( upper, 3UL );
      checkNonZeros( upper, 0UL );

      // Increasing the row capacity of the matrix
      upper.reserve( 0UL, 10UL );
      upper.reserve( 1UL, 15UL );
      upper.reserve( 2UL, 20UL );

      checkRows    ( upper,  3UL );
      checkColumns ( upper,  3UL );
      checkCapacity( upper, 45UL );
      checkCapacity( upper,  0UL, 10UL );
      checkCapacity( upper,  1UL, 15UL );
      checkCapacity( upper,  2UL, 20UL );

      // Trimming the matrix
      upper.trim();

      checkRows    ( upper,  3UL );
      checkColumns ( upper,  3UL );
      checkCapacity( upper, 45UL );
      checkCapacity( upper,  0UL, 0UL );
      checkCapacity( upper,  1UL, 0UL );
      checkCapacity( upper,  2UL, 0UL );
   }

   {
      test_ = "Row-major UpperMatrix::trim( size_t )";

      // Initialization check
      UT upper( 3UL, 3UL );

      checkRows    ( upper, 3UL );
      checkColumns ( upper, 3UL );
      checkNonZeros( upper, 0UL );

      // Increasing the row capacity of the matrix
      upper.reserve( 0UL, 10UL );
      upper.reserve( 1UL, 15UL );
      upper.reserve( 2UL, 20UL );

      checkRows    ( upper,  3UL );
      checkColumns ( upper,  3UL );
      checkCapacity( upper, 45UL );
      checkCapacity( upper,  0UL, 10UL );
      checkCapacity( upper,  1UL, 15UL );
      checkCapacity( upper,  2UL, 20UL );

      // Trimming the 0th row
      upper.trim( 0UL );

      checkRows    ( upper,  3UL );
      checkColumns ( upper,  3UL );
      checkCapacity( upper, 45UL );
      checkCapacity( upper,  0UL,  0UL );
      checkCapacity( upper,  1UL, 25UL );
      checkCapacity( upper,  2UL, 20UL );

      // Trimming the 1st row
      upper.trim( 1UL );

      checkRows    ( upper,  3UL );
      checkColumns ( upper,  3UL );
      checkCapacity( upper, 45UL );
      checkCapacity( upper,  0UL,  0UL );
      checkCapacity( upper,  1UL,  0UL );
      checkCapacity( upper,  2UL, 45UL );

      // Trimming the 2nd row
      upper.trim( 2UL );

      checkRows    ( upper,  3UL );
      checkColumns ( upper,  3UL );
      checkCapacity( upper, 45UL );
      checkCapacity( upper,  0UL, 0UL );
      checkCapacity( upper,  1UL, 0UL );
      checkCapacity( upper,  2UL, 0UL );
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major UpperMatrix::trim()";

      // Initialization check
      OUT upper( 3UL );

      checkRows    ( upper, 3UL );
      checkColumns ( upper, 3UL );
      checkNonZeros( upper, 0UL );

      // Increasing the row capacity of the matrix
      upper.reserve( 0UL, 10UL );
      upper.reserve( 1UL, 15UL );
      upper.reserve( 2UL, 20UL );

      checkRows    ( upper,  3UL );
      checkColumns ( upper,  3UL );
      checkCapacity( upper, 45UL );
      checkCapacity( upper,  0UL, 10UL );
      checkCapacity( upper,  1UL, 15UL );
      checkCapacity( upper,  2UL, 20UL );

      // Trimming the matrix
      upper.trim();

      checkRows    ( upper,  3UL );
      checkColumns ( upper,  3UL );
      checkCapacity( upper, 45UL );
      checkCapacity( upper,  0UL, 0UL );
      checkCapacity( upper,  1UL, 0UL );
      checkCapacity( upper,  2UL, 0UL );
   }

   {
      test_ = "Column-major UpperMatrix::trim( size_t )";

      // Initialization check
      OUT upper( 3UL, 3UL );

      checkRows    ( upper, 3UL );
      checkColumns ( upper, 3UL );
      checkNonZeros( upper, 0UL );

      // Increasing the column capacity of the matrix
      upper.reserve( 0UL, 10UL );
      upper.reserve( 1UL, 15UL );
      upper.reserve( 2UL, 20UL );

      checkRows    ( upper,  3UL );
      checkColumns ( upper,  3UL );
      checkCapacity( upper, 45UL );
      checkCapacity( upper,  0UL, 10UL );
      checkCapacity( upper,  1UL, 15UL );
      checkCapacity( upper,  2UL, 20UL );

      // Trimming the 0th column
      upper.trim( 0UL );

      checkRows    ( upper,  3UL );
      checkColumns ( upper,  3UL );
      checkCapacity( upper, 45UL );
      checkCapacity( upper,  0UL,  0UL );
      checkCapacity( upper,  1UL, 25UL );
      checkCapacity( upper,  2UL, 20UL );

      // Trimming the 1st column
      upper.trim( 1UL );

      checkRows    ( upper,  3UL );
      checkColumns ( upper,  3UL );
      checkCapacity( upper, 45UL );
      checkCapacity( upper,  0UL,  0UL );
      checkCapacity( upper,  1UL,  0UL );
      checkCapacity( upper,  2UL, 45UL );

      // Trimming the 2nd column
      upper.trim( 2UL );

      checkRows    ( upper,  3UL );
      checkColumns ( upper,  3UL );
      checkCapacity( upper, 45UL );
      checkCapacity( upper,  0UL, 0UL );
      checkCapacity( upper,  1UL, 0UL );
      checkCapacity( upper,  2UL, 0UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c shrinkToFit() member function of the UpperMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c shrinkToFit() member function of the UpperMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testShrinkToFit()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major UpperMatrix::shrinkToFit()";

      // Shrinking a matrix without excessive capacity
      {
         UT upper( 3UL, 6UL );
         upper(0,0) = 1;
         upper(0,1) = 2;
         upper(0,2) = 3;
         upper(1,1) = 4;
         upper(1,2) = 5;
         upper(2,2) = 6;

         upper.shrinkToFit();

         checkRows    ( upper, 3UL );
         checkColumns ( upper, 3UL );
         checkCapacity( upper, 6UL );
         checkNonZeros( upper, 6UL );
         checkNonZeros( upper, 0UL, 3UL );
         checkNonZeros( upper, 1UL, 2UL );
         checkNonZeros( upper, 2UL, 1UL );

         if( upper.capacity() != upper.nonZeros() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Shrinking the matrix failed\n"
                << " Details:\n"
                << "   Capacity         : " << upper.capacity() << "\n"
                << "   Expected capacity: " << upper.nonZeros() << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( upper(0,0) != 1 || upper(0,1) != 2 || upper(0,2) != 3 ||
             upper(1,0) != 0 || upper(1,1) != 4 || upper(1,2) != 5 ||
             upper(2,0) != 0 || upper(2,1) != 0 || upper(2,2) != 6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Shrinking the matrix failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 1 2 3 )\n( 0 4 5 )\n( 0 0 6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Shrinking a matrix with excessive capacity
      {
         UT upper( 3UL, 100UL );
         upper(0,0) = 1;
         upper(0,1) = 2;
         upper(0,2) = 3;
         upper(1,1) = 4;
         upper(1,2) = 5;
         upper(2,2) = 6;

         upper.shrinkToFit();

         checkRows    ( upper, 3UL );
         checkColumns ( upper, 3UL );
         checkCapacity( upper, 6UL );
         checkNonZeros( upper, 6UL );
         checkNonZeros( upper, 0UL, 3UL );
         checkNonZeros( upper, 1UL, 2UL );
         checkNonZeros( upper, 2UL, 1UL );

         if( upper.capacity() != upper.nonZeros() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Shrinking the matrix failed\n"
                << " Details:\n"
                << "   Capacity         : " << upper.capacity() << "\n"
                << "   Expected capacity: " << upper.nonZeros() << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( upper(0,0) != 1 || upper(0,1) != 2 || upper(0,2) != 3 ||
             upper(1,0) != 0 || upper(1,1) != 4 || upper(1,2) != 5 ||
             upper(2,0) != 0 || upper(2,1) != 0 || upper(2,2) != 6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Shrinking the matrix failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 1 2 3 )\n( 0 4 5 )\n( 0 0 6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major UpperMatrix::shrinkToFit()";

      // Shrinking a matrix without excessive capacity
      {
         OUT upper( 3UL, 6UL );
         upper(0,0) = 1;
         upper(0,1) = 2;
         upper(0,2) = 3;
         upper(1,1) = 4;
         upper(1,2) = 5;
         upper(2,2) = 6;

         upper.shrinkToFit();

         checkRows    ( upper, 3UL );
         checkColumns ( upper, 3UL );
         checkCapacity( upper, 6UL );
         checkNonZeros( upper, 6UL );
         checkNonZeros( upper, 0UL, 1UL );
         checkNonZeros( upper, 1UL, 2UL );
         checkNonZeros( upper, 2UL, 3UL );

         if( upper.capacity() != upper.nonZeros() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Shrinking the matrix failed\n"
                << " Details:\n"
                << "   Capacity         : " << upper.capacity() << "\n"
                << "   Expected capacity: " << upper.nonZeros() << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( upper(0,0) != 1 || upper(0,1) != 2 || upper(0,2) != 3 ||
             upper(1,0) != 0 || upper(1,1) != 4 || upper(1,2) != 5 ||
             upper(2,0) != 0 || upper(2,1) != 0 || upper(2,2) != 6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Shrinking the matrix failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 1 2 3 )\n( 0 4 5 )\n( 0 0 6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Shrinking a matrix with excessive capacity
      {
         OUT upper( 3UL, 100UL );
         upper(0,0) = 1;
         upper(0,1) = 2;
         upper(0,2) = 3;
         upper(1,1) = 4;
         upper(1,2) = 5;
         upper(2,2) = 6;

         upper.shrinkToFit();

         checkRows    ( upper, 3UL );
         checkColumns ( upper, 3UL );
         checkCapacity( upper, 6UL );
         checkNonZeros( upper, 6UL );
         checkNonZeros( upper, 0UL, 1UL );
         checkNonZeros( upper, 1UL, 2UL );
         checkNonZeros( upper, 2UL, 3UL );

         if( upper.capacity() != upper.nonZeros() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Shrinking the matrix failed\n"
                << " Details:\n"
                << "   Capacity         : " << upper.capacity() << "\n"
                << "   Expected capacity: " << upper.nonZeros() << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( upper(0,0) != 1 || upper(0,1) != 2 || upper(0,2) != 3 ||
             upper(1,0) != 0 || upper(1,1) != 4 || upper(1,2) != 5 ||
             upper(2,0) != 0 || upper(2,1) != 0 || upper(2,2) != 6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Shrinking the matrix failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 1 2 3 )\n( 0 4 5 )\n( 0 0 6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c swap() functionality of the UpperMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c swap() function of the UpperMatrix specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testSwap()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major UpperMatrix swap";

      UT upper1( 2UL );
      upper1(0,0) = 1;
      upper1(0,1) = 2;
      upper1(1,1) = 3;

      UT upper2( 2UL );
      upper2(0,0) = 4;
      upper2(0,1) = 5;
      upper2(1,1) = 0;

      swap( upper1, upper2 );

      checkRows    ( upper1, 2UL );
      checkColumns ( upper1, 2UL );
      checkCapacity( upper1, 2UL );
      checkNonZeros( upper1, 2UL );
      checkNonZeros( upper1, 0UL, 2UL );
      checkNonZeros( upper1, 1UL, 0UL );

      if( upper1(0,0) != 4 || upper1(0,1) != 5 || upper1(1,0) != 0 || upper1(1,1) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the first matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << upper1 << "\n"
             << "   Expected result:\n( 4 5 )\n( 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      checkRows    ( upper2, 2UL );
      checkColumns ( upper2, 2UL );
      checkCapacity( upper2, 4UL );
      checkNonZeros( upper2, 3UL );
      checkNonZeros( upper2, 0UL, 2UL );
      checkNonZeros( upper2, 1UL, 1UL );

      if( upper2(0,0) != 1 || upper2(0,1) != 2 || upper2(1,0) != 0 || upper2(1,1) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the second matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << upper2 << "\n"
             << "   Expected result:\n( 1 2 )\n( 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major UpperMatrix swap";

      OUT upper1( 2UL );
      upper1(0,0) = 1;
      upper1(0,1) = 2;
      upper1(1,1) = 3;

      OUT upper2( 2UL );
      upper2(0,0) = 4;
      upper2(0,1) = 5;
      upper2(1,1) = 0;

      swap( upper1, upper2 );

      checkRows    ( upper1, 2UL );
      checkColumns ( upper1, 2UL );
      checkCapacity( upper1, 2UL );
      checkNonZeros( upper1, 2UL );
      checkNonZeros( upper1, 0UL, 1UL );
      checkNonZeros( upper1, 1UL, 1UL );

      if( upper1(0,0) != 4 || upper1(0,1) != 5 || upper1(1,0) != 0 || upper1(1,1) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the first matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << upper1 << "\n"
             << "   Expected result:\n( 4 5 )\n( 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      checkRows    ( upper2, 2UL );
      checkColumns ( upper2, 2UL );
      checkCapacity( upper2, 4UL );
      checkNonZeros( upper2, 3UL );
      checkNonZeros( upper2, 0UL, 1UL );
      checkNonZeros( upper2, 1UL, 2UL );

      if( upper2(0,0) != 1 || upper2(0,1) != 2 || upper2(1,0) != 0 || upper2(1,1) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the second matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << upper2 << "\n"
             << "   Expected result:\n( 1 2 )\n( 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c set() member function of the UpperMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c set() member function of the UpperMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testSet()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major UpperMatrix::set()";

      using Iterator = UT::Iterator;

      // Initialization check
      UT upper( 4UL );

      checkRows    ( upper, 4UL );
      checkColumns ( upper, 4UL );
      checkNonZeros( upper, 0UL );
      checkNonZeros( upper, 0UL, 0UL );
      checkNonZeros( upper, 1UL, 0UL );
      checkNonZeros( upper, 2UL, 0UL );
      checkNonZeros( upper, 3UL, 0UL );

      // Setting a non-zero element
      {
         Iterator pos = upper.set( 1UL, 2UL, 1 );

         checkRows    ( upper, 4UL );
         checkColumns ( upper, 4UL );
         checkCapacity( upper, 1UL );
         checkNonZeros( upper, 1UL );
         checkNonZeros( upper, 0UL, 0UL );
         checkNonZeros( upper, 1UL, 1UL );
         checkNonZeros( upper, 2UL, 0UL );
         checkNonZeros( upper, 3UL, 0UL );

         if( pos->value() != 1 || pos->index() != 2UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 1\n"
                << "   Expected index: 2\n";
            throw std::runtime_error( oss.str() );
         }

         if( upper(1,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 1 0 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setting a second non-zero element
      {
         Iterator pos = upper.set( 1UL, 1UL, 2 );

         checkRows    ( upper, 4UL );
         checkColumns ( upper, 4UL );
         checkCapacity( upper, 2UL );
         checkNonZeros( upper, 2UL );
         checkNonZeros( upper, 0UL, 0UL );
         checkNonZeros( upper, 1UL, 2UL );
         checkNonZeros( upper, 2UL, 0UL );
         checkNonZeros( upper, 3UL, 0UL );

         if( pos->value() != 2 || pos->index() != 1UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 2\n"
                << "   Expected index: 1\n";
            throw std::runtime_error( oss.str() );
         }

         if( upper(1,1) != 2 || upper(1,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 2 1 0 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setting a third non-zero element
      {
         Iterator pos = upper.set( 1UL, 3UL, 3 );

         checkRows    ( upper, 4UL );
         checkColumns ( upper, 4UL );
         checkCapacity( upper, 3UL );
         checkNonZeros( upper, 3UL );
         checkNonZeros( upper, 0UL, 0UL );
         checkNonZeros( upper, 1UL, 3UL );
         checkNonZeros( upper, 2UL, 0UL );
         checkNonZeros( upper, 3UL, 0UL );

         if( pos->value() != 3 || pos->index() != 3UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 3\n"
                << "   Expected index: 3\n";
            throw std::runtime_error( oss.str() );
         }

         if( upper(1,1) != 2 || upper(1,2) != 1 || upper(1,3) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 2 1 3 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setting an already existing element
      {
         Iterator pos = upper.set( 1UL, 2UL, 4 );

         checkRows    ( upper, 4UL );
         checkColumns ( upper, 4UL );
         checkCapacity( upper, 3UL );
         checkNonZeros( upper, 3UL );
         checkNonZeros( upper, 0UL, 0UL );
         checkNonZeros( upper, 1UL, 3UL );
         checkNonZeros( upper, 2UL, 0UL );
         checkNonZeros( upper, 3UL, 0UL );

         if( pos->value() != 4 || pos->index() != 2UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 4\n"
                << "   Expected index: 2\n";
            throw std::runtime_error( oss.str() );
         }

         if( upper(1,1) != 2 || upper(1,2) != 4 || upper(1,3) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 2 4 3 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major UpperMatrix::set()";

      using Iterator = OUT::Iterator;

      // Initialization check
      OUT upper( 4UL );

      checkRows    ( upper, 4UL );
      checkColumns ( upper, 4UL );
      checkNonZeros( upper, 0UL );
      checkNonZeros( upper, 0UL, 0UL );
      checkNonZeros( upper, 1UL, 0UL );
      checkNonZeros( upper, 2UL, 0UL );
      checkNonZeros( upper, 3UL, 0UL );

      // Setting a non-zero element
      {
         Iterator pos = upper.set( 1UL, 2UL, 1 );

         checkRows    ( upper, 4UL );
         checkColumns ( upper, 4UL );
         checkCapacity( upper, 1UL );
         checkNonZeros( upper, 1UL );
         checkNonZeros( upper, 0UL, 0UL );
         checkNonZeros( upper, 1UL, 0UL );
         checkNonZeros( upper, 2UL, 1UL );
         checkNonZeros( upper, 3UL, 0UL );

         if( pos->value() != 1 || pos->index() != 1UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 1\n"
                << "   Expected index: 1\n";
            throw std::runtime_error( oss.str() );
         }

         if( upper(1,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 1 0 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setting a second non-zero element
      {
         Iterator pos = upper.set( 2UL, 2UL, 2 );

         checkRows    ( upper, 4UL );
         checkColumns ( upper, 4UL );
         checkCapacity( upper, 2UL );
         checkNonZeros( upper, 2UL );
         checkNonZeros( upper, 0UL, 0UL );
         checkNonZeros( upper, 1UL, 0UL );
         checkNonZeros( upper, 2UL, 2UL );
         checkNonZeros( upper, 3UL, 0UL );

         if( pos->value() != 2 || pos->index() != 2UL ) {
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

         if( upper(1,2) != 1 || upper(2,2) != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 1 0 )\n( 0 0 2 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setting a third non-zero element
      {
         Iterator pos = upper.set( 0UL, 2UL, 3 );

         checkRows    ( upper, 4UL );
         checkColumns ( upper, 4UL );
         checkCapacity( upper, 3UL );
         checkNonZeros( upper, 3UL );
         checkNonZeros( upper, 0UL, 0UL );
         checkNonZeros( upper, 1UL, 0UL );
         checkNonZeros( upper, 2UL, 3UL );
         checkNonZeros( upper, 3UL, 0UL );

         if( pos->value() != 3 || pos->index() != 0UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 3\n"
                << "   Expected index: 0\n";
            throw std::runtime_error( oss.str() );
         }

         if( upper(0,2) != 3 || upper(1,2) != 1 || upper(2,2) != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 0 0 3 0 )\n( 0 0 1 0 )\n( 0 0 2 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setting an already existing element
      {
         Iterator pos = upper.set( 1UL, 2UL, 4 );

         checkRows    ( upper, 4UL );
         checkColumns ( upper, 4UL );
         checkCapacity( upper, 3UL );
         checkNonZeros( upper, 3UL );
         checkNonZeros( upper, 0UL, 0UL );
         checkNonZeros( upper, 1UL, 0UL );
         checkNonZeros( upper, 2UL, 3UL );
         checkNonZeros( upper, 3UL, 0UL );

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

         if( upper(0,2) != 3 || upper(1,2) != 4 || upper(2,2) != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 0 0 3 0 )\n( 0 0 4 0 )\n( 0 0 2 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c insert() member function of the UpperMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c insert() member function of the UpperMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testInsert()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major UpperMatrix::insert()";

      using Iterator = UT::Iterator;

      // Initialization check
      UT upper( 4UL );

      checkRows    ( upper, 4UL );
      checkColumns ( upper, 4UL );
      checkNonZeros( upper, 0UL );
      checkNonZeros( upper, 0UL, 0UL );
      checkNonZeros( upper, 1UL, 0UL );
      checkNonZeros( upper, 2UL, 0UL );
      checkNonZeros( upper, 3UL, 0UL );

      // Inserting a non-zero element
      {
         Iterator pos = upper.insert( 1UL, 2UL, 1 );

         checkRows    ( upper, 4UL );
         checkColumns ( upper, 4UL );
         checkCapacity( upper, 1UL );
         checkNonZeros( upper, 1UL );
         checkNonZeros( upper, 0UL, 0UL );
         checkNonZeros( upper, 1UL, 1UL );
         checkNonZeros( upper, 2UL, 0UL );
         checkNonZeros( upper, 3UL, 0UL );

         if( pos->value() != 1 || pos->index() != 2UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 1\n"
                << "   Expected index: 2\n";
            throw std::runtime_error( oss.str() );
         }

         if( upper(1,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 1 0 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting a second non-zero element
      {
         Iterator pos = upper.insert( 1UL, 1UL, 2 );

         checkRows    ( upper, 4UL );
         checkColumns ( upper, 4UL );
         checkCapacity( upper, 2UL );
         checkNonZeros( upper, 2UL );
         checkNonZeros( upper, 0UL, 0UL );
         checkNonZeros( upper, 1UL, 2UL );
         checkNonZeros( upper, 2UL, 0UL );
         checkNonZeros( upper, 3UL, 0UL );

         if( pos->value() != 2 || pos->index() != 1UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 2\n"
                << "   Expected index: 1\n";
            throw std::runtime_error( oss.str() );
         }

         if( upper(1,1) != 2 || upper(1,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 2 1 0 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting a third non-zero element
      {
         Iterator pos = upper.insert( 1UL, 3UL, 3 );

         checkRows    ( upper, 4UL );
         checkColumns ( upper, 4UL );
         checkCapacity( upper, 3UL );
         checkNonZeros( upper, 3UL );
         checkNonZeros( upper, 0UL, 0UL );
         checkNonZeros( upper, 1UL, 3UL );
         checkNonZeros( upper, 2UL, 0UL );
         checkNonZeros( upper, 3UL, 0UL );

         if( pos->value() != 3 || pos->index() != 3UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 3\n"
                << "   Expected index: 3\n";
            throw std::runtime_error( oss.str() );
         }

         if( upper(1,1) != 2 || upper(1,2) != 1 || upper(1,3) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 2 1 3 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to insert an already existing element
      try {
         upper.insert( 2UL, 1UL, 4 );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Inserting an existing element succeeded\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 0 0 0 0 )\n( 0 2 1 3 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major UpperMatrix::insert()";

      using Iterator = OUT::Iterator;

      // Initialization check
      OUT upper( 4UL );

      checkRows    ( upper, 4UL );
      checkColumns ( upper, 4UL );
      checkNonZeros( upper, 0UL );
      checkNonZeros( upper, 0UL, 0UL );
      checkNonZeros( upper, 1UL, 0UL );
      checkNonZeros( upper, 2UL, 0UL );
      checkNonZeros( upper, 3UL, 0UL );

      // Inserting a non-zero element
      {
         Iterator pos = upper.insert( 1UL, 2UL, 1 );

         checkRows    ( upper, 4UL );
         checkColumns ( upper, 4UL );
         checkCapacity( upper, 1UL );
         checkNonZeros( upper, 1UL );
         checkNonZeros( upper, 0UL, 0UL );
         checkNonZeros( upper, 1UL, 0UL );
         checkNonZeros( upper, 2UL, 1UL );
         checkNonZeros( upper, 3UL, 0UL );

         if( pos->value() != 1 || pos->index() != 1UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 1\n"
                << "   Expected index: 1\n";
            throw std::runtime_error( oss.str() );
         }

         if( upper(1,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 1 0 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting a second non-zero element
      {
         Iterator pos = upper.insert( 2UL, 2UL, 2 );

         checkRows    ( upper, 4UL );
         checkColumns ( upper, 4UL );
         checkCapacity( upper, 2UL );
         checkNonZeros( upper, 2UL );
         checkNonZeros( upper, 0UL, 0UL );
         checkNonZeros( upper, 1UL, 0UL );
         checkNonZeros( upper, 2UL, 2UL );
         checkNonZeros( upper, 3UL, 0UL );

         if( pos->value() != 2 || pos->index() != 2UL ) {
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

         if( upper(1,2) != 1 || upper(2,2) != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 1 0 )\n( 0 0 2 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting a third non-zero element
      {
         Iterator pos = upper.insert( 0UL, 2UL, 3 );

         checkRows    ( upper, 4UL );
         checkColumns ( upper, 4UL );
         checkCapacity( upper, 3UL );
         checkNonZeros( upper, 3UL );
         checkNonZeros( upper, 0UL, 0UL );
         checkNonZeros( upper, 1UL, 0UL );
         checkNonZeros( upper, 2UL, 3UL );
         checkNonZeros( upper, 3UL, 0UL );

         if( pos->value() != 3 || pos->index() != 0UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 3\n"
                << "   Expected index: 0\n";
            throw std::runtime_error( oss.str() );
         }

         if( upper(0,2) != 3 || upper(1,2) != 1 || upper(2,2) != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 0 0 3 0 )\n( 0 0 1 0 )\n( 0 0 2 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to insert an already existing element
      try {
         upper.insert( 1UL, 2UL, 4 );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Inserting an existing element succeeded\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 0 0 3 0 )\n( 0 0 1 0 )\n( 0 0 2 0 )\n( 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c append() member function of the UpperMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c append() member function of the UpperMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testAppend()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major UpperMatrix::append()";

      // Appending with pre-allocation in each row
      {
         // Initialization check
         UT upper( 4UL, 5UL );
         upper.reserve( 0UL, 2UL );
         upper.reserve( 1UL, 1UL );
         upper.reserve( 2UL, 2UL );

         checkRows    ( upper, 4UL );
         checkColumns ( upper, 4UL );
         checkNonZeros( upper, 0UL );
         checkNonZeros( upper, 0UL, 0UL );
         checkNonZeros( upper, 1UL, 0UL );
         checkNonZeros( upper, 2UL, 0UL );
         checkNonZeros( upper, 3UL, 0UL );

         // Appending one non-zero element
         upper.append( 1UL, 2UL, 1 );

         checkRows    ( upper, 4UL );
         checkColumns ( upper, 4UL );
         checkCapacity( upper, 1UL );
         checkNonZeros( upper, 1UL );
         checkNonZeros( upper, 0UL, 0UL );
         checkNonZeros( upper, 1UL, 1UL );
         checkNonZeros( upper, 2UL, 0UL );
         checkNonZeros( upper, 3UL, 0UL );

         if( upper(1,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Appending operation failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 1 0 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Appending two more non-zero elements
         upper.append( 0UL, 0UL, 2 );
         upper.append( 2UL, 2UL, 3 );

         checkRows    ( upper, 4UL );
         checkColumns ( upper, 4UL );
         checkCapacity( upper, 3UL );
         checkNonZeros( upper, 3UL );
         checkNonZeros( upper, 0UL, 1UL );
         checkNonZeros( upper, 1UL, 1UL );
         checkNonZeros( upper, 2UL, 1UL );
         checkNonZeros( upper, 3UL, 0UL );

         if( upper(0,0) != 2 || upper(1,2) != 1 || upper(2,2) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 2 0 0 0 )\n( 0 0 1 0 )\n( 0 0 3 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Appending two more non-zero elements
         upper.append( 0UL, 3UL, 4 );
         upper.append( 2UL, 3UL, 5 );

         checkRows    ( upper, 4UL );
         checkColumns ( upper, 4UL );
         checkCapacity( upper, 5UL );
         checkNonZeros( upper, 5UL );
         checkNonZeros( upper, 0UL, 2UL );
         checkNonZeros( upper, 1UL, 1UL );
         checkNonZeros( upper, 2UL, 2UL );
         checkNonZeros( upper, 3UL, 0UL );

         if( upper(0,0) != 2 || upper(0,3) != 4 || upper(1,2) != 1 || upper(2,2) != 3 || upper(2,3) != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 2 0 0 4 )\n( 0 0 1 0 )\n( 0 0 3 5 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Appending with row finalization
      {
         // Initialization check
         UT upper( 4UL, 5UL );
         upper.reserve( 0UL, 1UL );
         upper.reserve( 1UL, 2UL );
         upper.reserve( 2UL, 2UL );

         // Appending one non-zero element
         upper.append( 0UL, 0UL, 1 );
         upper.finalize( 0UL );

         checkRows    ( upper, 4UL );
         checkColumns ( upper, 4UL );
         checkCapacity( upper, 1UL );
         checkNonZeros( upper, 1UL );
         checkNonZeros( upper, 0UL, 1UL );
         checkNonZeros( upper, 1UL, 0UL );
         checkNonZeros( upper, 2UL, 0UL );
         checkNonZeros( upper, 3UL, 0UL );

         if( upper(0,0) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Appending two more non-zero elements
         upper.append( 1UL, 1UL, 2 );
         upper.append( 1UL, 3UL, 3 );
         upper.finalize( 1UL );

         checkRows    ( upper, 4UL );
         checkColumns ( upper, 4UL );
         checkCapacity( upper, 3UL );
         checkNonZeros( upper, 3UL );
         checkNonZeros( upper, 0UL, 1UL );
         checkNonZeros( upper, 1UL, 2UL );
         checkNonZeros( upper, 2UL, 0UL );
         checkNonZeros( upper, 3UL, 0UL );

         if( upper(0,0) != 1 || upper(1,1) != 2 || upper(1,3) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 0 2 0 3 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Appending two more non-zero elements
         upper.append( 2UL, 2UL, 4 );
         upper.append( 2UL, 3UL, 5 );
         upper.finalize( 2UL );

         checkRows    ( upper, 4UL );
         checkColumns ( upper, 4UL );
         checkCapacity( upper, 5UL );
         checkNonZeros( upper, 5UL );
         checkNonZeros( upper, 0UL, 1UL );
         checkNonZeros( upper, 1UL, 2UL );
         checkNonZeros( upper, 2UL, 2UL );
         checkNonZeros( upper, 3UL, 0UL );

         if( upper(0,0) != 1 || upper(1,1) != 2 || upper(1,3) != 3 || upper(2,2) != 4 || upper(2,3) != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 0 2 0 3 )\n( 0 0 4 5 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major UpperMatrix::append()";

      // Appending with pre-allocation in each column
      {
         // Initialization check
         OUT upper( 4UL, 5UL );
         upper.reserve( 0UL, 1UL );
         upper.reserve( 2UL, 2UL );
         upper.reserve( 3UL, 2UL );

         checkRows    ( upper, 4UL );
         checkColumns ( upper, 4UL );
         checkNonZeros( upper, 0UL );
         checkNonZeros( upper, 0UL, 0UL );
         checkNonZeros( upper, 1UL, 0UL );
         checkNonZeros( upper, 2UL, 0UL );
         checkNonZeros( upper, 3UL, 0UL );

         // Appending one non-zero element
         upper.append( 1UL, 2UL, 1 );

         checkRows    ( upper, 4UL );
         checkColumns ( upper, 4UL );
         checkCapacity( upper, 1UL );
         checkNonZeros( upper, 1UL );
         checkNonZeros( upper, 0UL, 0UL );
         checkNonZeros( upper, 1UL, 0UL );
         checkNonZeros( upper, 2UL, 1UL );
         checkNonZeros( upper, 3UL, 0UL );

         if( upper(1,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 1 0 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Appending two more non-zero elements
         upper.append( 0UL, 0UL, 2 );
         upper.append( 0UL, 3UL, 3 );

         checkRows    ( upper, 4UL );
         checkColumns ( upper, 4UL );
         checkCapacity( upper, 3UL );
         checkNonZeros( upper, 3UL );
         checkNonZeros( upper, 0UL, 1UL );
         checkNonZeros( upper, 1UL, 0UL );
         checkNonZeros( upper, 2UL, 1UL );
         checkNonZeros( upper, 3UL, 1UL );

         if( upper(0,0) != 2 || upper(0,3) != 3 || upper(1,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 2 0 0 3 )\n( 0 0 1 0 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Appending two more non-zero elements
         upper.append( 2UL, 2UL, 4 );
         upper.append( 2UL, 3UL, 5 );

         checkRows    ( upper, 4UL );
         checkColumns ( upper, 4UL );
         checkCapacity( upper, 5UL );
         checkNonZeros( upper, 5UL );
         checkNonZeros( upper, 0UL, 1UL );
         checkNonZeros( upper, 1UL, 0UL );
         checkNonZeros( upper, 2UL, 2UL );
         checkNonZeros( upper, 3UL, 2UL );

         if( upper(0,0) != 2 || upper(0,3) != 3 || upper(1,2) != 1 || upper(2,2) != 4 || upper(2,3) != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 2 0 0 3 )\n( 0 0 1 0 )\n( 0 0 4 5 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Appending with column finalization
      {
         // Initialization check
         OUT upper( 4UL, 5UL );
         upper.reserve( 0UL, 1UL );
         upper.reserve( 2UL, 2UL );
         upper.reserve( 3UL, 2UL );

         // Appending one non-zero element
         upper.append( 0UL, 0UL, 1 );
         upper.finalize( 0UL );

         checkRows    ( upper, 4UL );
         checkColumns ( upper, 4UL );
         checkCapacity( upper, 1UL );
         checkNonZeros( upper, 1UL );
         checkNonZeros( upper, 0UL, 1UL );
         checkNonZeros( upper, 1UL, 0UL );
         checkNonZeros( upper, 2UL, 0UL );
         checkNonZeros( upper, 3UL, 0UL );

         if( upper(0,0) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Appending two more non-zero elements
         upper.append( 0UL, 2UL, 2 );
         upper.append( 1UL, 2UL, 3 );
         upper.finalize( 2UL );

         checkRows    ( upper, 4UL );
         checkColumns ( upper, 4UL );
         checkCapacity( upper, 3UL );
         checkNonZeros( upper, 3UL );
         checkNonZeros( upper, 0UL, 1UL );
         checkNonZeros( upper, 1UL, 0UL );
         checkNonZeros( upper, 2UL, 2UL );
         checkNonZeros( upper, 3UL, 0UL );

         if( upper(0,0) != 1 || upper(0,2) != 2 || upper(1,2) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 1 0 2 0 )\n( 0 0 3 0 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Appending two more non-zero elements
         upper.append( 0UL, 3UL, 4 );
         upper.append( 2UL, 3UL, 5 );
         upper.finalize( 3UL );

         checkRows    ( upper, 4UL );
         checkColumns ( upper, 4UL );
         checkCapacity( upper, 5UL );
         checkNonZeros( upper, 5UL );
         checkNonZeros( upper, 0UL, 1UL );
         checkNonZeros( upper, 1UL, 0UL );
         checkNonZeros( upper, 2UL, 2UL );
         checkNonZeros( upper, 3UL, 2UL );

         if( upper(0,0) != 1 || upper(0,2) != 2 || upper(0,3) != 4 || upper(1,2) != 3 || upper(2,3) != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 1 0 2 4 )\n( 0 0 3 0 )\n( 0 0 0 5 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c erase() member function of the UpperMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c erase() member function of the UpperMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testErase()
{
   //=====================================================================================
   // Row-major index-based erase function
   //=====================================================================================

   {
      test_ = "Row-major UpperMatrix::erase( size_t, size_t )";

      // Initialization check
      UT upper( 4UL, 8UL );
      upper(0,0) = 1;
      upper(0,2) = 2;
      upper(0,3) = 3;
      upper(1,2) = 4;
      upper(1,3) = 5;
      upper(2,2) = 6;
      upper(2,3) = 7;
      upper(3,3) = 8;

      checkRows    ( upper, 4UL );
      checkColumns ( upper, 4UL );
      checkCapacity( upper, 8UL );
      checkNonZeros( upper, 8UL );
      checkNonZeros( upper, 0UL, 3UL );
      checkNonZeros( upper, 1UL, 2UL );
      checkNonZeros( upper, 2UL, 2UL );
      checkNonZeros( upper, 3UL, 1UL );

      if( upper(0,0) != 1 || upper(0,2) != 2 || upper(0,3) != 3 ||
          upper(1,2) != 4 || upper(1,3) != 5 ||
          upper(2,2) != 6 || upper(2,3) != 7 ||
          upper(3,3) != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 0 2 3 )\n( 0 0 4 5 )\n( 0 0 6 7 )\n( 0 0 0 8 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element at (1,2)
      upper.erase( 1UL, 2UL );

      checkRows    ( upper, 4UL );
      checkColumns ( upper, 4UL );
      checkCapacity( upper, 8UL );
      checkNonZeros( upper, 7UL );
      checkNonZeros( upper, 0UL, 3UL );
      checkNonZeros( upper, 1UL, 1UL );
      checkNonZeros( upper, 2UL, 2UL );
      checkNonZeros( upper, 3UL, 1UL );

      if( upper(0,0) != 1 || upper(0,2) != 2 || upper(0,3) != 3 ||
          upper(1,3) != 5 ||
          upper(2,2) != 6 || upper(2,3) != 7 ||
          upper(3,3) != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 0 2 3 )\n( 0 0 0 5 )\n( 0 0 6 7 )\n( 0 0 0 8 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element at (2,3)
      upper.erase( 2UL, 3UL );

      checkRows    ( upper, 4UL );
      checkColumns ( upper, 4UL );
      checkCapacity( upper, 8UL );
      checkNonZeros( upper, 6UL );
      checkNonZeros( upper, 0UL, 3UL );
      checkNonZeros( upper, 1UL, 1UL );
      checkNonZeros( upper, 2UL, 1UL );
      checkNonZeros( upper, 3UL, 1UL );

      if( upper(0,0) != 1 || upper(0,2) != 2 || upper(0,3) != 3 ||
          upper(1,3) != 5 ||
          upper(2,2) != 6 ||
          upper(3,3) != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 0 2 3 )\n( 0 0 0 5 )\n( 0 0 6 0 )\n( 0 0 0 8 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element at (0,2)
      upper.erase( 0UL, 2UL );

      checkRows    ( upper, 4UL );
      checkColumns ( upper, 4UL );
      checkCapacity( upper, 8UL );
      checkNonZeros( upper, 5UL );
      checkNonZeros( upper, 0UL, 2UL );
      checkNonZeros( upper, 1UL, 1UL );
      checkNonZeros( upper, 2UL, 1UL );
      checkNonZeros( upper, 3UL, 1UL );

      if( upper(0,0) != 1 || upper(0,3) != 3 ||
          upper(1,3) != 5 ||
          upper(2,2) != 6 ||
          upper(3,3) != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 0 0 3 )\n( 0 0 0 5 )\n( 0 0 6 0 )\n( 0 0 0 8 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to erase a zero element
      upper.erase( 0UL, 1UL );

      checkRows    ( upper, 4UL );
      checkColumns ( upper, 4UL );
      checkCapacity( upper, 8UL );
      checkNonZeros( upper, 5UL );
      checkNonZeros( upper, 0UL, 2UL );
      checkNonZeros( upper, 1UL, 1UL );
      checkNonZeros( upper, 2UL, 1UL );
      checkNonZeros( upper, 3UL, 1UL );

      if( upper(0,0) != 1 || upper(0,3) != 3 ||
          upper(1,3) != 5 ||
          upper(2,2) != 6 ||
          upper(3,3) != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 0 0 3 )\n( 0 0 0 5 )\n( 0 0 6 0 )\n( 0 0 0 8 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major iterator-based erase function
   //=====================================================================================

   {
      test_ = "Row-major UpperMatrix::erase( size_t, Iterator )";

      using Iterator = UT::Iterator;

      // Initialization check
      UT upper( 4UL, 8UL );
      upper(0,0) = 1;
      upper(0,2) = 2;
      upper(0,3) = 3;
      upper(1,2) = 4;
      upper(1,3) = 5;
      upper(2,2) = 6;
      upper(2,3) = 7;
      upper(3,3) = 8;

      checkRows    ( upper, 4UL );
      checkColumns ( upper, 4UL );
      checkCapacity( upper, 8UL );
      checkNonZeros( upper, 8UL );
      checkNonZeros( upper, 0UL, 3UL );
      checkNonZeros( upper, 1UL, 2UL );
      checkNonZeros( upper, 2UL, 2UL );
      checkNonZeros( upper, 3UL, 1UL );

      if( upper(0,0) != 1 || upper(0,2) != 2 || upper(0,3) != 3 ||
          upper(1,2) != 4 || upper(1,3) != 5 ||
          upper(2,2) != 6 || upper(2,3) != 7 ||
          upper(3,3) != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 0 2 3 )\n( 0 0 4 5 )\n( 0 0 6 7 )\n( 0 0 0 8 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element at (1,2)
      {
         Iterator pos = upper.erase( 1UL, upper.find( 1UL, 2UL ) );

         checkRows    ( upper, 4UL );
         checkColumns ( upper, 4UL );
         checkCapacity( upper, 8UL );
         checkNonZeros( upper, 7UL );
         checkNonZeros( upper, 0UL, 3UL );
         checkNonZeros( upper, 1UL, 1UL );
         checkNonZeros( upper, 2UL, 2UL );
         checkNonZeros( upper, 3UL, 1UL );

         if( upper(0,0) != 1 || upper(0,2) != 2 || upper(0,3) != 3 ||
             upper(1,3) != 5 ||
             upper(2,2) != 6 || upper(2,3) != 7 ||
             upper(3,3) != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 1 0 2 3 )\n( 0 0 0 5 )\n( 0 0 6 7 )\n( 0 0 0 8 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos->value() != 5 || pos->index() != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 5\n"
                << "   Expected index: 3\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the element at (2,3)
      {
         Iterator pos = upper.erase( 2UL, upper.find( 2UL, 3UL ) );

         checkRows    ( upper, 4UL );
         checkColumns ( upper, 4UL );
         checkCapacity( upper, 8UL );
         checkNonZeros( upper, 6UL );
         checkNonZeros( upper, 0UL, 3UL );
         checkNonZeros( upper, 1UL, 1UL );
         checkNonZeros( upper, 2UL, 1UL );
         checkNonZeros( upper, 3UL, 1UL );

         if( upper(0,0) != 1 || upper(0,2) != 2 || upper(0,3) != 3 ||
             upper(1,3) != 5 ||
             upper(2,2) != 6 ||
             upper(3,3) != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 1 0 2 3 )\n( 0 0 0 5 )\n( 0 0 6 0 )\n( 0 0 0 8 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos != upper.end( 2UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the element at (0,2)
      {
         Iterator pos = upper.erase( 0UL, upper.find( 0UL, 2UL ) );

         checkRows    ( upper, 4UL );
         checkColumns ( upper, 4UL );
         checkCapacity( upper, 8UL );
         checkNonZeros( upper, 5UL );
         checkNonZeros( upper, 0UL, 2UL );
         checkNonZeros( upper, 1UL, 1UL );
         checkNonZeros( upper, 2UL, 1UL );
         checkNonZeros( upper, 3UL, 1UL );

         if( upper(0,0) != 1 || upper(0,3) != 3 ||
             upper(1,3) != 5 ||
             upper(2,2) != 6 ||
             upper(3,3) != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 1 0 0 3 )\n( 0 0 0 5 )\n( 0 0 6 0 )\n( 0 0 0 8 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos->value() != 3 || pos->index() != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 3\n"
                << "   Expected index: 3\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to erase a zero element
      {
         Iterator pos = upper.erase( 0UL, upper.find( 0UL, 1UL ) );

         checkRows    ( upper, 4UL );
         checkColumns ( upper, 4UL );
         checkCapacity( upper, 8UL );
         checkNonZeros( upper, 5UL );
         checkNonZeros( upper, 0UL, 2UL );
         checkNonZeros( upper, 1UL, 1UL );
         checkNonZeros( upper, 2UL, 1UL );
         checkNonZeros( upper, 3UL, 1UL );

         if( upper(0,0) != 1 || upper(0,3) != 3 ||
             upper(1,3) != 5 ||
             upper(2,2) != 6 ||
             upper(3,3) != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 1 0 0 3 )\n( 0 0 0 5 )\n( 0 0 6 0 )\n( 0 0 0 8 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos != upper.end( 0UL ) ) {
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
      test_ = "Row-major UpperMatrix::erase( size_t, Iterator, Iterator )";

      using Iterator = UT::Iterator;

      // Initialization check
      UT upper( 4UL, 8UL );
      upper(0,0) = 1;
      upper(0,2) = 2;
      upper(0,3) = 3;
      upper(1,2) = 4;
      upper(1,3) = 5;
      upper(2,2) = 6;
      upper(2,3) = 7;
      upper(3,3) = 8;

      checkRows    ( upper, 4UL );
      checkColumns ( upper, 4UL );
      checkCapacity( upper, 8UL );
      checkNonZeros( upper, 8UL );
      checkNonZeros( upper, 0UL, 3UL );
      checkNonZeros( upper, 1UL, 2UL );
      checkNonZeros( upper, 2UL, 2UL );
      checkNonZeros( upper, 3UL, 1UL );

      if( upper(0,0) != 1 || upper(0,2) != 2 || upper(0,3) != 3 ||
          upper(1,2) != 4 || upper(1,3) != 5 ||
          upper(2,2) != 6 || upper(2,3) != 7 ||
          upper(3,3) != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 0 2 3 )\n( 0 0 4 5 )\n( 0 0 6 7 )\n( 0 0 0 8 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the elements from the beginning of row 1 to (1,3)
      {
         Iterator pos = upper.erase( 1UL, upper.begin( 1UL ), upper.find( 1UL, 3UL ) );

         checkRows    ( upper, 4UL );
         checkColumns ( upper, 4UL );
         checkCapacity( upper, 8UL );
         checkNonZeros( upper, 7UL );
         checkNonZeros( upper, 0UL, 3UL );
         checkNonZeros( upper, 1UL, 1UL );
         checkNonZeros( upper, 2UL, 2UL );
         checkNonZeros( upper, 3UL, 1UL );

         if( upper(0,0) != 1 || upper(0,2) != 2 || upper(0,3) != 3 ||
             upper(1,3) != 5 ||
             upper(2,2) != 6 || upper(2,3) != 7 ||
             upper(3,3) != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a single-element range failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 1 0 2 3 )\n( 0 0 0 5 )\n( 0 0 6 7 )\n( 0 0 0 8 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos->value() != 5 || pos->index() != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 5\n"
                << "   Expected index: 3\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the elements from (2,3) to the row end
      {
         Iterator pos = upper.erase( 2UL, upper.find( 2UL, 3UL ), upper.end( 2UL ) );

         checkRows    ( upper, 4UL );
         checkColumns ( upper, 4UL );
         checkCapacity( upper, 8UL );
         checkNonZeros( upper, 6UL );
         checkNonZeros( upper, 0UL, 3UL );
         checkNonZeros( upper, 1UL, 1UL );
         checkNonZeros( upper, 2UL, 1UL );
         checkNonZeros( upper, 3UL, 1UL );

         if( upper(0,0) != 1 || upper(0,2) != 2 || upper(0,3) != 3 ||
             upper(1,3) != 5 ||
             upper(2,2) != 6 ||
             upper(3,3) != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a single-element range failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 1 0 2 3 )\n( 0 0 0 5 )\n( 0 0 6 0 )\n( 0 0 0 8 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos != upper.end( 2UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the elements from (0,0) to (0,3)
      {
         Iterator pos = upper.erase( 0UL, upper.find( 0UL, 0UL ), upper.find( 0UL, 3UL ) );

         checkRows    ( upper, 4UL );
         checkColumns ( upper, 4UL );
         checkCapacity( upper, 8UL );
         checkNonZeros( upper, 4UL );
         checkNonZeros( upper, 0UL, 1UL );
         checkNonZeros( upper, 1UL, 1UL );
         checkNonZeros( upper, 2UL, 1UL );
         checkNonZeros( upper, 3UL, 1UL );

         if( upper(0,3) != 3 ||
             upper(1,3) != 5 ||
             upper(2,2) != 6 ||
             upper(3,3) != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a multi-element range failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 0 0 0 3 )\n( 0 0 0 5 )\n( 0 0 6 0 )\n( 0 0 0 8 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos->value() != 3 || pos->index() != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 3\n"
                << "   Expected index: 3\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to erase an empty range
      {
         Iterator pos = upper.erase( 0UL, upper.find( 0UL, 3UL ), upper.find( 0UL, 3UL ) );

         checkRows    ( upper, 4UL );
         checkColumns ( upper, 4UL );
         checkCapacity( upper, 8UL );
         checkNonZeros( upper, 4UL );
         checkNonZeros( upper, 0UL, 1UL );
         checkNonZeros( upper, 1UL, 1UL );
         checkNonZeros( upper, 2UL, 1UL );
         checkNonZeros( upper, 3UL, 1UL );

         if( upper(0,3) != 3 ||
             upper(1,3) != 5 ||
             upper(2,2) != 6 ||
             upper(3,3) != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing an empty range failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 0 0 0 3 )\n( 0 0 0 5 )\n( 0 0 6 0 )\n( 0 0 0 8 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos->value() != 3 || pos->index() != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 3\n"
                << "   Expected index: 3\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Row-major erase function with predicate
   //=====================================================================================

   {
      test_ = "Row-major UpperMatrix::erase( Predicate )";

      // Initialization check
      UT upper( 4UL, 8UL );
      upper(0,0) = 1;
      upper(0,2) = 2;
      upper(0,3) = 3;
      upper(1,2) = 4;
      upper(1,3) = 5;
      upper(2,2) = 6;
      upper(2,3) = 7;
      upper(3,3) = 8;

      checkRows    ( upper, 4UL );
      checkColumns ( upper, 4UL );
      checkCapacity( upper, 8UL );
      checkNonZeros( upper, 8UL );
      checkNonZeros( upper, 0UL, 3UL );
      checkNonZeros( upper, 1UL, 2UL );
      checkNonZeros( upper, 2UL, 2UL );
      checkNonZeros( upper, 3UL, 1UL );

      if( upper(0,0) != 1 || upper(0,2) != 2 || upper(0,3) != 3 ||
          upper(1,2) != 4 || upper(1,3) != 5 ||
          upper(2,2) != 6 || upper(2,3) != 7 ||
          upper(3,3) != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 0 2 3 )\n( 0 0 4 5 )\n( 0 0 6 7 )\n( 0 0 0 8 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing a selection of elements
      upper.erase( []( int value ){ return value == 1 || value == 4 || value == 7; } );

      checkRows    ( upper, 4UL );
      checkColumns ( upper, 4UL );
      checkCapacity( upper, 8UL );
      checkNonZeros( upper, 5UL );
      checkNonZeros( upper, 0UL, 2UL );
      checkNonZeros( upper, 1UL, 1UL );
      checkNonZeros( upper, 2UL, 1UL );
      checkNonZeros( upper, 3UL, 1UL );

      if( upper(0,2) != 2 || upper(0,3) != 3 ||
          upper(1,3) != 5 ||
          upper(2,2) != 6 ||
          upper(3,3) != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a selection of elements failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 0 0 2 3 )\n( 0 0 0 5 )\n( 0 0 6 0 )\n( 0 0 0 8 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to erase all elements with value 1
      upper.erase( []( int value ){ return value == 1; } );

      checkRows    ( upper, 4UL );
      checkColumns ( upper, 4UL );
      checkCapacity( upper, 8UL );
      checkNonZeros( upper, 5UL );
      checkNonZeros( upper, 0UL, 2UL );
      checkNonZeros( upper, 1UL, 1UL );
      checkNonZeros( upper, 2UL, 1UL );
      checkNonZeros( upper, 3UL, 1UL );

      if( upper(0,2) != 2 || upper(0,3) != 3 ||
          upper(1,3) != 5 ||
          upper(2,2) != 6 ||
          upper(3,3) != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing all element with value 1 failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 0 0 2 3 )\n( 0 0 0 5 )\n( 0 0 6 0 )\n( 0 0 0 8 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major iterator-range-based erase function with predicate
   //=====================================================================================

   {
      test_ = "Row-major UpperMatrix::erase( size_t, Iterator, Iterator, Predicate )";

      // Initialization check
      UT upper( 4UL, 8UL );
      upper(0,0) = 1;
      upper(0,2) = 2;
      upper(0,3) = 3;
      upper(1,2) = 4;
      upper(1,3) = 5;
      upper(2,2) = 6;
      upper(2,3) = 7;
      upper(3,3) = 8;

      checkRows    ( upper, 4UL );
      checkColumns ( upper, 4UL );
      checkCapacity( upper, 8UL );
      checkNonZeros( upper, 8UL );
      checkNonZeros( upper, 0UL, 3UL );
      checkNonZeros( upper, 1UL, 2UL );
      checkNonZeros( upper, 2UL, 2UL );
      checkNonZeros( upper, 3UL, 1UL );

      if( upper(0,0) != 1 || upper(0,2) != 2 || upper(0,3) != 3 ||
          upper(1,2) != 4 || upper(1,3) != 5 ||
          upper(2,2) != 6 || upper(2,3) != 7 ||
          upper(3,3) != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 0 2 3 )\n( 0 0 4 5 )\n( 0 0 6 7 )\n( 0 0 0 8 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing a selection of elements
      upper.erase( 0UL, upper.begin( 0UL ), upper.find( 0UL, 3UL ),
                   []( int value ){ return value == 1 || value == 2; } );

      checkRows    ( upper, 4UL );
      checkColumns ( upper, 4UL );
      checkCapacity( upper, 8UL );
      checkNonZeros( upper, 6UL );
      checkNonZeros( upper, 0UL, 1UL );
      checkNonZeros( upper, 1UL, 2UL );
      checkNonZeros( upper, 2UL, 2UL );
      checkNonZeros( upper, 3UL, 1UL );

      if( upper(0,3) != 3 ||
          upper(1,2) != 4 || upper(1,3) != 5 ||
          upper(2,2) != 6 || upper(2,3) != 7 ||
          upper(3,3) != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a selection of elements failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 0 0 0 3 )\n( 0 0 4 5 )\n( 0 0 6 7 )\n( 0 0 0 8 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to erase from an empty range
      upper.erase( 0UL, upper.begin( 0UL ), upper.begin( 0UL ), []( int ){ return true; } );

      checkRows    ( upper, 4UL );
      checkColumns ( upper, 4UL );
      checkCapacity( upper, 8UL );
      checkNonZeros( upper, 6UL );
      checkNonZeros( upper, 0UL, 1UL );
      checkNonZeros( upper, 1UL, 2UL );
      checkNonZeros( upper, 2UL, 2UL );
      checkNonZeros( upper, 3UL, 1UL );

      if( upper(0,3) != 3 ||
          upper(1,2) != 4 || upper(1,3) != 5 ||
          upper(2,2) != 6 || upper(2,3) != 7 ||
          upper(3,3) != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing from an empty range failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 0 0 0 3 )\n( 0 0 4 5 )\n( 0 0 6 7 )\n( 0 0 0 8 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major index-based erase function
   //=====================================================================================

   {
      test_ = "Column-major UpperMatrix::erase( size_t, size_t )";

      // Initialization check
      OUT upper( 4UL, 8UL );
      upper(0,0) = 1;
      upper(0,1) = 2;
      upper(0,2) = 3;
      upper(0,3) = 4;
      upper(1,1) = 5;
      upper(1,2) = 6;
      upper(1,3) = 7;
      upper(3,3) = 8;

      checkRows    ( upper, 4UL );
      checkColumns ( upper, 4UL );
      checkCapacity( upper, 8UL );
      checkNonZeros( upper, 8UL );
      checkNonZeros( upper, 0UL, 1UL );
      checkNonZeros( upper, 1UL, 2UL );
      checkNonZeros( upper, 2UL, 2UL );
      checkNonZeros( upper, 3UL, 3UL );

      if( upper(0,0) != 1 || upper(0,1) != 2 || upper(0,2) != 3 || upper(0,3) != 4 ||
          upper(1,1) != 5 || upper(1,2) != 6 || upper(1,3) != 7 ||
          upper(3,3) != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 2 3 4 )\n( 0 5 6 7 )\n( 0 0 0 0 )\n( 0 0 0 8 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element at (0,1)
      upper.erase( 0UL, 1UL );

      checkRows    ( upper, 4UL );
      checkColumns ( upper, 4UL );
      checkCapacity( upper, 8UL );
      checkNonZeros( upper, 7UL );
      checkNonZeros( upper, 0UL, 1UL );
      checkNonZeros( upper, 1UL, 1UL );
      checkNonZeros( upper, 2UL, 2UL );
      checkNonZeros( upper, 3UL, 3UL );

      if( upper(0,0) != 1 || upper(0,2) != 3 || upper(0,3) != 4 ||
          upper(1,1) != 5 || upper(1,2) != 6 || upper(1,3) != 7 ||
          upper(3,3) != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 0 3 4 )\n( 0 5 6 7 )\n( 0 0 0 0 )\n( 0 0 0 8 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element at (1,2)
      upper.erase( 1UL, 2UL );

      checkRows    ( upper, 4UL );
      checkColumns ( upper, 4UL );
      checkCapacity( upper, 8UL );
      checkNonZeros( upper, 6UL );
      checkNonZeros( upper, 0UL, 1UL );
      checkNonZeros( upper, 1UL, 1UL );
      checkNonZeros( upper, 2UL, 1UL );
      checkNonZeros( upper, 3UL, 3UL );

      if( upper(0,0) != 1 || upper(0,2) != 3 || upper(0,3) != 4 ||
          upper(1,1) != 5 || upper(1,3) != 7 ||
          upper(3,3) != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 0 3 4 )\n( 0 5 0 7 )\n( 0 0 0 0 )\n( 0 0 0 8 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element at (1,3)
      upper.erase( 1UL, 3UL );

      checkRows    ( upper, 4UL );
      checkColumns ( upper, 4UL );
      checkCapacity( upper, 8UL );
      checkNonZeros( upper, 5UL );
      checkNonZeros( upper, 0UL, 1UL );
      checkNonZeros( upper, 1UL, 1UL );
      checkNonZeros( upper, 2UL, 1UL );
      checkNonZeros( upper, 3UL, 2UL );

      if( upper(0,0) != 1 || upper(0,2) != 3 || upper(0,3) != 4 ||
          upper(1,1) != 5 ||
          upper(3,3) != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 0 3 4 )\n( 0 5 0 0 )\n( 0 0 0 0 )\n( 0 0 0 8 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to erase a zero element
      upper.erase( 2UL, 3UL );

      checkRows    ( upper, 4UL );
      checkColumns ( upper, 4UL );
      checkCapacity( upper, 8UL );
      checkNonZeros( upper, 5UL );
      checkNonZeros( upper, 0UL, 1UL );
      checkNonZeros( upper, 1UL, 1UL );
      checkNonZeros( upper, 2UL, 1UL );
      checkNonZeros( upper, 3UL, 2UL );

      if( upper(0,0) != 1 || upper(0,2) != 3 || upper(0,3) != 4 ||
          upper(1,1) != 5 ||
          upper(3,3) != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 0 3 4 )\n( 0 5 0 0 )\n( 0 0 0 0 )\n( 0 0 0 8 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major iterator-based erase function
   //=====================================================================================

   {
      test_ = "Column-major UpperMatrix::erase( size_t, Iterator )";

      using Iterator = OUT::Iterator;

      // Initialization check
      OUT upper( 4UL, 8UL );
      upper(0,0) = 1;
      upper(0,1) = 2;
      upper(0,2) = 3;
      upper(0,3) = 4;
      upper(1,1) = 5;
      upper(1,2) = 6;
      upper(1,3) = 7;
      upper(3,3) = 8;

      checkRows    ( upper, 4UL );
      checkColumns ( upper, 4UL );
      checkCapacity( upper, 8UL );
      checkNonZeros( upper, 8UL );
      checkNonZeros( upper, 0UL, 1UL );
      checkNonZeros( upper, 1UL, 2UL );
      checkNonZeros( upper, 2UL, 2UL );
      checkNonZeros( upper, 3UL, 3UL );

      if( upper(0,0) != 1 || upper(0,1) != 2 || upper(0,2) != 3 || upper(0,3) != 4 ||
          upper(1,1) != 5 || upper(1,2) != 6 || upper(1,3) != 7 ||
          upper(3,3) != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 2 3 4 )\n( 0 5 6 7 )\n( 0 0 0 0 )\n( 0 0 0 8 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element at (0,1)
      {
         Iterator pos = upper.erase( 1UL, upper.find( 0UL, 1UL ) );

         checkRows    ( upper, 4UL );
         checkColumns ( upper, 4UL );
         checkCapacity( upper, 8UL );
         checkNonZeros( upper, 7UL );
         checkNonZeros( upper, 0UL, 1UL );
         checkNonZeros( upper, 1UL, 1UL );
         checkNonZeros( upper, 2UL, 2UL );
         checkNonZeros( upper, 3UL, 3UL );

         if( upper(0,0) != 1 || upper(0,2) != 3 || upper(0,3) != 4 ||
             upper(1,1) != 5 || upper(1,2) != 6 || upper(1,3) != 7 ||
             upper(3,3) != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 1 0 3 4 )\n( 0 5 6 7 )\n( 0 0 0 0 )\n( 0 0 0 8 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos->value() != 5 || pos->index() != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 5\n"
                << "   Expected index: 1\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the element at (1,2)
      {
         Iterator pos = upper.erase( 2UL, upper.find( 1UL, 2UL ) );

         checkRows    ( upper, 4UL );
         checkColumns ( upper, 4UL );
         checkCapacity( upper, 8UL );
         checkNonZeros( upper, 6UL );
         checkNonZeros( upper, 0UL, 1UL );
         checkNonZeros( upper, 1UL, 1UL );
         checkNonZeros( upper, 2UL, 1UL );
         checkNonZeros( upper, 3UL, 3UL );

         if( upper(0,0) != 1 || upper(0,2) != 3 || upper(0,3) != 4 ||
             upper(1,1) != 5 || upper(1,3) != 7 ||
             upper(3,3) != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 1 0 3 4 )\n( 0 5 0 7 )\n( 0 0 0 0 )\n( 0 0 0 8 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos != upper.end( 2UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the element at (1,3)
      {
         Iterator pos = upper.erase( 3UL, upper.find( 1UL, 3UL ) );

         checkRows    ( upper, 4UL );
         checkColumns ( upper, 4UL );
         checkCapacity( upper, 8UL );
         checkNonZeros( upper, 5UL );
         checkNonZeros( upper, 0UL, 1UL );
         checkNonZeros( upper, 1UL, 1UL );
         checkNonZeros( upper, 2UL, 1UL );
         checkNonZeros( upper, 3UL, 2UL );

         if( upper(0,0) != 1 || upper(0,2) != 3 || upper(0,3) != 4 ||
             upper(1,1) != 5 ||
             upper(3,3) != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 1 0 3 4 )\n( 0 5 0 0 )\n( 0 0 0 0 )\n( 0 0 0 8 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos->value() != 8 || pos->index() != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 8\n"
                << "   Expected index: 3\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to erase a zero element
      {
         Iterator pos = upper.erase( 3UL, upper.find( 2UL, 3UL ) );

         checkRows    ( upper, 4UL );
         checkColumns ( upper, 4UL );
         checkCapacity( upper, 8UL );
         checkNonZeros( upper, 5UL );
         checkNonZeros( upper, 0UL, 1UL );
         checkNonZeros( upper, 1UL, 1UL );
         checkNonZeros( upper, 2UL, 1UL );
         checkNonZeros( upper, 3UL, 2UL );

         if( upper(0,0) != 1 || upper(0,2) != 3 || upper(0,3) != 4 ||
             upper(1,1) != 5 ||
             upper(3,3) != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 1 0 3 4 )\n( 0 5 0 0 )\n( 0 0 0 0 )\n( 0 0 0 8 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos != upper.end( 3UL ) ) {
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
      test_ = "Column-major UpperMatrix::erase( size_t, Iterator, Iterator )";

      using Iterator = OUT::Iterator;

      // Initialization check
      OUT upper( 4UL, 8UL );
      upper(0,0) = 1;
      upper(0,1) = 2;
      upper(0,2) = 3;
      upper(0,3) = 4;
      upper(1,1) = 5;
      upper(1,2) = 6;
      upper(1,3) = 7;
      upper(3,3) = 8;

      checkRows    ( upper, 4UL );
      checkColumns ( upper, 4UL );
      checkCapacity( upper, 8UL );
      checkNonZeros( upper, 8UL );
      checkNonZeros( upper, 0UL, 1UL );
      checkNonZeros( upper, 1UL, 2UL );
      checkNonZeros( upper, 2UL, 2UL );
      checkNonZeros( upper, 3UL, 3UL );

      if( upper(0,0) != 1 || upper(0,1) != 2 || upper(0,2) != 3 || upper(0,3) != 4 ||
          upper(1,1) != 5 || upper(1,2) != 6 || upper(1,3) != 7 ||
          upper(3,3) != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 2 3 4 )\n( 0 5 6 7 )\n( 0 0 0 0 )\n( 0 0 0 8 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the elements from the beginning of column 1 to (1,1)
      {
         Iterator pos = upper.erase( 1UL, upper.begin( 1UL ), upper.find( 1UL, 1UL ) );

         checkRows    ( upper, 4UL );
         checkColumns ( upper, 4UL );
         checkCapacity( upper, 8UL );
         checkNonZeros( upper, 7UL );
         checkNonZeros( upper, 0UL, 1UL );
         checkNonZeros( upper, 1UL, 1UL );
         checkNonZeros( upper, 2UL, 2UL );
         checkNonZeros( upper, 3UL, 3UL );

         if( upper(0,0) != 1 || upper(0,2) != 3 || upper(0,3) != 4 ||
             upper(1,1) != 5 || upper(1,2) != 6 || upper(1,3) != 7 ||
             upper(3,3) != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a single-element range failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 1 0 3 4 )\n( 0 5 6 7 )\n( 0 0 0 0 )\n( 0 0 0 8 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos->value() != 5 || pos->index() != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 5\n"
                << "   Expected index: 1\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the elements from (1,2) to the column end
      {
         Iterator pos = upper.erase( 2UL, upper.find( 1UL, 2UL ), upper.end( 2UL ) );

         checkRows    ( upper, 4UL );
         checkColumns ( upper, 4UL );
         checkCapacity( upper, 8UL );
         checkNonZeros( upper, 6UL );
         checkNonZeros( upper, 0UL, 1UL );
         checkNonZeros( upper, 1UL, 1UL );
         checkNonZeros( upper, 2UL, 1UL );
         checkNonZeros( upper, 3UL, 3UL );

         if( upper(0,0) != 1 || upper(0,2) != 3 || upper(0,3) != 4 ||
             upper(1,1) != 5 || upper(1,3) != 7 ||
             upper(3,3) != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a single-element range failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 1 0 3 4 )\n( 0 5 0 7 )\n( 0 0 0 0 )\n( 0 0 0 8 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos != upper.end( 2UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the elements from (0,3) to (3,3)
      {
         Iterator pos = upper.erase( 3UL, upper.find( 0UL, 3UL ), upper.find( 3UL, 3UL ) );

         checkRows    ( upper, 4UL );
         checkColumns ( upper, 4UL );
         checkCapacity( upper, 8UL );
         checkNonZeros( upper, 4UL );
         checkNonZeros( upper, 0UL, 1UL );
         checkNonZeros( upper, 1UL, 1UL );
         checkNonZeros( upper, 2UL, 1UL );
         checkNonZeros( upper, 3UL, 1UL );

         if( upper(0,0) != 1 || upper(0,2) != 3 ||
             upper(1,1) != 5 ||
             upper(3,3) != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a multi-element range failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 1 0 3 0 )\n( 0 5 0 0 )\n( 0 0 0 0 )\n( 0 0 0 8 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos->value() != 8 || pos->index() != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 8\n"
                << "   Expected index: 3\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to erase an empty range
      {
         Iterator pos = upper.erase( 3UL, upper.find( 3UL, 3UL ), upper.find( 3UL, 3UL ) );

         checkRows    ( upper, 4UL );
         checkColumns ( upper, 4UL );
         checkCapacity( upper, 8UL );
         checkNonZeros( upper, 4UL );
         checkNonZeros( upper, 0UL, 1UL );
         checkNonZeros( upper, 1UL, 1UL );
         checkNonZeros( upper, 2UL, 1UL );
         checkNonZeros( upper, 3UL, 1UL );

         if( upper(0,0) != 1 || upper(0,2) != 3 ||
             upper(1,1) != 5 ||
             upper(3,3) != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing an empty range failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 1 0 3 0 )\n( 0 5 0 0 )\n( 0 0 0 0 )\n( 0 0 0 8 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos->value() != 8 || pos->index() != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 8\n"
                << "   Expected index: 3\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major erase function with predicate
   //=====================================================================================

   {
      test_ = "Column-major UpperMatrix::erase( Predicate )";

      // Initialization check
      OUT upper( 4UL, 8UL );
      upper(0,0) = 1;
      upper(0,2) = 2;
      upper(0,3) = 3;
      upper(1,2) = 4;
      upper(1,3) = 5;
      upper(2,2) = 6;
      upper(2,3) = 7;
      upper(3,3) = 8;

      checkRows    ( upper, 4UL );
      checkColumns ( upper, 4UL );
      checkCapacity( upper, 8UL );
      checkNonZeros( upper, 8UL );
      checkNonZeros( upper, 0UL, 1UL );
      checkNonZeros( upper, 1UL, 0UL );
      checkNonZeros( upper, 2UL, 3UL );
      checkNonZeros( upper, 3UL, 4UL );

      if( upper(0,0) != 1 || upper(0,2) != 2 || upper(0,3) != 3 ||
          upper(1,2) != 4 || upper(1,3) != 5 ||
          upper(2,2) != 6 || upper(2,3) != 7 ||
          upper(3,3) != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 0 2 3 )\n( 0 0 4 5 )\n( 0 0 6 7 )\n( 0 0 0 8 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing a selection of elements
      upper.erase( []( int value ){ return value == 1 || value == 4 || value == 7; } );

      checkRows    ( upper, 4UL );
      checkColumns ( upper, 4UL );
      checkCapacity( upper, 8UL );
      checkNonZeros( upper, 5UL );
      checkNonZeros( upper, 0UL, 0UL );
      checkNonZeros( upper, 1UL, 0UL );
      checkNonZeros( upper, 2UL, 2UL );
      checkNonZeros( upper, 3UL, 3UL );

      if( upper(0,2) != 2 || upper(0,3) != 3 ||
          upper(1,3) != 5 ||
          upper(2,2) != 6 ||
          upper(3,3) != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a selection of elements failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 0 0 2 3 )\n( 0 0 0 5 )\n( 0 0 6 0 )\n( 0 0 0 8 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to erase all elements with value 1
      upper.erase( []( int value ){ return value == 1; } );

      checkRows    ( upper, 4UL );
      checkColumns ( upper, 4UL );
      checkCapacity( upper, 8UL );
      checkNonZeros( upper, 5UL );
      checkNonZeros( upper, 0UL, 0UL );
      checkNonZeros( upper, 1UL, 0UL );
      checkNonZeros( upper, 2UL, 2UL );
      checkNonZeros( upper, 3UL, 3UL );

      if( upper(0,2) != 2 || upper(0,3) != 3 ||
          upper(1,3) != 5 ||
          upper(2,2) != 6 ||
          upper(3,3) != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing all element with value 1 failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 0 0 2 3 )\n( 0 0 0 5 )\n( 0 0 6 0 )\n( 0 0 0 8 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major iterator-range-based erase function with predicate
   //=====================================================================================

   {
      test_ = "Column-major UpperMatrix::erase( size_t, Iterator, Iterator, Predicate )";

      // Initialization check
      OUT upper( 4UL, 8UL );
      upper(0,0) = 1;
      upper(0,2) = 2;
      upper(0,3) = 3;
      upper(1,2) = 4;
      upper(1,3) = 5;
      upper(2,2) = 6;
      upper(2,3) = 7;
      upper(3,3) = 8;

      checkRows    ( upper, 4UL );
      checkColumns ( upper, 4UL );
      checkCapacity( upper, 8UL );
      checkNonZeros( upper, 8UL );
      checkNonZeros( upper, 0UL, 1UL );
      checkNonZeros( upper, 1UL, 0UL );
      checkNonZeros( upper, 2UL, 3UL );
      checkNonZeros( upper, 3UL, 4UL );

      if( upper(0,0) != 1 || upper(0,2) != 2 || upper(0,3) != 3 ||
          upper(1,2) != 4 || upper(1,3) != 5 ||
          upper(2,2) != 6 || upper(2,3) != 7 ||
          upper(3,3) != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 0 2 3 )\n( 0 0 4 5 )\n( 0 0 6 7 )\n( 0 0 0 8 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing a selection of elements
      upper.erase( 3UL, upper.begin( 3UL ), upper.find( 3UL, 3UL ),
                   []( int value ){ return value == 3 || value == 7; } );

      checkRows    ( upper, 4UL );
      checkColumns ( upper, 4UL );
      checkCapacity( upper, 8UL );
      checkNonZeros( upper, 6UL );
      checkNonZeros( upper, 0UL, 1UL );
      checkNonZeros( upper, 1UL, 0UL );
      checkNonZeros( upper, 2UL, 3UL );
      checkNonZeros( upper, 3UL, 2UL );

      if( upper(0,0) != 1 || upper(0,2) != 2 ||
          upper(1,2) != 4 || upper(1,3) != 5 ||
          upper(2,2) != 6 ||
          upper(3,3) != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a selection of elements failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 0 2 0 )\n( 0 0 4 5 )\n( 0 0 6 0 )\n( 0 0 0 8 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to erase from an empty range
      upper.erase( 3UL, upper.begin( 3UL ), upper.begin( 3UL ), []( int ){ return true; } );

      checkRows    ( upper, 4UL );
      checkColumns ( upper, 4UL );
      checkCapacity( upper, 8UL );
      checkNonZeros( upper, 6UL );
      checkNonZeros( upper, 0UL, 1UL );
      checkNonZeros( upper, 1UL, 0UL );
      checkNonZeros( upper, 2UL, 3UL );
      checkNonZeros( upper, 3UL, 2UL );

      if( upper(0,0) != 1 || upper(0,2) != 2 ||
          upper(1,2) != 4 || upper(1,3) != 5 ||
          upper(2,2) != 6 ||
          upper(3,3) != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing from an empty range failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 0 2 0 )\n( 0 0 4 5 )\n( 0 0 6 0 )\n( 0 0 0 8 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c find() member function of the UpperMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c find() member function of the UpperMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testFind()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major UpperMatrix::find()";

      using ConstIterator = UT::ConstIterator;

      // Initialization check
      UT upper( 8UL, 3UL );
      upper(1,2) = 1;
      upper(2,3) = 2;
      upper(5,6) = 3;

      checkRows    ( upper, 8UL );
      checkColumns ( upper, 8UL );
      checkCapacity( upper, 3UL );
      checkNonZeros( upper, 3UL );
      checkNonZeros( upper, 0UL, 0UL );
      checkNonZeros( upper, 1UL, 1UL );
      checkNonZeros( upper, 2UL, 1UL );
      checkNonZeros( upper, 3UL, 0UL );
      checkNonZeros( upper, 4UL, 0UL );
      checkNonZeros( upper, 5UL, 1UL );
      checkNonZeros( upper, 6UL, 0UL );
      checkNonZeros( upper, 7UL, 0UL );

      // Searching for the first element
      {
         ConstIterator pos( upper.find( 1UL, 2UL ) );

         if( pos == upper.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required position = (1,2)\n"
                << "   Current matrix:\n" << upper << "\n";
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
                << "   Current matrix:\n" << upper << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for the second element
      {
         ConstIterator pos( upper.find( 2UL, 3UL ) );

         if( pos == upper.end( 2UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required position = (2,3)\n"
                << "   Current matrix:\n" << upper << "\n";
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
                << "   Current matrix:\n" << upper << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for the third element
      {
         ConstIterator pos( upper.find( 5UL, 6UL ) );

         if( pos == upper.end( 5UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required position = (5,6)\n"
                << "   Current matrix:\n" << upper << "\n";
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
                << "   Current matrix:\n" << upper << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for a non-existing non-zero element
      {
         ConstIterator pos( upper.find( 0UL, 4UL ) );

         if( pos != upper.end( 0UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Non-existing element could be found\n"
                << " Details:\n"
                << "   Required index = 4\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 0\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << upper << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major UpperMatrix::find()";

      using ConstIterator = OUT::ConstIterator;

      // Initialization check
      OUT upper( 8UL, 3UL );
      upper(1,2) = 1;
      upper(2,3) = 2;
      upper(5,6) = 3;

      checkRows    ( upper, 8UL );
      checkColumns ( upper, 8UL );
      checkCapacity( upper, 3UL );
      checkNonZeros( upper, 3UL );
      checkNonZeros( upper, 0UL, 0UL );
      checkNonZeros( upper, 1UL, 0UL );
      checkNonZeros( upper, 2UL, 1UL );
      checkNonZeros( upper, 3UL, 1UL );
      checkNonZeros( upper, 4UL, 0UL );
      checkNonZeros( upper, 5UL, 0UL );
      checkNonZeros( upper, 6UL, 1UL );
      checkNonZeros( upper, 7UL, 0UL );

      // Searching for the first element
      {
         ConstIterator pos( upper.find( 1UL, 2UL ) );

         if( pos == upper.end( 2UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required position = (1,2)\n"
                << "   Current matrix:\n" << upper << "\n";
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
                << "   Current matrix:\n" << upper << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for the second element
      {
         ConstIterator pos( upper.find( 2UL, 3UL ) );

         if( pos == upper.end( 3UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required position = (2,3)\n"
                << "   Current matrix:\n" << upper << "\n";
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
                << "   Current matrix:\n" << upper << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for the third element
      {
         ConstIterator pos( upper.find( 5UL, 6UL ) );

         if( pos == upper.end( 6UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required position = (5,6)\n"
                << "   Current matrix:\n" << upper << "\n";
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
                << "   Current matrix:\n" << upper << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for a non-existing non-zero element
      {
         ConstIterator pos( upper.find( 0UL, 4UL ) );

         if( pos != upper.end( 4UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Non-existing element could be found\n"
                << " Details:\n"
                << "   Required index = 0\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 0\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << upper << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c lowerBound() member function of the UpperMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c lowerBound() member function of the UpperMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testLowerBound()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major UpperMatrix::lowerBound()";

      using ConstIterator = UT::ConstIterator;

      // Initialization check
      UT upper( 6UL, 2UL );
      upper(1,2) = 1;
      upper(1,4) = 2;

      checkRows    ( upper, 6UL );
      checkColumns ( upper, 6UL );
      checkCapacity( upper, 2UL );
      checkNonZeros( upper, 2UL );
      checkNonZeros( upper, 0UL, 0UL );
      checkNonZeros( upper, 1UL, 2UL );
      checkNonZeros( upper, 2UL, 0UL );
      checkNonZeros( upper, 3UL, 0UL );
      checkNonZeros( upper, 4UL, 0UL );
      checkNonZeros( upper, 5UL, 0UL );

      // Determining the lower bound for position (1,1)
      {
         ConstIterator pos( upper.lowerBound( 1UL, 1UL ) );

         if( pos == upper.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,1)\n"
                << "   Current matrix:\n" << upper << "\n";
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
                << "   Current matrix:\n" << upper << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for position (1,2)
      {
         ConstIterator pos( upper.lowerBound( 1UL, 2UL ) );

         if( pos == upper.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,2)\n"
                << "   Current matrix:\n" << upper << "\n";
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
                << "   Current matrix:\n" << upper << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for position (1,3)
      {
         ConstIterator pos( upper.lowerBound( 1UL, 3UL ) );

         if( pos == upper.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,3)\n"
                << "   Current matrix:\n" << upper << "\n";
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
                << "   Current matrix:\n" << upper << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for position (1,4)
      {
         ConstIterator pos( upper.lowerBound( 1UL, 4UL ) );

         if( pos == upper.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,4)\n"
                << "   Current matrix:\n" << upper << "\n";
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
                << "   Current matrix:\n" << upper << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for position (1,5)
      {
         ConstIterator pos( upper.lowerBound( 1UL, 5UL ) );

         if( pos != upper.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,5)\n"
                << "   Current matrix:\n" << upper << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major UpperMatrix::lowerBound()";

      using ConstIterator = OUT::ConstIterator;

      // Initialization check
      OUT upper( 6UL, 2UL );
      upper(1,4) = 1;
      upper(3,4) = 2;

      checkRows    ( upper, 6UL );
      checkColumns ( upper, 6UL );
      checkCapacity( upper, 2UL );
      checkNonZeros( upper, 2UL );
      checkNonZeros( upper, 0UL, 0UL );
      checkNonZeros( upper, 1UL, 0UL );
      checkNonZeros( upper, 2UL, 0UL );
      checkNonZeros( upper, 3UL, 0UL );
      checkNonZeros( upper, 4UL, 2UL );
      checkNonZeros( upper, 5UL, 0UL );

      // Determining the lower bound for position (0,4)
      {
         ConstIterator pos( upper.lowerBound( 0UL, 4UL ) );

         if( pos == upper.end( 4UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (0,4)\n"
                << "   Current matrix:\n" << upper << "\n";
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
                << "   Current matrix:\n" << upper << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for position (1,4)
      {
         ConstIterator pos( upper.lowerBound( 1UL, 4UL ) );

         if( pos == upper.end( 4UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,4)\n"
                << "   Current matrix:\n" << upper << "\n";
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
                << "   Current matrix:\n" << upper << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for position (2,4)
      {
         ConstIterator pos( upper.lowerBound( 2UL, 4UL ) );

         if( pos == upper.end( 4UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (2,4)\n"
                << "   Current matrix:\n" << upper << "\n";
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
                << "   Current matrix:\n" << upper << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for position (3,4)
      {
         ConstIterator pos( upper.lowerBound( 3UL, 4UL ) );

         if( pos == upper.end( 4UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (3,4)\n"
                << "   Current matrix:\n" << upper << "\n";
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
                << "   Current matrix:\n" << upper << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for position (4,4)
      {
         ConstIterator pos( upper.lowerBound( 4UL, 4UL ) );

         if( pos != upper.end( 4UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (4,4)\n"
                << "   Current matrix:\n" << upper << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c upperBound() member function of the UpperMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c upperBound() member function of the UpperMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testUpperBound()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major UpperMatrix::upperBound()";

      using ConstIterator = UT::ConstIterator;

      // Initialization check
      UT upper( 6UL, 2UL );
      upper(1,2) = 1;
      upper(1,4) = 2;

      checkRows    ( upper, 6UL );
      checkColumns ( upper, 6UL );
      checkCapacity( upper, 2UL );
      checkNonZeros( upper, 2UL );
      checkNonZeros( upper, 0UL, 0UL );
      checkNonZeros( upper, 1UL, 2UL );
      checkNonZeros( upper, 2UL, 0UL );
      checkNonZeros( upper, 3UL, 0UL );
      checkNonZeros( upper, 4UL, 0UL );
      checkNonZeros( upper, 5UL, 0UL );

      // Determining the upper bound for position (1,1)
      {
         ConstIterator pos( upper.upperBound( 1UL, 1UL ) );

         if( pos == upper.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,1)\n"
                << "   Current matrix:\n" << upper << "\n";
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
                << "   Current matrix:\n" << upper << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for position (1,2)
      {
         ConstIterator pos( upper.upperBound( 1UL, 2UL ) );

         if( pos == upper.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,2)\n"
                << "   Current matrix:\n" << upper << "\n";
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
                << "   Current matrix:\n" << upper << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for position (1,3)
      {
         ConstIterator pos( upper.upperBound( 1UL, 3UL ) );

         if( pos == upper.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,3)\n"
                << "   Current matrix:\n" << upper << "\n";
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
                << "   Current matrix:\n" << upper << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for position (1,4)
      {
         ConstIterator pos( upper.upperBound( 1UL, 4UL ) );

         if( pos != upper.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,4)\n"
                << "   Current matrix:\n" << upper << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for position (1,5)
      {
         ConstIterator pos( upper.upperBound( 1UL, 5UL ) );

         if( pos != upper.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,5)\n"
                << "   Current matrix:\n" << upper << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major UpperMatrix::upperBound()";

      using ConstIterator = OUT::ConstIterator;

      // Initialization check
      OUT upper( 6UL, 2UL );
      upper(1,4) = 1;
      upper(3,4) = 2;

      checkRows    ( upper, 6UL );
      checkColumns ( upper, 6UL );
      checkCapacity( upper, 2UL );
      checkNonZeros( upper, 2UL );
      checkNonZeros( upper, 0UL, 0UL );
      checkNonZeros( upper, 1UL, 0UL );
      checkNonZeros( upper, 2UL, 0UL );
      checkNonZeros( upper, 3UL, 0UL );
      checkNonZeros( upper, 4UL, 2UL );
      checkNonZeros( upper, 5UL, 0UL );

      // Determining the upper bound for position (0,4)
      {
         ConstIterator pos( upper.upperBound( 0UL, 4UL ) );

         if( pos == upper.end( 4UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (0,4)\n"
                << "   Current matrix:\n" << upper << "\n";
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
                << "   Current matrix:\n" << upper << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for position (1,4)
      {
         ConstIterator pos( upper.upperBound( 1UL, 4UL ) );

         if( pos == upper.end( 4UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,4)\n"
                << "   Current matrix:\n" << upper << "\n";
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
                << "   Current matrix:\n" << upper << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for position (2,4)
      {
         ConstIterator pos( upper.upperBound( 2UL, 4UL ) );

         if( pos == upper.end( 4UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (2,4)\n"
                << "   Current matrix:\n" << upper << "\n";
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
                << "   Current matrix:\n" << upper << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for position (3,4)
      {
         ConstIterator pos( upper.upperBound( 3UL, 4UL ) );

         if( pos != upper.end( 4UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (3,4)\n"
                << "   Current matrix:\n" << upper << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for position (4,4)
      {
         ConstIterator pos( upper.upperBound( 4UL, 4UL ) );

         if( pos != upper.end( 4UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (4,4)\n"
                << "   Current matrix:\n" << upper << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c isDefault() function with the UpperMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isDefault() function with the UpperMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testIsDefault()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major isDefault() function";

      // isDefault with 0x0 matrix
      {
         UT upper;

         if( isDefault( upper ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << upper << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with default matrix
      {
         UT upper( 3UL );

         if( isDefault( upper(0,1) ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix element:\n" << upper(0,1) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( upper ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << upper << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with non-default matrix
      {
         UT upper( 3UL );
         upper(0,1) = 1;

         if( isDefault( upper(0,1) ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix element:\n" << upper(0,1) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( upper ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << upper << "\n";
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
         OUT upper;

         if( isDefault( upper ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << upper << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with default matrix
      {
         OUT upper( 3UL );

         if( isDefault( upper(0,1) ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix element:\n" << upper(0,1) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( upper ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << upper << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with non-default matrix
      {
         OUT upper( 3UL );
         upper(0,1) = 1;

         if( isDefault( upper(0,1) ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix element:\n" << upper(0,1) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( upper ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << upper << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c submatrix() function with the UpperMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c submatrix() function with the UpperMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testSubmatrix()
{
   //=====================================================================================
   // Row-major general tests
   //=====================================================================================

   {
      test_ = "Row-major submatrix() function";

      using SMT = blaze::Submatrix<UT>;

      UT upper( 3UL );
      upper(0,0) =  1;
      upper(0,1) = -4;
      upper(0,2) =  7;
      upper(1,1) =  2;
      upper(2,2) =  3;

      SMT sm = submatrix( upper, 1UL, 1UL, 2UL, 2UL );

      if( sm(1,1) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator access failed\n"
             << " Details:\n"
             << "   Result: " << sm(1,1) << "\n"
             << "   Expected result: 3\n";
         throw std::runtime_error( oss.str() );
      }

      SMT::Iterator it = sm.begin(0UL);

      if( it == sm.end(0UL) || it->value() != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << it->value() << "\n"
             << "   Expected result: 2\n";
         throw std::runtime_error( oss.str() );
      }

      sm(0,1) = -5;

      if( sm(0,0) != 2 || sm(0,1) != -5 ||
          sm(1,0) != 0 || sm(1,1) !=  3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix access failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 2 -5 )\n( 0  3 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( upper(0,0) != 1 || upper(0,1) != -4 || upper(0,2) !=  7 ||
          upper(1,0) != 0 || upper(1,1) !=  2 || upper(1,2) != -5 ||
          upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) !=  3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix access failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 -4  7 )\n( 0  2 -5 )\n( 0  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }

      reset( sm );

      if( sm(0,0) != 0 || sm(0,1) != 0 ||
          sm(1,0) != 0 || sm(1,1) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix reset failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0 0 )\n( 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( upper(0,0) != 1 || upper(0,1) != -4 || upper(0,2) != 7 ||
          upper(1,0) != 0 || upper(1,1) !=  0 || upper(1,2) != 0 ||
          upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix reset failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 -4  7 )\n( 0  0  0 )\n( 0  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major general tests
   //=====================================================================================

   {
      test_ = "Column-major submatrix() function";

      using SMT = blaze::Submatrix<OUT>;

      OUT upper( 3UL );
      upper(0,0) =  1;
      upper(0,1) = -4;
      upper(0,2) =  7;
      upper(1,1) =  2;
      upper(2,2) =  3;

      SMT sm = submatrix( upper, 1UL, 1UL, 2UL, 2UL );

      if( sm(1,1) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator access failed\n"
             << " Details:\n"
             << "   Result: " << sm(1,1) << "\n"
             << "   Expected result: 3\n";
         throw std::runtime_error( oss.str() );
      }

      SMT::Iterator it = sm.begin(0UL);

      if( it == sm.end(0UL) || it->value() != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << it->value() << "\n"
             << "   Expected result: 2\n";
         throw std::runtime_error( oss.str() );
      }

      sm(0,1) = -5;

      if( sm(0,0) != 2 || sm(0,1) != -5 ||
          sm(1,0) != 0 || sm(1,1) !=  3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix access failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 2 -5 )\n( 0  3 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( upper(0,0) != 1 || upper(0,1) != -4 || upper(0,2) !=  7 ||
          upper(1,0) != 0 || upper(1,1) !=  2 || upper(1,2) != -5 ||
          upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) !=  3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix access failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 -4  7 )\n( 0  2 -5 )\n( 0  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }

      reset( sm );

      if( sm(0,0) != 0 || sm(0,1) != 0 ||
          sm(1,0) != 0 || sm(1,1) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix reset failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0 0 )\n( 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( upper(0,0) != 1 || upper(0,1) != -4 || upper(0,2) != 7 ||
          upper(1,0) != 0 || upper(1,1) !=  0 || upper(1,2) != 0 ||
          upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix reset failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 -4  7 )\n( 0  0  0 )\n( 0  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c row() function with the UpperMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c row() function with the UpperMatrix specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testRow()
{
   //=====================================================================================
   // Row-major general tests
   //=====================================================================================

   {
      test_ = "Row-major row() function";

      using RT = blaze::Row<UT>;

      UT upper( 3UL );
      upper(0,0) =  1;
      upper(0,1) = -4;
      upper(0,2) =  7;
      upper(1,1) =  2;
      upper(2,2) =  3;

      RT row1 = row( upper, 1UL );

      if( row1[1] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator access failed\n"
             << " Details:\n"
             << "   Result: " << row1[1] << "\n"
             << "   Expected result: 2\n";
         throw std::runtime_error( oss.str() );
      }

      RT::Iterator it( row1.begin() );

      if( it == row1.end() || it->value() != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << it->value() << "\n"
             << "   Expected result: 2\n";
         throw std::runtime_error( oss.str() );
      }

      row1[1] = -5;

      if( row1[0] != 0 || row1[1] != -5 || row1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row access failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( 0 -5  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( upper(0,0) != 1 || upper(0,1) != -4 || upper(0,2) != 7 ||
          upper(1,0) != 0 || upper(1,1) != -5 || upper(1,2) != 0 ||
          upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row access failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 -4  7 )\n( 0 -5  0 )\n( 0  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }

      reset( row1 );

      if( row1[0] != 0 || row1[1] != 0 || row1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row reset failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( upper(0,0) != 1 || upper(0,1) != -4 || upper(0,2) != 7 ||
          upper(1,0) != 0 || upper(1,1) !=  0 || upper(1,2) != 0 ||
          upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row reset failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 -4  7 )\n( 0  0  0 )\n( 0  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major general tests
   //=====================================================================================

   {
      test_ = "Column-major row() function";

      using RT = blaze::Row<OUT>;

      OUT upper( 3UL );
      upper(0,0) =  1;
      upper(0,1) = -4;
      upper(0,2) =  7;
      upper(1,1) =  2;
      upper(2,2) =  3;

      RT row1 = row( upper, 1UL );

      if( row1[1] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator access failed\n"
             << " Details:\n"
             << "   Result: " << row1[1] << "\n"
             << "   Expected result: 2\n";
         throw std::runtime_error( oss.str() );
      }

      RT::Iterator it( row1.begin() );

      if( it == row1.end() || it->value() != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << it->value() << "\n"
             << "   Expected result: 2\n";
         throw std::runtime_error( oss.str() );
      }

      row1[1] = -5;

      if( row1[0] != 0 || row1[1] != -5 || row1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row access failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( 0 -5  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( upper(0,0) != 1 || upper(0,1) != -4 || upper(0,2) != 7 ||
          upper(1,0) != 0 || upper(1,1) != -5 || upper(1,2) != 0 ||
          upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row access failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 -4  7 )\n( 0 -5  0 )\n( 0  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }

      reset( row1 );

      if( row1[0] != 0 || row1[1] != 0 || row1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row reset failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( upper(0,0) != 1 || upper(0,1) != -4 || upper(0,2) != 7 ||
          upper(1,0) != 0 || upper(1,1) !=  0 || upper(1,2) != 0 ||
          upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row reset failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 -4  7 )\n( 0  0  0 )\n( 0  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c column() function with the UpperMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c column() function with the UpperMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testColumn()
{
   //=====================================================================================
   // Row-major general tests
   //=====================================================================================

   {
      test_ = "Row-major column() function";

      using CT = blaze::Column<UT>;

      UT upper( 3UL );
      upper(0,0) =  1;
      upper(0,1) = -4;
      upper(0,2) =  7;
      upper(1,1) =  2;
      upper(2,2) =  3;

      CT col1 = column( upper, 1UL );

      if( col1[1] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator access failed\n"
             << " Details:\n"
             << "   Result: " << col1[1] << "\n"
             << "   Expected result: 2\n";
         throw std::runtime_error( oss.str() );
      }

      CT::Iterator it( col1.begin() );

      if( it == col1.end() || it->value() != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << it->value() << "\n"
             << "   Expected result: -4\n";
         throw std::runtime_error( oss.str() );
      }

      col1[1] = -5;

      if( col1[0] != -4 || col1[1] != -5 || col1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column access failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( -4 -5  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( upper(0,0) != 1 || upper(0,1) != -4 || upper(0,2) != 7 ||
          upper(1,0) != 0 || upper(1,1) != -5 || upper(1,2) != 0 ||
          upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column access failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 -4  7 )\n( 0 -5  0 )\n( 0  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }

      reset( col1 );

      if( col1[0] != 0 || col1[1] != 0 || col1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column reset failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( upper(0,0) != 1 || upper(0,1) != 0 || upper(0,2) != 7 ||
          upper(1,0) != 0 || upper(1,1) != 0 || upper(1,2) != 0 ||
          upper(2,0) != 0 || upper(2,1) != 0 || upper(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column reset failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 0 7 )\n( 0 0 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major general tests
   //=====================================================================================

   {
      test_ = "Column-major column() function";

      using CT = blaze::Column<OUT>;

      OUT upper( 3UL );
      upper(0,0) =  1;
      upper(0,1) = -4;
      upper(0,2) =  7;
      upper(1,1) =  2;
      upper(2,2) =  3;

      CT col1 = column( upper, 1UL );

      if( col1[1] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator access failed\n"
             << " Details:\n"
             << "   Result: " << col1[1] << "\n"
             << "   Expected result: 2\n";
         throw std::runtime_error( oss.str() );
      }

      CT::Iterator it( col1.begin() );

      if( it == col1.end() || it->value() != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << it->value() << "\n"
             << "   Expected result: -4\n";
         throw std::runtime_error( oss.str() );
      }

      col1[1] = -5;

      if( col1[0] != -4 || col1[1] != -5 || col1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column access failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( -4 -5  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( upper(0,0) != 1 || upper(0,1) != -4 || upper(0,2) != 7 ||
          upper(1,0) != 0 || upper(1,1) != -5 || upper(1,2) != 0 ||
          upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column access failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1  -4  7 )\n( 0 -5  0 )\n( 0  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }

      reset( col1 );

      if( col1[0] != 0 || col1[1] != 0 || col1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column reset failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( upper(0,0) != 1 || upper(0,1) != 0 || upper(0,2) != 7 ||
          upper(1,0) != 0 || upper(1,1) != 0 || upper(1,2) != 0 ||
          upper(2,0) != 0 || upper(2,1) != 0 || upper(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column reset failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 0 7 )\n( 0 0 0 )\n(  0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************

} // namespace uppermatrix


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
   std::cout << "   Running UpperMatrix sparse test (part 2)..." << std::endl;

   try
   {
      RUN_UPPERMATRIX_SPARSE_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during UpperMatrix sparse test (part 2):\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
