//=================================================================================================
/*!
//  \file src/mathtest/adaptors/uppermatrix/DenseTest2.cpp
//  \brief Source file for the UpperMatrix dense test (part 2)
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
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/Row.h>
#include <blaze/math/Submatrix.h>
#include <blaze/util/Complex.h>
#include <blazetest/mathtest/adaptors/uppermatrix/DenseTest.h>

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
/*!\brief Constructor for the UpperMatrix dense test.
//
// \exception std::runtime_error Operation error detected.
*/
DenseTest::DenseTest()
{
   testScaling();
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
/*!\brief Test of all UpperMatrix (self-)scaling operations.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all available ways to scale an instance of the UpperMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testScaling()
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
      checkCapacity( upper, 9UL );
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
      checkCapacity( upper, 9UL );
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
      checkCapacity( upper, 9UL );
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
      checkCapacity( upper, 9UL );
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
      checkCapacity( upper, 9UL );
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
      checkCapacity( upper, 9UL );
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
      checkCapacity( upper, 9UL );
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
      checkCapacity( upper, 9UL );
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

      blaze::UpperMatrix< blaze::DynamicMatrix<complex<float>,blaze::rowMajor> > upper( 2UL );
      upper(0,0) = complex<float>( 1.0F, 0.0F );
      upper(0,1) = complex<float>( 2.0F, 0.0F );
      upper(1,1) = complex<float>( 4.0F, 0.0F );

      upper.scale( complex<float>( 3.0F, 0.0F ) );

      checkRows    ( upper, 2UL );
      checkColumns ( upper, 2UL );
      checkCapacity( upper, 4UL );
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
      checkCapacity( upper, 9UL );
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
      checkCapacity( upper, 9UL );
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
      checkCapacity( upper, 9UL );
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
      checkCapacity( upper, 9UL );
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
      checkCapacity( upper, 9UL );
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
      checkCapacity( upper, 9UL );
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
      checkCapacity( upper, 9UL );
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
      checkCapacity( upper, 9UL );
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

      blaze::UpperMatrix< blaze::DynamicMatrix<complex<float>,blaze::columnMajor> > upper( 2UL );
      upper(0,0) = complex<float>( 1.0F, 0.0F );
      upper(0,1) = complex<float>( 2.0F, 0.0F );
      upper(1,1) = complex<float>( 4.0F, 0.0F );

      upper.scale( complex<float>( 3.0F, 0.0F ) );

      checkRows    ( upper, 2UL );
      checkColumns ( upper, 2UL );
      checkCapacity( upper, 4UL );
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
void DenseTest::testFunctionCall()
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
         checkCapacity( upper, 9UL );
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
         checkCapacity( upper, 9UL );
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

         // Writing the upper element (0,1)
         upper(0,1) = upper(1,2);

         checkRows    ( upper, 3UL );
         checkColumns ( upper, 3UL );
         checkCapacity( upper, 9UL );
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
         checkCapacity( upper, 9UL );
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
         checkCapacity( upper, 9UL );
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
         checkCapacity( upper, 9UL );
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
         checkCapacity( upper, 9UL );
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
         checkCapacity( upper, 9UL );
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
         checkCapacity( upper, 9UL );
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
         checkCapacity( upper, 9UL );
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
         checkCapacity( upper, 9UL );
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
         checkCapacity( upper, 9UL );
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
         checkCapacity( upper, 9UL );
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
         checkCapacity( upper, 9UL );
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
void DenseTest::testIterator()
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

         if( it == end( upper, 1UL ) || *it != 0 ) {
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

         const ptrdiff_t number( begin( upper, 0UL ) - end( upper, 0UL ) );

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

         const ptrdiff_t number( cend( upper, 1UL ) - cbegin( upper, 1UL ) );

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

         const ptrdiff_t number( cbegin( upper, 1UL ) - cend( upper, 1UL ) );

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

         ConstIterator it ( cbegin( upper, 0UL ) );
         ConstIterator end( cend( upper, 0UL ) );

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

      // Testing assignment to upper elements via Iterator
      {
         test_ = "Row-major assignment to upper elements via Iterator";

         int value = 0;

         for( Iterator it=begin( upper, 0UL ); it!=end( upper, 0UL ); ++it ) {
            *it = value++;
         }

         if( upper(0,0) != 0 || upper(0,1) !=  1 || upper(0,2) != 2 ||
             upper(1,0) != 0 || upper(1,1) != -2 || upper(1,2) != 0 ||
             upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 0  1  2 )\n( 0 -2  0 )\n( 0  0  4 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing assignment to lower elements via Iterator
      {
         test_ = "Row-major assignment to lower elements via Iterator";

         try {
            const Iterator it = begin( upper, 1UL );
            *it = 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to lower matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Testing addition assignment to upper elements via Iterator
      {
         test_ = "Row-major addition assignment to upper elements via Iterator";

         int value = -2;

         for( Iterator it=begin( upper, 0UL ); it!=end( upper, 0UL ); ++it ) {
            *it += value++;
         }

         if( upper(0,0) != -2 || upper(0,1) !=  0 || upper(0,2) != 2 ||
             upper(1,0) !=  0 || upper(1,1) != -2 || upper(1,2) != 0 ||
             upper(2,0) !=  0 || upper(2,1) !=  0 || upper(2,2) != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( -2  0  2 )\n(  0 -2  0 )\n(  0  0  4 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing addition assignment to lower elements via Iterator
      {
         test_ = "Row-major addition assignment to lower elements via Iterator";

         try {
            const Iterator it = begin( upper, 1UL );
            *it += 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to lower matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Testing subtraction assignment to upper elements via Iterator
      {
         test_ = "Row-major subtraction assignment to upper elements via Iterator";

         int value = -2;

         for( Iterator it=begin( upper, 0UL ); it!=end( upper, 0UL ); ++it ) {
            *it -= value++;
         }

         if( upper(0,0) != 0 || upper(0,1) !=  1 || upper(0,2) != 2 ||
             upper(1,0) != 0 || upper(1,1) != -2 || upper(1,2) != 0 ||
             upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 0  1  2 )\n( 0 -2  0 )\n( 0  0  4 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing subtraction assignment to lower elements via Iterator
      {
         test_ = "Row-major subtraction assignment to lower elements via Iterator";

         try {
            const Iterator it = begin( upper, 1UL );
            *it -= 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to lower matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Testing multiplication assignment to upper elements via Iterator
      {
         test_ = "Row-major multiplication assignment to upper elements via Iterator";

         int value = 1;

         for( Iterator it=begin( upper, 0UL ); it!=end( upper, 0UL ); ++it ) {
            *it *= value++;
         }

         if( upper(0,0) != 0 || upper(0,1) !=  2 || upper(0,2) != 6 ||
             upper(1,0) != 0 || upper(1,1) != -2 || upper(1,2) != 0 ||
             upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 0  2  6 )\n( 0 -2  0 )\n( 0  0  4 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing multiplication assignment to lower elements via Iterator
      {
         test_ = "Row-major multiplication assignment to lower elements via Iterator";

         try {
            const Iterator it = begin( upper, 1UL );
            *it *= 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to lower matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Testing division assignment to upper elements via Iterator
      {
         test_ = "Row-major division assignment to upper elements via Iterator";

         for( Iterator it=begin( upper, 0UL ); it!=end( upper, 0UL ); ++it ) {
            *it /= 2;
         }

         if( upper(0,0) != 0 || upper(0,1) !=  1 || upper(0,2) != 3 ||
             upper(1,0) != 0 || upper(1,1) != -2 || upper(1,2) != 0 ||
             upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 0  1  3 )\n( 0 -2  0 )\n( 0  0  4 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing division assignment to lower elements via Iterator
      {
         test_ = "Row-major division assignment to lower elements via Iterator";

         try {
            const Iterator it = begin( upper, 1UL );
            *it /= 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to lower matrix element succeeded\n"
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

         if( it == end( upper, 1UL ) || *it != 0 ) {
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

         const ptrdiff_t number( begin( upper, 0UL ) - end( upper, 0UL ) );

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

         const ptrdiff_t number( cend( upper, 1UL ) - cbegin( upper, 1UL ) );

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

         const ptrdiff_t number( cbegin( upper, 1UL ) - cend( upper, 1UL ) );

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

         ConstIterator it ( cbegin( upper, 2UL ) );
         ConstIterator end( cend( upper, 2UL ) );

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

         if( it == end || *it != 4 ) {
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

         if( it == end || *it != 4 ) {
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

      // Testing assignment to upper elements via Iterator
      {
         test_ = "Column-major assignment to upper elements via Iterator";

         int value = 0;

         for( Iterator it=begin( upper, 2UL ); it!=end( upper, 2UL ); ++it ) {
            *it = value++;
         }

         if( upper(0,0) != 1 || upper(0,1) !=  0 || upper(0,2) != 0 ||
             upper(1,0) != 0 || upper(1,1) != -2 || upper(1,2) != 1 ||
             upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 1  0  0 )\n( 0 -2  1 )\n( 0  0  2 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing assignment to lower elements via Iterator
      {
         test_ = "Column-major assignment to lower elements via Iterator";

         try {
            const Iterator it = begin( upper, 0UL ) + 1UL;
            *it = 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to lower matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Testing addition assignment to upper elements via Iterator
      {
         test_ = "Column-major addition assignment to upper elements via Iterator";

         int value = -2;

         for( Iterator it=begin( upper, 2UL ); it!=end( upper, 2UL ); ++it ) {
            *it += value++;
         }

         if( upper(0,0) != 1 || upper(0,1) !=  0 || upper(0,2) != -2 ||
             upper(1,0) != 0 || upper(1,1) != -2 || upper(1,2) !=  0 ||
             upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) !=  2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 1  0 -2 )\n( 0 -2  0 )\n( 0  0  2 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing addition assignment to lower elements via Iterator
      {
         test_ = "Column-major addition assignment to lower elements via Iterator";

         try {
            const Iterator it = begin( upper, 1UL ) + 2UL;
            *it += 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to lower matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Testing subtraction assignment to upper elements via Iterator
      {
         test_ = "Column-major subtraction assignment to upper elements via Iterator";

         int value = -2;

         for( Iterator it=begin( upper, 2UL ); it!=end( upper, 2UL ); ++it ) {
            *it -= value++;
         }

         if( upper(0,0) != 1 || upper(0,1) !=  0 || upper(0,2) != 0 ||
             upper(1,0) != 0 || upper(1,1) != -2 || upper(1,2) != 1 ||
             upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 1  0  0 )\n( 0 -2  1 )\n( 0  0  2 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing subtraction assignment to lower elements via Iterator
      {
         test_ = "Column-major subtraction assignment to lower elements via Iterator";

         try {
            const Iterator it = begin( upper, 0UL ) + 1UL;
            *it -= 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to lower matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Testing multiplication assignment to upper elements via Iterator
      {
         test_ = "Column-major multiplication assignment to upper elements via Iterator";

         int value = 1;

         for( Iterator it=begin( upper, 2UL ); it!=end( upper, 2UL ); ++it ) {
            *it *= value++;
         }

         if( upper(0,0) != 1 || upper(0,1) !=  0 || upper(0,2) != 0 ||
             upper(1,0) != 0 || upper(1,1) != -2 || upper(1,2) != 2 ||
             upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) != 6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 1  0  0 )\n( 0 -2  2 )\n( 0  0  6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing multiplication assignment to lower elements via Iterator
      {
         test_ = "Column-major multiplication assignment to lower elements via Iterator";

         try {
            const Iterator it = begin( upper, 1UL ) + 2UL;
            *it *= 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to lower matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Testing division assignment to upper elements via Iterator
      {
         test_ = "Column-major division assignment to upper elements via Iterator";

         for( Iterator it=begin( upper, 2UL ); it!=end( upper, 2UL ); ++it ) {
            *it /= 2;
         }

         if( upper(0,0) != 1 || upper(0,1) !=  0 || upper(0,2) != 0 ||
             upper(1,0) != 0 || upper(1,1) != -2 || upper(1,2) != 1 ||
             upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << upper << "\n"
                << "   Expected result:\n( 1  0  0 )\n( 0 -2  1 )\n( 0  0  3 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing division assignment to lower elements via Iterator
      {
         test_ = "Column-major division assignment to lower elements via Iterator";

         try {
            const Iterator it = begin( upper, 0UL ) + 1UL;
            *it /= 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to lower matrix element succeeded\n"
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
/*!\brief Test of the \c nonZeros() member function of the UpperMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c nonZeros() member function of the UpperMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testNonZeros()
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
         checkCapacity( upper, 9UL );
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
         checkCapacity( upper, 9UL );
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
         checkCapacity( upper, 9UL );
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
         checkCapacity( upper, 9UL );
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
         checkCapacity( upper, 9UL );
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
         checkCapacity( upper, 9UL );
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
void DenseTest::testReset()
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
      checkCapacity( upper, 9UL );
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
      checkCapacity( upper, 9UL );
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
      checkCapacity( upper, 9UL );
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
      checkCapacity( upper, 9UL );
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
      checkCapacity( upper, 9UL );
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
      checkCapacity( upper, 9UL );
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
      checkCapacity( upper, 9UL );
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
      checkCapacity( upper, 9UL );
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
      checkCapacity( upper, 9UL );
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
      checkCapacity( upper, 9UL );
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
void DenseTest::testClear()
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
      checkCapacity( upper, 9UL );
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
      checkCapacity( upper, 9UL );
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
      checkCapacity( upper, 9UL );
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
      checkCapacity( upper, 9UL );
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
      checkCapacity( upper, 9UL );
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
      checkCapacity( upper, 9UL );
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
void DenseTest::testResize()
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
      checkCapacity( upper, 4UL );

      if( upper(1,0) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( x x )\n( 0 x )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 4x4 and preserving the elements
      upper(0,0) = 1;
      upper(0,1) = 2;
      upper(1,1) = 3;
      upper.resize( 4UL, true );

      checkRows    ( upper,  4UL );
      checkColumns ( upper,  4UL );
      checkCapacity( upper, 16UL );

      if( upper(0,0) != 1 || upper(0,1) != 2 ||
          upper(1,0) != 0 || upper(1,1) != 3 ||
          upper(2,0) != 0 || upper(2,1) != 0 ||
          upper(3,0) != 0 || upper(3,1) != 0 || upper(3,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 2 x x )\n( 0 3 x x )\n( 0 0 x x )\n( 0 0 0 x )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 2x2
      upper(2,2) = 4;
      upper.resize( 2UL );

      checkRows    ( upper, 2UL );
      checkColumns ( upper, 2UL );
      checkCapacity( upper, 4UL );
      checkNonZeros( upper, 3UL );
      checkNonZeros( upper, 0UL, 2UL );
      checkNonZeros( upper, 1UL, 1UL );

      if( upper(0,0) != 1 || upper(0,1) != 2 ||
          upper(1,0) != 0 || upper(1,1) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 2 )\n( 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }

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
      checkCapacity( upper, 4UL );

      if( upper(1,0) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( x x )\n( 0 x )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 4x4 and preserving the elements
      upper(0,0) = 1;
      upper(0,1) = 2;
      upper(1,1) = 3;
      upper.resize( 4UL, true );

      checkRows    ( upper,  4UL );
      checkColumns ( upper,  4UL );
      checkCapacity( upper, 16UL );

      if( upper(0,0) != 1 || upper(0,1) != 2 ||
          upper(1,0) != 0 || upper(1,1) != 3 ||
          upper(2,0) != 0 || upper(2,1) != 0 ||
          upper(3,0) != 0 || upper(3,1) != 0 || upper(3,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 2 x x )\n( 0 3 x x )\n( 0 0 x x )\n( 0 0 0 x )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 2x2
      upper(2,2) = 4;
      upper.resize( 2UL );

      checkRows    ( upper, 2UL );
      checkColumns ( upper, 2UL );
      checkCapacity( upper, 4UL );
      checkNonZeros( upper, 3UL );
      checkNonZeros( upper, 0UL, 1UL );
      checkNonZeros( upper, 1UL, 2UL );

      if( upper(0,0) != 1 || upper(0,1) != 2 ||
          upper(1,0) != 0 || upper(1,1) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 2 )\n( 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 0x0
      upper.resize( 0UL );

      checkRows    ( upper, 0UL );
      checkColumns ( upper, 0UL );
      checkNonZeros( upper, 0UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c extend() member function of the UpperMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c extend() member function of the UpperMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testExtend()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major UpperMatrix::extend()";

      // Initialization check
      UT upper;

      checkRows    ( upper, 0UL );
      checkColumns ( upper, 0UL );
      checkNonZeros( upper, 0UL );

      // Extending the size of the matrix to 2x2
      upper.extend( 2UL );

      checkRows    ( upper, 2UL );
      checkColumns ( upper, 2UL );
      checkCapacity( upper, 4UL );

      if( upper(1,0) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Extending the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( x x )\n( 0 x )\n";
         throw std::runtime_error( oss.str() );
      }

      // Extending to 4x4 and preserving the elements
      upper(0,0) = 1;
      upper(0,1) = 2;
      upper(1,1) = 3;
      upper.extend( 2UL, true );

      checkRows    ( upper,  4UL );
      checkColumns ( upper,  4UL );
      checkCapacity( upper, 16UL );

      if( upper(0,0) != 1 || upper(0,1) != 2 ||
          upper(1,0) != 0 || upper(1,1) != 3 ||
          upper(2,0) != 0 || upper(2,1) != 0 ||
          upper(3,0) != 0 || upper(3,1) != 0 || upper(3,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Extending the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 2 x x )\n( 0 3 x x )\n( 0 0 x x )\n( 0 0 0 x )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major UpperMatrix::extend()";

      // Initialization check
      OUT upper;

      checkRows    ( upper, 0UL );
      checkColumns ( upper, 0UL );
      checkNonZeros( upper, 0UL );

      // Extending the size of the matrix to 2x2
      upper.extend( 2UL );

      checkRows    ( upper, 2UL );
      checkColumns ( upper, 2UL );
      checkCapacity( upper, 4UL );

      if( upper(1,0) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Extending the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( x x )\n( 0 x )\n";
         throw std::runtime_error( oss.str() );
      }

      // Extending to 4x4 and preserving the elements
      upper(0,0) = 1;
      upper(0,1) = 2;
      upper(1,1) = 3;
      upper.extend( 2UL, true );

      checkRows    ( upper,  4UL );
      checkColumns ( upper,  4UL );
      checkCapacity( upper, 16UL );

      if( upper(0,0) != 1 || upper(0,1) != 2 ||
          upper(1,0) != 0 || upper(1,1) != 3 ||
          upper(2,0) != 0 || upper(2,1) != 0 ||
          upper(3,0) != 0 || upper(3,1) != 0 || upper(3,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Extending the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 2 x x )\n( 0 3 x x )\n( 0 0 x x )\n( 0 0 0 x )\n";
         throw std::runtime_error( oss.str() );
      }
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
void DenseTest::testReserve()
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
/*!\brief Test of the \c shrinkToFit() member function of the UpperMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c shrinkToFit() member function of the UpperMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testShrinkToFit()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major UpperMatrix::shrinkToFit()";

      // Shrinking a matrix without excessive capacity
      {
         UT upper{ { 1, 2, 3 }, { 0, 4, 5 }, { 0, 0, 6 } };

         upper.shrinkToFit();

         checkRows    ( upper, 3UL );
         checkColumns ( upper, 3UL );
         checkCapacity( upper, 9UL );
         checkNonZeros( upper, 6UL );
         checkNonZeros( upper, 0UL, 3UL );
         checkNonZeros( upper, 1UL, 2UL );
         checkNonZeros( upper, 2UL, 1UL );

         if( upper.capacity() != upper.rows() * upper.spacing() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Shrinking the matrix failed\n"
                << " Details:\n"
                << "   Capacity         : " << upper.capacity() << "\n"
                << "   Expected capacity: " << ( upper.rows() * upper.spacing() ) << "\n";
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
         UT upper{ { 1, 2, 3 }, { 0, 4, 5 }, { 0, 0, 6 } };
         upper.reserve( 100UL );

         upper.shrinkToFit();

         checkRows    ( upper, 3UL );
         checkColumns ( upper, 3UL );
         checkCapacity( upper, 9UL );
         checkNonZeros( upper, 6UL );
         checkNonZeros( upper, 0UL, 3UL );
         checkNonZeros( upper, 1UL, 2UL );
         checkNonZeros( upper, 2UL, 1UL );

         if( upper.capacity() != upper.rows() * upper.spacing() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Shrinking the matrix failed\n"
                << " Details:\n"
                << "   Capacity         : " << upper.capacity() << "\n"
                << "   Expected capacity: " << ( upper.rows() * upper.spacing() ) << "\n";
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
         OUT upper{ { 1, 2, 3 }, { 0, 4, 5 }, { 0, 0, 6 } };

         upper.shrinkToFit();

         checkRows    ( upper, 3UL );
         checkColumns ( upper, 3UL );
         checkCapacity( upper, 9UL );
         checkNonZeros( upper, 6UL );
         checkNonZeros( upper, 0UL, 1UL );
         checkNonZeros( upper, 1UL, 2UL );
         checkNonZeros( upper, 2UL, 3UL );

         if( upper.capacity() != upper.spacing() * upper.columns() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Shrinking the matrix failed\n"
                << " Details:\n"
                << "   Capacity         : " << upper.capacity() << "\n"
                << "   Expected capacity: " << ( upper.spacing() * upper.columns() ) << "\n";
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
         OUT upper{ { 1, 2, 3 }, { 0, 4, 5 }, { 0, 0, 6 } };
         upper.reserve( 100UL );

         upper.shrinkToFit();

         checkRows    ( upper, 3UL );
         checkColumns ( upper, 3UL );
         checkCapacity( upper, 9UL );
         checkNonZeros( upper, 6UL );
         checkNonZeros( upper, 0UL, 1UL );
         checkNonZeros( upper, 1UL, 2UL );
         checkNonZeros( upper, 2UL, 3UL );

         if( upper.capacity() != upper.spacing() * upper.columns() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Shrinking the matrix failed\n"
                << " Details:\n"
                << "   Capacity         : " << upper.capacity() << "\n"
                << "   Expected capacity: " << ( upper.spacing() * upper.columns() ) << "\n";
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
void DenseTest::testSwap()
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
      checkCapacity( upper1, 4UL );
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
      checkCapacity( upper1, 4UL );
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
/*!\brief Test of the \c isDefault() function with the UpperMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isDefault() function with the UpperMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testIsDefault()
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
void DenseTest::testSubmatrix()
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

      if( it == sm.end(0UL) || *it != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << *it << "\n"
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
   // Row-major scalar assignment
   //=====================================================================================

   // ( 1 -4  7 -2 )      ( 1 12 12 -2 )
   // ( 0  2  0  0 )  =>  ( 0 12 12  0 )
   // ( 0  0  3  1 )      ( 0  0 12  1 )
   // ( 0  0  0  5 )      ( 0  0  0  5 )
   {
      test_ = "Row-major submatrix() function (scalar assignment test 1)";

      using SMT = blaze::Submatrix<UT>;

      UT upper( 4UL );
      upper(0,0) =  1;
      upper(0,1) = -4;
      upper(0,2) =  7;
      upper(0,3) = -2;
      upper(1,1) =  2;
      upper(2,2) =  3;
      upper(2,3) =  1;
      upper(3,3) =  5;

      SMT sm = submatrix( upper, 0UL, 1UL, 4UL, 2UL );
      sm = 12;

      checkRows    ( upper, 4UL );
      checkColumns ( upper, 4UL );
      checkNonZeros( upper, 9UL );
      checkNonZeros( upper, 0UL, 4UL );
      checkNonZeros( upper, 1UL, 2UL );
      checkNonZeros( upper, 2UL, 2UL );
      checkNonZeros( upper, 3UL, 1UL );

      if( sm(0,0) != 12 || sm(0,1) != 12 ||
          sm(1,0) != 12 || sm(1,1) != 12 ||
          sm(2,0) !=  0 || sm(2,1) != 12 ||
          sm(3,0) !=  0 || sm(3,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 12 12 )\n( 12 12 )\n(  0 12 )\n(  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( upper(0,0) != 1 || upper(0,1) != 12 || upper(0,2) != 12 || upper(0,3) != -2 ||
          upper(1,0) != 0 || upper(1,1) != 12 || upper(1,2) != 12 || upper(1,3) !=  0 ||
          upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) != 12 || upper(2,3) !=  1 ||
          upper(3,0) != 0 || upper(3,1) !=  0 || upper(3,2) !=  0 || upper(3,3) !=  5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 12 12 -2 )\n"
                                     "( 0 12 12  0 )\n"
                                     "( 0  0 12  1 )\n"
                                     "( 0  0  0  5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1 -4  7 -2 )      ( 1 -4  7 -2 )
   // ( 0  2  0  0 )  =>  ( 0 12 12 12 )
   // ( 0  0  3  1 )      ( 0  0 12 12 )
   // ( 0  0  0  5 )      ( 0  0  0  5 )
   {
      test_ = "Row-major submatrix() function (scalar assignment test 2)";

      using SMT = blaze::Submatrix<UT>;

      UT upper( 4UL );
      upper(0,0) =  1;
      upper(0,1) = -4;
      upper(0,2) =  7;
      upper(0,3) = -2;
      upper(1,1) =  2;
      upper(2,2) =  3;
      upper(2,3) =  1;
      upper(3,3) =  5;

      SMT sm = submatrix( upper, 1UL, 0UL, 2UL, 4UL );
      sm = 12;

      checkRows    ( upper,  4UL );
      checkColumns ( upper,  4UL );
      checkNonZeros( upper, 10UL );
      checkNonZeros( upper,  0UL, 4UL );
      checkNonZeros( upper,  1UL, 3UL );
      checkNonZeros( upper,  2UL, 2UL );
      checkNonZeros( upper,  3UL, 1UL );

      if( sm(0,0) != 0 || sm(0,1) != 12 || sm(0,2) != 12 || sm(0,3) != 12 ||
          sm(1,0) != 0 || sm(1,1) !=  0 || sm(1,2) != 12 || sm(1,3) != 12 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0 12 12 12 )\n( 0  0 12 12 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( upper(0,0) != 1 || upper(0,1) != -4 || upper(0,2) !=  7 || upper(0,3) != -2 ||
          upper(1,0) != 0 || upper(1,1) != 12 || upper(1,2) != 12 || upper(1,3) != 12 ||
          upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) != 12 || upper(2,3) != 12 ||
          upper(3,0) != 0 || upper(3,1) !=  0 || upper(3,2) !=  0 || upper(3,3) !=  5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 -4  7 -2 )\n"
                                     "( 0 12 12 12 )\n"
                                     "( 0  0 12 12 )\n"
                                     "( 0  0  0  5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1 -4  7 -2 )      ( 1 -4  7 -2 )
   // ( 0  2  0  0 )  =>  ( 0  2  0  0 )
   // ( 0  0  3  1 )      ( 0  0  3  1 )
   // ( 0  0  0  5 )      ( 0  0  0  5 )
   {
      test_ = "Row-major submatrix() function (scalar assignment test 3)";

      using SMT = blaze::Submatrix<UT>;

      UT upper( 4UL );
      upper(0,0) =  1;
      upper(0,1) = -4;
      upper(0,2) =  7;
      upper(0,3) = -2;
      upper(1,1) =  2;
      upper(2,2) =  3;
      upper(2,3) =  1;
      upper(3,3) =  5;

      SMT sm = submatrix( upper, 2UL, 0UL, 2UL, 2UL );
      sm = 12;

      checkRows    ( upper, 4UL );
      checkColumns ( upper, 4UL );
      checkNonZeros( upper, 8UL );
      checkNonZeros( upper, 0UL, 4UL );
      checkNonZeros( upper, 1UL, 1UL );
      checkNonZeros( upper, 2UL, 2UL );
      checkNonZeros( upper, 3UL, 1UL );

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

      if( upper(0,0) != 1 || upper(0,1) != -4 || upper(0,2) != 7 || upper(0,3) != -2 ||
          upper(1,0) != 0 || upper(1,1) !=  2 || upper(1,2) != 0 || upper(1,3) !=  0 ||
          upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) != 3 || upper(2,3) !=  1 ||
          upper(3,0) != 0 || upper(3,1) !=  0 || upper(3,2) != 0 || upper(3,3) !=  5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 -4  7 -2 )\n"
                                     "( 0  2  0  0 )\n"
                                     "( 0  0  3  1 )\n"
                                     "( 0  0  0  5 )\n";
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

      if( it == sm.end(0UL) || *it != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << *it << "\n"
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
   // Column-major scalar assignment
   //=====================================================================================

   // ( 1 -4  7 -2 )      ( 1 12 12 -2 )
   // ( 0  2  0  0 )  =>  ( 0 12 12  0 )
   // ( 0  0  3  1 )      ( 0  0 12  1 )
   // ( 0  0  0  5 )      ( 0  0  0  5 )
   {
      test_ = "Column-major submatrix() function (scalar assignment test 1)";

      using SMT = blaze::Submatrix<OUT>;

      OUT upper( 4UL );
      upper(0,0) =  1;
      upper(0,1) = -4;
      upper(0,2) =  7;
      upper(0,3) = -2;
      upper(1,1) =  2;
      upper(2,2) =  3;
      upper(2,3) =  1;
      upper(3,3) =  5;

      SMT sm = submatrix( upper, 0UL, 1UL, 4UL, 2UL );
      sm = 12;

      checkRows    ( upper, 4UL );
      checkColumns ( upper, 4UL );
      checkNonZeros( upper, 9UL );
      checkNonZeros( upper, 0UL, 1UL );
      checkNonZeros( upper, 1UL, 2UL );
      checkNonZeros( upper, 2UL, 3UL );
      checkNonZeros( upper, 3UL, 3UL );

      if( sm(0,0) != 12 || sm(0,1) != 12 ||
          sm(1,0) != 12 || sm(1,1) != 12 ||
          sm(2,0) !=  0 || sm(2,1) != 12 ||
          sm(3,0) !=  0 || sm(3,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 12 12 )\n( 12 12 )\n(  0 12 )\n(  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( upper(0,0) != 1 || upper(0,1) != 12 || upper(0,2) != 12 || upper(0,3) != -2 ||
          upper(1,0) != 0 || upper(1,1) != 12 || upper(1,2) != 12 || upper(1,3) !=  0 ||
          upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) != 12 || upper(2,3) !=  1 ||
          upper(3,0) != 0 || upper(3,1) !=  0 || upper(3,2) !=  0 || upper(3,3) !=  5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 12 12 -2 )\n"
                                     "( 0 12 12  0 )\n"
                                     "( 0  0 12  1 )\n"
                                     "( 0  0  0  5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1 -4  7 -2 )      ( 1 -4  7 -2 )
   // ( 0  2  0  0 )  =>  ( 0 12 12 12 )
   // ( 0  0  3  1 )      ( 0  0 12 12 )
   // ( 0  0  0  5 )      ( 0  0  0  5 )
   {
      test_ = "Column-major submatrix() function (scalar assignment test 2)";

      using SMT = blaze::Submatrix<OUT>;

      OUT upper( 4UL );
      upper(0,0) =  1;
      upper(0,1) = -4;
      upper(0,2) =  7;
      upper(0,3) = -2;
      upper(1,1) =  2;
      upper(2,2) =  3;
      upper(2,3) =  1;
      upper(3,3) =  5;

      SMT sm = submatrix( upper, 1UL, 0UL, 2UL, 4UL );
      sm = 12;

      checkRows    ( upper,  4UL );
      checkColumns ( upper,  4UL );
      checkNonZeros( upper, 10UL );
      checkNonZeros( upper,  0UL, 1UL );
      checkNonZeros( upper,  1UL, 2UL );
      checkNonZeros( upper,  2UL, 3UL );
      checkNonZeros( upper,  3UL, 4UL );

      if( sm(0,0) != 0 || sm(0,1) != 12 || sm(0,2) != 12 || sm(0,3) != 12 ||
          sm(1,0) != 0 || sm(1,1) !=  0 || sm(1,2) != 12 || sm(1,3) != 12 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0 12 12 12 )\n( 0  0 12 12 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( upper(0,0) != 1 || upper(0,1) != -4 || upper(0,2) !=  7 || upper(0,3) != -2 ||
          upper(1,0) != 0 || upper(1,1) != 12 || upper(1,2) != 12 || upper(1,3) != 12 ||
          upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) != 12 || upper(2,3) != 12 ||
          upper(3,0) != 0 || upper(3,1) !=  0 || upper(3,2) !=  0 || upper(3,3) !=  5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 -4  7 -2 )\n"
                                     "( 0 12 12 12 )\n"
                                     "( 0  0 12 12 )\n"
                                     "( 0  0  0  5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1 -4  7 -2 )      ( 1 -4  7 -2 )
   // ( 0  2  0  0 )  =>  ( 0  2  0  0 )
   // ( 0  0  3  1 )      ( 0  0  3  1 )
   // ( 0  0  0  5 )      ( 0  0  0  5 )
   {
      test_ = "Column-major submatrix() function (scalar assignment test 3)";

      using SMT = blaze::Submatrix<OUT>;

      OUT upper( 4UL );
      upper(0,0) =  1;
      upper(0,1) = -4;
      upper(0,2) =  7;
      upper(0,3) = -2;
      upper(1,1) =  2;
      upper(2,2) =  3;
      upper(2,3) =  1;
      upper(3,3) =  5;

      SMT sm = submatrix( upper, 2UL, 0UL, 2UL, 2UL );
      sm = 12;

      checkRows    ( upper, 4UL );
      checkColumns ( upper, 4UL );
      checkNonZeros( upper, 8UL );
      checkNonZeros( upper, 0UL, 1UL );
      checkNonZeros( upper, 1UL, 2UL );
      checkNonZeros( upper, 2UL, 2UL );
      checkNonZeros( upper, 3UL, 3UL );

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

      if( upper(0,0) != 1 || upper(0,1) != -4 || upper(0,2) != 7 || upper(0,3) != -2 ||
          upper(1,0) != 0 || upper(1,1) !=  2 || upper(1,2) != 0 || upper(1,3) !=  0 ||
          upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) != 3 || upper(2,3) !=  1 ||
          upper(3,0) != 0 || upper(3,1) !=  0 || upper(3,2) != 0 || upper(3,3) !=  5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 -4  7 -2 )\n"
                                     "( 0  2  0  0 )\n"
                                     "( 0  0  3  1 )\n"
                                     "( 0  0  0  5 )\n";
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
void DenseTest::testRow()
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

      if( it == row1.end() || *it != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << *it << "\n"
             << "   Expected result: 0\n";
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
   // Row-major scalar assignment
   //=====================================================================================

   {
      test_ = "Row-major row() function (scalar assignment test)";

      using RT = blaze::Row<UT>;

      UT upper( 3UL );
      upper(0,0) =  1;
      upper(0,1) = -4;
      upper(0,2) =  7;
      upper(1,1) =  2;
      upper(2,2) =  3;

      RT row1 = row( upper, 1UL );
      row1 = 8;

      checkRows    ( upper, 3UL );
      checkColumns ( upper, 3UL );
      checkNonZeros( upper, 6UL );
      checkNonZeros( upper, 0UL, 3UL );
      checkNonZeros( upper, 1UL, 2UL );
      checkNonZeros( upper, 2UL, 1UL );

      if( row1[0] != 0 || row1[1] != 8 || row1[2] != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row access failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( 0 8 8 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( upper(0,0) != 1 || upper(0,1) != -4 || upper(0,2) != 7 ||
          upper(1,0) != 0 || upper(1,1) !=  8 || upper(1,2) != 8 ||
          upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row access failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 -4  7 )\n( 0  8  8 )\n( 0  0  3 )\n";
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

      if( it == row1.end() || *it != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << *it << "\n"
             << "   Expected result: 0\n";
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
   // Column-major scalar assignment
   //=====================================================================================

   {
      test_ = "Column-major row() function (scalar assignment test)";

      using RT = blaze::Row<OUT>;

      OUT upper( 3UL );
      upper(0,0) =  1;
      upper(0,1) = -4;
      upper(0,2) =  7;
      upper(1,1) =  2;
      upper(2,2) =  3;

      RT row1 = row( upper, 1UL );
      row1 = 8;

      checkRows    ( upper, 3UL );
      checkColumns ( upper, 3UL );
      checkNonZeros( upper, 6UL );
      checkNonZeros( upper, 0UL, 1UL );
      checkNonZeros( upper, 1UL, 2UL );
      checkNonZeros( upper, 2UL, 3UL );

      if( row1[0] != 0 || row1[1] != 8 || row1[2] != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row access failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( 0 8 8 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( upper(0,0) != 1 || upper(0,1) != -4 || upper(0,2) != 7 ||
          upper(1,0) != 0 || upper(1,1) !=  8 || upper(1,2) != 8 ||
          upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row access failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 -4  7 )\n( 0  8  8 )\n( 0  0  3 )\n";
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
void DenseTest::testColumn()
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

      if( it == col1.end() || *it != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << *it << "\n"
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
   // Row-major scalar assignment
   //=====================================================================================

   {
      test_ = "Row-major column() function (scalar assignment test)";

      using CT = blaze::Column<UT>;

      UT upper( 3UL );
      upper(0,0) =  1;
      upper(0,1) = -4;
      upper(0,2) =  7;
      upper(1,1) =  2;
      upper(2,2) =  3;

      CT col1 = column( upper, 1UL );
      col1 = 8;

      checkRows    ( upper, 3UL );
      checkColumns ( upper, 3UL );
      checkNonZeros( upper, 5UL );
      checkNonZeros( upper, 0UL, 3UL );
      checkNonZeros( upper, 1UL, 1UL );
      checkNonZeros( upper, 2UL, 1UL );

      if( col1[0] != 8 || col1[1] != 8 || col1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column access failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( 8 8 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( upper(0,0) != 1 || upper(0,1) != 8 || upper(0,2) != 7 ||
          upper(1,0) != 0 || upper(1,1) != 8 || upper(1,2) != 0 ||
          upper(2,0) != 0 || upper(2,1) != 0 || upper(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column access failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 8 7 )\n( 0 8 0 )\n( 0 0 3 )\n";
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

      if( it == col1.end() || *it != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << *it << "\n"
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


   //=====================================================================================
   // Column-major scalar assignment
   //=====================================================================================

   {
      test_ = "Column-major column() function (scalar assignment test)";

      using CT = blaze::Column<OUT>;

      OUT upper( 3UL );
      upper(0,0) =  1;
      upper(0,1) = -4;
      upper(0,2) =  7;
      upper(1,1) =  2;
      upper(2,2) =  3;

      CT col1 = column( upper, 1UL );
      col1 = 8;

      checkRows    ( upper, 3UL );
      checkColumns ( upper, 3UL );
      checkNonZeros( upper, 5UL );
      checkNonZeros( upper, 0UL, 1UL );
      checkNonZeros( upper, 1UL, 2UL );
      checkNonZeros( upper, 2UL, 2UL );

      if( col1[0] != 8 || col1[1] != 8 || col1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column access failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( 8 8 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( upper(0,0) != 1 || upper(0,1) != 8 || upper(0,2) != 7 ||
          upper(1,0) != 0 || upper(1,1) != 8 || upper(1,2) != 0 ||
          upper(2,0) != 0 || upper(2,1) != 0 || upper(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column access failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n(  1 8 7 )\n( 0 8 0 )\n( 0 0 3 )\n";
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
   std::cout << "   Running UpperMatrix dense test (part 2)..." << std::endl;

   try
   {
      RUN_UPPERMATRIX_DENSE_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during UpperMatrix dense test (part 2):\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
