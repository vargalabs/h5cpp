//=================================================================================================
/*!
//  \file src/mathtest/adaptors/diagonalmatrix/SparseTest2.cpp
//  \brief Source file for the DiagonalMatrix sparse test (part 2)
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
#include <blazetest/mathtest/adaptors/diagonalmatrix/SparseTest.h>

#ifdef BLAZE_USE_HPX_THREADS
#  include <hpx/hpx_main.hpp>
#endif


namespace blazetest {

namespace mathtest {

namespace adaptors {

namespace diagonalmatrix {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the DiagonalMatrix sparse test.
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
/*!\brief Test of all DiagonalMatrix (self-)scaling operations.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all available ways to scale an instance of the DiagonalMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testScaling()
{
   //=====================================================================================
   // Row-major self-scaling (M*=s)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (M*=s)";

      DT diag( 3UL );
      diag(1,1) =  2;
      diag(2,2) = -3;

      diag *= 2;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 2UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 0UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) !=  0 ||
          diag(1,0) != 0 || diag(1,1) != 4 || diag(1,2) !=  0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 0  0  0 )\n( 0  4  0 )\n( 0  0 -6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major self-scaling (M=M*s)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (M=M*s)";

      DT diag( 3UL );
      diag(1,1) =  2;
      diag(2,2) = -3;

      diag = diag * 2;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 2UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 0UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) !=  0 ||
          diag(1,0) != 0 || diag(1,1) != 4 || diag(1,2) !=  0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 0  0  0 )\n( 0  4  0 )\n( 0  0 -6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major self-scaling (M=s*M)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (M=s*M)";

      DT diag( 3UL );
      diag(1,1) =  2;
      diag(2,2) = -3;

      diag = 2 * diag;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 2UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 0UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) !=  0 ||
          diag(1,0) != 0 || diag(1,1) != 4 || diag(1,2) !=  0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 0  0  0 )\n( 0  4  0 )\n( 0  0 -6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major self-scaling (M/=s)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (M/=s)";

      DT diag( 3UL );
      diag(1,1) =  4;
      diag(2,2) = -6;

      diag /= 2;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 2UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 0UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) !=  0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) !=  0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != -3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 0  0  0 )\n( 0  2  0 )\n( 0  0 -3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major self-scaling (M=M/s)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (M=M/s)";

      DT diag( 3UL );
      diag(1,1) =  4;
      diag(2,2) = -6;

      diag = diag / 2;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 2UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 0UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) !=  0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) !=  0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != -3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 0  0  0 )\n( 0  2  0 )\n( 0  0 -3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major DiagonalMatrix::scale()
   //=====================================================================================

   {
      test_ = "Row-major DiagonalMatrix::scale()";

      // Initialization check
      DT diag( 3UL );
      diag(1,1) =  2;
      diag(2,2) = -3;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 2UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 0UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) !=  0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) !=  0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != -3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 0  0  0 )\n( 0  2  0 )\n( 0  0 -3 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Integral scaling of the matrix
      diag.scale( 2 );

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 2UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 0UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) !=  0 ||
          diag(1,0) != 0 || diag(1,1) != 4 || diag(1,2) !=  0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 0  0  0 )\n( 0  4  0 )\n( 0  0 -6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Floating point scaling of the matrix
      diag.scale( 0.5 );

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 2UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 0UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) !=  0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) !=  0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != -3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 0  0  0 )\n( 0  2  0 )\n( 0  0 -3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major DiagonalMatrix::scale() (complex)";

      using blaze::complex;

      blaze::DiagonalMatrix< blaze::CompressedMatrix<complex<float>,blaze::rowMajor> > diag( 2UL );
      diag(0,0) = complex<float>( 1.0F, 0.0F );
      diag(1,1) = complex<float>( 2.0F, 0.0F );

      diag.scale( complex<float>( 3.0F, 0.0F ) );

      checkRows    ( diag, 2UL );
      checkColumns ( diag, 2UL );
      checkCapacity( diag, 2UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );

      if( diag(0,0) != complex<float>( 3.0F, 0.0F ) || diag(0,1) != complex<float>( 0.0F, 0.0F ) ||
          diag(1,0) != complex<float>( 0.0F, 0.0F ) || diag(1,1) != complex<float>( 6.0F, 0.0F ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( (3,0) (0,0)\n(0,0) (6,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (M*=s)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (M*=s)";

      ODT diag( 3UL );
      diag(1,1) =  2;
      diag(2,2) = -3;

      diag *= 2;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 2UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 0UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) !=  0 ||
          diag(1,0) != 0 || diag(1,1) != 4 || diag(1,2) !=  0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 0  0  0 )\n( 0  4  0 )\n( 0  0 -6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (M=M*s)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (M=M*s)";

      ODT diag( 3UL );
      diag(1,1) =  2;
      diag(2,2) = -3;

      diag = diag * 2;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 2UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 0UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) !=  0 ||
          diag(1,0) != 0 || diag(1,1) != 4 || diag(1,2) !=  0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 0  0  0 )\n( 0  4  0 )\n( 0  0 -6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (M=s*M)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (M=s*M)";

      ODT diag( 3UL );
      diag(1,1) =  2;
      diag(2,2) = -3;

      diag = 2 * diag;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 2UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 0UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) !=  0 ||
          diag(1,0) != 0 || diag(1,1) != 4 || diag(1,2) !=  0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 0  0  0 )\n( 0  4  0 )\n( 0  0 -6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (M/=s)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (M/=s)";

      ODT diag( 3UL );
      diag(1,1) =  4;
      diag(2,2) = -6;

      diag /= 2;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 2UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 0UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) !=  0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) !=  0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != -3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 0  0  0 )\n( 0  2  0 )\n( 0  0 -3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (M=M/s)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (M=M/s)";

      ODT diag( 3UL );
      diag(1,1) =  4;
      diag(2,2) = -6;

      diag = diag / 2;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 2UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 0UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) !=  0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) !=  0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != -3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 0  0  0 )\n( 0  2  0 )\n( 0  0 -3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major DiagonalMatrix::scale()
   //=====================================================================================

   {
      test_ = "Column-major DiagonalMatrix::scale()";

      // Initialization check
      ODT diag( 3UL );
      diag(1,1) =  2;
      diag(2,2) = -3;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 2UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 0UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) !=  0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) !=  0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != -3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 0  0  0 )\n( 0  2  0 )\n( 0  0 -3 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Integral scaling of the matrix
      diag.scale( 2 );

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 2UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 0UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) !=  0 ||
          diag(1,0) != 0 || diag(1,1) != 4 || diag(1,2) !=  0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 0  0  0 )\n( 0  4  0 )\n( 0  0 -6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Floating point scaling of the matrix
      diag.scale( 0.5 );

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 2UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 0UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) !=  0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) !=  0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != -3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 0  0  0 )\n( 0  2  0 )\n( 0  0 -3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major DiagonalMatrix::scale() (complex)";

      using blaze::complex;

      blaze::DiagonalMatrix< blaze::CompressedMatrix<complex<float>,blaze::columnMajor> > diag( 2UL );
      diag(0,0) = complex<float>( 1.0F, 0.0F );
      diag(1,1) = complex<float>( 2.0F, 0.0F );

      diag.scale( complex<float>( 3.0F, 0.0F ) );

      checkRows    ( diag, 2UL );
      checkColumns ( diag, 2UL );
      checkCapacity( diag, 2UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );

      if( diag(0,0) != complex<float>( 3.0F, 0.0F ) || diag(0,1) != complex<float>( 0.0F, 0.0F ) ||
          diag(1,0) != complex<float>( 0.0F, 0.0F ) || diag(1,1) != complex<float>( 6.0F, 0.0F ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( (3,0) (0,0)\n(0,0) (6,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the DiagonalMatrix function call operator.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of adding and accessing elements via the function call operator
// of the DiagonalMatrix specialization. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
void SparseTest::testFunctionCall()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major DiagonalMatrix::operator()";

      // Good cases
      {
         DT diag( 3UL );

         // Writing the element (1,1)
         diag(1,1) = 1;

         checkRows    ( diag, 3UL );
         checkColumns ( diag, 3UL );
         checkCapacity( diag, 1UL );
         checkNonZeros( diag, 1UL );
         checkNonZeros( diag, 0UL, 0UL );
         checkNonZeros( diag, 1UL, 1UL );
         checkNonZeros( diag, 2UL, 0UL );

         if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) != 0 ||
             diag(1,0) != 0 || diag(1,1) != 1 || diag(1,2) != 0 ||
             diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 0 0 0 )\n( 0 1 0 )\n( 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Writing the element (2,2)
         diag(2,2) = 2;

         checkRows    ( diag, 3UL );
         checkColumns ( diag, 3UL );
         checkCapacity( diag, 2UL );
         checkNonZeros( diag, 2UL );
         checkNonZeros( diag, 0UL, 0UL );
         checkNonZeros( diag, 1UL, 1UL );
         checkNonZeros( diag, 2UL, 1UL );

         if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) != 0 ||
             diag(1,0) != 0 || diag(1,1) != 1 || diag(1,2) != 0 ||
             diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 0 0 0 )\n( 0 1 0 )\n( 0 0 2 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Adding to the element (0,0)
         diag(0,0) += 3;

         checkRows    ( diag, 3UL );
         checkColumns ( diag, 3UL );
         checkCapacity( diag, 3UL );
         checkNonZeros( diag, 3UL );
         checkNonZeros( diag, 0UL, 1UL );
         checkNonZeros( diag, 1UL, 1UL );
         checkNonZeros( diag, 2UL, 1UL );

         if( diag(0,0) != 3 || diag(0,1) != 0 || diag(0,2) != 0 ||
             diag(1,0) != 0 || diag(1,1) != 1 || diag(1,2) != 0 ||
             diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 3 0 0 )\n( 0 1 0 )\n( 0 0 2 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Subtracting from the element (1,1)
         diag(1,1) -= 4;

         checkRows    ( diag, 3UL );
         checkColumns ( diag, 3UL );
         checkCapacity( diag, 3UL );
         checkNonZeros( diag, 3UL );
         checkNonZeros( diag, 0UL, 1UL );
         checkNonZeros( diag, 1UL, 1UL );
         checkNonZeros( diag, 2UL, 1UL );

         if( diag(0,0) != 3 || diag(0,1) !=  0 || diag(0,2) != 0 ||
             diag(1,0) != 0 || diag(1,1) != -3 || diag(1,2) != 0 ||
             diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 3  0  0 )\n( 0 -3  0 )\n( 0  0  2 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Multiplying the element (2,2)
         diag(2,2) *= -3;

         checkRows    ( diag, 3UL );
         checkColumns ( diag, 3UL );
         checkCapacity( diag, 3UL );
         checkNonZeros( diag, 3UL );
         checkNonZeros( diag, 0UL, 1UL );
         checkNonZeros( diag, 1UL, 1UL );
         checkNonZeros( diag, 2UL, 1UL );

         if( diag(0,0) != 3 || diag(0,1) !=  0 || diag(0,2) !=  0 ||
             diag(1,0) != 0 || diag(1,1) != -3 || diag(1,2) !=  0 ||
             diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != -6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 3  0  0 )\n( 0 -3  0 )\n( 0  0 -6 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Dividing the element (2,2)
         diag(2,2) /= 2;

         checkRows    ( diag, 3UL );
         checkColumns ( diag, 3UL );
         checkCapacity( diag, 3UL );
         checkNonZeros( diag, 3UL );
         checkNonZeros( diag, 0UL, 1UL );
         checkNonZeros( diag, 1UL, 1UL );
         checkNonZeros( diag, 2UL, 1UL );

         if( diag(0,0) != 3 || diag(0,1) !=  0 || diag(0,2) !=  0 ||
             diag(1,0) != 0 || diag(1,1) != -3 || diag(1,2) !=  0 ||
             diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != -3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 3  0  0 )\n( 0 -3  0 )\n( 0  0 -3 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Failure cases
      {
         DT diag( 3UL );

         // Trying to write the lower element (2,1)
         try {
            diag(2,1) = 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to lower matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to add to the lower element (2,1)
         try {
            diag(2,1) += 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Addition assignment to lower matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to subtract from the lower element (2,1)
         try {
            diag(2,1) -= 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subtraction assignment to lower matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to multiply the lower element (2,1)
         try {
            diag(2,1) *= 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Multiplication assignment to lower matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to divide the lower element (2,1)
         try {
            diag(2,1) /= 2;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Division assignment to lower matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to write the upper element (1,2)
         try {
            diag(1,2) = 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to upper matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to add to the upper element (1,2)
         try {
            diag(1,2) += 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Addition assignment to upper matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to subtract from the upper element (1,2)
         try {
            diag(1,2) -= 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subtraction assignment to upper matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to multiply the upper element (1,2)
         try {
            diag(1,2) *= 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Multiplication assignment to upper matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to divide the upper element (1,2)
         try {
            diag(1,2) /= 2;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Division assignment to upper matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major DiagonalMatrix::operator()";

      // Good cases
      {
         ODT diag( 3UL );

         // Writing the element (1,1)
         diag(1,1) = 1;

         checkRows    ( diag, 3UL );
         checkColumns ( diag, 3UL );
         checkCapacity( diag, 1UL );
         checkNonZeros( diag, 1UL );
         checkNonZeros( diag, 0UL, 0UL );
         checkNonZeros( diag, 1UL, 1UL );
         checkNonZeros( diag, 2UL, 0UL );

         if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) != 0 ||
             diag(1,0) != 0 || diag(1,1) != 1 || diag(1,2) != 0 ||
             diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 0 0 0 )\n( 0 1 0 )\n( 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Writing the element (2,2)
         diag(2,2) = 2;

         checkRows    ( diag, 3UL );
         checkColumns ( diag, 3UL );
         checkCapacity( diag, 2UL );
         checkNonZeros( diag, 2UL );
         checkNonZeros( diag, 0UL, 0UL );
         checkNonZeros( diag, 1UL, 1UL );
         checkNonZeros( diag, 2UL, 1UL );

         if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) != 0 ||
             diag(1,0) != 0 || diag(1,1) != 1 || diag(1,2) != 0 ||
             diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 0 0 0 )\n( 0 1 0 )\n( 0 0 2 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Adding to the element (0,0)
         diag(0,0) += 3;

         checkRows    ( diag, 3UL );
         checkColumns ( diag, 3UL );
         checkCapacity( diag, 3UL );
         checkNonZeros( diag, 3UL );
         checkNonZeros( diag, 0UL, 1UL );
         checkNonZeros( diag, 1UL, 1UL );
         checkNonZeros( diag, 2UL, 1UL );

         if( diag(0,0) != 3 || diag(0,1) != 0 || diag(0,2) != 0 ||
             diag(1,0) != 0 || diag(1,1) != 1 || diag(1,2) != 0 ||
             diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 3 0 0 )\n( 0 1 0 )\n( 0 0 2 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Subtracting from the element (1,1)
         diag(1,1) -= 4;

         checkRows    ( diag, 3UL );
         checkColumns ( diag, 3UL );
         checkCapacity( diag, 3UL );
         checkNonZeros( diag, 3UL );
         checkNonZeros( diag, 0UL, 1UL );
         checkNonZeros( diag, 1UL, 1UL );
         checkNonZeros( diag, 2UL, 1UL );

         if( diag(0,0) != 3 || diag(0,1) !=  0 || diag(0,2) != 0 ||
             diag(1,0) != 0 || diag(1,1) != -3 || diag(1,2) != 0 ||
             diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 3  0  0 )\n( 0 -3  0 )\n( 0  0  2 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Multiplying the element (2,2)
         diag(2,2) *= -3;

         checkRows    ( diag, 3UL );
         checkColumns ( diag, 3UL );
         checkCapacity( diag, 3UL );
         checkNonZeros( diag, 3UL );
         checkNonZeros( diag, 0UL, 1UL );
         checkNonZeros( diag, 1UL, 1UL );
         checkNonZeros( diag, 2UL, 1UL );

         if( diag(0,0) != 3 || diag(0,1) !=  0 || diag(0,2) !=  0 ||
             diag(1,0) != 0 || diag(1,1) != -3 || diag(1,2) !=  0 ||
             diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != -6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 3  0  0 )\n( 0 -3  0 )\n( 0  0 -6 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Dividing the element (2,2)
         diag(2,2) /= 2;

         checkRows    ( diag, 3UL );
         checkColumns ( diag, 3UL );
         checkCapacity( diag, 3UL );
         checkNonZeros( diag, 3UL );
         checkNonZeros( diag, 0UL, 1UL );
         checkNonZeros( diag, 1UL, 1UL );
         checkNonZeros( diag, 2UL, 1UL );

         if( diag(0,0) != 3 || diag(0,1) !=  0 || diag(0,2) !=  0 ||
             diag(1,0) != 0 || diag(1,1) != -3 || diag(1,2) !=  0 ||
             diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != -3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 3  0  0 )\n( 0 -3  0 )\n( 0  0 -3 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Failure cases
      {
         ODT diag( 3UL );

         // Trying to write the lower element (2,1)
         try {
            diag(2,1) = 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to lower matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to add to the lower element (2,1)
         try {
            diag(2,1) += 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Addition assignment to lower matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to subtract from the lower element (2,1)
         try {
            diag(2,1) -= 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subtraction assignment to lower matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to multiply the lower element (2,1)
         try {
            diag(2,1) *= 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Multiplication assignment to lower matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to divide the lower element (2,1)
         try {
            diag(2,1) /= 2;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Division assignment to lower matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to write the upper element (1,2)
         try {
            diag(1,2) = 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to upper matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to add to the upper element (1,2)
         try {
            diag(1,2) += 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Addition assignment to upper matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to subtract from the upper element (1,2)
         try {
            diag(1,2) -= 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subtraction assignment to upper matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to multiply the upper element (1,2)
         try {
            diag(1,2) *= 5;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Multiplication assignment to upper matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to divide the upper element (1,2)
         try {
            diag(1,2) /= 2;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Division assignment to upper matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the DiagonalMatrix iterator implementation.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the iterator implementation of the DiagonalMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testIterator()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      using Iterator      = DT::Iterator;
      using ConstIterator = DT::ConstIterator;

      DT diag( 3UL );
      diag(0,0) =  1;
      diag(1,1) = -2;
      diag(2,2) =  3;

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

         ConstIterator it( begin( diag, 1UL ) );

         if( it == end( diag, 1UL ) || it->value() != -2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Failed iterator conversion detected\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Counting the number of elements in 0th row via Iterator (end-begin)
      {
         test_ = "Row-major Iterator subtraction (end-begin)";

         const ptrdiff_t number( end( diag, 0UL ) - begin( diag, 0UL ) );

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

         const ptrdiff_t number( cend( diag, 1UL ) - cbegin( diag, 1UL ) );

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

         ConstIterator it ( cbegin( diag, 2UL ) );
         ConstIterator end( cend( diag, 2UL ) );

         if( it == end || it->value() != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid initial iterator detected\n";
            throw std::runtime_error( oss.str() );
         }

         ++it;

         if( it != end ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator pre-increment failed\n";
            throw std::runtime_error( oss.str() );
         }

         it = cbegin( diag, 2UL );
         ++it;

         if( it != end ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator post-increment failed\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing assignment to diagonal elements via Iterator
      {
         test_ = "Row-major assignment to diagonal elements via Iterator";

         const Iterator it = begin( diag, 0UL );
         *it = 4;

         if( diag(0,0) != 4 || diag(0,1) !=  0 || diag(0,2) != 0 ||
             diag(1,0) != 0 || diag(1,1) != -2 || diag(1,2) != 0 ||
             diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 4  0  0 )\n( 0 -2  0 )\n( 0  0  3 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing addition assignment to diagonal elements via Iterator
      {
         test_ = "Row-major addition assignment to diagonal elements via Iterator";

         const Iterator it = begin( diag, 1UL );
         *it += 3;

         if( diag(0,0) != 4 || diag(0,1) != 0 || diag(0,2) != 0 ||
             diag(1,0) != 0 || diag(1,1) != 1 || diag(1,2) != 0 ||
             diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 4 0 0 )\n( 0 1 0 )\n( 0 0 3 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing subtraction assignment to diagonal elements via Iterator
      {
         test_ = "Row-major subtraction assignment to diagonal elements via Iterator";

         const Iterator it = begin( diag, 2UL );
         *it -= 4;

         if( diag(0,0) != 4 || diag(0,1) != 0 || diag(0,2) !=  0 ||
             diag(1,0) != 0 || diag(1,1) != 1 || diag(1,2) !=  0 ||
             diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != -1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 4  0  0 )\n( 0  1  0 )\n( 0  0 -1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing multiplication assignment to diagonal elements via Iterator
      {
         test_ = "Row-major multiplication assignment to diagonal elements via Iterator";

         const Iterator it = begin( diag, 0UL );
         *it *= 2;

         if( diag(0,0) != 8 || diag(0,1) != 0 || diag(0,2) !=  0 ||
             diag(1,0) != 0 || diag(1,1) != 1 || diag(1,2) !=  0 ||
             diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != -1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 8  0  0 )\n( 0  1  0 )\n( 0  0 -1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing division assignment to diagonal elements via Iterator
      {
         test_ = "Row-major division assignment to diagonal elements via Iterator";

         const Iterator it = begin( diag, 0UL );
         *it /= 4;

         if( diag(0,0) != 2 || diag(0,1) != 0 || diag(0,2) !=  0 ||
             diag(1,0) != 0 || diag(1,1) != 1 || diag(1,2) !=  0 ||
             diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != -1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 2  0  0 )\n( 0  1  0 )\n( 0  0 -1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      using Iterator      = ODT::Iterator;
      using ConstIterator = ODT::ConstIterator;

      ODT diag( 3UL );
      diag(0,0) =  1;
      diag(1,1) = -2;
      diag(2,2) =  3;

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

         ConstIterator it( begin( diag, 1UL ) );

         if( it == end( diag, 1UL ) || it->value() != -2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Failed iterator conversion detected\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Counting the number of elements in 0th column via Iterator (end-begin)
      {
         test_ = "Column-major Iterator subtraction (end-begin)";

         const ptrdiff_t number( end( diag, 0UL ) - begin( diag, 0UL ) );

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

         const ptrdiff_t number( cend( diag, 1UL ) - cbegin( diag, 1UL ) );

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

         ConstIterator it ( cbegin( diag, 0UL ) );
         ConstIterator end( cend( diag, 0UL ) );

         if( it == end || it->value() != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid initial iterator detected\n";
            throw std::runtime_error( oss.str() );
         }

         ++it;

         if( it != end ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator pre-increment failed\n";
            throw std::runtime_error( oss.str() );
         }

         it = cbegin( diag, 0UL );
         it++;

         if( it != end ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator post-increment failed\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing assignment to diagonal elements via Iterator
      {
         test_ = "Column-major assignment to diagonal elements via Iterator";

         const Iterator it = begin( diag, 0UL );
         *it = 4;

         if( diag(0,0) != 4 || diag(0,1) !=  0 || diag(0,2) != 0 ||
             diag(1,0) != 0 || diag(1,1) != -2 || diag(1,2) != 0 ||
             diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 4  0  0 )\n( 0 -2  0 )\n( 0  0  3 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing addition assignment to diagonal elements via Iterator
      {
         test_ = "Column-major addition assignment to diagonal elements via Iterator";

         const Iterator it = begin( diag, 1UL );
         *it += 3;

         if( diag(0,0) != 4 || diag(0,1) != 0 || diag(0,2) != 0 ||
             diag(1,0) != 0 || diag(1,1) != 1 || diag(1,2) != 0 ||
             diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 4 0 0 )\n( 0 1 0 )\n( 0 0 3 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing subtraction assignment to diagonal elements via Iterator
      {
         test_ = "Column-major subtraction assignment to diagonal elements via Iterator";

         const Iterator it = begin( diag, 2UL );
         *it -= 4;

         if( diag(0,0) != 4 || diag(0,1) != 0 || diag(0,2) !=  0 ||
             diag(1,0) != 0 || diag(1,1) != 1 || diag(1,2) !=  0 ||
             diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != -1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 4  0  0 )\n( 0  1  0 )\n( 0  0 -1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing multiplication assignment to diagonal elements via Iterator
      {
         test_ = "Column-major multiplication assignment to diagonal elements via Iterator";

         const Iterator it = begin( diag, 0UL );
         *it *= 2;

         if( diag(0,0) != 8 || diag(0,1) != 0 || diag(0,2) !=  0 ||
             diag(1,0) != 0 || diag(1,1) != 1 || diag(1,2) !=  0 ||
             diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != -1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 8  0  0 )\n( 0  1  0 )\n( 0  0 -1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing division assignment to diagonal elements via Iterator
      {
         test_ = "Column-major division assignment to diagonal elements via Iterator";

         const Iterator it = begin( diag, 0UL );
         *it /= 4;

         if( diag(0,0) != 2 || diag(0,1) != 0 || diag(0,2) !=  0 ||
             diag(1,0) != 0 || diag(1,1) != 1 || diag(1,2) !=  0 ||
             diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != -1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 2  0  0 )\n( 0  1  0 )\n( 0  0 -1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c nonZeros() member function of the DiagonalMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c nonZeros() member function of the DiagonalMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testNonZeros()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major DiagonalMatrix::nonZeros()";

      // Empty matrix
      {
         DT diag( 3UL );

         checkRows    ( diag, 3UL );
         checkColumns ( diag, 3UL );
         checkNonZeros( diag, 0UL );
         checkNonZeros( diag, 0UL, 0UL );
         checkNonZeros( diag, 1UL, 0UL );
         checkNonZeros( diag, 2UL, 0UL );

         if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) != 0 ||
             diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
             diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 0 0 0 )\n( 0 0 0 )\n( 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Partially filled matrix
      {
         DT diag( 3UL );
         diag(0,0) =  1;
         diag(1,1) = -2;

         checkRows    ( diag, 3UL );
         checkColumns ( diag, 3UL );
         checkCapacity( diag, 2UL );
         checkNonZeros( diag, 2UL );
         checkNonZeros( diag, 0UL, 1UL );
         checkNonZeros( diag, 1UL, 1UL );
         checkNonZeros( diag, 2UL, 0UL );

         if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) != 0 ||
             diag(1,0) != 0 || diag(1,1) != -2 || diag(1,2) != 0 ||
             diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 1  0  0 )\n( 0 -2  0 )\n( 0  0  0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Fully filled matrix
      {
         DT diag( 3UL );
         diag(0,0) = -1;
         diag(1,1) =  2;
         diag(2,2) =  3;

         checkRows    ( diag, 3UL );
         checkColumns ( diag, 3UL );
         checkCapacity( diag, 3UL );
         checkNonZeros( diag, 3UL );
         checkNonZeros( diag, 0UL, 1UL );
         checkNonZeros( diag, 1UL, 1UL );
         checkNonZeros( diag, 2UL, 1UL );

         if( diag(0,0) != -1 || diag(0,1) != 0 || diag(0,2) != 0 ||
             diag(1,0) !=  0 || diag(1,1) != 2 || diag(1,2) != 0 ||
             diag(2,0) !=  0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( -1  0  0 )\n(  0  2  0 )\n(  0  0  3 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major DiagonalMatrix::nonZeros()";

      // Empty matrix
      {
         ODT diag( 3UL );

         checkRows    ( diag, 3UL );
         checkColumns ( diag, 3UL );
         checkNonZeros( diag, 0UL );
         checkNonZeros( diag, 0UL, 0UL );
         checkNonZeros( diag, 1UL, 0UL );
         checkNonZeros( diag, 2UL, 0UL );

         if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) != 0 ||
             diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
             diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 0 0 0 )\n( 0 0 0 )\n( 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Partially filled matrix
      {
         ODT diag( 3UL );
         diag(0,0) =  1;
         diag(1,1) = -2;

         checkRows    ( diag, 3UL );
         checkColumns ( diag, 3UL );
         checkCapacity( diag, 2UL );
         checkNonZeros( diag, 2UL );
         checkNonZeros( diag, 0UL, 1UL );
         checkNonZeros( diag, 1UL, 1UL );
         checkNonZeros( diag, 2UL, 0UL );

         if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) != 0 ||
             diag(1,0) != 0 || diag(1,1) != -2 || diag(1,2) != 0 ||
             diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 1  0  0 )\n( 0 -2  0 )\n( 0  0  0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Fully filled matrix
      {
         ODT diag( 3UL );
         diag(0,0) = -1;
         diag(1,1) =  2;
         diag(2,2) =  3;

         checkRows    ( diag, 3UL );
         checkColumns ( diag, 3UL );
         checkCapacity( diag, 3UL );
         checkNonZeros( diag, 3UL );
         checkNonZeros( diag, 0UL, 1UL );
         checkNonZeros( diag, 1UL, 1UL );
         checkNonZeros( diag, 2UL, 1UL );

         if( diag(0,0) != -1 || diag(0,1) != 0 || diag(0,2) != 0 ||
             diag(1,0) !=  0 || diag(1,1) != 2 || diag(1,2) != 0 ||
             diag(2,0) !=  0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( -1  0  0 )\n(  0  2  0 )\n(  0  0  3 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c reset() member function of the DiagonalMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c reset() member function of the DiagonalMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testReset()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major DiagonalMatrix::reset()";

      // Initialization check
      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 3UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting a diagonal element
      reset( diag(1,1) );

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 3UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting a lower element
      reset( diag(1,0) );

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 3UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting an upper element
      reset( diag(0,1) );

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 3UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting row 2
      reset( diag, 2UL );

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 3UL );
      checkNonZeros( diag, 1UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 0UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting the entire matrix
      reset( diag );

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 3UL );
      checkNonZeros( diag, 0UL );
      checkNonZeros( diag, 0UL, 0UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 0UL );

      if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 0 0 0 )\n( 0 0 0 )\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major DiagonalMatrix::reset()";

      // Initialization check
      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 3UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting a diagonal element
      reset( diag(1,1) );

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 3UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting a lower element
      reset( diag(1,0) );

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 3UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting an upper element
      reset( diag(0,1) );

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 3UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting row 2
      reset( diag, 2UL );

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 3UL );
      checkNonZeros( diag, 1UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 0UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting the entire matrix
      reset( diag );

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 3UL );
      checkNonZeros( diag, 0UL );
      checkNonZeros( diag, 0UL, 0UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 0UL );

      if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 0 0 0 )\n( 0 0 0 )\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c clear() member function of the DiagonalMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c clear() member function of the DiagonalMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testClear()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major DiagonalMatrix::clear()";

      // Initialization check
      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 3UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing a diagonal element
      clear( diag(1,1) );

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 3UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing a lower element
      clear( diag(1,0) );

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 3UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing an upper element
      clear( diag(0,1) );

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 3UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing the matrix
      clear( diag );

      checkRows    ( diag, 0UL );
      checkColumns ( diag, 0UL );
      checkNonZeros( diag, 0UL );
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major DiagonalMatrix::clear()";

      // Initialization check
      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 3UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing a diagonal element
      clear( diag(1,1) );

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 3UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing a lower element
      clear( diag(1,0) );

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 3UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing an upper element
      clear( diag(0,1) );

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 3UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing the matrix
      clear( diag );

      checkRows    ( diag, 0UL );
      checkColumns ( diag, 0UL );
      checkNonZeros( diag, 0UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c resize() member function of the DiagonalMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c resize() member function of the DiagonalMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testResize()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major DiagonalMatrix::resize()";

      // Initialization check
      DT diag;

      checkRows    ( diag, 0UL );
      checkColumns ( diag, 0UL );
      checkNonZeros( diag, 0UL );

      // Resizing to 2x2
      diag.resize( 2UL );

      checkRows    ( diag, 2UL );
      checkColumns ( diag, 2UL );
      checkNonZeros( diag, 0UL );
      checkNonZeros( diag, 0UL, 0UL );
      checkNonZeros( diag, 1UL, 0UL );

      // Resizing to 4x4 and preserving the elements
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag.resize( 4UL, true );

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkCapacity( diag, 2UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 0UL );
      checkNonZeros( diag, 3UL, 0UL );

      // Resizing to 2x2
      diag(2,2) = 3;
      diag.resize( 2UL );

      checkRows    ( diag, 2UL );
      checkColumns ( diag, 2UL );
      checkCapacity( diag, 2UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );

      // Resizing to 0x0
      diag.resize( 0UL );

      checkRows    ( diag, 0UL );
      checkColumns ( diag, 0UL );
      checkNonZeros( diag, 0UL );
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major DiagonalMatrix::resize()";

      // Initialization check
      ODT diag;

      checkRows    ( diag, 0UL );
      checkColumns ( diag, 0UL );
      checkNonZeros( diag, 0UL );

      // Resizing to 2x2
      diag.resize( 2UL );

      checkRows    ( diag, 2UL );
      checkColumns ( diag, 2UL );
      checkNonZeros( diag, 0UL );
      checkNonZeros( diag, 0UL, 0UL );
      checkNonZeros( diag, 1UL, 0UL );

      // Resizing to 4x4 and preserving the elements
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag.resize( 4UL, true );

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkCapacity( diag, 2UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 0UL );
      checkNonZeros( diag, 3UL, 0UL );

      // Resizing to 2x2
      diag(2,2) = 3;
      diag.resize( 2UL );

      checkRows    ( diag, 2UL );
      checkColumns ( diag, 2UL );
      checkCapacity( diag, 2UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );

      // Resizing to 0x0
      diag.resize( 0UL );

      checkRows    ( diag, 0UL );
      checkColumns ( diag, 0UL );
      checkNonZeros( diag, 0UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c reserve() member function of the DiagonalMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c reserve() member function of the DiagonalMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testReserve()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major DiagonalMatrix::reserve()";

      // Initialization check
      DT diag;

      checkRows    ( diag, 0UL );
      checkColumns ( diag, 0UL );
      checkNonZeros( diag, 0UL );

      // Increasing the capacity of the matrix
      diag.reserve( 10UL );

      checkRows    ( diag,  0UL );
      checkColumns ( diag,  0UL );
      checkCapacity( diag, 10UL );
      checkNonZeros( diag,  0UL );

      // Further increasing the capacity of the matrix
      diag.reserve( 20UL );

      checkRows    ( diag,  0UL );
      checkColumns ( diag,  0UL );
      checkCapacity( diag, 20UL );
      checkNonZeros( diag,  0UL );
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major DiagonalMatrix::reserve()";

      // Initialization check
      ODT diag;

      checkRows    ( diag, 0UL );
      checkColumns ( diag, 0UL );
      checkNonZeros( diag, 0UL );

      // Increasing the capacity of the matrix
      diag.reserve( 10UL );

      checkRows    ( diag,  0UL );
      checkColumns ( diag,  0UL );
      checkCapacity( diag, 10UL );
      checkNonZeros( diag,  0UL );

      // Further increasing the capacity of the matrix
      diag.reserve( 20UL );

      checkRows    ( diag,  0UL );
      checkColumns ( diag,  0UL );
      checkCapacity( diag, 20UL );
      checkNonZeros( diag,  0UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c trim() member function of the DiagonalMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c trim() member function of the DiagonalMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testTrim()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major DiagonalMatrix::trim()";

      // Initialization check
      DT diag( 3UL );

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkNonZeros( diag, 0UL );

      // Increasing the row capacity of the matrix
      diag.reserve( 0UL, 10UL );
      diag.reserve( 1UL, 15UL );
      diag.reserve( 2UL, 20UL );

      checkRows    ( diag,  3UL );
      checkColumns ( diag,  3UL );
      checkCapacity( diag, 45UL );
      checkCapacity( diag,  0UL, 10UL );
      checkCapacity( diag,  1UL, 15UL );
      checkCapacity( diag,  2UL, 20UL );

      // Trimming the matrix
      diag.trim();

      checkRows    ( diag,  3UL );
      checkColumns ( diag,  3UL );
      checkCapacity( diag, 45UL );
      checkCapacity( diag,  0UL, 0UL );
      checkCapacity( diag,  1UL, 0UL );
      checkCapacity( diag,  2UL, 0UL );
   }

   {
      test_ = "Row-major DiagonalMatrix::trim( size_t )";

      // Initialization check
      DT diag( 3UL, 3UL );

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkNonZeros( diag, 0UL );

      // Increasing the row capacity of the matrix
      diag.reserve( 0UL, 10UL );
      diag.reserve( 1UL, 15UL );
      diag.reserve( 2UL, 20UL );

      checkRows    ( diag,  3UL );
      checkColumns ( diag,  3UL );
      checkCapacity( diag, 45UL );
      checkCapacity( diag,  0UL, 10UL );
      checkCapacity( diag,  1UL, 15UL );
      checkCapacity( diag,  2UL, 20UL );

      // Trimming the 0th row
      diag.trim( 0UL );

      checkRows    ( diag,  3UL );
      checkColumns ( diag,  3UL );
      checkCapacity( diag, 45UL );
      checkCapacity( diag,  0UL,  0UL );
      checkCapacity( diag,  1UL, 25UL );
      checkCapacity( diag,  2UL, 20UL );

      // Trimming the 1st row
      diag.trim( 1UL );

      checkRows    ( diag,  3UL );
      checkColumns ( diag,  3UL );
      checkCapacity( diag, 45UL );
      checkCapacity( diag,  0UL,  0UL );
      checkCapacity( diag,  1UL,  0UL );
      checkCapacity( diag,  2UL, 45UL );

      // Trimming the 2nd row
      diag.trim( 2UL );

      checkRows    ( diag,  3UL );
      checkColumns ( diag,  3UL );
      checkCapacity( diag, 45UL );
      checkCapacity( diag,  0UL, 0UL );
      checkCapacity( diag,  1UL, 0UL );
      checkCapacity( diag,  2UL, 0UL );
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major DiagonalMatrix::trim()";

      // Initialization check
      ODT diag( 3UL );

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkNonZeros( diag, 0UL );

      // Increasing the row capacity of the matrix
      diag.reserve( 0UL, 10UL );
      diag.reserve( 1UL, 15UL );
      diag.reserve( 2UL, 20UL );

      checkRows    ( diag,  3UL );
      checkColumns ( diag,  3UL );
      checkCapacity( diag, 45UL );
      checkCapacity( diag,  0UL, 10UL );
      checkCapacity( diag,  1UL, 15UL );
      checkCapacity( diag,  2UL, 20UL );

      // Trimming the matrix
      diag.trim();

      checkRows    ( diag,  3UL );
      checkColumns ( diag,  3UL );
      checkCapacity( diag, 45UL );
      checkCapacity( diag,  0UL, 0UL );
      checkCapacity( diag,  1UL, 0UL );
      checkCapacity( diag,  2UL, 0UL );
   }

   {
      test_ = "Column-major DiagonalMatrix::trim( size_t )";

      // Initialization check
      ODT diag( 3UL, 3UL );

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkNonZeros( diag, 0UL );

      // Increasing the column capacity of the matrix
      diag.reserve( 0UL, 10UL );
      diag.reserve( 1UL, 15UL );
      diag.reserve( 2UL, 20UL );

      checkRows    ( diag,  3UL );
      checkColumns ( diag,  3UL );
      checkCapacity( diag, 45UL );
      checkCapacity( diag,  0UL, 10UL );
      checkCapacity( diag,  1UL, 15UL );
      checkCapacity( diag,  2UL, 20UL );

      // Trimming the 0th column
      diag.trim( 0UL );

      checkRows    ( diag,  3UL );
      checkColumns ( diag,  3UL );
      checkCapacity( diag, 45UL );
      checkCapacity( diag,  0UL,  0UL );
      checkCapacity( diag,  1UL, 25UL );
      checkCapacity( diag,  2UL, 20UL );

      // Trimming the 1st column
      diag.trim( 1UL );

      checkRows    ( diag,  3UL );
      checkColumns ( diag,  3UL );
      checkCapacity( diag, 45UL );
      checkCapacity( diag,  0UL,  0UL );
      checkCapacity( diag,  1UL,  0UL );
      checkCapacity( diag,  2UL, 45UL );

      // Trimming the 2nd column
      diag.trim( 2UL );

      checkRows    ( diag,  3UL );
      checkColumns ( diag,  3UL );
      checkCapacity( diag, 45UL );
      checkCapacity( diag,  0UL, 0UL );
      checkCapacity( diag,  1UL, 0UL );
      checkCapacity( diag,  2UL, 0UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c shrinkToFit() member function of the DiagonalMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c shrinkToFit() member function of the DiagonalMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testShrinkToFit()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major DiagonalMatrix::shrinkToFit()";

      // Shrinking a matrix without excessive capacity
      {
         DT diag( 3UL, 3UL );
         diag(0,0) = 1;
         diag(1,1) = 2;
         diag(2,2) = 3;

         diag.shrinkToFit();

         checkRows    ( diag, 3UL );
         checkColumns ( diag, 3UL );
         checkCapacity( diag, 3UL );
         checkNonZeros( diag, 3UL );
         checkNonZeros( diag, 0UL, 1UL );
         checkNonZeros( diag, 1UL, 1UL );
         checkNonZeros( diag, 2UL, 1UL );

         if( diag.capacity() != diag.nonZeros() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Shrinking the matrix failed\n"
                << " Details:\n"
                << "   Capacity         : " << diag.capacity() << "\n"
                << "   Expected capacity: " << diag.nonZeros() << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
             diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 ||
             diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Shrinking the matrix failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Shrinking a matrix with excessive capacity
      {
         DT diag( 3UL, 100UL );
         diag(0,0) = 1;
         diag(1,1) = 2;
         diag(2,2) = 3;

         diag.shrinkToFit();

         checkRows    ( diag, 3UL );
         checkColumns ( diag, 3UL );
         checkCapacity( diag, 3UL );
         checkNonZeros( diag, 3UL );
         checkNonZeros( diag, 0UL, 1UL );
         checkNonZeros( diag, 1UL, 1UL );
         checkNonZeros( diag, 2UL, 1UL );

         if( diag.capacity() != diag.nonZeros() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Shrinking the matrix failed\n"
                << " Details:\n"
                << "   Capacity         : " << diag.capacity() << "\n"
                << "   Expected capacity: " << diag.nonZeros() << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
             diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 ||
             diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Shrinking the matrix failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major DiagonalMatrix::shrinkToFit()";

      // Shrinking a matrix without excessive capacity
      {
         ODT diag( 3UL, 3UL );
         diag(0,0) = 1;
         diag(1,1) = 2;
         diag(2,2) = 3;

         diag.shrinkToFit();

         checkRows    ( diag, 3UL );
         checkColumns ( diag, 3UL );
         checkCapacity( diag, 3UL );
         checkNonZeros( diag, 3UL );
         checkNonZeros( diag, 0UL, 1UL );
         checkNonZeros( diag, 1UL, 1UL );
         checkNonZeros( diag, 2UL, 1UL );

         if( diag.capacity() != diag.nonZeros() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Shrinking the matrix failed\n"
                << " Details:\n"
                << "   Capacity         : " << diag.capacity() << "\n"
                << "   Expected capacity: " << diag.nonZeros() << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
             diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 ||
             diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Shrinking the matrix failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Shrinking a matrix with excessive capacity
      {
         ODT diag( 3UL, 100UL );
         diag(0,0) = 1;
         diag(1,1) = 2;
         diag(2,2) = 3;

         diag.shrinkToFit();

         checkRows    ( diag, 3UL );
         checkColumns ( diag, 3UL );
         checkCapacity( diag, 3UL );
         checkNonZeros( diag, 3UL );
         checkNonZeros( diag, 0UL, 1UL );
         checkNonZeros( diag, 1UL, 1UL );
         checkNonZeros( diag, 2UL, 1UL );

         if( diag.capacity() != diag.nonZeros() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Shrinking the matrix failed\n"
                << " Details:\n"
                << "   Capacity         : " << diag.capacity() << "\n"
                << "   Expected capacity: " << diag.nonZeros() << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
             diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 ||
             diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Shrinking the matrix failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c swap() functionality of the DiagonalMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c swap() function of the DiagonalMatrix specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testSwap()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major DiagonalMatrix swap";

      DT diag1( 2UL );
      diag1(0,0) = 1;
      diag1(1,1) = 2;

      DT diag2( 3UL );
      diag2(0,0) = 3;
      diag2(1,1) = 4;
      diag2(2,2) = 5;

      swap( diag1, diag2 );

      checkRows    ( diag1, 3UL );
      checkColumns ( diag1, 3UL );
      checkCapacity( diag1, 3UL );
      checkNonZeros( diag1, 3UL );
      checkNonZeros( diag1, 0UL, 1UL );
      checkNonZeros( diag1, 1UL, 1UL );
      checkNonZeros( diag1, 2UL, 1UL );

      if( diag1(0,0) != 3 || diag1(0,1) != 0 || diag1(0,2) != 0 ||
          diag1(1,0) != 0 || diag1(1,1) != 4 || diag1(1,2) != 0 ||
          diag1(2,0) != 0 || diag1(2,1) != 0 || diag1(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the first matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag1 << "\n"
             << "   Expected result:\n( 3 0 0 )\n( 0 4 0 )\n( 0 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }

      checkRows    ( diag2, 2UL );
      checkColumns ( diag2, 2UL );
      checkCapacity( diag2, 2UL );
      checkNonZeros( diag2, 2UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(1,0) != 0 || diag2(1,1) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the second matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 )\n( 0 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major DiagonalMatrix swap";

      ODT diag1( 2UL );
      diag1(0,0) = 1;
      diag1(1,1) = 2;

      ODT diag2( 3UL );
      diag2(0,0) = 3;
      diag2(1,1) = 4;
      diag2(2,2) = 5;

      swap( diag1, diag2 );

      checkRows    ( diag1, 3UL );
      checkColumns ( diag1, 3UL );
      checkCapacity( diag1, 3UL );
      checkNonZeros( diag1, 3UL );
      checkNonZeros( diag1, 0UL, 1UL );
      checkNonZeros( diag1, 1UL, 1UL );
      checkNonZeros( diag1, 2UL, 1UL );

      if( diag1(0,0) != 3 || diag1(0,1) != 0 || diag1(0,2) != 0 ||
          diag1(1,0) != 0 || diag1(1,1) != 4 || diag1(1,2) != 0 ||
          diag1(2,0) != 0 || diag1(2,1) != 0 || diag1(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the first matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag1 << "\n"
             << "   Expected result:\n( 3 0 0 )\n( 0 4 0 )\n( 0 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }

      checkRows    ( diag2, 2UL );
      checkColumns ( diag2, 2UL );
      checkCapacity( diag2, 2UL );
      checkNonZeros( diag2, 2UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(1,0) != 0 || diag2(1,1) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the second matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 )\n( 0 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c set() member function of the DiagonalMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c set() member function of the DiagonalMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testSet()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major DiagonalMatrix::set()";

      using Iterator = DT::Iterator;

      // Initialization check
      DT diag( 4UL );

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkNonZeros( diag, 0UL );
      checkNonZeros( diag, 0UL, 0UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 0UL );
      checkNonZeros( diag, 3UL, 0UL );

      // Setting a non-zero element
      {
         Iterator pos = diag.set( 2UL, 2UL, 1 );

         checkRows    ( diag, 4UL );
         checkColumns ( diag, 4UL );
         checkCapacity( diag, 1UL );
         checkNonZeros( diag, 1UL );
         checkNonZeros( diag, 0UL, 0UL );
         checkNonZeros( diag, 1UL, 0UL );
         checkNonZeros( diag, 2UL, 1UL );
         checkNonZeros( diag, 3UL, 0UL );

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

         if( diag(2,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 0 1 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setting a second non-zero element
      {
         Iterator pos = diag.set( 3UL, 3UL, 2 );

         checkRows    ( diag, 4UL );
         checkColumns ( diag, 4UL );
         checkCapacity( diag, 2UL );
         checkNonZeros( diag, 2UL );
         checkNonZeros( diag, 0UL, 0UL );
         checkNonZeros( diag, 1UL, 0UL );
         checkNonZeros( diag, 2UL, 1UL );
         checkNonZeros( diag, 3UL, 1UL );

         if( pos->value() != 2 || pos->index() != 3UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 2\n"
                << "   Expected index: 3\n";
            throw std::runtime_error( oss.str() );
         }

         if( diag(2,2) != 1 || diag(3,3) != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 0 1 0 )\n( 0 0 0 2 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setting a third non-zero element
      {
         Iterator pos = diag.set( 0UL, 0UL, 3 );

         checkRows    ( diag, 4UL );
         checkColumns ( diag, 4UL );
         checkCapacity( diag, 3UL );
         checkNonZeros( diag, 3UL );
         checkNonZeros( diag, 0UL, 1UL );
         checkNonZeros( diag, 1UL, 0UL );
         checkNonZeros( diag, 2UL, 1UL );
         checkNonZeros( diag, 3UL, 1UL );

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

         if( diag(0,0) != 3 || diag(2,2) != 1 || diag(3,3) != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 3 0 0 0 )\n( 0 0 0 0 )\n( 0 0 1 0 )\n( 0 0 0 2 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setting an already existing element
      {
         Iterator pos = diag.set( 2UL, 2UL, 4 );

         checkRows    ( diag, 4UL );
         checkColumns ( diag, 4UL );
         checkCapacity( diag, 3UL );
         checkNonZeros( diag, 3UL );
         checkNonZeros( diag, 0UL, 1UL );
         checkNonZeros( diag, 1UL, 0UL );
         checkNonZeros( diag, 2UL, 1UL );
         checkNonZeros( diag, 3UL, 1UL );

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

         if( diag(0,0) != 3 || diag(2,2) != 4 || diag(3,3) != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 3 0 0 0 )\n( 0 0 0 0 )\n( 0 0 4 0 )\n( 0 0 0 2 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major DiagonalMatrix::set()";

      using Iterator = ODT::Iterator;

      // Initialization check
      ODT diag( 4UL );

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkNonZeros( diag, 0UL );
      checkNonZeros( diag, 0UL, 0UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 0UL );
      checkNonZeros( diag, 3UL, 0UL );

      // Setting a non-zero element
      {
         Iterator pos = diag.set( 2UL, 2UL, 1 );

         checkRows    ( diag, 4UL );
         checkColumns ( diag, 4UL );
         checkCapacity( diag, 1UL );
         checkNonZeros( diag, 1UL );
         checkNonZeros( diag, 0UL, 0UL );
         checkNonZeros( diag, 1UL, 0UL );
         checkNonZeros( diag, 2UL, 1UL );
         checkNonZeros( diag, 3UL, 0UL );

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

         if( diag(2,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 0 1 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setting a second non-zero element
      {
         Iterator pos = diag.set( 3UL, 3UL, 2 );

         checkRows    ( diag, 4UL );
         checkColumns ( diag, 4UL );
         checkCapacity( diag, 2UL );
         checkNonZeros( diag, 2UL );
         checkNonZeros( diag, 0UL, 0UL );
         checkNonZeros( diag, 1UL, 0UL );
         checkNonZeros( diag, 2UL, 1UL );
         checkNonZeros( diag, 3UL, 1UL );

         if( pos->value() != 2 || pos->index() != 3UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 2\n"
                << "   Expected index: 3\n";
            throw std::runtime_error( oss.str() );
         }

         if( diag(2,2) != 1 || diag(3,3) != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 0 1 0 )\n( 0 0 0 2 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setting a third non-zero element
      {
         Iterator pos = diag.set( 0UL, 0UL, 3 );

         checkRows    ( diag, 4UL );
         checkColumns ( diag, 4UL );
         checkCapacity( diag, 3UL );
         checkNonZeros( diag, 3UL );
         checkNonZeros( diag, 0UL, 1UL );
         checkNonZeros( diag, 1UL, 0UL );
         checkNonZeros( diag, 2UL, 1UL );
         checkNonZeros( diag, 3UL, 1UL );

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

         if( diag(0,0) != 3 || diag(2,2) != 1 || diag(3,3) != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 3 0 0 0 )\n( 0 0 0 0 )\n( 0 0 1 0 )\n( 0 0 0 2 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setting an already existing element
      {
         Iterator pos = diag.set( 2UL, 2UL, 4 );

         checkRows    ( diag, 4UL );
         checkColumns ( diag, 4UL );
         checkCapacity( diag, 3UL );
         checkNonZeros( diag, 3UL );
         checkNonZeros( diag, 0UL, 1UL );
         checkNonZeros( diag, 1UL, 0UL );
         checkNonZeros( diag, 2UL, 1UL );
         checkNonZeros( diag, 3UL, 1UL );

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

         if( diag(0,0) != 3 || diag(2,2) != 4 || diag(3,3) != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 3 0 0 0 )\n( 0 0 0 0 )\n( 0 0 4 0 )\n( 0 0 0 2 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c insert() member function of the DiagonalMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c insert() member function of the DiagonalMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testInsert()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major DiagonalMatrix::insert()";

      using Iterator = DT::Iterator;

      // Initialization check
      DT diag( 4UL );

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkNonZeros( diag, 0UL );
      checkNonZeros( diag, 0UL, 0UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 0UL );
      checkNonZeros( diag, 3UL, 0UL );

      // Inserting a non-zero element
      {
         Iterator pos = diag.insert( 2UL, 2UL, 1 );

         checkRows    ( diag, 4UL );
         checkColumns ( diag, 4UL );
         checkCapacity( diag, 1UL );
         checkNonZeros( diag, 1UL );
         checkNonZeros( diag, 0UL, 0UL );
         checkNonZeros( diag, 1UL, 0UL );
         checkNonZeros( diag, 2UL, 1UL );
         checkNonZeros( diag, 3UL, 0UL );

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

         if( diag(2,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 0 1 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting a second non-zero element
      {
         Iterator pos = diag.insert( 3UL, 3UL, 2 );

         checkRows    ( diag, 4UL );
         checkColumns ( diag, 4UL );
         checkCapacity( diag, 2UL );
         checkNonZeros( diag, 2UL );
         checkNonZeros( diag, 0UL, 0UL );
         checkNonZeros( diag, 1UL, 0UL );
         checkNonZeros( diag, 2UL, 1UL );
         checkNonZeros( diag, 3UL, 1UL );

         if( pos->value() != 2 || pos->index() != 3UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 2\n"
                << "   Expected index: 3\n";
            throw std::runtime_error( oss.str() );
         }

         if( diag(2,2) != 1 || diag(3,3) != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 0 1 0 )\n( 0 0 0 2 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting a third non-zero element
      {
         Iterator pos = diag.insert( 0UL, 0UL, 3 );

         checkRows    ( diag, 4UL );
         checkColumns ( diag, 4UL );
         checkCapacity( diag, 3UL );
         checkNonZeros( diag, 3UL );
         checkNonZeros( diag, 0UL, 1UL );
         checkNonZeros( diag, 1UL, 0UL );
         checkNonZeros( diag, 2UL, 1UL );
         checkNonZeros( diag, 3UL, 1UL );

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

         if( diag(0,0) != 3 || diag(2,2) != 1 || diag(3,3) != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 3 0 0 0 )\n( 0 0 0 0 )\n( 0 0 1 0 )\n( 0 0 0 2 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to insert an already existing element
      try {
         diag.insert( 2UL, 2UL, 4 );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Inserting an existing element succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 3 0 0 0 )\n( 0 0 0 0 )\n( 0 0 1 0 )\n( 0 0 0 2 )\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major DiagonalMatrix::insert()";

      using Iterator = ODT::Iterator;

      // Initialization check
      ODT diag( 4UL );

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkNonZeros( diag, 0UL );
      checkNonZeros( diag, 0UL, 0UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 0UL );
      checkNonZeros( diag, 3UL, 0UL );

      // Inserting a non-zero element
      {
         Iterator pos = diag.insert( 2UL, 2UL, 1 );

         checkRows    ( diag, 4UL );
         checkColumns ( diag, 4UL );
         checkCapacity( diag, 1UL );
         checkNonZeros( diag, 1UL );
         checkNonZeros( diag, 0UL, 0UL );
         checkNonZeros( diag, 1UL, 0UL );
         checkNonZeros( diag, 2UL, 1UL );
         checkNonZeros( diag, 3UL, 0UL );

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

         if( diag(2,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 0 1 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting a second non-zero element
      {
         Iterator pos = diag.insert( 3UL, 3UL, 2 );

         checkRows    ( diag, 4UL );
         checkColumns ( diag, 4UL );
         checkCapacity( diag, 2UL );
         checkNonZeros( diag, 2UL );
         checkNonZeros( diag, 0UL, 0UL );
         checkNonZeros( diag, 1UL, 0UL );
         checkNonZeros( diag, 2UL, 1UL );
         checkNonZeros( diag, 3UL, 1UL );

         if( pos->value() != 2 || pos->index() != 3UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 2\n"
                << "   Expected index: 3\n";
            throw std::runtime_error( oss.str() );
         }

         if( diag(2,2) != 1 || diag(3,3) != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 0 1 0 )\n( 0 0 0 2 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting a third non-zero element
      {
         Iterator pos = diag.insert( 0UL, 0UL, 3 );

         checkRows    ( diag, 4UL );
         checkColumns ( diag, 4UL );
         checkCapacity( diag, 3UL );
         checkNonZeros( diag, 3UL );
         checkNonZeros( diag, 0UL, 1UL );
         checkNonZeros( diag, 1UL, 0UL );
         checkNonZeros( diag, 2UL, 1UL );
         checkNonZeros( diag, 3UL, 1UL );

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

         if( diag(0,0) != 3 || diag(2,2) != 1 || diag(3,3) != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 3 0 0 0 )\n( 0 0 0 0 )\n( 0 0 1 0 )\n( 0 0 0 2 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to insert an already existing element
      try {
         diag.insert( 2UL, 2UL, 4 );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Inserting an existing element succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 3 0 0 0 )\n( 0 0 0 0 )\n( 0 0 1 0 )\n( 0 0 0 2 )\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c append() member function of the DiagonalMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c append() member function of the DiagonalMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testAppend()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major DiagonalMatrix::append()";

      // Appending with pre-allocation in each row
      {
         // Initialization check
         DT diag( 4UL, 3UL );
         diag.reserve( 0UL, 1UL );
         diag.reserve( 2UL, 1UL );
         diag.reserve( 3UL, 1UL );

         checkRows    ( diag, 4UL );
         checkColumns ( diag, 4UL );
         checkNonZeros( diag, 0UL );
         checkNonZeros( diag, 0UL, 0UL );
         checkNonZeros( diag, 1UL, 0UL );
         checkNonZeros( diag, 2UL, 0UL );
         checkNonZeros( diag, 3UL, 0UL );

         // Appending one non-zero element
         diag.append( 2UL, 2UL, 1 );

         checkRows    ( diag, 4UL );
         checkColumns ( diag, 4UL );
         checkCapacity( diag, 1UL );
         checkNonZeros( diag, 1UL );
         checkNonZeros( diag, 0UL, 0UL );
         checkNonZeros( diag, 1UL, 0UL );
         checkNonZeros( diag, 2UL, 1UL );
         checkNonZeros( diag, 3UL, 0UL );

         if( diag(2,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 0 2 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Appending two more non-zero elements
         diag.append( 3UL, 3UL, 2 );
         diag.append( 0UL, 0UL, 3 );

         checkRows    ( diag, 4UL );
         checkColumns ( diag, 4UL );
         checkCapacity( diag, 3UL );
         checkNonZeros( diag, 3UL );
         checkNonZeros( diag, 0UL, 1UL );
         checkNonZeros( diag, 1UL, 0UL );
         checkNonZeros( diag, 2UL, 1UL );
         checkNonZeros( diag, 3UL, 1UL );

         if( diag(0,0) != 3 || diag(2,2) != 1 || diag(3,3) != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 3 0 0 0 )\n( 0 0 0 0 )\n( 0 0 1 0 )\n( 0 0 0 2 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Appending with row finalization
      {
         // Initialization check
         DT diag( 4UL, 3UL );

         // Appending one non-zero element
         diag.append( 0UL, 0UL, 1 );
         diag.finalize( 0UL );

         checkRows    ( diag, 4UL );
         checkColumns ( diag, 4UL );
         checkCapacity( diag, 1UL );
         checkNonZeros( diag, 1UL );
         checkNonZeros( diag, 0UL, 1UL );
         checkNonZeros( diag, 1UL, 0UL );
         checkNonZeros( diag, 2UL, 0UL );
         checkNonZeros( diag, 3UL, 0UL );

         if( diag(0,0) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Appending two more non-zero elements
         diag.finalize( 1UL );
         diag.append( 2UL, 2UL, 2 );
         diag.finalize( 2UL );
         diag.append( 3UL, 3UL, 3 );
         diag.finalize( 3UL );

         checkRows    ( diag, 4UL );
         checkColumns ( diag, 4UL );
         checkCapacity( diag, 3UL );
         checkNonZeros( diag, 3UL );
         checkNonZeros( diag, 0UL, 1UL );
         checkNonZeros( diag, 1UL, 0UL );
         checkNonZeros( diag, 2UL, 1UL );
         checkNonZeros( diag, 3UL, 1UL );

         if( diag(0,0) != 1 || diag(2,2) != 2 || diag(3,3) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 0 0 0 0 )\n( 0 0 2 0 )\n( 0 0 0 3 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major DiagonalMatrix::append()";

      // Appending with pre-allocation in each row
      {
         // Initialization check
         ODT diag( 4UL, 3UL );
         diag.reserve( 0UL, 1UL );
         diag.reserve( 2UL, 1UL );
         diag.reserve( 3UL, 1UL );

         checkRows    ( diag, 4UL );
         checkColumns ( diag, 4UL );
         checkNonZeros( diag, 0UL );
         checkNonZeros( diag, 0UL, 0UL );
         checkNonZeros( diag, 1UL, 0UL );
         checkNonZeros( diag, 2UL, 0UL );
         checkNonZeros( diag, 3UL, 0UL );

         // Appending one non-zero element
         diag.append( 2UL, 2UL, 1 );

         checkRows    ( diag, 4UL );
         checkColumns ( diag, 4UL );
         checkCapacity( diag, 1UL );
         checkNonZeros( diag, 1UL );
         checkNonZeros( diag, 0UL, 0UL );
         checkNonZeros( diag, 1UL, 0UL );
         checkNonZeros( diag, 2UL, 1UL );
         checkNonZeros( diag, 3UL, 0UL );

         if( diag(2,2) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 0 2 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Appending two more non-zero elements
         diag.append( 3UL, 3UL, 2 );
         diag.append( 0UL, 0UL, 3 );

         checkRows    ( diag, 4UL );
         checkColumns ( diag, 4UL );
         checkCapacity( diag, 3UL );
         checkNonZeros( diag, 3UL );
         checkNonZeros( diag, 0UL, 1UL );
         checkNonZeros( diag, 1UL, 0UL );
         checkNonZeros( diag, 2UL, 1UL );
         checkNonZeros( diag, 3UL, 1UL );

         if( diag(0,0) != 3 || diag(2,2) != 1 || diag(3,3) != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 3 0 0 0 )\n( 0 0 0 0 )\n( 0 0 1 0 )\n( 0 0 0 2 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Appending with row finalization
      {
         // Initialization check
         DT diag( 4UL, 3UL );

         // Appending one non-zero element
         diag.append( 0UL, 0UL, 1 );
         diag.finalize( 0UL );

         checkRows    ( diag, 4UL );
         checkColumns ( diag, 4UL );
         checkCapacity( diag, 1UL );
         checkNonZeros( diag, 1UL );
         checkNonZeros( diag, 0UL, 1UL );
         checkNonZeros( diag, 1UL, 0UL );
         checkNonZeros( diag, 2UL, 0UL );
         checkNonZeros( diag, 3UL, 0UL );

         if( diag(0,0) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Appending two more non-zero elements
         diag.finalize( 1UL );
         diag.append( 2UL, 2UL, 2 );
         diag.finalize( 2UL );
         diag.append( 3UL, 3UL, 3 );
         diag.finalize( 3UL );

         checkRows    ( diag, 4UL );
         checkColumns ( diag, 4UL );
         checkCapacity( diag, 3UL );
         checkNonZeros( diag, 3UL );
         checkNonZeros( diag, 0UL, 1UL );
         checkNonZeros( diag, 1UL, 0UL );
         checkNonZeros( diag, 2UL, 1UL );
         checkNonZeros( diag, 3UL, 1UL );

         if( diag(0,0) != 1 || diag(2,2) != 2 || diag(3,3) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 0 0 0 0 )\n( 0 0 2 0 )\n( 0 0 0 3 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c erase() member function of the DiagonalMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c erase() member function of the DiagonalMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testErase()
{
   //=====================================================================================
   // Row-major index-based erase function
   //=====================================================================================

   {
      test_ = "Row-major DiagonalMatrix::erase( size_t, size_t )";

      // Initialization check
      DT diag( 4UL, 4UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;
      diag(3,3) = 4;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkCapacity( diag, 4UL );
      checkNonZeros( diag, 4UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );
      checkNonZeros( diag, 3UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) != 0 || diag(3,2) != 0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 0 2 0 0 )\n( 0 0 3 0 )\n( 0 0 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element at (1,1)
      diag.erase( 1UL, 1UL );

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkCapacity( diag, 4UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );
      checkNonZeros( diag, 3UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) != 0 || diag(3,2) != 0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 0 0 0 0 )\n( 0 0 3 0 )\n( 0 0 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element at (3,3)
      diag.erase( 3UL, 3UL );

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkCapacity( diag, 4UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );
      checkNonZeros( diag, 3UL, 0UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) != 0 || diag(3,2) != 0 || diag(3,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 0 0 0 0 )\n( 0 0 3 0 )\n( 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element at (0,0)
      diag.erase( 0UL, size_t(0) );

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkCapacity( diag, 4UL );
      checkNonZeros( diag, 1UL );
      checkNonZeros( diag, 0UL, 0UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );
      checkNonZeros( diag, 3UL, 0UL );

      if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) != 0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) != 0 || diag(3,2) != 0 || diag(3,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 0 3 0 )\n( 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to erase a zero diagonal element
      diag.erase( 1UL, 1UL );

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkCapacity( diag, 4UL );
      checkNonZeros( diag, 1UL );
      checkNonZeros( diag, 0UL, 0UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );
      checkNonZeros( diag, 3UL, 0UL );

      if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) != 0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) != 0 || diag(3,2) != 0 || diag(3,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 0 3 0 )\n( 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to erase a zero lower element
      diag.erase( 2UL, 1UL );

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkCapacity( diag, 4UL );
      checkNonZeros( diag, 1UL );
      checkNonZeros( diag, 0UL, 0UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );
      checkNonZeros( diag, 3UL, 0UL );

      if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) != 0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) != 0 || diag(3,2) != 0 || diag(3,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 0 3 0 )\n( 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to erase a zero upper element
      diag.erase( 1UL, 2UL );

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkCapacity( diag, 4UL );
      checkNonZeros( diag, 1UL );
      checkNonZeros( diag, 0UL, 0UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );
      checkNonZeros( diag, 3UL, 0UL );

      if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) != 0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) != 0 || diag(3,2) != 0 || diag(3,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 0 3 0 )\n( 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major iterator-based erase function
   //=====================================================================================

   {
      test_ = "Row-major DiagonalMatrix::erase( size_t, Iterator )";

      using Iterator = DT::Iterator;

      // Initialization check
      DT diag( 4UL, 4UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;
      diag(3,3) = 4;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkCapacity( diag, 4UL );
      checkNonZeros( diag, 4UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );
      checkNonZeros( diag, 3UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) != 0 || diag(3,2) != 0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 0 2 0 0 )\n( 0 0 3 0 )\n( 0 0 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element at (1,1)
      {
         Iterator pos = diag.erase( 1UL, diag.find( 1UL, 1UL ) );

         checkRows    ( diag, 4UL );
         checkColumns ( diag, 4UL );
         checkCapacity( diag, 4UL );
         checkNonZeros( diag, 3UL );
         checkNonZeros( diag, 0UL, 1UL );
         checkNonZeros( diag, 1UL, 0UL );
         checkNonZeros( diag, 2UL, 1UL );
         checkNonZeros( diag, 3UL, 1UL );

         if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 || diag(0,3) != 0 ||
             diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 || diag(1,3) != 0 ||
             diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 || diag(2,3) != 0 ||
             diag(3,0) != 0 || diag(3,1) != 0 || diag(3,2) != 0 || diag(3,3) != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 0 0 0 0 )\n( 0 0 3 0 )\n( 0 0 0 4 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos != diag.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the element at (3,3)
      {
         Iterator pos = diag.erase( 3UL, diag.find( 3UL, 3UL ) );

         checkRows    ( diag, 4UL );
         checkColumns ( diag, 4UL );
         checkCapacity( diag, 4UL );
         checkNonZeros( diag, 2UL );
         checkNonZeros( diag, 0UL, 1UL );
         checkNonZeros( diag, 1UL, 0UL );
         checkNonZeros( diag, 2UL, 1UL );
         checkNonZeros( diag, 3UL, 0UL );

         if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 || diag(0,3) != 0 ||
             diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 || diag(1,3) != 0 ||
             diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 || diag(2,3) != 0 ||
             diag(3,0) != 0 || diag(3,1) != 0 || diag(3,2) != 0 || diag(3,3) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 0 0 0 0 )\n( 0 0 3 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos != diag.end( 3UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the element at (0,0)
      {
         Iterator pos = diag.erase( 0UL, diag.find( 0UL, 0UL ) );

         checkRows    ( diag, 4UL );
         checkColumns ( diag, 4UL );
         checkCapacity( diag, 4UL );
         checkNonZeros( diag, 1UL );
         checkNonZeros( diag, 0UL, 0UL );
         checkNonZeros( diag, 1UL, 0UL );
         checkNonZeros( diag, 2UL, 1UL );
         checkNonZeros( diag, 3UL, 0UL );

         if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) != 0 || diag(0,3) != 0 ||
             diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 || diag(1,3) != 0 ||
             diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 || diag(2,3) != 0 ||
             diag(3,0) != 0 || diag(3,1) != 0 || diag(3,2) != 0 || diag(3,3) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 0 3 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos != diag.end( 0UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to erase a zero diagonal element
      {
         Iterator pos = diag.erase( 1UL, diag.find( 1UL, 1UL ) );

         checkRows    ( diag, 4UL );
         checkColumns ( diag, 4UL );
         checkCapacity( diag, 4UL );
         checkNonZeros( diag, 1UL );
         checkNonZeros( diag, 0UL, 0UL );
         checkNonZeros( diag, 1UL, 0UL );
         checkNonZeros( diag, 2UL, 1UL );
         checkNonZeros( diag, 3UL, 0UL );

         if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) != 0 || diag(0,3) != 0 ||
             diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 || diag(1,3) != 0 ||
             diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 || diag(2,3) != 0 ||
             diag(3,0) != 0 || diag(3,1) != 0 || diag(3,2) != 0 || diag(3,3) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 0 3 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos != diag.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to erase a zero lower element
      {
         Iterator pos = diag.erase( 2UL, diag.find( 2UL, 1UL ) );

         checkRows    ( diag, 4UL );
         checkColumns ( diag, 4UL );
         checkCapacity( diag, 4UL );
         checkNonZeros( diag, 1UL );
         checkNonZeros( diag, 0UL, 0UL );
         checkNonZeros( diag, 1UL, 0UL );
         checkNonZeros( diag, 2UL, 1UL );
         checkNonZeros( diag, 3UL, 0UL );

         if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) != 0 || diag(0,3) != 0 ||
             diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 || diag(1,3) != 0 ||
             diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 || diag(2,3) != 0 ||
             diag(3,0) != 0 || diag(3,1) != 0 || diag(3,2) != 0 || diag(3,3) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 0 3 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos != diag.end( 2UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to erase a zero upper element
      {
         Iterator pos = diag.erase( 1UL, diag.find( 1UL, 2UL ) );

         checkRows    ( diag, 4UL );
         checkColumns ( diag, 4UL );
         checkCapacity( diag, 4UL );
         checkNonZeros( diag, 1UL );
         checkNonZeros( diag, 0UL, 0UL );
         checkNonZeros( diag, 1UL, 0UL );
         checkNonZeros( diag, 2UL, 1UL );
         checkNonZeros( diag, 3UL, 0UL );

         if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) != 0 || diag(0,3) != 0 ||
             diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 || diag(1,3) != 0 ||
             diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 || diag(2,3) != 0 ||
             diag(3,0) != 0 || diag(3,1) != 0 || diag(3,2) != 0 || diag(3,3) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 0 3 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos != diag.end( 1UL ) ) {
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
      test_ = "Row-major DiagonalMatrix::erase( size_t, Iterator, Iterator )";

      using Iterator = DT::Iterator;

      // Initialization check
      DT diag( 4UL, 4UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;
      diag(3,3) = 4;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkCapacity( diag, 4UL );
      checkNonZeros( diag, 4UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );
      checkNonZeros( diag, 3UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) != 0 || diag(3,2) != 0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 0 2 0 0 )\n( 0 0 3 0 )\n( 0 0 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the elements from the beginning of row 1 to the row end
      {
         Iterator pos = diag.erase( 1UL, diag.begin( 1UL ), diag.end( 1UL ) );

         checkRows    ( diag, 4UL );
         checkColumns ( diag, 4UL );
         checkCapacity( diag, 4UL );
         checkNonZeros( diag, 3UL );
         checkNonZeros( diag, 0UL, 1UL );
         checkNonZeros( diag, 1UL, 0UL );
         checkNonZeros( diag, 2UL, 1UL );
         checkNonZeros( diag, 3UL, 1UL );

         if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 || diag(0,3) != 0 ||
             diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 || diag(1,3) != 0 ||
             diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 || diag(2,3) != 0 ||
             diag(3,0) != 0 || diag(3,1) != 0 || diag(3,2) != 0 || diag(3,3) != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a single-element range failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 0 0 0 0 )\n( 0 0 3 0 )\n( 0 0 0 4 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos != diag.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the elements from (2,2) to the row end
      {
         Iterator pos = diag.erase( 2UL, diag.find( 2UL, 2UL ), diag.end( 2UL ) );

         checkRows    ( diag, 4UL );
         checkColumns ( diag, 4UL );
         checkCapacity( diag, 4UL );
         checkNonZeros( diag, 2UL );
         checkNonZeros( diag, 0UL, 1UL );
         checkNonZeros( diag, 1UL, 0UL );
         checkNonZeros( diag, 2UL, 0UL );
         checkNonZeros( diag, 3UL, 1UL );

         if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 || diag(0,3) != 0 ||
             diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 || diag(1,3) != 0 ||
             diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 0 || diag(2,3) != 0 ||
             diag(3,0) != 0 || diag(3,1) != 0 || diag(3,2) != 0 || diag(3,3) != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a single-element range failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 0 0 4 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos != diag.end( 2UL ) ) {
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
         Iterator pos = diag.erase( 3UL, diag.find( 3UL, 3UL ), diag.find( 3UL, 3UL ) );

         checkRows    ( diag, 4UL );
         checkColumns ( diag, 4UL );
         checkCapacity( diag, 4UL );
         checkNonZeros( diag, 2UL );
         checkNonZeros( diag, 0UL, 1UL );
         checkNonZeros( diag, 1UL, 0UL );
         checkNonZeros( diag, 2UL, 0UL );
         checkNonZeros( diag, 3UL, 1UL );

         if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 || diag(0,3) != 0 ||
             diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 || diag(1,3) != 0 ||
             diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 0 || diag(2,3) != 0 ||
             diag(3,0) != 0 || diag(3,1) != 0 || diag(3,2) != 0 || diag(3,3) != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing an empty range failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 0 0 4 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos->value() != 4 || pos->index() != 3 ) {
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
      }
   }


   //=====================================================================================
   // Row-major erase function with predicate
   //=====================================================================================

   {
      test_ = "Row-major DiagonalMatrix::erase( Predicate )";

      // Initialization check
      DT diag( 4UL, 4UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;
      diag(3,3) = 4;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkCapacity( diag, 4UL );
      checkNonZeros( diag, 4UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );
      checkNonZeros( diag, 3UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) != 0 || diag(3,2) != 0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 0 2 0 0 )\n( 0 0 3 0 )\n( 0 0 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing a selection of elements
      diag.erase( []( int value ){ return value == 1 || value == 3; } );

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkCapacity( diag, 4UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 0UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 0UL );
      checkNonZeros( diag, 3UL, 1UL );

      if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) != 0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 0 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) != 0 || diag(3,2) != 0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a selection of elements failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 0 0 0 0 )\n( 0 2 0 0 )\n( 0 0 0 0 )\n( 0 0 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to erase all elements with value 1
      diag.erase( []( int value ){ return value == 1; } );

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkCapacity( diag, 4UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 0UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 0UL );
      checkNonZeros( diag, 3UL, 1UL );

      if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) != 0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 0 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) != 0 || diag(3,2) != 0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing all elements with value 1 failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 0 0 0 0 )\n( 0 2 0 0 )\n( 0 0 0 0 )\n( 0 0 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major iterator-range-based erase function with predicate
   //=====================================================================================

   {
      test_ = "Row-major DiagonalMatrix::erase( size_t, Iterator, Iterator, Predicate )";

      // Initialization check
      DT diag( 4UL, 4UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;
      diag(3,3) = 4;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkCapacity( diag, 4UL );
      checkNonZeros( diag, 4UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );
      checkNonZeros( diag, 3UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) != 0 || diag(3,2) != 0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 0 2 0 0 )\n( 0 0 3 0 )\n( 0 0 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing a selection of elements
      diag.erase( 1UL, diag.begin( 1UL ), diag.end( 1UL ), []( int value ){ return value == 2; } );

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkCapacity( diag, 4UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );
      checkNonZeros( diag, 3UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) != 0 || diag(3,2) != 0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a selection of elements failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 0 0 0 0 )\n( 0 0 3 0 )\n( 0 0 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to erase from an empty range
      diag.erase( 2UL, diag.begin( 2UL ), diag.begin( 2UL ), []( int ){ return true; } );

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkCapacity( diag, 4UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );
      checkNonZeros( diag, 3UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) != 0 || diag(3,2) != 0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing from an empty range failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 0 0 0 0 )\n( 0 0 3 0 )\n( 0 0 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major index-based erase function
   //=====================================================================================

   {
      test_ = "Column-major DiagonalMatrix::erase( size_t, size_t )";

      // Initialization check
      ODT diag( 4UL, 4UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;
      diag(3,3) = 4;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkCapacity( diag, 4UL );
      checkNonZeros( diag, 4UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );
      checkNonZeros( diag, 3UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) != 0 || diag(3,2) != 0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 0 2 0 0 )\n( 0 0 3 0 )\n( 0 0 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element at (1,1)
      diag.erase( 1UL, 1UL );

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkCapacity( diag, 4UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );
      checkNonZeros( diag, 3UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) != 0 || diag(3,2) != 0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 0 0 0 0 )\n( 0 0 3 0 )\n( 0 0 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element at (3,3)
      diag.erase( 3UL, 3UL );

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkCapacity( diag, 4UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );
      checkNonZeros( diag, 3UL, 0UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) != 0 || diag(3,2) != 0 || diag(3,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 0 0 0 0 )\n( 0 0 3 0 )\n( 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element at (0,0)
      diag.erase( 0UL, size_t(0) );

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkCapacity( diag, 4UL );
      checkNonZeros( diag, 1UL );
      checkNonZeros( diag, 0UL, 0UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );
      checkNonZeros( diag, 3UL, 0UL );

      if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) != 0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) != 0 || diag(3,2) != 0 || diag(3,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 0 3 0 )\n( 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to erase a zero diagonal element
      diag.erase( 1UL, 1UL );

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkCapacity( diag, 4UL );
      checkNonZeros( diag, 1UL );
      checkNonZeros( diag, 0UL, 0UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );
      checkNonZeros( diag, 3UL, 0UL );

      if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) != 0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) != 0 || diag(3,2) != 0 || diag(3,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 0 3 0 )\n( 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to erase a zero lower element
      diag.erase( 2UL, 1UL );

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkCapacity( diag, 4UL );
      checkNonZeros( diag, 1UL );
      checkNonZeros( diag, 0UL, 0UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );
      checkNonZeros( diag, 3UL, 0UL );

      if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) != 0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) != 0 || diag(3,2) != 0 || diag(3,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 0 3 0 )\n( 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to erase a zero upper element
      diag.erase( 1UL, 2UL );

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkCapacity( diag, 4UL );
      checkNonZeros( diag, 1UL );
      checkNonZeros( diag, 0UL, 0UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );
      checkNonZeros( diag, 3UL, 0UL );

      if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) != 0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) != 0 || diag(3,2) != 0 || diag(3,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 0 3 0 )\n( 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major iterator-based erase function
   //=====================================================================================

   {
      test_ = "Column-major DiagonalMatrix::erase( size_t, Iterator )";

      using Iterator = ODT::Iterator;

      // Initialization check
      ODT diag( 4UL, 4UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;
      diag(3,3) = 4;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkCapacity( diag, 4UL );
      checkNonZeros( diag, 4UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );
      checkNonZeros( diag, 3UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) != 0 || diag(3,2) != 0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 0 2 0 0 )\n( 0 0 3 0 )\n( 0 0 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element at (1,1)
      {
         Iterator pos = diag.erase( 1UL, diag.find( 1UL, 1UL ) );

         checkRows    ( diag, 4UL );
         checkColumns ( diag, 4UL );
         checkCapacity( diag, 4UL );
         checkNonZeros( diag, 3UL );
         checkNonZeros( diag, 0UL, 1UL );
         checkNonZeros( diag, 1UL, 0UL );
         checkNonZeros( diag, 2UL, 1UL );
         checkNonZeros( diag, 3UL, 1UL );

         if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 || diag(0,3) != 0 ||
             diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 || diag(1,3) != 0 ||
             diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 || diag(2,3) != 0 ||
             diag(3,0) != 0 || diag(3,1) != 0 || diag(3,2) != 0 || diag(3,3) != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 0 0 0 0 )\n( 0 0 3 0 )\n( 0 0 0 4 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos != diag.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the element at (3,3)
      {
         Iterator pos = diag.erase( 3UL, diag.find( 3UL, 3UL ) );

         checkRows    ( diag, 4UL );
         checkColumns ( diag, 4UL );
         checkCapacity( diag, 4UL );
         checkNonZeros( diag, 2UL );
         checkNonZeros( diag, 0UL, 1UL );
         checkNonZeros( diag, 1UL, 0UL );
         checkNonZeros( diag, 2UL, 1UL );
         checkNonZeros( diag, 3UL, 0UL );

         if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 || diag(0,3) != 0 ||
             diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 || diag(1,3) != 0 ||
             diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 || diag(2,3) != 0 ||
             diag(3,0) != 0 || diag(3,1) != 0 || diag(3,2) != 0 || diag(3,3) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 0 0 0 0 )\n( 0 0 3 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos != diag.end( 3UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the element at (0,0)
      {
         Iterator pos = diag.erase( 0UL, diag.find( 0UL, 0UL ) );

         checkRows    ( diag, 4UL );
         checkColumns ( diag, 4UL );
         checkCapacity( diag, 4UL );
         checkNonZeros( diag, 1UL );
         checkNonZeros( diag, 0UL, 0UL );
         checkNonZeros( diag, 1UL, 0UL );
         checkNonZeros( diag, 2UL, 1UL );
         checkNonZeros( diag, 3UL, 0UL );

         if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) != 0 || diag(0,3) != 0 ||
             diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 || diag(1,3) != 0 ||
             diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 || diag(2,3) != 0 ||
             diag(3,0) != 0 || diag(3,1) != 0 || diag(3,2) != 0 || diag(3,3) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 0 3 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos != diag.end( 0UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to erase a zero diagonal element
      {
         Iterator pos = diag.erase( 1UL, diag.find( 1UL, 1UL ) );

         checkRows    ( diag, 4UL );
         checkColumns ( diag, 4UL );
         checkCapacity( diag, 4UL );
         checkNonZeros( diag, 1UL );
         checkNonZeros( diag, 0UL, 0UL );
         checkNonZeros( diag, 1UL, 0UL );
         checkNonZeros( diag, 2UL, 1UL );
         checkNonZeros( diag, 3UL, 0UL );

         if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) != 0 || diag(0,3) != 0 ||
             diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 || diag(1,3) != 0 ||
             diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 || diag(2,3) != 0 ||
             diag(3,0) != 0 || diag(3,1) != 0 || diag(3,2) != 0 || diag(3,3) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 0 3 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos != diag.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to erase a zero lower element
      {
         Iterator pos = diag.erase( 1UL, diag.find( 2UL, 1UL ) );

         checkRows    ( diag, 4UL );
         checkColumns ( diag, 4UL );
         checkCapacity( diag, 4UL );
         checkNonZeros( diag, 1UL );
         checkNonZeros( diag, 0UL, 0UL );
         checkNonZeros( diag, 1UL, 0UL );
         checkNonZeros( diag, 2UL, 1UL );
         checkNonZeros( diag, 3UL, 0UL );

         if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) != 0 || diag(0,3) != 0 ||
             diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 || diag(1,3) != 0 ||
             diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 || diag(2,3) != 0 ||
             diag(3,0) != 0 || diag(3,1) != 0 || diag(3,2) != 0 || diag(3,3) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 0 3 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos != diag.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to erase a zero upper element
      {
         Iterator pos = diag.erase( 2UL, diag.find( 1UL, 2UL ) );

         checkRows    ( diag, 4UL );
         checkColumns ( diag, 4UL );
         checkCapacity( diag, 4UL );
         checkNonZeros( diag, 1UL );
         checkNonZeros( diag, 0UL, 0UL );
         checkNonZeros( diag, 1UL, 0UL );
         checkNonZeros( diag, 2UL, 1UL );
         checkNonZeros( diag, 3UL, 0UL );

         if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) != 0 || diag(0,3) != 0 ||
             diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 || diag(1,3) != 0 ||
             diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 || diag(2,3) != 0 ||
             diag(3,0) != 0 || diag(3,1) != 0 || diag(3,2) != 0 || diag(3,3) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 0 3 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos != diag.end( 2UL ) ) {
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
      test_ = "Column-major DiagonalMatrix::erase( size_t, Iterator, Iterator )";

      using Iterator = ODT::Iterator;

      // Initialization check
      ODT diag( 4UL, 4UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;
      diag(3,3) = 4;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkCapacity( diag, 4UL );
      checkNonZeros( diag, 4UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );
      checkNonZeros( diag, 3UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) != 0 || diag(3,2) != 0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 0 2 0 0 )\n( 0 0 3 0 )\n( 0 0 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the elements from the beginning of row 1 to the row end
      {
         Iterator pos = diag.erase( 1UL, diag.begin( 1UL ), diag.end( 1UL ) );

         checkRows    ( diag, 4UL );
         checkColumns ( diag, 4UL );
         checkCapacity( diag, 4UL );
         checkNonZeros( diag, 3UL );
         checkNonZeros( diag, 0UL, 1UL );
         checkNonZeros( diag, 1UL, 0UL );
         checkNonZeros( diag, 2UL, 1UL );
         checkNonZeros( diag, 3UL, 1UL );

         if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 || diag(0,3) != 0 ||
             diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 || diag(1,3) != 0 ||
             diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 || diag(2,3) != 0 ||
             diag(3,0) != 0 || diag(3,1) != 0 || diag(3,2) != 0 || diag(3,3) != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a single-element range failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 0 0 0 0 )\n( 0 0 3 0 )\n( 0 0 0 4 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos != diag.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the elements from (2,2) to the row end
      {
         Iterator pos = diag.erase( 2UL, diag.find( 2UL, 2UL ), diag.end( 2UL ) );

         checkRows    ( diag, 4UL );
         checkColumns ( diag, 4UL );
         checkCapacity( diag, 4UL );
         checkNonZeros( diag, 2UL );
         checkNonZeros( diag, 0UL, 1UL );
         checkNonZeros( diag, 1UL, 0UL );
         checkNonZeros( diag, 2UL, 0UL );
         checkNonZeros( diag, 3UL, 1UL );

         if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 || diag(0,3) != 0 ||
             diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 || diag(1,3) != 0 ||
             diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 0 || diag(2,3) != 0 ||
             diag(3,0) != 0 || diag(3,1) != 0 || diag(3,2) != 0 || diag(3,3) != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a single-element range failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 0 0 4 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos != diag.end( 2UL ) ) {
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
         Iterator pos = diag.erase( 3UL, diag.find( 3UL, 3UL ), diag.find( 3UL, 3UL ) );

         checkRows    ( diag, 4UL );
         checkColumns ( diag, 4UL );
         checkCapacity( diag, 4UL );
         checkNonZeros( diag, 2UL );
         checkNonZeros( diag, 0UL, 1UL );
         checkNonZeros( diag, 1UL, 0UL );
         checkNonZeros( diag, 2UL, 0UL );
         checkNonZeros( diag, 3UL, 1UL );

         if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 || diag(0,3) != 0 ||
             diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 || diag(1,3) != 0 ||
             diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 0 || diag(2,3) != 0 ||
             diag(3,0) != 0 || diag(3,1) != 0 || diag(3,2) != 0 || diag(3,3) != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing an empty range failed\n"
                << " Details:\n"
                << "   Result:\n" << diag << "\n"
                << "   Expected result:\n( 1 0 0 0 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n( 0 0 0 4 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos->value() != 4 || pos->index() != 3 ) {
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
      }
   }


   //=====================================================================================
   // Column-major erase function with predicate
   //=====================================================================================

   {
      test_ = "Column-major DiagonalMatrix::erase( Predicate )";

      // Initialization check
      ODT diag( 4UL, 4UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;
      diag(3,3) = 4;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkCapacity( diag, 4UL );
      checkNonZeros( diag, 4UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );
      checkNonZeros( diag, 3UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) != 0 || diag(3,2) != 0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 0 2 0 0 )\n( 0 0 3 0 )\n( 0 0 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing a selection of elements
      diag.erase( []( int value ){ return value == 1 || value == 3; } );

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkCapacity( diag, 4UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 0UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 0UL );
      checkNonZeros( diag, 3UL, 1UL );

      if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) != 0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 0 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) != 0 || diag(3,2) != 0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a selection of elements failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 0 0 0 0 )\n( 0 2 0 0 )\n( 0 0 0 0 )\n( 0 0 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to erase all elements with value 1
      diag.erase( []( int value ){ return value == 1; } );

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkCapacity( diag, 4UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 0UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 0UL );
      checkNonZeros( diag, 3UL, 1UL );

      if( diag(0,0) != 0 || diag(0,1) != 0 || diag(0,2) != 0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 0 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) != 0 || diag(3,2) != 0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing all elements with value 1 failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 0 0 0 0 )\n( 0 2 0 0 )\n( 0 0 0 0 )\n( 0 0 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major iterator-range-based erase function with predicate
   //=====================================================================================

   {
      test_ = "Column-major DiagonalMatrix::erase( size_t, Iterator, Iterator, Predicate )";

      // Initialization check
      ODT diag( 4UL, 4UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;
      diag(3,3) = 4;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkCapacity( diag, 4UL );
      checkNonZeros( diag, 4UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );
      checkNonZeros( diag, 3UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) != 0 || diag(3,2) != 0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 0 2 0 0 )\n( 0 0 3 0 )\n( 0 0 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing a selection of elements
      diag.erase( 1UL, diag.begin( 1UL ), diag.end( 1UL ), []( int value ){ return value == 2; } );

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkCapacity( diag, 4UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );
      checkNonZeros( diag, 3UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) != 0 || diag(3,2) != 0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a selection of elements failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 0 0 0 0 )\n( 0 0 3 0 )\n( 0 0 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to erase from an empty range
      diag.erase( 2UL, diag.begin( 2UL ), diag.begin( 2UL ), []( int ){ return true; } );

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkCapacity( diag, 4UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );
      checkNonZeros( diag, 3UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) != 0 || diag(3,2) != 0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing from an empty range failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 0 0 0 0 )\n( 0 0 3 0 )\n( 0 0 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c find() member function of the DiagonalMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c find() member function of the DiagonalMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testFind()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major DiagonalMatrix::find()";

      using ConstIterator = DT::ConstIterator;

      // Initialization check
      DT diag( 8UL, 3UL );
      diag(2,2) = 1;
      diag(3,3) = 2;
      diag(6,6) = 3;

      checkRows    ( diag, 8UL );
      checkColumns ( diag, 8UL );
      checkCapacity( diag, 3UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 0UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );
      checkNonZeros( diag, 3UL, 1UL );
      checkNonZeros( diag, 4UL, 0UL );
      checkNonZeros( diag, 5UL, 0UL );
      checkNonZeros( diag, 6UL, 1UL );
      checkNonZeros( diag, 7UL, 0UL );

      // Searching for the first element
      {
         ConstIterator pos( diag.find( 2UL, 2UL ) );

         if( pos == diag.end( 2UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required position = (2,2)\n"
                << "   Current matrix:\n" << diag << "\n";
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
                << "   Current matrix:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for the second element
      {
         ConstIterator pos( diag.find( 3UL, 3UL ) );

         if( pos == diag.end( 3UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required position = (3,3)\n"
                << "   Current matrix:\n" << diag << "\n";
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
                << "   Current matrix:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for the third element
      {
         ConstIterator pos( diag.find( 6UL, 6UL ) );

         if( pos == diag.end( 6UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required position = (6,6)\n"
                << "   Current matrix:\n" << diag << "\n";
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
                << "   Current matrix:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for a non-existing non-zero element
      {
         ConstIterator pos( diag.find( 4UL, 0UL ) );

         if( pos != diag.end( 4UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Non-existing element could be found\n"
                << " Details:\n"
                << "   Required index = 0\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 0\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major DiagonalMatrix::find()";

      using ConstIterator = ODT::ConstIterator;

      // Initialization check
      ODT diag( 8UL, 3UL );
      diag(2,2) = 1;
      diag(3,3) = 2;
      diag(6,6) = 3;

      checkRows    ( diag, 8UL );
      checkColumns ( diag, 8UL );
      checkCapacity( diag, 3UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 0UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );
      checkNonZeros( diag, 3UL, 1UL );
      checkNonZeros( diag, 4UL, 0UL );
      checkNonZeros( diag, 5UL, 0UL );
      checkNonZeros( diag, 6UL, 1UL );
      checkNonZeros( diag, 7UL, 0UL );

      // Searching for the first element
      {
         ConstIterator pos( diag.find( 2UL, 2UL ) );

         if( pos == diag.end( 2UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required position = (2,2)\n"
                << "   Current matrix:\n" << diag << "\n";
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
                << "   Current matrix:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for the second element
      {
         ConstIterator pos( diag.find( 3UL, 3UL ) );

         if( pos == diag.end( 3UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required position = (3,3)\n"
                << "   Current matrix:\n" << diag << "\n";
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
                << "   Current matrix:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for the third element
      {
         ConstIterator pos( diag.find( 6UL, 6UL ) );

         if( pos == diag.end( 6UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required position = (6,6)\n"
                << "   Current matrix:\n" << diag << "\n";
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
                << "   Current matrix:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for a non-existing non-zero element
      {
         ConstIterator pos( diag.find( 4UL, 0UL ) );

         if( pos != diag.end( 0UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Non-existing element could be found\n"
                << " Details:\n"
                << "   Required index = 4\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 0\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c lowerBound() member function of the DiagonalMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c lowerBound() member function of the DiagonalMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testLowerBound()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major DiagonalMatrix::lowerBound()";

      using ConstIterator = DT::ConstIterator;

      // Initialization check
      DT diag( 3UL, 1UL );
      diag(1,1) = 1;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 1UL );
      checkNonZeros( diag, 1UL );
      checkNonZeros( diag, 0UL, 0UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 0UL );

      // Determining the lower bound for position (1,0)
      {
         ConstIterator pos( diag.lowerBound( 1UL, 0UL ) );

         if( pos == diag.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,0)\n"
                << "   Current matrix:\n" << diag << "\n";
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
                << "   Current matrix:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for position (1,1)
      {
         ConstIterator pos( diag.lowerBound( 1UL, 1UL ) );

         if( pos == diag.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,1)\n"
                << "   Current matrix:\n" << diag << "\n";
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
                << "   Current matrix:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for position (1,2)
      {
         ConstIterator pos( diag.lowerBound( 1UL, 2UL ) );

         if( pos != diag.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,2)\n"
                << "   Current matrix:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major DiagonalMatrix::lowerBound()";

      using ConstIterator = ODT::ConstIterator;

      // Initialization check
      ODT diag( 3UL, 1UL );
      diag(1,1) = 1;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 1UL );
      checkNonZeros( diag, 1UL );
      checkNonZeros( diag, 0UL, 0UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 0UL );

      // Determining the lower bound for position (0,1)
      {
         ConstIterator pos( diag.lowerBound( 0UL, 1UL ) );

         if( pos == diag.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (0,1)\n"
                << "   Current matrix:\n" << diag << "\n";
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
                << "   Current matrix:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for position (1,1)
      {
         ConstIterator pos( diag.lowerBound( 1UL, 1UL ) );

         if( pos == diag.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,1)\n"
                << "   Current matrix:\n" << diag << "\n";
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
                << "   Current matrix:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for position (2,1)
      {
         ConstIterator pos( diag.lowerBound( 2UL, 1UL ) );

         if( pos != diag.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (2,1)\n"
                << "   Current matrix:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c upperBound() member function of the DiagonalMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c upperBound() member function of the DiagonalMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testUpperBound()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major DiagonalMatrix::upperBound()";

      using ConstIterator = DT::ConstIterator;

      // Initialization check
      DT diag( 3UL, 1UL );
      diag(1,1) = 1;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 1UL );
      checkNonZeros( diag, 1UL );
      checkNonZeros( diag, 0UL, 0UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 0UL );

      // Determining the upper bound for position (1,0)
      {
         ConstIterator pos( diag.upperBound( 1UL, 0UL ) );

         if( pos == diag.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,0)\n"
                << "   Current matrix:\n" << diag << "\n";
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
                << "   Current matrix:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for position (1,1)
      {
         ConstIterator pos( diag.upperBound( 1UL, 1UL ) );

         if( pos != diag.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,1)\n"
                << "   Current matrix:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for position (1,2)
      {
         ConstIterator pos( diag.upperBound( 1UL, 2UL ) );

         if( pos != diag.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,2)\n"
                << "   Current matrix:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major DiagonalMatrix::lowerBound()";

      using ConstIterator = ODT::ConstIterator;

      // Initialization check
      ODT diag( 3UL, 1UL );
      diag(1,1) = 1;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 1UL );
      checkNonZeros( diag, 1UL );
      checkNonZeros( diag, 0UL, 0UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 0UL );

      // Determining the upper bound for position (0,1)
      {
         ConstIterator pos( diag.upperBound( 0UL, 1UL ) );

         if( pos == diag.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (0,1)\n"
                << "   Current matrix:\n" << diag << "\n";
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
                << "   Current matrix:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for position (1,1)
      {
         ConstIterator pos( diag.upperBound( 1UL, 1UL ) );

         if( pos != diag.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,1)\n"
                << "   Current matrix:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for position (2,1)
      {
         ConstIterator pos( diag.upperBound( 2UL, 1UL ) );

         if( pos != diag.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (2,1)\n"
                << "   Current matrix:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c isDefault() function with the DiagonalMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isDefault() function with the DiagonalMatrix
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
         DT diag;

         if( isDefault( diag ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with default matrix
      {
         DT diag( 3UL );

         if( isDefault( diag(1,1) ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix element:\n" << diag(1,1) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( diag ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with non-default matrix
      {
         DT diag( 3UL );
         diag(1,1) = 1;

         if( isDefault( diag(1,1) ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix element:\n" << diag(1,1) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( diag ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << diag << "\n";
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
         ODT diag;

         if( isDefault( diag ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with default matrix
      {
         ODT diag( 3UL );

         if( isDefault( diag(1,1) ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix element:\n" << diag(1,1) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( diag ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with non-default matrix
      {
         ODT diag( 3UL );
         diag(1,1) = 1;

         if( isDefault( diag(1,1) ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix element:\n" << diag(1,1) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( diag ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << diag << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c submatrix() function with the DiagonalMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c submatrix() function with the DiagonalMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testSubmatrix()
{
   //=====================================================================================
   // Row-major general tests
   //=====================================================================================

   {
      test_ = "Row-major submatrix() function";

      using SMT = blaze::Submatrix<DT>;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      SMT sm = submatrix( diag, 1UL, 1UL, 2UL, 2UL );

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

      sm(0,0) = -5;

      if( sm(0,0) != -5 || sm(0,1) != 0 ||
          sm(1,0) !=  0 || sm(1,1) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix access failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( -5  0 )\n(  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != -5 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix access failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0 )\n( 0 -5  0 )\n( 0  0  3 )\n";
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

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix reset failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major general tests
   //=====================================================================================

   {
      test_ = "Column-major submatrix() function";

      using SMT = blaze::Submatrix<ODT>;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      SMT sm = submatrix( diag, 1UL, 1UL, 2UL, 2UL );

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

      sm(0,0) = -5;

      if( sm(0,0) != -5 || sm(0,1) != 0 ||
          sm(1,0) !=  0 || sm(1,1) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix access failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( -5  0 )\n(  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != -5 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix access failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0 )\n( 0 -5  0 )\n( 0  0  3 )\n";
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

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix reset failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c row() function with the DiagonalMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c row() function with the DiagonalMatrix specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testRow()
{
   //=====================================================================================
   // Row-major general tests
   //=====================================================================================

   {
      test_ = "Row-major row() function";

      using RT = blaze::Row<DT>;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      RT row1 = row( diag, 1UL );

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

      if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != -5 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row access failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n(  1  0  0 )\n( -4 -5  0 )\n(  7  0  3 )\n";
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

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row reset failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major general tests
   //=====================================================================================

   {
      test_ = "Column-major row() function";

      using RT = blaze::Row<ODT>;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      RT row1 = row( diag, 1UL );

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
             << "   Expected result:\n( -4 -5  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != -5 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row access failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0 )\n( 0 -5  0 )\n( 0  0  3 )\n";
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

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row reset failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c column() function with the DiagonalMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c column() function with the DiagonalMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testColumn()
{
   //=====================================================================================
   // Row-major general tests
   //=====================================================================================

   {
      test_ = "Row-major column() function";

      using CT = blaze::Column<DT>;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      CT col1 = column( diag, 1UL );

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

      if( it == col1.end() || it->value() != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << it->value()<< "\n"
             << "   Expected result: 2\n";
         throw std::runtime_error( oss.str() );
      }

      col1[1] = -5;

      if( col1[0] != 0 || col1[1] != -5 || col1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column access failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( 0 -5  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != -5 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column access failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0 )\n( 0 -5  0 )\n( 0  0  3 )\n";
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

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column reset failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major general tests
   //=====================================================================================

   {
      test_ = "Column-major column() function";

      using CT = blaze::Column<ODT>;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      CT col1 = column( diag, 1UL );

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

      if( it == col1.end() || it->value()!= 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << it->value()<< "\n"
             << "   Expected result: 2\n";
         throw std::runtime_error( oss.str() );
      }

      col1[1] = -5;

      if( col1[0] != 0 || col1[1] != -5 || col1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column access failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( 0 -5  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != -5 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column access failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0 )\n( 0 -5  0 )\n( 0  0  3 )\n";
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

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column reset failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************

} // namespace diagonalmatrix

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
   std::cout << "   Running DiagonalMatrix sparse test (part 2)..." << std::endl;

   try
   {
      RUN_DIAGONALMATRIX_SPARSE_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during DiagonalMatrix sparse test (part 2):\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
