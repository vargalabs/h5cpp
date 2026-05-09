//=================================================================================================
/*!
//  \file src/mathtest/adaptors/symmetricmatrix/SparseScalarTest2.cpp
//  \brief Source file for the SymmetricMatrix sparse scalar test (part 2)
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
#include <blazetest/mathtest/adaptors/symmetricmatrix/SparseScalarTest.h>

#ifdef BLAZE_USE_HPX_THREADS
#  include <hpx/hpx_main.hpp>
#endif


namespace blazetest {

namespace mathtest {

namespace adaptors {

namespace symmetricmatrix {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the SymmetricMatrix sparse scalar test.
//
// \exception std::runtime_error Operation error detected.
*/
SparseScalarTest::SparseScalarTest()
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
   testTranspose();
   testCTranspose();
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
/*!\brief Test of all SymmetricMatrix (self-)scaling operations.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all available ways to scale an instance of the SymmetricMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseScalarTest::testScaling()
{
   //=====================================================================================
   // Row-major self-scaling (M*=s)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (M*=s)";

      ST sym( 3UL );
      sym(1,2) =  1;
      sym(2,0) = -2;
      sym(2,2) =  3;

      sym *= 2;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 5UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 1UL );
      checkNonZeros( sym, 1UL, 1UL );
      checkNonZeros( sym, 2UL, 3UL );

      if( sym(0,0) !=  0 || sym(0,1) != 0 || sym(0,2) != -4 ||
          sym(1,0) !=  0 || sym(1,1) != 0 || sym(1,2) !=  2 ||
          sym(2,0) != -4 || sym(2,1) != 2 || sym(2,2) !=  6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  0 0 -4 )\n(  0 0  2 )\n( -4 2  6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major self-scaling (M=M*s)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (M=M*s)";

      ST sym( 3UL );
      sym(1,2) =  1;
      sym(2,0) = -2;
      sym(2,2) =  3;

      sym = sym * 2;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 5UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 1UL );
      checkNonZeros( sym, 1UL, 1UL );
      checkNonZeros( sym, 2UL, 3UL );

      if( sym(0,0) !=  0 || sym(0,1) != 0 || sym(0,2) != -4 ||
          sym(1,0) !=  0 || sym(1,1) != 0 || sym(1,2) !=  2 ||
          sym(2,0) != -4 || sym(2,1) != 2 || sym(2,2) !=  6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  0 0 -4 )\n(  0 0  2 )\n( -4 2  6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major self-scaling (M=s*M)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (M=s*M)";

      ST sym( 3UL );
      sym(1,2) =  1;
      sym(2,0) = -2;
      sym(2,2) =  3;

      sym = 2 * sym;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 5UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 1UL );
      checkNonZeros( sym, 1UL, 1UL );
      checkNonZeros( sym, 2UL, 3UL );

      if( sym(0,0) !=  0 || sym(0,1) != 0 || sym(0,2) != -4 ||
          sym(1,0) !=  0 || sym(1,1) != 0 || sym(1,2) !=  2 ||
          sym(2,0) != -4 || sym(2,1) != 2 || sym(2,2) !=  6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  0 0 -4 )\n(  0 0  2 )\n( -4 2  6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major self-scaling (M/=s)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (M/=s)";

      ST sym( 3UL );
      sym(1,2) =  2;
      sym(2,0) = -4;
      sym(2,2) =  6;

      sym /= 2;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 5UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 1UL );
      checkNonZeros( sym, 1UL, 1UL );
      checkNonZeros( sym, 2UL, 3UL );

      if( sym(0,0) !=  0 || sym(0,1) != 0 || sym(0,2) != -2 ||
          sym(1,0) !=  0 || sym(1,1) != 0 || sym(1,2) !=  1 ||
          sym(2,0) != -2 || sym(2,1) != 1 || sym(2,2) !=  3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  0 0 -2 )\n(  0 0  1 )\n( -2 1  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major self-scaling (M=M/s)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (M=M/s)";

      ST sym( 3UL );
      sym(1,2) =  2;
      sym(2,0) = -4;
      sym(2,2) =  6;

      sym = sym / 2;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 5UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 1UL );
      checkNonZeros( sym, 1UL, 1UL );
      checkNonZeros( sym, 2UL, 3UL );

      if( sym(0,0) !=  0 || sym(0,1) != 0 || sym(0,2) != -2 ||
          sym(1,0) !=  0 || sym(1,1) != 0 || sym(1,2) !=  1 ||
          sym(2,0) != -2 || sym(2,1) != 1 || sym(2,2) !=  3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  0 0 -2 )\n(  0 0  1 )\n( -2 1  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major SymmetricMatrix::scale()
   //=====================================================================================

   {
      test_ = "Row-major SymmetricMatrix::scale()";

      // Initialization check
      ST sym( 3UL );
      sym(1,2) =  1;
      sym(2,0) = -2;
      sym(2,2) =  3;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 5UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 1UL );
      checkNonZeros( sym, 1UL, 1UL );
      checkNonZeros( sym, 2UL, 3UL );

      if( sym(0,0) !=  0 || sym(0,1) != 0 || sym(0,2) != -2 ||
          sym(1,0) !=  0 || sym(1,1) != 0 || sym(1,2) !=  1 ||
          sym(2,0) != -2 || sym(2,1) != 1 || sym(2,2) !=  3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  0 0 -2 )\n(  0 0 1 )\n( -2 1 3 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Integral scaling of the matrix
      sym.scale( 2 );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 5UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 1UL );
      checkNonZeros( sym, 1UL, 1UL );
      checkNonZeros( sym, 2UL, 3UL );

      if( sym(0,0) !=  0 || sym(0,1) != 0 || sym(0,2) != -4 ||
          sym(1,0) !=  0 || sym(1,1) != 0 || sym(1,2) !=  2 ||
          sym(2,0) != -4 || sym(2,1) != 2 || sym(2,2) !=  6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  0 0 -4 )\n(  0 0 2 )\n( -4 2 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Floating point scaling of the matrix
      sym.scale( 0.5 );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 5UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 1UL );
      checkNonZeros( sym, 1UL, 1UL );
      checkNonZeros( sym, 2UL, 3UL );

      if( sym(0,0) !=  0 || sym(0,1) != 0 || sym(0,2) != -2 ||
          sym(1,0) !=  0 || sym(1,1) != 0 || sym(1,2) !=  1 ||
          sym(2,0) != -2 || sym(2,1) != 1 || sym(2,2) !=  3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  0 0 -2 )\n(  0 0 1 )\n( -2 1 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major SymmetricMatrix::scale() (complex)";

      using blaze::complex;

      blaze::SymmetricMatrix< blaze::CompressedMatrix<complex<float>,blaze::rowMajor> > sym( 2UL );
      sym(0,0) = complex<float>( 1.0F, 0.0F );
      sym(0,1) = complex<float>( 2.0F, 0.0F );
      sym(1,1) = complex<float>( 4.0F, 0.0F );

      sym.scale( complex<float>( 3.0F, 0.0F ) );

      checkRows    ( sym, 2UL );
      checkColumns ( sym, 2UL );
      checkCapacity( sym, 4UL );
      checkNonZeros( sym, 4UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );

      if( sym(0,0) != complex<float>( 3.0F, 0.0F ) || sym(0,1) != complex<float>(  6.0F, 0.0F ) ||
          sym(1,0) != complex<float>( 6.0F, 0.0F ) || sym(1,1) != complex<float>( 12.0F, 0.0F ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( ( 3,0) ( 6,0)\n( 6,0) (12,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (M*=s)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (M*=s)";

      OST sym( 3UL );
      sym(1,2) =  1;
      sym(2,0) = -2;
      sym(2,2) =  3;

      sym *= 2;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 5UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 1UL );
      checkNonZeros( sym, 1UL, 1UL );
      checkNonZeros( sym, 2UL, 3UL );

      if( sym(0,0) !=  0 || sym(0,1) != 0 || sym(0,2) != -4 ||
          sym(1,0) !=  0 || sym(1,1) != 0 || sym(1,2) !=  2 ||
          sym(2,0) != -4 || sym(2,1) != 2 || sym(2,2) !=  6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  0 0 -4 )\n(  0 0  2 )\n( -4 2  6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (M=M*s)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (M=M*s)";

      OST sym( 3UL );
      sym(1,2) =  1;
      sym(2,0) = -2;
      sym(2,2) =  3;

      sym = sym * 2;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 5UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 1UL );
      checkNonZeros( sym, 1UL, 1UL );
      checkNonZeros( sym, 2UL, 3UL );

      if( sym(0,0) !=  0 || sym(0,1) != 0 || sym(0,2) != -4 ||
          sym(1,0) !=  0 || sym(1,1) != 0 || sym(1,2) !=  2 ||
          sym(2,0) != -4 || sym(2,1) != 2 || sym(2,2) !=  6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  0 0 -4 )\n(  0 0  2 )\n( -4 2  6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (M=s*M)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (M=s*M)";

      OST sym( 3UL );
      sym(1,2) =  1;
      sym(2,0) = -2;
      sym(2,2) =  3;

      sym = 2 * sym;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 5UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 1UL );
      checkNonZeros( sym, 1UL, 1UL );
      checkNonZeros( sym, 2UL, 3UL );

      if( sym(0,0) !=  0 || sym(0,1) != 0 || sym(0,2) != -4 ||
          sym(1,0) !=  0 || sym(1,1) != 0 || sym(1,2) !=  2 ||
          sym(2,0) != -4 || sym(2,1) != 2 || sym(2,2) !=  6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  0 0 -4 )\n(  0 0  2 )\n( -4 2  6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (M/=s)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (M/=s)";

      OST sym( 3UL );
      sym(1,2) =  2;
      sym(2,0) = -4;
      sym(2,2) =  6;

      sym /= 2;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 5UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 1UL );
      checkNonZeros( sym, 1UL, 1UL );
      checkNonZeros( sym, 2UL, 3UL );

      if( sym(0,0) !=  0 || sym(0,1) != 0 || sym(0,2) != -2 ||
          sym(1,0) !=  0 || sym(1,1) != 0 || sym(1,2) !=  1 ||
          sym(2,0) != -2 || sym(2,1) != 1 || sym(2,2) !=  3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  0 0 -2 )\n(  0 0  1 )\n( -2 1  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (M=M/s)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (M=M/s)";

      OST sym( 3UL );
      sym(1,2) =  2;
      sym(2,0) = -4;
      sym(2,2) =  6;

      sym = sym / 2;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 5UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 1UL );
      checkNonZeros( sym, 1UL, 1UL );
      checkNonZeros( sym, 2UL, 3UL );

      if( sym(0,0) !=  0 || sym(0,1) != 0 || sym(0,2) != -2 ||
          sym(1,0) !=  0 || sym(1,1) != 0 || sym(1,2) !=  1 ||
          sym(2,0) != -2 || sym(2,1) != 1 || sym(2,2) !=  3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  0 0 -2 )\n(  0 0  1 )\n( -2 1  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major SymmetricMatrix::scale()
   //=====================================================================================

   {
      test_ = "Row-major SymmetricMatrix::scale()";

      // Initialization check
      OST sym( 3UL );
      sym(1,2) =  1;
      sym(2,0) = -2;
      sym(2,2) =  3;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 5UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 1UL );
      checkNonZeros( sym, 1UL, 1UL );
      checkNonZeros( sym, 2UL, 3UL );

      if( sym(0,0) !=  0 || sym(0,1) != 0 || sym(0,2) != -2 ||
          sym(1,0) !=  0 || sym(1,1) != 0 || sym(1,2) !=  1 ||
          sym(2,0) != -2 || sym(2,1) != 1 || sym(2,2) !=  3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  0 0 -2 )\n(  0 0 1 )\n( -2 1 3 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Integral scaling of the matrix
      sym.scale( 2 );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 5UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 1UL );
      checkNonZeros( sym, 1UL, 1UL );
      checkNonZeros( sym, 2UL, 3UL );

      if( sym(0,0) !=  0 || sym(0,1) != 0 || sym(0,2) != -4 ||
          sym(1,0) !=  0 || sym(1,1) != 0 || sym(1,2) !=  2 ||
          sym(2,0) != -4 || sym(2,1) != 2 || sym(2,2) !=  6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  0 0 -4 )\n(  0 0 2 )\n( -4 2 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Floating point scaling of the matrix
      sym.scale( 0.5 );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 5UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 1UL );
      checkNonZeros( sym, 1UL, 1UL );
      checkNonZeros( sym, 2UL, 3UL );

      if( sym(0,0) !=  0 || sym(0,1) != 0 || sym(0,2) != -2 ||
          sym(1,0) !=  0 || sym(1,1) != 0 || sym(1,2) !=  1 ||
          sym(2,0) != -2 || sym(2,1) != 1 || sym(2,2) !=  3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  0 0 -2 )\n(  0 0 1 )\n( -2 1 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major SymmetricMatrix::scale() (complex)";

      using blaze::complex;

      blaze::SymmetricMatrix< blaze::CompressedMatrix<complex<float>,blaze::columnMajor> > sym( 2UL );
      sym(0,0) = complex<float>( 1.0F, 0.0F );
      sym(0,1) = complex<float>( 2.0F, 0.0F );
      sym(1,1) = complex<float>( 4.0F, 0.0F );

      sym.scale( complex<float>( 3.0F, 0.0F ) );

      checkRows    ( sym, 2UL );
      checkColumns ( sym, 2UL );
      checkCapacity( sym, 4UL );
      checkNonZeros( sym, 4UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );

      if( sym(0,0) != complex<float>( 3.0F, 0.0F ) || sym(0,1) != complex<float>(  6.0F, 0.0F ) ||
          sym(1,0) != complex<float>( 6.0F, 0.0F ) || sym(1,1) != complex<float>( 12.0F, 0.0F ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( ( 3,0) ( 6,0)\n( 6,0) (12,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the SymmetricMatrix function call operator.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of adding and accessing elements via the function call operator
// of the SymmetricMatrix specialization. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
void SparseScalarTest::testFunctionCall()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major SymmetricMatrix::operator()";

      ST sym( 3UL );

      // Writing the element (1,1)
      sym(1,1) = 1;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 1UL );
      checkNonZeros( sym, 1UL );
      checkNonZeros( sym, 0UL, 0UL );
      checkNonZeros( sym, 1UL, 1UL );
      checkNonZeros( sym, 2UL, 0UL );

      if( sym(0,0) != 0 || sym(0,1) != 0 || sym(0,2) != 0 ||
          sym(1,0) != 0 || sym(1,1) != 1 || sym(1,2) != 0 ||
          sym(2,0) != 0 || sym(2,1) != 0 || sym(2,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 0 0 0 )\n( 0 1 0 )\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Writing the elements (2,1) and (1,2)
      sym(2,1) = 2;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 3UL );
      checkNonZeros( sym, 3UL );
      checkNonZeros( sym, 0UL, 0UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 1UL );

      if( sym(0,0) != 0 || sym(0,1) != 0 || sym(0,2) != 0 ||
          sym(1,0) != 0 || sym(1,1) != 1 || sym(1,2) != 2 ||
          sym(2,0) != 0 || sym(2,1) != 2 || sym(2,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 0 0 0 )\n( 0 1 2 )\n( 0 2 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Writing the elements (0,2) and (2,0)
      sym(0,2) = sym(1,2);

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 5UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 1UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 2UL );

      if( sym(0,0) != 0 || sym(0,1) != 0 || sym(0,2) != 2 ||
          sym(1,0) != 0 || sym(1,1) != 1 || sym(1,2) != 2 ||
          sym(2,0) != 2 || sym(2,1) != 2 || sym(2,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 0 0 2 )\n( 0 1 2 )\n( 2 2 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Adding to the elements (1,2) and (2,1)
      sym(1,2) += 3;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 5UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 1UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 2UL );

      if( sym(0,0) != 0 || sym(0,1) != 0 || sym(0,2) != 2 ||
          sym(1,0) != 0 || sym(1,1) != 1 || sym(1,2) != 5 ||
          sym(2,0) != 2 || sym(2,1) != 5 || sym(2,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 0 0 2 )\n( 0 1 5 )\n( 2 5 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Subtracting from the elements (0,1) and (1,0)
      sym(0,1) -= 4;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 7UL );
      checkNonZeros( sym, 7UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 3UL );
      checkNonZeros( sym, 2UL, 2UL );

      if( sym(0,0) !=  0 || sym(0,1) != -4 || sym(0,2) != 2 ||
          sym(1,0) != -4 || sym(1,1) !=  1 || sym(1,2) != 5 ||
          sym(2,0) !=  2 || sym(2,1) !=  5 || sym(2,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  0 -4  2 )\n( -4  1  5 )\n(  2  5  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Multiplying the element (1,1)
      sym(2,0) *= -3;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 7UL );
      checkNonZeros( sym, 7UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 3UL );
      checkNonZeros( sym, 2UL, 2UL );

      if( sym(0,0) !=  0 || sym(0,1) != -4 || sym(0,2) != -6 ||
          sym(1,0) != -4 || sym(1,1) !=  1 || sym(1,2) !=  5 ||
          sym(2,0) != -6 || sym(2,1) !=  5 || sym(2,2) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  0 -4 -6 )\n( -4  1  5 )\n( -6  5  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Dividing the elements (0,2) and (2,0)
      sym(1,0) /= 2;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 7UL );
      checkNonZeros( sym, 7UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 3UL );
      checkNonZeros( sym, 2UL, 2UL );

      if( sym(0,0) !=  0 || sym(0,1) != -2 || sym(0,2) != -6 ||
          sym(1,0) != -2 || sym(1,1) !=  1 || sym(1,2) !=  5 ||
          sym(2,0) != -6 || sym(2,1) !=  5 || sym(2,2) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  0 -2 -6 )\n( -2  1  5 )\n( -6  5  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      // Testing assignment to non-synced elements
      {
         test_ = "Row-major SymmetricMatrix::operator() (assignment to non-synced element)";

         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 6UL );
         mat(0,0) = 1;
         mat(0,2) = 7;
         mat(1,1) = 2;
         mat(2,0) = 7;
         mat(2,2) = 3;
         mat.insert( 1UL, 2UL, 0 );

         ST sym;
         sym = mat;

         sym(1,2) = 9;

         checkRows    ( sym, 3UL );
         checkColumns ( sym, 3UL );
         checkCapacity( sym, 7UL );
         checkNonZeros( sym, 7UL );
         checkNonZeros( sym, 0UL, 2UL );
         checkNonZeros( sym, 1UL, 2UL );
         checkNonZeros( sym, 2UL, 3UL );

         if( sym(0,0) != 1 || sym(0,1) != 0 || sym(0,2) != 7 ||
             sym(1,0) != 0 || sym(1,1) != 2 || sym(1,2) != 9 ||
             sym(2,0) != 7 || sym(2,1) != 9 || sym(2,2) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to non-synced element failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 1  0  7 )\n( 0  2  9 )\n( 7  9  3 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing addition assignment to non-synced elements
      {
         test_ = "Row-major SymmetricMatrix::operator() (addition assignment to non-synced element)";

         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 6UL );
         mat(0,0) = 1;
         mat(0,2) = 7;
         mat(1,1) = 2;
         mat(2,0) = 7;
         mat(2,2) = 3;
         mat.insert( 1UL, 2UL, 0 );

         ST sym;
         sym = mat;

         sym(1,2) += 9;

         checkRows    ( sym, 3UL );
         checkColumns ( sym, 3UL );
         checkCapacity( sym, 7UL );
         checkNonZeros( sym, 7UL );
         checkNonZeros( sym, 0UL, 2UL );
         checkNonZeros( sym, 1UL, 2UL );
         checkNonZeros( sym, 2UL, 3UL );

         if( sym(0,0) != 1 || sym(0,1) != 0 || sym(0,2) != 7 ||
             sym(1,0) != 0 || sym(1,1) != 2 || sym(1,2) != 9 ||
             sym(2,0) != 7 || sym(2,1) != 9 || sym(2,2) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to non-synced element failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 1  0  7 )\n( 0  2  9 )\n( 7  9  3 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing subtraction assignment to non-synced elements
      {
         test_ = "Row-major SymmetricMatrix::operator() (subtraction assignment to non-synced element)";

         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 6UL );
         mat(0,0) = 1;
         mat(0,2) = 7;
         mat(1,1) = 2;
         mat(2,0) = 7;
         mat(2,2) = 3;
         mat.insert( 1UL, 2UL, 0 );

         ST sym;
         sym = mat;

         sym(1,2) -= -9;

         checkRows    ( sym, 3UL );
         checkColumns ( sym, 3UL );
         checkCapacity( sym, 7UL );
         checkNonZeros( sym, 7UL );
         checkNonZeros( sym, 0UL, 2UL );
         checkNonZeros( sym, 1UL, 2UL );
         checkNonZeros( sym, 2UL, 3UL );

         if( sym(0,0) != 1 || sym(0,1) != 0 || sym(0,2) != 7 ||
             sym(1,0) != 0 || sym(1,1) != 2 || sym(1,2) != 9 ||
             sym(2,0) != 7 || sym(2,1) != 9 || sym(2,2) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to non-synced element failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 1  0  7 )\n( 0  2  9 )\n( 7  9  3 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing multiplication assignment to non-synced elements
      {
         test_ = "Row-major SymmetricMatrix::operator() (multiplication assignment to non-synced element)";

         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 6UL );
         mat(0,0) = 1;
         mat(0,2) = 7;
         mat(1,1) = 2;
         mat(2,0) = 7;
         mat(2,2) = 3;
         mat.insert( 1UL, 2UL, 0 );

         ST sym;
         sym = mat;

         sym(1,2) *= -9;

         checkRows    ( sym, 3UL );
         checkColumns ( sym, 3UL );
         checkCapacity( sym, 5UL );
         checkNonZeros( sym, 5UL );
         checkNonZeros( sym, 0UL, 2UL );
         checkNonZeros( sym, 1UL, 1UL );
         checkNonZeros( sym, 2UL, 2UL );

         if( sym(0,0) != 1 || sym(0,1) != 0 || sym(0,2) != 7 ||
             sym(1,0) != 0 || sym(1,1) != 2 || sym(1,2) != 0 ||
             sym(2,0) != 7 || sym(2,1) != 0 || sym(2,2) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to non-synced element failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 1  0  7 )\n( 0  2  0 )\n( 7  0  3 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing division assignment to non-synced elements
      {
         test_ = "Row-major SymmetricMatrix::operator() (division assignment to non-synced element)";

         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 6UL );
         mat(0,0) = 1;
         mat(0,2) = 7;
         mat(1,1) = 2;
         mat(2,0) = 7;
         mat(2,2) = 3;
         mat.insert( 1UL, 2UL, 0 );

         ST sym;
         sym = mat;

         sym(1,2) /= -9;

         checkRows    ( sym, 3UL );
         checkColumns ( sym, 3UL );
         checkCapacity( sym, 5UL );
         checkNonZeros( sym, 5UL );
         checkNonZeros( sym, 0UL, 2UL );
         checkNonZeros( sym, 1UL, 1UL );
         checkNonZeros( sym, 2UL, 2UL );

         if( sym(0,0) != 1 || sym(0,1) != 0 || sym(0,2) != 7 ||
             sym(1,0) != 0 || sym(1,1) != 2 || sym(1,2) != 0 ||
             sym(2,0) != 7 || sym(2,1) != 0 || sym(2,2) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to non-synced element failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 1  0  7 )\n( 0  2  0 )\n( 7  0  3 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major SymmetricMatrix::operator()";

      OST sym( 3UL );

      // Writing the element (1,1)
      sym(1,1) = 1;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 1UL );
      checkNonZeros( sym, 1UL );
      checkNonZeros( sym, 0UL, 0UL );
      checkNonZeros( sym, 1UL, 1UL );
      checkNonZeros( sym, 2UL, 0UL );

      if( sym(0,0) != 0 || sym(0,1) != 0 || sym(0,2) != 0 ||
          sym(1,0) != 0 || sym(1,1) != 1 || sym(1,2) != 0 ||
          sym(2,0) != 0 || sym(2,1) != 0 || sym(2,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 0 0 0 )\n( 0 1 0 )\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Writing the elements (2,1) and (1,2)
      sym(2,1) = 2;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 3UL );
      checkNonZeros( sym, 3UL );
      checkNonZeros( sym, 0UL, 0UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 1UL );

      if( sym(0,0) != 0 || sym(0,1) != 0 || sym(0,2) != 0 ||
          sym(1,0) != 0 || sym(1,1) != 1 || sym(1,2) != 2 ||
          sym(2,0) != 0 || sym(2,1) != 2 || sym(2,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 0 0 0 )\n( 0 1 2 )\n( 0 2 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Writing the elements (0,2) and (2,0)
      sym(0,2) = sym(1,2);

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 5UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 1UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 2UL );

      if( sym(0,0) != 0 || sym(0,1) != 0 || sym(0,2) != 2 ||
          sym(1,0) != 0 || sym(1,1) != 1 || sym(1,2) != 2 ||
          sym(2,0) != 2 || sym(2,1) != 2 || sym(2,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 0 0 2 )\n( 0 1 2 )\n( 2 2 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Adding to the elements (1,2) and (2,1)
      sym(1,2) += 3;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 5UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 1UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 2UL );

      if( sym(0,0) != 0 || sym(0,1) != 0 || sym(0,2) != 2 ||
          sym(1,0) != 0 || sym(1,1) != 1 || sym(1,2) != 5 ||
          sym(2,0) != 2 || sym(2,1) != 5 || sym(2,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 0 0 2 )\n( 0 1 5 )\n( 2 5 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Subtracting from the elements (0,1) and (1,0)
      sym(0,1) -= 4;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 7UL );
      checkNonZeros( sym, 7UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 3UL );
      checkNonZeros( sym, 2UL, 2UL );

      if( sym(0,0) !=  0 || sym(0,1) != -4 || sym(0,2) != 2 ||
          sym(1,0) != -4 || sym(1,1) !=  1 || sym(1,2) != 5 ||
          sym(2,0) !=  2 || sym(2,1) !=  5 || sym(2,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  0 -4  2 )\n( -4  1  5 )\n(  2  5  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Multiplying the element (1,1)
      sym(2,0) *= -3;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 7UL );
      checkNonZeros( sym, 7UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 3UL );
      checkNonZeros( sym, 2UL, 2UL );

      if( sym(0,0) !=  0 || sym(0,1) != -4 || sym(0,2) != -6 ||
          sym(1,0) != -4 || sym(1,1) !=  1 || sym(1,2) !=  5 ||
          sym(2,0) != -6 || sym(2,1) !=  5 || sym(2,2) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  0 -4 -6 )\n( -4  1  5 )\n( -6  5  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Dividing the elements (0,2) and (2,0)
      sym(1,0) /= 2;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 7UL );
      checkNonZeros( sym, 7UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 3UL );
      checkNonZeros( sym, 2UL, 2UL );

      if( sym(0,0) !=  0 || sym(0,1) != -2 || sym(0,2) != -6 ||
          sym(1,0) != -2 || sym(1,1) !=  1 || sym(1,2) !=  5 ||
          sym(2,0) != -6 || sym(2,1) !=  5 || sym(2,2) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  0 -2 -6 )\n( -2  1  5 )\n( -6  5  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      // Testing assignment to non-synced elements
      {
         test_ = "Column-major SymmetricMatrix::operator() (assignment to non-synced element)";

         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 6UL );
         mat(0,0) = 1;
         mat(0,2) = 7;
         mat(1,1) = 2;
         mat(2,0) = 7;
         mat(2,2) = 3;
         mat.insert( 1UL, 2UL, 0 );

         OST sym;
         sym = mat;

         sym(1,2) = 9;

         checkRows    ( sym, 3UL );
         checkColumns ( sym, 3UL );
         checkCapacity( sym, 7UL );
         checkNonZeros( sym, 7UL );
         checkNonZeros( sym, 0UL, 2UL );
         checkNonZeros( sym, 1UL, 2UL );
         checkNonZeros( sym, 2UL, 3UL );

         if( sym(0,0) != 1 || sym(0,1) != 0 || sym(0,2) != 7 ||
             sym(1,0) != 0 || sym(1,1) != 2 || sym(1,2) != 9 ||
             sym(2,0) != 7 || sym(2,1) != 9 || sym(2,2) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to non-synced element failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 1  0  7 )\n( 0  2  9 )\n( 7  9  3 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing addition assignment to non-synced elements
      {
         test_ = "Column-major SymmetricMatrix::operator() (addition assignment to non-synced element)";

         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 6UL );
         mat(0,0) = 1;
         mat(0,2) = 7;
         mat(1,1) = 2;
         mat(2,0) = 7;
         mat(2,2) = 3;
         mat.insert( 1UL, 2UL, 0 );

         OST sym;
         sym = mat;

         sym(1,2) += 9;

         checkRows    ( sym, 3UL );
         checkColumns ( sym, 3UL );
         checkCapacity( sym, 7UL );
         checkNonZeros( sym, 7UL );
         checkNonZeros( sym, 0UL, 2UL );
         checkNonZeros( sym, 1UL, 2UL );
         checkNonZeros( sym, 2UL, 3UL );

         if( sym(0,0) != 1 || sym(0,1) != 0 || sym(0,2) != 7 ||
             sym(1,0) != 0 || sym(1,1) != 2 || sym(1,2) != 9 ||
             sym(2,0) != 7 || sym(2,1) != 9 || sym(2,2) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to non-synced element failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 1  0  7 )\n( 0  2  9 )\n( 7  9  3 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing subtraction assignment to non-synced elements
      {
         test_ = "Column-major SymmetricMatrix::operator() (subtraction assignment to non-synced element)";

         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 6UL );
         mat(0,0) = 1;
         mat(0,2) = 7;
         mat(1,1) = 2;
         mat(2,0) = 7;
         mat(2,2) = 3;
         mat.insert( 1UL, 2UL, 0 );

         OST sym;
         sym = mat;

         sym(1,2) -= -9;

         checkRows    ( sym, 3UL );
         checkColumns ( sym, 3UL );
         checkCapacity( sym, 7UL );
         checkNonZeros( sym, 7UL );
         checkNonZeros( sym, 0UL, 2UL );
         checkNonZeros( sym, 1UL, 2UL );
         checkNonZeros( sym, 2UL, 3UL );

         if( sym(0,0) != 1 || sym(0,1) != 0 || sym(0,2) != 7 ||
             sym(1,0) != 0 || sym(1,1) != 2 || sym(1,2) != 9 ||
             sym(2,0) != 7 || sym(2,1) != 9 || sym(2,2) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to non-synced element failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 1  0  7 )\n( 0  2  9 )\n( 7  9  3 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing multiplication assignment to non-synced elements
      {
         test_ = "Column-major SymmetricMatrix::operator() (multiplication assignment to non-synced element)";

         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 6UL );
         mat(0,0) = 1;
         mat(0,2) = 7;
         mat(1,1) = 2;
         mat(2,0) = 7;
         mat(2,2) = 3;
         mat.insert( 1UL, 2UL, 0 );

         OST sym;
         sym = mat;

         sym(1,2) *= -9;

         checkRows    ( sym, 3UL );
         checkColumns ( sym, 3UL );
         checkCapacity( sym, 5UL );
         checkNonZeros( sym, 5UL );
         checkNonZeros( sym, 0UL, 2UL );
         checkNonZeros( sym, 1UL, 1UL );
         checkNonZeros( sym, 2UL, 2UL );

         if( sym(0,0) != 1 || sym(0,1) != 0 || sym(0,2) != 7 ||
             sym(1,0) != 0 || sym(1,1) != 2 || sym(1,2) != 0 ||
             sym(2,0) != 7 || sym(2,1) != 0 || sym(2,2) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to non-synced element failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 1  0  7 )\n( 0  2  0 )\n( 7  0  3 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing division assignment to non-synced elements
      {
         test_ = "Column-major SymmetricMatrix::operator() (division assignment to non-synced element)";

         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 6UL );
         mat(0,0) = 1;
         mat(0,2) = 7;
         mat(1,1) = 2;
         mat(2,0) = 7;
         mat(2,2) = 3;
         mat.insert( 1UL, 2UL, 0 );

         OST sym;
         sym = mat;

         sym(1,2) /= -9;

         checkRows    ( sym, 3UL );
         checkColumns ( sym, 3UL );
         checkCapacity( sym, 5UL );
         checkNonZeros( sym, 5UL );
         checkNonZeros( sym, 0UL, 2UL );
         checkNonZeros( sym, 1UL, 1UL );
         checkNonZeros( sym, 2UL, 2UL );

         if( sym(0,0) != 1 || sym(0,1) != 0 || sym(0,2) != 7 ||
             sym(1,0) != 0 || sym(1,1) != 2 || sym(1,2) != 0 ||
             sym(2,0) != 7 || sym(2,1) != 0 || sym(2,2) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to non-synced element failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 1  0  7 )\n( 0  2  0 )\n( 7  0  3 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the SymmetricMatrix iterator implementation.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the iterator implementation of the SymmetricMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseScalarTest::testIterator()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      using Iterator      = ST::Iterator;
      using ConstIterator = ST::ConstIterator;

      ST sym( 3UL );
      sym(0,1) =  1;
      sym(1,2) = -2;
      sym(2,2) =  3;

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

         ConstIterator it( begin( sym, 1UL ) );

         if( it == end( sym, 1UL ) || it->value() != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Failed iterator conversion detected\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Counting the number of elements in 0th row via Iterator (end-begin)
      {
         test_ = "Row-major Iterator subtraction (end-begin)";

         const ptrdiff_t number( end( sym, 0UL ) - begin( sym, 0UL ) );

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

         const ptrdiff_t number( cend( sym, 1UL ) - cbegin( sym, 1UL ) );

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

         ConstIterator it ( cbegin( sym, 2UL ) );
         ConstIterator end( cend( sym, 2UL ) );

         if( it == end || it->value() != -2 ) {
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

      // Testing assignment via Iterator
      {
         test_ = "Row-major assignment via Iterator";

         int value = 7;

         for( Iterator it=begin( sym, 2UL ); it!=end( sym, 2UL ); ++it ) {
            *it = value++;
         }

         if( sym(0,0) != 0 || sym(0,1) != 1 || sym(0,2) != 0 ||
             sym(1,0) != 1 || sym(1,1) != 0 || sym(1,2) != 7 ||
             sym(2,0) != 0 || sym(2,1) != 7 || sym(2,2) != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 0 1 0 )\n( 1 0 7 )\n( 0 7 8 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing addition assignment via Iterator
      {
         test_ = "Row-major addition assignment via Iterator";

         int value = 4;

         for( Iterator it=begin( sym, 1UL ); it!=end( sym, 1UL ); ++it ) {
            *it += value++;
         }

         if( sym(0,0) != 0 || sym(0,1) !=  5 || sym(0,2) !=  0 ||
             sym(1,0) != 5 || sym(1,1) !=  0 || sym(1,2) != 12 ||
             sym(2,0) != 0 || sym(2,1) != 12 || sym(2,2) !=  8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 0  5  0 )\n( 5  0 12 )\n( 0 12  8 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing subtraction assignment via Iterator
      {
         test_ = "Row-major subtraction assignment via Iterator";

         int value = 4;

         for( Iterator it=begin( sym, 1UL ); it!=end( sym, 1UL ); ++it ) {
            *it -= value++;
         }

         if( sym(0,0) != 0 || sym(0,1) != 1 || sym(0,2) != 0 ||
             sym(1,0) != 1 || sym(1,1) != 0 || sym(1,2) != 7 ||
             sym(2,0) != 0 || sym(2,1) != 7 || sym(2,2) != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 0 1 0 )\n( 1 0 7 )\n( 0 7 8 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing multiplication assignment via Iterator
      {
         test_ = "Row-major multiplication assignment via Iterator";

         for( Iterator it=begin( sym, 1UL ); it!=end( sym, 1UL ); ++it ) {
            *it *= 2;
         }

         if( sym(0,0) != 0 || sym(0,1) !=  2 || sym(0,2) !=  0 ||
             sym(1,0) != 2 || sym(1,1) !=  0 || sym(1,2) != 14 ||
             sym(2,0) != 0 || sym(2,1) != 14 || sym(2,2) !=  8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 0  2  0 )\n( 1  0 14 )\n( 0 14  8 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing division assignment via Iterator
      {
         test_ = "Row-major division assignment via Iterator";

         for( Iterator it=begin( sym, 1UL ); it!=end( sym, 1UL ); ++it ) {
            *it /= 2;
         }

         if( sym(0,0) != 0 || sym(0,1) != 1 || sym(0,2) != 0 ||
             sym(1,0) != 1 || sym(1,1) != 0 || sym(1,2) != 7 ||
             sym(2,0) != 0 || sym(2,1) != 7 || sym(2,2) != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 0 1 0 )\n( 1 0 7 )\n( 0 7 8 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   {
      using Iterator = ST::Iterator;

      // Testing assignment to via Iterator non-synced elements
      {
         test_ = "Row-major assignment via Iterator to non-synced elements";

         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 6UL );
         mat(0,0) = 1;
         mat(0,2) = 7;
         mat(1,1) = 2;
         mat(2,0) = 7;
         mat(2,2) = 3;
         mat.insert( 1UL, 2UL, 0 );

         ST sym;
         sym = mat;

         Iterator it = sym.begin( 1UL );
         ++it;
         it->value() = 9;

         checkRows    ( sym, 3UL );
         checkColumns ( sym, 3UL );
         checkCapacity( sym, 7UL );
         checkNonZeros( sym, 7UL );
         checkNonZeros( sym, 0UL, 2UL );
         checkNonZeros( sym, 1UL, 2UL );
         checkNonZeros( sym, 2UL, 3UL );

         if( sym(0,0) != 1 || sym(0,1) != 0 || sym(0,2) != 7 ||
             sym(1,0) != 0 || sym(1,1) != 2 || sym(1,2) != 9 ||
             sym(2,0) != 7 || sym(2,1) != 9 || sym(2,2) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to non-synced element failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 1  0  7 )\n( 0  2  9 )\n( 7  9  3 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing addition assignment to via Iterator non-synced elements
      {
         test_ = "Row-major addition assignment via Iterator to non-synced elements";

         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 6UL );
         mat(0,0) = 1;
         mat(0,2) = 7;
         mat(1,1) = 2;
         mat(2,0) = 7;
         mat(2,2) = 3;
         mat.insert( 1UL, 2UL, 0 );

         ST sym;
         sym = mat;

         Iterator it = sym.begin( 1UL );
         ++it;
         it->value() += 9;

         checkRows    ( sym, 3UL );
         checkColumns ( sym, 3UL );
         checkCapacity( sym, 7UL );
         checkNonZeros( sym, 7UL );
         checkNonZeros( sym, 0UL, 2UL );
         checkNonZeros( sym, 1UL, 2UL );
         checkNonZeros( sym, 2UL, 3UL );

         if( sym(0,0) != 1 || sym(0,1) != 0 || sym(0,2) != 7 ||
             sym(1,0) != 0 || sym(1,1) != 2 || sym(1,2) != 9 ||
             sym(2,0) != 7 || sym(2,1) != 9 || sym(2,2) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to non-synced element failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 1  0  7 )\n( 0  2  9 )\n( 7  9  3 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing subtraction assignment to via Iterator non-synced elements
      {
         test_ = "Row-major subtraction assignment via Iterator to non-synced elements";

         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 6UL );
         mat(0,0) = 1;
         mat(0,2) = 7;
         mat(1,1) = 2;
         mat(2,0) = 7;
         mat(2,2) = 3;
         mat.insert( 1UL, 2UL, 0 );

         ST sym;
         sym = mat;

         Iterator it = sym.begin( 1UL );
         ++it;
         it->value() -= -9;

         checkRows    ( sym, 3UL );
         checkColumns ( sym, 3UL );
         checkCapacity( sym, 7UL );
         checkNonZeros( sym, 7UL );
         checkNonZeros( sym, 0UL, 2UL );
         checkNonZeros( sym, 1UL, 2UL );
         checkNonZeros( sym, 2UL, 3UL );

         if( sym(0,0) != 1 || sym(0,1) != 0 || sym(0,2) != 7 ||
             sym(1,0) != 0 || sym(1,1) != 2 || sym(1,2) != 9 ||
             sym(2,0) != 7 || sym(2,1) != 9 || sym(2,2) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to non-synced element failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 1  0  7 )\n( 0  2  9 )\n( 7  9  3 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing multiplication assignment via Iterator to non-synced elements
      {
         test_ = "Row-major multiplication assignment via Iterator to non-synced elements";

         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 6UL );
         mat(0,0) = 1;
         mat(0,2) = 7;
         mat(1,1) = 2;
         mat(2,0) = 7;
         mat(2,2) = 3;
         mat.insert( 1UL, 2UL, 0 );

         ST sym;
         sym = mat;

         Iterator it = sym.begin( 1UL );
         ++it;
         it->value() *= 9;

         checkRows    ( sym, 3UL );
         checkColumns ( sym, 3UL );
         checkCapacity( sym, 6UL );
         checkNonZeros( sym, 6UL );
         checkNonZeros( sym, 0UL, 2UL );
         checkNonZeros( sym, 1UL, 2UL );
         checkNonZeros( sym, 2UL, 2UL );

         if( sym(0,0) != 1 || sym(0,1) != 0 || sym(0,2) != 7 ||
             sym(1,0) != 0 || sym(1,1) != 2 || sym(1,2) != 0 ||
             sym(2,0) != 7 || sym(2,1) != 0 || sym(2,2) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to non-synced element failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 1  0  7 )\n( 0  2  0 )\n( 7  0  3 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing division assignment to via Iterator non-synced elements
      {
         test_ = "Row-major division assignment to via Iterator non-synced elements";

         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 6UL );
         mat(0,0) = 1;
         mat(0,2) = 7;
         mat(1,1) = 2;
         mat(2,0) = 7;
         mat(2,2) = 3;
         mat.insert( 1UL, 2UL, 0 );

         ST sym;
         sym = mat;

         Iterator it = sym.begin( 1UL );
         ++it;
         it->value() /= 9;

         checkRows    ( sym, 3UL );
         checkColumns ( sym, 3UL );
         checkCapacity( sym, 6UL );
         checkNonZeros( sym, 6UL );
         checkNonZeros( sym, 0UL, 2UL );
         checkNonZeros( sym, 1UL, 2UL );
         checkNonZeros( sym, 2UL, 2UL );

         if( sym(0,0) != 1 || sym(0,1) != 0 || sym(0,2) != 7 ||
             sym(1,0) != 0 || sym(1,1) != 2 || sym(1,2) != 0 ||
             sym(2,0) != 7 || sym(2,1) != 0 || sym(2,2) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to non-synced element failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 1  0  7 )\n( 0  2  0 )\n( 7  0  3 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      using Iterator      = OST::Iterator;
      using ConstIterator = OST::ConstIterator;

      OST sym( 3UL );
      sym(0,1) =  1;
      sym(1,2) = -2;
      sym(2,2) =  3;

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

         ConstIterator it( begin( sym, 1UL ) );

         if( it == end( sym, 1UL ) || it->value() != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Failed iterator conversion detected\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Counting the number of elements in 0th row via Iterator (end-begin)
      {
         test_ = "Row-major Iterator subtraction (end-begin)";

         const ptrdiff_t number( end( sym, 0UL ) - begin( sym, 0UL ) );

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

         const ptrdiff_t number( cend( sym, 1UL ) - cbegin( sym, 1UL ) );

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

         ConstIterator it ( cbegin( sym, 2UL ) );
         ConstIterator end( cend( sym, 2UL ) );

         if( it == end || it->value() != -2 ) {
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

      // Testing assignment via Iterator
      {
         test_ = "Row-major assignment via Iterator";

         int value = 7;

         for( Iterator it=begin( sym, 2UL ); it!=end( sym, 2UL ); ++it ) {
            *it = value++;
         }

         if( sym(0,0) != 0 || sym(0,1) != 1 || sym(0,2) != 0 ||
             sym(1,0) != 1 || sym(1,1) != 0 || sym(1,2) != 7 ||
             sym(2,0) != 0 || sym(2,1) != 7 || sym(2,2) != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 0 1 0 )\n( 1 0 7 )\n( 0 7 8 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing addition assignment via Iterator
      {
         test_ = "Row-major addition assignment via Iterator";

         int value = 4;

         for( Iterator it=begin( sym, 1UL ); it!=end( sym, 1UL ); ++it ) {
            *it += value++;
         }

         if( sym(0,0) != 0 || sym(0,1) !=  5 || sym(0,2) !=  0 ||
             sym(1,0) != 5 || sym(1,1) !=  0 || sym(1,2) != 12 ||
             sym(2,0) != 0 || sym(2,1) != 12 || sym(2,2) !=  8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 0  5  0 )\n( 5  0 12 )\n( 0 12  8 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing subtraction assignment via Iterator
      {
         test_ = "Row-major subtraction assignment via Iterator";

         int value = 4;

         for( Iterator it=begin( sym, 1UL ); it!=end( sym, 1UL ); ++it ) {
            *it -= value++;
         }

         if( sym(0,0) != 0 || sym(0,1) != 1 || sym(0,2) != 0 ||
             sym(1,0) != 1 || sym(1,1) != 0 || sym(1,2) != 7 ||
             sym(2,0) != 0 || sym(2,1) != 7 || sym(2,2) != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 0 1 0 )\n( 1 0 7 )\n( 0 7 8 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing multiplication assignment via Iterator
      {
         test_ = "Row-major multiplication assignment via Iterator";

         for( Iterator it=begin( sym, 1UL ); it!=end( sym, 1UL ); ++it ) {
            *it *= 2;
         }

         if( sym(0,0) != 0 || sym(0,1) !=  2 || sym(0,2) !=  0 ||
             sym(1,0) != 2 || sym(1,1) !=  0 || sym(1,2) != 14 ||
             sym(2,0) != 0 || sym(2,1) != 14 || sym(2,2) !=  8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 0  2  0 )\n( 1  0 14 )\n( 0 14  8 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing division assignment via Iterator
      {
         test_ = "Row-major division assignment via Iterator";

         for( Iterator it=begin( sym, 1UL ); it!=end( sym, 1UL ); ++it ) {
            *it /= 2;
         }

         if( sym(0,0) != 0 || sym(0,1) != 1 || sym(0,2) != 0 ||
             sym(1,0) != 1 || sym(1,1) != 0 || sym(1,2) != 7 ||
             sym(2,0) != 0 || sym(2,1) != 7 || sym(2,2) != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 0 1 0 )\n( 1 0 7 )\n( 0 7 8 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   {
      using Iterator = OST::Iterator;

      // Testing assignment to via Iterator non-synced elements
      {
         test_ = "Row-major assignment to via Iterator non-synced elements";

         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 6UL );
         mat(0,0) = 1;
         mat(0,2) = 7;
         mat(1,1) = 2;
         mat(2,0) = 7;
         mat(2,2) = 3;
         mat.insert( 1UL, 2UL, 0 );

         OST sym;
         sym = mat;

         Iterator it = sym.begin( 2UL );
         ++it;
         it->value() = 9;

         checkRows    ( sym, 3UL );
         checkColumns ( sym, 3UL );
         checkCapacity( sym, 7UL );
         checkNonZeros( sym, 7UL );
         checkNonZeros( sym, 0UL, 2UL );
         checkNonZeros( sym, 1UL, 2UL );
         checkNonZeros( sym, 2UL, 3UL );

         if( sym(0,0) != 1 || sym(0,1) != 0 || sym(0,2) != 7 ||
             sym(1,0) != 0 || sym(1,1) != 2 || sym(1,2) != 9 ||
             sym(2,0) != 7 || sym(2,1) != 9 || sym(2,2) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to non-synced element failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 1  0  7 )\n( 0  2  9 )\n( 7  9  3 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing addition assignment to via Iterator non-synced elements
      {
         test_ = "Row-major addition assignment to via Iterator non-synced elements";

         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 6UL );
         mat(0,0) = 1;
         mat(0,2) = 7;
         mat(1,1) = 2;
         mat(2,0) = 7;
         mat(2,2) = 3;
         mat.insert( 1UL, 2UL, 0 );

         OST sym;
         sym = mat;

         Iterator it = sym.begin( 2UL );
         ++it;
         it->value() += 9;

         checkRows    ( sym, 3UL );
         checkColumns ( sym, 3UL );
         checkCapacity( sym, 7UL );
         checkNonZeros( sym, 7UL );
         checkNonZeros( sym, 0UL, 2UL );
         checkNonZeros( sym, 1UL, 2UL );
         checkNonZeros( sym, 2UL, 3UL );

         if( sym(0,0) != 1 || sym(0,1) != 0 || sym(0,2) != 7 ||
             sym(1,0) != 0 || sym(1,1) != 2 || sym(1,2) != 9 ||
             sym(2,0) != 7 || sym(2,1) != 9 || sym(2,2) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to non-synced element failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 1  0  7 )\n( 0  2  9 )\n( 7  9  3 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing subtraction assignment to via Iterator non-synced elements
      {
         test_ = "Row-major subtraction assignment to via Iterator non-synced elements";

         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 6UL );
         mat(0,0) = 1;
         mat(0,2) = 7;
         mat(1,1) = 2;
         mat(2,0) = 7;
         mat(2,2) = 3;
         mat.insert( 1UL, 2UL, 0 );

         OST sym;
         sym = mat;

         Iterator it = sym.begin( 2UL );
         ++it;
         it->value() -= -9;

         checkRows    ( sym, 3UL );
         checkColumns ( sym, 3UL );
         checkCapacity( sym, 7UL );
         checkNonZeros( sym, 7UL );
         checkNonZeros( sym, 0UL, 2UL );
         checkNonZeros( sym, 1UL, 2UL );
         checkNonZeros( sym, 2UL, 3UL );

         if( sym(0,0) != 1 || sym(0,1) != 0 || sym(0,2) != 7 ||
             sym(1,0) != 0 || sym(1,1) != 2 || sym(1,2) != 9 ||
             sym(2,0) != 7 || sym(2,1) != 9 || sym(2,2) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to non-synced element failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 1  0  7 )\n( 0  2  9 )\n( 7  9  3 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing multiplication assignment to via Iterator non-synced elements
      {
         test_ = "Row-major multiplication assignment to via Iterator non-synced elements";

         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 6UL );
         mat(0,0) = 1;
         mat(0,2) = 7;
         mat(1,1) = 2;
         mat(2,0) = 7;
         mat(2,2) = 3;
         mat.insert( 1UL, 2UL, 0 );

         OST sym;
         sym = mat;

         Iterator it = sym.begin( 2UL );
         ++it;
         it->value() *= 9;

         checkRows    ( sym, 3UL );
         checkColumns ( sym, 3UL );
         checkCapacity( sym, 6UL );
         checkNonZeros( sym, 6UL );
         checkNonZeros( sym, 0UL, 2UL );
         checkNonZeros( sym, 1UL, 1UL );
         checkNonZeros( sym, 2UL, 3UL );

         if( sym(0,0) != 1 || sym(0,1) != 0 || sym(0,2) != 7 ||
             sym(1,0) != 0 || sym(1,1) != 2 || sym(1,2) != 0 ||
             sym(2,0) != 7 || sym(2,1) != 0 || sym(2,2) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to non-synced element failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 1  0  7 )\n( 0  2  0 )\n( 7  0  3 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing division assignment to via Iterator non-synced elements
      {
         test_ = "Row-major division assignment to via Iterator non-synced elements";

         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 6UL );
         mat(0,0) = 1;
         mat(0,2) = 7;
         mat(1,1) = 2;
         mat(2,0) = 7;
         mat(2,2) = 3;
         mat.insert( 1UL, 2UL, 0 );

         OST sym;
         sym = mat;

         Iterator it = sym.begin( 2UL );
         ++it;
         it->value() /= 9;

         checkRows    ( sym, 3UL );
         checkColumns ( sym, 3UL );
         checkCapacity( sym, 6UL );
         checkNonZeros( sym, 6UL );
         checkNonZeros( sym, 0UL, 2UL );
         checkNonZeros( sym, 1UL, 1UL );
         checkNonZeros( sym, 2UL, 3UL );

         if( sym(0,0) != 1 || sym(0,1) != 0 || sym(0,2) != 7 ||
             sym(1,0) != 0 || sym(1,1) != 2 || sym(1,2) != 0 ||
             sym(2,0) != 7 || sym(2,1) != 0 || sym(2,2) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to non-synced element failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 1  0  7 )\n( 0  2  0 )\n( 7  0  3 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c nonZeros() member function of the SymmetricMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c nonZeros() member function of the SymmetricMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseScalarTest::testNonZeros()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major SymmetricMatrix::nonZeros()";

      // Empty matrix
      {
         ST sym( 3UL );

         checkRows    ( sym, 3UL );
         checkColumns ( sym, 3UL );
         checkNonZeros( sym, 0UL );
         checkNonZeros( sym, 0UL, 0UL );
         checkNonZeros( sym, 1UL, 0UL );
         checkNonZeros( sym, 2UL, 0UL );

         if( sym(0,0) != 0 || sym(0,1) != 0 || sym(0,2) != 0 ||
             sym(1,0) != 0 || sym(1,1) != 0 || sym(1,2) != 0 ||
             sym(2,0) != 0 || sym(2,1) != 0 || sym(2,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 0 0 0 )\n( 0 0 0 )\n( 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Partially filled matrix
      {
         ST sym( 3UL );
         sym(0,0) =  1;
         sym(1,2) = -2;
         sym(2,0) =  0;
         sym(2,2) =  3;

         checkRows    ( sym, 3UL );
         checkColumns ( sym, 3UL );
         checkCapacity( sym, 4UL );
         checkNonZeros( sym, 4UL );
         checkNonZeros( sym, 0UL, 1UL );
         checkNonZeros( sym, 1UL, 1UL );
         checkNonZeros( sym, 2UL, 2UL );

         if( sym(0,0) != 1 || sym(0,1) !=  0 || sym(0,2) !=  0 ||
             sym(1,0) != 0 || sym(1,1) !=  0 || sym(1,2) != -2 ||
             sym(2,0) != 0 || sym(2,1) != -2 || sym(2,2) !=  3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 1  0  0 )\n( 0  0 -2 )\n( 0 -2  3 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Fully filled matrix
      {
         ST sym( 3UL );
         sym(0,0) = -1;
         sym(0,1) =  2;
         sym(0,2) = -3;
         sym(1,1) =  4;
         sym(1,2) = -5;
         sym(2,2) =  6;

         checkRows    ( sym, 3UL );
         checkColumns ( sym, 3UL );
         checkCapacity( sym, 9UL );
         checkNonZeros( sym, 9UL );
         checkNonZeros( sym, 0UL, 3UL );
         checkNonZeros( sym, 1UL, 3UL );
         checkNonZeros( sym, 2UL, 3UL );

         if( sym(0,0) != -1 || sym(0,1) !=  2 || sym(0,2) != -3 ||
             sym(1,0) !=  2 || sym(1,1) !=  4 || sym(1,2) != -5 ||
             sym(2,0) != -3 || sym(2,1) != -5 || sym(2,2) !=  6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( -1  2 -3 )\n(  2  4 -5 )\n( -3 -5  6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major SymmetricMatrix::nonZeros()";

      // Empty matrix
      {
         OST sym( 3UL );

         checkRows    ( sym, 3UL );
         checkColumns ( sym, 3UL );
         checkNonZeros( sym, 0UL );
         checkNonZeros( sym, 0UL, 0UL );
         checkNonZeros( sym, 1UL, 0UL );
         checkNonZeros( sym, 2UL, 0UL );

         if( sym(0,0) != 0 || sym(0,1) != 0 || sym(0,2) != 0 ||
             sym(1,0) != 0 || sym(1,1) != 0 || sym(1,2) != 0 ||
             sym(2,0) != 0 || sym(2,1) != 0 || sym(2,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 0 0 0 )\n( 0 0 0 )\n( 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Partially filled matrix
      {
         OST sym( 3UL );
         sym(0,0) =  1;
         sym(1,2) = -2;
         sym(2,0) =  0;
         sym(2,2) =  3;

         checkRows    ( sym, 3UL );
         checkColumns ( sym, 3UL );
         checkCapacity( sym, 4UL );
         checkNonZeros( sym, 4UL );
         checkNonZeros( sym, 0UL, 1UL );
         checkNonZeros( sym, 1UL, 1UL );
         checkNonZeros( sym, 2UL, 2UL );

         if( sym(0,0) != 1 || sym(0,1) !=  0 || sym(0,2) !=  0 ||
             sym(1,0) != 0 || sym(1,1) !=  0 || sym(1,2) != -2 ||
             sym(2,0) != 0 || sym(2,1) != -2 || sym(2,2) !=  3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 1  0  0 )\n( 0  0 -2 )\n( 0 -2  3 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Fully filled matrix
      {
         OST sym( 3UL );
         sym(0,0) = -1;
         sym(0,1) =  2;
         sym(0,2) = -3;
         sym(1,1) =  4;
         sym(1,2) = -5;
         sym(2,2) =  6;

         checkRows    ( sym, 3UL );
         checkColumns ( sym, 3UL );
         checkCapacity( sym, 9UL );
         checkNonZeros( sym, 9UL );
         checkNonZeros( sym, 0UL, 3UL );
         checkNonZeros( sym, 1UL, 3UL );
         checkNonZeros( sym, 2UL, 3UL );

         if( sym(0,0) != -1 || sym(0,1) !=  2 || sym(0,2) != -3 ||
             sym(1,0) !=  2 || sym(1,1) !=  4 || sym(1,2) != -5 ||
             sym(2,0) != -3 || sym(2,1) != -5 || sym(2,2) !=  6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( -1  2 -3 )\n(  2  4 -5 )\n( -3 -5  6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c reset() member function of the SymmetricMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c reset() member function of the SymmetricMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseScalarTest::testReset()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major SymmetricMatrix::reset()";

      // Initialization check
      ST sym( 3UL );
      sym(0,0) = 1;
      sym(0,1) = 2;
      sym(0,2) = 3;
      sym(1,1) = 4;
      sym(1,2) = 5;
      sym(2,2) = 6;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 9UL );
      checkNonZeros( sym, 0UL, 3UL );
      checkNonZeros( sym, 1UL, 3UL );
      checkNonZeros( sym, 2UL, 3UL );

      if( sym(0,0) != 1 || sym(0,1) != 2 || sym(0,2) != 3 ||
          sym(1,0) != 2 || sym(1,1) != 4 || sym(1,2) != 5 ||
          sym(2,0) != 3 || sym(2,1) != 5 || sym(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 2 4 5 )\n( 3 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting a single element
      reset( sym(0,1) );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 7UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 3UL );

      if( sym(0,0) != 1 || sym(0,1) != 0 || sym(0,2) != 3 ||
          sym(1,0) != 0 || sym(1,1) != 4 || sym(1,2) != 5 ||
          sym(2,0) != 3 || sym(2,1) != 5 || sym(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 1 0 3 )\n( 0 4 5 )\n( 3 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting row 1
      reset( sym, 1UL );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 4UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 0UL );
      checkNonZeros( sym, 2UL, 2UL );

      if( sym(0,0) != 1 || sym(0,1) != 0 || sym(0,2) != 3 ||
          sym(1,0) != 0 || sym(1,1) != 0 || sym(1,2) != 0 ||
          sym(2,0) != 3 || sym(2,1) != 0 || sym(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 1 0 3 )\n( 0 0 0 )\n( 3 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting the entire matrix
      reset( sym );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 0UL );
      checkNonZeros( sym, 0UL, 0UL );
      checkNonZeros( sym, 1UL, 0UL );
      checkNonZeros( sym, 2UL, 0UL );

      if( sym(0,0) != 0 || sym(0,1) != 0 || sym(0,2) != 0 ||
          sym(1,0) != 0 || sym(1,1) != 0 || sym(1,2) != 0 ||
          sym(2,0) != 0 || sym(2,1) != 0 || sym(2,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 0 0 0 )\n( 0 0 0 )\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major SymmetricMatrix::reset()";

      // Initialization check
      OST sym( 3UL );
      sym(0,0) = 1;
      sym(0,1) = 2;
      sym(0,2) = 3;
      sym(1,1) = 4;
      sym(1,2) = 5;
      sym(2,2) = 6;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 9UL );
      checkNonZeros( sym, 0UL, 3UL );
      checkNonZeros( sym, 1UL, 3UL );
      checkNonZeros( sym, 2UL, 3UL );

      if( sym(0,0) != 1 || sym(0,1) != 2 || sym(0,2) != 3 ||
          sym(1,0) != 2 || sym(1,1) != 4 || sym(1,2) != 5 ||
          sym(2,0) != 3 || sym(2,1) != 5 || sym(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 2 4 5 )\n( 3 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting a single element
      reset( sym(0,1) );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 7UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 3UL );

      if( sym(0,0) != 1 || sym(0,1) != 0 || sym(0,2) != 3 ||
          sym(1,0) != 0 || sym(1,1) != 4 || sym(1,2) != 5 ||
          sym(2,0) != 3 || sym(2,1) != 5 || sym(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 1 0 3 )\n( 0 4 5 )\n( 3 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting column 1
      reset( sym, 1UL );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 4UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 0UL );
      checkNonZeros( sym, 2UL, 2UL );

      if( sym(0,0) != 1 || sym(0,1) != 0 || sym(0,2) != 3 ||
          sym(1,0) != 0 || sym(1,1) != 0 || sym(1,2) != 0 ||
          sym(2,0) != 3 || sym(2,1) != 0 || sym(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 1 0 3 )\n( 0 0 0 )\n( 3 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting the entire matrix
      reset( sym );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 0UL );
      checkNonZeros( sym, 0UL, 0UL );
      checkNonZeros( sym, 1UL, 0UL );
      checkNonZeros( sym, 2UL, 0UL );

      if( sym(0,0) != 0 || sym(0,1) != 0 || sym(0,2) != 0 ||
          sym(1,0) != 0 || sym(1,1) != 0 || sym(1,2) != 0 ||
          sym(2,0) != 0 || sym(2,1) != 0 || sym(2,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 0 0 0 )\n( 0 0 0 )\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c clear() member function of the SymmetricMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c clear() member function of the SymmetricMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseScalarTest::testClear()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major SymmetricMatrix::clear()";

      // Initialization check
      ST sym( 3UL );
      sym(0,0) = 1;
      sym(0,1) = 2;
      sym(0,2) = 3;
      sym(1,1) = 4;
      sym(1,2) = 5;
      sym(2,2) = 6;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 9UL );
      checkNonZeros( sym, 0UL, 3UL );
      checkNonZeros( sym, 1UL, 3UL );
      checkNonZeros( sym, 2UL, 3UL );

      if( sym(0,0) != 1 || sym(0,1) != 2 || sym(0,2) != 3 ||
          sym(1,0) != 2 || sym(1,1) != 4 || sym(1,2) != 5 ||
          sym(2,0) != 3 || sym(2,1) != 5 || sym(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 2 4 5 )\n( 3 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing a single element
      clear( sym(0,1) );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 7UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 3UL );

      if( sym(0,0) != 1 || sym(0,1) != 0 || sym(0,2) != 3 ||
          sym(1,0) != 0 || sym(1,1) != 4 || sym(1,2) != 5 ||
          sym(2,0) != 3 || sym(2,1) != 5 || sym(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 1 0 3 )\n( 0 4 5 )\n( 3 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing the matrix
      clear( sym );

      checkRows    ( sym, 0UL );
      checkColumns ( sym, 0UL );
      checkNonZeros( sym, 0UL );
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major SymmetricMatrix::clear()";

      // Initialization check
      OST sym( 3UL );
      sym(0,0) = 1;
      sym(0,1) = 2;
      sym(0,2) = 3;
      sym(1,1) = 4;
      sym(1,2) = 5;
      sym(2,2) = 6;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 9UL );
      checkNonZeros( sym, 0UL, 3UL );
      checkNonZeros( sym, 1UL, 3UL );
      checkNonZeros( sym, 2UL, 3UL );

      if( sym(0,0) != 1 || sym(0,1) != 2 || sym(0,2) != 3 ||
          sym(1,0) != 2 || sym(1,1) != 4 || sym(1,2) != 5 ||
          sym(2,0) != 3 || sym(2,1) != 5 || sym(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 2 4 5 )\n( 3 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing a single element
      clear( sym(0,1) );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 7UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 3UL );

      if( sym(0,0) != 1 || sym(0,1) != 0 || sym(0,2) != 3 ||
          sym(1,0) != 0 || sym(1,1) != 4 || sym(1,2) != 5 ||
          sym(2,0) != 3 || sym(2,1) != 5 || sym(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 1 0 3 )\n( 0 4 5 )\n( 3 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing the matrix
      clear( sym );

      checkRows    ( sym, 0UL );
      checkColumns ( sym, 0UL );
      checkNonZeros( sym, 0UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c resize() member function of the SymmetricMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c resize() member function of the SymmetricMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseScalarTest::testResize()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major SymmetricMatrix::resize()";

      // Initialization check
      ST sym;

      checkRows    ( sym, 0UL );
      checkColumns ( sym, 0UL );
      checkNonZeros( sym, 0UL );

      // Resizing to 2x2
      sym.resize( 2UL );

      checkRows    ( sym, 2UL );
      checkColumns ( sym, 2UL );
      checkNonZeros( sym, 0UL );
      checkNonZeros( sym, 0UL, 0UL );
      checkNonZeros( sym, 1UL, 0UL );

      if( sym(0,0) != 0 || sym(0,1) != 0 ||
          sym(1,0) != 0 || sym(1,1) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 0 0 )\n( 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 4x4 and preserving the elements
      sym(0,1) = 1;
      sym(1,1) = 2;
      sym.resize( 4UL, true );

      checkRows    ( sym, 4UL );
      checkColumns ( sym, 4UL );
      checkCapacity( sym, 3UL );
      checkNonZeros( sym, 3UL );
      checkNonZeros( sym, 0UL, 1UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 0UL );
      checkNonZeros( sym, 3UL, 0UL );

      if( sym(0,0) != 0 || sym(0,1) != 1 || sym(0,2) != 0 || sym(0,3) != 0 ||
          sym(1,0) != 1 || sym(1,1) != 2 || sym(1,2) != 0 || sym(1,3) != 0 ||
          sym(2,0) != 0 || sym(2,1) != 0 || sym(2,2) != 0 || sym(2,3) != 0 ||
          sym(3,0) != 0 || sym(3,1) != 0 || sym(3,2) != 0 || sym(3,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 0 1 0 0 )\n( 1 2 0 0 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 2x2
      sym(2,2) = 3;
      sym.resize( 2UL );

      checkRows    ( sym, 2UL );
      checkColumns ( sym, 2UL );
      checkCapacity( sym, 3UL );
      checkNonZeros( sym, 3UL );
      checkNonZeros( sym, 0UL, 1UL );
      checkNonZeros( sym, 1UL, 2UL );

      if( sym(0,0) != 0 || sym(0,1) != 1 ||
          sym(1,0) != 1 || sym(1,1) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 0 1 )\n( 1 2 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 0x0
      sym.resize( 0UL );

      checkRows    ( sym, 0UL );
      checkColumns ( sym, 0UL );
      checkNonZeros( sym, 0UL );
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major SymmetricMatrix::resize()";

      // Initialization check
      OST sym;

      checkRows    ( sym, 0UL );
      checkColumns ( sym, 0UL );
      checkNonZeros( sym, 0UL );

      // Resizing to 2x2
      sym.resize( 2UL );

      checkRows    ( sym, 2UL );
      checkColumns ( sym, 2UL );
      checkNonZeros( sym, 0UL );
      checkNonZeros( sym, 0UL, 0UL );
      checkNonZeros( sym, 1UL, 0UL );

      if( sym(0,0) != 0 || sym(0,1) != 0 ||
          sym(1,0) != 0 || sym(1,1) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 0 0 )\n( 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 4x4 and preserving the elements
      sym(0,1) = 1;
      sym(1,1) = 2;
      sym.resize( 4UL );

      checkRows    ( sym, 4UL );
      checkColumns ( sym, 4UL );
      checkCapacity( sym, 3UL );
      checkNonZeros( sym, 3UL );
      checkNonZeros( sym, 0UL, 1UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 0UL );
      checkNonZeros( sym, 3UL, 0UL );

      if( sym(0,0) != 0 || sym(0,1) != 1 || sym(0,2) != 0 || sym(0,3) != 0 ||
          sym(1,0) != 1 || sym(1,1) != 2 || sym(1,2) != 0 || sym(1,3) != 0 ||
          sym(2,0) != 0 || sym(2,1) != 0 || sym(2,2) != 0 || sym(2,3) != 0 ||
          sym(3,0) != 0 || sym(3,1) != 0 || sym(3,2) != 0 || sym(3,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 0 1 0 0 )\n( 1 2 0 0 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 2x2
      sym(2,2) = 2;
      sym.resize( 2UL );

      checkRows    ( sym, 2UL );
      checkColumns ( sym, 2UL );
      checkCapacity( sym, 3UL );
      checkNonZeros( sym, 3UL );
      checkNonZeros( sym, 0UL, 1UL );
      checkNonZeros( sym, 1UL, 2UL );

      if( sym(0,0) != 0 || sym(0,1) != 1 ||
          sym(1,0) != 1 || sym(1,1) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 0 1 )\n( 1 2 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 0x0
      sym.resize( 0UL );

      checkRows    ( sym, 0UL );
      checkColumns ( sym, 0UL );
      checkNonZeros( sym, 0UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c reserve() member function of the SymmetricMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c reserve() member function of the SymmetricMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseScalarTest::testReserve()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major SymmetricMatrix::reserve()";

      // Initialization check
      ST sym;

      checkRows    ( sym, 0UL );
      checkColumns ( sym, 0UL );
      checkNonZeros( sym, 0UL );

      // Increasing the capacity of the matrix
      sym.reserve( 10UL );

      checkRows    ( sym,  0UL );
      checkColumns ( sym,  0UL );
      checkCapacity( sym, 10UL );
      checkNonZeros( sym,  0UL );

      // Further increasing the capacity of the matrix
      sym.reserve( 20UL );

      checkRows    ( sym,  0UL );
      checkColumns ( sym,  0UL );
      checkCapacity( sym, 20UL );
      checkNonZeros( sym,  0UL );
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major SymmetricMatrix::reserve()";

      // Initialization check
      OST sym;

      checkRows    ( sym, 0UL );
      checkColumns ( sym, 0UL );
      checkNonZeros( sym, 0UL );

      // Increasing the capacity of the matrix
      sym.reserve( 10UL );

      checkRows    ( sym,  0UL );
      checkColumns ( sym,  0UL );
      checkCapacity( sym, 10UL );
      checkNonZeros( sym,  0UL );

      // Further increasing the capacity of the matrix
      sym.reserve( 20UL );

      checkRows    ( sym,  0UL );
      checkColumns ( sym,  0UL );
      checkCapacity( sym, 20UL );
      checkNonZeros( sym,  0UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c trim() member function of the SymmetricMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c trim() member function of the SymmetricMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseScalarTest::testTrim()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major SymmetricMatrix::trim()";

      // Initialization check
      ST sym( 3UL );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 0UL );

      // Increasing the row capacity of the matrix
      sym.reserve( 0UL, 10UL );
      sym.reserve( 1UL, 15UL );
      sym.reserve( 2UL, 20UL );

      checkRows    ( sym,  3UL );
      checkColumns ( sym,  3UL );
      checkCapacity( sym, 45UL );
      checkCapacity( sym,  0UL, 10UL );
      checkCapacity( sym,  1UL, 15UL );
      checkCapacity( sym,  2UL, 20UL );

      // Trimming the matrix
      sym.trim();

      checkRows    ( sym,  3UL );
      checkColumns ( sym,  3UL );
      checkCapacity( sym, 45UL );
      checkCapacity( sym,  0UL, 0UL );
      checkCapacity( sym,  1UL, 0UL );
      checkCapacity( sym,  2UL, 0UL );
   }

   {
      test_ = "Row-major SymmetricMatrix::trim( size_t )";

      // Initialization check
      ST sym( 3UL, 3UL );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 0UL );

      // Increasing the row capacity of the matrix
      sym.reserve( 0UL, 10UL );
      sym.reserve( 1UL, 15UL );
      sym.reserve( 2UL, 20UL );

      checkRows    ( sym,  3UL );
      checkColumns ( sym,  3UL );
      checkCapacity( sym, 45UL );
      checkCapacity( sym,  0UL, 10UL );
      checkCapacity( sym,  1UL, 15UL );
      checkCapacity( sym,  2UL, 20UL );

      // Trimming the 0th row
      sym.trim( 0UL );

      checkRows    ( sym,  3UL );
      checkColumns ( sym,  3UL );
      checkCapacity( sym, 45UL );
      checkCapacity( sym,  0UL,  0UL );
      checkCapacity( sym,  1UL, 25UL );
      checkCapacity( sym,  2UL, 20UL );

      // Trimming the 1st row
      sym.trim( 1UL );

      checkRows    ( sym,  3UL );
      checkColumns ( sym,  3UL );
      checkCapacity( sym, 45UL );
      checkCapacity( sym,  0UL,  0UL );
      checkCapacity( sym,  1UL,  0UL );
      checkCapacity( sym,  2UL, 45UL );

      // Trimming the 2nd row
      sym.trim( 2UL );

      checkRows    ( sym,  3UL );
      checkColumns ( sym,  3UL );
      checkCapacity( sym, 45UL );
      checkCapacity( sym,  0UL, 0UL );
      checkCapacity( sym,  1UL, 0UL );
      checkCapacity( sym,  2UL, 0UL );
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major SymmetricMatrix::trim()";

      // Initialization check
      OST sym( 3UL );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 0UL );

      // Increasing the row capacity of the matrix
      sym.reserve( 0UL, 10UL );
      sym.reserve( 1UL, 15UL );
      sym.reserve( 2UL, 20UL );

      checkRows    ( sym,  3UL );
      checkColumns ( sym,  3UL );
      checkCapacity( sym, 45UL );
      checkCapacity( sym,  0UL, 10UL );
      checkCapacity( sym,  1UL, 15UL );
      checkCapacity( sym,  2UL, 20UL );

      // Trimming the matrix
      sym.trim();

      checkRows    ( sym,  3UL );
      checkColumns ( sym,  3UL );
      checkCapacity( sym, 45UL );
      checkCapacity( sym,  0UL, 0UL );
      checkCapacity( sym,  1UL, 0UL );
      checkCapacity( sym,  2UL, 0UL );
   }

   {
      test_ = "Column-major SymmetricMatrix::trim( size_t )";

      // Initialization check
      OST sym( 3UL, 3UL );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 0UL );

      // Increasing the column capacity of the matrix
      sym.reserve( 0UL, 10UL );
      sym.reserve( 1UL, 15UL );
      sym.reserve( 2UL, 20UL );

      checkRows    ( sym,  3UL );
      checkColumns ( sym,  3UL );
      checkCapacity( sym, 45UL );
      checkCapacity( sym,  0UL, 10UL );
      checkCapacity( sym,  1UL, 15UL );
      checkCapacity( sym,  2UL, 20UL );

      // Trimming the 0th column
      sym.trim( 0UL );

      checkRows    ( sym,  3UL );
      checkColumns ( sym,  3UL );
      checkCapacity( sym, 45UL );
      checkCapacity( sym,  0UL,  0UL );
      checkCapacity( sym,  1UL, 25UL );
      checkCapacity( sym,  2UL, 20UL );

      // Trimming the 1st column
      sym.trim( 1UL );

      checkRows    ( sym,  3UL );
      checkColumns ( sym,  3UL );
      checkCapacity( sym, 45UL );
      checkCapacity( sym,  0UL,  0UL );
      checkCapacity( sym,  1UL,  0UL );
      checkCapacity( sym,  2UL, 45UL );

      // Trimming the 2nd column
      sym.trim( 2UL );

      checkRows    ( sym,  3UL );
      checkColumns ( sym,  3UL );
      checkCapacity( sym, 45UL );
      checkCapacity( sym,  0UL, 0UL );
      checkCapacity( sym,  1UL, 0UL );
      checkCapacity( sym,  2UL, 0UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c shrinkToFit() member function of the SymmetricMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c shrinkToFit() member function of the SymmetricMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseScalarTest::testShrinkToFit()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major SymmetricMatrix::shrinkToFit()";

      // Shrinking a matrix without excessive capacity
      {
         ST sym( 3UL, 5UL );
         sym(0,0) = 1;
         sym(0,2) = 2;
         sym(1,1) = 3;
         sym(2,2) = 4;

         sym.shrinkToFit();

         checkRows    ( sym, 3UL );
         checkColumns ( sym, 3UL );
         checkCapacity( sym, 5UL );
         checkNonZeros( sym, 5UL );
         checkNonZeros( sym, 0UL, 2UL );
         checkNonZeros( sym, 1UL, 1UL );
         checkNonZeros( sym, 2UL, 2UL );

         if( sym.capacity() != sym.nonZeros() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Shrinking the matrix failed\n"
                << " Details:\n"
                << "   Capacity         : " << sym.capacity() << "\n"
                << "   Expected capacity: " << sym.nonZeros() << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) != 1 || sym(0,1) != 0 || sym(0,2) != 2 ||
             sym(1,0) != 0 || sym(1,1) != 3 || sym(1,2) != 0 ||
             sym(2,0) != 2 || sym(2,1) != 0 || sym(2,2) != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Shrinking the matrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 1 0 2 )\n( 0 3 0 )\n( 2 0 4 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Shrinking a matrix with excessive capacity
      {
         ST sym{ 3UL, 100UL };
         sym(0,0) = 1;
         sym(0,2) = 2;
         sym(1,1) = 3;
         sym(2,2) = 4;

         sym.shrinkToFit();

         checkRows    ( sym, 3UL );
         checkColumns ( sym, 3UL );
         checkCapacity( sym, 5UL );
         checkNonZeros( sym, 5UL );
         checkNonZeros( sym, 0UL, 2UL );
         checkNonZeros( sym, 1UL, 1UL );
         checkNonZeros( sym, 2UL, 2UL );

         if( sym.capacity() != sym.nonZeros() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Shrinking the matrix failed\n"
                << " Details:\n"
                << "   Capacity         : " << sym.capacity() << "\n"
                << "   Expected capacity: " << sym.nonZeros() << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) != 1 || sym(0,1) != 0 || sym(0,2) != 2 ||
             sym(1,0) != 0 || sym(1,1) != 3 || sym(1,2) != 0 ||
             sym(2,0) != 2 || sym(2,1) != 0 || sym(2,2) != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Shrinking the matrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 1 0 2 )\n( 0 3 0 )\n( 2 0 4 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major SymmetricMatrix::shrinkToFit()";

      // Shrinking a matrix without excessive capacity
      {
         OST sym( 3UL, 5UL );
         sym(0,0) = 1;
         sym(0,2) = 2;
         sym(1,1) = 3;
         sym(2,2) = 4;

         sym.shrinkToFit();

         checkRows    ( sym, 3UL );
         checkColumns ( sym, 3UL );
         checkCapacity( sym, 5UL );
         checkNonZeros( sym, 5UL );
         checkNonZeros( sym, 0UL, 2UL );
         checkNonZeros( sym, 1UL, 1UL );
         checkNonZeros( sym, 2UL, 2UL );

         if( sym.capacity() != sym.nonZeros() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Shrinking the matrix failed\n"
                << " Details:\n"
                << "   Capacity         : " << sym.capacity() << "\n"
                << "   Expected capacity: " << sym.nonZeros() << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) != 1 || sym(0,1) != 0 || sym(0,2) != 2 ||
             sym(1,0) != 0 || sym(1,1) != 3 || sym(1,2) != 0 ||
             sym(2,0) != 2 || sym(2,1) != 0 || sym(2,2) != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Shrinking the matrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 1 0 2 )\n( 0 3 0 )\n( 2 0 4 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Shrinking a matrix with excessive capacity
      {
         OST sym{ 3UL, 100UL };
         sym(0,0) = 1;
         sym(0,2) = 2;
         sym(1,1) = 3;
         sym(2,2) = 4;

         sym.shrinkToFit();

         checkRows    ( sym, 3UL );
         checkColumns ( sym, 3UL );
         checkCapacity( sym, 5UL );
         checkNonZeros( sym, 5UL );
         checkNonZeros( sym, 0UL, 2UL );
         checkNonZeros( sym, 1UL, 1UL );
         checkNonZeros( sym, 2UL, 2UL );

         if( sym.capacity() != sym.nonZeros() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Shrinking the matrix failed\n"
                << " Details:\n"
                << "   Capacity         : " << sym.capacity() << "\n"
                << "   Expected capacity: " << sym.nonZeros() << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) != 1 || sym(0,1) != 0 || sym(0,2) != 2 ||
             sym(1,0) != 0 || sym(1,1) != 3 || sym(1,2) != 0 ||
             sym(2,0) != 2 || sym(2,1) != 0 || sym(2,2) != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Shrinking the matrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 1 0 2 )\n( 0 3 0 )\n( 2 0 4 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c swap() functionality of the SymmetricMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c swap() function of the SymmetricMatrix specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseScalarTest::testSwap()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major SymmetricMatrix swap";

      ST sym1( 2UL );
      sym1(0,0) = 1;
      sym1(0,1) = 2;
      sym1(1,1) = 3;

      ST sym2( 2UL );
      sym2(0,0) = 4;
      sym2(0,1) = 5;

      swap( sym1, sym2 );

      checkRows    ( sym1, 2UL );
      checkColumns ( sym1, 2UL );
      checkCapacity( sym1, 4UL );
      checkNonZeros( sym1, 3UL );
      checkNonZeros( sym1, 0UL, 2UL );
      checkNonZeros( sym1, 1UL, 1UL );

      if( sym1(0,0) != 4 || sym1(0,1) != 5 || sym1(1,0) != 5 || sym1(1,1) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the first matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sym1 << "\n"
             << "   Expected result:\n( 4 5 )\n( 5 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      checkRows    ( sym2, 2UL );
      checkColumns ( sym2, 2UL );
      checkCapacity( sym2, 4UL );
      checkNonZeros( sym2, 4UL );
      checkNonZeros( sym2, 0UL, 2UL );
      checkNonZeros( sym2, 1UL, 2UL );

      if( sym2(0,0) != 1 || sym2(0,1) != 2 || sym2(1,0) != 2 || sym2(1,1) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the second matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( 1 2 )\n( 2 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major SymmetricMatrix swap";

      OST sym1( 2UL );
      sym1(0,0) = 1;
      sym1(0,1) = 2;
      sym1(1,1) = 3;

      OST sym2( 2UL );
      sym2(0,0) = 4;
      sym2(0,1) = 5;

      swap( sym1, sym2 );

      checkRows    ( sym1, 2UL );
      checkColumns ( sym1, 2UL );
      checkCapacity( sym1, 4UL );
      checkNonZeros( sym1, 3UL );
      checkNonZeros( sym1, 0UL, 2UL );
      checkNonZeros( sym1, 1UL, 1UL );

      if( sym1(0,0) != 4 || sym1(0,1) != 5 || sym1(1,0) != 5 || sym1(1,1) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the first matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sym1 << "\n"
             << "   Expected result:\n( 4 5 )\n( 5 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      checkRows    ( sym2, 2UL );
      checkColumns ( sym2, 2UL );
      checkCapacity( sym2, 4UL );
      checkNonZeros( sym2, 4UL );
      checkNonZeros( sym2, 0UL, 2UL );
      checkNonZeros( sym2, 1UL, 2UL );

      if( sym2(0,0) != 1 || sym2(0,1) != 2 || sym2(1,0) != 2 || sym2(1,1) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the second matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( 1 2 )\n( 2 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c set() member function of the SymmetricMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c set() member function of the SymmetricMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseScalarTest::testSet()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major SymmetricMatrix::set()";

      using Iterator = ST::Iterator;

      // Initialization check
      ST sym( 4UL );

      checkRows    ( sym, 4UL );
      checkColumns ( sym, 4UL );
      checkNonZeros( sym, 0UL );
      checkNonZeros( sym, 0UL, 0UL );
      checkNonZeros( sym, 1UL, 0UL );
      checkNonZeros( sym, 2UL, 0UL );
      checkNonZeros( sym, 3UL, 0UL );

      // Setting a non-zero element
      {
         Iterator pos = sym.set( 2UL, 1UL, 1 );

         checkRows    ( sym, 4UL );
         checkColumns ( sym, 4UL );
         checkCapacity( sym, 2UL );
         checkNonZeros( sym, 2UL );
         checkNonZeros( sym, 0UL, 0UL );
         checkNonZeros( sym, 1UL, 1UL );
         checkNonZeros( sym, 2UL, 1UL );
         checkNonZeros( sym, 3UL, 0UL );

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

         if( sym(1,2) != 1 || sym(2,1) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 1 0 )\n( 0 1 0 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setting a second non-zero element
      {
         Iterator pos = sym.set( 2UL, 2UL, 2 );

         checkRows    ( sym, 4UL );
         checkColumns ( sym, 4UL );
         checkCapacity( sym, 3UL );
         checkNonZeros( sym, 3UL );
         checkNonZeros( sym, 0UL, 0UL );
         checkNonZeros( sym, 1UL, 1UL );
         checkNonZeros( sym, 2UL, 2UL );
         checkNonZeros( sym, 3UL, 0UL );

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

         if( sym(1,2) != 1 || sym(2,1) != 1 || sym(2,2) != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 1 0 )\n( 0 1 2 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setting a third non-zero element
      {
         Iterator pos = sym.set( 2UL, 0UL, 3 );

         checkRows    ( sym, 4UL );
         checkColumns ( sym, 4UL );
         checkCapacity( sym, 5UL );
         checkNonZeros( sym, 5UL );
         checkNonZeros( sym, 0UL, 1UL );
         checkNonZeros( sym, 1UL, 1UL );
         checkNonZeros( sym, 2UL, 3UL );
         checkNonZeros( sym, 3UL, 0UL );

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

         if( sym(0,2) != 3 || sym(1,2) != 1 || sym(2,0) != 3 || sym(2,1) != 1 || sym(2,2) != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 0 0 3 0 )\n( 0 0 1 0 )\n( 3 1 2 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setting an already existing element
      {
         Iterator pos = sym.set( 1UL, 2UL, 4 );

         checkRows    ( sym, 4UL );
         checkColumns ( sym, 4UL );
         checkCapacity( sym, 5UL );
         checkNonZeros( sym, 5UL );
         checkNonZeros( sym, 0UL, 1UL );
         checkNonZeros( sym, 1UL, 1UL );
         checkNonZeros( sym, 2UL, 3UL );
         checkNonZeros( sym, 3UL, 0UL );

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

         if( sym(0,2) != 3 || sym(1,2) != 4 || sym(2,0) != 3 || sym(2,1) != 4 || sym(2,2) != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 0 0 3 0 )\n( 0 0 4 0 )\n( 3 4 2 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major SymmetricMatrix::set()";

      using Iterator = OST::Iterator;

      // Initialization check
      OST sym( 4UL );

      checkRows    ( sym, 4UL );
      checkColumns ( sym, 4UL );
      checkNonZeros( sym, 0UL );
      checkNonZeros( sym, 0UL, 0UL );
      checkNonZeros( sym, 1UL, 0UL );
      checkNonZeros( sym, 2UL, 0UL );
      checkNonZeros( sym, 3UL, 0UL );

      // Setting a non-zero element
      {
         Iterator pos = sym.set( 1UL, 2UL, 1 );

         checkRows    ( sym, 4UL );
         checkColumns ( sym, 4UL );
         checkCapacity( sym, 2UL );
         checkNonZeros( sym, 2UL );
         checkNonZeros( sym, 0UL, 0UL );
         checkNonZeros( sym, 1UL, 1UL );
         checkNonZeros( sym, 2UL, 1UL );
         checkNonZeros( sym, 3UL, 0UL );

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

         if( sym(1,2) != 1 || sym(2,1) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 1 0 )\n( 0 1 0 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setting a second non-zero element
      {
         Iterator pos = sym.set( 2UL, 2UL, 2 );

         checkRows    ( sym, 4UL );
         checkColumns ( sym, 4UL );
         checkCapacity( sym, 3UL );
         checkNonZeros( sym, 3UL );
         checkNonZeros( sym, 0UL, 0UL );
         checkNonZeros( sym, 1UL, 1UL );
         checkNonZeros( sym, 2UL, 2UL );
         checkNonZeros( sym, 3UL, 0UL );

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

         if( sym(1,2) != 1 || sym(2,1) != 1 || sym(2,2) != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 1 0 )\n( 0 1 2 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setting a third non-zero element
      {
         Iterator pos = sym.set( 0UL, 2UL, 3 );

         checkRows    ( sym, 4UL );
         checkColumns ( sym, 4UL );
         checkCapacity( sym, 5UL );
         checkNonZeros( sym, 5UL );
         checkNonZeros( sym, 0UL, 1UL );
         checkNonZeros( sym, 1UL, 1UL );
         checkNonZeros( sym, 2UL, 3UL );
         checkNonZeros( sym, 3UL, 0UL );

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

         if( sym(0,2) != 3 || sym(1,2) != 1 || sym(2,0) != 3 || sym(2,1) != 1 || sym(2,2) != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 0 0 3 0 )\n( 0 0 1 0 )\n( 3 1 2 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setting an already existing element
      {
         Iterator pos = sym.set( 2UL, 1UL, 4 );

         checkRows    ( sym, 4UL );
         checkColumns ( sym, 4UL );
         checkCapacity( sym, 5UL );
         checkNonZeros( sym, 5UL );
         checkNonZeros( sym, 0UL, 1UL );
         checkNonZeros( sym, 1UL, 1UL );
         checkNonZeros( sym, 2UL, 3UL );
         checkNonZeros( sym, 3UL, 0UL );

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

         if( sym(0,2) != 3 || sym(1,2) != 4 || sym(2,0) != 3 || sym(2,1) != 4 || sym(2,2) != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 0 0 3 0 )\n( 0 0 4 0 )\n( 3 4 2 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c insert() member function of the SymmetricMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c insert() member function of the SymmetricMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseScalarTest::testInsert()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major SymmetricMatrix::insert()";

      using Iterator = ST::Iterator;

      // Initialization check
      ST sym( 4UL );

      checkRows    ( sym, 4UL );
      checkColumns ( sym, 4UL );
      checkNonZeros( sym, 0UL );
      checkNonZeros( sym, 0UL, 0UL );
      checkNonZeros( sym, 1UL, 0UL );
      checkNonZeros( sym, 2UL, 0UL );
      checkNonZeros( sym, 3UL, 0UL );

      // Inserting a non-zero element
      {
         Iterator pos = sym.insert( 2UL, 1UL, 1 );

         checkRows    ( sym, 4UL );
         checkColumns ( sym, 4UL );
         checkCapacity( sym, 2UL );
         checkNonZeros( sym, 2UL );
         checkNonZeros( sym, 0UL, 0UL );
         checkNonZeros( sym, 1UL, 1UL );
         checkNonZeros( sym, 2UL, 1UL );
         checkNonZeros( sym, 3UL, 0UL );

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

         if( sym(1,2) != 1 || sym(2,1) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 1 0 )\n( 0 1 0 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting a second non-zero element
      {
         Iterator pos = sym.insert( 2UL, 2UL, 2 );

         checkRows    ( sym, 4UL );
         checkColumns ( sym, 4UL );
         checkCapacity( sym, 3UL );
         checkNonZeros( sym, 3UL );
         checkNonZeros( sym, 0UL, 0UL );
         checkNonZeros( sym, 1UL, 1UL );
         checkNonZeros( sym, 2UL, 2UL );
         checkNonZeros( sym, 3UL, 0UL );

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

         if( sym(1,2) != 1 || sym(2,1) != 1 || sym(2,2) != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 1 0 )\n( 0 1 2 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting a third non-zero element
      {
         Iterator pos = sym.insert( 2UL, 0UL, 3 );

         checkRows    ( sym, 4UL );
         checkColumns ( sym, 4UL );
         checkCapacity( sym, 5UL );
         checkNonZeros( sym, 5UL );
         checkNonZeros( sym, 0UL, 1UL );
         checkNonZeros( sym, 1UL, 1UL );
         checkNonZeros( sym, 2UL, 3UL );
         checkNonZeros( sym, 3UL, 0UL );

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

         if( sym(0,2) != 3 || sym(1,2) != 1 || sym(2,0) != 3 || sym(2,1) != 1 || sym(2,2) != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 0 0 3 0 )\n( 0 0 1 0 )\n( 3 1 2 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to insert an already existing element
      try {
         sym.insert( 1UL, 2UL, 4 );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Inserting an existing element succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 0 0 3 0 )\n( 0 0 1 0 )\n( 3 1 2 0 )\n( 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major SymmetricMatrix::insert()";

      using Iterator = OST::Iterator;

      // Initialization check
      OST sym( 4UL );

      checkRows    ( sym, 4UL );
      checkColumns ( sym, 4UL );
      checkNonZeros( sym, 0UL );
      checkNonZeros( sym, 0UL, 0UL );
      checkNonZeros( sym, 1UL, 0UL );
      checkNonZeros( sym, 2UL, 0UL );
      checkNonZeros( sym, 3UL, 0UL );

      // Inserting a non-zero element
      {
         Iterator pos = sym.insert( 1UL, 2UL, 1 );

         checkRows    ( sym, 4UL );
         checkColumns ( sym, 4UL );
         checkCapacity( sym, 2UL );
         checkNonZeros( sym, 2UL );
         checkNonZeros( sym, 0UL, 0UL );
         checkNonZeros( sym, 1UL, 1UL );
         checkNonZeros( sym, 2UL, 1UL );
         checkNonZeros( sym, 3UL, 0UL );

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

         if( sym(1,2) != 1 || sym(2,1) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 1 0 )\n( 0 1 0 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting a second non-zero element
      {
         Iterator pos = sym.insert( 2UL, 2UL, 2 );

         checkRows    ( sym, 4UL );
         checkColumns ( sym, 4UL );
         checkCapacity( sym, 3UL );
         checkNonZeros( sym, 3UL );
         checkNonZeros( sym, 0UL, 0UL );
         checkNonZeros( sym, 1UL, 1UL );
         checkNonZeros( sym, 2UL, 2UL );
         checkNonZeros( sym, 3UL, 0UL );

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

         if( sym(1,2) != 1 || sym(2,1) != 1 || sym(2,2) != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 1 0 )\n( 0 1 2 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting a third non-zero element
      {
         Iterator pos = sym.insert( 0UL, 2UL, 3 );

         checkRows    ( sym, 4UL );
         checkColumns ( sym, 4UL );
         checkCapacity( sym, 5UL );
         checkNonZeros( sym, 5UL );
         checkNonZeros( sym, 0UL, 1UL );
         checkNonZeros( sym, 1UL, 1UL );
         checkNonZeros( sym, 2UL, 3UL );
         checkNonZeros( sym, 3UL, 0UL );

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

         if( sym(0,2) != 3 || sym(1,2) != 1 || sym(2,0) != 3 || sym(2,1) != 1 || sym(2,2) != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 0 0 3 0 )\n( 0 0 1 0 )\n( 3 1 2 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to insert an already existing element
      try {
         sym.insert( 2UL, 1UL, 4 );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Inserting an existing element succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 0 0 3 0 )\n( 0 0 1 0 )\n( 3 1 2 0 )\n( 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c append() member function of the SymmetricMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c append() member function of the SymmetricMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseScalarTest::testAppend()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major SymmetricMatrix::append()";

      // Appending with pre-allocation in each row
      {
         // Initialization check
         ST sym( 4UL, 9UL );
         sym.reserve( 0UL, 2UL );
         sym.reserve( 1UL, 2UL );
         sym.reserve( 2UL, 2UL );
         sym.reserve( 3UL, 3UL );

         checkRows    ( sym, 4UL );
         checkColumns ( sym, 4UL );
         checkCapacity( sym, 9UL );
         checkNonZeros( sym, 0UL );
         checkNonZeros( sym, 0UL, 0UL );
         checkNonZeros( sym, 1UL, 0UL );
         checkNonZeros( sym, 2UL, 0UL );
         checkNonZeros( sym, 3UL, 0UL );

         // Appending one non-zero element
         sym.append( 2UL, 1UL, 1 );

         checkRows    ( sym, 4UL );
         checkColumns ( sym, 4UL );
         checkCapacity( sym, 9UL );
         checkNonZeros( sym, 2UL );
         checkNonZeros( sym, 0UL, 0UL );
         checkNonZeros( sym, 1UL, 1UL );
         checkNonZeros( sym, 2UL, 1UL );
         checkNonZeros( sym, 3UL, 0UL );

         if( sym(1,2) != 1 || sym(2,1) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 1 0 )\n( 0 1 0 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Appending two more non-zero elements
         sym.append( 0UL, 0UL, 2 );
         sym.append( 0UL, 3UL, 3 );

         checkRows    ( sym, 4UL );
         checkColumns ( sym, 4UL );
         checkCapacity( sym, 9UL );
         checkNonZeros( sym, 5UL );
         checkNonZeros( sym, 0UL, 2UL );
         checkNonZeros( sym, 1UL, 1UL );
         checkNonZeros( sym, 2UL, 1UL );
         checkNonZeros( sym, 3UL, 1UL );

         if( sym(0,0) != 2 || sym(0,3) != 3 ||
             sym(1,2) != 1 ||
             sym(2,1) != 1 ||
             sym(3,0) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 2 0 0 3 )\n( 0 0 1 0 )\n( 0 1 0 0 )\n( 3 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Appending two more non-zero elements
         sym.append( 3UL, 1UL, 4 );
         sym.append( 3UL, 2UL, 5 );

         checkRows    ( sym, 4UL );
         checkColumns ( sym, 4UL );
         checkCapacity( sym, 9UL );
         checkNonZeros( sym, 9UL );
         checkNonZeros( sym, 0UL, 2UL );
         checkNonZeros( sym, 1UL, 2UL );
         checkNonZeros( sym, 2UL, 2UL );
         checkNonZeros( sym, 3UL, 3UL );

         if( sym(0,0) != 2 || sym(0,3) != 3 ||
             sym(1,2) != 1 || sym(1,3) != 4 ||
             sym(2,1) != 1 || sym(2,3) != 5 ||
             sym(3,0) != 3 || sym(3,1) != 4 || sym(3,2) != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 2 0 0 3 )\n( 0 0 1 4 )\n( 0 1 0 5 )\n( 3 4 5 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Appending with row finalization
      {
         // Initialization check
         ST sym( 4UL, 9UL );
         sym.reserve( 0UL, 2UL );
         sym.reserve( 1UL, 4UL );
         sym.reserve( 2UL, 1UL );
         sym.reserve( 3UL, 2UL );

         // Appending one non-zero element
         sym.append( 0UL, 1UL, 1 );
         sym.finalize( 0UL );

         checkRows    ( sym, 4UL );
         checkColumns ( sym, 4UL );
         checkCapacity( sym, 9UL );
         checkNonZeros( sym, 2UL );
         checkNonZeros( sym, 0UL, 1UL );
         checkNonZeros( sym, 1UL, 1UL );
         checkNonZeros( sym, 2UL, 0UL );
         checkNonZeros( sym, 3UL, 0UL );

         if( sym(0,1) != 1 || sym(1,0) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 0 1 0 0 )\n( 1 0 0 0 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Appending two more non-zero elements
         sym.append( 1UL, 1UL, 2 );
         sym.append( 1UL, 2UL, 3 );
         sym.finalize( 1UL );

         checkRows    ( sym, 4UL );
         checkColumns ( sym, 4UL );
         checkCapacity( sym, 9UL );
         checkNonZeros( sym, 5UL );
         checkNonZeros( sym, 0UL, 1UL );
         checkNonZeros( sym, 1UL, 3UL );
         checkNonZeros( sym, 2UL, 1UL );
         checkNonZeros( sym, 3UL, 0UL );

         if( sym(0,1) != 1 ||
             sym(1,0) != 1 || sym(1,1) != 2 || sym(1,2) != 3 ||
             sym(2,1) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 0 1 0 0 )\n( 1 2 3 0 )\n( 0 3 0 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Appending two more non-zero elements
         sym.append( 3UL, 0UL, 4 );
         sym.append( 3UL, 1UL, 5 );
         sym.finalize( 3UL );

         checkRows    ( sym, 4UL );
         checkColumns ( sym, 4UL );
         checkCapacity( sym, 9UL );
         checkNonZeros( sym, 9UL );
         checkNonZeros( sym, 0UL, 2UL );
         checkNonZeros( sym, 1UL, 4UL );
         checkNonZeros( sym, 2UL, 1UL );
         checkNonZeros( sym, 3UL, 2UL );

         if( sym(0,1) != 1 || sym(0,3) != 4 ||
             sym(1,0) != 1 || sym(1,1) != 2 || sym(1,2) != 3 || sym(1,3) != 5 ||
             sym(2,1) != 3 ||
             sym(3,0) != 4 || sym(3,1) != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 0 1 0 4 )\n( 1 2 3 5 )\n( 0 3 0 0 )\n( 4 5 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major SymmetricMatrix::append()";

      // Appending with pre-allocation in each column
      {
         // Initialization check
         OST sym( 4UL, 9UL );
         sym.reserve( 0UL, 2UL );
         sym.reserve( 1UL, 2UL );
         sym.reserve( 2UL, 2UL );
         sym.reserve( 3UL, 3UL );

         checkRows    ( sym, 4UL );
         checkColumns ( sym, 4UL );
         checkCapacity( sym, 9UL );
         checkNonZeros( sym, 0UL );
         checkNonZeros( sym, 0UL, 0UL );
         checkNonZeros( sym, 1UL, 0UL );
         checkNonZeros( sym, 2UL, 0UL );
         checkNonZeros( sym, 3UL, 0UL );

         // Appending one non-zero element
         sym.append( 1UL, 2UL, 1 );

         checkRows    ( sym, 4UL );
         checkColumns ( sym, 4UL );
         checkCapacity( sym, 9UL );
         checkNonZeros( sym, 2UL );
         checkNonZeros( sym, 0UL, 0UL );
         checkNonZeros( sym, 1UL, 1UL );
         checkNonZeros( sym, 2UL, 1UL );
         checkNonZeros( sym, 3UL, 0UL );

         if( sym(1,2) != 1 || sym(2,1) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 0 0 0 0 )\n( 0 0 1 0 )\n( 0 1 0 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Appending two more non-zero elements
         sym.append( 0UL, 0UL, 2 );
         sym.append( 3UL, 0UL, 3 );

         checkRows    ( sym, 4UL );
         checkColumns ( sym, 4UL );
         checkCapacity( sym, 9UL );
         checkNonZeros( sym, 5UL );
         checkNonZeros( sym, 0UL, 2UL );
         checkNonZeros( sym, 1UL, 1UL );
         checkNonZeros( sym, 2UL, 1UL );
         checkNonZeros( sym, 3UL, 1UL );

         if( sym(0,0) != 2 || sym(0,3) != 3 ||
             sym(1,2) != 1 ||
             sym(2,1) != 1 ||
             sym(3,0) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 2 0 0 3 )\n( 0 0 1 0 )\n( 0 1 0 0 )\n( 3 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Appending two more non-zero elements
         sym.append( 1UL, 3UL, 4 );
         sym.append( 2UL, 3UL, 5 );

         checkRows    ( sym, 4UL );
         checkColumns ( sym, 4UL );
         checkCapacity( sym, 9UL );
         checkNonZeros( sym, 9UL );
         checkNonZeros( sym, 0UL, 2UL );
         checkNonZeros( sym, 1UL, 2UL );
         checkNonZeros( sym, 2UL, 2UL );
         checkNonZeros( sym, 3UL, 3UL );

         if( sym(0,0) != 2 || sym(0,3) != 3 ||
             sym(1,2) != 1 || sym(1,3) != 4 ||
             sym(2,1) != 1 || sym(2,3) != 5 ||
             sym(3,0) != 3 || sym(3,1) != 4 || sym(3,2) != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 2 0 0 3 )\n( 0 0 1 4 )\n( 0 1 0 5 )\n( 3 4 5 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Appending with column finalization
      {
         // Initialization check
         OST sym( 4UL, 9UL );
         sym.reserve( 0UL, 2UL );
         sym.reserve( 1UL, 4UL );
         sym.reserve( 2UL, 1UL );
         sym.reserve( 3UL, 2UL );

         // Appending one non-zero element
         sym.append( 1UL, 0UL, 1 );
         sym.finalize( 0UL );

         checkRows    ( sym, 4UL );
         checkColumns ( sym, 4UL );
         checkCapacity( sym, 9UL );
         checkNonZeros( sym, 2UL );
         checkNonZeros( sym, 0UL, 1UL );
         checkNonZeros( sym, 1UL, 1UL );
         checkNonZeros( sym, 2UL, 0UL );
         checkNonZeros( sym, 3UL, 0UL );

         if( sym(0,1) != 1 || sym(1,0) != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 0 1 0 0 )\n( 1 0 0 0 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Appending two more non-zero elements
         sym.append( 1UL, 1UL, 2 );
         sym.append( 2UL, 1UL, 3 );
         sym.finalize( 1UL );

         checkRows    ( sym, 4UL );
         checkColumns ( sym, 4UL );
         checkCapacity( sym, 9UL );
         checkNonZeros( sym, 5UL );
         checkNonZeros( sym, 0UL, 1UL );
         checkNonZeros( sym, 1UL, 3UL );
         checkNonZeros( sym, 2UL, 1UL );
         checkNonZeros( sym, 3UL, 0UL );

         if( sym(0,1) != 1 ||
             sym(1,0) != 1 || sym(1,1) != 2 || sym(1,2) != 3 ||
             sym(2,1) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 0 1 0 0 )\n( 1 2 3 0 )\n( 0 3 0 0 )\n( 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Appending two more non-zero elements
         sym.append( 0UL, 3UL, 4 );
         sym.append( 1UL, 3UL, 5 );
         sym.finalize( 3UL );

         checkRows    ( sym, 4UL );
         checkColumns ( sym, 4UL );
         checkCapacity( sym, 9UL );
         checkNonZeros( sym, 9UL );
         checkNonZeros( sym, 0UL, 2UL );
         checkNonZeros( sym, 1UL, 4UL );
         checkNonZeros( sym, 2UL, 1UL );
         checkNonZeros( sym, 3UL, 2UL );

         if( sym(0,1) != 1 || sym(0,3) != 4 ||
             sym(1,0) != 1 || sym(1,1) != 2 || sym(1,2) != 3 || sym(1,3) != 5 ||
             sym(2,1) != 3 ||
             sym(3,0) != 4 || sym(3,1) != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 0 1 0 4 )\n( 1 2 3 5 )\n( 0 3 0 0 )\n( 4 5 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c erase() member function of the SymmetricMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c erase() member function of the SymmetricMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseScalarTest::testErase()
{
   //=====================================================================================
   // Row-major index-based erase function
   //=====================================================================================

   {
      test_ = "Row-major SymmetricMatrix::erase( size_t, size_t )";

      // Initialization check
      ST sym( 4UL );
      sym(0,0) = 1;
      sym(0,2) = 2;
      sym(0,3) = 3;
      sym(1,1) = 4;
      sym(1,2) = 5;
      sym(2,2) = 6;
      sym(2,3) = 7;

      checkRows    ( sym,  4UL );
      checkColumns ( sym,  4UL );
      checkCapacity( sym, 11UL );
      checkNonZeros( sym, 11UL );
      checkNonZeros( sym, 0UL, 3UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 4UL );
      checkNonZeros( sym, 3UL, 2UL );

      if( sym(0,0) != 1 || sym(0,2) != 2 || sym(0,3) != 3 ||
          sym(1,1) != 4 || sym(1,2) != 5 ||
          sym(2,0) != 2 || sym(2,1) != 5 || sym(2,2) != 6 || sym(2,3) != 7 ||
          sym(3,0) != 3 || sym(3,2) != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 1 0 2 3 )\n( 0 4 5 0 )\n( 2 5 6 7 )\n( 3 0 7 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element at (0,0)
      sym.erase( 0UL, 0UL );

      checkRows    ( sym,  4UL );
      checkColumns ( sym,  4UL );
      checkCapacity( sym, 11UL );
      checkNonZeros( sym, 10UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 4UL );
      checkNonZeros( sym, 3UL, 2UL );

      if( sym(0,2) != 2 || sym(0,3) != 3 ||
          sym(1,1) != 4 || sym(1,2) != 5 ||
          sym(2,0) != 2 || sym(2,1) != 5 || sym(2,2) != 6 || sym(2,3) != 7 ||
          sym(3,0) != 3 || sym(3,2) != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 0 0 2 3 )\n( 0 4 5 0 )\n( 2 5 6 7 )\n( 3 0 7 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element at (1,2)
      sym.erase( 1UL, 2UL );

      checkRows    ( sym,  4UL );
      checkColumns ( sym,  4UL );
      checkCapacity( sym, 11UL );
      checkNonZeros( sym,  8UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 1UL );
      checkNonZeros( sym, 2UL, 3UL );
      checkNonZeros( sym, 3UL, 2UL );

      if( sym(0,2) != 2 || sym(0,3) != 3 ||
          sym(1,1) != 4 ||
          sym(2,0) != 2 || sym(2,2) != 6 || sym(2,3) != 7 ||
          sym(3,0) != 3 || sym(3,2) != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 0 0 2 3 )\n( 0 4 0 0 )\n( 2 0 6 7 )\n( 3 0 7 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element at (0,2)
      sym.erase( 0UL, 2UL );

      checkRows    ( sym,  4UL );
      checkColumns ( sym,  4UL );
      checkCapacity( sym, 11UL );
      checkNonZeros( sym,  6UL );
      checkNonZeros( sym, 0UL, 1UL );
      checkNonZeros( sym, 1UL, 1UL );
      checkNonZeros( sym, 2UL, 2UL );
      checkNonZeros( sym, 3UL, 2UL );

      if( sym(0,3) != 3 ||
          sym(1,1) != 4 ||
          sym(2,2) != 6 || sym(2,3) != 7 ||
          sym(3,0) != 3 || sym(3,2) != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 0 0 0 3 )\n( 0 4 0 0 )\n( 0 0 6 7 )\n( 3 0 7 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to erase a zero element
      sym.erase( 0UL, 1UL );

      checkRows    ( sym,  4UL );
      checkColumns ( sym,  4UL );
      checkCapacity( sym, 11UL );
      checkNonZeros( sym,  6UL );
      checkNonZeros( sym, 0UL, 1UL );
      checkNonZeros( sym, 1UL, 1UL );
      checkNonZeros( sym, 2UL, 2UL );
      checkNonZeros( sym, 3UL, 2UL );

      if( sym(0,3) != 3 ||
          sym(1,1) != 4 ||
          sym(2,2) != 6 || sym(2,3) != 7 ||
          sym(3,0) != 3 || sym(3,2) != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 0 0 0 3 )\n( 0 4 0 0 )\n( 0 0 6 7 )\n( 3 0 7 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major iterator-based erase function
   //=====================================================================================

   {
      test_ = "Row-major SymmetricMatrix::erase( size_t, Iterator )";

      using Iterator = ST::Iterator;

      // Initialization check
      ST sym( 4UL );
      sym(0,0) = 1;
      sym(0,2) = 2;
      sym(0,3) = 3;
      sym(1,1) = 4;
      sym(1,2) = 5;
      sym(2,2) = 6;
      sym(2,3) = 7;

      checkRows    ( sym,  4UL );
      checkColumns ( sym,  4UL );
      checkCapacity( sym, 11UL );
      checkNonZeros( sym, 11UL );
      checkNonZeros( sym, 0UL, 3UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 4UL );
      checkNonZeros( sym, 3UL, 2UL );

      if( sym(0,0) != 1 || sym(0,2) != 2 || sym(0,3) != 3 ||
          sym(1,1) != 4 || sym(1,2) != 5 ||
          sym(2,0) != 2 || sym(2,1) != 5 || sym(2,2) != 6 || sym(2,3) != 7 ||
          sym(3,0) != 3 || sym(3,2) != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 1 0 2 3 )\n( 0 4 5 0 )\n( 2 5 6 7 )\n( 3 0 7 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element at (0,0)
      {
         Iterator pos = sym.erase( 0UL, sym.find( 0UL, 0UL ) );

         checkRows    ( sym,  4UL );
         checkColumns ( sym,  4UL );
         checkCapacity( sym, 11UL );
         checkNonZeros( sym, 10UL );
         checkNonZeros( sym, 0UL, 2UL );
         checkNonZeros( sym, 1UL, 2UL );
         checkNonZeros( sym, 2UL, 4UL );
         checkNonZeros( sym, 3UL, 2UL );

         if( sym(0,2) != 2 || sym(0,3) != 3 ||
             sym(1,1) != 4 || sym(1,2) != 5 ||
             sym(2,0) != 2 || sym(2,1) != 5 || sym(2,2) != 6 || sym(2,3) != 7 ||
             sym(3,0) != 3 || sym(3,2) != 7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 0 0 2 3 )\n( 0 4 5 0 )\n( 2 5 6 7 )\n( 3 0 7 0 )\n";
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
         Iterator pos = sym.erase( 1UL, sym.find( 1UL, 2UL ) );

         checkRows    ( sym,  4UL );
         checkColumns ( sym,  4UL );
         checkCapacity( sym, 11UL );
         checkNonZeros( sym,  8UL );
         checkNonZeros( sym, 0UL, 2UL );
         checkNonZeros( sym, 1UL, 1UL );
         checkNonZeros( sym, 2UL, 3UL );
         checkNonZeros( sym, 3UL, 2UL );

         if( sym(0,2) != 2 || sym(0,3) != 3 ||
             sym(1,1) != 4 ||
             sym(2,0) != 2 || sym(2,2) != 6 || sym(2,3) != 7 ||
             sym(3,0) != 3 || sym(3,2) != 7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 0 0 2 3 )\n( 0 4 0 0 )\n( 2 0 6 7 )\n( 3 0 7 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos != sym.end( 1UL ) ) {
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
         Iterator pos = sym.erase( 0UL, sym.find( 0UL, 2UL ) );

         checkRows    ( sym,  4UL );
         checkColumns ( sym,  4UL );
         checkCapacity( sym, 11UL );
         checkNonZeros( sym,  6UL );
         checkNonZeros( sym, 0UL, 1UL );
         checkNonZeros( sym, 1UL, 1UL );
         checkNonZeros( sym, 2UL, 2UL );
         checkNonZeros( sym, 3UL, 2UL );

         if( sym(0,3) != 3 ||
             sym(1,1) != 4 ||
             sym(2,2) != 6 || sym(2,3) != 7 ||
             sym(3,0) != 3 || sym(3,2) != 7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 0 0 0 3 )\n( 0 4 0 0 )\n( 0 0 6 7 )\n( 3 0 7 0 )\n";
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
         Iterator pos = sym.erase( 0UL, sym.find( 0UL, 1UL ) );

         checkRows    ( sym,  4UL );
         checkColumns ( sym,  4UL );
         checkCapacity( sym, 11UL );
         checkNonZeros( sym,  6UL );
         checkNonZeros( sym, 0UL, 1UL );
         checkNonZeros( sym, 1UL, 1UL );
         checkNonZeros( sym, 2UL, 2UL );
         checkNonZeros( sym, 3UL, 2UL );

         if( sym(0,3) != 3 ||
             sym(1,1) != 4 ||
             sym(2,2) != 6 || sym(2,3) != 7 ||
             sym(3,0) != 3 || sym(3,2) != 7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 0 0 0 3 )\n( 0 4 0 0 )\n( 0 0 6 7 )\n( 3 0 7 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos != sym.end( 0UL ) ) {
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
      test_ = "Row-major SymmetricMatrix::erase( size_t, Iterator, Iterator )";

      using Iterator = ST::Iterator;

      // Initialization check
      ST sym( 4UL );
      sym(0,0) = 1;
      sym(0,2) = 2;
      sym(0,3) = 3;
      sym(1,1) = 4;
      sym(1,2) = 5;
      sym(2,2) = 6;
      sym(2,3) = 7;

      checkRows    ( sym,  4UL );
      checkColumns ( sym,  4UL );
      checkCapacity( sym, 11UL );
      checkNonZeros( sym, 11UL );
      checkNonZeros( sym, 0UL, 3UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 4UL );
      checkNonZeros( sym, 3UL, 2UL );

      if( sym(0,0) != 1 || sym(0,2) != 2 || sym(0,3) != 3 ||
          sym(1,1) != 4 || sym(1,2) != 5 ||
          sym(2,0) != 2 || sym(2,1) != 5 || sym(2,2) != 6 || sym(2,3) != 7 ||
          sym(3,0) != 3 || sym(3,2) != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 1 0 2 3 )\n( 0 4 5 0 )\n( 2 5 6 7 )\n( 3 0 7 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element from (0,0) to (0,2)
      {
         Iterator pos = sym.erase( 0UL, sym.find( 0UL, 0UL ), sym.find( 0UL, 2UL ) );

         checkRows    ( sym,  4UL );
         checkColumns ( sym,  4UL );
         checkCapacity( sym, 11UL );
         checkNonZeros( sym, 10UL );
         checkNonZeros( sym, 0UL, 2UL );
         checkNonZeros( sym, 1UL, 2UL );
         checkNonZeros( sym, 2UL, 4UL );
         checkNonZeros( sym, 3UL, 2UL );

         if( sym(0,2) != 2 || sym(0,3) != 3 ||
             sym(1,1) != 4 || sym(1,2) != 5 ||
             sym(2,0) != 2 || sym(2,1) != 5 || sym(2,2) != 6 || sym(2,3) != 7 ||
             sym(3,0) != 3 || sym(3,2) != 7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a single-element range failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 0 0 2 3 )\n( 0 4 5 0 )\n( 2 5 6 7 )\n( 3 0 7 0 )\n";
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

      // Erasing the element from (2,1) to (2,3)
      {
         Iterator pos = sym.erase( 2UL, sym.find( 2UL, 1UL ), sym.find( 2UL, 3UL ) );

         checkRows    ( sym,  4UL );
         checkColumns ( sym,  4UL );
         checkCapacity( sym, 11UL );
         checkNonZeros( sym,  7UL );
         checkNonZeros( sym, 0UL, 2UL );
         checkNonZeros( sym, 1UL, 1UL );
         checkNonZeros( sym, 2UL, 2UL );
         checkNonZeros( sym, 3UL, 2UL );

         if( sym(0,2) != 2 || sym(0,3) != 3 ||
             sym(1,1) != 4 ||
             sym(2,0) != 2 || sym(2,3) != 7 ||
             sym(3,0) != 3 || sym(3,2) != 7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a single-element range failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 0 0 2 3 )\n( 0 4 0 0 )\n( 2 0 0 7 )\n( 3 0 7 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos->value() != 7 || pos->index() != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 7\n"
                << "   Expected index: 3\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the element from (3,2) to the row end
      {
         Iterator pos = sym.erase( 3UL, sym.find( 3UL, 2UL ), sym.end( 3UL ) );

         checkRows    ( sym,  4UL );
         checkColumns ( sym,  4UL );
         checkCapacity( sym, 11UL );
         checkNonZeros( sym,  5UL );
         checkNonZeros( sym, 0UL, 2UL );
         checkNonZeros( sym, 1UL, 1UL );
         checkNonZeros( sym, 2UL, 1UL );
         checkNonZeros( sym, 3UL, 1UL );

         if( sym(0,2) != 2 || sym(0,3) != 3 ||
             sym(1,1) != 4 ||
             sym(2,0) != 2 ||
             sym(3,0) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a single-element range failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 0 0 2 3 )\n( 0 4 0 0 )\n( 2 0 0 0 )\n( 3 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos != sym.end( 3UL ) ) {
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
         Iterator pos = sym.erase( 2UL, sym.find( 2UL, 0UL ), sym.find( 2UL, 0UL ) );

         checkRows    ( sym,  4UL );
         checkColumns ( sym,  4UL );
         checkCapacity( sym, 11UL );
         checkNonZeros( sym,  5UL );
         checkNonZeros( sym, 0UL, 2UL );
         checkNonZeros( sym, 1UL, 1UL );
         checkNonZeros( sym, 2UL, 1UL );
         checkNonZeros( sym, 3UL, 1UL );

         if( sym(0,2) != 2 || sym(0,3) != 3 ||
             sym(1,1) != 4 ||
             sym(2,0) != 2 ||
             sym(3,0) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing an empty range failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 0 0 2 3 )\n( 0 4 0 0 )\n( 2 0 0 0 )\n( 3 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos != sym.find( 2UL, 0UL ) ) {
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
      test_ = "Row-major SymmetricMatrix::erase( Predicate )";

      // Initialization check
      ST sym( 4UL );
      sym(0,0) = 1;
      sym(0,2) = 2;
      sym(0,3) = 3;
      sym(1,1) = 4;
      sym(1,2) = 5;
      sym(2,2) = 6;
      sym(2,3) = 7;

      checkRows    ( sym,  4UL );
      checkColumns ( sym,  4UL );
      checkCapacity( sym, 11UL );
      checkNonZeros( sym, 11UL );
      checkNonZeros( sym, 0UL, 3UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 4UL );
      checkNonZeros( sym, 3UL, 2UL );

      if( sym(0,0) != 1 || sym(0,2) != 2 || sym(0,3) != 3 ||
          sym(1,1) != 4 || sym(1,2) != 5 ||
          sym(2,0) != 2 || sym(2,1) != 5 || sym(2,2) != 6 || sym(2,3) != 7 ||
          sym(3,0) != 3 || sym(3,2) != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 1 0 2 3 )\n( 0 4 5 0 )\n( 2 5 6 7 )\n( 3 0 7 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing a selection of elements
      sym.erase( []( int value ){ return value == 1 || value == 5 || value == 6; } );

      checkRows    ( sym,  4UL );
      checkColumns ( sym,  4UL );
      checkCapacity( sym, 11UL );
      checkNonZeros( sym,  7UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 1UL );
      checkNonZeros( sym, 2UL, 2UL );
      checkNonZeros( sym, 3UL, 2UL );

      if( sym(0,2) != 2 || sym(0,3) != 3 ||
          sym(1,1) != 4 ||
          sym(2,0) != 2 || sym(2,3) != 7 ||
          sym(3,0) != 3 || sym(3,2) != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a selection of elements failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 0 0 2 3 )\n( 0 4 0 0 )\n( 2 0 0 7 )\n( 3 0 7 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to erase all elements with value 1
      sym.erase( []( int value ){ return value == 1; } );

      checkRows    ( sym,  4UL );
      checkColumns ( sym,  4UL );
      checkCapacity( sym, 11UL );
      checkNonZeros( sym,  7UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 1UL );
      checkNonZeros( sym, 2UL, 2UL );
      checkNonZeros( sym, 3UL, 2UL );

      if( sym(0,2) != 2 || sym(0,3) != 3 ||
          sym(1,1) != 4 ||
          sym(2,0) != 2 || sym(2,3) != 7 ||
          sym(3,0) != 3 || sym(3,2) != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing all elements with value 1 failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 0 0 2 3 )\n( 0 4 0 0 )\n( 2 0 0 7 )\n( 3 0 7 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major iterator-range-based erase function with predicate
   //=====================================================================================

   {
      test_ = "Row-major SymmetricMatrix::erase( size_t, Iterator, Iterator, Predicate )";

      // Initialization check
      ST sym( 4UL );
      sym(0,0) = 1;
      sym(0,2) = 2;
      sym(0,3) = 3;
      sym(1,1) = 4;
      sym(1,2) = 5;
      sym(2,2) = 6;
      sym(2,3) = 7;

      checkRows    ( sym,  4UL );
      checkColumns ( sym,  4UL );
      checkCapacity( sym, 11UL );
      checkNonZeros( sym, 11UL );
      checkNonZeros( sym, 0UL, 3UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 4UL );
      checkNonZeros( sym, 3UL, 2UL );

      if( sym(0,0) != 1 || sym(0,2) != 2 || sym(0,3) != 3 ||
          sym(1,1) != 4 || sym(1,2) != 5 ||
          sym(2,0) != 2 || sym(2,1) != 5 || sym(2,2) != 6 || sym(2,3) != 7 ||
          sym(3,0) != 3 || sym(3,2) != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 1 0 2 3 )\n( 0 4 5 0 )\n( 2 5 6 7 )\n( 3 0 7 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing a selection of elements
      sym.erase( 2UL, sym.begin( 2UL ), sym.find( 2UL, 3UL ),
                 []( int value ){ return value == 2 || value == 6; } );

      checkRows    ( sym,  4UL );
      checkColumns ( sym,  4UL );
      checkCapacity( sym, 11UL );
      checkNonZeros( sym,  8UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 2UL );
      checkNonZeros( sym, 3UL, 2UL );

      if( sym(0,0) != 1 || sym(0,3) != 3 ||
          sym(1,1) != 4 || sym(1,2) != 5 ||
          sym(2,1) != 5 || sym(2,3) != 7 ||
          sym(3,0) != 3 || sym(3,2) != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a selection of elements failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 1 0 0 3 )\n( 0 4 5 0 )\n( 0 5 0 7 )\n( 3 0 7 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to erase from an empty range
      sym.erase( 1UL, sym.begin( 1UL ), sym.begin( 1UL ), []( int ){ return true; } );

      checkRows    ( sym,  4UL );
      checkColumns ( sym,  4UL );
      checkCapacity( sym, 11UL );
      checkNonZeros( sym,  8UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 2UL );
      checkNonZeros( sym, 3UL, 2UL );

      if( sym(0,0) != 1 || sym(0,3) != 3 ||
          sym(1,1) != 4 || sym(1,2) != 5 ||
          sym(2,1) != 5 || sym(2,3) != 7 ||
          sym(3,0) != 3 || sym(3,2) != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing from an empty range failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 1 0 0 3 )\n( 0 4 5 0 )\n( 0 5 0 7 )\n( 3 0 7 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major index-based erase function
   //=====================================================================================

   {
      test_ = "Column-major SymmetricMatrix::erase( size_t, size_t )";

      // Initialization check
      OST sym( 4UL );
      sym(0,0) = 1;
      sym(0,2) = 2;
      sym(0,3) = 3;
      sym(1,1) = 4;
      sym(1,2) = 5;
      sym(2,2) = 6;
      sym(2,3) = 7;

      checkRows    ( sym,  4UL );
      checkColumns ( sym,  4UL );
      checkCapacity( sym, 11UL );
      checkNonZeros( sym, 11UL );
      checkNonZeros( sym, 0UL, 3UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 4UL );
      checkNonZeros( sym, 3UL, 2UL );

      if( sym(0,0) != 1 || sym(0,2) != 2 || sym(0,3) != 3 ||
          sym(1,1) != 4 || sym(1,2) != 5 ||
          sym(2,0) != 2 || sym(2,1) != 5 || sym(2,2) != 6 || sym(2,3) != 7 ||
          sym(3,0) != 3 || sym(3,2) != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 1 0 2 3 )\n( 0 4 5 0 )\n( 2 5 6 7 )\n( 3 0 7 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element at (0,0)
      sym.erase( 0UL, 0UL );

      checkRows    ( sym,  4UL );
      checkColumns ( sym,  4UL );
      checkCapacity( sym, 11UL );
      checkNonZeros( sym, 10UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 4UL );
      checkNonZeros( sym, 3UL, 2UL );

      if( sym(0,2) != 2 || sym(0,3) != 3 ||
          sym(1,1) != 4 || sym(1,2) != 5 ||
          sym(2,0) != 2 || sym(2,1) != 5 || sym(2,2) != 6 || sym(2,3) != 7 ||
          sym(3,0) != 3 || sym(3,2) != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 0 0 2 3 )\n( 0 4 5 0 )\n( 2 5 6 7 )\n( 3 0 7 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element at (2,1)
      sym.erase( 2UL, 1UL );

      checkRows    ( sym,  4UL );
      checkColumns ( sym,  4UL );
      checkCapacity( sym, 11UL );
      checkNonZeros( sym,  8UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 1UL );
      checkNonZeros( sym, 2UL, 3UL );
      checkNonZeros( sym, 3UL, 2UL );

      if( sym(0,2) != 2 || sym(0,3) != 3 ||
          sym(1,1) != 4 ||
          sym(2,0) != 2 || sym(2,2) != 6 || sym(2,3) != 7 ||
          sym(3,0) != 3 || sym(3,2) != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 0 0 2 3 )\n( 0 4 0 0 )\n( 2 0 6 7 )\n( 3 0 7 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element at (2,0)
      sym.erase( 2UL, 0UL );

      checkRows    ( sym,  4UL );
      checkColumns ( sym,  4UL );
      checkCapacity( sym, 11UL );
      checkNonZeros( sym,  6UL );
      checkNonZeros( sym, 0UL, 1UL );
      checkNonZeros( sym, 1UL, 1UL );
      checkNonZeros( sym, 2UL, 2UL );
      checkNonZeros( sym, 3UL, 2UL );

      if( sym(0,3) != 3 ||
          sym(1,1) != 4 ||
          sym(2,2) != 6 || sym(2,3) != 7 ||
          sym(3,0) != 3 || sym(3,2) != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 0 0 0 3 )\n( 0 4 0 0 )\n( 0 0 6 7 )\n( 3 0 7 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to erase a zero element
      sym.erase( 1UL, 0UL );

      checkRows    ( sym,  4UL );
      checkColumns ( sym,  4UL );
      checkCapacity( sym, 11UL );
      checkNonZeros( sym,  6UL );
      checkNonZeros( sym, 0UL, 1UL );
      checkNonZeros( sym, 1UL, 1UL );
      checkNonZeros( sym, 2UL, 2UL );
      checkNonZeros( sym, 3UL, 2UL );

      if( sym(0,3) != 3 ||
          sym(1,1) != 4 ||
          sym(2,2) != 6 || sym(2,3) != 7 ||
          sym(3,0) != 3 || sym(3,2) != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 0 0 0 3 )\n( 0 4 0 0 )\n( 0 0 6 7 )\n( 3 0 7 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major iterator-based erase function
   //=====================================================================================

   {
      test_ = "Column-major SymmetricMatrix::erase( size_t, Iterator )";

      using Iterator = OST::Iterator;

      // Initialization check
      OST sym( 4UL );
      sym(0,0) = 1;
      sym(0,2) = 2;
      sym(0,3) = 3;
      sym(1,1) = 4;
      sym(1,2) = 5;
      sym(2,2) = 6;
      sym(2,3) = 7;

      checkRows    ( sym,  4UL );
      checkColumns ( sym,  4UL );
      checkCapacity( sym, 11UL );
      checkNonZeros( sym, 11UL );
      checkNonZeros( sym, 0UL, 3UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 4UL );
      checkNonZeros( sym, 3UL, 2UL );

      if( sym(0,0) != 1 || sym(0,2) != 2 || sym(0,3) != 3 ||
          sym(1,1) != 4 || sym(1,2) != 5 ||
          sym(2,0) != 2 || sym(2,1) != 5 || sym(2,2) != 6 || sym(2,3) != 7 ||
          sym(3,0) != 3 || sym(3,2) != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 1 0 2 3 )\n( 0 4 5 0 )\n( 2 5 6 7 )\n( 3 0 7 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element at (0,0)
      {
         Iterator pos = sym.erase( 0UL, sym.find( 0UL, 0UL ) );

         checkRows    ( sym,  4UL );
         checkColumns ( sym,  4UL );
         checkCapacity( sym, 11UL );
         checkNonZeros( sym, 10UL );
         checkNonZeros( sym, 0UL, 2UL );
         checkNonZeros( sym, 1UL, 2UL );
         checkNonZeros( sym, 2UL, 4UL );
         checkNonZeros( sym, 3UL, 2UL );

         if( sym(0,2) != 2 || sym(0,3) != 3 ||
             sym(1,1) != 4 || sym(1,2) != 5 ||
             sym(2,0) != 2 || sym(2,1) != 5 || sym(2,2) != 6 || sym(2,3) != 7 ||
             sym(3,0) != 3 || sym(3,2) != 7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 0 0 2 3 )\n( 0 4 5 0 )\n( 2 5 6 7 )\n( 3 0 7 0 )\n";
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
         Iterator pos = sym.erase( 1UL, sym.find( 2UL, 1UL ) );

         checkRows    ( sym,  4UL );
         checkColumns ( sym,  4UL );
         checkCapacity( sym, 11UL );
         checkNonZeros( sym,  8UL );
         checkNonZeros( sym, 0UL, 2UL );
         checkNonZeros( sym, 1UL, 1UL );
         checkNonZeros( sym, 2UL, 3UL );
         checkNonZeros( sym, 3UL, 2UL );

         if( sym(0,2) != 2 || sym(0,3) != 3 ||
             sym(1,1) != 4 ||
             sym(2,0) != 2 || sym(2,2) != 6 || sym(2,3) != 7 ||
             sym(3,0) != 3 || sym(3,2) != 7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 0 0 2 3 )\n( 0 4 0 0 )\n( 2 0 6 7 )\n( 3 0 7 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos != sym.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the element at (2,0)
      {
         Iterator pos = sym.erase( 0UL, sym.find( 2UL, 0UL ) );

         checkRows    ( sym,  4UL );
         checkColumns ( sym,  4UL );
         checkCapacity( sym, 11UL );
         checkNonZeros( sym,  6UL );
         checkNonZeros( sym, 0UL, 1UL );
         checkNonZeros( sym, 1UL, 1UL );
         checkNonZeros( sym, 2UL, 2UL );
         checkNonZeros( sym, 3UL, 2UL );

         if( sym(0,3) != 3 ||
             sym(1,1) != 4 ||
             sym(2,2) != 6 || sym(2,3) != 7 ||
             sym(3,0) != 3 || sym(3,2) != 7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 0 0 0 3 )\n( 0 4 0 0 )\n( 0 0 6 7 )\n( 3 0 7 0 )\n";
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
         Iterator pos = sym.erase( 0UL, sym.find( 1UL, 0UL ) );

         checkRows    ( sym,  4UL );
         checkColumns ( sym,  4UL );
         checkCapacity( sym, 11UL );
         checkNonZeros( sym,  6UL );
         checkNonZeros( sym, 0UL, 1UL );
         checkNonZeros( sym, 1UL, 1UL );
         checkNonZeros( sym, 2UL, 2UL );
         checkNonZeros( sym, 3UL, 2UL );

         if( sym(0,3) != 3 ||
             sym(1,1) != 4 ||
             sym(2,2) != 6 || sym(2,3) != 7 ||
             sym(3,0) != 3 || sym(3,2) != 7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 0 0 0 3 )\n( 0 4 0 0 )\n( 0 0 6 7 )\n( 3 0 7 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos != sym.end( 0UL ) ) {
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
      test_ = "Column-major SymmetricMatrix::erase( size_t, Iterator, Iterator )";

      using Iterator = OST::Iterator;

      // Initialization check
      OST sym( 4UL );
      sym(0,0) = 1;
      sym(0,2) = 2;
      sym(0,3) = 3;
      sym(1,1) = 4;
      sym(1,2) = 5;
      sym(2,2) = 6;
      sym(2,3) = 7;

      checkRows    ( sym,  4UL );
      checkColumns ( sym,  4UL );
      checkCapacity( sym, 11UL );
      checkNonZeros( sym, 11UL );
      checkNonZeros( sym, 0UL, 3UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 4UL );
      checkNonZeros( sym, 3UL, 2UL );

      if( sym(0,0) != 1 || sym(0,2) != 2 || sym(0,3) != 3 ||
          sym(1,1) != 4 || sym(1,2) != 5 ||
          sym(2,0) != 2 || sym(2,1) != 5 || sym(2,2) != 6 || sym(2,3) != 7 ||
          sym(3,0) != 3 || sym(3,2) != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 1 0 2 3 )\n( 0 4 5 0 )\n( 2 5 6 7 )\n( 3 0 7 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element from (0,0) to (2,0)
      {
         Iterator pos = sym.erase( 0UL, sym.find( 0UL, 0UL ), sym.find( 2UL, 0UL ) );

         checkRows    ( sym,  4UL );
         checkColumns ( sym,  4UL );
         checkCapacity( sym, 11UL );
         checkNonZeros( sym, 10UL );
         checkNonZeros( sym, 0UL, 2UL );
         checkNonZeros( sym, 1UL, 2UL );
         checkNonZeros( sym, 2UL, 4UL );
         checkNonZeros( sym, 3UL, 2UL );

         if( sym(0,2) != 2 || sym(0,3) != 3 ||
             sym(1,1) != 4 || sym(1,2) != 5 ||
             sym(2,0) != 2 || sym(2,1) != 5 || sym(2,2) != 6 || sym(2,3) != 7 ||
             sym(3,0) != 3 || sym(3,2) != 7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a single-element range failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 0 0 2 3 )\n( 0 4 5 0 )\n( 2 5 6 7 )\n( 3 0 7 0 )\n";
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

      // Erasing the element from (1,2) to (3,2)
      {
         Iterator pos = sym.erase( 2UL, sym.find( 1UL, 2UL ), sym.find( 3UL, 2UL ) );

         checkRows    ( sym,  4UL );
         checkColumns ( sym,  4UL );
         checkCapacity( sym, 11UL );
         checkNonZeros( sym,  7UL );
         checkNonZeros( sym, 0UL, 2UL );
         checkNonZeros( sym, 1UL, 1UL );
         checkNonZeros( sym, 2UL, 2UL );
         checkNonZeros( sym, 3UL, 2UL );

         if( sym(0,2) != 2 || sym(0,3) != 3 ||
             sym(1,1) != 4 ||
             sym(2,0) != 2 || sym(2,3) != 7 ||
             sym(3,0) != 3 || sym(3,2) != 7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a single-element range failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 0 0 2 3 )\n( 0 4 0 0 )\n( 2 0 0 7 )\n( 3 0 7 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos->value() != 7 || pos->index() != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: 7\n"
                << "   Expected index: 3\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the element from (2,3) to the column end
      {
         Iterator pos = sym.erase( 3UL, sym.find( 2UL, 3UL ), sym.end( 3UL ) );

         checkRows    ( sym,  4UL );
         checkColumns ( sym,  4UL );
         checkCapacity( sym, 11UL );
         checkNonZeros( sym,  5UL );
         checkNonZeros( sym, 0UL, 2UL );
         checkNonZeros( sym, 1UL, 1UL );
         checkNonZeros( sym, 2UL, 1UL );
         checkNonZeros( sym, 3UL, 1UL );

         if( sym(0,2) != 2 || sym(0,3) != 3 ||
             sym(1,1) != 4 ||
             sym(2,0) != 2 ||
             sym(3,0) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a single-element range failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 0 0 2 3 )\n( 0 4 0 0 )\n( 2 0 0 0 )\n( 3 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos != sym.end( 3UL ) ) {
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
         Iterator pos = sym.erase( 2UL, sym.find( 0UL, 2UL ), sym.find( 0UL, 2UL ) );

         checkRows    ( sym,  4UL );
         checkColumns ( sym,  4UL );
         checkCapacity( sym, 11UL );
         checkNonZeros( sym,  5UL );
         checkNonZeros( sym, 0UL, 2UL );
         checkNonZeros( sym, 1UL, 1UL );
         checkNonZeros( sym, 2UL, 1UL );
         checkNonZeros( sym, 3UL, 1UL );

         if( sym(0,2) != 2 || sym(0,3) != 3 ||
             sym(1,1) != 4 ||
             sym(2,0) != 2 ||
             sym(3,0) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing an empty range failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 0 0 2 3 )\n( 0 4 0 0 )\n( 2 0 0 0 )\n( 3 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos != sym.find( 0UL, 2UL ) ) {
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
   // Column-major erase function with predicate
   //=====================================================================================

   {
      test_ = "Column-major SymmetricMatrix::erase( Predicate )";

      // Initialization check
      OST sym( 4UL );
      sym(0,0) = 1;
      sym(0,2) = 2;
      sym(0,3) = 3;
      sym(1,1) = 4;
      sym(1,2) = 5;
      sym(2,2) = 6;
      sym(2,3) = 7;

      checkRows    ( sym,  4UL );
      checkColumns ( sym,  4UL );
      checkCapacity( sym, 11UL );
      checkNonZeros( sym, 11UL );
      checkNonZeros( sym, 0UL, 3UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 4UL );
      checkNonZeros( sym, 3UL, 2UL );

      if( sym(0,0) != 1 || sym(0,2) != 2 || sym(0,3) != 3 ||
          sym(1,1) != 4 || sym(1,2) != 5 ||
          sym(2,0) != 2 || sym(2,1) != 5 || sym(2,2) != 6 || sym(2,3) != 7 ||
          sym(3,0) != 3 || sym(3,2) != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 1 0 2 3 )\n( 0 4 5 0 )\n( 2 5 6 7 )\n( 3 0 7 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing a selection of elements
      sym.erase( []( int value ){ return value == 1 || value == 5 || value == 6; } );

      checkRows    ( sym,  4UL );
      checkColumns ( sym,  4UL );
      checkCapacity( sym, 11UL );
      checkNonZeros( sym,  7UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 1UL );
      checkNonZeros( sym, 2UL, 2UL );
      checkNonZeros( sym, 3UL, 2UL );

      if( sym(0,2) != 2 || sym(0,3) != 3 ||
          sym(1,1) != 4 ||
          sym(2,0) != 2 || sym(2,3) != 7 ||
          sym(3,0) != 3 || sym(3,2) != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a selection of elements failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 0 0 2 3 )\n( 0 4 0 0 )\n( 2 0 0 7 )\n( 3 0 7 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to erase all elements with value 1
      sym.erase( []( int value ){ return value == 1; } );

      checkRows    ( sym,  4UL );
      checkColumns ( sym,  4UL );
      checkCapacity( sym, 11UL );
      checkNonZeros( sym,  7UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 1UL );
      checkNonZeros( sym, 2UL, 2UL );
      checkNonZeros( sym, 3UL, 2UL );

      if( sym(0,2) != 2 || sym(0,3) != 3 ||
          sym(1,1) != 4 ||
          sym(2,0) != 2 || sym(2,3) != 7 ||
          sym(3,0) != 3 || sym(3,2) != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing all elements with value 1 failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 0 0 2 3 )\n( 0 4 0 0 )\n( 2 0 0 7 )\n( 3 0 7 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major iterator-range-based erase function with predicate
   //=====================================================================================

   {
      test_ = "Column-major SymmetricMatrix::erase( size_t, Iterator, Iterator, Predicate )";

      // Initialization check
      OST sym( 4UL );
      sym(0,0) = 1;
      sym(0,2) = 2;
      sym(0,3) = 3;
      sym(1,1) = 4;
      sym(1,2) = 5;
      sym(2,2) = 6;
      sym(2,3) = 7;

      checkRows    ( sym,  4UL );
      checkColumns ( sym,  4UL );
      checkCapacity( sym, 11UL );
      checkNonZeros( sym, 11UL );
      checkNonZeros( sym, 0UL, 3UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 4UL );
      checkNonZeros( sym, 3UL, 2UL );

      if( sym(0,0) != 1 || sym(0,2) != 2 || sym(0,3) != 3 ||
          sym(1,1) != 4 || sym(1,2) != 5 ||
          sym(2,0) != 2 || sym(2,1) != 5 || sym(2,2) != 6 || sym(2,3) != 7 ||
          sym(3,0) != 3 || sym(3,2) != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 1 0 2 3 )\n( 0 4 5 0 )\n( 2 5 6 7 )\n( 3 0 7 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing a selection of elements
      sym.erase( 2UL, sym.begin( 2UL ), sym.find( 3UL, 2UL ),
                 []( int value ){ return value == 2 || value == 6; } );

      checkRows    ( sym,  4UL );
      checkColumns ( sym,  4UL );
      checkCapacity( sym, 11UL );
      checkNonZeros( sym,  8UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 2UL );
      checkNonZeros( sym, 3UL, 2UL );

      if( sym(0,0) != 1 || sym(0,3) != 3 ||
          sym(1,1) != 4 || sym(1,2) != 5 ||
          sym(2,1) != 5 || sym(2,3) != 7 ||
          sym(3,0) != 3 || sym(3,2) != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a selection of elements failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 1 0 0 3 )\n( 0 4 5 0 )\n( 0 5 0 7 )\n( 3 0 7 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to erase from an empty range
      sym.erase( 1UL, sym.begin( 1UL ), sym.begin( 1UL ), []( int ){ return true; } );

      checkRows    ( sym,  4UL );
      checkColumns ( sym,  4UL );
      checkCapacity( sym, 11UL );
      checkNonZeros( sym,  8UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 2UL );
      checkNonZeros( sym, 3UL, 2UL );

      if( sym(0,0) != 1 || sym(0,3) != 3 ||
          sym(1,1) != 4 || sym(1,2) != 5 ||
          sym(2,1) != 5 || sym(2,3) != 7 ||
          sym(3,0) != 3 || sym(3,2) != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing from an empty range failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 1 0 0 3 )\n( 0 4 5 0 )\n( 0 5 0 7 )\n( 3 0 7 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c find() member function of the SymmetricMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c find() member function of the SymmetricMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseScalarTest::testFind()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major SymmetricMatrix::find()";

      using ConstIterator = ST::ConstIterator;

      // Initialization check
      ST sym( 8UL, 3UL );
      sym(1,2) = 1;
      sym(2,3) = 2;
      sym(6,5) = 3;

      checkRows    ( sym, 8UL );
      checkColumns ( sym, 8UL );
      checkCapacity( sym, 3UL );
      checkNonZeros( sym, 6UL );
      checkNonZeros( sym, 0UL, 0UL );
      checkNonZeros( sym, 1UL, 1UL );
      checkNonZeros( sym, 2UL, 2UL );
      checkNonZeros( sym, 3UL, 1UL );
      checkNonZeros( sym, 4UL, 0UL );
      checkNonZeros( sym, 5UL, 1UL );
      checkNonZeros( sym, 6UL, 1UL );
      checkNonZeros( sym, 7UL, 0UL );

      // Searching for the first element
      {
         ConstIterator pos( sym.find( 1UL, 2UL ) );

         if( pos == sym.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required position = (1,2)\n"
                << "   Current matrix:\n" << sym << "\n";
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
                << "   Current matrix:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for the second element
      {
         ConstIterator pos( sym.find( 2UL, 3UL ) );

         if( pos == sym.end( 2UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required position = (2,3)\n"
                << "   Current matrix:\n" << sym << "\n";
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
                << "   Current matrix:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for the third element
      {
         ConstIterator pos( sym.find( 6UL, 5UL ) );

         if( pos == sym.end( 6UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required position = (6,5)\n"
                << "   Current matrix:\n" << sym << "\n";
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
                << "   Current matrix:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for a non-existing non-zero element
      {
         ConstIterator pos( sym.find( 4UL, 0UL ) );

         if( pos != sym.end( 4UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Non-existing element could be found\n"
                << " Details:\n"
                << "   Required index = 4\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 0\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major SymmetricMatrix::find()";

      using ConstIterator = OST::ConstIterator;

      // Initialization check
      OST sym( 8UL, 3UL );
      sym(2,1) = 1;
      sym(3,2) = 2;
      sym(5,6) = 3;

      checkRows    ( sym, 8UL );
      checkColumns ( sym, 8UL );
      checkCapacity( sym, 3UL );
      checkNonZeros( sym, 6UL );
      checkNonZeros( sym, 0UL, 0UL );
      checkNonZeros( sym, 1UL, 1UL );
      checkNonZeros( sym, 2UL, 2UL );
      checkNonZeros( sym, 3UL, 1UL );
      checkNonZeros( sym, 4UL, 0UL );
      checkNonZeros( sym, 5UL, 1UL );
      checkNonZeros( sym, 6UL, 1UL );
      checkNonZeros( sym, 7UL, 0UL );

      // Searching for the first element
      {
         ConstIterator pos( sym.find( 2UL, 1UL ) );

         if( pos == sym.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required position = (2,1)\n"
                << "   Current matrix:\n" << sym << "\n";
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
                << "   Current matrix:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for the second element
      {
         ConstIterator pos( sym.find( 3UL, 2UL ) );

         if( pos == sym.end( 2UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required position = (3,2)\n"
                << "   Current matrix:\n" << sym << "\n";
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
                << "   Current matrix:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for the third element
      {
         ConstIterator pos( sym.find( 5UL, 6UL ) );

         if( pos == sym.end( 6UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required position = (5,6)\n"
                << "   Current matrix:\n" << sym << "\n";
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
                << "   Current matrix:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for a non-existing non-zero element
      {
         ConstIterator pos( sym.find( 0UL, 4UL ) );

         if( pos != sym.end( 4UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Non-existing element could be found\n"
                << " Details:\n"
                << "   Required index = 4\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 0\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c lowerBound() member function of the SymmetricMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c lowerBound() member function of the SymmetricMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseScalarTest::testLowerBound()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major SymmetricMatrix::lowerBound()";

      using ConstIterator = ST::ConstIterator;

      // Initialization check
      ST sym( 6UL, 3UL );
      sym(1,2) = 1;
      sym(1,4) = 2;

      checkRows    ( sym, 6UL );
      checkColumns ( sym, 6UL );
      checkCapacity( sym, 4UL );
      checkNonZeros( sym, 4UL );
      checkNonZeros( sym, 0UL, 0UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 1UL );
      checkNonZeros( sym, 3UL, 0UL );
      checkNonZeros( sym, 4UL, 1UL );
      checkNonZeros( sym, 5UL, 0UL );

      // Determining the lower bound for position (1,1)
      {
         ConstIterator pos( sym.lowerBound( 1UL, 1UL ) );

         if( pos == sym.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,1)\n"
                << "   Current matrix:\n" << sym << "\n";
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
                << "   Current matrix:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for position (1,2)
      {
         ConstIterator pos( sym.lowerBound( 1UL, 2UL ) );

         if( pos == sym.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,2)\n"
                << "   Current matrix:\n" << sym << "\n";
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
                << "   Current matrix:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for position (1,3)
      {
         ConstIterator pos( sym.lowerBound( 1UL, 3UL ) );

         if( pos == sym.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,3)\n"
                << "   Current matrix:\n" << sym << "\n";
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
                << "   Current matrix:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for position (1,4)
      {
         ConstIterator pos( sym.lowerBound( 1UL, 4UL ) );

         if( pos == sym.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,4)\n"
                << "   Current matrix:\n" << sym << "\n";
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
                << "   Current matrix:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for position (1,5)
      {
         ConstIterator pos( sym.lowerBound( 1UL, 5UL ) );

         if( pos != sym.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,5)\n"
                << "   Current matrix:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major SymmetricMatrix::lowerBound()";

      using ConstIterator = OST::ConstIterator;

      // Initialization check
      OST sym( 6UL, 3UL );
      sym(2,1) = 1;
      sym(4,1) = 2;

      checkRows    ( sym, 6UL );
      checkColumns ( sym, 6UL );
      checkCapacity( sym, 4UL );
      checkNonZeros( sym, 4UL );
      checkNonZeros( sym, 0UL, 0UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 1UL );
      checkNonZeros( sym, 3UL, 0UL );
      checkNonZeros( sym, 4UL, 1UL );
      checkNonZeros( sym, 5UL, 0UL );

      // Determining the lower bound for position (1,1)
      {
         ConstIterator pos( sym.lowerBound( 1UL, 1UL ) );

         if( pos == sym.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,1)\n"
                << "   Current matrix:\n" << sym << "\n";
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
                << "   Current matrix:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for position (2,1)
      {
         ConstIterator pos( sym.lowerBound( 2UL, 1UL ) );

         if( pos == sym.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (2,1)\n"
                << "   Current matrix:\n" << sym << "\n";
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
                << "   Current matrix:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for position (3,1)
      {
         ConstIterator pos( sym.lowerBound( 3UL, 1UL ) );

         if( pos == sym.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (3,1)\n"
                << "   Current matrix:\n" << sym << "\n";
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
                << "   Current matrix:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for position (4,1)
      {
         ConstIterator pos( sym.lowerBound( 4UL, 1UL ) );

         if( pos == sym.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (4,1)\n"
                << "   Current matrix:\n" << sym << "\n";
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
                << "   Current matrix:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for position (5,1)
      {
         ConstIterator pos( sym.lowerBound( 5UL, 1UL ) );

         if( pos != sym.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (5,1)\n"
                << "   Current matrix:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c upperBound() member function of the SymmetricMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c upperBound() member function of the SymmetricMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseScalarTest::testUpperBound()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major SymmetricMatrix::upperBound()";

      using ConstIterator = ST::ConstIterator;

      // Initialization check
      ST sym( 6UL, 3UL );
      sym(1,2) = 1;
      sym(1,4) = 2;

      checkRows    ( sym, 6UL );
      checkColumns ( sym, 6UL );
      checkCapacity( sym, 4UL );
      checkNonZeros( sym, 4UL );
      checkNonZeros( sym, 0UL, 0UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 1UL );
      checkNonZeros( sym, 3UL, 0UL );
      checkNonZeros( sym, 4UL, 1UL );
      checkNonZeros( sym, 5UL, 0UL );

      // Determining the upper bound for position (1,1)
      {
         ConstIterator pos( sym.upperBound( 1UL, 1UL ) );

         if( pos == sym.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,1)\n"
                << "   Current matrix:\n" << sym << "\n";
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
                << "   Current matrix:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for position (1,2)
      {
         ConstIterator pos( sym.upperBound( 1UL, 2UL ) );

         if( pos == sym.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,2)\n"
                << "   Current matrix:\n" << sym << "\n";
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
                << "   Current matrix:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for position (1,3)
      {
         ConstIterator pos( sym.upperBound( 1UL, 3UL ) );

         if( pos == sym.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,3)\n"
                << "   Current matrix:\n" << sym << "\n";
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
                << "   Current matrix:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for position (1,4)
      {
         ConstIterator pos( sym.upperBound( 1UL, 4UL ) );

         if( pos != sym.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,4)\n"
                << "   Current matrix:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for position (1,5)
      {
         ConstIterator pos( sym.upperBound( 1UL, 5UL ) );

         if( pos != sym.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,5)\n"
                << "   Current matrix:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major SymmetricMatrix::upperBound()";

      using ConstIterator = OST::ConstIterator;

      // Initialization check
      OST sym( 6UL, 3UL );
      sym(2,1) = 1;
      sym(4,1) = 2;

      checkRows    ( sym, 6UL );
      checkColumns ( sym, 6UL );
      checkCapacity( sym, 4UL );
      checkNonZeros( sym, 4UL );
      checkNonZeros( sym, 0UL, 0UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 1UL );
      checkNonZeros( sym, 3UL, 0UL );
      checkNonZeros( sym, 4UL, 1UL );
      checkNonZeros( sym, 5UL, 0UL );

      // Determining the upper bound for position (1,1)
      {
         ConstIterator pos( sym.upperBound( 1UL, 1UL ) );

         if( pos == sym.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,1)\n"
                << "   Current matrix:\n" << sym << "\n";
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
                << "   Current matrix:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for position (2,1)
      {
         ConstIterator pos( sym.upperBound( 2UL, 1UL ) );

         if( pos == sym.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (2,1)\n"
                << "   Current matrix:\n" << sym << "\n";
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
                << "   Current matrix:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for position (3,1)
      {
         ConstIterator pos( sym.upperBound( 3UL, 1UL ) );

         if( pos == sym.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (3,1)\n"
                << "   Current matrix:\n" << sym << "\n";
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
                << "   Current matrix:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for position (4,1)
      {
         ConstIterator pos( sym.upperBound( 4UL, 1UL ) );

         if( pos != sym.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (4,1)\n"
                << "   Current matrix:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for position (5,1)
      {
         ConstIterator pos( sym.upperBound( 5UL, 1UL ) );

         if( pos != sym.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (5,1)\n"
                << "   Current matrix:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c transpose() member function of the SymmetricMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c transpose() member function of the SymmetricMatrix
// specialization. Additionally, it performs a test of self-transpose via the \c trans()
// function. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseScalarTest::testTranspose()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major self-transpose via transpose()";

      ST sym( 4UL );
      sym(0,0) = 1;
      sym(0,2) = 2;
      sym(0,3) = 3;
      sym(1,1) = 4;
      sym(1,3) = 5;
      sym(2,2) = 6;
      sym(2,3) = 7;

      transpose( sym );

      checkRows    ( sym,  4UL );
      checkColumns ( sym,  4UL );
      checkCapacity( sym, 11UL );
      checkNonZeros( sym, 11UL );
      checkNonZeros( sym,  0UL, 3UL );
      checkNonZeros( sym,  1UL, 2UL );
      checkNonZeros( sym,  2UL, 3UL );
      checkNonZeros( sym,  3UL, 3UL );

      if( sym(0,0) != 1 || sym(0,1) != 0 || sym(0,2) != 2 || sym(0,3) != 3 ||
          sym(1,0) != 0 || sym(1,1) != 4 || sym(1,2) != 0 || sym(1,3) != 5 ||
          sym(2,0) != 2 || sym(2,1) != 0 || sym(2,2) != 6 || sym(2,3) != 7 ||
          sym(3,0) != 3 || sym(3,1) != 5 || sym(3,2) != 7 || sym(3,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Transpose operation failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 1 0 2 3 )\n( 0 4 0 5 )\n( 2 0 6 7 )\n( 3 5 7 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major self-transpose via trans()";

      ST sym( 4UL );
      sym(0,0) = 1;
      sym(0,2) = 2;
      sym(0,3) = 3;
      sym(1,1) = 4;
      sym(1,3) = 5;
      sym(2,2) = 6;
      sym(2,3) = 7;

      sym = trans( sym );

      checkRows    ( sym,  4UL );
      checkColumns ( sym,  4UL );
      checkCapacity( sym, 11UL );
      checkNonZeros( sym, 11UL );
      checkNonZeros( sym,  0UL, 3UL );
      checkNonZeros( sym,  1UL, 2UL );
      checkNonZeros( sym,  2UL, 3UL );
      checkNonZeros( sym,  3UL, 3UL );

      if( sym(0,0) != 1 || sym(0,1) != 0 || sym(0,2) != 2 || sym(0,3) != 3 ||
          sym(1,0) != 0 || sym(1,1) != 4 || sym(1,2) != 0 || sym(1,3) != 5 ||
          sym(2,0) != 2 || sym(2,1) != 0 || sym(2,2) != 6 || sym(2,3) != 7 ||
          sym(3,0) != 3 || sym(3,1) != 5 || sym(3,2) != 7 || sym(3,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Transpose operation failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 1 0 2 3 )\n( 0 4 0 5 )\n( 2 0 6 7 )\n( 3 5 7 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major self-transpose via transpose()";

      OST sym( 4UL );
      sym(0,0) = 1;
      sym(0,2) = 2;
      sym(0,3) = 3;
      sym(1,1) = 4;
      sym(1,3) = 5;
      sym(2,2) = 6;
      sym(2,3) = 7;

      transpose( sym );

      checkRows    ( sym,  4UL );
      checkColumns ( sym,  4UL );
      checkCapacity( sym, 11UL );
      checkNonZeros( sym, 11UL );
      checkNonZeros( sym,  0UL, 3UL );
      checkNonZeros( sym,  1UL, 2UL );
      checkNonZeros( sym,  2UL, 3UL );
      checkNonZeros( sym,  3UL, 3UL );

      if( sym(0,0) != 1 || sym(0,1) != 0 || sym(0,2) != 2 || sym(0,3) != 3 ||
          sym(1,0) != 0 || sym(1,1) != 4 || sym(1,2) != 0 || sym(1,3) != 5 ||
          sym(2,0) != 2 || sym(2,1) != 0 || sym(2,2) != 6 || sym(2,3) != 7 ||
          sym(3,0) != 3 || sym(3,1) != 5 || sym(3,2) != 7 || sym(3,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Transpose operation failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 1 0 2 3 )\n( 0 4 0 5 )\n( 2 0 6 7 )\n( 3 5 7 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major self-transpose via trans()";

      OST sym( 4UL );
      sym(0,0) = 1;
      sym(0,2) = 2;
      sym(0,3) = 3;
      sym(1,1) = 4;
      sym(1,3) = 5;
      sym(2,2) = 6;
      sym(2,3) = 7;

      sym = trans( sym );

      checkRows    ( sym,  4UL );
      checkColumns ( sym,  4UL );
      checkCapacity( sym, 11UL );
      checkNonZeros( sym, 11UL );
      checkNonZeros( sym,  0UL, 3UL );
      checkNonZeros( sym,  1UL, 2UL );
      checkNonZeros( sym,  2UL, 3UL );
      checkNonZeros( sym,  3UL, 3UL );

      if( sym(0,0) != 1 || sym(0,1) != 0 || sym(0,2) != 2 || sym(0,3) != 3 ||
          sym(1,0) != 0 || sym(1,1) != 4 || sym(1,2) != 0 || sym(1,3) != 5 ||
          sym(2,0) != 2 || sym(2,1) != 0 || sym(2,2) != 6 || sym(2,3) != 7 ||
          sym(3,0) != 3 || sym(3,1) != 5 || sym(3,2) != 7 || sym(3,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Transpose operation failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 1 0 2 3 )\n( 0 4 0 5 )\n( 2 0 6 7 )\n( 3 5 7 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c ctranspose() member function of the SymmetricMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c ctranspose() member function of the SymmetricMatrix
// specialization. Additionally, it performs a test of self-transpose via the \c ctrans()
// function. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseScalarTest::testCTranspose()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major self-transpose via ctranspose()";

      ST sym( 4UL );
      sym(0,0) = 1;
      sym(0,2) = 2;
      sym(0,3) = 3;
      sym(1,1) = 4;
      sym(1,3) = 5;
      sym(2,2) = 6;
      sym(2,3) = 7;

      ctranspose( sym );

      checkRows    ( sym,  4UL );
      checkColumns ( sym,  4UL );
      checkCapacity( sym, 11UL );
      checkNonZeros( sym, 11UL );
      checkNonZeros( sym,  0UL, 3UL );
      checkNonZeros( sym,  1UL, 2UL );
      checkNonZeros( sym,  2UL, 3UL );
      checkNonZeros( sym,  3UL, 3UL );

      if( sym(0,0) != 1 || sym(0,1) != 0 || sym(0,2) != 2 || sym(0,3) != 3 ||
          sym(1,0) != 0 || sym(1,1) != 4 || sym(1,2) != 0 || sym(1,3) != 5 ||
          sym(2,0) != 2 || sym(2,1) != 0 || sym(2,2) != 6 || sym(2,3) != 7 ||
          sym(3,0) != 3 || sym(3,1) != 5 || sym(3,2) != 7 || sym(3,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Transpose operation failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 1 0 2 3 )\n( 0 4 0 5 )\n( 2 0 6 7 )\n( 3 5 7 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major self-transpose via ctrans()";

      ST sym( 4UL );
      sym(0,0) = 1;
      sym(0,2) = 2;
      sym(0,3) = 3;
      sym(1,1) = 4;
      sym(1,3) = 5;
      sym(2,2) = 6;
      sym(2,3) = 7;

      sym = ctrans( sym );

      checkRows    ( sym,  4UL );
      checkColumns ( sym,  4UL );
      checkCapacity( sym, 11UL );
      checkNonZeros( sym, 11UL );
      checkNonZeros( sym,  0UL, 3UL );
      checkNonZeros( sym,  1UL, 2UL );
      checkNonZeros( sym,  2UL, 3UL );
      checkNonZeros( sym,  3UL, 3UL );

      if( sym(0,0) != 1 || sym(0,1) != 0 || sym(0,2) != 2 || sym(0,3) != 3 ||
          sym(1,0) != 0 || sym(1,1) != 4 || sym(1,2) != 0 || sym(1,3) != 5 ||
          sym(2,0) != 2 || sym(2,1) != 0 || sym(2,2) != 6 || sym(2,3) != 7 ||
          sym(3,0) != 3 || sym(3,1) != 5 || sym(3,2) != 7 || sym(3,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Transpose operation failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 1 0 2 3 )\n( 0 4 0 5 )\n( 2 0 6 7 )\n( 3 5 7 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major self-transpose via ctranspose()";

      OST sym( 4UL );
      sym(0,0) = 1;
      sym(0,2) = 2;
      sym(0,3) = 3;
      sym(1,1) = 4;
      sym(1,3) = 5;
      sym(2,2) = 6;
      sym(2,3) = 7;

      ctranspose( sym );

      checkRows    ( sym,  4UL );
      checkColumns ( sym,  4UL );
      checkCapacity( sym, 11UL );
      checkNonZeros( sym, 11UL );
      checkNonZeros( sym,  0UL, 3UL );
      checkNonZeros( sym,  1UL, 2UL );
      checkNonZeros( sym,  2UL, 3UL );
      checkNonZeros( sym,  3UL, 3UL );

      if( sym(0,0) != 1 || sym(0,1) != 0 || sym(0,2) != 2 || sym(0,3) != 3 ||
          sym(1,0) != 0 || sym(1,1) != 4 || sym(1,2) != 0 || sym(1,3) != 5 ||
          sym(2,0) != 2 || sym(2,1) != 0 || sym(2,2) != 6 || sym(2,3) != 7 ||
          sym(3,0) != 3 || sym(3,1) != 5 || sym(3,2) != 7 || sym(3,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Transpose operation failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 1 0 2 3 )\n( 0 4 0 5 )\n( 2 0 6 7 )\n( 3 5 7 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major self-transpose via ctrans()";

      OST sym( 4UL );
      sym(0,0) = 1;
      sym(0,2) = 2;
      sym(0,3) = 3;
      sym(1,1) = 4;
      sym(1,3) = 5;
      sym(2,2) = 6;
      sym(2,3) = 7;

      sym = ctrans( sym );

      checkRows    ( sym,  4UL );
      checkColumns ( sym,  4UL );
      checkCapacity( sym, 11UL );
      checkNonZeros( sym, 11UL );
      checkNonZeros( sym,  0UL, 3UL );
      checkNonZeros( sym,  1UL, 2UL );
      checkNonZeros( sym,  2UL, 3UL );
      checkNonZeros( sym,  3UL, 3UL );

      if( sym(0,0) != 1 || sym(0,1) != 0 || sym(0,2) != 2 || sym(0,3) != 3 ||
          sym(1,0) != 0 || sym(1,1) != 4 || sym(1,2) != 0 || sym(1,3) != 5 ||
          sym(2,0) != 2 || sym(2,1) != 0 || sym(2,2) != 6 || sym(2,3) != 7 ||
          sym(3,0) != 3 || sym(3,1) != 5 || sym(3,2) != 7 || sym(3,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Transpose operation failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 1 0 2 3 )\n( 0 4 0 5 )\n( 2 0 6 7 )\n( 3 5 7 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c isDefault() function with the SymmetricMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isDefault() function with the SymmetricMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseScalarTest::testIsDefault()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major isDefault() function";

      // isDefault with 0x0 matrix
      {
         ST sym;

         if( isDefault( sym ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with default matrix
      {
         ST sym( 3UL );

         if( isDefault( sym(0,1) ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix element:\n" << sym(0,1) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( sym ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with non-default matrix
      {
         ST sym( 3UL );
         sym(0,1) = 1;

         if( isDefault( sym(0,1) ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix element:\n" << sym(0,1) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( sym ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << sym << "\n";
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
         OST sym;

         if( isDefault( sym ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with default matrix
      {
         OST sym( 3UL );

         if( isDefault( sym(0,1) ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix element:\n" << sym(0,1) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( sym ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with non-default matrix
      {
         OST sym( 3UL );
         sym(1,0) = 1;

         if( isDefault( sym(0,1) ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix element:\n" << sym(0,1) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( sym ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c submatrix() function with the SymmetricMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c submatrix() function with the SymmetricMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseScalarTest::testSubmatrix()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major submatrix() function";

      using SMT = blaze::Submatrix<ST>;

      ST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      SMT sm = submatrix( sym, 0UL, 1UL, 2UL, 2UL );

      if( sm(0,1) != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator access failed\n"
             << " Details:\n"
             << "   Result: " << sm(0,1) << "\n"
             << "   Expected result: 7\n";
         throw std::runtime_error( oss.str() );
      }

      SMT::Iterator it = sm.begin(0UL);

      if( it == sm.end(0UL) || it->value() != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << it->value() << "\n"
             << "   Expected result: -4\n";
         throw std::runtime_error( oss.str() );
      }

      sm(1,1) = -5;

      if( sm(0,0) != -4 || sm(0,1) !=  7 ||
          sm(1,0) !=  2 || sm(1,1) != -5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix access failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( -4  7 )\n(  2 -5 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( sym(0,0) !=  1 || sym(0,1) != -4 || sym(0,2) !=  7 ||
          sym(1,0) != -4 || sym(1,1) !=  2 || sym(1,2) != -5 ||
          sym(2,0) !=  7 || sym(2,1) != -5 || sym(2,2) !=  3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix access failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  1 -4  7 )\n( -4  2 -5 )\n(  7 -5  3 )\n";
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

      if( sym(0,0) != 1 || sym(0,1) != 0 || sym(0,2) != 0 ||
          sym(1,0) != 0 || sym(1,1) != 0 || sym(1,2) != 0 ||
          sym(2,0) != 0 || sym(2,1) != 0 || sym(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix reset failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major submatrix() function";

      using SMT = blaze::Submatrix<OST>;

      OST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      SMT sm = submatrix( sym, 0UL, 1UL, 2UL, 2UL );

      if( sm(0,1) != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator access failed\n"
             << " Details:\n"
             << "   Result: " << sm(0,1) << "\n"
             << "   Expected result: 7\n";
         throw std::runtime_error( oss.str() );
      }

      SMT::Iterator it = sm.begin(0UL);

      if( it == sm.end(0UL) || it->value() != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << it->value() << "\n"
             << "   Expected result: -4\n";
         throw std::runtime_error( oss.str() );
      }

      sm(1,1) = -5;

      if( sm(0,0) != -4 || sm(0,1) !=  7 ||
          sm(1,0) !=  2 || sm(1,1) != -5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix access failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( -4  7 )\n(  2 -5 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( sym(0,0) !=  1 || sym(0,1) != -4 || sym(0,2) !=  7 ||
          sym(1,0) != -4 || sym(1,1) !=  2 || sym(1,2) != -5 ||
          sym(2,0) !=  7 || sym(2,1) != -5 || sym(2,2) !=  3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix access failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  1 -4  7 )\n( -4  2 -5 )\n(  7 -5  3 )\n";
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

      if( sym(0,0) != 1 || sym(0,1) != 0 || sym(0,2) != 0 ||
          sym(1,0) != 0 || sym(1,1) != 0 || sym(1,2) != 0 ||
          sym(2,0) != 0 || sym(2,1) != 0 || sym(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix reset failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c row() function with the SymmetricMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c row() function with the SymmetricMatrix specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseScalarTest::testRow()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major row() function";

      using RT = blaze::Row<ST>;

      ST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      RT row1 = row( sym, 1UL );

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

      if( it == row1.end() || it->value() != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << it->value() << "\n"
             << "   Expected result: -4\n";
         throw std::runtime_error( oss.str() );
      }

      row1[2] = -5;

      if( row1[0] != -4 || row1[1] != 2 || row1[2] != -5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row access failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( -4 2 -5 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( sym(0,0) !=  1 || sym(0,1) != -4 || sym(0,2) !=  7 ||
          sym(1,0) != -4 || sym(1,1) !=  2 || sym(1,2) != -5 ||
          sym(2,0) !=  7 || sym(2,1) != -5 || sym(2,2) !=  3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row access failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  1 -4  7 )\n( -4  2 -5 )\n(  7 -5  3 )\n";
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

      if( sym(0,0) != 1 || sym(0,1) != 0 || sym(0,2) != 7 ||
          sym(1,0) != 0 || sym(1,1) != 0 || sym(1,2) != 0 ||
          sym(2,0) != 7 || sym(2,1) != 0 || sym(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row reset failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 1 0 7 )\n( 0 0 0 )\n( 7 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major row() function";

      using RT = blaze::Row<OST>;

      OST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      RT row1 = row( sym, 1UL );

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

      if( it == row1.end() || it->value() != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << it->value() << "\n"
             << "   Expected result: -4\n";
         throw std::runtime_error( oss.str() );
      }

      row1[2] = -5;

      if( row1[0] != -4 || row1[1] != 2 || row1[2] != -5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row access failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( -4 2 -5 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( sym(0,0) !=  1 || sym(0,1) != -4 || sym(0,2) !=  7 ||
          sym(1,0) != -4 || sym(1,1) !=  2 || sym(1,2) != -5 ||
          sym(2,0) !=  7 || sym(2,1) != -5 || sym(2,2) !=  3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row access failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  1 -4  7 )\n( -4  2 -5 )\n(  7 -5  3 )\n";
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

      if( sym(0,0) != 1 || sym(0,1) != 0 || sym(0,2) != 7 ||
          sym(1,0) != 0 || sym(1,1) != 0 || sym(1,2) != 0 ||
          sym(2,0) != 7 || sym(2,1) != 0 || sym(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row reset failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 1 0 7 )\n( 0 0 0 )\n( 7 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c column() function with the SymmetricMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c column() function with the SymmetricMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseScalarTest::testColumn()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major column() function";

      using CT = blaze::Column<ST>;

      ST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      CT col1 = column( sym, 1UL );

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

      col1[2] = -5;

      if( col1[0] != -4 || col1[1] != 2 || col1[2] != -5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column access failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( -4 2 -5 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( sym(0,0) !=  1 || sym(0,1) != -4 || sym(0,2) !=  7 ||
          sym(1,0) != -4 || sym(1,1) !=  2 || sym(1,2) != -5 ||
          sym(2,0) !=  7 || sym(2,1) != -5 || sym(2,2) !=  3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column access failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  1 -4  7 )\n( -4  2 -5 )\n(  7 -5  3 )\n";
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

      if( sym(0,0) != 1 || sym(0,1) != 0 || sym(0,2) != 7 ||
          sym(1,0) != 0 || sym(1,1) != 0 || sym(1,2) != 0 ||
          sym(2,0) != 7 || sym(2,1) != 0 || sym(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column reset failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 1 0 7 )\n( 0 0 0 )\n( 7 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major column() function";

      using CT = blaze::Column<OST>;

      OST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      CT col1 = column( sym, 1UL );

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

      col1[2] = -5;

      if( col1[0] != -4 || col1[1] != 2 || col1[2] != -5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column access failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( -4 2 -5 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( sym(0,0) !=  1 || sym(0,1) != -4 || sym(0,2) !=  7 ||
          sym(1,0) != -4 || sym(1,1) !=  2 || sym(1,2) != -5 ||
          sym(2,0) !=  7 || sym(2,1) != -5 || sym(2,2) !=  3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column access failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  1 -4  7 )\n( -4  2 -5 )\n(  7 -5  3 )\n";
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

      if( sym(0,0) != 1 || sym(0,1) != 0 || sym(0,2) != 7 ||
          sym(1,0) != 0 || sym(1,1) != 0 || sym(1,2) != 0 ||
          sym(2,0) != 7 || sym(2,1) != 0 || sym(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column reset failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 1 0 7 )\n( 0 0 0 )\n( 7 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************

} // namespace symmetricmatrix

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
   std::cout << "   Running SymmetricMatrix sparse scalar test (part 2)..." << std::endl;

   try
   {
      RUN_SYMMETRICMATRIX_SPARSESCALAR_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during SymmetricMatrix sparse scalar test (part 2):\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
