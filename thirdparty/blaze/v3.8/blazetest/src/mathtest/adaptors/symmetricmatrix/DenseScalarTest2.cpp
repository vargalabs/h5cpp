//=================================================================================================
/*!
//  \file src/mathtest/adaptors/symmetricmatrix/DenseScalarTest2.cpp
//  \brief Source file for the SymmetricMatrix dense scalar test (part 2)
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
#include <blazetest/mathtest/adaptors/symmetricmatrix/DenseScalarTest.h>

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
/*!\brief Constructor for the SymmetricMatrix dense scalar test.
//
// \exception std::runtime_error Operation error detected.
*/
DenseScalarTest::DenseScalarTest()
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
void DenseScalarTest::testScaling()
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
      checkCapacity( sym, 9UL );
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
      checkCapacity( sym, 9UL );
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
      checkCapacity( sym, 9UL );
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
      checkCapacity( sym, 9UL );
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
      checkCapacity( sym, 9UL );
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
      checkCapacity( sym, 9UL );
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
      checkCapacity( sym, 9UL );
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
      checkCapacity( sym, 9UL );
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

      blaze::SymmetricMatrix< blaze::DynamicMatrix<complex<float>,blaze::rowMajor> > sym( 2UL );
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
      checkCapacity( sym, 9UL );
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
      checkCapacity( sym, 9UL );
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
      checkCapacity( sym, 9UL );
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
      checkCapacity( sym, 9UL );
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
      checkCapacity( sym, 9UL );
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
      checkCapacity( sym, 9UL );
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
      checkCapacity( sym, 9UL );
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
      checkCapacity( sym, 9UL );
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

      blaze::SymmetricMatrix< blaze::DynamicMatrix<complex<float>,blaze::columnMajor> > sym( 2UL );
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
void DenseScalarTest::testFunctionCall()
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
      checkCapacity( sym, 9UL );
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
      checkCapacity( sym, 9UL );
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
      checkCapacity( sym, 9UL );
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
      checkCapacity( sym, 9UL );
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
      checkCapacity( sym, 9UL );
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
      checkCapacity( sym, 9UL );
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
      checkCapacity( sym, 9UL );
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
      checkCapacity( sym, 9UL );
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
      checkCapacity( sym, 9UL );
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
      checkCapacity( sym, 9UL );
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
      checkCapacity( sym, 9UL );
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
      checkCapacity( sym, 9UL );
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
      checkCapacity( sym, 9UL );
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
      checkCapacity( sym, 9UL );
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
void DenseScalarTest::testIterator()
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

         if( it == end( sym, 1UL ) || *it != 1 ) {
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

         const ptrdiff_t number( begin( sym, 0UL ) - end( sym, 0UL ) );

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

         const ptrdiff_t number( cend( sym, 1UL ) - cbegin( sym, 1UL ) );

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

         const ptrdiff_t number( cbegin( sym, 1UL ) - cend( sym, 1UL ) );

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

         ConstIterator it ( cbegin( sym, 2UL ) );
         ConstIterator end( cend( sym, 2UL ) );

         if( it == end || *it != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid initial iterator detected\n";
            throw std::runtime_error( oss.str() );
         }

         ++it;

         if( it == end || *it != -2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator pre-increment failed\n";
            throw std::runtime_error( oss.str() );
         }

         --it;

         if( it == end || *it != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator pre-decrement failed\n";
            throw std::runtime_error( oss.str() );
         }

         it++;

         if( it == end || *it != -2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator post-increment failed\n";
            throw std::runtime_error( oss.str() );
         }

         it--;

         if( it == end || *it != 0 ) {
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

         if( it == end || *it != 0 ) {
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

         if( it == end || *it != 0 ) {
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

      // Testing assignment via Iterator
      {
         test_ = "Row-major assignment via Iterator";

         int value = 7;

         for( Iterator it=begin( sym, 2UL ); it!=end( sym, 2UL ); ++it ) {
            *it = value++;
         }

         if( sym(0,0) != 0 || sym(0,1) != 1 || sym(0,2) != 7 ||
             sym(1,0) != 1 || sym(1,1) != 0 || sym(1,2) != 8 ||
             sym(2,0) != 7 || sym(2,1) != 8 || sym(2,2) != 9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 0 1 7 )\n( 1 0 8 )\n( 7 8 9 )\n";
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

         if( sym(0,0) != 0 || sym(0,1) !=  5 || sym(0,2) !=  7 ||
             sym(1,0) != 5 || sym(1,1) !=  5 || sym(1,2) != 14 ||
             sym(2,0) != 7 || sym(2,1) != 14 || sym(2,2) !=  9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 0  5  7 )\n( 5  5 14 )\n( 7 14  9 )\n";
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

         if( sym(0,0) != 0 || sym(0,1) != 1 || sym(0,2) != 7 ||
             sym(1,0) != 1 || sym(1,1) != 0 || sym(1,2) != 8 ||
             sym(2,0) != 7 || sym(2,1) != 8 || sym(2,2) != 9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 0 1 7 )\n( 1 0 8 )\n( 7 8 9 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing multiplication assignment via Iterator
      {
         test_ = "Row-major multiplication assignment via Iterator";

         int value = 2;

         for( Iterator it=begin( sym, 1UL ); it!=end( sym, 1UL ); ++it ) {
            *it *= value++;
         }

         if( sym(0,0) != 0 || sym(0,1) !=  2 || sym(0,2) !=  7 ||
             sym(1,0) != 2 || sym(1,1) !=  0 || sym(1,2) != 32 ||
             sym(2,0) != 7 || sym(2,1) != 32 || sym(2,2) !=  9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 0  2  7 )\n( 2  0 32 )\n( 7 32  9 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing division assignment via Iterator
      {
         test_ = "Row-major division assignment via Iterator";

         for( Iterator it=begin( sym, 1UL ); it!=end( sym, 1UL ); ++it ) {
            *it /= 2;
         }

         if( sym(0,0) != 0 || sym(0,1) !=  1 || sym(0,2) !=  7 ||
             sym(1,0) != 1 || sym(1,1) !=  0 || sym(1,2) != 16 ||
             sym(2,0) != 7 || sym(2,1) != 16 || sym(2,2) !=  9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 0  1  7 )\n( 1  0 16 )\n( 7 16  9 )\n";
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

         ConstIterator it( begin( sym, 1UL ) );

         if( it == end( sym, 1UL ) || *it != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Failed iterator conversion detected\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Counting the number of elements in 0th column via Iterator (end-begin)
      {
         test_ = "Column-major Iterator subtraction (end-begin)";

         const ptrdiff_t number( end( sym, 0UL ) - begin( sym, 0UL ) );

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

         const ptrdiff_t number( begin( sym, 0UL ) - end( sym, 0UL ) );

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

         const ptrdiff_t number( cend( sym, 1UL ) - cbegin( sym, 1UL ) );

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

         const ptrdiff_t number( cbegin( sym, 1UL ) - cend( sym, 1UL ) );

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

         ConstIterator it ( cbegin( sym, 2UL ) );
         ConstIterator end( cend( sym, 2UL ) );

         if( it == end || *it != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid initial iterator detected\n";
            throw std::runtime_error( oss.str() );
         }

         ++it;

         if( it == end || *it != -2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator pre-increment failed\n";
            throw std::runtime_error( oss.str() );
         }

         --it;

         if( it == end || *it != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator pre-decrement failed\n";
            throw std::runtime_error( oss.str() );
         }

         it++;

         if( it == end || *it != -2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator post-increment failed\n";
            throw std::runtime_error( oss.str() );
         }

         it--;

         if( it == end || *it != 0 ) {
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

         if( it == end || *it != 0 ) {
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

         if( it == end || *it != 0 ) {
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

      // Testing assignment via Iterator
      {
         test_ = "Column-major assignment via Iterator";

         int value = 7;

         for( Iterator it=begin( sym, 2UL ); it!=end( sym, 2UL ); ++it ) {
            *it = value++;
         }

         if( sym(0,0) != 0 || sym(0,1) != 1 || sym(0,2) != 7 ||
             sym(1,0) != 1 || sym(1,1) != 0 || sym(1,2) != 8 ||
             sym(2,0) != 7 || sym(2,1) != 8 || sym(2,2) != 9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 0 1 7 )\n( 1 0 8 )\n( 7 8 9 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing addition assignment via Iterator
      {
         test_ = "Column-major addition assignment via Iterator";

         int value = 4;

         for( Iterator it=begin( sym, 1UL ); it!=end( sym, 1UL ); ++it ) {
            *it += value++;
         }

         if( sym(0,0) != 0 || sym(0,1) !=  5 || sym(0,2) !=  7 ||
             sym(1,0) != 5 || sym(1,1) !=  5 || sym(1,2) != 14 ||
             sym(2,0) != 7 || sym(2,1) != 14 || sym(2,2) !=  9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 0  5  7 )\n( 5  5 14 )\n( 7 14  9 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing subtraction assignment via Iterator
      {
         test_ = "Column-major subtraction assignment via Iterator";

         int value = 4;

         for( Iterator it=begin( sym, 1UL ); it!=end( sym, 1UL ); ++it ) {
            *it -= value++;
         }

         if( sym(0,0) != 0 || sym(0,1) != 1 || sym(0,2) != 7 ||
             sym(1,0) != 1 || sym(1,1) != 0 || sym(1,2) != 8 ||
             sym(2,0) != 7 || sym(2,1) != 8 || sym(2,2) != 9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 0 1 7 )\n( 1 0 8 )\n( 7 8 9 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing multiplication assignment via Iterator
      {
         test_ = "Column-major multiplication assignment via Iterator";

         int value = 2;

         for( Iterator it=begin( sym, 1UL ); it!=end( sym, 1UL ); ++it ) {
            *it *= value++;
         }

         if( sym(0,0) != 0 || sym(0,1) !=  2 || sym(0,2) !=  7 ||
             sym(1,0) != 2 || sym(1,1) !=  0 || sym(1,2) != 32 ||
             sym(2,0) != 7 || sym(2,1) != 32 || sym(2,2) !=  9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 0  2  7 )\n( 2  0 32 )\n( 7 32  9 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing division assignment via Iterator
      {
         test_ = "Column-major division assignment via Iterator";

         for( Iterator it=begin( sym, 1UL ); it!=end( sym, 1UL ); ++it ) {
            *it /= 2;
         }

         if( sym(0,0) != 0 || sym(0,1) !=  1 || sym(0,2) !=  7 ||
             sym(1,0) != 1 || sym(1,1) !=  0 || sym(1,2) != 16 ||
             sym(2,0) != 7 || sym(2,1) != 16 || sym(2,2) !=  9 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 0  1  7 )\n( 1  0 16 )\n( 7 16  9 )\n";
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
void DenseScalarTest::testNonZeros()
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
         checkCapacity( sym, 9UL );
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
         checkCapacity( sym, 9UL );
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
void DenseScalarTest::testReset()
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
void DenseScalarTest::testClear()
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
void DenseScalarTest::testResize()
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
      checkCapacity( sym, 4UL );
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

      checkRows    ( sym,  4UL );
      checkColumns ( sym,  4UL );
      checkCapacity( sym, 16UL );
      checkNonZeros( sym,  3UL );
      checkNonZeros( sym,  0UL, 1UL );
      checkNonZeros( sym,  1UL, 2UL );
      checkNonZeros( sym,  2UL, 0UL );
      checkNonZeros( sym,  3UL, 0UL );

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
      checkCapacity( sym, 4UL );
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
      checkCapacity( sym, 4UL );
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

      checkRows    ( sym,  4UL );
      checkColumns ( sym,  4UL );
      checkCapacity( sym, 16UL );
      checkNonZeros( sym,  3UL );
      checkNonZeros( sym,  0UL, 1UL );
      checkNonZeros( sym,  1UL, 2UL );
      checkNonZeros( sym,  2UL, 0UL );
      checkNonZeros( sym,  3UL, 0UL );

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
      checkCapacity( sym, 4UL );
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
/*!\brief Test of the \c extend() member function of the SymmetricMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c extend() member function of the SymmetricMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseScalarTest::testExtend()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major SymmetricMatrix::extend()";

      // Initialization check
      ST sym;

      checkRows    ( sym, 0UL );
      checkColumns ( sym, 0UL );
      checkNonZeros( sym, 0UL );

      // Extending the size of the matrix to 2x2
      sym.extend( 2UL );

      checkRows    ( sym, 2UL );
      checkColumns ( sym, 2UL );
      checkCapacity( sym, 4UL );
      checkNonZeros( sym, 0UL );
      checkNonZeros( sym, 0UL, 0UL );
      checkNonZeros( sym, 1UL, 0UL );

      if( sym(0,0) != 0 || sym(0,1) != 0 ||
          sym(1,0) != 0 || sym(1,1) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Extending the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 0 0 )\n( 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Extending to 4x4 and preserving the elements
      sym(0,1) = 1;
      sym(1,1) = 2;
      sym.extend( 2UL, true );

      checkRows    ( sym,  4UL );
      checkColumns ( sym,  4UL );
      checkCapacity( sym, 16UL );
      checkNonZeros( sym,  3UL );
      checkNonZeros( sym,  0UL, 1UL );
      checkNonZeros( sym,  1UL, 2UL );
      checkNonZeros( sym,  2UL, 0UL );
      checkNonZeros( sym,  3UL, 0UL );

      if( sym(0,0) != 0 || sym(0,1) != 1 || sym(0,2) != 0 || sym(0,3) != 0 ||
          sym(1,0) != 1 || sym(1,1) != 2 || sym(1,2) != 0 || sym(1,3) != 0 ||
          sym(2,0) != 0 || sym(2,1) != 0 || sym(2,2) != 0 || sym(2,3) != 0 ||
          sym(3,0) != 0 || sym(3,1) != 0 || sym(3,2) != 0 || sym(3,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Extending the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 0 1 0 0 )\n( 1 2 0 0 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major SymmetricMatrix::extend()";

      // Initialization check
      OST sym;

      checkRows    ( sym, 0UL );
      checkColumns ( sym, 0UL );
      checkNonZeros( sym, 0UL );

      // Extending the size of the matrix to 2x2
      sym.extend( 2UL );

      checkRows    ( sym, 2UL );
      checkColumns ( sym, 2UL );
      checkCapacity( sym, 4UL );
      checkNonZeros( sym, 0UL );
      checkNonZeros( sym, 0UL, 0UL );
      checkNonZeros( sym, 1UL, 0UL );

      if( sym(0,0) != 0 || sym(0,1) != 0 ||
          sym(1,0) != 0 || sym(1,1) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Extending the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 0 0 )\n( 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Extending to 4x4 and preserving the elements
      sym(0,1) = 1;
      sym(1,1) = 2;
      sym.extend( 2UL, true );

      checkRows    ( sym,  4UL );
      checkColumns ( sym,  4UL );
      checkCapacity( sym, 16UL );
      checkNonZeros( sym,  3UL );
      checkNonZeros( sym,  0UL, 1UL );
      checkNonZeros( sym,  1UL, 2UL );
      checkNonZeros( sym,  2UL, 0UL );
      checkNonZeros( sym,  3UL, 0UL );

      if( sym(0,0) != 0 || sym(0,1) != 1 || sym(0,2) != 0 || sym(0,3) != 0 ||
          sym(1,0) != 1 || sym(1,1) != 2 || sym(1,2) != 0 || sym(1,3) != 0 ||
          sym(2,0) != 0 || sym(2,1) != 0 || sym(2,2) != 0 || sym(2,3) != 0 ||
          sym(3,0) != 0 || sym(3,1) != 0 || sym(3,2) != 0 || sym(3,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Extending the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 0 1 0 0 )\n( 1 2 0 0 )\n( 0 0 0 0 )\n( 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
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
void DenseScalarTest::testReserve()
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
/*!\brief Test of the \c shrinkToFit() member function of the SymmetricMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c shrinkToFit() member function of the SymmetricMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseScalarTest::testShrinkToFit()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major SymmetricMatrix::shrinkToFit()";

      // Shrinking a matrix without excessive capacity
      {
         ST sym{ { 1, 2, 3 }, { 2, 4, 5 }, { 3, 5, 6 } };

         sym.shrinkToFit();

         checkRows    ( sym, 3UL );
         checkColumns ( sym, 3UL );
         checkCapacity( sym, 9UL );
         checkNonZeros( sym, 9UL );
         checkNonZeros( sym, 0UL, 3UL );
         checkNonZeros( sym, 1UL, 3UL );
         checkNonZeros( sym, 2UL, 3UL );

         if( sym.capacity() != sym.rows() * sym.spacing() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Shrinking the matrix failed\n"
                << " Details:\n"
                << "   Capacity         : " << sym.capacity() << "\n"
                << "   Expected capacity: " << ( sym.rows() * sym.spacing() ) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) != 1 || sym(0,1) != 2 || sym(0,2) != 3 ||
             sym(1,0) != 2 || sym(1,1) != 4 || sym(1,2) != 5 ||
             sym(2,0) != 3 || sym(2,1) != 5 || sym(2,2) != 6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Shrinking the matrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 1 2 3 )\n( 2 4 5 )\n( 3 5 6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Shrinking a matrix with excessive capacity
      {
         ST sym{ { 1, 2, 3 }, { 2, 4, 5 }, { 3, 5, 6 } };
         sym.reserve( 100UL );

         sym.shrinkToFit();

         checkRows    ( sym, 3UL );
         checkColumns ( sym, 3UL );
         checkCapacity( sym, 9UL );
         checkNonZeros( sym, 9UL );
         checkNonZeros( sym, 0UL, 3UL );
         checkNonZeros( sym, 1UL, 3UL );
         checkNonZeros( sym, 2UL, 3UL );

         if( sym.capacity() != sym.rows() * sym.spacing() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Shrinking the matrix failed\n"
                << " Details:\n"
                << "   Capacity         : " << sym.capacity() << "\n"
                << "   Expected capacity: " << ( sym.rows() * sym.spacing() ) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) != 1 || sym(0,1) != 2 || sym(0,2) != 3 ||
             sym(1,0) != 2 || sym(1,1) != 4 || sym(1,2) != 5 ||
             sym(2,0) != 3 || sym(2,1) != 5 || sym(2,2) != 6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Shrinking the matrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 1 2 3 )\n( 2 4 5 )\n( 3 5 6 )\n";
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
         OST sym{ { 1, 2, 3 }, { 2, 4, 5 }, { 3, 5, 6 } };

         sym.shrinkToFit();

         checkRows    ( sym, 3UL );
         checkColumns ( sym, 3UL );
         checkCapacity( sym, 9UL );
         checkNonZeros( sym, 9UL );
         checkNonZeros( sym, 0UL, 3UL );
         checkNonZeros( sym, 1UL, 3UL );
         checkNonZeros( sym, 2UL, 3UL );

         if( sym.capacity() != sym.spacing() * sym.columns() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Shrinking the matrix failed\n"
                << " Details:\n"
                << "   Capacity         : " << sym.capacity() << "\n"
                << "   Expected capacity: " << ( sym.spacing() * sym.columns() ) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) != 1 || sym(0,1) != 2 || sym(0,2) != 3 ||
             sym(1,0) != 2 || sym(1,1) != 4 || sym(1,2) != 5 ||
             sym(2,0) != 3 || sym(2,1) != 5 || sym(2,2) != 6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Shrinking the matrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 1 2 3 )\n( 2 4 5 )\n( 3 5 6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Shrinking a matrix with excessive capacity
      {
         OST sym{ { 1, 2, 3 }, { 2, 4, 5 }, { 3, 5, 6 } };
         sym.reserve( 100UL );

         sym.shrinkToFit();

         checkRows    ( sym, 3UL );
         checkColumns ( sym, 3UL );
         checkCapacity( sym, 9UL );
         checkNonZeros( sym, 9UL );
         checkNonZeros( sym, 0UL, 3UL );
         checkNonZeros( sym, 1UL, 3UL );
         checkNonZeros( sym, 2UL, 3UL );

         if( sym.capacity() != sym.spacing() * sym.columns() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Shrinking the matrix failed\n"
                << " Details:\n"
                << "   Capacity         : " << sym.capacity() << "\n"
                << "   Expected capacity: " << ( sym.spacing() * sym.columns() ) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) != 1 || sym(0,1) != 2 || sym(0,2) != 3 ||
             sym(1,0) != 2 || sym(1,1) != 4 || sym(1,2) != 5 ||
             sym(2,0) != 3 || sym(2,1) != 5 || sym(2,2) != 6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Shrinking the matrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 1 2 3 )\n( 2 4 5 )\n( 3 5 6 )\n";
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
void DenseScalarTest::testSwap()
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
/*!\brief Test of the \c transpose() member function of the SymmetricMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c transpose() member function of the SymmetricMatrix
// specialization. Additionally, it performs a test of self-transpose via the \c trans()
// function. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseScalarTest::testTranspose()
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
      checkCapacity( sym, 16UL );
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
      checkCapacity( sym, 16UL );
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
      checkCapacity( sym, 16UL );
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
      checkCapacity( sym, 16UL );
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
void DenseScalarTest::testCTranspose()
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
      checkCapacity( sym, 16UL );
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
      checkCapacity( sym, 16UL );
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
      checkCapacity( sym, 16UL );
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
      checkCapacity( sym, 16UL );
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
void DenseScalarTest::testIsDefault()
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

         if( isDefault( sym(1,0) ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix element:\n" << sym(1,0) << "\n";
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

         if( isDefault( sym(1,0) ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix element:\n" << sym(1,0) << "\n";
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
void DenseScalarTest::testSubmatrix()
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

      if( it == sm.end(0UL) || *it != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << *it << "\n"
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

      if( it == sm.end(0UL) || *it != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << *it << "\n"
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
void DenseScalarTest::testRow()
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

      if( it == row1.end() || *it != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << *it << "\n"
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

      if( it == row1.end() || *it != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << *it << "\n"
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
void DenseScalarTest::testColumn()
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

      if( it == col1.end() || *it != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << *it << "\n"
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

      if( it == col1.end() || *it != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << *it << "\n"
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
   std::cout << "   Running SymmetricMatrix dense scalar test (part 2)..." << std::endl;

   try
   {
      RUN_SYMMETRICMATRIX_DENSESCALAR_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during SymmetricMatrix dense scalar test (part 2):\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
