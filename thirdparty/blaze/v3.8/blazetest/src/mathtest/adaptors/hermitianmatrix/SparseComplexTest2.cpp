//=================================================================================================
/*!
//  \file src/mathtest/adaptors/hermitianmatrix/SparseComplexTest2.cpp
//  \brief Source file for the HermitianMatrix sparse complex test (part 2)
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
#include <blazetest/mathtest/adaptors/hermitianmatrix/SparseComplexTest.h>

#ifdef BLAZE_USE_HPX_THREADS
#  include <hpx/hpx_main.hpp>
#endif


namespace blazetest {

namespace mathtest {

namespace adaptors {

namespace hermitianmatrix {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the HermitianMatrix sparse complex test.
//
// \exception std::runtime_error Operation error detected.
*/
SparseComplexTest::SparseComplexTest()
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
/*!\brief Test of all HermitianMatrix (self-)scaling operations.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all available ways to scale an instance of the HermitianMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseComplexTest::testScaling()
{
   //=====================================================================================
   // Row-major self-scaling (M*=s)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (M*=s)";

      HT herm( 3UL );
      herm(1,2) = cplx( 1,-2);
      herm(2,0) = cplx(-2, 0);
      herm(2,2) = cplx( 3, 0);

      herm *= 2;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 5UL );
      checkNonZeros( herm, 5UL );
      checkNonZeros( herm, 0UL, 1UL );
      checkNonZeros( herm, 1UL, 1UL );
      checkNonZeros( herm, 2UL, 3UL );

      if( herm(0,0) != cplx( 0,0) || herm(0,1) != cplx(0,0) || herm(0,2) != cplx(-4, 0) ||
          herm(1,0) != cplx( 0,0) || herm(1,1) != cplx(0,0) || herm(1,2) != cplx( 2,-4) ||
          herm(2,0) != cplx(-4,0) || herm(2,1) != cplx(2,4) || herm(2,2) != cplx( 6, 0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( ( 0,0) (0,0) (-4, 0) )\n"
                                     "( ( 0,0) (0,0) ( 2,-4) )\n"
                                     "( (-4,0) (2,4) ( 6, 0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major self-scaling (M=M*s)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (M=M*s)";

      HT herm( 3UL );
      herm(1,2) = cplx( 1,-2);
      herm(2,0) = cplx(-2, 0);
      herm(2,2) = cplx( 3, 0);

      herm = herm * 2;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 5UL );
      checkNonZeros( herm, 5UL );
      checkNonZeros( herm, 0UL, 1UL );
      checkNonZeros( herm, 1UL, 1UL );
      checkNonZeros( herm, 2UL, 3UL );

      if( herm(0,0) != cplx( 0,0) || herm(0,1) != cplx(0,0) || herm(0,2) != cplx(-4, 0) ||
          herm(1,0) != cplx( 0,0) || herm(1,1) != cplx(0,0) || herm(1,2) != cplx( 2,-4) ||
          herm(2,0) != cplx(-4,0) || herm(2,1) != cplx(2,4) || herm(2,2) != cplx( 6, 0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( ( 0,0) (0,0) (-4, 0) )\n"
                                     "( ( 0,0) (0,0) ( 2,-4) )\n"
                                     "( (-4,0) (2,4) ( 6, 0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major self-scaling (M=s*M)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (M=s*M)";

      HT herm( 3UL );
      herm(1,2) = cplx( 1,-2);
      herm(2,0) = cplx(-2, 0);
      herm(2,2) = cplx( 3, 0);

      herm = 2 * herm;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 5UL );
      checkNonZeros( herm, 5UL );
      checkNonZeros( herm, 0UL, 1UL );
      checkNonZeros( herm, 1UL, 1UL );
      checkNonZeros( herm, 2UL, 3UL );

      if( herm(0,0) != cplx( 0,0) || herm(0,1) != cplx(0,0) || herm(0,2) != cplx(-4, 0) ||
          herm(1,0) != cplx( 0,0) || herm(1,1) != cplx(0,0) || herm(1,2) != cplx( 2,-4) ||
          herm(2,0) != cplx(-4,0) || herm(2,1) != cplx(2,4) || herm(2,2) != cplx( 6, 0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( ( 0,0) (0,0) (-4, 0) )\n"
                                     "( ( 0,0) (0,0) ( 2,-4) )\n"
                                     "( (-4,0) (2,4) ( 6, 0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major self-scaling (M/=s)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (M/=s)";

      HT herm( 3UL );
      herm(1,2) = cplx( 2,-4);
      herm(2,0) = cplx(-4, 0);
      herm(2,2) = cplx( 6, 0);

      herm /= 2;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 5UL );
      checkNonZeros( herm, 5UL );
      checkNonZeros( herm, 0UL, 1UL );
      checkNonZeros( herm, 1UL, 1UL );
      checkNonZeros( herm, 2UL, 3UL );

      if( herm(0,0) != cplx( 0,0) || herm(0,1) != cplx(0,0) || herm(0,2) != cplx(-2, 0) ||
          herm(1,0) != cplx( 0,0) || herm(1,1) != cplx(0,0) || herm(1,2) != cplx( 1,-2) ||
          herm(2,0) != cplx(-2,0) || herm(2,1) != cplx(1,2) || herm(2,2) != cplx( 3, 0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( ( 0,0) (0,0) (-2, 0) )\n"
                                     "( ( 0,0) (0,0) ( 1,-2) )\n"
                                     "( (-2,0) (1,2) ( 3, 0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major self-scaling (M=M/s)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (M=M/s)";

      HT herm( 3UL );
      herm(1,2) = cplx( 2,-4);
      herm(2,0) = cplx(-4, 0);
      herm(2,2) = cplx( 6, 0);

      herm = herm / 2;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 5UL );
      checkNonZeros( herm, 5UL );
      checkNonZeros( herm, 0UL, 1UL );
      checkNonZeros( herm, 1UL, 1UL );
      checkNonZeros( herm, 2UL, 3UL );

      if( herm(0,0) != cplx( 0,0) || herm(0,1) != cplx(0,0) || herm(0,2) != cplx(-2, 0) ||
          herm(1,0) != cplx( 0,0) || herm(1,1) != cplx(0,0) || herm(1,2) != cplx( 1,-2) ||
          herm(2,0) != cplx(-2,0) || herm(2,1) != cplx(1,2) || herm(2,2) != cplx( 3, 0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( ( 0,0) (0,0) (-2, 0) )\n"
                                     "( ( 0,0) (0,0) ( 1,-2) )\n"
                                     "( (-2,0) (1,2) ( 3, 0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major HermitianMatrix::scale()
   //=====================================================================================

   {
      test_ = "Row-major HermitianMatrix::scale()";

      // Initialization check
      HT herm( 3UL );
      herm(1,2) = cplx( 1,-2);
      herm(2,0) = cplx(-2, 0);
      herm(2,2) = cplx( 3, 0);

      herm.scale( 2 );

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 5UL );
      checkNonZeros( herm, 5UL );
      checkNonZeros( herm, 0UL, 1UL );
      checkNonZeros( herm, 1UL, 1UL );
      checkNonZeros( herm, 2UL, 3UL );

      if( herm(0,0) != cplx( 0,0) || herm(0,1) != cplx(0,0) || herm(0,2) != cplx(-4, 0) ||
          herm(1,0) != cplx( 0,0) || herm(1,1) != cplx(0,0) || herm(1,2) != cplx( 2,-4) ||
          herm(2,0) != cplx(-4,0) || herm(2,1) != cplx(2,4) || herm(2,2) != cplx( 6, 0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( ( 0,0) (0,0) (-4, 0) )\n"
                                     "( ( 0,0) (0,0) ( 2,-4) )\n"
                                     "( (-4,0) (2,4) ( 6, 0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major HermitianMatrix::scale() (complex)";

      // Initialization check
      HT herm( 3UL );
      herm(1,2) = cplx( 1,-2);
      herm(2,0) = cplx(-2, 0);
      herm(2,2) = cplx( 3, 0);

      herm.scale( cplx(3,0) );

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 5UL );
      checkNonZeros( herm, 5UL );
      checkNonZeros( herm, 0UL, 1UL );
      checkNonZeros( herm, 1UL, 1UL );
      checkNonZeros( herm, 2UL, 3UL );

      if( herm(0,0) != cplx( 0,0) || herm(0,1) != cplx(0,0) || herm(0,2) != cplx(-6, 0) ||
          herm(1,0) != cplx( 0,0) || herm(1,1) != cplx(0,0) || herm(1,2) != cplx( 3,-6) ||
          herm(2,0) != cplx(-6,0) || herm(2,1) != cplx(3,6) || herm(2,2) != cplx( 9, 0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( ( 0,0) (0,0) (-6, 0) )\n"
                                     "( ( 0,0) (0,0) ( 3,-6) )\n"
                                     "( (-6,0) (3,6) ( 9, 0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (M*=s)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (M*=s)";

      OHT herm( 3UL );
      herm(1,2) = cplx( 1,-2);
      herm(2,0) = cplx(-2, 0);
      herm(2,2) = cplx( 3, 0);

      herm *= 2;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 5UL );
      checkNonZeros( herm, 5UL );
      checkNonZeros( herm, 0UL, 1UL );
      checkNonZeros( herm, 1UL, 1UL );
      checkNonZeros( herm, 2UL, 3UL );

      if( herm(0,0) != cplx( 0,0) || herm(0,1) != cplx(0,0) || herm(0,2) != cplx(-4, 0) ||
          herm(1,0) != cplx( 0,0) || herm(1,1) != cplx(0,0) || herm(1,2) != cplx( 2,-4) ||
          herm(2,0) != cplx(-4,0) || herm(2,1) != cplx(2,4) || herm(2,2) != cplx( 6, 0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( ( 0,0) (0,0) (-4, 0) )\n"
                                     "( ( 0,0) (0,0) ( 2,-4) )\n"
                                     "( (-4,0) (2,4) ( 6, 0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (M=M*s)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (M=M*s)";

      OHT herm( 3UL );
      herm(1,2) = cplx( 1,-2);
      herm(2,0) = cplx(-2, 0);
      herm(2,2) = cplx( 3, 0);

      herm = herm * 2;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 5UL );
      checkNonZeros( herm, 5UL );
      checkNonZeros( herm, 0UL, 1UL );
      checkNonZeros( herm, 1UL, 1UL );
      checkNonZeros( herm, 2UL, 3UL );

      if( herm(0,0) != cplx( 0,0) || herm(0,1) != cplx(0,0) || herm(0,2) != cplx(-4, 0) ||
          herm(1,0) != cplx( 0,0) || herm(1,1) != cplx(0,0) || herm(1,2) != cplx( 2,-4) ||
          herm(2,0) != cplx(-4,0) || herm(2,1) != cplx(2,4) || herm(2,2) != cplx( 6, 0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( ( 0,0) (0,0) (-4, 0) )\n"
                                     "( ( 0,0) (0,0) ( 2,-4) )\n"
                                     "( (-4,0) (2,4) ( 6, 0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (M=s*M)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (M=s*M)";

      OHT herm( 3UL );
      herm(1,2) = cplx( 1,-2);
      herm(2,0) = cplx(-2, 0);
      herm(2,2) = cplx( 3, 0);

      herm = 2 * herm;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 5UL );
      checkNonZeros( herm, 5UL );
      checkNonZeros( herm, 0UL, 1UL );
      checkNonZeros( herm, 1UL, 1UL );
      checkNonZeros( herm, 2UL, 3UL );

      if( herm(0,0) != cplx( 0,0) || herm(0,1) != cplx(0,0) || herm(0,2) != cplx(-4, 0) ||
          herm(1,0) != cplx( 0,0) || herm(1,1) != cplx(0,0) || herm(1,2) != cplx( 2,-4) ||
          herm(2,0) != cplx(-4,0) || herm(2,1) != cplx(2,4) || herm(2,2) != cplx( 6, 0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( ( 0,0) (0,0) (-4, 0) )\n"
                                     "( ( 0,0) (0,0) ( 2,-4) )\n"
                                     "( (-4,0) (2,4) ( 6, 0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (M/=s)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (M/=s)";

      OHT herm( 3UL );
      herm(1,2) = cplx( 2,-4);
      herm(2,0) = cplx(-4, 0);
      herm(2,2) = cplx( 6, 0);

      herm /= 2;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 5UL );
      checkNonZeros( herm, 5UL );
      checkNonZeros( herm, 0UL, 1UL );
      checkNonZeros( herm, 1UL, 1UL );
      checkNonZeros( herm, 2UL, 3UL );

      if( herm(0,0) != cplx( 0,0) || herm(0,1) != cplx(0,0) || herm(0,2) != cplx(-2, 0) ||
          herm(1,0) != cplx( 0,0) || herm(1,1) != cplx(0,0) || herm(1,2) != cplx( 1,-2) ||
          herm(2,0) != cplx(-2,0) || herm(2,1) != cplx(1,2) || herm(2,2) != cplx( 3, 0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( ( 0,0) (0,0) (-2, 0) )\n"
                                     "( ( 0,0) (0,0) ( 1,-2) )\n"
                                     "( (-2,0) (1,2) ( 3, 0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (M=M/s)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (M=M/s)";

      OHT herm( 3UL );
      herm(1,2) = cplx( 2,-4);
      herm(2,0) = cplx(-4, 0);
      herm(2,2) = cplx( 6, 0);

      herm = herm / 2;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 5UL );
      checkNonZeros( herm, 5UL );
      checkNonZeros( herm, 0UL, 1UL );
      checkNonZeros( herm, 1UL, 1UL );
      checkNonZeros( herm, 2UL, 3UL );

      if( herm(0,0) != cplx( 0,0) || herm(0,1) != cplx(0,0) || herm(0,2) != cplx(-2, 0) ||
          herm(1,0) != cplx( 0,0) || herm(1,1) != cplx(0,0) || herm(1,2) != cplx( 1,-2) ||
          herm(2,0) != cplx(-2,0) || herm(2,1) != cplx(1,2) || herm(2,2) != cplx( 3, 0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( ( 0,0) (0,0) (-2, 0) )\n"
                                     "( ( 0,0) (0,0) ( 1,-2) )\n"
                                     "( (-2,0) (1,2) ( 3, 0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major HermitianMatrix::scale()
   //=====================================================================================

   {
      test_ = "Column-major HermitianMatrix::scale()";

      // Initialization check
      OHT herm( 3UL );
      herm(1,2) = cplx( 1,-2);
      herm(2,0) = cplx(-2, 0);
      herm(2,2) = cplx( 3, 0);

      herm.scale( 2 );

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 5UL );
      checkNonZeros( herm, 5UL );
      checkNonZeros( herm, 0UL, 1UL );
      checkNonZeros( herm, 1UL, 1UL );
      checkNonZeros( herm, 2UL, 3UL );

      if( herm(0,0) != cplx( 0,0) || herm(0,1) != cplx(0,0) || herm(0,2) != cplx(-4, 0) ||
          herm(1,0) != cplx( 0,0) || herm(1,1) != cplx(0,0) || herm(1,2) != cplx( 2,-4) ||
          herm(2,0) != cplx(-4,0) || herm(2,1) != cplx(2,4) || herm(2,2) != cplx( 6, 0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( ( 0,0) (0,0) (-4, 0) )\n"
                                     "( ( 0,0) (0,0) ( 2,-4) )\n"
                                     "( (-4,0) (2,4) ( 6, 0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major HermitianMatrix::scale() (complex)";

      // Initialization check
      OHT herm( 3UL );
      herm(1,2) = cplx( 1,-2);
      herm(2,0) = cplx(-2, 0);
      herm(2,2) = cplx( 3, 0);

      herm.scale( cplx(3,0) );

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 5UL );
      checkNonZeros( herm, 5UL );
      checkNonZeros( herm, 0UL, 1UL );
      checkNonZeros( herm, 1UL, 1UL );
      checkNonZeros( herm, 2UL, 3UL );

      if( herm(0,0) != cplx( 0,0) || herm(0,1) != cplx(0,0) || herm(0,2) != cplx(-6, 0) ||
          herm(1,0) != cplx( 0,0) || herm(1,1) != cplx(0,0) || herm(1,2) != cplx( 3,-6) ||
          herm(2,0) != cplx(-6,0) || herm(2,1) != cplx(3,6) || herm(2,2) != cplx( 9, 0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( ( 0,0) (0,0) (-6, 0) )\n"
                                     "( ( 0,0) (0,0) ( 3,-6) )\n"
                                     "( (-6,0) (3,6) ( 9, 0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the HermitianMatrix function call operator.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of adding and accessing elements via the function call operator
// of the HermitianMatrix specialization. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
void SparseComplexTest::testFunctionCall()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major HermitianMatrix::operator()";

      // Good cases
      {
         HT herm( 3UL );

         // Writing the element (1,1)
         herm(1,1) = cplx(1,0);

         checkRows    ( herm, 3UL );
         checkColumns ( herm, 3UL );
         checkCapacity( herm, 1UL );
         checkNonZeros( herm, 1UL );
         checkNonZeros( herm, 0UL, 0UL );
         checkNonZeros( herm, 1UL, 1UL );
         checkNonZeros( herm, 2UL, 0UL );

         if( herm(0,0) != cplx(0,0) || herm(0,1) != cplx(0,0) || herm(0,2) != cplx(0,0) ||
             herm(1,0) != cplx(0,0) || herm(1,1) != cplx(1,0) || herm(1,2) != cplx(0,0) ||
             herm(2,0) != cplx(0,0) || herm(2,1) != cplx(0,0) || herm(2,2) != cplx(0,0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (0,0) (0,0) (0,0) )\n"
                                        "( (0,0) (1,0) (0,0) )\n"
                                        "( (0,0) (0,0) (0,0) )\n";
            throw std::runtime_error( oss.str() );
         }

         // Writing the elements (2,1) and (1,2)
         herm(2,1) = cplx(2,2);

         checkRows    ( herm, 3UL );
         checkColumns ( herm, 3UL );
         checkCapacity( herm, 3UL );
         checkNonZeros( herm, 3UL );
         checkNonZeros( herm, 0UL, 0UL );
         checkNonZeros( herm, 1UL, 2UL );
         checkNonZeros( herm, 2UL, 1UL );

         if( herm(0,0) != cplx(0,0) || herm(0,1) != cplx(0,0) || herm(0,2) != cplx(0, 0) ||
             herm(1,0) != cplx(0,0) || herm(1,1) != cplx(1,0) || herm(1,2) != cplx(2,-2) ||
             herm(2,0) != cplx(0,0) || herm(2,1) != cplx(2,2) || herm(2,2) != cplx(0, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (0,0) (0,0) (0, 0) )\n"
                                        "( (0,0) (1,0) (2,-2) )\n"
                                        "( (0,0) (2,2) (0, 0) )\n";
            throw std::runtime_error( oss.str() );
         }

         // Writing the elements (0,2) and (2,0)
         herm(0,2) = herm(1,2);

         checkRows    ( herm, 3UL );
         checkColumns ( herm, 3UL );
         checkCapacity( herm, 5UL );
         checkNonZeros( herm, 5UL );
         checkNonZeros( herm, 0UL, 1UL );
         checkNonZeros( herm, 1UL, 2UL );
         checkNonZeros( herm, 2UL, 2UL );

         if( herm(0,0) != cplx(0,0) || herm(0,1) != cplx(0,0) || herm(0,2) != cplx(2,-2) ||
             herm(1,0) != cplx(0,0) || herm(1,1) != cplx(1,0) || herm(1,2) != cplx(2,-2) ||
             herm(2,0) != cplx(2,2) || herm(2,1) != cplx(2,2) || herm(2,2) != cplx(0, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (0,0) (0,0) (2,-2) )\n"
                                        "( (0,0) (1,0) (2,-2) )\n"
                                        "( (2,2) (2,2) (0, 0) )\n";
            throw std::runtime_error( oss.str() );
         }

         // Adding to the elements (1,2) and (2,1)
         herm(1,2) += cplx(3,3);

         checkRows    ( herm, 3UL );
         checkColumns ( herm, 3UL );
         checkCapacity( herm, 5UL );
         checkNonZeros( herm, 5UL );
         checkNonZeros( herm, 0UL, 1UL );
         checkNonZeros( herm, 1UL, 2UL );
         checkNonZeros( herm, 2UL, 2UL );

         if( herm(0,0) != cplx(0,0) || herm(0,1) != cplx(0, 0) || herm(0,2) != cplx(2,-2) ||
             herm(1,0) != cplx(0,0) || herm(1,1) != cplx(1, 0) || herm(1,2) != cplx(5, 1) ||
             herm(2,0) != cplx(2,2) || herm(2,1) != cplx(5,-1) || herm(2,2) != cplx(0, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (0,0) (0, 0) (2,-2) )\n"
                                        "( (0,0) (1, 0) (5, 1) )\n"
                                        "( (2,2) (5,-1) (0, 0) )\n";
            throw std::runtime_error( oss.str() );
         }

         // Subtracting from the elements (0,1) and (1,0)
         herm(0,1) -= cplx(4,4);

         checkRows    ( herm, 3UL );
         checkColumns ( herm, 3UL );
         checkCapacity( herm, 7UL );
         checkNonZeros( herm, 7UL );
         checkNonZeros( herm, 0UL, 2UL );
         checkNonZeros( herm, 1UL, 3UL );
         checkNonZeros( herm, 2UL, 2UL );

         if( herm(0,0) != cplx( 0,0) || herm(0,1) != cplx(-4,-4) || herm(0,2) != cplx(2,-2) ||
             herm(1,0) != cplx(-4,4) || herm(1,1) != cplx( 1, 0) || herm(1,2) != cplx(5, 1) ||
             herm(2,0) != cplx( 2,2) || herm(2,1) != cplx( 5,-1) || herm(2,2) != cplx(0, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( ( 0,0) (-4,-4) (2,-2) )\n"
                                        "( (-4,4) ( 1, 0) (5, 1) )\n"
                                        "( ( 2,2) ( 5,-1) (0, 0) )\n";
            throw std::runtime_error( oss.str() );
         }

         // Multiplying the element (1,1)
         herm(2,0) *= cplx(-3,1);

         checkRows    ( herm, 3UL );
         checkColumns ( herm, 3UL );
         checkCapacity( herm, 7UL );
         checkNonZeros( herm, 7UL );
         checkNonZeros( herm, 0UL, 2UL );
         checkNonZeros( herm, 1UL, 3UL );
         checkNonZeros( herm, 2UL, 2UL );

         if( herm(0,0) != cplx( 0, 0) || herm(0,1) != cplx(-4,-4) || herm(0,2) != cplx(-8,4) ||
             herm(1,0) != cplx(-4, 4) || herm(1,1) != cplx( 1, 0) || herm(1,2) != cplx( 5,1) ||
             herm(2,0) != cplx(-8,-4) || herm(2,1) != cplx( 5,-1) || herm(2,2) != cplx( 0,0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( ( 0, 0) (-4,-4) (-8,4) )\n"
                                        "( (-4, 4) ( 1, 0) ( 5,1) )\n"
                                        "( (-8,-4) ( 5,-1) ( 0,0) )\n";
            throw std::runtime_error( oss.str() );
         }

         // Dividing the elements (0,2) and (2,0)
         herm(1,0) /= 2;

         checkRows    ( herm, 3UL );
         checkColumns ( herm, 3UL );
         checkCapacity( herm, 7UL );
         checkNonZeros( herm, 7UL );
         checkNonZeros( herm, 0UL, 2UL );
         checkNonZeros( herm, 1UL, 3UL );
         checkNonZeros( herm, 2UL, 2UL );

         if( herm(0,0) != cplx( 0, 0) || herm(0,1) != cplx(-2,-2) || herm(0,2) != cplx(-8,4) ||
             herm(1,0) != cplx(-2, 2) || herm(1,1) != cplx( 1, 0) || herm(1,2) != cplx( 5,1) ||
             herm(2,0) != cplx(-8,-4) || herm(2,1) != cplx( 5,-1) || herm(2,2) != cplx( 0,0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( ( 0, 0) (-2,-2) (-8,4) )\n"
                                        "( (-2, 2) ( 1, 0) ( 5,1) )\n"
                                        "( (-8,-4) ( 5,-1) ( 0,0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Failure cases
      {
         HT herm( 3UL );

         // Trying to write the diagonal element (0,0)
         try {
            herm(0,0) = cplx(5,5);

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to add to the diagonal element (1,1)
         try {
            herm(1,1) += cplx(5,5);

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Addition assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to subtract from the diagonal element (2,2)
         try {
            herm(2,2) -= cplx(5,5);

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subtraction assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to multiply the diagonal element (1,1)
         try {
            herm(1,1) *= cplx(5,5);

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Multiplication assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to divide the diagonal element (1,1)
         try {
            herm(1,1) /= cplx(5,5);

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Division assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major HermitianMatrix::operator()";

      // Good cases
      {
         OHT herm( 3UL );

         // Writing the element (1,1)
         herm(1,1) = cplx(1,0);

         checkRows    ( herm, 3UL );
         checkColumns ( herm, 3UL );
         checkCapacity( herm, 1UL );
         checkNonZeros( herm, 1UL );
         checkNonZeros( herm, 0UL, 0UL );
         checkNonZeros( herm, 1UL, 1UL );
         checkNonZeros( herm, 2UL, 0UL );

         if( herm(0,0) != cplx(0,0) || herm(0,1) != cplx(0,0) || herm(0,2) != cplx(0,0) ||
             herm(1,0) != cplx(0,0) || herm(1,1) != cplx(1,0) || herm(1,2) != cplx(0,0) ||
             herm(2,0) != cplx(0,0) || herm(2,1) != cplx(0,0) || herm(2,2) != cplx(0,0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (0,0) (0,0) (0,0) )\n"
                                        "( (0,0) (1,0) (0,0) )\n"
                                        "( (0,0) (0,0) (0,0) )\n";
            throw std::runtime_error( oss.str() );
         }

         // Writing the elements (2,1) and (1,2)
         herm(2,1) = cplx(2,2);

         checkRows    ( herm, 3UL );
         checkColumns ( herm, 3UL );
         checkCapacity( herm, 3UL );
         checkNonZeros( herm, 3UL );
         checkNonZeros( herm, 0UL, 0UL );
         checkNonZeros( herm, 1UL, 2UL );
         checkNonZeros( herm, 2UL, 1UL );

         if( herm(0,0) != cplx(0,0) || herm(0,1) != cplx(0,0) || herm(0,2) != cplx(0, 0) ||
             herm(1,0) != cplx(0,0) || herm(1,1) != cplx(1,0) || herm(1,2) != cplx(2,-2) ||
             herm(2,0) != cplx(0,0) || herm(2,1) != cplx(2,2) || herm(2,2) != cplx(0, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (0,0) (0,0) (0, 0) )\n"
                                        "( (0,0) (1,0) (2,-2) )\n"
                                        "( (0,0) (2,2) (0, 0) )\n";
            throw std::runtime_error( oss.str() );
         }

         // Writing the elements (0,2) and (2,0)
         herm(0,2) = herm(1,2);

         checkRows    ( herm, 3UL );
         checkColumns ( herm, 3UL );
         checkCapacity( herm, 5UL );
         checkNonZeros( herm, 5UL );
         checkNonZeros( herm, 0UL, 1UL );
         checkNonZeros( herm, 1UL, 2UL );
         checkNonZeros( herm, 2UL, 2UL );

         if( herm(0,0) != cplx(0,0) || herm(0,1) != cplx(0,0) || herm(0,2) != cplx(2,-2) ||
             herm(1,0) != cplx(0,0) || herm(1,1) != cplx(1,0) || herm(1,2) != cplx(2,-2) ||
             herm(2,0) != cplx(2,2) || herm(2,1) != cplx(2,2) || herm(2,2) != cplx(0, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (0,0) (0,0) (2,-2) )\n"
                                        "( (0,0) (1,0) (2,-2) )\n"
                                        "( (2,2) (2,2) (0, 0) )\n";
            throw std::runtime_error( oss.str() );
         }

         // Adding to the elements (1,2) and (2,1)
         herm(1,2) += cplx(3,3);

         checkRows    ( herm, 3UL );
         checkColumns ( herm, 3UL );
         checkCapacity( herm, 5UL );
         checkNonZeros( herm, 5UL );
         checkNonZeros( herm, 0UL, 1UL );
         checkNonZeros( herm, 1UL, 2UL );
         checkNonZeros( herm, 2UL, 2UL );

         if( herm(0,0) != cplx(0,0) || herm(0,1) != cplx(0, 0) || herm(0,2) != cplx(2,-2) ||
             herm(1,0) != cplx(0,0) || herm(1,1) != cplx(1, 0) || herm(1,2) != cplx(5, 1) ||
             herm(2,0) != cplx(2,2) || herm(2,1) != cplx(5,-1) || herm(2,2) != cplx(0, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (0,0) (0, 0) (2,-2) )\n"
                                        "( (0,0) (1, 0) (5, 1) )\n"
                                        "( (2,2) (5,-1) (0, 0) )\n";
            throw std::runtime_error( oss.str() );
         }

         // Subtracting from the elements (0,1) and (1,0)
         herm(0,1) -= cplx(4,4);

         checkRows    ( herm, 3UL );
         checkColumns ( herm, 3UL );
         checkCapacity( herm, 7UL );
         checkNonZeros( herm, 7UL );
         checkNonZeros( herm, 0UL, 2UL );
         checkNonZeros( herm, 1UL, 3UL );
         checkNonZeros( herm, 2UL, 2UL );

         if( herm(0,0) != cplx( 0,0) || herm(0,1) != cplx(-4,-4) || herm(0,2) != cplx(2,-2) ||
             herm(1,0) != cplx(-4,4) || herm(1,1) != cplx( 1, 0) || herm(1,2) != cplx(5, 1) ||
             herm(2,0) != cplx( 2,2) || herm(2,1) != cplx( 5,-1) || herm(2,2) != cplx(0, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( ( 0,0) (-4,-4) (2,-2) )\n"
                                        "( (-4,4) ( 1, 0) (5, 1) )\n"
                                        "( ( 2,2) ( 5,-1) (0, 0) )\n";
            throw std::runtime_error( oss.str() );
         }

         // Multiplying the element (1,1)
         herm(2,0) *= cplx(-3,1);

         checkRows    ( herm, 3UL );
         checkColumns ( herm, 3UL );
         checkCapacity( herm, 7UL );
         checkNonZeros( herm, 7UL );
         checkNonZeros( herm, 0UL, 2UL );
         checkNonZeros( herm, 1UL, 3UL );
         checkNonZeros( herm, 2UL, 2UL );

         if( herm(0,0) != cplx( 0, 0) || herm(0,1) != cplx(-4,-4) || herm(0,2) != cplx(-8,4) ||
             herm(1,0) != cplx(-4, 4) || herm(1,1) != cplx( 1, 0) || herm(1,2) != cplx( 5,1) ||
             herm(2,0) != cplx(-8,-4) || herm(2,1) != cplx( 5,-1) || herm(2,2) != cplx( 0,0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( ( 0, 0) (-4,-4) (-8,4) )\n"
                                        "( (-4, 4) ( 1, 0) ( 5,1) )\n"
                                        "( (-8,-4) ( 5,-1) ( 0,0) )\n";
            throw std::runtime_error( oss.str() );
         }

         // Dividing the elements (0,2) and (2,0)
         herm(1,0) /= 2;

         checkRows    ( herm, 3UL );
         checkColumns ( herm, 3UL );
         checkCapacity( herm, 7UL );
         checkNonZeros( herm, 7UL );
         checkNonZeros( herm, 0UL, 2UL );
         checkNonZeros( herm, 1UL, 3UL );
         checkNonZeros( herm, 2UL, 2UL );

         if( herm(0,0) != cplx( 0, 0) || herm(0,1) != cplx(-2,-2) || herm(0,2) != cplx(-8,4) ||
             herm(1,0) != cplx(-2, 2) || herm(1,1) != cplx( 1, 0) || herm(1,2) != cplx( 5,1) ||
             herm(2,0) != cplx(-8,-4) || herm(2,1) != cplx( 5,-1) || herm(2,2) != cplx( 0,0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( ( 0, 0) (-2,-2) (-8,4) )\n"
                                        "( (-2, 2) ( 1, 0) ( 5,1) )\n"
                                        "( (-8,-4) ( 5,-1) ( 0,0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Failure cases
      {
         OHT herm( 3UL );

         // Trying to write the diagonal element (0,0)
         try {
            herm(0,0) = cplx(5,5);

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to add to the diagonal element (1,1)
         try {
            herm(1,1) += cplx(5,5);

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Addition assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to subtract from the diagonal element (2,2)
         try {
            herm(2,2) -= cplx(5,5);

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subtraction assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to multiply the diagonal element (1,1)
         try {
            herm(1,1) *= cplx(5,5);

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Multiplication assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}

         // Trying to divide the diagonal element (1,1)
         try {
            herm(1,1) /= cplx(5,5);

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Division assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the HermitianMatrix iterator implementation.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the iterator implementation of the HermitianMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseComplexTest::testIterator()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      using Iterator      = HT::Iterator;
      using ConstIterator = HT::ConstIterator;

      HT herm( 3UL );
      herm(0,0) = cplx( 4, 0);
      herm(0,1) = cplx( 1,-2);
      herm(1,2) = cplx(-2, 0);
      herm(2,2) = cplx( 3, 0);

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

         ConstIterator it( begin( herm, 1UL ) );

         if( it == end( herm, 1UL ) || it->value() != cplx(1,2) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Failed iterator conversion detected\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Counting the number of elements in 0th row via Iterator (end-begin)
      {
         test_ = "Row-major Iterator subtraction (end-begin)";

         const ptrdiff_t number( end( herm, 0UL ) - begin( herm, 0UL ) );

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

         const ptrdiff_t number( cend( herm, 1UL ) - cbegin( herm, 1UL ) );

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

         ConstIterator it ( cbegin( herm, 2UL ) );
         ConstIterator end( cend( herm, 2UL ) );

         if( it == end || it->value() != cplx(-2,0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid initial iterator detected\n";
            throw std::runtime_error( oss.str() );
         }

         ++it;

         if( it == end || it->value() != cplx(3,0) ) {
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

         Iterator it=begin( herm, 2UL );
         *it = cplx(2,-3);
         ++it;
         *it = cplx(-3,0);

         if( herm(0,0) != cplx(4,0) || herm(0,1) != cplx(1,-2) || herm(0,2) != cplx( 0,0) ||
             herm(1,0) != cplx(1,2) || herm(1,1) != cplx(0, 0) || herm(1,2) != cplx( 2,3) ||
             herm(2,0) != cplx(0,0) || herm(2,1) != cplx(2,-3) || herm(2,2) != cplx(-3,0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (4,0) (1,-2) ( 0,0) )\n"
                                        "( (1,2) (0, 0) ( 2,3) )\n"
                                        "( (0,0) (2,-3) (-3,0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing assignment to diagonal element via Iterator
      {
         test_ = "Row-major assignment to diagonal element via Iterator";

         try {
            const Iterator it = begin( herm, 0UL );
            *it = cplx(5,5);

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Testing addition assignment via Iterator
      {
         test_ = "Row-major addition assignment via Iterator";

         Iterator it=begin( herm, 2UL );
         *it += cplx(2,-3);
         ++it;
         *it += cplx(-3,0);

         if( herm(0,0) != cplx(4,0) || herm(0,1) != cplx(1,-2) || herm(0,2) != cplx( 0,0) ||
             herm(1,0) != cplx(1,2) || herm(1,1) != cplx(0, 0) || herm(1,2) != cplx( 4,6) ||
             herm(2,0) != cplx(0,0) || herm(2,1) != cplx(4,-6) || herm(2,2) != cplx(-6,0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (4,0) (1,-2) ( 0,0) )\n"
                                        "( (1,2) (0, 0) ( 4,6) )\n"
                                        "( (0,0) (4,-6) (-6,0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing addition assignment to diagonal element via Iterator
      {
         test_ = "Row-major addition assignment to diagonal element via Iterator";

         try {
            const Iterator it = begin( herm, 0UL );
            *it += cplx(5,5);

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Testing subtraction assignment via Iterator
      {
         test_ = "Row-major subtraction assignment via Iterator";

         Iterator it=begin( herm, 2UL );
         *it -= cplx(2,-3);
         ++it;
         *it -= cplx(-3,0);

         if( herm(0,0) != cplx(4,0) || herm(0,1) != cplx(1,-2) || herm(0,2) != cplx( 0,0) ||
             herm(1,0) != cplx(1,2) || herm(1,1) != cplx(0, 0) || herm(1,2) != cplx( 2,3) ||
             herm(2,0) != cplx(0,0) || herm(2,1) != cplx(2,-3) || herm(2,2) != cplx(-3,0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (4,0) (1,-2) ( 0,0) )\n"
                                        "( (1,2) (0, 0) ( 2,3) )\n"
                                        "( (0,0) (2,-3) (-3,0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing subtraction assignment to diagonal element via Iterator
      {
         test_ = "Row-major subtraction assignment to diagonal element via Iterator";

         try {
            const Iterator it = begin( herm, 0UL );
            *it -= cplx(5,5);

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Testing multiplication assignment via Iterator
      {
         test_ = "Row-major multiplication assignment via Iterator";

         Iterator it=begin( herm, 2UL );
         *it *= 2;
         ++it;
         *it *= 2;

         if( herm(0,0) != cplx(4,0) || herm(0,1) != cplx(1,-2) || herm(0,2) != cplx( 0,0) ||
             herm(1,0) != cplx(1,2) || herm(1,1) != cplx(0, 0) || herm(1,2) != cplx( 4,6) ||
             herm(2,0) != cplx(0,0) || herm(2,1) != cplx(4,-6) || herm(2,2) != cplx(-6,0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (4,0) (1,-2) ( 0,0) )\n"
                                        "( (1,2) (0, 0) ( 4,6) )\n"
                                        "( (0,0) (4,-6) (-6,0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing multiplication assignment to diagonal element via Iterator
      {
         test_ = "Row-major multiplication assignment to diagonal element via Iterator";

         try {
            const Iterator it = begin( herm, 0UL );
            *it *= cplx(5,5);

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Testing division assignment via Iterator
      {
         test_ = "Row-major division assignment via Iterator";

         Iterator it=begin( herm, 2UL );
         *it /= 2;
         ++it;
         *it /= 2;

         if( herm(0,0) != cplx(4,0) || herm(0,1) != cplx(1,-2) || herm(0,2) != cplx( 0,0) ||
             herm(1,0) != cplx(1,2) || herm(1,1) != cplx(0, 0) || herm(1,2) != cplx( 2,3) ||
             herm(2,0) != cplx(0,0) || herm(2,1) != cplx(2,-3) || herm(2,2) != cplx(-3,0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (4,0) (1,-2) ( 0,0) )\n"
                                        "( (1,2) (0, 0) ( 2,3) )\n"
                                        "( (0,0) (2,-3) (-3,0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing division assignment to diagonal element via Iterator
      {
         test_ = "Row-major division assignment to diagonal element via Iterator";

         try {
            const Iterator it = begin( herm, 0UL );
            *it /= cplx(5,5);

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      using Iterator      = OHT::Iterator;
      using ConstIterator = OHT::ConstIterator;

      OHT herm( 3UL );
      herm(0,0) = cplx( 4, 0);
      herm(0,1) = cplx( 1,-2);
      herm(1,2) = cplx(-2, 0);
      herm(2,2) = cplx( 3, 0);

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

         ConstIterator it( begin( herm, 1UL ) );

         if( it == end( herm, 1UL ) || it->value() != cplx(1,-2) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Failed iterator conversion detected\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Counting the number of elements in 0th row via Iterator (end-begin)
      {
         test_ = "Column-major Iterator subtraction (end-begin)";

         const ptrdiff_t number( end( herm, 0UL ) - begin( herm, 0UL ) );

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
         test_ = "Column-major ConstIterator subtraction (end-begin)";

         const ptrdiff_t number( cend( herm, 1UL ) - cbegin( herm, 1UL ) );

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

         ConstIterator it ( cbegin( herm, 2UL ) );
         ConstIterator end( cend( herm, 2UL ) );

         if( it == end || it->value() != cplx(-2,0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid initial iterator detected\n";
            throw std::runtime_error( oss.str() );
         }

         ++it;

         if( it == end || it->value() != cplx(3,0) ) {
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
         test_ = "Column-major assignment via Iterator";

         Iterator it=begin( herm, 2UL );
         *it = cplx(2,3);
         ++it;
         *it = cplx(-3,0);

         if( herm(0,0) != cplx(4,0) || herm(0,1) != cplx(1,-2) || herm(0,2) != cplx( 0,0) ||
             herm(1,0) != cplx(1,2) || herm(1,1) != cplx(0, 0) || herm(1,2) != cplx( 2,3) ||
             herm(2,0) != cplx(0,0) || herm(2,1) != cplx(2,-3) || herm(2,2) != cplx(-3,0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (4,0) (1,-2) ( 0,0) )\n"
                                        "( (1,2) (0, 0) ( 2,3) )\n"
                                        "( (0,0) (2,-3) (-3,0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing assignment to diagonal element via Iterator
      {
         test_ = "Column-major assignment to diagonal element via Iterator";

         try {
            const Iterator it = begin( herm, 0UL );
            *it = cplx(5,5);

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Testing addition assignment via Iterator
      {
         test_ = "Column-major addition assignment via Iterator";

         Iterator it=begin( herm, 2UL );
         *it += cplx(2,3);
         ++it;
         *it += cplx(-3,0);

         if( herm(0,0) != cplx(4,0) || herm(0,1) != cplx(1,-2) || herm(0,2) != cplx( 0,0) ||
             herm(1,0) != cplx(1,2) || herm(1,1) != cplx(0, 0) || herm(1,2) != cplx( 4,6) ||
             herm(2,0) != cplx(0,0) || herm(2,1) != cplx(4,-6) || herm(2,2) != cplx(-6,0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (4,0) (1,-2) ( 0,0) )\n"
                                        "( (1,2) (0, 0) ( 4,6) )\n"
                                        "( (0,0) (4,-6) (-6,0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing addition assignment to diagonal element via Iterator
      {
         test_ = "Column-major addition assignment to diagonal element via Iterator";

         try {
            const Iterator it = begin( herm, 0UL );
            *it += cplx(5,5);

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Testing subtraction assignment via Iterator
      {
         test_ = "Column-major subtraction assignment via Iterator";

         Iterator it=begin( herm, 2UL );
         *it -= cplx(2,3);
         ++it;
         *it -= cplx(-3,0);

         if( herm(0,0) != cplx(4,0) || herm(0,1) != cplx(1,-2) || herm(0,2) != cplx( 0,0) ||
             herm(1,0) != cplx(1,2) || herm(1,1) != cplx(0, 0) || herm(1,2) != cplx( 2,3) ||
             herm(2,0) != cplx(0,0) || herm(2,1) != cplx(2,-3) || herm(2,2) != cplx(-3,0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (4,0) (1,-2) ( 0,0) )\n"
                                        "( (1,2) (0, 0) ( 2,3) )\n"
                                        "( (0,0) (2,-3) (-3,0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing subtraction assignment to diagonal element via Iterator
      {
         test_ = "Column-major subtraction assignment to diagonal element via Iterator";

         try {
            const Iterator it = begin( herm, 0UL );
            *it -= cplx(5,5);

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Testing multiplication assignment via Iterator
      {
         test_ = "Column-major multiplication assignment via Iterator";

         Iterator it=begin( herm, 2UL );
         *it *= 2;
         ++it;
         *it *= 2;

         if( herm(0,0) != cplx(4,0) || herm(0,1) != cplx(1,-2) || herm(0,2) != cplx( 0,0) ||
             herm(1,0) != cplx(1,2) || herm(1,1) != cplx(0, 0) || herm(1,2) != cplx( 4,6) ||
             herm(2,0) != cplx(0,0) || herm(2,1) != cplx(4,-6) || herm(2,2) != cplx(-6,0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (4,0) (1,-2) ( 0,0) )\n"
                                        "( (1,2) (0, 0) ( 4,6) )\n"
                                        "( (0,0) (4,-6) (-6,0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing multiplication assignment to diagonal element via Iterator
      {
         test_ = "Column-major multiplication assignment to diagonal element via Iterator";

         try {
            const Iterator it = begin( herm, 0UL );
            *it *= cplx(5,5);

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Testing division assignment via Iterator
      {
         test_ = "Column-major division assignment via Iterator";

         Iterator it=begin( herm, 2UL );
         *it /= 2;
         ++it;
         *it /= 2;

         if( herm(0,0) != cplx(4,0) || herm(0,1) != cplx(1,-2) || herm(0,2) != cplx( 0,0) ||
             herm(1,0) != cplx(1,2) || herm(1,1) != cplx(0, 0) || herm(1,2) != cplx( 2,3) ||
             herm(2,0) != cplx(0,0) || herm(2,1) != cplx(2,-3) || herm(2,2) != cplx(-3,0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (4,0) (1,-2) ( 0,0) )\n"
                                        "( (1,2) (0, 0) ( 2,3) )\n"
                                        "( (0,0) (2,-3) (-3,0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Testing division assignment to diagonal element via Iterator
      {
         test_ = "Column-major division assignment to diagonal element via Iterator";

         try {
            const Iterator it = begin( herm, 0UL );
            *it /= cplx(5,5);

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to diagonal matrix element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c nonZeros() member function of the HermitianMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c nonZeros() member function of the HermitianMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseComplexTest::testNonZeros()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major HermitianMatrix::nonZeros()";

      // Empty matrix
      {
         HT herm( 3UL );

         checkRows    ( herm, 3UL );
         checkColumns ( herm, 3UL );
         checkNonZeros( herm, 0UL );
         checkNonZeros( herm, 0UL, 0UL );
         checkNonZeros( herm, 1UL, 0UL );
         checkNonZeros( herm, 2UL, 0UL );

         if( herm(0,0) != cplx(0,0) || herm(0,1) != cplx(0,0) || herm(0,2) != cplx(0,0) ||
             herm(1,0) != cplx(0,0) || herm(1,1) != cplx(0,0) || herm(1,2) != cplx(0,0) ||
             herm(2,0) != cplx(0,0) || herm(2,1) != cplx(0,0) || herm(2,2) != cplx(0,0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (0,0) (0,0) (0,0) )\n"
                                        "( (0,0) (0,0) (0,0) )\n"
                                        "( (0,0) (0,0) (0,0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Partially filled matrix
      {
         HT herm( 3UL );
         herm(0,0) = cplx( 1, 0);
         herm(1,2) = cplx(-2,-3);
         herm(2,0) = cplx( 0, 0);
         herm(2,2) = cplx( 3, 0);

         checkRows    ( herm, 3UL );
         checkColumns ( herm, 3UL );
         checkCapacity( herm, 4UL );
         checkNonZeros( herm, 4UL );
         checkNonZeros( herm, 0UL, 1UL );
         checkNonZeros( herm, 1UL, 1UL );
         checkNonZeros( herm, 2UL, 2UL );

         if( herm(0,0) != cplx(1,0) || herm(0,1) != cplx( 0,0) || herm(0,2) != cplx( 0, 0) ||
             herm(1,0) != cplx(0,0) || herm(1,1) != cplx( 0,0) || herm(1,2) != cplx(-2,-3) ||
             herm(2,0) != cplx(0,0) || herm(2,1) != cplx(-2,3) || herm(2,2) != cplx( 3, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (1,0) ( 0,0) ( 0, 0) )\n"
                                        "( (0,0) ( 0,0) (-2,-3) )\n"
                                        "( (0,0) (-2,3) ( 3, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Fully filled matrix
      {
         HT herm( 3UL );
         herm(0,0) = cplx(-1, 0);
         herm(0,1) = cplx( 2, 1);
         herm(0,2) = cplx(-3,-2);
         herm(1,1) = cplx( 4, 0);
         herm(1,2) = cplx(-5,-1);
         herm(2,2) = cplx( 6, 0);

         checkRows    ( herm, 3UL );
         checkColumns ( herm, 3UL );
         checkCapacity( herm, 9UL );
         checkNonZeros( herm, 9UL );
         checkNonZeros( herm, 0UL, 3UL );
         checkNonZeros( herm, 1UL, 3UL );
         checkNonZeros( herm, 2UL, 3UL );

         if( herm(0,0) != cplx(-1, 0) || herm(0,1) != cplx( 2,1) || herm(0,2) != cplx(-3,-2) ||
             herm(1,0) != cplx( 2,-1) || herm(1,1) != cplx( 4,0) || herm(1,2) != cplx(-5,-1) ||
             herm(2,0) != cplx(-3, 2) || herm(2,1) != cplx(-5,1) || herm(2,2) != cplx( 6, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (-1, 0) ( 2,1) (-3,-2) )\n"
                                        "( ( 2,-1) ( 4,0) (-5,-1) )\n"
                                        "( (-3, 2) (-5,1) ( 6, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major HermitianMatrix::nonZeros()";

      // Empty matrix
      {
         OHT herm( 3UL );

         checkRows    ( herm, 3UL );
         checkColumns ( herm, 3UL );
         checkNonZeros( herm, 0UL );
         checkNonZeros( herm, 0UL, 0UL );
         checkNonZeros( herm, 1UL, 0UL );
         checkNonZeros( herm, 2UL, 0UL );

         if( herm(0,0) != cplx(0,0) || herm(0,1) != cplx(0,0) || herm(0,2) != cplx(0,0) ||
             herm(1,0) != cplx(0,0) || herm(1,1) != cplx(0,0) || herm(1,2) != cplx(0,0) ||
             herm(2,0) != cplx(0,0) || herm(2,1) != cplx(0,0) || herm(2,2) != cplx(0,0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (0,0) (0,0) (0,0) )\n"
                                        "( (0,0) (0,0) (0,0) )\n"
                                        "( (0,0) (0,0) (0,0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Partially filled matrix
      {
         HT herm( 3UL );
         herm(0,0) = cplx( 1, 0);
         herm(1,2) = cplx(-2,-3);
         herm(2,0) = cplx( 0, 0);
         herm(2,2) = cplx( 3, 0);

         checkRows    ( herm, 3UL );
         checkColumns ( herm, 3UL );
         checkCapacity( herm, 4UL );
         checkNonZeros( herm, 4UL );
         checkNonZeros( herm, 0UL, 1UL );
         checkNonZeros( herm, 1UL, 1UL );
         checkNonZeros( herm, 2UL, 2UL );

         if( herm(0,0) != cplx(1,0) || herm(0,1) != cplx( 0,0) || herm(0,2) != cplx( 0, 0) ||
             herm(1,0) != cplx(0,0) || herm(1,1) != cplx( 0,0) || herm(1,2) != cplx(-2,-3) ||
             herm(2,0) != cplx(0,0) || herm(2,1) != cplx(-2,3) || herm(2,2) != cplx( 3, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (1,0) ( 0,0) ( 0, 0) )\n"
                                        "( (0,0) ( 0,0) (-2,-3) )\n"
                                        "( (0,0) (-2,3) ( 3, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Fully filled matrix
      {
         HT herm( 3UL );
         herm(0,0) = cplx(-1, 0);
         herm(0,1) = cplx( 2, 1);
         herm(0,2) = cplx(-3,-2);
         herm(1,1) = cplx( 4, 0);
         herm(1,2) = cplx(-5,-1);
         herm(2,2) = cplx( 6, 0);

         checkRows    ( herm, 3UL );
         checkColumns ( herm, 3UL );
         checkCapacity( herm, 9UL );
         checkNonZeros( herm, 9UL );
         checkNonZeros( herm, 0UL, 3UL );
         checkNonZeros( herm, 1UL, 3UL );
         checkNonZeros( herm, 2UL, 3UL );

         if( herm(0,0) != cplx(-1, 0) || herm(0,1) != cplx( 2,1) || herm(0,2) != cplx(-3,-2) ||
             herm(1,0) != cplx( 2,-1) || herm(1,1) != cplx( 4,0) || herm(1,2) != cplx(-5,-1) ||
             herm(2,0) != cplx(-3, 2) || herm(2,1) != cplx(-5,1) || herm(2,2) != cplx( 6, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (-1, 0) ( 2,1) (-3,-2) )\n"
                                        "( ( 2,-1) ( 4,0) (-5,-1) )\n"
                                        "( (-3, 2) (-5,1) ( 6, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c reset() member function of the HermitianMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c reset() member function of the HermitianMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseComplexTest::testReset()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major HermitianMatrix::reset()";

      // Initialization check
      HT herm( 3UL );
      herm(0,0) = cplx(-1, 0);
      herm(0,1) = cplx( 2, 1);
      herm(0,2) = cplx(-3,-2);
      herm(1,1) = cplx( 4, 0);
      herm(1,2) = cplx(-5,-1);
      herm(2,2) = cplx( 6, 0);

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 9UL );
      checkNonZeros( herm, 9UL );
      checkNonZeros( herm, 0UL, 3UL );
      checkNonZeros( herm, 1UL, 3UL );
      checkNonZeros( herm, 2UL, 3UL );

      if( herm(0,0) != cplx(-1, 0) || herm(0,1) != cplx( 2, 1) || herm(0,2) != cplx(-3,-2) ||
          herm(1,0) != cplx( 2,-1) || herm(1,1) != cplx( 4, 0) || herm(1,2) != cplx(-5,-1) ||
          herm(2,0) != cplx(-3, 2) || herm(2,1) != cplx(-5, 1) || herm(2,2) != cplx( 6, 0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (-1, 0) ( 2,1) (-3,-2) )\n"
                                     "( ( 2,-1) ( 4,0) (-5,-1) )\n"
                                     "( (-3, 2) (-5,1) ( 6, 0) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting a single element
      reset( herm(0,1) );

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 9UL );
      checkNonZeros( herm, 7UL );
      checkNonZeros( herm, 0UL, 2UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 3UL );

      if( herm(0,0) != cplx(-1,0) || herm(0,1) != cplx( 0,0) || herm(0,2) != cplx(-3,-2) ||
          herm(1,0) != cplx( 0,0) || herm(1,1) != cplx( 4,0) || herm(1,2) != cplx(-5,-1) ||
          herm(2,0) != cplx(-3,2) || herm(2,1) != cplx(-5,1) || herm(2,2) != cplx( 6, 0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (-1,0) ( 0,0) (-3,-2) )\n"
                                     "( ( 0,0) ( 4,0) (-5,-1) )\n"
                                     "( (-3,2) (-5,1) ( 6, 0) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting row 1
      reset( herm, 1UL );

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 9UL );
      checkNonZeros( herm, 4UL );
      checkNonZeros( herm, 0UL, 2UL );
      checkNonZeros( herm, 1UL, 0UL );
      checkNonZeros( herm, 2UL, 2UL );

      if( herm(0,0) != cplx(-1,0) || herm(0,1) != cplx(0,0) || herm(0,2) != cplx(-3,-2) ||
          herm(1,0) != cplx( 0,0) || herm(1,1) != cplx(0,0) || herm(1,2) != cplx( 0, 0) ||
          herm(2,0) != cplx(-3,2) || herm(2,1) != cplx(0,0) || herm(2,2) != cplx( 6, 0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (-1,0) (0,0) (-3,-2) )\n"
                                     "( ( 0,0) (0,0) ( 0, 0) )\n"
                                     "( (-3,2) (0,0) ( 6, 0) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting the entire matrix
      reset( herm );

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 9UL );
      checkNonZeros( herm, 0UL );
      checkNonZeros( herm, 0UL, 0UL );
      checkNonZeros( herm, 1UL, 0UL );
      checkNonZeros( herm, 2UL, 0UL );

      if( herm(0,0) != cplx(0,0) || herm(0,1) != cplx(0,0) || herm(0,2) != cplx(0,0) ||
          herm(1,0) != cplx(0,0) || herm(1,1) != cplx(0,0) || herm(1,2) != cplx(0,0) ||
          herm(2,0) != cplx(0,0) || herm(2,1) != cplx(0,0) || herm(2,2) != cplx(0,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (0,0) (0,0) (0,0) )\n"
                                     "( (0,0) (0,0) (0,0) )\n"
                                     "( (0,0) (0,0) (0,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major HermitianMatrix::reset()";

      // Initialization check
      OHT herm( 3UL );
      herm(0,0) = cplx(-1, 0);
      herm(0,1) = cplx( 2, 1);
      herm(0,2) = cplx(-3,-2);
      herm(1,1) = cplx( 4, 0);
      herm(1,2) = cplx(-5,-1);
      herm(2,2) = cplx( 6, 0);

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 9UL );
      checkNonZeros( herm, 9UL );
      checkNonZeros( herm, 0UL, 3UL );
      checkNonZeros( herm, 1UL, 3UL );
      checkNonZeros( herm, 2UL, 3UL );

      if( herm(0,0) != cplx(-1, 0) || herm(0,1) != cplx( 2, 1) || herm(0,2) != cplx(-3,-2) ||
          herm(1,0) != cplx( 2,-1) || herm(1,1) != cplx( 4, 0) || herm(1,2) != cplx(-5,-1) ||
          herm(2,0) != cplx(-3, 2) || herm(2,1) != cplx(-5, 1) || herm(2,2) != cplx( 6, 0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (-1, 0) ( 2,1) (-3,-2) )\n"
                                     "( ( 2,-1) ( 4,0) (-5,-1) )\n"
                                     "( (-3, 2) (-5,1) ( 6, 0) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting a single element
      reset( herm(0,1) );

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 9UL );
      checkNonZeros( herm, 7UL );
      checkNonZeros( herm, 0UL, 2UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 3UL );

      if( herm(0,0) != cplx(-1,0) || herm(0,1) != cplx( 0,0) || herm(0,2) != cplx(-3,-2) ||
          herm(1,0) != cplx( 0,0) || herm(1,1) != cplx( 4,0) || herm(1,2) != cplx(-5,-1) ||
          herm(2,0) != cplx(-3,2) || herm(2,1) != cplx(-5,1) || herm(2,2) != cplx( 6, 0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (-1,0) ( 0,0) (-3,-2) )\n"
                                     "( ( 0,0) ( 4,0) (-5,-1) )\n"
                                     "( (-3,2) (-5,1) ( 6, 0) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting row 1
      reset( herm, 1UL );

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 9UL );
      checkNonZeros( herm, 4UL );
      checkNonZeros( herm, 0UL, 2UL );
      checkNonZeros( herm, 1UL, 0UL );
      checkNonZeros( herm, 2UL, 2UL );

      if( herm(0,0) != cplx(-1,0) || herm(0,1) != cplx(0,0) || herm(0,2) != cplx(-3,-2) ||
          herm(1,0) != cplx( 0,0) || herm(1,1) != cplx(0,0) || herm(1,2) != cplx( 0, 0) ||
          herm(2,0) != cplx(-3,2) || herm(2,1) != cplx(0,0) || herm(2,2) != cplx( 6, 0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (-1,0) (0,0) (-3,-2) )\n"
                                     "( ( 0,0) (0,0) ( 0, 0) )\n"
                                     "( (-3,2) (0,0) ( 6, 0) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting the entire matrix
      reset( herm );

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 9UL );
      checkNonZeros( herm, 0UL );
      checkNonZeros( herm, 0UL, 0UL );
      checkNonZeros( herm, 1UL, 0UL );
      checkNonZeros( herm, 2UL, 0UL );

      if( herm(0,0) != cplx(0,0) || herm(0,1) != cplx(0,0) || herm(0,2) != cplx(0,0) ||
          herm(1,0) != cplx(0,0) || herm(1,1) != cplx(0,0) || herm(1,2) != cplx(0,0) ||
          herm(2,0) != cplx(0,0) || herm(2,1) != cplx(0,0) || herm(2,2) != cplx(0,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (0,0) (0,0) (0,0) )\n"
                                     "( (0,0) (0,0) (0,0) )\n"
                                     "( (0,0) (0,0) (0,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c clear() member function of the HermitianMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c clear() member function of the HermitianMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseComplexTest::testClear()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major HermitianMatrix::clear()";

      // Initialization check
      HT herm( 3UL );
      herm(0,0) = cplx(-1, 0);
      herm(0,1) = cplx( 2, 1);
      herm(0,2) = cplx(-3,-2);
      herm(1,1) = cplx( 4, 0);
      herm(1,2) = cplx(-5,-1);
      herm(2,2) = cplx( 6, 0);

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 9UL );
      checkNonZeros( herm, 9UL );
      checkNonZeros( herm, 0UL, 3UL );
      checkNonZeros( herm, 1UL, 3UL );
      checkNonZeros( herm, 2UL, 3UL );

      if( herm(0,0) != cplx(-1, 0) || herm(0,1) != cplx( 2, 1) || herm(0,2) != cplx(-3,-2) ||
          herm(1,0) != cplx( 2,-1) || herm(1,1) != cplx( 4, 0) || herm(1,2) != cplx(-5,-1) ||
          herm(2,0) != cplx(-3, 2) || herm(2,1) != cplx(-5, 1) || herm(2,2) != cplx( 6, 0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (-1, 0) ( 2,1) (-3,-2) )\n"
                                     "( ( 2,-1) ( 4,0) (-5,-1) )\n"
                                     "( (-3, 2) (-5,1) ( 6, 0) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing a single element
      clear( herm(0,1) );

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 9UL );
      checkNonZeros( herm, 7UL );
      checkNonZeros( herm, 0UL, 2UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 3UL );

      if( herm(0,0) != cplx(-1,0) || herm(0,1) != cplx( 0,0) || herm(0,2) != cplx(-3,-2) ||
          herm(1,0) != cplx( 0,0) || herm(1,1) != cplx( 4,0) || herm(1,2) != cplx(-5,-1) ||
          herm(2,0) != cplx(-3,2) || herm(2,1) != cplx(-5,1) || herm(2,2) != cplx( 6, 0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (-1,0) ( 0,0) (-3,-2) )\n"
                                     "( ( 0,0) ( 4,0) (-5,-1) )\n"
                                     "( (-3,2) (-5,1) ( 6, 0) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing the matrix
      clear( herm );

      checkRows    ( herm, 0UL );
      checkColumns ( herm, 0UL );
      checkNonZeros( herm, 0UL );
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major HermitianMatrix::clear()";

      // Initialization check
      OHT herm( 3UL );
      herm(0,0) = cplx(-1, 0);
      herm(0,1) = cplx( 2, 1);
      herm(0,2) = cplx(-3,-2);
      herm(1,1) = cplx( 4, 0);
      herm(1,2) = cplx(-5,-1);
      herm(2,2) = cplx( 6, 0);

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 9UL );
      checkNonZeros( herm, 9UL );
      checkNonZeros( herm, 0UL, 3UL );
      checkNonZeros( herm, 1UL, 3UL );
      checkNonZeros( herm, 2UL, 3UL );

      if( herm(0,0) != cplx(-1, 0) || herm(0,1) != cplx( 2, 1) || herm(0,2) != cplx(-3,-2) ||
          herm(1,0) != cplx( 2,-1) || herm(1,1) != cplx( 4, 0) || herm(1,2) != cplx(-5,-1) ||
          herm(2,0) != cplx(-3, 2) || herm(2,1) != cplx(-5, 1) || herm(2,2) != cplx( 6, 0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (-1, 0) ( 2,1) (-3,-2) )\n"
                                     "( ( 2,-1) ( 4,0) (-5,-1) )\n"
                                     "( (-3, 2) (-5,1) ( 6, 0) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing a single element
      clear( herm(0,1) );

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 9UL );
      checkNonZeros( herm, 7UL );
      checkNonZeros( herm, 0UL, 2UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 3UL );

      if( herm(0,0) != cplx(-1,0) || herm(0,1) != cplx( 0,0) || herm(0,2) != cplx(-3,-2) ||
          herm(1,0) != cplx( 0,0) || herm(1,1) != cplx( 4,0) || herm(1,2) != cplx(-5,-1) ||
          herm(2,0) != cplx(-3,2) || herm(2,1) != cplx(-5,1) || herm(2,2) != cplx( 6, 0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (-1,0) ( 0,0) (-3,-2) )\n"
                                     "( ( 0,0) ( 4,0) (-5,-1) )\n"
                                     "( (-3,2) (-5,1) ( 6, 0) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing the matrix
      clear( herm );

      checkRows    ( herm, 0UL );
      checkColumns ( herm, 0UL );
      checkNonZeros( herm, 0UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c resize() member function of the HermitianMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c resize() member function of the HermitianMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseComplexTest::testResize()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major HermitianMatrix::resize()";

      // Initialization check
      HT herm;

      checkRows    ( herm, 0UL );
      checkColumns ( herm, 0UL );
      checkNonZeros( herm, 0UL );

      // Resizing to 2x2
      herm.resize( 2UL );

      checkRows    ( herm, 2UL );
      checkColumns ( herm, 2UL );
      checkNonZeros( herm, 0UL );
      checkNonZeros( herm, 0UL, 0UL );
      checkNonZeros( herm, 1UL, 0UL );

      if( herm(0,0) != cplx(0,0) || herm(0,1) != cplx(0,0) ||
          herm(1,0) != cplx(0,0) || herm(1,1) != cplx(0,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (0,0) (0,0) )\n"
                                     "( (0,0) (0,0) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 4x4 and preserving the elements
      herm(0,1) = cplx(1,-1);
      herm(1,1) = cplx(2, 0);
      herm.resize( 4UL, true );

      checkRows    ( herm, 4UL );
      checkColumns ( herm, 4UL );
      checkCapacity( herm, 3UL );
      checkNonZeros( herm, 3UL );
      checkNonZeros( herm, 0UL, 1UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 0UL );
      checkNonZeros( herm, 3UL, 0UL );

      if( herm(0,0) != cplx(0,0) || herm(0,1) != cplx(1,-1) || herm(0,2) != cplx(0,0) || herm(0,3) != cplx(0,0) ||
          herm(1,0) != cplx(1,1) || herm(1,1) != cplx(2, 0) || herm(1,2) != cplx(0,0) || herm(1,3) != cplx(0,0) ||
          herm(2,0) != cplx(0,0) || herm(2,1) != cplx(0, 0) || herm(2,2) != cplx(0,0) || herm(2,3) != cplx(0,0) ||
          herm(3,0) != cplx(0,0) || herm(3,1) != cplx(0, 0) || herm(3,2) != cplx(0,0) || herm(3,3) != cplx(0,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (0,0) (1,-1) (0,0) (0,0) )\n"
                                     "( (1,1) (2, 0) (0,0) (0,0) )\n"
                                     "( (0,0) (0, 0) (0,0) (0,0) )\n"
                                     "( (0,0) (0, 0) (0,0) (0,0) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 2x2
      herm(2,2) = cplx(3,0);
      herm.resize( 2UL );

      checkRows    ( herm, 2UL );
      checkColumns ( herm, 2UL );
      checkCapacity( herm, 3UL );
      checkNonZeros( herm, 3UL );
      checkNonZeros( herm, 0UL, 1UL );
      checkNonZeros( herm, 1UL, 2UL );

      if( herm(0,0) != cplx(0,0) || herm(0,1) != cplx(1,-1) ||
          herm(1,0) != cplx(1,1) || herm(1,1) != cplx(2,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (0,0) (1,-1) )\n"
                                     "( (1,1) (2,0) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 0x0
      herm.resize( 0UL );

      checkRows    ( herm, 0UL );
      checkColumns ( herm, 0UL );
      checkNonZeros( herm, 0UL );
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major HermitianMatrix::resize()";

      // Initialization check
      OHT herm;

      checkRows    ( herm, 0UL );
      checkColumns ( herm, 0UL );
      checkNonZeros( herm, 0UL );

      // Resizing to 2x2
      herm.resize( 2UL );

      checkRows    ( herm, 2UL );
      checkColumns ( herm, 2UL );
      checkNonZeros( herm, 0UL );
      checkNonZeros( herm, 0UL, 0UL );
      checkNonZeros( herm, 1UL, 0UL );

      if( herm(0,0) != cplx(0,0) || herm(0,1) != cplx(0,0) ||
          herm(1,0) != cplx(0,0) || herm(1,1) != cplx(0,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (0,0) (0,0) )\n"
                                     "( (0,0) (0,0) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 4x4 and preserving the elements
      herm(0,1) = cplx(1,-1);
      herm(1,1) = cplx(2, 0);
      herm.resize( 4UL, true );

      checkRows    ( herm, 4UL );
      checkColumns ( herm, 4UL );
      checkCapacity( herm, 3UL );
      checkNonZeros( herm, 3UL );
      checkNonZeros( herm, 0UL, 1UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 0UL );
      checkNonZeros( herm, 3UL, 0UL );

      if( herm(0,0) != cplx(0,0) || herm(0,1) != cplx(1,-1) || herm(0,2) != cplx(0,0) || herm(0,3) != cplx(0,0) ||
          herm(1,0) != cplx(1,1) || herm(1,1) != cplx(2, 0) || herm(1,2) != cplx(0,0) || herm(1,3) != cplx(0,0) ||
          herm(2,0) != cplx(0,0) || herm(2,1) != cplx(0, 0) || herm(2,2) != cplx(0,0) || herm(2,3) != cplx(0,0) ||
          herm(3,0) != cplx(0,0) || herm(3,1) != cplx(0, 0) || herm(3,2) != cplx(0,0) || herm(3,3) != cplx(0,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (0,0) (1,-1) (0,0) (0,0) )\n"
                                     "( (1,1) (2, 0) (0,0) (0,0) )\n"
                                     "( (0,0) (0, 0) (0,0) (0,0) )\n"
                                     "( (0,0) (0, 0) (0,0) (0,0) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 2x2
      herm(2,2) = cplx(3,0);
      herm.resize( 2UL );

      checkRows    ( herm, 2UL );
      checkColumns ( herm, 2UL );
      checkCapacity( herm, 3UL );
      checkNonZeros( herm, 3UL );
      checkNonZeros( herm, 0UL, 1UL );
      checkNonZeros( herm, 1UL, 2UL );

      if( herm(0,0) != cplx(0,0) || herm(0,1) != cplx(1,-1) ||
          herm(1,0) != cplx(1,1) || herm(1,1) != cplx(2,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (0,0) (1,-1) )\n"
                                     "( (1,1) (2,0) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 0x0
      herm.resize( 0UL );

      checkRows    ( herm, 0UL );
      checkColumns ( herm, 0UL );
      checkNonZeros( herm, 0UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c reserve() member function of the HermitianMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c reserve() member function of the HermitianMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseComplexTest::testReserve()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major HermitianMatrix::reserve()";

      // Initialization check
      HT herm;

      checkRows    ( herm, 0UL );
      checkColumns ( herm, 0UL );
      checkNonZeros( herm, 0UL );

      // Increasing the capacity of the matrix
      herm.reserve( 10UL );

      checkRows    ( herm,  0UL );
      checkColumns ( herm,  0UL );
      checkCapacity( herm, 10UL );
      checkNonZeros( herm,  0UL );

      // Further increasing the capacity of the matrix
      herm.reserve( 20UL );

      checkRows    ( herm,  0UL );
      checkColumns ( herm,  0UL );
      checkCapacity( herm, 20UL );
      checkNonZeros( herm,  0UL );
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major HermitianMatrix::reserve()";

      // Initialization check
      OHT herm;

      checkRows    ( herm, 0UL );
      checkColumns ( herm, 0UL );
      checkNonZeros( herm, 0UL );

      // Increasing the capacity of the matrix
      herm.reserve( 10UL );

      checkRows    ( herm,  0UL );
      checkColumns ( herm,  0UL );
      checkCapacity( herm, 10UL );
      checkNonZeros( herm,  0UL );

      // Further increasing the capacity of the matrix
      herm.reserve( 20UL );

      checkRows    ( herm,  0UL );
      checkColumns ( herm,  0UL );
      checkCapacity( herm, 20UL );
      checkNonZeros( herm,  0UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c trim() member function of the HermitianMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c trim() member function of the HermitianMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseComplexTest::testTrim()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major HermitianMatrix::trim()";

      // Initialization check
      HT herm( 3UL );

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkNonZeros( herm, 0UL );

      // Increasing the row capacity of the matrix
      herm.reserve( 0UL, 10UL );
      herm.reserve( 1UL, 15UL );
      herm.reserve( 2UL, 20UL );

      checkRows    ( herm,  3UL );
      checkColumns ( herm,  3UL );
      checkCapacity( herm, 45UL );
      checkCapacity( herm,  0UL, 10UL );
      checkCapacity( herm,  1UL, 15UL );
      checkCapacity( herm,  2UL, 20UL );

      // Trimming the matrix
      herm.trim();

      checkRows    ( herm,  3UL );
      checkColumns ( herm,  3UL );
      checkCapacity( herm, 45UL );
      checkCapacity( herm,  0UL, 0UL );
      checkCapacity( herm,  1UL, 0UL );
      checkCapacity( herm,  2UL, 0UL );
   }

   {
      test_ = "Row-major HermitianMatrix::trim( size_t )";

      // Initialization check
      HT herm( 3UL, 3UL );

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkNonZeros( herm, 0UL );

      // Increasing the row capacity of the matrix
      herm.reserve( 0UL, 10UL );
      herm.reserve( 1UL, 15UL );
      herm.reserve( 2UL, 20UL );

      checkRows    ( herm,  3UL );
      checkColumns ( herm,  3UL );
      checkCapacity( herm, 45UL );
      checkCapacity( herm,  0UL, 10UL );
      checkCapacity( herm,  1UL, 15UL );
      checkCapacity( herm,  2UL, 20UL );

      // Trimming the 0th row
      herm.trim( 0UL );

      checkRows    ( herm,  3UL );
      checkColumns ( herm,  3UL );
      checkCapacity( herm, 45UL );
      checkCapacity( herm,  0UL,  0UL );
      checkCapacity( herm,  1UL, 25UL );
      checkCapacity( herm,  2UL, 20UL );

      // Trimming the 1st row
      herm.trim( 1UL );

      checkRows    ( herm,  3UL );
      checkColumns ( herm,  3UL );
      checkCapacity( herm, 45UL );
      checkCapacity( herm,  0UL,  0UL );
      checkCapacity( herm,  1UL,  0UL );
      checkCapacity( herm,  2UL, 45UL );

      // Trimming the 2nd row
      herm.trim( 2UL );

      checkRows    ( herm,  3UL );
      checkColumns ( herm,  3UL );
      checkCapacity( herm, 45UL );
      checkCapacity( herm,  0UL, 0UL );
      checkCapacity( herm,  1UL, 0UL );
      checkCapacity( herm,  2UL, 0UL );
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major HermitianMatrix::trim()";

      // Initialization check
      OHT herm( 3UL );

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkNonZeros( herm, 0UL );

      // Increasing the row capacity of the matrix
      herm.reserve( 0UL, 10UL );
      herm.reserve( 1UL, 15UL );
      herm.reserve( 2UL, 20UL );

      checkRows    ( herm,  3UL );
      checkColumns ( herm,  3UL );
      checkCapacity( herm, 45UL );
      checkCapacity( herm,  0UL, 10UL );
      checkCapacity( herm,  1UL, 15UL );
      checkCapacity( herm,  2UL, 20UL );

      // Trimming the matrix
      herm.trim();

      checkRows    ( herm,  3UL );
      checkColumns ( herm,  3UL );
      checkCapacity( herm, 45UL );
      checkCapacity( herm,  0UL, 0UL );
      checkCapacity( herm,  1UL, 0UL );
      checkCapacity( herm,  2UL, 0UL );
   }

   {
      test_ = "Column-major HermitianMatrix::trim( size_t )";

      // Initialization check
      OHT herm( 3UL, 3UL );

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkNonZeros( herm, 0UL );

      // Increasing the column capacity of the matrix
      herm.reserve( 0UL, 10UL );
      herm.reserve( 1UL, 15UL );
      herm.reserve( 2UL, 20UL );

      checkRows    ( herm,  3UL );
      checkColumns ( herm,  3UL );
      checkCapacity( herm, 45UL );
      checkCapacity( herm,  0UL, 10UL );
      checkCapacity( herm,  1UL, 15UL );
      checkCapacity( herm,  2UL, 20UL );

      // Trimming the 0th column
      herm.trim( 0UL );

      checkRows    ( herm,  3UL );
      checkColumns ( herm,  3UL );
      checkCapacity( herm, 45UL );
      checkCapacity( herm,  0UL,  0UL );
      checkCapacity( herm,  1UL, 25UL );
      checkCapacity( herm,  2UL, 20UL );

      // Trimming the 1st column
      herm.trim( 1UL );

      checkRows    ( herm,  3UL );
      checkColumns ( herm,  3UL );
      checkCapacity( herm, 45UL );
      checkCapacity( herm,  0UL,  0UL );
      checkCapacity( herm,  1UL,  0UL );
      checkCapacity( herm,  2UL, 45UL );

      // Trimming the 2nd column
      herm.trim( 2UL );

      checkRows    ( herm,  3UL );
      checkColumns ( herm,  3UL );
      checkCapacity( herm, 45UL );
      checkCapacity( herm,  0UL, 0UL );
      checkCapacity( herm,  1UL, 0UL );
      checkCapacity( herm,  2UL, 0UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c shrinkToFit() member function of the HermitianMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c shrinkToFit() member function of the HermitianMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseComplexTest::testShrinkToFit()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major HermitianMatrix::shrinkToFit()";

      // Shrinking a matrix without excessive capacity
      {
         HT herm( 3UL, 5UL );
         herm(0,0) = cplx(1,0);
         herm(0,2) = cplx(2,2);
         herm(1,1) = cplx(3,0);
         herm(2,2) = cplx(4,0);

         herm.shrinkToFit();

         checkRows    ( herm, 3UL );
         checkColumns ( herm, 3UL );
         checkCapacity( herm, 5UL );
         checkNonZeros( herm, 5UL );
         checkNonZeros( herm, 0UL, 2UL );
         checkNonZeros( herm, 1UL, 1UL );
         checkNonZeros( herm, 2UL, 2UL );

         if( herm.capacity() != herm.nonZeros() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Shrinking the matrix failed\n"
                << " Details:\n"
                << "   Capacity         : " << herm.capacity() << "\n"
                << "   Expected capacity: " << herm.nonZeros() << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx(1, 0) || herm(0,1) != cplx(0,0) || herm(0,2) != cplx(2,2) ||
             herm(1,0) != cplx(0, 0) || herm(1,1) != cplx(3,0) || herm(1,2) != cplx(0,0) ||
             herm(2,0) != cplx(2,-2) || herm(2,1) != cplx(0,0) || herm(2,2) != cplx(4,0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Shrinking the matrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (1, 0) (0,0) (2,2) )\n"
                                        "( (0, 0) (3,0) (0,0) )\n"
                                        "( (2,-2) (0,0) (4,0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Shrinking a matrix with excessive capacity
      {
         HT herm( 3UL, 100UL );
         herm(0,0) = cplx(1,0);
         herm(0,2) = cplx(2,2);
         herm(1,1) = cplx(3,0);
         herm(2,2) = cplx(4,0);

         herm.shrinkToFit();

         checkRows    ( herm, 3UL );
         checkColumns ( herm, 3UL );
         checkCapacity( herm, 5UL );
         checkNonZeros( herm, 5UL );
         checkNonZeros( herm, 0UL, 2UL );
         checkNonZeros( herm, 1UL, 1UL );
         checkNonZeros( herm, 2UL, 2UL );

         if( herm.capacity() != herm.nonZeros() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Shrinking the matrix failed\n"
                << " Details:\n"
                << "   Capacity         : " << herm.capacity() << "\n"
                << "   Expected capacity: " << herm.nonZeros() << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx(1, 0) || herm(0,1) != cplx(0,0) || herm(0,2) != cplx(2,2) ||
             herm(1,0) != cplx(0, 0) || herm(1,1) != cplx(3,0) || herm(1,2) != cplx(0,0) ||
             herm(2,0) != cplx(2,-2) || herm(2,1) != cplx(0,0) || herm(2,2) != cplx(4,0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Shrinking the matrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (1, 0) (0,0) (2,2) )\n"
                                        "( (0, 0) (3,0) (0,0) )\n"
                                        "( (2,-2) (0,0) (4,0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major HermitianMatrix::shrinkToFit()";

      // Shrinking a matrix without excessive capacity
      {
         OHT herm( 3UL, 5UL );
         herm(0,0) = cplx(1,0);
         herm(0,2) = cplx(2,2);
         herm(1,1) = cplx(3,0);
         herm(2,2) = cplx(4,0);

         herm.shrinkToFit();

         checkRows    ( herm, 3UL );
         checkColumns ( herm, 3UL );
         checkCapacity( herm, 5UL );
         checkNonZeros( herm, 5UL );
         checkNonZeros( herm, 0UL, 2UL );
         checkNonZeros( herm, 1UL, 1UL );
         checkNonZeros( herm, 2UL, 2UL );

         if( herm.capacity() != herm.nonZeros() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Shrinking the matrix failed\n"
                << " Details:\n"
                << "   Capacity         : " << herm.capacity() << "\n"
                << "   Expected capacity: " << herm.nonZeros() << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx(1, 0) || herm(0,1) != cplx(0,0) || herm(0,2) != cplx(2,2) ||
             herm(1,0) != cplx(0, 0) || herm(1,1) != cplx(3,0) || herm(1,2) != cplx(0,0) ||
             herm(2,0) != cplx(2,-2) || herm(2,1) != cplx(0,0) || herm(2,2) != cplx(4,0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Shrinking the matrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (1, 0) (0,0) (2,2) )\n"
                                        "( (0, 0) (3,0) (0,0) )\n"
                                        "( (2,-2) (0,0) (4,0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Shrinking a matrix with excessive capacity
      {
         OHT herm( 3UL, 100UL );
         herm(0,0) = cplx(1,0);
         herm(0,2) = cplx(2,2);
         herm(1,1) = cplx(3,0);
         herm(2,2) = cplx(4,0);

         herm.shrinkToFit();

         checkRows    ( herm, 3UL );
         checkColumns ( herm, 3UL );
         checkCapacity( herm, 5UL );
         checkNonZeros( herm, 5UL );
         checkNonZeros( herm, 0UL, 2UL );
         checkNonZeros( herm, 1UL, 1UL );
         checkNonZeros( herm, 2UL, 2UL );

         if( herm.capacity() != herm.nonZeros() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Shrinking the matrix failed\n"
                << " Details:\n"
                << "   Capacity         : " << herm.capacity() << "\n"
                << "   Expected capacity: " << herm.nonZeros() << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx(1, 0) || herm(0,1) != cplx(0,0) || herm(0,2) != cplx(2,2) ||
             herm(1,0) != cplx(0, 0) || herm(1,1) != cplx(3,0) || herm(1,2) != cplx(0,0) ||
             herm(2,0) != cplx(2,-2) || herm(2,1) != cplx(0,0) || herm(2,2) != cplx(4,0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Shrinking the matrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (1, 0) (0,0) (2,2) )\n"
                                        "( (0, 0) (3,0) (0,0) )\n"
                                        "( (2,-2) (0,0) (4,0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c swap() functionality of the HermitianMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c swap() function of the HermitianMatrix specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseComplexTest::testSwap()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major HermitianMatrix swap";

      HT herm1( 2UL );
      herm1(0,0) = cplx(1,0);
      herm1(0,1) = cplx(2,1);
      herm1(1,1) = cplx(3,0);

      HT herm2( 2UL );
      herm2(0,0) = cplx(4,0);
      herm2(0,1) = cplx(5,1);

      swap( herm1, herm2 );

      checkRows    ( herm1, 2UL );
      checkColumns ( herm1, 2UL );
      checkCapacity( herm1, 4UL );
      checkNonZeros( herm1, 3UL );
      checkNonZeros( herm1, 0UL, 2UL );
      checkNonZeros( herm1, 1UL, 1UL );

      if( herm1(0,0) != cplx(4, 0) || herm1(0,1) != cplx(5,1) ||
          herm1(1,0) != cplx(5,-1) || herm1(1,1) != cplx(0,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the first matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << herm1 << "\n"
             << "   Expected result:\n( (4, 0) (5,1) )\n"
                                     "( (5,-1) (0,0) )\n";
         throw std::runtime_error( oss.str() );
      }

      checkRows    ( herm2, 2UL );
      checkColumns ( herm2, 2UL );
      checkCapacity( herm2, 4UL );
      checkNonZeros( herm2, 4UL );
      checkNonZeros( herm2, 0UL, 2UL );
      checkNonZeros( herm2, 1UL, 2UL );

      if( herm2(0,0) != cplx(1, 0) || herm2(0,1) != cplx(2,1) ||
          herm2(1,0) != cplx(2,-1) || herm2(1,1) != cplx(3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the second matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n( (1, 0) (2,1) )\n"
                                     "( (2,-1) (3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major HermitianMatrix swap";

      OHT herm1( 2UL );
      herm1(0,0) = cplx(1,0);
      herm1(0,1) = cplx(2,1);
      herm1(1,1) = cplx(3,0);

      OHT herm2( 2UL );
      herm2(0,0) = cplx(4,0);
      herm2(0,1) = cplx(5,1);

      swap( herm1, herm2 );

      checkRows    ( herm1, 2UL );
      checkColumns ( herm1, 2UL );
      checkCapacity( herm1, 4UL );
      checkNonZeros( herm1, 3UL );
      checkNonZeros( herm1, 0UL, 2UL );
      checkNonZeros( herm1, 1UL, 1UL );

      if( herm1(0,0) != cplx(4, 0) || herm1(0,1) != cplx(5,1) ||
          herm1(1,0) != cplx(5,-1) || herm1(1,1) != cplx(0,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the first matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << herm1 << "\n"
             << "   Expected result:\n( (4, 0) (5,1) )\n"
                                     "( (5,-1) (0,0) )\n";
         throw std::runtime_error( oss.str() );
      }

      checkRows    ( herm2, 2UL );
      checkColumns ( herm2, 2UL );
      checkCapacity( herm2, 4UL );
      checkNonZeros( herm2, 4UL );
      checkNonZeros( herm2, 0UL, 2UL );
      checkNonZeros( herm2, 1UL, 2UL );

      if( herm2(0,0) != cplx(1, 0) || herm2(0,1) != cplx(2,1) ||
          herm2(1,0) != cplx(2,-1) || herm2(1,1) != cplx(3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the second matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n( (1, 0) (2,1) )\n"
                                     "( (2,-1) (3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c set() member function of the HermitianMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c set() member function of the HermitianMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseComplexTest::testSet()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major HermitianMatrix::set()";

      using Iterator = HT::Iterator;

      // Initialization check
      HT herm( 4UL );

      checkRows    ( herm, 4UL );
      checkColumns ( herm, 4UL );
      checkNonZeros( herm, 0UL );
      checkNonZeros( herm, 0UL, 0UL );
      checkNonZeros( herm, 1UL, 0UL );
      checkNonZeros( herm, 2UL, 0UL );
      checkNonZeros( herm, 3UL, 0UL );

      // Setting a non-zero element
      {
         Iterator pos = herm.set( 2UL, 1UL, cplx(1,1) );

         checkRows    ( herm, 4UL );
         checkColumns ( herm, 4UL );
         checkCapacity( herm, 2UL );
         checkNonZeros( herm, 2UL );
         checkNonZeros( herm, 0UL, 0UL );
         checkNonZeros( herm, 1UL, 1UL );
         checkNonZeros( herm, 2UL, 1UL );
         checkNonZeros( herm, 3UL, 0UL );

         if( pos->value() != cplx(1,1) || pos->index() != 1UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: (1,1)\n"
                << "   Expected index: 1\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(1,2) != cplx(1,-1) || herm(2,1) != cplx(1,1) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (0,0) (0,0) (0, 0) (0,0) )\n"
                                        "( (0,0) (0,0) (1,-1) (0,0) )\n"
                                        "( (0,0) (1,1) (0, 0) (0,0) )\n"
                                        "( (0,0) (0,0) (0, 0) (0,0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setting a second non-zero element
      {
         Iterator pos = herm.set( 2UL, 2UL, cplx(2,0) );

         checkRows    ( herm, 4UL );
         checkColumns ( herm, 4UL );
         checkCapacity( herm, 3UL );
         checkNonZeros( herm, 3UL );
         checkNonZeros( herm, 0UL, 0UL );
         checkNonZeros( herm, 1UL, 1UL );
         checkNonZeros( herm, 2UL, 2UL );
         checkNonZeros( herm, 3UL, 0UL );

         if( pos->value() != cplx(2,0) || pos->index() != 2UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: (2,0)\n"
                << "   Expected index: 2\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(1,2) != cplx(1,-1) || herm(2,1) != cplx(1,1) || herm(2,2) != cplx(2,0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (0,0) (0,0) (0, 0) (0,0) )\n"
                                        "( (0,0) (0,0) (1,-1) (0,0) )\n"
                                        "( (0,0) (1,1) (2, 0) (0,0) )\n"
                                        "( (0,0) (0,0) (0, 0) (0,0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setting a third non-zero element
      {
         Iterator pos = herm.set( 2UL, 0UL, cplx(3,3) );

         checkRows    ( herm, 4UL );
         checkColumns ( herm, 4UL );
         checkCapacity( herm, 5UL );
         checkNonZeros( herm, 5UL );
         checkNonZeros( herm, 0UL, 1UL );
         checkNonZeros( herm, 1UL, 1UL );
         checkNonZeros( herm, 2UL, 3UL );
         checkNonZeros( herm, 3UL, 0UL );

         if( pos->value() != cplx(3,3) || pos->index() != 0UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: (3,3)\n"
                << "   Expected index: 0\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,2) != cplx(3,-3) || herm(1,2) != cplx(1,-1) ||
             herm(2,0) != cplx(3,3) || herm(2,1) != cplx(1,1) || herm(2,2) != cplx(2,0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (0,0) (0,0) (3,-3) (0,0) )\n"
                                        "( (0,0) (0,0) (1,-1) (0,0) )\n"
                                        "( (3,3) (1,1) (2, 0) (0,0) )\n"
                                        "( (0,0) (0,0) (0, 0) (0,0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setting an already existing element
      {
         Iterator pos = herm.set( 1UL, 2UL, cplx(4,4) );

         checkRows    ( herm, 4UL );
         checkColumns ( herm, 4UL );
         checkCapacity( herm, 5UL );
         checkNonZeros( herm, 5UL );
         checkNonZeros( herm, 0UL, 1UL );
         checkNonZeros( herm, 1UL, 1UL );
         checkNonZeros( herm, 2UL, 3UL );
         checkNonZeros( herm, 3UL, 0UL );

         if( pos->value() != cplx(4,4) || pos->index() != 2UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: (4,4)\n"
                << "   Expected index: 2\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,2) != cplx(3,-3) || herm(1,2) != cplx(4,4) ||
             herm(2,0) != cplx(3,3) || herm(2,1) != cplx(4,-4) || herm(2,2) != cplx(2,0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (0,0) (0, 0) (3,-3) (0,0) )\n"
                                        "( (0,0) (0, 0) (4, 4) (0,0) )\n"
                                        "( (3,3) (4,-4) (2, 0) (0,0) )\n"
                                        "( (0,0) (0, 0) (0, 0) (0,0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to set an invalid diagonal element
      try {
         herm.set( 1UL, 1UL, cplx(5,5) );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setting an invalid diagonal element succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (0,0) (0,0) (3,-3) (0,0) )\n"
                                     "( (0,0) (0,0) (1,-1) (0,0) )\n"
                                     "( (3,3) (1,1) (2, 0) (0,0) )\n"
                                     "( (0,0) (0,0) (0, 0) (0,0) )\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major HermitianMatrix::set()";

      using Iterator = OHT::Iterator;

      // Initialization check
      OHT herm( 4UL );

      checkRows    ( herm, 4UL );
      checkColumns ( herm, 4UL );
      checkNonZeros( herm, 0UL );
      checkNonZeros( herm, 0UL, 0UL );
      checkNonZeros( herm, 1UL, 0UL );
      checkNonZeros( herm, 2UL, 0UL );
      checkNonZeros( herm, 3UL, 0UL );

      // Setting a non-zero element
      {
         Iterator pos = herm.set( 2UL, 1UL, cplx(1,1) );

         checkRows    ( herm, 4UL );
         checkColumns ( herm, 4UL );
         checkCapacity( herm, 2UL );
         checkNonZeros( herm, 2UL );
         checkNonZeros( herm, 0UL, 0UL );
         checkNonZeros( herm, 1UL, 1UL );
         checkNonZeros( herm, 2UL, 1UL );
         checkNonZeros( herm, 3UL, 0UL );

         if( pos->value() != cplx(1,1) || pos->index() != 2UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: (1,1)\n"
                << "   Expected index: 2\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(1,2) != cplx(1,-1) || herm(2,1) != cplx(1,1) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (0,0) (0,0) (0, 0) (0,0) )\n"
                                        "( (0,0) (0,0) (1,-1) (0,0) )\n"
                                        "( (0,0) (1,1) (0, 0) (0,0) )\n"
                                        "( (0,0) (0,0) (0, 0) (0,0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setting a second non-zero element
      {
         Iterator pos = herm.set( 2UL, 2UL, cplx(2,0) );

         checkRows    ( herm, 4UL );
         checkColumns ( herm, 4UL );
         checkCapacity( herm, 3UL );
         checkNonZeros( herm, 3UL );
         checkNonZeros( herm, 0UL, 0UL );
         checkNonZeros( herm, 1UL, 1UL );
         checkNonZeros( herm, 2UL, 2UL );
         checkNonZeros( herm, 3UL, 0UL );

         if( pos->value() != cplx(2,0) || pos->index() != 2UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: (2,0)\n"
                << "   Expected index: 2\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(1,2) != cplx(1,-1) || herm(2,1) != cplx(1,1) || herm(2,2) != cplx(2,0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (0,0) (0,0) (0, 0) (0,0) )\n"
                                        "( (0,0) (0,0) (1,-1) (0,0) )\n"
                                        "( (0,0) (1,1) (2, 0) (0,0) )\n"
                                        "( (0,0) (0,0) (0, 0) (0,0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setting a third non-zero element
      {
         Iterator pos = herm.set( 2UL, 0UL, cplx(3,3) );

         checkRows    ( herm, 4UL );
         checkColumns ( herm, 4UL );
         checkCapacity( herm, 5UL );
         checkNonZeros( herm, 5UL );
         checkNonZeros( herm, 0UL, 1UL );
         checkNonZeros( herm, 1UL, 1UL );
         checkNonZeros( herm, 2UL, 3UL );
         checkNonZeros( herm, 3UL, 0UL );

         if( pos->value() != cplx(3,3) || pos->index() != 2UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: (3,3)\n"
                << "   Expected index: 2\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,2) != cplx(3,-3) || herm(1,2) != cplx(1,-1) ||
             herm(2,0) != cplx(3,3) || herm(2,1) != cplx(1,1) || herm(2,2) != cplx(2,0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (0,0) (0,0) (3,-3) (0,0) )\n"
                                        "( (0,0) (0,0) (1,-1) (0,0) )\n"
                                        "( (3,3) (1,1) (2, 0) (0,0) )\n"
                                        "( (0,0) (0,0) (0, 0) (0,0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setting an already existing element
      {
         Iterator pos = herm.set( 1UL, 2UL, cplx(4,4) );

         checkRows    ( herm, 4UL );
         checkColumns ( herm, 4UL );
         checkCapacity( herm, 5UL );
         checkNonZeros( herm, 5UL );
         checkNonZeros( herm, 0UL, 1UL );
         checkNonZeros( herm, 1UL, 1UL );
         checkNonZeros( herm, 2UL, 3UL );
         checkNonZeros( herm, 3UL, 0UL );

         if( pos->value() != cplx(4,4) || pos->index() != 1UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: (4,4)\n"
                << "   Expected index: 1\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,2) != cplx(3,-3) || herm(1,2) != cplx(4,4) ||
             herm(2,0) != cplx(3,3) || herm(2,1) != cplx(4,-4) || herm(2,2) != cplx(2,0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (0,0) (0, 0) (3,-3) (0,0) )\n"
                                        "( (0,0) (0, 0) (4, 4) (0,0) )\n"
                                        "( (3,3) (4,-4) (2, 0) (0,0) )\n"
                                        "( (0,0) (0, 0) (0, 0) (0,0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to set an invalid diagonal element
      try {
         herm.set( 1UL, 1UL, cplx(5,5) );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setting an invalid diagonal element succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (0,0) (0,0) (3,-3) (0,0) )\n"
                                     "( (0,0) (0,0) (1,-1) (0,0) )\n"
                                     "( (3,3) (1,1) (2, 0) (0,0) )\n"
                                     "( (0,0) (0,0) (0, 0) (0,0) )\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c insert() member function of the HermitianMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c insert() member function of the HermitianMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseComplexTest::testInsert()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major HermitianMatrix::insert()";

      using Iterator = HT::Iterator;

      // Initialization check
      HT herm( 4UL );

      checkRows    ( herm, 4UL );
      checkColumns ( herm, 4UL );
      checkNonZeros( herm, 0UL );
      checkNonZeros( herm, 0UL, 0UL );
      checkNonZeros( herm, 1UL, 0UL );
      checkNonZeros( herm, 2UL, 0UL );
      checkNonZeros( herm, 3UL, 0UL );

      // Inserting a non-zero element
      {
         Iterator pos = herm.insert( 2UL, 1UL, cplx(1,1) );

         checkRows    ( herm, 4UL );
         checkColumns ( herm, 4UL );
         checkCapacity( herm, 2UL );
         checkNonZeros( herm, 2UL );
         checkNonZeros( herm, 0UL, 0UL );
         checkNonZeros( herm, 1UL, 1UL );
         checkNonZeros( herm, 2UL, 1UL );
         checkNonZeros( herm, 3UL, 0UL );

         if( pos->value() != cplx(1,1) || pos->index() != 1UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: (1,1)\n"
                << "   Expected index: 1\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(1,2) != cplx(1,-1) || herm(2,1) != cplx(1,1) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (0,0) (0,0) (0, 0) (0,0) )\n"
                                        "( (0,0) (0,0) (1,-1) (0,0) )\n"
                                        "( (0,0) (1,1) (0, 0) (0,0) )\n"
                                        "( (0,0) (0,0) (0, 0) (0,0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting a second non-zero element
      {
         Iterator pos = herm.insert( 2UL, 2UL, cplx(2,0) );

         checkRows    ( herm, 4UL );
         checkColumns ( herm, 4UL );
         checkCapacity( herm, 3UL );
         checkNonZeros( herm, 3UL );
         checkNonZeros( herm, 0UL, 0UL );
         checkNonZeros( herm, 1UL, 1UL );
         checkNonZeros( herm, 2UL, 2UL );
         checkNonZeros( herm, 3UL, 0UL );

         if( pos->value() != cplx(2,0) || pos->index() != 2UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: (2,0)\n"
                << "   Expected index: 2\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(1,2) != cplx(1,-1) || herm(2,1) != cplx(1,1) || herm(2,2) != cplx(2,0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (0,0) (0,0) (0, 0) (0,0) )\n"
                                        "( (0,0) (0,0) (1,-1) (0,0) )\n"
                                        "( (0,0) (1,1) (2, 0) (0,0) )\n"
                                        "( (0,0) (0,0) (0, 0) (0,0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting a third non-zero element
      {
         Iterator pos = herm.insert( 2UL, 0UL, cplx(3,3) );

         checkRows    ( herm, 4UL );
         checkColumns ( herm, 4UL );
         checkCapacity( herm, 5UL );
         checkNonZeros( herm, 5UL );
         checkNonZeros( herm, 0UL, 1UL );
         checkNonZeros( herm, 1UL, 1UL );
         checkNonZeros( herm, 2UL, 3UL );
         checkNonZeros( herm, 3UL, 0UL );

         if( pos->value() != cplx(3,3) || pos->index() != 0UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: cplx(3,3)\n"
                << "   Expected index: 0\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,2) != cplx(3,-3) || herm(1,2) != cplx(1,-1) ||
             herm(2,0) != cplx(3,3) || herm(2,1) != cplx(1,1) || herm(2,2) != cplx(2,0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (0,0) (0,0) (3,-3) (0,0) )\n"
                                        "( (0,0) (0,0) (1,-1) (0,0) )\n"
                                        "( (3,3) (1,1) (2, 0) (0,0) )\n"
                                        "( (0,0) (0,0) (0, 0) (0,0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to insert an already existing element
      try {
         herm.insert( 1UL, 2UL, cplx(4,4) );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Inserting an existing element succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (0,0) (0,0) (3,-3) (0,0) )\n"
                                     "( (0,0) (0,0) (1,-1) (0,0) )\n"
                                     "( (3,3) (1,1) (2, 0) (0,0) )\n"
                                     "( (0,0) (0,0) (0, 0) (0,0) )\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      // Trying to insert an invalid diagonal element
      try {
         herm.insert( 1UL, 1UL, cplx(5,5) );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Inserting an invalid diagonal element succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (0,0) (0,0) (3,-3) (0,0) )\n"
                                     "( (0,0) (0,0) (1,-1) (0,0) )\n"
                                     "( (3,3) (1,1) (2, 0) (0,0) )\n"
                                     "( (0,0) (0,0) (0, 0) (0,0) )\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major HermitianMatrix::insert()";

      using Iterator = OHT::Iterator;

      // Initialization check
      OHT herm( 4UL );

      checkRows    ( herm, 4UL );
      checkColumns ( herm, 4UL );
      checkNonZeros( herm, 0UL );
      checkNonZeros( herm, 0UL, 0UL );
      checkNonZeros( herm, 1UL, 0UL );
      checkNonZeros( herm, 2UL, 0UL );
      checkNonZeros( herm, 3UL, 0UL );

      // Inserting a non-zero element
      {
         Iterator pos = herm.insert( 2UL, 1UL, cplx(1,1) );

         checkRows    ( herm, 4UL );
         checkColumns ( herm, 4UL );
         checkCapacity( herm, 2UL );
         checkNonZeros( herm, 2UL );
         checkNonZeros( herm, 0UL, 0UL );
         checkNonZeros( herm, 1UL, 1UL );
         checkNonZeros( herm, 2UL, 1UL );
         checkNonZeros( herm, 3UL, 0UL );

         if( pos->value() != cplx(1,1) || pos->index() != 2UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: (1,1)\n"
                << "   Expected index: 2\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(1,2) != cplx(1,-1) || herm(2,1) != cplx(1,1) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (0,0) (0,0) (0, 0) (0,0) )\n"
                                        "( (0,0) (0,0) (1,-1) (0,0) )\n"
                                        "( (0,0) (1,1) (0, 0) (0,0) )\n"
                                        "( (0,0) (0,0) (0, 0) (0,0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting a second non-zero element
      {
         Iterator pos = herm.insert( 2UL, 2UL, cplx(2,0) );

         checkRows    ( herm, 4UL );
         checkColumns ( herm, 4UL );
         checkCapacity( herm, 3UL );
         checkNonZeros( herm, 3UL );
         checkNonZeros( herm, 0UL, 0UL );
         checkNonZeros( herm, 1UL, 1UL );
         checkNonZeros( herm, 2UL, 2UL );
         checkNonZeros( herm, 3UL, 0UL );

         if( pos->value() != cplx(2,0) || pos->index() != 2UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: (2,0)\n"
                << "   Expected index: 2\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(1,2) != cplx(1,-1) || herm(2,1) != cplx(1,1) || herm(2,2) != cplx(2,0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (0,0) (0,0) (0, 0) (0,0) )\n"
                                        "( (0,0) (0,0) (1,-1) (0,0) )\n"
                                        "( (0,0) (1,1) (2, 0) (0,0) )\n"
                                        "( (0,0) (0,0) (0, 0) (0,0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting a third non-zero element
      {
         Iterator pos = herm.insert( 2UL, 0UL, cplx(3,3) );

         checkRows    ( herm, 4UL );
         checkColumns ( herm, 4UL );
         checkCapacity( herm, 5UL );
         checkNonZeros( herm, 5UL );
         checkNonZeros( herm, 0UL, 1UL );
         checkNonZeros( herm, 1UL, 1UL );
         checkNonZeros( herm, 2UL, 3UL );
         checkNonZeros( herm, 3UL, 0UL );

         if( pos->value() != cplx(3,3) || pos->index() != 2UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: cplx(3,3)\n"
                << "   Expected index: 2\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,2) != cplx(3,-3) || herm(1,2) != cplx(1,-1) ||
             herm(2,0) != cplx(3,3) || herm(2,1) != cplx(1,1) || herm(2,2) != cplx(2,0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (0,0) (0,0) (3,-3) (0,0) )\n"
                                        "( (0,0) (0,0) (1,-1) (0,0) )\n"
                                        "( (3,3) (1,1) (2, 0) (0,0) )\n"
                                        "( (0,0) (0,0) (0, 0) (0,0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to insert an already existing element
      try {
         herm.insert( 1UL, 2UL, cplx(4,4) );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Inserting an existing element succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (0,0) (0,0) (3,-3) (0,0) )\n"
                                     "( (0,0) (0,0) (1,-1) (0,0) )\n"
                                     "( (3,3) (1,1) (2, 0) (0,0) )\n"
                                     "( (0,0) (0,0) (0, 0) (0,0) )\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      // Trying to insert an invalid diagonal element
      try {
         herm.insert( 1UL, 1UL, cplx(5,5) );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Inserting an invalid diagonal element succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (0,0) (0,0) (3,-3) (0,0) )\n"
                                     "( (0,0) (0,0) (1,-1) (0,0) )\n"
                                     "( (3,3) (1,1) (2, 0) (0,0) )\n"
                                     "( (0,0) (0,0) (0, 0) (0,0) )\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c append() member function of the HermitianMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c append() member function of the HermitianMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseComplexTest::testAppend()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major HermitianMatrix::append()";

      // Appending with pre-allocation in each row
      {
         // Initialization check
         HT herm( 4UL, 10UL );
         herm.reserve( 0UL, 2UL );
         herm.reserve( 1UL, 2UL );
         herm.reserve( 2UL, 2UL );
         herm.reserve( 3UL, 4UL );

         checkRows    ( herm, 4UL );
         checkColumns ( herm, 4UL );
         checkCapacity( herm, 9UL );
         checkNonZeros( herm, 0UL );
         checkNonZeros( herm, 0UL, 0UL );
         checkNonZeros( herm, 1UL, 0UL );
         checkNonZeros( herm, 2UL, 0UL );
         checkNonZeros( herm, 3UL, 0UL );

         // Appending one non-zero element
         herm.append( 2UL, 1UL, cplx(1,1) );

         checkRows    ( herm, 4UL );
         checkColumns ( herm, 4UL );
         checkCapacity( herm, 9UL );
         checkNonZeros( herm, 2UL );
         checkNonZeros( herm, 0UL, 0UL );
         checkNonZeros( herm, 1UL, 1UL );
         checkNonZeros( herm, 2UL, 1UL );
         checkNonZeros( herm, 3UL, 0UL );

         if( herm(1,2) != cplx(1,-1) || herm(2,1) != cplx(1,1) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (0,0) (0,0) (0, 0) (0,0) )\n"
                                        "( (0,0) (0,0) (1,-1) (0,0) )\n"
                                        "( (0,0) (1,1) (0, 0) (0,0) )\n"
                                        "( (0,0) (0,0) (0, 0) (0,0) )\n";
            throw std::runtime_error( oss.str() );
         }

         // Appending two more non-zero elements
         herm.append( 0UL, 0UL, cplx(2,0) );
         herm.append( 0UL, 3UL, cplx(3,3) );

         checkRows    ( herm, 4UL );
         checkColumns ( herm, 4UL );
         checkCapacity( herm, 9UL );
         checkNonZeros( herm, 5UL );
         checkNonZeros( herm, 0UL, 2UL );
         checkNonZeros( herm, 1UL, 1UL );
         checkNonZeros( herm, 2UL, 1UL );
         checkNonZeros( herm, 3UL, 1UL );

         if( herm(0,0) != cplx(2, 0) || herm(0,3) != cplx(3,3) ||
             herm(1,2) != cplx(1,-1) ||
             herm(2,1) != cplx(1, 1) ||
             herm(3,0) != cplx(3,-3) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (2, 0) (0,0) (0, 0) (3,3) )\n"
                                        "( (0, 0) (0,0) (1,-1) (0,0) )\n"
                                        "( (0, 0) (1,1) (0, 0) (0,0) )\n"
                                        "( (3,-3) (0,0) (0, 0) (0,0) )\n";
            throw std::runtime_error( oss.str() );
         }

         // Appending two more non-zero elements
         herm.append( 3UL, 1UL, cplx(4,4) );
         herm.append( 3UL, 2UL, cplx(5,5) );

         checkRows    ( herm, 4UL );
         checkColumns ( herm, 4UL );
         checkCapacity( herm, 9UL );
         checkNonZeros( herm, 9UL );
         checkNonZeros( herm, 0UL, 2UL );
         checkNonZeros( herm, 1UL, 2UL );
         checkNonZeros( herm, 2UL, 2UL );
         checkNonZeros( herm, 3UL, 3UL );

         if( herm(0,0) != cplx(2, 0) || herm(0,3) != cplx(3, 3) ||
             herm(1,2) != cplx(1,-1) || herm(1,3) != cplx(4,-4) ||
             herm(2,1) != cplx(1, 1) || herm(2,3) != cplx(5,-5) ||
             herm(3,0) != cplx(3,-3) || herm(3,1) != cplx(4, 4) || herm(3,2) != cplx(5,5) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (2, 0) (0,0) (0, 0) (3, 3) )\n"
                                        "( (0, 0) (0,0) (1,-1) (4,-4) )\n"
                                        "( (0, 0) (1,1) (0, 0) (5,-5) )\n"
                                        "( (3,-3) (4,4) (5, 5) (0, 0) )\n";
            throw std::runtime_error( oss.str() );
         }

         // Trying to append an invalid diagonal element
         try {
            herm.append( 3UL, 3UL, cplx(6,6) );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Appending an invalid diagonal element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (2, 0) (0,0) (0, 0) (3, 3) )\n"
                                        "( (0, 0) (0,0) (1,-1) (4,-4) )\n"
                                        "( (0, 0) (1,1) (0, 0) (5,-5) )\n"
                                        "( (3,-3) (4,4) (5, 5) (0, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Appending with row finalization
      {
         // Initialization check
         HT herm( 4UL, 9UL );
         herm.reserve( 0UL, 2UL );
         herm.reserve( 1UL, 4UL );
         herm.reserve( 2UL, 1UL );
         herm.reserve( 3UL, 2UL );

         // Appending one non-zero element
         herm.append( 0UL, 1UL, cplx(1,1) );
         herm.finalize( 0UL );

         checkRows    ( herm, 4UL );
         checkColumns ( herm, 4UL );
         checkCapacity( herm, 9UL );
         checkNonZeros( herm, 2UL );
         checkNonZeros( herm, 0UL, 1UL );
         checkNonZeros( herm, 1UL, 1UL );
         checkNonZeros( herm, 2UL, 0UL );
         checkNonZeros( herm, 3UL, 0UL );

         if( herm(0,1) != cplx(1,1) || herm(1,0) != cplx(1,-1) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (0, 0) (1,1) (0,0) (0,0) )\n"
                                        "( (1,-1) (0,0) (0,0) (0,0) )\n"
                                        "( (0, 0) (0,0) (0,0) (0,0) )\n"
                                        "( (0, 0) (0,0) (0,0) (0,0) )\n";
            throw std::runtime_error( oss.str() );
         }

         // Appending two more non-zero elements
         herm.append( 1UL, 1UL, cplx(2,0) );
         herm.append( 1UL, 2UL, cplx(3,3) );
         herm.finalize( 1UL );

         checkRows    ( herm, 4UL );
         checkColumns ( herm, 4UL );
         checkCapacity( herm, 9UL );
         checkNonZeros( herm, 5UL );
         checkNonZeros( herm, 0UL, 1UL );
         checkNonZeros( herm, 1UL, 3UL );
         checkNonZeros( herm, 2UL, 1UL );
         checkNonZeros( herm, 3UL, 0UL );

         if( herm(0,1) != cplx(1, 1) ||
             herm(1,0) != cplx(1,-1) || herm(1,1) != cplx(2,0) || herm(1,2) != cplx(3,3) ||
             herm(2,1) != cplx(3,-3) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (0, 0) (1, 1) (0,0) (0,0) )\n"
                                        "( (1,-1) (2, 0) (3,3) (0,0) )\n"
                                        "( (0, 0) (3,-3) (0,0) (0,0) )\n"
                                        "( (0, 0) (0, 0) (0,0) (0,0) )\n";
            throw std::runtime_error( oss.str() );
         }

         // Appending two more non-zero elements
         herm.append( 3UL, 0UL, cplx(4,4) );
         herm.append( 3UL, 1UL, cplx(5,5) );
         herm.finalize( 3UL );

         checkRows    ( herm, 4UL );
         checkColumns ( herm, 4UL );
         checkCapacity( herm, 9UL );
         checkNonZeros( herm, 9UL );
         checkNonZeros( herm, 0UL, 2UL );
         checkNonZeros( herm, 1UL, 4UL );
         checkNonZeros( herm, 2UL, 1UL );
         checkNonZeros( herm, 3UL, 2UL );

         if( herm(0,1) != cplx(1, 1) || herm(0,3) != cplx(4,-4) ||
             herm(1,0) != cplx(1,-1) || herm(1,1) != cplx(2, 0) || herm(1,2) != cplx(3,3) || herm(1,3) != cplx(5,-5) ||
             herm(2,1) != cplx(3,-3) ||
             herm(3,0) != cplx(4, 4) || herm(3,1) != cplx(5, 5) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (0, 0) (1, 1) (0,0) (4,-4) )\n"
                                        "( (1,-1) (2, 0) (3,3) (5,-5) )\n"
                                        "( (0, 0) (3,-3) (0,0) (0, 0) )\n"
                                        "( (4, 4) (5, 5) (0,0) (0, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major HermitianMatrix::append()";

      // Appending with pre-allocation in each column
      {
         // Initialization check
         OHT herm( 4UL, 10UL );
         herm.reserve( 0UL, 2UL );
         herm.reserve( 1UL, 2UL );
         herm.reserve( 2UL, 2UL );
         herm.reserve( 3UL, 4UL );

         checkRows    ( herm, 4UL );
         checkColumns ( herm, 4UL );
         checkCapacity( herm, 9UL );
         checkNonZeros( herm, 0UL );
         checkNonZeros( herm, 0UL, 0UL );
         checkNonZeros( herm, 1UL, 0UL );
         checkNonZeros( herm, 2UL, 0UL );
         checkNonZeros( herm, 3UL, 0UL );

         // Appending one non-zero element
         herm.append( 1UL, 2UL, cplx(1,1) );

         checkRows    ( herm, 4UL );
         checkColumns ( herm, 4UL );
         checkCapacity( herm, 9UL );
         checkNonZeros( herm, 2UL );
         checkNonZeros( herm, 0UL, 0UL );
         checkNonZeros( herm, 1UL, 1UL );
         checkNonZeros( herm, 2UL, 1UL );
         checkNonZeros( herm, 3UL, 0UL );

         if( herm(1,2) != cplx(1,1) || herm(2,1) != cplx(1,-1) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (0,0) (0, 0) (0,0) (0,0) )\n"
                                        "( (0,0) (0, 0) (1,1) (0,0) )\n"
                                        "( (0,0) (1,-1) (0,0) (0,0) )\n"
                                        "( (0,0) (0, 0) (0,0) (0,0) )\n";
            throw std::runtime_error( oss.str() );
         }

         // Appending two more non-zero elements
         herm.append( 0UL, 0UL, cplx(2,0) );
         herm.append( 3UL, 0UL, cplx(3,3) );

         checkRows    ( herm, 4UL );
         checkColumns ( herm, 4UL );
         checkCapacity( herm, 9UL );
         checkNonZeros( herm, 5UL );
         checkNonZeros( herm, 0UL, 2UL );
         checkNonZeros( herm, 1UL, 1UL );
         checkNonZeros( herm, 2UL, 1UL );
         checkNonZeros( herm, 3UL, 1UL );

         if( herm(0,0) != cplx(2, 0) || herm(0,3) != cplx(3,-3) ||
             herm(1,2) != cplx(1, 1) ||
             herm(2,1) != cplx(1,-1) ||
             herm(3,0) != cplx(3, 3) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (2,0) (0, 0) (0,0) (3,-3) )\n"
                                        "( (0,0) (0, 0) (1,1) (0, 0) )\n"
                                        "( (0,0) (1,-1) (0,0) (0, 0) )\n"
                                        "( (3,3) (0, 0) (0,0) (0, 0) )\n";
            throw std::runtime_error( oss.str() );
         }

         // Appending two more non-zero elements
         herm.append( 1UL, 3UL, cplx(4,4) );
         herm.append( 2UL, 3UL, cplx(5,5) );

         checkRows    ( herm, 4UL );
         checkColumns ( herm, 4UL );
         checkCapacity( herm, 9UL );
         checkNonZeros( herm, 9UL );
         checkNonZeros( herm, 0UL, 2UL );
         checkNonZeros( herm, 1UL, 2UL );
         checkNonZeros( herm, 2UL, 2UL );
         checkNonZeros( herm, 3UL, 3UL );

         if( herm(0,0) != cplx(2, 0) || herm(0,3) != cplx(3,-3) ||
             herm(1,2) != cplx(1, 1) || herm(1,3) != cplx(4, 4) ||
             herm(2,1) != cplx(1,-1) || herm(2,3) != cplx(5, 5) ||
             herm(3,0) != cplx(3, 3) || herm(3,1) != cplx(4,-4) || herm(3,2) != cplx(5,-5) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (2,0) (0, 0) (0, 0) (3,-3) )\n"
                                        "( (0,0) (0, 0) (1, 1) (4, 4) )\n"
                                        "( (0,0) (1,-1) (0, 0) (5, 5) )\n"
                                        "( (3,3) (4,-4) (5,-5) (0, 0) )\n";
            throw std::runtime_error( oss.str() );
         }

         // Trying to append an invalid diagonal element
         try {
            herm.append( 3UL, 3UL, cplx(6,6) );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Appending an invalid diagonal element succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (2,0) (0, 0) (0, 0) (3,-3) )\n"
                                        "( (0,0) (0, 0) (1, 1) (4, 4) )\n"
                                        "( (0,0) (1,-1) (0, 0) (5, 5) )\n"
                                        "( (3,3) (4,-4) (5,-5) (0, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Appending with column finalization
      {
         // Initialization check
         OHT herm( 4UL, 9UL );
         herm.reserve( 0UL, 2UL );
         herm.reserve( 1UL, 4UL );
         herm.reserve( 2UL, 1UL );
         herm.reserve( 3UL, 2UL );

         // Appending one non-zero element
         herm.append( 1UL, 0UL, cplx(1,1) );
         herm.finalize( 0UL );

         checkRows    ( herm, 4UL );
         checkColumns ( herm, 4UL );
         checkCapacity( herm, 9UL );
         checkNonZeros( herm, 2UL );
         checkNonZeros( herm, 0UL, 1UL );
         checkNonZeros( herm, 1UL, 1UL );
         checkNonZeros( herm, 2UL, 0UL );
         checkNonZeros( herm, 3UL, 0UL );

         if( herm(0,1) != cplx(1,-1) || herm(1,0) != cplx(1,1) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (0,0) (1,-1) (0,0) (0,0) )\n"
                                        "( (1,1) (0, 0) (0,0) (0,0) )\n"
                                        "( (0,0) (0, 0) (0,0) (0,0) )\n"
                                        "( (0,0) (0, 0) (0,0) (0,0) )\n";
            throw std::runtime_error( oss.str() );
         }

         // Appending two more non-zero elements
         herm.append( 1UL, 1UL, cplx(2,0) );
         herm.append( 2UL, 1UL, cplx(3,3) );
         herm.finalize( 1UL );

         checkRows    ( herm, 4UL );
         checkColumns ( herm, 4UL );
         checkCapacity( herm, 9UL );
         checkNonZeros( herm, 5UL );
         checkNonZeros( herm, 0UL, 1UL );
         checkNonZeros( herm, 1UL, 3UL );
         checkNonZeros( herm, 2UL, 1UL );
         checkNonZeros( herm, 3UL, 0UL );

         if( herm(0,1) != cplx(1,-1) ||
             herm(1,0) != cplx(1, 1) || herm(1,1) != cplx(2,0) || herm(1,2) != cplx(3,-3) ||
             herm(2,1) != cplx(3, 3) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (0,0) (1,-1) (0, 0) (0,0) )\n"
                                        "( (1,1) (2, 0) (3,-3) (0,0) )\n"
                                        "( (0,0) (3, 3) (0, 0) (0,0) )\n"
                                        "( (0,0) (0, 0) (0, 0) (0,0) )\n";
            throw std::runtime_error( oss.str() );
         }

         // Appending two more non-zero elements
         herm.append( 0UL, 3UL, cplx(4,4) );
         herm.append( 1UL, 3UL, cplx(5,5) );
         herm.finalize( 3UL );

         checkRows    ( herm, 4UL );
         checkColumns ( herm, 4UL );
         checkCapacity( herm, 9UL );
         checkNonZeros( herm, 9UL );
         checkNonZeros( herm, 0UL, 2UL );
         checkNonZeros( herm, 1UL, 4UL );
         checkNonZeros( herm, 2UL, 1UL );
         checkNonZeros( herm, 3UL, 2UL );

         if( herm(0,1) != cplx(1,-1) || herm(0,3) != cplx(4, 4) ||
             herm(1,0) != cplx(1, 1) || herm(1,1) != cplx(2, 0) || herm(1,2) != cplx(3,-3) || herm(1,3) != cplx(5,5) ||
             herm(2,1) != cplx(3, 3) ||
             herm(3,0) != cplx(4,-4) || herm(3,1) != cplx(5,-5) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Append operation failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (0, 0) (1,-1) (0, 0) (4,4) )\n"
                                        "( (1, 1) (2, 0) (3,-3) (5,5) )\n"
                                        "( (0, 0) (3, 3) (0, 0) (0,0) )\n"
                                        "( (4,-4) (5,-5) (0, 0) (0,0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c erase() member function of the HermitianMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c erase() member function of the HermitianMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseComplexTest::testErase()
{
   //=====================================================================================
   // Row-major index-based erase function
   //=====================================================================================

   {
      test_ = "Row-major HermitianMatrix::erase( size_t, size_t )";

      // Initialization check
      HT herm( 4UL );
      herm(0,0) = cplx(1,0);
      herm(0,2) = cplx(2,2);
      herm(0,3) = cplx(3,3);
      herm(1,1) = cplx(4,0);
      herm(1,2) = cplx(5,5);
      herm(2,2) = cplx(6,0);
      herm(2,3) = cplx(7,7);

      checkRows    ( herm,  4UL );
      checkColumns ( herm,  4UL );
      checkCapacity( herm, 11UL );
      checkNonZeros( herm, 11UL );
      checkNonZeros( herm, 0UL, 3UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 4UL );
      checkNonZeros( herm, 3UL, 2UL );

      if( herm(0,0) != cplx(1, 0) || herm(0,2) != cplx(2, 2) || herm(0,3) != cplx(3,3) ||
          herm(1,1) != cplx(4, 0) || herm(1,2) != cplx(5, 5) ||
          herm(2,0) != cplx(2,-2) || herm(2,1) != cplx(5,-5) || herm(2,2) != cplx(6,0) || herm(2,3) != cplx(7,7) ||
          herm(3,0) != cplx(3,-3) || herm(3,2) != cplx(7,-7) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (1, 0) (0, 0) (2, 2) (3,3) )\n"
                                     "( (0, 0) (4, 0) (5, 5) (0,0) )\n"
                                     "( (2,-2) (5,-5) (6, 0) (7,7) )\n"
                                     "( (3,-3) (0, 0) (7,-7) (0,0) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element at (0,0)
      herm.erase( 0UL, 0UL );

      checkRows    ( herm,  4UL );
      checkColumns ( herm,  4UL );
      checkCapacity( herm, 11UL );
      checkNonZeros( herm, 10UL );
      checkNonZeros( herm, 0UL, 2UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 4UL );
      checkNonZeros( herm, 3UL, 2UL );

      if( herm(0,2) != cplx(2, 2) || herm(0,3) != cplx(3, 3) ||
          herm(1,1) != cplx(4, 0) || herm(1,2) != cplx(5, 5) ||
          herm(2,0) != cplx(2,-2) || herm(2,1) != cplx(5,-5) || herm(2,2) != cplx(6,0) || herm(2,3) != cplx(7,7) ||
          herm(3,0) != cplx(3,-3) || herm(3,2) != cplx(7,-7) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (0, 0) (0, 0) (2, 2) (3,3) )\n"
                                     "( (0, 0) (4, 0) (5, 5) (0,0) )\n"
                                     "( (2,-2) (5,-5) (6, 0) (7,7) )\n"
                                     "( (3,-3) (0, 0) (7,-7) (0,0) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element at (1,2)
      herm.erase( 1UL, 2UL );

      checkRows    ( herm,  4UL );
      checkColumns ( herm,  4UL );
      checkCapacity( herm, 11UL );
      checkNonZeros( herm,  8UL );
      checkNonZeros( herm, 0UL, 2UL );
      checkNonZeros( herm, 1UL, 1UL );
      checkNonZeros( herm, 2UL, 3UL );
      checkNonZeros( herm, 3UL, 2UL );

      if( herm(0,2) != cplx(2, 2) || herm(0,3) != cplx(3, 3) ||
          herm(1,1) != cplx(4, 0) ||
          herm(2,0) != cplx(2,-2) || herm(2,2) != cplx(6, 0) || herm(2,3) != cplx(7,7) ||
          herm(3,0) != cplx(3,-3) || herm(3,2) != cplx(7,-7) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (0, 0) (0,0) (2, 2) (3,3) )\n"
                                     "( (0, 0) (4,0) (0, 0) (0,0) )\n"
                                     "( (2,-2) (0,0) (6, 0) (7,7) )\n"
                                     "( (3,-3) (0,0) (7,-7) (0,0) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element at (0,2)
      herm.erase( 0UL, 2UL );

      checkRows    ( herm,  4UL );
      checkColumns ( herm,  4UL );
      checkCapacity( herm, 11UL );
      checkNonZeros( herm,  6UL );
      checkNonZeros( herm, 0UL, 1UL );
      checkNonZeros( herm, 1UL, 1UL );
      checkNonZeros( herm, 2UL, 2UL );
      checkNonZeros( herm, 3UL, 2UL );

      if( herm(0,3) != cplx(3, 3) ||
          herm(1,1) != cplx(4, 0) ||
          herm(2,2) != cplx(6, 0) || herm(2,3) != cplx(7,7) ||
          herm(3,0) != cplx(3,-3) || herm(3,2) != cplx(7,-7) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (0, 0) (0,0) (0, 0) (3,3) )\n"
                                     "( (0, 0) (4,0) (0, 0) (0,0) )\n"
                                     "( (0, 0) (0,0) (6, 0) (7,7) )\n"
                                     "( (3,-3) (0,0) (7,-7) (0,0) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to erase a zero element
      herm.erase( 0UL, 1UL );

      checkRows    ( herm,  4UL );
      checkColumns ( herm,  4UL );
      checkCapacity( herm, 11UL );
      checkNonZeros( herm,  6UL );
      checkNonZeros( herm, 0UL, 1UL );
      checkNonZeros( herm, 1UL, 1UL );
      checkNonZeros( herm, 2UL, 2UL );
      checkNonZeros( herm, 3UL, 2UL );

      if( herm(0,3) != cplx(3, 3) ||
          herm(1,1) != cplx(4, 0) ||
          herm(2,2) != cplx(6, 0) || herm(2,3) != cplx(7,7) ||
          herm(3,0) != cplx(3,-3) || herm(3,2) != cplx(7,-7) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (0, 0) (0,0) (0, 0) (3,3) )\n"
                                     "( (0, 0) (4,0) (0, 0) (0,0) )\n"
                                     "( (0, 0) (0,0) (6, 0) (7,7) )\n"
                                     "( (3,-3) (0,0) (7,-7) (0,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major iterator-based erase function
   //=====================================================================================

   {
      test_ = "Row-major HermitianMatrix::erase( size_t, Iterator )";

      using Iterator = HT::Iterator;

      // Initialization check
      HT herm( 4UL );
      herm(0,0) = cplx(1,0);
      herm(0,2) = cplx(2,2);
      herm(0,3) = cplx(3,3);
      herm(1,1) = cplx(4,0);
      herm(1,2) = cplx(5,5);
      herm(2,2) = cplx(6,0);
      herm(2,3) = cplx(7,7);

      checkRows    ( herm,  4UL );
      checkColumns ( herm,  4UL );
      checkCapacity( herm, 11UL );
      checkNonZeros( herm, 11UL );
      checkNonZeros( herm, 0UL, 3UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 4UL );
      checkNonZeros( herm, 3UL, 2UL );

      if( herm(0,0) != cplx(1, 0) || herm(0,2) != cplx(2, 2) || herm(0,3) != cplx(3,3) ||
          herm(1,1) != cplx(4, 0) || herm(1,2) != cplx(5, 5) ||
          herm(2,0) != cplx(2,-2) || herm(2,1) != cplx(5,-5) || herm(2,2) != cplx(6,0) || herm(2,3) != cplx(7,7) ||
          herm(3,0) != cplx(3,-3) || herm(3,2) != cplx(7,-7) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (1, 0) (0, 0) (2, 2) (3,3) )\n"
                                     "( (0, 0) (4, 0) (5, 5) (0,0) )\n"
                                     "( (2,-2) (5,-5) (6, 0) (7,7) )\n"
                                     "( (3,-3) (0, 0) (7,-7) (0,0) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element at (0,0)
      {
         Iterator pos = herm.erase( 0UL, herm.find( 0UL, 0UL ) );

         checkRows    ( herm,  4UL );
         checkColumns ( herm,  4UL );
         checkCapacity( herm, 11UL );
         checkNonZeros( herm, 10UL );
         checkNonZeros( herm, 0UL, 2UL );
         checkNonZeros( herm, 1UL, 2UL );
         checkNonZeros( herm, 2UL, 4UL );
         checkNonZeros( herm, 3UL, 2UL );

         if( herm(0,2) != cplx(2, 2) || herm(0,3) != cplx(3, 3) ||
             herm(1,1) != cplx(4, 0) || herm(1,2) != cplx(5, 5) ||
             herm(2,0) != cplx(2,-2) || herm(2,1) != cplx(5,-5) || herm(2,2) != cplx(6,0) || herm(2,3) != cplx(7,7) ||
             herm(3,0) != cplx(3,-3) || herm(3,2) != cplx(7,-7) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (0, 0) (0, 0) (2, 2) (3,3) )\n"
                                        "( (0, 0) (4, 0) (5, 5) (0,0) )\n"
                                        "( (2,-2) (5,-5) (6, 0) (7,7) )\n"
                                        "( (3,-3) (0, 0) (7,-7) (0,0) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos->value() != cplx(2,2) || pos->index() != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: (2,2)\n"
                << "   Expected index: 2\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the element at (1,2)
      {
         Iterator pos = herm.erase( 1UL, herm.find( 1UL, 2UL ) );

         checkRows    ( herm,  4UL );
         checkColumns ( herm,  4UL );
         checkCapacity( herm, 11UL );
         checkNonZeros( herm,  8UL );
         checkNonZeros( herm, 0UL, 2UL );
         checkNonZeros( herm, 1UL, 1UL );
         checkNonZeros( herm, 2UL, 3UL );
         checkNonZeros( herm, 3UL, 2UL );

         if( herm(0,2) != cplx(2, 2) || herm(0,3) != cplx(3, 3) ||
             herm(1,1) != cplx(4, 0) ||
             herm(2,0) != cplx(2,-2) || herm(2,2) != cplx(6, 0) || herm(2,3) != cplx(7,7) ||
             herm(3,0) != cplx(3,-3) || herm(3,2) != cplx(7,-7) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (0, 0) (0,0) (2, 2) (3,3) )\n"
                                        "( (0, 0) (4,0) (0, 0) (0,0) )\n"
                                        "( (2,-2) (0,0) (6, 0) (7,7) )\n"
                                        "( (3,-3) (0,0) (7,-7) (0,0) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos != herm.end( 1UL ) ) {
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
         Iterator pos = herm.erase( 0UL, herm.find( 0UL, 2UL ) );

         checkRows    ( herm,  4UL );
         checkColumns ( herm,  4UL );
         checkCapacity( herm, 11UL );
         checkNonZeros( herm,  6UL );
         checkNonZeros( herm, 0UL, 1UL );
         checkNonZeros( herm, 1UL, 1UL );
         checkNonZeros( herm, 2UL, 2UL );
         checkNonZeros( herm, 3UL, 2UL );

         if( herm(0,3) != cplx(3, 3) ||
             herm(1,1) != cplx(4, 0) ||
             herm(2,2) != cplx(6, 0) || herm(2,3) != cplx(7,7) ||
             herm(3,0) != cplx(3,-3) || herm(3,2) != cplx(7,-7) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (0, 0) (0,0) (0, 0) (3,3) )\n"
                                        "( (0, 0) (4,0) (0, 0) (0,0) )\n"
                                        "( (0, 0) (0,0) (6, 0) (7,7) )\n"
                                        "( (3,-3) (0,0) (7,-7) (0,0) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos->value() != cplx(3,3) || pos->index() != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: (3,3)\n"
                << "   Expected index: 3\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to erase a zero element
      {
         Iterator pos = herm.erase( 0UL, herm.find( 0UL, 1UL ) );

         checkRows    ( herm,  4UL );
         checkColumns ( herm,  4UL );
         checkCapacity( herm, 11UL );
         checkNonZeros( herm,  6UL );
         checkNonZeros( herm, 0UL, 1UL );
         checkNonZeros( herm, 1UL, 1UL );
         checkNonZeros( herm, 2UL, 2UL );
         checkNonZeros( herm, 3UL, 2UL );

         if( herm(0,3) != cplx(3, 3) ||
             herm(1,1) != cplx(4, 0) ||
             herm(2,2) != cplx(6, 0) || herm(2,3) != cplx(7,7) ||
             herm(3,0) != cplx(3,-3) || herm(3,2) != cplx(7,-7) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (0, 0) (0,0) (0, 0) (3,3) )\n"
                                        "( (0, 0) (4,0) (0, 0) (0,0) )\n"
                                        "( (0, 0) (0,0) (6, 0) (7,7) )\n"
                                        "( (3,-3) (0,0) (7,-7) (0,0) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos != herm.end( 0UL ) ) {
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
      test_ = "Row-major HermitianMatrix::erase( size_t, Iterator, Iterator )";

      using Iterator = HT::Iterator;

      // Initialization check
      HT herm( 4UL );
      herm(0,0) = cplx(1,0);
      herm(0,2) = cplx(2,2);
      herm(0,3) = cplx(3,3);
      herm(1,1) = cplx(4,0);
      herm(1,2) = cplx(5,5);
      herm(2,2) = cplx(6,0);
      herm(2,3) = cplx(7,7);

      checkRows    ( herm,  4UL );
      checkColumns ( herm,  4UL );
      checkCapacity( herm, 11UL );
      checkNonZeros( herm, 11UL );
      checkNonZeros( herm, 0UL, 3UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 4UL );
      checkNonZeros( herm, 3UL, 2UL );

      if( herm(0,0) != cplx(1, 0) || herm(0,2) != cplx(2, 2) || herm(0,3) != cplx(3,3) ||
          herm(1,1) != cplx(4, 0) || herm(1,2) != cplx(5, 5) ||
          herm(2,0) != cplx(2,-2) || herm(2,1) != cplx(5,-5) || herm(2,2) != cplx(6,0) || herm(2,3) != cplx(7,7) ||
          herm(3,0) != cplx(3,-3) || herm(3,2) != cplx(7,-7) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (1, 0) (0, 0) (2, 2) (3,3) )\n"
                                     "( (0, 0) (4, 0) (5, 5) (0,0) )\n"
                                     "( (2,-2) (5,-5) (6, 0) (7,7) )\n"
                                     "( (3,-3) (0, 0) (7,-7) (0,0) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element from (0,0) to (0,2)
      {
         Iterator pos = herm.erase( 0UL, herm.find( 0UL, 0UL ), herm.find( 0UL, 2UL ) );

         checkRows    ( herm,  4UL );
         checkColumns ( herm,  4UL );
         checkCapacity( herm, 11UL );
         checkNonZeros( herm, 10UL );
         checkNonZeros( herm, 0UL, 2UL );
         checkNonZeros( herm, 1UL, 2UL );
         checkNonZeros( herm, 2UL, 4UL );
         checkNonZeros( herm, 3UL, 2UL );

         if( herm(0,2) != cplx(2, 2) || herm(0,3) != cplx(3, 3) ||
             herm(1,1) != cplx(4, 0) || herm(1,2) != cplx(5, 5) ||
             herm(2,0) != cplx(2,-2) || herm(2,1) != cplx(5,-5) || herm(2,2) != cplx(6,0) || herm(2,3) != cplx(7,7) ||
             herm(3,0) != cplx(3,-3) || herm(3,2) != cplx(7,-7) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a single-element range failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (0, 0) (0, 0) (2, 2) (3,3) )\n"
                                        "( (0, 0) (4, 0) (5, 5) (0,0) )\n"
                                        "( (2,-2) (5,-5) (6, 0) (7,7) )\n"
                                        "( (3,-3) (0, 0) (7,-7) (0,0) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos->value() != cplx(2,2) || pos->index() != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: (2,2)\n"
                << "   Expected index: 2\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the element from (2,1) to (2,3)
      {
         Iterator pos = herm.erase( 2UL, herm.find( 2UL, 1UL ), herm.find( 2UL, 3UL ) );

         checkRows    ( herm,  4UL );
         checkColumns ( herm,  4UL );
         checkCapacity( herm, 11UL );
         checkNonZeros( herm,  7UL );
         checkNonZeros( herm, 0UL, 2UL );
         checkNonZeros( herm, 1UL, 1UL );
         checkNonZeros( herm, 2UL, 2UL );
         checkNonZeros( herm, 3UL, 2UL );

         if( herm(0,2) != cplx(2, 2) || herm(0,3) != cplx(3, 3) ||
             herm(1,1) != cplx(4, 0) ||
             herm(2,0) != cplx(2,-2) || herm(2,3) != cplx(7, 7) ||
             herm(3,0) != cplx(3,-3) || herm(3,2) != cplx(7,-7) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a single-element range failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (0, 0) (0,0) (2, 2) (3,3) )\n"
                                        "( (0, 0) (4,0) (0, 0) (0,0) )\n"
                                        "( (2,-2) (0,0) (0, 0) (7,7) )\n"
                                        "( (3,-3) (0,0) (7,-7) (0,0) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos->value() != cplx(7,7) || pos->index() != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: (7,7)\n"
                << "   Expected index: 3\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the element from (3,2) to the row end
      {
         Iterator pos = herm.erase( 3UL, herm.find( 3UL, 2UL ), herm.end( 3UL ) );

         checkRows    ( herm,  4UL );
         checkColumns ( herm,  4UL );
         checkCapacity( herm, 11UL );
         checkNonZeros( herm,  5UL );
         checkNonZeros( herm, 0UL, 2UL );
         checkNonZeros( herm, 1UL, 1UL );
         checkNonZeros( herm, 2UL, 1UL );
         checkNonZeros( herm, 3UL, 1UL );

         if( herm(0,2) != cplx(2, 2) || herm(0,3) != cplx(3, 3) ||
             herm(1,1) != cplx(4, 0) ||
             herm(2,0) != cplx(2,-2) ||
             herm(3,0) != cplx(3,-3) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a single-element range failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (0, 0) (0,0) (2,2) (3,3) )\n"
                                        "( (0, 0) (4,0) (0,0) (0,0) )\n"
                                        "( (2,-2) (0,0) (0,0) (0,0) )\n"
                                        "( (3,-3) (0,0) (0,0) (0,0) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos != herm.end( 3UL ) ) {
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
         Iterator pos = herm.erase( 2UL, herm.find( 2UL, 0UL ), herm.find( 2UL, 0UL ) );

         checkRows    ( herm,  4UL );
         checkColumns ( herm,  4UL );
         checkCapacity( herm, 11UL );
         checkNonZeros( herm,  5UL );
         checkNonZeros( herm, 0UL, 2UL );
         checkNonZeros( herm, 1UL, 1UL );
         checkNonZeros( herm, 2UL, 1UL );
         checkNonZeros( herm, 3UL, 1UL );

         if( herm(0,2) != cplx(2, 2) || herm(0,3) != cplx(3, 3) ||
             herm(1,1) != cplx(4, 0) ||
             herm(2,0) != cplx(2,-2) ||
             herm(3,0) != cplx(3,-3) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing an empty range failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (0, 0) (0,0) (2,2) (3,3) )\n"
                                        "( (0, 0) (4,0) (0,0) (0,0) )\n"
                                        "( (2,-2) (0,0) (0,0) (0,0) )\n"
                                        "( (3,-3) (0,0) (0,0) (0,0) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos != herm.find( 2UL, 0UL ) ) {
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
      test_ = "Row-major HermitianMatrix::erase( Predicate )";

      // Initialization check
      HT herm( 4UL );
      herm(0,0) = cplx(1,0);
      herm(0,2) = cplx(2,2);
      herm(0,3) = cplx(3,3);
      herm(1,1) = cplx(4,0);
      herm(1,2) = cplx(5,5);
      herm(2,2) = cplx(6,0);
      herm(2,3) = cplx(7,7);

      checkRows    ( herm,  4UL );
      checkColumns ( herm,  4UL );
      checkCapacity( herm, 11UL );
      checkNonZeros( herm, 11UL );
      checkNonZeros( herm, 0UL, 3UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 4UL );
      checkNonZeros( herm, 3UL, 2UL );

      if( herm(0,0) != cplx(1, 0) || herm(0,2) != cplx(2, 2) || herm(0,3) != cplx(3,3) ||
          herm(1,1) != cplx(4, 0) || herm(1,2) != cplx(5, 5) ||
          herm(2,0) != cplx(2,-2) || herm(2,1) != cplx(5,-5) || herm(2,2) != cplx(6,0) || herm(2,3) != cplx(7,7) ||
          herm(3,0) != cplx(3,-3) || herm(3,2) != cplx(7,-7) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (1, 0) (0, 0) (2, 2) (3,3) )\n"
                                     "( (0, 0) (4, 0) (5, 5) (0,0) )\n"
                                     "( (2,-2) (5,-5) (6, 0) (7,7) )\n"
                                     "( (3,-3) (0, 0) (7,-7) (0,0) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing a selection of elements
      herm.erase( []( const cplx& value ) {
                     return value == cplx(1,0) || value == cplx(5,5) || value == cplx(6,0);
                  } );

      checkRows    ( herm,  4UL );
      checkColumns ( herm,  4UL );
      checkCapacity( herm, 11UL );
      checkNonZeros( herm,  7UL );
      checkNonZeros( herm, 0UL, 2UL );
      checkNonZeros( herm, 1UL, 1UL );
      checkNonZeros( herm, 2UL, 2UL );
      checkNonZeros( herm, 3UL, 2UL );

      if( herm(0,2) != cplx(2, 2) || herm(0,3) != cplx(3, 3) ||
          herm(1,1) != cplx(4, 0) ||
          herm(2,0) != cplx(2,-2) || herm(2,3) != cplx(7, 7) ||
          herm(3,0) != cplx(3,-3) || herm(3,2) != cplx(7,-7) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a selection of elements failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (0, 0) (0, 0) (2, 2) (3,3) )\n"
                                     "( (0, 0) (4, 0) (0, 0) (0,0) )\n"
                                     "( (2,-2) (0, 0) (0, 0) (7,7) )\n"
                                     "( (3,-3) (0, 0) (7,-7) (0,0) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to erase all elements with value (1,0)
      herm.erase( []( const cplx& value ){ return value == cplx(1,0); } );

      checkRows    ( herm,  4UL );
      checkColumns ( herm,  4UL );
      checkCapacity( herm, 11UL );
      checkNonZeros( herm,  7UL );
      checkNonZeros( herm, 0UL, 2UL );
      checkNonZeros( herm, 1UL, 1UL );
      checkNonZeros( herm, 2UL, 2UL );
      checkNonZeros( herm, 3UL, 2UL );

      if( herm(0,2) != cplx(2, 2) || herm(0,3) != cplx(3, 3) ||
          herm(1,1) != cplx(4, 0) ||
          herm(2,0) != cplx(2,-2) || herm(2,3) != cplx(7, 7) ||
          herm(3,0) != cplx(3,-3) || herm(3,2) != cplx(7,-7) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing all elements with value (1,0) failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (0, 0) (0, 0) (2, 2) (3,3) )\n"
                                     "( (0, 0) (4, 0) (0, 0) (0,0) )\n"
                                     "( (2,-2) (0, 0) (0, 0) (7,7) )\n"
                                     "( (3,-3) (0, 0) (7,-7) (0,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major iterator-range-based erase function with predicate
   //=====================================================================================

   {
      test_ = "Row-major HermitianMatrix::erase( size_t, Iterator, Iterator, Predicate )";

      // Initialization check
      HT herm( 4UL );
      herm(0,0) = cplx(1,0);
      herm(0,2) = cplx(2,2);
      herm(0,3) = cplx(3,3);
      herm(1,1) = cplx(4,0);
      herm(1,2) = cplx(5,5);
      herm(2,2) = cplx(6,0);
      herm(2,3) = cplx(7,7);

      checkRows    ( herm,  4UL );
      checkColumns ( herm,  4UL );
      checkCapacity( herm, 11UL );
      checkNonZeros( herm, 11UL );
      checkNonZeros( herm, 0UL, 3UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 4UL );
      checkNonZeros( herm, 3UL, 2UL );

      if( herm(0,0) != cplx(1, 0) || herm(0,2) != cplx(2, 2) || herm(0,3) != cplx(3,3) ||
          herm(1,1) != cplx(4, 0) || herm(1,2) != cplx(5, 5) ||
          herm(2,0) != cplx(2,-2) || herm(2,1) != cplx(5,-5) || herm(2,2) != cplx(6,0) || herm(2,3) != cplx(7,7) ||
          herm(3,0) != cplx(3,-3) || herm(3,2) != cplx(7,-7) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (1, 0) (0, 0) (2, 2) (3,3) )\n"
                                     "( (0, 0) (4, 0) (5, 5) (0,0) )\n"
                                     "( (2,-2) (5,-5) (6, 0) (7,7) )\n"
                                     "( (3,-3) (0, 0) (7,-7) (0,0) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing a selection of elements
      herm.erase( 2UL, herm.begin( 2UL ), herm.find( 2UL, 3UL ),
                  []( const cplx& value ){ return value == cplx(2,-2) || value == cplx(6,0); } );

      checkRows    ( herm,  4UL );
      checkColumns ( herm,  4UL );
      checkCapacity( herm, 11UL );
      checkNonZeros( herm,  8UL );
      checkNonZeros( herm, 0UL, 2UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 2UL );
      checkNonZeros( herm, 3UL, 2UL );

      if( herm(0,0) != cplx(1, 0) || herm(0,3) != cplx(3, 3) ||
          herm(1,1) != cplx(4, 0) || herm(1,2) != cplx(5, 5) ||
          herm(2,1) != cplx(5,-5) || herm(2,3) != cplx(7, 7) ||
          herm(3,0) != cplx(3,-3) || herm(3,2) != cplx(7,-7) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a selection of elements failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (0, 0) (0, 0) (2, 2) (3,3) )\n"
                                     "( (0, 0) (4, 0) (0, 0) (0,0) )\n"
                                     "( (2,-2) (0, 0) (0, 0) (7,7) )\n"
                                     "( (3,-3) (0, 0) (7,-7) (0,0) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to erase from an empty range
      herm.erase( 1UL, herm.begin( 1UL ), herm.begin( 1UL ), []( const cplx& ){ return true; } );

      checkRows    ( herm,  4UL );
      checkColumns ( herm,  4UL );
      checkCapacity( herm, 11UL );
      checkNonZeros( herm,  8UL );
      checkNonZeros( herm, 0UL, 2UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 2UL );
      checkNonZeros( herm, 3UL, 2UL );

      if( herm(0,0) != cplx(1, 0) || herm(0,3) != cplx(3, 3) ||
          herm(1,1) != cplx(4, 0) || herm(1,2) != cplx(5, 5) ||
          herm(2,1) != cplx(5,-5) || herm(2,3) != cplx(7, 7) ||
          herm(3,0) != cplx(3,-3) || herm(3,2) != cplx(7,-7) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing from an empty range failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (0, 0) (0, 0) (2, 2) (3,3) )\n"
                                     "( (0, 0) (4, 0) (0, 0) (0,0) )\n"
                                     "( (2,-2) (0, 0) (0, 0) (7,7) )\n"
                                     "( (3,-3) (0, 0) (7,-7) (0,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major index-based erase function
   //=====================================================================================

   {
      test_ = "Column-major HermitianMatrix::erase( size_t, size_t )";

      // Initialization check
      OHT herm( 4UL );
      herm(0,0) = cplx(1,0);
      herm(0,2) = cplx(2,2);
      herm(0,3) = cplx(3,3);
      herm(1,1) = cplx(4,0);
      herm(1,2) = cplx(5,5);
      herm(2,2) = cplx(6,0);
      herm(2,3) = cplx(7,7);

      checkRows    ( herm,  4UL );
      checkColumns ( herm,  4UL );
      checkCapacity( herm, 11UL );
      checkNonZeros( herm, 11UL );
      checkNonZeros( herm, 0UL, 3UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 4UL );
      checkNonZeros( herm, 3UL, 2UL );

      if( herm(0,0) != cplx(1, 0) || herm(0,2) != cplx(2, 2) || herm(0,3) != cplx(3,3) ||
          herm(1,1) != cplx(4, 0) || herm(1,2) != cplx(5, 5) ||
          herm(2,0) != cplx(2,-2) || herm(2,1) != cplx(5,-5) || herm(2,2) != cplx(6,0) || herm(2,3) != cplx(7,7) ||
          herm(3,0) != cplx(3,-3) || herm(3,2) != cplx(7,-7) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (1, 0) (0, 0) (2, 2) (3,3) )\n"
                                     "( (0, 0) (4, 0) (5, 5) (0,0) )\n"
                                     "( (2,-2) (5,-5) (6, 0) (7,7) )\n"
                                     "( (3,-3) (0, 0) (7,-7) (0,0) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element at (0,0)
      herm.erase( 0UL, 0UL );

      checkRows    ( herm,  4UL );
      checkColumns ( herm,  4UL );
      checkCapacity( herm, 11UL );
      checkNonZeros( herm, 10UL );
      checkNonZeros( herm, 0UL, 2UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 4UL );
      checkNonZeros( herm, 3UL, 2UL );

      if( herm(0,2) != cplx(2, 2) || herm(0,3) != cplx(3, 3) ||
          herm(1,1) != cplx(4, 0) || herm(1,2) != cplx(5, 5) ||
          herm(2,0) != cplx(2,-2) || herm(2,1) != cplx(5,-5) || herm(2,2) != cplx(6,0) || herm(2,3) != cplx(7,7) ||
          herm(3,0) != cplx(3,-3) || herm(3,2) != cplx(7,-7) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (0, 0) (0, 0) (2, 2) (3,3) )\n"
                                     "( (0, 0) (4, 0) (5, 5) (0,0) )\n"
                                     "( (2,-2) (5,-5) (6, 0) (7,7) )\n"
                                     "( (3,-3) (0, 0) (7,-7) (0,0) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element at (2,1)
      herm.erase( 2UL, 1UL );

      checkRows    ( herm,  4UL );
      checkColumns ( herm,  4UL );
      checkCapacity( herm, 11UL );
      checkNonZeros( herm,  8UL );
      checkNonZeros( herm, 0UL, 2UL );
      checkNonZeros( herm, 1UL, 1UL );
      checkNonZeros( herm, 2UL, 3UL );
      checkNonZeros( herm, 3UL, 2UL );

      if( herm(0,2) != cplx(2, 2) || herm(0,3) != cplx(3, 3) ||
          herm(1,1) != cplx(4, 0) ||
          herm(2,0) != cplx(2,-2) || herm(2,2) != cplx(6, 0) || herm(2,3) != cplx(7,7) ||
          herm(3,0) != cplx(3,-3) || herm(3,2) != cplx(7,-7) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (0, 0) (0,0) (2, 2) (3,3) )\n"
                                     "( (0, 0) (4,0) (0, 0) (0,0) )\n"
                                     "( (2,-2) (0,0) (6, 0) (7,7) )\n"
                                     "( (3,-3) (0,0) (7,-7) (0,0) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element at (2,0)
      herm.erase( 2UL, 0UL );

      checkRows    ( herm,  4UL );
      checkColumns ( herm,  4UL );
      checkCapacity( herm, 11UL );
      checkNonZeros( herm,  6UL );
      checkNonZeros( herm, 0UL, 1UL );
      checkNonZeros( herm, 1UL, 1UL );
      checkNonZeros( herm, 2UL, 2UL );
      checkNonZeros( herm, 3UL, 2UL );

      if( herm(0,3) != cplx(3, 3) ||
          herm(1,1) != cplx(4, 0) ||
          herm(2,2) != cplx(6, 0) || herm(2,3) != cplx(7,7) ||
          herm(3,0) != cplx(3,-3) || herm(3,2) != cplx(7,-7) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (0, 0) (0,0) (0, 0) (3,3) )\n"
                                     "( (0, 0) (4,0) (0, 0) (0,0) )\n"
                                     "( (0, 0) (0,0) (6, 0) (7,7) )\n"
                                     "( (3,-3) (0,0) (7,-7) (0,0) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to erase a zero element
      herm.erase( 1UL, 0UL );

      checkRows    ( herm,  4UL );
      checkColumns ( herm,  4UL );
      checkCapacity( herm, 11UL );
      checkNonZeros( herm,  6UL );
      checkNonZeros( herm, 0UL, 1UL );
      checkNonZeros( herm, 1UL, 1UL );
      checkNonZeros( herm, 2UL, 2UL );
      checkNonZeros( herm, 3UL, 2UL );

      if( herm(0,3) != cplx(3, 3) ||
          herm(1,1) != cplx(4, 0) ||
          herm(2,2) != cplx(6, 0) || herm(2,3) != cplx(7,7) ||
          herm(3,0) != cplx(3,-3) || herm(3,2) != cplx(7,-7) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (0, 0) (0,0) (0, 0) (3,3) )\n"
                                     "( (0, 0) (4,0) (0, 0) (0,0) )\n"
                                     "( (0, 0) (0,0) (6, 0) (7,7) )\n"
                                     "( (3,-3) (0,0) (7,-7) (0,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major iterator-based erase function
   //=====================================================================================

   {
      test_ = "Column-major HermitianMatrix::erase( size_t, Iterator )";

      using Iterator = OHT::Iterator;

      // Initialization check
      OHT herm( 4UL );
      herm(0,0) = cplx(1,0);
      herm(0,2) = cplx(2,2);
      herm(0,3) = cplx(3,3);
      herm(1,1) = cplx(4,0);
      herm(1,2) = cplx(5,5);
      herm(2,2) = cplx(6,0);
      herm(2,3) = cplx(7,7);

      checkRows    ( herm,  4UL );
      checkColumns ( herm,  4UL );
      checkCapacity( herm, 11UL );
      checkNonZeros( herm, 11UL );
      checkNonZeros( herm, 0UL, 3UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 4UL );
      checkNonZeros( herm, 3UL, 2UL );

      if( herm(0,0) != cplx(1, 0) || herm(0,2) != cplx(2, 2) || herm(0,3) != cplx(3,3) ||
          herm(1,1) != cplx(4, 0) || herm(1,2) != cplx(5, 5) ||
          herm(2,0) != cplx(2,-2) || herm(2,1) != cplx(5,-5) || herm(2,2) != cplx(6,0) || herm(2,3) != cplx(7,7) ||
          herm(3,0) != cplx(3,-3) || herm(3,2) != cplx(7,-7) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (1, 0) (0, 0) (2, 2) (3,3) )\n"
                                     "( (0, 0) (4, 0) (5, 5) (0,0) )\n"
                                     "( (2,-2) (5,-5) (6, 0) (7,7) )\n"
                                     "( (3,-3) (0, 0) (7,-7) (0,0) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element at (0,0)
      {
         Iterator pos = herm.erase( 0UL, herm.find( 0UL, 0UL ) );

         checkRows    ( herm,  4UL );
         checkColumns ( herm,  4UL );
         checkCapacity( herm, 11UL );
         checkNonZeros( herm, 10UL );
         checkNonZeros( herm, 0UL, 2UL );
         checkNonZeros( herm, 1UL, 2UL );
         checkNonZeros( herm, 2UL, 4UL );
         checkNonZeros( herm, 3UL, 2UL );

         if( herm(0,2) != cplx(2, 2) || herm(0,3) != cplx(3, 3) ||
             herm(1,1) != cplx(4, 0) || herm(1,2) != cplx(5, 5) ||
             herm(2,0) != cplx(2,-2) || herm(2,1) != cplx(5,-5) || herm(2,2) != cplx(6,0) || herm(2,3) != cplx(7,7) ||
             herm(3,0) != cplx(3,-3) || herm(3,2) != cplx(7,-7) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (0, 0) (0, 0) (2, 2) (3,3) )\n"
                                        "( (0, 0) (4, 0) (5, 5) (0,0) )\n"
                                        "( (2,-2) (5,-5) (6, 0) (7,7) )\n"
                                        "( (3,-3) (0, 0) (7,-7) (0,0) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos->value() != cplx(2,-2) || pos->index() != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: (2,-2)\n"
                << "   Expected index: 2\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the element at (2,1)
      {
         Iterator pos = herm.erase( 1UL, herm.find( 2UL, 1UL ) );

         checkRows    ( herm,  4UL );
         checkColumns ( herm,  4UL );
         checkCapacity( herm, 11UL );
         checkNonZeros( herm,  8UL );
         checkNonZeros( herm, 0UL, 2UL );
         checkNonZeros( herm, 1UL, 1UL );
         checkNonZeros( herm, 2UL, 3UL );
         checkNonZeros( herm, 3UL, 2UL );

         if( herm(0,2) != cplx(2, 2) || herm(0,3) != cplx(3, 3) ||
             herm(1,1) != cplx(4, 0) ||
             herm(2,0) != cplx(2,-2) || herm(2,2) != cplx(6, 0) || herm(2,3) != cplx(7,7) ||
             herm(3,0) != cplx(3,-3) || herm(3,2) != cplx(7,-7) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (0, 0) (0,0) (2, 2) (3,3) )\n"
                                        "( (0, 0) (4,0) (0, 0) (0,0) )\n"
                                        "( (2,-2) (0,0) (6, 0) (7,7) )\n"
                                        "( (3,-3) (0,0) (7,-7) (0,0) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos != herm.end( 1UL ) ) {
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
         Iterator pos = herm.erase( 0UL, herm.find( 2UL, 0UL ) );

         checkRows    ( herm,  4UL );
         checkColumns ( herm,  4UL );
         checkCapacity( herm, 11UL );
         checkNonZeros( herm,  6UL );
         checkNonZeros( herm, 0UL, 1UL );
         checkNonZeros( herm, 1UL, 1UL );
         checkNonZeros( herm, 2UL, 2UL );
         checkNonZeros( herm, 3UL, 2UL );

         if( herm(0,3) != cplx(3, 3) ||
             herm(1,1) != cplx(4, 0) ||
             herm(2,2) != cplx(6, 0) || herm(2,3) != cplx(7,7) ||
             herm(3,0) != cplx(3,-3) || herm(3,2) != cplx(7,-7) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (0, 0) (0,0) (0, 0) (3,3) )\n"
                                        "( (0, 0) (4,0) (0, 0) (0,0) )\n"
                                        "( (0, 0) (0,0) (6, 0) (7,7) )\n"
                                        "( (3,-3) (0,0) (7,-7) (0,0) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos->value() != cplx(3,-3) || pos->index() != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: (3,-3)\n"
                << "   Expected index: 3\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to erase a zero element
      {
         Iterator pos = herm.erase( 0UL, herm.find( 1UL, 0UL ) );

         checkRows    ( herm,  4UL );
         checkColumns ( herm,  4UL );
         checkCapacity( herm, 11UL );
         checkNonZeros( herm,  6UL );
         checkNonZeros( herm, 0UL, 1UL );
         checkNonZeros( herm, 1UL, 1UL );
         checkNonZeros( herm, 2UL, 2UL );
         checkNonZeros( herm, 3UL, 2UL );

         if( herm(0,3) != cplx(3, 3) ||
             herm(1,1) != cplx(4, 0) ||
             herm(2,2) != cplx(6, 0) || herm(2,3) != cplx(7,7) ||
             herm(3,0) != cplx(3,-3) || herm(3,2) != cplx(7,-7) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (0, 0) (0,0) (0, 0) (3,3) )\n"
                                        "( (0, 0) (4,0) (0, 0) (0,0) )\n"
                                        "( (0, 0) (0,0) (6, 0) (7,7) )\n"
                                        "( (3,-3) (0,0) (7,-7) (0,0) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos != herm.end( 0UL ) ) {
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
      test_ = "Column-major HermitianMatrix::erase( size_t, Iterator, Iterator )";

      using Iterator = OHT::Iterator;

      // Initialization check
      OHT herm( 4UL );
      herm(0,0) = cplx(1,0);
      herm(0,2) = cplx(2,2);
      herm(0,3) = cplx(3,3);
      herm(1,1) = cplx(4,0);
      herm(1,2) = cplx(5,5);
      herm(2,2) = cplx(6,0);
      herm(2,3) = cplx(7,7);

      checkRows    ( herm,  4UL );
      checkColumns ( herm,  4UL );
      checkCapacity( herm, 11UL );
      checkNonZeros( herm, 11UL );
      checkNonZeros( herm, 0UL, 3UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 4UL );
      checkNonZeros( herm, 3UL, 2UL );

      if( herm(0,0) != cplx(1, 0) || herm(0,2) != cplx(2, 2) || herm(0,3) != cplx(3,3) ||
          herm(1,1) != cplx(4, 0) || herm(1,2) != cplx(5, 5) ||
          herm(2,0) != cplx(2,-2) || herm(2,1) != cplx(5,-5) || herm(2,2) != cplx(6,0) || herm(2,3) != cplx(7,7) ||
          herm(3,0) != cplx(3,-3) || herm(3,2) != cplx(7,-7) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (1, 0) (0, 0) (2, 2) (3,3) )\n"
                                     "( (0, 0) (4, 0) (5, 5) (0,0) )\n"
                                     "( (2,-2) (5,-5) (6, 0) (7,7) )\n"
                                     "( (3,-3) (0, 0) (7,-7) (0,0) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the element from (0,0) to (2,0)
      {
         Iterator pos = herm.erase( 0UL, herm.find( 0UL, 0UL ), herm.find( 2UL, 0UL ) );

         checkRows    ( herm,  4UL );
         checkColumns ( herm,  4UL );
         checkCapacity( herm, 11UL );
         checkNonZeros( herm, 10UL );
         checkNonZeros( herm, 0UL, 2UL );
         checkNonZeros( herm, 1UL, 2UL );
         checkNonZeros( herm, 2UL, 4UL );
         checkNonZeros( herm, 3UL, 2UL );

         if( herm(0,2) != cplx(2, 2) || herm(0,3) != cplx(3, 3) ||
             herm(1,1) != cplx(4, 0) || herm(1,2) != cplx(5, 5) ||
             herm(2,0) != cplx(2,-2) || herm(2,1) != cplx(5,-5) || herm(2,2) != cplx(6,0) || herm(2,3) != cplx(7,7) ||
             herm(3,0) != cplx(3,-3) || herm(3,2) != cplx(7,-7) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a single-element range failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (0, 0) (0, 0) (2, 2) (3,3) )\n"
                                        "( (0, 0) (4, 0) (5, 5) (0,0) )\n"
                                        "( (2,-2) (5,-5) (6, 0) (7,7) )\n"
                                        "( (3,-3) (0, 0) (7,-7) (0,0) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos->value() != cplx(2,-2) || pos->index() != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: (2,-2)\n"
                << "   Expected index: 2\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the element from (1,2) to (3,2)
      {
         Iterator pos = herm.erase( 2UL, herm.find( 1UL, 2UL ), herm.find( 3UL, 2UL ) );

         checkRows    ( herm,  4UL );
         checkColumns ( herm,  4UL );
         checkCapacity( herm, 11UL );
         checkNonZeros( herm,  7UL );
         checkNonZeros( herm, 0UL, 2UL );
         checkNonZeros( herm, 1UL, 1UL );
         checkNonZeros( herm, 2UL, 2UL );
         checkNonZeros( herm, 3UL, 2UL );

         if( herm(0,2) != cplx(2, 2) || herm(0,3) != cplx(3, 3) ||
             herm(1,1) != cplx(4, 0) ||
             herm(2,0) != cplx(2,-2) || herm(2,3) != cplx(7, 7) ||
             herm(3,0) != cplx(3,-3) || herm(3,2) != cplx(7,-7) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a single-element range failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (0, 0) (0,0) (2, 2) (3,3) )\n"
                                        "( (0, 0) (4,0) (0, 0) (0,0) )\n"
                                        "( (2,-2) (0,0) (0, 0) (7,7) )\n"
                                        "( (3,-3) (0,0) (7,-7) (0,0) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos->value() != cplx(7,-7) || pos->index() != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: (7,-7)\n"
                << "   Expected index: 3\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the element from (2,3) to the column end
      {
         Iterator pos = herm.erase( 3UL, herm.find( 2UL, 3UL ), herm.end( 3UL ) );

         checkRows    ( herm,  4UL );
         checkColumns ( herm,  4UL );
         checkCapacity( herm, 11UL );
         checkNonZeros( herm,  5UL );
         checkNonZeros( herm, 0UL, 2UL );
         checkNonZeros( herm, 1UL, 1UL );
         checkNonZeros( herm, 2UL, 1UL );
         checkNonZeros( herm, 3UL, 1UL );

         if( herm(0,2) != cplx(2, 2) || herm(0,3) != cplx(3, 3) ||
             herm(1,1) != cplx(4, 0) ||
             herm(2,0) != cplx(2,-2) ||
             herm(3,0) != cplx(3,-3) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a single-element range failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (0, 0) (0,0) (2,2) (3,3) )\n"
                                        "( (0, 0) (4,0) (0,0) (0,0) )\n"
                                        "( (2,-2) (0,0) (0,0) (0,0) )\n"
                                        "( (3,-3) (0,0) (0,0) (0,0) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos != herm.end( 3UL ) ) {
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
         Iterator pos = herm.erase( 2UL, herm.find( 0UL, 2UL ), herm.find( 0UL, 2UL ) );

         checkRows    ( herm,  4UL );
         checkColumns ( herm,  4UL );
         checkCapacity( herm, 11UL );
         checkNonZeros( herm,  5UL );
         checkNonZeros( herm, 0UL, 2UL );
         checkNonZeros( herm, 1UL, 1UL );
         checkNonZeros( herm, 2UL, 1UL );
         checkNonZeros( herm, 3UL, 1UL );

         if( herm(0,2) != cplx(2, 2) || herm(0,3) != cplx(3, 3) ||
             herm(1,1) != cplx(4, 0) ||
             herm(2,0) != cplx(2,-2) ||
             herm(3,0) != cplx(3,-3) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing an empty range failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (0, 0) (0,0) (2,2) (3,3) )\n"
                                        "( (0, 0) (4,0) (0,0) (0,0) )\n"
                                        "( (2,-2) (0,0) (0,0) (0,0) )\n"
                                        "( (3,-3) (0,0) (0,0) (0,0) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( pos != herm.find( 0UL, 2UL ) ) {
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
      test_ = "Column-major HermitianMatrix::erase( Predicate )";

      // Initialization check
      OHT herm( 4UL );
      herm(0,0) = cplx(1,0);
      herm(0,2) = cplx(2,2);
      herm(0,3) = cplx(3,3);
      herm(1,1) = cplx(4,0);
      herm(1,2) = cplx(5,5);
      herm(2,2) = cplx(6,0);
      herm(2,3) = cplx(7,7);

      checkRows    ( herm,  4UL );
      checkColumns ( herm,  4UL );
      checkCapacity( herm, 11UL );
      checkNonZeros( herm, 11UL );
      checkNonZeros( herm, 0UL, 3UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 4UL );
      checkNonZeros( herm, 3UL, 2UL );

      if( herm(0,0) != cplx(1, 0) || herm(0,2) != cplx(2, 2) || herm(0,3) != cplx(3,3) ||
          herm(1,1) != cplx(4, 0) || herm(1,2) != cplx(5, 5) ||
          herm(2,0) != cplx(2,-2) || herm(2,1) != cplx(5,-5) || herm(2,2) != cplx(6,0) || herm(2,3) != cplx(7,7) ||
          herm(3,0) != cplx(3,-3) || herm(3,2) != cplx(7,-7) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (1, 0) (0, 0) (2, 2) (3,3) )\n"
                                     "( (0, 0) (4, 0) (5, 5) (0,0) )\n"
                                     "( (2,-2) (5,-5) (6, 0) (7,7) )\n"
                                     "( (3,-3) (0, 0) (7,-7) (0,0) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing a selection of elements
      herm.erase( []( const cplx& value ) {
                     return value == cplx(1,0) || value == cplx(5,5) || value == cplx(6,0);
                  } );

      checkRows    ( herm,  4UL );
      checkColumns ( herm,  4UL );
      checkCapacity( herm, 11UL );
      checkNonZeros( herm,  7UL );
      checkNonZeros( herm, 0UL, 2UL );
      checkNonZeros( herm, 1UL, 1UL );
      checkNonZeros( herm, 2UL, 2UL );
      checkNonZeros( herm, 3UL, 2UL );

      if( herm(0,2) != cplx(2, 2) || herm(0,3) != cplx(3, 3) ||
          herm(1,1) != cplx(4, 0) ||
          herm(2,0) != cplx(2,-2) || herm(2,3) != cplx(7, 7) ||
          herm(3,0) != cplx(3,-3) || herm(3,2) != cplx(7,-7) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a selection of elements failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (0, 0) (0, 0) (2, 2) (3,3) )\n"
                                     "( (0, 0) (4, 0) (0, 0) (0,0) )\n"
                                     "( (2,-2) (0, 0) (0, 0) (7,7) )\n"
                                     "( (3,-3) (0, 0) (7,-7) (0,0) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to erase all elements with value (1,0)
      herm.erase( []( const cplx& value ){ return value == cplx(1,0); } );

      checkRows    ( herm,  4UL );
      checkColumns ( herm,  4UL );
      checkCapacity( herm, 11UL );
      checkNonZeros( herm,  7UL );
      checkNonZeros( herm, 0UL, 2UL );
      checkNonZeros( herm, 1UL, 1UL );
      checkNonZeros( herm, 2UL, 2UL );
      checkNonZeros( herm, 3UL, 2UL );

      if( herm(0,2) != cplx(2, 2) || herm(0,3) != cplx(3, 3) ||
          herm(1,1) != cplx(4, 0) ||
          herm(2,0) != cplx(2,-2) || herm(2,3) != cplx(7, 7) ||
          herm(3,0) != cplx(3,-3) || herm(3,2) != cplx(7,-7) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing all elements with value (1,0) failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (0, 0) (0, 0) (2, 2) (3,3) )\n"
                                     "( (0, 0) (4, 0) (0, 0) (0,0) )\n"
                                     "( (2,-2) (0, 0) (0, 0) (7,7) )\n"
                                     "( (3,-3) (0, 0) (7,-7) (0,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major iterator-range-based erase function with predicate
   //=====================================================================================

   {
      test_ = "Column-major HermitianMatrix::erase( size_t, Iterator, Iterator, Predicate )";

      // Initialization check
      OHT herm( 4UL );
      herm(0,0) = cplx(1,0);
      herm(0,2) = cplx(2,2);
      herm(0,3) = cplx(3,3);
      herm(1,1) = cplx(4,0);
      herm(1,2) = cplx(5,5);
      herm(2,2) = cplx(6,0);
      herm(2,3) = cplx(7,7);

      checkRows    ( herm,  4UL );
      checkColumns ( herm,  4UL );
      checkCapacity( herm, 11UL );
      checkNonZeros( herm, 11UL );
      checkNonZeros( herm, 0UL, 3UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 4UL );
      checkNonZeros( herm, 3UL, 2UL );

      if( herm(0,0) != cplx(1, 0) || herm(0,2) != cplx(2, 2) || herm(0,3) != cplx(3,3) ||
          herm(1,1) != cplx(4, 0) || herm(1,2) != cplx(5, 5) ||
          herm(2,0) != cplx(2,-2) || herm(2,1) != cplx(5,-5) || herm(2,2) != cplx(6,0) || herm(2,3) != cplx(7,7) ||
          herm(3,0) != cplx(3,-3) || herm(3,2) != cplx(7,-7) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (1, 0) (0, 0) (2, 2) (3,3) )\n"
                                     "( (0, 0) (4, 0) (5, 5) (0,0) )\n"
                                     "( (2,-2) (5,-5) (6, 0) (7,7) )\n"
                                     "( (3,-3) (0, 0) (7,-7) (0,0) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing a selection of elements
      herm.erase( 2UL, herm.begin( 2UL ), herm.find( 3UL, 2UL ),
                  []( const cplx& value ){ return value == cplx(2,2) || value == cplx(6,0); } );

      checkRows    ( herm,  4UL );
      checkColumns ( herm,  4UL );
      checkCapacity( herm, 11UL );
      checkNonZeros( herm,  8UL );
      checkNonZeros( herm, 0UL, 2UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 2UL );
      checkNonZeros( herm, 3UL, 2UL );

      if( herm(0,0) != cplx(1, 0) || herm(0,3) != cplx(3, 3) ||
          herm(1,1) != cplx(4, 0) || herm(1,2) != cplx(5, 5) ||
          herm(2,1) != cplx(5,-5) || herm(2,3) != cplx(7, 7) ||
          herm(3,0) != cplx(3,-3) || herm(3,2) != cplx(7,-7) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a selection of elements failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (0, 0) (0, 0) (2, 2) (3,3) )\n"
                                     "( (0, 0) (4, 0) (0, 0) (0,0) )\n"
                                     "( (2,-2) (0, 0) (0, 0) (7,7) )\n"
                                     "( (3,-3) (0, 0) (7,-7) (0,0) )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to erase from an empty range
      herm.erase( 1UL, herm.begin( 1UL ), herm.begin( 1UL ), []( const cplx& ){ return true; } );

      checkRows    ( herm,  4UL );
      checkColumns ( herm,  4UL );
      checkCapacity( herm, 11UL );
      checkNonZeros( herm,  8UL );
      checkNonZeros( herm, 0UL, 2UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 2UL );
      checkNonZeros( herm, 3UL, 2UL );

      if( herm(0,0) != cplx(1, 0) || herm(0,3) != cplx(3, 3) ||
          herm(1,1) != cplx(4, 0) || herm(1,2) != cplx(5, 5) ||
          herm(2,1) != cplx(5,-5) || herm(2,3) != cplx(7, 7) ||
          herm(3,0) != cplx(3,-3) || herm(3,2) != cplx(7,-7) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing from an empty range failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (0, 0) (0, 0) (2, 2) (3,3) )\n"
                                     "( (0, 0) (4, 0) (0, 0) (0,0) )\n"
                                     "( (2,-2) (0, 0) (0, 0) (7,7) )\n"
                                     "( (3,-3) (0, 0) (7,-7) (0,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c find() member function of the HermitianMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c find() member function of the HermitianMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseComplexTest::testFind()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major HermitianMatrix::find()";

      using ConstIterator = HT::ConstIterator;

      // Initialization check
      HT herm( 8UL, 3UL );
      herm(1,2) = cplx(1,1);
      herm(2,3) = cplx(2,2);
      herm(6,5) = cplx(3,3);

      checkRows    ( herm, 8UL );
      checkColumns ( herm, 8UL );
      checkCapacity( herm, 3UL );
      checkNonZeros( herm, 6UL );
      checkNonZeros( herm, 0UL, 0UL );
      checkNonZeros( herm, 1UL, 1UL );
      checkNonZeros( herm, 2UL, 2UL );
      checkNonZeros( herm, 3UL, 1UL );
      checkNonZeros( herm, 4UL, 0UL );
      checkNonZeros( herm, 5UL, 1UL );
      checkNonZeros( herm, 6UL, 1UL );
      checkNonZeros( herm, 7UL, 0UL );

      // Searching for the first element
      {
         ConstIterator pos( herm.find( 1UL, 2UL ) );

         if( pos == herm.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required position = (1,2)\n"
                << "   Current matrix:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 2 || pos->value() != cplx(1,1) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 2\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = (1,1)\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for the second element
      {
         ConstIterator pos( herm.find( 2UL, 3UL ) );

         if( pos == herm.end( 2UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required position = (2,3)\n"
                << "   Current matrix:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 3 || pos->value() != cplx(2,2) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 3\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = (2,2)\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for the third element
      {
         ConstIterator pos( herm.find( 6UL, 5UL ) );

         if( pos == herm.end( 6UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required position = (6,5)\n"
                << "   Current matrix:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 5 || pos->value() != cplx(3,3) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 5\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = (3,3)\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for a non-existing non-zero element
      {
         ConstIterator pos( herm.find( 4UL, 0UL ) );

         if( pos != herm.end( 4UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Non-existing element could be found\n"
                << " Details:\n"
                << "   Required index = 4\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 0\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major HermitianMatrix::find()";

      using ConstIterator = OHT::ConstIterator;

      // Initialization check
      OHT herm( 8UL, 3UL );
      herm(2,1) = cplx(1,1);
      herm(3,2) = cplx(2,2);
      herm(5,6) = cplx(3,3);

      checkRows    ( herm, 8UL );
      checkColumns ( herm, 8UL );
      checkCapacity( herm, 3UL );
      checkNonZeros( herm, 6UL );
      checkNonZeros( herm, 0UL, 0UL );
      checkNonZeros( herm, 1UL, 1UL );
      checkNonZeros( herm, 2UL, 2UL );
      checkNonZeros( herm, 3UL, 1UL );
      checkNonZeros( herm, 4UL, 0UL );
      checkNonZeros( herm, 5UL, 1UL );
      checkNonZeros( herm, 6UL, 1UL );
      checkNonZeros( herm, 7UL, 0UL );

      // Searching for the first element
      {
         ConstIterator pos( herm.find( 2UL, 1UL ) );

         if( pos == herm.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required position = (2,1)\n"
                << "   Current matrix:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 2 || pos->value() != cplx(1,1) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 2\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = (1,1)\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for the second element
      {
         ConstIterator pos( herm.find( 3UL, 2UL ) );

         if( pos == herm.end( 2UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required position = (3,2)\n"
                << "   Current matrix:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 3 || pos->value() != cplx(2,2) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 3\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = (2,2)\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for the third element
      {
         ConstIterator pos( herm.find( 5UL, 6UL ) );

         if( pos == herm.end( 6UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Element could not be found\n"
                << " Details:\n"
                << "   Required position = (5,6)\n"
                << "   Current matrix:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 5 || pos->value() != cplx(3,3) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 5\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = (3,3)\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Searching for a non-existing non-zero element
      {
         ConstIterator pos( herm.find( 0UL, 4UL ) );

         if( pos != herm.end( 4UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Non-existing element could be found\n"
                << " Details:\n"
                << "   Required index = 4\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = 0\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c lowerBound() member function of the HermitianMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c lowerBound() member function of the HermitianMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseComplexTest::testLowerBound()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major HermitianMatrix::lowerBound()";

      using ConstIterator = HT::ConstIterator;

      // Initialization check
      HT herm( 6UL, 3UL );
      herm(1,2) = cplx(1,1);
      herm(1,4) = cplx(2,2);

      checkRows    ( herm, 6UL );
      checkColumns ( herm, 6UL );
      checkCapacity( herm, 4UL );
      checkNonZeros( herm, 4UL );
      checkNonZeros( herm, 0UL, 0UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 1UL );
      checkNonZeros( herm, 3UL, 0UL );
      checkNonZeros( herm, 4UL, 1UL );
      checkNonZeros( herm, 5UL, 0UL );

      // Determining the lower bound for position (1,1)
      {
         ConstIterator pos( herm.lowerBound( 1UL, 1UL ) );

         if( pos == herm.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,1)\n"
                << "   Current matrix:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 2 || pos->value() != cplx(1,1) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 2\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = (1,1)\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for position (1,2)
      {
         ConstIterator pos( herm.lowerBound( 1UL, 2UL ) );

         if( pos == herm.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,2)\n"
                << "   Current matrix:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 2 || pos->value() != cplx(1,1) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 2\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = (1,1)\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for position (1,3)
      {
         ConstIterator pos( herm.lowerBound( 1UL, 3UL ) );

         if( pos == herm.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,3)\n"
                << "   Current matrix:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 4 || pos->value() != cplx(2,2) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 4\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = cplx(2,2)\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for position (1,4)
      {
         ConstIterator pos( herm.lowerBound( 1UL, 4UL ) );

         if( pos == herm.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,4)\n"
                << "   Current matrix:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 4 || pos->value() != cplx(2,2) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 4\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = cplx(2,2)\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for position (1,5)
      {
         ConstIterator pos( herm.lowerBound( 1UL, 5UL ) );

         if( pos != herm.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,5)\n"
                << "   Current matrix:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major HermitianMatrix::lowerBound()";

      using ConstIterator = OHT::ConstIterator;

      // Initialization check
      OHT herm( 6UL, 3UL );
      herm(2,1) = cplx(1,1);
      herm(4,1) = cplx(2,2);

      checkRows    ( herm, 6UL );
      checkColumns ( herm, 6UL );
      checkCapacity( herm, 4UL );
      checkNonZeros( herm, 4UL );
      checkNonZeros( herm, 0UL, 0UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 1UL );
      checkNonZeros( herm, 3UL, 0UL );
      checkNonZeros( herm, 4UL, 1UL );
      checkNonZeros( herm, 5UL, 0UL );

      // Determining the lower bound for position (1,1)
      {
         ConstIterator pos( herm.lowerBound( 1UL, 1UL ) );

         if( pos == herm.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,1)\n"
                << "   Current matrix:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 2 || pos->value() != cplx(1,1) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 2\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = (1,1)\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for position (2,1)
      {
         ConstIterator pos( herm.lowerBound( 2UL, 1UL ) );

         if( pos == herm.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (2,1)\n"
                << "   Current matrix:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 2 || pos->value() != cplx(1,1) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 2\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = (1,1)\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for position (3,1)
      {
         ConstIterator pos( herm.lowerBound( 3UL, 1UL ) );

         if( pos == herm.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (3,1)\n"
                << "   Current matrix:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 4 || pos->value() != cplx(2,2) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 4\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = (2,2)\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for position (4,1)
      {
         ConstIterator pos( herm.lowerBound( 4UL, 1UL ) );

         if( pos == herm.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (4,1)\n"
                << "   Current matrix:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 4 || pos->value() != cplx(2,2) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 4\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = (2,2)\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the lower bound for position (5,1)
      {
         ConstIterator pos( herm.lowerBound( 5UL, 1UL ) );

         if( pos != herm.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Lower bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (5,1)\n"
                << "   Current matrix:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c upperBound() member function of the HermitianMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c upperBound() member function of the HermitianMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseComplexTest::testUpperBound()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major HermitianMatrix::upperBound()";

      using ConstIterator = HT::ConstIterator;

      // Initialization check
      HT herm( 6UL, 3UL );
      herm(1,2) = cplx(1,1);
      herm(1,4) = cplx(2,2);

      checkRows    ( herm, 6UL );
      checkColumns ( herm, 6UL );
      checkCapacity( herm, 4UL );
      checkNonZeros( herm, 4UL );
      checkNonZeros( herm, 0UL, 0UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 1UL );
      checkNonZeros( herm, 3UL, 0UL );
      checkNonZeros( herm, 4UL, 1UL );
      checkNonZeros( herm, 5UL, 0UL );

      // Determining the upper bound for position (1,1)
      {
         ConstIterator pos( herm.upperBound( 1UL, 1UL ) );

         if( pos == herm.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,1)\n"
                << "   Current matrix:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 2 || pos->value() != cplx(1,1) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 2\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = (1,1)\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for position (1,2)
      {
         ConstIterator pos( herm.upperBound( 1UL, 2UL ) );

         if( pos == herm.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,2)\n"
                << "   Current matrix:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 4 || pos->value() != cplx(2,2) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 4\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = (2,2)\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for position (1,3)
      {
         ConstIterator pos( herm.upperBound( 1UL, 3UL ) );

         if( pos == herm.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,3)\n"
                << "   Current matrix:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 4 || pos->value() != cplx(2,2) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 4\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = (2,2)\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for position (1,4)
      {
         ConstIterator pos( herm.upperBound( 1UL, 4UL ) );

         if( pos != herm.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,4)\n"
                << "   Current matrix:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for position (1,5)
      {
         ConstIterator pos( herm.upperBound( 1UL, 5UL ) );

         if( pos != herm.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,5)\n"
                << "   Current matrix:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major HermitianMatrix::upperBound()";

      using ConstIterator = OHT::ConstIterator;

      // Initialization check
      OHT herm( 6UL, 3UL );
      herm(2,1) = cplx(1,1);
      herm(4,1) = cplx(2,2);

      checkRows    ( herm, 6UL );
      checkColumns ( herm, 6UL );
      checkCapacity( herm, 4UL );
      checkNonZeros( herm, 4UL );
      checkNonZeros( herm, 0UL, 0UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 1UL );
      checkNonZeros( herm, 3UL, 0UL );
      checkNonZeros( herm, 4UL, 1UL );
      checkNonZeros( herm, 5UL, 0UL );

      // Determining the upper bound for position (1,1)
      {
         ConstIterator pos( herm.upperBound( 1UL, 1UL ) );

         if( pos == herm.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (1,1)\n"
                << "   Current matrix:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 2 || pos->value() != cplx(1,1) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 2\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = (1,1)\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for position (2,1)
      {
         ConstIterator pos( herm.upperBound( 2UL, 1UL ) );

         if( pos == herm.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (2,1)\n"
                << "   Current matrix:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 4 || pos->value() != cplx(2,2) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 4\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = (2,2)\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for position (3,1)
      {
         ConstIterator pos( herm.upperBound( 3UL, 1UL ) );

         if( pos == herm.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (3,1)\n"
                << "   Current matrix:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         else if( pos->index() != 4 || pos->value() != cplx(2,2) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Wrong element found\n"
                << " Details:\n"
                << "   Required index = 4\n"
                << "   Found index    = " << pos->index() << "\n"
                << "   Expected value = (2,2)\n"
                << "   Value at index = " << pos->value() << "\n"
                << "   Current matrix:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for position (4,1)
      {
         ConstIterator pos( herm.upperBound( 4UL, 1UL ) );

         if( pos != herm.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (4,1)\n"
                << "   Current matrix:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Determining the upper bound for position (5,1)
      {
         ConstIterator pos( herm.upperBound( 5UL, 1UL ) );

         if( pos != herm.end( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Upper bound could not be determined\n"
                << " Details:\n"
                << "   Required position = (5,1)\n"
                << "   Current matrix:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c transpose() member function of the HermitianMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c transpose() member function of the HermitianMatrix
// specialization. Additionally, it performs a test of self-transpose via the \c trans()
// function. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseComplexTest::testTranspose()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major self-transpose via transpose()";

      HT herm( 4UL );
      herm(0,0) = cplx(1, 0);
      herm(0,2) = cplx(2,-1);
      herm(0,3) = cplx(3, 2);
      herm(1,1) = cplx(4, 0);
      herm(1,3) = cplx(5,-3);
      herm(2,2) = cplx(6, 0);
      herm(2,3) = cplx(7, 1);

      transpose( herm );

      checkRows    ( herm,  4UL );
      checkColumns ( herm,  4UL );
      checkCapacity( herm, 11UL );
      checkNonZeros( herm, 11UL );
      checkNonZeros( herm,  0UL, 3UL );
      checkNonZeros( herm,  1UL, 2UL );
      checkNonZeros( herm,  2UL, 3UL );
      checkNonZeros( herm,  3UL, 3UL );

      if( herm(0,0) != cplx(1, 0) || herm(0,1) != cplx(0, 0) || herm(0,2) != cplx(2, 1) || herm(0,3) != cplx(3,-2) ||
          herm(1,0) != cplx(0, 0) || herm(1,1) != cplx(4, 0) || herm(1,2) != cplx(0, 0) || herm(1,3) != cplx(5, 3) ||
          herm(2,0) != cplx(2,-1) || herm(2,1) != cplx(0, 0) || herm(2,2) != cplx(6, 0) || herm(2,3) != cplx(7,-1) ||
          herm(3,0) != cplx(3, 2) || herm(3,1) != cplx(5,-3) || herm(3,2) != cplx(7, 1) || herm(3,3) != cplx(0, 0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Transpose operation failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (1, 0) (0, 0) (2,1) (3,-2) )\n"
                                     "( (0, 0) (4, 0) (0,0) (5, 3) )\n"
                                     "( (2,-1) (0, 0) (6,0) (7,-1) )\n"
                                     "( (3, 2) (5,-3) (7,1) (0, 0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major self-transpose via trans()";

      HT herm( 4UL );
      herm(0,0) = cplx(1, 0);
      herm(0,2) = cplx(2,-1);
      herm(0,3) = cplx(3, 2);
      herm(1,1) = cplx(4, 0);
      herm(1,3) = cplx(5,-3);
      herm(2,2) = cplx(6, 0);
      herm(2,3) = cplx(7, 1);

      herm = trans( herm );

      checkRows    ( herm,  4UL );
      checkColumns ( herm,  4UL );
      checkCapacity( herm, 11UL );
      checkNonZeros( herm, 11UL );
      checkNonZeros( herm,  0UL, 3UL );
      checkNonZeros( herm,  1UL, 2UL );
      checkNonZeros( herm,  2UL, 3UL );
      checkNonZeros( herm,  3UL, 3UL );

      if( herm(0,0) != cplx(1, 0) || herm(0,1) != cplx(0, 0) || herm(0,2) != cplx(2, 1) || herm(0,3) != cplx(3,-2) ||
          herm(1,0) != cplx(0, 0) || herm(1,1) != cplx(4, 0) || herm(1,2) != cplx(0, 0) || herm(1,3) != cplx(5, 3) ||
          herm(2,0) != cplx(2,-1) || herm(2,1) != cplx(0, 0) || herm(2,2) != cplx(6, 0) || herm(2,3) != cplx(7,-1) ||
          herm(3,0) != cplx(3, 2) || herm(3,1) != cplx(5,-3) || herm(3,2) != cplx(7, 1) || herm(3,3) != cplx(0, 0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Transpose operation failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (1, 0) (0, 0) (2,1) (3,-2) )\n"
                                     "( (0, 0) (4, 0) (0,0) (5, 3) )\n"
                                     "( (2,-1) (0, 0) (6,0) (7,-1) )\n"
                                     "( (3, 2) (5,-3) (7,1) (0, 0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major self-transpose via transpose()";

      OHT herm( 4UL );
      herm(0,0) = cplx(1, 0);
      herm(0,2) = cplx(2,-1);
      herm(0,3) = cplx(3, 2);
      herm(1,1) = cplx(4, 0);
      herm(1,3) = cplx(5,-3);
      herm(2,2) = cplx(6, 0);
      herm(2,3) = cplx(7, 1);

      transpose( herm );

      checkRows    ( herm,  4UL );
      checkColumns ( herm,  4UL );
      checkCapacity( herm, 11UL );
      checkNonZeros( herm, 11UL );
      checkNonZeros( herm,  0UL, 3UL );
      checkNonZeros( herm,  1UL, 2UL );
      checkNonZeros( herm,  2UL, 3UL );
      checkNonZeros( herm,  3UL, 3UL );

      if( herm(0,0) != cplx(1, 0) || herm(0,1) != cplx(0, 0) || herm(0,2) != cplx(2, 1) || herm(0,3) != cplx(3,-2) ||
          herm(1,0) != cplx(0, 0) || herm(1,1) != cplx(4, 0) || herm(1,2) != cplx(0, 0) || herm(1,3) != cplx(5, 3) ||
          herm(2,0) != cplx(2,-1) || herm(2,1) != cplx(0, 0) || herm(2,2) != cplx(6, 0) || herm(2,3) != cplx(7,-1) ||
          herm(3,0) != cplx(3, 2) || herm(3,1) != cplx(5,-3) || herm(3,2) != cplx(7, 1) || herm(3,3) != cplx(0, 0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Transpose operation failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (1, 0) (0, 0) (2,1) (3,-2) )\n"
                                     "( (0, 0) (4, 0) (0,0) (5, 3) )\n"
                                     "( (2,-1) (0, 0) (6,0) (7,-1) )\n"
                                     "( (3, 2) (5,-3) (7,1) (0, 0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major self-transpose via trans()";

      OHT herm( 4UL );
      herm(0,0) = cplx(1, 0);
      herm(0,2) = cplx(2,-1);
      herm(0,3) = cplx(3, 2);
      herm(1,1) = cplx(4, 0);
      herm(1,3) = cplx(5,-3);
      herm(2,2) = cplx(6, 0);
      herm(2,3) = cplx(7, 1);

      herm = trans( herm );

      checkRows    ( herm,  4UL );
      checkColumns ( herm,  4UL );
      checkCapacity( herm, 11UL );
      checkNonZeros( herm, 11UL );
      checkNonZeros( herm,  0UL, 3UL );
      checkNonZeros( herm,  1UL, 2UL );
      checkNonZeros( herm,  2UL, 3UL );
      checkNonZeros( herm,  3UL, 3UL );

      if( herm(0,0) != cplx(1, 0) || herm(0,1) != cplx(0, 0) || herm(0,2) != cplx(2, 1) || herm(0,3) != cplx(3,-2) ||
          herm(1,0) != cplx(0, 0) || herm(1,1) != cplx(4, 0) || herm(1,2) != cplx(0, 0) || herm(1,3) != cplx(5, 3) ||
          herm(2,0) != cplx(2,-1) || herm(2,1) != cplx(0, 0) || herm(2,2) != cplx(6, 0) || herm(2,3) != cplx(7,-1) ||
          herm(3,0) != cplx(3, 2) || herm(3,1) != cplx(5,-3) || herm(3,2) != cplx(7, 1) || herm(3,3) != cplx(0, 0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Transpose operation failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (1, 0) (0, 0) (2,1) (3,-2) )\n"
                                     "( (0, 0) (4, 0) (0,0) (5, 3) )\n"
                                     "( (2,-1) (0, 0) (6,0) (7,-1) )\n"
                                     "( (3, 2) (5,-3) (7,1) (0, 0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c ctranspose() member function of the HermitianMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c ctranspose() member function of the HermitianMatrix
// specialization. Additionally, it performs a test of self-transpose via the \c ctrans()
// function. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseComplexTest::testCTranspose()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major self-transpose via ctranspose()";

      HT herm( 4UL );
      herm(0,0) = cplx(1, 0);
      herm(0,2) = cplx(2,-1);
      herm(0,3) = cplx(3, 2);
      herm(1,1) = cplx(4, 0);
      herm(1,3) = cplx(5,-3);
      herm(2,2) = cplx(6, 0);
      herm(2,3) = cplx(7, 1);

      ctranspose( herm );

      checkRows    ( herm,  4UL );
      checkColumns ( herm,  4UL );
      checkCapacity( herm, 11UL );
      checkNonZeros( herm, 11UL );
      checkNonZeros( herm,  0UL, 3UL );
      checkNonZeros( herm,  1UL, 2UL );
      checkNonZeros( herm,  2UL, 3UL );
      checkNonZeros( herm,  3UL, 3UL );

      if( herm(0,0) != cplx(1, 0) || herm(0,1) != cplx(0, 0) || herm(0,2) != cplx(2,-1) || herm(0,3) != cplx(3, 2) ||
          herm(1,0) != cplx(0, 0) || herm(1,1) != cplx(4, 0) || herm(1,2) != cplx(0, 0) || herm(1,3) != cplx(5,-3) ||
          herm(2,0) != cplx(2, 1) || herm(2,1) != cplx(0, 0) || herm(2,2) != cplx(6, 0) || herm(2,3) != cplx(7, 1) ||
          herm(3,0) != cplx(3,-2) || herm(3,1) != cplx(5, 3) || herm(3,2) != cplx(7,-1) || herm(3,3) != cplx(0, 0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Transpose operation failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (1, 0) (0, 0) (2,-1) (3, 2) )\n"
                                     "( (0, 0) (4, 0) (0, 0) (5,-3) )\n"
                                     "( (2, 1) (0, 0) (6, 0) (7, 1) )\n"
                                     "( (3,-2) (5, 3) (7,-1) (0, 0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major self-transpose via ctrans()";

      HT herm( 4UL );
      herm(0,0) = cplx(1, 0);
      herm(0,2) = cplx(2,-1);
      herm(0,3) = cplx(3, 2);
      herm(1,1) = cplx(4, 0);
      herm(1,3) = cplx(5,-3);
      herm(2,2) = cplx(6, 0);
      herm(2,3) = cplx(7, 1);

      herm = ctrans( herm );

      checkRows    ( herm,  4UL );
      checkColumns ( herm,  4UL );
      checkCapacity( herm, 11UL );
      checkNonZeros( herm, 11UL );
      checkNonZeros( herm,  0UL, 3UL );
      checkNonZeros( herm,  1UL, 2UL );
      checkNonZeros( herm,  2UL, 3UL );
      checkNonZeros( herm,  3UL, 3UL );

      if( herm(0,0) != cplx(1, 0) || herm(0,1) != cplx(0, 0) || herm(0,2) != cplx(2,-1) || herm(0,3) != cplx(3, 2) ||
          herm(1,0) != cplx(0, 0) || herm(1,1) != cplx(4, 0) || herm(1,2) != cplx(0, 0) || herm(1,3) != cplx(5,-3) ||
          herm(2,0) != cplx(2, 1) || herm(2,1) != cplx(0, 0) || herm(2,2) != cplx(6, 0) || herm(2,3) != cplx(7, 1) ||
          herm(3,0) != cplx(3,-2) || herm(3,1) != cplx(5, 3) || herm(3,2) != cplx(7,-1) || herm(3,3) != cplx(0, 0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Transpose operation failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (1, 0) (0, 0) (2,-1) (3, 2) )\n"
                                     "( (0, 0) (4, 0) (0, 0) (5,-3) )\n"
                                     "( (2, 1) (0, 0) (6, 0) (7, 1) )\n"
                                     "( (3,-2) (5, 3) (7,-1) (0, 0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major self-transpose via ctranspose()";

      OHT herm( 4UL );
      herm(0,0) = cplx(1, 0);
      herm(0,2) = cplx(2,-1);
      herm(0,3) = cplx(3, 2);
      herm(1,1) = cplx(4, 0);
      herm(1,3) = cplx(5,-3);
      herm(2,2) = cplx(6, 0);
      herm(2,3) = cplx(7, 1);

      ctranspose( herm );

      checkRows    ( herm,  4UL );
      checkColumns ( herm,  4UL );
      checkCapacity( herm, 11UL );
      checkNonZeros( herm, 11UL );
      checkNonZeros( herm,  0UL, 3UL );
      checkNonZeros( herm,  1UL, 2UL );
      checkNonZeros( herm,  2UL, 3UL );
      checkNonZeros( herm,  3UL, 3UL );

      if( herm(0,0) != cplx(1, 0) || herm(0,1) != cplx(0, 0) || herm(0,2) != cplx(2,-1) || herm(0,3) != cplx(3, 2) ||
          herm(1,0) != cplx(0, 0) || herm(1,1) != cplx(4, 0) || herm(1,2) != cplx(0, 0) || herm(1,3) != cplx(5,-3) ||
          herm(2,0) != cplx(2, 1) || herm(2,1) != cplx(0, 0) || herm(2,2) != cplx(6, 0) || herm(2,3) != cplx(7, 1) ||
          herm(3,0) != cplx(3,-2) || herm(3,1) != cplx(5, 3) || herm(3,2) != cplx(7,-1) || herm(3,3) != cplx(0, 0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Transpose operation failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (1, 0) (0, 0) (2,-1) (3, 2) )\n"
                                     "( (0, 0) (4, 0) (0, 0) (5,-3) )\n"
                                     "( (2, 1) (0, 0) (6, 0) (7, 1) )\n"
                                     "( (3,-2) (5, 3) (7,-1) (0, 0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major self-transpose via ctrans()";

      OHT herm( 4UL );
      herm(0,0) = cplx(1, 0);
      herm(0,2) = cplx(2,-1);
      herm(0,3) = cplx(3, 2);
      herm(1,1) = cplx(4, 0);
      herm(1,3) = cplx(5,-3);
      herm(2,2) = cplx(6, 0);
      herm(2,3) = cplx(7, 1);

      herm = ctrans( herm );

      checkRows    ( herm,  4UL );
      checkColumns ( herm,  4UL );
      checkCapacity( herm, 11UL );
      checkNonZeros( herm, 11UL );
      checkNonZeros( herm,  0UL, 3UL );
      checkNonZeros( herm,  1UL, 2UL );
      checkNonZeros( herm,  2UL, 3UL );
      checkNonZeros( herm,  3UL, 3UL );

      if( herm(0,0) != cplx(1, 0) || herm(0,1) != cplx(0, 0) || herm(0,2) != cplx(2,-1) || herm(0,3) != cplx(3, 2) ||
          herm(1,0) != cplx(0, 0) || herm(1,1) != cplx(4, 0) || herm(1,2) != cplx(0, 0) || herm(1,3) != cplx(5,-3) ||
          herm(2,0) != cplx(2, 1) || herm(2,1) != cplx(0, 0) || herm(2,2) != cplx(6, 0) || herm(2,3) != cplx(7, 1) ||
          herm(3,0) != cplx(3,-2) || herm(3,1) != cplx(5, 3) || herm(3,2) != cplx(7,-1) || herm(3,3) != cplx(0, 0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Transpose operation failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (1, 0) (0, 0) (2,-1) (3, 2) )\n"
                                     "( (0, 0) (4, 0) (0, 0) (5,-3) )\n"
                                     "( (2, 1) (0, 0) (6, 0) (7, 1) )\n"
                                     "( (3,-2) (5, 3) (7,-1) (0, 0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c isDefault() function with the HermitianMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isDefault() function with the HermitianMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseComplexTest::testIsDefault()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major isDefault() function";

      // isDefault with 0x0 matrix
      {
         HT herm;

         if( isDefault( herm ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with default matrix
      {
         HT herm( 3UL );

         if( isDefault( herm(0,1) ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix element:\n" << herm(0,1) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( herm ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with non-default matrix
      {
         HT herm( 3UL );
         herm(0,1) = cplx(1,1);

         if( isDefault( herm(0,1) ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix element:\n" << herm(0,1) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( herm ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << herm << "\n";
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
         OHT herm;

         if( isDefault( herm ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with default matrix
      {
         OHT herm( 3UL );

         if( isDefault( herm(0,1) ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix element:\n" << herm(0,1) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( herm ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with non-default matrix
      {
         OHT herm( 3UL );
         herm(0,1) = cplx(1,1);

         if( isDefault( herm(0,1) ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix element:\n" << herm(0,1) << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( isDefault( herm ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c submatrix() function with the HermitianMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c submatrix() function with the HermitianMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseComplexTest::testSubmatrix()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major submatrix() function";

      using SMT = blaze::Submatrix<HT>;

      HT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      SMT sm = submatrix( herm, 0UL, 1UL, 2UL, 2UL );

      if( sm(0,1) != cplx(7,3) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator access failed\n"
             << " Details:\n"
             << "   Result: " << sm(0,1) << "\n"
             << "   Expected result: (7,3)\n";
         throw std::runtime_error( oss.str() );
      }

      SMT::Iterator it = sm.begin(0UL);

      if( it == sm.end(0UL) || it->value() != cplx(-4,-1) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << it->value() << "\n"
             << "   Expected result: (-4,-1)\n";
         throw std::runtime_error( oss.str() );
      }

      sm(1,1) = cplx(-5,2);

      if( sm(0,0) != cplx(-4,-1) || sm(0,1) != cplx( 7,3) ||
          sm(1,0) != cplx( 2, 0) || sm(1,1) != cplx(-5,2) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix access failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( (-4,-1) ( 7,3) )\n"
                                     "( ( 2, 0) (-5,2) )\n";
         throw std::runtime_error( oss.str() );
      }

      if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-4,-1) || herm(0,2) != cplx( 7,3) ||
          herm(1,0) != cplx(-4, 1) || herm(1,1) != cplx( 2, 0) || herm(1,2) != cplx(-5,2) ||
          herm(2,0) != cplx( 7,-3) || herm(2,1) != cplx(-5,-2) || herm(2,2) != cplx( 3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix access failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( ( 1, 0) (-4,-1) ( 7,3) )\n"
                                     "( (-4, 1) ( 2, 0) (-5,2) )\n"
                                     "( ( 7,-3) (-5,-2) ( 3,0) )\n";
         throw std::runtime_error( oss.str() );
      }

      reset( sm );

      if( sm(0,0) != cplx(0,0) || sm(0,1) != cplx(0,0) ||
          sm(1,0) != cplx(0,0) || sm(1,1) != cplx(0,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix reset failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( (0,0) (0,0) )\n"
                                     "( (0,0) (0,0) )\n";
         throw std::runtime_error( oss.str() );
      }

      if( herm(0,0) != cplx(1,0) || herm(0,1) != cplx(0,0) || herm(0,2) != cplx(0,0) ||
          herm(1,0) != cplx(0,0) || herm(1,1) != cplx(0,0) || herm(1,2) != cplx(0,0) ||
          herm(2,0) != cplx(0,0) || herm(2,1) != cplx(0,0) || herm(2,2) != cplx(3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix reset failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (1,0) (0,0) (0,0) )\n"
                                     "( (0,0) (0,0) (0,0) )\n"
                                     "( (0,0) (0,0) (3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major submatrix() function";

      using SMT = blaze::Submatrix<OHT>;

      OHT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      SMT sm = submatrix( herm, 0UL, 1UL, 2UL, 2UL );

      if( sm(0,1) != cplx(7,3) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator access failed\n"
             << " Details:\n"
             << "   Result: " << sm(0,1) << "\n"
             << "   Expected result: (7,3)\n";
         throw std::runtime_error( oss.str() );
      }

      SMT::Iterator it = sm.begin(0UL);

      if( it == sm.end(0UL) || it->value() != cplx(-4,-1) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << it->value() << "\n"
             << "   Expected result: (-4,-1)\n";
         throw std::runtime_error( oss.str() );
      }

      sm(1,1) = cplx(-5,2);

      if( sm(0,0) != cplx(-4,-1) || sm(0,1) != cplx( 7,3) ||
          sm(1,0) != cplx( 2, 0) || sm(1,1) != cplx(-5,2) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix access failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( (-4,-1) ( 7,3) )\n"
                                     "( ( 2, 0) (-5,2) )\n";
         throw std::runtime_error( oss.str() );
      }

      if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-4,-1) || herm(0,2) != cplx( 7,3) ||
          herm(1,0) != cplx(-4, 1) || herm(1,1) != cplx( 2, 0) || herm(1,2) != cplx(-5,2) ||
          herm(2,0) != cplx( 7,-3) || herm(2,1) != cplx(-5,-2) || herm(2,2) != cplx( 3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix access failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( ( 1, 0) (-4,-1) ( 7,3) )\n"
                                     "( (-4, 1) ( 2, 0) (-5,2) )\n"
                                     "( ( 7,-3) (-5,-2) ( 3,0) )\n";
         throw std::runtime_error( oss.str() );
      }

      reset( sm );

      if( sm(0,0) != cplx(0,0) || sm(0,1) != cplx(0,0) ||
          sm(1,0) != cplx(0,0) || sm(1,1) != cplx(0,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix reset failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( (0,0) (0,0) )\n"
                                     "( (0,0) (0,0) )\n";
         throw std::runtime_error( oss.str() );
      }

      if( herm(0,0) != cplx(1,0) || herm(0,1) != cplx(0,0) || herm(0,2) != cplx(0,0) ||
          herm(1,0) != cplx(0,0) || herm(1,1) != cplx(0,0) || herm(1,2) != cplx(0,0) ||
          herm(2,0) != cplx(0,0) || herm(2,1) != cplx(0,0) || herm(2,2) != cplx(3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix reset failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (1,0) (0,0) (0,0) )\n"
                                     "( (0,0) (0,0) (0,0) )\n"
                                     "( (0,0) (0,0) (3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c row() function with the HermitianMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c row() function with the HermitianMatrix specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseComplexTest::testRow()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major row() function";

      using RT = blaze::Row<HT>;

      HT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      RT row1 = row( herm, 1UL );

      if( row1[1] != cplx(2,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator access failed\n"
             << " Details:\n"
             << "   Result: " << row1[1] << "\n"
             << "   Expected result: (2,0)\n";
         throw std::runtime_error( oss.str() );
      }

      RT::Iterator it( row1.begin() );

      if( it == row1.end() || it->value() != cplx(-4,1) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << it->value() << "\n"
             << "   Expected result: (-4,1)\n";
         throw std::runtime_error( oss.str() );
      }

      row1[2] = cplx(-5,2);

      if( row1[0] != cplx(-4,1) || row1[1] != cplx(2,0) || row1[2] != cplx(-5,2) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row access failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( (-4,1) (2,0) (-5,2) )\n";
         throw std::runtime_error( oss.str() );
      }

      if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-4,-1) || herm(0,2) != cplx( 7,3) ||
          herm(1,0) != cplx(-4, 1) || herm(1,1) != cplx( 2, 0) || herm(1,2) != cplx(-5,2) ||
          herm(2,0) != cplx( 7,-3) || herm(2,1) != cplx(-5,-2) || herm(2,2) != cplx( 3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row access failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( ( 1, 0) (-4,-1) ( 7,3) )\n"
                                     "( (-4, 1) ( 2, 0) (-5,2) )\n"
                                     "( ( 7,-3) (-5, 0) ( 3,0) )\n";
         throw std::runtime_error( oss.str() );
      }

      reset( row1 );

      if( row1[0] != cplx(0,0) || row1[1] != cplx(0,0) || row1[2] != cplx(0,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row reset failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( (0,0) (0,0) (0,0) )\n";
         throw std::runtime_error( oss.str() );
      }

      if( herm(0,0) != cplx(1, 0) || herm(0,1) != cplx(0,0) || herm(0,2) != cplx(7,3) ||
          herm(1,0) != cplx(0, 0) || herm(1,1) != cplx(0,0) || herm(1,2) != cplx(0,0) ||
          herm(2,0) != cplx(7,-3) || herm(2,1) != cplx(0,0) || herm(2,2) != cplx(3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row reset failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (1, 0) (0,0) (7,3) )\n"
                                     "( (0, 0) (0,0) (0,0) )\n"
                                     "( (7,-3) (0,0) (3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major row() function";

      using RT = blaze::Row<OHT>;

      OHT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      RT row1 = row( herm, 1UL );

      if( row1[1] != cplx(2,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator access failed\n"
             << " Details:\n"
             << "   Result: " << row1[1] << "\n"
             << "   Expected result: (2,0)\n";
         throw std::runtime_error( oss.str() );
      }

      RT::Iterator it( row1.begin() );

      if( it == row1.end() || it->value() != cplx(-4,1) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << it->value() << "\n"
             << "   Expected result: (-4,1)\n";
         throw std::runtime_error( oss.str() );
      }

      row1[2] = cplx(-5,2);

      if( row1[0] != cplx(-4,1) || row1[1] != cplx(2,0) || row1[2] != cplx(-5,2) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row access failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( (-4,1) (2,0) (-5,2) )\n";
         throw std::runtime_error( oss.str() );
      }

      if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-4,-1) || herm(0,2) != cplx( 7,3) ||
          herm(1,0) != cplx(-4, 1) || herm(1,1) != cplx( 2, 0) || herm(1,2) != cplx(-5,2) ||
          herm(2,0) != cplx( 7,-3) || herm(2,1) != cplx(-5,-2) || herm(2,2) != cplx( 3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row access failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( ( 1, 0) (-4,-1) ( 7,3) )\n"
                                     "( (-4, 1) ( 2, 0) (-5,2) )\n"
                                     "( ( 7,-3) (-5, 0) ( 3,0) )\n";
         throw std::runtime_error( oss.str() );
      }

      reset( row1 );

      if( row1[0] != cplx(0,0) || row1[1] != cplx(0,0) || row1[2] != cplx(0,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row reset failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( (0,0) (0,0) (0,0) )\n";
         throw std::runtime_error( oss.str() );
      }

      if( herm(0,0) != cplx(1, 0) || herm(0,1) != cplx(0,0) || herm(0,2) != cplx(7,3) ||
          herm(1,0) != cplx(0, 0) || herm(1,1) != cplx(0,0) || herm(1,2) != cplx(0,0) ||
          herm(2,0) != cplx(7,-3) || herm(2,1) != cplx(0,0) || herm(2,2) != cplx(3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row reset failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (1, 0) (0,0) (7,3) )\n"
                                     "( (0, 0) (0,0) (0,0) )\n"
                                     "( (7,-3) (0,0) (3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c column() function with the HermitianMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c column() function with the HermitianMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseComplexTest::testColumn()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major column() function";

      using CT = blaze::Column<HT>;

      HT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      CT col1 = column( herm, 1UL );

      if( col1[1] != cplx(2,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator access failed\n"
             << " Details:\n"
             << "   Result: " << col1[1] << "\n"
             << "   Expected result: (2,0)\n";
         throw std::runtime_error( oss.str() );
      }

      CT::Iterator it( col1.begin() );

      if( it == col1.end() || it->value() != cplx(-4,-1) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << it->value() << "\n"
             << "   Expected result: (-4,-1)\n";
         throw std::runtime_error( oss.str() );
      }

      col1[2] = cplx(-5,-2);

      if( col1[0] != cplx(-4,-1) || col1[1] != cplx(2,0) || col1[2] != cplx(-5,-2) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column access failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( (-4,-1) ( 2,0) (-5,-2) )\n";
         throw std::runtime_error( oss.str() );
      }

      if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-4,-1) || herm(0,2) != cplx( 7,3) ||
          herm(1,0) != cplx(-4, 1) || herm(1,1) != cplx( 2, 0) || herm(1,2) != cplx(-5,2) ||
          herm(2,0) != cplx( 7,-3) || herm(2,1) != cplx(-5,-2) || herm(2,2) != cplx( 3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column access failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( ( 1, 0) (-4,-1) ( 7,3) )\n"
                                     "( (-4, 1) ( 2, 0) (-5,2) )\n"
                                     "( ( 7,-3) (-5,-2) ( 3,0) )\n";
         throw std::runtime_error( oss.str() );
      }

      reset( col1 );

      if( col1[0] != cplx(0,0) || col1[1] != cplx(0,0) || col1[2] != cplx(0,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column reset failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( (0,0) (0,0) (0,0) )\n";
         throw std::runtime_error( oss.str() );
      }

      if( herm(0,0) != cplx(1, 0) || herm(0,1) != cplx(0,0) || herm(0,2) != cplx(7,3) ||
          herm(1,0) != cplx(0, 0) || herm(1,1) != cplx(0,0) || herm(1,2) != cplx(0,0) ||
          herm(2,0) != cplx(7,-3) || herm(2,1) != cplx(0,0) || herm(2,2) != cplx(3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column reset failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (1, 0) (0,0) (7,3) )\n"
                                     "( (0, 0) (0,0) (0,0) )\n"
                                     "( (7,-3) (0,0) (3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major column() function";

      using CT = blaze::Column<OHT>;

      OHT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      CT col1 = column( herm, 1UL );

      if( col1[1] != cplx(2,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator access failed\n"
             << " Details:\n"
             << "   Result: " << col1[1] << "\n"
             << "   Expected result: (2,0)\n";
         throw std::runtime_error( oss.str() );
      }

      CT::Iterator it( col1.begin() );

      if( it == col1.end() || it->value() != cplx(-4,-1) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << it->value() << "\n"
             << "   Expected result: (-4,-1)\n";
         throw std::runtime_error( oss.str() );
      }

      col1[2] = cplx(-5,-2);

      if( col1[0] != cplx(-4,-1) || col1[1] != cplx(2,0) || col1[2] != cplx(-5,-2) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column access failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( (-4,-1) ( 2,0) (-5,-2) )\n";
         throw std::runtime_error( oss.str() );
      }

      if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-4,-1) || herm(0,2) != cplx( 7,3) ||
          herm(1,0) != cplx(-4, 1) || herm(1,1) != cplx( 2, 0) || herm(1,2) != cplx(-5,2) ||
          herm(2,0) != cplx( 7,-3) || herm(2,1) != cplx(-5,-2) || herm(2,2) != cplx( 3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column access failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( ( 1, 0) (-4,-1) ( 7,3) )\n"
                                     "( (-4, 1) ( 2, 0) (-5,2) )\n"
                                     "( ( 7,-3) (-5,-2) ( 3,0) )\n";
         throw std::runtime_error( oss.str() );
      }

      reset( col1 );

      if( col1[0] != cplx(0,0) || col1[1] != cplx(0,0) || col1[2] != cplx(0,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column reset failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( (0,0) (0,0) (0,0) )\n";
         throw std::runtime_error( oss.str() );
      }

      if( herm(0,0) != cplx(1, 0) || herm(0,1) != cplx(0,0) || herm(0,2) != cplx(7,3) ||
          herm(1,0) != cplx(0, 0) || herm(1,1) != cplx(0,0) || herm(1,2) != cplx(0,0) ||
          herm(2,0) != cplx(7,-3) || herm(2,1) != cplx(0,0) || herm(2,2) != cplx(3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column reset failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (1, 0) (0,0) (7,3) )\n"
                                     "( (0, 0) (0,0) (0,0) )\n"
                                     "( (7,-3) (0,0) (3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************

} // namespace hermitianmatrix

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
   std::cout << "   Running HermitianMatrix sparse complex test (part 2)..." << std::endl;

   try
   {
      RUN_HERMITIANMATRIX_SPARSECOMPLEX_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during HermitianMatrix sparse complex test (part 2):\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
