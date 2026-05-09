//=================================================================================================
/*!
//  \file src/mathtest/adaptors/lowermatrix/DenseTest2.cpp
//  \brief Source file for the LowerMatrix dense test (part 2)
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
#include <blaze/math/CompressedMatrix.h>
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/Row.h>
#include <blaze/math/Submatrix.h>
#include <blaze/util/Complex.h>
#include <blazetest/mathtest/adaptors/lowermatrix/DenseTest.h>

#ifdef BLAZE_USE_HPX_THREADS
#  include <hpx/hpx_main.hpp>
#endif


namespace blazetest {

namespace mathtest {

namespace adaptors {

namespace lowermatrix {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the LowerMatrix dense test.
//
// \exception std::runtime_error Operation error detected.
*/
DenseTest::DenseTest()
{
   testMultAssign();
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
/*!\brief Test of the LowerMatrix multiplication assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the multiplication assignment operators of the LowerMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testMultAssign()
{
   //=====================================================================================
   // Row-major dense matrix multiplication assignment
   //=====================================================================================

   // Row-major/row-major dense matrix multiplication assignment (lower)
   {
      test_ = "Row-major/row-major LowerMatrix dense matrix multiplication assignment (lower)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(0,0) = 2;
      mat(1,1) = 2;
      mat(2,2) = 2;

      LT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

      lower *= mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 2UL );

      if( lower(0,0) !=  2 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -8 || lower(1,1) != 4 || lower(1,2) != 0 ||
          lower(2,0) != 14 || lower(2,1) != 0 || lower(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  2 0 0 )\n( -8 4 0 )\n( 14 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix multiplication assignment (lower)
   {
      test_ = "Row-major/column-major LowerMatrix dense matrix multiplication assignment (lower)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(0,0) = 2;
      mat(1,1) = 2;
      mat(2,2) = 2;

      LT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

      lower *= mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 2UL );

      if( lower(0,0) !=  2 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -8 || lower(1,1) != 4 || lower(1,2) != 0 ||
          lower(2,0) != 14 || lower(2,1) != 0 || lower(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  2 0 0 )\n( -8 4 0 )\n( 14 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major dense matrix multiplication assignment (non-lower)
   {
      test_ = "Row-major/row-major LowerMatrix dense matrix multiplication assignment (non-lower)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      LT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

      try {
         lower *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of non-lower row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major dense matrix multiplication assignment (non-lower)
   {
      test_ = "Row-major/column-major LowerMatrix dense matrix multiplication assignment (non-lower)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      LT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

      try {
         lower *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of non-lower column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major dense matrix multiplication assignment (LowerMatrix)
   {
      test_ = "Row-major/row-major LowerMatrix dense matrix multiplication assignment (LowerMatrix)";

      LT lower1( 3UL );
      lower1(0,0) = 2;
      lower1(1,1) = 2;
      lower1(2,2) = 2;

      LT lower2( 3UL );
      lower2(0,0) =  1;
      lower2(1,0) = -4;
      lower2(1,1) =  2;
      lower2(2,0) =  7;
      lower2(2,2) =  3;

      lower2 *= lower1;

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkCapacity( lower2, 9UL );
      checkNonZeros( lower2, 5UL );
      checkNonZeros( lower2, 0UL, 1UL );
      checkNonZeros( lower2, 1UL, 2UL );
      checkNonZeros( lower2, 2UL, 2UL );

      if( lower2(0,0) !=  2 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -8 || lower2(1,1) != 4 || lower2(1,2) != 0 ||
          lower2(2,0) != 14 || lower2(2,1) != 0 || lower2(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  2 0 0 )\n( -8 4 0 )\n( 14 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix multiplication assignment (LowerMatrix)
   {
      test_ = "Row-major/column-major LowerMatrix dense matrix multiplication assignment (LowerMatrix)";

      OLT lower1( 3UL );
      lower1(0,0) = 2;
      lower1(1,1) = 2;
      lower1(2,2) = 2;

      LT lower2( 3UL );
      lower2(0,0) =  1;
      lower2(1,0) = -4;
      lower2(1,1) =  2;
      lower2(2,0) =  7;
      lower2(2,2) =  3;

      lower2 *= lower1;

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkCapacity( lower2, 9UL );
      checkNonZeros( lower2, 5UL );
      checkNonZeros( lower2, 0UL, 1UL );
      checkNonZeros( lower2, 1UL, 2UL );
      checkNonZeros( lower2, 2UL, 2UL );

      if( lower2(0,0) !=  2 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -8 || lower2(1,1) != 4 || lower2(1,2) != 0 ||
          lower2(2,0) != 14 || lower2(2,1) != 0 || lower2(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  2 0 0 )\n( -8 4 0 )\n( 14 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix multiplication assignment
   //=====================================================================================

   // Row-major/row-major sparse matrix multiplication assignment (lower)
   {
      test_ = "Row-major/row-major LowerMatrix sparse matrix multiplication assignment (lower)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 4UL );
      mat(0,0) = 2;
      mat(1,1) = 2;
      mat(2,2) = 2;
      mat.insert( 1UL, 2UL, 0 );

      LT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

      lower *= mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 2UL );

      if( lower(0,0) !=  2 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -8 || lower(1,1) != 4 || lower(1,2) != 0 ||
          lower(2,0) != 14 || lower(2,1) != 0 || lower(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  2 0 0 )\n( -8 4 0 )\n( 14 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major sparse matrix multiplication assignment (lower)
   {
      test_ = "Row-major/column-major LowerMatrix sparse matrix multiplication assignment (lower)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 4UL );
      mat(0,0) = 2;
      mat(1,1) = 2;
      mat(2,2) = 2;
      mat.insert( 1UL, 2UL, 0 );

      LT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

      lower *= mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 2UL );

      if( lower(0,0) !=  2 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -8 || lower(1,1) != 4 || lower(1,2) != 0 ||
          lower(2,0) != 14 || lower(2,1) != 0 || lower(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  2 0 0 )\n( -8 4 0 )\n( 14 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major sparse matrix multiplication assignment (non-lower)
   {
      test_ = "Row-major/row-major LowerMatrix sparse matrix multiplication assignment (non-lower)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 4UL );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      LT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

      try {
         lower *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of non-lower row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major sparse matrix multiplication assignment (non-lower)
   {
      test_ = "Row-major/column-major LowerMatrix sparse matrix multiplication assignment (non-lower)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 4UL );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      LT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

      try {
         lower *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of non-lower column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major sparse matrix multiplication assignment (LowerMatrix)
   {
      test_ = "Row-major/row-major LowerMatrix sparse matrix multiplication assignment (LowerMatrix)";

      blaze::LowerMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > lower1( 3UL, 3UL );
      lower1(0,0) = 2;
      lower1(1,1) = 2;
      lower1(2,2) = 2;

      LT lower2( 3UL );
      lower2(0,0) =  1;
      lower2(1,0) = -4;
      lower2(1,1) =  2;
      lower2(2,0) =  7;
      lower2(2,2) =  3;

      lower2 *= lower1;

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkCapacity( lower2, 9UL );
      checkNonZeros( lower2, 5UL );
      checkNonZeros( lower2, 0UL, 1UL );
      checkNonZeros( lower2, 1UL, 2UL );
      checkNonZeros( lower2, 2UL, 2UL );

      if( lower2(0,0) !=  2 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -8 || lower2(1,1) != 4 || lower2(1,2) != 0 ||
          lower2(2,0) != 14 || lower2(2,1) != 0 || lower2(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  2 0 0 )\n( -8 4 0 )\n( 14 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major sparse matrix multiplication assignment (LowerMatrix)
   {
      test_ = "Row-major/column-major LowerMatrix sparse matrix multiplication assignment (LowerMatrix)";

      blaze::LowerMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > lower1( 3UL, 3UL );
      lower1(0,0) = 2;
      lower1(1,1) = 2;
      lower1(2,2) = 2;

      LT lower2( 3UL );
      lower2(0,0) =  1;
      lower2(1,0) = -4;
      lower2(1,1) =  2;
      lower2(2,0) =  7;
      lower2(2,2) =  3;

      lower2 *= lower1;

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkCapacity( lower2, 9UL );
      checkNonZeros( lower2, 5UL );
      checkNonZeros( lower2, 0UL, 1UL );
      checkNonZeros( lower2, 1UL, 2UL );
      checkNonZeros( lower2, 2UL, 2UL );

      if( lower2(0,0) !=  2 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -8 || lower2(1,1) != 4 || lower2(1,2) != 0 ||
          lower2(2,0) != 14 || lower2(2,1) != 0 || lower2(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  2 0 0 )\n( -8 4 0 )\n( 14 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix multiplication assignment
   //=====================================================================================

   // Column-major/row-major dense matrix multiplication assignment (lower)
   {
      test_ = "Column-major/row-major LowerMatrix dense matrix multiplication assignment (lower)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(0,0) = 2;
      mat(1,1) = 2;
      mat(2,2) = 2;

      OLT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

      lower *= mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 3UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) !=  2 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -8 || lower(1,1) != 4 || lower(1,2) != 0 ||
          lower(2,0) != 14 || lower(2,1) != 0 || lower(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  2 0 0 )\n( -8 4 0 )\n( 14 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix multiplication assignment (lower)
   {
      test_ = "Column-major/column-major LowerMatrix dense matrix multiplication assignment (lower)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(0,0) = 2;
      mat(1,1) = 2;
      mat(2,2) = 2;

      OLT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

      lower *= mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 3UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) !=  2 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -8 || lower(1,1) != 4 || lower(1,2) != 0 ||
          lower(2,0) != 14 || lower(2,1) != 0 || lower(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  2 0 0 )\n( -8 4 0 )\n( 14 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major dense matrix multiplication assignment (non-lower)
   {
      test_ = "Column-major/row-major LowerMatrix dense matrix multiplication assignment (non-lower)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      OLT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

      try {
         lower *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of non-lower row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major dense matrix multiplication assignment (non-lower)
   {
      test_ = "Column-major/column-major LowerMatrix dense matrix multiplication assignment (non-lower)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      OLT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

      try {
         lower *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of non-lower column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major dense matrix multiplication assignment (LowerMatrix)
   {
      test_ = "Column-major/row-major LowerMatrix dense matrix multiplication assignment (LowerMatrix)";

      LT lower1( 3UL );
      lower1(0,0) = 2;
      lower1(1,1) = 2;
      lower1(2,2) = 2;

      OLT lower2( 3UL );
      lower2(0,0) =  1;
      lower2(1,0) = -4;
      lower2(1,1) =  2;
      lower2(2,0) =  7;
      lower2(2,2) =  3;

      lower2 *= lower1;

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkCapacity( lower2, 9UL );
      checkNonZeros( lower2, 5UL );
      checkNonZeros( lower2, 0UL, 3UL );
      checkNonZeros( lower2, 1UL, 1UL );
      checkNonZeros( lower2, 2UL, 1UL );

      if( lower2(0,0) !=  2 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -8 || lower2(1,1) != 4 || lower2(1,2) != 0 ||
          lower2(2,0) != 14 || lower2(2,1) != 0 || lower2(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  2 0 0 )\n( -8 4 0 )\n( 14 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix multiplication assignment (LowerMatrix)
   {
      test_ = "Column-major/column-major LowerMatrix dense matrix multiplication assignment (LowerMatrix)";

      OLT lower1( 3UL );
      lower1(0,0) = 2;
      lower1(1,1) = 2;
      lower1(2,2) = 2;

      OLT lower2( 3UL );
      lower2(0,0) =  1;
      lower2(1,0) = -4;
      lower2(1,1) =  2;
      lower2(2,0) =  7;
      lower2(2,2) =  3;

      lower2 *= lower1;

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkCapacity( lower2, 9UL );
      checkNonZeros( lower2, 5UL );
      checkNonZeros( lower2, 0UL, 3UL );
      checkNonZeros( lower2, 1UL, 1UL );
      checkNonZeros( lower2, 2UL, 1UL );

      if( lower2(0,0) !=  2 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -8 || lower2(1,1) != 4 || lower2(1,2) != 0 ||
          lower2(2,0) != 14 || lower2(2,1) != 0 || lower2(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  2 0 0 )\n( -8 4 0 )\n( 14 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse matrix multiplication assignment
   //=====================================================================================

   // Column-major/row-major sparse matrix multiplication assignment (lower)
   {
      test_ = "Column-major/row-major LowerMatrix sparse matrix multiplication assignment (lower)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 4UL );
      mat(0,0) = 2;
      mat(1,1) = 2;
      mat(2,2) = 2;
      mat.insert( 1UL, 2UL, 0 );

      OLT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

      lower *= mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 3UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) !=  2 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -8 || lower(1,1) != 4 || lower(1,2) != 0 ||
          lower(2,0) != 14 || lower(2,1) != 0 || lower(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  2 0 0 )\n( -8 4 0 )\n( 14 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major sparse matrix multiplication assignment (lower)
   {
      test_ = "Column-major/column-major LowerMatrix sparse matrix multiplication assignment (lower)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 4UL );
      mat(0,0) = 2;
      mat(1,1) = 2;
      mat(2,2) = 2;
      mat.insert( 1UL, 2UL, 0 );

      OLT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

      lower *= mat;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 3UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) !=  2 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -8 || lower(1,1) != 4 || lower(1,2) != 0 ||
          lower(2,0) != 14 || lower(2,1) != 0 || lower(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  2 0 0 )\n( -8 4 0 )\n( 14 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major sparse matrix multiplication assignment (non-lower)
   {
      test_ = "Column-major/row-major LowerMatrix sparse matrix multiplication assignment (non-lower)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 4UL );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      OLT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

      try {
         lower *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of non-lower row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major sparse matrix multiplication assignment (non-lower)
   {
      test_ = "Column-major/column-major LowerMatrix sparse matrix multiplication assignment (non-lower)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 4UL );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      OLT lower( 3UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

      try {
         lower *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of non-lower column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major sparse matrix multiplication assignment (LowerMatrix)
   {
      test_ = "Column-major/row-major LowerMatrix sparse matrix multiplication assignment (LowerMatrix)";

      blaze::LowerMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > lower1( 3UL, 3UL );
      lower1(0,0) = 2;
      lower1(1,1) = 2;
      lower1(2,2) = 2;

      OLT lower2( 3UL );
      lower2(0,0) =  1;
      lower2(1,0) = -4;
      lower2(1,1) =  2;
      lower2(2,0) =  7;
      lower2(2,2) =  3;

      lower2 *= lower1;

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkCapacity( lower2, 9UL );
      checkNonZeros( lower2, 5UL );
      checkNonZeros( lower2, 0UL, 3UL );
      checkNonZeros( lower2, 1UL, 1UL );
      checkNonZeros( lower2, 2UL, 1UL );

      if( lower2(0,0) !=  2 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -8 || lower2(1,1) != 4 || lower2(1,2) != 0 ||
          lower2(2,0) != 14 || lower2(2,1) != 0 || lower2(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  2 0 0 )\n( -8 4 0 )\n( 14 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major sparse matrix multiplication assignment (LowerMatrix)
   {
      test_ = "Column-major/column-major LowerMatrix sparse matrix multiplication assignment (LowerMatrix)";

      blaze::LowerMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > lower1( 3UL, 3UL );
      lower1(0,0) = 2;
      lower1(1,1) = 2;
      lower1(2,2) = 2;

      OLT lower2( 3UL );
      lower2(0,0) =  1;
      lower2(1,0) = -4;
      lower2(1,1) =  2;
      lower2(2,0) =  7;
      lower2(2,2) =  3;

      lower2 *= lower1;

      checkRows    ( lower2, 3UL );
      checkColumns ( lower2, 3UL );
      checkCapacity( lower2, 9UL );
      checkNonZeros( lower2, 5UL );
      checkNonZeros( lower2, 0UL, 3UL );
      checkNonZeros( lower2, 1UL, 1UL );
      checkNonZeros( lower2, 2UL, 1UL );

      if( lower2(0,0) !=  2 || lower2(0,1) != 0 || lower2(0,2) != 0 ||
          lower2(1,0) != -8 || lower2(1,1) != 4 || lower2(1,2) != 0 ||
          lower2(2,0) != 14 || lower2(2,1) != 0 || lower2(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n(  2 0 0 )\n( -8 4 0 )\n( 14 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of all LowerMatrix (self-)scaling operations.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all available ways to scale an instance of the LowerMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testScaling()
{
   //=====================================================================================
   // Row-major self-scaling (M*=s)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (M*=s)";

      LT lower( 3UL );
      lower(1,0) =  1;
      lower(2,0) = -2;
      lower(2,1) =  3;
      lower(2,2) = -4;

      lower *= 2;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 4UL );
      checkNonZeros( lower, 0UL, 0UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 3UL );

      if( lower(0,0) !=  0 || lower(0,1) != 0 || lower(0,2) !=  0 ||
          lower(1,0) !=  2 || lower(1,1) != 0 || lower(1,2) !=  0 ||
          lower(2,0) != -4 || lower(2,1) != 6 || lower(2,2) != -8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  0  0  0 )\n(  2  0  0 )\n( -4  6 -8 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major self-scaling (M=M*s)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (M=M*s)";

      LT lower( 3UL );
      lower(1,0) =  1;
      lower(2,0) = -2;
      lower(2,1) =  3;
      lower(2,2) = -4;

      lower = lower * 2;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 4UL );
      checkNonZeros( lower, 0UL, 0UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 3UL );

      if( lower(0,0) !=  0 || lower(0,1) != 0 || lower(0,2) !=  0 ||
          lower(1,0) !=  2 || lower(1,1) != 0 || lower(1,2) !=  0 ||
          lower(2,0) != -4 || lower(2,1) != 6 || lower(2,2) != -8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  0  0  0 )\n(  2  0  0 )\n( -4  6 -8 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major self-scaling (M=s*M)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (M=s*M)";

      LT lower( 3UL );
      lower(1,0) =  1;
      lower(2,0) = -2;
      lower(2,1) =  3;
      lower(2,2) = -4;

      lower = 2 * lower;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 4UL );
      checkNonZeros( lower, 0UL, 0UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 3UL );

      if( lower(0,0) !=  0 || lower(0,1) != 0 || lower(0,2) !=  0 ||
          lower(1,0) !=  2 || lower(1,1) != 0 || lower(1,2) !=  0 ||
          lower(2,0) != -4 || lower(2,1) != 6 || lower(2,2) != -8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  0  0  0 )\n(  2  0  0 )\n( -4  6 -8 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major self-scaling (M/=s)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (M/=s)";

      LT lower( 3UL );
      lower(1,0) =  2;
      lower(2,0) = -4;
      lower(2,1) =  6;
      lower(2,2) = -8;

      lower /= 2;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 4UL );
      checkNonZeros( lower, 0UL, 0UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 3UL );

      if( lower(0,0) !=  0 || lower(0,1) != 0 || lower(0,2) !=  0 ||
          lower(1,0) !=  1 || lower(1,1) != 0 || lower(1,2) !=  0 ||
          lower(2,0) != -2 || lower(2,1) != 3 || lower(2,2) != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  0  0  0 )\n(  1  0  0 )\n( -2  3 -4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major self-scaling (M=M/s)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (M=M/s)";

      LT lower( 3UL );
      lower(1,0) =  2;
      lower(2,0) = -4;
      lower(2,1) =  6;
      lower(2,2) = -8;

      lower = lower / 2;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 4UL );
      checkNonZeros( lower, 0UL, 0UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 3UL );

      if( lower(0,0) !=  0 || lower(0,1) != 0 || lower(0,2) !=  0 ||
          lower(1,0) !=  1 || lower(1,1) != 0 || lower(1,2) !=  0 ||
          lower(2,0) != -2 || lower(2,1) != 3 || lower(2,2) != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  0  0  0 )\n(  1  0  0 )\n( -2  3 -4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major LowerMatrix::scale()
   //=====================================================================================

   {
      test_ = "Row-major LowerMatrix::scale()";

      // Initialization check
      LT lower( 3UL );
      lower(1,0) =  1;
      lower(2,0) = -2;
      lower(2,1) =  3;
      lower(2,2) = -4;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 4UL );
      checkNonZeros( lower, 0UL, 0UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 3UL );

      if( lower(0,0) !=  0 || lower(0,1) != 0 || lower(0,2) !=  0 ||
          lower(1,0) !=  1 || lower(1,1) != 0 || lower(1,2) !=  0 ||
          lower(2,0) != -2 || lower(2,1) != 3 || lower(2,2) != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  0  0  0 )\n(  1  0  0 )\n( -2  3 -4 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Integral scaling of the matrix
      lower.scale( 2 );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 4UL );
      checkNonZeros( lower, 0UL, 0UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 3UL );

      if( lower(0,0) !=  0 || lower(0,1) != 0 || lower(0,2) !=  0 ||
          lower(1,0) !=  2 || lower(1,1) != 0 || lower(1,2) !=  0 ||
          lower(2,0) != -4 || lower(2,1) != 6 || lower(2,2) != -8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  0  0  0 )\n(  2  0  0 )\n( -4  6 -8 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Floating point scaling of the matrix
      lower.scale( 0.5 );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 4UL );
      checkNonZeros( lower, 0UL, 0UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 3UL );

      if( lower(0,0) !=  0 || lower(0,1) != 0 || lower(0,2) !=  0 ||
          lower(1,0) !=  1 || lower(1,1) != 0 || lower(1,2) !=  0 ||
          lower(2,0) != -2 || lower(2,1) != 3 || lower(2,2) != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  0  0  0 )\n(  1  0  0 )\n( -2  3 -4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major LowerMatrix::scale() (complex)";

      using blaze::complex;

      blaze::LowerMatrix< blaze::DynamicMatrix<complex<float>,blaze::rowMajor> > lower( 2UL );
      lower(0,0) = complex<float>( 1.0F, 0.0F );
      lower(1,0) = complex<float>( 2.0F, 0.0F );
      lower(1,1) = complex<float>( 4.0F, 0.0F );

      lower.scale( complex<float>( 3.0F, 0.0F ) );

      checkRows    ( lower, 2UL );
      checkColumns ( lower, 2UL );
      checkCapacity( lower, 4UL );
      checkNonZeros( lower, 3UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );

      if( lower(0,0) != complex<float>( 3.0F, 0.0F ) || lower(0,1) != complex<float>(  0.0F, 0.0F ) ||
          lower(1,0) != complex<float>( 6.0F, 0.0F ) || lower(1,1) != complex<float>( 12.0F, 0.0F ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( ( 3,0) ( 0,0)\n( 6,0) (12,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (M*=s)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (M*=s)";

      OLT lower( 3UL );
      lower(1,0) =  1;
      lower(2,0) = -2;
      lower(2,1) =  3;
      lower(2,2) = -4;

      lower *= 2;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 4UL );
      checkNonZeros( lower, 0UL, 2UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) !=  0 || lower(0,1) != 0 || lower(0,2) !=  0 ||
          lower(1,0) !=  2 || lower(1,1) != 0 || lower(1,2) !=  0 ||
          lower(2,0) != -4 || lower(2,1) != 6 || lower(2,2) != -8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  0  0  0 )\n(  2  0  0 )\n( -4  6 -8 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (M=M*s)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (M=M*s)";

      OLT lower( 3UL );
      lower(1,0) =  1;
      lower(2,0) = -2;
      lower(2,1) =  3;
      lower(2,2) = -4;

      lower = lower * 2;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 4UL );
      checkNonZeros( lower, 0UL, 2UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) !=  0 || lower(0,1) != 0 || lower(0,2) !=  0 ||
          lower(1,0) !=  2 || lower(1,1) != 0 || lower(1,2) !=  0 ||
          lower(2,0) != -4 || lower(2,1) != 6 || lower(2,2) != -8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  0  0  0 )\n(  2  0  0 )\n( -4  6 -8 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (M=s*M)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (M=s*M)";

      OLT lower( 3UL );
      lower(1,0) =  1;
      lower(2,0) = -2;
      lower(2,1) =  3;
      lower(2,2) = -4;

      lower = 2 * lower;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 4UL );
      checkNonZeros( lower, 0UL, 2UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) !=  0 || lower(0,1) != 0 || lower(0,2) !=  0 ||
          lower(1,0) !=  2 || lower(1,1) != 0 || lower(1,2) !=  0 ||
          lower(2,0) != -4 || lower(2,1) != 6 || lower(2,2) != -8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  0  0  0 )\n(  2  0  0 )\n( -4  6 -8 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (M/=s)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (M/=s)";

      OLT lower( 3UL );
      lower(1,0) =  2;
      lower(2,0) = -4;
      lower(2,1) =  6;
      lower(2,2) = -8;

      lower /= 2;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 4UL );
      checkNonZeros( lower, 0UL, 2UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) !=  0 || lower(0,1) != 0 || lower(0,2) !=  0 ||
          lower(1,0) !=  1 || lower(1,1) != 0 || lower(1,2) !=  0 ||
          lower(2,0) != -2 || lower(2,1) != 3 || lower(2,2) != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  0  0  0 )\n(  1  0  0 )\n( -2  3 -4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (M=M/s)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (M=M/s)";

      OLT lower( 3UL );
      lower(1,0) =  2;
      lower(2,0) = -4;
      lower(2,1) =  6;
      lower(2,2) = -8;

      lower = lower / 2;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 4UL );
      checkNonZeros( lower, 0UL, 2UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) !=  0 || lower(0,1) != 0 || lower(0,2) !=  0 ||
          lower(1,0) !=  1 || lower(1,1) != 0 || lower(1,2) !=  0 ||
          lower(2,0) != -2 || lower(2,1) != 3 || lower(2,2) != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  0  0  0 )\n(  1  0  0 )\n( -2  3 -4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major LowerMatrix::scale()
   //=====================================================================================

   {
      test_ = "Column-major LowerMatrix::scale()";

      // Initialization check
      OLT lower( 3UL );
      lower(1,0) =  1;
      lower(2,0) = -2;
      lower(2,1) =  3;
      lower(2,2) = -4;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 4UL );
      checkNonZeros( lower, 0UL, 2UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) !=  0 || lower(0,1) != 0 || lower(0,2) !=  0 ||
          lower(1,0) !=  1 || lower(1,1) != 0 || lower(1,2) !=  0 ||
          lower(2,0) != -2 || lower(2,1) != 3 || lower(2,2) != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  0  0  0 )\n(  1  0  0 )\n( -2  3 -4 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Integral scaling of the matrix
      lower.scale( 2 );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 4UL );
      checkNonZeros( lower, 0UL, 2UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) !=  0 || lower(0,1) != 0 || lower(0,2) !=  0 ||
          lower(1,0) !=  2 || lower(1,1) != 0 || lower(1,2) !=  0 ||
          lower(2,0) != -4 || lower(2,1) != 6 || lower(2,2) != -8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  0  0  0 )\n(  2  0  0 )\n( -4  6 -8 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Floating point scaling of the matrix
      lower.scale( 0.5 );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 4UL );
      checkNonZeros( lower, 0UL, 2UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) !=  0 || lower(0,1) != 0 || lower(0,2) !=  0 ||
          lower(1,0) !=  1 || lower(1,1) != 0 || lower(1,2) !=  0 ||
          lower(2,0) != -2 || lower(2,1) != 3 || lower(2,2) != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  0  0  0 )\n(  1  0  0 )\n( -2  3 -4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major LowerMatrix::scale() (complex)";

      using blaze::complex;

      blaze::LowerMatrix< blaze::DynamicMatrix<complex<float>,blaze::columnMajor> > lower( 2UL );
      lower(0,0) = complex<float>( 1.0F, 0.0F );
      lower(1,0) = complex<float>( 2.0F, 0.0F );
      lower(1,1) = complex<float>( 4.0F, 0.0F );

      lower.scale( complex<float>( 3.0F, 0.0F ) );

      checkRows    ( lower, 2UL );
      checkColumns ( lower, 2UL );
      checkCapacity( lower, 4UL );
      checkNonZeros( lower, 3UL );
      checkNonZeros( lower, 0UL, 2UL );
      checkNonZeros( lower, 1UL, 1UL );

      if( lower(0,0) != complex<float>( 3.0F, 0.0F ) || lower(0,1) != complex<float>(  0.0F, 0.0F ) ||
          lower(1,0) != complex<float>( 6.0F, 0.0F ) || lower(1,1) != complex<float>( 12.0F, 0.0F ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( ( 3,0) ( 0,0)\n( 6,0) (12,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the LowerMatrix function call operator.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of adding and accessing elements via the function call operator
// of the LowerMatrix specialization. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
void DenseTest::testFunctionCall()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major LowerMatrix::operator()";

      // Good cases
      {
         LT lower( 3UL );

         // Writing the diagonal element (1,1)
         lower(1,1) = 1;

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 9UL );
         checkNonZeros( lower, 1UL );
         checkNonZeros( lower, 0UL, 0UL );
         checkNonZeros( lower, 1UL, 1UL );
         checkNonZeros( lower, 2UL, 0UL );

         if( lower(0,0) != 0 || lower(0,1) != 0 || lower(0,2) != 0 ||
             lower(1,0) != 0 || lower(1,1) != 1 || lower(1,2) != 0 ||
             lower(2,0) != 0 || lower(2,1) != 0 || lower(2,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 0 0 0 )\n( 0 1 0 )\n( 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Writing the lower element (2,1)
         lower(2,1) = 2;

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 9UL );
         checkNonZeros( lower, 2UL );
         checkNonZeros( lower, 0UL, 0UL );
         checkNonZeros( lower, 1UL, 1UL );
         checkNonZeros( lower, 2UL, 1UL );

         if( lower(0,0) != 0 || lower(0,1) != 0 || lower(0,2) != 0 ||
             lower(1,0) != 0 || lower(1,1) != 1 || lower(1,2) != 0 ||
             lower(2,0) != 0 || lower(2,1) != 2 || lower(2,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 0 0 0 )\n( 0 1 0 )\n( 0 2 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Writing the lower element (1,0)
         lower(1,0) = lower(2,1);

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 9UL );
         checkNonZeros( lower, 3UL );
         checkNonZeros( lower, 0UL, 0UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 1UL );

         if( lower(0,0) != 0 || lower(0,1) != 0 || lower(0,2) != 0 ||
             lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 ||
             lower(2,0) != 0 || lower(2,1) != 2 || lower(2,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 0 0 0 )\n( 2 1 0 )\n( 0 2 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Adding to the lower element (2,0)
         lower(2,0) += 3;

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 9UL );
         checkNonZeros( lower, 4UL );
         checkNonZeros( lower, 0UL, 0UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 2UL );

         if( lower(0,0) != 0 || lower(0,1) != 0 || lower(0,2) != 0 ||
             lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 ||
             lower(2,0) != 3 || lower(2,1) != 2 || lower(2,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 0 0 0 )\n( 2 1 0 )\n( 3 2 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Subtracting from the lower element (1,0)
         lower(1,0) -= 4;

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 9UL );
         checkNonZeros( lower, 4UL );
         checkNonZeros( lower, 0UL, 0UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 2UL );

         if( lower(0,0) !=  0 || lower(0,1) != 0 || lower(0,2) != 0 ||
             lower(1,0) != -2 || lower(1,1) != 1 || lower(1,2) != 0 ||
             lower(2,0) !=  3 || lower(2,1) != 2 || lower(2,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n(  0 0 0 )\n( -2 1 0 )\n(  3 2 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Multiplying the lower element (2,1)
         lower(2,1) *= -3;

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 9UL );
         checkNonZeros( lower, 4UL );
         checkNonZeros( lower, 0UL, 0UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 2UL );

         if( lower(0,0) !=  0 || lower(0,1) !=  0 || lower(0,2) != 0 ||
             lower(1,0) != -2 || lower(1,1) !=  1 || lower(1,2) != 0 ||
             lower(2,0) !=  3 || lower(2,1) != -6 || lower(2,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n(  0  0  0 )\n( -2  1  0 )\n(  3 -6  0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Dividing the lower element (2,1)
         lower(2,1) /= 2;

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 9UL );
         checkNonZeros( lower, 4UL );
         checkNonZeros( lower, 0UL, 0UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 2UL );

         if( lower(0,0) !=  0 || lower(0,1) !=  0 || lower(0,2) != 0 ||
             lower(1,0) != -2 || lower(1,1) !=  1 || lower(1,2) != 0 ||
             lower(2,0) !=  3 || lower(2,1) != -3 || lower(2,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n(  0  0  0 )\n( -2  1  0 )\n(  3 -3  0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Failure cases
      {
         LT lower( 3UL );

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
      test_ = "Column-major LowerMatrix::operator()";

      // Good cases
      {
         OLT lower( 3UL );

         // Writing the diagonal element (1,1)
         lower(1,1) = 1;

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 9UL );
         checkNonZeros( lower, 1UL );
         checkNonZeros( lower, 0UL, 0UL );
         checkNonZeros( lower, 1UL, 1UL );
         checkNonZeros( lower, 2UL, 0UL );

         if( lower(0,0) != 0 || lower(0,1) != 0 || lower(0,2) != 0 ||
             lower(1,0) != 0 || lower(1,1) != 1 || lower(1,2) != 0 ||
             lower(2,0) != 0 || lower(2,1) != 0 || lower(2,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 0 0 0 )\n( 0 1 0 )\n( 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Writing the lower element (2,1)
         lower(2,1) = 2;

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 9UL );
         checkNonZeros( lower, 2UL );
         checkNonZeros( lower, 0UL, 0UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 0UL );

         if( lower(0,0) != 0 || lower(0,1) != 0 || lower(0,2) != 0 ||
             lower(1,0) != 0 || lower(1,1) != 1 || lower(1,2) != 0 ||
             lower(2,0) != 0 || lower(2,1) != 2 || lower(2,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 0 0 0 )\n( 0 1 0 )\n( 0 2 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Writing the lower element (1,0)
         lower(1,0) = lower(2,1);

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 9UL );
         checkNonZeros( lower, 3UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 0UL );

         if( lower(0,0) != 0 || lower(0,1) != 0 || lower(0,2) != 0 ||
             lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 ||
             lower(2,0) != 0 || lower(2,1) != 2 || lower(2,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 0 0 0 )\n( 2 1 0 )\n( 0 2 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Adding to the lower element (2,0)
         lower(2,0) += 3;

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 9UL );
         checkNonZeros( lower, 4UL );
         checkNonZeros( lower, 0UL, 2UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 0UL );

         if( lower(0,0) != 0 || lower(0,1) != 0 || lower(0,2) != 0 ||
             lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 ||
             lower(2,0) != 3 || lower(2,1) != 2 || lower(2,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 0 0 0 )\n( 2 1 0 )\n( 3 2 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Subtracting from the lower element (1,0)
         lower(1,0) -= 4;

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 9UL );
         checkNonZeros( lower, 4UL );
         checkNonZeros( lower, 0UL, 2UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 0UL );

         if( lower(0,0) !=  0 || lower(0,1) != 0 || lower(0,2) != 0 ||
             lower(1,0) != -2 || lower(1,1) != 1 || lower(1,2) != 0 ||
             lower(2,0) !=  3 || lower(2,1) != 2 || lower(2,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n(  0 0 0 )\n( -2 1 0 )\n(  3 2 0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Multiplying the lower element (2,1)
         lower(2,1) *= -3;

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 9UL );
         checkNonZeros( lower, 4UL );
         checkNonZeros( lower, 0UL, 2UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 0UL );

         if( lower(0,0) !=  0 || lower(0,1) !=  0 || lower(0,2) != 0 ||
             lower(1,0) != -2 || lower(1,1) !=  1 || lower(1,2) != 0 ||
             lower(2,0) !=  3 || lower(2,1) != -6 || lower(2,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n(  0  0  0 )\n( -2  1  0 )\n(  3 -6  0 )\n";
            throw std::runtime_error( oss.str() );
         }

         // Dividing the lower element (2,1)
         lower(2,1) /= 2;

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 9UL );
         checkNonZeros( lower, 4UL );
         checkNonZeros( lower, 0UL, 2UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 0UL );

         if( lower(0,0) !=  0 || lower(0,1) !=  0 || lower(0,2) != 0 ||
             lower(1,0) != -2 || lower(1,1) !=  1 || lower(1,2) != 0 ||
             lower(2,0) !=  3 || lower(2,1) != -3 || lower(2,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Function call operator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n(  0  0  0 )\n( -2  1  0 )\n(  3 -3  0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Failure cases
      {
         OLT lower( 3UL );

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
/*!\brief Test of the LowerMatrix iterator implementation.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the iterator implementation of the LowerMatrix
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
      lower(0,0) =  1;
      lower(1,1) = -2;
      lower(2,0) =  3;
      lower(2,2) =  4;

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

      // Testing assignment to lower elements via Iterator
      {
         test_ = "Row-major assignment to lower elements via Iterator";

         int value = 0;

         for( Iterator it=begin( lower, 2UL ); it!=end( lower, 2UL ); ++it ) {
            *it = value++;
         }

         if( lower(0,0) != 1 || lower(0,1) !=  0 || lower(0,2) != 0 ||
             lower(1,0) != 0 || lower(1,1) != -2 || lower(1,2) != 0 ||
             lower(2,0) != 0 || lower(2,1) !=  1 || lower(2,2) != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1  0  0 )\n( 0 -2  0 )\n( 0  1  2 )\n";
            throw std::runtime_error( oss.str() );
         }
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

         int value = -2;

         for( Iterator it=begin( lower, 2UL ); it!=end( lower, 2UL ); ++it ) {
            *it += value++;
         }

         if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) != 0 ||
             lower(1,0) !=  0 || lower(1,1) != -2 || lower(1,2) != 0 ||
             lower(2,0) != -2 || lower(2,1) !=  0 || lower(2,2) != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n(  1  0  0 )\n(  0 -2  0 )\n( -2  0  2 )\n";
            throw std::runtime_error( oss.str() );
         }
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

         int value = -2;

         for( Iterator it=begin( lower, 2UL ); it!=end( lower, 2UL ); ++it ) {
            *it -= value++;
         }

         if( lower(0,0) != 1 || lower(0,1) !=  0 || lower(0,2) != 0 ||
             lower(1,0) != 0 || lower(1,1) != -2 || lower(1,2) != 0 ||
             lower(2,0) != 0 || lower(2,1) !=  1 || lower(2,2) != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1  0  0 )\n( 0 -2  0 )\n( 0  1  2 )\n";
            throw std::runtime_error( oss.str() );
         }
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

         int value = 1;

         for( Iterator it=begin( lower, 2UL ); it!=end( lower, 2UL ); ++it ) {
            *it *= value++;
         }

         if( lower(0,0) != 1 || lower(0,1) !=  0 || lower(0,2) != 0 ||
             lower(1,0) != 0 || lower(1,1) != -2 || lower(1,2) != 0 ||
             lower(2,0) != 0 || lower(2,1) !=  2 || lower(2,2) != 6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1  0  0 )\n( 0 -2  0 )\n( 0  2  6 )\n";
            throw std::runtime_error( oss.str() );
         }
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

         for( Iterator it=begin( lower, 2UL ); it!=end( lower, 2UL ); ++it ) {
            *it /= 2;
         }

         if( lower(0,0) != 1 || lower(0,1) !=  0 || lower(0,2) != 0 ||
             lower(1,0) != 0 || lower(1,1) != -2 || lower(1,2) != 0 ||
             lower(2,0) != 0 || lower(2,1) !=  1 || lower(2,2) != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1  0  0 )\n( 0 -2  0 )\n( 0  1  3 )\n";
            throw std::runtime_error( oss.str() );
         }
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
      lower(0,0) =  1;
      lower(1,1) = -2;
      lower(2,0) =  3;
      lower(2,2) =  4;

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

         int value = 0;

         for( Iterator it=begin( lower, 0UL ); it!=end( lower, 0UL ); ++it ) {
            *it = value++;
         }

         if( lower(0,0) != 0 || lower(0,1) !=  0 || lower(0,2) != 0 ||
             lower(1,0) != 1 || lower(1,1) != -2 || lower(1,2) != 0 ||
             lower(2,0) != 2 || lower(2,1) !=  0 || lower(2,2) != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 0  0  0 )\n( 1 -2  0 )\n( 2  0  4 )\n";
            throw std::runtime_error( oss.str() );
         }
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

         int value = -2;

         for( Iterator it=begin( lower, 0UL ); it!=end( lower, 0UL ); ++it ) {
            *it += value++;
         }

         if( lower(0,0) != -2 || lower(0,1) !=  0 || lower(0,2) != 0 ||
             lower(1,0) !=  0 || lower(1,1) != -2 || lower(1,2) != 0 ||
             lower(2,0) !=  2 || lower(2,1) !=  0 || lower(2,2) != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( -2  0  0 )\n( -1 -2  0 )\n(  0  0  4 )\n";
            throw std::runtime_error( oss.str() );
         }
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

         int value = -2;

         for( Iterator it=begin( lower, 0UL ); it!=end( lower, 0UL ); ++it ) {
            *it -= value++;
         }

         if( lower(0,0) != 0 || lower(0,1) !=  0 || lower(0,2) != 0 ||
             lower(1,0) != 1 || lower(1,1) != -2 || lower(1,2) != 0 ||
             lower(2,0) != 2 || lower(2,1) !=  0 || lower(2,2) != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 0  0  0 )\n( 1 -2  0 )\n( 2  0  4 )\n";
            throw std::runtime_error( oss.str() );
         }
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

         int value = 1;

         for( Iterator it=begin( lower, 0UL ); it!=end( lower, 0UL ); ++it ) {
            *it *= value++;
         }

         if( lower(0,0) != 0 || lower(0,1) !=  0 || lower(0,2) != 0 ||
             lower(1,0) != 2 || lower(1,1) != -2 || lower(1,2) != 0 ||
             lower(2,0) != 6 || lower(2,1) !=  0 || lower(2,2) != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 0  0  0 )\n( 2 -2  0 )\n( 6  0  4 )\n";
            throw std::runtime_error( oss.str() );
         }
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

         for( Iterator it=begin( lower, 0UL ); it!=end( lower, 0UL ); ++it ) {
            *it /= 2;
         }

         if( lower(0,0) != 0 || lower(0,1) !=  0 || lower(0,2) != 0 ||
             lower(1,0) != 1 || lower(1,1) != -2 || lower(1,2) != 0 ||
             lower(2,0) != 3 || lower(2,1) !=  0 || lower(2,2) != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment via iterator failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 0  0  0 )\n( 1 -2  0 )\n( 3  0  4 )\n";
            throw std::runtime_error( oss.str() );
         }
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
/*!\brief Test of the \c nonZeros() member function of the LowerMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c nonZeros() member function of the LowerMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testNonZeros()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major LowerMatrix::nonZeros()";

      // Empty matrix
      {
         LT lower( 3UL );

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 9UL );
         checkNonZeros( lower, 0UL );
         checkNonZeros( lower, 0UL, 0UL );
         checkNonZeros( lower, 1UL, 0UL );
         checkNonZeros( lower, 2UL, 0UL );

         if( lower(0,0) != 0 || lower(0,1) != 0 || lower(0,2) != 0 ||
             lower(1,0) != 0 || lower(1,1) != 0 || lower(1,2) != 0 ||
             lower(2,0) != 0 || lower(2,1) != 0 || lower(2,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 0 0 0 )\n( 0 0 0 )\n( 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Partially filled matrix
      {
         LT lower( 3UL );
         lower(0,0) =  1;
         lower(1,1) = -2;
         lower(2,1) =  3;
         lower(2,2) = -4;

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 9UL );
         checkNonZeros( lower, 4UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 1UL );
         checkNonZeros( lower, 2UL, 2UL );

         if( lower(0,0) != 1 || lower(0,1) !=  0 || lower(0,2) !=  0 ||
             lower(1,0) != 0 || lower(1,1) != -2 || lower(1,2) !=  0 ||
             lower(2,0) != 0 || lower(2,1) !=  3 || lower(2,2) != -4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1  0  0 )\n( 0 -2  0 )\n( 0  3 -4 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Fully filled matrix
      {
         LT lower( 3UL );
         lower(0,0) = -1;
         lower(1,0) =  2;
         lower(1,1) =  3;
         lower(2,0) = -4;
         lower(2,1) = -5;
         lower(2,2) =  6;

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 9UL );
         checkNonZeros( lower, 6UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 3UL );

         if( lower(0,0) != -1 || lower(0,1) !=  0 || lower(0,2) != 0 ||
             lower(1,0) !=  2 || lower(1,1) !=  3 || lower(1,2) != 0 ||
             lower(2,0) != -4 || lower(2,1) != -5 || lower(2,2) != 6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( -1  0  0 )\n(  2  3  0 )\n( -4 -5  6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major LowerMatrix::nonZeros()";

      // Empty matrix
      {
         OLT lower( 3UL );

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 9UL );
         checkNonZeros( lower, 0UL );
         checkNonZeros( lower, 0UL, 0UL );
         checkNonZeros( lower, 1UL, 0UL );
         checkNonZeros( lower, 2UL, 0UL );

         if( lower(0,0) != 0 || lower(0,1) != 0 || lower(0,2) != 0 ||
             lower(1,0) != 0 || lower(1,1) != 0 || lower(1,2) != 0 ||
             lower(2,0) != 0 || lower(2,1) != 0 || lower(2,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 0 0 0 )\n( 0 0 0 )\n( 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Partially filled matrix
      {
         OLT lower( 3UL );
         lower(0,0) =  1;
         lower(1,1) = -2;
         lower(2,1) =  3;
         lower(2,2) = -4;

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 9UL );
         checkNonZeros( lower, 4UL );
         checkNonZeros( lower, 0UL, 1UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 1UL );

         if( lower(0,0) != 1 || lower(0,1) !=  0 || lower(0,2) !=  0 ||
             lower(1,0) != 0 || lower(1,1) != -2 || lower(1,2) !=  0 ||
             lower(2,0) != 0 || lower(2,1) !=  3 || lower(2,2) != -4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1  0  0 )\n( 0 -2  0 )\n( 0  3 -4 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Fully filled matrix
      {
         OLT lower( 3UL );
         lower(0,0) = -1;
         lower(1,0) =  2;
         lower(1,1) =  3;
         lower(2,0) = -4;
         lower(2,1) = -5;
         lower(2,2) =  6;

         checkRows    ( lower, 3UL );
         checkColumns ( lower, 3UL );
         checkCapacity( lower, 9UL );
         checkNonZeros( lower, 6UL );
         checkNonZeros( lower, 0UL, 3UL );
         checkNonZeros( lower, 1UL, 2UL );
         checkNonZeros( lower, 2UL, 1UL );

         if( lower(0,0) != -1 || lower(0,1) !=  0 || lower(0,2) != 0 ||
             lower(1,0) !=  2 || lower(1,1) !=  3 || lower(1,2) != 0 ||
             lower(2,0) != -4 || lower(2,1) != -5 || lower(2,2) != 6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( -1  0  0 )\n(  2  3  0 )\n( -4 -5  6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c reset() member function of the LowerMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c reset() member function of the LowerMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testReset()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major LowerMatrix::reset()";

      // Initialization check
      LT lower( 3UL );
      lower(0,0) = 1;
      lower(1,0) = 2;
      lower(1,1) = 3;
      lower(2,0) = 4;
      lower(2,1) = 5;
      lower(2,2) = 6;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 6UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 3UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 2 || lower(1,1) != 3 || lower(1,2) != 0 ||
          lower(2,0) != 4 || lower(2,1) != 5 || lower(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 2 3 0 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting a lower element
      reset( lower(1,0) );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 3UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 0 || lower(1,1) != 3 || lower(1,2) != 0 ||
          lower(2,0) != 4 || lower(2,1) != 5 || lower(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 3 0 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting an upper element
      reset( lower(0,1) );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 3UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 0 || lower(1,1) != 3 || lower(1,2) != 0 ||
          lower(2,0) != 4 || lower(2,1) != 5 || lower(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 3 0 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting row 1
      reset( lower, 1UL );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 4UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 0UL );
      checkNonZeros( lower, 2UL, 3UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 0 || lower(1,1) != 0 || lower(1,2) != 0 ||
          lower(2,0) != 4 || lower(2,1) != 5 || lower(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting the entire matrix
      reset( lower );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 0UL );
      checkNonZeros( lower, 0UL, 0UL );
      checkNonZeros( lower, 1UL, 0UL );
      checkNonZeros( lower, 2UL, 0UL );

      if( lower(0,0) != 0 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 0 || lower(1,1) != 0 || lower(1,2) != 0 ||
          lower(2,0) != 0 || lower(2,1) != 0 || lower(2,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 0 0 0 )\n( 0 0 0 )\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major LowerMatrix::reset()";

      // Initialization check
      OLT lower( 3UL );
      lower(0,0) = 1;
      lower(1,0) = 2;
      lower(1,1) = 3;
      lower(2,0) = 4;
      lower(2,1) = 5;
      lower(2,2) = 6;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 6UL );
      checkNonZeros( lower, 0UL, 3UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 2 || lower(1,1) != 3 || lower(1,2) != 0 ||
          lower(2,0) != 4 || lower(2,1) != 5 || lower(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 2 3 0 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting a lower element
      reset( lower(1,0) );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 2UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 0 || lower(1,1) != 3 || lower(1,2) != 0 ||
          lower(2,0) != 4 || lower(2,1) != 5 || lower(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 3 0 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting an upper element
      reset( lower(0,1) );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 2UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 0 || lower(1,1) != 3 || lower(1,2) != 0 ||
          lower(2,0) != 4 || lower(2,1) != 5 || lower(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 3 0 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting column 1
      reset( lower, 1UL );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 3UL );
      checkNonZeros( lower, 0UL, 2UL );
      checkNonZeros( lower, 1UL, 0UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 0 || lower(1,1) != 0 || lower(1,2) != 0 ||
          lower(2,0) != 4 || lower(2,1) != 0 || lower(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting the entire matrix
      reset( lower );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 0UL );
      checkNonZeros( lower, 0UL, 0UL );
      checkNonZeros( lower, 1UL, 0UL );
      checkNonZeros( lower, 2UL, 0UL );

      if( lower(0,0) != 0 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 0 || lower(1,1) != 0 || lower(1,2) != 0 ||
          lower(2,0) != 0 || lower(2,1) != 0 || lower(2,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 0 0 0 )\n( 0 0 0 )\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c clear() member function of the LowerMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c clear() member function of the LowerMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testClear()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major LowerMatrix::clear()";

      // Initialization check
      LT lower( 3UL );
      lower(0,0) = 1;
      lower(1,0) = 2;
      lower(1,1) = 3;
      lower(2,0) = 4;
      lower(2,1) = 5;
      lower(2,2) = 6;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 6UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 3UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 2 || lower(1,1) != 3 || lower(1,2) != 0 ||
          lower(2,0) != 4 || lower(2,1) != 5 || lower(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 2 3 0 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing a lower element
      clear( lower(1,0) );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 3UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 0 || lower(1,1) != 3 || lower(1,2) != 0 ||
          lower(2,0) != 4 || lower(2,1) != 5 || lower(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 3 0 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing an upper element
      clear( lower(0,1) );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 3UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 0 || lower(1,1) != 3 || lower(1,2) != 0 ||
          lower(2,0) != 4 || lower(2,1) != 5 || lower(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 3 0 )\n( 4 5 6 )\n";
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
      test_ = "Column-major LowerMatrix::clear()";

      // Initialization check
      OLT lower( 3UL );
      lower(0,0) = 1;
      lower(1,0) = 2;
      lower(1,1) = 3;
      lower(2,0) = 4;
      lower(2,1) = 5;
      lower(2,2) = 6;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 6UL );
      checkNonZeros( lower, 0UL, 3UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 2 || lower(1,1) != 3 || lower(1,2) != 0 ||
          lower(2,0) != 4 || lower(2,1) != 5 || lower(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 2 3 0 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing a lower element
      clear( lower(1,0) );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 2UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 0 || lower(1,1) != 3 || lower(1,2) != 0 ||
          lower(2,0) != 4 || lower(2,1) != 5 || lower(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 3 0 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing an upper element
      clear( lower(0,1) );

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkCapacity( lower, 9UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 2UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 0 || lower(1,1) != 3 || lower(1,2) != 0 ||
          lower(2,0) != 4 || lower(2,1) != 5 || lower(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 3 0 )\n( 4 5 6 )\n";
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
/*!\brief Test of the \c resize() member function of the LowerMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c resize() member function of the LowerMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testResize()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major LowerMatrix::resize()";

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

      if( lower(0,1) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( x 0 )\n( x x )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 4x4 and preserving the elements
      lower(0,0) = 1;
      lower(1,0) = 2;
      lower(1,1) = 3;
      lower.resize( 4UL, true );

      checkRows    ( lower,  4UL );
      checkColumns ( lower,  4UL );
      checkCapacity( lower, 16UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 || lower(0,3) != 0 ||
          lower(1,0) != 2 || lower(1,1) != 3 || lower(1,2) != 0 || lower(1,3) != 0 ||
                                                                   lower(2,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 2 3 0 0 )\n( x x x 0 )\n( x x x x )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 2x2
      lower(2,2) = 4;
      lower.resize( 2UL );

      checkRows    ( lower, 2UL );
      checkColumns ( lower, 2UL );
      checkCapacity( lower, 4UL );
      checkNonZeros( lower, 3UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 ||
          lower(1,0) != 2 || lower(1,1) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 )\n( 2 3 )\n";
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
      test_ = "Column-major LowerMatrix::resize()";

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

      if( lower(0,1) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( x 0 )\n( x x )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 4x4 and preserving the elements
      lower(0,0) = 1;
      lower(1,0) = 2;
      lower(1,1) = 3;
      lower.resize( 4UL, true );

      checkRows    ( lower,  4UL );
      checkColumns ( lower,  4UL );
      checkCapacity( lower, 16UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 || lower(0,3) != 0 ||
          lower(1,0) != 2 || lower(1,1) != 3 || lower(1,2) != 0 || lower(1,3) != 0 ||
                                                                   lower(2,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 2 3 0 0 )\n( x x x 0 )\n( x x x x )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 2x2
      lower(2,2) = 4;
      lower.resize( 2UL );

      checkRows    ( lower, 2UL );
      checkColumns ( lower, 2UL );
      checkCapacity( lower, 4UL );
      checkNonZeros( lower, 3UL );
      checkNonZeros( lower, 0UL, 2UL );
      checkNonZeros( lower, 1UL, 1UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 ||
          lower(1,0) != 2 || lower(1,1) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 )\n( 2 3 )\n";
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
/*!\brief Test of the \c extend() member function of the LowerMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c extend() member function of the LowerMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testExtend()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major LowerMatrix::extend()";

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

      if( lower(0,1) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Extending the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( x 0 )\n( x x )\n";
         throw std::runtime_error( oss.str() );
      }

      // Extending to 4x4 and preserving the elements
      lower(0,0) = 1;
      lower(1,0) = 2;
      lower(1,1) = 3;
      lower.extend( 2UL, true );

      checkRows    ( lower,  4UL );
      checkColumns ( lower,  4UL );
      checkCapacity( lower, 16UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 || lower(0,3) != 0 ||
          lower(1,0) != 2 || lower(1,1) != 3 || lower(1,2) != 0 || lower(1,3) != 0 ||
                                                                   lower(2,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Extending the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 2 3 0 0 )\n( x x x 0 )\n( x x x x )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major LowerMatrix::extend()";

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

      if( lower(0,1) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Extending the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( x 0 )\n( x x )\n";
         throw std::runtime_error( oss.str() );
      }

      // Extending to 4x4 and preserving the elements
      lower(0,0) = 1;
      lower(1,0) = 2;
      lower(1,1) = 3;
      lower.extend( 2UL, true );

      checkRows    ( lower,  4UL );
      checkColumns ( lower,  4UL );
      checkCapacity( lower, 16UL );

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 || lower(0,3) != 0 ||
          lower(1,0) != 2 || lower(1,1) != 3 || lower(1,2) != 0 || lower(1,3) != 0 ||
                                                                   lower(2,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Extending the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 0 )\n( 2 3 0 0 )\n( x x x 0 )\n( x x x x )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c reserve() member function of the LowerMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c reserve() member function of the LowerMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testReserve()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major LowerMatrix::reserve()";

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
      test_ = "Column-major LowerMatrix::reserve()";

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
/*!\brief Test of the \c shrinkToFit() member function of the LowerMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c shrinkToFit() member function of the LowerMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testShrinkToFit()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major LowerMatrix::shrinkToFit()";

      // Shrinking a matrix without excessive capacity
      {
         LT lower{ { 1, 0, 0 }, { 2, 3, 0 }, { 4, 5, 6 } };

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
             lower(1,0) != 2 || lower(1,1) != 3 || lower(1,2) != 0 ||
             lower(2,0) != 4 || lower(2,1) != 5 || lower(2,2) != 6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Shrinking the matrix failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 )\n( 2 3 0 )\n( 4 5 6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Shrinking a matrix with excessive capacity
      {
         LT lower{ { 1, 0, 0 }, { 2, 3, 0 }, { 4, 5, 6 } };
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
             lower(1,0) != 2 || lower(1,1) != 3 || lower(1,2) != 0 ||
             lower(2,0) != 4 || lower(2,1) != 5 || lower(2,2) != 6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Shrinking the matrix failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 )\n( 2 3 0 )\n( 4 5 6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major LowerMatrix::shrinkToFit()";

      // Shrinking a matrix without excessive capacity
      {
         OLT lower{ { 1, 0, 0 }, { 2, 3, 0 }, { 4, 5, 6 } };

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
             lower(1,0) != 2 || lower(1,1) != 3 || lower(1,2) != 0 ||
             lower(2,0) != 4 || lower(2,1) != 5 || lower(2,2) != 6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Shrinking the matrix failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 )\n( 2 3 0 )\n( 4 5 6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Shrinking a matrix with excessive capacity
      {
         OLT lower{ { 1, 0, 0 }, { 2, 3, 0 }, { 4, 5, 6 } };
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
             lower(1,0) != 2 || lower(1,1) != 3 || lower(1,2) != 0 ||
             lower(2,0) != 4 || lower(2,1) != 5 || lower(2,2) != 6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Shrinking the matrix failed\n"
                << " Details:\n"
                << "   Result:\n" << lower << "\n"
                << "   Expected result:\n( 1 0 0 )\n( 2 3 0 )\n( 4 5 6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c swap() functionality of the LowerMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c swap() function of the LowerMatrix specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testSwap()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major LowerMatrix swap";

      LT lower1( 2UL );
      lower1(0,0) = 1;
      lower1(1,0) = 2;
      lower1(1,1) = 3;

      LT lower2( 2UL );
      lower2(0,0) = 4;
      lower2(1,0) = 5;
      lower2(1,1) = 0;

      swap( lower1, lower2 );

      checkRows    ( lower1, 2UL );
      checkColumns ( lower1, 2UL );
      checkCapacity( lower1, 4UL );
      checkNonZeros( lower1, 2UL );
      checkNonZeros( lower1, 0UL, 1UL );
      checkNonZeros( lower1, 1UL, 1UL );

      if( lower1(0,0) != 4 || lower1(0,1) != 0 || lower1(1,0) != 5 || lower1(1,1) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the first matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower1 << "\n"
             << "   Expected result:\n( 4 0 )\n( 5 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      checkRows    ( lower2, 2UL );
      checkColumns ( lower2, 2UL );
      checkCapacity( lower2, 4UL );
      checkNonZeros( lower2, 3UL );
      checkNonZeros( lower2, 0UL, 1UL );
      checkNonZeros( lower2, 1UL, 2UL );

      if( lower2(0,0) != 1 || lower2(0,1) != 0 || lower2(1,0) != 2 || lower2(1,1) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the second matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n( 1 0 )\n( 2 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major LowerMatrix swap";

      OLT lower1( 2UL );
      lower1(0,0) = 1;
      lower1(1,0) = 2;
      lower1(1,1) = 3;

      OLT lower2( 2UL );
      lower2(0,0) = 4;
      lower2(1,0) = 5;
      lower2(1,1) = 0;

      swap( lower1, lower2 );

      checkRows    ( lower1, 2UL );
      checkColumns ( lower1, 2UL );
      checkCapacity( lower1, 4UL );
      checkNonZeros( lower1, 2UL );
      checkNonZeros( lower1, 0UL, 2UL );
      checkNonZeros( lower1, 1UL, 0UL );

      if( lower1(0,0) != 4 || lower1(0,1) != 0 || lower1(1,0) != 5 || lower1(1,1) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the first matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower1 << "\n"
             << "   Expected result:\n( 4 0 )\n( 5 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      checkRows    ( lower2, 2UL );
      checkColumns ( lower2, 2UL );
      checkCapacity( lower2, 4UL );
      checkNonZeros( lower2, 3UL );
      checkNonZeros( lower2, 0UL, 2UL );
      checkNonZeros( lower2, 1UL, 1UL );

      if( lower2(0,0) != 1 || lower2(0,1) != 0 || lower2(1,0) != 2 || lower2(1,1) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the second matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower2 << "\n"
             << "   Expected result:\n( 1 0 )\n( 2 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c isDefault() function with the LowerMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isDefault() function with the LowerMatrix
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
      test_ = "Column-major isDefault() function";

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
/*!\brief Test of the assignment to submatrices of the LowerMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the assignment to submatrices of the LowerMatrix
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
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

      SMT sm = submatrix( lower, 1UL, 1UL, 2UL, 2UL );

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

      sm(1,0) = -5;

      if( sm(0,0) !=  2 || sm(0,1) != 0 ||
          sm(1,0) != -5 || sm(1,1) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix access failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  2  0 )\n( -5  3 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) != 0 ||
          lower(1,0) != -4 || lower(1,1) !=  2 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != -5 || lower(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix access failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0 )\n( -4  2  0 )\n(  7 -5  3 )\n";
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

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -4 || lower(1,1) != 0 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 0 || lower(2,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix reset failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 0 0 )\n(  7 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major scalar assignment
   //=====================================================================================

   // (  1  0  0  0 )      (  1  0  0  0 )
   // ( -4  2  0  0 )  =>  ( -4 12  0  0 )
   // (  7  0  3  0 )      (  7 12 12  0 )
   // ( -2  0  1  5 )      ( -2 12 12  5 )
   {
      test_ = "Row-major submatrix() function (scalar assignment test 1)";

      using SMT = blaze::Submatrix<LT>;

      LT lower( 4UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;
      lower(3,0) = -2;
      lower(3,2) =  1;
      lower(3,3) =  5;

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
          sm(1,0) != 12 || sm(1,1) !=  0 ||
          sm(2,0) != 12 || sm(2,1) != 12 ||
          sm(3,0) != 12 || sm(3,1) != 12 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0  0 )\n( 12  0 )\n( 12 12 )\n( 12 12 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) !=  0 || lower(0,3) != 0 ||
          lower(1,0) != -4 || lower(1,1) != 12 || lower(1,2) !=  0 || lower(1,3) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 12 || lower(2,2) != 12 || lower(2,3) != 0 ||
          lower(3,0) != -2 || lower(3,1) != 12 || lower(3,2) != 12 || lower(3,3) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0  0 )\n"
                                     "( -4 12  0  0 )\n"
                                     "(  7 12 12  0 )\n"
                                     "( -2 12 12  5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // (  1  0  0  0 )      (  1  0  0  0 )
   // ( -4  2  0  0 )  =>  ( 12 12  0  0 )
   // (  7  0  3  0 )      ( 12 12 12  0 )
   // ( -2  0  1  5 )      ( -2  0  1  5 )
   {
      test_ = "Row-major submatrix() function (scalar assignment test 2)";

      using SMT = blaze::Submatrix<LT>;

      LT lower( 4UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;
      lower(3,0) = -2;
      lower(3,2) =  1;
      lower(3,3) =  5;

      SMT sm = submatrix( lower, 1UL, 0UL, 2UL, 4UL );
      sm = 12;

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkNonZeros( lower, 9UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 3UL );
      checkNonZeros( lower, 3UL, 3UL );

      if( sm(0,0) != 12 || sm(0,1) != 12 || sm(0,2) !=  0 || sm(0,3) != 0 ||
          sm(1,0) != 12 || sm(1,1) != 12 || sm(1,2) != 12 || sm(1,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 12 12  0  0 )\n( 12 12 12  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) !=  0 || lower(0,3) != 0 ||
          lower(1,0) != 12 || lower(1,1) != 12 || lower(1,2) !=  0 || lower(1,3) != 0 ||
          lower(2,0) != 12 || lower(2,1) != 12 || lower(2,2) != 12 || lower(2,3) != 0 ||
          lower(3,0) != -2 || lower(3,1) !=  0 || lower(3,2) !=  1 || lower(3,3) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0  0 )\n"
                                     "( 12 12  0  0 )\n"
                                     "( 12 12 12  0 )\n"
                                     "( -2  0  1  5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // (  1  0  0  0 )      (  1  0  0  0 )
   // ( -4  2  0  0 )  =>  ( -4  2  0  0 )
   // (  7  0  3  0 )      (  7  0  3  0 )
   // ( -2  0  1  5 )      ( -2  0  1  5 )
   {
      test_ = "Row-major submatrix() function (scalar assignment test 3)";

      using SMT = blaze::Submatrix<LT>;

      LT lower( 4UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;
      lower(3,0) = -2;
      lower(3,2) =  1;
      lower(3,3) =  5;

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
          lower(1,0) != -4 || lower(1,1) != 2 || lower(1,2) != 0 || lower(1,3) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 0 || lower(2,2) != 3 || lower(2,3) != 0 ||
          lower(3,0) != -2 || lower(3,1) != 0 || lower(3,2) != 1 || lower(3,3) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0  0 )\n"
                                     "( -4  2  0  0 )\n"
                                     "(  7  0  3  0 )\n"
                                     "( -2  0  1  5 )\n";
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
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

      SMT sm = submatrix( lower, 1UL, 1UL, 2UL, 2UL );

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

      sm(1,0) = -5;

      if( sm(0,0) !=  2 || sm(0,1) != 0 ||
          sm(1,0) != -5 || sm(1,1) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix access failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  2  0 )\n( -5  3 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) != 0 ||
          lower(1,0) != -4 || lower(1,1) !=  2 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != -5 || lower(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix access failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0 )\n( -4  2  0 )\n(  7 -5  3 )\n";
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

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -4 || lower(1,1) != 0 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 0 || lower(2,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix reset failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 0 0 )\n(  7 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major scalar assignment
   //=====================================================================================

   // (  1  0  0  0 )      (  1  0  0  0 )
   // ( -4  2  0  0 )  =>  ( -4 12  0  0 )
   // (  7  0  3  0 )      (  7 12 12  0 )
   // ( -2  0  1  5 )      ( -2 12 12  5 )
   {
      test_ = "Column-major submatrix() function (scalar assignment test 1)";

      using SMT = blaze::Submatrix<OLT>;

      OLT lower( 4UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;
      lower(3,0) = -2;
      lower(3,2) =  1;
      lower(3,3) =  5;

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
          sm(1,0) != 12 || sm(1,1) !=  0 ||
          sm(2,0) != 12 || sm(2,1) != 12 ||
          sm(3,0) != 12 || sm(3,1) != 12 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0  0 )\n( 12  0 )\n( 12 12 )\n( 12 12 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) !=  0 || lower(0,3) != 0 ||
          lower(1,0) != -4 || lower(1,1) != 12 || lower(1,2) !=  0 || lower(1,3) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 12 || lower(2,2) != 12 || lower(2,3) != 0 ||
          lower(3,0) != -2 || lower(3,1) != 12 || lower(3,2) != 12 || lower(3,3) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0  0 )\n"
                                     "( -4 12  0  0 )\n"
                                     "(  7 12 12  0 )\n"
                                     "( -2 12 12  5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // (  1  0  0  0 )      (  1  0  0  0 )
   // ( -4  2  0  0 )  =>  ( 12 12  0  0 )
   // (  7  0  3  0 )      ( 12 12 12  0 )
   // ( -2  0  1  5 )      ( -2  0  1  5 )
   {
      test_ = "Column-major submatrix() function (scalar assignment test 2)";

      using SMT = blaze::Submatrix<OLT>;

      OLT lower( 4UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;
      lower(3,0) = -2;
      lower(3,2) =  1;
      lower(3,3) =  5;

      SMT sm = submatrix( lower, 1UL, 0UL, 2UL, 4UL );
      sm = 12;

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkNonZeros( lower, 9UL );
      checkNonZeros( lower, 0UL, 4UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 2UL );
      checkNonZeros( lower, 3UL, 1UL );

      if( sm(0,0) != 12 || sm(0,1) != 12 || sm(0,2) !=  0 || sm(0,3) != 0 ||
          sm(1,0) != 12 || sm(1,1) != 12 || sm(1,2) != 12 || sm(1,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 12 12  0  0 )\n( 12 12 12  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) !=  0 || lower(0,3) != 0 ||
          lower(1,0) != 12 || lower(1,1) != 12 || lower(1,2) !=  0 || lower(1,3) != 0 ||
          lower(2,0) != 12 || lower(2,1) != 12 || lower(2,2) != 12 || lower(2,3) != 0 ||
          lower(3,0) != -2 || lower(3,1) !=  0 || lower(3,2) !=  1 || lower(3,3) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0  0 )\n"
                                     "( 12 12  0  0 )\n"
                                     "( 12 12 12  0 )\n"
                                     "( -2  0  1  5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // (  1  0  0  0 )      (  1  0  0  0 )
   // ( -4  2  0  0 )  =>  ( -4  2  0  0 )
   // (  7  0  3  0 )      (  7  0  3  0 )
   // ( -2  0  1  5 )      ( -2  0  1  5 )
   {
      test_ = "Column-major submatrix() function (scalar assignment test 3)";

      using SMT = blaze::Submatrix<OLT>;

      OLT lower( 4UL );
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;
      lower(3,0) = -2;
      lower(3,2) =  1;
      lower(3,3) =  5;

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
          lower(1,0) != -4 || lower(1,1) != 2 || lower(1,2) != 0 || lower(1,3) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 0 || lower(2,2) != 3 || lower(2,3) != 0 ||
          lower(3,0) != -2 || lower(3,1) != 0 || lower(3,2) != 1 || lower(3,3) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0  0 )\n"
                                     "( -4  2  0  0 )\n"
                                     "(  7  0  3  0 )\n"
                                     "( -2  0  1  5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the assignment to rows of the LowerMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the assignment to rows of the LowerMatrix specialization.
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
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

      RT row1 = row( lower, 1UL );

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

      row1[1] = -5;

      if( row1[0] != -4 || row1[1] != -5 || row1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row access failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( -4 -5  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) != 0 ||
          lower(1,0) != -4 || lower(1,1) != -5 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) !=  0 || lower(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row access failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
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

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 0 || lower(1,1) != 0 || lower(1,2) != 0 ||
          lower(2,0) != 7 || lower(2,1) != 0 || lower(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row reset failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 7 0 3 )\n";
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
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

      RT row1 = row( lower, 1UL );
      row1 = 8;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 2UL );

      if( row1[0] != 8 || row1[1] != 8 || row1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row access failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( 8 8 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 8 || lower(1,1) != 8 || lower(1,2) != 0 ||
          lower(2,0) != 7 || lower(2,1) != 0 || lower(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row access failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 8 8 0 )\n( 7 0 3 )\n";
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
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

      RT row1 = row( lower, 1UL );

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

      row1[1] = -5;

      if( row1[0] != -4 || row1[1] != -5 || row1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row access failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( -4 -5  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) != 0 ||
          lower(1,0) != -4 || lower(1,1) != -5 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) !=  0 || lower(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row access failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
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

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 0 || lower(1,1) != 0 || lower(1,2) != 0 ||
          lower(2,0) != 7 || lower(2,1) != 0 || lower(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row reset failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 7 0 3 )\n";
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
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

      RT row1 = row( lower, 1UL );
      row1 = 8;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkNonZeros( lower, 5UL );
      checkNonZeros( lower, 0UL, 3UL );
      checkNonZeros( lower, 1UL, 1UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( row1[0] != 8 || row1[1] != 8 || row1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row access failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( 8 8 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 8 || lower(1,1) != 8 || lower(1,2) != 0 ||
          lower(2,0) != 7 || lower(2,1) != 0 || lower(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row access failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 8 8 0 )\n( 7 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the assignment to columns of the LowerMatrix specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the assignment to columns of the LowerMatrix specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
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
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

      CT col1 = column( lower, 1UL );

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

      if( it == col1.end() || *it != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << *it << "\n"
             << "   Expected result: 0\n";
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

      if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) != 0 ||
          lower(1,0) != -4 || lower(1,1) != -5 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) !=  0 || lower(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column access failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0 )\n( -4 -5  0 )\n(  7  0  3 )\n";
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

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -4 || lower(1,1) != 0 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 0 || lower(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column reset failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 0 0 )\n(  7 0 3 )\n";
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
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

      CT col1 = column( lower, 1UL );
      col1 = 8;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkNonZeros( lower, 6UL );
      checkNonZeros( lower, 0UL, 1UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 3UL );

      if( col1[0] != 0 || col1[1] != 8 || col1[2] != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column access failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( 0 8 8 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -4 || lower(1,1) != 8 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 8 || lower(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column access failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 8 0 )\n(  7 8 3 )\n";
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
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

      CT col1 = column( lower, 1UL );

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

      if( it == col1.end() || *it != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << *it << "\n"
             << "   Expected result: 0\n";
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

      if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) != 0 ||
          lower(1,0) != -4 || lower(1,1) != -5 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) !=  0 || lower(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column access failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0 )\n( -4 -5  0 )\n(  7  0  3 )\n";
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

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -4 || lower(1,1) != 0 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 0 || lower(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column reset failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 0 0 )\n(  7 0 3 )\n";
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
      lower(0,0) =  1;
      lower(1,0) = -4;
      lower(1,1) =  2;
      lower(2,0) =  7;
      lower(2,2) =  3;

      CT col1 = column( lower, 1UL );
      col1 = 8;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkNonZeros( lower, 6UL );
      checkNonZeros( lower, 0UL, 3UL );
      checkNonZeros( lower, 1UL, 2UL );
      checkNonZeros( lower, 2UL, 1UL );

      if( col1[0] != 0 || col1[1] != 8 || col1[2] != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column access failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( 0 8 8 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -4 || lower(1,1) != 8 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 8 || lower(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Column access failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -4 8 0 )\n(  7 8 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************

} // namespace lowermatrix

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
   std::cout << "   Running LowerMatrix dense test (part 2)..." << std::endl;

   try
   {
      RUN_LOWERMATRIX_DENSE_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during LowerMatrix dense test (part 2):\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
