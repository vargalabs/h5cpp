//=================================================================================================
/*!
//  \file src/mathtest/matrices/uniformmatrix/ClassTest2.cpp
//  \brief Source file for the UniformMatrix class test (part 2)
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
#include <blaze/math/CompressedMatrix.h>
#include <blaze/math/CustomMatrix.h>
#include <blaze/math/DynamicMatrix.h>
#include <blaze/util/Complex.h>
#include <blaze/util/Memory.h>
#include <blaze/util/policies/Deallocate.h>
#include <blaze/util/Random.h>
#include <blazetest/mathtest/matrices/uniformmatrix/ClassTest.h>
#include <blazetest/mathtest/RandomMaximum.h>
#include <blazetest/mathtest/RandomMinimum.h>

#ifdef BLAZE_USE_HPX_THREADS
#  include <hpx/hpx_main.hpp>
#endif


namespace blazetest {

namespace mathtest {

namespace matrices {

namespace uniformmatrix {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the UniformMatrix class test.
//
// \exception std::runtime_error Operation error detected.
*/
ClassTest::ClassTest()
{
   testSchurAssign();
   testMultAssign();
   testScaling();
   testFunctionCall();
   testAt();
   testIterator();
   testNonZeros();
   testReset();
   testClear();
   testResize();
   testExtend();
   testSwap();
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
/*!\brief Test of the UniformMatrix Schur product assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the Schur product assignment operators of the UniformMatrix
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testSchurAssign()
{
   //=====================================================================================
   // Row-major dense matrix Schur product assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major UniformMatrix dense matrix Schur product assignment (mixed type)";

      blaze::UniformMatrix<short,blaze::rowMajor> mat1( 2UL, 3UL, 2 );

      blaze::UniformMatrix<int,blaze::rowMajor> mat2( 2UL, 3UL, 1 );

      mat2 %= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 3UL );
      checkNonZeros( mat2, 1UL, 3UL );

      if( mat2(0,0) != 2 || mat2(0,1) != 2 || mat2(0,2) != 2 ||
          mat2(1,0) != 2 || mat2(1,1) != 2 || mat2(1,2) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major UniformMatrix dense matrix Schur product assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::rowMajor;

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,rowMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 32UL ) );
      AlignedPadded mat1( memory.get(), 2UL, 3UL, 16UL );
      mat1 = 2;

      blaze::UniformMatrix<int,blaze::rowMajor> mat2( 2UL, 3UL, 1 );

      mat2 %= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 3UL );
      checkNonZeros( mat2, 1UL, 3UL );

      if( mat2(0,0) != 2 || mat2(0,1) != 2 || mat2(0,2) != 2 ||
          mat2(1,0) != 2 || mat2(1,1) != 2 || mat2(1,2) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major UniformMatrix dense matrix Schur product assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowMajor;

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>;
      std::unique_ptr<int[]> memory( new int[7UL] );
      UnalignedUnpadded mat1( memory.get()+1UL, 2UL, 3UL );
      mat1 = 2;

      blaze::UniformMatrix<int,blaze::rowMajor> mat2( 2UL, 3UL, 1 );

      mat2 %= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 3UL );
      checkNonZeros( mat2, 1UL, 3UL );

      if( mat2(0,0) != 2 || mat2(0,1) != 2 || mat2(0,2) != 2 ||
          mat2(1,0) != 2 || mat2(1,1) != 2 || mat2(1,2) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major UniformMatrix dense matrix Schur product assignment (non-uniform)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat1{ { 2, 2, 2 }, { 2, 0, 2 } };

      try {
         blaze::UniformMatrix<int,blaze::rowMajor> mat2( 2UL, 3UL, 1 );
         mat2 %= mat1;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-uniform dense matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   {
      test_ = "Row-major/column-major UniformMatrix dense matrix Schur product assignment (mixed type)";

      blaze::UniformMatrix<short,blaze::columnMajor> mat1( 2UL, 3UL, 2 );

      blaze::UniformMatrix<int,blaze::rowMajor> mat2( 2UL, 3UL, 1 );

      mat2 %= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 3UL );
      checkNonZeros( mat2, 1UL, 3UL );

      if( mat2(0,0) != 2 || mat2(0,1) != 2 || mat2(0,2) != 2 ||
          mat2(1,0) != 2 || mat2(1,1) != 2 || mat2(1,2) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major UniformMatrix dense matrix Schur product assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::columnMajor;

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,columnMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      AlignedPadded mat1( memory.get(), 2UL, 3UL, 16UL );
      mat1 = 2;

      blaze::UniformMatrix<int,blaze::rowMajor> mat2( 2UL, 3UL, 1 );

      mat2 %= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 3UL );
      checkNonZeros( mat2, 1UL, 3UL );

      if( mat2(0,0) != 2 || mat2(0,1) != 2 || mat2(0,2) != 2 ||
          mat2(1,0) != 2 || mat2(1,1) != 2 || mat2(1,2) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major UniformMatrix dense matrix Schur product assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::columnMajor;

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>;
      std::unique_ptr<int[]> memory( new int[7UL] );
      UnalignedUnpadded mat1( memory.get()+1UL, 2UL, 3UL );
      mat1 = 2;

      blaze::UniformMatrix<int,blaze::rowMajor> mat2( 2UL, 3UL, 1 );

      mat2 %= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 3UL );
      checkNonZeros( mat2, 1UL, 3UL );

      if( mat2(0,0) != 2 || mat2(0,1) != 2 || mat2(0,2) != 2 ||
          mat2(1,0) != 2 || mat2(1,1) != 2 || mat2(1,2) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major UniformMatrix dense matrix Schur product assignment (non-uniform)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat1{ { 2, 2, 2 }, { 2, 0, 2 } };

      try {
         blaze::UniformMatrix<int,blaze::rowMajor> mat2( 2UL, 3UL, 1 );
         mat2 %= mat1;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-uniform dense matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Row-major sparse matrix Schur product assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major UniformMatrix sparse matrix Schur product assignment";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1{ { 2, 2, 2 }, { 2, 2, 2 } };

      blaze::UniformMatrix<int,blaze::rowMajor> mat2( 2UL, 3UL, 1 );

      mat2 %= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 3UL );
      checkNonZeros( mat2, 1UL, 3UL );

      if( mat2(0,0) != 2 || mat2(0,1) != 2 || mat2(0,2) != 2 ||
          mat2(1,0) != 2 || mat2(1,1) != 2 || mat2(1,2) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major UniformMatrix sparse matrix Schur product assignment (non-uniform)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1{ { 2, 2, 2 }, { 2, 0, 2 } };

      try {
         blaze::UniformMatrix<int,blaze::rowMajor> mat2( 2UL, 3UL, 1 );
         mat2 %= mat1;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-uniform sparse matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   {
      test_ = "Row-major/column-major UniformMatrix sparse matrix Schur product assignment";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1{ { 2, 2, 2 }, { 2, 2, 2 } };

      blaze::UniformMatrix<int,blaze::rowMajor> mat2( 2UL, 3UL, 1 );

      mat2 %= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 3UL );
      checkNonZeros( mat2, 1UL, 3UL );

      if( mat2(0,0) != 2 || mat2(0,1) != 2 || mat2(0,2) != 2 ||
          mat2(1,0) != 2 || mat2(1,1) != 2 || mat2(1,2) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major UniformMatrix sparse matrix Schur product assignment (non-uniform)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1{ { 2, 2, 2 }, { 2, 0, 2 } };

      try {
         blaze::UniformMatrix<int,blaze::rowMajor> mat2( 2UL, 3UL, 1 );
         mat2 %= mat1;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-uniform sparse matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Column-major dense matrix Schur product assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major UniformMatrix dense matrix Schur product assignment (mixed type)";

      blaze::UniformMatrix<short,blaze::rowMajor> mat1( 2UL, 3UL, 2 );

      blaze::UniformMatrix<int,blaze::columnMajor> mat2( 2UL, 3UL, 1 );

      mat2 %= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 2 || mat2(0,1) != 2 || mat2(0,2) != 2 ||
          mat2(1,0) != 2 || mat2(1,1) != 2 || mat2(1,2) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major UniformMatrix dense matrix Schur product assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::rowMajor;

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,rowMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 32UL ) );
      AlignedPadded mat1( memory.get(), 2UL, 3UL, 16UL );
      mat1 = 2;

      blaze::UniformMatrix<int,blaze::columnMajor> mat2( 2UL, 3UL, 1 );

      mat2 %= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 2 || mat2(0,1) != 2 || mat2(0,2) != 2 ||
          mat2(1,0) != 2 || mat2(1,1) != 2 || mat2(1,2) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major UniformMatrix dense matrix Schur product assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowMajor;

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>;
      std::unique_ptr<int[]> memory( new int[7UL] );
      UnalignedUnpadded mat1( memory.get()+1UL, 2UL, 3UL );
      mat1 = 2;

      blaze::UniformMatrix<int,blaze::columnMajor> mat2( 2UL, 3UL, 1 );

      mat2 %= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 2 || mat2(0,1) != 2 || mat2(0,2) != 2 ||
          mat2(1,0) != 2 || mat2(1,1) != 2 || mat2(1,2) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major UniformMatrix dense matrix Schur product assignment (non-uniform)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat1{ { 2, 2, 2 }, { 2, 0, 2 } };

      try {
         blaze::UniformMatrix<int,blaze::columnMajor> mat2( 2UL, 3UL, 1 );
         mat2 %= mat1;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-uniform dense matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   {
      test_ = "Column-major/column-major UniformMatrix dense matrix Schur product assignment (mixed type)";

      blaze::UniformMatrix<short,blaze::columnMajor> mat1( 2UL, 3UL, 2 );

      blaze::UniformMatrix<int,blaze::columnMajor> mat2( 2UL, 3UL, 1 );

      mat2 %= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 2 || mat2(0,1) != 2 || mat2(0,2) != 2 ||
          mat2(1,0) != 2 || mat2(1,1) != 2 || mat2(1,2) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major UniformMatrix dense matrix Schur product assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::columnMajor;

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,columnMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      AlignedPadded mat1( memory.get(), 2UL, 3UL, 16UL );
      mat1 = 2;

      blaze::UniformMatrix<int,blaze::columnMajor> mat2( 2UL, 3UL, 1 );

      mat2 %= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 2 || mat2(0,1) != 2 || mat2(0,2) != 2 ||
          mat2(1,0) != 2 || mat2(1,1) != 2 || mat2(1,2) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major UniformMatrix dense matrix Schur product assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::columnMajor;

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>;
      std::unique_ptr<int[]> memory( new int[7UL] );
      UnalignedUnpadded mat1( memory.get()+1UL, 2UL, 3UL );
      mat1 = 2;

      blaze::UniformMatrix<int,blaze::columnMajor> mat2( 2UL, 3UL, 1 );

      mat2 %= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 2 || mat2(0,1) != 2 || mat2(0,2) != 2 ||
          mat2(1,0) != 2 || mat2(1,1) != 2 || mat2(1,2) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major UniformMatrix dense matrix Schur product assignment (non-uniform)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat1{ { 2, 2, 2 }, { 2, 0, 2 } };

      try {
         blaze::UniformMatrix<int,blaze::rowMajor> mat2( 2UL, 3UL, 1 );
         mat2 %= mat1;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-uniform dense matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Column-major sparse matrix Schur product assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major UniformMatrix sparse matrix Schur product assignment";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1{ { 2, 2, 2 }, { 2, 2, 2 } };

      blaze::UniformMatrix<int,blaze::columnMajor> mat2( 2UL, 3UL, 1 );

      mat2 %= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 2 || mat2(0,1) != 2 || mat2(0,2) != 2 ||
          mat2(1,0) != 2 || mat2(1,1) != 2 || mat2(1,2) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major UniformMatrix sparse matrix Schur product assignment (non-uniform)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1{ { 2, 2, 2 }, { 2, 0, 2 } };

      try {
         blaze::UniformMatrix<int,blaze::columnMajor> mat2( 2UL, 3UL, 1 );
         mat2 %= mat1;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-uniform sparse matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   {
      test_ = "Column-major/column-major UniformMatrix sparse matrix Schur product assignment";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1{ { 2, 2, 2 }, { 2, 2, 2 } };

      blaze::UniformMatrix<int,blaze::columnMajor> mat2( 2UL, 3UL, 1 );

      mat2 %= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 2 || mat2(0,1) != 2 || mat2(0,2) != 2 ||
          mat2(1,0) != 2 || mat2(1,1) != 2 || mat2(1,2) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major UniformMatrix sparse matrix Schur product assignment (non-uniform)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1{ { 2, 2, 2 }, { 2, 0, 2 } };

      try {
         blaze::UniformMatrix<int,blaze::columnMajor> mat2( 2UL, 3UL, 1 );
         mat2 %= mat1;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-uniform sparse matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the UniformMatrix multiplication assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the multiplication assignment operators of the UniformMatrix
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testMultAssign()
{
   //=====================================================================================
   // Row-major dense matrix multiplication assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major UniformMatrix dense matrix multiplication assignment (mixed type)";

      blaze::UniformMatrix<short,blaze::rowMajor> mat1( 3U, 4UL, 2 );
      blaze::UniformMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL, 1 );

      mat2 *= mat1;

      checkRows    ( mat2,  3UL );
      checkColumns ( mat2,  4UL );
      checkNonZeros( mat2, 12UL );
      checkNonZeros( mat2, 0UL, 4UL );
      checkNonZeros( mat2, 1UL, 4UL );
      checkNonZeros( mat2, 2UL, 4UL );

      if( mat2(0,0) != 6 || mat2(0,1) != 6 || mat2(0,2) != 6 || mat2(0,3) != 6 ||
          mat2(1,0) != 6 || mat2(1,1) != 6 || mat2(1,2) != 6 || mat2(1,3) != 6 ||
          mat2(2,0) != 6 || mat2(2,1) != 6 || mat2(2,2) != 6 || mat2(2,3) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 6 6 6 6 )\n( 6 6 6 6 )\n( 6 6 6 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major UniformMatrix dense matrix multiplication assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::rowMajor;

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,rowMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      AlignedPadded mat1( memory.get(), 3UL, 4UL, 16UL );
      mat1 = 2;

      blaze::UniformMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL, 1 );

      mat2 *= mat1;

      checkRows    ( mat2,  3UL );
      checkColumns ( mat2,  4UL );
      checkNonZeros( mat2, 12UL );
      checkNonZeros( mat2, 0UL, 4UL );
      checkNonZeros( mat2, 1UL, 4UL );
      checkNonZeros( mat2, 2UL, 4UL );

      if( mat2(0,0) != 6 || mat2(0,1) != 6 || mat2(0,2) != 6 || mat2(0,3) != 6 ||
          mat2(1,0) != 6 || mat2(1,1) != 6 || mat2(1,2) != 6 || mat2(1,3) != 6 ||
          mat2(2,0) != 6 || mat2(2,1) != 6 || mat2(2,2) != 6 || mat2(2,3) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 6 6 6 6 )\n( 6 6 6 6 )\n( 6 6 6 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major UniformMatrix dense matrix multiplication assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowMajor;

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>;
      std::unique_ptr<int[]> memory( new int[13UL] );
      UnalignedUnpadded mat1( memory.get()+1UL, 3UL, 4UL );
      mat1 = 2;

      blaze::UniformMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL, 1 );

      mat2 *= mat1;

      checkRows    ( mat2,  3UL );
      checkColumns ( mat2,  4UL );
      checkNonZeros( mat2, 12UL );
      checkNonZeros( mat2, 0UL, 4UL );
      checkNonZeros( mat2, 1UL, 4UL );
      checkNonZeros( mat2, 2UL, 4UL );

      if( mat2(0,0) != 6 || mat2(0,1) != 6 || mat2(0,2) != 6 || mat2(0,3) != 6 ||
          mat2(1,0) != 6 || mat2(1,1) != 6 || mat2(1,2) != 6 || mat2(1,3) != 6 ||
          mat2(2,0) != 6 || mat2(2,1) != 6 || mat2(2,2) != 6 || mat2(2,3) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 6 6 6 6 )\n( 6 6 6 6 )\n( 6 6 6 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major UniformMatrix dense matrix multiplication assignment (non-uniform)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat1{ { 2, 2, 2 }, { 2, 0, 2 }, { 2, 2, 2 } };

      try {
         blaze::UniformMatrix<int,blaze::rowMajor> mat2( 2UL, 3UL, 1 );
         mat2 *= mat1;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-uniform dense matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   {
      test_ = "Row-major/column-major UniformMatrix dense matrix multiplication assignment (mixed type)";

      blaze::UniformMatrix<short,blaze::columnMajor> mat1( 3U, 4UL, 2 );
      blaze::UniformMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL, 1 );

      mat2 *= mat1;

      checkRows    ( mat2,  3UL );
      checkColumns ( mat2,  4UL );
      checkNonZeros( mat2, 12UL );
      checkNonZeros( mat2, 0UL, 4UL );
      checkNonZeros( mat2, 1UL, 4UL );
      checkNonZeros( mat2, 2UL, 4UL );

      if( mat2(0,0) != 6 || mat2(0,1) != 6 || mat2(0,2) != 6 || mat2(0,3) != 6 ||
          mat2(1,0) != 6 || mat2(1,1) != 6 || mat2(1,2) != 6 || mat2(1,3) != 6 ||
          mat2(2,0) != 6 || mat2(2,1) != 6 || mat2(2,2) != 6 || mat2(2,3) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 6 6 6 6 )\n( 6 6 6 6 )\n( 6 6 6 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major UniformMatrix dense matrix multiplication assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::columnMajor;

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,columnMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 64UL ) );
      AlignedPadded mat1( memory.get(), 3UL, 4UL, 16UL );
      mat1 = 2;

      blaze::UniformMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL, 1 );

      mat2 *= mat1;

      checkRows    ( mat2,  3UL );
      checkColumns ( mat2,  4UL );
      checkNonZeros( mat2, 12UL );
      checkNonZeros( mat2, 0UL, 4UL );
      checkNonZeros( mat2, 1UL, 4UL );
      checkNonZeros( mat2, 2UL, 4UL );

      if( mat2(0,0) != 6 || mat2(0,1) != 6 || mat2(0,2) != 6 || mat2(0,3) != 6 ||
          mat2(1,0) != 6 || mat2(1,1) != 6 || mat2(1,2) != 6 || mat2(1,3) != 6 ||
          mat2(2,0) != 6 || mat2(2,1) != 6 || mat2(2,2) != 6 || mat2(2,3) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 6 6 6 6 )\n( 6 6 6 6 )\n( 6 6 6 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major UniformMatrix dense matrix multiplication assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::columnMajor;

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>;
      std::unique_ptr<int[]> memory( new int[13UL] );
      UnalignedUnpadded mat1( memory.get()+1UL, 3UL, 4UL );
      mat1 = 2;

      blaze::UniformMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL, 1 );

      mat2 *= mat1;

      checkRows    ( mat2,  3UL );
      checkColumns ( mat2,  4UL );
      checkNonZeros( mat2, 12UL );
      checkNonZeros( mat2, 0UL, 4UL );
      checkNonZeros( mat2, 1UL, 4UL );
      checkNonZeros( mat2, 2UL, 4UL );

      if( mat2(0,0) != 6 || mat2(0,1) != 6 || mat2(0,2) != 6 || mat2(0,3) != 6 ||
          mat2(1,0) != 6 || mat2(1,1) != 6 || mat2(1,2) != 6 || mat2(1,3) != 6 ||
          mat2(2,0) != 6 || mat2(2,1) != 6 || mat2(2,2) != 6 || mat2(2,3) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 6 6 6 6 )\n( 6 6 6 6 )\n( 6 6 6 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major UniformMatrix dense matrix multiplication assignment (non-uniform)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat1{ { 2, 2, 2 }, { 2, 0, 2 }, { 2, 2, 2 } };

      try {
         blaze::UniformMatrix<int,blaze::rowMajor> mat2( 2UL, 3UL, 1 );
         mat2 *= mat1;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-uniform dense matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Row-major sparse matrix multiplication assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major UniformMatrix sparse matrix multiplication assignment";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1{ { 2, 2, 2, 2 }
                                                       , { 2, 2, 2, 2 }
                                                       , { 2, 2, 2, 2 } };

      blaze::UniformMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL, 1 );

      mat2 *= mat1;

      checkRows    ( mat2,  3UL );
      checkColumns ( mat2,  4UL );
      checkNonZeros( mat2, 12UL );
      checkNonZeros( mat2, 0UL, 4UL );
      checkNonZeros( mat2, 1UL, 4UL );
      checkNonZeros( mat2, 2UL, 4UL );

      if( mat2(0,0) != 6 || mat2(0,1) != 6 || mat2(0,2) != 6 || mat2(0,3) != 6 ||
          mat2(1,0) != 6 || mat2(1,1) != 6 || mat2(1,2) != 6 || mat2(1,3) != 6 ||
          mat2(2,0) != 6 || mat2(2,1) != 6 || mat2(2,2) != 6 || mat2(2,3) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 6 6 6 6 )\n( 6 6 6 6 )\n( 6 6 6 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major UniformMatrix sparse matrix multiplication assignment (non-uniform)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1{ { 2, 2, 2, 2 }
                                                       , { 2, 0, 2, 2 }
                                                       , { 2, 2, 2, 2 } };

      try {
         blaze::UniformMatrix<int,blaze::rowMajor> mat2( 2UL, 3UL, 1 );
         mat2 *= mat1;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-uniform sparse matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   {
      test_ = "Row-major/column-major UniformMatrix sparse matrix multiplication assignment";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1{ { 2, 2, 2, 2 }
                                                          , { 2, 2, 2, 2 }
                                                          , { 2, 2, 2, 2 } };

      blaze::UniformMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL, 1 );

      mat2 *= mat1;

      checkRows    ( mat2,  3UL );
      checkColumns ( mat2,  4UL );
      checkNonZeros( mat2, 12UL );
      checkNonZeros( mat2, 0UL, 4UL );
      checkNonZeros( mat2, 1UL, 4UL );
      checkNonZeros( mat2, 2UL, 4UL );

      if( mat2(0,0) != 6 || mat2(0,1) != 6 || mat2(0,2) != 6 || mat2(0,3) != 6 ||
          mat2(1,0) != 6 || mat2(1,1) != 6 || mat2(1,2) != 6 || mat2(1,3) != 6 ||
          mat2(2,0) != 6 || mat2(2,1) != 6 || mat2(2,2) != 6 || mat2(2,3) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 6 6 6 6 )\n( 6 6 6 6 )\n( 6 6 6 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major UniformMatrix sparse matrix multiplication assignment (non-uniform)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1{ { 2, 2, 2, 2 }
                                                          , { 2, 0, 2, 2 }
                                                          , { 2, 2, 2, 2 } };

      try {
         blaze::UniformMatrix<int,blaze::rowMajor> mat2( 2UL, 3UL, 1 );
         mat2 *= mat1;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-uniform sparse matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Column-major dense matrix multiplication assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major UniformMatrix dense matrix multiplication assignment (mixed type)";

      blaze::UniformMatrix<short,blaze::rowMajor> mat1( 3U, 4UL, 2 );
      blaze::UniformMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL, 1 );

      mat2 *= mat1;

      checkRows    ( mat2,  3UL );
      checkColumns ( mat2,  4UL );
      checkNonZeros( mat2, 12UL );
      checkNonZeros( mat2, 0UL, 3UL );
      checkNonZeros( mat2, 1UL, 3UL );
      checkNonZeros( mat2, 2UL, 3UL );
      checkNonZeros( mat2, 3UL, 3UL );

      if( mat2(0,0) != 6 || mat2(0,1) != 6 || mat2(0,2) != 6 || mat2(0,3) != 6 ||
          mat2(1,0) != 6 || mat2(1,1) != 6 || mat2(1,2) != 6 || mat2(1,3) != 6 ||
          mat2(2,0) != 6 || mat2(2,1) != 6 || mat2(2,2) != 6 || mat2(2,3) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 6 6 6 6 )\n( 6 6 6 6 )\n( 6 6 6 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major UniformMatrix dense matrix multiplication assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::rowMajor;

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,rowMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      AlignedPadded mat1( memory.get(), 3UL, 4UL, 16UL );
      mat1 = 2;

      blaze::UniformMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL, 1 );

      mat2 *= mat1;

      checkRows    ( mat2,  3UL );
      checkColumns ( mat2,  4UL );
      checkNonZeros( mat2, 12UL );
      checkNonZeros( mat2, 0UL, 3UL );
      checkNonZeros( mat2, 1UL, 3UL );
      checkNonZeros( mat2, 2UL, 3UL );
      checkNonZeros( mat2, 3UL, 3UL );

      if( mat2(0,0) != 6 || mat2(0,1) != 6 || mat2(0,2) != 6 || mat2(0,3) != 6 ||
          mat2(1,0) != 6 || mat2(1,1) != 6 || mat2(1,2) != 6 || mat2(1,3) != 6 ||
          mat2(2,0) != 6 || mat2(2,1) != 6 || mat2(2,2) != 6 || mat2(2,3) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 6 6 6 6 )\n( 6 6 6 6 )\n( 6 6 6 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major UniformMatrix dense matrix multiplication assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowMajor;

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>;
      std::unique_ptr<int[]> memory( new int[13UL] );
      UnalignedUnpadded mat1( memory.get()+1UL, 3UL, 4UL );
      mat1 = 2;

      blaze::UniformMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL, 1 );

      mat2 *= mat1;

      checkRows    ( mat2,  3UL );
      checkColumns ( mat2,  4UL );
      checkNonZeros( mat2, 12UL );
      checkNonZeros( mat2, 0UL, 3UL );
      checkNonZeros( mat2, 1UL, 3UL );
      checkNonZeros( mat2, 2UL, 3UL );
      checkNonZeros( mat2, 3UL, 3UL );

      if( mat2(0,0) != 6 || mat2(0,1) != 6 || mat2(0,2) != 6 || mat2(0,3) != 6 ||
          mat2(1,0) != 6 || mat2(1,1) != 6 || mat2(1,2) != 6 || mat2(1,3) != 6 ||
          mat2(2,0) != 6 || mat2(2,1) != 6 || mat2(2,2) != 6 || mat2(2,3) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 6 6 6 6 )\n( 6 6 6 6 )\n( 6 6 6 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major UniformMatrix dense matrix multiplication assignment (non-uniform)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat1{ { 2, 2, 2 }, { 2, 0, 2 }, { 2, 2, 2 } };

      try {
         blaze::UniformMatrix<int,blaze::columnMajor> mat2( 2UL, 3UL, 1 );
         mat2 *= mat1;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-uniform dense matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   {
      test_ = "Column-major/column-major UniformMatrix dense matrix multiplication assignment (mixed type)";

      blaze::UniformMatrix<short,blaze::columnMajor> mat1( 3U, 4UL, 2 );
      blaze::UniformMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL, 1 );

      mat2 *= mat1;

      checkRows    ( mat2,  3UL );
      checkColumns ( mat2,  4UL );
      checkNonZeros( mat2, 12UL );
      checkNonZeros( mat2, 0UL, 3UL );
      checkNonZeros( mat2, 1UL, 3UL );
      checkNonZeros( mat2, 2UL, 3UL );
      checkNonZeros( mat2, 3UL, 3UL );

      if( mat2(0,0) != 6 || mat2(0,1) != 6 || mat2(0,2) != 6 || mat2(0,3) != 6 ||
          mat2(1,0) != 6 || mat2(1,1) != 6 || mat2(1,2) != 6 || mat2(1,3) != 6 ||
          mat2(2,0) != 6 || mat2(2,1) != 6 || mat2(2,2) != 6 || mat2(2,3) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 6 6 6 6 )\n( 6 6 6 6 )\n( 6 6 6 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major UniformMatrix dense matrix multiplication assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::columnMajor;

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,columnMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 64UL ) );
      AlignedPadded mat1( memory.get(), 3UL, 4UL, 16UL );
      mat1 = 2;

      blaze::UniformMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL, 1 );

      mat2 *= mat1;

      checkRows    ( mat2,  3UL );
      checkColumns ( mat2,  4UL );
      checkNonZeros( mat2, 12UL );
      checkNonZeros( mat2, 0UL, 3UL );
      checkNonZeros( mat2, 1UL, 3UL );
      checkNonZeros( mat2, 2UL, 3UL );
      checkNonZeros( mat2, 3UL, 3UL );

      if( mat2(0,0) != 6 || mat2(0,1) != 6 || mat2(0,2) != 6 || mat2(0,3) != 6 ||
          mat2(1,0) != 6 || mat2(1,1) != 6 || mat2(1,2) != 6 || mat2(1,3) != 6 ||
          mat2(2,0) != 6 || mat2(2,1) != 6 || mat2(2,2) != 6 || mat2(2,3) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 6 6 6 6 )\n( 6 6 6 6 )\n( 6 6 6 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major UniformMatrix dense matrix multiplication assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::columnMajor;

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>;
      std::unique_ptr<int[]> memory( new int[13UL] );
      UnalignedUnpadded mat1( memory.get()+1UL, 3UL, 4UL );
      mat1 = 2;

      blaze::UniformMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL, 1 );

      mat2 *= mat1;

      checkRows    ( mat2,  3UL );
      checkColumns ( mat2,  4UL );
      checkNonZeros( mat2, 12UL );
      checkNonZeros( mat2, 0UL, 3UL );
      checkNonZeros( mat2, 1UL, 3UL );
      checkNonZeros( mat2, 2UL, 3UL );
      checkNonZeros( mat2, 3UL, 3UL );

      if( mat2(0,0) != 6 || mat2(0,1) != 6 || mat2(0,2) != 6 || mat2(0,3) != 6 ||
          mat2(1,0) != 6 || mat2(1,1) != 6 || mat2(1,2) != 6 || mat2(1,3) != 6 ||
          mat2(2,0) != 6 || mat2(2,1) != 6 || mat2(2,2) != 6 || mat2(2,3) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 6 6 6 6 )\n( 6 6 6 6 )\n( 6 6 6 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major UniformMatrix dense matrix multiplication assignment (non-uniform)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat1{ { 2, 2, 2 }, { 2, 0, 2 }, { 2, 2, 2 } };

      try {
         blaze::UniformMatrix<int,blaze::columnMajor> mat2( 2UL, 3UL, 1 );
         mat2 *= mat1;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-uniform dense matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Column-major sparse matrix multiplication assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major UniformMatrix sparse matrix multiplication assignment";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1{ { 2, 2, 2, 2 }
                                                       , { 2, 2, 2, 2 }
                                                       , { 2, 2, 2, 2 } };

      blaze::UniformMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL, 1 );

      mat2 *= mat1;

      checkRows    ( mat2,  3UL );
      checkColumns ( mat2,  4UL );
      checkNonZeros( mat2, 12UL );
      checkNonZeros( mat2, 0UL, 3UL );
      checkNonZeros( mat2, 1UL, 3UL );
      checkNonZeros( mat2, 2UL, 3UL );
      checkNonZeros( mat2, 3UL, 3UL );

      if( mat2(0,0) != 6 || mat2(0,1) != 6 || mat2(0,2) != 6 || mat2(0,3) != 6 ||
          mat2(1,0) != 6 || mat2(1,1) != 6 || mat2(1,2) != 6 || mat2(1,3) != 6 ||
          mat2(2,0) != 6 || mat2(2,1) != 6 || mat2(2,2) != 6 || mat2(2,3) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 6 6 6 6 )\n( 6 6 6 6 )\n( 6 6 6 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major UniformMatrix sparse matrix multiplication assignment (non-uniform)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1{ { 2, 2, 2, 2 }
                                                       , { 2, 0, 2, 2 }
                                                       , { 2, 2, 2, 2 } };

      try {
         blaze::UniformMatrix<int,blaze::columnMajor> mat2( 2UL, 3UL, 1 );
         mat2 *= mat1;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-uniform sparse matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   {
      test_ = "Column-major/column-major UniformMatrix sparse matrix multiplication assignment";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1{ { 2, 2, 2, 2 }
                                                          , { 2, 2, 2, 2 }
                                                          , { 2, 2, 2, 2 } };

      blaze::UniformMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL, 1 );

      mat2 *= mat1;

      checkRows    ( mat2,  3UL );
      checkColumns ( mat2,  4UL );
      checkNonZeros( mat2, 12UL );
      checkNonZeros( mat2, 0UL, 3UL );
      checkNonZeros( mat2, 1UL, 3UL );
      checkNonZeros( mat2, 2UL, 3UL );
      checkNonZeros( mat2, 3UL, 3UL );

      if( mat2(0,0) != 6 || mat2(0,1) != 6 || mat2(0,2) != 6 || mat2(0,3) != 6 ||
          mat2(1,0) != 6 || mat2(1,1) != 6 || mat2(1,2) != 6 || mat2(1,3) != 6 ||
          mat2(2,0) != 6 || mat2(2,1) != 6 || mat2(2,2) != 6 || mat2(2,3) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 6 6 6 6 )\n( 6 6 6 6 )\n( 6 6 6 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major UniformMatrix sparse matrix multiplication assignment (non-uniform)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1{ { 2, 2, 2, 2 }
                                                          , { 2, 0, 2, 2 }
                                                          , { 2, 2, 2, 2 } };

      try {
         blaze::UniformMatrix<int,blaze::columnMajor> mat2( 2UL, 3UL, 1 );
         mat2 *= mat1;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-uniform sparse matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of all UniformMatrix (self-)scaling operations.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all available ways to scale an instance of the UniformMatrix
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testScaling()
{
   //=====================================================================================
   // Row-major self-scaling (M*=s)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (M*=s)";

      blaze::UniformMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 2 );

      mat *= 2;

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 3UL );
      checkNonZeros( mat, 9UL );
      checkNonZeros( mat, 0UL, 3UL );
      checkNonZeros( mat, 1UL, 3UL );
      checkNonZeros( mat, 2UL, 3UL );

      if( mat(0,0) != 4 || mat(0,1) != 4 || mat(0,2) != 4 ||
          mat(1,0) != 4 || mat(1,1) != 4 || mat(1,2) != 4 ||
          mat(2,0) != 4 || mat(2,1) != 4 || mat(2,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 4 4 4 )\n( 4 4 4 )\n( 4 4 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major self-scaling (M=M*s)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (M=M*s)";

      blaze::UniformMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 2 );

      mat = mat * 2;

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 3UL );
      checkNonZeros( mat, 9UL );
      checkNonZeros( mat, 0UL, 3UL );
      checkNonZeros( mat, 1UL, 3UL );
      checkNonZeros( mat, 2UL, 3UL );

      if( mat(0,0) != 4 || mat(0,1) != 4 || mat(0,2) != 4 ||
          mat(1,0) != 4 || mat(1,1) != 4 || mat(1,2) != 4 ||
          mat(2,0) != 4 || mat(2,1) != 4 || mat(2,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 4 4 4 )\n( 4 4 4 )\n( 4 4 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major self-scaling (M=s*M)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (M=s*M)";

      blaze::UniformMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 2 );

      mat = 2 * mat;

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 3UL );
      checkNonZeros( mat, 9UL );
      checkNonZeros( mat, 0UL, 3UL );
      checkNonZeros( mat, 1UL, 3UL );
      checkNonZeros( mat, 2UL, 3UL );

      if( mat(0,0) != 4 || mat(0,1) != 4 || mat(0,2) != 4 ||
          mat(1,0) != 4 || mat(1,1) != 4 || mat(1,2) != 4 ||
          mat(2,0) != 4 || mat(2,1) != 4 || mat(2,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 4 4 4 )\n( 4 4 4 )\n( 4 4 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major self-scaling (M/=s)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (M/=s)";

      blaze::UniformMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 4 );

      mat /= 2;

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 3UL );
      checkNonZeros( mat, 9UL );
      checkNonZeros( mat, 0UL, 3UL );
      checkNonZeros( mat, 1UL, 3UL );
      checkNonZeros( mat, 2UL, 3UL );

      if( mat(0,0) != 2 || mat(0,1) != 2 || mat(0,2) != 2 ||
          mat(1,0) != 2 || mat(1,1) != 2 || mat(1,2) != 2 ||
          mat(2,0) != 2 || mat(2,1) != 2 || mat(2,2) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major self-scaling (M=M/s)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (M=M/s)";

      blaze::UniformMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 4 );

      mat = mat / 2;

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 3UL );
      checkNonZeros( mat, 9UL );
      checkNonZeros( mat, 0UL, 3UL );
      checkNonZeros( mat, 1UL, 3UL );
      checkNonZeros( mat, 2UL, 3UL );

      if( mat(0,0) != 2 || mat(0,1) != 2 || mat(0,2) != 2 ||
          mat(1,0) != 2 || mat(1,1) != 2 || mat(1,2) != 2 ||
          mat(2,0) != 2 || mat(2,1) != 2 || mat(2,2) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major UniformMatrix::scale()
   //=====================================================================================

   {
      test_ = "Row-major UniformMatrix::scale() (int)";

      // Initialization check
      blaze::UniformMatrix<int,blaze::rowMajor> mat( 3UL, 2UL, 2 );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 6UL );
      checkNonZeros( mat, 0UL, 2UL );
      checkNonZeros( mat, 1UL, 2UL );
      checkNonZeros( mat, 2UL, 2UL );

      if( mat(0,0) != 2 || mat(0,1) != 2 ||
          mat(1,0) != 2 || mat(1,1) != 2 ||
          mat(2,0) != 2 || mat(2,1) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 2 2 )\n( 2 2 )\n( 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Integral scaling of the matrix
      mat.scale( 2 );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 6UL );
      checkNonZeros( mat, 0UL, 2UL );
      checkNonZeros( mat, 1UL, 2UL );
      checkNonZeros( mat, 2UL, 2UL );

      if( mat(0,0) != 4 || mat(0,1) != 4 ||
          mat(1,0) != 4 || mat(1,1) != 4 ||
          mat(2,0) != 4 || mat(2,1) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 4 4 )\n( 4 4 )\n( 4 4 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Floating point scaling of the matrix
      mat.scale( 0.5 );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 6UL );
      checkNonZeros( mat, 0UL, 2UL );
      checkNonZeros( mat, 1UL, 2UL );
      checkNonZeros( mat, 2UL, 2UL );

      if( mat(0,0) != 2 || mat(0,1) != 2 ||
          mat(1,0) != 2 || mat(1,1) != 2 ||
          mat(2,0) != 2 || mat(2,1) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 2 2 )\n( 2 2 )\n( 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major UniformMatrix::scale() (complex)";

      using blaze::complex;

      blaze::UniformMatrix<complex<float>,blaze::rowMajor> mat( 2UL, 2UL, complex<float>( 2.0F, 0.0F ) );
      mat.scale( complex<float>( 3.0F, 0.0F ) );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 4UL );
      checkNonZeros( mat, 4UL );
      checkNonZeros( mat, 0UL, 2UL );
      checkNonZeros( mat, 1UL, 2UL );

      if( mat(0,0) != complex<float>( 6.0F, 0.0F ) || mat(0,1) != complex<float>( 6.0F, 0.0F ) ||
          mat(1,0) != complex<float>( 6.0F, 0.0F ) || mat(1,1) != complex<float>( 6.0F, 0.0F ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( (6,0) (6,0)\n(6,0) (6,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (M*=s)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (M*=s)";

      blaze::UniformMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 2 );

      mat *= 2;

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 3UL );
      checkNonZeros( mat, 9UL );
      checkNonZeros( mat, 0UL, 3UL );
      checkNonZeros( mat, 1UL, 3UL );
      checkNonZeros( mat, 2UL, 3UL );

      if( mat(0,0) != 4 || mat(0,1) != 4 || mat(0,2) != 4 ||
          mat(1,0) != 4 || mat(1,1) != 4 || mat(1,2) != 4 ||
          mat(2,0) != 4 || mat(2,1) != 4 || mat(2,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 4 4 4 )\n( 4 4 4 )\n( 4 4 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (M=M*s)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (M=M*s)";

      blaze::UniformMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 2 );

      mat = mat * 2;

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 3UL );
      checkNonZeros( mat, 9UL );
      checkNonZeros( mat, 0UL, 3UL );
      checkNonZeros( mat, 1UL, 3UL );
      checkNonZeros( mat, 2UL, 3UL );

      if( mat(0,0) != 4 || mat(0,1) != 4 || mat(0,2) != 4 ||
          mat(1,0) != 4 || mat(1,1) != 4 || mat(1,2) != 4 ||
          mat(2,0) != 4 || mat(2,1) != 4 || mat(2,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 4 4 4 )\n( 4 4 4 )\n( 4 4 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (M=s*M)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (M=s*M)";

      blaze::UniformMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 2 );

      mat = 2 * mat;

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 3UL );
      checkNonZeros( mat, 9UL );
      checkNonZeros( mat, 0UL, 3UL );
      checkNonZeros( mat, 1UL, 3UL );
      checkNonZeros( mat, 2UL, 3UL );

      if( mat(0,0) != 4 || mat(0,1) != 4 || mat(0,2) != 4 ||
          mat(1,0) != 4 || mat(1,1) != 4 || mat(1,2) != 4 ||
          mat(2,0) != 4 || mat(2,1) != 4 || mat(2,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 4 4 4 )\n( 4 4 4 )\n( 4 4 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (M/=s)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (M/=s)";

      blaze::UniformMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 4 );

      mat /= 2;

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 3UL );
      checkNonZeros( mat, 9UL );
      checkNonZeros( mat, 0UL, 3UL );
      checkNonZeros( mat, 1UL, 3UL );
      checkNonZeros( mat, 2UL, 3UL );

      if( mat(0,0) != 2 || mat(0,1) != 2 || mat(0,2) != 2 ||
          mat(1,0) != 2 || mat(1,1) != 2 || mat(1,2) != 2 ||
          mat(2,0) != 2 || mat(2,1) != 2 || mat(2,2) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (M=M/s)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (M=M/s)";

      blaze::UniformMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 4 );

      mat = mat / 2;

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 3UL );
      checkNonZeros( mat, 9UL );
      checkNonZeros( mat, 0UL, 3UL );
      checkNonZeros( mat, 1UL, 3UL );
      checkNonZeros( mat, 2UL, 3UL );

      if( mat(0,0) != 2 || mat(0,1) != 2 || mat(0,2) != 2 ||
          mat(1,0) != 2 || mat(1,1) != 2 || mat(1,2) != 2 ||
          mat(2,0) != 2 || mat(2,1) != 2 || mat(2,2) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major UniformMatrix::scale()
   //=====================================================================================

   {
      test_ = "Column-major UniformMatrix::scale() (int)";

      // Initialization check
      blaze::UniformMatrix<int,blaze::columnMajor> mat( 3UL, 2UL, 2 );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 6UL );
      checkNonZeros( mat, 0UL, 3UL );
      checkNonZeros( mat, 1UL, 3UL );

      if( mat(0,0) != 2 || mat(0,1) != 2 ||
          mat(1,0) != 2 || mat(1,1) != 2 ||
          mat(2,0) != 2 || mat(2,1) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 2 2 )\n( 2 2 )\n( 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Integral scaling of the matrix
      mat.scale( 2 );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 6UL );
      checkNonZeros( mat, 0UL, 3UL );
      checkNonZeros( mat, 1UL, 3UL );

      if( mat(0,0) != 4 || mat(0,1) != 4 ||
          mat(1,0) != 4 || mat(1,1) != 4 ||
          mat(2,0) != 4 || mat(2,1) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 4 4 )\n( 4 4 )\n( 4 4 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Floating point scaling of the matrix
      mat.scale( 0.5 );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 6UL );
      checkNonZeros( mat, 0UL, 3UL );
      checkNonZeros( mat, 1UL, 3UL );

      if( mat(0,0) != 2 || mat(0,1) != 2 ||
          mat(1,0) != 2 || mat(1,1) != 2 ||
          mat(2,0) != 2 || mat(2,1) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 2 2 )\n( 2 2 )\n( 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major UniformMatrix::scale() (complex)";

      using blaze::complex;

      blaze::UniformMatrix<complex<float>,blaze::columnMajor> mat( 2UL, 2UL, complex<float>( 2.0F, 0.0F ) );
      mat.scale( complex<float>( 3.0F, 0.0F ) );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 4UL );
      checkNonZeros( mat, 4UL );
      checkNonZeros( mat, 0UL, 2UL );
      checkNonZeros( mat, 1UL, 2UL );

      if( mat(0,0) != complex<float>( 6.0F, 0.0F ) || mat(0,1) != complex<float>( 6.0F, 0.0F ) ||
          mat(1,0) != complex<float>( 6.0F, 0.0F ) || mat(1,1) != complex<float>( 6.0F, 0.0F ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( (6,0) (6,0)\n(6,0) (6,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the UniformMatrix function call operator.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of adding and accessing elements via the function call operator
// of the UniformMatrix class template. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
void ClassTest::testFunctionCall()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major UniformMatrix::operator()";

      blaze::UniformMatrix<int,blaze::rowMajor> mat( 3UL, 5UL, 2 );

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  5UL );
      checkCapacity( mat, 15UL );
      checkNonZeros( mat, 15UL );
      checkNonZeros( mat, 0UL, 5UL );
      checkNonZeros( mat, 1UL, 5UL );
      checkNonZeros( mat, 2UL, 5UL );

      // Accessing all elements
      if( mat(0,0) != 2 || mat(0,1) != 2 || mat(0,2) != 2 || mat(0,3) != 2 || mat(0,4) != 2 ||
          mat(1,0) != 2 || mat(1,1) != 2 || mat(1,2) != 2 || mat(1,3) != 2 || mat(1,4) != 2 ||
          mat(2,0) != 2 || mat(2,1) != 2 || mat(2,2) != 2 || mat(2,3) != 2 || mat(2,4) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 2 2 2 2 2 )\n( 2 2 2 2 2 )\n( 2 2 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major UniformMatrix::operator()";

      blaze::UniformMatrix<int,blaze::columnMajor> mat( 3UL, 5UL, 2 );

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  5UL );
      checkCapacity( mat, 15UL );
      checkNonZeros( mat, 15UL );
      checkNonZeros( mat, 0UL, 3UL );
      checkNonZeros( mat, 1UL, 3UL );
      checkNonZeros( mat, 2UL, 3UL );
      checkNonZeros( mat, 3UL, 3UL );
      checkNonZeros( mat, 4UL, 3UL );

      // Accessing all elements
      if( mat(0,0) != 2 || mat(0,1) != 2 || mat(0,2) != 2 || mat(0,3) != 2 || mat(0,4) != 2 ||
          mat(1,0) != 2 || mat(1,1) != 2 || mat(1,2) != 2 || mat(1,3) != 2 || mat(1,4) != 2 ||
          mat(2,0) != 2 || mat(2,1) != 2 || mat(2,2) != 2 || mat(2,3) != 2 || mat(2,4) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 2 2 2 2 2 )\n( 2 2 2 2 2 )\n( 2 2 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c at() member function of the UniformMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of adding and accessing elements via the \c at() member function
// of the UniformMatrix class template. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
void ClassTest::testAt()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major UniformMatrix::at()";

      blaze::UniformMatrix<int,blaze::rowMajor> mat( 3UL, 5UL, 2 );

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  5UL );
      checkCapacity( mat, 15UL );
      checkNonZeros( mat, 15UL );
      checkNonZeros( mat, 0UL, 5UL );
      checkNonZeros( mat, 1UL, 5UL );
      checkNonZeros( mat, 2UL, 5UL );

      // Accessing all elements
      if( mat.at(0,0) != 2 || mat.at(0,1) != 2 || mat.at(0,2) != 2 || mat.at(0,3) != 2 || mat.at(0,4) != 2 ||
          mat.at(1,0) != 2 || mat.at(1,1) != 2 || mat.at(1,2) != 2 || mat.at(1,3) != 2 || mat.at(1,4) != 2 ||
          mat.at(2,0) != 2 || mat.at(2,1) != 2 || mat.at(2,2) != 2 || mat.at(2,3) != 2 || mat.at(2,4) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 2 2 2 2 2 )\n( 2 2 2 2 2 )\n( 2 2 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Attempt to assign to the element (3,0)
      try {
         mat.at(3,0);

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Out-of-bound access succeeded\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 2 2 2 2 2 )\n( 2 2 2 2 2 )\n( 2 2 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::out_of_range& ) {}

      // Attempt to assign to the element (0,5)
      try {
         mat.at(0,5);

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Out-of-bound access succeeded\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 2 2 2 2 2 )\n( 2 2 2 2 2 )\n( 2 2 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::out_of_range& ) {}
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major UniformMatrix::at()";

      blaze::UniformMatrix<int,blaze::columnMajor> mat( 3UL, 5UL, 2 );

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  5UL );
      checkCapacity( mat, 15UL );
      checkNonZeros( mat, 15UL );
      checkNonZeros( mat, 0UL, 3UL );
      checkNonZeros( mat, 1UL, 3UL );
      checkNonZeros( mat, 2UL, 3UL );
      checkNonZeros( mat, 3UL, 3UL );
      checkNonZeros( mat, 4UL, 3UL );

      // Accessing all elements
      if( mat.at(0,0) != 2 || mat.at(0,1) != 2 || mat.at(0,2) != 2 || mat.at(0,3) != 2 || mat.at(0,4) != 2 ||
          mat.at(1,0) != 2 || mat.at(1,1) != 2 || mat.at(1,2) != 2 || mat.at(1,3) != 2 || mat.at(1,4) != 2 ||
          mat.at(2,0) != 2 || mat.at(2,1) != 2 || mat.at(2,2) != 2 || mat.at(2,3) != 2 || mat.at(2,4) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Function call operator failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 2 2 2 2 2 )\n( 2 2 2 2 2 )\n( 2 2 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Attempt to assign to the element (3,0)
      try {
         mat.at(3,0);

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Out-of-bound access succeeded\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 2 2 2 2 2 )\n( 2 2 2 2 2 )\n( 2 2 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::out_of_range& ) {}

      // Attempt to assign to the element (0,5)
      try {
         mat.at(0,5);

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Out-of-bound access succeeded\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 2 2 2 2 2 )\n( 2 2 2 2 2 )\n( 2 2 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::out_of_range& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the UniformMatrix iterator implementation.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the iterator implementation of the UniformMatrix class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testIterator()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      using MatrixType    = blaze::UniformMatrix<int,blaze::rowMajor>;
      using ConstIterator = MatrixType::ConstIterator;

      MatrixType mat( 3UL, 3UL, 2 );

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

      // Counting the number of elements in 1st row via ConstIterator (end-begin)
      {
         test_ = "Row-major ConstIterator subtraction (end-begin)";

         const ptrdiff_t number( cend( mat, 1UL ) - cbegin( mat, 1UL ) );

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

         const ptrdiff_t number( cbegin( mat, 1UL ) - cend( mat, 1UL ) );

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

         ConstIterator it ( cbegin( mat, 2UL ) );
         ConstIterator end( cend( mat, 2UL ) );

         if( it == end || *it != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid initial iterator detected\n";
            throw std::runtime_error( oss.str() );
         }

         ++it;

         if( it == end || *it != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator pre-increment failed\n";
            throw std::runtime_error( oss.str() );
         }

         --it;

         if( it == end || *it != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator pre-decrement failed\n";
            throw std::runtime_error( oss.str() );
         }

         it++;

         if( it == end || *it != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator post-increment failed\n";
            throw std::runtime_error( oss.str() );
         }

         it--;

         if( it == end || *it != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator post-decrement failed\n";
            throw std::runtime_error( oss.str() );
         }

         it += 2UL;

         if( it == end || *it != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator addition assignment failed\n";
            throw std::runtime_error( oss.str() );
         }

         it -= 2UL;

         if( it == end || *it != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator subtraction assignment failed\n";
            throw std::runtime_error( oss.str() );
         }

         it = it + 2UL;

         if( it == end || *it != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator/scalar addition failed\n";
            throw std::runtime_error( oss.str() );
         }

         it = it - 2UL;

         if( it == end || *it != 2 ) {
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
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      using MatrixType    = blaze::UniformMatrix<int,blaze::columnMajor>;
      using ConstIterator = MatrixType::ConstIterator;

      MatrixType mat( 3UL, 3UL, 2 );

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

      // Counting the number of elements in 1st row via ConstIterator (end-begin)
      {
         test_ = "Column-major ConstIterator subtraction (end-begin)";

         const ptrdiff_t number( cend( mat, 1UL ) - cbegin( mat, 1UL ) );

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
         test_ = "Column-major ConstIterator subtraction (begin-end)";

         const ptrdiff_t number( cbegin( mat, 1UL ) - cend( mat, 1UL ) );

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

         ConstIterator it ( cbegin( mat, 2UL ) );
         ConstIterator end( cend( mat, 2UL ) );

         if( it == end || *it != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid initial iterator detected\n";
            throw std::runtime_error( oss.str() );
         }

         ++it;

         if( it == end || *it != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator pre-increment failed\n";
            throw std::runtime_error( oss.str() );
         }

         --it;

         if( it == end || *it != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator pre-decrement failed\n";
            throw std::runtime_error( oss.str() );
         }

         it++;

         if( it == end || *it != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator post-increment failed\n";
            throw std::runtime_error( oss.str() );
         }

         it--;

         if( it == end || *it != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator post-decrement failed\n";
            throw std::runtime_error( oss.str() );
         }

         it += 2UL;

         if( it == end || *it != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator addition assignment failed\n";
            throw std::runtime_error( oss.str() );
         }

         it -= 2UL;

         if( it == end || *it != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator subtraction assignment failed\n";
            throw std::runtime_error( oss.str() );
         }

         it = it + 2UL;

         if( it == end || *it != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator/scalar addition failed\n";
            throw std::runtime_error( oss.str() );
         }

         it = it - 2UL;

         if( it == end || *it != 2 ) {
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
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c nonZeros() member function of the UniformMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c nonZeros() member function of the UniformMatrix class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testNonZeros()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major UniformMatrix::nonZeros()";

      {
         blaze::UniformMatrix<int,blaze::rowMajor> mat( 2UL, 3UL, 0 );

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 6UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );

         if( mat(0,0) != 0 || mat(0,1) != 0 || mat(0,2) != 0 ||
             mat(1,0) != 0 || mat(1,1) != 0 || mat(1,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 0 0 )\n( 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::UniformMatrix<int,blaze::rowMajor> mat( 2UL, 3UL, 2 );

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 6UL );
         checkNonZeros( mat, 6UL );
         checkNonZeros( mat, 0UL, 3UL );
         checkNonZeros( mat, 1UL, 3UL );

         if( mat(0,0) != 2 || mat(0,1) != 2 || mat(0,2) != 2 ||
             mat(1,0) != 2 || mat(1,1) != 2 || mat(1,2) != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major UniformMatrix::nonZeros()";

      {
         blaze::UniformMatrix<int,blaze::columnMajor> mat( 2UL, 3UL, 0 );

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 6UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 0UL );

         if( mat(0,0) != 0 || mat(0,1) != 0 || mat(0,2) != 0 ||
             mat(1,0) != 0 || mat(1,1) != 0 || mat(1,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 0 0 )\n( 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::UniformMatrix<int,blaze::columnMajor> mat( 2UL, 3UL, 2 );

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 6UL );
         checkNonZeros( mat, 6UL );
         checkNonZeros( mat, 0UL, 2UL );
         checkNonZeros( mat, 1UL, 2UL );
         checkNonZeros( mat, 2UL, 2UL );

         if( mat(0,0) != 2 || mat(0,1) != 2 || mat(0,2) != 2 ||
             mat(1,0) != 2 || mat(1,1) != 2 || mat(1,2) != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c reset() member function of the UniformMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c reset() member function of the UniformMatrix class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testReset()
{
   using blaze::reset;


   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major UniformMatrix::reset()";

      // Resetting a default initialized matrix
      {
         blaze::UniformMatrix<int,blaze::rowMajor> mat;

         reset( mat );

         checkRows    ( mat, 0UL );
         checkColumns ( mat, 0UL );
         checkNonZeros( mat, 0UL );
      }

      // Resetting an initialized matrix
      {
         blaze::UniformMatrix<int,blaze::rowMajor> mat( 2UL, 3UL, 2 );

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 6UL );
         checkNonZeros( mat, 6UL );
         checkNonZeros( mat, 0UL, 3UL );
         checkNonZeros( mat, 1UL, 3UL );

         if( mat(0,0) != 2 || mat(0,1) != 2 || mat(0,2) != 2 ||
             mat(1,0) != 2 || mat(1,1) != 2 || mat(1,2) != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n";
            throw std::runtime_error( oss.str() );
         }

         reset( mat );

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 6UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );

         if( mat(0,0) != 0 || mat(0,1) != 0 || mat(0,2) != 0 ||
             mat(1,0) != 0 || mat(1,1) != 0 || mat(1,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Reset operation failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 0 0 )\n( 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major UniformMatrix::reset()";

      // Resetting a default initialized matrix
      {
         blaze::UniformMatrix<int,blaze::columnMajor> mat;

         reset( mat );

         checkRows    ( mat, 0UL );
         checkColumns ( mat, 0UL );
         checkNonZeros( mat, 0UL );
      }

      // Resetting an initialized matrix
      {
         blaze::UniformMatrix<int,blaze::columnMajor> mat( 2UL, 3UL, 2 );

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 6UL );
         checkNonZeros( mat, 6UL );
         checkNonZeros( mat, 0UL, 2UL );
         checkNonZeros( mat, 1UL, 2UL );
         checkNonZeros( mat, 2UL, 2UL );

         if( mat(0,0) != 2 || mat(0,1) != 2 || mat(0,2) != 2 ||
             mat(1,0) != 2 || mat(1,1) != 2 || mat(1,2) != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n";
            throw std::runtime_error( oss.str() );
         }

         reset( mat );

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 6UL );
         checkNonZeros( mat, 0UL );
         checkNonZeros( mat, 0UL, 0UL );
         checkNonZeros( mat, 1UL, 0UL );
         checkNonZeros( mat, 2UL, 0UL );

         if( mat(0,0) != 0 || mat(0,1) != 0 || mat(0,2) != 0 ||
             mat(1,0) != 0 || mat(1,1) != 0 || mat(1,2) != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Reset operation failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 0 0 0 )\n( 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c clear() member function of the UniformMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c clear() member function of the UniformMatrix class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testClear()
{
   using blaze::clear;


   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major UniformMatrix::clear()";

      // Clearing a default constructed matrix
      {
         blaze::UniformMatrix<int,blaze::rowMajor> mat;

         clear( mat );

         checkRows    ( mat, 0UL );
         checkColumns ( mat, 0UL );
         checkNonZeros( mat, 0UL );
      }

      // Clearing an initialized matrix
      {
         blaze::UniformMatrix<int,blaze::rowMajor> mat( 2UL, 3UL, 2 );

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 6UL );
         checkNonZeros( mat, 6UL );
         checkNonZeros( mat, 0UL, 3UL );
         checkNonZeros( mat, 1UL, 3UL );

         if( mat(0,0) != 2 || mat(0,1) != 2 || mat(0,2) != 2 ||
             mat(1,0) != 2 || mat(1,1) != 2 || mat(1,2) != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n";
            throw std::runtime_error( oss.str() );
         }

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
      test_ = "Column-major UniformMatrix::clear()";

      // Clearing a default constructed matrix
      {
         blaze::UniformMatrix<int,blaze::columnMajor> mat;

         clear( mat );

         checkRows    ( mat, 0UL );
         checkColumns ( mat, 0UL );
         checkNonZeros( mat, 0UL );
      }

      // Clearing an initialized matrix
      {
         blaze::UniformMatrix<int,blaze::columnMajor> mat( 2UL, 3UL, 2 );

         checkRows    ( mat, 2UL );
         checkColumns ( mat, 3UL );
         checkCapacity( mat, 6UL );
         checkNonZeros( mat, 6UL );
         checkNonZeros( mat, 0UL, 2UL );
         checkNonZeros( mat, 1UL, 2UL );
         checkNonZeros( mat, 2UL, 2UL );

         if( mat(0,0) != 2 || mat(0,1) != 2 || mat(0,2) != 2 ||
             mat(1,0) != 2 || mat(1,1) != 2 || mat(1,2) != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Initialization failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n";
            throw std::runtime_error( oss.str() );
         }

         clear( mat );

         checkRows    ( mat, 0UL );
         checkColumns ( mat, 0UL );
         checkNonZeros( mat, 0UL );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c resize() member function of the UniformMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c resize() member function of the UniformMatrix class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testResize()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major UniformMatrix::resize()";

      // Initialization check
      blaze::UniformMatrix<int,blaze::rowMajor> mat;

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

      // Resizing to 2x1
      mat.resize( 2UL, 1UL );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 1UL );
      checkCapacity( mat, 2UL );

      if( mat(0,0) != 0 || mat(1,0) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 )\n( 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 3x2 and preserving the elements
      mat = 5;
      mat.resize( 3UL, 2UL, true );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 6UL );

      if( mat(0,0) != 5 || mat(0,1) != 5 ||
          mat(1,0) != 5 || mat(1,1) != 5 ||
          mat(2,0) != 5 || mat(2,1) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 5 5 )\n( 5 5 )\n( 5 5 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 2x2 and preserving the elements
      mat.resize( 2UL, 2UL, true );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 4UL );
      checkNonZeros( mat, 4UL );
      checkNonZeros( mat, 0UL, 2UL );
      checkNonZeros( mat, 1UL, 2UL );

      if( mat(0,0) != 5 || mat(0,1) != 5 ||
          mat(1,0) != 5 || mat(1,1) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 5 5 )\n( 5 5 )\n";
         throw std::runtime_error( oss.str() );
      }

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
      test_ = "Column-major UniformMatrix::resize()";

      // Initialization check
      blaze::UniformMatrix<int,blaze::columnMajor> mat;

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

      // Resizing to 2x1
      mat.resize( 2UL, 1UL );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 1UL );
      checkCapacity( mat, 2UL );

      if( mat(0,0) != 0 || mat(1,0) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 )\n( 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 3x2 and preserving the elements
      mat = 5;
      mat.resize( 3UL, 2UL, true );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 6UL );

      if( mat(0,0) != 5 || mat(0,1) != 5 ||
          mat(1,0) != 5 || mat(1,1) != 5 ||
          mat(2,0) != 5 || mat(2,1) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 5 5 )\n( 5 5 )\n( 5 5 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 2x2 and preserving the elements
      mat.resize( 2UL, 2UL, true );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 4UL );
      checkNonZeros( mat, 4UL );
      checkNonZeros( mat, 0UL, 2UL );
      checkNonZeros( mat, 1UL, 2UL );

      if( mat(0,0) != 5 || mat(0,1) != 5 ||
          mat(1,0) != 5 || mat(1,1) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 5 5 )\n( 5 5 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 0x0
      mat.resize( 0UL, 0UL );

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c extend() member function of the UniformMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c extend() member function of the UniformMatrix class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testExtend()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major UniformMatrix::extend()";

      // Initialization check
      blaze::UniformMatrix<int,blaze::rowMajor> mat;

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );

      // Increasing the size of the matrix
      mat.extend( 2UL, 2UL );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 3UL );

      if( mat(0,0) != 0 || mat(0,1) != 0 ||
          mat(1,0) != 0 || mat(1,1) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 )\n( 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Further increasing the size of the matrix and preserving the elements
      mat.extend( 1UL, 1UL, true );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 9UL );

      if( mat(0,0) != 0 || mat(0,1) != 0 || mat(0,2) != 0 ||
          mat(1,0) != 0 || mat(1,1) != 0 || mat(1,2) != 0 ||
          mat(2,0) != 0 || mat(2,1) != 0 || mat(2,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 0 )\n( 0 0 0 )\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Further increasing the size of the matrix
      mat.extend( 4UL, 10UL, false );

      checkRows    ( mat,  7UL );
      checkColumns ( mat, 13UL );
      checkCapacity( mat, 91UL );
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major UniformMatrix::extend()";

      // Initialization check
      blaze::UniformMatrix<int,blaze::columnMajor> mat;

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );

      // Increasing the size of the matrix
      mat.extend( 2UL, 2UL );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 3UL );

      if( mat(0,0) != 0 || mat(0,1) != 0 ||
          mat(1,0) != 0 || mat(1,1) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 )\n( 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Further increasing the size of the matrix and preserving the elements
      mat.extend( 1UL, 1UL, true );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 9UL );

      if( mat(0,0) != 0 || mat(0,1) != 0 || mat(0,2) != 0 ||
          mat(1,0) != 0 || mat(1,1) != 0 || mat(1,2) != 0 ||
          mat(2,0) != 0 || mat(2,1) != 0 || mat(2,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resizing the matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 0 )\n( 0 0 0 )\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Further increasing the size of the matrix
      mat.extend( 4UL, 10UL, false );

      checkRows    ( mat,  7UL );
      checkColumns ( mat, 13UL );
      checkCapacity( mat, 91UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c swap() functionality of the UniformMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c swap() function of the UniformMatrix class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testSwap()
{
   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major UniformMatrix swap";

      blaze::UniformMatrix<int,blaze::rowMajor> mat1( 3UL, 2UL, 2 );
      blaze::UniformMatrix<int,blaze::rowMajor> mat2( 2UL, 3UL, 5 );

      swap( mat1, mat2 );

      checkRows    ( mat1, 2UL );
      checkColumns ( mat1, 3UL );
      checkCapacity( mat1, 6UL );
      checkNonZeros( mat1, 6UL );
      checkNonZeros( mat1, 0UL, 3UL );
      checkNonZeros( mat1, 1UL, 3UL );

      if( mat1(0,0) != 5 || mat1(0,1) != 5 || mat1(0,2) != 5 ||
          mat1(1,0) != 5 || mat1(1,1) != 5 || mat1(1,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the first matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n( 5 5 5 )\n( 5 5 5 )\n";
         throw std::runtime_error( oss.str() );
      }

      checkRows    ( mat2, 3UL );
      checkColumns ( mat2, 2UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 2 || mat2(0,1) != 2 ||
          mat2(1,0) != 2 || mat2(1,1) != 2 ||
          mat2(2,0) != 2 || mat2(2,1) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the second matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 2 2 )\n( 2 2 )\n( 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major UniformMatrix swap";

      blaze::UniformMatrix<int,blaze::columnMajor> mat1( 3UL, 2UL, 2 );
      blaze::UniformMatrix<int,blaze::columnMajor> mat2( 2UL, 3UL, 5 );

      swap( mat1, mat2 );

      checkRows    ( mat1, 2UL );
      checkColumns ( mat1, 3UL );
      checkCapacity( mat1, 6UL );
      checkNonZeros( mat1, 6UL );
      checkNonZeros( mat1, 0UL, 2UL );
      checkNonZeros( mat1, 1UL, 2UL );
      checkNonZeros( mat1, 2UL, 2UL );

      if( mat1(0,0) != 5 || mat1(0,1) != 5 || mat1(0,2) != 5 ||
          mat1(1,0) != 5 || mat1(1,1) != 5 || mat1(1,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the first matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n( 5 5 5 )\n( 5 5 5 )\n";
         throw std::runtime_error( oss.str() );
      }

      checkRows    ( mat2, 3UL );
      checkColumns ( mat2, 2UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 3UL );
      checkNonZeros( mat2, 1UL, 3UL );

      if( mat2(0,0) != 2 || mat2(0,1) != 2 ||
          mat2(1,0) != 2 || mat2(1,1) != 2 ||
          mat2(2,0) != 2 || mat2(2,1) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the second matrix failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 2 2 )\n( 2 2 )\n( 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c transpose() member function of the UniformMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c transpose() member function of the UniformMatrix
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
         blaze::UniformMatrix<int,blaze::rowMajor> mat( 3UL, 5UL, 2 );

         transpose( mat );

         checkRows    ( mat,  5UL );
         checkColumns ( mat,  3UL );
         checkCapacity( mat, 15UL );
         checkNonZeros( mat, 15UL );
         checkNonZeros( mat,  0UL, 3UL );
         checkNonZeros( mat,  1UL, 3UL );
         checkNonZeros( mat,  2UL, 3UL );
         checkNonZeros( mat,  3UL, 3UL );
         checkNonZeros( mat,  4UL, 3UL );

         if( mat(0,0) != 2 || mat(0,1) != 2 || mat(0,2) != 2 ||
             mat(1,0) != 2 || mat(1,1) != 2 || mat(1,2) != 2 ||
             mat(2,0) != 2 || mat(2,1) != 2 || mat(2,2) != 2 ||
             mat(3,0) != 2 || mat(3,1) != 2 || mat(3,2) != 2 ||
             mat(4,0) != 2 || mat(4,1) != 2 || mat(4,2) != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Transpose operation failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n( 2 2 2 )\n( 2 2 2 )\n( 2 2 2 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Self-transpose of a 5x3 matrix
      {
         blaze::UniformMatrix<int,blaze::rowMajor> mat( 5UL, 3UL, 2 );

         transpose( mat );

         checkRows    ( mat,  3UL );
         checkColumns ( mat,  5UL );
         checkCapacity( mat, 15UL );
         checkNonZeros( mat, 15UL );
         checkNonZeros( mat,  0UL, 5UL );
         checkNonZeros( mat,  1UL, 5UL );
         checkNonZeros( mat,  2UL, 5UL );

         if( mat(0,0) != 2 || mat(0,1) != 2 || mat(0,2) != 2 || mat(0,3) != 2 || mat(0,4) != 2 ||
             mat(1,0) != 2 || mat(1,1) != 2 || mat(1,2) != 2 || mat(1,3) != 2 || mat(1,4) != 2 ||
             mat(2,0) != 2 || mat(2,1) != 2 || mat(2,2) != 2 || mat(2,3) != 2 || mat(2,4) != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Transpose operation failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 2 2 2 2 2 )\n( 2 2 2 2 2 )\n( 2 2 2 2 2 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   {
      test_ = "Row-major self-transpose via trans()";

      // Self-transpose of a 3x5 matrix
      {
         blaze::UniformMatrix<int,blaze::rowMajor> mat( 3UL, 5UL, 2 );

         mat = trans( mat );

         checkRows    ( mat,  5UL );
         checkColumns ( mat,  3UL );
         checkCapacity( mat, 15UL );
         checkNonZeros( mat, 15UL );
         checkNonZeros( mat,  0UL, 3UL );
         checkNonZeros( mat,  1UL, 3UL );
         checkNonZeros( mat,  2UL, 3UL );
         checkNonZeros( mat,  3UL, 3UL );
         checkNonZeros( mat,  4UL, 3UL );

         if( mat(0,0) != 2 || mat(0,1) != 2 || mat(0,2) != 2 ||
             mat(1,0) != 2 || mat(1,1) != 2 || mat(1,2) != 2 ||
             mat(2,0) != 2 || mat(2,1) != 2 || mat(2,2) != 2 ||
             mat(3,0) != 2 || mat(3,1) != 2 || mat(3,2) != 2 ||
             mat(4,0) != 2 || mat(4,1) != 2 || mat(4,2) != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Transpose operation failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n( 2 2 2 )\n( 2 2 2 )\n( 2 2 2 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Self-transpose of a 5x3 matrix
      {
         blaze::UniformMatrix<int,blaze::rowMajor> mat( 5UL, 3UL, 2 );

         mat = trans( mat );

         checkRows    ( mat,  3UL );
         checkColumns ( mat,  5UL );
         checkCapacity( mat, 15UL );
         checkNonZeros( mat, 15UL );
         checkNonZeros( mat,  0UL, 5UL );
         checkNonZeros( mat,  1UL, 5UL );
         checkNonZeros( mat,  2UL, 5UL );

         if( mat(0,0) != 2 || mat(0,1) != 2 || mat(0,2) != 2 || mat(0,3) != 2 || mat(0,4) != 2 ||
             mat(1,0) != 2 || mat(1,1) != 2 || mat(1,2) != 2 || mat(1,3) != 2 || mat(1,4) != 2 ||
             mat(2,0) != 2 || mat(2,1) != 2 || mat(2,2) != 2 || mat(2,3) != 2 || mat(2,4) != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Transpose operation failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 2 2 2 2 2 )\n( 2 2 2 2 2 )\n( 2 2 2 2 2 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   {
      test_ = "Row-major self-transpose (stress test)";

      const size_t m( blaze::rand<size_t>( 0UL, 100UL ) );
      const size_t n( blaze::rand<size_t>( 0UL, 100UL ) );

      blaze::UniformMatrix<int,blaze::rowMajor> mat1( m, n );
      randomize( mat1 );
      blaze::UniformMatrix<int,blaze::rowMajor> mat2( mat1 );

      transpose( mat1 );

      if( mat1 != trans( mat2 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Transpose operation failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << trans( mat2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major self-transpose via transpose()";

      // Self-transpose of a 3x5 matrix
      {
         blaze::UniformMatrix<int,blaze::columnMajor> mat( 3UL, 5UL, 2 );

         transpose( mat );

         checkRows    ( mat,  5UL );
         checkColumns ( mat,  3UL );
         checkCapacity( mat, 15UL );
         checkNonZeros( mat, 15UL );
         checkNonZeros( mat,  0UL, 5UL );
         checkNonZeros( mat,  1UL, 5UL );
         checkNonZeros( mat,  2UL, 5UL );

         if( mat(0,0) != 2 || mat(0,1) != 2 || mat(0,2) != 2 ||
             mat(1,0) != 2 || mat(1,1) != 2 || mat(1,2) != 2 ||
             mat(2,0) != 2 || mat(2,1) != 2 || mat(2,2) != 2 ||
             mat(3,0) != 2 || mat(3,1) != 2 || mat(3,2) != 2 ||
             mat(4,0) != 2 || mat(4,1) != 2 || mat(4,2) != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Transpose operation failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n( 2 2 2 )\n( 2 2 2 )\n( 2 2 2 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Self-transpose of a 5x3 matrix
      {
         blaze::UniformMatrix<int,blaze::columnMajor> mat( 5UL, 3UL, 2 );

         transpose( mat );

         checkRows    ( mat,  3UL );
         checkColumns ( mat,  5UL );
         checkCapacity( mat, 15UL );
         checkNonZeros( mat, 15UL );
         checkNonZeros( mat,  0UL, 3UL );
         checkNonZeros( mat,  1UL, 3UL );
         checkNonZeros( mat,  2UL, 3UL );
         checkNonZeros( mat,  2UL, 3UL );
         checkNonZeros( mat,  2UL, 3UL );

         if( mat(0,0) != 2 || mat(0,1) != 2 || mat(0,2) != 2 || mat(0,3) != 2 || mat(0,4) != 2 ||
             mat(1,0) != 2 || mat(1,1) != 2 || mat(1,2) != 2 || mat(1,3) != 2 || mat(1,4) != 2 ||
             mat(2,0) != 2 || mat(2,1) != 2 || mat(2,2) != 2 || mat(2,3) != 2 || mat(2,4) != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Transpose operation failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 2 2 2 2 2 )\n( 2 2 2 2 2 )\n( 2 2 2 2 2 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   {
      test_ = "Column-major self-transpose via trans()";

      // Self-transpose of a 3x5 matrix
      {
         blaze::UniformMatrix<int,blaze::columnMajor> mat( 3UL, 5UL, 2 );

         mat = trans( mat );

         checkRows    ( mat,  5UL );
         checkColumns ( mat,  3UL );
         checkCapacity( mat, 15UL );
         checkNonZeros( mat, 15UL );
         checkNonZeros( mat,  0UL, 5UL );
         checkNonZeros( mat,  1UL, 5UL );
         checkNonZeros( mat,  2UL, 5UL );

         if( mat(0,0) != 2 || mat(0,1) != 2 || mat(0,2) != 2 ||
             mat(1,0) != 2 || mat(1,1) != 2 || mat(1,2) != 2 ||
             mat(2,0) != 2 || mat(2,1) != 2 || mat(2,2) != 2 ||
             mat(3,0) != 2 || mat(3,1) != 2 || mat(3,2) != 2 ||
             mat(4,0) != 2 || mat(4,1) != 2 || mat(4,2) != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Transpose operation failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n( 2 2 2 )\n( 2 2 2 )\n( 2 2 2 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Self-transpose of a 5x3 matrix
      {
         blaze::UniformMatrix<int,blaze::columnMajor> mat( 5UL, 3UL, 2 );

         mat = trans( mat );

         checkRows    ( mat,  3UL );
         checkColumns ( mat,  5UL );
         checkCapacity( mat, 15UL );
         checkNonZeros( mat, 15UL );
         checkNonZeros( mat,  0UL, 3UL );
         checkNonZeros( mat,  1UL, 3UL );
         checkNonZeros( mat,  2UL, 3UL );
         checkNonZeros( mat,  2UL, 3UL );
         checkNonZeros( mat,  2UL, 3UL );

         if( mat(0,0) != 2 || mat(0,1) != 2 || mat(0,2) != 2 || mat(0,3) != 2 || mat(0,4) != 2 ||
             mat(1,0) != 2 || mat(1,1) != 2 || mat(1,2) != 2 || mat(1,3) != 2 || mat(1,4) != 2 ||
             mat(2,0) != 2 || mat(2,1) != 2 || mat(2,2) != 2 || mat(2,3) != 2 || mat(2,4) != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Transpose operation failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( 2 2 2 2 2 )\n( 2 2 2 2 2 )\n( 2 2 2 2 2 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   {
      test_ = "Column-major self-transpose (stress test)";

      const size_t m( blaze::rand<size_t>( 0UL, 100UL ) );
      const size_t n( blaze::rand<size_t>( 0UL, 100UL ) );

      blaze::UniformMatrix<int,blaze::columnMajor> mat1( m, n );
      randomize( mat1 );
      blaze::UniformMatrix<int,blaze::columnMajor> mat2( mat1 );

      transpose( mat1 );

      if( mat1 != trans( mat2 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Transpose operation failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << trans( mat2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c ctranspose() member function of the UniformMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c ctranspose() member function of the UniformMatrix
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

      // Self-transpose of a 4x4 matrix
      {
         blaze::UniformMatrix<cplx,blaze::rowMajor> mat( 4UL, 4UL, cplx(2,-2) );

         ctranspose( mat );

         checkRows    ( mat,  4UL );
         checkColumns ( mat,  4UL );
         checkCapacity( mat, 16UL );
         checkNonZeros( mat, 16UL );
         checkNonZeros( mat,  0UL, 4UL );
         checkNonZeros( mat,  1UL, 4UL );
         checkNonZeros( mat,  2UL, 4UL );
         checkNonZeros( mat,  3UL, 4UL );

         if( mat(0,0) != cplx(2,2) || mat(0,1) != cplx(2,2) || mat(0,2) != cplx(2,2) || mat(0,3) != cplx(2,2) ||
             mat(1,0) != cplx(2,2) || mat(1,1) != cplx(2,2) || mat(1,2) != cplx(2,2) || mat(1,3) != cplx(2,2) ||
             mat(2,0) != cplx(2,2) || mat(2,1) != cplx(2,2) || mat(2,2) != cplx(2,2) || mat(2,3) != cplx(2,2) ||
             mat(3,0) != cplx(2,2) || mat(3,1) != cplx(2,2) || mat(3,2) != cplx(2,2) || mat(3,3) != cplx(2,2) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Transpose operation failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( (2,2) (2,2) (2,2) (2,2) )\n"
                                        "( (2,2) (2,2) (2,2) (2,2) )\n"
                                        "( (2,2) (2,2) (2,2) (2,2) )\n"
                                        "( (2,2) (2,2) (2,2) (2,2) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Self-transpose of a 3x5 matrix
      {
         blaze::UniformMatrix<cplx,blaze::rowMajor> mat( 3UL, 5UL, cplx(2,-2) );

         ctranspose( mat );

         checkRows    ( mat,  5UL );
         checkColumns ( mat,  3UL );
         checkCapacity( mat, 15UL );
         checkNonZeros( mat, 15UL );
         checkNonZeros( mat,  0UL, 3UL );
         checkNonZeros( mat,  1UL, 3UL );
         checkNonZeros( mat,  2UL, 3UL );
         checkNonZeros( mat,  3UL, 3UL );
         checkNonZeros( mat,  4UL, 3UL );

         if( mat(0,0) != cplx(2,2) || mat(0,1) != cplx(2,2) || mat(0,2) != cplx(2,2) ||
             mat(1,0) != cplx(2,2) || mat(1,1) != cplx(2,2) || mat(1,2) != cplx(2,2) ||
             mat(2,0) != cplx(2,2) || mat(2,1) != cplx(2,2) || mat(2,2) != cplx(2,2) ||
             mat(3,0) != cplx(2,2) || mat(3,1) != cplx(2,2) || mat(3,2) != cplx(2,2) ||
             mat(4,0) != cplx(2,2) || mat(4,1) != cplx(2,2) || mat(4,2) != cplx(2,2) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Transpose operation failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( (2,2) (2,2) (2,2) )\n"
                                        "( (2,2) (2,2) (2,2) )\n"
                                        "( (2,2) (2,2) (2,2) )\n"
                                        "( (2,2) (2,2) (2,2) )\n"
                                        "( (2,2) (2,2) (2,2) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Self-transpose of a 5x3 matrix
      {
         blaze::UniformMatrix<cplx,blaze::rowMajor> mat( 5UL, 3UL, cplx(2,-2) );

         ctranspose( mat );

         checkRows    ( mat,  3UL );
         checkColumns ( mat,  5UL );
         checkCapacity( mat, 15UL );
         checkNonZeros( mat, 15UL );
         checkNonZeros( mat,  0UL, 5UL );
         checkNonZeros( mat,  1UL, 5UL );
         checkNonZeros( mat,  2UL, 5UL );

         if( mat(0,0) != cplx(2,2) || mat(0,1) != cplx(2,2) || mat(0,2) != cplx(2,2) || mat(0,3) != cplx(2,2) || mat(0,4) != cplx(2,2) ||
             mat(1,0) != cplx(2,2) || mat(1,1) != cplx(2,2) || mat(1,2) != cplx(2,2) || mat(1,3) != cplx(2,2) || mat(1,4) != cplx(2,2) ||
             mat(2,0) != cplx(2,2) || mat(2,1) != cplx(2,2) || mat(2,2) != cplx(2,2) || mat(2,3) != cplx(2,2) || mat(2,4) != cplx(2,2) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Transpose operation failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( (2,2) (2,2) (2,2) (2,2) (2,2) )\n"
                                        "( (2,2) (2,2) (2,2) (2,2) (2,2) )\n"
                                        "( (2,2) (2,2) (2,2) (2,2) (2,2) )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   {
      test_ = "Row-major self-transpose via ctranspose() (stress test)";

      using cplx = blaze::complex<int>;

      const size_t m( blaze::rand<size_t>( 0UL, 100UL ) );
      const size_t n( blaze::rand<size_t>( 0UL, 100UL ) );

      blaze::UniformMatrix<cplx,blaze::rowMajor> mat1( m, n );
      randomize( mat1 );
      blaze::UniformMatrix<cplx,blaze::rowMajor> mat2( mat1 );

      ctranspose( mat1 );

      if( mat1 != ctrans( mat2 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Transpose operation failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << ctrans( mat2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major self-transpose via ctrans()";

      using cplx = blaze::complex<int>;

      // Self-transpose of a 4x4 matrix
      {
         blaze::UniformMatrix<cplx,blaze::rowMajor> mat( 4UL, 4UL, cplx(2,-2) );

         mat = ctrans( mat );

         checkRows    ( mat,  4UL );
         checkColumns ( mat,  4UL );
         checkCapacity( mat, 16UL );
         checkNonZeros( mat, 16UL );
         checkNonZeros( mat,  0UL, 4UL );
         checkNonZeros( mat,  1UL, 4UL );
         checkNonZeros( mat,  2UL, 4UL );
         checkNonZeros( mat,  3UL, 4UL );

         if( mat(0,0) != cplx(2,2) || mat(0,1) != cplx(2,2) || mat(0,2) != cplx(2,2) || mat(0,3) != cplx(2,2) ||
             mat(1,0) != cplx(2,2) || mat(1,1) != cplx(2,2) || mat(1,2) != cplx(2,2) || mat(1,3) != cplx(2,2) ||
             mat(2,0) != cplx(2,2) || mat(2,1) != cplx(2,2) || mat(2,2) != cplx(2,2) || mat(2,3) != cplx(2,2) ||
             mat(3,0) != cplx(2,2) || mat(3,1) != cplx(2,2) || mat(3,2) != cplx(2,2) || mat(3,3) != cplx(2,2) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Transpose operation failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( (2,2) (2,2) (2,2) (2,2) )\n"
                                        "( (2,2) (2,2) (2,2) (2,2) )\n"
                                        "( (2,2) (2,2) (2,2) (2,2) )\n"
                                        "( (2,2) (2,2) (2,2) (2,2) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Self-transpose of a 3x5 matrix
      {
         blaze::UniformMatrix<cplx,blaze::rowMajor> mat( 3UL, 5UL, cplx(2,-2) );

         mat = ctrans( mat );

         checkRows    ( mat,  5UL );
         checkColumns ( mat,  3UL );
         checkCapacity( mat, 15UL );
         checkNonZeros( mat, 15UL );
         checkNonZeros( mat,  0UL, 3UL );
         checkNonZeros( mat,  1UL, 3UL );
         checkNonZeros( mat,  2UL, 3UL );
         checkNonZeros( mat,  3UL, 3UL );
         checkNonZeros( mat,  4UL, 3UL );

         if( mat(0,0) != cplx(2,2) || mat(0,1) != cplx(2,2) || mat(0,2) != cplx(2,2) ||
             mat(1,0) != cplx(2,2) || mat(1,1) != cplx(2,2) || mat(1,2) != cplx(2,2) ||
             mat(2,0) != cplx(2,2) || mat(2,1) != cplx(2,2) || mat(2,2) != cplx(2,2) ||
             mat(3,0) != cplx(2,2) || mat(3,1) != cplx(2,2) || mat(3,2) != cplx(2,2) ||
             mat(4,0) != cplx(2,2) || mat(4,1) != cplx(2,2) || mat(4,2) != cplx(2,2) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Transpose operation failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( (2,2) (2,2) (2,2) )\n"
                                        "( (2,2) (2,2) (2,2) )\n"
                                        "( (2,2) (2,2) (2,2) )\n"
                                        "( (2,2) (2,2) (2,2) )\n"
                                        "( (2,2) (2,2) (2,2) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Self-transpose of a 5x3 matrix
      {
         blaze::UniformMatrix<cplx,blaze::rowMajor> mat( 5UL, 3UL, cplx(2,-2) );

         mat = ctrans( mat );

         checkRows    ( mat,  3UL );
         checkColumns ( mat,  5UL );
         checkCapacity( mat, 15UL );
         checkNonZeros( mat, 15UL );
         checkNonZeros( mat,  0UL, 5UL );
         checkNonZeros( mat,  1UL, 5UL );
         checkNonZeros( mat,  2UL, 5UL );

         if( mat(0,0) != cplx(2,2) || mat(0,1) != cplx(2,2) || mat(0,2) != cplx(2,2) || mat(0,3) != cplx(2,2) || mat(0,4) != cplx(2,2) ||
             mat(1,0) != cplx(2,2) || mat(1,1) != cplx(2,2) || mat(1,2) != cplx(2,2) || mat(1,3) != cplx(2,2) || mat(1,4) != cplx(2,2) ||
             mat(2,0) != cplx(2,2) || mat(2,1) != cplx(2,2) || mat(2,2) != cplx(2,2) || mat(2,3) != cplx(2,2) || mat(2,4) != cplx(2,2) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Transpose operation failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( (2,2) (2,2) (2,2) (2,2) (2,2) )\n"
                                        "( (2,2) (2,2) (2,2) (2,2) (2,2) )\n"
                                        "( (2,2) (2,2) (2,2) (2,2) (2,2) )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   {
      test_ = "Row-major self-transpose via ctrans() (stress test)";

      using cplx = blaze::complex<int>;

      const size_t m( blaze::rand<size_t>( 0UL, 100UL ) );
      const size_t n( blaze::rand<size_t>( 0UL, 100UL ) );

      blaze::UniformMatrix<cplx,blaze::rowMajor> mat1( m, n );
      randomize( mat1 );
      blaze::UniformMatrix<cplx,blaze::rowMajor> mat2( mat1 );

      mat1 = ctrans( mat1 );

      if( mat1 != ctrans( mat2 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Transpose operation failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << ctrans( mat2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      test_ = "Column-major self-transpose via ctranspose()";

      using cplx = blaze::complex<int>;

      // Self-transpose of a 4x4 matrix
      {
         blaze::UniformMatrix<cplx,blaze::columnMajor> mat( 4UL, 4UL, cplx(2,-2) );

         ctranspose( mat );

         checkRows    ( mat,  4UL );
         checkColumns ( mat,  4UL );
         checkCapacity( mat, 16UL );
         checkNonZeros( mat, 16UL );
         checkNonZeros( mat,  0UL, 4UL );
         checkNonZeros( mat,  1UL, 4UL );
         checkNonZeros( mat,  2UL, 4UL );
         checkNonZeros( mat,  3UL, 4UL );

         if( mat(0,0) != cplx(2,2) || mat(0,1) != cplx(2,2) || mat(0,2) != cplx(2,2) || mat(0,3) != cplx(2,2) ||
             mat(1,0) != cplx(2,2) || mat(1,1) != cplx(2,2) || mat(1,2) != cplx(2,2) || mat(1,3) != cplx(2,2) ||
             mat(2,0) != cplx(2,2) || mat(2,1) != cplx(2,2) || mat(2,2) != cplx(2,2) || mat(2,3) != cplx(2,2) ||
             mat(3,0) != cplx(2,2) || mat(3,1) != cplx(2,2) || mat(3,2) != cplx(2,2) || mat(3,3) != cplx(2,2) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Transpose operation failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( (2,2) (2,2) (2,2) (2,2) )\n"
                                        "( (2,2) (2,2) (2,2) (2,2) )\n"
                                        "( (2,2) (2,2) (2,2) (2,2) )\n"
                                        "( (2,2) (2,2) (2,2) (2,2) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Self-transpose of a 3x5 matrix
      {
         blaze::UniformMatrix<cplx,blaze::columnMajor> mat( 3UL, 5UL, cplx(2,-2) );

         ctranspose( mat );

         checkRows    ( mat,  5UL );
         checkColumns ( mat,  3UL );
         checkCapacity( mat, 15UL );
         checkNonZeros( mat, 15UL );
         checkNonZeros( mat,  0UL, 5UL );
         checkNonZeros( mat,  1UL, 5UL );
         checkNonZeros( mat,  2UL, 5UL );

         if( mat(0,0) != cplx(2,2) || mat(0,1) != cplx(2,2) || mat(0,2) != cplx(2,2) ||
             mat(1,0) != cplx(2,2) || mat(1,1) != cplx(2,2) || mat(1,2) != cplx(2,2) ||
             mat(2,0) != cplx(2,2) || mat(2,1) != cplx(2,2) || mat(2,2) != cplx(2,2) ||
             mat(3,0) != cplx(2,2) || mat(3,1) != cplx(2,2) || mat(3,2) != cplx(2,2) ||
             mat(4,0) != cplx(2,2) || mat(4,1) != cplx(2,2) || mat(4,2) != cplx(2,2) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Transpose operation failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( (2,2) (2,2) (2,2) )\n"
                                        "( (2,2) (2,2) (2,2) )\n"
                                        "( (2,2) (2,2) (2,2) )\n"
                                        "( (2,2) (2,2) (2,2) )\n"
                                        "( (2,2) (2,2) (2,2) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Self-transpose of a 5x3 matrix
      {
         blaze::UniformMatrix<cplx,blaze::columnMajor> mat( 5UL, 3UL, cplx(2,-2) );

         ctranspose( mat );

         checkRows    ( mat,  3UL );
         checkColumns ( mat,  5UL );
         checkCapacity( mat, 15UL );
         checkNonZeros( mat, 15UL );
         checkNonZeros( mat,  0UL, 3UL );
         checkNonZeros( mat,  1UL, 3UL );
         checkNonZeros( mat,  2UL, 3UL );
         checkNonZeros( mat,  3UL, 3UL );
         checkNonZeros( mat,  4UL, 3UL );

         if( mat(0,0) != cplx(2,2) || mat(0,1) != cplx(2,2) || mat(0,2) != cplx(2,2) || mat(0,3) != cplx(2,2) || mat(0,4) != cplx(2,2) ||
             mat(1,0) != cplx(2,2) || mat(1,1) != cplx(2,2) || mat(1,2) != cplx(2,2) || mat(1,3) != cplx(2,2) || mat(1,4) != cplx(2,2) ||
             mat(2,0) != cplx(2,2) || mat(2,1) != cplx(2,2) || mat(2,2) != cplx(2,2) || mat(2,3) != cplx(2,2) || mat(2,4) != cplx(2,2) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Transpose operation failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( (2,2) (2,2) (2,2) (2,2) (2,2) )\n"
                                        "( (2,2) (2,2) (2,2) (2,2) (2,2) )\n"
                                        "( (2,2) (2,2) (2,2) (2,2) (2,2) )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   {
      test_ = "Column-major self-transpose via ctranspose() (stress test)";

      using cplx = blaze::complex<int>;

      const size_t m( blaze::rand<size_t>( 0UL, 100UL ) );
      const size_t n( blaze::rand<size_t>( 0UL, 100UL ) );

      blaze::UniformMatrix<cplx,blaze::columnMajor> mat1( m, n );
      randomize( mat1 );
      blaze::UniformMatrix<cplx,blaze::columnMajor> mat2( mat1 );

      ctranspose( mat1 );

      if( mat1 != ctrans( mat2 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Transpose operation failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << ctrans( mat2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major self-transpose via ctrans()";

      using cplx = blaze::complex<int>;

      // Self-transpose of a 4x4 matrix
      {
         blaze::UniformMatrix<cplx,blaze::columnMajor> mat( 4UL, 4UL, cplx(2,-2) );

         mat = ctrans( mat );

         checkRows    ( mat,  4UL );
         checkColumns ( mat,  4UL );
         checkCapacity( mat, 16UL );
         checkNonZeros( mat, 16UL );
         checkNonZeros( mat,  0UL, 4UL );
         checkNonZeros( mat,  1UL, 4UL );
         checkNonZeros( mat,  2UL, 4UL );
         checkNonZeros( mat,  3UL, 4UL );

         if( mat(0,0) != cplx(2,2) || mat(0,1) != cplx(2,2) || mat(0,2) != cplx(2,2) || mat(0,3) != cplx(2,2) ||
             mat(1,0) != cplx(2,2) || mat(1,1) != cplx(2,2) || mat(1,2) != cplx(2,2) || mat(1,3) != cplx(2,2) ||
             mat(2,0) != cplx(2,2) || mat(2,1) != cplx(2,2) || mat(2,2) != cplx(2,2) || mat(2,3) != cplx(2,2) ||
             mat(3,0) != cplx(2,2) || mat(3,1) != cplx(2,2) || mat(3,2) != cplx(2,2) || mat(3,3) != cplx(2,2) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Transpose operation failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( (2,2) (2,2) (2,2) (2,2) )\n"
                                        "( (2,2) (2,2) (2,2) (2,2) )\n"
                                        "( (2,2) (2,2) (2,2) (2,2) )\n"
                                        "( (2,2) (2,2) (2,2) (2,2) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Self-transpose of a 3x5 matrix
      {
         blaze::UniformMatrix<cplx,blaze::columnMajor> mat( 3UL, 5UL, cplx(2,-2) );

         mat = ctrans( mat );

         checkRows    ( mat,  5UL );
         checkColumns ( mat,  3UL );
         checkCapacity( mat, 15UL );
         checkNonZeros( mat, 15UL );
         checkNonZeros( mat,  0UL, 5UL );
         checkNonZeros( mat,  1UL, 5UL );
         checkNonZeros( mat,  2UL, 5UL );

         if( mat(0,0) != cplx(2,2) || mat(0,1) != cplx(2,2) || mat(0,2) != cplx(2,2) ||
             mat(1,0) != cplx(2,2) || mat(1,1) != cplx(2,2) || mat(1,2) != cplx(2,2) ||
             mat(2,0) != cplx(2,2) || mat(2,1) != cplx(2,2) || mat(2,2) != cplx(2,2) ||
             mat(3,0) != cplx(2,2) || mat(3,1) != cplx(2,2) || mat(3,2) != cplx(2,2) ||
             mat(4,0) != cplx(2,2) || mat(4,1) != cplx(2,2) || mat(4,2) != cplx(2,2) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Transpose operation failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( (2,2) (2,2) (2,2) )\n"
                                        "( (2,2) (2,2) (2,2) )\n"
                                        "( (2,2) (2,2) (2,2) )\n"
                                        "( (2,2) (2,2) (2,2) )\n"
                                        "( (2,2) (2,2) (2,2) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Self-transpose of a 5x3 matrix
      {
         blaze::UniformMatrix<cplx,blaze::columnMajor> mat( 5UL, 3UL, cplx(2,-2) );

         mat = ctrans( mat );

         checkRows    ( mat,  3UL );
         checkColumns ( mat,  5UL );
         checkCapacity( mat, 15UL );
         checkNonZeros( mat, 15UL );
         checkNonZeros( mat,  0UL, 3UL );
         checkNonZeros( mat,  1UL, 3UL );
         checkNonZeros( mat,  2UL, 3UL );
         checkNonZeros( mat,  3UL, 3UL );
         checkNonZeros( mat,  4UL, 3UL );

         if( mat(0,0) != cplx(2,2) || mat(0,1) != cplx(2,2) || mat(0,2) != cplx(2,2) || mat(0,3) != cplx(2,2) || mat(0,4) != cplx(2,2) ||
             mat(1,0) != cplx(2,2) || mat(1,1) != cplx(2,2) || mat(1,2) != cplx(2,2) || mat(1,3) != cplx(2,2) || mat(1,4) != cplx(2,2) ||
             mat(2,0) != cplx(2,2) || mat(2,1) != cplx(2,2) || mat(2,2) != cplx(2,2) || mat(2,3) != cplx(2,2) || mat(2,4) != cplx(2,2) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Transpose operation failed\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n"
                << "   Expected result:\n( (2,2) (2,2) (2,2) (2,2) (2,2) )\n"
                                        "( (2,2) (2,2) (2,2) (2,2) (2,2) )\n"
                                        "( (2,2) (2,2) (2,2) (2,2) (2,2) )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   {
      test_ = "Column-major self-transpose via ctrans() (stress test)";

      using cplx = blaze::complex<int>;

      const size_t m( blaze::rand<size_t>( 0UL, 100UL ) );
      const size_t n( blaze::rand<size_t>( 0UL, 100UL ) );

      blaze::UniformMatrix<cplx,blaze::columnMajor> mat1( m, n );
      randomize( mat1 );
      blaze::UniformMatrix<cplx,blaze::columnMajor> mat2( mat1 );

      mat1 = ctrans( mat1 );

      if( mat1 != ctrans( mat2 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Transpose operation failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << ctrans( mat2 ) << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c isDefault() function with the UniformMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isDefault() function with the UniformMatrix class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testIsDefault()
{
   using blaze::isDefault;


   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      test_ = "Row-major isDefault() function";

      // isDefault with 0x0 matrix (default)
      {
         blaze::UniformMatrix<int,blaze::rowMajor> mat;

         if( isDefault( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with 0x3 matrix (non-default)
      {
         blaze::UniformMatrix<int,blaze::rowMajor> mat( 0UL, 3UL );

         if( isDefault( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with 2x0 matrix (non-default)
      {
         blaze::UniformMatrix<int,blaze::rowMajor> mat( 2UL, 0UL );

         if( isDefault( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with 2x3 matrix (non-default)
      {
         blaze::UniformMatrix<int,blaze::rowMajor> mat( 2UL, 3UL, 0 );

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

      // isDefault with 3x2 matrix (non-default)
      {
         blaze::UniformMatrix<int,blaze::rowMajor> mat( 3UL, 2UL, 1 );

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

      // isDefault with 0x0 matrix (default)
      {
         blaze::UniformMatrix<int,blaze::columnMajor> mat;

         if( isDefault( mat ) != true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with 0x3 matrix (non-default)
      {
         blaze::UniformMatrix<int,blaze::columnMajor> mat( 0UL, 3UL );

         if( isDefault( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with 2x0 matrix (non-default)
      {
         blaze::UniformMatrix<int,blaze::columnMajor> mat( 2UL, 0UL );

         if( isDefault( mat ) != false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isDefault evaluation\n"
                << " Details:\n"
                << "   Matrix:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isDefault with 2x3 matrix (non-default)
      {
         blaze::UniformMatrix<int,blaze::columnMajor> mat( 2UL, 3UL, 0 );

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

      // isDefault with 3x2 matrix (non-default)
      {
         blaze::UniformMatrix<int,blaze::columnMajor> mat( 3UL, 2UL, 1 );

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

} // namespace uniformmatrix

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
   std::cout << "   Running UniformMatrix class test (part 2)..." << std::endl;

   try
   {
      RUN_UNIFORMMATRIX_CLASS_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during UniformMatrix class test (part 2):\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
