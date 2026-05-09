//=================================================================================================
/*!
//  \file src/mathtest/matrices/uniformmatrix/ClassTest1.cpp
//  \brief Source file for the UniformMatrix class test (part 1)
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
   testConstructors();
   testAssignment();
   testAddAssign();
   testSubAssign();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the UniformMatrix constructors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all constructors of the UniformMatrix class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testConstructors()
{
   //=====================================================================================
   // Row-major default constructor
   //=====================================================================================

   // Default constructor
   {
      test_ = "Row-major UniformMatrix default constructor";

      blaze::UniformMatrix<int,blaze::rowMajor> mat;

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );
   }


   //=====================================================================================
   // Row-major size constructor
   //=====================================================================================

   {
      test_ = "Row-major UniformMatrix size constructor (0x0)";

      blaze::UniformMatrix<int,blaze::rowMajor> mat( 0UL, 0UL );

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );
   }

   {
      test_ = "Row-major UniformMatrix size constructor (0x4)";

      blaze::UniformMatrix<int,blaze::rowMajor> mat( 0UL, 4UL );

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 4UL );
      checkNonZeros( mat, 0UL );
   }

   {
      test_ = "Row-major UniformMatrix size constructor (3x0)";

      blaze::UniformMatrix<int,blaze::rowMajor> mat( 3UL, 0UL );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );
   }

   {
      test_ = "Row-major UniformMatrix size constructor (3x4)";

      blaze::UniformMatrix<int,blaze::rowMajor> mat( 3UL, 4UL );

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  4UL );
      checkCapacity( mat, 12UL );
   }


   //=====================================================================================
   // Row-major homogeneous initialization
   //=====================================================================================

   {
      test_ = "Row-major UniformMatrix homogeneous initialization constructor (0x0)";

      blaze::UniformMatrix<int,blaze::rowMajor> mat( 0UL, 0UL, 2 );

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );
   }

   {
      test_ = "Row-major UniformMatrix homogeneous initialization constructor (0x4)";

      blaze::UniformMatrix<int,blaze::rowMajor> mat( 0UL, 4UL, 2 );

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 4UL );
      checkNonZeros( mat, 0UL );
   }

   {
      test_ = "Row-major UniformMatrix homogeneous initialization constructor (3x0)";

      blaze::UniformMatrix<int,blaze::rowMajor> mat( 3UL, 0UL, 2 );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );
   }

   {
      test_ = "Row-major UniformMatrix homogeneous initialization constructor (3x4)";

      blaze::UniformMatrix<int,blaze::rowMajor> mat( 3UL, 4UL, 2 );

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  4UL );
      checkCapacity( mat, 12UL );
      checkNonZeros( mat, 12UL );
      checkNonZeros( mat,  0UL, 4UL );
      checkNonZeros( mat,  1UL, 4UL );
      checkNonZeros( mat,  2UL, 4UL );

      if( mat(0,0) != 2 || mat(0,1) != 2 || mat(0,2) != 2 || mat(0,3) != 2 ||
          mat(1,0) != 2 || mat(1,1) != 2 || mat(1,2) != 2 || mat(1,3) != 2 ||
          mat(2,0) != 2 || mat(2,1) != 2 || mat(2,2) != 2 || mat(2,3) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 2 2 2 2 )\n( 2 2 2 2 )\n( 2 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major copy constructor
   //=====================================================================================

   {
      test_ = "Row-major UniformMatrix copy constructor (0x0)";

      blaze::UniformMatrix<int,blaze::rowMajor> mat1( 0UL, 0UL );
      blaze::UniformMatrix<int,blaze::rowMajor> mat2( mat1 );

      checkRows    ( mat2, 0UL );
      checkColumns ( mat2, 0UL );
      checkNonZeros( mat2, 0UL );
   }

   {
      test_ = "Row-major UniformMatrix copy constructor (0x3)";

      blaze::UniformMatrix<int,blaze::rowMajor> mat1( 0UL, 3UL );
      blaze::UniformMatrix<int,blaze::rowMajor> mat2( mat1 );

      checkRows    ( mat2, 0UL );
      checkColumns ( mat2, 3UL );
      checkNonZeros( mat2, 0UL );
   }

   {
      test_ = "Row-major UniformMatrix copy constructor (2x0)";

      blaze::UniformMatrix<int,blaze::rowMajor> mat1( 2UL, 0UL );
      blaze::UniformMatrix<int,blaze::rowMajor> mat2( mat1 );

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 0UL );
      checkNonZeros( mat2, 0UL );
   }

   {
      test_ = "Row-major UniformMatrix copy constructor (2x3)";

      blaze::UniformMatrix<int,blaze::rowMajor> mat1( 2UL, 3UL, 2 );
      blaze::UniformMatrix<int,blaze::rowMajor> mat2( mat1 );

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
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major move constructor
   //=====================================================================================

   {
      test_ = "Row-major UniformMatrix move constructor (0x0)";

      blaze::UniformMatrix<int,blaze::rowMajor> mat1( 0UL, 0UL );
      blaze::UniformMatrix<int,blaze::rowMajor> mat2( std::move( mat1 ) );

      checkRows    ( mat2, 0UL );
      checkColumns ( mat2, 0UL );
      checkNonZeros( mat2, 0UL );
   }

   {
      test_ = "Row-major UniformMatrix move constructor (0x3)";

      blaze::UniformMatrix<int,blaze::rowMajor> mat1( 0UL, 3UL );
      blaze::UniformMatrix<int,blaze::rowMajor> mat2( std::move( mat1 ) );

      checkRows    ( mat2, 0UL );
      checkColumns ( mat2, 3UL );
      checkNonZeros( mat2, 0UL );
   }

   {
      test_ = "Row-major UniformMatrix move constructor (2x0)";

      blaze::UniformMatrix<int,blaze::rowMajor> mat1( 2UL, 0UL );
      blaze::UniformMatrix<int,blaze::rowMajor> mat2( std::move( mat1 ) );

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 0UL );
      checkNonZeros( mat2, 0UL );
   }

   {
      test_ = "Row-major UniformMatrix copy constructor (2x3)";

      blaze::UniformMatrix<int,blaze::rowMajor> mat1( 2UL, 3UL, 2 );
      blaze::UniformMatrix<int,blaze::rowMajor> mat2( std::move( mat1 ) );

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
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major dense matrix constructor
   //=====================================================================================

   {
      test_ = "Row-major/row-major UniformMatrix dense matrix constructor (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::rowMajor;

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,rowMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 32UL ) );
      AlignedPadded mat1( memory.get(), 2UL, 3UL, 16UL );
      mat1(0,0) = 2;
      mat1(0,1) = 2;
      mat1(0,2) = 2;
      mat1(1,0) = 2;
      mat1(1,1) = 2;
      mat1(1,2) = 2;

      blaze::UniformMatrix<int,blaze::rowMajor> mat2( mat1 );

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
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major UniformMatrix dense matrix constructor (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowMajor;

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>;
      std::unique_ptr<int[]> memory( new int[7UL] );
      UnalignedUnpadded mat1( memory.get()+1UL, 2UL, 3UL );
      mat1(0,0) = 2;
      mat1(0,1) = 2;
      mat1(0,2) = 2;
      mat1(1,0) = 2;
      mat1(1,1) = 2;
      mat1(1,2) = 2;

      blaze::UniformMatrix<int,blaze::rowMajor> mat2( mat1 );

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
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major UniformMatrix dense matrix constructor (non-uniform)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat1{ { 2, 2, 2 }, { 2, 0, 2 } };

      try {
         blaze::UniformMatrix<int,blaze::rowMajor> mat2( mat1 );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of non-uniform UniformMatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   {
      test_ = "Row-major/column-major UniformMatrix dense matrix constructor (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::columnMajor;

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,columnMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      AlignedPadded mat1( memory.get(), 2UL, 3UL, 16UL );
      mat1(0,0) = 2;
      mat1(0,1) = 2;
      mat1(0,2) = 2;
      mat1(1,0) = 2;
      mat1(1,1) = 2;
      mat1(1,2) = 2;

      blaze::UniformMatrix<int,blaze::rowMajor> mat2( mat1 );

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
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major UniformMatrix dense matrix constructor (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::columnMajor;

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>;
      std::unique_ptr<int[]> memory( new int[7UL] );
      UnalignedUnpadded mat1( memory.get()+1UL, 2UL, 3UL );
      mat1(0,0) = 2;
      mat1(0,1) = 2;
      mat1(0,2) = 2;
      mat1(1,0) = 2;
      mat1(1,1) = 2;
      mat1(1,2) = 2;

      blaze::UniformMatrix<int,blaze::rowMajor> mat2( mat1 );

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
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major UniformMatrix dense matrix constructor (non-uniform)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat1{ { 2, 2, 2 }, { 2, 0, 2 } };

      try {
         blaze::UniformMatrix<int,blaze::rowMajor> mat2( mat1 );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of non-uniform UniformMatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Row-major sparse matrix constructor
   //=====================================================================================

   {
      test_ = "Row-major/row-major UniformMatrix sparse matrix constructor (uniform)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1{ { 2, 2, 2 }, { 2, 2, 2 } };
      blaze::UniformMatrix<int,blaze::rowMajor> mat2( mat1 );

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
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major UniformMatrix sparse matrix constructor (non-uniform)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1{ { 2, 2, 2 }, { 2, 0, 2 } };

      try {
         blaze::UniformMatrix<int,blaze::rowMajor> mat2( mat1 );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of non-uniform UniformMatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   {
      test_ = "Row-major/column-major UniformMatrix sparse matrix constructor (uniform)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1{ { 2, 2, 2 }, { 2, 2, 2 } };
      blaze::UniformMatrix<int,blaze::rowMajor> mat2( mat1 );

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
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major UniformMatrix sparse matrix constructor (non-uniform)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1{ { 2, 2, 2 }, { 2, 0, 2 } };

      try {
         blaze::UniformMatrix<int,blaze::rowMajor> mat2( mat1 );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of non-uniform UniformMatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Column-major default constructor
   //=====================================================================================

   {
      test_ = "Column-major UniformMatrix default constructor";

      blaze::UniformMatrix<int,blaze::columnMajor> mat;

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );
   }


   //=====================================================================================
   // Column-major size constructor
   //=====================================================================================

   {
      test_ = "Column-major UniformMatrix size constructor (0x0)";

      blaze::UniformMatrix<int,blaze::columnMajor> mat( 0UL, 0UL );

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );
   }

   {
      test_ = "Column-major UniformMatrix size constructor (0x4)";

      blaze::UniformMatrix<int,blaze::columnMajor> mat( 0UL, 4UL );

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 4UL );
      checkNonZeros( mat, 0UL );
   }

   {
      test_ = "Column-major UniformMatrix size constructor (3x0)";

      blaze::UniformMatrix<int,blaze::columnMajor> mat( 3UL, 0UL );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );
   }

   {
      test_ = "Column-major UniformMatrix size constructor (3x4)";

      blaze::UniformMatrix<int,blaze::columnMajor> mat( 3UL, 4UL );

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  4UL );
      checkCapacity( mat, 12UL );
   }


   //=====================================================================================
   // Column-major homogeneous initialization
   //=====================================================================================

   {
      test_ = "Column-major UniformMatrix homogeneous initialization constructor (0x0)";

      blaze::UniformMatrix<int,blaze::columnMajor> mat( 0UL, 0UL, 2 );

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );
   }

   {
      test_ = "Column-major UniformMatrix homogeneous initialization constructor (0x4)";

      blaze::UniformMatrix<int,blaze::columnMajor> mat( 0UL, 4UL, 2 );

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 4UL );
      checkNonZeros( mat, 0UL );
   }

   {
      test_ = "Column-major UniformMatrix homogeneous initialization constructor (3x0)";

      blaze::UniformMatrix<int,blaze::columnMajor> mat( 3UL, 0UL, 2 );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );
   }

   {
      test_ = "Column-major UniformMatrix homogeneous initialization constructor (3x4)";

      blaze::UniformMatrix<int,blaze::columnMajor> mat( 3UL, 4UL, 2 );

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  4UL );
      checkCapacity( mat, 12UL );
      checkNonZeros( mat, 12UL );
      checkNonZeros( mat,  0UL, 3UL );
      checkNonZeros( mat,  1UL, 3UL );
      checkNonZeros( mat,  2UL, 3UL );
      checkNonZeros( mat,  3UL, 3UL );

      if( mat(0,0) != 2 || mat(0,1) != 2 || mat(0,2) != 2 || mat(0,3) != 2 ||
          mat(1,0) != 2 || mat(1,1) != 2 || mat(1,2) != 2 || mat(1,3) != 2 ||
          mat(2,0) != 2 || mat(2,1) != 2 || mat(2,2) != 2 || mat(2,3) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 2 2 2 2 )\n( 2 2 2 2 )\n( 2 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major copy constructor
   //=====================================================================================

   {
      test_ = "Column-major UniformMatrix copy constructor (0x0)";

      blaze::UniformMatrix<int,blaze::columnMajor> mat1( 0UL, 0UL );
      blaze::UniformMatrix<int,blaze::columnMajor> mat2( mat1 );

      checkRows    ( mat2, 0UL );
      checkColumns ( mat2, 0UL );
      checkNonZeros( mat2, 0UL );
   }

   {
      test_ = "Column-major UniformMatrix copy constructor (0x3)";

      blaze::UniformMatrix<int,blaze::columnMajor> mat1( 0UL, 3UL );
      blaze::UniformMatrix<int,blaze::columnMajor> mat2( mat1 );

      checkRows    ( mat2, 0UL );
      checkColumns ( mat2, 3UL );
      checkNonZeros( mat2, 0UL );
   }

   {
      test_ = "Column-major UniformMatrix copy constructor (2x0)";

      blaze::UniformMatrix<int,blaze::columnMajor> mat1( 2UL, 0UL );
      blaze::UniformMatrix<int,blaze::columnMajor> mat2( mat1 );

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 0UL );
      checkNonZeros( mat2, 0UL );
   }

   {
      test_ = "Column-major UniformMatrix copy constructor (2x3)";

      blaze::UniformMatrix<int,blaze::columnMajor> mat1( 2UL, 3UL, 2 );
      blaze::UniformMatrix<int,blaze::columnMajor> mat2( mat1 );

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
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major move constructor
   //=====================================================================================

   {
      test_ = "Column-major UniformMatrix move constructor (0x0)";

      blaze::UniformMatrix<int,blaze::columnMajor> mat1( 0UL, 0UL );
      blaze::UniformMatrix<int,blaze::columnMajor> mat2( std::move( mat1 ) );

      checkRows    ( mat2, 0UL );
      checkColumns ( mat2, 0UL );
      checkNonZeros( mat2, 0UL );
   }

   {
      test_ = "Column-major UniformMatrix move constructor (0x3)";

      blaze::UniformMatrix<int,blaze::columnMajor> mat1( 0UL, 3UL );
      blaze::UniformMatrix<int,blaze::columnMajor> mat2( std::move( mat1 ) );

      checkRows    ( mat2, 0UL );
      checkColumns ( mat2, 3UL );
      checkNonZeros( mat2, 0UL );
   }

   {
      test_ = "Column-major UniformMatrix move constructor (2x0)";

      blaze::UniformMatrix<int,blaze::columnMajor> mat1( 2UL, 0UL );
      blaze::UniformMatrix<int,blaze::columnMajor> mat2( std::move( mat1 ) );

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 0UL );
      checkNonZeros( mat2, 0UL );
   }

   {
      test_ = "Column-major UniformMatrix move constructor (2x3)";

      blaze::UniformMatrix<int,blaze::columnMajor> mat1( 2UL, 3UL, 2 );
      blaze::UniformMatrix<int,blaze::columnMajor> mat2( std::move( mat1 ) );

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
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix constructor
   //=====================================================================================

   {
      test_ = "Column-major/row-major UniformMatrix dense matrix constructor (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::rowMajor;

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,rowMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 32UL ) );
      AlignedPadded mat1( memory.get(), 2UL, 3UL, 16UL );
      mat1(0,0) = 2;
      mat1(0,1) = 2;
      mat1(0,2) = 2;
      mat1(1,0) = 2;
      mat1(1,1) = 2;
      mat1(1,2) = 2;

      blaze::UniformMatrix<int,blaze::columnMajor> mat2( mat1 );

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
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major UniformMatrix dense matrix constructor (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowMajor;

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>;
      std::unique_ptr<int[]> memory( new int[7UL] );
      UnalignedUnpadded mat1( memory.get()+1UL, 2UL, 3UL );
      mat1(0,0) = 2;
      mat1(0,1) = 2;
      mat1(0,2) = 2;
      mat1(1,0) = 2;
      mat1(1,1) = 2;
      mat1(1,2) = 2;

      blaze::UniformMatrix<int,blaze::columnMajor> mat2( mat1 );

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
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major UniformMatrix dense matrix constructor (non-uniform)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat1{ { 2, 2, 2 }, { 2, 0, 2 } };

      try {
         blaze::UniformMatrix<int,blaze::columnMajor> mat2( mat1 );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of non-uniform UniformMatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   {
      test_ = "Column-major/column-major UniformMatrix dense matrix constructor (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::columnMajor;

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,columnMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      AlignedPadded mat1( memory.get(), 2UL, 3UL, 16UL );
      mat1(0,0) = 2;
      mat1(0,1) = 2;
      mat1(0,2) = 2;
      mat1(1,0) = 2;
      mat1(1,1) = 2;
      mat1(1,2) = 2;

      blaze::UniformMatrix<int,blaze::columnMajor> mat2( mat1 );

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
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major UniformMatrix dense matrix constructor (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::columnMajor;

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>;
      std::unique_ptr<int[]> memory( new int[7UL] );
      UnalignedUnpadded mat1( memory.get()+1UL, 2UL, 3UL );
      mat1(0,0) = 2;
      mat1(0,1) = 2;
      mat1(0,2) = 2;
      mat1(1,0) = 2;
      mat1(1,1) = 2;
      mat1(1,2) = 2;

      blaze::UniformMatrix<int,blaze::columnMajor> mat2( mat1 );

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
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major UniformMatrix dense matrix constructor (non-uniform)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat1{ { 2, 2, 2 }, { 2, 0, 2 } };

      try {
         blaze::UniformMatrix<int,blaze::columnMajor> mat2( mat1 );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of non-uniform UniformMatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Column-major sparse matrix constructor
   //=====================================================================================

   {
      test_ = "Column-major/row-major UniformMatrix sparse matrix constructor (uniform)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1{ { 2, 2, 2 }, { 2, 2, 2 } };
      blaze::UniformMatrix<int,blaze::columnMajor> mat2( mat1 );

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
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major UniformMatrix sparse matrix constructor (non-uniform)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1{ { 2, 2, 2 }, { 2, 0, 2 } };

      try {
         blaze::UniformMatrix<int,blaze::columnMajor> mat2( mat1 );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of non-uniform UniformMatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   {
      test_ = "Column-major/column-major UniformMatrix sparse matrix constructor (uniform)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1{ { 2, 2, 2 }, { 2, 2, 2 } };
      blaze::UniformMatrix<int,blaze::columnMajor> mat2( mat1 );

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
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major UniformMatrix sparse matrix constructor (non-uniform)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1{ { 2, 2, 2 }, { 2, 0, 2 } };

      try {
         blaze::UniformMatrix<int,blaze::columnMajor> mat2( mat1 );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of non-uniform UniformMatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the UniformMatrix assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all assignment operators of the UniformMatrix class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testAssignment()
{
   //=====================================================================================
   // Row-major homogeneous assignment
   //=====================================================================================

   {
      test_ = "Row-major UniformMatrix homogeneous assignment";

      blaze::UniformMatrix<int,blaze::rowMajor> mat( 3UL, 4UL );
      mat = 2;

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  4UL );
      checkCapacity( mat, 12UL );
      checkNonZeros( mat, 12UL );
      checkNonZeros( mat,  0UL, 4UL );
      checkNonZeros( mat,  1UL, 4UL );
      checkNonZeros( mat,  2UL, 4UL );

      if( mat(0,0) != 2 || mat(0,1) != 2 || mat(0,2) != 2 || mat(0,3) != 2 ||
          mat(1,0) != 2 || mat(1,1) != 2 || mat(1,2) != 2 || mat(1,3) != 2 ||
          mat(2,0) != 2 || mat(2,1) != 2 || mat(2,2) != 2 || mat(2,3) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 2 2 2 2 )\n( 2 2 2 2 )\n( 2 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major copy assignment
   //=====================================================================================

   {
      test_ = "Row-major UniformMatrix copy assignment";

      blaze::UniformMatrix<int,blaze::rowMajor> mat1( 2UL, 3UL, 2 );
      blaze::UniformMatrix<int,blaze::rowMajor> mat2;
      mat2 = mat1;

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
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major UniformMatrix copy assignment stress test";

      using RandomMatrixType = blaze::UniformMatrix<int,blaze::rowMajor>;

      blaze::UniformMatrix<int,blaze::rowMajor> mat1;
      const int min( randmin );
      const int max( randmax );

      for( size_t i=0UL; i<100UL; ++i )
      {
         const size_t rows   ( blaze::rand<size_t>( 0UL, 10UL ) );
         const size_t columns( blaze::rand<size_t>( 0UL, 10UL ) );
         const RandomMatrixType mat2( blaze::rand<RandomMatrixType>( rows, columns, min, max ) );

         mat1 = mat2;

         if( mat1 != mat2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << mat1 << "\n"
                << "   Expected result:\n" << mat2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Row-major move assignment
   //=====================================================================================

   {
      test_ = "Row-major UniformMatrix move assignment";

      blaze::UniformMatrix<int,blaze::rowMajor> mat1( 2UL, 3UL,  2 );
      blaze::UniformMatrix<int,blaze::rowMajor> mat2( 4UL, 1UL, 11 );

      mat2 = std::move( mat1 );

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
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major dense matrix assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major UniformMatrix dense matrix assignment (mixed type)";

      blaze::UniformMatrix<short,blaze::rowMajor> mat1( 2UL, 3UL, 2 );
      blaze::UniformMatrix<int,blaze::rowMajor> mat2;
      mat2 = mat1;

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
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major UniformMatrix dense matrix assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::rowMajor;

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,rowMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 32UL ) );
      AlignedPadded mat1( memory.get(), 2UL, 3UL, 16UL );
      mat1(0,0) = 2;
      mat1(0,1) = 2;
      mat1(0,2) = 2;
      mat1(1,0) = 2;
      mat1(1,1) = 2;
      mat1(1,2) = 2;

      blaze::UniformMatrix<int,blaze::rowMajor> mat2;
      mat2 = mat1;

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
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major UniformMatrix dense matrix assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowMajor;

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>;
      std::unique_ptr<int[]> memory( new int[7UL] );
      UnalignedUnpadded mat1( memory.get()+1UL, 2UL, 3UL );
      mat1(0,0) = 2;
      mat1(0,1) = 2;
      mat1(0,2) = 2;
      mat1(1,0) = 2;
      mat1(1,1) = 2;
      mat1(1,2) = 2;

      blaze::UniformMatrix<int,blaze::rowMajor> mat2;
      mat2 = mat1;

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
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major UniformMatrix dense matrix assignment (non-uniform)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat1{ { 2, 2, 2 }, { 2, 0, 2 } };

      try {
         blaze::UniformMatrix<int,blaze::rowMajor> mat2;
         mat2 = mat1;

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
      test_ = "Row-major/column-major UniformMatrix dense matrix assignment (mixed type)";

      blaze::UniformMatrix<short,blaze::columnMajor> mat1( 2UL, 3UL, 2 );
      blaze::UniformMatrix<int,blaze::rowMajor> mat2;
      mat2 = mat1;

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
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major UniformMatrix dense matrix assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::columnMajor;

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,columnMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      AlignedPadded mat1( memory.get(), 2UL, 3UL, 16UL );
      mat1(0,0) = 2;
      mat1(0,1) = 2;
      mat1(0,2) = 2;
      mat1(1,0) = 2;
      mat1(1,1) = 2;
      mat1(1,2) = 2;

      blaze::UniformMatrix<int,blaze::rowMajor> mat2;
      mat2 = mat1;

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
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major UniformMatrix dense matrix assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::columnMajor;

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>;
      std::unique_ptr<int[]> memory( new int[7UL] );
      UnalignedUnpadded mat1( memory.get()+1UL, 2UL, 3UL );
      mat1(0,0) = 2;
      mat1(0,1) = 2;
      mat1(0,2) = 2;
      mat1(1,0) = 2;
      mat1(1,1) = 2;
      mat1(1,2) = 2;

      blaze::UniformMatrix<int,blaze::rowMajor> mat2;
      mat2 = mat1;

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
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major UniformMatrix dense matrix assignment (non-uniform)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat1{ { 2, 2, 2 }, { 2, 0, 2 } };

      try {
         blaze::UniformMatrix<int,blaze::rowMajor> mat2;
         mat2 = mat1;

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
   // Row-major sparse matrix assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major UniformMatrix sparse matrix assignment (uniform)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1{ { 2, 2, 2 }, { 2, 2, 2 } };
      blaze::UniformMatrix<int,blaze::rowMajor> mat2;
      mat2 = mat1;

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
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major UniformMatrix sparse matrix assignment (non-uniform)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1{ { 2, 2, 2 }, { 2, 0, 2 } };

      try {
         blaze::UniformMatrix<int,blaze::rowMajor> mat2;
         mat2 = mat1;

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
      test_ = "Row-major/column-major UniformMatrix sparse matrix assignment (uniform)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1{ { 2, 2, 2 }, { 2, 2, 2 } };
      blaze::UniformMatrix<int,blaze::rowMajor> mat2;
      mat2 = mat1;

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
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major UniformMatrix sparse matrix assignment (non-uniform)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1{ { 2, 2, 2 }, { 2, 0, 2 } };

      try {
         blaze::UniformMatrix<int,blaze::rowMajor> mat2;
         mat2 = mat1;

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
   // Column-major homogeneous assignment
   //=====================================================================================

   {
      test_ = "Column-major UniformMatrix homogeneous assigment";

      blaze::UniformMatrix<int,blaze::columnMajor> mat( 3UL, 4UL );
      mat = 2;

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  4UL );
      checkCapacity( mat, 12UL );
      checkNonZeros( mat, 12UL );
      checkNonZeros( mat,  0UL, 3UL );
      checkNonZeros( mat,  1UL, 3UL );
      checkNonZeros( mat,  2UL, 3UL );
      checkNonZeros( mat,  3UL, 3UL );

      if( mat(0,0) != 2 || mat(0,1) != 2 || mat(0,2) != 2 || mat(0,3) != 2 ||
          mat(1,0) != 2 || mat(1,1) != 2 || mat(1,2) != 2 || mat(1,3) != 2 ||
          mat(2,0) != 2 || mat(2,1) != 2 || mat(2,2) != 2 || mat(2,3) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 2 2 2 2 )\n( 2 2 2 2 )\n( 2 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major copy assignment
   //=====================================================================================

   {
      test_ = "Column-major UniformMatrix copy assignment";

      blaze::UniformMatrix<int,blaze::columnMajor> mat1( 2UL, 3UL, 2 );
      blaze::UniformMatrix<int,blaze::columnMajor> mat2;
      mat2 = mat1;

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
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major UniformMatrix copy assignment stress test";

      using RandomMatrixType = blaze::UniformMatrix<int,blaze::columnMajor>;

      blaze::UniformMatrix<int,blaze::columnMajor> mat1;
      const int min( randmin );
      const int max( randmax );

      for( size_t i=0UL; i<100UL; ++i )
      {
         const size_t rows   ( blaze::rand<size_t>( 0UL, 10UL ) );
         const size_t columns( blaze::rand<size_t>( 0UL, 10UL ) );
         const RandomMatrixType mat2( blaze::rand<RandomMatrixType>( rows, columns, min, max ) );

         mat1 = mat2;

         if( mat1 != mat2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << mat1 << "\n"
                << "   Expected result:\n" << mat2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major move assignment
   //=====================================================================================

   {
      test_ = "Column-major UniformMatrix move assignment";

      blaze::UniformMatrix<int,blaze::columnMajor> mat1( 2UL, 3UL,  2 );
      blaze::UniformMatrix<int,blaze::columnMajor> mat2( 4UL, 1UL, 11 );

      mat2 = std::move( mat1 );

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
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major UniformMatrix dense matrix assignment (mixed type)";

      blaze::UniformMatrix<short,blaze::rowMajor> mat1( 2UL, 3UL, 2 );
      blaze::UniformMatrix<int,blaze::columnMajor> mat2;
      mat2 = mat1;

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
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major UniformMatrix dense matrix assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::rowMajor;

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,rowMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 32UL ) );
      AlignedPadded mat1( memory.get(), 2UL, 3UL, 16UL );
      mat1(0,0) = 2;
      mat1(0,1) = 2;
      mat1(0,2) = 2;
      mat1(1,0) = 2;
      mat1(1,1) = 2;
      mat1(1,2) = 2;

      blaze::UniformMatrix<int,blaze::columnMajor> mat2;
      mat2 = mat1;

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
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major UniformMatrix dense matrix assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowMajor;

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>;
      std::unique_ptr<int[]> memory( new int[7UL] );
      UnalignedUnpadded mat1( memory.get()+1UL, 2UL, 3UL );
      mat1(0,0) = 2;
      mat1(0,1) = 2;
      mat1(0,2) = 2;
      mat1(1,0) = 2;
      mat1(1,1) = 2;
      mat1(1,2) = 2;

      blaze::UniformMatrix<int,blaze::columnMajor> mat2;
      mat2 = mat1;

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
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major UniformMatrix dense matrix assignment (non-uniform)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat1{ { 2, 2, 2 }, { 2, 0, 2 } };

      try {
         blaze::UniformMatrix<int,blaze::columnMajor> mat2;
         mat2 = mat1;

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
      test_ = "Column-major/column-major UniformMatrix dense matrix assignment (mixed type)";

      blaze::UniformMatrix<short,blaze::columnMajor> mat1( 2UL, 3UL, 2 );
      blaze::UniformMatrix<int,blaze::columnMajor> mat2;
      mat2 = mat1;

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
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major UniformMatrix dense matrix assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::columnMajor;

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,columnMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      AlignedPadded mat1( memory.get(), 2UL, 3UL, 16UL );
      mat1(0,0) = 2;
      mat1(0,1) = 2;
      mat1(0,2) = 2;
      mat1(1,0) = 2;
      mat1(1,1) = 2;
      mat1(1,2) = 2;

      blaze::UniformMatrix<int,blaze::columnMajor> mat2;
      mat2 = mat1;

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
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major UniformMatrix dense matrix assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::columnMajor;

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>;
      std::unique_ptr<int[]> memory( new int[7UL] );
      UnalignedUnpadded mat1( memory.get()+1UL, 2UL, 3UL );
      mat1(0,0) = 2;
      mat1(0,1) = 2;
      mat1(0,2) = 2;
      mat1(1,0) = 2;
      mat1(1,1) = 2;
      mat1(1,2) = 2;

      blaze::UniformMatrix<int,blaze::columnMajor> mat2;
      mat2 = mat1;

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
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major UniformMatrix dense matrix assignment (non-uniform)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat1{ { 2, 2, 2 }, { 2, 0, 2 } };

      try {
         blaze::UniformMatrix<int,blaze::columnMajor> mat2;
         mat2 = mat1;

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
   // Column-major sparse matrix assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major UniformMatrix sparse matrix assignment (uniform)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1{ { 2, 2, 2 }, { 2, 2, 2 } };
      blaze::UniformMatrix<int,blaze::columnMajor> mat2;
      mat2 = mat1;

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
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major UniformMatrix sparse matrix assignment (non-uniform)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1{ { 2, 2, 2 }, { 2, 0, 2 } };

      try {
         blaze::UniformMatrix<int,blaze::columnMajor> mat2;
         mat2 = mat1;

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
      test_ = "Column-major/column-major UniformMatrix sparse matrix assignment (uniform)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1{ { 2, 2, 2 }, { 2, 2, 2 } };
      blaze::UniformMatrix<int,blaze::columnMajor> mat2;
      mat2 = mat1;

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
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 2 2 2 )\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major UniformMatrix sparse matrix assignment (non-uniform)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1{ { 2, 2, 2 }, { 2, 0, 2 } };

      try {
         blaze::UniformMatrix<int,blaze::columnMajor> mat2;
         mat2 = mat1;

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
/*!\brief Test of the UniformMatrix addition assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the addition assignment operators of the UniformMatrix class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testAddAssign()
{
   //=====================================================================================
   // Row-major dense matrix addition assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major UniformMatrix dense matrix addition assignment (mixed type)";

      blaze::UniformMatrix<short,blaze::rowMajor> mat1( 2UL, 3UL, 2 );

      blaze::UniformMatrix<int,blaze::rowMajor> mat2( 2UL, 3UL, 1 );

      mat2 += mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 3UL );
      checkNonZeros( mat2, 1UL, 3UL );

      if( mat2(0,0) != 3 || mat2(0,1) != 3 || mat2(0,2) != 3 ||
          mat2(1,0) != 3 || mat2(1,1) != 3 || mat2(1,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 3 3 3 )\n( 3 3 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major UniformMatrix dense matrix addition assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::rowMajor;

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,rowMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 32UL ) );
      AlignedPadded mat1( memory.get(), 2UL, 3UL, 16UL );
      mat1 = 2;

      blaze::UniformMatrix<int,blaze::rowMajor> mat2( 2UL, 3UL, 1 );

      mat2 += mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 3UL );
      checkNonZeros( mat2, 1UL, 3UL );

      if( mat2(0,0) != 3 || mat2(0,1) != 3 || mat2(0,2) != 3 ||
          mat2(1,0) != 3 || mat2(1,1) != 3 || mat2(1,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 3 3 3 )\n( 3 3 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major UniformMatrix dense matrix addition assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowMajor;

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>;
      std::unique_ptr<int[]> memory( new int[7UL] );
      UnalignedUnpadded mat1( memory.get()+1UL, 2UL, 3UL );
      mat1 = 2;

      blaze::UniformMatrix<int,blaze::rowMajor> mat2( 2UL, 3UL, 1 );

      mat2 += mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 3UL );
      checkNonZeros( mat2, 1UL, 3UL );

      if( mat2(0,0) != 3 || mat2(0,1) != 3 || mat2(0,2) != 3 ||
          mat2(1,0) != 3 || mat2(1,1) != 3 || mat2(1,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 3 3 3 )\n( 3 3 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major UniformMatrix dense matrix addition assignment (non-uniform)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat1{ { 2, 2, 2 }, { 2, 0, 2 } };

      try {
         blaze::UniformMatrix<int,blaze::rowMajor> mat2( 2UL, 3UL, 1 );
         mat2 += mat1;

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
      test_ = "Row-major/column-major UniformMatrix dense matrix addition assignment (mixed type)";

      blaze::UniformMatrix<short,blaze::columnMajor> mat1( 2UL, 3UL, 2 );

      blaze::UniformMatrix<int,blaze::rowMajor> mat2( 2UL, 3UL, 1 );

      mat2 += mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 3UL );
      checkNonZeros( mat2, 1UL, 3UL );

      if( mat2(0,0) != 3 || mat2(0,1) != 3 || mat2(0,2) != 3 ||
          mat2(1,0) != 3 || mat2(1,1) != 3 || mat2(1,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 3 3 3 )\n( 3 3 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major UniformMatrix dense matrix addition assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::columnMajor;

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,columnMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      AlignedPadded mat1( memory.get(), 2UL, 3UL, 16UL );
      mat1 = 2;

      blaze::UniformMatrix<int,blaze::rowMajor> mat2( 2UL, 3UL, 1 );

      mat2 += mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 3UL );
      checkNonZeros( mat2, 1UL, 3UL );

      if( mat2(0,0) != 3 || mat2(0,1) != 3 || mat2(0,2) != 3 ||
          mat2(1,0) != 3 || mat2(1,1) != 3 || mat2(1,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 3 3 3 )\n( 3 3 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major UniformMatrix dense matrix addition assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::columnMajor;

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>;
      std::unique_ptr<int[]> memory( new int[7UL] );
      UnalignedUnpadded mat1( memory.get()+1UL, 2UL, 3UL );
      mat1 = 2;

      blaze::UniformMatrix<int,blaze::rowMajor> mat2( 2UL, 3UL, 1 );

      mat2 += mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 3UL );
      checkNonZeros( mat2, 1UL, 3UL );

      if( mat2(0,0) != 3 || mat2(0,1) != 3 || mat2(0,2) != 3 ||
          mat2(1,0) != 3 || mat2(1,1) != 3 || mat2(1,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 3 3 3 )\n( 3 3 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major UniformMatrix dense matrix addition assignment (non-uniform)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat1{ { 2, 2, 2 }, { 2, 0, 2 } };

      try {
         blaze::UniformMatrix<int,blaze::rowMajor> mat2( 2UL, 3UL, 1 );
         mat2 += mat1;

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
   // Row-major sparse matrix addition assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major UniformMatrix sparse matrix addition assignment";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1{ { 2, 2, 2 }, { 2, 2, 2 } };

      blaze::UniformMatrix<int,blaze::rowMajor> mat2( 2UL, 3UL, 1 );

      mat2 += mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 3UL );
      checkNonZeros( mat2, 1UL, 3UL );

      if( mat2(0,0) != 3 || mat2(0,1) != 3 || mat2(0,2) != 3 ||
          mat2(1,0) != 3 || mat2(1,1) != 3 || mat2(1,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 3 3 3 )\n( 3 3 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major UniformMatrix sparse matrix addition assignment (non-uniform)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1{ { 2, 2, 2 }, { 2, 0, 2 } };

      try {
         blaze::UniformMatrix<int,blaze::rowMajor> mat2( 2UL, 3UL, 1 );
         mat2 += mat1;

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
      test_ = "Row-major/column-major UniformMatrix sparse matrix addition assignment";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1{ { 2, 2, 2 }, { 2, 2, 2 } };

      blaze::UniformMatrix<int,blaze::rowMajor> mat2( 2UL, 3UL, 1 );

      mat2 += mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 3UL );
      checkNonZeros( mat2, 1UL, 3UL );

      if( mat2(0,0) != 3 || mat2(0,1) != 3 || mat2(0,2) != 3 ||
          mat2(1,0) != 3 || mat2(1,1) != 3 || mat2(1,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 3 3 3 )\n( 3 3 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major UniformMatrix sparse matrix addition assignment (non-uniform)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1{ { 2, 2, 2 }, { 2, 0, 2 } };

      try {
         blaze::UniformMatrix<int,blaze::rowMajor> mat2( 2UL, 3UL, 1 );
         mat2 += mat1;

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
   // Column-major dense matrix addition assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major UniformMatrix dense matrix addition assignment (mixed type)";

      blaze::UniformMatrix<short,blaze::rowMajor> mat1( 2UL, 3UL, 2 );

      blaze::UniformMatrix<int,blaze::columnMajor> mat2( 2UL, 3UL, 1 );

      mat2 += mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 3 || mat2(0,1) != 3 || mat2(0,2) != 3 ||
          mat2(1,0) != 3 || mat2(1,1) != 3 || mat2(1,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 3 3 3 )\n( 3 3 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major UniformMatrix dense matrix addition assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::rowMajor;

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,rowMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 32UL ) );
      AlignedPadded mat1( memory.get(), 2UL, 3UL, 16UL );
      mat1 = 2;

      blaze::UniformMatrix<int,blaze::columnMajor> mat2( 2UL, 3UL, 1 );

      mat2 += mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 3 || mat2(0,1) != 3 || mat2(0,2) != 3 ||
          mat2(1,0) != 3 || mat2(1,1) != 3 || mat2(1,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 3 3 3 )\n( 3 3 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major UniformMatrix dense matrix addition assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowMajor;

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>;
      std::unique_ptr<int[]> memory( new int[7UL] );
      UnalignedUnpadded mat1( memory.get()+1UL, 2UL, 3UL );
      mat1 = 2;

      blaze::UniformMatrix<int,blaze::columnMajor> mat2( 2UL, 3UL, 1 );

      mat2 += mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 3 || mat2(0,1) != 3 || mat2(0,2) != 3 ||
          mat2(1,0) != 3 || mat2(1,1) != 3 || mat2(1,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 3 3 3 )\n( 3 3 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major UniformMatrix dense matrix addition assignment (non-uniform)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat1{ { 2, 2, 2 }, { 2, 0, 2 } };

      try {
         blaze::UniformMatrix<int,blaze::columnMajor> mat2( 2UL, 3UL, 1 );
         mat2 += mat1;

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
      test_ = "Column-major/column-major UniformMatrix dense matrix addition assignment (mixed type)";

      blaze::UniformMatrix<short,blaze::columnMajor> mat1( 2UL, 3UL, 2 );

      blaze::UniformMatrix<int,blaze::columnMajor> mat2( 2UL, 3UL, 1 );

      mat2 += mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 3 || mat2(0,1) != 3 || mat2(0,2) != 3 ||
          mat2(1,0) != 3 || mat2(1,1) != 3 || mat2(1,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 3 3 3 )\n( 3 3 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major UniformMatrix dense matrix addition assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::columnMajor;

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,columnMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      AlignedPadded mat1( memory.get(), 2UL, 3UL, 16UL );
      mat1 = 2;

      blaze::UniformMatrix<int,blaze::columnMajor> mat2( 2UL, 3UL, 1 );

      mat2 += mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 3 || mat2(0,1) != 3 || mat2(0,2) != 3 ||
          mat2(1,0) != 3 || mat2(1,1) != 3 || mat2(1,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 3 3 3 )\n( 3 3 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major UniformMatrix dense matrix addition assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::columnMajor;

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>;
      std::unique_ptr<int[]> memory( new int[7UL] );
      UnalignedUnpadded mat1( memory.get()+1UL, 2UL, 3UL );
      mat1 = 2;

      blaze::UniformMatrix<int,blaze::columnMajor> mat2( 2UL, 3UL, 1 );

      mat2 += mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 3 || mat2(0,1) != 3 || mat2(0,2) != 3 ||
          mat2(1,0) != 3 || mat2(1,1) != 3 || mat2(1,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 3 3 3 )\n( 3 3 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major UniformMatrix dense matrix addition assignment (non-uniform)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat1{ { 2, 2, 2 }, { 2, 0, 2 } };

      try {
         blaze::UniformMatrix<int,blaze::rowMajor> mat2( 2UL, 3UL, 1 );
         mat2 += mat1;

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
   // Column-major sparse matrix addition assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major UniformMatrix sparse matrix addition assignment";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1{ { 2, 2, 2 }, { 2, 2, 2 } };

      blaze::UniformMatrix<int,blaze::columnMajor> mat2( 2UL, 3UL, 1 );

      mat2 += mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 3 || mat2(0,1) != 3 || mat2(0,2) != 3 ||
          mat2(1,0) != 3 || mat2(1,1) != 3 || mat2(1,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 3 3 3 )\n( 3 3 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major UniformMatrix sparse matrix addition assignment (non-uniform)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1{ { 2, 2, 2 }, { 2, 0, 2 } };

      try {
         blaze::UniformMatrix<int,blaze::columnMajor> mat2( 2UL, 3UL, 1 );
         mat2 += mat1;

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
      test_ = "Column-major/column-major UniformMatrix sparse matrix addition assignment";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1{ { 2, 2, 2 }, { 2, 2, 2 } };

      blaze::UniformMatrix<int,blaze::columnMajor> mat2( 2UL, 3UL, 1 );

      mat2 += mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 3 || mat2(0,1) != 3 || mat2(0,2) != 3 ||
          mat2(1,0) != 3 || mat2(1,1) != 3 || mat2(1,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 3 3 3 )\n( 3 3 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major UniformMatrix sparse matrix addition assignment (non-uniform)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1{ { 2, 2, 2 }, { 2, 0, 2 } };

      try {
         blaze::UniformMatrix<int,blaze::columnMajor> mat2( 2UL, 3UL, 1 );
         mat2 += mat1;

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
/*!\brief Test of the UniformMatrix subtraction assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the subtraction assignment operators of the UniformMatrix
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testSubAssign()
{
   //=====================================================================================
   // Row-major dense matrix subtraction assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major UniformMatrix dense matrix subtraction assignment (mixed type)";

      blaze::UniformMatrix<short,blaze::rowMajor> mat1( 2UL, 3UL, 2 );

      blaze::UniformMatrix<int,blaze::rowMajor> mat2( 2UL, 3UL, 1 );

      mat2 -= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 3UL );
      checkNonZeros( mat2, 1UL, 3UL );

      if( mat2(0,0) != -1 || mat2(0,1) != -1 || mat2(0,2) != -1 ||
          mat2(1,0) != -1 || mat2(1,1) != -1 || mat2(1,2) != -1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( -1 -1 -1 )\n( -1 -1 -1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major UniformMatrix dense matrix subtraction assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::rowMajor;

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,rowMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 32UL ) );
      AlignedPadded mat1( memory.get(), 2UL, 3UL, 16UL );
      mat1 = 2;

      blaze::UniformMatrix<int,blaze::rowMajor> mat2( 2UL, 3UL, 1 );

      mat2 -= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 3UL );
      checkNonZeros( mat2, 1UL, 3UL );

      if( mat2(0,0) != -1 || mat2(0,1) != -1 || mat2(0,2) != -1 ||
          mat2(1,0) != -1 || mat2(1,1) != -1 || mat2(1,2) != -1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( -1 -1 -1 )\n( -1 -1 -1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major UniformMatrix dense matrix subtraction assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowMajor;

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>;
      std::unique_ptr<int[]> memory( new int[7UL] );
      UnalignedUnpadded mat1( memory.get()+1UL, 2UL, 3UL );
      mat1 = 2;

      blaze::UniformMatrix<int,blaze::rowMajor> mat2( 2UL, 3UL, 1 );

      mat2 -= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 3UL );
      checkNonZeros( mat2, 1UL, 3UL );

      if( mat2(0,0) != -1 || mat2(0,1) != -1 || mat2(0,2) != -1 ||
          mat2(1,0) != -1 || mat2(1,1) != -1 || mat2(1,2) != -1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( -1 -1 -1 )\n( -1 -1 -1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major UniformMatrix dense matrix subtraction assignment (non-uniform)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat1{ { 2, 2, 2 }, { 2, 0, 2 } };

      try {
         blaze::UniformMatrix<int,blaze::rowMajor> mat2( 2UL, 3UL, 1 );
         mat2 -= mat1;

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
      test_ = "Row-major/column-major UniformMatrix dense matrix subtraction assignment (mixed type)";

      blaze::UniformMatrix<short,blaze::columnMajor> mat1( 2UL, 3UL, 2 );

      blaze::UniformMatrix<int,blaze::rowMajor> mat2( 2UL, 3UL, 1 );

      mat2 -= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 3UL );
      checkNonZeros( mat2, 1UL, 3UL );

      if( mat2(0,0) != -1 || mat2(0,1) != -1 || mat2(0,2) != -1 ||
          mat2(1,0) != -1 || mat2(1,1) != -1 || mat2(1,2) != -1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( -1 -1 -1 )\n( -1 -1 -1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major UniformMatrix dense matrix subtraction assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::columnMajor;

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,columnMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      AlignedPadded mat1( memory.get(), 2UL, 3UL, 16UL );
      mat1 = 2;

      blaze::UniformMatrix<int,blaze::rowMajor> mat2( 2UL, 3UL, 1 );

      mat2 -= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 3UL );
      checkNonZeros( mat2, 1UL, 3UL );

      if( mat2(0,0) != -1 || mat2(0,1) != -1 || mat2(0,2) != -1 ||
          mat2(1,0) != -1 || mat2(1,1) != -1 || mat2(1,2) != -1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( -1 -1 -1 )\n( -1 -1 -1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major UniformMatrix dense matrix subtraction assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::columnMajor;

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>;
      std::unique_ptr<int[]> memory( new int[7UL] );
      UnalignedUnpadded mat1( memory.get()+1UL, 2UL, 3UL );
      mat1 = 2;

      blaze::UniformMatrix<int,blaze::rowMajor> mat2( 2UL, 3UL, 1 );

      mat2 -= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 3UL );
      checkNonZeros( mat2, 1UL, 3UL );

      if( mat2(0,0) != -1 || mat2(0,1) != -1 || mat2(0,2) != -1 ||
          mat2(1,0) != -1 || mat2(1,1) != -1 || mat2(1,2) != -1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( -1 -1 -1 )\n( -1 -1 -1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major UniformMatrix dense matrix subtraction assignment (non-uniform)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat1{ { 2, 2, 2 }, { 2, 0, 2 } };

      try {
         blaze::UniformMatrix<int,blaze::rowMajor> mat2( 2UL, 3UL, 1 );
         mat2 -= mat1;

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
   // Row-major sparse matrix subtraction assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major UniformMatrix sparse matrix subtraction assignment";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1{ { 2, 2, 2 }, { 2, 2, 2 } };

      blaze::UniformMatrix<int,blaze::rowMajor> mat2( 2UL, 3UL, 1 );

      mat2 -= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 3UL );
      checkNonZeros( mat2, 1UL, 3UL );

      if( mat2(0,0) != -1 || mat2(0,1) != -1 || mat2(0,2) != -1 ||
          mat2(1,0) != -1 || mat2(1,1) != -1 || mat2(1,2) != -1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( -1 -1 -1 )\n( -1 -1 -1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major UniformMatrix sparse matrix subtraction assignment (non-uniform)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1{ { 2, 2, 2 }, { 2, 0, 2 } };

      try {
         blaze::UniformMatrix<int,blaze::rowMajor> mat2( 2UL, 3UL, 1 );
         mat2 -= mat1;

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
      test_ = "Row-major/column-major UniformMatrix sparse matrix subtraction assignment";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1{ { 2, 2, 2 }, { 2, 2, 2 } };

      blaze::UniformMatrix<int,blaze::rowMajor> mat2( 2UL, 3UL, 1 );

      mat2 -= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 3UL );
      checkNonZeros( mat2, 1UL, 3UL );

      if( mat2(0,0) != -1 || mat2(0,1) != -1 || mat2(0,2) != -1 ||
          mat2(1,0) != -1 || mat2(1,1) != -1 || mat2(1,2) != -1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( -1 -1 -1 )\n( -1 -1 -1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major UniformMatrix sparse matrix subtraction assignment (non-uniform)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1{ { 2, 2, 2 }, { 2, 0, 2 } };

      try {
         blaze::UniformMatrix<int,blaze::rowMajor> mat2( 2UL, 3UL, 1 );
         mat2 -= mat1;

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
   // Column-major dense matrix subtraction assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major UniformMatrix dense matrix subtraction assignment (mixed type)";

      blaze::UniformMatrix<short,blaze::rowMajor> mat1( 2UL, 3UL, 2 );

      blaze::UniformMatrix<int,blaze::columnMajor> mat2( 2UL, 3UL, 1 );

      mat2 -= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != -1 || mat2(0,1) != -1 || mat2(0,2) != -1 ||
          mat2(1,0) != -1 || mat2(1,1) != -1 || mat2(1,2) != -1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( -1 -1 -1 )\n( -1 -1 -1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major UniformMatrix dense matrix subtraction assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::rowMajor;

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,rowMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 32UL ) );
      AlignedPadded mat1( memory.get(), 2UL, 3UL, 16UL );
      mat1 = 2;

      blaze::UniformMatrix<int,blaze::columnMajor> mat2( 2UL, 3UL, 1 );

      mat2 -= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != -1 || mat2(0,1) != -1 || mat2(0,2) != -1 ||
          mat2(1,0) != -1 || mat2(1,1) != -1 || mat2(1,2) != -1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( -1 -1 -1 )\n( -1 -1 -1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major UniformMatrix dense matrix subtraction assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowMajor;

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>;
      std::unique_ptr<int[]> memory( new int[7UL] );
      UnalignedUnpadded mat1( memory.get()+1UL, 2UL, 3UL );
      mat1 = 2;

      blaze::UniformMatrix<int,blaze::columnMajor> mat2( 2UL, 3UL, 1 );

      mat2 -= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != -1 || mat2(0,1) != -1 || mat2(0,2) != -1 ||
          mat2(1,0) != -1 || mat2(1,1) != -1 || mat2(1,2) != -1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( -1 -1 -1 )\n( -1 -1 -1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major UniformMatrix dense matrix subtraction assignment (non-uniform)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat1{ { 2, 2, 2 }, { 2, 0, 2 } };

      try {
         blaze::UniformMatrix<int,blaze::columnMajor> mat2( 2UL, 3UL, 1 );
         mat2 -= mat1;

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
      test_ = "Column-major/column-major UniformMatrix dense matrix subtraction assignment (mixed type)";

      blaze::UniformMatrix<short,blaze::columnMajor> mat1( 2UL, 3UL, 2 );

      blaze::UniformMatrix<int,blaze::columnMajor> mat2( 2UL, 3UL, 1 );

      mat2 -= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != -1 || mat2(0,1) != -1 || mat2(0,2) != -1 ||
          mat2(1,0) != -1 || mat2(1,1) != -1 || mat2(1,2) != -1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( -1 -1 -1 )\n( -1 -1 -1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major UniformMatrix dense matrix subtraction assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::columnMajor;

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,columnMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      AlignedPadded mat1( memory.get(), 2UL, 3UL, 16UL );
      mat1 = 2;

      blaze::UniformMatrix<int,blaze::columnMajor> mat2( 2UL, 3UL, 1 );

      mat2 -= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != -1 || mat2(0,1) != -1 || mat2(0,2) != -1 ||
          mat2(1,0) != -1 || mat2(1,1) != -1 || mat2(1,2) != -1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( -1 -1 -1 )\n( -1 -1 -1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major UniformMatrix dense matrix subtraction assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::columnMajor;

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>;
      std::unique_ptr<int[]> memory( new int[7UL] );
      UnalignedUnpadded mat1( memory.get()+1UL, 2UL, 3UL );
      mat1 = 2;

      blaze::UniformMatrix<int,blaze::columnMajor> mat2( 2UL, 3UL, 1 );

      mat2 -= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != -1 || mat2(0,1) != -1 || mat2(0,2) != -1 ||
          mat2(1,0) != -1 || mat2(1,1) != -1 || mat2(1,2) != -1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( -1 -1 -1 )\n( -1 -1 -1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major UniformMatrix dense matrix subtraction assignment (non-uniform)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat1{ { 2, 2, 2 }, { 2, 0, 2 } };

      try {
         blaze::UniformMatrix<int,blaze::rowMajor> mat2( 2UL, 3UL, 1 );
         mat2 -= mat1;

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
   // Column-major sparse matrix subtraction assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major UniformMatrix sparse matrix subtraction assignment";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1{ { 2, 2, 2 }, { 2, 2, 2 } };

      blaze::UniformMatrix<int,blaze::columnMajor> mat2( 2UL, 3UL, 1 );

      mat2 -= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != -1 || mat2(0,1) != -1 || mat2(0,2) != -1 ||
          mat2(1,0) != -1 || mat2(1,1) != -1 || mat2(1,2) != -1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( -1 -1 -1 )\n( -1 -1 -1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major UniformMatrix sparse matrix subtraction assignment (non-uniform)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1{ { 2, 2, 2 }, { 2, 0, 2 } };

      try {
         blaze::UniformMatrix<int,blaze::columnMajor> mat2( 2UL, 3UL, 1 );
         mat2 -= mat1;

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
      test_ = "Column-major/column-major UniformMatrix sparse matrix subtraction assignment";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1{ { 2, 2, 2 }, { 2, 2, 2 } };

      blaze::UniformMatrix<int,blaze::columnMajor> mat2( 2UL, 3UL, 1 );

      mat2 -= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != -1 || mat2(0,1) != -1 || mat2(0,2) != -1 ||
          mat2(1,0) != -1 || mat2(1,1) != -1 || mat2(1,2) != -1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( -1 -1 -1 )\n( -1 -1 -1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major UniformMatrix sparse matrix subtraction assignment (non-uniform)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1{ { 2, 2, 2 }, { 2, 0, 2 } };

      try {
         blaze::UniformMatrix<int,blaze::columnMajor> mat2( 2UL, 3UL, 1 );
         mat2 -= mat1;

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
   std::cout << "   Running UniformMatrix class test (part 1)..." << std::endl;

   try
   {
      RUN_UNIFORMMATRIX_CLASS_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during UniformMatrix class test (part 1):\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
