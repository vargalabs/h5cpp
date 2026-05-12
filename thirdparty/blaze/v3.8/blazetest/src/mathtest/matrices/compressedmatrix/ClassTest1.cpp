//=================================================================================================
/*!
//  \file src/mathtest/matrices/compressedmatrix/ClassTest1.cpp
//  \brief Source file for the CompressedMatrix class test (part 1)
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
#include <blaze/math/DiagonalMatrix.h>
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/LowerMatrix.h>
#include <blaze/math/UpperMatrix.h>
#include <blaze/util/Complex.h>
#include <blaze/util/Random.h>
#include <blazetest/mathtest/matrices/compressedmatrix/ClassTest.h>
#include <blazetest/mathtest/RandomMaximum.h>
#include <blazetest/mathtest/RandomMinimum.h>

#ifdef BLAZE_USE_HPX_THREADS
#  include <hpx/hpx_main.hpp>
#endif


namespace blazetest {

namespace mathtest {

namespace matrices {

namespace compressedmatrix {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the CompressedMatrix class test.
//
// \exception std::runtime_error Operation error detected.
*/
ClassTest::ClassTest()
{
   testConstructors();
   testAssignment();
   testAddAssign();
   testSubAssign();
   testSchurAssign();
   testMultAssign();
   testScaling();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the CompressedMatrix constructors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all constructors of the CompressedMatrix class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testConstructors()
{
   //=====================================================================================
   // Row-major default constructor
   //=====================================================================================

   {
      test_ = "Row-major CompressedMatrix default constructor";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat;

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );
   }


   //=====================================================================================
   // Row-major size constructor
   //=====================================================================================

   {
      test_ = "Row-major CompressedMatrix size constructor (0x0)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 0UL, 0UL );

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );
   }

   {
      test_ = "Row-major CompressedMatrix size constructor (0x4)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 0UL, 4UL );

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 4UL );
      checkNonZeros( mat, 0UL );
   }

   {
      test_ = "Row-major CompressedMatrix size constructor (3x0)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 0UL );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );
   }

   {
      test_ = "Row-major CompressedMatrix size constructor (3x4)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 4UL );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 4UL );
      checkNonZeros( mat, 0UL );
      checkNonZeros( mat, 0UL, 0UL );
      checkNonZeros( mat, 1UL, 0UL );
      checkNonZeros( mat, 2UL, 0UL );
   }


   //=====================================================================================
   // Row-major size/non-zeros constructor
   //=====================================================================================

   {
      test_ = "Row-major CompressedMatrix size/non-zeros constructor (0x0)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 0UL, 0UL, 5UL );

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 0UL );
      checkCapacity( mat, 5UL );
      checkNonZeros( mat, 0UL );
   }

   {
      test_ = "Row-major CompressedMatrix size/non-zeros constructor (0x4)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 0UL, 4UL, 5UL );

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 4UL );
      checkCapacity( mat, 5UL );
      checkNonZeros( mat, 0UL );
   }

   {
      test_ = "Row-major CompressedMatrix size/non-zeros constructor (3x0)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 0UL, 5UL );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 0UL );
      checkCapacity( mat, 5UL );
      checkNonZeros( mat, 0UL );
   }

   {
      test_ = "Row-major CompressedMatrix size/non-zeros constructor (3x4)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 4UL, 5UL );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 4UL );
      checkCapacity( mat, 5UL );
      checkNonZeros( mat, 0UL );
      checkNonZeros( mat, 0UL, 0UL );
      checkNonZeros( mat, 1UL, 0UL );
      checkNonZeros( mat, 2UL, 0UL );
   }


   //=====================================================================================
   // Row-major size/vector constructor
   //=====================================================================================

   {
      test_ = "Row-major CompressedMatrix size/vector constructor (0x0)";

      std::vector<size_t> nonzeros;
      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 0UL, 0UL, nonzeros );

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 0UL );
      checkCapacity( mat, 0UL );
      checkNonZeros( mat, 0UL );
   }

   {
      test_ = "Row-major CompressedMatrix size/vector constructor (0x5)";

      std::vector<size_t> nonzeros;
      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 0UL, 5UL, nonzeros );

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 5UL );
      checkCapacity( mat, 0UL );
      checkNonZeros( mat, 0UL );
   }

   {
      test_ = "Row-major CompressedMatrix size/vector constructor (3x4)";

      std::vector<size_t> nonzeros( 3UL );
      nonzeros[0] = 2UL;
      nonzeros[1] = 1UL;
      nonzeros[2] = 3UL;
      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 4UL, nonzeros );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 4UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 0UL );
      checkNonZeros( mat, 0UL, 0UL );
      checkNonZeros( mat, 1UL, 0UL );
      checkNonZeros( mat, 2UL, 0UL );

      mat.append( 0UL, 0UL, 1 );
      mat.append( 0UL, 1UL, 2 );
      mat.append( 1UL, 0UL, 3 );
      mat.append( 2UL, 0UL, 4 );
      mat.append( 2UL, 1UL, 5 );
      mat.append( 2UL, 2UL, 6 );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 4UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 6UL );
      checkNonZeros( mat, 0UL, 2UL );
      checkNonZeros( mat, 1UL, 1UL );
      checkNonZeros( mat, 2UL, 3UL );
   }


   //=====================================================================================
   // Row-major list initialization
   //=====================================================================================

   {
      test_ = "Row-major CompressedMatrix initializer list constructor (complete list)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat{ { 0, 2, 0 }, { 4, 0, 6 } };

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 3UL );
      checkNonZeros( mat, 3UL );
      checkNonZeros( mat, 0UL, 1UL );
      checkNonZeros( mat, 1UL, 2UL );

      if( mat(0,0) != 0 || mat(0,1) != 2 || mat(0,2) != 0 ||
          mat(1,0) != 4 || mat(1,1) != 0 || mat(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 2 0 )\n( 4 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major CompressedMatrix initializer list constructor (incomplete list)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat{ { 1 }, { 4, 0, 6 } };

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 3UL );
      checkNonZeros( mat, 3UL );
      checkNonZeros( mat, 0UL, 1UL );
      checkNonZeros( mat, 1UL, 2UL );

      if( mat(0,0) != 1 || mat(0,1) != 0 || mat(0,2) != 0 ||
          mat(1,0) != 4 || mat(1,1) != 0 || mat(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 4 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major copy constructor
   //=====================================================================================

   {
      test_ = "Row-major CompressedMatrix copy constructor (0x0)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1( 0UL, 0UL, 3UL );
      blaze::CompressedMatrix<int,blaze::rowMajor> mat2( mat1 );

      checkRows    ( mat2, 0UL );
      checkColumns ( mat2, 0UL );
      checkNonZeros( mat2, 0UL );
   }

   {
      test_ = "Row-major CompressedMatrix copy constructor (0x3)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1( 0UL, 3UL, 3UL );
      blaze::CompressedMatrix<int,blaze::rowMajor> mat2( mat1 );

      checkRows    ( mat2, 0UL );
      checkColumns ( mat2, 3UL );
      checkNonZeros( mat2, 0UL );
   }

   {
      test_ = "Row-major CompressedMatrix copy constructor (2x0)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1( 2UL, 0UL, 3UL );
      blaze::CompressedMatrix<int,blaze::rowMajor> mat2( mat1 );

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 0UL );
      checkNonZeros( mat2, 0UL );
   }

   {
      test_ = "Row-major CompressedMatrix copy constructor (2x3)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1{ { 1, 0, 2 },
                                                         { 0, 3, 0 } };

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2( mat1 );

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 3UL );
      checkNonZeros( mat2, 3UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 1UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 2 ||
          mat2(1,0) != 0 || mat2(1,1) != 3 || mat2(1,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 2 )\n( 0 3 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major move constructor
   //=====================================================================================

   {
      test_ = "Row-major CompressedMatrix move constructor (0x0)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1( 0UL, 0UL, 3UL );
      blaze::CompressedMatrix<int,blaze::rowMajor> mat2( std::move( mat1 ) );

      checkRows    ( mat2, 0UL );
      checkColumns ( mat2, 0UL );
      checkNonZeros( mat2, 0UL );
   }

   {
      test_ = "Row-major CompressedMatrix move constructor (0x3)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1( 0UL, 3UL, 3UL );
      blaze::CompressedMatrix<int,blaze::rowMajor> mat2( std::move( mat1 ) );

      checkRows    ( mat2, 0UL );
      checkColumns ( mat2, 3UL );
      checkNonZeros( mat2, 0UL );
   }

   {
      test_ = "Row-major CompressedMatrix move constructor (2x0)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1( 2UL, 0UL, 3UL );
      blaze::CompressedMatrix<int,blaze::rowMajor> mat2( std::move( mat1 ) );

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 0UL );
      checkNonZeros( mat2, 0UL );
   }

   {
      test_ = "Row-major CompressedMatrix move constructor (2x3)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1{ { 1, 0, 2 },
                                                         { 0, 3, 0 } };

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2( std::move( mat1 ) );

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 3UL );
      checkNonZeros( mat2, 3UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 1UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 2 ||
          mat2(1,0) != 0 || mat2(1,1) != 3 || mat2(1,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 2 )\n( 0 3 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major dense matrix constructor
   //=====================================================================================

   {
      test_ = "Row-major/row-major CompressedMatrix dense matrix constructor";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat1{ { 0, 1 },
                                                      { 0, 2 },
                                                      { 3, 4 } };

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2( mat1 );

      checkRows    ( mat2, 3UL );
      checkColumns ( mat2, 2UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 1UL );
      checkNonZeros( mat2, 1UL, 1UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 0 || mat2(0,1) != 1 ||
          mat2(1,0) != 0 || mat2(1,1) != 2 ||
          mat2(2,0) != 3 || mat2(2,1) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 0 1 )\n( 0 2 )\n( 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CompressedMatrix dense matrix constructor";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat1{ { 0, 1 },
                                                         { 0, 2 },
                                                         { 3, 4 } };

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2( mat1 );

      checkRows    ( mat2, 3UL );
      checkColumns ( mat2, 2UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 1UL );
      checkNonZeros( mat2, 1UL, 1UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 0 || mat2(0,1) != 1 ||
          mat2(1,0) != 0 || mat2(1,1) != 2 ||
          mat2(2,0) != 3 || mat2(2,1) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 0 1 )\n( 0 2 )\n( 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix constructor
   //=====================================================================================

   {
      test_ = "Row-major/row-major CompressedMatrix sparse matrix constructor";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1{ { 0, 0, 1 },
                                                            { 2, 3, 4 } };

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2( trans( mat1 ) );

      checkRows    ( mat2, 3UL );
      checkColumns ( mat2, 2UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 1UL );
      checkNonZeros( mat2, 1UL, 1UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 0 || mat2(0,1) != 2 ||
          mat2(1,0) != 0 || mat2(1,1) != 3 ||
          mat2(2,0) != 1 || mat2(2,1) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 0 2 )\n( 0 3 )\n( 1 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CompressedMatrix sparse matrix constructor";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1{ { 1, 0, 2 },
                                                            { 0, 3, 4 } };

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2( mat1 );

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 2 ||
          mat2(1,0) != 0 || mat2(1,1) != 3 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 2 )\n( 0 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major default constructor
   //=====================================================================================

   {
      test_ = "Column-major CompressedMatrix default constructor";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat;

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );
   }


   //=====================================================================================
   // Column-major size constructor
   //=====================================================================================

   {
      test_ = "Column-major CompressedMatrix size constructor (0x0)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 0UL, 0UL );

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );
   }

   {
      test_ = "Column-major CompressedMatrix size constructor (0x4)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 0UL, 4UL );

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 4UL );
      checkNonZeros( mat, 0UL );
   }

   {
      test_ = "Column-major CompressedMatrix size constructor (3x0)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 0UL );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );
   }

   {
      test_ = "Column-major CompressedMatrix size constructor (3x4)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 4UL );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 4UL );
      checkNonZeros( mat, 0UL );
      checkNonZeros( mat, 0UL, 0UL );
      checkNonZeros( mat, 1UL, 0UL );
      checkNonZeros( mat, 2UL, 0UL );
      checkNonZeros( mat, 3UL, 0UL );
   }


   //=====================================================================================
   // Column-major size/non-zeros constructor
   //=====================================================================================

   {
      test_ = "Column-major CompressedMatrix size/non-zeros constructor (0x0)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 0UL, 0UL, 5UL );

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 0UL );
      checkCapacity( mat, 5UL );
      checkNonZeros( mat, 0UL );
   }

   {
      test_ = "Column-major CompressedMatrix size/non-zeros constructor (0x4)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 0UL, 4UL, 5UL );

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 4UL );
      checkCapacity( mat, 5UL );
      checkNonZeros( mat, 0UL );
   }

   {
      test_ = "Column-major CompressedMatrix size/non-zeros constructor (3x0)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 0UL, 5UL );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 0UL );
      checkCapacity( mat, 5UL );
      checkNonZeros( mat, 0UL );
   }

   {
      test_ = "Column-major CompressedMatrix size/non-zeros constructor (3x4)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 4UL, 5UL );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 4UL );
      checkCapacity( mat, 5UL );
      checkNonZeros( mat, 0UL );
      checkNonZeros( mat, 0UL, 0UL );
      checkNonZeros( mat, 1UL, 0UL );
      checkNonZeros( mat, 2UL, 0UL );
      checkNonZeros( mat, 3UL, 0UL );
   }


   //=====================================================================================
   // Column-major size/vector constructor
   //=====================================================================================

   {
      test_ = "Column-major CompressedMatrix size/vector constructor (0x0)";

      std::vector<size_t> nonzeros;
      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 0UL, 0UL, nonzeros );

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 0UL );
      checkCapacity( mat, 0UL );
      checkNonZeros( mat, 0UL );
   }

   {
      test_ = "Column-major CompressedMatrix size/vector constructor (5x0)";

      std::vector<size_t> nonzeros;
      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 5UL, 0UL, nonzeros );

      checkRows    ( mat, 5UL );
      checkColumns ( mat, 0UL );
      checkCapacity( mat, 0UL );
      checkNonZeros( mat, 0UL );
   }

   {
      test_ = "Column-major CompressedMatrix size/vector constructor (4x3)";

      std::vector<size_t> nonzeros( 3UL );
      nonzeros[0] = 2UL;
      nonzeros[1] = 1UL;
      nonzeros[2] = 3UL;
      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 4UL, 3UL, nonzeros );

      checkRows    ( mat, 4UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 0UL );
      checkNonZeros( mat, 0UL, 0UL );
      checkNonZeros( mat, 1UL, 0UL );
      checkNonZeros( mat, 2UL, 0UL );

      mat.append( 0UL, 0UL, 1 );
      mat.append( 1UL, 0UL, 2 );
      mat.append( 0UL, 1UL, 3 );
      mat.append( 0UL, 2UL, 4 );
      mat.append( 1UL, 2UL, 5 );
      mat.append( 2UL, 2UL, 6 );

      checkRows    ( mat, 4UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 6UL );
      checkNonZeros( mat, 0UL, 2UL );
      checkNonZeros( mat, 1UL, 1UL );
      checkNonZeros( mat, 2UL, 3UL );
   }


   //=====================================================================================
   // Column-major list initialization
   //=====================================================================================

   {
      test_ = "Column-major CompressedMatrix initializer list constructor (complete list)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat{ { 0, 2, 0 }, { 4, 0, 6 } };

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 3UL );
      checkNonZeros( mat, 3UL );
      checkNonZeros( mat, 0UL, 1UL );
      checkNonZeros( mat, 1UL, 1UL );
      checkNonZeros( mat, 2UL, 1UL );

      if( mat(0,0) != 0 || mat(0,1) != 2 || mat(0,2) != 0 ||
          mat(1,0) != 4 || mat(1,1) != 0 || mat(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 2 0 )\n( 4 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major CompressedMatrix initializer list constructor (incomplete list)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat{ { 1 }, { 4, 0, 6 } };

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 3UL );
      checkNonZeros( mat, 3UL );
      checkNonZeros( mat, 0UL, 2UL );
      checkNonZeros( mat, 1UL, 0UL );
      checkNonZeros( mat, 2UL, 1UL );

      if( mat(0,0) != 1 || mat(0,1) != 0 || mat(0,2) != 0 ||
          mat(1,0) != 4 || mat(1,1) != 0 || mat(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 4 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major copy constructor
   //=====================================================================================

   {
      test_ = "Column-major CompressedMatrix copy constructor (0x0)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1( 0UL, 0UL, 3UL );
      blaze::CompressedMatrix<int,blaze::columnMajor> mat2( mat1 );

      checkRows    ( mat2, 0UL );
      checkColumns ( mat2, 0UL );
      checkNonZeros( mat2, 0UL );
   }

   {
      test_ = "Column-major CompressedMatrix copy constructor (0x3)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1( 0UL, 3UL, 3UL );
      blaze::CompressedMatrix<int,blaze::columnMajor> mat2( mat1 );

      checkRows    ( mat2, 0UL );
      checkColumns ( mat2, 3UL );
      checkNonZeros( mat2, 0UL );
   }

   {
      test_ = "Column-major CompressedMatrix copy constructor (2x0)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1( 2UL, 0UL, 3UL );
      blaze::CompressedMatrix<int,blaze::columnMajor> mat2( mat1 );

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 0UL );
      checkNonZeros( mat2, 0UL );
   }

   {
      test_ = "Column-major CompressedMatrix copy constructor (2x3)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1{ { 1, 0, 2 },
                                                            { 0, 3, 0 } };

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2( mat1 );

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 3UL );
      checkNonZeros( mat2, 3UL );
      checkNonZeros( mat2, 0UL, 1UL );
      checkNonZeros( mat2, 1UL, 1UL );
      checkNonZeros( mat2, 2UL, 1UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 2 ||
          mat2(1,0) != 0 || mat2(1,1) != 3 || mat2(1,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 2 )\n( 0 3 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major move constructor
   //=====================================================================================

   {
      test_ = "Column-major CompressedMatrix move constructor (0x0)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1( 0UL, 0UL, 3UL );
      blaze::CompressedMatrix<int,blaze::columnMajor> mat2( std::move( mat1 ) );

      checkRows    ( mat2, 0UL );
      checkColumns ( mat2, 0UL );
      checkNonZeros( mat2, 0UL );
   }

   {
      test_ = "Column-major CompressedMatrix move constructor (0x3)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1( 0UL, 3UL, 3UL );
      blaze::CompressedMatrix<int,blaze::columnMajor> mat2( std::move( mat1 ) );

      checkRows    ( mat2, 0UL );
      checkColumns ( mat2, 3UL );
      checkNonZeros( mat2, 0UL );
   }

   {
      test_ = "Column-major CompressedMatrix move constructor (2x0)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1( 2UL, 0UL, 3UL );
      blaze::CompressedMatrix<int,blaze::columnMajor> mat2( std::move( mat1 ) );

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 0UL );
      checkNonZeros( mat2, 0UL );
   }

   {
      test_ = "Column-major CompressedMatrix move constructor (2x3)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1{ { 1, 0, 2 },
                                                            { 0, 3, 0 } };

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2( std::move( mat1 ) );

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 3UL );
      checkNonZeros( mat2, 3UL );
      checkNonZeros( mat2, 0UL, 1UL );
      checkNonZeros( mat2, 1UL, 1UL );
      checkNonZeros( mat2, 2UL, 1UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 2 ||
          mat2(1,0) != 0 || mat2(1,1) != 3 || mat2(1,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 2 )\n( 0 3 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix constructor
   //=====================================================================================

   {
      test_ = "Column-major/row-major CompressedMatrix dense matrix constructor";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat1{ { 0, 1 },
                                                      { 0, 2 },
                                                      { 3, 4 } };

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2( mat1 );

      checkRows    ( mat2, 3UL );
      checkColumns ( mat2, 2UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 1UL );
      checkNonZeros( mat2, 1UL, 3UL );

      if( mat2(0,0) != 0 || mat2(0,1) != 1 ||
          mat2(1,0) != 0 || mat2(1,1) != 2 ||
          mat2(2,0) != 3 || mat2(2,1) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 0 1 )\n( 0 2 )\n( 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major CompressedMatrix dense matrix constructor";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat1{ { 0, 1 },
                                                         { 0, 2 },
                                                         { 3, 4 } };

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2( mat1 );

      checkRows    ( mat2, 3UL );
      checkColumns ( mat2, 2UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 1UL );
      checkNonZeros( mat2, 1UL, 3UL );

      if( mat2(0,0) != 0 || mat2(0,1) != 1 ||
          mat2(1,0) != 0 || mat2(1,1) != 2 ||
          mat2(2,0) != 3 || mat2(2,1) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 0 1 )\n( 0 2 )\n( 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse matrix constructor
   //=====================================================================================

   {
      test_ = "Column-major/row-major CompressedMatrix sparse matrix constructor";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1{ { 1, 0, 2 },
                                                         { 0, 3, 4 } };

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2( mat1 );

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 1UL );
      checkNonZeros( mat2, 1UL, 1UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 2 ||
          mat2(1,0) != 0 || mat2(1,1) != 3 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 2 )\n( 0 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/Column-major CompressedMatrix sparse matrix constructor";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1{ { 0, 0, 1 },
                                                         { 2, 3, 4 } };

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2( trans( mat1 ) );

      checkRows    ( mat2, 3UL );
      checkColumns ( mat2, 2UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 1UL );
      checkNonZeros( mat2, 1UL, 3UL );

      if( mat2(0,0) != 0 || mat2(0,1) != 2 ||
          mat2(1,0) != 0 || mat2(1,1) != 3 ||
          mat2(2,0) != 1 || mat2(2,1) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 0 2 )\n( 0 3 )\n( 1 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the CompressedMatrix assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all assignment operators of the CompressedMatrix class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testAssignment()
{
   //=====================================================================================
   // Row-major list assignment
   //=====================================================================================

   {
      test_ = "Row-major CompressedMatrix initializer list assignment (complete list)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat;
      mat = { { 0, 2, 0 }, { 4, 0, 6 } };

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 3UL );
      checkNonZeros( mat, 3UL );
      checkNonZeros( mat, 0UL, 1UL );
      checkNonZeros( mat, 1UL, 2UL );

      if( mat(0,0) != 0 || mat(0,1) != 2 || mat(0,2) != 0 ||
          mat(1,0) != 4 || mat(1,1) != 0 || mat(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 2 0 )\n( 4 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major CompressedMatrix initializer list assignment (incomplete list)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat;
      mat = { { 1 }, { 4, 0, 6 } };

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 3UL );
      checkNonZeros( mat, 3UL );
      checkNonZeros( mat, 0UL, 1UL );
      checkNonZeros( mat, 1UL, 2UL );

      if( mat(0,0) != 1 || mat(0,1) != 0 || mat(0,2) != 0 ||
          mat(1,0) != 4 || mat(1,1) != 0 || mat(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 4 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major copy assignment
   //=====================================================================================

   {
      test_ = "Row-major CompressedMatrix copy assignment";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1{ { 1, 0, 2 },
                                                         { 0, 3, 0 } };

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2;
      mat2 = mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 3UL );
      checkNonZeros( mat2, 3UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 1UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 2 ||
          mat2(1,0) != 0 || mat2(1,1) != 3 || mat2(1,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 2 )\n( 0 3 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major CompressedMatrix copy assignment stress test";

      using RandomMatrixType = blaze::CompressedMatrix<int,blaze::rowMajor>;

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1;
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
      test_ = "Row-major CompressedMatrix move assignment";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1{ { 1, 0, 2 },
                                                         { 0, 3, 0 } };

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2{ {}, { 11 }, {}, { 12 } };

      mat2 = std::move( mat1 );

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 3UL );
      checkNonZeros( mat2, 3UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 1UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 2 ||
          mat2(1,0) != 0 || mat2(1,1) != 3 || mat2(1,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 2 )\n( 0 3 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major dense matrix assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major CompressedMatrix dense matrix assignment";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat1{ { 0, 1 },
                                                      { 0, 2 },
                                                      { 3, 4 } };

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2( 2UL, 3UL );
      mat2 = mat1;

      checkRows    ( mat2, 3UL );
      checkColumns ( mat2, 2UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 1UL );
      checkNonZeros( mat2, 1UL, 1UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 0 || mat2(0,1) != 1 ||
          mat2(1,0) != 0 || mat2(1,1) != 2 ||
          mat2(2,0) != 3 || mat2(2,1) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 0 1 )\n( 0 2 )\n( 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major CompressedMatrix dense matrix assignment stress test";

      using RandomMatrixType = blaze::DynamicMatrix<int,blaze::rowMajor>;

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1;
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

   {
      test_ = "Row-major/column-major CompressedMatrix dense matrix assignment";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat1{ { 0, 1 },
                                                         { 0, 2 },
                                                         { 3, 4 } };

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2( 2UL, 3UL );
      mat2 = mat1;

      checkRows    ( mat2, 3UL );
      checkColumns ( mat2, 2UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 1UL );
      checkNonZeros( mat2, 1UL, 1UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 0 || mat2(0,1) != 1 ||
          mat2(1,0) != 0 || mat2(1,1) != 2 ||
          mat2(2,0) != 3 || mat2(2,1) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 0 1 )\n( 0 2 )\n( 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CompressedMatrix dense matrix assignment stress test";

      using RandomMatrixType = blaze::DynamicMatrix<int,blaze::columnMajor>;

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1;
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

   {
      test_ = "Row-major/row-major CompressedMatrix dense matrix assignment (lower)";

      blaze::LowerMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL );
      randomize( mat2 );

      mat2 = mat1;

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

   {
      test_ = "Row-major/column-major CompressedMatrix dense matrix assignment (lower)";

      blaze::LowerMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL );
      randomize( mat2 );

      mat2 = mat1;

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

   {
      test_ = "Row-major/row-major CompressedMatrix dense matrix assignment (upper)";

      blaze::UpperMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL );
      randomize( mat2 );

      mat2 = mat1;

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

   {
      test_ = "Row-major/column-major CompressedMatrix dense matrix assignment (upper)";

      blaze::UpperMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL );
      randomize( mat2 );

      mat2 = mat1;

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

   {
      test_ = "Row-major/row-major CompressedMatrix dense matrix assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL );
      randomize( mat2 );

      mat2 = mat1;

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

   {
      test_ = "Row-major/column-major CompressedMatrix dense matrix assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL );
      randomize( mat2 );

      mat2 = mat1;

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


   //=====================================================================================
   // Row-major sparse matrix assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major CompressedMatrix sparse matrix assignment";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1{ { 0, 0, 1 },
                                                            { 2, 3, 4 } };

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2( 2UL, 3UL );
      mat2 = trans( mat1 );

      checkRows    ( mat2, 3UL );
      checkColumns ( mat2, 2UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 1UL );
      checkNonZeros( mat2, 1UL, 1UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 0 || mat2(0,1) != 2 ||
          mat2(1,0) != 0 || mat2(1,1) != 3 ||
          mat2(2,0) != 1 || mat2(2,1) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 0 2 )\n( 0 3 )\n( 1 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major CompressedMatrix sparse matrix assignment stress test";

      using RandomMatrixType = blaze::CompressedMatrix<short,blaze::rowMajor>;

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1;
      const short min( randmin );
      const short max( randmax );

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

   {
      test_ = "Row-major/column-major CompressedMatrix sparse matrix assignment";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1{ { 1, 0, 2 },
                                                            { 0, 3, 4 } };

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2( 3UL, 2UL );
      mat2 = mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 2 ||
          mat2(1,0) != 0 || mat2(1,1) != 3 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 2 )\n( 0 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CompressedMatrix sparse matrix assignment stress test";

      using RandomMatrixType = blaze::CompressedMatrix<short,blaze::columnMajor>;

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1;
      const short min( randmin );
      const short max( randmax );

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

   {
      test_ = "Row-major/row-major CompressedMatrix sparse matrix assignment (lower)";

      blaze::LowerMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL );
      randomize( mat2 );

      mat2 = mat1;

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

   {
      test_ = "Row-major/column-major CompressedMatrix sparse matrix assignment (lower)";

      blaze::LowerMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL );
      randomize( mat2 );

      mat2 = mat1;

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

   {
      test_ = "Row-major/row-major CompressedMatrix sparse matrix assignment (upper)";

      blaze::UpperMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL );
      randomize( mat2 );

      mat2 = mat1;

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

   {
      test_ = "Row-major/column-major CompressedMatrix sparse matrix assignment (upper)";

      blaze::UpperMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL );
      randomize( mat2 );

      mat2 = mat1;

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

   {
      test_ = "Row-major/row-major CompressedMatrix sparse matrix assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL );
      randomize( mat2 );

      mat2 = mat1;

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

   {
      test_ = "Row-major/column-major CompressedMatrix sparse matrix assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL );
      randomize( mat2 );

      mat2 = mat1;

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


   //=====================================================================================
   // Column-major list assignment
   //=====================================================================================

   {
      test_ = "Column-major CompressedMatrix initializer list assignment (complete list)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat;
      mat = { { 0, 2, 0 }, { 4, 0, 6 } };

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 3UL );
      checkNonZeros( mat, 3UL );
      checkNonZeros( mat, 0UL, 1UL );
      checkNonZeros( mat, 1UL, 1UL );
      checkNonZeros( mat, 2UL, 1UL );

      if( mat(0,0) != 0 || mat(0,1) != 2 || mat(0,2) != 0 ||
          mat(1,0) != 4 || mat(1,1) != 0 || mat(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 2 0 )\n( 4 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major CompressedMatrix initializer list assignment (incomplete list)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat;
      mat = { { 1 }, { 4, 0, 6 } };

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 3UL );
      checkNonZeros( mat, 3UL );
      checkNonZeros( mat, 0UL, 2UL );
      checkNonZeros( mat, 1UL, 0UL );
      checkNonZeros( mat, 2UL, 1UL );

      if( mat(0,0) != 1 || mat(0,1) != 0 || mat(0,2) != 0 ||
          mat(1,0) != 4 || mat(1,1) != 0 || mat(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 4 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major copy assignment
   //=====================================================================================

   {
      test_ = "Column-major CompressedMatrix copy assignment";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1{ { 1, 0, 2 },
                                                            { 0, 3, 0 } };

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2;
      mat2 = mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 3UL );
      checkNonZeros( mat2, 3UL );
      checkNonZeros( mat2, 0UL, 1UL );
      checkNonZeros( mat2, 1UL, 1UL );
      checkNonZeros( mat2, 2UL, 1UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 2 ||
          mat2(1,0) != 0 || mat2(1,1) != 3 || mat2(1,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 2 )\n( 0 3 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major CompressedMatrix copy assignment stress test";

      using RandomMatrixType = blaze::CompressedMatrix<int,blaze::columnMajor>;

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1;
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
      test_ = "Column-major CompressedMatrix move assignment";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1{ { 1, 0, 2 },
                                                            { 0, 3, 0 } };

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2{ {}, { 11 }, {}, { 12 } };

      mat2 = std::move( mat1 );

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 3UL );
      checkNonZeros( mat2, 3UL );
      checkNonZeros( mat2, 0UL, 1UL );
      checkNonZeros( mat2, 1UL, 1UL );
      checkNonZeros( mat2, 2UL, 1UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 2 ||
          mat2(1,0) != 0 || mat2(1,1) != 3 || mat2(1,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 2 )\n( 0 3 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major CompressedMatrix dense matrix assignment";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat1{ { 0, 1 },
                                                      { 0, 2 },
                                                      { 3, 4 } };

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2( 2UL, 3UL );
      mat2 = mat1;

      checkRows    ( mat2, 3UL );
      checkColumns ( mat2, 2UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 1UL );
      checkNonZeros( mat2, 1UL, 3UL );

      if( mat2(0,0) != 0 || mat2(0,1) != 1 ||
          mat2(1,0) != 0 || mat2(1,1) != 2 ||
          mat2(2,0) != 3 || mat2(2,1) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 0 1 )\n( 0 2 )\n( 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major CompressedMatrix dense matrix assignment stress test";

      using RandomMatrixType = blaze::DynamicMatrix<int,blaze::rowMajor>;

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1;
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

   {
      test_ = "Column-major/column-major CompressedMatrix dense matrix assignment";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat1{ { 0, 1 },
                                                         { 0, 2 },
                                                         { 3, 4 } };

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2( 2UL, 3UL );
      mat2 = mat1;

      checkRows    ( mat2, 3UL );
      checkColumns ( mat2, 2UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 1UL );
      checkNonZeros( mat2, 1UL, 3UL );

      if( mat2(0,0) != 0 || mat2(0,1) != 1 ||
          mat2(1,0) != 0 || mat2(1,1) != 2 ||
          mat2(2,0) != 3 || mat2(2,1) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 0 1 )\n( 0 2 )\n( 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major CompressedMatrix dense matrix assignment stress test";

      using RandomMatrixType = blaze::DynamicMatrix<int,blaze::columnMajor>;

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1;
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

   {
      test_ = "Column-major/row-major CompressedMatrix dense matrix assignment (lower)";

      blaze::LowerMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL );
      randomize( mat2 );

      mat2 = mat1;

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

   {
      test_ = "Column-major/column-major CompressedMatrix dense matrix assignment (lower)";

      blaze::LowerMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL );
      randomize( mat2 );

      mat2 = mat1;

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

   {
      test_ = "Column-major/row-major CompressedMatrix dense matrix assignment (upper)";

      blaze::UpperMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL );
      randomize( mat2 );

      mat2 = mat1;

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

   {
      test_ = "Column-major/column-major CompressedMatrix dense matrix assignment (upper)";

      blaze::UpperMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL );
      randomize( mat2 );

      mat2 = mat1;

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

   {
      test_ = "Column-major/row-major CompressedMatrix dense matrix assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL );
      randomize( mat2 );

      mat2 = mat1;

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

   {
      test_ = "Column-major/column-major CompressedMatrix dense matrix assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL );
      randomize( mat2 );

      mat2 = mat1;

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


   //=====================================================================================
   // Column-major sparse matrix assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major CompressedMatrix sparse matrix assignment";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1{ { 1, 0, 2 },
                                                         { 0, 3, 4 } };

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2( 3UL, 2UL );
      mat2 = mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 1UL );
      checkNonZeros( mat2, 1UL, 1UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 2 ||
          mat2(1,0) != 0 || mat2(1,1) != 3 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 2 )\n( 0 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major CompressedMatrix sparse matrix assignment stress test";

      using RandomMatrixType = blaze::CompressedMatrix<short,blaze::rowMajor>;

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1;
      const short min( randmin );
      const short max( randmax );

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

   {
      test_ = "Column-major/Column-major CompressedMatrix sparse matrix assignment";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1{ { 0, 0, 1 },
                                                         { 2, 3, 4 } };

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2( 2UL, 3UL );
      mat2 = trans( mat1 );

      checkRows    ( mat2, 3UL );
      checkColumns ( mat2, 2UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 1UL );
      checkNonZeros( mat2, 1UL, 3UL );

      if( mat2(0,0) != 0 || mat2(0,1) != 2 ||
          mat2(1,0) != 0 || mat2(1,1) != 3 ||
          mat2(2,0) != 1 || mat2(2,1) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 0 2 )\n( 0 3 )\n( 1 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major CompressedMatrix sparse matrix assignment stress test";

      using RandomMatrixType = blaze::CompressedMatrix<short,blaze::columnMajor>;

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1;
      const short min( randmin );
      const short max( randmax );

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

   {
      test_ = "Column-major/row-major CompressedMatrix sparse matrix assignment (lower)";

      blaze::LowerMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL );
      randomize( mat2 );

      mat2 = mat1;

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

   {
      test_ = "Column-major/column-major CompressedMatrix sparse matrix assignment (lower)";

      blaze::LowerMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL );
      randomize( mat2 );

      mat2 = mat1;

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

   {
      test_ = "Column-major/row-major CompressedMatrix sparse matrix assignment (upper)";

      blaze::UpperMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL );
      randomize( mat2 );

      mat2 = mat1;

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

   {
      test_ = "Column-major/column-major CompressedMatrix sparse matrix assignment (upper)";

      blaze::UpperMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL );
      randomize( mat2 );

      mat2 = mat1;

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

   {
      test_ = "Column-major/row-major CompressedMatrix sparse matrix assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL );
      randomize( mat2 );

      mat2 = mat1;

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

   {
      test_ = "Column-major/column-major CompressedMatrix sparse matrix assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL );
      randomize( mat2 );

      mat2 = mat1;

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
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the CompressedMatrix addition assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the addition assignment operators of the CompressedMatrix
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testAddAssign()
{
   //=====================================================================================
   // Row-major dense matrix addition assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major CompressedMatrix dense matrix addition assignment";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat1{ { 0, 0, 1 },
                                                      { 2, 3, 4 } };

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2{ { 5, 0, 6 },
                                                         { 0, 7, 0 } };

      mat2 += mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkNonZeros( mat2, 5UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 3UL );

      if( mat2(0,0) != 5 || mat2(0,1) !=  0 || mat2(0,2) != 7 ||
          mat2(1,0) != 2 || mat2(1,1) != 10 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 5  0 7 )\n( 2 10 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CompressedMatrix dense matrix addition assignment";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat1{ { 0, 0, 1 },
                                                         { 2, 3, 4 } };

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2{ { 5, 0, 6 },
                                                         { 0, 7, 0 } };

      mat2 += mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkNonZeros( mat2, 5UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 3UL );

      if( mat2(0,0) != 5 || mat2(0,1) !=  0 || mat2(0,2) != 7 ||
          mat2(1,0) != 2 || mat2(1,1) != 10 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 5  0 7 )\n( 2 10 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major CompressedMatrix dense matrix addition assignment (lower)";

      blaze::LowerMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL );

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CompressedMatrix dense matrix addition assignment (lower)";

      blaze::LowerMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL );

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major CompressedMatrix dense matrix addition assignment (upper)";

      blaze::UpperMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL );

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CompressedMatrix dense matrix addition assignment (upper)";

      blaze::UpperMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL );

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major CompressedMatrix dense matrix addition assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL );

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CompressedMatrix dense matrix addition assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL );

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix addition assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major CompressedMatrix sparse matrix addition assignment";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1{ { 0, 0, 1 },
                                                         { 2, 3, 4 } };

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2{ { 5, 0, 6 },
                                                         { 0, 7, 0 } };

      mat2 += mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkNonZeros( mat2, 5UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 3UL );

      if( mat2(0,0) != 5 || mat2(0,1) !=  0 || mat2(0,2) != 7 ||
          mat2(1,0) != 2 || mat2(1,1) != 10 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 5  0 7 )\n( 2 10 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CompressedMatrix sparse matrix addition assignment";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1{ { 0, 0, 1 },
                                                            { 2, 3, 4 } };

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2{ { 5, 0, 6 },
                                                         { 0, 7, 0 } };

      mat2 += mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkNonZeros( mat2, 5UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 3UL );

      if( mat2(0,0) != 5 || mat2(0,1) !=  0 || mat2(0,2) != 7 ||
          mat2(1,0) != 2 || mat2(1,1) != 10 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 5  0 7 )\n( 2 10 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major CompressedMatrix sparse matrix addition assignment (lower)";

      blaze::LowerMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL );

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CompressedMatrix sparse matrix addition assignment (lower)";

      blaze::LowerMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL );

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major CompressedMatrix sparse matrix addition assignment (upper)";

      blaze::UpperMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL );

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CompressedMatrix sparse matrix addition assignment (upper)";

      blaze::UpperMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL );

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major CompressedMatrix sparse matrix addition assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL );

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CompressedMatrix sparse matrix addition assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL );

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix addition assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major CompressedMatrix dense matrix addition assignment";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat1{ { 0, 0, 1 },
                                                      { 2, 3, 4 } };

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2{ { 5, 0, 6 },
                                                            { 0, 7, 0 } };

      mat2 += mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkNonZeros( mat2, 5UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 1UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 5 || mat2(0,1) !=  0 || mat2(0,2) != 7 ||
          mat2(1,0) != 2 || mat2(1,1) != 10 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 5  0 7 )\n( 2 10 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CompressedMatrix dense matrix addition assignment";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat1{ { 0, 0, 1 },
                                                         { 2, 3, 4 } };

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2{ { 5, 0, 6 },
                                                            { 0, 7, 0 } };

      mat2 += mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkNonZeros( mat2, 5UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 1UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 5 || mat2(0,1) !=  0 || mat2(0,2) != 7 ||
          mat2(1,0) != 2 || mat2(1,1) != 10 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 5  0 7 )\n( 2 10 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major CompressedMatrix dense matrix addition assignment (lower)";

      blaze::LowerMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL );

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major CompressedMatrix dense matrix addition assignment (lower)";

      blaze::LowerMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL );

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major CompressedMatrix dense matrix addition assignment (upper)";

      blaze::UpperMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL );

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major CompressedMatrix dense matrix addition assignment (upper)";

      blaze::UpperMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL );

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major CompressedMatrix dense matrix addition assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL );

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major CompressedMatrix dense matrix addition assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL );

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse matrix addition assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major CompressedMatrix sparse matrix addition assignment";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1{ { 0, 0, 1 },
                                                         { 2, 3, 4 } };

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2{ { 5, 0, 6 },
                                                            { 0, 7, 0 } };

      mat2 += mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkNonZeros( mat2, 5UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 1UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 5 || mat2(0,1) !=  0 || mat2(0,2) != 7 ||
          mat2(1,0) != 2 || mat2(1,1) != 10 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 5  0 7 )\n( 2 10 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major CompressedMatrix sparse matrix addition assignment";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1{ { 0, 0, 1 },
                                                            { 2, 3, 4 } };

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2{ { 5, 0, 6 },
                                                            { 0, 7, 0 } };

      mat2 += mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkNonZeros( mat2, 5UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 1UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 5 || mat2(0,1) !=  0 || mat2(0,2) != 7 ||
          mat2(1,0) != 2 || mat2(1,1) != 10 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 5  0 7 )\n( 2 10 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major CompressedMatrix sparse matrix addition assignment (lower)";

      blaze::LowerMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL );

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CompressedMatrix sparse matrix addition assignment (lower)";

      blaze::LowerMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL );

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major CompressedMatrix sparse matrix addition assignment (upper)";

      blaze::UpperMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL );

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CompressedMatrix sparse matrix addition assignment (upper)";

      blaze::UpperMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL );

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major CompressedMatrix sparse matrix addition assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL );

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CompressedMatrix sparse matrix addition assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL );

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the CompressedMatrix subtraction assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the subtraction assignment operators of the CompressedMatrix
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testSubAssign()
{
   //=====================================================================================
   // Row-major dense matrix subtraction assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major CompressedMatrix dense matrix subtraction assignment";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat1{ { 0, 0, 1 },
                                                      { 2, 3, 4 } };

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2{ { 5, 0, 6 },
                                                         { 0, 7, 0 } };

      mat2 -= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkNonZeros( mat2, 5UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 3UL );

      if( mat2(0,0) !=  5 || mat2(0,1) != 0 || mat2(0,2) !=  5 ||
          mat2(1,0) != -2 || mat2(1,1) != 4 || mat2(1,2) != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n(  5 0  5 )\n( -2 4 -4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CompressedMatrix dense matrix subtraction assignment";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat1{ { 0, 0, 1 },
                                                         { 2, 3, 4 } };

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2{ { 5, 0, 6 },
                                                         { 0, 7, 0 } };

      mat2 -= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkNonZeros( mat2, 5UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 3UL );

      if( mat2(0,0) !=  5 || mat2(0,1) != 0 || mat2(0,2) !=  5 ||
          mat2(1,0) != -2 || mat2(1,1) != 4 || mat2(1,2) != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n(  5 0  5 )\n( -2 4 -4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major CompressedMatrix dense matrix subtraction assignment (lower)";

      blaze::LowerMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL );

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CompressedMatrix dense matrix subtraction assignment (lower)";

      blaze::LowerMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL );

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major CompressedMatrix dense matrix subtraction assignment (upper)";

      blaze::UpperMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL );

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CompressedMatrix dense matrix subtraction assignment (upper)";

      blaze::UpperMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL );

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major CompressedMatrix dense matrix subtraction assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL );

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CompressedMatrix dense matrix subtraction assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL );

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix subtraction assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major CompressedMatrix sparse matrix subtraction assignment";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1{ { 0, 0, 1 },
                                                         { 2, 3, 4 } };

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2{ { 5, 0, 6 },
                                                         { 0, 7, 0 } };

      mat2 -= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkNonZeros( mat2, 5UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 3UL );

      if( mat2(0,0) !=  5 || mat2(0,1) != 0 || mat2(0,2) !=  5 ||
          mat2(1,0) != -2 || mat2(1,1) != 4 || mat2(1,2) != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n(  5 0  5 )\n( -2 4 -4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CompressedMatrix sparse matrix subtraction assignment";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1{ { 0, 0, 1 },
                                                            { 2, 3, 4 } };

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2{ { 5, 0, 6 },
                                                         { 0, 7, 0 } };

      mat2 -= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkNonZeros( mat2, 5UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 3UL );

      if( mat2(0,0) !=  5 || mat2(0,1) != 0 || mat2(0,2) !=  5 ||
          mat2(1,0) != -2 || mat2(1,1) != 4 || mat2(1,2) != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n(  5 0  5 )\n( -2 4 -4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major CompressedMatrix sparse matrix subtraction assignment (lower)";

      blaze::LowerMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL );

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CompressedMatrix sparse matrix subtraction assignment (lower)";

      blaze::LowerMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL );

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major CompressedMatrix sparse matrix subtraction assignment (upper)";

      blaze::UpperMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL );

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CompressedMatrix sparse matrix subtraction assignment (upper)";

      blaze::UpperMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL );

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major CompressedMatrix sparse matrix subtraction assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL );

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CompressedMatrix sparse matrix subtraction assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL );

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix subtraction assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major CompressedMatrix dense matrix subtraction assignment";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat1{ { 0, 0, 1 },
                                                      { 2, 3, 4 } };

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2{ { 5, 0, 6 },
                                                            { 0, 7, 0 } };

      mat2 -= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkNonZeros( mat2, 5UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 1UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) !=  5 || mat2(0,1) != 0 || mat2(0,2) !=  5 ||
          mat2(1,0) != -2 || mat2(1,1) != 4 || mat2(1,2) != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n(  5 0  5 )\n( -2 4 -4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CompressedMatrix dense matrix subtraction assignment";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat1{ { 0, 0, 1 },
                                                         { 2, 3, 4 } };

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2{ { 5, 0, 6 },
                                                            { 0, 7, 0 } };

      mat2 -= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkNonZeros( mat2, 5UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 1UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) !=  5 || mat2(0,1) != 0 || mat2(0,2) !=  5 ||
          mat2(1,0) != -2 || mat2(1,1) != 4 || mat2(1,2) != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n(  5 0  5 )\n( -2 4 -4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major CompressedMatrix dense matrix subtraction assignment (lower)";

      blaze::LowerMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL );

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major CompressedMatrix dense matrix subtraction assignment (lower)";

      blaze::LowerMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL );

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major CompressedMatrix dense matrix subtraction assignment (upper)";

      blaze::UpperMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL );

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major CompressedMatrix dense matrix subtraction assignment (upper)";

      blaze::UpperMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL );

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major CompressedMatrix dense matrix subtraction assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL );

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major CompressedMatrix dense matrix subtraction assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL );

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse matrix subtraction assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major CompressedMatrix sparse matrix subtraction assignment";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1{ { 0, 0, 1 },
                                                         { 2, 3, 4 } };

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2{ { 5, 0, 6 },
                                                            { 0, 7, 0 } };

      mat2 -= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkNonZeros( mat2, 5UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 1UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) !=  5 || mat2(0,1) != 0 || mat2(0,2) !=  5 ||
          mat2(1,0) != -2 || mat2(1,1) != 4 || mat2(1,2) != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n(  5 0  5 )\n( -2 4 -4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major CompressedMatrix sparse matrix subtraction assignment";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1{ { 0, 0, 1 },
                                                            { 2, 3, 4 } };

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2{ { 5, 0, 6 },
                                                            { 0, 7, 0 } };

      mat2 -= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkNonZeros( mat2, 5UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 1UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) !=  5 || mat2(0,1) != 0 || mat2(0,2) !=  5 ||
          mat2(1,0) != -2 || mat2(1,1) != 4 || mat2(1,2) != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n(  5 0  5 )\n( -2 4 -4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major CompressedMatrix sparse matrix subtraction assignment (lower)";

      blaze::LowerMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL );

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major CompressedMatrix sparse matrix subtraction assignment (lower)";

      blaze::LowerMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL );

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major CompressedMatrix sparse matrix subtraction assignment (upper)";

      blaze::UpperMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL );

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major CompressedMatrix sparse matrix subtraction assignment (upper)";

      blaze::UpperMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL );

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major CompressedMatrix sparse matrix subtraction assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL );

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major CompressedMatrix sparse matrix subtraction assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL );

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the CompressedMatrix Schur product assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the Schur product assignment operators of the CompressedMatrix
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testSchurAssign()
{
   //=====================================================================================
   // Row-major dense matrix Schur product assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major CompressedMatrix dense matrix Schur product assignment";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat1{ { 0, 0, 1 },
                                                      { 2, 3, 4 } };

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2{ { 5, 0, 6 },
                                                         { 0, 7, 0 } };

      mat2 %= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkNonZeros( mat2, 3UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 1UL );

      if( mat2(0,0) != 0 || mat2(0,1) !=  0 || mat2(0,2) != 6 ||
          mat2(1,0) != 0 || mat2(1,1) != 21 || mat2(1,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 0  0 6 )\n( 0 21 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CompressedMatrix dense matrix Schur product assignment";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat1{ { 0, 0, 1 },
                                                         { 2, 3, 4 } };

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2{ { 5, 0, 6 },
                                                         { 0, 7, 0 } };

      mat2 %= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkNonZeros( mat2, 3UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 1UL );

      if( mat2(0,0) != 0 || mat2(0,1) !=  0 || mat2(0,2) != 6 ||
          mat2(1,0) != 0 || mat2(1,1) != 21 || mat2(1,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 0  0 6 )\n( 0 21 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major CompressedMatrix dense matrix Schur product assignment (lower)";

      blaze::LowerMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2( blaze::DynamicMatrix<int>( 3UL, 3UL, 1 ) );

      mat2 %= mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CompressedMatrix dense matrix Schur product assignment (lower)";

      blaze::LowerMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2( blaze::DynamicMatrix<int>( 3UL, 3UL, 1 ) );

      mat2 %= mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major CompressedMatrix dense matrix Schur product assignment (upper)";

      blaze::UpperMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2( blaze::DynamicMatrix<int>( 3UL, 3UL, 1 ) );

      mat2 %= mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CompressedMatrix dense matrix Schur product assignment (upper)";

      blaze::UpperMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2( blaze::DynamicMatrix<int>( 3UL, 3UL, 1 ) );

      mat2 %= mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major CompressedMatrix dense matrix Schur product assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2( blaze::DynamicMatrix<int>( 3UL, 3UL, 1 ) );

      mat2 %= mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CompressedMatrix dense matrix Schur product assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2( blaze::DynamicMatrix<int>( 3UL, 3UL, 1 ) );

      mat2 %= mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix Schur product assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major CompressedMatrix sparse matrix Schur product assignment";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1{ { 0, 0, 1 },
                                                         { 2, 3, 4 } };

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2{ { 5, 0, 6 },
                                                         { 0, 7, 0 } };

      mat2 %= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkNonZeros( mat2, 2UL );
      checkNonZeros( mat2, 0UL, 1UL );
      checkNonZeros( mat2, 1UL, 1UL );

      if( mat2(0,0) != 0 || mat2(0,1) !=  0 || mat2(0,2) != 6 ||
          mat2(1,0) != 0 || mat2(1,1) != 21 || mat2(1,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 0  0 6 )\n( 0 21 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CompressedMatrix sparse matrix Schur product assignment";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1{ { 0, 0, 1 },
                                                            { 2, 3, 4 } };

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2{ { 5, 0, 6 },
                                                         { 0, 7, 0 } };

      mat2 %= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkNonZeros( mat2, 2UL );
      checkNonZeros( mat2, 0UL, 1UL );
      checkNonZeros( mat2, 1UL, 1UL );

      if( mat2(0,0) != 0 || mat2(0,1) !=  0 || mat2(0,2) != 6 ||
          mat2(1,0) != 0 || mat2(1,1) != 21 || mat2(1,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 0  0 6 )\n( 0 21 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major CompressedMatrix sparse matrix Schur product assignment (lower)";

      blaze::LowerMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2( blaze::DynamicMatrix<int>( 3UL, 3UL, 1 ) );

      mat2 %= mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CompressedMatrix sparse matrix Schur product assignment (lower)";

      blaze::LowerMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2( blaze::DynamicMatrix<int>( 3UL, 3UL, 1 ) );

      mat2 %= mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major CompressedMatrix sparse matrix Schur product assignment (upper)";

      blaze::UpperMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2( blaze::DynamicMatrix<int>( 3UL, 3UL, 1 ) );

      mat2 %= mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CompressedMatrix sparse matrix Schur product assignment (upper)";

      blaze::UpperMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2( blaze::DynamicMatrix<int>( 3UL, 3UL, 1 ) );

      mat2 %= mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major CompressedMatrix sparse matrix Schur product assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2( blaze::DynamicMatrix<int>( 3UL, 3UL, 1 ) );

      mat2 %= mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CompressedMatrix sparse matrix Schur product assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2( blaze::DynamicMatrix<int>( 3UL, 3UL, 1 ) );

      mat2 %= mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix Schur product assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major CompressedMatrix dense matrix Schur product assignment";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat1{ { 0, 0, 1 },
                                                      { 2, 3, 4 } };

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2{ { 5, 0, 6 },
                                                            { 0, 7, 0 } };

      mat2 %= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkNonZeros( mat2, 3UL );
      checkNonZeros( mat2, 0UL, 1UL );
      checkNonZeros( mat2, 1UL, 1UL );
      checkNonZeros( mat2, 2UL, 1UL );

      if( mat2(0,0) != 0 || mat2(0,1) !=  0 || mat2(0,2) != 6 ||
          mat2(1,0) != 0 || mat2(1,1) != 21 || mat2(1,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 0  0 6 )\n( 0 21 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CompressedMatrix dense matrix Schur product assignment";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat1{ { 0, 0, 1 },
                                                         { 2, 3, 4 } };

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2{ { 5, 0, 6 },
                                                            { 0, 7, 0 } };

      mat2 %= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkNonZeros( mat2, 3UL );
      checkNonZeros( mat2, 0UL, 1UL );
      checkNonZeros( mat2, 1UL, 1UL );
      checkNonZeros( mat2, 2UL, 1UL );

      if( mat2(0,0) != 0 || mat2(0,1) !=  0 || mat2(0,2) != 6 ||
          mat2(1,0) != 0 || mat2(1,1) != 21 || mat2(1,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 0  0 6 )\n( 0 21 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major CompressedMatrix dense matrix Schur product assignment (lower)";

      blaze::LowerMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2( blaze::DynamicMatrix<int>( 3UL, 3UL, 1 ) );

      mat2 %= mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major CompressedMatrix dense matrix Schur product assignment (lower)";

      blaze::LowerMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2( blaze::DynamicMatrix<int>( 3UL, 3UL, 1 ) );

      mat2 %= mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major CompressedMatrix dense matrix Schur product assignment (upper)";

      blaze::UpperMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2( blaze::DynamicMatrix<int>( 3UL, 3UL, 1 ) );

      mat2 %= mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major CompressedMatrix dense matrix Schur product assignment (upper)";

      blaze::UpperMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2( blaze::DynamicMatrix<int>( 3UL, 3UL, 1 ) );

      mat2 %= mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major CompressedMatrix dense matrix Schur product assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2( blaze::DynamicMatrix<int>( 3UL, 3UL, 1 ) );

      mat2 %= mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major CompressedMatrix dense matrix Schur product assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2( blaze::DynamicMatrix<int>( 3UL, 3UL, 1 ) );

      mat2 %= mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse matrix Schur product assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major CompressedMatrix sparse matrix Schur product assignment";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1{ { 0, 0, 1 },
                                                         { 2, 3, 4 } };

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2{ { 5, 0, 6 },
                                                            { 0, 7, 0 } };

      mat2 %= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkNonZeros( mat2, 2UL );
      checkNonZeros( mat2, 0UL, 0UL );
      checkNonZeros( mat2, 1UL, 1UL );
      checkNonZeros( mat2, 2UL, 1UL );

      if( mat2(0,0) != 0 || mat2(0,1) !=  0 || mat2(0,2) != 6 ||
          mat2(1,0) != 0 || mat2(1,1) != 21 || mat2(1,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 0  0 6 )\n( 0 21 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major CompressedMatrix sparse matrix Schur product assignment";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1{ { 0, 0, 1 },
                                                            { 2, 3, 4 } };

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2{ { 5, 0, 6 },
                                                            { 0, 7, 0 } };

      mat2 %= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkNonZeros( mat2, 2UL );
      checkNonZeros( mat2, 0UL, 0UL );
      checkNonZeros( mat2, 1UL, 1UL );
      checkNonZeros( mat2, 2UL, 1UL );

      if( mat2(0,0) != 0 || mat2(0,1) !=  0 || mat2(0,2) != 6 ||
          mat2(1,0) != 0 || mat2(1,1) != 21 || mat2(1,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 0  0 6 )\n( 0 21 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major CompressedMatrix sparse matrix Schur product assignment (lower)";

      blaze::LowerMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2( blaze::DynamicMatrix<int>( 3UL, 3UL, 1 ) );

      mat2 %= mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CompressedMatrix sparse matrix Schur product assignment (lower)";

      blaze::LowerMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2( blaze::DynamicMatrix<int>( 3UL, 3UL, 1 ) );

      mat2 %= mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major CompressedMatrix sparse matrix Schur product assignment (upper)";

      blaze::UpperMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2( blaze::DynamicMatrix<int>( 3UL, 3UL, 1 ) );

      mat2 %= mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CompressedMatrix sparse matrix Schur product assignment (upper)";

      blaze::UpperMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2( blaze::DynamicMatrix<int>( 3UL, 3UL, 1 ) );

      mat2 %= mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major CompressedMatrix sparse matrix Schur product assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2( blaze::DynamicMatrix<int>( 3UL, 3UL, 1 ) );

      mat2 %= mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CompressedMatrix sparse matrix Schur product assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2( blaze::DynamicMatrix<int>( 3UL, 3UL, 1 ) );

      mat2 %= mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the CompressedMatrix multiplication assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the multiplication assignment operators of the CompressedMatrix
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testMultAssign()
{
   //=====================================================================================
   // Row-major dense matrix multiplication assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major CompressedMatrix dense matrix multiplication assignment";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat1{ { 0, 2, 0, 0 },
                                                      { 1, 3, 0, 4 },
                                                      { 0, 0, 0, 5 } };

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2{ { 1, 0, 2 },
                                                         { 0, 3, 0 },
                                                         { 4, 0, 5 } };

      mat2 *= mat1;

      checkRows    ( mat2, 3UL );
      checkColumns ( mat2, 4UL );
      checkNonZeros( mat2, 7UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 3UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 0 || mat2(0,1) != 2 || mat2(0,2) != 0 || mat2(0,3) != 10 ||
          mat2(1,0) != 3 || mat2(1,1) != 9 || mat2(1,2) != 0 || mat2(1,3) != 12 ||
          mat2(2,0) != 0 || mat2(2,1) != 8 || mat2(2,2) != 0 || mat2(2,3) != 25 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 0 2 0 10 )\n( 3 9 0 12 )\n( 0 8 0 25 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CompressedMatrix dense matrix multiplication assignment";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat1{ { 0, 2, 0, 0 },
                                                         { 1, 3, 0, 4 },
                                                         { 0, 0, 0, 5 } };

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2{ { 1, 0, 2 },
                                                         { 0, 3, 0 },
                                                         { 4, 0, 5 } };

      mat2 *= mat1;

      checkRows    ( mat2, 3UL );
      checkColumns ( mat2, 4UL );
      checkNonZeros( mat2, 7UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 3UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 0 || mat2(0,1) != 2 || mat2(0,2) != 0 || mat2(0,3) != 10 ||
          mat2(1,0) != 3 || mat2(1,1) != 9 || mat2(1,2) != 0 || mat2(1,3) != 12 ||
          mat2(2,0) != 0 || mat2(2,1) != 8 || mat2(2,2) != 0 || mat2(2,3) != 25 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 0 2 0 10 )\n( 3 9 0 12 )\n( 0 8 0 25 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix multiplication assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major CompressedMatrix sparse matrix multiplication assignment";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1{ { 0, 2, 0, 0 },
                                                         { 1, 3, 0, 4 },
                                                         { 0, 0, 0, 5 } };

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2{ { 1, 0, 2 },
                                                         { 0, 3, 0 },
                                                         { 4, 0, 5 } };

      mat2 *= mat1;

      checkRows    ( mat2, 3UL );
      checkColumns ( mat2, 4UL );
      checkNonZeros( mat2, 7UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 3UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 0 || mat2(0,1) != 2 || mat2(0,2) != 0 || mat2(0,3) != 10 ||
          mat2(1,0) != 3 || mat2(1,1) != 9 || mat2(1,2) != 0 || mat2(1,3) != 12 ||
          mat2(2,0) != 0 || mat2(2,1) != 8 || mat2(2,2) != 0 || mat2(2,3) != 25 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 0 2 0 10 )\n( 3 9 0 12 )\n( 0 8 0 25 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CompressedMatrix sparse matrix multiplication assignment";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1{ { 0, 2, 0, 0 },
                                                            { 1, 3, 0, 4 },
                                                            { 0, 0, 0, 5 } };

      blaze::CompressedMatrix<int,blaze::rowMajor> mat2{ { 1, 0, 2 },
                                                         { 0, 3, 0 },
                                                         { 4, 0, 5 } };

      mat2 *= mat1;

      checkRows    ( mat2, 3UL );
      checkColumns ( mat2, 4UL );
      checkNonZeros( mat2, 7UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 3UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 0 || mat2(0,1) != 2 || mat2(0,2) != 0 || mat2(0,3) != 10 ||
          mat2(1,0) != 3 || mat2(1,1) != 9 || mat2(1,2) != 0 || mat2(1,3) != 12 ||
          mat2(2,0) != 0 || mat2(2,1) != 8 || mat2(2,2) != 0 || mat2(2,3) != 25 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 0 2 0 10 )\n( 3 9 0 12 )\n( 0 8 0 25 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix multiplication assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major CompressedMatrix dense matrix multiplication assignment";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat1{ { 0, 2, 0, 0 },
                                                      { 1, 3, 0, 4 },
                                                      { 0, 0, 0, 5 } };

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2{ { 1, 0, 2 },
                                                            { 0, 3, 0 },
                                                            { 4, 0, 5 } };

      mat2 *= mat1;

      checkRows    ( mat2, 3UL );
      checkColumns ( mat2, 4UL );
      checkNonZeros( mat2, 7UL );
      checkNonZeros( mat2, 0UL, 1UL );
      checkNonZeros( mat2, 1UL, 3UL );
      checkNonZeros( mat2, 2UL, 0UL );
      checkNonZeros( mat2, 3UL, 3UL );

      if( mat2(0,0) != 0 || mat2(0,1) != 2 || mat2(0,2) != 0 || mat2(0,3) != 10 ||
          mat2(1,0) != 3 || mat2(1,1) != 9 || mat2(1,2) != 0 || mat2(1,3) != 12 ||
          mat2(2,0) != 0 || mat2(2,1) != 8 || mat2(2,2) != 0 || mat2(2,3) != 25 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 0 2 0 10 )\n( 3 9 0 12 )\n( 0 8 0 25 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CompressedMatrix dense matrix multiplication assignment";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat1{ { 0, 2, 0, 0 },
                                                         { 1, 3, 0, 4 },
                                                         { 0, 0, 0, 5 } };

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2{ { 1, 0, 2 },
                                                            { 0, 3, 0 },
                                                            { 4, 0, 5 } };

      mat2 *= mat1;

      checkRows    ( mat2, 3UL );
      checkColumns ( mat2, 4UL );
      checkNonZeros( mat2, 7UL );
      checkNonZeros( mat2, 0UL, 1UL );
      checkNonZeros( mat2, 1UL, 3UL );
      checkNonZeros( mat2, 2UL, 0UL );
      checkNonZeros( mat2, 3UL, 3UL );

      if( mat2(0,0) != 0 || mat2(0,1) != 2 || mat2(0,2) != 0 || mat2(0,3) != 10 ||
          mat2(1,0) != 3 || mat2(1,1) != 9 || mat2(1,2) != 0 || mat2(1,3) != 12 ||
          mat2(2,0) != 0 || mat2(2,1) != 8 || mat2(2,2) != 0 || mat2(2,3) != 25 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 0 2 0 10 )\n( 3 9 0 12 )\n( 0 8 0 25 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse matrix multiplication assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major CompressedMatrix sparse matrix multiplication assignment";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1{ { 0, 2, 0, 0 },
                                                         { 1, 3, 0, 4 },
                                                         { 0, 0, 0, 5 } };

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2{ { 1, 0, 2 },
                                                            { 0, 3, 0 },
                                                            { 4, 0, 5 } };

      mat2 *= mat1;

      checkRows    ( mat2, 3UL );
      checkColumns ( mat2, 4UL );
      checkNonZeros( mat2, 7UL );
      checkNonZeros( mat2, 0UL, 1UL );
      checkNonZeros( mat2, 1UL, 3UL );
      checkNonZeros( mat2, 2UL, 0UL );
      checkNonZeros( mat2, 3UL, 3UL );

      if( mat2(0,0) != 0 || mat2(0,1) != 2 || mat2(0,2) != 0 || mat2(0,3) != 10 ||
          mat2(1,0) != 3 || mat2(1,1) != 9 || mat2(1,2) != 0 || mat2(1,3) != 12 ||
          mat2(2,0) != 0 || mat2(2,1) != 8 || mat2(2,2) != 0 || mat2(2,3) != 25 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 0 2 0 10 )\n( 3 9 0 12 )\n( 0 8 0 25 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CompressedMatrix sparse matrix multiplication assignment";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1{ { 0, 2, 0, 0 },
                                                            { 1, 3, 0, 4 },
                                                            { 0, 0, 0, 5 } };

      blaze::CompressedMatrix<int,blaze::columnMajor> mat2{ { 1, 0, 2 },
                                                            { 0, 3, 0 },
                                                            { 4, 0, 5 } };

      mat2 *= mat1;

      checkRows    ( mat2, 3UL );
      checkColumns ( mat2, 4UL );
      checkNonZeros( mat2, 7UL );
      checkNonZeros( mat2, 0UL, 1UL );
      checkNonZeros( mat2, 1UL, 3UL );
      checkNonZeros( mat2, 2UL, 0UL );
      checkNonZeros( mat2, 3UL, 3UL );

      if( mat2(0,0) != 0 || mat2(0,1) != 2 || mat2(0,2) != 0 || mat2(0,3) != 10 ||
          mat2(1,0) != 3 || mat2(1,1) != 9 || mat2(1,2) != 0 || mat2(1,3) != 12 ||
          mat2(2,0) != 0 || mat2(2,1) != 8 || mat2(2,2) != 0 || mat2(2,3) != 25 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 0 2 0 10 )\n( 3 9 0 12 )\n( 0 8 0 25 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of all CompressedMatrix (self-)scaling operations.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all available ways to scale an instance of the CompressedMatrix
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testScaling()
{
   //=====================================================================================
   // Row-major self-scaling (M*=s)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (M*=s)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat{ {  0, 0, 0 },
                                                        {  0, 0, 1 },
                                                        { -2, 0, 3 } };

      mat *= 2;

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 3UL );
      checkNonZeros( mat, 3UL );
      checkNonZeros( mat, 0UL, 0UL );
      checkNonZeros( mat, 1UL, 1UL );
      checkNonZeros( mat, 2UL, 2UL );

      if( mat(0,0) !=  0 || mat(0,1) != 0 || mat(0,2) != 0 ||
          mat(1,0) !=  0 || mat(1,1) != 0 || mat(1,2) != 2 ||
          mat(2,0) != -4 || mat(2,1) != 0 || mat(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  0 0 0 )\n(  0 0 2 )\n( -4 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major self-scaling (M=M*s)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (M=M*s)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat{ {  0, 0, 0 },
                                                        {  0, 0, 1 },
                                                        { -2, 0, 3 } };

      mat = mat * 2;

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 3UL );
      checkNonZeros( mat, 3UL );
      checkNonZeros( mat, 0UL, 0UL );
      checkNonZeros( mat, 1UL, 1UL );
      checkNonZeros( mat, 2UL, 2UL );

      if( mat(0,0) !=  0 || mat(0,1) != 0 || mat(0,2) != 0 ||
          mat(1,0) !=  0 || mat(1,1) != 0 || mat(1,2) != 2 ||
          mat(2,0) != -4 || mat(2,1) != 0 || mat(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  0 0 0 )\n(  0 0 2 )\n( -4 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major self-scaling (M=s*M)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (M=s*M)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat{ {  0, 0, 0 },
                                                        {  0, 0, 1 },
                                                        { -2, 0, 3 } };

      mat = 2 * mat;

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 3UL );
      checkNonZeros( mat, 3UL );
      checkNonZeros( mat, 0UL, 0UL );
      checkNonZeros( mat, 1UL, 1UL );
      checkNonZeros( mat, 2UL, 2UL );

      if( mat(0,0) !=  0 || mat(0,1) != 0 || mat(0,2) != 0 ||
          mat(1,0) !=  0 || mat(1,1) != 0 || mat(1,2) != 2 ||
          mat(2,0) != -4 || mat(2,1) != 0 || mat(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  0 0 0 )\n(  0 0 2 )\n( -4 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major self-scaling (M/=s)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (M/=s)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat{ {  0, 0, 0 },
                                                        {  0, 0, 2 },
                                                        { -4, 0, 6 } };

      mat /= 2;

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 3UL );
      checkNonZeros( mat, 3UL );
      checkNonZeros( mat, 0UL, 0UL );
      checkNonZeros( mat, 1UL, 1UL );
      checkNonZeros( mat, 2UL, 2UL );

      if( mat(0,0) !=  0 || mat(0,1) != 0 || mat(0,2) != 0 ||
          mat(1,0) !=  0 || mat(1,1) != 0 || mat(1,2) != 1 ||
          mat(2,0) != -2 || mat(2,1) != 0 || mat(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  0 0 0 )\n(  0 0 1 )\n( -2 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major self-scaling (M=M/s)
   //=====================================================================================

   {
      test_ = "Row-major self-scaling (M=M/s)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat{ {  0, 0, 0 },
                                                        {  0, 0, 2 },
                                                        { -4, 0, 6 } };

      mat = mat / 2;

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 3UL );
      checkNonZeros( mat, 3UL );
      checkNonZeros( mat, 0UL, 0UL );
      checkNonZeros( mat, 1UL, 1UL );
      checkNonZeros( mat, 2UL, 2UL );

      if( mat(0,0) !=  0 || mat(0,1) != 0 || mat(0,2) != 0 ||
          mat(1,0) !=  0 || mat(1,1) != 0 || mat(1,2) != 1 ||
          mat(2,0) != -2 || mat(2,1) != 0 || mat(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  0 0 0 )\n(  0 0 1 )\n( -2 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major CompressedMatrix::scale()
   //=====================================================================================

   {
      test_ = "Row-major CompressedMatrix::scale() (int)";

      // Initialization check
      blaze::CompressedMatrix<int,blaze::rowMajor> mat{ { 1, 0 },
                                                        { 2, 0 },
                                                        { 0, 3 } };

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 3UL );
      checkNonZeros( mat, 3UL );
      checkNonZeros( mat, 0UL, 1UL );
      checkNonZeros( mat, 1UL, 1UL );
      checkNonZeros( mat, 2UL, 1UL );

      if( mat(0,0) != 1 || mat(1,0) != 2 || mat(2,1) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 0 )\n( 2 0 )\n( 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Integral scaling of the matrix
      mat.scale( 2 );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 3UL );
      checkNonZeros( mat, 3UL );
      checkNonZeros( mat, 0UL, 1UL );
      checkNonZeros( mat, 1UL, 1UL );
      checkNonZeros( mat, 2UL, 1UL );

      if( mat(0,0) != 2 || mat(1,0) != 4 || mat(2,1) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 2 0 )\n( 4 0 )\n( 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Floating point scaling of the matrix
      mat.scale( 0.5 );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 3UL );
      checkNonZeros( mat, 3UL );
      checkNonZeros( mat, 0UL, 1UL );
      checkNonZeros( mat, 1UL, 1UL );
      checkNonZeros( mat, 2UL, 1UL );

      if( mat(0,0) != 1 || mat(1,0) != 2 || mat(2,1) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 0 )\n( 2 0 )\n( 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major CompressedMatrix::scale() (complex)";

      using blaze::complex;

      blaze::CompressedMatrix<complex<float>,blaze::rowMajor> mat( 2UL, 2UL );
      mat(0,0) = complex<float>( 1.0F, 0.0F );
      mat(1,1) = complex<float>( 2.0F, 0.0F );
      mat.scale( complex<float>( 3.0F, 0.0F ) );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 2UL );
      checkNonZeros( mat, 2UL );
      checkNonZeros( mat, 0UL, 1UL );
      checkNonZeros( mat, 1UL, 1UL );

      if( mat(0,0) != complex<float>( 3.0F, 0.0F ) || mat(1,1) != complex<float>( 6.0F, 0.0F ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( (3,0) (0,0)\n(0,0) (6,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (M*=s)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (M*=s)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat{ {  0, 0, 0 },
                                                           {  0, 0, 1 },
                                                           { -2, 0, 3 } };

      mat *= 2;

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 3UL );
      checkNonZeros( mat, 3UL );
      checkNonZeros( mat, 0UL, 1UL );
      checkNonZeros( mat, 1UL, 0UL );
      checkNonZeros( mat, 2UL, 2UL );

      if( mat(0,0) !=  0 || mat(0,1) != 0 || mat(0,2) != 0 ||
          mat(1,0) !=  0 || mat(1,1) != 0 || mat(1,2) != 2 ||
          mat(2,0) != -4 || mat(2,1) != 0 || mat(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  0 0 0 )\n(  0 0 2 )\n( -4 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (M=M*s)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (M=M*s)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat{ {  0, 0, 0 },
                                                           {  0, 0, 1 },
                                                           { -2, 0, 3 } };

      mat = mat * 2;

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 3UL );
      checkNonZeros( mat, 3UL );
      checkNonZeros( mat, 0UL, 1UL );
      checkNonZeros( mat, 1UL, 0UL );
      checkNonZeros( mat, 2UL, 2UL );

      if( mat(0,0) !=  0 || mat(0,1) != 0 || mat(0,2) != 0 ||
          mat(1,0) !=  0 || mat(1,1) != 0 || mat(1,2) != 2 ||
          mat(2,0) != -4 || mat(2,1) != 0 || mat(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  0 0 0 )\n(  0 0 2 )\n( -4 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (M=s*M)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (M=s*M)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat{ {  0, 0, 0 },
                                                           {  0, 0, 1 },
                                                           { -2, 0, 3 } };

      mat = 2 * mat;

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 3UL );
      checkNonZeros( mat, 3UL );
      checkNonZeros( mat, 0UL, 1UL );
      checkNonZeros( mat, 1UL, 0UL );
      checkNonZeros( mat, 2UL, 2UL );

      if( mat(0,0) !=  0 || mat(0,1) != 0 || mat(0,2) != 0 ||
          mat(1,0) !=  0 || mat(1,1) != 0 || mat(1,2) != 2 ||
          mat(2,0) != -4 || mat(2,1) != 0 || mat(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  0 0 0 )\n(  0 0 2 )\n( -4 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (M/=s)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (M/=s)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat{ {  0, 0, 0 },
                                                           {  0, 0, 2 },
                                                           { -4, 0, 6 } };

      mat /= 2;

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 3UL );
      checkNonZeros( mat, 3UL );
      checkNonZeros( mat, 0UL, 1UL );
      checkNonZeros( mat, 1UL, 0UL );
      checkNonZeros( mat, 2UL, 2UL );

      if( mat(0,0) !=  0 || mat(0,1) != 0 || mat(0,2) != 0 ||
          mat(1,0) !=  0 || mat(1,1) != 0 || mat(1,2) != 1 ||
          mat(2,0) != -2 || mat(2,1) != 0 || mat(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  0 0 0 )\n(  0 0 1 )\n( -2 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major self-scaling (M=M/s)
   //=====================================================================================

   {
      test_ = "Column-major self-scaling (M=M/s)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat{ {  0, 0, 0 },
                                                           {  0, 0, 2 },
                                                           { -4, 0, 6 } };

      mat = mat / 2;

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 3UL );
      checkNonZeros( mat, 3UL );
      checkNonZeros( mat, 0UL, 1UL );
      checkNonZeros( mat, 1UL, 0UL );
      checkNonZeros( mat, 2UL, 2UL );

      if( mat(0,0) !=  0 || mat(0,1) != 0 || mat(0,2) != 0 ||
          mat(1,0) !=  0 || mat(1,1) != 0 || mat(1,2) != 1 ||
          mat(2,0) != -2 || mat(2,1) != 0 || mat(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(  0 0 0 )\n(  0 0 1 )\n( -2 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major CompressedMatrix::scale()
   //=====================================================================================

   {
      test_ = "Column-major CompressedMatrix::scale() (int)";

      // Initialization check
      blaze::CompressedMatrix<int,blaze::columnMajor> mat{ { 1, 0 },
                                                           { 2, 0 },
                                                           { 0, 3 } };

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 3UL );
      checkNonZeros( mat, 3UL );
      checkNonZeros( mat, 0UL, 2UL );
      checkNonZeros( mat, 1UL, 1UL );

      if( mat(0,0) != 1 || mat(1,0) != 2 || mat(2,1) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 0 )\n( 2 0 )\n( 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Integral scaling of the matrix
      mat.scale( 2 );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 3UL );
      checkNonZeros( mat, 3UL );
      checkNonZeros( mat, 0UL, 2UL );
      checkNonZeros( mat, 1UL, 1UL );

      if( mat(0,0) != 2 || mat(1,0) != 4 || mat(2,1) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 2 0 )\n( 4 0 )\n( 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Floating point scaling of the matrix
      mat.scale( 0.5 );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 3UL );
      checkNonZeros( mat, 3UL );
      checkNonZeros( mat, 0UL, 2UL );
      checkNonZeros( mat, 1UL, 1UL );

      if( mat(0,0) != 1 || mat(1,0) != 2 || mat(2,1) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 0 )\n( 2 0 )\n( 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major CompressedMatrix::scale() (complex)";

      using blaze::complex;

      blaze::CompressedMatrix<complex<float>,blaze::columnMajor> mat( 2UL, 2UL );
      mat(0,0) = complex<float>( 1.0F, 0.0F );
      mat(1,1) = complex<float>( 2.0F, 0.0F );
      mat.scale( complex<float>( 3.0F, 0.0F ) );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 2UL );
      checkCapacity( mat, 2UL );
      checkNonZeros( mat, 2UL );
      checkNonZeros( mat, 0UL, 1UL );
      checkNonZeros( mat, 1UL, 1UL );

      if( mat(0,0) != complex<float>( 3.0F, 0.0F ) || mat(1,1) != complex<float>( 6.0F, 0.0F ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( (3,0) (0,0)\n(0,0) (6,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************

} // namespace compressedmatrix

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
   std::cout << "   Running CompressedMatrix class test (part 1)..." << std::endl;

   try
   {
      RUN_COMPRESSEDMATRIX_CLASS_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during CompressedMatrix class test (part 1):\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
