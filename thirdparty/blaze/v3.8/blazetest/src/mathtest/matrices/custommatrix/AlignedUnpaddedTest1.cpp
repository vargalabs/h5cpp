//=================================================================================================
/*!
//  \file src/mathtest/matrices/custommatrix/AlignedUnpaddedTest1.cpp
//  \brief Source file for the aligned/unpadded CustomMatrix class test (part 1)
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

#include <array>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <blaze/math/CompressedMatrix.h>
#include <blaze/math/DiagonalMatrix.h>
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/LowerMatrix.h>
#include <blaze/math/shims/NextMultiple.h>
#include <blaze/math/UpperMatrix.h>
#include <blaze/util/Complex.h>
#include <blaze/util/Memory.h>
#include <blaze/util/policies/Deallocate.h>
#include <blaze/util/Random.h>
#include <blaze/util/typetraits/AlignmentOf.h>
#include <blazetest/mathtest/matrices/custommatrix/AlignedUnpaddedTest.h>
#include <blazetest/mathtest/RandomMaximum.h>
#include <blazetest/mathtest/RandomMinimum.h>

#ifdef BLAZE_USE_HPX_THREADS
#  include <hpx/hpx_main.hpp>
#endif


namespace blazetest {

namespace mathtest {

namespace matrices {

namespace custommatrix {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the CustomMatrix class test.
//
// \exception std::runtime_error Operation error detected.
*/
AlignedUnpaddedTest::AlignedUnpaddedTest()
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
/*!\brief Test of the CustomMatrix constructors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all constructors of the CustomMatrix class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void AlignedUnpaddedTest::testConstructors()
{
   //=====================================================================================
   // Row-major default constructor
   //=====================================================================================

   {
      test_ = "Row-major CustomMatrix default constructor";

      MT mat;

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );
   }


   //=====================================================================================
   // Row-major constructor ( Type*, size_t, size_t )
   //=====================================================================================

   {
      test_ = "Row-major CustomMatrix constructor ( Type*, size_t, size_t )";

      // Constructor a 2x16 custom matrix
      {
         std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 32UL ) );
         MT mat( memory.get(), 2UL, 16UL );

         checkRows    ( mat,  2UL );
         checkColumns ( mat, 16UL );
         checkCapacity( mat, 32UL );
      }

      // Trying to construct a custom matrix with invalid array of elements
      try {
         MT mat( nullptr, 0UL, 0UL );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Constructing a custom matrix with a nullptr succeeded\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      // Trying to construct a custom matrix with invalid alignment
      if( blaze::AlignmentOf<int>::value > sizeof(int) )
      {
         try {
            std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 33UL ) );
            MT mat( memory.get()+1UL, 2UL, 16UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Constructing a custom matrix with invalid alignment succeeded\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Trying to construct a custom matrix with invalid row alignment
      if( blaze::AlignmentOf<int>::value > sizeof(int) )
      {
         try {
            std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 30UL ) );
            MT mat( memory.get(), 2UL, 15UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Constructing a custom matrix with invalid row alignment succeeded\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }


   //=====================================================================================
   // Row-major constructor ( Type*, size_t, size_t, size_t )
   //=====================================================================================

   {
      test_ = "Row-major CustomMatrix constructor ( Type*, size_t, size_t, size_t )";

      // Constructor a 2x3 custom matrix
      {
         std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 32UL ) );
         MT mat( memory.get(), 2UL, 3UL, 16UL );

         checkRows    ( mat,  2UL );
         checkColumns ( mat,  3UL );
         checkCapacity( mat, 32UL );
      }

      // Trying to construct a custom matrix with invalid array of elements
      try {
         MT mat( nullptr, 0UL, 0UL, 0UL );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Constructing a custom matrix with a nullptr succeeded\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      // Trying to construct a custom matrix with invalid alignment
      if( blaze::AlignmentOf<int>::value > sizeof(int) )
      {
         try {
            std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 33UL ) );
            MT mat( memory.get()+1UL, 2UL, 2UL, 16UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Constructing a custom matrix with invalid alignment succeeded\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Trying to construct a custom matrix with invalid row alignment
      if( blaze::AlignmentOf<int>::value > sizeof(int) )
      {
         try {
            std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 30UL ) );
            MT mat( memory.get(), 2UL, 2UL, 15UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Constructing a custom matrix with invalid row alignment succeeded\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }


   //=====================================================================================
   // Row-major copy constructor
   //=====================================================================================

   {
      test_ = "Row-major CustomMatrix copy constructor (0x0)";

      MT mat1;
      MT mat2( mat1 );

      checkRows    ( mat2, 0UL );
      checkColumns ( mat2, 0UL );
      checkNonZeros( mat2, 0UL );
   }

   {
      test_ = "Row-major CustomMatrix copy constructor (0x3)";

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 16UL ) );
      MT mat1( memory.get(), 0UL, 3UL, 16UL );
      MT mat2( mat1 );

      checkRows    ( mat2, 0UL );
      checkColumns ( mat2, 3UL );
      checkNonZeros( mat2, 0UL );
   }

   {
      test_ = "Row-major CustomMatrix copy constructor (2x0)";

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 10UL ) );
      MT mat1( memory.get(), 2UL, 0UL, 0UL );
      MT mat2( mat1 );

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 0UL );
      checkNonZeros( mat2, 0UL );
   }

   {
      test_ = "Row-major CustomMatrix copy constructor (2x3)";

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 32UL ) );
      MT mat1( memory.get(), 2UL, 3UL, 16UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(0,2) = 3;
      mat1(1,0) = 4;
      mat1(1,1) = 5;
      mat1(1,2) = 6;

      MT mat2( mat1 );

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 32UL );
      checkNonZeros( mat2,  6UL );
      checkNonZeros( mat2,  0UL, 3UL );
      checkNonZeros( mat2,  1UL, 3UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 3 ||
          mat2(1,0) != 4 || mat2(1,1) != 5 || mat2(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major move constructor
   //=====================================================================================

   {
      test_ = "Row-major CustomMatrix move constructor (0x0)";

      MT mat1;
      MT mat2( std::move( mat1 ) );

      checkRows    ( mat2, 0UL );
      checkColumns ( mat2, 0UL );
      checkNonZeros( mat2, 0UL );
   }

   {
      test_ = "Row-major CustomMatrix move constructor (0x3)";

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 16UL ) );
      MT mat1( memory.get(), 0UL, 3UL, 16UL );
      MT mat2( std::move( mat1 ) );

      checkRows    ( mat2, 0UL );
      checkColumns ( mat2, 3UL );
      checkNonZeros( mat2, 0UL );
   }

   {
      test_ = "Row-major CustomMatrix move constructor (2x0)";

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 10UL ) );
      MT mat1( memory.get(), 2UL, 0UL, 0UL );
      MT mat2( std::move( mat1 ) );

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 0UL );
      checkNonZeros( mat2, 0UL );
   }

   {
      test_ = "Row-major CustomMatrix move constructor (2x3)";

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 32UL ) );
      MT mat1( memory.get(), 2UL, 3UL, 16UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(0,2) = 3;
      mat1(1,0) = 4;
      mat1(1,1) = 5;
      mat1(1,2) = 6;

      MT mat2( std::move( mat1 ) );

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 32UL );
      checkNonZeros( mat2,  6UL );
      checkNonZeros( mat2,  0UL, 3UL );
      checkNonZeros( mat2,  1UL, 3UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 3 ||
          mat2(1,0) != 4 || mat2(1,1) != 5 || mat2(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major default constructor
   //=====================================================================================

   {
      test_ = "Column-major CustomMatrix default constructor";

      OMT mat;

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );
   }


   //=====================================================================================
   // Column-major constructor ( Type*, size_t, size_t )
   //=====================================================================================

   {
      test_ = "Column-major CustomMatrix constructor ( Type*, size_t, size_t )";

      // Constructor a 2x3 custom matrix
      {
         std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
         OMT mat( memory.get(), 16UL, 3UL );

         checkRows    ( mat, 16UL );
         checkColumns ( mat,  3UL );
         checkCapacity( mat, 48UL );
      }

      // Trying to construct a custom matrix with invalid array of elements
      try {
         OMT mat( nullptr, 0UL, 0UL );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Constructing a custom matrix with a nullptr succeeded\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      // Trying to construct a custom matrix with invalid alignment
      if( blaze::AlignmentOf<int>::value > sizeof(int) )
      {
         try {
            std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 33UL ) );
            OMT mat( memory.get()+1UL, 16UL, 2UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Constructing a custom matrix with invalid alignment succeeded\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Trying to construct a custom matrix with invalid column alignment
      if( blaze::AlignmentOf<int>::value > sizeof(int) )
      {
         try {
            std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 30UL ) );
            OMT mat( memory.get(), 15UL, 2UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Constructing a custom matrix with invalid column alignment succeeded\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }


   //=====================================================================================
   // Column-major constructor ( Type*, size_t, size_t, size_t )
   //=====================================================================================

   {
      test_ = "Column-major CustomMatrix constructor ( Type*, size_t, size_t, size_t )";

      // Constructor a 2x3 custom matrix
      {
         std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
         OMT mat( memory.get(), 2UL, 3UL, 16UL );

         checkRows    ( mat,  2UL );
         checkColumns ( mat,  3UL );
         checkCapacity( mat, 48UL );
      }

      // Trying to construct a custom matrix with invalid array of elements
      try {
         OMT mat( nullptr, 0UL, 0UL, 0UL );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Constructing a custom matrix with a nullptr succeeded\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      // Trying to construct a custom matrix with invalid alignment
      if( blaze::AlignmentOf<int>::value > sizeof(int) )
      {
         try {
            std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 33UL ) );
            OMT mat( memory.get()+1UL, 2UL, 2UL, 16UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Constructing a custom matrix with invalid alignment succeeded\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      // Trying to construct a custom matrix with invalid column alignment
      if( blaze::AlignmentOf<int>::value > sizeof(int) )
      {
         try {
            std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 30UL ) );
            OMT mat( memory.get(), 2UL, 2UL, 15UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Constructing a custom matrix with invalid column alignment succeeded\n"
                << " Details:\n"
                << "   Result:\n" << mat << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }


   //=====================================================================================
   // Column-major copy constructor
   //=====================================================================================

   {
      test_ = "Column-major CustomMatrix copy constructor (0x0)";

      OMT mat1;
      OMT mat2( mat1 );

      checkRows    ( mat2, 0UL );
      checkColumns ( mat2, 0UL );
      checkNonZeros( mat2, 0UL );
   }

   {
      test_ = "Column-major CustomMatrix copy constructor (0x3)";

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 10UL ) );
      OMT mat1( memory.get(), 0UL, 3UL, 0UL );
      OMT mat2( mat1 );

      checkRows    ( mat2, 0UL );
      checkColumns ( mat2, 3UL );
      checkNonZeros( mat2, 0UL );
   }

   {
      test_ = "Column-major CustomMatrix copy constructor (2x0)";

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 16UL ) );
      OMT mat1( memory.get(), 2UL, 0UL, 16UL );
      OMT mat2( mat1 );

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 0UL );
      checkNonZeros( mat2, 0UL );
   }

   {
      test_ = "Column-major CustomMatrix copy constructor (2x3)";

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      OMT mat1( memory.get(), 2UL, 3UL, 16UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(0,2) = 3;
      mat1(1,0) = 4;
      mat1(1,1) = 5;
      mat1(1,2) = 6;

      OMT mat2( mat1 );

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 48UL );
      checkNonZeros( mat2,  6UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 2UL );
      checkNonZeros( mat2,  2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 3 ||
          mat2(1,0) != 4 || mat2(1,1) != 5 || mat2(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major move constructor
   //=====================================================================================

   {
      test_ = "Column-major CustomMatrix move constructor (0x0)";

      OMT mat1;
      OMT mat2( std::move( mat1 ) );

      checkRows    ( mat2, 0UL );
      checkColumns ( mat2, 0UL );
      checkNonZeros( mat2, 0UL );
   }

   {
      test_ = "Column-major CustomMatrix move constructor (0x3)";

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 10UL ) );
      OMT mat1( memory.get(), 0UL, 3UL, 0UL );
      OMT mat2( std::move( mat1 ) );

      checkRows    ( mat2, 0UL );
      checkColumns ( mat2, 3UL );
      checkNonZeros( mat2, 0UL );
   }

   {
      test_ = "Column-major CustomMatrix move constructor (2x0)";

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 16UL ) );
      OMT mat1( memory.get(), 2UL, 0UL, 16UL );
      OMT mat2( std::move( mat1 ) );

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 0UL );
      checkNonZeros( mat2, 0UL );
   }

   {
      test_ = "Column-major CustomMatrix move constructor (2x3)";

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      OMT mat1( memory.get(), 2UL, 3UL, 16UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(0,2) = 3;
      mat1(1,0) = 4;
      mat1(1,1) = 5;
      mat1(1,2) = 6;

      OMT mat2( std::move( mat1 ) );

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 48UL );
      checkNonZeros( mat2,  6UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 2UL );
      checkNonZeros( mat2,  2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 3 ||
          mat2(1,0) != 4 || mat2(1,1) != 5 || mat2(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the CustomMatrix assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all assignment operators of the CustomMatrix class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void AlignedUnpaddedTest::testAssignment()
{
   //=====================================================================================
   // Row-major homogeneous assignment
   //=====================================================================================

   {
      test_ = "Row-major CustomMatrix homogeneous assignment";

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      MT mat( memory.get(), 3UL, 4UL, 16UL );
      mat = 2;

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  4UL );
      checkCapacity( mat, 48UL );
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
   // Row-major list assignment
   //=====================================================================================

   {
      test_ = "Row-major CustomMatrix initializer list assignment (complete list)";

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 32UL ) );
      MT mat( memory.get(), 2UL, 3UL, 16UL );
      mat = { { 1, 2, 3 }, { 4, 5, 6 } };

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 6UL );
      checkNonZeros( mat, 0UL, 3UL );
      checkNonZeros( mat, 1UL, 3UL );

      if( mat(0,0) != 1 || mat(0,1) != 2 || mat(0,2) != 3 ||
          mat(1,0) != 4 || mat(1,1) != 5 || mat(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major StaticMatrix initializer list assignment (incomplete list)";

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 32UL ) );
      MT mat( memory.get(), 2UL, 3UL, 16UL );
      mat = { { 1 }, { 4, 5, 6 } };

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 4UL );
      checkNonZeros( mat, 0UL, 1UL );
      checkNonZeros( mat, 1UL, 3UL );

      if( mat(0,0) != 1 || mat(0,1) != 0 || mat(0,2) != 0 ||
          mat(1,0) != 4 || mat(1,1) != 5 || mat(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major array assignment
   //=====================================================================================

   {
      test_ = "Row-major CustomMatrix static array assignment";

      const int array[2UL][3UL] = { { 1, 2, 3 }, { 4, 5, 6 } };
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 32UL ) );
      MT mat( memory.get(), 2UL, 3UL, 16UL );
      mat = array;

      checkRows    ( mat,  2UL );
      checkColumns ( mat,  3UL );
      checkCapacity( mat, 32UL );
      checkNonZeros( mat,  6UL );
      checkNonZeros( mat,  0UL, 3UL );
      checkNonZeros( mat,  1UL, 3UL );

      if( mat(0,0) != 1 || mat(0,1) != 2 || mat(0,2) != 3 ||
          mat(1,0) != 4 || mat(1,1) != 5 || mat(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major CustomMatrix std::array assignment";

      const std::array<std::array<int,3UL>,2UL> array{ { { 1, 2, 3 }, { 4, 5, 6 } } };
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 32UL ) );
      MT mat( memory.get(), 2UL, 3UL, 16UL );
      mat = array;

      checkRows    ( mat,  2UL );
      checkColumns ( mat,  3UL );
      checkCapacity( mat, 32UL );
      checkNonZeros( mat,  6UL );
      checkNonZeros( mat,  0UL, 3UL );
      checkNonZeros( mat,  1UL, 3UL );

      if( mat(0,0) != 1 || mat(0,1) != 2 || mat(0,2) != 3 ||
          mat(1,0) != 4 || mat(1,1) != 5 || mat(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major copy assignment
   //=====================================================================================

   {
      test_ = "Row-major CustomMatrix copy assignment";

      std::unique_ptr<int[],blaze::Deallocate> memory1( blaze::allocate<int>( 32UL ) );
      MT mat1( memory1.get(), 2UL, 3UL, 16UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(0,2) = 3;
      mat1(1,0) = 4;
      mat1(1,1) = 5;
      mat1(1,2) = 6;

      std::unique_ptr<int[],blaze::Deallocate> memory2( blaze::allocate<int>( 32UL ) );
      MT mat2( memory2.get(), 2UL, 3UL, 16UL );
      mat2 = mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 32UL );
      checkNonZeros( mat2,  6UL );
      checkNonZeros( mat2,  0UL, 3UL );
      checkNonZeros( mat2,  1UL, 3UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 3 ||
          mat2(1,0) != 4 || mat2(1,1) != 5 || mat2(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major move assignment
   //=====================================================================================

   {
      test_ = "Row-major CustomMatrix move assignment";

      std::unique_ptr<int[],blaze::Deallocate> memory1( blaze::allocate<int>( 32UL ) );
      MT mat1( memory1.get(), 2UL, 3UL, 16UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(0,2) = 3;
      mat1(1,0) = 4;
      mat1(1,1) = 5;
      mat1(1,2) = 6;

      std::unique_ptr<int[],blaze::Deallocate> memory2( blaze::allocate<int>( 32UL ) );
      MT mat2( memory2.get(), 2UL, 3UL, 16UL );
      mat2 = std::move( mat1 );

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 32UL );
      checkNonZeros( mat2,  6UL );
      checkNonZeros( mat2,  0UL, 3UL );
      checkNonZeros( mat2,  1UL, 3UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 3 ||
          mat2(1,0) != 4 || mat2(1,1) != 5 || mat2(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major dense matrix assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major CustomMatrix dense matrix assignment (mixed type)";

      using blaze::aligned;
      using blaze::unpadded;
      using blaze::rowMajor;

      using AlignedUnpadded = blaze::CustomMatrix<short,aligned,unpadded,rowMajor>;
      std::unique_ptr<short[],blaze::Deallocate> memory1( blaze::allocate<short>( 64UL ) );
      AlignedUnpadded mat1( memory1.get(), 2UL, 3UL, 32UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(0,2) = 3;
      mat1(1,0) = 4;
      mat1(1,1) = 5;
      mat1(1,2) = 6;

      std::unique_ptr<int[],blaze::Deallocate> memory2( blaze::allocate<int>( 32UL ) );
      MT mat2( memory2.get(), 2UL, 3UL, 16UL );
      mat2 = mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 32UL );
      checkNonZeros( mat2,  6UL );
      checkNonZeros( mat2,  0UL, 3UL );
      checkNonZeros( mat2,  1UL, 3UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 3 ||
          mat2(1,0) != 4 || mat2(1,1) != 5 || mat2(1,2) != 6 ) {
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
      test_ = "Row-major/row-major CustomMatrix dense matrix assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::rowMajor;

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,rowMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory1( blaze::allocate<int>( 32UL ) );
      AlignedPadded mat1( memory1.get(), 2UL, 3UL, 16UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(0,2) = 3;
      mat1(1,0) = 4;
      mat1(1,1) = 5;
      mat1(1,2) = 6;

      std::unique_ptr<int[],blaze::Deallocate> memory2( blaze::allocate<int>( 32UL ) );
      MT mat2( memory2.get(), 2UL, 3UL, 16UL );
      mat2 = mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 32UL );
      checkNonZeros( mat2,  6UL );
      checkNonZeros( mat2,  0UL, 3UL );
      checkNonZeros( mat2,  1UL, 3UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 3 ||
          mat2(1,0) != 4 || mat2(1,1) != 5 || mat2(1,2) != 6 ) {
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
      test_ = "Row-major/row-major CustomMatrix dense matrix assignment stress test (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::rowMajor;

      const int min( randmin );
      const int max( randmax );

      for( size_t i=0UL; i<10UL; ++i )
      {
         const size_t rows   ( blaze::rand<size_t>( 0UL, 16UL ) );
         const size_t columns( blaze::rand<size_t>( 0UL, 16UL ) );
         const size_t spacing( blaze::nextMultiple<size_t>( columns, 16UL ) );

         using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,rowMajor>;
         std::unique_ptr<int[],blaze::Deallocate> memory1( blaze::allocate<int>( rows*spacing ) );
         AlignedPadded mat1( memory1.get(), rows, columns, spacing );
         randomize( mat1, min, max );

         std::unique_ptr<int[],blaze::Deallocate> memory2( blaze::allocate<int>( rows*spacing ) );
         MT mat2( memory2.get(), rows, columns, spacing );
         mat2 = mat1;

         if( mat1 != mat2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << mat2 << "\n"
                << "   Expected result:\n" << mat1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   {
      test_ = "Row-major/row-major CustomMatrix dense matrix assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowMajor;

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>;
      std::unique_ptr<int[]> memory1( new int[7UL] );
      UnalignedUnpadded mat1( memory1.get()+1UL, 2UL, 3UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(0,2) = 3;
      mat1(1,0) = 4;
      mat1(1,1) = 5;
      mat1(1,2) = 6;

      std::unique_ptr<int[],blaze::Deallocate> memory2( blaze::allocate<int>( 32UL ) );
      MT mat2( memory2.get(), 2UL, 3UL, 16UL );
      mat2 = mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 32UL );
      checkNonZeros( mat2,  6UL );
      checkNonZeros( mat2,  0UL, 3UL );
      checkNonZeros( mat2,  1UL, 3UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 3 ||
          mat2(1,0) != 4 || mat2(1,1) != 5 || mat2(1,2) != 6 ) {
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
      test_ = "Row-major/row-major CustomMatrix dense matrix assignment stress test (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowMajor;

      const int min( randmin );
      const int max( randmax );

      for( size_t i=0UL; i<10UL; ++i )
      {
         const size_t rows   ( blaze::rand<size_t>( 0UL, 16UL ) );
         const size_t columns( blaze::rand<size_t>( 0UL, 16UL ) );
         const size_t spacing( blaze::nextMultiple<size_t>( columns, 16UL ) );

         using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>;
         std::unique_ptr<int[]> memory1( new int[rows*columns+1UL] );
         UnalignedUnpadded mat1( memory1.get()+1UL, rows, columns );
         randomize( mat1, min, max );

         std::unique_ptr<int[],blaze::Deallocate> memory2( blaze::allocate<int>( rows*spacing ) );
         MT mat2( memory2.get(), rows, columns, spacing );
         mat2 = mat1;

         if( mat1 != mat2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << mat2 << "\n"
                << "   Expected result:\n" << mat1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   {
      test_ = "Row-major/column-major CustomMatrix dense matrix assignment (mixed type)";

      using blaze::aligned;
      using blaze::unpadded;
      using blaze::columnMajor;

      using AlignedUnpadded = blaze::CustomMatrix<short,aligned,unpadded,columnMajor>;
      std::unique_ptr<short[],blaze::Deallocate> memory1( blaze::allocate<short>( 96UL ) );
      AlignedUnpadded mat1( memory1.get(), 2UL, 3UL, 32UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(0,2) = 3;
      mat1(1,0) = 4;
      mat1(1,1) = 5;
      mat1(1,2) = 6;

      std::unique_ptr<int[],blaze::Deallocate> memory2( blaze::allocate<int>( 32UL ) );
      MT mat2( memory2.get(), 2UL, 3UL, 16UL );
      mat2 = mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 32UL );
      checkNonZeros( mat2,  6UL );
      checkNonZeros( mat2,  0UL, 3UL );
      checkNonZeros( mat2,  1UL, 3UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 3 ||
          mat2(1,0) != 4 || mat2(1,1) != 5 || mat2(1,2) != 6 ) {
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
      test_ = "Row-major/column-major CustomMatrix dense matrix assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::columnMajor;

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,columnMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory1( blaze::allocate<int>( 48UL ) );
      AlignedPadded mat1( memory1.get(), 2UL, 3UL, 16UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(0,2) = 3;
      mat1(1,0) = 4;
      mat1(1,1) = 5;
      mat1(1,2) = 6;

      std::unique_ptr<int[],blaze::Deallocate> memory2( blaze::allocate<int>( 32UL ) );
      MT mat2( memory2.get(), 2UL, 3UL, 16UL );
      mat2 = mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 32UL );
      checkNonZeros( mat2,  6UL );
      checkNonZeros( mat2,  0UL, 3UL );
      checkNonZeros( mat2,  1UL, 3UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 3 ||
          mat2(1,0) != 4 || mat2(1,1) != 5 || mat2(1,2) != 6 ) {
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
      test_ = "Row-major/column-major CustomMatrix dense matrix assignment stress test (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::columnMajor;

      const int min( randmin );
      const int max( randmax );

      for( size_t i=0UL; i<10UL; ++i )
      {
         const size_t rows    ( blaze::rand<size_t>( 0UL, 16UL ) );
         const size_t columns ( blaze::rand<size_t>( 0UL, 16UL ) );
         const size_t spacing1( blaze::nextMultiple<size_t>( rows, 16UL ) );
         const size_t spacing2( blaze::nextMultiple<size_t>( columns, 16UL ) );

         using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,columnMajor>;
         std::unique_ptr<int[],blaze::Deallocate> memory1( blaze::allocate<int>( spacing1*columns ) );
         AlignedPadded mat1( memory1.get(), rows, columns, spacing1 );
         randomize( mat1, min, max );

         std::unique_ptr<int[],blaze::Deallocate> memory2( blaze::allocate<int>( rows*spacing2 ) );
         MT mat2( memory2.get(), rows, columns, spacing2 );
         mat2 = mat1;

         if( mat1 != mat2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << mat2 << "\n"
                << "   Expected result:\n" << mat1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   {
      test_ = "Row-major/column-major CustomMatrix dense matrix assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::columnMajor;

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>;
      std::unique_ptr<int[]> memory1( new int[7UL] );
      UnalignedUnpadded mat1( memory1.get()+1UL, 2UL, 3UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(0,2) = 3;
      mat1(1,0) = 4;
      mat1(1,1) = 5;
      mat1(1,2) = 6;

      std::unique_ptr<int[],blaze::Deallocate> memory2( blaze::allocate<int>( 32UL ) );
      MT mat2( memory2.get(), 2UL, 3UL, 16UL );
      mat2 = mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 32UL );
      checkNonZeros( mat2,  6UL );
      checkNonZeros( mat2,  0UL, 3UL );
      checkNonZeros( mat2,  1UL, 3UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 3 ||
          mat2(1,0) != 4 || mat2(1,1) != 5 || mat2(1,2) != 6 ) {
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
      test_ = "Row-major/column-major CustomMatrix dense matrix assignment stress test (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::columnMajor;

      const int min( randmin );
      const int max( randmax );

      for( size_t i=0UL; i<10UL; ++i )
      {
         const size_t rows   ( blaze::rand<size_t>( 0UL, 16UL ) );
         const size_t columns( blaze::rand<size_t>( 0UL, 16UL ) );
         const size_t spacing( blaze::nextMultiple<size_t>( columns, 16UL ) );

         using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>;
         std::unique_ptr<int[]> memory1( new int[rows*columns+1UL] );
         UnalignedUnpadded mat1( memory1.get()+1UL, rows, columns );
         randomize( mat1, min, max );

         std::unique_ptr<int[],blaze::Deallocate> memory2( blaze::allocate<int>( rows*spacing ) );
         MT mat2( memory2.get(), rows, columns, spacing );
         mat2 = mat1;

         if( mat1 != mat2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << mat2 << "\n"
                << "   Expected result:\n" << mat1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   {
      test_ = "Row-major/row-major CustomMatrix dense matrix assignment (lower)";

      blaze::LowerMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      MT mat2( memory.get(), 3UL, 3UL, 16UL );
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
      test_ = "Row-major/column-major CustomMatrix dense matrix assignment (lower)";

      blaze::LowerMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      MT mat2( memory.get(), 3UL, 3UL, 16UL );
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
      test_ = "Row-major/row-major CustomMatrix dense matrix assignment (upper)";

      blaze::UpperMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      MT mat2( memory.get(), 3UL, 3UL, 16UL );
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
      test_ = "Row-major/column-major CustomMatrix dense matrix assignment (upper)";

      blaze::UpperMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      MT mat2( memory.get(), 3UL, 3UL, 16UL );
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
      test_ = "Row-major/row-major CustomMatrix dense matrix assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      MT mat2( memory.get(), 3UL, 3UL, 16UL );
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
      test_ = "Row-major/column-major CustomMatrix dense matrix assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      MT mat2( memory.get(), 3UL, 3UL, 16UL );
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
      test_ = "Row-major/row-major CustomMatrix sparse matrix assignment";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1( 2UL, 3UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(1,0) = 3;
      mat1(1,2) = 4;

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 32UL ) );
      MT mat2( memory.get(), 2UL, 3UL, 16UL );
      mat2 = mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 32UL );
      checkNonZeros( mat2,  4UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 0 ||
          mat2(1,0) != 3 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 0 )\n( 3 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major CustomMatrix sparse matrix assignment stress test";

      const int min( randmin );
      const int max( randmax );

      for( size_t i=0UL; i<10UL; ++i )
      {
         const size_t rows   ( blaze::rand<size_t>( 0UL, 16UL ) );
         const size_t columns( blaze::rand<size_t>( 0UL, 16UL ) );
         const size_t spacing( blaze::nextMultiple<size_t>( columns, 16UL ) );

         blaze::CompressedMatrix<int,blaze::rowMajor> mat1( rows, columns );
         randomize( mat1, min, max );

         std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( rows*spacing ) );
         MT mat2( memory.get(), rows, columns, spacing );
         mat2 = mat1;

         if( mat1 != mat2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << mat2 << "\n"
                << "   Expected result:\n" << mat1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   {
      test_ = "Row-major/column-major CustomMatrix sparse matrix assignment";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1( 2UL, 3UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(1,0) = 3;
      mat1(1,2) = 4;

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 32UL ) );
      MT mat2( memory.get(), 2UL, 3UL, 16UL );
      mat2 = mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 32UL );
      checkNonZeros( mat2,  4UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 0 ||
          mat2(1,0) != 3 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 0 )\n( 3 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CustomMatrix sparse matrix assignment stress test";

      const int min( randmin );
      const int max( randmax );

      for( size_t i=0UL; i<10UL; ++i )
      {
         const size_t rows   ( blaze::rand<size_t>( 0UL, 16UL ) );
         const size_t columns( blaze::rand<size_t>( 0UL, 16UL ) );
         const size_t spacing( blaze::nextMultiple<size_t>( columns, 16UL ) );

         blaze::CompressedMatrix<int,blaze::columnMajor> mat1( rows, columns );
         randomize( mat1, min, max );

         std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( rows*spacing ) );
         MT mat2( memory.get(), rows, columns, spacing );
         mat2 = mat1;

         if( mat1 != mat2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << mat2 << "\n"
                << "   Expected result:\n" << mat1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   {
      test_ = "Row-major/row-major CustomMatrix sparse matrix assignment (lower)";

      blaze::LowerMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      MT mat2( memory.get(), 3UL, 3UL, 16UL );
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
      test_ = "Row-major/column-major CustomMatrix sparse matrix assignment (lower)";

      blaze::LowerMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      MT mat2( memory.get(), 3UL, 3UL, 16UL );
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
      test_ = "Row-major/row-major CustomMatrix sparse matrix assignment (upper)";

      blaze::UpperMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      MT mat2( memory.get(), 3UL, 3UL, 16UL );
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
      test_ = "Row-major/column-major CustomMatrix sparse matrix assignment (upper)";

      blaze::UpperMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      MT mat2( memory.get(), 3UL, 3UL, 16UL );
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
      test_ = "Row-major/row-major CustomMatrix sparse matrix assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      MT mat2( memory.get(), 3UL, 3UL, 16UL );
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
      test_ = "Row-major/column-major CustomMatrix sparse matrix assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      MT mat2( memory.get(), 3UL, 3UL, 16UL );
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
   // Column-major homogeneous assignment
   //=====================================================================================

   {
      test_ = "Column-major CustomMatrix homogeneous assigment";

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 64UL ) );
      OMT mat( memory.get(), 3UL, 4UL, 16UL );
      mat = 2;

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  4UL );
      checkCapacity( mat, 64UL );
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
   // Column-major list assignment
   //=====================================================================================

   {
      test_ = "Column-major CustomMatrix initializer list assignment (complete list)";

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      OMT mat( memory.get(), 2UL, 3UL, 16UL );
      mat = { { 1, 2, 3 }, { 4, 5, 6 } };

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 6UL );
      checkNonZeros( mat, 0UL, 2UL );
      checkNonZeros( mat, 1UL, 2UL );
      checkNonZeros( mat, 2UL, 2UL );

      if( mat(0,0) != 1 || mat(0,1) != 2 || mat(0,2) != 3 ||
          mat(1,0) != 4 || mat(1,1) != 5 || mat(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major CustomMatrix initializer list assignment (incomplete list)";

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      OMT mat( memory.get(), 2UL, 3UL, 16UL );
      mat = { { 1 }, { 4, 5, 6 } };

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 4UL );
      checkNonZeros( mat, 0UL, 2UL );
      checkNonZeros( mat, 1UL, 1UL );
      checkNonZeros( mat, 2UL, 1UL );

      if( mat(0,0) != 1 || mat(0,1) != 0 || mat(0,2) != 0 ||
          mat(1,0) != 4 || mat(1,1) != 5 || mat(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major array assignment
   //=====================================================================================

   {
      test_ = "Column-major CustomMatrix static array assignment";

      const int array[2][3] = { { 1, 2, 3 }, { 4, 5, 6 } };
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      OMT mat( memory.get(), 2UL, 3UL, 16UL );
      mat = array;

      checkRows    ( mat,  2UL );
      checkColumns ( mat,  3UL );
      checkCapacity( mat, 48UL );
      checkNonZeros( mat,  6UL );
      checkNonZeros( mat,  0UL, 2UL );
      checkNonZeros( mat,  1UL, 2UL );
      checkNonZeros( mat,  2UL, 2UL );

      if( mat(0,0) != 1 || mat(0,1) != 2 || mat(0,2) != 3 ||
          mat(1,0) != 4 || mat(1,1) != 5 || mat(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major CustomMatrix std::array assignment";

      const std::array<std::array<int,3UL>,2UL> array{ { { 1, 2, 3 }, { 4, 5, 6 } } };
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      OMT mat( memory.get(), 2UL, 3UL, 16UL );
      mat = array;

      checkRows    ( mat,  2UL );
      checkColumns ( mat,  3UL );
      checkCapacity( mat, 48UL );
      checkNonZeros( mat,  6UL );
      checkNonZeros( mat,  0UL, 2UL );
      checkNonZeros( mat,  1UL, 2UL );
      checkNonZeros( mat,  2UL, 2UL );

      if( mat(0,0) != 1 || mat(0,1) != 2 || mat(0,2) != 3 ||
          mat(1,0) != 4 || mat(1,1) != 5 || mat(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major copy assignment
   //=====================================================================================

   {
      test_ = "Column-major CustomMatrix copy assignment";

      std::unique_ptr<int[],blaze::Deallocate> memory1( blaze::allocate<int>( 48UL ) );
      OMT mat1( memory1.get(), 2UL, 3UL, 16UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(0,2) = 3;
      mat1(1,0) = 4;
      mat1(1,1) = 5;
      mat1(1,2) = 6;

      std::unique_ptr<int[],blaze::Deallocate> memory2( blaze::allocate<int>( 48UL ) );
      OMT mat2( memory2.get(), 2UL, 3UL, 16UL );
      mat2 = mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 48UL );
      checkNonZeros( mat2,  6UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 2UL );
      checkNonZeros( mat2,  2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 3 ||
          mat2(1,0) != 4 || mat2(1,1) != 5 || mat2(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major move assignment
   //=====================================================================================

   {
      test_ = "Column-major CustomMatrix move assignment";

      std::unique_ptr<int[],blaze::Deallocate> memory1( blaze::allocate<int>( 48UL ) );
      OMT mat1( memory1.get(), 2UL, 3UL, 16UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(0,2) = 3;
      mat1(1,0) = 4;
      mat1(1,1) = 5;
      mat1(1,2) = 6;

      std::unique_ptr<int[],blaze::Deallocate> memory2( blaze::allocate<int>( 48UL ) );
      OMT mat2( memory2.get(), 2UL, 3UL, 16UL );
      mat2 = std::move( mat1 );

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 48UL );
      checkNonZeros( mat2,  6UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 2UL );
      checkNonZeros( mat2,  2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 3 ||
          mat2(1,0) != 4 || mat2(1,1) != 5 || mat2(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major CustomMatrix dense matrix assignment (mixed type)";

      using blaze::aligned;
      using blaze::unpadded;
      using blaze::rowMajor;

      using AlignedUnpadded = blaze::CustomMatrix<short,aligned,unpadded,rowMajor>;
      std::unique_ptr<short[],blaze::Deallocate> memory1( blaze::allocate<short>( 64UL ) );
      AlignedUnpadded mat1( memory1.get(), 2UL, 3UL, 32UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(0,2) = 3;
      mat1(1,0) = 4;
      mat1(1,1) = 5;
      mat1(1,2) = 6;

      std::unique_ptr<int[],blaze::Deallocate> memory2( blaze::allocate<int>( 48UL ) );
      OMT mat2( memory2.get(), 2UL, 3UL, 16UL );
      mat2 = mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 48UL );
      checkNonZeros( mat2,  6UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 2UL );
      checkNonZeros( mat2,  2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 3 ||
          mat2(1,0) != 4 || mat2(1,1) != 5 || mat2(1,2) != 6 ) {
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
      test_ = "Column-major/row-major CustomMatrix dense matrix assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::rowMajor;

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,rowMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory1( blaze::allocate<int>( 32UL ) );
      AlignedPadded mat1( memory1.get(), 2UL, 3UL, 16UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(0,2) = 3;
      mat1(1,0) = 4;
      mat1(1,1) = 5;
      mat1(1,2) = 6;

      std::unique_ptr<int[],blaze::Deallocate> memory2( blaze::allocate<int>( 48UL ) );
      OMT mat2( memory2.get(), 2UL, 3UL, 16UL );
      mat2 = mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 48UL );
      checkNonZeros( mat2,  6UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 2UL );
      checkNonZeros( mat2,  2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 3 ||
          mat2(1,0) != 4 || mat2(1,1) != 5 || mat2(1,2) != 6 ) {
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
      test_ = "Column-major/row-major CustomMatrix dense matrix assignment stress test (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::rowMajor;

      const int min( randmin );
      const int max( randmax );

      for( size_t i=0UL; i<10UL; ++i )
      {
         const size_t rows    ( blaze::rand<size_t>( 0UL, 16UL ) );
         const size_t columns ( blaze::rand<size_t>( 0UL, 16UL ) );
         const size_t spacing1( blaze::nextMultiple<size_t>( columns, 16UL ) );
         const size_t spacing2( blaze::nextMultiple<size_t>( rows, 16UL ) );

         using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,rowMajor>;
         std::unique_ptr<int[],blaze::Deallocate> memory1( blaze::allocate<int>( rows*spacing1 ) );
         AlignedPadded mat1( memory1.get(), rows, columns, spacing1 );
         randomize( mat1, min, max );

         std::unique_ptr<int[],blaze::Deallocate> memory2( blaze::allocate<int>( spacing2*columns ) );
         OMT mat2( memory2.get(), rows, columns, spacing2 );
         mat2 = mat1;

         if( mat1 != mat2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << mat2 << "\n"
                << "   Expected result:\n" << mat1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   {
      test_ = "Column-major/row-major CustomMatrix dense matrix assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowMajor;

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>;
      std::unique_ptr<int[]> memory1( new int[7UL] );
      UnalignedUnpadded mat1( memory1.get()+1UL, 2UL, 3UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(0,2) = 3;
      mat1(1,0) = 4;
      mat1(1,1) = 5;
      mat1(1,2) = 6;

      std::unique_ptr<int[],blaze::Deallocate> memory2( blaze::allocate<int>( 48UL ) );
      OMT mat2( memory2.get(), 2UL, 3UL, 16UL );
      mat2 = mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 48UL );
      checkNonZeros( mat2,  6UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 2UL );
      checkNonZeros( mat2,  2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 3 ||
          mat2(1,0) != 4 || mat2(1,1) != 5 || mat2(1,2) != 6 ) {
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
      test_ = "Column-major/row-major CustomMatrix dense matrix assignment stress test (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowMajor;

      const int min( randmin );
      const int max( randmax );

      for( size_t i=0UL; i<10UL; ++i )
      {
         const size_t rows   ( blaze::rand<size_t>( 0UL, 16UL ) );
         const size_t columns( blaze::rand<size_t>( 0UL, 16UL ) );
         const size_t spacing( blaze::nextMultiple<size_t>( rows, 16UL ) );

         using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>;
         std::unique_ptr<int[]> memory1( new int[rows*columns+1UL] );
         UnalignedUnpadded mat1( memory1.get()+1UL, rows, columns );
         randomize( mat1, min, max );

         std::unique_ptr<int[],blaze::Deallocate> memory2( blaze::allocate<int>( spacing*columns ) );
         OMT mat2( memory2.get(), rows, columns, spacing );
         mat2 = mat1;

         if( mat1 != mat2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << mat2 << "\n"
                << "   Expected result:\n" << mat1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   {
      test_ = "Column-major/column-major CustomMatrix dense matrix assignment (mixed type)";

      using blaze::aligned;
      using blaze::unpadded;
      using blaze::columnMajor;

      using AlignedUnpadded = blaze::CustomMatrix<short,aligned,unpadded,columnMajor>;
      std::unique_ptr<short[],blaze::Deallocate> memory1( blaze::allocate<short>( 96UL ) );
      AlignedUnpadded mat1( memory1.get(), 2UL, 3UL, 32UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(0,2) = 3;
      mat1(1,0) = 4;
      mat1(1,1) = 5;
      mat1(1,2) = 6;

      std::unique_ptr<int[],blaze::Deallocate> memory2( blaze::allocate<int>( 48UL ) );
      OMT mat2( memory2.get(), 2UL, 3UL, 16UL );
      mat2 = mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 48UL );
      checkNonZeros( mat2,  6UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 2UL );
      checkNonZeros( mat2,  2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 3 ||
          mat2(1,0) != 4 || mat2(1,1) != 5 || mat2(1,2) != 6 ) {
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
      test_ = "Column-major/column-major CustomMatrix dense matrix assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::columnMajor;

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,columnMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory1( blaze::allocate<int>( 48UL ) );
      AlignedPadded mat1( memory1.get(), 2UL, 3UL, 16UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(0,2) = 3;
      mat1(1,0) = 4;
      mat1(1,1) = 5;
      mat1(1,2) = 6;

      std::unique_ptr<int[],blaze::Deallocate> memory2( blaze::allocate<int>( 48UL ) );
      OMT mat2( memory2.get(), 2UL, 3UL, 16UL );
      mat2 = mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 48UL );
      checkNonZeros( mat2,  6UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 2UL );
      checkNonZeros( mat2,  2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 3 ||
          mat2(1,0) != 4 || mat2(1,1) != 5 || mat2(1,2) != 6 ) {
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
      test_ = "Column-major/column-major CustomMatrix dense matrix assignment stress test (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::columnMajor;

      const int min( randmin );
      const int max( randmax );

      for( size_t i=0UL; i<10UL; ++i )
      {
         const size_t rows   ( blaze::rand<size_t>( 0UL, 16UL ) );
         const size_t columns( blaze::rand<size_t>( 0UL, 16UL ) );
         const size_t spacing( blaze::nextMultiple<size_t>( rows, 16UL ) );

         using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,columnMajor>;
         std::unique_ptr<int[],blaze::Deallocate> memory1( blaze::allocate<int>( spacing*columns ) );
         AlignedPadded mat1( memory1.get(), rows, columns, spacing );
         randomize( mat1, min, max );

         std::unique_ptr<int[],blaze::Deallocate> memory2( blaze::allocate<int>( spacing*columns ) );
         OMT mat2( memory2.get(), rows, columns, spacing );
         mat2 = mat1;

         if( mat1 != mat2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << mat2 << "\n"
                << "   Expected result:\n" << mat1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   {
      test_ = "Column-major/column-major CustomMatrix dense matrix assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::columnMajor;

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>;
      std::unique_ptr<int[]> memory1( new int[7UL] );
      UnalignedUnpadded mat1( memory1.get()+1UL, 2UL, 3UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(0,2) = 3;
      mat1(1,0) = 4;
      mat1(1,1) = 5;
      mat1(1,2) = 6;

      std::unique_ptr<int[],blaze::Deallocate> memory2( blaze::allocate<int>( 48UL ) );
      OMT mat2( memory2.get(), 2UL, 3UL, 16UL );
      mat2 = mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 48UL );
      checkNonZeros( mat2,  6UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 2UL );
      checkNonZeros( mat2,  2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 3 ||
          mat2(1,0) != 4 || mat2(1,1) != 5 || mat2(1,2) != 6 ) {
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
      test_ = "Column-major/column-major CustomMatrix dense matrix assignment stress test (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::columnMajor;

      const int min( randmin );
      const int max( randmax );

      for( size_t i=0UL; i<10UL; ++i )
      {
         const size_t rows   ( blaze::rand<size_t>( 0UL, 16UL ) );
         const size_t columns( blaze::rand<size_t>( 0UL, 16UL ) );
         const size_t spacing( blaze::nextMultiple<size_t>( rows, 16UL ) );

         using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>;
         std::unique_ptr<int[]> memory1( new int[rows*columns+1UL] );
         UnalignedUnpadded mat1( memory1.get()+1UL, rows, columns );
         randomize( mat1, min, max );

         std::unique_ptr<int[],blaze::Deallocate> memory2( blaze::allocate<int>( spacing*columns ) );
         OMT mat2( memory2.get(), rows, columns, spacing );
         mat2 = mat1;

         if( mat1 != mat2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << mat2 << "\n"
                << "   Expected result:\n" << mat1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   {
      test_ = "Column-major/row-major CustomMatrix dense matrix assignment (lower)";

      blaze::LowerMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      OMT mat2( memory.get(), 3UL, 3UL, 16UL );
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
      test_ = "Column-major/column-major CustomMatrix dense matrix assignment (lower)";

      blaze::LowerMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      OMT mat2( memory.get(), 3UL, 3UL, 16UL );
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
      test_ = "Column-major/row-major CustomMatrix dense matrix assignment (upper)";

      blaze::UpperMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      OMT mat2( memory.get(), 3UL, 3UL, 16UL );
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
      test_ = "Column-major/column-major CustomMatrix dense matrix assignment (upper)";

      blaze::UpperMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      OMT mat2( memory.get(), 3UL, 3UL, 16UL );
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
      test_ = "Column-major/row-major CustomMatrix dense matrix assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      OMT mat2( memory.get(), 3UL, 3UL, 16UL );
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
      test_ = "Column-major/column-major CustomMatrix dense matrix assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      OMT mat2( memory.get(), 3UL, 3UL, 16UL );
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
      test_ = "Column-major/row-major CustomMatrix sparse matrix assignment";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1( 2UL, 3UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(1,0) = 3;
      mat1(1,2) = 4;

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      OMT mat2( memory.get(), 2UL, 3UL, 16UL );
      mat2 = mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 48UL );
      checkNonZeros( mat2,  4UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 1UL );
      checkNonZeros( mat2,  2UL, 1UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 0 ||
          mat2(1,0) != 3 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 0 )\n( 3 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major CustomMatrix sparse matrix assignment stress test";

      const int min( randmin );
      const int max( randmax );

      for( size_t i=0UL; i<10UL; ++i )
      {
         const size_t rows   ( blaze::rand<size_t>( 0UL, 16UL ) );
         const size_t columns( blaze::rand<size_t>( 0UL, 16UL ) );
         const size_t spacing( blaze::nextMultiple<size_t>( rows, 16UL ) );

         blaze::CompressedMatrix<int,blaze::rowMajor> mat1( rows, columns );
         randomize( mat1, min, max );

         std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( spacing*columns ) );
         OMT mat2( memory.get(), rows, columns, spacing );
         mat2 = mat1;

         if( mat1 != mat2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << mat2 << "\n"
                << "   Expected result:\n" << mat1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   {
      test_ = "Column-major/column-major CustomMatrix sparse matrix assignment";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1( 2UL, 3UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(1,0) = 3;
      mat1(1,2) = 4;

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      OMT mat2( memory.get(), 2UL, 3UL, 16UL );
      mat2 = mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 48UL );
      checkNonZeros( mat2,  4UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 1UL );
      checkNonZeros( mat2,  2UL, 1UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 0 ||
          mat2(1,0) != 3 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 0 )\n( 3 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major CustomMatrix sparse matrix assignment stress test";

      const int min( randmin );
      const int max( randmax );

      for( size_t i=0UL; i<10UL; ++i )
      {
         const size_t rows   ( blaze::rand<size_t>( 0UL, 16UL ) );
         const size_t columns( blaze::rand<size_t>( 0UL, 16UL ) );
         const size_t spacing( blaze::nextMultiple<size_t>( rows, 16UL ) );

         blaze::CompressedMatrix<int,blaze::columnMajor> mat1( rows, columns );
         randomize( mat1, min, max );

         std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( spacing*columns ) );
         OMT mat2( memory.get(), rows, columns, spacing );
         mat2 = mat1;

         if( mat1 != mat2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << mat2 << "\n"
                << "   Expected result:\n" << mat1 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   {
      test_ = "Column-major/row-major CustomMatrix sparse matrix assignment (lower)";

      blaze::LowerMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      OMT mat2( memory.get(), 3UL, 3UL, 16UL );
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
      test_ = "Column-major/column-major CustomMatrix sparse matrix assignment (lower)";

      blaze::LowerMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      OMT mat2( memory.get(), 3UL, 3UL, 16UL );
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
      test_ = "Column-major/row-major CustomMatrix sparse matrix assignment (upper)";

      blaze::UpperMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      OMT mat2( memory.get(), 3UL, 3UL, 16UL );
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
      test_ = "Column-major/column-major CustomMatrix sparse matrix assignment (upper)";

      blaze::UpperMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      OMT mat2( memory.get(), 3UL, 3UL, 16UL );
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
      test_ = "Column-major/row-major CustomMatrix sparse matrix assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      OMT mat2( memory.get(), 3UL, 3UL, 16UL );
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
      test_ = "Column-major/column-major CustomMatrix sparse matrix assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      OMT mat2( memory.get(), 3UL, 3UL, 16UL );
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
/*!\brief Test of the CustomMatrix addition assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the addition assignment operators of the CustomMatrix class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void AlignedUnpaddedTest::testAddAssign()
{
   //=====================================================================================
   // Row-major dense matrix addition assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major CustomMatrix dense matrix addition assignment (mixed type)";

      using blaze::aligned;
      using blaze::unpadded;
      using blaze::rowMajor;

      using AlignedUnpadded = blaze::CustomMatrix<short,aligned,unpadded,rowMajor>;
      std::unique_ptr<short[],blaze::Deallocate> memory1( blaze::allocate<short>( 64UL ) );
      AlignedUnpadded mat1( memory1.get(), 2UL, 3UL, 32UL );
      mat1 = 0;
      mat1(0,0) =  1;
      mat1(0,1) =  2;
      mat1(1,0) = -3;
      mat1(1,2) =  4;

      std::unique_ptr<int[],blaze::Deallocate> memory2( blaze::allocate<int>( 32UL ) );
      MT mat2( memory2.get(), 2UL, 3UL, 16UL );
      mat2 = 0;
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 += mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 32UL );
      checkNonZeros( mat2,  4UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major CustomMatrix dense matrix addition assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::rowMajor;

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,rowMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory1( blaze::allocate<int>( 32UL ) );
      AlignedPadded mat1( memory1.get(), 2UL, 3UL, 16UL );
      mat1 = 0;
      mat1(0,0) =  1;
      mat1(0,1) =  2;
      mat1(1,0) = -3;
      mat1(1,2) =  4;

      std::unique_ptr<int[],blaze::Deallocate> memory2( blaze::allocate<int>( 32UL ) );
      MT mat2( memory2.get(), 2UL, 3UL, 16UL );
      mat2 = 0;
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 += mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 32UL );
      checkNonZeros( mat2,  4UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major CustomMatrix dense matrix addition assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowMajor;

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>;
      std::unique_ptr<int[]> memory1( new int[7UL] );
      UnalignedUnpadded mat1( memory1.get()+1UL, 2UL, 3UL );
      mat1 = 0;
      mat1(0,0) =  1;
      mat1(0,1) =  2;
      mat1(1,0) = -3;
      mat1(1,2) =  4;

      std::unique_ptr<int[],blaze::Deallocate> memory2( blaze::allocate<int>( 32UL ) );
      MT mat2( memory2.get(), 2UL, 3UL, 16UL );
      mat2 = 0;
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 += mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 32UL );
      checkNonZeros( mat2,  4UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CustomMatrix dense matrix addition assignment (mixed)";

      using blaze::aligned;
      using blaze::unpadded;
      using blaze::columnMajor;

      using AlignedUnpadded = blaze::CustomMatrix<short,aligned,unpadded,columnMajor>;
      std::unique_ptr<short[],blaze::Deallocate> memory1( blaze::allocate<short>( 96UL ) );
      AlignedUnpadded mat1( memory1.get(), 2UL, 3UL, 32UL );
      mat1 = 0;
      mat1(0,0) =  1;
      mat1(0,1) =  2;
      mat1(1,0) = -3;
      mat1(1,2) =  4;

      std::unique_ptr<int[],blaze::Deallocate> memory2( blaze::allocate<int>( 32UL ) );
      MT mat2( memory2.get(), 2UL, 3UL, 16UL );
      mat2 = 0;
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 += mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 32UL );
      checkNonZeros( mat2,  4UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CustomMatrix dense matrix addition assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::columnMajor;

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,columnMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory1( blaze::allocate<int>( 48UL ) );
      AlignedPadded mat1( memory1.get(), 2UL, 3UL, 16UL );
      mat1 = 0;
      mat1(0,0) =  1;
      mat1(0,1) =  2;
      mat1(1,0) = -3;
      mat1(1,2) =  4;

      std::unique_ptr<int[],blaze::Deallocate> memory2( blaze::allocate<int>( 32UL ) );
      MT mat2( memory2.get(), 2UL, 3UL, 16UL );
      mat2 = 0;
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 += mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 32UL );
      checkNonZeros( mat2,  4UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CustomMatrix dense matrix addition assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::columnMajor;

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>;
      std::unique_ptr<int[]> memory1( new int[7UL] );
      UnalignedUnpadded mat1( memory1.get()+1UL, 2UL, 3UL );
      mat1 = 0;
      mat1(0,0) =  1;
      mat1(0,1) =  2;
      mat1(1,0) = -3;
      mat1(1,2) =  4;

      std::unique_ptr<int[],blaze::Deallocate> memory2( blaze::allocate<int>( 32UL ) );
      MT mat2( memory2.get(), 2UL, 3UL, 16UL );
      mat2 = 0;
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 += mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 32UL );
      checkNonZeros( mat2,  4UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major CustomMatrix dense matrix addition assignment (lower)";

      blaze::LowerMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      MT mat2( memory.get(), 3UL, 3UL, 16UL );
      mat2 = 0;

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
      test_ = "Row-major/column-major CustomMatrix dense matrix addition assignment (lower)";

      blaze::LowerMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      MT mat2( memory.get(), 3UL, 3UL, 16UL );
      mat2 = 0;

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
      test_ = "Row-major/row-major CustomMatrix dense matrix addition assignment (upper)";

      blaze::UpperMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      MT mat2( memory.get(), 3UL, 3UL, 16UL );
      mat2 = 0;

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
      test_ = "Row-major/column-major CustomMatrix dense matrix addition assignment (upper)";

      blaze::UpperMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      MT mat2( memory.get(), 3UL, 3UL, 16UL );
      mat2 = 0;

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
      test_ = "Row-major/row-major CustomMatrix dense matrix addition assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      MT mat2( memory.get(), 3UL, 3UL, 16UL );
      mat2 = 0;

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
      test_ = "Row-major/column-major CustomMatrix dense matrix addition assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      MT mat2( memory.get(), 3UL, 3UL, 16UL );
      mat2 = 0;

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
      test_ = "Row-major/row-major CustomMatrix sparse matrix addition assignment";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1( 2UL, 3UL, 4UL );
      mat1(0,0) =  1;
      mat1(0,1) =  2;
      mat1(1,0) = -3;
      mat1(1,2) =  4;

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 32UL ) );
      MT mat2( memory.get(), 2UL, 3UL, 16UL );
      mat2 = 0;
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 += mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 32UL );
      checkNonZeros( mat2,  4UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CustomMatrix sparse matrix addition assignment";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1( 2UL, 3UL, 4UL );
      mat1(0,0) =  1;
      mat1(0,1) =  2;
      mat1(1,0) = -3;
      mat1(1,2) =  4;

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 32UL ) );
      MT mat2( memory.get(), 2UL, 3UL, 16UL );
      mat2 = 0;
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 += mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 32UL );
      checkNonZeros( mat2,  4UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major CustomMatrix sparse matrix addition assignment (lower)";

      blaze::LowerMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      MT mat2( memory.get(), 3UL, 3UL, 16UL );
      mat2 = 0;

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
      test_ = "Row-major/column-major CustomMatrix sparse matrix addition assignment (lower)";

      blaze::LowerMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      MT mat2( memory.get(), 3UL, 3UL, 16UL );
      mat2 = 0;

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
      test_ = "Row-major/row-major CustomMatrix sparse matrix addition assignment (upper)";

      blaze::UpperMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      MT mat2( memory.get(), 3UL, 3UL, 16UL );
      mat2 = 0;

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
      test_ = "Row-major/column-major CustomMatrix sparse matrix addition assignment (upper)";

      blaze::UpperMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      MT mat2( memory.get(), 3UL, 3UL, 16UL );
      mat2 = 0;

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
      test_ = "Row-major/row-major CustomMatrix sparse matrix addition assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      MT mat2( memory.get(), 3UL, 3UL, 16UL );
      mat2 = 0;

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
      test_ = "Row-major/column-major CustomMatrix sparse matrix addition assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      MT mat2( memory.get(), 3UL, 3UL, 16UL );
      mat2 = 0;

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
      test_ = "Column-major/row-major CustomMatrix dense matrix addition assignment (mixed type)";

      using blaze::aligned;
      using blaze::unpadded;
      using blaze::rowMajor;

      using AlignedUnpadded = blaze::CustomMatrix<short,aligned,unpadded,rowMajor>;
      std::unique_ptr<short[],blaze::Deallocate> memory1( blaze::allocate<short>( 64UL ) );
      AlignedUnpadded mat1( memory1.get(), 2UL, 3UL, 32UL );
      mat1 = 0;
      mat1(0,0) =  1;
      mat1(0,1) =  2;
      mat1(1,0) = -3;
      mat1(1,2) =  4;

      std::unique_ptr<int[],blaze::Deallocate> memory2( blaze::allocate<int>( 48UL ) );
      OMT mat2( memory2.get(), 2UL, 3UL, 16UL );
      mat2 = 0;
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 += mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 48UL );
      checkNonZeros( mat2,  4UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 0UL );
      checkNonZeros( mat2,  2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major CustomMatrix dense matrix addition assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::rowMajor;

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,rowMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory1( blaze::allocate<int>( 32UL ) );
      AlignedPadded mat1( memory1.get(), 2UL, 3UL, 16UL );
      mat1 = 0;
      mat1(0,0) =  1;
      mat1(0,1) =  2;
      mat1(1,0) = -3;
      mat1(1,2) =  4;

      std::unique_ptr<int[],blaze::Deallocate> memory2( blaze::allocate<int>( 48UL ) );
      OMT mat2( memory2.get(), 2UL, 3UL, 16UL );
      mat2 = 0;
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 += mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 48UL );
      checkNonZeros( mat2,  4UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 0UL );
      checkNonZeros( mat2,  2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major CustomMatrix dense matrix addition assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowMajor;

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>;
      std::unique_ptr<int[]> memory1( new int[7UL] );
      UnalignedUnpadded mat1( memory1.get()+1UL, 2UL, 3UL );
      mat1 = 0;
      mat1(0,0) =  1;
      mat1(0,1) =  2;
      mat1(1,0) = -3;
      mat1(1,2) =  4;

      std::unique_ptr<int[],blaze::Deallocate> memory2( blaze::allocate<int>( 48UL ) );
      OMT mat2( memory2.get(), 2UL, 3UL, 16UL );
      mat2 = 0;
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 += mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 48UL );
      checkNonZeros( mat2,  4UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 0UL );
      checkNonZeros( mat2,  2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major CustomMatrix dense matrix addition assignment (mixed type)";

      using blaze::aligned;
      using blaze::unpadded;
      using blaze::columnMajor;

      using AlignedUnpadded = blaze::CustomMatrix<short,aligned,unpadded,columnMajor>;
      std::unique_ptr<short[],blaze::Deallocate> memory1( blaze::allocate<short>( 96UL ) );
      AlignedUnpadded mat1( memory1.get(), 2UL, 3UL, 32UL );
      mat1 = 0;
      mat1(0,0) =  1;
      mat1(0,1) =  2;
      mat1(1,0) = -3;
      mat1(1,2) =  4;

      std::unique_ptr<int[],blaze::Deallocate> memory2( blaze::allocate<int>( 48UL ) );
      OMT mat2( memory2.get(), 2UL, 3UL, 16UL );
      mat2 = 0;
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 += mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 48UL );
      checkNonZeros( mat2,  4UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 0UL );
      checkNonZeros( mat2,  2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major CustomMatrix dense matrix addition assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::columnMajor;

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,columnMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory1( blaze::allocate<int>( 48UL ) );
      AlignedPadded mat1( memory1.get(), 2UL, 3UL, 16UL );
      mat1 = 0;
      mat1(0,0) =  1;
      mat1(0,1) =  2;
      mat1(1,0) = -3;
      mat1(1,2) =  4;

      std::unique_ptr<int[],blaze::Deallocate> memory2( blaze::allocate<int>( 48UL ) );
      OMT mat2( memory2.get(), 2UL, 3UL, 16UL );
      mat2 = 0;
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 += mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 48UL );
      checkNonZeros( mat2,  4UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 0UL );
      checkNonZeros( mat2,  2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major CustomMatrix dense matrix addition assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::columnMajor;

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>;
      std::unique_ptr<int[]> memory1( new int[7UL] );
      UnalignedUnpadded mat1( memory1.get()+1UL, 2UL, 3UL );
      mat1 = 0;
      mat1(0,0) =  1;
      mat1(0,1) =  2;
      mat1(1,0) = -3;
      mat1(1,2) =  4;

      std::unique_ptr<int[],blaze::Deallocate> memory2( blaze::allocate<int>( 48UL ) );
      OMT mat2( memory2.get(), 2UL, 3UL, 16UL );
      mat2 = 0;
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 += mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 48UL );
      checkNonZeros( mat2,  4UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 0UL );
      checkNonZeros( mat2,  2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major CustomMatrix dense matrix addition assignment (lower)";

      blaze::LowerMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      OMT mat2( memory.get(), 3UL, 3UL, 16UL );
      mat2 = 0;

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
      test_ = "Column-major/column-major CustomMatrix dense matrix addition assignment (lower)";

      blaze::LowerMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      OMT mat2( memory.get(), 3UL, 3UL, 16UL );
      mat2 = 0;

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
      test_ = "Column-major/row-major CustomMatrix dense matrix addition assignment (upper)";

      blaze::UpperMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      OMT mat2( memory.get(), 3UL, 3UL, 16UL );
      mat2 = 0;

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
      test_ = "Column-major/column-major CustomMatrix dense matrix addition assignment (upper)";

      blaze::UpperMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      OMT mat2( memory.get(), 3UL, 3UL, 16UL );
      mat2 = 0;

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
      test_ = "Column-major/row-major CustomMatrix dense matrix addition assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      OMT mat2( memory.get(), 3UL, 3UL, 16UL );
      mat2 = 0;

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
      test_ = "Column-major/column-major CustomMatrix dense matrix addition assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      OMT mat2( memory.get(), 3UL, 3UL, 16UL );
      mat2 = 0;

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
      test_ = "Column-major/row-major CustomMatrix sparse matrix addition assignment";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1( 2UL, 3UL, 4UL );
      mat1(0,0) =  1;
      mat1(0,1) =  2;
      mat1(1,0) = -3;
      mat1(1,2) =  4;

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      OMT mat2( memory.get(), 2UL, 3UL, 16UL );
      mat2 = 0;
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 += mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 48UL );
      checkNonZeros( mat2,  4UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 0UL );
      checkNonZeros( mat2,  2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major CustomMatrix sparse matrix addition assignment";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1( 2UL, 3UL, 4UL );
      mat1(0,0) =  1;
      mat1(0,1) =  2;
      mat1(1,0) = -3;
      mat1(1,2) =  4;

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      OMT mat2( memory.get(), 2UL, 3UL, 16UL );
      mat2 = 0;
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 += mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 48UL );
      checkNonZeros( mat2,  4UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 0UL );
      checkNonZeros( mat2,  2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major CustomMatrix sparse matrix addition assignment (lower)";

      blaze::LowerMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      OMT mat2( memory.get(), 3UL, 3UL, 16UL );
      mat2 = 0;

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
      test_ = "Row-major/column-major CustomMatrix sparse matrix addition assignment (lower)";

      blaze::LowerMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      OMT mat2( memory.get(), 3UL, 3UL, 16UL );
      mat2 = 0;

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
      test_ = "Row-major/row-major CustomMatrix sparse matrix addition assignment (upper)";

      blaze::UpperMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      OMT mat2( memory.get(), 3UL, 3UL, 16UL );
      mat2 = 0;

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
      test_ = "Row-major/column-major CustomMatrix sparse matrix addition assignment (upper)";

      blaze::UpperMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      OMT mat2( memory.get(), 3UL, 3UL, 16UL );
      mat2 = 0;

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
      test_ = "Row-major/row-major CustomMatrix sparse matrix addition assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      OMT mat2( memory.get(), 3UL, 3UL, 16UL );
      mat2 = 0;

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
      test_ = "Row-major/column-major CustomMatrix sparse matrix addition assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      OMT mat2( memory.get(), 3UL, 3UL, 16UL );
      mat2 = 0;

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
/*!\brief Test of the CustomMatrix subtraction assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the subtraction assignment operators of the CustomMatrix
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void AlignedUnpaddedTest::testSubAssign()
{
   //=====================================================================================
   // Row-major dense matrix subtraction assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major CustomMatrix dense matrix subtraction assignment (mixed type)";

      using blaze::aligned;
      using blaze::unpadded;
      using blaze::rowMajor;

      using AlignedUnpadded = blaze::CustomMatrix<short,aligned,unpadded,rowMajor>;
      std::unique_ptr<short[],blaze::Deallocate> memory1( blaze::allocate<short>( 64UL ) );
      AlignedUnpadded mat1( memory1.get(), 2UL, 3UL, 32UL );
      mat1 = 0;
      mat1(0,0) = -1;
      mat1(0,1) = -2;
      mat1(1,0) =  3;
      mat1(1,2) = -4;

      std::unique_ptr<int[],blaze::Deallocate> memory2( blaze::allocate<int>( 32UL ) );
      MT mat2( memory2.get(), 2UL, 3UL, 16UL );
      mat2 = 0;
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 -= mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 32UL );
      checkNonZeros( mat2,  4UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major CustomMatrix dense matrix subtraction assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::rowMajor;

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,rowMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory1( blaze::allocate<int>( 32UL ) );
      AlignedPadded mat1( memory1.get(), 2UL, 3UL, 16UL );
      mat1 = 0;
      mat1(0,0) = -1;
      mat1(0,1) = -2;
      mat1(1,0) =  3;
      mat1(1,2) = -4;

      std::unique_ptr<int[],blaze::Deallocate> memory2( blaze::allocate<int>( 32UL ) );
      MT mat2( memory2.get(), 2UL, 3UL, 16UL );
      mat2 = 0;
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 -= mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 32UL );
      checkNonZeros( mat2,  4UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major CustomMatrix dense matrix subtraction assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowMajor;

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>;
      std::unique_ptr<int[]> memory1( new int[7UL] );
      UnalignedUnpadded mat1( memory1.get()+1UL, 2UL, 3UL );
      mat1 = 0;
      mat1(0,0) = -1;
      mat1(0,1) = -2;
      mat1(1,0) =  3;
      mat1(1,2) = -4;

      std::unique_ptr<int[],blaze::Deallocate> memory2( blaze::allocate<int>( 32UL ) );
      MT mat2( memory2.get(), 2UL, 3UL, 16UL );
      mat2 = 0;
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 -= mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 32UL );
      checkNonZeros( mat2,  4UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CustomMatrix dense matrix subtraction assignment (mixed type)";

      using blaze::aligned;
      using blaze::unpadded;
      using blaze::columnMajor;

      using AlignedUnpadded = blaze::CustomMatrix<short,aligned,unpadded,columnMajor>;
      std::unique_ptr<short[],blaze::Deallocate> memory1( blaze::allocate<short>( 96UL ) );
      AlignedUnpadded mat1( memory1.get(), 2UL, 3UL, 32UL );
      mat1 = 0;
      mat1(0,0) = -1;
      mat1(0,1) = -2;
      mat1(1,0) =  3;
      mat1(1,2) = -4;

      std::unique_ptr<int[],blaze::Deallocate> memory2( blaze::allocate<int>( 32UL ) );
      MT mat2( memory2.get(), 2UL, 3UL, 16UL );
      mat2 = 0;
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 -= mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 32UL );
      checkNonZeros( mat2,  4UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CustomMatrix dense matrix subtraction assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::columnMajor;

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,columnMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory1( blaze::allocate<int>( 48UL ) );
      AlignedPadded mat1( memory1.get(), 2UL, 3UL, 16UL );
      mat1 = 0;
      mat1(0,0) = -1;
      mat1(0,1) = -2;
      mat1(1,0) =  3;
      mat1(1,2) = -4;

      std::unique_ptr<int[],blaze::Deallocate> memory2( blaze::allocate<int>( 32UL ) );
      MT mat2( memory2.get(), 2UL, 3UL, 16UL );
      mat2 = 0;
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 -= mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 32UL );
      checkNonZeros( mat2,  4UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CustomMatrix dense matrix subtraction assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::columnMajor;

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>;
      std::unique_ptr<int[]> memory1( new int[7UL] );
      UnalignedUnpadded mat1( memory1.get()+1UL, 2UL, 3UL );
      mat1 = 0;
      mat1(0,0) = -1;
      mat1(0,1) = -2;
      mat1(1,0) =  3;
      mat1(1,2) = -4;

      std::unique_ptr<int[],blaze::Deallocate> memory2( blaze::allocate<int>( 32UL ) );
      MT mat2( memory2.get(), 2UL, 3UL, 16UL );
      mat2 = 0;
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 -= mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 32UL );
      checkNonZeros( mat2,  4UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major CustomMatrix dense matrix subtraction assignment (lower)";

      blaze::LowerMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      MT mat2( memory.get(), 3UL, 3UL, 16UL );
      mat2 = 0;

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
      test_ = "Row-major/column-major CustomMatrix dense matrix subtraction assignment (lower)";

      blaze::LowerMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      MT mat2( memory.get(), 3UL, 3UL, 16UL );
      mat2 = 0;

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
      test_ = "Row-major/row-major CustomMatrix dense matrix subtraction assignment (upper)";

      blaze::UpperMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      MT mat2( memory.get(), 3UL, 3UL, 16UL );
      mat2 = 0;

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
      test_ = "Row-major/column-major CustomMatrix dense matrix subtraction assignment (upper)";

      blaze::UpperMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      MT mat2( memory.get(), 3UL, 3UL, 16UL );
      mat2 = 0;

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
      test_ = "Row-major/row-major CustomMatrix dense matrix subtraction assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      MT mat2( memory.get(), 3UL, 3UL, 16UL );
      mat2 = 0;

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
      test_ = "Row-major/column-major CustomMatrix dense matrix subtraction assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      MT mat2( memory.get(), 3UL, 3UL, 16UL );
      mat2 = 0;

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
      test_ = "Row-major/row-major CustomMatrix sparse matrix subtraction assignment";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1( 2UL, 3UL, 4UL );
      mat1(0,0) = -1;
      mat1(0,1) = -2;
      mat1(1,0) =  3;
      mat1(1,2) = -4;

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 32UL ) );
      MT mat2( memory.get(), 2UL, 3UL, 16UL );
      mat2 = 0;
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 -= mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 32UL );
      checkNonZeros( mat2,  4UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major CustomMatrix sparse matrix subtraction assignment";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1( 2UL, 3UL, 4UL );
      mat1(0,0) = -1;
      mat1(0,1) = -2;
      mat1(1,0) =  3;
      mat1(1,2) = -4;

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 32UL ) );
      MT mat2( memory.get(), 2UL, 3UL, 16UL );
      mat2 = 0;
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 -= mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 32UL );
      checkNonZeros( mat2,  4UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major CustomMatrix sparse matrix subtraction assignment (lower)";

      blaze::LowerMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      MT mat2( memory.get(), 3UL, 3UL, 16UL );
      mat2 = 0;

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
      test_ = "Row-major/column-major CustomMatrix sparse matrix subtraction assignment (lower)";

      blaze::LowerMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      MT mat2( memory.get(), 3UL, 3UL, 16UL );
      mat2 = 0;

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
      test_ = "Row-major/row-major CustomMatrix sparse matrix subtraction assignment (upper)";

      blaze::UpperMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      MT mat2( memory.get(), 3UL, 3UL, 16UL );
      mat2 = 0;

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
      test_ = "Row-major/column-major CustomMatrix sparse matrix subtraction assignment (upper)";

      blaze::UpperMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      MT mat2( memory.get(), 3UL, 3UL, 16UL );
      mat2 = 0;

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
      test_ = "Row-major/row-major CustomMatrix sparse matrix subtraction assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      MT mat2( memory.get(), 3UL, 3UL, 16UL );
      mat2 = 0;

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
      test_ = "Row-major/column-major CustomMatrix sparse matrix subtraction assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      MT mat2( memory.get(), 3UL, 3UL, 16UL );
      mat2 = 0;

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
      test_ = "Column-major/row-major CustomMatrix dense matrix subtraction assignment (mixed type)";

      using blaze::aligned;
      using blaze::unpadded;
      using blaze::rowMajor;

      using AlignedUnpadded = blaze::CustomMatrix<short,aligned,unpadded,rowMajor>;
      std::unique_ptr<short[],blaze::Deallocate> memory1( blaze::allocate<short>( 64UL ) );
      AlignedUnpadded mat1( memory1.get(), 2UL, 3UL, 32UL );
      mat1 = 0;
      mat1(0,0) = -1;
      mat1(0,1) = -2;
      mat1(1,0) =  3;
      mat1(1,2) = -4;

      std::unique_ptr<int[],blaze::Deallocate> memory2( blaze::allocate<int>( 48UL ) );
      OMT mat2( memory2.get(), 2UL, 3UL, 16UL );
      mat2 = 0;
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 -= mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 48UL );
      checkNonZeros( mat2,  4UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 0UL );
      checkNonZeros( mat2,  2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major CustomMatrix dense matrix subtraction assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::rowMajor;

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,rowMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory1( blaze::allocate<int>( 32UL ) );
      AlignedPadded mat1( memory1.get(), 2UL, 3UL, 16UL );
      mat1 = 0;
      mat1(0,0) = -1;
      mat1(0,1) = -2;
      mat1(1,0) =  3;
      mat1(1,2) = -4;

      std::unique_ptr<int[],blaze::Deallocate> memory2( blaze::allocate<int>( 48UL ) );
      OMT mat2( memory2.get(), 2UL, 3UL, 16UL );
      mat2 = 0;
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 -= mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 48UL );
      checkNonZeros( mat2,  4UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 0UL );
      checkNonZeros( mat2,  2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major CustomMatrix dense matrix subtraction assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowMajor;

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>;
      std::unique_ptr<int[]> memory1( new int[7UL] );
      UnalignedUnpadded mat1( memory1.get()+1UL, 2UL, 3UL );
      mat1 = 0;
      mat1(0,0) = -1;
      mat1(0,1) = -2;
      mat1(1,0) =  3;
      mat1(1,2) = -4;

      std::unique_ptr<int[],blaze::Deallocate> memory2( blaze::allocate<int>( 48UL ) );
      OMT mat2( memory2.get(), 2UL, 3UL, 16UL );
      mat2 = 0;
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 -= mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 48UL );
      checkNonZeros( mat2,  4UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 0UL );
      checkNonZeros( mat2,  2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major CustomMatrix dense matrix subtraction assignment (mixed type)";

      using blaze::aligned;
      using blaze::unpadded;
      using blaze::columnMajor;

      using AlignedUnpadded = blaze::CustomMatrix<short,aligned,unpadded,columnMajor>;
      std::unique_ptr<short[],blaze::Deallocate> memory1( blaze::allocate<short>( 96UL ) );
      AlignedUnpadded mat1( memory1.get(), 2UL, 3UL, 32UL );
      mat1 = 0;
      mat1(0,0) = -1;
      mat1(0,1) = -2;
      mat1(1,0) =  3;
      mat1(1,2) = -4;

      std::unique_ptr<int[],blaze::Deallocate> memory2( blaze::allocate<int>( 48UL ) );
      OMT mat2( memory2.get(), 2UL, 3UL, 16UL );
      mat2 = 0;
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 -= mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 48UL );
      checkNonZeros( mat2,  4UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 0UL );
      checkNonZeros( mat2,  2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major CustomMatrix dense matrix subtraction assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::columnMajor;

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,columnMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory1( blaze::allocate<int>( 48UL ) );
      AlignedPadded mat1( memory1.get(), 2UL, 3UL, 16UL );
      mat1 = 0;
      mat1(0,0) = -1;
      mat1(0,1) = -2;
      mat1(1,0) =  3;
      mat1(1,2) = -4;

      std::unique_ptr<int[],blaze::Deallocate> memory2( blaze::allocate<int>( 48UL ) );
      OMT mat2( memory2.get(), 2UL, 3UL, 16UL );
      mat2 = 0;
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 -= mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 48UL );
      checkNonZeros( mat2,  4UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 0UL );
      checkNonZeros( mat2,  2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major CustomMatrix dense matrix subtraction assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::columnMajor;

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>;
      std::unique_ptr<int[]> memory1( new int[7UL] );
      UnalignedUnpadded mat1( memory1.get()+1UL, 2UL, 3UL );
      mat1 = 0;
      mat1(0,0) = -1;
      mat1(0,1) = -2;
      mat1(1,0) =  3;
      mat1(1,2) = -4;

      std::unique_ptr<int[],blaze::Deallocate> memory2( blaze::allocate<int>( 48UL ) );
      OMT mat2( memory2.get(), 2UL, 3UL, 16UL );
      mat2 = 0;
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 -= mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 48UL );
      checkNonZeros( mat2,  4UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 0UL );
      checkNonZeros( mat2,  2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major CustomMatrix dense matrix subtraction assignment (lower)";

      blaze::LowerMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      OMT mat2( memory.get(), 3UL, 3UL, 16UL );
      mat2 = 0;

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
      test_ = "Column-major/column-major CustomMatrix dense matrix subtraction assignment (lower)";

      blaze::LowerMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      OMT mat2( memory.get(), 3UL, 3UL, 16UL );
      mat2 = 0;

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
      test_ = "Column-major/row-major CustomMatrix dense matrix subtraction assignment (upper)";

      blaze::UpperMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      OMT mat2( memory.get(), 3UL, 3UL, 16UL );
      mat2 = 0;

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
      test_ = "Column-major/column-major CustomMatrix dense matrix subtraction assignment (upper)";

      blaze::UpperMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      OMT mat2( memory.get(), 3UL, 3UL, 16UL );
      mat2 = 0;

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
      test_ = "Column-major/row-major CustomMatrix dense matrix subtraction assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      OMT mat2( memory.get(), 3UL, 3UL, 16UL );
      mat2 = 0;

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
      test_ = "Column-major/column-major CustomMatrix dense matrix subtraction assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      OMT mat2( memory.get(), 3UL, 3UL, 16UL );
      mat2 = 0;

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
      test_ = "Column-major/row-major CustomMatrix sparse matrix subtraction assignment";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1( 2UL, 3UL, 4UL );
      mat1(0,0) = -1;
      mat1(0,1) = -2;
      mat1(1,0) =  3;
      mat1(1,2) = -4;

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      OMT mat2( memory.get(), 2UL, 3UL, 16UL );
      mat2 = 0;
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 -= mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 48UL );
      checkNonZeros( mat2,  4UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 0UL );
      checkNonZeros( mat2,  2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major CustomMatrix sparse matrix subtraction assignment";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1( 2UL, 3UL, 4UL );
      mat1(0,0) = -1;
      mat1(0,1) = -2;
      mat1(1,0) =  3;
      mat1(1,2) = -4;

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      OMT mat2( memory.get(), 2UL, 3UL, 16UL );
      mat2 = 0;
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 -= mat1;

      checkRows    ( mat2,  2UL );
      checkColumns ( mat2,  3UL );
      checkCapacity( mat2, 48UL );
      checkNonZeros( mat2,  4UL );
      checkNonZeros( mat2,  0UL, 2UL );
      checkNonZeros( mat2,  1UL, 0UL );
      checkNonZeros( mat2,  2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major CustomMatrix sparse matrix subtraction assignment (lower)";

      blaze::LowerMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      OMT mat2( memory.get(), 3UL, 3UL, 16UL );
      mat2 = 0;

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
      test_ = "Column-major/column-major CustomMatrix sparse matrix subtraction assignment (lower)";

      blaze::LowerMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      OMT mat2( memory.get(), 3UL, 3UL, 16UL );
      mat2 = 0;

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
      test_ = "Column-major/row-major CustomMatrix sparse matrix subtraction assignment (upper)";

      blaze::UpperMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      OMT mat2( memory.get(), 3UL, 3UL, 16UL );
      mat2 = 0;

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
      test_ = "Column-major/column-major CustomMatrix sparse matrix subtraction assignment (upper)";

      blaze::UpperMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      OMT mat2( memory.get(), 3UL, 3UL, 16UL );
      mat2 = 0;

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
      test_ = "Column-major/row-major CustomMatrix sparse matrix subtraction assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      OMT mat2( memory.get(), 3UL, 3UL, 16UL );
      mat2 = 0;

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
      test_ = "Column-major/column-major CustomMatrix sparse matrix subtraction assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      OMT mat2( memory.get(), 3UL, 3UL, 16UL );
      mat2 = 0;

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

} // namespace custommatrix

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
   std::cout << "   Running aligned/unpadded CustomMatrix class test (part 1)..." << std::endl;

   try
   {
      RUN_CUSTOMMATRIX_ALIGNED_UNPADDED_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during aligned/unpadded CustomMatrix class test (part 1):\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
