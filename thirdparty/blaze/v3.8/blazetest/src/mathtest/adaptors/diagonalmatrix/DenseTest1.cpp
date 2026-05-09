//=================================================================================================
/*!
//  \file src/mathtest/adaptors/diagonalmatrix/DenseTest1.cpp
//  \brief Source file for the DiagonalMatrix dense test (part 1)
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
#include <blaze/math/HybridMatrix.h>
#include <blaze/math/StaticMatrix.h>
#include <blaze/util/policies/ArrayDelete.h>
#include <blazetest/mathtest/adaptors/diagonalmatrix/DenseTest.h>

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
/*!\brief Constructor for the DiagonalMatrix dense test.
//
// \exception std::runtime_error Operation error detected.
*/
DenseTest::DenseTest()
{
   testConstructors();
   testAssignment();
   testAddAssign();
   testSubAssign();
   testSchurAssign();
   testMultAssign();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the DiagonalMatrix constructors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all constructors of the DiagonalMatrix specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testConstructors()
{
   //=====================================================================================
   // Row-major default constructor
   //=====================================================================================

   // Default constructor (StaticMatrix)
   {
      test_ = "Row-major DiagonalMatrix default constructor (StaticMatrix)";

      const blaze::DiagonalMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> > diag;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkNonZeros( diag, 0UL );
   }

   // Default constructor (HybridMatrix)
   {
      test_ = "Row-major DiagonalMatrix default constructor (HybridMatrix)";

      const blaze::DiagonalMatrix< blaze::HybridMatrix<int,3UL,3UL,blaze::rowMajor> > diag;

      checkRows    ( diag, 0UL );
      checkColumns ( diag, 0UL );
      checkNonZeros( diag, 0UL );
   }

   // Default constructor (DynamicMatrix)
   {
      test_ = "Row-major DiagonalMatrix default constructor (DynamicMatrix)";

      const DT diag;

      checkRows    ( diag, 0UL );
      checkColumns ( diag, 0UL );
      checkNonZeros( diag, 0UL );
   }


   //=====================================================================================
   // Row-major single argument constructor
   //=====================================================================================

   // Single argument constructor (StaticMatrix)
   {
      test_ = "Row-major DiagonalMatrix single argument constructor (StaticMatrix)";

      const blaze::DiagonalMatrix< blaze::StaticMatrix<int,2UL,2UL,blaze::rowMajor> > diag( 5 );

      checkRows    ( diag, 2UL );
      checkColumns ( diag, 2UL );
      checkCapacity( diag, 4UL );
      checkNonZeros( diag, 2UL );

      if( diag(0,0) != 5 || diag(0,1) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 5 0 )\n( 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Single argument constructor (HybridMatrix)
   {
      test_ = "Row-major DiagonalMatrix single argument constructor (HybridMatrix)";

      const blaze::DiagonalMatrix< blaze::HybridMatrix<int,3UL,3UL,blaze::rowMajor> > diag( 2UL );

      checkRows    ( diag, 2UL );
      checkColumns ( diag, 2UL );
      checkCapacity( diag, 4UL );
      checkNonZeros( diag, 0UL );

      if( diag(0,0) != 0 || diag(0,1) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 0 0 )\n( 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Single argument constructor (DynamicMatrix)
   {
      test_ = "Row-major DiagonalMatrix single argument constructor (DynamicMatrix)";

      const DT diag( 2UL );

      checkRows    ( diag, 2UL );
      checkColumns ( diag, 2UL );
      checkCapacity( diag, 4UL );
      checkNonZeros( diag, 0UL );

      if( diag(0,0) != 0 || diag(0,1) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 0 0 )\n( 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Single argument constructor (0x0)
   {
      test_ = "Row-major DiagonalMatrix single argument constructor (0x0)";

      const blaze::DynamicMatrix<int,blaze::rowMajor> mat;
      const DT diag( mat );

      checkRows    ( diag, 0UL );
      checkColumns ( diag, 0UL );
      checkNonZeros( diag, 0UL );
   }

   // Single argument constructor (diagonal)
   {
      test_ = "Row-major DiagonalMatrix single argument constructor (diagonal)";

      blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> mat;
      mat(0,0) = 1;
      mat(1,1) = 2;
      mat(2,2) = 3;

      const DT diag( mat );

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Single argument constructor (lower)
   {
      test_ = "Row-major DiagonalMatrix single argument constructor (lower)";

      blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> mat;
      mat(0,0) = 1;
      mat(1,1) = 2;
      mat(2,0) = 5;
      mat(2,2) = 3;

      try {
         const DT diag( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of non-diagonal DiagonalMatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Single argument constructor (upper)
   {
      test_ = "Row-major DiagonalMatrix single argument constructor (upper)";

      blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> mat;
      mat(0,0) = 1;
      mat(0,2) = 5;
      mat(1,1) = 2;
      mat(2,2) = 3;

      try {
         const DT diag( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of non-diagonal DiagonalMatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Single argument constructor (DiagonalMatrix)
   {
      test_ = "Row-major DiagonalMatrix single argument constructor (DiagonalMatrix)";

      blaze::DiagonalMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> > diag1;
      diag1(0,0) = 1;
      diag1(1,1) = 2;
      diag1(2,2) = 3;

      const DT diag2( diag1 );

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 3UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 1UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 2 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major two argument constructor
   //=====================================================================================

   // Two argument constructor (HybridMatrix)
   {
      test_ = "Row-major DiagonalMatrix two argument constructor (HybridMatrix)";

      const blaze::DiagonalMatrix< blaze::HybridMatrix<int,3UL,3UL,blaze::rowMajor> > diag( 2UL, 5 );

      checkRows    ( diag, 2UL );
      checkColumns ( diag, 2UL );
      checkCapacity( diag, 4UL );
      checkNonZeros( diag, 2UL );

      if( diag(0,0) != 5 || diag(0,1) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 5 0 )\n( 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Two argument constructor (DynamicMatrix)
   {
      test_ = "Row-major DiagonalMatrix two argument constructor (DynamicMatrix)";

      const DT diag( 2UL, 5 );

      checkRows    ( diag, 2UL );
      checkColumns ( diag, 2UL );
      checkCapacity( diag, 4UL );
      checkNonZeros( diag, 2UL );

      if( diag(0,0) != 5 || diag(0,1) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 5 0 )\n( 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major list initialization
   //=====================================================================================

   // Complete initializer list
   {
      test_ = "Row-major DiagonalMatrix initializer list constructor (complete list)";

      const DT diag{ { 1, 0, 0 }, { 0, 2, 0 }, { 0, 0, 3 } };

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 3UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Incomplete initializer list
   {
      test_ = "Row-major DiagonalMatrix initializer list constructor (incomplete list)";

      const DT diag{ { 1 }, { 0, 2 }, { 0, 0, 3 } };

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 3UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major array initialization
   //=====================================================================================

   // Dynamic array initialization constructor
   {
      test_ = "Row-major DiagonalMatrix dynamic array initialization constructor";

      std::unique_ptr<int[]> array( new int[9] );
      array[0] = 1;
      array[1] = 0;
      array[2] = 0;
      array[3] = 0;
      array[4] = 2;
      array[5] = 0;
      array[6] = 0;
      array[7] = 0;
      array[8] = 3;
      const DT diag( 3UL, array.get() );

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 3UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Static array initialization constructor
   {
      test_ = "Row-major DiagonalMatrix static array initialization constructor";

      const int array[3][3] = { { 1, 0, 0 }, { 0, 2, 0 }, { 0, 0, 3 } };
      const DT diag( array );

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 3UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major custom matrix constructors
   //=====================================================================================

   // Custom matrix constructor (ElementType*, size_t)
   {
      test_ = "Row-major DiagonalMatrix custom matrix constructor (ElementType*, size_t)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowMajor;

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>;
      std::unique_ptr<int[]> memory( new int[5UL] );
      memory[1] = 1;
      memory[2] = 0;
      memory[3] = 0;
      memory[4] = 2;
      const blaze::DiagonalMatrix<UnalignedUnpadded> diag( memory.get()+1UL, 2UL );

      checkRows    ( diag, 2UL );
      checkColumns ( diag, 2UL );
      checkCapacity( diag, 4UL );
      checkNonZeros( diag, 2UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 )\n( 0 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Custom matrix constructor (ElementType*, size_t, size_t)
   {
      test_ = "Row-major DiagonalMatrix custom matrix constructor (ElementType*, size_t, size_t)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowMajor;

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>;
      std::unique_ptr<int[]> memory( new int[11UL] );
      memory[1] = 1;
      memory[2] = 0;
      memory[6] = 0;
      memory[7] = 2;
      const blaze::DiagonalMatrix<UnalignedUnpadded> diag( memory.get()+1UL, 2UL, 5UL );

      checkRows    ( diag, 2UL );
      checkColumns ( diag, 2UL );
      checkCapacity( diag, 4UL );
      checkNonZeros( diag, 2UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 )\n( 0 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major copy constructor
   //=====================================================================================

   // Copy constructor (0x0)
   {
      test_ = "Row-major DiagonalMatrix copy constructor (0x0)";

      const DT diag1;
      const DT diag2( diag1 );

      checkRows    ( diag2, 0UL );
      checkColumns ( diag2, 0UL );
      checkNonZeros( diag2, 0UL );
   }

   // Copy constructor (3x3)
   {
      test_ = "Row-major DiagonalMatrix copy constructor (3x3)";

      DT diag1( 3UL );
      diag1(0,0) = 1;
      diag1(1,1) = 2;
      diag1(2,2) = 3;

      const DT diag2( diag1 );

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 3UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 1UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 2 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major move constructor
   //=====================================================================================

   // Move constructor (0x0)
   {
      test_ = "Row-major DiagonalMatrix move constructor (0x0)";

      DT diag1;
      DT diag2( std::move( diag1 ) );

      checkRows    ( diag2, 0UL );
      checkColumns ( diag2, 0UL );
      checkNonZeros( diag2, 0UL );
   }

   // Move constructor (3x3)
   {
      test_ = "Row-major DiagonalMatrix move constructor (3x3)";

      DT diag1( 3UL );
      diag1(0,0) = 1;
      diag1(1,1) = 2;
      diag1(2,2) = 3;

      DT diag2( std::move( diag1 ) );

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 3UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 1UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 2 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major default constructor
   //=====================================================================================

   // Default constructor (StaticMatrix)
   {
      test_ = "Column-major DiagonalMatrix default constructor (StaticMatrix)";

      blaze::DiagonalMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> > diag;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkNonZeros( diag, 0UL );
   }

   // Default constructor (HybridMatrix)
   {
      test_ = "Column-major DiagonalMatrix default constructor (HybridMatrix)";

      blaze::DiagonalMatrix< blaze::HybridMatrix<int,3UL,3UL,blaze::columnMajor> > diag;

      checkRows    ( diag, 0UL );
      checkColumns ( diag, 0UL );
      checkNonZeros( diag, 0UL );
   }

   // Default constructor (DynamicMatrix)
   {
      test_ = "Column-major DiagonalMatrix default constructor (DynamicMatrix)";

      ODT diag;

      checkRows    ( diag, 0UL );
      checkColumns ( diag, 0UL );
      checkNonZeros( diag, 0UL );
   }


   //=====================================================================================
   // Column-major single argument constructor
   //=====================================================================================

   // Single argument constructor (StaticMatrix)
   {
      test_ = "Column-major DiagonalMatrix single argument constructor (StaticMatrix)";

      const blaze::DiagonalMatrix< blaze::StaticMatrix<int,2UL,2UL,blaze::columnMajor> > diag( 5 );

      checkRows    ( diag, 2UL );
      checkColumns ( diag, 2UL );
      checkCapacity( diag, 4UL );
      checkNonZeros( diag, 2UL );

      if( diag(0,0) != 5 || diag(0,1) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 5 0 )\n( 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Single argument constructor (HybridMatrix)
   {
      test_ = "Column-major DiagonalMatrix single argument constructor (HybridMatrix)";

      const blaze::DiagonalMatrix< blaze::HybridMatrix<int,3UL,3UL,blaze::columnMajor> > diag( 2UL );

      checkRows    ( diag, 2UL );
      checkColumns ( diag, 2UL );
      checkCapacity( diag, 4UL );
      checkNonZeros( diag, 0UL );

      if( diag(0,0) != 0 || diag(0,1) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 0 0 )\n( 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Single argument constructor (DynamicMatrix)
   {
      test_ = "Column-major DiagonalMatrix single argument constructor (DynamicMatrix)";

      const ODT diag( 2UL );

      checkRows    ( diag, 2UL );
      checkColumns ( diag, 2UL );
      checkCapacity( diag, 4UL );
      checkNonZeros( diag, 0UL );

      if( diag(0,0) != 0 || diag(0,1) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 0 0 )\n( 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Single argument constructor (0x0)
   {
      test_ = "Column-major DiagonalMatrix single argument constructor (0x0)";

      const blaze::DynamicMatrix<int,blaze::columnMajor> mat;
      const ODT diag( mat );

      checkRows    ( diag, 0UL );
      checkColumns ( diag, 0UL );
      checkNonZeros( diag, 0UL );
   }

   // Single argument constructor (diagonal)
   {
      test_ = "Column-major DiagonalMatrix single argument constructor (diagonal)";

      blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> mat;
      mat(0,0) = 1;
      mat(1,1) = 2;
      mat(2,2) = 3;

      const ODT diag( mat );

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Single argument constructor (lower)
   {
      test_ = "Column-major DiagonalMatrix single argument constructor (lower)";

      blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> mat;
      mat(0,0) = 1;
      mat(1,1) = 2;
      mat(2,0) = 5;
      mat(2,2) = 3;

      try {
         const ODT diag( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of non-diagonal DiagonalMatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Single argument constructor (upper)
   {
      test_ = "Column-major DiagonalMatrix single argument constructor (upper)";

      blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> mat;
      mat(0,0) = 1;
      mat(0,2) = 5;
      mat(1,1) = 2;
      mat(2,2) = 3;

      try {
         const ODT diag( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of non-diagonal DiagonalMatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Single argument constructor (DiagonalMatrix)
   {
      test_ = "Column-major DiagonalMatrix single argument constructor (DiagonalMatrix)";

      blaze::DiagonalMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> > diag1;
      diag1(0,0) = 1;
      diag1(1,1) = 2;
      diag1(2,2) = 3;

      const ODT diag2( diag1 );

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 3UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 1UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 2 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major two argument constructor
   //=====================================================================================

   // Two argument constructor (HybridMatrix)
   {
      test_ = "Column-major DiagonalMatrix two argument constructor (HybridMatrix)";

      const blaze::DiagonalMatrix< blaze::HybridMatrix<int,3UL,3UL,blaze::columnMajor> > diag( 2UL, 5 );

      checkRows    ( diag, 2UL );
      checkColumns ( diag, 2UL );
      checkCapacity( diag, 4UL );
      checkNonZeros( diag, 2UL );

      if( diag(0,0) != 5 || diag(0,1) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 5 0 )\n( 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Two argument constructor (DynamicMatrix)
   {
      test_ = "Column-major DiagonalMatrix two argument constructor (DynamicMatrix)";

      const ODT diag( 2UL, 5 );

      checkRows    ( diag, 2UL );
      checkColumns ( diag, 2UL );
      checkCapacity( diag, 4UL );
      checkNonZeros( diag, 2UL );

      if( diag(0,0) != 5 || diag(0,1) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 5 0 )\n( 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major list initialization
   //=====================================================================================

   // Complete initializer list
   {
      test_ = "Column-major DiagonalMatrix initializer list constructor (complete list)";

      const ODT diag{ { 1, 0, 0 }, { 0, 2, 0 }, { 0, 0, 3 } };

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 3UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Incomplete initializer list
   {
      test_ = "Column-major DiagonalMatrix initializer list constructor (incomplete list)";

      const ODT diag{ { 1 }, { 0, 2 }, { 0, 0, 3 } };

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 3UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major array initialization
   //=====================================================================================

   // Dynamic array initialization constructor
   {
      test_ = "Column-major DiagonalMatrix dynamic array initialization constructor";

      std::unique_ptr<int[]> array( new int[9] );
      array[0] = 1;
      array[1] = 0;
      array[2] = 0;
      array[3] = 0;
      array[4] = 2;
      array[5] = 0;
      array[6] = 0;
      array[7] = 0;
      array[8] = 3;
      const ODT diag( 3UL, array.get() );

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 3UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Static array initialization constructor
   {
      test_ = "Column-major DiagonalMatrix static array initialization constructor";

      const int array[3][3] = { { 1, 0, 0 }, { 0, 2, 0 }, { 0, 0, 3 } };
      const ODT diag( array );

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 3UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major custom matrix constructors
   //=====================================================================================

   // Custom matrix constructor (ElementType*, size_t)
   {
      test_ = "Column-major DiagonalMatrix custom matrix constructor (ElementType*, size_t)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::columnMajor;

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>;
      std::unique_ptr<int[]> memory( new int[5UL] );
      memory[1] = 1;
      memory[2] = 0;
      memory[3] = 0;
      memory[4] = 2;
      const blaze::DiagonalMatrix<UnalignedUnpadded> diag( memory.get()+1UL, 2UL );

      checkRows    ( diag, 2UL );
      checkColumns ( diag, 2UL );
      checkCapacity( diag, 4UL );
      checkNonZeros( diag, 2UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 )\n( 0 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Custom matrix constructor (ElementType*, size_t, size_t)
   {
      test_ = "Column-major DiagonalMatrix custom matrix constructor (ElementType*, size_t, size_t)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::columnMajor;

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>;
      std::unique_ptr<int[]> memory( new int[11UL] );
      memory[1] = 1;
      memory[2] = 0;
      memory[6] = 0;
      memory[7] = 2;
      const blaze::DiagonalMatrix<UnalignedUnpadded> diag( memory.get()+1UL, 2UL, 5UL );

      checkRows    ( diag, 2UL );
      checkColumns ( diag, 2UL );
      checkCapacity( diag, 4UL );
      checkNonZeros( diag, 2UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 )\n( 0 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major copy constructor
   //=====================================================================================

   // Copy constructor (0x0)
   {
      test_ = "Column-major DiagonalMatrix copy constructor (0x0)";

      const ODT diag1;
      const ODT diag2( diag1 );

      checkRows    ( diag2, 0UL );
      checkColumns ( diag2, 0UL );
      checkNonZeros( diag2, 0UL );
   }

   // Copy constructor (3x3)
   {
      test_ = "Column-major DiagonalMatrix copy constructor (3x3)";

      ODT diag1( 3UL );
      diag1(0,0) = 1;
      diag1(1,1) = 2;
      diag1(2,2) = 3;

      const ODT diag2( diag1 );

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 3UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 1UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 2 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major move constructor
   //=====================================================================================

   // Move constructor (0x0)
   {
      test_ = "Column-major DiagonalMatrix move constructor (0x0)";

      ODT diag1;
      ODT diag2( std::move( diag1 ) );

      checkRows    ( diag2, 0UL );
      checkColumns ( diag2, 0UL );
      checkNonZeros( diag2, 0UL );
   }

   // Move constructor (3x3)
   {
      test_ = "Column-major DiagonalMatrix move constructor (3x3)";

      ODT diag1( 3UL );
      diag1(0,0) = 1;
      diag1(1,1) = 2;
      diag1(2,2) = 3;

      ODT diag2( std::move( diag1 ) );

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 3UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 1UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 2 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the DiagonalMatrix assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all assignment operators of the DiagonalMatrix specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testAssignment()
{
   //=====================================================================================
   // Row-major homogeneous assignment
   //=====================================================================================

   // Homogeneous assignment (3x3)
   {
      test_ = "Row-major DiagonalMatrix homogeneous assignment (3x3)";

      DT diag( 3UL );
      diag = 2;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 2 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 2 0 0 )\n( 0 2 0 )\n( 0 0 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major list assignment
   //=====================================================================================

   // Complete initializer list
   {
      test_ = "Row-major DiagonalMatrix initializer list assignment (complete list)";

      DT diag;
      diag = { { 1, 0, 0 }, { 0, 2, 0 }, { 0, 0, 3 } };

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Incomplete initializer list
   {
      test_ = "Row-major DiagonalMatrix initializer list assignment (incomplete list)";

      DT diag;
      diag = { { 1 }, { 0, 2 }, { 0, 0, 3 } };

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major array assignment
   //=====================================================================================

   // Array assignment
   {
      test_ = "Row-major DiagonalMatrix array assignment";

      const int array[3][3] = { { 1, 0, 0 }, { 0, 2, 0 }, { 0, 0, 3 } };
      DT diag;
      diag = array;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major copy assignment
   //=====================================================================================

   // Copy assignment (0x0)
   {
      test_ = "Row-major DiagonalMatrix copy assignment (0x0)";

      DT diag1, diag2;

      diag2 = diag1;

      checkRows    ( diag2, 0UL );
      checkColumns ( diag2, 0UL );
      checkNonZeros( diag2, 0UL );
   }

   // Copy assignment (3x3)
   {
      test_ = "Row-major DiagonalMatrix copy assignment (3x3)";

      DT diag1( 3UL );
      diag1(0,0) = 1;
      diag1(1,1) = 2;
      diag1(2,2) = 3;

      DT diag2;
      diag2 = diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkNonZeros( diag2, 3UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 1UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 2 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major move assignment
   //=====================================================================================

   // Move assignment (0x0)
   {
      test_ = "Row-major DiagonalMatrix move assignment (0x0)";

      DT diag1, diag2;

      diag2 = std::move( diag1 );

      checkRows    ( diag2, 0UL );
      checkColumns ( diag2, 0UL );
      checkNonZeros( diag2, 0UL );
   }

   // Move assignment (3x3)
   {
      test_ = "Row-major DiagonalMatrix move assignment (3x3)";

      DT diag1( 3UL );
      diag1(0,0) = 1;
      diag1(1,1) = 2;
      diag1(2,2) = 3;

      DT diag2;
      diag2 = std::move( diag1 );

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkNonZeros( diag2, 3UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 1UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 2 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major dense matrix assignment
   //=====================================================================================

   // Conversion assignment (0x0)
   {
      test_ = "Row-major DiagonalMatrix dense matrix assignment (0x0)";

      const blaze::DynamicMatrix<int,blaze::rowMajor> mat;

      DT diag;
      diag = mat;

      checkRows    ( diag, 0UL );
      checkColumns ( diag, 0UL );
      checkNonZeros( diag, 0UL );
   }

   // Row-major/row-major dense matrix assignment (diagonal)
   {
      test_ = "Row-major/row-major DiagonalMatrix dense matrix assignment (diagonal)";

      blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> mat;
      mat(0,0) = 1;
      mat(1,1) = 2;
      mat(2,2) = 3;

      DT diag;
      diag = mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix assignment (diagonal)
   {
      test_ = "Row-major/column-major DiagonalMatrix dense matrix assignment (diagonal)";

      blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> mat;
      mat(0,0) = 1;
      mat(1,1) = 2;
      mat(2,2) = 3;

      DT diag;
      diag = mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major dense matrix assignment (lower)
   {
      test_ = "Row-major/row-major DiagonalMatrix dense matrix assignment (lower)";

      blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> mat;
      mat(0,0) = 1;
      mat(1,1) = 2;
      mat(2,0) = 5;
      mat(2,2) = 3;

      try {
         DT diag;
         diag = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of lower row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major dense matrix assignment (lower)
   {
      test_ = "Row-major/column-major DiagonalMatrix dense matrix assignment (lower)";

      blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> mat;
      mat(0,0) = 1;
      mat(1,1) = 2;
      mat(2,0) = 5;
      mat(2,2) = 3;

      try {
         DT diag;
         diag = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of lower column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major dense matrix assignment (upper)
   {
      test_ = "Row-major/row-major DiagonalMatrix dense matrix assignment (upper)";

      blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> mat;
      mat(0,0) = 1;
      mat(0,2) = 5;
      mat(1,1) = 2;
      mat(2,2) = 3;

      try {
         DT diag;
         diag = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of upper row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major dense matrix assignment (upper)
   {
      test_ = "Row-major/column-major DiagonalMatrix dense matrix assignment (upper)";

      blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> mat;
      mat(0,0) = 1;
      mat(0,2) = 5;
      mat(1,1) = 2;
      mat(2,2) = 3;

      try {
         DT diag;
         diag = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of upper column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major dense matrix assignment (DiagonalMatrix)
   {
      test_ = "Row-major/row-major DiagonalMatrix dense matrix assignment (DiagonalMatrix)";

      blaze::DiagonalMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> > diag1;
      diag1(0,0) = 1;
      diag1(1,1) = 2;
      diag1(2,2) = 3;

      DT diag2;
      diag2 = diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkNonZeros( diag2, 3UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 1UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 2 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix assignment (DiagonalMatrix)
   {
      test_ = "Row-major/column-major DiagonalMatrix dense matrix assignment (DiagonalMatrix)";

      blaze::DiagonalMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> > diag1;
      diag1(0,0) = 1;
      diag1(1,1) = 2;
      diag1(2,2) = 3;

      DT diag2;
      diag2 = diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkNonZeros( diag2, 3UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 1UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 2 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix assignment
   //=====================================================================================

   // Conversion assignment (0x0)
   {
      test_ = "Row-major DiagonalMatrix sparse matrix assignment (0x0)";

      const blaze::CompressedMatrix<int,blaze::rowMajor> mat;

      DT diag;
      diag = mat;

      checkRows    ( diag, 0UL );
      checkColumns ( diag, 0UL );
      checkNonZeros( diag, 0UL );
   }

   // Row-major/row-major sparse matrix assignment (diagonal)
   {
      test_ = "Row-major/row-major DiagonalMatrix sparse matrix assignment (diagonal)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 5UL );
      mat(0,0) = 1;
      mat(1,1) = 2;
      mat(2,2) = 3;
      mat.insert( 1UL, 2UL, 0 );
      mat.insert( 2UL, 1UL, 0 );

      DT diag;
      diag = mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major sparse matrix assignment (diagonal)
   {
      test_ = "Row-major/column-major DiagonalMatrix sparse matrix assignment (diagonal)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 5UL );
      mat(0,0) = 1;
      mat(1,1) = 2;
      mat(2,2) = 3;
      mat.insert( 1UL, 2UL, 0 );
      mat.insert( 2UL, 1UL, 0 );

      DT diag;
      diag = mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major sparse matrix assignment (lower)
   {
      test_ = "Row-major/row-major DiagonalMatrix sparse matrix assignment (lower)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 4UL );
      mat(0,0) = 1;
      mat(1,1) = 2;
      mat(2,0) = 5;
      mat(2,2) = 3;

      try {
         DT diag;
         diag = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of lower row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major sparse matrix assignment (lower)
   {
      test_ = "Row-major/column-major DiagonalMatrix sparse matrix assignment (lower)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 4UL );
      mat(0,0) = 1;
      mat(1,1) = 2;
      mat(2,0) = 5;
      mat(2,2) = 3;

      try {
         DT diag;
         diag = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of lower column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major sparse matrix assignment (upper)
   {
      test_ = "Row-major/row-major DiagonalMatrix sparse matrix assignment (upper)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 4UL );
      mat(0,0) = 1;
      mat(0,2) = 5;
      mat(1,1) = 2;
      mat(2,2) = 3;

      try {
         DT diag;
         diag = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of upper row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major sparse matrix assignment (upper)
   {
      test_ = "Row-major/column-major DiagonalMatrix sparse matrix assignment (upper)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 4UL );
      mat(0,0) = 1;
      mat(0,2) = 5;
      mat(1,1) = 2;
      mat(2,2) = 3;

      try {
         DT diag;
         diag = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of upper column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major sparse matrix assignment (DiagonalMatrix)
   {
      test_ = "Row-major/row-major DiagonalMatrix sparse matrix assignment (DiagonalMatrix)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > diag1( 3UL, 3UL );
      diag1(0,0) = 1;
      diag1(1,1) = 2;
      diag1(2,2) = 3;

      DT diag2;
      diag2 = diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkNonZeros( diag2, 3UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 1UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 2 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major sparse matrix assignment (DiagonalMatrix)
   {
      test_ = "Row-major/column-major DiagonalMatrix sparse matrix assignment (DiagonalMatrix)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > diag1( 3UL, 3UL );
      diag1(0,0) = 1;
      diag1(1,1) = 2;
      diag1(2,2) = 3;

      DT diag2;
      diag2 = diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkNonZeros( diag2, 3UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 1UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 2 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major homogeneous assignment
   //=====================================================================================

   // Homogeneous assignment (3x3)
   {
      test_ = "Column-major DiagonalMatrix homogeneous assignment (3x3)";

      ODT diag( 3UL );
      diag = 2;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 2 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 2 0 0 )\n( 0 2 0 )\n( 0 0 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major list assignment
   //=====================================================================================

   // Complete initializer list
   {
      test_ = "Column-major DiagonalMatrix initializer list assignment (complete list)";

      ODT diag;
      diag = { { 1, 0, 0 }, { 0, 2, 0 }, { 0, 0, 3 } };

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Incomplete initializer list
   {
      test_ = "Column-major DiagonalMatrix initializer list assignment (incomplete list)";

      ODT diag;
      diag = { { 1 }, { 0, 2 }, { 0, 0, 3 } };

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major array assignment
   //=====================================================================================

   // Array assignment
   {
      test_ = "Column-major DiagonalMatrix array assignment";

      const int array[3][3] = { { 1, 0, 0 }, { 0, 2, 0 }, { 0, 0, 3 } };
      ODT diag;
      diag = array;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major copy assignment
   //=====================================================================================

   // Copy assignment (0x0)
   {
      test_ = "Column-major DiagonalMatrix copy assignment (0x0)";

      ODT diag1, diag2;

      diag2 = diag1;

      checkRows    ( diag2, 0UL );
      checkColumns ( diag2, 0UL );
      checkNonZeros( diag2, 0UL );
   }

   // Copy assignment (3x3)
   {
      test_ = "Column-major DiagonalMatrix copy assignment (3x3)";

      ODT diag1( 3UL );
      diag1(0,0) = 1;
      diag1(1,1) = 2;
      diag1(2,2) = 3;

      ODT diag2;
      diag2 = diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkNonZeros( diag2, 3UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 1UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 2 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major move assignment
   //=====================================================================================

   // Move assignment (0x0)
   {
      test_ = "Column-major DiagonalMatrix move assignment (0x0)";

      ODT diag1, diag2;

      diag2 = std::move( diag1 );

      checkRows    ( diag2, 0UL );
      checkColumns ( diag2, 0UL );
      checkNonZeros( diag2, 0UL );
   }

   // Move assignment (3x3)
   {
      test_ = "Column-major DiagonalMatrix move assignment (3x3)";

      ODT diag1( 3UL );
      diag1(0,0) = 1;
      diag1(1,1) = 2;
      diag1(2,2) = 3;

      ODT diag2;
      diag2 = std::move( diag1 );

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkNonZeros( diag2, 3UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 1UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 2 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix assignment
   //=====================================================================================

   // Conversion assignment (0x0)
   {
      test_ = "Column-major DiagonalMatrix dense matrix assignment (0x0)";

      const blaze::DynamicMatrix<int,blaze::rowMajor> mat;

      ODT diag;
      diag = mat;

      checkRows    ( diag, 0UL );
      checkColumns ( diag, 0UL );
      checkNonZeros( diag, 0UL );
   }

   // Column-major/row-major dense matrix assignment (diagonal)
   {
      test_ = "Column-major/row-major DiagonalMatrix dense matrix assignment (diagonal)";

      blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> mat;
      mat(0,0) = 1;
      mat(1,1) = 2;
      mat(2,2) = 3;

      ODT diag;
      diag = mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix assignment (diagonal)
   {
      test_ = "Column-major/column-major DiagonalMatrix dense matrix assignment (diagonal)";

      blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> mat;
      mat(0,0) = 1;
      mat(1,1) = 2;
      mat(2,2) = 3;

      ODT diag;
      diag = mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major dense matrix assignment (lower)
   {
      test_ = "Column-major/row-major DiagonalMatrix dense matrix assignment (lower)";

      blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> mat;
      mat(0,0) = 1;
      mat(1,1) = 2;
      mat(2,0) = 5;
      mat(2,2) = 3;

      try {
         ODT diag;
         diag = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of lower row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major dense matrix assignment (lower)
   {
      test_ = "Column-major/column-major DiagonalMatrix dense matrix assignment (lower)";

      blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> mat;
      mat(0,0) = 1;
      mat(1,1) = 2;
      mat(2,0) = 5;
      mat(2,2) = 3;

      try {
         ODT diag;
         diag = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of lower column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major dense matrix assignment (upper)
   {
      test_ = "Column-major/row-major DiagonalMatrix dense matrix assignment (upper)";

      blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> mat;
      mat(0,0) = 1;
      mat(0,2) = 5;
      mat(1,1) = 2;
      mat(2,2) = 3;

      try {
         ODT diag;
         diag = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of upper row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major dense matrix assignment (upper)
   {
      test_ = "Column-major/column-major DiagonalMatrix dense matrix assignment (upper)";

      blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> mat;
      mat(0,0) = 1;
      mat(0,2) = 5;
      mat(1,1) = 2;
      mat(2,2) = 3;

      try {
         ODT diag;
         diag = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of upper column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major dense matrix assignment (DiagonalMatrix)
   {
      test_ = "Column-major/row-major DiagonalMatrix dense matrix assignment (DiagonalMatrix)";

      blaze::DiagonalMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> > diag1;
      diag1(0,0) = 1;
      diag1(1,1) = 2;
      diag1(2,2) = 3;

      ODT diag2;
      diag2 = diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkNonZeros( diag2, 3UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 1UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 2 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix assignment (DiagonalMatrix)
   {
      test_ = "Column-major/column-major DiagonalMatrix dense matrix assignment (DiagonalMatrix)";

      blaze::DiagonalMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> > diag1;
      diag1(0,0) = 1;
      diag1(1,1) = 2;
      diag1(2,2) = 3;

      ODT diag2;
      diag2 = diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkNonZeros( diag2, 3UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 1UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 2 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse matrix assignment
   //=====================================================================================

   // Conversion assignment (0x0)
   {
      test_ = "Column-major DiagonalMatrix sparse matrix assignment (0x0)";

      const blaze::CompressedMatrix<int,blaze::rowMajor> mat;

      ODT diag;
      diag = mat;

      checkRows    ( diag, 0UL );
      checkColumns ( diag, 0UL );
      checkNonZeros( diag, 0UL );
   }

   // Column-major/row-major sparse matrix assignment (diagonal)
   {
      test_ = "Column-major/row-major DiagonalMatrix sparse matrix assignment (diagonal)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 5UL );
      mat(0,0) = 1;
      mat(1,1) = 2;
      mat(2,2) = 3;
      mat.insert( 1UL, 2UL, 0 );
      mat.insert( 2UL, 1UL, 0 );

      ODT diag;
      diag = mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major sparse matrix assignment (diagonal)
   {
      test_ = "Column-major/column-major DiagonalMatrix sparse matrix assignment (diagonal)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 5UL );
      mat(0,0) = 1;
      mat(1,1) = 2;
      mat(2,2) = 3;
      mat.insert( 1UL, 2UL, 0 );
      mat.insert( 2UL, 1UL, 0 );

      ODT diag;
      diag = mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 2 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major sparse matrix assignment (lower)
   {
      test_ = "Column-major/row-major DiagonalMatrix sparse matrix assignment (lower)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 4UL );
      mat(0,0) = 1;
      mat(1,1) = 2;
      mat(2,0) = 5;
      mat(2,2) = 3;

      try {
         ODT diag;
         diag = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of lower row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major sparse matrix assignment (lower)
   {
      test_ = "Column-major/column-major DiagonalMatrix sparse matrix assignment (lower)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 4UL );
      mat(0,0) = 1;
      mat(1,1) = 2;
      mat(2,0) = 5;
      mat(2,2) = 3;

      try {
         ODT diag;
         diag = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of lower column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major sparse matrix assignment (upper)
   {
      test_ = "Column-major/row-major DiagonalMatrix sparse matrix assignment (upper)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 4UL );
      mat(0,0) = 1;
      mat(0,2) = 5;
      mat(1,1) = 2;
      mat(2,2) = 3;

      try {
         ODT diag;
         diag = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of upper row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major sparse matrix assignment (upper)
   {
      test_ = "Column-major/column-major DiagonalMatrix sparse matrix assignment (upper)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 4UL );
      mat(0,0) = 1;
      mat(0,2) = 5;
      mat(1,1) = 2;
      mat(2,2) = 3;

      try {
         ODT diag;
         diag = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of upper column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major sparse matrix assignment (DiagonalMatrix)
   {
      test_ = "Column-major/row-major DiagonalMatrix sparse matrix assignment (DiagonalMatrix)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > diag1( 3UL, 3UL );
      diag1(0,0) = 1;
      diag1(1,1) = 2;
      diag1(2,2) = 3;

      ODT diag2;
      diag2 = diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkNonZeros( diag2, 3UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 1UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 2 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major sparse matrix assignment (DiagonalMatrix)
   {
      test_ = "Column-major/column-major DiagonalMatrix sparse matrix assignment (DiagonalMatrix)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > diag1( 3UL, 3UL );
      diag1(0,0) = 1;
      diag1(1,1) = 2;
      diag1(2,2) = 3;

      ODT diag2;
      diag2 = diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkNonZeros( diag2, 3UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 1UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 2 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 2 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the DiagonalMatrix addition assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the addition assignment operators of the DiagonalMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testAddAssign()
{
   //=====================================================================================
   // Row-major dense matrix addition assignment
   //=====================================================================================

   // Row-major/row-major dense matrix addition assignment (diagonal)
   {
      test_ = "Row-major/row-major DiagonalMatrix dense matrix addition assignment (diagonal)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(1,1) = -2;
      mat(2,2) =  2;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      diag += mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix addition assignment (diagonal)
   {
      test_ = "Row-major/column-major DiagonalMatrix dense matrix addition assignment (diagonal)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(1,1) = -2;
      mat(2,2) =  2;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      diag += mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major dense matrix addition assignment (lower)
   {
      test_ = "Row-major/row-major DiagonalMatrix dense matrix addition assignment (lower)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(2,0) = 5;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of lower row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major dense matrix addition assignment (lower)
   {
      test_ = "Row-major/column-major DiagonalMatrix dense matrix addition assignment (lower)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(2,0) = 5;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of lower column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major dense matrix addition assignment (upper)
   {
      test_ = "Row-major/row-major DiagonalMatrix dense matrix addition assignment (upper)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(0,2) = 5;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of upper row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major dense matrix addition assignment (upper)
   {
      test_ = "Row-major/column-major DiagonalMatrix dense matrix addition assignment (upper)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(0,2) = 5;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of upper column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major dense matrix addition assignment (DiagonalMatrix)
   {
      test_ = "Row-major/row-major DiagonalMatrix dense matrix addition assignment (DiagonalMatrix)";

      DT diag1( 3UL );
      diag1(1,1) = -2;
      diag1(2,2) =  2;

      DT diag2( 3UL );
      diag2(0,0) = 1;
      diag2(1,1) = 2;
      diag2(2,2) = 3;

      diag2 += diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 2UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 0UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 0 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix addition assignment (DiagonalMatrix)
   {
      test_ = "Row-major/column-major DiagonalMatrix dense matrix addition assignment (DiagonalMatrix)";

      ODT diag1( 3UL );
      diag1(1,1) = -2;
      diag1(2,2) =  2;

      DT diag2( 3UL );
      diag2(0,0) = 1;
      diag2(1,1) = 2;
      diag2(2,2) = 3;

      diag2 += diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 2UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 0UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 0 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix addition assignment
   //=====================================================================================

   // Row-major/row-major sparse matrix addition assignment (diagonal)
   {
      test_ = "Row-major/row-major DiagonalMatrix sparse matrix addition assignment (diagonal)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 4UL );
      mat(1,1) = -2;
      mat(2,2) =  2;
      mat.insert( 1UL, 2UL, 0 );
      mat.insert( 2UL, 1UL, 0 );

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      diag += mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major sparse matrix addition assignment (diagonal)
   {
      test_ = "Row-major/column-major DiagonalMatrix sparse matrix addition assignment (diagonal)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 4UL );
      mat(1,1) = -2;
      mat(2,2) =  2;
      mat.insert( 1UL, 2UL, 0 );
      mat.insert( 2UL, 1UL, 0 );

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      diag += mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major sparse matrix addition assignment (lower)
   {
      test_ = "Row-major/row-major DiagonalMatrix sparse matrix addition assignment (lower)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 1UL );
      mat(2,0) = 5;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of lower row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major sparse matrix addition assignment (lower)
   {
      test_ = "Row-major/column-major DiagonalMatrix sparse matrix addition assignment (lower)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 1UL );
      mat(2,0) = 5;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of lower column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major sparse matrix addition assignment (upper)
   {
      test_ = "Row-major/row-major DiagonalMatrix sparse matrix addition assignment (upper)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 1UL );
      mat(0,2) = 5;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of upper row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major sparse matrix addition assignment (upper)
   {
      test_ = "Row-major/column-major DiagonalMatrix sparse matrix addition assignment (upper)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 1UL );
      mat(0,2) = 5;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of upper column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major sparse matrix addition assignment (DiagonalMatrix)
   {
      test_ = "Row-major/row-major DiagonalMatrix sparse matrix addition assignment (DiagonalMatrix)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > diag1( 3UL, 2UL );
      diag1(1,1) = -2;
      diag1(2,2) =  2;

      DT diag2( 3UL );
      diag2(0,0) = 1;
      diag2(1,1) = 2;
      diag2(2,2) = 3;

      diag2 += diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 2UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 0UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 0 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major sparse matrix addition assignment (DiagonalMatrix)
   {
      test_ = "Row-major/column-major DiagonalMatrix sparse matrix addition assignment (DiagonalMatrix)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > diag1( 3UL, 2UL );
      diag1(1,1) = -2;
      diag1(2,2) =  2;

      DT diag2( 3UL );
      diag2(0,0) = 1;
      diag2(1,1) = 2;
      diag2(2,2) = 3;

      diag2 += diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 2UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 0UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 0 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix addition assignment
   //=====================================================================================

   // Column-major/row-major dense matrix addition assignment (diagonal)
   {
      test_ = "Column-major/row-major DiagonalMatrix dense matrix addition assignment (diagonal)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(1,1) = -2;
      mat(2,2) =  2;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      diag += mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix addition assignment (diagonal)
   {
      test_ = "Column-major/column-major DiagonalMatrix dense matrix addition assignment (diagonal)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(1,1) = -2;
      mat(2,2) =  2;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      diag += mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major dense matrix addition assignment (lower)
   {
      test_ = "Column-major/row-major DiagonalMatrix dense matrix addition assignment (lower)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(2,0) = 5;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of lower row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major dense matrix addition assignment (lower)
   {
      test_ = "Column-major/column-major DiagonalMatrix dense matrix addition assignment (lower)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(2,0) = 5;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of lower column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major dense matrix addition assignment (upper)
   {
      test_ = "Column-major/row-major DiagonalMatrix dense matrix addition assignment (upper)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(0,2) = 5;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of upper row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major dense matrix addition assignment (upper)
   {
      test_ = "Column-major/column-major DiagonalMatrix dense matrix addition assignment (upper)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(0,2) = 5;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of upper column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major dense matrix addition assignment (DiagonalMatrix)
   {
      test_ = "Column-major/row-major DiagonalMatrix dense matrix addition assignment (DiagonalMatrix)";

      DT diag1( 3UL );
      diag1(1,1) = -2;
      diag1(2,2) =  2;

      ODT diag2( 3UL );
      diag2(0,0) = 1;
      diag2(1,1) = 2;
      diag2(2,2) = 3;

      diag2 += diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 2UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 0UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 0 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix addition assignment (DiagonalMatrix)
   {
      test_ = "Column-major/column-major DiagonalMatrix dense matrix addition assignment (DiagonalMatrix)";

      ODT diag1( 3UL );
      diag1(1,1) = -2;
      diag1(2,2) =  2;

      ODT diag2( 3UL );
      diag2(0,0) = 1;
      diag2(1,1) = 2;
      diag2(2,2) = 3;

      diag2 += diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 2UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 0UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 0 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix addition assignment
   //=====================================================================================

   // Column-major/row-major sparse matrix addition assignment (diagonal)
   {
      test_ = "Column-major/row-major DiagonalMatrix sparse matrix addition assignment (diagonal)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 4UL );
      mat(1,1) = -2;
      mat(2,2) =  2;
      mat.insert( 1UL, 2UL, 0 );
      mat.insert( 2UL, 1UL, 0 );

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      diag += mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major sparse matrix addition assignment (diagonal)
   {
      test_ = "Column-major/column-major DiagonalMatrix sparse matrix addition assignment (diagonal)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 4UL );
      mat(1,1) = -2;
      mat(2,2) =  2;
      mat.insert( 1UL, 2UL, 0 );
      mat.insert( 2UL, 1UL, 0 );

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      diag += mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major sparse matrix addition assignment (lower)
   {
      test_ = "Column-major/row-major DiagonalMatrix sparse matrix addition assignment (lower)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 1UL );
      mat(2,0) = 5;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of lower row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major sparse matrix addition assignment (lower)
   {
      test_ = "Column-major/column-major DiagonalMatrix sparse matrix addition assignment (lower)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 1UL );
      mat(2,0) = 5;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of lower column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major sparse matrix addition assignment (upper)
   {
      test_ = "Column-major/row-major DiagonalMatrix sparse matrix addition assignment (upper)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 1UL );
      mat(0,2) = 5;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of upper row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major sparse matrix addition assignment (upper)
   {
      test_ = "Column-major/column-major DiagonalMatrix sparse matrix addition assignment (upper)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 1UL );
      mat(0,2) = 5;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of upper column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major sparse matrix addition assignment (DiagonalMatrix)
   {
      test_ = "Column-major/row-major DiagonalMatrix sparse matrix addition assignment (DiagonalMatrix)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > diag1( 3UL, 2UL );
      diag1(1,1) = -2;
      diag1(2,2) =  2;

      ODT diag2( 3UL );
      diag2(0,0) = 1;
      diag2(1,1) = 2;
      diag2(2,2) = 3;

      diag2 += diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 2UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 0UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 0 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major sparse matrix addition assignment (DiagonalMatrix)
   {
      test_ = "Column-major/column-major DiagonalMatrix sparse matrix addition assignment (DiagonalMatrix)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > diag1( 3UL, 2UL );
      diag1(1,1) = -2;
      diag1(2,2) =  2;

      ODT diag2( 3UL );
      diag2(0,0) = 1;
      diag2(1,1) = 2;
      diag2(2,2) = 3;

      diag2 += diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 2UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 0UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 0 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the DiagonalMatrix subtraction assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the subtraction assignment operators of the DiagonalMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testSubAssign()
{
   //=====================================================================================
   // Row-major dense matrix subtraction assignment
   //=====================================================================================

   // Row-major/row-major dense matrix subtraction assignment (diagonal)
   {
      test_ = "Row-major/row-major DiagonalMatrix dense matrix subtraction assignment (diagonal)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(1,1) =  2;
      mat(2,2) = -2;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      diag -= mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix subtraction assignment (diagonal)
   {
      test_ = "Row-major/column-major DiagonalMatrix dense matrix subtraction assignment (diagonal)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(1,1) =  2;
      mat(2,2) = -2;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      diag -= mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major dense matrix subtraction assignment (lower)
   {
      test_ = "Row-major/row-major DiagonalMatrix dense matrix subtraction assignment (lower)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(2,0) = 5;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of lower row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major dense matrix subtraction assignment (lower)
   {
      test_ = "Row-major/column-major DiagonalMatrix dense matrix subtraction assignment (lower)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(2,0) = 5;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of lower column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major dense matrix subtraction assignment (upper)
   {
      test_ = "Row-major/row-major DiagonalMatrix dense matrix subtraction assignment (upper)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(0,2) = 5;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of upper row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major dense matrix subtraction assignment (upper)
   {
      test_ = "Row-major/column-major DiagonalMatrix dense matrix subtraction assignment (upper)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(0,2) = 5;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of upper column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major dense matrix subtraction assignment (DiagonalMatrix)
   {
      test_ = "Row-major/row-major DiagonalMatrix dense matrix subtraction assignment (DiagonalMatrix)";

      DT diag1( 3UL );
      diag1(1,1) =  2;
      diag1(2,2) = -2;

      DT diag2( 3UL );
      diag2(0,0) = 1;
      diag2(1,1) = 2;
      diag2(2,2) = 3;

      diag2 -= diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 2UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 0UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 0 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix subtraction assignment (DiagonalMatrix)
   {
      test_ = "Row-major/column-major DiagonalMatrix dense matrix subtraction assignment (DiagonalMatrix)";

      ODT diag1( 3UL );
      diag1(1,1) =  2;
      diag1(2,2) = -2;

      DT diag2( 3UL );
      diag2(0,0) = 1;
      diag2(1,1) = 2;
      diag2(2,2) = 3;

      diag2 -= diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 2UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 0UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 0 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix subtraction assignment
   //=====================================================================================

   // Row-major/row-major sparse matrix subtraction assignment (diagonal)
   {
      test_ = "Row-major/row-major DiagonalMatrix sparse matrix subtraction assignment (diagonal)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 4UL );
      mat(1,1) =  2;
      mat(2,2) = -2;
      mat.insert( 1UL, 2UL, 0 );
      mat.insert( 2UL, 1UL, 0 );

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      diag -= mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major sparse matrix subtraction assignment (diagonal)
   {
      test_ = "Row-major/column-major DiagonalMatrix sparse matrix subtraction assignment (diagonal)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 4UL );
      mat(1,1) =  2;
      mat(2,2) = -2;
      mat.insert( 1UL, 2UL, 0 );
      mat.insert( 2UL, 1UL, 0 );

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      diag -= mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major sparse matrix subtraction assignment (lower)
   {
      test_ = "Row-major/row-major DiagonalMatrix sparse matrix subtraction assignment (lower)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 1UL );
      mat(2,0) = 5;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of lower row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major sparse matrix subtraction assignment (lower)
   {
      test_ = "Row-major/column-major DiagonalMatrix sparse matrix subtraction assignment (lower)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 1UL );
      mat(2,0) = 5;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of lower column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major sparse matrix subtraction assignment (upper)
   {
      test_ = "Row-major/row-major DiagonalMatrix sparse matrix subtraction assignment (upper)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 1UL );
      mat(0,2) = 5;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of upper row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major sparse matrix subtraction assignment (upper)
   {
      test_ = "Row-major/column-major DiagonalMatrix sparse matrix subtraction assignment (upper)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 1UL );
      mat(0,2) = 5;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of upper column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major sparse matrix subtraction assignment (DiagonalMatrix)
   {
      test_ = "Row-major/row-major DiagonalMatrix sparse matrix subtraction assignment (DiagonalMatrix)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > diag1( 3UL, 2UL );
      diag1(1,1) =  2;
      diag1(2,2) = -2;

      DT diag2( 3UL );
      diag2(0,0) = 1;
      diag2(1,1) = 2;
      diag2(2,2) = 3;

      diag2 -= diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 2UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 0UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 0 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major sparse matrix subtraction assignment (DiagonalMatrix)
   {
      test_ = "Row-major/column-major DiagonalMatrix sparse matrix subtraction assignment (DiagonalMatrix)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > diag1( 3UL, 2UL );
      diag1(1,1) =  2;
      diag1(2,2) = -2;

      DT diag2( 3UL );
      diag2(0,0) = 1;
      diag2(1,1) = 2;
      diag2(2,2) = 3;

      diag2 -= diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 2UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 0UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 0 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix subtraction assignment
   //=====================================================================================

   // Column-major/row-major dense matrix subtraction assignment (diagonal)
   {
      test_ = "Row-major/row-major DiagonalMatrix dense matrix subtraction assignment (diagonal)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(1,1) =  2;
      mat(2,2) = -2;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      diag -= mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix subtraction assignment (diagonal)
   {
      test_ = "Column-major/column-major DiagonalMatrix dense matrix subtraction assignment (diagonal)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(1,1) =  2;
      mat(2,2) = -2;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      diag -= mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major dense matrix subtraction assignment (lower)
   {
      test_ = "Column-major/row-major DiagonalMatrix dense matrix subtraction assignment (lower)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(2,0) = 5;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of lower row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major dense matrix subtraction assignment (lower)
   {
      test_ = "Column-major/column-major DiagonalMatrix dense matrix subtraction assignment (lower)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(2,0) = 5;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of lower column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major dense matrix subtraction assignment (upper)
   {
      test_ = "Column-major/row-major DiagonalMatrix dense matrix subtraction assignment (upper)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(0,2) = 5;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of upper row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major dense matrix subtraction assignment (upper)
   {
      test_ = "Column-major/column-major DiagonalMatrix dense matrix subtraction assignment (upper)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(0,2) = 5;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of upper column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major dense matrix subtraction assignment (DiagonalMatrix)
   {
      test_ = "Column-major/row-major DiagonalMatrix dense matrix subtraction assignment (DiagonalMatrix)";

      DT diag1( 3UL );
      diag1(1,1) =  2;
      diag1(2,2) = -2;

      ODT diag2( 3UL );
      diag2(0,0) = 1;
      diag2(1,1) = 2;
      diag2(2,2) = 3;

      diag2 -= diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 2UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 0UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 0 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix subtraction assignment (DiagonalMatrix)
   {
      test_ = "Column-major/column-major DiagonalMatrix dense matrix subtraction assignment (DiagonalMatrix)";

      ODT diag1( 3UL );
      diag1(1,1) =  2;
      diag1(2,2) = -2;

      ODT diag2( 3UL );
      diag2(0,0) = 1;
      diag2(1,1) = 2;
      diag2(2,2) = 3;

      diag2 -= diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 2UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 0UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 0 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse matrix subtraction assignment
   //=====================================================================================

   // Column-major/row-major sparse matrix subtraction assignment (diagonal)
   {
      test_ = "Column-major/row-major DiagonalMatrix sparse matrix subtraction assignment (diagonal)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 4UL );
      mat(1,1) =  2;
      mat(2,2) = -2;
      mat.insert( 1UL, 2UL, 0 );
      mat.insert( 2UL, 1UL, 0 );

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      diag -= mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major sparse matrix subtraction assignment (diagonal)
   {
      test_ = "Column-major/column-major DiagonalMatrix sparse matrix subtraction assignment (diagonal)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 4UL );
      mat(1,1) =  2;
      mat(2,2) = -2;
      mat.insert( 1UL, 2UL, 0 );
      mat.insert( 2UL, 1UL, 0 );

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      diag -= mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major sparse matrix subtraction assignment (lower)
   {
      test_ = "Column-major/row-major DiagonalMatrix sparse matrix subtraction assignment (lower)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 1UL );
      mat(2,0) = 5;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of lower row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major sparse matrix subtraction assignment (lower)
   {
      test_ = "Column-major/column-major DiagonalMatrix sparse matrix subtraction assignment (lower)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 1UL );
      mat(2,0) = 5;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of lower column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major sparse matrix subtraction assignment (upper)
   {
      test_ = "Column-major/row-major DiagonalMatrix sparse matrix subtraction assignment (upper)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 1UL );
      mat(0,2) = 5;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of upper row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major sparse matrix subtraction assignment (upper)
   {
      test_ = "Column-major/column-major DiagonalMatrix sparse matrix subtraction assignment (upper)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 1UL );
      mat(0,2) = 5;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of upper column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major sparse matrix subtraction assignment (DiagonalMatrix)
   {
      test_ = "Column-major/row-major DiagonalMatrix sparse matrix subtraction assignment (DiagonalMatrix)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > diag1( 3UL, 2UL );
      diag1(1,1) =  2;
      diag1(2,2) = -2;

      ODT diag2( 3UL );
      diag2(0,0) = 1;
      diag2(1,1) = 2;
      diag2(2,2) = 3;

      diag2 -= diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 2UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 0UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 0 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major sparse matrix subtraction assignment (DiagonalMatrix)
   {
      test_ = "Column-major/column-major DiagonalMatrix sparse matrix subtraction assignment (DiagonalMatrix)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > diag1( 3UL, 2UL );
      diag1(1,1) =  2;
      diag1(2,2) = -2;

      ODT diag2( 3UL );
      diag2(0,0) = 1;
      diag2(1,1) = 2;
      diag2(2,2) = 3;

      diag2 -= diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 2UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 0UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 0 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the DiagonalMatrix Schur product assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the Schur product assignment operators of the DiagonalMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testSchurAssign()
{
   //=====================================================================================
   // Row-major dense matrix Schur product assignment
   //=====================================================================================

   // Row-major/row-major dense matrix Schur product assignment (general)
   {
      test_ = "Row-major/row-major DiagonalMatrix dense matrix Schur product assignment (general)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat{ { 1, 0, 9 }, { 0, 0, 0 }, { 9, 0, 3 } };

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      diag %= mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix Schur product assignment (general)
   {
      test_ = "Row-major/column-major DiagonalMatrix dense matrix Schur product assignment (general)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat{ { 1, 9, 9 }, { 9, 0, 9 }, { 9, 9, 3 } };

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      diag %= mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major dense matrix Schur product assignment (DiagonalMatrix)
   {
      test_ = "Row-major/row-major DiagonalMatrix dense matrix Schur product assignment (DiagonalMatrix)";

      DT diag1( 3UL );
      diag1(0,0) = 1;
      diag1(2,2) = 3;

      DT diag2( 3UL );
      diag2(0,0) = 1;
      diag2(1,1) = 2;
      diag2(2,2) = 3;

      diag2 %= diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 2UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 0UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 0 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix Schur product assignment (DiagonalMatrix)
   {
      test_ = "Row-major/column-major DiagonalMatrix dense matrix Schur product assignment (DiagonalMatrix)";

      ODT diag1( 3UL );
      diag1(0,0) = 1;
      diag1(2,2) = 3;

      DT diag2( 3UL );
      diag2(0,0) = 1;
      diag2(1,1) = 2;
      diag2(2,2) = 3;

      diag2 %= diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 2UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 0UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 0 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix Schur product assignment
   //=====================================================================================

   // Row-major/row-major sparse matrix Schur product assignment (general)
   {
      test_ = "Row-major/row-major DiagonalMatrix sparse matrix Schur product assignment (general)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 6UL );
      mat(0,0) = 1;
      mat(0,2) = 9;
      mat(2,0) = 9;
      mat(2,2) = 3;
      mat.insert( 1UL, 2UL, 0 );
      mat.insert( 2UL, 1UL, 0 );

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      diag %= mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major sparse matrix Schur product assignment (general)
   {
      test_ = "Row-major/column-major DiagonalMatrix sparse matrix Schur product assignment (general)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 4UL );
      mat(0,0) = 1;
      mat(0,2) = 9;
      mat(2,0) = 9;
      mat(2,2) = 3;
      mat.insert( 1UL, 2UL, 0 );
      mat.insert( 2UL, 1UL, 0 );

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      diag %= mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major sparse matrix Schur product assignment (DiagonalMatrix)
   {
      test_ = "Row-major/row-major DiagonalMatrix sparse matrix Schur product assignment (DiagonalMatrix)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > diag1( 3UL, 2UL );
      diag1(0,0) = 1;
      diag1(2,2) = 3;

      DT diag2( 3UL );
      diag2(0,0) = 1;
      diag2(1,1) = 2;
      diag2(2,2) = 3;

      diag2 %= diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 2UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 0UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 0 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major sparse matrix Schur product assignment (DiagonalMatrix)
   {
      test_ = "Row-major/column-major DiagonalMatrix sparse matrix Schur product assignment (DiagonalMatrix)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > diag1( 3UL, 2UL );
      diag1(0,0) = 1;
      diag1(2,2) = 3;

      DT diag2( 3UL );
      diag2(0,0) = 1;
      diag2(1,1) = 2;
      diag2(2,2) = 3;

      diag2 %= diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 2UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 0UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 0 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix Schur product assignment
   //=====================================================================================

   // Column-major/row-major dense matrix Schur product assignment (general)
   {
      test_ = "Column-major/row-major DiagonalMatrix dense matrix Schur product assignment (general)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat{ { 1, 0, 9 }, { 0, 0, 0 }, { 9, 0, 3 } };

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      diag %= mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix Schur product assignment (general)
   {
      test_ = "Column-major/column-major DiagonalMatrix dense matrix Schur product assignment (general)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat{ { 1, 0, 9 }, { 0, 0, 0 }, { 9, 0, 3 } };

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      diag %= mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major dense matrix Schur product assignment (DiagonalMatrix)
   {
      test_ = "Column-major/row-major DiagonalMatrix dense matrix Schur product assignment (DiagonalMatrix)";

      DT diag1( 3UL );
      diag1(0,0) = 1;
      diag1(2,2) = 3;

      ODT diag2( 3UL );
      diag2(0,0) = 1;
      diag2(1,1) = 2;
      diag2(2,2) = 3;

      diag2 %= diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 2UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 0UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 0 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix Schur product assignment (DiagonalMatrix)
   {
      test_ = "Column-major/column-major DiagonalMatrix dense matrix Schur product assignment (DiagonalMatrix)";

      ODT diag1( 3UL );
      diag1(0,0) = 1;
      diag1(2,2) = 3;

      ODT diag2( 3UL );
      diag2(0,0) = 1;
      diag2(1,1) = 2;
      diag2(2,2) = 3;

      diag2 %= diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 2UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 0UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 0 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix Schur product assignment
   //=====================================================================================

   // Column-major/row-major sparse matrix Schur product assignment (general)
   {
      test_ = "Column-major/row-major DiagonalMatrix sparse matrix Schur product assignment (general)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 6UL );
      mat(0,0) = 1;
      mat(0,2) = 9;
      mat(2,0) = 9;
      mat(2,2) = 3;
      mat.insert( 1UL, 2UL, 0 );
      mat.insert( 2UL, 1UL, 0 );

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      diag %= mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major sparse matrix Schur product assignment (general)
   {
      test_ = "Column-major/column-major DiagonalMatrix sparse matrix Schur product assignment (general)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 6UL );
      mat(0,0) = 1;
      mat(0,2) = 9;
      mat(2,0) = 9;
      mat(2,2) = 3;
      mat.insert( 1UL, 2UL, 0 );
      mat.insert( 2UL, 1UL, 0 );

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      diag %= mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 2UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 0UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 0 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major sparse matrix Schur product assignment (DiagonalMatrix)
   {
      test_ = "Column-major/row-major DiagonalMatrix sparse matrix Schur product assignment (DiagonalMatrix)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > diag1( 3UL, 2UL );
      diag1(0,0) = 1;
      diag1(2,2) = 3;

      ODT diag2( 3UL );
      diag2(0,0) = 1;
      diag2(1,1) = 2;
      diag2(2,2) = 3;

      diag2 %= diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 2UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 0UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 0 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major sparse matrix Schur product assignment (DiagonalMatrix)
   {
      test_ = "Column-major/column-major DiagonalMatrix sparse matrix Schur product assignment (DiagonalMatrix)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > diag1( 3UL, 2UL );
      diag1(0,0) = 1;
      diag1(2,2) = 3;

      ODT diag2( 3UL );
      diag2(0,0) = 1;
      diag2(1,1) = 2;
      diag2(2,2) = 3;

      diag2 %= diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 2UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 0UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 1 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 0 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 0 0 )\n( 0 0 9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the DiagonalMatrix multiplication assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the multiplication assignment operators of the DiagonalMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testMultAssign()
{
   //=====================================================================================
   // Row-major dense matrix multiplication assignment
   //=====================================================================================

   // Row-major/row-major dense matrix multiplication assignment (diagonal)
   {
      test_ = "Row-major/row-major DiagonalMatrix dense matrix multiplication assignment (diagonal)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(0,0) = 2;
      mat(1,1) = 2;
      mat(2,2) = 2;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      diag *= mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 2 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 4 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 2 0 0 )\n( 0 4 0 )\n( 0 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix multiplication assignment (diagonal)
   {
      test_ = "Row-major/column-major DiagonalMatrix dense matrix multiplication assignment (diagonal)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(0,0) = 2;
      mat(1,1) = 2;
      mat(2,2) = 2;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      diag *= mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 2 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 4 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 2 0 0 )\n( 0 4 0 )\n( 0 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major dense matrix multiplication assignment (lower)
   {
      test_ = "Row-major/row-major DiagonalMatrix dense matrix multiplication assignment (lower)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(2,0) = 5;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of lower row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major dense matrix multiplication assignment (lower)
   {
      test_ = "Row-major/column-major DiagonalMatrix dense matrix multiplication assignment (lower)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(2,0) = 5;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of lower column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major dense matrix multiplication assignment (upper)
   {
      test_ = "Row-major/row-major DiagonalMatrix dense matrix multiplication assignment (upper)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(0,2) = 5;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of upper row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major dense matrix multiplication assignment (upper)
   {
      test_ = "Row-major/column-major DiagonalMatrix dense matrix multiplication assignment (upper)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(0,2) = 5;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of upper column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major dense matrix multiplication assignment (DiagonalMatrix)
   {
      test_ = "Row-major/row-major DiagonalMatrix dense matrix multiplication assignment (DiagonalMatrix)";

      DT diag1( 3UL );
      diag1(0,0) = 2;
      diag1(1,1) = 2;
      diag1(2,2) = 2;

      DT diag2( 3UL );
      diag2(0,0) = 1;
      diag2(1,1) = 2;
      diag2(2,2) = 3;

      diag2 *= diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 3UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 1UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 2 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 4 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 2 0 0 )\n( 0 4 0 )\n( 0 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix multiplication assignment (DiagonalMatrix)
   {
      test_ = "Row-major/column-major DiagonalMatrix dense matrix multiplication assignment (DiagonalMatrix)";

      ODT diag1( 3UL );
      diag1(0,0) = 2;
      diag1(1,1) = 2;
      diag1(2,2) = 2;

      DT diag2( 3UL );
      diag2(0,0) = 1;
      diag2(1,1) = 2;
      diag2(2,2) = 3;

      diag2 *= diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 3UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 1UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 2 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 4 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 2 0 0 )\n( 0 4 0 )\n( 0 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix multiplication assignment
   //=====================================================================================

   // Row-major/row-major sparse matrix multiplication assignment (diagonal)
   {
      test_ = "Row-major/row-major DiagonalMatrix sparse matrix multiplication assignment (diagonal)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 5UL );
      mat(0,0) = 2;
      mat(1,1) = 2;
      mat(2,2) = 2;
      mat.insert( 1UL, 2UL, 0 );
      mat.insert( 2UL, 1UL, 0 );

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      diag *= mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 2 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 4 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 2 0 0 )\n( 0 4 0 )\n( 0 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major sparse matrix multiplication assignment (diagonal)
   {
      test_ = "Row-major/column-major DiagonalMatrix sparse matrix multiplication assignment (diagonal)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 5UL );
      mat(0,0) = 2;
      mat(1,1) = 2;
      mat(2,2) = 2;
      mat.insert( 1UL, 2UL, 0 );
      mat.insert( 2UL, 1UL, 0 );

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      diag *= mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 2 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 4 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 2 0 0 )\n( 0 4 0 )\n( 0 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major sparse matrix multiplication assignment (lower)
   {
      test_ = "Row-major/row-major DiagonalMatrix sparse matrix multiplication assignment (lower)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 1UL );
      mat(2,0) = 5;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of lower row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major sparse matrix multiplication assignment (lower)
   {
      test_ = "Row-major/column-major DiagonalMatrix sparse matrix multiplication assignment (lower)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 1UL );
      mat(2,0) = 5;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of lower column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major sparse matrix multiplication assignment (upper)
   {
      test_ = "Row-major/row-major DiagonalMatrix sparse matrix multiplication assignment (upper)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 1UL );
      mat(0,2) = 5;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of upper row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major sparse matrix multiplication assignment (upper)
   {
      test_ = "Row-major/column-major DiagonalMatrix sparse matrix multiplication assignment (upper)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 1UL );
      mat(0,2) = 5;

      DT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of upper column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major sparse matrix multiplication assignment (DiagonalMatrix)
   {
      test_ = "Row-major/row-major DiagonalMatrix sparse matrix multiplication assignment (DiagonalMatrix)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > diag1( 3UL, 3UL );
      diag1(0,0) = 2;
      diag1(1,1) = 2;
      diag1(2,2) = 2;

      DT diag2( 3UL );
      diag2(0,0) = 1;
      diag2(1,1) = 2;
      diag2(2,2) = 3;

      diag2 *= diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 3UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 1UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 2 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 4 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 2 0 0 )\n( 0 4 0 )\n( 0 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major sparse matrix multiplication assignment (DiagonalMatrix)
   {
      test_ = "Row-major/column-major DiagonalMatrix sparse matrix multiplication assignment (DiagonalMatrix)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > diag1( 3UL, 3UL );
      diag1(0,0) = 2;
      diag1(1,1) = 2;
      diag1(2,2) = 2;

      DT diag2( 3UL );
      diag2(0,0) = 1;
      diag2(1,1) = 2;
      diag2(2,2) = 3;

      diag2 *= diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 3UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 1UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 2 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 4 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 2 0 0 )\n( 0 4 0 )\n( 0 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix multiplication assignment
   //=====================================================================================

   // Column-major/row-major dense matrix multiplication assignment (diagonal)
   {
      test_ = "Column-major/row-major DiagonalMatrix dense matrix multiplication assignment (diagonal)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(0,0) = 2;
      mat(1,1) = 2;
      mat(2,2) = 2;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      diag *= mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 2 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 4 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 2 0 0 )\n( 0 4 0 )\n( 0 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix multiplication assignment (diagonal)
   {
      test_ = "Column-major/column-major DiagonalMatrix dense matrix multiplication assignment (diagonal)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(0,0) = 2;
      mat(1,1) = 2;
      mat(2,2) = 2;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      diag *= mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 2 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 4 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 2 0 0 )\n( 0 4 0 )\n( 0 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major dense matrix multiplication assignment (lower)
   {
      test_ = "Column-major/row-major DiagonalMatrix dense matrix multiplication assignment (lower)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(2,0) = 5;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of lower row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major dense matrix multiplication assignment (lower)
   {
      test_ = "Column-major/column-major DiagonalMatrix dense matrix multiplication assignment (lower)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(2,0) = 5;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of lower column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major dense matrix multiplication assignment (upper)
   {
      test_ = "Column";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(0,2) = 5;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of upper row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major dense matrix multiplication assignment (upper)
   {
      test_ = "Column-major/column-major DiagonalMatrix dense matrix multiplication assignment (upper)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(0,2) = 5;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of upper column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major dense matrix multiplication assignment (DiagonalMatrix)
   {
      test_ = "Column-major/row-major DiagonalMatrix dense matrix multiplication assignment (DiagonalMatrix)";

      ODT diag1( 3UL );
      diag1(0,0) = 2;
      diag1(1,1) = 2;
      diag1(2,2) = 2;

      DT diag2( 3UL );
      diag2(0,0) = 1;
      diag2(1,1) = 2;
      diag2(2,2) = 3;

      diag2 *= diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 3UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 1UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 2 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 4 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 2 0 0 )\n( 0 4 0 )\n( 0 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix multiplication assignment (DiagonalMatrix)
   {
      test_ = "Column-major/column-major DiagonalMatrix dense matrix multiplication assignment (DiagonalMatrix)";

      ODT diag1( 3UL );
      diag1(0,0) = 2;
      diag1(1,1) = 2;
      diag1(2,2) = 2;

      ODT diag2( 3UL );
      diag2(0,0) = 1;
      diag2(1,1) = 2;
      diag2(2,2) = 3;

      diag2 *= diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 3UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 1UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 2 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 4 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 2 0 0 )\n( 0 4 0 )\n( 0 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse matrix multiplication assignment
   //=====================================================================================

   // Column-major/row-major sparse matrix multiplication assignment (diagonal)
   {
      test_ = "Column-major/row-major DiagonalMatrix sparse matrix multiplication assignment (diagonal)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 5UL );
      mat(0,0) = 2;
      mat(1,1) = 2;
      mat(2,2) = 2;
      mat.insert( 1UL, 2UL, 0 );
      mat.insert( 2UL, 1UL, 0 );

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      diag *= mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 2 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 4 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 2 0 0 )\n( 0 4 0 )\n( 0 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major sparse matrix multiplication assignment (diagonal)
   {
      test_ = "Column-major/column-major DiagonalMatrix sparse matrix multiplication assignment (diagonal)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 5UL );
      mat(0,0) = 2;
      mat(1,1) = 2;
      mat(2,2) = 2;
      mat.insert( 1UL, 2UL, 0 );
      mat.insert( 2UL, 1UL, 0 );

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      diag *= mat;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkCapacity( diag, 9UL );
      checkNonZeros( diag, 3UL );
      checkNonZeros( diag, 0UL, 1UL );
      checkNonZeros( diag, 1UL, 1UL );
      checkNonZeros( diag, 2UL, 1UL );

      if( diag(0,0) != 2 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 4 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 2 0 0 )\n( 0 4 0 )\n( 0 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major sparse matrix multiplication assignment (lower)
   {
      test_ = "Column-major/row-major DiagonalMatrix sparse matrix multiplication assignment (lower)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 1UL );
      mat(2,0) = 5;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of lower row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major sparse matrix multiplication assignment (lower)
   {
      test_ = "Column-major/column-major DiagonalMatrix sparse matrix multiplication assignment (lower)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 1UL );
      mat(2,0) = 5;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of lower column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major sparse matrix multiplication assignment (upper)
   {
      test_ = "Column-major/row-major DiagonalMatrix sparse matrix multiplication assignment (upper)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 1UL );
      mat(0,2) = 5;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of upper row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major sparse matrix multiplication assignment (upper)
   {
      test_ = "Column-major/column-major DiagonalMatrix sparse matrix multiplication assignment (upper)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 1UL );
      mat(0,2) = 5;

      ODT diag( 3UL );
      diag(0,0) = 1;
      diag(1,1) = 2;
      diag(2,2) = 3;

      try {
         diag *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of upper column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major sparse matrix multiplication assignment (DiagonalMatrix)
   {
      test_ = "Column-major/row-major DiagonalMatrix sparse matrix multiplication assignment (DiagonalMatrix)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > diag1( 3UL, 3UL );
      diag1(0,0) = 2;
      diag1(1,1) = 2;
      diag1(2,2) = 2;

      ODT diag2( 3UL );
      diag2(0,0) = 1;
      diag2(1,1) = 2;
      diag2(2,2) = 3;

      diag2 *= diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 3UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 1UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 2 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 4 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 2 0 0 )\n( 0 4 0 )\n( 0 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major sparse matrix multiplication assignment (DiagonalMatrix)
   {
      test_ = "Column-major/column-major DiagonalMatrix sparse matrix multiplication assignment (DiagonalMatrix)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > diag1( 3UL, 3UL );
      diag1(0,0) = 2;
      diag1(1,1) = 2;
      diag1(2,2) = 2;

      ODT diag2( 3UL );
      diag2(0,0) = 1;
      diag2(1,1) = 2;
      diag2(2,2) = 3;

      diag2 *= diag1;

      checkRows    ( diag2, 3UL );
      checkColumns ( diag2, 3UL );
      checkCapacity( diag2, 9UL );
      checkNonZeros( diag2, 3UL );
      checkNonZeros( diag2, 0UL, 1UL );
      checkNonZeros( diag2, 1UL, 1UL );
      checkNonZeros( diag2, 2UL, 1UL );

      if( diag2(0,0) != 2 || diag2(0,1) != 0 || diag2(0,2) != 0 ||
          diag2(1,0) != 0 || diag2(1,1) != 4 || diag2(1,2) != 0 ||
          diag2(2,0) != 0 || diag2(2,1) != 0 || diag2(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << diag2 << "\n"
             << "   Expected result:\n( 2 0 0 )\n( 0 4 0 )\n( 0 0 6 )\n";
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
   std::cout << "   Running DiagonalMatrix dense test (part 1)..." << std::endl;

   try
   {
      RUN_DIAGONALMATRIX_DENSE_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during DiagonalMatrix dense test (part 1):\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
